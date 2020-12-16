function c2 = threeComponentProcessing(c2, cstation, n, J, K, dt, width,...
    adjustAlignment, gainControl)

% threeComponentProcessing: create three-component objects for each of the
%   records in the network-correlation object, applies rotation and a
%   polarization filter to enhance body wave arrivals.
% Input:
%   cstation : cell array with the names of the chosen stations
%   n, J, K, dt, width: are the polarization-filter parameters described by
%       Montabetti & Kanasewich, 1970. Enhancement of teleseismic body wav
%       phases with a polarisation analysis, Geophys. J. R. Astron. Soc.
%   adjustAlignment: boolean, whether to adjust seismogram alignment again
%   gainControl: boolean, whether to apply automatic gain control
% Output:
%   c2 : network correlation object with substructure: Channel, cstation,
%       correlation/catalog object. 

% Can be easily parallelized across events, by setting event loop to
% "parfor", but function is still not the quickest for larger datasets.
% c2 is a structure with the 

SECPERDAY = 60*60*24;
DAYPERSEC = 1/SECPERDAY;

c2.Zp = c2.Z;
c2.R = c2.Z;
c2.Rp = c2.Z;
c2.T = c2.Z;
c2.Tp = c2.Z;
c2.TC = c2.Z;


%for each station
for k=1:1:length(cstation)
    c2.TC.(cstation{k}) = repmat(threecomp,...
        size(get(c2.Z.(cstation{k}).corr,'waveforms')));
    %check that all events have three components at station
    if c2.Z.(cstation{k}).cat.numberOfEvents ~=...
            c2.E.(cstation{k}).cat.numberOfEvents ||...
            c2.Z.(cstation{k}).cat.numberOfEvents ~=...
            c2.N.(cstation{k}).cat.numberOfEvents
        error('Error: missing waveform on one or more components')
    end
        
    Z = get(c2.Z.(cstation{k}).corr,'waveforms');
    N = get(c2.N.(cstation{k}).corr,'waveforms');
    E = get(c2.E.(cstation{k}).corr,'waveforms');
    Ztrigger = get(c2.Z.(cstation{k}).corr,'trig');
    nw = length(Z);
    
    % check that all waveforms have same length across components
    for j=1:1:nw
        lZ = length(get(Z(j),'data'));
        lN = length(get(N(j),'data'));
        lE = length(get(E(j),'data'));
        if lZ ~= lN || lZ ~= lE
            warning(['Components for station ', cstation{k}, ' of event ',...
                datestr(c2.Z.(cstation{k}).cat.otime(j),'yyyymmddHHMM'),...
                ' are not the same length, cutting all three to the ',...
                'shortest waveform.'])
            [shortest, shortestI] = min([lZ, lN, lE]);
            if shortestI == 1
                starttime = get(Z(j),'start');
                endtime = get(Z(j),'end');
            elseif shortestI == 2
                starttime = get(N(j),'start');
                endtime = get(N(j),'end');
            elseif shortestI == 3
                starttime = get(E(j),'start');
                endtime = get(E(j),'end');
            end
            Z(j) = extract(Z(j),'TIME', starttime, endtime);
            Z(j) = set(Z(j),'start',starttime);
            N(j) = extract(N(j),'TIME', starttime, endtime);
            N(j) = set(N(j),'start',starttime);
            E(j) = extract(E(j),'TIME', starttime, endtime);
            E(j) = set(E(j),'start',starttime);
        end     
    end
    

    % BackAzimuth from great circle arc between station and
    % hypoDD-location
    evcat = c2.Z.(cstation{k}).cat;
    backazimuth = backazmiuthStationToEvent(cstation{k}, evcat);
        
    R = Z;
    T = Z;
    Zp = Z;
    Rp = Z;
    Tp = Z;
    
    % create dummy array of three-comp objects
    zrttc = threecomp([Z(1),N(1),E(1)],0,0,[0 0 0 90 90 90]);
    zrttc = repmat(zrttc,numel(Z),1);
    % for each event
    for j=1:1:nw
        % Check if any channel is NaN or all zero, then don't rotate
        if any(isnan(get(Z(j),'data'))) || any(isnan(get(N(j),'data'))) ||...
                any(isnan(get(E(j),'data'))) ||...
                sum(get(Z(j),'data'))==0 || sum(get(N(j),'data'))==0 ||...
                sum(get(E(j),'data'))==0
            warning(['One component contains NaN or all zero, cannot ',...
                'rotate (setting R and T to NaN'])
            % fill all non-rotated/filtered components with zeros
            Zp(j) = Z(j);
            Zp(j) = set(Zp(j),'data',zeros(size((get(Zp(j),'data')))));
            
            R(j) = Z(j);
            chn = get(R(j),'ChannelTag');
            chn.channel = [chn.channel(1:end-1), 'R'];
            R(j) = set(R(j),'ChannelTag',chn);
            R(j) = set(R(j),'data',zeros(size((get(R(j),'data')))));           
            
            T(j) = Z(j);
            chn = get(T(j),'ChannelTag');
            chn.channel = [chn.channel(1:end-1), 'T'];
            T(j) = set(T(j),'ChannelTag',chn);
            T(j) = set(T(j),'data',zeros(size((get(T(j),'data')))));
            
            Rp(j) = set(R(j),'data',zeros(size((get(R(j),'data')))));
            Tp(j) = set(T(j),'data',zeros(size((get(T(j),'data')))));
        else
            W = [Z(j), N(j), E(j)];
            % Check for sub-sample time differences of channels
            st = get(W,'start');
            if st(1) ~= st(2) || st(1) ~= st(3)
                % Check that the difference is less than 10 % of one sample
                if st(1) - st(2) <= 0.1*1/get(W(1)*DAYPERSEC,'freq') &&...
                        st(1) - st(3) <= 0.1*1/get(W(1)*DAYPERSEC,'freq')
                    starttim = get(W(1),'start');
                    W(2) = set(W(2), 'start', starttim);
                    W(3) = set(W(3), 'start', starttim);
                end
            end 
            
            % now fill the required information on backazimuth, trigger
            % etc.
            trigger = [Ztrigger(j), Ztrigger(j), Ztrigger(j)]; 
                      
            orientation = [0 0 0 90 90 90];
            baz = backazimuth(j);
            tc = threecomp(W, trigger, baz, orientation);
            % set rotated records
            zrt = rotate(tc);

            % dt here as the sampling frequency. and width as twice the 
            % longest period
            % (or maybe better once to twice the DOMINANT period? -->)
            % looks like 1.5 times the longest period is a good compromise!
            ww = get(zrt,'waveform');
            R(j) = ww(2);
            T(j) = ww(3);
            
            % set polarization-filtered records
            zrt = polarizationfilter(zrt, n, J, K, dt, width);
            ww = get(zrt,'waveform');
            Zp(j) = ww(1);
            Rp(j) = ww(2);
            Tp(j) = ww(3);
        end
    end
    
    if gainControl
        comp = fieldnames(c2);
        for nc = 1:1:length(comps)-1
            c2.(comp{nc}).(cstation{k}).corr = agc(c2.(comp{nc}).(...
                cstation{k}).corr,1);
        end
    end
    
    c2.R.(cstation{k}).corr = set(c2.R.(cstation{k}).corr,'waveforms',R);
    c2.T.(cstation{k}).corr = set(c2.T.(cstation{k}).corr,'waveforms',T);
    
    c2.Zp.(cstation{k}).corr = set(c2.Zp.(cstation{k}).corr,'waveforms',Zp);
    c2.Rp.(cstation{k}).corr = set(c2.Rp.(cstation{k}).corr,'waveforms',Rp);
    c2.Tp.(cstation{k}).corr = set(c2.Tp.(cstation{k}).corr,'waveforms',Tp);
    
       
    % Fill NaN-value gaps % don't think this is necessary.. let's see
    if adjustAlignment
        % Adjust triggers according to the polarizationfiltered-vertical
        % Component.
        % 1. run xcorr on first ~1.5 s of Zp-component
        c2.Zp.(cstation{k}).corr = xcorr(c2.Zp.(cstation{k}).corr,...
            [-width/3 width]);
        if length(get(c2.Zp.(cstation{k}).corr,'trig'))>1
            c2.Zp.(cstation{k}).corr = linkage(c2.Zp.(cstation{k}).corr);
        % 2. Correct all components with the same adjustments from Zp-
        %    component
            Zplags = get(c2.Zp.(cstation{k}).corr,'lag');
            comp = fieldnames(c2);
            for j=1:1:numel(comp)
                % don't try this on the three-component object
                if ~strcmp(comp{j},'TC')
                    c2.(comp{j}).(cstation{k}).corr = set(...
                        c2.(comp{j}).(cstation{k}).corr,'lag', Zplags);
                    c2.(comp{j}).(cstation{k}).corr = adjusttrig(...
                        c2.(comp{j}).(cstation{k}).corr, 'MMM', width/3);
                    %0.06
            % 3. run xcorr again to update xcorr/lag fields etc.
                    c2.(comp{j}).(cstation{k}).corr = xcorr(...
                        c2.(comp{j}).(cstation{k}).corr, [-0.2 10]);
                    c2.(comp{j}).(cstation{k}).corr = linkage(...
                        c2.(comp{j}).(cstation{k}).corr);
                end
            end
        end
    end
    
    disp(['Three-component processing: ', num2str(...
        k/length(cstation)*100), ' % completed']);
end

% Sample code to plot waveforms and rectilinearity/weighting functions
% www = [get(zrt,'waveform'), get(zrt,'rectilinearity')];
% fig = figure();
% for j=1:1:length(www)
%     ax(j) = subplot(length(www),1,j);
%     plot(www(j),'axeshandle',ax(j))
% end
