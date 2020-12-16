
function [correlationCatalog] = buildCorrelationCatalog(seisCatalog,...
    eventI, cstation, corcomp, targetSamplingRate)
    %function to build a complete correlatable catalog for many stations, many
    %channels, and many events. those can then be individually accessed within
    %that structure through a "." (dot)-notation.
    %seisCatalog is a catalog object of earthquakes with waveforms 
    % cstation is a cell array with the names of the desired
    % stations
    % corcomp is a cell array with the names of the desired
    % channels
    % eventI are the indices of the events to be processed
   
    %don't know what I wanted this paramter for
    deleteArrivals = false;

	secPerDay = 60*60*24;
    
    %for each channel
    for ncc=1:1:numel(corcomp)
        c.(corcomp{ncc})=struct();
    end
    
    if size(eventI,1) > size(eventI,2)
        eventI = eventI';
    end
    
    %for each event
    for j = eventI
        w = seisCatalog.waveforms{j};
        %check if there is no waveform at all for that event
        if isempty(w)
            continue
        end
        wstation = get(w,'station');
        wchannel = get(w,'channel');
        % I don't understand why this happens, but sometimes the above two
        % get statements return a X*X cell of station names for a
        % waveform object with only X waveforms
        if size(wstation,1) == size(wstation,2) && size(wstation,1) > 1
            wstation = wstation(:,1);
        end
        if size(wchannel,1) == size(wchannel,2) && size(wchannel,1) > 1
            wchannel = wchannel(:,1);
        end
        wchannel = char(wchannel);
        wcomponent = cellstr(wchannel(:,end));
        %for each station
        for k=1:1:length(cstation)
            %find the trigger
            trigSta = {seisCatalog.arrivals{j}.stacode};
            ptriggers = {seisCatalog.arrivals{j}.p_time};
            striggers = {seisCatalog.arrivals{j}.s_time};
            trigStaI = find(strcmp(trigSta,cstation{k}));
            %set the trigger to the P-pick (or maybe an assumed value)
            if ~isempty(trigStaI) && ~isempty(ptriggers{trigStaI})
%                     if there are multiple picks, take the earlier one
                if numel(trigStaI) > 1
                    [earliest, earliestI] = min(ptriggers{trigStaI});
                    trigStaI = trigStaI(earliestI);
                end
                PTrig = ptriggers{trigStaI};
            else
                % For now: if P-trigger is not known, then skip this record
                PTrig = NaN;
                continue
            end

            %in case i may need the s-pick at some time
            if ~isempty(trigStaI) && ~isempty(striggers{trigStaI})
                setSTrig = striggers{trigStaI};
                % if there is an S-pick but no P-pick, then set
                % the P-pick to an assumed value:
%                             if isnan(setPTrig)
%                                 %calculate P-arrival from S-arrival and S/P
% %                                 seisCatalog.arrivals{j}(trigStaI).s_dis
%                                 setPTrig = NaN;
%                             end
            else
                setSTrig = NaN;
            end
            
            %for each channel
            for ncc=1:1:numel(corcomp)

                clear prevWav prevTrig;
                if ~isfield(c.(corcomp{ncc}),cstation{k})
                    c.(corcomp{ncc}).(cstation{k}).corr = correlation();
                    c.(corcomp{ncc}).(cstation{k}).cat = seisCatalog.subset([]);
                else
                    prevWav = get(c.(corcomp{ncc}).(cstation{k}).corr,'waveform');
                    prevTrig = get(c.(corcomp{ncc}).(cstation{k}).corr,'trig');
                end
                
                if exist('prevTrig','var') && ~isnan(PTrig)
                    setPTrig = [prevTrig; PTrig];
                else
                    setPTrig = PTrig;
                end
                
                if ~isnan(PTrig)
                    %find the waveform
                    staWavI = find(strcmp(wstation,cstation{k}));
                    staCompWavI = find(strcmp(wstation,cstation{k}) &...
                        strcmp(wcomponent,corcomp{ncc}));
                    if ~isempty(staCompWavI)
                        % if there are duplicate waveforms, then select the
                        % the longer, or the more densely sampled one, for 
                        % processing
                        if numel(staCompWavI) > 1
                            wavedata = get(seisCatalog.waveforms{j}(...
                                staCompWavI),'data');
                            wavefreq = get(seisCatalog.waveforms{j}(...
                                staCompWavI),'freq');
                            waveNsamples = cellfun(@numel,wavedata);
                            waveDuration = waveNsamples./wavefreq;
                            [longestwave, longI] = max(waveDuration);
                            staCompWavI = staCompWavI(longI);
                        end

                        testWav = seisCatalog.waveforms{j}(staCompWavI);
    %                     testWav = testWav(1);
                        testfreq = get(testWav,'freq');
                        testdata = get(testWav,'data');
                        if testfreq > 9.99 
                            % & sum(testdata)~=0
                            %resample waveform to target Sampling rate
                            crunchFactor = get(testWav,'freq') / targetSamplingRate;
                            crunchFactor = round(crunchFactor,4);
                            [Q, P] = numden(sym(crunchFactor));
                            Q = double(Q);
                            P = double(P);

                            D = double(testWav);
                            ResampleD = resample(D,P,Q); 
                            % put back into waveform, but don't forget to 
                            % update the frequency
                            testWav = set(testWav,'data',ResampleD,...
                                'Freq', targetSamplingRate);
                            % and for good measure... update the waveform's history
%                               testWav = testWav.addHistory('Resampled data outside of waveform object');
                        else
                            % Set a NAN-waveform
                            testWav = set(testWav,'data',NaN(...
                                targetSamplingRate/testfreq * length(get(testWav,'data')),1));
                            testWav = set(testWav,'freq',targetSamplingRate);
                        end
                    elseif ~isempty(staWavI) 
                        % if there are waveforms for the other channels of 
                        % that station, then set a dummy NaN-waveform
                        % set a NaN-waveform
                        testWav = waveform();
                        ct = get(seisCatalog.waveforms{j}(staWavI(1)),'ChannelTag');
                        ct = ct(1);
                        ct.channel = [ct.channel(1:end-1), corcomp{ncc}];
                        testWav = set(testWav,'ChannelTag', ct);
                        testWav = set(testWav,'data', NaN(120*targetSamplingRate,1));
                        testWav = set(testWav,'freq',targetSamplingRate);
                        testWav = set(testWav,'start',...
                            get(seisCatalog.waveforms{j}(staWavI(1)),'start'));
                    else
                        continue
                    end
                    
                    % Cut waveform to the time from the trigger and some
                    testWav = extract(testWav, 'TIME', PTrig-15/secPerDay,...
                        PTrig + 30/secPerDay);

                    %set the waveform for the correlation event
                    setWav = seisCatalog.waveforms{j}(staCompWavI);
                    if exist('prevWav', 'var')
%                                 setWav = [prevWav;setWav];
                        setWav = [prevWav;testWav];
                    else
                        setWav = testWav;
                    end
                    %set the event details for each waveform within the
                    %correlation event.  THIS IS SLOOOOOOWWWWW. 
                    c.(corcomp{ncc}).(cstation{k}).cat.table = ...
                        [c.(corcomp{ncc}).(cstation{k}).cat.table; ...
                        seisCatalog.subset(j).table(:,1:end)];
                    %delete arrivals from new per-station-event-catalog
                    c.(corcomp{ncc}).(cstation{k}).cat.table.arrivals{end,1}(:)=[];
                    % sort the arrivals for that very station for this event to
                    % the event-at-station-catalog
                    c.(corcomp{ncc}).(cstation{k}).cat.table.arrivals{end,1} =...
                        seisCatalog.table.arrivals{j}(trigStaI);

                    %set triggers, waveforms, and catalog
                    c.(corcomp{ncc}).(cstation{k}).corr =...
                        set(c.(corcomp{ncc}).(cstation{k}).corr,'trig',setPTrig);
                    c.(corcomp{ncc}).(cstation{k}).corr =...
                        set(c.(corcomp{ncc}).(cstation{k}).corr,'waveform',setWav);
                end
            end
            % make threecomp-objects here
        end
    end
    
    correlationCatalog = c;
end
