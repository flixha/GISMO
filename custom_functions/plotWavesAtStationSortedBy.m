function [pickOffset] = plotWavesAtStationSortedBy(c, comp, sta, varargin)

% varargin: method, plottype, linewidth, scale, 
% plottype is 'sha' or 'wig' 
SECPERDAY = 24*60*60;


p = inputParser;

defaultMethod = 'DISTANCEFROMTRENCH';
validMethods = {'LAT', 'LATITUDE', 'LON', 'LONGITUDE', 'X',...
    'DISTFROMTRENCH', 'DISTANCEFROMTRENCH', 'Y', 'DISTALONGSTRIKE',...
    'Z', 'DEPTH', 'DISTFROMSLABTOP', 'DISTANCEFROMSLABTOP', 'SLABTOP',...
    'DFST', 'S-P','SMP','SMINUSP'};
checkMethod = @(x) any(validatestring(x,validMethods));

defaultPlotType = 'WIG';
validPlotType = {'SHA', 'WIG', 'WIGBYY', 'CWIG', 'BWIG'};
checkPlotType = @(x) any(validatestring(x,validPlotType));

defaultLineWidth = 0.5;
defaultScale = 1;

defaultMarkArrivalTimes = false;

addRequired(p,'c',@isstruct);
addRequired(p,'comp',@ischar);
addRequired(p,'sta',@ischar);
addOptional(p,'method', defaultMethod, checkMethod);
addOptional(p,'plottype', defaultPlotType, checkPlotType);
addParameter(p,'linewidth', defaultLineWidth, @isnumeric)
addParameter(p,'scale', defaultScale, @isnumeric)
addParameter(p,'markArrivalTimes', defaultMarkArrivalTimes, @islogical)

defaultArrivals = struct();
addParameter(p,'arrivals', defaultArrivals, @isstruct)

defaultLabelArrivals = false;
addOptional(p,'labelArrivals', defaultLabelArrivals, @islogical);

defaultMaxConversions = 1;
addOptional(p,'maxConversions', defaultMaxConversions, @isnumeric);

defaultMaxReflections = 1;
addOptional(p,'maxReflections', defaultMaxReflections, @isnumeric);

defaultMaxInteractions = 1;
addOptional(p,'maxInteractions', defaultMaxInteractions, @isnumeric);

defaultConvAndRefl = 0;
addOptional(p,'maxConvAndRefl', defaultConvAndRefl, @isnumeric);

defaultAvoidPtoS = true;
addOptional(p,'avoidPtoS', defaultAvoidPtoS, @islogical);

defaultPlotPMS = true;
addOptional(p,'plotPMS', defaultPlotPMS, @islogical);

defaultSV_vs_P = 0;
addOptional(p,'SV_vs_P', defaultSV_vs_P, @isnumeric);

defaultSH_vs_P = 0;
addOptional(p,'SH_vs_P', defaultSH_vs_P, @isnumeric);

defaultMaxPhaseLength = 5;
addParameter(p,'maxPhaseLength', defaultMaxPhaseLength, @isnumeric)

defaultPlotEnvelope = false;
addParameter(p,'plotEnvelope', defaultPlotEnvelope, @islogical)

defaultMinTafterDirP = 0.4;
addParameter(p,'minTafterDirP', defaultMinTafterDirP, @isnumeric)

defaultColorFul = false;
addParameter(p, 'colorful', defaultColorFul, @islogical);

defaultSeismogramType = 'velocity';
validSeismogramType = {'DISPLACEMENT', 'VELOCITY'};
checkSeismogramType = @(x) any(validatestring(upper(x),validSeismogramType));
addOptional(p,'seismogramtype', defaultSeismogramType, checkSeismogramType);

parse(p, c, comp, sta, varargin{:});

res = p.Results;
method = res.method;
plottype = res.plottype;
plotEnvelope = res.plotEnvelope;
linewidth = res.linewidth;
arrivals = res.arrivals;
markArrivalTimes = res.markArrivalTimes;
maxConversions = res.maxConversions;
labelArrivals = res.labelArrivals;
maxReflections  = res.maxReflections;
maxInteractions = res.maxInteractions;
maxConvAndRefl = res.maxConvAndRefl;
avoidPtoS = res.avoidPtoS;
plotPMS = res.plotPMS;
SV_vs_P = res.SV_vs_P;
SH_vs_P = res.SH_vs_P;
maxPhaseLength = res.maxPhaseLength;
seismogramtype = res.seismogramtype;
colorful = res.colorful;
minTafterDirP = res.minTafterDirP;


invertplot = false;
isStackedCatalog = false;

if ismember('stddepth', c.(comp).(sta).cat.table.Properties.VariableNames)
    isStackedCatalog = true;
end

if maxInteractions < maxConversions + maxReflections
    maxInteractions = max(maxConversions, maxReflections);
end

% Plot not the waveforms but their envelopes/Hilbert transforms
plusEnv = 0;
if plotEnvelope
    c.(comp).(sta).corr = hilbert(c.(comp).(sta).corr);
    %w = get(c.(comp).(sta).corr,'waveform');
    % w = power(w,3/4);
    %c.(comp).(sta).corr = set(c.(comp).(sta).corr,'waveform',w);
    plusEnv = -0.5;
    
    % Reduce the sampling by factor of two
    
end

if strcmpi(seismogramtype,'DISPLACEMENT')
    c.(comp).(sta).corr = integrate(c.(comp).(sta).corr);
    c.(comp).(sta).corr = demean(c.(comp).(sta).corr);
    c.(comp).(sta).corr = detrend(c.(comp).(sta).corr); 
end

% if muteFirstArrivals
%     wavs = get(c.(comp).(sta).corr,'waveform');
%     ptriggers = c.(comp).(sta).corr,'trig');
% need to find a way to store the striggers in the correlation object...
%     striggers = 
%     for j=1:1:numel(wavs)
%         preP = taper(wavs(j),[-4, ptriggers(j)-0.1]);
%         postP = taper(wavs(j),[ptriggers(j)+0.2 striggers(j)-0.1]);
%         postS = taper(wavs(j),[striggers(j)+0.2 20]);
%     end
% end

switch upper(method)
    case {'LAT','LATITUDE'}
        %sort by lat
        [sortedLats,sortLatIndices] = sort(c.(comp).(sta).cat.lat);
        traceOrder = sortLatIndices;
        sortedValues = sortedLats;
        titlename = 'latitude';
        ylabelname = titlename;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdlat(traceOrder);
        end
    case {'LON','LONGITUDE'}
        %sort by lat
        [sortedLons,sortLonIndices] = sort(c.(comp).(sta).cat.lon);
        traceOrder = sortLonIndices;
        sortedValues = sortedLons;
        titlename = 'longitude';
        ylabelname = titlename;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdlon(traceOrder);
        end        
    case {'X','DISTFROMTRENCH','DISTANCEFROMTRENCH'}
        %sort by distance from trench greekevents.table.x
        %- this sorting might be preferrable, conversion waveforms align quite well
        [sortedXs,sortXIndices] = sort(c.(comp).(sta).cat.table.x);
        traceOrder = sortXIndices;
        sortedValues = sortedXs;
        titlename = 'distance from trench';
        ylabelname = [titlename, ' (km)'];
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdx(traceOrder);
        end
    case {'Y','DISTALONGSTRIKE'}
        %sort by distance along strike (along trench) greekevents.table.x
        [sortedYs,sortYIndices] = sort(c.(comp).(sta).cat.table.y);
        traceOrder = sortYIndices;
        sortedValues = sortedYs;
        titlename = 'distance along strike';
        ylabelname = [titlename, ' (km)'];
        invertplot = true;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdy(traceOrder);
        end
    case {'Z','DEPTH'}
        %sort by depth greekevents.table.z
        [sortedDepths, sortDepthIndices] = sort(c.(comp).(sta).cat.table.depth);
        traceOrder = sortDepthIndices;
        sortedValues = sortedDepths;
        titlename = 'depth';
        ylabelname = [titlename, ' (km)'];
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stddepth(traceOrder);
        end
    case {'DISTFROMSLABTOP','DISTANCEFROMSLABTOP','SLABTOP','DFST'}
        %sort by distance from slab Top
        [sortedDFSTs,sortDFSTIndices] = sort(...
            c.(comp).(sta).cat.table.DistFromSlabTop,'ascend');
        traceOrder = sortDFSTIndices;
        sortedValues = sortedDFSTs;
        titlename = 'distance from slab top';
        ylabelname = [titlename, ' (km)'];
        invertplot = true;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdDistFromSlabTop(traceOrder);
        end
    case {'S-P','SMP','SMINUSP'}
        %sort by S-P
        nsmp = size(c.(comp).(sta).cat.arrivals,1);
        SmPtime = zeros(nsmp,1);
        for j=1:1:nsmp
            %if there are any arrivals at all for that station and event
            if size((c.(comp).(sta).cat.arrivals{j}),2) > 0
                %if there are both p- and s- arrivals
                if ~isempty(c.(comp).(sta).cat.arrivals{j}.s_time) &&...
                        ~isempty(c.(comp).(sta).cat.arrivals{j}.p_time)
                    SmPtime(j) = c.(comp).(sta).cat.arrivals{j}.s_time-...
                        c.(comp).(sta).cat.arrivals{j}.p_time;
                end
            else
                SmPtime.(sta)(j) = NaN;
            end
        end
        SmPtime = SmPtime * SECPERDAY;
        [sortedSmP,sortSmPIndices] = sort(SmPtime);
        traceOrder = sortSmPIndices;
        sortedValues = sortedSmP;
        titlename = 'S-minus-P arrival time';
        ylabelname = [titlename, ' (s)'];
    case {'MAG','M','MAGNITUDE'}
        %sort by distance from slab Top
        [sortedMags,sortMagIndices] = sort(...
            c.(comp).(sta).cat.mag,'descend');
        traceOrder = sortMagIndices;
        sortedValues = sortedMags;
        titlename = 'Magnitude';
        ylabelname = titlename;
        %invertplot = true;
end


if ~any(contains(p.UsingDefaults,'scale'))
    scale = p.Results.scale;
    if strcmpi(plottype,'WIGBYY')
        markScale = abs(sortedValues(end) - sortedValues(1))/numel(...
            traceOrder)*1;
    else
        markScale = 1;
    end
end


% amplitude factor that needs to be met to NOT plot P-to-S or S-to-P
tolerance = 2;
% initiate the Amplitude ratio variable for later checks
if length(SV_vs_P) ~= length(traceOrder)
    SV_vs_P = NaN(length(traceOrder),1);
else
    avoidPtoS = false;
end

if length(SH_vs_P) ~= length(traceOrder)
    SH_vs_P = 10 * ones(length(traceOrder),1);
else
    avoidPtoS = false;
end

% plot(c.(comp).(sta).corr,'wig',1,traceOrder);
% plot(c.(comp).(sta).corr,'sha',1,traceOrder);
% plot(c.(comp).(sta).corr, plottype, 1, traceOrder);
% To plot .....
% plot(c.(comp).(sta).corr,'wig',0.05,traceOrder,c.(comp).(sta).cat.table.x,'linewidth',0.1);

% Switch around waveforms and plot when the chosen sort value should go
% from large at the top to small at the bottom (that applies only to the 
% distance from slab top for now)
if invertplot
    if contains(upper(plottype),'WI')
        w = get(c.(comp).(sta).corr,'waveform');
        w = times(w,-1);
        c.(comp).(sta).corr = set(c.(comp).(sta).corr,'waveform',w);
    end
    ydirection = 'normal';
    plusEnv = -1 * plusEnv;
else
    ydirection = 'reverse';
end
            
switch upper(plottype)
    case {'WIG', 'SHA', 'CWIG', 'BWIG'}
        % Find a suitable value for the scale of wiggles if none provided:
        if ~exist('scale', 'var')
            scale = 1;
            markScale = 1;
        end
        plot(c.(comp).(sta).corr, plottype, scale, traceOrder);
        
        ylabel(ylabelname,'fontsize',12)
        labelValues = sortedValues(get(gca,'YTick'))';
        if size(labelValues,1) < size(labelValues,2)
            labelValues = labelValues';
        end
        
        % If this is a stacked object, then display number of traces in
        % each stack and standard deviations
        nlabels = numel(labelValues);
        if isStackedCatalog
            newLabels = [ num2str(labelValues, '%6.2f'),...
                repmat('\pm: ', nlabels, 1),...
                num2str(stdValues, '%3.1f'),...
                repmat(' n: ',nlabels,1),...
                num2str(c.(comp).(sta).cat.table.nStackedTraces(traceOrder),...
                '%2.0f')];
        else
            newLabels = num2str(labelValues, '%6.2f');
        end
        
        set(gca,'YTickLabel', newLabels);
    case 'WIGBYY'
        if ~exist('scale', 'var')
           scale = abs(sortedValues(end) - sortedValues(1))/numel(traceOrder)*1;
           markScale = scale;
        end
        % Set y-axis to the sorted parameter, i.e. plot wiggle against that
        % sort-parameter
        plot(c.(comp).(sta).corr, 'wig', scale, traceOrder, sortedValues,...
           'linewidth', linewidth);
%         plot(c.(comp).(sta).corr, 'wig', scale, sortedValues,...
%             'linewidth', linewidth);
        
        addsome = abs(sortedValues(end) - sortedValues(1))*0.02;      
        ylabel(ylabelname,'fontsize',12)
        
        ylim([min(sortedValues)-addsome, max(sortedValues)+addsome]);
        set(gca,'YTickMode','auto')
        set(gca,'YTickLabelMode','auto')
        set(gca,'YTickLabelMode','auto')
        set(gca,'YMinorTick','on')
end

%logical array to store whether arrivals were labelled on the last few
%lines or no
nLines = numel(traceOrder);
minLinesWithoutLabel = round(nLines/8);
lastXLines = false(minLinesWithoutLabel,1);

ax1 = gca;
hold(ax1, 'on')

if markArrivalTimes
%     load('arrivals.mat');
%     load('arrivals_VgridsDepthCorrectlyInverted.mat');
%     load('arrivals_wNconv.mat');
    % check if arrivals were sent to function, else load from file
    if isempty(fields(arrivals))
        load('arrivals_0.02degGrid.mat','arrivals');
    end
    
    prevEventArrivalPhases = cell(0,0);
    
    % Go from the bottom up if the plot is inverted at the end so that the
    % phase-text appears at the top
    if invertplot
        traceOrder = flip(traceOrder);
        sortedValues = flip(sortedValues);
    end
    
    pickOffset = zeros(nLines,1);
    
    p = 0;
    for j=traceOrder'
        pickColor = zeros(0, 3);
        p = p+1;
        evID = c.(comp).(sta).cat.table.EventID(j);
        
        if arrivals(evID).EventID ~= evID
            warning('event IDs do not match')
            continue
        end
            
            
        thisEventArrivalPhases = cell(0,0);
        thisEventArrivalTimes = zeros(0,0);
        %Check that there
%         if traceOrder(j) > length(arrivals)
%             continue
%         end
                
        % get the right line for plotting
        [arrivalLine, arrIdx] = find(traceOrder==j);
        if invertplot
            arrivalLine = length(traceOrder) - arrivalLine + 1;
        end
        
        if strcmpi(plottype,'WIGBYY')
            arrivalLine = sortedValues(arrIdx);
        end

        % plot arrivals 
        % check that theere are arrivals for that event
        if isempty(arrivals(evID).arrivals)
            continue
        end
        % check that there are arrivals for that station 
        if isempty(arrivals(evID).arrivals.(sta).time)
            continue
        end
        
        narrivals = numel(arrivals(evID).arrivals.(sta).time);

        % Two ways to find the first arrival: from the FM3D paths
        % or from the picks (triggers)
        firstArrival = getFirstArrivalTimeFromFM3D(c, arrivals,...
            sta, comp, traceOrder, evID, arrIdx);

        triggers = get(c.(comp).(sta).corr,'trig');
        firstArrivalFromPicks = (triggers(j) - c.(comp).(sta...
            ).cat.table.otime(j)) * SECPERDAY;
        pickOffset(j) = firstArrival - firstArrivalFromPicks;
        
        firstArrivalS = getFirstArrivalTimeFromFM3D(c, arrivals,...
            sta, comp, traceOrder, evID, arrIdx, 'Phase', 'S');

        %noe plot every arrival
        q = 0;
        for k=1:1:narrivals
            nConv = arrivals(evID).arrivals.(sta).nConv(k);
            nRefl = arrivals(evID).arrivals.(sta).nReflections(k);
            time = arrivals(evID).arrivals.(sta).time(k);
            % Check the number of conversions along a path, and
            % only allow plotting when it is less than
            % maxConversions
            thisPhase = char(arrivals(evID).arrivals.(sta).phase(k));
%                     plotPhase = true; nConv = 0;
%                     pIdx = 1:3:length(thisPhase);
%                     for p=2:1:length(pIdx)
%                         if ~strcmp(thisPhase(pIdx(p)),...
%                                 thisPhase(pIdx(p-1)))
%                             nConv = nConv + 1;
%                         end
%                     end
            plotPhase = true;

            % do not plot when phase contains too many conversions
            if nConv > maxConversions
                plotPhase = false;
            % don't plot if too many reflections
            elseif nRefl >...
                    maxReflections
                plotPhase = false;
                
            % avoid S-to-P if too little energy from both SH and SV
            elseif SV_vs_P(p) < 1/tolerance && SH_vs_P(p) < 1*tolerance &&...
                strcmp(thisPhase(end),'P') && nConv >= 1 && nRefl == 0 
                plotPhase = false;
            %avoid P-to-S if too little energy from P
            elseif SV_vs_P(p) > 1*tolerance  &&...
                strcmp(thisPhase(end),'S') && nConv >= 1 && nRefl == 0 
                plotPhase = false;
                
            % don't plot if it is a P-to -S conversion
            elseif avoidPtoS && nConv >= 1 && strcmp(thisPhase(end),'S')...
                    && ~(plotPMS && strcmp(thisPhase,'P_MS'))
                % exception from the rule: PMS-arrival, or if this
                % station receiver little S-energy compared to P-energy,
                % and this is a conversion without reflection anywhere
                plotPhase = false;
                                
            % if polarizationfiltered: don't plot the arrivals that
            % shouldn't be visible on that channel
            %         elseif strcmp(comp,'Zp') && strcmp(thisPhase(end),'S')
            %             plotPhase = false;
            %         elseif (contains(comp,'R') || contains(comp,'T')) &&...
            %                 strcmp(thisPhase(end),'P')
            %             plotPhase = false;
            elseif isThisPhaseTooLong(thisPhase, maxPhaseLength,...
                    c.(comp).(sta).cat.table.DistFromSlabTop(j))
                plotPhase = false;
            % Don't plot if it is very shortly after P or S (i.e. for 
            % interface EQs)
            elseif (nConv >= 1 || nRefl >= 1) &&...
                    (time - firstArrival <= minTafterDirP ||...
                    abs(time - firstArrivalS) <= minTafterDirP)
                plotPhase = false;                
            % don't plot if too many interactions
            elseif  nConv + nRefl > maxInteractions
                plotPhase = false;
                % still plot if reflection and conversion was the
                % same interaction (but this can be turned off)
                % need to be no more than x interactions and x
                % reflections
                if nConv==maxConversions && nRefl == maxReflections
                    if nConv + nRefl <= maxConvAndRefl &&...
                            are_reflec_and_convers_oneInteract(thisPhase)
                        plotPhase = true;
                    end
                end
            end

            if plotPhase
                q = q + 1;
                thisEventArrivalPhases{q,1} = char(arrivals(...
                    evID).arrivals.(sta).phase(k));
                thisEventArrivalTimes(q,1) = arrivals(evID).arrivals.(...
                    sta).time(k) - firstArrival;
            end


            if plotPhase    
                %plot(ax1, arrivals(evID).arrivals.(sta).time(k) -...
                %    firstArrival,arrivalLine, '+k');
                x = arrivals(evID).arrivals.(sta).time(k) - firstArrival;
                if x < 0
                    warning(['Arrival earlier than first arrival,',...
                        'something went wrong']);
                end
                y = arrivalLine;
                if strcmpi(plottype,'wigbyy')
                    y = sortedValues(j);
                end
                
                % chose color for picks
                if strcmpi(plottype,'CWIG')
                    pickColor = [pickColor; 1 0.64 0];
                elseif strcmpi(plottype, 'WIG') ||...
                    strcmpi(plottype, 'WIGBYY') ||...
                        strcmpi(plottype,'BWIG')
                    if colorful
                        % reflections: red
                        if nRefl >= 1 && nConv == 0
                            pickColor = [pickColor; 1 0 0];
                        %conversions: blue
                        elseif nConv >= 1 && nRefl == 0
                            pickColor = [pickColor; 0 0 1];
                        %reflected & converted: lila
                        elseif nRefl >= 1 && nConv >= 1
                            %pickColor = [pickColor; 1 0.4 1];
                            % strong lila:
                            %pickColor = [pickColor; 0.6 0.2 0.85]; 
                            % weak lila:
                            pickColor = [pickColor; 0.8 0.55 1]; 
                        else
                        % direct: green
                            pickColor = [pickColor; 0.15 0.7 0.25];
                        end
                    else
                        pickColor = [pickColor; 1 0 0];
                    end
                else
                    pickColor = [pickColor; 0 0 0];
                end

                plot(ax1, [x, x], [y-markScale/2+plusEnv...
                    y+markScale/2+plusEnv],...
                    'color', pickColor(q,:), 'LineWidth', 0.5);
            end
        end

        
        if ~labelArrivals
            continue
        end
        % Mark the phases with text only if the
        % previous event did not have the same arrivals
        plotText = true;
        % thisEventArrivalPhases = arrivals(evID).arrivals.(sta).phase;
        if length(thisEventArrivalPhases) == length(prevEventArrivalPhases)
            phasesEqual = strcmp(prevEventArrivalPhases,...
                thisEventArrivalPhases);
            if all(phasesEqual)
                plotText = false;
            end
        end
        %mark arrivals only for the first MW, first IF, first crustal,first
        % slabmantle EQ; don't label when the last 10 lines weren't empty
        if any(lastXLines)
            plotText = false;
        end

        if plotText
            %Y = repmat(y, length(thisEventArrivalTimes),1);
            %text(ax1, thisEventArrivalTimes, Y,...
            %            thisEventArrivalPhases,...
            %            'color', pickColor,'Rotation', 90,...
            %            'FontSize', 10);
            [sortedEventArrivalTimes, sortI] = sort(thisEventArrivalTimes);
            sortedEventPhases = thisEventArrivalPhases(sortI,:);
            sortedPickColors = pickColor(sortI,:);
            for k=1:1:length(thisEventArrivalTimes)
                %if k>1
                %    pickDist = sortedEventArrivalTimes(k) -...
                %        sortedEventArrivalTimes(k-1);
                %else
                %    pickDist = 99;
                %end
                %get a shortened phase string
                cat = subset(c.(comp).(sta).cat, j);
                shortPhase = shorten_phase_descriptor(cat,...
                    sortedEventPhases{k}, sortedEventArrivalTimes(k));
                %if there is long enough between arrivals, then print
                % with defautl alignment
                if ~strcmp(shortPhase,'')
                    if mod(k,2) || length(traceOrder) < 15
                        text(ax1, sortedEventArrivalTimes(k),...
                            y + scale/2, shortPhase,...
                            'color', sortedPickColors(k,:),'Rotation', 90,...
                            'FontSize', 10, 'HorizontalAlignment','left');
                    else
                        text(ax1, sortedEventArrivalTimes(k),...
                            y - scale/2, shortPhase,...
                            'color', sortedPickColors(k,:),'Rotation', 90,...
                            'FontSize', 10, 'HorizontalAlignment','right');
                    end
                end

                %     if pickDist > 1
                %         text(ax1, sortedEventArrivalTimes(k), y,...
                %                 sortedEventPhases(k),...
                %                 'color', pickColor,'Rotation', 90,...
                %                 'FontSize', 10, 'HorizontalAlignment','right');
                %     else
                %         % print offset to the left/bottom
                % 
                % if mod(k,2)
                %     text(ax1, sortedEventArrivalTimes(k), y,...
                %             sortedEventPhases(k),...
                %             'color', pickColor,'Rotation', 90,...
                %             'FontSize', 10, 'HorizontalAlignment','left');
                % else
                %     text(ax1, sortedEventArrivalTimes(k), y,...
                %         sortedEventPhases(k),...
                %         'color', pickColor,'Rotation', 90,...
                %         'FontSize', 10, 'HorizontalAlignment','right');
                % end
            end
        end
        lastXLines(1:end-1,1) = lastXLines(2:end,1);
        lastXLines(end,1) = plotText;
        
        prevEventArrivalPhases = thisEventArrivalPhases;
    end
end

%change axis scale for shaded plotting
if strcmpi(plottype,'SHA') 
    if scale~=1
        caxis(caxis./scale)
    end
    if plotEnvelope
        if invertplot
            colormap(flipud(hot))
        else
            colormap(hot)
        end
    end
end

% some plot format improvements 
set(gca, 'YDir', ydirection);
ax1 = gca;


title([sta, '.', comp, ' sorted by ', titlename],'fontsize',14)
ax1.Title.Position(2) = ax1.Title.Position(2) * 1.05;

% or fix the title relative to y-axis values if "wigbyy"
if strcmpi(plottype, 'wigbyy')
    bounds = get(gca,'YLim');
    if invertplot
        ax1.Title.Position(2) = bounds(2) + diff(bounds)*0.06;
    else
        ax1.Title.Position(2) = bounds(1) - diff(bounds)*0.06;
    end
end

nYTicks = length(ax1.YTickLabel);
maxNYTicks = 50;
if nYTicks > maxNYTicks
    yTickStep = round(nYTicks/maxNYTicks);
    newYTickLabels = ax1.YTickLabel(1:yTickStep:end,:);
    ax1.YTick = ax1.YTick(1:yTickStep:end);
    ax1.YTickLabel = newYTickLabels;
end
    
ax1.Position = [0.1 0.03 0.89 0.9];
ax1.XMinorTick = 'on';
ax1.XGrid = 'on';
xlim([-4 19])

fig1 = gcf;
oldpos = fig1.Position;
fig1.Position(4) = oldpos(4);
