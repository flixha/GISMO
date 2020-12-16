function [traceOrder, sortedValues, titlename, ylabelname, stdValues,...
    invertplot] = getTraceOrderAndPlotMethodProps(c, method, comp, sta,...
    isStackedCatalog)
% getTraceOrderAndPlotMethodProps sorts traces by a chosen method and
%   returns the sorted order of the traces.
% Input:
%   c       : network-correlation object
%   method  : one of 'LAT' or 'LATITUDE',
%                   'LON' or 'LONGITUDE',
%                   'X' or 'DISTFROMTRENCH' or 'DISTANCEFROMTRENCH'
%                   'Y' or 'DISTALONGSTRIKE', 
%                   'Y' or 'DISTALONGSTRIKE',
%                   'Z' or 'DEPTH',
%                   'DISTFROMSLABTOP' or 'DISTANCEFROMSLABTOP' or 
%                       'SLABTOP','DFST',
%                   'S-P' or 'SMP' or 'SMINUSP',
%                   'MAG' or 'M' or 'MAGNITUDE'.
%   sta     : cell array with the names of the chosen stations
%   comp    : cell array with the names of the chosen channels
%   isStackedCatalog: indicated whether there are stacked events in the
%             catalog.
% Output:
%   traceOrder  : array of the sorted trace indices
%   sortedValues: array of the sorted values by which the method sorted
%   titlename   : str of the name by which the traces were sorted
%   ylabelname  : str of the name by which the traces were sorted
%   stdValues   : array of standard deviation values if more than one event
%               are considered for the sorting of a single trace ( is the
%               case for stacked traces).
%   invertplot  : whether the order is increasing or decreasing

    stdValues = zeros(height(c.(comp).(sta).cat.table),1);
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