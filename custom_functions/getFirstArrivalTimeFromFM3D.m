function firstArrival = getFirstArrivalTimeFromFM3D(c, arrivals,...
    sta, comp, traceOrder, evID, arrivalLine, varargin)

% getFirstArrivalTimeFromFM3D retrieves the first arrival time for a 
%   station-event pair from the a set of arrival times (here originally
%   supplied by the FM3D software by Rawlinson, Sambridge, deKool.
% Input:
%   c        : network-correlation object
%   sta      : cell array with the names of the chosen stations
%   comp     : cell array with the names of the chosen channels
%   traceOrder: array of indices of the traces
%   evID     : ID of the event
%   arrivalLine: the line in the plot where the arrival will need to be
%       plotted into to fit the trace.
% Output:
%   firstArrival : arrival object for the first arrival

    p = inputParser;

%     addRequired(p,'c',@isstruct);
%     addRequired(p,'arrivals',@isstruct);
%     addRequired(p,'sta',@ischar);
%     addRequired(p,'comp',@ischar);
%     addRequired(p,'traceOrder',@isnumeric);
%     addRequired(p,'evID',@isnumeric);
%     addRequired(p,'arrivalLine',@isnumeric);
    
    defaultPhase = 'P';
    validPhases = {'P','S'};
    checkPhases = @(x) any(validatestring(x,validPhases));
    addOptional(p,'Phase', defaultPhase, checkPhases);
    
    % parse(p, c, arrivals, sta, comp, traceOrder, evID, arrivalLine, varargin{:});
    parse(p, varargin{:});
    res = p.Results;
    Phase = res.Phase;
    
    [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase, arrivals, evID, sta);
        
    firstArrivalPhases = firstArrivalPhase(1:3:end)';
    if strcmp(firstArrivalPhases,repmat(Phase,...
            length(firstArrivalPhases),1))
        firstArrival = arrivals(evID).arrivals.(sta).time(idx);
    else
        % find an alternative first arrival from neighboring
        % events
        warning(['Earliest arrival is NOT all P-phases,',... 
            'did FM3D find no path?']);

        jm1 = arrivalLine;
        firstArrM1 = nan;
        while jm1 > 1
            jm1 = jm1 - 1;
            evIDm1 = c.(comp).(sta).cat.table.EventID(traceOrder(jm1));

            if ~isempty(arrivals(evIDm1).arrivals)
                if ~isempty(arrivals(evIDm1).arrivals.(sta).time)
                    [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase,...
                        arrivals, evIDm1, sta);
                    firstArrivalPhases = firstArrivalPhase(1:3:end)';
                    if strcmp(firstArrivalPhases,repmat(Phase,...
                            length(firstArrivalPhases),1))
                        firstArrM1 = arrivals(evIDm1).arrivals.(sta).time(idx);
                        continue
                    end
                end
            end
        end

        jp1 = arrivalLine;
        firstArrP1 = nan;
        while jp1 < c.(comp).(sta).cat.numberOfEvents - 1
            jp1 = jp1 + 1;
            evIDp1 = c.(comp).(sta).cat.table.EventID(traceOrder(jp1));

            if ~isempty(arrivals(evIDp1).arrivals)
                if ~isempty(arrivals(evIDp1).arrivals.(sta).time)
                    [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase,...
                        arrivals, evIDp1, sta);
                    firstArrivalPhases = firstArrivalPhase(1:3:end)';
                    if strcmp(firstArrivalPhases,repmat(Phase,...
                            length(firstArrivalPhases),1))
                        firstArrP1 = arrivals(evIDp1).arrivals.(sta).time(idx);
                        continue
                    end
                end
            end
        end                                      

        firstArrival = nanmean([firstArrM1, firstArrP1]);
    end

function [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase, arrivals,...
    evID, sta)
% getFirstArrivalPhase return the first arriving phase of a chosen type
%   (e.g., 'P' or 'S') for one event and a selected station.
% Input:
%   Phase   : phase type, e.g., 'P' or 'S'
%   arrivals: set of arrivals
%   evID    : ID of the event
%   sta     : station name
% Output:
% firstArrivalPhase : the first arriving phase
% idx               : the index of the first arriving phase in the input
%                       set of arrivals

    if strcmp(Phase,'P')
        firstArrivalPhase = char(arrivals(evID).arrivals.(sta).phase(1));
    elseif strcmp(Phase,'S')
        noP = ~contains(arrivals(evID).arrivals.(sta).phase,'P');
        phasesWithoutP = arrivals(evID).arrivals.(sta).phase(noP);
        firstArrivalPhase = char(phasesWithoutP(1));
    end
    idx = find(ismember(arrivals(evID).arrivals.(sta).phase,...
        firstArrivalPhase));
    
    