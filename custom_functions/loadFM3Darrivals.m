function arrivals = loadFM3Darrivals(fm3Dpath, varargin)

% loadFM3Darrivals loads a set of arrivals from theoretical calculations in
%   FM3D (fast marching 3D program by Rawlinson, Sambridge, deKool).
% Input:
%   fm3Dpath: path to the FM3D-run
% Output:
%   arrivals: GISMO arrivals object

%check if it should read from an folder-per-event structure
folderPerEvent = false;

if length(varargin) >= 1
    jEvents = varargin{1};
    folderPerEvent = true;
    if size(jEvents,1) > size(jEvents,2)
        jEvents = jEvents';
    end
else
    jEvents = 1;
end

%loop through the requested events  (folders)
arrivals = struct();
kev = 0;
for jev = jEvents
    kev = kev + 1;
    disp(['Reading event ', num2str(kev), ' which is event with ID ',...
        num2str(jev)])
    
    
    if kev==1
        arrivals(kev,1) = struct();
        arrivals(kev,1).arrivals = struct();
    else
        % copy structure of the first
        arrivals(kev,1) = arrivals(1,1);
        %delete all arrivals
        arrivals(kev,1).arrivals = struct();
    end
    
    if folderPerEvent
        eventFolder = ['/Event_', sprintf('%03d', jev),'/'];
    else
        eventFolder = '';
    end
    
    fm3Dsubpath = [fm3Dpath, eventFolder];

      
    recIn = dlmread([fm3Dsubpath,'/receivers.in']);
    stations = readtable([fm3Dsubpath,'/stations.dat']);
    stations.Properties.VariableNames = {'name' 'lat' 'lon' 'elev'};

    sourceData = dlmread([fm3Dsubpath,'/sources.in']);

    % read in source
    sourceLat = sourceData(3,2);
    sourceLon = sourceData(3,3);
    sourceDep = sourceData(3,1);
    nSourcepaths = sourceData(4,1);
    
    arrivals(kev).lat = sourceLat;
    arrivals(kev).lon = sourceLon;
    arrivals(kev).depth = sourceDep;
    arrivals(kev).EventID = jev;

    arrivalFile = [fm3Dsubpath,'/arrivals.dat'];
    if exist(arrivalFile, 'file') == 2
        try
            arrivalData = readtable(arrivalFile);
            arrivalData.Properties.VariableNames = {'receiver' 'source'...
                'path', 'refl' 'time' 'isDiffracted' 'isHeadwave'};
        catch 
            disp(['Could not read ', arrivalFile])
            arrivals(kev,1).arrivals = zeros(0,0);
            continue
        end
    else
        arrivals(kev,1).arrivals = zeros(0,0);
        continue
    end
    
    phases = struct();
    phases.name = '';
    phases.time = zeros(0,0);
    phaseNames = cell(0,0);
    phases.nReflections = zeros(0,0);
    phases.nConversions = zeros(0,0);

    % create phase names like 'PtPMP' / or P5P4P
    sourcePath = cell(nSourcepaths,1);
    sourcePhases = cell(nSourcepaths,1);
    for j=1:1:nSourcepaths
        nSegmentsPerPath = sourceData(3*(j-1) + 5,1);
        sourcePath{j,1} = sourceData(3*(j-1) + 6,:);
        sourcePhases{j,1} = sourceData(3*(j-1) + 7,:);

        %initialise empty phase name
        phaseName = repmat(' ', 1, 2*nSegmentsPerPath - 1);
        for k=1:1:nSegmentsPerPath 
            if sourcePath{j,1}(k*2) ~= 0
                % call it 'P' or 'S'
                switch sourcePhases{j,1}(k)
                    case 1
                        phaseNow = 'P';
                    case 2
                        phaseNow = 'S';
                end

                % call it by a number of by the interface
                switch sourcePath{j}((k*2))
                    case 2
                        boundary = '_T';
                    case 3
                        boundary = '_M';
                    case 4
                        boundary = '_t';
                    case 5
                        boundary = '_m';
                end
                g = 3*k - 2;
                phaseName(g:g+2) = [phaseNow, boundary];

            end
        end
        phaseName = phaseName(1:end-2);
        phaseNames{j,1} = phaseName;

    end

    % read in receivers

    nrec = recIn(1,1);
    receivers.lon = zeros(nrec,1);
    receivers.lat = zeros(nrec,1);
    receivers.dep = zeros(nrec,1);
    receivers.name = strings(nrec,1);

    for j=1:1:nrec
        receivers.lat(j) = recIn(2 + (j-1) * 4, 2);
        receivers.lon(j) = recIn(2 + (j-1) * 4, 3);
        receivers.dep(j) = recIn(2 + (j-1) * 4, 1);
    end

    % find receiver name in station-list by their location
    for j=1:1:nrec
        stIdx = find ( stations.lat == receivers.lat(j) &...
            stations.lon == receivers.lon(j));
        receivers.name(j) = stations.name(stIdx(1));
    end

    narrivals = height(arrivalData);

    for j=1:1:nrec
        staName = receivers.name(j);
        arrivals(kev).arrivals.(staName) = struct();
        arrivals(kev).arrivals.(staName).phase = cell(0,0);
        arrivals(kev).arrivals.(staName).time = zeros(0,0);
        arrivals(kev).arrivals.(staName).nReflections = zeros(0,0);
    end

    % Find the valid arrivals that are not "duplicates" (i.e. same arrival
    % time, probably because P or S-changed in a zero-thickness layer)
    for j=1:1:nrec
        staName = receivers.name(j);
        idx = find(arrivalData.receiver == j);
        arrivalsAtThisSta = arrivalData(idx, :);
        nPathsToRec = height(arrivalsAtThisSta);
        p = 1;

        %for k=1:1:nSourcepaths
        for k=1:1:nPathsToRec
            if arrivalsAtThisSta.time(k) > 0
                simultaneousArrival = false;

                nPrevArrivals = length(arrivals(kev).arrivals.(staName).time);
                newArrivalTime = arrivalsAtThisSta.time(k);
                for m=1:1:nPrevArrivals
                    if arrivals(kev).arrivals.(staName).time(m) ==...
                            newArrivalTime
                        simultaneousArrival = true;
                    end
                end
                if ~simultaneousArrival
                    arrivals(kev).arrivals.(staName).phase(p,1) =...
                        cellstr(phaseNames(arrivalsAtThisSta.path(k),:));
                    arrivals(kev).arrivals.(staName).time(p,1) = newArrivalTime;

                    thisPhase = char(...
                        arrivals(kev).arrivals.(staName).phase(p,1));                    
                    % if this phase descirptor is longer than the first
                    % phase of that event, then call it a reflected arrival
                    lengthDiff = length(thisPhase) -...
                        length(arrivals(kev).arrivals.(staName).phase{1});
                    %return 1 if lengthDiff is 3, return 2, 3... if
                    %lengthDiff is 9, 15, 21 etc.
                    nRefl = round(lengthDiff/3);
                    if lengthDiff > 0
                        nRefl = 1;
                    end                  
                    
                    arrivals(kev).arrivals.(staName).nReflections(p,1) =...
                        nRefl;
                    
                    % Check the number of conversions along a path, and
                    % only allow plotting when it is less than
                    % maxConversions
                    nConv = 0;
                    pIdx = 1:3:length(thisPhase);
                    for m=2:1:length(pIdx)
                        if ~strcmp(thisPhase(pIdx(m)),...
                                thisPhase(pIdx(m-1)))
                            nConv = nConv + 1;
                        end
                    end
                    arrivals(kev).arrivals.(staName).nConv(p,1) = nConv;
                    
                    p = p+1;
                end
            end
        end

    end
end
