function baz = backazmiuthStationToEvent(stationstr, catalog)

% returns back azimuth between a station specified by stationstr and
% events in GISMO catalog format

if exist('stations.mat', 'file') == 2
%     load('stations.mat');
    load stations.mat stations;
else
    stations = readtable('station.dat');
    stations.Properties.VariableNames = {'name' 'lat' 'lon'};
end

if isempty(stations.name) || isempty(stationstr)
    baz = [];
    return
end

[requStationID, stationSort] = ismember(stations.name, stationstr);

requestedStations = stations(requStationID,:);
stationSortOrder = stationSort(requStationID);
%[fi] = find(strcmp(stations.name, stationstr));
%[fi] = find(strcmp(requestedStations.name, stationstr));

[ism inversSortOrder] = ismember([1:1:numel(stationSortOrder)],...
    stationSortOrder);
reqStationsSorted = requestedStations(inversSortOrder,:);
lon1 = reqStationsSorted.lon;
lat1 = reqStationsSorted.lat;

%lon1 = stations.lon(fi);
%lat1 = stations.lat(fi);

lon2 = catalog.lon;
lat2 = catalog.lat;

% either multiple stations or multiple events
if length(lon1) >= 1 && numel(stationstr) >= 1
    lon2 = repmat(lon2,size(lon1));
    lat2 = repmat(lat2,size(lon1));
elseif length(lon2) >= 1 && numel(stationstr) == 1
    lon1 = repmat(lon1,size(lon2));
    lat1 = repmat(lat1,size(lon2));
else 
    error('cannot calculate backazimuth for multiple sources and stations')
end

baz = simpleAzimuth(lon1, lat1, lon2, lat2);