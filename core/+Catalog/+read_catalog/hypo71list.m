function self = hypo71list(zmapdata)
%readEvents.zmap Translate a hypo71 inputlist into a Catalog object
% Hypo71 has xx columns: longitude, latitude, decimal year, month, day,
% magnitude, depth, hour, minute, second

hypo71list










%     lon = zmapdata(:,1);
%     lat = zmapdata(:,2);
%     time = datenum( floor(zmapdata(:,3)), zmapdata(:,4), zmapdata(:,5), ...
%         zmapdata(:,8), zmapdata(:,9), zmapdata(:,10) );
%     mag = zmapdata(:,6);
%     depth = zmapdata(:,7);
%     self = Catalog(time, lon, lat, depth, mag, {}, {});
% create Catalog object
    self = Catalog(dnum', lon', lat', depth', mag', magtype', etype', ...
        'request', request, 'arrivals', arrivals');
end
