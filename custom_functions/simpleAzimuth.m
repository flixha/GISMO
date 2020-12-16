function az = simpleAzimuth(lon1, lat1, lon2, lat2)

% simpleAzimuth returns the azimuth between two locations on Earth.
% Input:
%   lon1: longitude at location 1
%   lat1: latittude at location 1
%   lon2: longitude at location 2
%   lat2: latittude at location 2
% Output:
%   az  : azimuth from location 1 to location 2

    phi1 = deg2rad(lat1);
    phi2 = deg2rad(lat2);
    lambda1 = deg2rad(lon1);
    lambda2 = deg2rad(lon2);
    dphi = deg2rad(lat2 - lat1);
    dlambda = deg2rad(lon2 - lon1);

    % bearing:
    y = sin(lambda2 - lambda1) .* cos(phi2);
    x = cos(phi1).*sin(phi2) - sin(phi1).*cos(phi2).*cos(lambda2 - lambda1);
    az = atan2(y, x);

    az = rad2deg(az);