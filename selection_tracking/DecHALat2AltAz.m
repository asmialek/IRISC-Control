function [ Alt, Az ] = DecHALat2AltAz( Dec, HA, Lat )
% converts Declination, Hour Angle, Latitude into Alt, Az
% coordinates (output in degrees)

% convert input data from decimal degree to radians 
Dec = Dec*pi/180;
HA = HA*pi/180;
Lat = Lat*pi/180;

Alt = asin(sin(Dec)*sin(Lat)+cos(Dec)*cos(Lat)*cos(HA));

Az = acos((sin(Dec) - sin(Alt)*sin(Lat)) / (cos(Alt)*cos(Lat)));

if(sin(HA) > 0)
    Az = 2*pi - Az;
end
    
Alt = Alt * 180/pi;
Az = Az *180/pi;

end

