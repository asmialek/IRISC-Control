% test tracking

%% input 1
Dec = 41.2667;
HA = 101.5602;
Lat = 67.8557;

%% input 2
Dec = 36.466667;
HA = 54.382617;
Lat = 52.5;


%% test scipt DecHALat2AltAz
[Alt, Az] = DecHALat2AltAz( Dec, HA, Lat )
