% Selection and Tracking


%% Input
RADec_Star_Tracker = [25.4, 30.0, 15.0];  % RA (deg), Dec, Roll position of Star Tracker (decimal)
Angle_Enc = [5.6, 27.8, 10.4];            % offset from zero position by the encoders (decimal)

% GPS position (latitude, longitude (decimal))
GPS_pos = [67.85572, 20.22513];            % test values: Kiruna
% time (hours, minutes, seconds)
RPi_time = [10, 35, 22];
% date (year, month, day)
RPi_date = [2019, 07, 15];



%% Constants
%target_list_RADec = zeros(10,5);  % list of all celestial targets (index, RA, Dec, Mag, type_prioritisation)
%target_list_RADec = importdata('input/target_list.txt', ',');
load('target_list.mat')
FoV_gimbal = 180;         % operational FoV of the gimbal (symmetrical)

%% Conversion of variables
target_list_RADec(:,2) = target_list_RADec(:,2) * 15;   % convert RA from hours to degrees

UT = ((RPi_time(3)/60)+RPi_time(2))/60 + RPi_time(1);                  % UT from RPi_time
d_J2000 =  date2J2000(RPi_date(1), RPi_date(2), RPi_date(3), RPi_time(1), RPi_time(2), RPi_time(3));         % format J2000, calculated by RPi_time and RPi_date

LST = 100.46 + 0.985647 * d_J2000 + GPS_pos(2) + 15*UT; % Local Sidereal Time -> correct formula, accuracy?
                    % bring LST in format [0,360]
LST = mod(LST,360);
target_HA = LST - target_list_RADec(:,2);   % Hour Angle derived from LST and RA of each target


%% Calculation of gondola attitude

RADec_gon = RADec_Star_Tracker - Angle_Enc;


%% Current position of the celestial targets

% target_list_AltAz             % Alt, Az of each target for current
% time and location

l = length(target_list_RADec(:,1));
target_list_AltAz = zeros(l,2);
for i=1:l   % calculate Alt and Az for every celestial target
    [target_list_AltAz(i,1),target_list_AltAz(i,2)] = DecHALat2AltAz(target_list_RADec(i,3), target_HA(i), GPS_pos(1));
end



%% Target selection

% init
exp_prio_list = [1,2,3,4,4,4,3,2,1,1,1];

target_prio = zeros(l,6); % no. exposures, Magnitude, position parameter, exposure parameter, type prioritisation, total prio parameter
target_prio(:,2) = target_list_RADec(:,4); % Magnitude of target
target_prio(:,5) = target_list_RADec(:,5); % prioritisation parameter based on type of target (Nebula, Galaxy, ...)

for i=1:l
    
    % position parameter
    tar_pos_gondola = RADec_gon(1) - target_list_RADec(i,2);      % position relative to gondola (z-axis)
    if abs(tar_pos_gondola) < FoV_gimbal/2  
        target_prio(i,3) = (FoV_gimbal/2 - abs(tar_pos_gondola));
    else
        target_prio(i,3) = 0;
    end
    
    % exposure parameter
    if target_prio(i,1) > 10
        target_prio(i,4) = 0;
    else
        target_prio(i,4) = exp_prio_list(target_prio(i,1)+1);
    end
    
    % total prioritisation parameter
    target_prio(i,6) = target_prio(i,2) * target_prio(i,3) * target_prio(i,4) * target_prio(i,4);
        
end

[prio,index] = max(target_prio(:,6));
target = target_list_RADec(index,:);              % name, RA, Dec, Magnitude of selected target


%% Target tracking

% init
% for testing: just add a value in the order of msec/sec/minutes/... to
% d_J2000, UT to simulate exposure
timestep = 1/3600;          % 1s
T = timestep;               % sample time tracking
roll_ST = RADec_Star_Tracker(3);       % roll starts at initial position of star tracker

exp_time = 1/T;             % equiv. to 1h exposure
yaw = zeros(exp_time,1);
pitch = zeros(exp_time,1);
roll = zeros(exp_time,1);
% do before & during exposure (for loop as placeholder)
for i=1:exp_time

% call current time and date
%RPi_time            % time (format?)
%RPi_date            % date (format?)


% calculate current position of target
d_J2000 = d_J2000 + T;             % format J2000, calculated by GPS_time and GPS_date
UT = UT + T;                  % UT from GPS_time

LST = 100.46 + 0.985647 * d_J2000 + GPS_pos(2) + 15*UT; % Local Sidereal Time -> correct formula, accuracy?
LST = mod(LST,360);                    % bring LST in format [0,360]
                    
target_HA = LST - target(2);   % Hour Angle derived from LST and RA of target

% yaw                 % az
% pitch               % alt
[yaw(i), pitch(i)] = DecHALat2AltAz(target(3), target_HA, GPS_pos(1));

% roll 
roll(i) = roll_ST + timestep*i*360/23.9333;     % approx. 15deg/h

end

time_tr = timestep:timestep:timestep*exp_time;
figure();
plot(time_tr, yaw);
xlabel('time since start of exposure (in h)');
ylabel('Angle (in deg)');
hold on;

plot(time_tr, pitch);
plot(time_tr, roll);
legend('yaw', 'pitch', 'roll');

