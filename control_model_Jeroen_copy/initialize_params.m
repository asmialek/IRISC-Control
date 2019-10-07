% routine to initialize all parameters

% simulation parameters
dt_sim = .01;
t_end = 1000;

% plant parameters
I_el = 5;  % inertia [kgm2]

% motor parameters
Tmax = .1;  % [Nm]

% controller parameters
wn = .003*2*pi;        % desired bandwidth [rad/sec]
zeta = 1;           % desired damping ratio [-]

Kp = I_el*wn^2;
Kd = I_el*2*zeta*wn;

T = 3/(zeta*wn);
Ki = I_el*(wn^2/T);

dt_control = 1/50;

% star tracker parameters
dt_str = 5;
f_str = 1/dt_str;
PSD_str = (.1/180*pi)^2*f_str;

% gyro parameters
f_gyr = 50;
dt_gyr = 1/f_gyr;
PSD_gyr = (.15*pi/180/60)^2;    % arw [rad/rt-sec]
gyro_bias_0 = 10/180*pi/3600;      % gyro bias on/off repeatability [rad/s]
gyr_delay = 2*dt_gyr;

