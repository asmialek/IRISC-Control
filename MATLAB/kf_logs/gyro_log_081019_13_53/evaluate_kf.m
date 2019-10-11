% results Kalman filter test with real gyroscope data
% gyro frequency: 100Hz
% initialisation time: 300s
% total run time: 1000s

%% Import data from log files
clear all;
close all;
clc;

%% gyro.log
opts = detectImportOptions('gyro.log');
opts = setvartype(opts, 1, 'datetime');
opts = setvartype(opts, 2:4, 'double');
opts = setvaropts(opts, 1, 'InputFormat', 'HH:mm:ss.SSS');
opts = setvaropts(opts, 2:4, 'DecimalSeparator', '.');

[time, gyro_data_x, gyro_data_y, gyro_data_z] = readvars('gyro.log', opts);

%% nu_next.log
opts = detectImportOptions('nu_next.log');
opts = setvartype(opts, 2, 'double');
opts = setvaropts(opts, 2, 'DecimalSeparator', '.');

[~,nu] = readvars('nu_next.log', opts);

%% p_upd.log
opts = detectImportOptions('p_upd.log');
opts = setvartype(opts, 2:5, 'double');
opts = setvaropts(opts, 2:5, 'DecimalSeparator', '.');

[~,P11, P12, P21, P22] = readvars('p_upd.log', opts);

%% s_next.log
opts = detectImportOptions('s_next.log');
opts = setvartype(opts, 2, 'double');
opts = setvaropts(opts, 2, 'DecimalSeparator', '.');

[~,S] = readvars('s_next.log', opts);

%% x_upd.log
opts = detectImportOptions('x_upd.log');
opts = setvartype(opts, 2:3, 'double');
opts = setvaropts(opts, 2:3, 'DecimalSeparator', '.');

[~,z_est, beta_est] = readvars('x_upd.log', opts);

clear opts;
save('log_files.mat');

%% Adjust imported data
% total number of measurements
N = length(time);

% lengthen nu and S to fit length(time)
vec = zeros(N,1);
nu = [nu; vec];
nu = nu(1:N);

vec = zeros(N,1);
S = [S; vec];
S = S(1:N);

clear vec;

% shorten variables to the same length (and correct gyro data log error - no gyro data for the first 1.04s (104 steps))
gyro_error_steps = 96;
gyro_error = zeros(gyro_error_steps, 1);

gyro_data_x = [gyro_error; gyro_data_x(1:(N-gyro_error_steps),1)];
gyro_data_y = [gyro_error; gyro_data_y(1:(N-gyro_error_steps),1)];
gyro_data_z = [gyro_error; gyro_data_z(1:(N-gyro_error_steps),1)];

gyro_data_x = gyro_data_x(1:N,1);
gyro_data_y = gyro_data_y(1:N,1);
gyro_data_z = gyro_data_z(1:N,1);
z_est = z_est(1:N,1) - 1;
beta_est = beta_est(1:N,1);
P11 = P11(1:N,1);
P12 = P12(1:N,1);
P21 = P21(1:N,1);
P22 = P22(1:N,1);

% re-format P
P = zeros(N,2,2);
P(:,1,1) = P11;
P(:,1,2) = P12;
P(:,2,1) = P21;
P(:,2,2) = P22;
P = P(1:N,:,:);

%% Run Kalman filter for imported data
dt = 1/100; 
init_time = 300;
ref_angle = 0;
gyro_bias_0 = 10;

ARW = 0.15;
RRW = 1;
gyro_bias_0 = gyro_bias_0/180*pi/3600;
ref_angle_noise = 0.1;

x_est_kf = kalman_filter_w_init(gyro_data_x', dt, init_time, ARW, RRW, gyro_bias_0, ref_angle, ref_angle_noise);
y_est_kf = kalman_filter_w_init(gyro_data_y', dt, init_time, ARW, RRW, gyro_bias_0, ref_angle, ref_angle_noise);
z_est_kf = kalman_filter_w_init(gyro_data_z', dt, init_time, ARW, RRW, gyro_bias_0, ref_angle, ref_angle_noise);

x_est_kf = x_est_kf';
y_est_kf = y_est_kf';
z_est_kf = z_est_kf';

x_est_no_kf = cumsum(gyro_data_x*dt);
y_est_no_kf = cumsum(gyro_data_y*dt);
z_est_no_kf = cumsum(gyro_data_z*dt);

save('log_files_r2plot');


%% Display results
figure()

plot(time, gyro_data_x);
hold on
plot(time, gyro_data_y);
plot(time, gyro_data_z);
legend('gyro data x', 'gyro data y', 'gyro data z');
xlabel ('time');
ylabel('angular rate [^\circ/s]')
grid on

figure()
hold on
plot(time, x_est_no_kf);
plot(time, x_est_kf(:,1));
legend('estimated angle x-axis without Kalman filter', 'estimated angle x-axis with post-processed Kalman filter');
xlabel ('time');
ylabel('estimated angle [^\circ]')
grid on

figure()
hold on
plot(time, y_est_no_kf);
plot(time, y_est_kf(:,1));
legend('estimated angle y-axis without Kalman filter', 'estimated angle y-axis with post-processed Kalman filter');
xlabel ('time');
ylabel('estimated angle [^\circ]')
grid on

figure()
plot(time, z_est);
hold on
plot(time, z_est_no_kf);
plot(time, z_est_kf(:,1));
plot(time, z_est - z_est_kf(:,1));
plot(time,(P11.^.5)*180/pi,'k')
plot(time,-(P11.^.5)*180/pi,'k')
legend('estimated angle z-axis  with Kalman filter', 'estimated angle z-axis without Kalman filter', 'estimated angle z-axis with post-processed Kalman filter', 'difference between OBSW and post-processed Kalman filter', '', '');
xlabel ('time');
ylabel('estimated angle [^\circ]')
grid on


figure()
plot(time, z_est);
hold on
plot(time, x_est_no_kf);
plot(time, y_est_no_kf);
plot(time, z_est_no_kf);
legend('estimated angle z-axis  with Kalman filter', 'estimated angle x-axis without Kalman filter', 'estimated angle y-axis without Kalman filter', 'estimated angle z-axis without Kalman filter');
xlabel ('time');
ylabel('estimated angle [^\circ]')
grid on





