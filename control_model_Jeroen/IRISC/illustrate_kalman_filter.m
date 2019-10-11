% routine to illustrate Kalman filter

clear all
close all
clc

% simulation parameters
dt = 30;
t_end = 3600;
time = [0:dt:t_end];
N = length(time);

% gyro parameters
ARW = (.15/180*pi/60)^2;                                % angular random walk acc. to data sheet [deg/rt(h)]
gyro_bias_0 = 10/180*pi/3600;                           % initial gyro bias [deg/h]
RRW = (1/180*pi/3600/sqrt(3600))^2;

% star tracker parameters
dt_str = 5;
s_str = 0.1/180*pi;

% generate reference profile
w_real = zeros(1,N);
t_real = zeros(1,N);

% generate gyro measurements
gyro_noise = randn(1,N)*sqrt(ARW)/sqrt(dt);             % gaussian gyro noise with standard deviation of ARW
bias_step = gyro_bias_0*ones(1,N);
%bias_step(1:300) = 0;                                   % just for illustration
bias_drift = cumsum(randn(1,N)*sqrt(RRW)*sqrt(dt));     % models time-variable random bias drift 
gyro_bias = randn(1,1)*gyro_bias_0*ones(1,N) + 0*bias_step + bias_drift;  % models total gyro bias (random constant gyro bias + time-variable random bias drift + whatever other gyro errors one can think of)
w_meas = w_real + gyro_noise + gyro_bias;

% generate star tracker measurements
str_noise = randn(1,N)*s_str;
t_meas = t_real + str_noise;

% run the kalman filter
x_est = zeros(2,N);
x_est(:,1) = [t_meas(:,1);0];

% propagation matrices
Phi = [1 -dt;0 1];
Gamma = [dt;0];
Upsilon = [dt;0];
Upsilon2 = [0;dt];

% measurement matrix
H = [1 0];

% tuning values
Q = [ARW];
Q2 = [RRW];
R = [s_str^2];

nu_kf = zeros(1,N);
S_kf = zeros(1,1,N);
P = zeros(2,2,N);
P(:,:,1) = diag([s_str^2 gyro_bias_0^2]);
for k=1:N-1
    % grab previous state vector
    x_prev = x_est(:,k);
    P_prev = P(:,:,k);
    
    % grab gyro measurement
    w_meas_k = w_meas(:,k);
    
    % propagate state vector
    x_est_next = Phi * x_prev + Gamma * w_meas_k;
    
    % propagate
    P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
    
    
    % innovation
    nu_next = t_meas(k+1) - H*x_est_next;
    S_next = H*P_next*H' + R;
    
    % compute the Kalman gain
    K = P_next*H'/S_next;
    
    % update of states and covariance
    x_upd = x_est_next + K*nu_next;
    P_upd = (eye(2)-K*H)*P_next;
    
    % include bias step
%     if k==300
%         P_upd(2,2) = P_upd(2,2) + gyro_bias_0^2;
%     end
    
    % save estimates
    x_est(:,k+1) = x_upd;
    P(:,:,k+1) = P_upd;
    nu_kf(:,k+1) = nu_next;
    S_kf(:,:,k+1) = S_next;
    
end
    
P11 = P(1,1,:);
P11 = P11(:);

P22 = P(2,2,:);
P22 = P22(:);

S11 = S_kf(1,1,:);
S11 = S11(:);

% plot results
figure
subplot(211)
plot(time,t_real*180/pi)
hold all
plot(time,t_meas*180/pi)
xlabel('Time [sec]')
ylabel('\theta [^\circ]')
subplot(212)
plot(time,w_real*180/pi*3600)
hold all
plot(time,w_meas*180/pi*3600)
xlabel('Time [sec]')
ylabel('\omega [^\circ/h]')

figure
subplot(211)
plot(time,(x_est(1,:)-t_real)*180/pi)
hold all
plot(time,(P11.^.5)*180/pi,'k--')
plot(time,-(P11.^.5)*180/pi,'k--')
xlabel('Time [sec]')
ylabel('est. error \theta [^\circ]')
subplot(212)
plot(time,(x_est(2,:)-gyro_bias)*180/pi*3600)
hold all
plot(time,(P22.^.5)*180/pi*3600,'k--')
plot(time,-(P22.^.5)*180/pi*3600,'k--')
xlabel('Time [sec]')
ylabel('est. error bias [^\circ/h]')

figure
plot(time,nu_kf(1,:)*180/pi)
hold all
plot(time,(S11.^.5)*180/pi,'k--')
plot(time,-(S11.^.5)*180/pi,'k--')
xlabel('Time [sec]')
ylabel('innovation [^\circ]')

