% routine to illustrate Kalman filter

clear all
close all
clc
rng(17); % set seed for random numbers to achieve the same results for different runs - only interesting for trying out e.g. different sample times. Disable if different noise pattern are needed!

% simulation parameters
dt = 0.1;                                                 % sample time of the system, e.g. 1/gyro_freq if gyro is sensor with highest sample rate
t_end = 100;
%% 
time = [0:dt:t_end];
%time = 0;
N = length(time);

% gyro parameters
ARW = (.15/180*pi/60)^2;                                % angular random walk acc. to data sheet [deg/rt(h)]
gyro_bias_0 = 10/180*pi/3600;                           % initial gyro bias [deg/h]
RRW = (1/180*pi/3600/sqrt(3600))^2;

% star tracker parameters
s_init = 0.1/180*pi;
init_time = 30;                                         % initialisation time of the star tracker (in s) 
init_steps = init_time/dt;                          % number of initialisation steps
init_steps = floor(init_steps);

% generate reference profile
w_real = 1/1800*pi*ones(1,N);                            % angular rotational speed profile [deg/s]
%w_real = 5/1800*pi*cos(time*sqrt(9.81/100));
w_real = [zeros(1, init_steps), w_real];
w_real = w_real(:,1:N);                                 % include motionless initialisation state in the beginning
w_real(:,ceil(3*N/4):N) = 0;
w_real(:,1:N) = 0;
t_real = cumsum(dt*w_real);                             % angle profile [deg]
t_real = cumsum(dt*w_real) + 1;

% generate gyro measurements
gyro_noise = randn(1,N)*sqrt(ARW)/sqrt(dt);             % gaussian gyro noise with standard deviation of ARW                              
gyro_bias = randn(1,1)*gyro_bias_0*ones(1,N);           % models total gyro bias (random constant gyro bias)
w_meas = w_real + gyro_noise + gyro_bias;
w_meas(:,1:N) = 0.001;

% generate star tracker measurements
%str_noise = randn(1,N)*s_init;
str_noise = 0;
t_meas = t_real + str_noise;

% run the kalman filter
x_est = zeros(2,N);
x_est(:,1) = [t_meas(:,1);0];

% propagation matrices
Phi = [1 -dt;0 1];
Gamma = [dt;0];
Upsilon = [dt;0];                           % for gyro noise (angular rate)
Upsilon2 = [0;dt];                          % for gyro bias (angle)

% measurement matrix
H = [1 0];

% tuning values
Q = [ARW/dt];
Q2 = [RRW/dt];
R = [s_init^2];

nu_kf = zeros(1,N);
x_est_next_kf = zeros(2,N);
S_kf = zeros(1,1,N);
P = zeros(2,2,N);
P1 = zeros(2,2,N);
l = 1;
P(:,:,1) = diag([s_init^2 gyro_bias_0^2]);

for k=1:N-1
    % grab previous state vector
    x_prev = x_est(:,k);
    P_prev = P(:,:,k);
    
    % grab gyro measurement
    w_meas_k = w_meas(:,k);
    
    % propagate state vector
    x_est_next = Phi * x_prev + Gamma * w_meas_k;
    x_est_next_kf(:,k) = x_est_next;
    
    % propagate
    P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
    P1(:,:,k) = P_next;
    
    if l <= init_steps % update only for the first measurements where pointing is known (no motion)
        
        % innovation
        nu_next = t_meas(k+1) - H*x_est_next;                       % 'propagation error' (difference between estimated state and measured state)
        S_next = H*P_next*H' + R;
        
        % compute the Kalman gain
        K = P_next*H'/S_next;
        
        % update of states and covariance
        x_upd = x_est_next + K*nu_next;
        P_upd = (eye(2)-K*H)*P_next;
        
        % save estimates
        x_est(:,k+1) = x_upd;
        P(:,:,k+1) = P_upd;
        nu_kf(:,k+1) = nu_next;
        S_kf(:,:,k+1) = S_next;
        l = l+1;
        
        % propagate with initialisation measurement
        P_prev = P_upd;
        P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
        
        
    else % initialisation completed - only propagate from this point
        nu_next = t_meas(k+1) - H*x_est_next;
        
        % save estimates
        x_est(:,k+1) = x_est_next;              % only propagate, don't update
        P(:,:,k+1) = P_next;                    % no updated values
        nu_kf(:,k+1) = nu_next;                 % update
        S_kf(:,:,k+1) = S_kf(:,:,k);            % no updated values
        
        k
        x_est_next
        P_next
    
    end
    
    
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
plot(time,x_est(1,:)*180/pi)
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
plot(time,(P11.^.5)*180/pi,'k')
plot(time,-(P11.^.5)*180/pi,'k')
xlabel('Time [sec]')
ylabel('est. error \theta [^\circ]')
subplot(212)
plot(time,(x_est(2,:)-gyro_bias)*180/pi*3600)
hold all
plot(time,(P22.^.5)*180/pi*3600,'k')
plot(time,-(P22.^.5)*180/pi*3600,'k')
xlabel('Time [sec]')
ylabel('est. error bias [^\circ/h]')


