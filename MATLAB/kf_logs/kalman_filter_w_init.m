function [angle_est] = kalman_filter_w_init(gyro_data_input, dt, init_time, ARW, RRW, gyro_bias_0, ref_angle, ref_angle_noise)
%kalman_filter_w_init kalman filter with initialisation by motionless state
%   Kalman filter for gyroscope with sampling time dt assuming a motionless
%   state for init_time at the beginning, pointing at ref_angle. Noise
%   parameters for the gyroscope are ARW, RRW, gyro_bias_0.
%   gyro data: input in degrees/s
%   ARW: angle random walk, input in deg/root(h) (gyro noise)
%   RRW: rate random walk, input in  (bias drift)
%   gyro_bias_0: initial gyro bias drift in deg/h
%   ref_angle_noise: estimated noise of the reference angle during
%   initialisation in degrees
%   initialisation (ideally 0, tuning parameter)

N = length(gyro_data_input);
%t_end = N*dt;
%time = [0:dt:t_end];

ARW = (ARW/180*pi/60)^2;
RRW = (RRW/180*pi/3600/sqrt(3600))^2;

ref_angle_noise = ref_angle_noise*pi/180;

init_steps = init_time/dt;                          % number of initialisation steps
init_steps = floor(init_steps);


% run the kalman filter
angle_est = zeros(2,N);
angle_est(:,1) = [ref_angle;0];

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
R = [ref_angle_noise^2];

nu_kf = zeros(1,N);
angle_est_next_kf = zeros(2,N);
S_kf = zeros(1,1,N);
P = zeros(2,2,N);
P1 = zeros(2,2,N);
l = 1;
P(:,:,1) = diag([ref_angle_noise^2 gyro_bias_0^2]);

for k=1:N-1
    % grab previous state vector
    x_prev = angle_est(:,k);
    P_prev = P(:,:,k);
    
    % grab gyro measurement
    w_meas_k = gyro_data_input(:,k);
    
    % propagate state vector
    angle_est_next = Phi * x_prev + Gamma * w_meas_k;
    angle_est_next_kf(:,k) = angle_est_next;
    
    % propagate
    P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
    P1(:,:,k) = P_next;
    
    if l <= init_steps % update only for the first measurements where pointing is known (no motion)
        
        % innovation
        nu_next = ref_angle - H*angle_est_next;                       % 'propagation error' (difference between estimated state and measured state)
        S_next = H*P_next*H' + R;
        
        % compute the Kalman gain
        K = P_next*H'/S_next;
        
        % update of states and covariance
        angle_upd = angle_est_next + K*nu_next;
        P_upd = (eye(2)-K*H)*P_next;
        
        % save estimates
        angle_est(:,k+1) = angle_upd;
        P(:,:,k+1) = P_upd;
        nu_kf(:,k+1) = nu_next;
        S_kf(:,:,k+1) = S_next;
        l = l+1;
        
        % propagate with initialisation measurement
        P_prev = P_upd;
        P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
        
        
    else % initialisation completed - only propagate from this point
        %nu_next = t_meas(k+1) - H*angle_est_next;
        
        % save estimates
        angle_est(:,k+1) = angle_est_next;              % only propagate, don't update
        P(:,:,k+1) = P_next;                    % no updated values
        nu_kf(:,k+1) = nu_next;                 % update
        S_kf(:,:,k+1) = S_kf(:,:,k);            % no updated values
    
    end
    
    
end




end

