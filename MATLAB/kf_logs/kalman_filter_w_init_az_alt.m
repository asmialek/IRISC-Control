function [angle_est] = kalman_filter_w_init_az_alt(gyro_data_input, dt, init_time, ARW, RRW, gyro_bias_0, ref_angle, ref_angle_noise)
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
angle_est = zeros(2,2,N);
angle_est(:,1,1) = [ref_angle(1);0];
angle_est(:,2,1) = [ref_angle(2);0];

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

nu_kf = zeros(1,2,N);
angle_est_next_kf = zeros(2,2,N);
S_kf = zeros(1,1,2,N);
P = zeros(2,2,2,N);
P1 = zeros(2,2,2,N);
l = 1;
P(:,:,1,1) = diag([ref_angle_noise^2 gyro_bias_0^2]);
P(:,:,2,1) = diag([ref_angle_noise^2 gyro_bias_0^2]);

for k=1:N-1
    if(k/100 - ceil(k/100) ==0)
        k
    end
    % altitude axis
    % grab previous state vector
    x_prev(:,1) = angle_est(:,1,k);
    P_prev(:,:,1) = P(:,:,1,k);
    
    % grab gyro measurement
    w_meas_k(1) = gyro_data_input(3,k); % take z axis for altitude angle
    
    [angle_est_next(:,1), P_next(:,:,1)] = kf_propagate(Phi, x_prev(:,1), P_prev(:,:,1), Gamma, w_meas_k(1), Upsilon, Q, Upsilon2, Q2);
    % propagate state vector
%     angle_est_next(:,1) = Phi * x_prev(:,1) + Gamma * w_meas_k(1);
%     
%     % propagate
%     P_next(:,:,1) = Phi*P_prev(:,:,1)*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
    
    
    angle_est_next_kf(:,1,k) = angle_est_next(:,1);
    P1(:,:,1,k) = P_next(:,:,1);
    
    if l <= init_steps % update only for the first measurements where pointing is known (no motion)
        
        [nu_next(:,1),S_next(:,:,1), angle_upd(:,1), P_upd(:,:,1)] = kf_innovate(ref_angle(1), H, angle_est_next(:,1), P_next(:,:,1), R);
        % innovation
%         nu_next(:,1) = ref_angle(1) - H*angle_est_next(:,1);                       % 'propagation error' (difference between estimated state and measured state)
%         S_next(:,:,1) = H*P_next(:,:,1)*H' + R;
%         
%         % compute the Kalman gain
%         K = P_next(:,:,1)*H'/S_next(:,:,1);
%         
%         % update of states and covariance
%         angle_upd(:,1) = angle_est_next(:,1) + K*nu_next(:,1);
%         P_upd(:,:,1) = (eye(2)-K*H)*P_next(:,:,1);
        
        
        % save estimates
        angle_est(:,1,k+1) = angle_upd(:,1);
        P(:,:,1,k+1) = P_upd(:,:,1);
        nu_kf(:,1,k+1) = nu_next(:,1);
        S_kf(:,:,1,k+1) = S_next(:,:,1);
        % l = l+1;
        
        [P_next(:,:,1)] = kf_prop_upd(P_upd(:,:,1), Phi, Upsilon, Q, Upsilon2, Q2);
        % propagate with initialisation measurement
%         P_prev(:,:,1) = P_upd(:,:,1);
%         P_next(:,:,1) = Phi*P_prev(:,:,1)*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
        
        
    else % initialisation completed - only propagate from this point
        %nu_next = t_meas(k+1) - H*angle_est_next;
        
        % save estimates
        angle_est(:,1,k+1) = angle_est_next(:,1);              % only propagate, don't update
        P(:,:,1,k+1) = P_next(:,:,1);                    % no updated values
        nu_kf(:,1,k+1) = nu_next(:,1);                 % update
        S_kf(:,:,1,k+1) = S_kf(:,:,1,k);            % no updated values
        
    end % if init
    
    
    
    % azimuth
    % grab previous state vector
    x_prev(:,2) = angle_est(:,2,k);
    P_prev(:,:,2) = P(:,:,2,k);
    
    % grab gyro measurement
    % TODO merge x and y axis
    w_meas_k(2) = gyro_data_input(2,k)*sin(angle_est(1,1,k+1)) - gyro_data_input(1,k)*cos(angle_est(1,1,k+1)); % merge x and y axis for projected gyro input
    
    [angle_est_next(:,2), P_next(:,:,2)] = kf_propagate(Phi, x_prev(:,2), P_prev(:,:,2), Gamma, w_meas_k(2), Upsilon, Q, Upsilon2, Q2);
    % propagate state vector
%     angle_est_next(:,2) = Phi * x_prev(:,2) + Gamma * w_meas_k(2);
%     
%     % propagate
%     P_next(:,:,2) = Phi*P_prev(:,:,2)*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
%     
%     angle_est_next_kf(:,2,k) = angle_est_next(:,2);
%     P1(:,:,k,2) = P_next(:,:,2);
    
    if l <= init_steps % update only for the first measurements where pointing is known (no motion)
        
        [nu_next(:,2),S_next(:,:,2), angle_upd(:,2), P_upd(:,:,2)] = kf_innovate(ref_angle(2), H, angle_est_next(:,2), P_next(:,:,2), R);
%         % innovation
%         nu_next(:,2) = ref_angle(2) - H*angle_est_next(:,2);                       % 'propagation error' (difference between estimated state and measured state)
%         S_next(:,:,2) = H*P_next(:,:,2)*H' + R;
%         
%         % compute the Kalman gain
%         K = P_next(:,:,2)*H'/S_next(:,:,2);
%         
%         % update of states and covariance
%         angle_upd(:,2) = angle_est_next(:,2) + K*nu_next(:,2);
%         P_upd(:,:,2) = (eye(2)-K*H)*P_next(:,:,2);
        
        
        % save estimates
        angle_est(:,2,k+1) = angle_upd(:,2);
        P(:,:,2,k+1) = P_upd(:,:,2);
        nu_kf(:,2,k+1) = nu_next(:,2);
        S_kf(:,:,2,k+1) = S_next(:,:,2);
        l = l+1;
        
        [P_next(:,:,2)] = kf_prop_upd(P_upd(:,:,2), Phi, Upsilon, Q, Upsilon2, Q2);
        % propagate with initialisation measurement
%         P_prev(:,:,2) = P_upd(:,:,2);
%         P_next(:,:,2) = Phi*P_prev(:,:,2)*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
        
        
    else % initialisation completed - only propagate from this point
        %nu_next = t_meas(k+1) - H*angle_est_next;
        
        % save estimates
        angle_est(:,2,k+1) = angle_est_next(:,2);              % only propagate, don't update
        P(:,:,2,k+1) = P_next(:,:,2);                    % no updated values
        nu_kf(:,2,k+1) = nu_next(:,2);                 % update
        S_kf(:,:,2,k+1) = S_kf(:,:,2,k);            % no updated values
        
    end % if init
    
    
    
    
end % for loop




end

