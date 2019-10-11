function [angle_est_next, P_next] = kf_propagate(Phi, x_prev, P_prev, Gamma, w_meas_k, Upsilon, Q, Upsilon2, Q2)
%kf_propagate propagation Kalman filter
%   propagation of estimated state and P matrix

    % propagate state vector
    angle_est_next = Phi * x_prev + Gamma * w_meas_k;
    
    % propagate
    P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';

end

