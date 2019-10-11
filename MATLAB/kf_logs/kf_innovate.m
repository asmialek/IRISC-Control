function [nu_next,S_next, angle_upd, P_upd] = kf_innovate(ref_angle, H, angle_est_next, P_next, R)
%kf_innovate innovation of Kalman filter
%   innovation of Kalman filter and update of states and covariance

 % innovation
        nu_next = ref_angle - H*angle_est_next;                       % 'propagation error' (difference between estimated state and measured state)
        S_next = H*P_next*H' + R;
        
        % compute the Kalman gain
        K = P_next*H'/S_next;
        
        % update of states and covariance
        angle_upd = angle_est_next + K*nu_next;
        P_upd = (eye(2)-K*H)*P_next;
end

