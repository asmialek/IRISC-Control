function [P_next] = kf_prop_upd(P_upd, Phi, Upsilon, Q, Upsilon2, Q2)
%kf_prop_upd propagate after update step

% propagate with initialisation measurement
        P_prev = P_upd;
        P_next = Phi*P_prev*Phi' + Upsilon*Q*Upsilon' + Upsilon2*Q2*Upsilon2';
end

