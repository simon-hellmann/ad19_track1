function J = costWLS_sub(theta_sub, keep_idx, thetaFixed, ...
                         y_meas_long, t_meas_long, out_idx, scale_long, ...
                         t_events, u_segment, x0, odeFunc, measFunc, odeOpts)
% WLS cost for PI #2: optimises only the identifiable parameter subset.
% Non-identifiable parameters are held fixed at thetaFixed.
% Embeds theta_sub back into the full parameter vector, then calls costWLS.

    theta_full           = thetaFixed;
    theta_full(keep_idx) = theta_sub;
    J = costWLS(theta_full, y_meas_long, t_meas_long, out_idx, scale_long, ...
                t_events, u_segment, x0, odeFunc, measFunc, odeOpts);
end
