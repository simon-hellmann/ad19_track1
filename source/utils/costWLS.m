function J = costWLS(theta, y_meas_long, t_meas_long, out_idx, ...
                     scale_long, t_events, u_segments, xi_segments, x0, ...
                     odeFunc, measFunc, odeOpts, J_penalty)
% Weighted least-squares cost function for parameter identification.
% Thin wrapper around simulateLong; all simulation logic lives there.
%
%   J = sum_i (1/n_i) * sum_t [ (y_sim_i(t,theta) - y_meas_i(t)) / sigma_i ]^2
%
% where the normalisation by n_i and sigma_i is encoded in scale_long.
%
% If the ODE solver fails (non-finite output), J_penalty is returned so
% fmincon can backtrack rather than propagating NaN.

    y_sim = simulateLong(theta, t_meas_long, out_idx, t_events, ...
                         u_segments, xi_segments, x0, odeFunc, measFunc, odeOpts);

    % guard for failed cost function evaluations:
    if any(~isfinite(y_sim))
        J = J_penalty;
        return;
    end

    r_scaled = (y_sim - y_meas_long) ./ scale_long;
    J        = r_scaled' * r_scaled;
end
