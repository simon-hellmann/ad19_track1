%% costWLS.m
% Weighted least-squares cost for parameter identification.
% Thin wrapper around simulateLong; all simulation logic lives there.
%
%   J = sum_i sum_t [ (y_sim_i(t,theta) - y_meas_i(t)) / (sigma_i * sqrt(n_i)) ]^2
%
% where the normalisation by sigma_i and n_i is encoded in scale_long.
% Returns J_penalty when the ODE solver fails so fmincon can backtrack
% instead of propagating NaN into the gradient.
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   J:              scalar total WLS cost
%   J_ch:           (1 x n_out) cost split by output channel
%   r_scaled_cell:  (1 x n_out) cell — scaled residuals per channel
%   y_sim:          (N_long x 1) simulated outputs in long-vector ordering
%
%% Input
%
%   theta:          (n_theta x 1) parameter vector
%   y_meas_long:    (N_long x 1) stacked measured values (time-sorted)
%   t_meas_long:    (N_long x 1) corresponding measurement times        [d]
%   out_idx:        (N_long x 1) output channel index per entry         [1..n_out]
%   scale_long:     (N_long x 1) per-entry scaling  sigma_k * sqrt(n_k)
%   t_events:       (n_ev x 1) feed-event time grid                     [d]
%   u_segments:     (n_ev-1 x 1) feed rate per segment                  [m^3/d]
%   xi_segments:    (n_ev-1 x n_xi) inlet composition per segment
%   x0:             (n_states x 1) initial state vector
%   odeFunc:        @(x, u, xi, theta) ODE right-hand side
%   measFunc:       @(x, theta) output equation
%   odeOpts:        odeset options struct
%   J_penalty:      scalar cost returned on ODE solver failure
%

function [J, J_ch, r_scaled_cell, y_sim] = costWLS(theta, y_meas_long, ...
        t_meas_long, out_idx, scale_long, t_events, u_segments, ...
        xi_segments, x0, odeFunc, measFunc, odeOpts, J_penalty)

y_sim = simulateLong(theta, t_meas_long, out_idx, t_events, ...
                     u_segments, xi_segments, x0, odeFunc, measFunc, odeOpts);

if any(~isfinite(y_sim)) || any(isnan(y_sim))
    J = J_penalty;
    return;
end

r_scaled = (y_sim - y_meas_long) ./ scale_long;
J        = r_scaled' * r_scaled;

% split total cost into output channels:
n_out         = max(out_idx);
J_ch          = nan(1, n_out);
r_scaled_cell = cell(1, n_out);
y_sim_cell    = cell(1, n_out);

for ch_k = 1:n_out
    mask             = out_idx == ch_k;
    r_ch             = r_scaled(mask);
    J_ch(ch_k)       = r_ch' * r_ch;
    r_scaled_cell{ch_k} = r_ch;
    y_sim_cell{ch_k}    = y_sim(mask);
end % for

end % fun
