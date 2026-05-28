%% costWLS_sub.m
% WLS cost for PI #2: embeds the identifiable parameter subset back into
% the full parameter vector and delegates to costWLS.
% Non-identifiable parameters are held fixed at thetaFixed.
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   J:              scalar WLS cost
%
%% Input
%
%   theta_sub:      (n_sub x 1) identifiable parameter subset in physical units
%   keep_idx:       index vector selecting identifiable entries in the full vector
%   thetaFixed:     (n_theta x 1) full parameter vector; non-identifiable entries fixed
%   y_meas_long:    (N_long x 1) stacked measured values (time-sorted)
%   t_meas_long:    (N_long x 1) corresponding measurement times            [d]
%   out_idx:        (N_long x 1) output channel index per entry             [1..n_out]
%   scale_long:     (N_long x 1) per-entry scaling  sigma_k * sqrt(n_k)
%   t_events:       (n_ev x 1) feed-event time grid                         [d]
%   u_segments:     (n_ev-1 x 1) feed rate per segment                      [m^3/d]
%   xi_segments:    (n_ev-1 x n_xi) inlet composition per segment
%   x0:             (n_states x 1) initial state vector
%   odeFunc:        @(x, u, xi, theta) ODE right-hand side
%   measFunc:       @(x, theta) output equation
%   odeOpts:        odeset options struct
%   J_penalty:      scalar cost returned on ODE solver failure
%

function J = costWLS_sub(theta_sub, keep_idx, thetaFixed, ...
                         y_meas_long, t_meas_long, out_idx, scale_long, ...
                         t_events, u_segments, xi_segments, x0, ...
                         odeFunc, measFunc, odeOpts, J_penalty)

theta_full           = thetaFixed;
theta_full(keep_idx) = theta_sub;
J = costWLS(theta_full, y_meas_long, t_meas_long, out_idx, scale_long, ...
            t_events, u_segments, xi_segments, x0, odeFunc, measFunc, odeOpts, J_penalty);

end % fun
