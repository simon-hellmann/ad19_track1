%% simulateLong.m
% Piecewise ODE integration across feeding event boundaries, then assemble
% simulated outputs at measurement times and — optionally — on an
% equidistant fine time grid for smooth plotting over the full dataset window.
%
% Author: Simon Hellmann. Created: 2026/05/15. Version: Matlab R2022b, Update 6
%
%% Output
%
%   y_sim_long:     simulated outputs at measurement times (N_long x 1),
%                   same row ordering as t_meas_long / out_idx
%   x_end:          terminal state at t_events(end) (n_states x 1);
%                   use as x0_CV when data_cross follows data_auto in time
%   t_fine:         equidistant fine time grid from t0 to tf (n_fine x 1) [d];
%                   only computed when nargout >= 3
%   x_fine:         state trajectories at t_fine (n_fine x n_states);
%                   only computed when nargout >= 3
%   y_fine:         all model outputs at t_fine (n_fine x n_out);
%                   only computed when nargout >= 3
%
%% Input
%
%   theta:          parameter vector (n_params x 1)
%   t_meas_long:    sorted measurement times (N_long x 1) [d]
%   out_idx:        output index for each long-vector entry (N_long x 1)
%   t_events:       sorted boundary times incl. t0 and tf (n_seg+1 x 1) [d]
%   u_segments:     constant feed flow for each segment (n_seg x 1) [m^3/d]
%   xi_segments:    inlet composition for each segment (n_seg x n_xi)
%   x0:             initial state (n_states x 1)
%   odeFunc:        handle: @(x, u, xi, theta) -> f  (RHS, first output only)
%   measFunc:       handle: @(x, theta) -> g         (n_out x 1)
%   odeOpts:        odeset options struct
%   dt_fine:        fine grid resolution [d] (optional, default: 10/1440 = 10 min)
%

function [y_sim_long, x_end, t_fine, x_fine, y_fine] = simulateLong( ...
    theta, t_meas_long, out_idx, t_events, u_segments, xi_segments, ...
    x0, odeFunc, measFunc, odeOpts, dt_fine)

%% optional argument default

if nargin < 11 || isempty(dt_fine)
    dt_fine = 10/1440;  % [d] = 10 min
end

%% piecewise ODE integration

tSol = [];
xSol = [];
x_k  = x0(:); % enforce column vector

for k = 1:numel(t_events) - 1
    tSeg = [t_events(k), t_events(k+1)];
    [tSeg_out, xSeg_out] = ode15s(@(t,x) odeFunc(x, u_segments(k), xi_segments(k,:)', ...
        theta), tSeg, x_k, odeOpts);
    if k == 1
        tSol = tSeg_out;
        xSol = xSeg_out;
    else
        tSol = [tSol; tSeg_out(2:end)];        %#ok<AGROW>
        xSol = [xSol; xSeg_out(2:end, :)];     %#ok<AGROW>
    end
    x_k = xSeg_out(end, :)';
end % for
x_end = x_k;

%% assemble long output vector at measurement times

t_unique   = unique(t_meas_long);
xAtUnique  = interp1(tSol, xSol, t_unique);    % (n_uniq x n_states)
N_long     = numel(t_meas_long);
y_sim_long = zeros(N_long, 1);

for k = 1:numel(t_unique)
    mask             = (t_meas_long == t_unique(k));
    gVec             = measFunc(xAtUnique(k,:)', theta);   % (n_out x 1)
    y_sim_long(mask) = gVec(out_idx(mask));
end

%% fine equidistant grid (only computed when caller requests it)

if nargout >= 3
    n_fine = ceil((t_events(end) - t_events(1)) / dt_fine) + 1;
    t_fine = linspace(t_events(1), t_events(end), n_fine)';    % (n_fine x 1) [d]
    x_fine = interp1(tSol, xSol, t_fine);                      % (n_fine x n_states)

    g0     = measFunc(x_fine(1,:)', theta);
    n_out  = numel(g0);
    y_fine = zeros(n_fine, n_out);
    for k = 1:n_fine
        y_fine(k,:) = measFunc(x_fine(k,:)', theta)';
    end % for
end % if

end % fun
