%% computeX0.m
% Dynamic state initialisation in two steps:
%
%   Step 1 -- Steady-state pre-simulation (skipped when t_ss = 0)
%     Run the model for t_ss days under constant average feed rate and
%     volume-weighted average inlet composition.  Drives the state to a
%     representative operating point regardless of x0_rough.
%     If measFunc is provided, a sour-SS guard checks pH output 4 < 6.5
%     after the pre-simulation and falls back to x0_init if triggered.
%
%   Step 2 -- Dynamic swing-up over data_init
%     Simulate piecewise across the data_init window using its actual
%     feeding schedule and per-event inlet compositions, starting from the
%     steady state found in step 1.  The terminal state is returned as x0.
%
% If flag_plot is true, two diagnostic figures are opened:
%   Fig 1 -- all state trajectories during SS pre-simulation (2-col grid,
%            omitted when t_ss = 0)
%   Fig 2 -- model output trajectories during swing-up overlaid with
%            data_init measurements; feed events as red dashed lines
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   x0:              warm-started initial state for the training window (n_states x 1)
%
%% Input
%
%   theta:           current parameter vector                            (n_params x 1)
%   data_init:       init dataset struct (output of prepare_data.m) with fields
%                    .t_feed_start  -- (n_ev x 1) feed start times      [d]
%                    .feed_mass     -- (n_ev x 1) total mass per event   [kg]
%                    .xi_feed       -- (n_ev x n_xi) inlet compositions
%                    .tMeas         -- {n_out x 1} measurement times     [d]
%                    .yMeas         -- {n_out x 1} measured values
%                    .t0, .tf       -- window start / end                [d]
%   t_ss:            SS pre-simulation duration; 0 skips step 1 and    [d]
%                   uses x0_init directly as the swing-up warm start
%   x0_init:         rough initial state (n_states x 1)
%   odeFunc:         function handle  @(x, u, xi, theta) -> dx/dt
%   odeOpts:         odeset options struct
%   feeding_duration: feeding event duration                            [d]
%   flag_plot:       true → open diagnostic figures (default: false)
%   measFunc:        function handle  @(x, theta) -> y
%                    required when flag_plot = true; optional when
%                    t_ss > 0: enables sour-SS guard (pH < 6.5 triggers
%                    fallback to x0_init)

function x0 = computeX0(theta, data_init, t_ss, x0_init, odeFunc, odeOpts, ...
    feeding_duration, flag_plot, measFunc)

if nargin < 8
    flag_plot = false;
end
if nargin < 9
    measFunc = [];
end
if flag_plot && isempty(measFunc)
    error('computeX0: measFunc is required when flag_plot = true.');
end

rho_feed = 1000;   % [kg/m^3] substrate density
feed_volumes = data_init.feed_mass(:) / rho_feed;   % [m^3] per event

%% Step 1: SS pre-simulation with average feed rate

if t_ss > 0
    total_volume = sum(feed_volumes);

    u_avg  = total_volume / data_init.tf;                             % [m^3/d] avg flow
    xi_avg = (feed_volumes' * data_init.xi_feed)' / total_volume;    % (n_xi x 1) vol-weighted

    [t_sol_ss, x_sol_ss] = ode15s(@(t,x) odeFunc(x, u_avg, xi_avg, theta), ...
                                   [0, t_ss], x0_init(:), odeOpts);
    x_ss = x_sol_ss(end,:)';

    % Sour-SS guard: if pH < 6.5 after pre-simulation, fall back to x0_init
    if ~isempty(measFunc)
        y_ss  = measFunc(x_ss, theta);
        if y_ss(4) < 6.5
            warning("computeX0: sour SS detected after pre-simulation " + ...
                    "(pH = %.2f). Falling back to x0_init as warm start.", y_ss(4));
            x_ss = x0_init(:);
        end
    end
else
    % t_ss = 0: skip SS pre-simulation, use provided warm start directly
    t_sol_ss = [];
    x_sol_ss = [];
    x_ss     = x0_init(:);
end

%% Step 2: dynamic swing-up across the data_init window

t_feed_end_init   = data_init.t_feed_start(:) + feeding_duration;   % [d]
u_feed_value_init = feed_volumes / feeding_duration;                 % [m^3/d]

t_events_init  = unique([data_init.t0; ...
                          data_init.t_feed_start(:); ...
                          t_feed_end_init; ...
                          data_init.tf]);
t_mid_init     = (t_events_init(1:end-1) + t_events_init(2:end)) / 2;
n_seg          = numel(t_events_init) - 1;
n_xi           = size(data_init.xi_feed, 2);
u_segment_init  = zeros(n_seg, 1);
xi_segment_init = zeros(n_seg, n_xi);

for i = 1:numel(data_init.t_feed_start)
    active = t_mid_init >= data_init.t_feed_start(i) & ...
             t_mid_init <= t_feed_end_init(i);
    u_segment_init(active)    = u_feed_value_init(i);
    xi_segment_init(active,:) = repmat(data_init.xi_feed(i,:), sum(active), 1);
end % for

if flag_plot
    y_test  = measFunc(x_ss, theta);
    n_out   = numel(y_test);
    t_swing = [];
    y_swing = [];
end

x_k = x_ss;
for k = 1:n_seg
    t_seg = [t_events_init(k), t_events_init(k+1)];
    [t_seg_out, x_seg_out] = ode15s( ...
        @(t,x) odeFunc(x, u_segment_init(k), xi_segment_init(k,:)', theta), ...
        t_seg, x_k, odeOpts);
    if flag_plot
        n_pts   = numel(t_seg_out);
        y_seg   = zeros(n_pts, n_out);
        for j = 1:n_pts
            y_seg(j,:) = measFunc(x_seg_out(j,:)', theta)';
        end % for
        t_swing = [t_swing; t_seg_out];   %#ok<AGROW>
        y_swing = [y_swing; y_seg];       %#ok<AGROW>
    end
    x_k = x_seg_out(end,:)';
end % for

x0 = x_k;

%% Diagnostic plots

if ~flag_plot
    return
end

n_states = numel(x0);
n_cols   = 2;
n_rows_s = ceil(n_states / n_cols);
n_rows_o = ceil(n_out   / n_cols);

% --- Fig 1: SS pre-simulation (states) -----------------------------------
if ~isempty(x_sol_ss)
    figure('Name','computeX0 — SS pre-simulation');
    for i_s = 1:n_states
        subplot(n_rows_s, n_cols, i_s);
        plot(t_sol_ss, x_sol_ss(:, i_s));
        xlabel('t [d]');
        ylabel(sprintf('x_{%d}', i_s));
        grid on
    end % for
    sgtitle('SS pre-simulation');
end

% --- Fig 2: swing-up (outputs + measurements) ----------------------------
figure('Name','computeX0 — dynamic swing-up');
for i_o = 1:n_out
    subplot(n_rows_o, n_cols, i_o);
    plot(t_swing, y_swing(:, i_o), 'b-');
    hold on
    if ~isempty(data_init.tMeas{i_o})
        plot(data_init.tMeas{i_o}, data_init.yMeas{i_o}, 'k.', 'MarkerSize',8);
    end
    xline(data_init.t_feed_start, 'r--');
    hold off
    xlabel('t [d]');
    ylabel(sprintf('y_{%d}', i_o));
    grid on
end % for
sgtitle('Dynamic swing-up (data\_init)');

end % fun
