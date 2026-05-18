function x0 = computeX0(theta, data_init, t_ss, x0_rough, odeFunc, odeOpts)
% Dynamic state initialisation in two steps:
%
%   Step 1 -- Steady-state pre-simulation
%     Run the model for t_ss days under constant average feed rate and
%     volume-weighted average inlet composition.  This drives the state to
%     a representative operating point regardless of x0_rough.
%
%   Step 2 -- Dynamic swing-up over data_init
%     Simulate piecewise across the data_init window using its actual
%     feeding schedule and per-event inlet compositions, starting from the
%     steady state found in step 1.  The terminal state is returned as x0.
%
% Inputs:
%   theta      -- current parameter vector (n_params x 1)
%   data_init  -- init dataset struct from split_data_pss (output of mess2pssData)
%   t_ss       -- steady-state pre-simulation duration [d], e.g. 500
%   x0_rough   -- rough initial state (n_states x 1), e.g. from literature
%   odeFunc    -- handle: @(x, u, xi, theta) -> f
%   odeOpts    -- odeset options struct
%
% Output:
%   x0         -- warm-started initial state for the training window (n_states x 1)

    % ------------------------------------------------------------------
    %  Step 1: steady-state pre-simulation with average feed rate
    % ------------------------------------------------------------------

    % Volume delivered by each feeding event [m^3]:
    feed_volumes   = data_init.u_feed_value(:) .* ...
                     (data_init.t_feed_end(:) - data_init.t_feed_start(:));
    total_volume   = sum(feed_volumes);

    % Average volumetric feed flow [m^3/d]:
    u_avg  = total_volume / data_init.tf;

    % Volume-weighted average inlet composition (n_xi x 1):
    xi_avg = (feed_volumes' * data_init.xi_feed)' / total_volume;

    [~, xSol_ss] = ode15s(@(t,x) odeFunc(x, u_avg, xi_avg, theta), ...
                           [0, t_ss], x0_rough(:), odeOpts);
    x_ss = xSol_ss(end,:)';

    % ------------------------------------------------------------------
    %  Step 2: dynamic swing-up across the data_init window
    % ------------------------------------------------------------------

    % Build piecewise feeding grid for data_init
    t_events_init  = unique([data_init.t0; ...
                              data_init.t_feed_start(:); ...
                              data_init.t_feed_end(:); ...
                              data_init.tf]);
    t_mid_init     = (t_events_init(1:end-1) + t_events_init(2:end)) / 2;
    n_seg          = numel(t_events_init) - 1;
    n_xi           = size(data_init.xi_feed, 2);
    u_segment_init  = zeros(n_seg, 1);
    xi_segment_init = zeros(n_seg, n_xi);

    for i = 1:numel(data_init.t_feed_start)
        active = t_mid_init >= data_init.t_feed_start(i) & ...
                 t_mid_init <= data_init.t_feed_end(i);
        u_segment_init(active)    = data_init.u_feed_value(i);
        xi_segment_init(active,:) = repmat(data_init.xi_feed(i,:), sum(active), 1);
    end

    % Piecewise integration across data_init
    x_k = x_ss;
    for k = 1:n_seg
        tSeg = [t_events_init(k), t_events_init(k+1)];
        [~, xSeg_out] = ode15s( ...
            @(t,x) odeFunc(x, u_segment_init(k), xi_segment_init(k,:)', theta), ...
            tSeg, x_k, odeOpts);
        x_k = xSeg_out(end,:)';
    end

    x0 = x_k;   % terminal state of the swing-up
end
