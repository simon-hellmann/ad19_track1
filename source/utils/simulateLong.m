function [y_sim_long, x_end] = simulateLong(theta, t_meas_long, out_idx, ...
                                             t_events, u_segments, xi_segments, x0, ...
                                             odeFunc, measFunc, odeOpts)
% Piecewise ODE integration across feeding event boundaries, then assemble
% the simulated output long vector at measurement times.
%
% Inputs:
%   theta        -- parameter vector (n_params x 1)
%   t_meas_long  -- sorted measurement times (N_long x 1)
%   out_idx      -- output index for each long-vector entry (N_long x 1)
%   t_events     -- sorted boundary times incl. t0 and tf (n_seg+1 x 1)
%   u_segments   -- constant feed flow for each segment (n_seg x 1)
%   xi_segments  -- inlet composition for each segment (n_seg x n_xi)
%   x0           -- initial state (n_states x 1)
%   odeFunc      -- handle: @(x, u, xi, theta) -> f  (RHS, first output only)
%   measFunc     -- handle: @(x, theta) -> g         (n_out x 1)
%   odeOpts      -- odeset options struct
%
% Outputs:
%   y_sim_long   -- simulated outputs at measurement times (N_long x 1),
%                   same row ordering as t_meas_long / out_idx
%   x_end        -- terminal state at t_events(end) (n_states x 1);
%                   use as x0_CV when data_cross follows data_auto in time

    tSol = [];
    xSol = [];
    x_k  = x0(:); % enforce column vector

    % get state trajectories by iterating over all feeding events:
    for k = 1:numel(t_events) - 1
        tSeg = [t_events(k), t_events(k+1)];
        [tSeg_out, xSeg_out] = ode15s( ...
            @(t,x) odeFunc(x, u_segments(k), xi_segments(k,:)', theta), tSeg, x_k, odeOpts);
        if k == 1
            tSol = tSeg_out;
            xSol = xSeg_out;
        else
            tSol = [tSol; tSeg_out(2:end)];        %#ok<AGROW>
            xSol = [xSol; xSeg_out(2:end, :)];     %#ok<AGROW>
        end
        x_k = xSeg_out(end, :)'; % update for next iteration
    end
    x_end = x_k;
    
    % interpolate to measurement time grid: 
    t_unique   = unique(t_meas_long);
    xAtUnique  = interp1(tSol, xSol, t_unique);    % (n_uniq x n_states)
    N_long     = numel(t_meas_long);
    y_sim_long = zeros(N_long, 1);

    % get output trajectories:
    for k = 1:numel(t_unique)
        mask             = (t_meas_long == t_unique(k));
        gVec             = measFunc(xAtUnique(k,:)', theta);   % (n_out x 1)
        y_sim_long(mask) = gVec(out_idx(mask));
    end
end
