function y_sim_long = simulateLong(theta, t_meas_long, out_idx, ...
                                    t_events, u_segment, x0, ...
                                    odeFunc, measFunc, odeOpts)
% Piecewise ODE integration across feeding event boundaries, then assemble
% the simulated output long vector at measurement times.
%
% Inputs:
%   theta        -- parameter vector (n_params x 1)
%   t_meas_long  -- sorted measurement times (N_long x 1)
%   out_idx      -- output index for each long-vector entry (N_long x 1)
%   t_events     -- sorted boundary times incl. t0 and tf (n_seg+1 x 1)
%   u_segment    -- constant input value for each segment (n_seg x 1)
%   x0           -- initial state (n_states x 1)
%   odeFunc      -- handle: @(x, u, theta) -> f  (RHS, first output only)
%   measFunc     -- handle: @(x, theta) -> g     (n_out x 1)
%   odeOpts      -- odeset options struct
%
% Output:
%   y_sim_long   -- simulated outputs at measurement times (N_long x 1),
%                   same row ordering as t_meas_long / out_idx

    tSol = [];
    xSol = [];
    x_k  = x0(:);

    for k = 1:numel(t_events) - 1
        tSeg = [t_events(k), t_events(k+1)];
        [tSeg_out, xSeg_out] = ode15s( ...
            @(t,x) odeFunc(x, u_segment(k), theta), tSeg, x_k, odeOpts);
        if k == 1
            tSol = tSeg_out;
            xSol = xSeg_out;
        else
            tSol = [tSol; tSeg_out(2:end)];        %#ok<AGROW>
            xSol = [xSol; xSeg_out(2:end, :)];     %#ok<AGROW>
        end
        x_k = xSeg_out(end, :)';
    end

    t_unique   = unique(t_meas_long);
    xAtUniq    = interp1(tSol, xSol, t_unique);    % (n_uniq x n_states)
    N_long     = numel(t_meas_long);
    y_sim_long = zeros(N_long, 1);

    for k = 1:numel(t_unique)
        mask             = (t_meas_long == t_unique(k));
        gVec             = measFunc(xAtUniq(k,:)', theta);   % (n_out x 1)
        y_sim_long(mask) = gVec(out_idx(mask));
    end
end
