%% preprocessData.m
% Remove gas-production measurements (output channel 1) that fall within
% exclusion windows around specified disturbance events, or that are below
% a minimum physically plausible gas production value.
%
% Gas flow is disturbed by feeding and by manual sampling for IN (ch 5) and
% AC (ch 6) concentrations.  Affected samples are deleted from tMeas{1} and
% yMeas{1} so they never enter the WLS cost or the scaling normalisation.
%
% Author: Simon Hellmann. Created: 2026/05/23. Version: Matlab R2022b, Update 6
%
%% Output
%
%   data_out:   data struct — same fields as data_in; tMeas{1} and yMeas{1}
%               have disturbance-affected samples removed
%
%% Input
%
%   data_in:    data struct with fields
%               .tMeas          -- {n_out x 1} cell, measurement times per channel [d]
%               .yMeas          -- {n_out x 1} cell, measured values per channel
%               .t_feed_start   -- (n_ev x 1) feed event start times               [d]
%               (all other fields are passed through unchanged)
%
%   opts:       exclusion options struct with fields
%               .flag_filter_feed   -- exclude gas data around feed events
%               .flag_filter_IN     -- exclude gas data around IN-sample events (ch 5)
%               .flag_filter_AC     -- exclude gas data around AC-sample events (ch 6)
%               .dt_feed_before     -- exclusion window before feed start           [d]
%               .dt_feed_after      -- exclusion window after feed end              [d]
%               .dt_IN_before       -- exclusion window before each IN sample       [d]
%               .dt_IN_after        -- exclusion window after each IN sample        [d]
%               .dt_AC_before       -- exclusion window before each AC sample       [d]
%               .dt_AC_after        -- exclusion window after each AC sample        [d]
%               .q_gas_min          -- minimum valid gas production value [m^3/d]
%                                      samples below this are always excluded
%
%   feeding_duration  -- scalar, feeding event duration [d]; used to compute
%                        t_feed_end = t_feed_start + feeding_duration
%

function data_out = preprocessData(data_in, opts, feeding_duration)

data_out = data_in;   % pass all fields through; overwrite gas channel below

t_feed_end = data_in.t_feed_start + feeding_duration;   % [d], local only

% --- Collect exclusion intervals [t_lo, t_hi] ----------------------------
intervals = zeros(0, 2);   % (n_intervals x 2)

if opts.flag_filter_feed
    n_ev = numel(data_in.t_feed_start);
    for ev_k = 1:n_ev
        t_lo = data_in.t_feed_start(ev_k) - opts.dt_feed_before;
        t_hi = t_feed_end(ev_k)            + opts.dt_feed_after;
        intervals(end+1, :) = [t_lo, t_hi];
    end % for
end

if opts.flag_filter_IN
    t_ev = data_in.tMeas{5};
    for ev_k = 1:numel(t_ev)
        t_lo = t_ev(ev_k) - opts.dt_IN_before;
        t_hi = t_ev(ev_k) + opts.dt_IN_after;
        intervals(end+1, :) = [t_lo, t_hi];
    end % for
end

if opts.flag_filter_AC
    t_ev = data_in.tMeas{6};
    for ev_k = 1:numel(t_ev)
        t_lo = t_ev(ev_k) - opts.dt_AC_before;
        t_hi = t_ev(ev_k) + opts.dt_AC_after;
        intervals(end+1, :) = [t_lo, t_hi];
    end % for
end

% --- Apply exclusion mask to gas channel (output 1) ----------------------
t_gas     = data_in.tMeas{1};
y_gas     = data_in.yMeas{1};
flag_excl = false(numel(t_gas), 1);

% OR-accumulate: a sample is excluded if it falls inside ANY exclusion window.
for intv_k = 1:size(intervals, 1)
    flag_excl = flag_excl | ...
        (t_gas >= intervals(intv_k, 1) & t_gas <= intervals(intv_k, 2));
end % for

% Exclude samples below the minimum gas production threshold (sensor noise,
% reactor not yet running, or post-sample depression artifacts).
flag_excl = flag_excl | (y_gas < opts.q_gas_min);

data_out.tMeas{1} = t_gas(~flag_excl);
data_out.yMeas{1} = y_gas(~flag_excl);

end % fun
