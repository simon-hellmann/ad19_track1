%% cutPSS.m
% Cut a PSS data struct to the time window [T_start, T_end) and re-zero
% all relative timestamps to T_start.
%
% Both measurement channels and feed events are trimmed using a half-open
% interval so adjacent windows share no boundary samples.
%
% Author: Simon Hellmann. Created: 2026/05/27. Version: Matlab R2022b, Update 6
%
%% Output
%
%   data_out:   PSS struct — same fields as data_in; tMeas/yMeas and feed
%               arrays trimmed and re-zeroed to T_start
%
%% Input
%
%   data_in:    PSS struct with fields
%               .tMeas          -- {n_out x 1} cell, measurement times  [d]
%               .yMeas          -- {n_out x 1} cell, measured values
%               .t_feed_start   -- (n_ev x 1) feed event start times    [d]
%               .feed_mass      -- (n_ev x 1) mass per feed event       [kg]
%               .xi_feed        -- (n_ev x n_xi) substrate composition
%               .T_start        -- absolute datetime origin of data_in
%               (all other fields are passed through unchanged)
%
%   T_start:    window start (datetime scalar)
%
%   T_end:      window end   (datetime scalar)

function data_out = cutPSS(data_in, T_start, T_end)

if T_start >= T_end
    error('cutPSS: T_start must be strictly before T_end.');
end

% Align timezones so datetime arithmetic is unambiguous
ref_tz = data_in.T_start.TimeZone;
if ~isempty(ref_tz)
    if isempty(T_start.TimeZone)
        T_start = datetime(T_start, 'TimeZone',ref_tz);
    end
    if isempty(T_end.TimeZone)
        T_end = datetime(T_end, 'TimeZone',ref_tz);
    end
end

t_lo = days(T_start - data_in.T_start);   % [d] window start relative to data origin
t_hi = days(T_end   - data_in.T_start);   % [d] window end   relative to data origin

data_out = data_in;

% --- Measurement channels ------------------------------------------------
n_ch = numel(data_in.tMeas);
for i_ch = 1:n_ch
    t_i                   = data_in.tMeas{i_ch};
    mask                  = t_i >= t_lo & t_i < t_hi;
    data_out.tMeas{i_ch}  = t_i(mask) - t_lo;   % re-zero to window start
    data_out.yMeas{i_ch}  = data_in.yMeas{i_ch}(mask);
end % for

% --- Feed events ----------------------------------------------------------
t_f    = data_in.t_feed_start;
mask_f = t_f >= t_lo & t_f < t_hi;
data_out.t_feed_start = t_f(mask_f) - t_lo;   % re-zero to window start
data_out.feed_mass    = data_in.feed_mass(mask_f);
if isfield(data_in, 'xi_feed') && ~isempty(data_in.xi_feed)
    data_out.xi_feed = data_in.xi_feed(mask_f, :);
end

% --- Window metadata ------------------------------------------------------
data_out.T_start = T_start;
data_out.t0      = 0;
data_out.tf      = days(T_end - T_start);   % [d] window duration

end % fun
