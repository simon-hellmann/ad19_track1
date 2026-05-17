function MESS_cut = cutMESS(MESS, T_start, T_end)
% Cut a MESS struct to the time window [T_start, T_end] and compute
% relative time vectors in days from T_start for every sub-struct.
%
% Mirrors the logic of cut_MESS_4_obs from the ad_para_ident project,
% adapted for the PSS workflow.
%
% Inputs:
%   MESS     -- full measurement struct (datetime-indexed)
%   T_start  -- window start (datetime scalar)
%   T_end    -- window end   (datetime scalar)
%
% Output:
%   MESS_cut -- struct with same fields, cut to [T_start, T_end] and
%               enriched with .rel_time fields in days

    if T_start >= T_end
        error('cutMESS: T_start must be strictly before T_end.');
    end

    MESS_cut = MESS;

    % --- online gas (gasflow, methane, co2 share one time vector) ---------
    mask = MESS.online_gas.time < T_start | MESS.online_gas.time > T_end;
    MESS_cut.online_gas.time(mask)    = [];
    MESS_cut.online_gas.gasflow(mask) = [];
    MESS_cut.online_gas.methane(mask) = [];
    MESS_cut.online_gas.co2(mask)     = [];
    MESS_cut.online_gas.rel_time_gasflow = ...
        days(MESS_cut.online_gas.time - T_start);

    % --- online liquid (pH) -----------------------------------------------
    mask = MESS.online_liq.time < T_start | MESS.online_liq.time > T_end;
    MESS_cut.online_liq.time(mask) = [];
    MESS_cut.online_liq.pH(mask)   = [];
    MESS_cut.online_liq.rel_time   = days(MESS_cut.online_liq.time - T_start);

    % --- offline acetate --------------------------------------------------
    mask = MESS.offline_ac_eq.time < T_start | MESS.offline_ac_eq.time > T_end;
    MESS_cut.offline_ac_eq.time(mask)  = [];
    MESS_cut.offline_ac_eq.ac_eq(mask) = [];
    MESS_cut.offline_ac_eq.rel_time    = days(MESS_cut.offline_ac_eq.time - T_start);

    % --- offline inorganic nitrogen ---------------------------------------
    mask = MESS.offline_nh3.time < T_start | MESS.offline_nh3.time > T_end;
    MESS_cut.offline_nh3.time(mask) = [];
    MESS_cut.offline_nh3.sin(mask)  = [];
    MESS_cut.offline_nh3.rel_time   = days(MESS_cut.offline_nh3.time - T_start);

    % --- feeding events ---------------------------------------------------
    mask = MESS.feed.time < T_start | MESS.feed.time > T_end;
    MESS_cut.feed.time(mask)        = [];
    MESS_cut.feed.mass(mask)        = [];
    MESS_cut.feed.charac(mask, :)   = [];
    MESS_cut.feed.rel_time          = days(MESS_cut.feed.time - T_start);

    % --- window metadata --------------------------------------------------
    MESS_cut.T_start = T_start;
    MESS_cut.T_end   = T_end;
    MESS_cut.t_span  = days(T_end - T_start);   % total duration [d]
end
