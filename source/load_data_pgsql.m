%% load_data_pgsql.m
% Load measurement data from the PostgreSQL database and assemble the full
% PSS dataset (data_full) in the format expected by main_pss_workflow.m.
%
% Fills tMeas / yMeas for all 6 output channels and the feed-event fields.
% Currently implemented: channels 1-6 and feeding events from SQL tables.
%
% Author: Simon Hellmann. Created: 2026/05/24. Version: Matlab R2022b, Update 6
%
% no inputs and outputs (main file)

clc
clear
close all

% change working directory to the script's directory
[scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
cd(scriptPath);

addpath('utils');

%% User settings

full_path_raw = '../data/raw/automated_feeder/data_raw_auto_feeder.mat';

reactor_id = 1;

% absolute time window of the full dataset (adjust to match the DB records)
T_start = datetime('29-Apr-2026 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
T_end   = datetime('27-May-2026 12:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');

% feeding event settings
delta_bundle_min = 15;    % max gap for bundling nearby events [min]
rho_water        = 1000;  % [kg/m³] density of water feeds; used for xi volume-weighting

% path to the JSON database config file (relative to this script)
db_config_file = 'config/db_config.json';

%% Connect to PostgreSQL

db = jsondecode(fileread(db_config_file)).database;

fprintf("Connecting to '%s' on %s...\n", db.database, db.host);
try
    conn = postgresql(db.user, db.password, 'Server',db.host, ...
        'DatabaseName',db.database, 'PortNumber',db.port);
catch ME
    error("Database connection failed: %s", ME.message);
end
fprintf("Connected.\n");

%% Query online measurements

% Include UTC offset in the SQL strings so PostgreSQL treats the bounds as
% TIMESTAMPTZ literals and compares them in the correct timezone (no DST
% ambiguity).  'Z' in MATLAB datetime format gives the offset, e.g. +0200.
T_start_str = string(T_start, "yyyy-MM-dd HH:mm:ss Z");
T_end_str   = string(T_end,   "yyyy-MM-dd HH:mm:ss Z");

sql_online = sprintf( ...
    "SELECT time, v_dot_gas, pch4, pco2, ph " + ...
    "FROM public.online_measurements " + ...
    "WHERE reactor_id = %d " + ...
    "  AND time >= '%s' " + ...
    "  AND time <= '%s' " + ...
    "ORDER BY time;", ...
    reactor_id, T_start_str, T_end_str);

fprintf("Querying online measurements...\n");
tbl_online = fetch(conn, sql_online);
fprintf("  %d rows returned.\n", height(tbl_online));

%% Build tMeas / yMeas from online measurements (with NaN handling)

% all online channels share the same time column; NULLs come back as NaN.
% MATLAB returns TIMESTAMPTZ columns as timezone-aware datetime objects;
% subtracting T_start (also timezone-aware) gives correct relative times
% across DST transitions.
t_abs_online = tbl_online.time;                  % datetime (TIMESTAMPTZ, timezone-aware, in UTC)
t_rel_online = days(t_abs_online - T_start);     % [d] relative to T_start

data_full.tMeas = cell(6, 1);
data_full.yMeas = cell(6, 1);

% channel 1: v_dot_gas — gasflow [m³/d]
mask_gas           = ~isnan(tbl_online.v_dot_gas);
data_full.tMeas{1} = t_rel_online(mask_gas);
data_full.yMeas{1} = tbl_online.v_dot_gas(mask_gas);

% channel 2: pch4 — partial pressure CH4 [bar]
mask_ch4           = ~isnan(tbl_online.pch4);
data_full.tMeas{2} = t_rel_online(mask_ch4);
data_full.yMeas{2} = tbl_online.pch4(mask_ch4);

% channel 3: pco2 — partial pressure CO2 [bar]
mask_co2           = ~isnan(tbl_online.pco2);
data_full.tMeas{3} = t_rel_online(mask_co2);
data_full.yMeas{3} = tbl_online.pco2(mask_co2);

% channel 4: ph [-]
mask_ph            = ~isnan(tbl_online.ph);
data_full.tMeas{4} = t_rel_online(mask_ph);
data_full.yMeas{4} = tbl_online.ph(mask_ph);

% channels 5-6: filled below from offline_measurements + offline_measurements_returns

%% Query offline measurements — channel 5: S_IN (type = 'NH4-N')

% sample time (om.time) is used for tMeas — this is when the reactor state was observed.
% if a measurement has multiple return entries, the latest (highest omr.id) is taken
% to pick up re-processed values over the original ones.
sql_sin = sprintf( ...
    "SELECT DISTINCT ON (om.id) om.time, omr.value " + ...
    "FROM public.offline_measurements om " + ...
    "JOIN public.offline_measurements_returns omr ON omr.measurement_id = om.id " + ...
    "WHERE om.reactor_id = %d " + ...
    "  AND om.type = 'NH4-N' " + ...
    "  AND om.flag_valid != 0 " + ...
    "  AND om.time >= '%s' " + ...
    "  AND om.time <= '%s' " + ...
    "ORDER BY om.id, omr.id DESC;", ...
    reactor_id, T_start_str, T_end_str);

fprintf("Querying offline S_IN (NH4-N) measurements...\n");
tbl_sin = fetch(conn, sql_sin);
fprintf("  %d rows returned.\n", height(tbl_sin));

t_rel_sin          = days(tbl_sin.time - T_start);   % [d] relative to T_start
data_full.tMeas{5} = t_rel_sin;
data_full.yMeas{5} = tbl_sin.value;

%% Query offline measurements — channel 6: S_ac (type = 'GC')

sql_sac = sprintf( ...
    "SELECT DISTINCT ON (om.id) om.time, omr.value " + ...
    "FROM public.offline_measurements om " + ...
    "JOIN public.offline_measurements_returns omr ON omr.measurement_id = om.id " + ...
    "WHERE om.reactor_id = %d " + ...
    "  AND om.type = 'GC' " + ...
    "  AND om.flag_valid != 0 " + ...
    "  AND om.time >= '%s' " + ...
    "  AND om.time <= '%s' " + ...
    "ORDER BY om.id, omr.id DESC;", ...
    reactor_id, T_start_str, T_end_str);

fprintf("Querying offline S_ac (GC) measurements...\n");
tbl_sac = fetch(conn, sql_sac);
fprintf("  %d rows returned.\n", height(tbl_sac));

t_rel_sac          = days(tbl_sac.time - T_start);   % [d] relative to T_start
data_full.tMeas{6} = t_rel_sac;
data_full.yMeas{6} = tbl_sac.value;

%% Query feeding events

% feedings table holds one row per substrate per event; multiple substrates
% fed at the same timestamp are mixed into a single event below.
sql_feedings = sprintf( ...
    "SELECT feedings.time, feedings.amount, " + ...
    "COALESCE(substrates.BK_number, 'water') AS bk_number " + ...
    "FROM feedings " + ...
    "LEFT JOIN substrates ON feedings.substrate_id = substrates.id " + ...
    "WHERE feedings.reactor_id = %d " + ...
    "  AND feedings.time >= '%s' " + ...
    "  AND feedings.time <= '%s' " + ...
    "  AND feedings.flag_valid = true " + ...
    "ORDER BY feedings.time;", ...
    reactor_id, T_start_str, T_end_str);

fprintf("Querying feeding events...\n");
tbl_feedings = fetch(conn, sql_feedings);
fprintf("  %d feeding entries returned.\n", height(tbl_feedings));

%% Fetch xi-vectors for each unique substrate (BK_number)

% for each unique substrate, get the most recent xi-vector from xi_vectors;
% water feeds get a zero xi-vector and rho = 1000 kg/m³.
bk_names = string(unique(tbl_feedings.bk_number));
n_bk     = numel(bk_names);
xi_vecs  = cell(n_bk, 1);    % row vector (1 x n_xi) per substrate
xi_rhos  = nan(n_bk, 1);     % [kg/m³]

for i_bk = 1:n_bk
    bk = bk_names(i_bk);
    if strcmp(bk, "") % water
        sql_xi = "SELECT xi_vector FROM xi_vectors " + ...
                 "ORDER BY time DESC LIMIT 1;";
        tbl_xi        = fetch(conn, sql_xi);
        n_xi_len      = numel(jsondecode(tbl_xi.xi_vector));
        xi_vecs{i_bk} = zeros(1, n_xi_len);
        xi_rhos(i_bk) = rho_water;
    else
        sql_xi = sprintf( ...
            "SELECT xi_vector, rho FROM xi_vectors " + ...
            "LEFT JOIN substrates ON substrate_id = substrates.id " + ...
            "WHERE substrates.BK_number = '%s' " + ...
            "ORDER BY xi_vectors.time DESC LIMIT 1;", bk);
        tbl_xi        = fetch(conn, sql_xi);
        xi_vecs{i_bk} = jsondecode(tbl_xi.xi_vector)';   % column → row vector
        xi_rhos(i_bk) = tbl_xi.rho;
    end
end % for

%% Mix substrates at each unique event timestamp

t_abs_feedings = tbl_feedings.time;                   % datetime (TIMESTAMPTZ, timezone-aware)
t_rel_feedings = days(t_abs_feedings - T_start);      % [d] relative to T_start
bk_feedings    = string(tbl_feedings.bk_number);

event_times = unique(t_rel_feedings);
n_events    = numel(event_times);
n_xi        = numel(xi_vecs{1});

% intermediate arrays: one entry per unique raw timestamp
mass_raw   = nan(n_events, 1);     % [kg] total mass fed
vol_raw    = nan(n_events, 1);     % [m³] total volume (for xi volume-weighting only)
xi_raw_mat = nan(n_events, n_xi);  % volume-weighted xi

for i_ev = 1:n_events
    t_k      = event_times(i_ev);
    mask_k   = t_rel_feedings == t_k;
    amount_k = tbl_feedings.amount(mask_k);   % [kg] per substrate at this event
    bk_k     = bk_feedings(mask_k);

    [~, xi_idx] = ismember(bk_k, bk_names);
    rho_k       = xi_rhos(xi_idx);             % (n_subs x 1) [kg/m³]
    xi_k        = cell2mat(xi_vecs(xi_idx));   % (n_subs x n_xi)

    vol_k       = amount_k ./ rho_k;           % [m³] per substrate
    total_vol_k = sum(vol_k);                  % [m³] total

    mass_raw(i_ev)        = sum(amount_k);                             % [kg]
    vol_raw(i_ev)         = total_vol_k;
    xi_raw_mat(i_ev, :)   = (vol_k' * xi_k) / total_vol_k;            % volume-weighted mix
end % for

%% Bundle nearby events within delta_bundle_min

% a new group starts whenever the gap to the previous event reaches the threshold;
% all events within the same group are collapsed into a single feed pulse.
delta_bundle_days = delta_bundle_min / (24*60);   % [d]

group_id      = ones(n_events, 1);
current_group = 1;
for i_ev = 2:n_events
    if event_times(i_ev) - event_times(i_ev-1) >= delta_bundle_days
        current_group = current_group + 1;
    end
    group_id(i_ev) = current_group;
end % for
n_groups = current_group;

t_feed_start_arr = nan(n_groups, 1);
feed_mass_arr    = nan(n_groups, 1);
xi_feed_arr      = nan(n_groups, n_xi);

for i_g = 1:n_groups
    mask_g      = group_id == i_g;
    vols_g      = vol_raw(mask_g);
    total_vol_g = sum(vols_g);

    t_feed_start_arr(i_g) = min(event_times(mask_g));
    feed_mass_arr(i_g)    = sum(mass_raw(mask_g));                             % [kg]
    xi_feed_arr(i_g, :)   = (vols_g' * xi_raw_mat(mask_g, :)) / total_vol_g;  % volume-weighted
end % for

data_full.t_feed_start = t_feed_start_arr;
data_full.feed_mass    = feed_mass_arr;    % [kg]; t_feed_end and u_feed_value derived in main
data_full.xi_feed      = xi_feed_arr;

%% Timing fields

data_full.t0      = 0;
data_full.tf      = days(T_end - T_start);   % [d] total window duration
data_full.T_start = T_start;

%% Close database connection

close(conn);
fprintf("Database connection closed.\n");

%% Sanity check — sample counts per channel

output_names = {'v_dot_gas [m³/d]', 'pch4 [bar]', 'pco2 [bar]', 'ph [-]', ...
                'S_IN [g/L]', 'S_ac [g/L]'};
fprintf("\nSample counts per output channel:\n");
for i = 1:6
    fprintf("  ch %d  %-22s  %4d samples\n", i, output_names{i}, numel(data_full.tMeas{i}));
end
fprintf("  feed events: %d bundled from %d raw timestamps (%d substrate entries)\n", ...
    n_groups, n_events, height(tbl_feedings));

%% Save

fprintf("\nSaving data_full to %s...\n", full_path_raw);
save(full_path_raw, 'data_full');
fprintf("Saved.\n");
