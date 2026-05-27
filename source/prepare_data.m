%% prepare_data.m
% Load, preprocess, and split raw experiment data into three non-overlapping
% time windows for use in the PI/PSS workflow (main_pss_workflow.m):
%
%   data_init  -- initialisation window (state warm-up, not used in fitting)
%   data_auto  -- auto-validation window (training for PI #1, PSS, PI #2)
%   data_cross -- cross-validation window (independent evaluation)
%
% Set the dataset flag in the USER SETTINGS block to switch experiments.
% All file paths and names are derived from that flag — nothing else changes.
%
% Reads from  : data/raw/<dataset>/
% Writes to   : data/processed/<dataset>/
%
% Author: Simon Hellmann. Created: 2026/05/27. Version: Matlab R2022b, Update 6
%
% no inputs and outputs (main script)

clc; clear; close all

[script_path, ~, ~] = fileparts(mfilename('fullpath'));
cd(script_path);
addpath('utils');

%% User settings

dataset = 'intensiv';   % 'intensiv' | 'automated_feeder'

% --- Preprocessing options (gas channel 1; applied to the full dataset) --
feeding_duration              = 15/(24*60);  % [d] feeding pulse length

% --- Time window boundaries (absolute datetimes; adjust per dataset) -----
switch dataset
    case 'intensiv'
        T_init_start  = datetime('01-Jun-2022 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss');
        T_init_end    = datetime('20-Jun-2022 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss');
        T_auto_start  = datetime('20-Jun-2022 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss');
        T_auto_end    = datetime('15-Jul-2022 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss');
        T_cross_start = datetime('15-Jul-2022 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss');
        T_cross_end   = datetime('05-Aug-2022 15:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss');

        preproc_opts.flag_filter_feed = true;   % true: remove gas data around feedings
        preproc_opts.flag_filter_IN   = true;   % true: remove gas data around IN samples
        preproc_opts.flag_filter_AC   = true;   % true: remove gas data around AC samples
        preproc_opts.dt_feed_before   = 1/24;   % [d] 
        preproc_opts.dt_feed_after    = 2/24;   % [d] 
        preproc_opts.dt_IN_before     = 0.5/24; % [d] 
        preproc_opts.dt_IN_after      = 1/24;   % [d]
        preproc_opts.dt_AC_before     = 0.5/24; % [d]
        preproc_opts.dt_AC_after      = 1/24;   % [d]
        preproc_opts.q_gas_min        = 0.002;  % [m³/d] below this: exclude
    case 'automated_feeder'
        % Adjust these boundaries to match your desired analysis period
        T_init_start  = datetime('29-Apr-2026 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
        T_init_end    = datetime('04-May-2026 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
        T_auto_start  = datetime('04-May-2026 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
        T_auto_end    = datetime('18-May-2026 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
        T_cross_start = datetime('18-May-2026 00:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
        T_cross_end   = datetime('27-May-2026 12:00:00', 'InputFormat','dd-MMM-yyyy HH:mm:ss', 'TimeZone','Europe/Berlin');
        
        preproc_opts.flag_filter_feed = false;   % true: remove gas data around feedings
        preproc_opts.flag_filter_IN   = false;   % true: remove gas data around IN samples
        preproc_opts.flag_filter_AC   = false;   % true: remove gas data around AC samples
        preproc_opts.q_gas_min        = 0.002;  % [m³/d] below this: exclude
    otherwise
        error('prepare_data: unknown dataset ''%s''.', dataset);
end

%% Paths  (derived from dataset flag — no hardcoded names beyond here)

raw_dir       = fullfile('..', 'data', 'raw',       dataset);
processed_dir = fullfile('..', 'data', 'processed', dataset);

%% Load and convert to unified PSS format
% Both sources produce an identical PSS struct (data_raw) after this block.
% Source-specific variable names are confined to their case branch.

fprintf("Loading %s...\n", dataset);
switch dataset
    case 'intensiv'
        load(fullfile(raw_dir, 'MESS_struct_IntBePro_R36_mod_gasflow.mat'), 'MESS');
        MESS_full = cutMESS(MESS, T_init_start, T_cross_end);
        data_raw  = mess2pssData(MESS_full);
    case 'automated_feeder'
        load(fullfile(raw_dir, 'data_raw_auto_feeder.mat'), 'data_full');
        data_raw = data_full;
end
fprintf("  Loaded.\n");

%% Bundle feeding events

fprintf("Bundling feed events...\n");
data_raw = bundleFeedings(data_raw, feeding_duration);
fprintf("  Done.\n");

%% Preprocess (full dataset, before splitting)

fprintf("Preprocessing...\n");
data_raw = preprocessData(data_raw, preproc_opts, feeding_duration);
fprintf("  Done.\n");

%% Split

fprintf("Splitting time windows...\n");
data_init  = cutPSS(data_raw, T_init_start,  T_init_end);
data_auto  = cutPSS(data_raw, T_auto_start,  T_auto_end);
data_cross = cutPSS(data_raw, T_cross_start, T_cross_end);

t0_init  = 0;
t0_auto  = days(T_auto_start  - T_init_start);   % [d] abs. offset from t = 0
t0_cross = days(T_cross_start - T_init_start);   % [d]

fprintf("  Init:  %2d feed events, abs. t = [%5.1f, %5.1f] d\n", ...
    numel(data_init.t_feed_start),  t0_init,  t0_init  + data_init.tf);
fprintf("  Auto:  %2d feed events, abs. t = [%5.1f, %5.1f] d\n", ...
    numel(data_auto.t_feed_start),  t0_auto,  t0_auto  + data_auto.tf);
fprintf("  Cross: %2d feed events, abs. t = [%5.1f, %5.1f] d\n", ...
    numel(data_cross.t_feed_start), t0_cross, t0_cross + data_cross.tf);

%% Sanity checks

output_names = {'gasflow','p_CH4','p_CO2','pH','S_IN','S_ac'};
for i_ds = 1:3
    switch i_ds
        case 1; d = data_init;  label = 'init';
        case 2; d = data_auto;  label = 'auto';
        case 3; d = data_cross; label = 'cross';
    end
    for i_ch = 1:6
        if isempty(d.tMeas{i_ch})
            warning('prepare_data: output %s has NO measurements in %s window.', ...
                output_names{i_ch}, label);
        end
    end
end

if T_init_end > T_auto_start
    warning('prepare_data: init and auto windows overlap.');
end
if T_auto_end > T_cross_start
    warning('prepare_data: auto and cross windows overlap.');
end

%% Save

fprintf("Saving to %s...\n", processed_dir);
save(fullfile(processed_dir, 'data_init.mat'),       'data_init');
save(fullfile(processed_dir, 'data_auto.mat'),       'data_auto');
save(fullfile(processed_dir, 'data_cross.mat'),      'data_cross');
save(fullfile(processed_dir, 'feeding_duration.mat'), 'feeding_duration');
fprintf("Saved: data_init.mat, data_auto.mat, data_cross.mat\n");
fprintf("Done.\n");
