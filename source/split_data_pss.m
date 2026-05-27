%% split_data_pss.m
% Split a preprocessed MESS struct into three non-overlapping time windows:
%
%   data_init  -- initialization window (state warm-up, not used in fitting)
%   data_auto  -- auto-validation window (training for PI #1 and PI #2)
%   data_cross -- cross-validation window (independent evaluation)
%
% Expects the preprocessed MESS struct in data/raw/.
% Saves data_init, data_auto, data_cross to data/processed/.
%
% Compatibility: MATLAB R2022b
% Author: Simon Hellmann
% Created: 2026-05-17

clear; clc;

% change working directory to this script's location so relative paths work
[scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
cd(scriptPath);

%% paths
addpath('utils');
dataRawPath       = '../data/raw/';
dataProcessedPath = '../data/processed/';

%% -----------------------------------------------------------------------
%  USER SETTINGS -- adjust for each dataset
% -----------------------------------------------------------------------
flag_raw_file = 'feeder'; % 'feeder' | 'intensiv'

% --- time window boundaries ---------------------------------------------
% Adjust these datetimes to match the dataset in use.

% Initialization window: used to spin up the ODE to a realistic state.
T_init_start  = datetime('01-Jun-2022 00:00:00', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
T_init_end    = datetime('20-Jun-2022 00:00:00', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

% Auto-validation (training) window: PI #1, PSS, PI #2 all run on this.
T_auto_start  = datetime('20-Jun-2022 00:00:00', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
T_auto_end    = datetime('15-Jul-2022 00:00:00', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

% Cross-validation window: never used in fitting, only for evaluation.
T_cross_start = datetime('15-Jul-2022 00:00:00', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
T_cross_end   = datetime('05-Aug-2022 15:00:00', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

%% -----------------------------------------------------------------------
%  LOAD RAW DATA
% -----------------------------------------------------------------------

% --- raw data file (place in data/raw/) ---------------------------------
switch flag_raw_file
    case 'intensiv'
        rawFileName = 'MESS_struct_IntBePro_R36_mod_gasflow.mat'; 
        full_path_raw = fullfile(dataRawPath, rawFileName);
        fprintf('Loading raw data from %s...\n', full_path_raw);
        load(fullfile(dataRawPath, rawFileName), 'MESS');
    case 'feeder'
        rawFileName = 'data_raw_auto_feeder.mat';     
        full_path_raw = fullfile(dataRawPath, rawFileName);
        fprintf('Loading raw data from %s...\n', full_path_raw);
        load(fullfile(dataRawPath, rawFileName), 'data_full');
        MESS = data_full; % rename so compatible with the remaining script
    otherwise 
        error('Invalid raw file source!')
end
fprintf('Done.\n');

%% -----------------------------------------------------------------------
%  CUT TIME WINDOWS
% -----------------------------------------------------------------------

fprintf('Cutting time windows...\n');
data_init  = cutMESS(MESS, T_init_start,  T_init_end);
data_auto  = cutMESS(MESS, T_auto_start,  T_auto_end);
data_cross = cutMESS(MESS, T_cross_start, T_cross_end);

t0_init  = 0;
t0_auto  = days(T_auto_start  - T_init_start);
t0_cross = days(T_cross_start - T_init_start);
fprintf('  Init:  %2d feed events, abs. t = [%5.1f, %5.1f] d\n', ...
    numel(data_init.feed.rel_time),  t0_init,  t0_init  + data_init.t_span);
fprintf('  Auto:  %2d feed events, abs. t = [%5.1f, %5.1f] d\n', ...
    numel(data_auto.feed.rel_time),  t0_auto,  t0_auto  + data_auto.t_span);
fprintf('  Cross: %2d feed events, abs. t = [%5.1f, %5.1f] d\n', ...
    numel(data_cross.feed.rel_time), t0_cross, t0_cross + data_cross.t_span);

%% -----------------------------------------------------------------------
%  CONVERT TO PSS DATA FORMAT
% -----------------------------------------------------------------------

fprintf('Converting to PSS data format...\n');
data_init  = mess2pssData(data_init);
data_auto  = mess2pssData(data_auto);
data_cross = mess2pssData(data_cross);

%% -----------------------------------------------------------------------
%  SANITY CHECKS
% -----------------------------------------------------------------------

% Warn if any output channel has no measurements in a window
outputNames = {'gasflow','p_CH4','p_CO2','pH','S_IN','S_ac'};
for dataset_k = 1:3
    switch dataset_k
        case 1; d = data_init;  label = 'init';
        case 2; d = data_auto;  label = 'auto';
        case 3; d = data_cross; label = 'cross';
    end
    for i = 1:6
        if isempty(d.tMeas{i})
            warning('split_data_pss: output %s has NO measurements in %s window.', ...
                outputNames{i}, label);
        end
    end
end

% Check that windows do not overlap (windows are half-open [start, end))
if T_init_end > T_auto_start
    warning('split_data_pss: init and auto windows overlap.');
end
if T_auto_end > T_cross_start
    warning('split_data_pss: auto and cross windows overlap.');
end

%% -----------------------------------------------------------------------
%  SAVE
% -----------------------------------------------------------------------

fprintf('Saving to %s...\n', dataProcessedPath);
save(fullfile(dataProcessedPath, 'data_init.mat'),  'data_init');
save(fullfile(dataProcessedPath, 'data_auto.mat'),  'data_auto');
save(fullfile(dataProcessedPath, 'data_cross.mat'), 'data_cross');

fprintf('Saved: data_init.mat, data_auto.mat, data_cross.mat\n');
fprintf('Split complete.\n');
