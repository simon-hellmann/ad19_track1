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

% --- raw data file (place in data/raw/) ---------------------------------
rawFileName = 'MESS_struct.mat';    % update to actual filename

% --- feeding event settings ---------------------------------------------
delta_feed_min  = 5;               % feeding duration [min]
delta_feed_days = delta_feed_min / (24*60);   % convert to days

rho_substrate   = 1000;            % substrate density [kg/m^3]
                                   % used to convert feed mass -> volume flow

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

fprintf('Loading raw data from %s...\n', fullfile(dataRawPath, rawFileName));
load(fullfile(dataRawPath, rawFileName), 'MESS');
fprintf('Done.\n');

%% -----------------------------------------------------------------------
%  CUT TIME WINDOWS
% -----------------------------------------------------------------------

fprintf('Cutting time windows...\n');
MESS_init  = cutMESS(MESS, T_init_start,  T_init_end);
MESS_auto  = cutMESS(MESS, T_auto_start,  T_auto_end);
MESS_cross = cutMESS(MESS, T_cross_start, T_cross_end);

fprintf('  Init:  %d feed events, t = [0, %.1f] d\n', ...
    numel(MESS_init.feed.rel_time),  MESS_init.t_span);
fprintf('  Auto:  %d feed events, t = [0, %.1f] d\n', ...
    numel(MESS_auto.feed.rel_time),  MESS_auto.t_span);
fprintf('  Cross: %d feed events, t = [0, %.1f] d\n', ...
    numel(MESS_cross.feed.rel_time), MESS_cross.t_span);

%% -----------------------------------------------------------------------
%  CONVERT TO PSS DATA FORMAT
% -----------------------------------------------------------------------

fprintf('Converting to PSS data format...\n');
data_init  = mess2pssData(MESS_init,  delta_feed_days, rho_substrate);
data_auto  = mess2pssData(MESS_auto,  delta_feed_days, rho_substrate);
data_cross = mess2pssData(MESS_cross, delta_feed_days, rho_substrate);

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

% Check that windows do not overlap
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
