%% test_estimate_sigma_online.m
% Estimates measurement noise sigmas (σ) for the three online model outputs
%   ch 1: gas_flow,  ch 2: p_CH4,  ch 4: pH
% from the preprocessed automated_feeder init dataset.
%
% Method: robust MAD-of-first-differences on quiet samples (option 2).
% Quiet = outside a ±1 h exclusion window around every bundled feeding event.
%   d      = diff(y_quiet)
%   sigma  = median(|d - median(d)|) / (0.6745 * sqrt(2))
%
% Results are printed for direct use as sigmaY entries in main_pss_workflow.
%
% Author: Simon Hellmann. Created: 2026/05/27. Version: Matlab R2022b, Update 6
%
% no inputs and outputs (main script)

clc; clear; close all

[script_path, ~, ~] = fileparts(mfilename('fullpath'));
cd(script_path);

%% Load data

data_file = fullfile('..', 'data', 'processed', 'automated_feeder', 'data_init.mat');
load(data_file, 'data_init');

%% Settings

dt_mask    = 1/24;       % [d] exclusion half-width around each feeding event (±1 h)
channels   = [1, 2, 4];  % output indices: gas_flow, p_CH4, pH
chan_names  = {"gas_flow", "p_CH4", "pH"};

%% Estimate sigma per channel

n_ch   = numel(channels);
sigma  = nan(n_ch, 1);
t_feed = data_init.t_feed_start;   % [d] bundled feeding event start times

fprintf("Estimating sigmas from %d feed events (±1 h mask)...\n\n", numel(t_feed));

for i_ch = 1:n_ch

    ch    = channels(i_ch);
    t_raw = data_init.tMeas{ch};
    y_raw = data_init.yMeas{ch};

    % Build quiet mask: exclude samples within ±dt_mask of any feeding event
    quiet_mask = true(numel(t_raw), 1);
    for i_ev = 1:numel(t_feed)
        in_window  = t_raw >= t_feed(i_ev) - dt_mask & ...
                     t_raw <= t_feed(i_ev) + dt_mask;
        quiet_mask = quiet_mask & ~in_window;
    end % for

    y_quiet = y_raw(quiet_mask);
    n_quiet = numel(y_quiet);

    if n_quiet < 4
        warning("test_estimate_sigma_online: ch %d has only %d quiet samples — sigma unreliable.", ...
            ch, n_quiet);
        continue
    end

    d           = diff(y_quiet);
    sigma(i_ch) = median(abs(d - median(d))) / (0.6745 * sqrt(2));

    fprintf("  ch %d (%s):  %4d / %4d samples quiet,  sigma = %.4g\n", ...
        ch, chan_names{i_ch}, n_quiet, numel(t_raw), sigma(i_ch));

end % for

%% Summary

fprintf("\n=== Sigma estimates for sigmaY in main_pss_workflow ===\n");
fprintf("  gas_flow  (ch 1):  %.4g\n", sigma(1));
fprintf("  p_CH4     (ch 2):  %.4g\n", sigma(2));
fprintf("  pH        (ch 4):  %.4g\n", sigma(3));
fprintf("\nCurrent sigmaY (intensiv):  [4e-4, 1.78e-2, _, 2e-2, _, _]\n");
