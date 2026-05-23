%% main_pss_workflow.m
% Full workflow: Parameter Identification -> PSS -> Parameter Identification
%
% Goal: Show that PSS yields comparable fit quality while drastically
%       reducing parameter uncertainty (wider confidence intervals before
%       PSS, tight ones after).
%
% Compatibility: MATLAB R2022b
%
% Sections:
%   1.  Setup & paths
%   2.  Model definition & nominal parameters
%   3.  Load datasets (output of split_data_pss.m)
%   4.  State initialisation for PI #1 -- simulate data_init with theta0
%   5.  PI #1 -- WLS on full parameter vector (fmincon)
%   6.  Post-PI #1 analysis: fit, Fisher uncertainty, cross-validation
%   7.  PSS -- compute scaled sensitivities, SVD, QRP, select subset
%   8.  State initialisation for PI #2 -- re-simulate data_init with thetaHat1
%   9.  PI #2 -- WLS on identifiable parameter subset (fmincon)
%   10. Post-PI #2 analysis: fit, Fisher uncertainty, cross-validation
%   11. Comparison plots: fit overlay, confidence interval shrinkage
%
% Note on state initialisation (sections 4 and 8):
%   x0 cannot be precomputed because it depends on the current parameter
%   vector.  Simulating data_init forward propagates a rough starting state
%   into a physically realistic one; the final state of that simulation
%   becomes the initial condition for the autovalidation (and CV) windows.
%
% Author: Simon Hellmann
% Created: 2026-05-17

clear; clc; close all;

%% -----------------------------------------------------------------------
%  1. SETUP & PATHS
% -----------------------------------------------------------------------

% add subfolders (model files, helper functions, PSS routine, etc.)
addpath('model_files');
addpath('utils');  

% random seed for reproducibility (only matters if using stochastic optimiser)
rng(42);

%% -----------------------------------------------------------------------
%  2. MODEL DEFINITION & NOMINAL PARAMETERS
% -----------------------------------------------------------------------

% --- Reactor dimensions --------------------------------------------------
V_liq = 0.012;   % liquid volume [m^3]
V_gas = 0.003;   % gas headspace volume [m^3]

% --- Load ADM1 physico-chemical constants --------------------------------
addpath('model_files/model_data');
load('ADM1_parameters.mat', 'parameters');
parameters_r3 = parameters.ADM1_R3;

% --- Assemble model (ADM1-R3) --------------------------------------------
% Returns c, a, odeFunc, measFunc, theta0, bounds, and metadata struct p.
% To switch models, replace this call (e.g. setup_ADM1_R3x1).
[c, a, odeFunc, measFunc, theta0, thetaLB, thetaUB, p] = ...
    setup_ADM1_R3(parameters_r3, V_liq, V_gas);

% --- ODE solver options --------------------------------------------------
% Loose tolerances inside fmincon (~3-5x faster, negligible effect on optimum).
% NonNegative prevents ode15s driving biological concentration states below zero.
non_negative_state_idx = 1:14;
odeOptsOpt  = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 0.5, ...
                     'NonNegative', non_negative_state_idx);
% Tight tolerances for post-processing (FD sensitivity, CV, plots):
odeOptsPost = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.5, ...
                     'NonNegative', non_negative_state_idx);

% --- Output measurement noise standard deviations -----------------------
% Units: [m^3/d, bar, bar, -, g/L, g/L]
sigmaY = [4e-4, 1.78e-2, 2.68e-2, 2e-2, 0.12, 5e-2];

%% -----------------------------------------------------------------------
%  3. LOAD & PREPARE DATASETS
% -----------------------------------------------------------------------

% Note: Run split_data_pss.m first to generate these files from the raw MESS struct.
%
load('../data/processed/data_init.mat',  'data_init');
load('../data/processed/data_auto.mat',  'data_auto');
load('../data/processed/data_cross.mat', 'data_cross');

% Unpack training (auto-validation) dataset:
%   tMeas{i}        -- (n_i x 1) measurement times per output [d]
%   yMeas{i}        -- (n_i x 1) measured values per output
%   t_feed_start    -- (n_events x 1) feed start times [d]
%   t_feed_end      -- (n_events x 1) feed end times [d]
%   u_feed_value    -- (n_events x 1) feed volume flow [m^3/d]
%   t0, tf          -- simulation window [d]
tMeas        = data_auto.tMeas;
yMeas        = data_auto.yMeas;
t_feed_start = data_auto.t_feed_start;
t_feed_end   = data_auto.t_feed_end;
u_feed_value = data_auto.u_feed_value;
t0           = data_auto.t0;
tf           = data_auto.tf;
n_out        = 6; % # model output channels

% Unpack cross-validation dataset (same fields, _CV suffix):
tMeasCV        = data_cross.tMeas;
yMeasCV        = data_cross.yMeas;
t_feed_start_CV = data_cross.t_feed_start;
t_feed_end_CV   = data_cross.t_feed_end;
u_feed_value_CV = data_cross.u_feed_value;
t0_CV           = data_cross.t0;
tf_CV           = data_cross.tf;

% --- Build the long observation vector (autovalidation) -----------------------
% Entries are ordered by time, then by output index within each time point.
% At each measurement time only the outputs actually sampled there appear.
%
%   y_meas_long: (N_long x 1)  -- stacked measured values
%   t_meas_long: (N_long x 1)  -- corresponding time points
%   out_idx:     (N_long x 1)  -- output index (1..n_out) for each entry

rows_raw = [];
for out_k = 1:n_out % iterate over output channels
    n_samples_k   = numel(tMeas{out_k});
    rows_raw = [rows_raw; tMeas{out_k}(:), out_k*ones(n_samples_k,1), yMeas{out_k}(:)];
end
rows        = sortrows(rows_raw, [1, 2]);  % primary: time, secondary: output index
t_meas_long = rows(:,1);
out_idx     = rows(:,2);
y_meas_long = rows(:,3);

% --- Scaling vector for the autovalidation long residual vector ---------
% Each residual r_j is divided by sigma_k * sqrt(n_samples_k), where k = 
% out_idx(j) and n_samples_k is the total number of samples for output k 
% across the dataset.
nSamples   = cellfun(@numel, tMeas);                        % (n_out x 1)
scale_long = sigmaY(out_idx)' .* sqrt(nSamples(out_idx));   % (N_long x 1)

% --- Feeding event grid for autovalidation ------------------------------
t_events   = unique([t0; t_feed_start(:); t_feed_end(:); tf]);
t_mid      = (t_events(1:end-1) + t_events(2:end)) / 2;
n_xi        = size(data_auto.xi_feed, 2);
u_segments  = zeros(numel(t_events)-1, 1);
xi_segments = zeros(numel(t_events)-1, n_xi);
for event_k = 1:numel(t_feed_start)
    active              = (t_mid >= t_feed_start(event_k)) & (t_mid <= t_feed_end(event_k));
    u_segments(active)  = u_feed_value(event_k);
    xi_segments(active,:) = repmat(data_auto.xi_feed(event_k,:), sum(active), 1);
end

% --- Build the long observation vector for cross-validation -------------
% Identical construction; uses CV measurement cell arrays.

rowsCV_raw        = [];
for out_k = 1:n_out
    n_samples_k   = numel(tMeasCV{out_k});
    rowsCV_raw    = [rowsCV_raw; tMeasCV{out_k}(:), out_k*ones(n_samples_k,1), yMeasCV{out_k}(:)];
end
rowsCV           = sortrows(rowsCV_raw, [1, 2]);
t_meas_long_CV   = rowsCV(:,1);
out_idx_CV       = rowsCV(:,2);
y_meas_long_CV   = rowsCV(:,3);

% --- Scaling vector for the CV long residual vector ---------------------
% n_samples_k is taken from the CV dataset, not the autovalidation dataset, so 
% the normalisation reflects the actual CV sample density per output.

nSamplesCV    = cellfun(@numel, tMeasCV);                            % (n_out x 1)
scale_long_CV = sigmaY(out_idx_CV)' .* sqrt(nSamplesCV(out_idx_CV)); % (N_long_CV x 1)

% --- Feeding event grid for cross-validation ----------------------------
t_events_CV  = unique([t0_CV; t_feed_start_CV(:); t_feed_end_CV(:); tf_CV]);
t_mid_CV     = (t_events_CV(1:end-1) + t_events_CV(2:end)) / 2;
u_segments_CV  = zeros(numel(t_events_CV)-1, 1);
xi_segments_CV = zeros(numel(t_events_CV)-1, n_xi);
for event_k = 1:numel(t_feed_start_CV)
    active = (t_mid_CV >= t_feed_start_CV(event_k)) & (t_mid_CV <= t_feed_end_CV(event_k));
    u_segments_CV(active)   = u_feed_value_CV(event_k);
    xi_segments_CV(active,:) = repmat(data_cross.xi_feed(event_k,:), sum(active), 1);
end

%% -----------------------------------------------------------------------
%  4. STATE INITIALISATION FOR PI #1
% -----------------------------------------------------------------------

% Two-step dynamic initialisation (see utils/computeX0.m):
%
%   Step 1 -- 500 d steady-state pre-simulation under average feed rate
%     Drives an arbitrary rough state to a representative operating point.
%     The average feed rate is computed from data_init automatically.
%
%   Step 2 -- Piecewise swing-up across data_init
%     Simulates the actual (discrete) feeding schedule of data_init,
%     starting from the steady state.  The terminal state is x0.

% rough initial state (n_states x 1), cf. Hellmann et al. (2026), DOI: 10.1016/j.jprocont.2026.103703
x0_init = [0.049; % S_ac 
            0.012; % S_ch4
            4.975; % S_IC
            0.964; % S_IN
            2.962; % X_ch
            0.949; % X_pr
            0.412; % X_li
            1.926; % X_bac
            0.522; % X_ac
            0.049; % S_ac-
            4.546; % S_hco3-
            0.022; % S_nh3
            0.358; % S_ch4_gas
            0.660];% S_co2_gas
t_ss     = 500;     % [d] pre-simulation duration for steady-state
x0 = computeX0(theta0, data_init, t_ss, x0_init, odeFunc, odeOptsPost);

%% -----------------------------------------------------------------------
%  4b. LATIN HYPERCUBE SAMPLING -- GLOBAL PRE-SEARCH FOR PI #1
% -----------------------------------------------------------------------
%
%  Draws N_LHS parameter vectors via Latin Hypercube Sampling over the full
%  parameter bounds, evaluates the cost function at each, and selects the
%  best candidate as the starting point for PI #1.
%
%  Sampling scale:
%    - Parameters 1-7 and 9 (strictly positive bounds): log10-scale.
%    - Parameter 8 (Delta_S_ion, LB = -1e-2 < 0): linear scale.

N_LHS     = 50;
n_theta   = numel(theta0);
J_penalty = 1e10;   % also used in section 5

log_idx = setdiff(1:n_theta, 8);   % all parameters except Delta_S_ion
lin_idx = 8;                       % Delta_S_ion: bounds straddle zero -> linear

% LHS design in [0,1]^n_theta
lhs_unit = lhsdesign(N_LHS, n_theta);   % (N_LHS x n_theta)

% Map unit hypercube to physical parameter space
theta_lhs = nan(N_LHS, n_theta);
for k = log_idx
    lo = log10(thetaLB(k));
    hi = log10(thetaUB(k));
    theta_lhs(:, k) = 10.^(lo + lhs_unit(:, k) .* (hi - lo));
end
theta_lhs(:, lin_idx) = thetaLB(lin_idx) + ...
    lhs_unit(:, lin_idx) .* (thetaUB(lin_idx) - thetaLB(lin_idx));

% Evaluate cost at each LHS sample (optimization-grade ODE tolerances)
fprintf('Evaluating cost at %d LHS samples (parfor)...\n', N_LHS);
J_lhs = nan(N_LHS, 1);
maxNumCores = feature('numCores');
if isempty(gcp('nocreate'))
    parpool(maxNumCores, 'IdleTimeout', Inf); % pool that never loses connection
end
tic_lhs = tic;
parfor i_lhs = 1:N_LHS
    J_lhs(i_lhs) = costWLS(theta_lhs(i_lhs, :)', y_meas_long, t_meas_long, ...
        out_idx, scale_long, t_events, u_segments, xi_segments, x0, ...
        odeFunc, measFunc, odeOptsOpt, J_penalty);
end
computing_time_lhs = toc(tic_lhs); % [s]
fprintf('LHS done.  Cost range: [%.4g, %.4g]  (%.1f s)\n', min(J_lhs), max(J_lhs), computing_time_lhs);
delete(gcp('nocreate'));

% Best LHS candidate becomes the starting point for PI #1
[J0, lhs_best_idx] = min(J_lhs);
theta_lhs_best    = theta_lhs(lhs_best_idx, :)';
fprintf('Best LHS candidate: sample #%d,  J = %.4g\n', lhs_best_idx, min(J_lhs));

%% -----------------------------------------------------------------------
%  5. PI #1 -- MAX LIKELIHOOD ESTIMATION ON FULL PARAMETER VECTOR (fmincon)
% -----------------------------------------------------------------------

% Core simulation helper -- shared by costWLS and the FD sensitivity routine:
%
% See utils/simulateLong.m -- piecewise ODE integration then long-vector assembly.
% See utils/costWLS.m      -- thin wrapper: residuals + J = r'*r.

% Evaluate cost at theta0 once (with optimization tolerances) to normalize J.
% Dividing by J0 brings the cost to O(1) at the start, which keeps the
% gradient components in a numerically tractable range for fmincon's FD.
% The optimum location is unchanged; only the gradient magnitude is rescaled.
J0 = costWLS(theta0, y_meas_long, t_meas_long, out_idx, ...
             scale_long, t_events, u_segments, xi_segments, x0, ...
             odeFunc, measFunc, odeOptsOpt, J_penalty);
fprintf('J(theta0) = %.4g  (normalization factor)\n', J0);

objFun1 = @(theta) costWLS(theta, y_meas_long, t_meas_long, out_idx, ...
              scale_long, t_events, u_segments, xi_segments, x0, ...
              odeFunc, measFunc, odeOptsOpt, J_penalty) / J0;

% --- fmincon options ----------------------------------------------------
fd_stepSize = odeOptsOpt.RelTol^(1/3); % FD step size ≈ 0.046 for RelTol=1e-4
opts1 = optimoptions('fmincon', ...
    'Display',                  'iter-detailed', ...
    'Algorithm',                'interior-point', ...
    'MaxFunctionEvaluations',   500, ...
    'TypicalX',                 theta0, ...
    'FiniteDifferenceStepSize', fd_stepSize, ... 
    'FiniteDifferenceType',     'central');

% --- Run PI #1 (starting from best LHS candidate) -----------------------
disp("Running PI1 with fmincon...")
tic_pi1 = tic;
[thetaHat1, fval1_norm, exitflag1, output1] = fmincon(objFun1, theta_lhs_best, ...
    [], [], [], [], thetaLB, thetaUB, [], opts1);
computing_time_pi1 = toc(tic_pi1); % [s]

% thetaHat1: estimated parameter vector after PI #1

%% -----------------------------------------------------------------------
%  6. POST-PI #1 ANALYSIS
% -----------------------------------------------------------------------

% ---- 5a. Simulate with thetaHat1 and compute fit quality ---------------
% ySimLong1 has the same row ordering as y_meas_long (time-sorted).

[ySimLong1, x0_CV] = simulateLong(thetaHat1, t_meas_long, out_idx, t_events, ...
                         u_segments, xi_segments, x0, odeFunc, measFunc, odeOptsPost);
r_scaled1   = (ySimLong1 - y_meas_long) ./ scale_long; % scaled residuals
RMSE1_auto = sqrt(mean(r_scaled1.^2));
plotFit(t_meas_long, y_meas_long, ySimLong1, out_idx, p, 'PI #1 --autovalidation fit');

% ---- 5b. Scaled output sensitivity matrix via central finite differences -
%
% Two sensitivity matrices are built from dydth1_raw (see Step 2 below):
%   dydth1_os(j,k) = (d y_j/d theta_k) / scale_long(j)
%     Output-scaled only.  FIM = dydth1_os'*dydth1_os gives Fisher information
%     in physical parameter units (consistent with the WLS cost).
%   dydth1(j,k) = dydth1_os(j,k) * thetaHat1(k)
%     Also parameter-scaled (dimensionless).  Used for PSS only.
%     Must NOT be used for FIM: the thetaHat scaling introduces D_theta on
%     both sides of the FIM (FIM_scaled = D_theta*FIM_wls*D_theta), giving
%     std devs in the normalized space rather than physical units.
%
% Step 1: compute raw (unscaled) one-at-a-time FD sensitivities
%   dydth1_raw(j,k) = d y_j / d theta_k   [output units / parameter units]
%   rows follow the same ordering as y_meas_long (time-sorted, then out_idx)

h_rel        = odeOptsPost.RelTol^(1/3);   % ≈ 0.01 for RelTol=1e-6; matched to ODE noise floor
N_long       = numel(t_meas_long);
dydth1_raw   = nan(N_long, p.nParameters); 

disp("Computing FD sensitivities (parfor)...")
if isempty(gcp('nocreate'))
    parpool(maxNumCores,'IdleTimeout',Inf); % pool that never loses connection
end

parfor k = 1:p.nParameters
    fprintf("Computing FD sensitivities of param %i of %i...\n", k, p.nParameters); 
    delta_k   = h_rel * abs(thetaHat1(k));
    if delta_k == 0;  delta_k = h_rel;  end

    theta_fwd = thetaHat1;  theta_fwd(k) = thetaHat1(k) + delta_k;
    theta_bwd = thetaHat1;  theta_bwd(k) = thetaHat1(k) - delta_k;

    y_fwd = simulateLong(theta_fwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0, odeFunc, measFunc, odeOptsPost);
    y_bwd = simulateLong(theta_bwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0, odeFunc, measFunc, odeOptsPost);

    dydth1_raw(:,k) = (y_fwd - y_bwd) ./ (2*delta_k);
end
delete(gcp('nocreate')) % Shut down parallel pool

% Step 2: two separate scaled matrices with distinct purposes
%
%   dydth1_os  -- output-scaled only (for FIM / Cramer-Rao)
%     Row j divided by scale_long(j) = sigma_i * sqrt(n_i).
%     FIM = dydth1_os' * dydth1_os gives the Fisher information in the
%     original physical parameter units, consistent with the WLS cost.
%
%   dydth1_ops -- output- AND parameter-scaled (for PSS only)
%     Additionally multiply column k by thetaHat1(k) to make all columns
%     dimensionless and comparable across parameters of different magnitudes.
%     This matrix must NOT be used for FIM: it introduces D_theta on both
%     sides of the FIM (FIM_scaled = D_theta * FIM_wls * D_theta), which
%     would give std devs in normalized space, not physical units.

dydth1_os = dydth1_raw ./ scale_long;       % (N_long x nP): output-scaled
dydth1_ops = dydth1_os  .* thetaHat1(:)';   % (N_long x nP): also param-scaled (PSS only)

FIM1      = dydth1_os' * dydth1_os;         % Fisher Information Matrix (physical units)
C_theta1  = inv(FIM1);                      % Cramer-Rao lower bound
stdTheta1 = sqrt(diag(C_theta1));           % std dev in physical parameter units
plotUncertainty(thetaHat1, stdTheta1, p, 'PI #1 --parameter uncertainty');

% ---- 5c. Cross-validation with independent dataset ---------------------
ySimLong1_CV = simulateLong(thetaHat1, t_meas_long_CV, out_idx_CV, ...
                   t_events_CV, u_segments_CV, xi_segments_CV, x0_CV, odeFunc, measFunc, odeOptsPost);
r_scaled1_CV = (ySimLong1_CV - y_meas_long_CV) ./ scale_long_CV;
RMSE1_cv     = sqrt(mean(r_scaled1_CV.^2));
plotFit(t_meas_long_CV, y_meas_long_CV, ySimLong1_CV, out_idx_CV, p, ...
        'PI #1 --cross-validation');

%% -----------------------------------------------------------------------
%  7. PARAMETER SUBSET SELECTION (PSS)
% -----------------------------------------------------------------------

% Uses the scaled sensitivity matrix from PI #1 as input.
% Applies SVD + QRP decomposition to identify the practically identifiable
% parameter subset (see subsetSelection.m -- Lopez et al., 2015).

% --- PSS thresholds ------------------------------
kappa_max = 500; % maximum condition number
gamma_max = 15;   % maximum collinearity index

% --- Run PSS ------------------------------------------------------------
% dydth1 is doubly-scaled (rows / scale_long, cols * thetaHat1).
% C_pp returned is in physical parameter units (re-transformed inside).
[si, keep_idx, C_pp, pi_decomp, epsilon, kappa, gamma] = param_subset_selection(dydth1_ops, ...
    kappa_max, gamma_max, p, thetaHat1);

% keep_idx: indices of identifiable parameters
% The complement can be fixed at thetaHat1 or at literature values.

% --- Summarise PSS result -----------------------------------------------
fprintf('\nPSS result: %d / %d parameters are identifiable.\n', ...
        numel(keep_idx), p.nParameters);
disp('Identifiable parameters:');
disp(p.names(keep_idx));

plotSingularValues(si, epsilon, p, 'PSS --singular value spectrum');
plotVarianceDecomposition(pi_decomp, p, 'PSS --variance decomposition');

%% -----------------------------------------------------------------------
%  8. STATE INITIALISATION FOR PI #2
% -----------------------------------------------------------------------

% thetaHat1 may differ substantially from theta0, so x0 from section 4
% is no longer the best initial condition.  Repeat the same two-step
% procedure with thetaHat1 to obtain x0_2.

x0_2 = computeX0(thetaHat1, data_init, t_ss, x0_init, odeFunc, odeOptsPost);

% For cross-validation, data_cross immediately follows data_auto in time,
% so x0_CV is the terminal state of the auto simulation (captured in §6).
% x0_CV_2 (self-consistent with thetaHat2) is computed in §10 after PI #2.

%% -----------------------------------------------------------------------
%  9. PI #2 -- MAX LIKELIHOOD ESTIMATION ON IDENTIFIABLE PARAMETER SUBSET
% -----------------------------------------------------------------------

% Non-identifiable parameters are fixed at thetaHat1 (their PI #1 values).
% The same long-vector cost function is reused; only the free variables change.

% --- Initialise fixed and free parameter vectors ------------------------
thetaFixed  = thetaHat1;          % full vector; non-id entries stay fixed
theta1_sub  = thetaHat1(keep_idx);
thetaLB_sub = thetaLB(keep_idx);
thetaUB_sub = thetaUB(keep_idx);

% --- Objective: embed the free sub-vector back into the full vector -----
% See utils/costWLS_sub.m

objFun2 = @(theta_sub) costWLS_sub(theta_sub, keep_idx, thetaFixed, ...
              y_meas_long, t_meas_long, out_idx, scale_long, ...
              t_events, u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsOpt, J_penalty) / J0;

% --- fmincon options for PI #2 (TypicalX uses the subset estimate) ------
opts2 = optimoptions(opts1, 'TypicalX',theta1_sub);

% --- Run PI #2 ----------------------------------------------------------
disp("Running PI2 with fmincon...")
tic_pi2 = tic; 
[thetaSub2, fval2, exitflag2, output2] = fmincon(objFun2, theta1_sub, ...
    [], [], [], [], thetaLB_sub, thetaUB_sub, [], opts2);
computing_time_pi2 = toc(tic_pi2); % [s]

% --- Reconstruct the full parameter vector ------------------------------
thetaHat2           = thetaFixed;
thetaHat2(keep_idx) = thetaSub2;

%% -----------------------------------------------------------------------
%  10. POST-PI #2 ANALYSIS
% -----------------------------------------------------------------------

% ---- 10a. Fit quality ---------------------------------------------------
[ySimLong2, x0_CV_2] = simulateLong(thetaHat2, t_meas_long, out_idx, t_events, ...
                           u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsPost);
r_scaled2   = (ySimLong2 - y_meas_long) ./ scale_long;
RMSE2_auto = sqrt(mean(r_scaled2.^2));
plotFit(t_meas_long, y_meas_long, ySimLong2, out_idx, p, 'PI #2 --autovalidation fit');

% ---- 10b. Scaled output sensitivity matrix for the identifiable subset --
% Same FD procedure as 5b, but columns restricted to keep_idx.
% dydth2 has the same row structure as y_meas_long; it has numel(keep_idx)
% columns rather than p.nParameters.
%
% Step 1: raw FD sensitivities, columns restricted to keep_idx only

dydth2_raw = nan(N_long, numel(keep_idx));

disp("Computing FD sensitivities (parfor)...")
if isempty(gcp('nocreate'))
    parpool(maxNumCores,'IdleTimeout',Inf); % pool that never loses connection
end

parfor ki = 1:numel(keep_idx)
    k         = keep_idx(ki);
    delta_k   = h_rel * abs(thetaHat2(k));
    if delta_k == 0;  delta_k = h_rel;  end

    theta_fwd = thetaHat2;  theta_fwd(k) = thetaHat2(k) + delta_k;
    theta_bwd = thetaHat2;  theta_bwd(k) = thetaHat2(k) - delta_k;

    y_fwd = simulateLong(theta_fwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsPost);
    y_bwd = simulateLong(theta_bwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsPost);

    dydth2_raw(:,ki) = (y_fwd - y_bwd) ./ (2*delta_k);
end
delete(gcp('nocreate')) % Shut down parallel pool

% Step 2: same two-matrix scheme as §6b
%   dydth2_os  -- output-scaled only (for FIM / Cramer-Rao, physical units)
%   dydth2     -- also param-scaled by thetaHat2(keep_idx) (for reference/PSS)

dydth2_os = dydth2_raw ./ scale_long;
dydth2    = dydth2_os  .* thetaHat2(keep_idx)';

FIM2          = dydth2_os' * dydth2_os;
C_theta2_sub  = inv(FIM2);
stdTheta2_sub = sqrt(diag(C_theta2_sub));
stdTheta2_full             = NaN(p.nParameters, 1);
stdTheta2_full(keep_idx)   = stdTheta2_sub;
plotUncertainty(thetaHat2, stdTheta2_full, p, ...
                'PI #2 --parameter uncertainty (identifiable subset)');

% ---- 10c. Cross-validation ----------------------------------------------
ySimLong2_CV = simulateLong(thetaHat2, t_meas_long_CV, out_idx_CV, ...
                   t_events_CV, u_segments_CV, xi_segments, x0_CV_2, odeFunc, measFunc, odeOptsPost);
r_scaled2_CV = (ySimLong2_CV - y_meas_long_CV) ./ scale_long_CV;
RMSE2_cv     = sqrt(mean(r_scaled2_CV.^2));
plotFit(t_meas_long_CV, y_meas_long_CV, ySimLong2_CV, out_idx_CV, p, ...
        'PI #2 --cross-validation');

%% -----------------------------------------------------------------------
%  11. COMPARISON: BEFORE vs. AFTER PSS
% -----------------------------------------------------------------------

% ---- 11a. Fit overlay (autovalidation + CV) -----------------------------------
% Show that ySimLong1 ~= ySimLong2 (fit quality preserved after PSS)
plotFitComparison(t_meas_long, y_meas_long, ySimLong1, ySimLong2, out_idx, p, ...
    'Fit comparison: PI #1 vs PI #2 (autovalidation)');
plotFitComparison(t_meas_long_CV, y_meas_long_CV, ySimLong1_CV, ySimLong2_CV, ...
    out_idx_CV, p, 'Fit comparison: PI #1 vs PI #2 (cross-validation)');

% ---- 11b. Confidence interval shrinkage ---------------------------------
% Compare stdTheta1 (all params, PI #1) vs stdTheta2_sub (subset, PI #2)
% Expected result: much smaller std devs after PSS because the ill-posed
%                  directions have been removed.
plotCIComparison(thetaHat1, stdTheta1, thetaHat2, stdTheta2_sub, ...
                 keep_idx, p, 'Parameter uncertainty: PI #1 vs PI #2');

% ---- 11c. Print summary table -------------------------------------------
fprintf('\n=== RMSE Summary ===\n');
fprintf('%-30s  Train  |  CV\n', 'Configuration');
fprintf('%-30s  %.4f | %.4f\n', 'PI #1 (full set)',     RMSE1_auto, RMSE1_cv);
fprintf('%-30s  %.4f | %.4f\n', 'PI #2 (PSS subset)',   RMSE2_auto, RMSE2_cv);
