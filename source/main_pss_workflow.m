%% main_pss_workflow.m
% AD19 workshop 'AD modeling toolbox': part 2, parallel session 3
% Parameter Identification -> PSS -> Parameter Identification
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

run_id = 3; % document in list of runs (user's responsibility)

% add subfolders (model files, helper functions, PSS routine, etc.)
addpath('model_files');
addpath('utils');  
addpath('model_files/model_data');

% random seed for reproducibility (only matters if using stochastic optimiser)
rng(42);

% --- Workflow flags -------------------------------------------------------
dataset          = 'automated_feeder';  % 'intensiv' | 'automated_feeder'
model_name       = "ADM1-R3-x1";% model variant: "ADM1-R3" | "ADM1-R3-x1" | "ADM1-R3-x2"
flag_skip_lhs    = true;        % true → skip LHS pre-search, start PI #1 from theta0
flag_omit_co2    = true;        % true → drop p_CO2 (output 3) from PI; q_gas + p_CH4 are enough
flag_plot_x0     = true;       % true → open SS + swing-up diagnostic plots from computeX0
dt_fine = 15/1440;              % [d] fine time-grid resolution for smooth output plots

%% -----------------------------------------------------------------------
%  2. MODEL DEFINITION & NOMINAL PARAMETERS
% -----------------------------------------------------------------------

% --- Reactor dimensions --------------------------------------------------
switch dataset
    case 'intensiv'
        V_liq = 0.012;   % liquid volume [m^3]
        V_gas = 0.003;   % gas headspace volume [m^3]
    case 'automated_feeder'
        V_liq = 0.011;   % liquid volume [m^3]
        V_gas = 0.004;   % gas headspace volume [m^3]
end

% --- Output measurement noise standard deviations -----------------------
% Units: [m^3/d, bar, bar, -, g/L, g/L]
switch dataset
    case 'intensiv'
        sigmaY = [4e-4, 1.78e-2, 2.68e-2, 2e-2, 0.12, 5e-2];
    case 'automated_feeder'
        sigmaY = [3e-4, 1e-3, NaN, 5e-3, 0.12, 5e-2]; % reduced pH noise, no CO2 measurement
end 

% --- Load ADM1 physico-chemical constants --------------------------------
load('ADM1_parameters.mat', 'parameters');
parameters_r3 = parameters.ADM1_R3;

% --- Assemble model -------------------------------------------------------
% Switch model_name in section 1. Add a case here for each new variant.
switch model_name
    case "ADM1-R3-x1"
        [c, a, odeFunc, measFunc, theta0, thetaLB, thetaUB, p] = ...
            setup_ADM1_R3_x1(parameters_r3, V_liq, V_gas);
    case "ADM1-R3-x2"
        [c, a, odeFunc, measFunc, theta0, thetaLB, thetaUB, p] = ...
            setup_ADM1_R3_x2(parameters_r3, V_liq, V_gas);
    case "ADM1-R3"
        [c, a, odeFunc, measFunc, theta0, thetaLB, thetaUB, p] = ...
            setup_ADM1_R3(parameters_r3, V_liq, V_gas);
    otherwise
        error("Unknown model_name '%s'. Add a case for new variants here.", model_name);
end % switch

% --- Log10 reparametrisation indices and phi-space bounds ----------------
% Strictly positive parameters are optimised in log10-space (phi_k = log10(theta_k)).
% This normalises all phi ~ O(1), improves L-BFGS Hessian conditioning, and
% makes the FD step a uniform ~1 % relative perturbation on theta regardless
% of parameter magnitude.  Delta_S_ion (param 8) has LB < 0, so it stays linear.
n_theta = numel(theta0);
log_idx = setdiff(1:n_theta, 8);   % strictly positive bounds → log10-scale
lin_idx = 8;                       % Delta_S_ion: LB < 0 → linear scale

phiLB = thetaLB;  phiLB(log_idx) = log10(thetaLB(log_idx));   % bounds in phi-space
phiUB = thetaUB;  phiUB(log_idx) = log10(thetaUB(log_idx));

% --- ODE solver options --------------------------------------------------
% Relaxed tolerances inside fmincon (10x looser than post-processing) to speed
% up cost evaluations. NonNegative prevents ode15s driving biological
% concentration states below zero.
non_negative_state_idx = 1:14;
odeOptsOpt  = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', 0.5/24, ...
                     'NonNegative', non_negative_state_idx);
% Tight tolerances for post-processing (FD sensitivity, CV, plots):
odeOptsPost = odeset('RelTol', 1e-8, 'AbsTol', 1e-9, 'MaxStep', 0.5/24, ...
                     'NonNegative', non_negative_state_idx);

% --- LHS pre-screening (only valid if flag_skip_lhs=false) --------------
N_LHS = 100;

%% -----------------------------------------------------------------------
%  3. LOAD & PREPROCESS DATASETS
% -----------------------------------------------------------------------

% Note: Run prepare_data.m first (with the same dataset flag) to generate these files.
%
proc_dir = fullfile('..', 'data', 'processed', dataset);
load(fullfile(proc_dir, 'data_init.mat'),        'data_init');
load(fullfile(proc_dir, 'data_auto.mat'),        'data_auto');
load(fullfile(proc_dir, 'data_cross.mat'),       'data_cross');
load(fullfile(proc_dir, 'feeding_duration.mat'), 'feeding_duration');

% --- Derive feed flow rates -----------------------------------------------
rho_feed = 1000;   % [kg/m^3] substrate density

% Unpack training (auto-validation) dataset:
%   tMeas{i}        -- (n_i x 1) measurement times per output [d]
%   yMeas{i}        -- (n_i x 1) measured values per output
%   t_feed_start    -- (n_events x 1) feed start times [d]
%   t_feed_end      -- derived: t_feed_start + feeding_duration [d]
%   u_feed_value    -- derived: feed_mass / (rho_feed * feeding_duration) [m^3/d]
%   t0, tf          -- simulation window [d]
tMeas        = data_auto.tMeas;
yMeas        = data_auto.yMeas;
t_feed_start = data_auto.t_feed_start;
t_feed_end   = t_feed_start + feeding_duration;
u_feed_value = data_auto.feed_mass ./ (rho_feed * feeding_duration);
t0           = data_auto.t0;
tf           = data_auto.tf;
n_out        = 6; % # model output channels

% Unpack cross-validation dataset (same derivation, _CV suffix):
tMeasCV         = data_cross.tMeas;
yMeasCV         = data_cross.yMeas;
t_feed_start_CV = data_cross.t_feed_start;
t_feed_end_CV   = t_feed_start_CV + feeding_duration;
u_feed_value_CV = data_cross.feed_mass ./ (rho_feed * feeding_duration);
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
nSamples   = cellfun(@numel, tMeas(:));                          % (n_out x 1) always column
scale_long = sigmaY(out_idx(:))' .* sqrt(nSamples(out_idx(:)));  % (N_long x 1)

% --- Feeding event grid for autovalidation ------------------------------
t_events   = unique([t0; t_feed_start(:); t_feed_end(:); tf]);
t_mid      = (t_events(1:end-1) + t_events(2:end)) / 2;
n_xi        = size(data_auto.xi_feed, 2);
u_segments  = zeros(numel(t_events)-1, 1);
xi_segments = zeros(numel(t_events)-1, n_xi);
for event_k = 1:numel(t_feed_start)
    mask_active = (t_mid >= t_feed_start(event_k)) & (t_mid <= t_feed_end(event_k));
    u_segments(mask_active) = u_feed_value(event_k);
    xi_segments(mask_active,:) = repmat(data_auto.xi_feed(event_k,:), sum(mask_active), 1);
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

nSamplesCV    = cellfun(@numel, tMeasCV(:));                               % (n_out x 1) always column
scale_long_CV = sigmaY(out_idx_CV(:))' .* sqrt(nSamplesCV(out_idx_CV(:)));   % (N_long_CV x 1)

% --- Feeding event grid for cross-validation ----------------------------
t_events_CV  = unique([t0_CV; t_feed_start_CV(:); t_feed_end_CV(:); tf_CV]);
t_mid_CV     = (t_events_CV(1:end-1) + t_events_CV(2:end)) / 2;
u_segments_CV  = zeros(numel(t_events_CV)-1, 1);
xi_segments_CV = zeros(numel(t_events_CV)-1, n_xi);
for event_k = 1:numel(t_feed_start_CV)
    mask_active = (t_mid_CV >= t_feed_start_CV(event_k)) & (t_mid_CV <= t_feed_end_CV(event_k));
    u_segments_CV(mask_active)   = u_feed_value_CV(event_k);
    xi_segments_CV(mask_active,:) = repmat(data_cross.xi_feed(event_k,:), sum(mask_active), 1);
end

% --- Optionally drop p_CO2 (output 3) from both PI datasets --------------
% q_gas, p_CH4, p_CO2 are linearly dependent via the gas-phase mole balance,
% so all three together over-constrain the problem; two suffice.
idx_co2 = 3;
if flag_omit_co2
    keep          = out_idx    ~= idx_co2;
    t_meas_long   = t_meas_long(keep);
    out_idx       = out_idx(keep);
    y_meas_long   = y_meas_long(keep);
    scale_long    = scale_long(keep);

    keep_CV         = out_idx_CV ~= idx_co2;
    t_meas_long_CV  = t_meas_long_CV(keep_CV);
    out_idx_CV      = out_idx_CV(keep_CV);
    y_meas_long_CV  = y_meas_long_CV(keep_CV);
    scale_long_CV   = scale_long_CV(keep_CV);

    fprintf("flag_omit_co2 = true: p_CO2 removed. " + ...
            "%d auto-val samples, %d CV samples remain.\n", ...
            numel(t_meas_long), numel(t_meas_long_CV));
end

plotMeasurements(t_meas_long,    y_meas_long,    out_idx,    p, t_events,    u_segments,    'Measurement data -- autovalidation');
plotMeasurements(t_meas_long_CV, y_meas_long_CV, out_idx_CV, p, t_events_CV, u_segments_CV, 'Measurement data -- cross-validation');

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
x0 = computeX0(theta0, data_init, t_ss, x0_init, odeFunc, odeOptsOpt, ...
    feeding_duration, flag_plot_x0, measFunc);

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

J_penalty   = 1e10;
maxNumCores = feature('numCores');

if flag_skip_lhs
    disp("LHS skipped — will start PI #1 from theta0.")
else
    % LHS design in [0,1]^n_theta
    lhs_unit  = lhsdesign(N_LHS, n_theta);   % (N_LHS x n_theta)

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
    fprintf("Evaluating cost at %d LHS samples (parfor)...\n", N_LHS);
    J_lhs = nan(N_LHS, 1);
    if isempty(gcp('nocreate'))
        parpool(maxNumCores, 'IdleTimeout',Inf); % pool that never loses connection
    end
    tic_lhs = tic;
    parfor i_lhs = 1:N_LHS
        J_lhs(i_lhs) = costWLS(theta_lhs(i_lhs, :)', y_meas_long, t_meas_long, ...
            out_idx, scale_long, t_events, u_segments, xi_segments, x0, ...
            odeFunc, measFunc, odeOptsOpt, J_penalty);
    end
    computing_time_lhs = toc(tic_lhs); % [s]
    fprintf("LHS done.  Cost range: [%.4g, %.4g]  (%.1f s)\n", min(J_lhs), max(J_lhs), computing_time_lhs);
    delete(gcp('nocreate'));

    % Best LHS candidate becomes the starting point for PI #1
    [J0_lhs, lhs_best_idx] = min(J_lhs);
    theta_lhs_best    = theta_lhs(lhs_best_idx, :)';
    fprintf("Best LHS candidate: sample #%d,  J = %.4g\n", lhs_best_idx, J0_lhs);
    theta0 = theta_lhs_best;
end % if flag_skip_lhs

% Starting point for fmincon in phi-space (all phi ~ O(1) after log10 transform)
phi0          = theta0;
phi0(log_idx) = log10(theta0(log_idx));

%% -----------------------------------------------------------------------
% sanity-checks of initial param
% -----------------------------------------------------------------------

% Evaluate cost at theta0 once (with optimization tolerances) to normalize J.
% Dividing by J0 brings the cost to O(1) at the start, which keeps the
% gradient components in a numerically tractable range for fmincon's FD.
% The optimum location is unchanged; only the gradient magnitude is rescaled.
[J0, J_ch0, r_scaled_cell0] = costWLS(theta0, y_meas_long, ...
        t_meas_long, out_idx, scale_long, t_events, u_segments, ...
        xi_segments, x0, odeFunc, measFunc, odeOptsOpt, J_penalty);
fprintf('J(theta0) = %.4g  (normalization factor)\n', J0);

[~, ~, t_fine0, ~, y_fine0] = simulateLong(theta0, t_meas_long, out_idx, t_events, ...
    u_segments, xi_segments, x0, odeFunc, measFunc, odeOptsPost, dt_fine);

plotCostChannels(J_ch0, p, 'Cost channels at init param theta0')
plotResidualBoxplot(r_scaled_cell0, p, 'Scaled residual at init param theta0')
plotFit(t_meas_long, y_meas_long, t_fine0, y_fine0, out_idx, p, t_events, ...
    u_segments, 'Fit at init param theta0')

%% -----------------------------------------------------------------------
%  5. PI #1 -- MAX LIKELIHOOD ESTIMATION ON FULL PARAMETER VECTOR (fmincon)
% -----------------------------------------------------------------------

% Core simulation helper -- shared by costWLS and the FD sensitivity routine:
%
% See utils/simulateLong.m -- piecewise ODE integration then long-vector assembly.
% See utils/costWLS.m      -- thin wrapper: residuals + J = r'*r.

objFun1 = @(phi) costWLS(phi2theta(phi, log_idx), y_meas_long, t_meas_long, out_idx, ...
              scale_long, t_events, u_segments, xi_segments, x0, ...
              odeFunc, measFunc, odeOptsOpt, J_penalty) / J0;

% --- fmincon options ----------------------------------------------------
% FD step in phi-space: h ~ ODE_noise^(1/3) balances truncation error O(h^2)
% against gradient noise O(ODE_noise/h).  All phi ~ O(1) so a uniform absolute
% step is appropriate; no per-parameter relative scaling needed.  RelTol=1e-7 → h ≈ 4.6e-3.
fd_stepSize = odeOptsOpt.RelTol^(1/3);  % optimal FD step size
gradient_noise_floor = odeOptsOpt.RelTol/fd_stepSize; 
opts1 = optimoptions('fmincon', ...
    'Display',                  'iter-detailed', ...
    'Algorithm',                'interior-point', ...
    'UseParallel',              true, ...
    'FiniteDifferenceType',     'central', ...
    'FiniteDifferenceStepSize', fd_stepSize, ...
    'HessianApproximation',     'lbfgs', ...    % default 10-pair memory; 4 pairs caused cost
    ...                                     % to increase in early iterations (Hessian too coarse)
    'StepTolerance',            1e-10, ...  % MATLAB default: larger values fire SOONER
    ...                                     % (step must shrink BELOW threshold to trigger)
    'OptimalityTolerance',      20*gradient_noise_floor, ...   % above gradient noise floor:
    ...                                     %   only stop when slope is genuinely flat
    'TypicalX',                 ones(n_theta, 1), ...   % phi ~ O(1) for all params
    'MaxFunctionEvaluations',   5000, ...
    'MaxIterations',            800);

% --- Run PI #1 (starting from best LHS candidate) -----------------------
disp("Running PI1 with fmincon...")
tic_pi1 = tic;
[phiHat1, fval1_norm, exitflag1, output1] = fmincon(objFun1, phi0, ...
    [], [], [], [], phiLB, phiUB, [], opts1);
computing_time_pi1 = toc(tic_pi1); % [s]
thetaHat1 = phi2theta(phiHat1, log_idx);   % convert phi-space result to physical params

%% -----------------------------------------------------------------------
%  6. POST-PI #1 ANALYSIS
% -----------------------------------------------------------------------

% ---- 5a. Simulate with thetaHat1 and compute fit quality ---------------
% ySimLong1 has the same row ordering as y_meas_long (time-sorted).

[ySimLong1, x0_CV, t_fine1, ~, y_fine1] = simulateLong(thetaHat1, t_meas_long, ...
    out_idx, t_events, u_segments, xi_segments, x0, odeFunc, measFunc, ...
    odeOptsPost, dt_fine);
r_scaled1  = (ySimLong1 - y_meas_long) ./ scale_long;
RMSE1_auto = sqrt(mean(r_scaled1.^2));
plotFit(t_meas_long, y_meas_long, t_fine1, y_fine1, out_idx, p, t_events, ...
    u_segments, 'PI #1 --autovalidation fit');

[J1, J_ch1, r_scaled_cell1] = costWLS(thetaHat1, y_meas_long, ...
        t_meas_long, out_idx, scale_long, t_events, u_segments, ...
        xi_segments, x0, odeFunc, measFunc, odeOptsOpt, J_penalty);
plotCostChannels(J_ch1, p, 'Cost channels PI #1')
plotResidualBoxplot(r_scaled_cell1, p, 'Scaled residual PI #1')

% ---- 5b. Scaled output sensitivity matrix via central finite differences in phi-space --
%
% A uniform absolute step h_phi is used for all parameters.
% For log-scaled params (log_idx): step h_phi in phi = log10(theta) is equivalent
% to multiplying theta_k by 10^(±h_phi), i.e. a ~0.5 % relative step on theta
% regardless of parameter magnitude.
% For the linear param (lin_idx): the same h_phi is an absolute step on theta_8 directly.
%
% dydphi1_os(j,k) = (d y_j / d phi_k) / scale_long(j)   [dimensionless]
%
% This single matrix serves both FIM and PSS — no second matrix needed.
% Relation to the old theta-space matrices for log-scaled columns:
%   dydth1_ops(:,k)  = (dy/dtheta_k) * thetaHat1(k) / scale_j   [old, param-scaled]
%   dydphi1_os(:,k)  = dydth1_ops(:,k) * ln(10)                  [phi-space, log_idx]
%   dydphi1_os(:,k)  = dydth1_os(:,k)                            [phi-space, lin_idx]
% The factor ln(10) is a uniform column scale that does not change PSS collinearity
% detection but makes the FIM live in the correct phi-space (log10-unit) coordinates.

h_phi  = odeOptsPost.RelTol^(1/3);   % ≈ 2.15e-3 for RelTol=1e-8; uniform in phi-space
N_long = numel(t_meas_long);
dydphi1_raw = nan(N_long, n_theta);

disp("Computing FD sensitivities in phi-space (parfor)...")
if isempty(gcp('nocreate'))
    parpool(maxNumCores,'IdleTimeout',Inf); % pool that never loses connection
end

parfor k = 1:n_theta
    fprintf("Computing FD sensitivities of param %i of %i...\n", k, n_theta);
    phi_fwd    = phiHat1;  phi_fwd(k)  = phiHat1(k) + h_phi;
    phi_bwd    = phiHat1;  phi_bwd(k)  = phiHat1(k) - h_phi;
    theta_fwd  = phi2theta(phi_fwd, log_idx);
    theta_bwd  = phi2theta(phi_bwd, log_idx);

    y_fwd = simulateLong(theta_fwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0, odeFunc, measFunc, odeOptsPost);
    y_bwd = simulateLong(theta_bwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0, odeFunc, measFunc, odeOptsPost);

    dydphi1_raw(:,k) = (y_fwd - y_bwd) ./ (2*h_phi);
end
delete(gcp('nocreate'))

dydphi1_os = dydphi1_raw ./ scale_long;   % (N_long x n_theta): output-scaled phi-space sensitivity

FIM1_phi = dydphi1_os' * dydphi1_os;      % Fisher information in phi-space
C_phi1   = inv(FIM1_phi);
std_phi1 = sqrt(diag(C_phi1));            % std dev of phiHat1 [log10 units for log-params]

% ±1-sigma confidence interval bounds in theta-space
% log-params: multiplicative  theta * 10^(±std_phi)
% lin-param:  additive        theta ± std_phi
theta1_lo              = thetaHat1;
theta1_hi              = thetaHat1;
theta1_lo(log_idx)     = thetaHat1(log_idx) .* 10.^(-std_phi1(log_idx));
theta1_hi(log_idx)     = thetaHat1(log_idx) .* 10.^(+std_phi1(log_idx));
theta1_lo(lin_idx)     = thetaHat1(lin_idx) - std_phi1(lin_idx);
theta1_hi(lin_idx)     = thetaHat1(lin_idx) + std_phi1(lin_idx);

% delta-method std in theta-space for plotting (first-order approx, valid when std_phi < 0.5)
std_theta1             = thetaHat1 .* log(10) .* std_phi1;
std_theta1(lin_idx)    = std_phi1(lin_idx);
plotUncertainty(thetaHat1, std_theta1, p, 'PI #1 --parameter uncertainty');

% ---- 5c. Cross-validation with independent dataset ---------------------
[ySimLong1_CV, ~, t_fine1_CV, ~, y_fine1_CV] = simulateLong(thetaHat1, t_meas_long_CV, ...
    out_idx_CV, t_events_CV, u_segments_CV, xi_segments_CV, x0_CV, odeFunc, measFunc, odeOptsPost, dt_fine);
r_scaled1_CV = (ySimLong1_CV - y_meas_long_CV) ./ scale_long_CV;
RMSE1_cv     = sqrt(mean(r_scaled1_CV.^2));
plotFit(t_meas_long_CV, y_meas_long_CV, t_fine1_CV, y_fine1_CV, out_idx_CV, p, ...
    t_events_CV, u_segments_CV, 'PI #1 --cross-validation');

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
% dydphi1_os is output-scaled in phi-space; equivalent to output-and-
% parameter-scaled dydth1_ops up to a factor ln(10) on log-param columns.
% That uniform column scale does not change the SVD collinearity structure.
% C_pp from PSS is suppressed: the function assumes dydth_ops (double-scaled 
% by outputs and thetaHat), but we pass dydphi1_os (scaled by thetaHat*ln(10)
% for log-params), so its D_theta back-transformation produces wrong units. 
% C_phi1 (§6b) is the correct covariance.
[si, keep_idx, ~, pi_decomp, epsilon, kappa, gamma] = param_subset_selection(...
    dydphi1_os, kappa_max, gamma_max, p, thetaHat1);

% keep_idx: indices of identifiable parameters
% The complement should be fixed at thetaHat1 or at literature values.

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
% is no longer the best initial condition.  Use x0 (PI #1 warm start, on
% the normal operating branch) with t_ss = 0 to skip the SS pre-simulation.
% This avoids the sour-SS bifurcation that can occur when the 500 d
% pre-simulation is run with thetaHat1 under the automated-feeder feed rate.

x0_2 = computeX0(thetaHat1, data_init, 0, x0_CV, odeFunc, odeOptsOpt, ...
    feeding_duration, flag_plot_x0, measFunc);

% For cross-validation, data_cross immediately follows data_auto in time,
% so x0_CV is the terminal state of the auto simulation (captured in §6).
% x0_CV_2 (self-consistent with thetaHat2) is computed in §10 after PI #2.

%% -----------------------------------------------------------------------
%  9. PI #2 -- MAX LIKELIHOOD ESTIMATION ON IDENTIFIABLE PARAMETER SUBSET
% -----------------------------------------------------------------------

% Non-identifiable parameters are fixed at thetaHat1 (their PI #1 values).
% The same long-vector cost function is reused; only the free variables change.

% --- Initialise fixed and free parameter vectors in phi-space -----------
thetaFixed    = thetaHat1;             % full vector; non-id params stay fixed at PI#1 result
keep_log_mask = ismember(keep_idx, log_idx);   % logical: which subset entries are log-scaled

phi1_sub  = phiHat1(keep_idx);         % phi starting point restricted to identifiable subset
phiLB_sub = phiLB(keep_idx);
phiUB_sub = phiUB(keep_idx);

% --- Objective: phi_sub -> theta_sub -> embed in thetaFixed -> WLS cost -
% See utils/costWLS_sub.m and utils/phi2theta.m
objFun2 = @(phi_sub) costWLS_sub(phi2theta(phi_sub, keep_log_mask), keep_idx, thetaFixed, ...
              y_meas_long, t_meas_long, out_idx, scale_long, ...
              t_events, u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsOpt, J_penalty) / J0;

% --- fmincon options for PI #2 (all phi_sub ~ O(1)) --------------------
opts2 = optimoptions(opts1, 'TypicalX',ones(numel(keep_idx), 1));

% --- Run PI #2 ----------------------------------------------------------
disp("Running PI2 with fmincon...")
tic_pi2 = tic;
[phiSub2, fval2, exitflag2, output2] = fmincon(objFun2, phi1_sub, ...
    [], [], [], [], phiLB_sub, phiUB_sub, [], opts2);
computing_time_pi2 = toc(tic_pi2); % [s]

% --- Reconstruct full parameter vector in theta-space ------------------
thetaSub2           = phi2theta(phiSub2, keep_log_mask);
thetaHat2           = thetaFixed;
thetaHat2(keep_idx) = thetaSub2;

%% -----------------------------------------------------------------------
%  10. POST-PI #2 ANALYSIS
% -----------------------------------------------------------------------

% ---- 10a. Fit quality ---------------------------------------------------
[ySimLong2, x0_CV_2, t_fine2, ~, y_fine2] = simulateLong(thetaHat2, t_meas_long, out_idx, ...
    t_events, u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsPost, dt_fine);
r_scaled2  = (ySimLong2 - y_meas_long) ./ scale_long;
RMSE2_auto = sqrt(mean(r_scaled2.^2));
plotFit(t_meas_long, y_meas_long, t_fine2, y_fine2, out_idx, p, t_events, u_segments, 'PI #2 --autovalidation fit');

% ---- 10b. Scaled output sensitivity in phi-space for the identifiable subset --
% Same structure as §6b, columns restricted to keep_idx.
% keep_log_mask identifies which subset entries are log-scaled.

dydphi2_raw = nan(N_long, numel(keep_idx));

disp("Computing FD sensitivities in phi-space (parfor)...")
if isempty(gcp('nocreate'))
    parpool(maxNumCores,'IdleTimeout',Inf); % pool that never loses connection
end

parfor ki = 1:numel(keep_idx)
    fprintf("Computing FD sensitivities of param %i of %i...\n", ki, numel(keep_idx));
    phi_fwd_sub   = phiSub2;  phi_fwd_sub(ki) = phiSub2(ki) + h_phi;
    phi_bwd_sub   = phiSub2;  phi_bwd_sub(ki) = phiSub2(ki) - h_phi;
    theta_fwd_sub = phi2theta(phi_fwd_sub, keep_log_mask);
    theta_bwd_sub = phi2theta(phi_bwd_sub, keep_log_mask);

    theta_fwd = thetaFixed;  theta_fwd(keep_idx) = theta_fwd_sub;
    theta_bwd = thetaFixed;  theta_bwd(keep_idx) = theta_bwd_sub;

    y_fwd = simulateLong(theta_fwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsPost);
    y_bwd = simulateLong(theta_bwd, t_meas_long, out_idx, t_events, ...
                u_segments, xi_segments, x0_2, odeFunc, measFunc, odeOptsPost);

    dydphi2_raw(:,ki) = (y_fwd - y_bwd) ./ (2*h_phi);
end
delete(gcp('nocreate'))

dydphi2_os = dydphi2_raw ./ scale_long;   % (N_long x numel(keep_idx)): output-scaled phi-space

FIM2_phi = dydphi2_os' * dydphi2_os;
C_phi2   = inv(FIM2_phi);
std_phi2 = sqrt(diag(C_phi2));            % std dev of phiSub2 [log10 units for log-params]

% ±1-sigma confidence interval bounds in theta-space
thetaSub2_lo              = thetaSub2;
thetaSub2_hi              = thetaSub2;
thetaSub2_lo(keep_log_mask)  = thetaSub2(keep_log_mask)  .* 10.^(-std_phi2(keep_log_mask));
thetaSub2_hi(keep_log_mask)  = thetaSub2(keep_log_mask)  .* 10.^(+std_phi2(keep_log_mask));
thetaSub2_lo(~keep_log_mask) = thetaSub2(~keep_log_mask) - std_phi2(~keep_log_mask);
thetaSub2_hi(~keep_log_mask) = thetaSub2(~keep_log_mask) + std_phi2(~keep_log_mask);

% delta-method std in theta-space for plotting
std_theta2_sub                  = thetaSub2 .* log(10) .* std_phi2;
std_theta2_sub(~keep_log_mask)  = std_phi2(~keep_log_mask);
std_theta2_full                 = NaN(n_theta, 1);
std_theta2_full(keep_idx)       = std_theta2_sub;
plotUncertainty(thetaHat2, std_theta2_full, p, ...
                'PI #2 --parameter uncertainty (identifiable subset)');

% ---- 10c. Cross-validation ----------------------------------------------
[ySimLong2_CV, ~, t_fine2_CV, ~, y_fine2_CV] = simulateLong(thetaHat2, t_meas_long_CV, ...
    out_idx_CV, t_events_CV, u_segments_CV, xi_segments_CV, x0_CV_2, odeFunc, measFunc, odeOptsPost, dt_fine);
r_scaled2_CV = (ySimLong2_CV - y_meas_long_CV) ./ scale_long_CV;
RMSE2_cv     = sqrt(mean(r_scaled2_CV.^2));
plotFit(t_meas_long_CV, y_meas_long_CV, t_fine2_CV, y_fine2_CV, out_idx_CV, p, ...
    t_events_CV, u_segments_CV, 'PI #2 --cross-validation');

%% -----------------------------------------------------------------------
%  11. COMPARISON: BEFORE vs. AFTER PSS
% -----------------------------------------------------------------------

% ---- 11a. Fit overlay (autovalidation + CV) -----------------------------------
% Show that ySimLong1 ~= ySimLong2 (fit quality preserved after PSS)
plotFitComparison(t_meas_long, y_meas_long, t_fine1, y_fine1, y_fine2, out_idx, p, ...
    t_events, u_segments, 'Fit comparison: PI #1 vs PI #2 (autovalidation)');
plotFitComparison(t_meas_long_CV, y_meas_long_CV, t_fine1_CV, y_fine1_CV, y_fine2_CV, ...
    out_idx_CV, p, t_events_CV, u_segments_CV, 'Fit comparison: PI #1 vs PI #2 (cross-validation)');

% ---- 11b. Confidence interval shrinkage ---------------------------------
% Compare std_theta1 (all params, PI #1) vs std_theta2_sub (subset, PI #2)
% Expected result: much smaller std devs after PSS because the ill-posed
%                  directions have been removed.
plotCIComparison(thetaHat1, std_theta1, thetaHat2, std_theta2_sub, ...
                 keep_idx, p, 'Parameter uncertainty: PI #1 vs PI #2');

% ---- 11c. Print summary table -------------------------------------------
fprintf('\n=== RMSE Summary ===\n');
fprintf('%-30s  Train  |  CV\n', 'Configuration');
fprintf('%-30s  %.4f | %.4f\n', 'PI #1 (full set)',     RMSE1_auto, RMSE1_cv);
fprintf('%-30s  %.4f | %.4f\n', 'PI #2 (PSS subset)',   RMSE2_auto, RMSE2_cv);

% save workspace: 
save(sprintf(fullfile('..', 'data', 'generated', 'workspace_run%i.mat'), run_id))
 