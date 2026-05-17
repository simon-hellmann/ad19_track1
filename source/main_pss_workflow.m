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
%   becomes the initial condition for the training (and CV) windows.
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

% --- Reactor geometry ----------------------------------------------------
V_liq = 0.012;   % liquid volume [m^3]
V_gas = 0.003;   % gas headspace volume [m^3]

% --- Load ADM1 fixed physico-chemical constants -------------------------
% Sourced from Soeren's model data (ad_para_ident reference project).
% Update this path if the reference project moves.
addpath('/Users/simonhellmann/Documents/GIT/ad_para_ident/source/Parameter Identification/Helper Functions/soerens_files/model_data');
load('ADM1_parameters.mat', 'parameters');
parameters_r3 = parameters.ADM1_R3;

p_atm    = 1.01325;        % atmospheric pressure [bar]
T_K      = 0 + 273.15;    % standard temperature [K]
p_h2o   = 0;               % dry biogas assumption [bar]

K_H_ch4  = parameters_r3{1,:};   K_H_co2  = parameters_r3{2,:};
K_S_IN   = parameters_r3{3,:};   K_I_nh3  = parameters_r3{4,:};
K_a_IN   = parameters_r3{5,:};   K_a_ac   = parameters_r3{6,:};
K_a_co2  = parameters_r3{7,:};   K_S_ac   = parameters_r3{8,:};
K_w      = parameters_r3{9,:};   R_gas    = parameters_r3{10,:};
k_AB_IN  = parameters_r3{12,:};  k_AB_ac  = parameters_r3{13,:};
k_AB_co2 = parameters_r3{14,:};  k_La     = parameters_r3{15,:};
k_p      = parameters_r3{20,:};
pK_l_ac  = parameters_r3{22,:};  pK_u_ac  = parameters_r3{23,:};

M_ch4 = 16;   M_co2 = 44;
n_ac  = 3 / (pK_u_ac - pK_l_ac);

% --- Assemble constant vector c (31 x 1) --------------------------------
c = nan(31, 1);
c(1)  = 1 / V_liq;
c(2)  = n_ac;
c(3)  = 10^(-(3/2) * (pK_u_ac + pK_l_ac) / (pK_u_ac - pK_l_ac));
c(4)  = 4 * K_w;
c(5)  = k_La;
c(6)  = k_La * K_H_ch4 * R_gas * T_K;
c(7)  = k_La * K_H_co2 * R_gas * T_K;
c(8)  = K_S_IN;
c(9)  = k_AB_ac;
c(10) = k_AB_co2;
c(11) = k_AB_IN;
c(12) = k_La * V_liq / V_gas;
c(13) = k_p / p_atm * (R_gas * T_K / M_ch4)^2;
c(14) = 2 * k_p / p_atm * (R_gas * T_K)^2 / M_ch4 / M_co2;
c(15) = k_p / p_atm * (R_gas * T_K / M_co2)^2;
c(16) = k_p / p_atm * R_gas * T_K / M_ch4 * (2*p_h2o - p_atm);
c(17) = k_p / p_atm * R_gas * T_K / M_co2 * (2*p_h2o - p_atm);
c(18) = k_p / p_atm * (p_h2o - p_atm) * p_h2o;
c(19) = R_gas * T_K / M_ch4;
c(20) = R_gas * T_K / M_co2;
c(21) = -k_p / V_gas / p_atm * (R_gas * T_K / M_ch4)^2;
c(22) = -2 * k_p / V_gas / p_atm * (R_gas * T_K)^2 / M_ch4 / M_co2;
c(23) = -k_p / V_gas / p_atm * (R_gas * T_K / M_co2)^2;
c(24) = -k_p / V_gas / p_atm * (R_gas * T_K / M_ch4) * (2*p_h2o - p_atm);
c(25) = -k_p / V_gas / p_atm * (R_gas * T_K / M_co2) * (2*p_h2o - p_atm);
c(26) = -k_La * V_liq / V_gas * K_H_ch4 * R_gas * T_K ...
        - k_p / V_gas / p_atm * (p_h2o - p_atm) * p_h2o;
c(27) = -k_La * V_liq / V_gas * K_H_co2 * R_gas * T_K ...
        - k_p / V_gas / p_atm * (p_h2o - p_atm) * p_h2o;
c(28) = k_AB_ac  * K_a_ac;
c(29) = k_AB_co2 * K_a_co2;
c(30) = k_AB_IN  * K_a_IN;
c(31) = V_liq / V_gas;

% --- Petersen stoichiometry matrix a (14 states x 11 reactions) ---------
% Rows: states; columns: reactions (ch, pr, li, dec, X_ac-dec, X_su-dec,
%   X_aa-dec, ac-methanogenesis, [rest gas/liq]).
% c(31) = V_liq/V_gas appears in the gas-phase rows.
a = [  0.6555,  0.081837,  0.2245,  -0.016932, -1,      0,      0,      0.11246, 0,  0,  0,  0,  0,       0;
       0.9947,  0.069636,  0.10291,  0.17456,   0,     -1,      0,      0.13486, 0,  0,  0,  0,  0,       0;
       1.7651,  0.19133,  -0.64716, -0.024406,  0,      0,     -1,      0.1621,  0,  0,  0,  0,  0,       0;
      -26.5447, 6.7367,   18.4808,  -0.15056,   0,      0,      0,      0,       1,  0,  0,  0,  0,       0;
       0,       0,         0,        0,          0.18,  0.77,   0.05,  -1,       0,  0,  0,  0,  0,       0;
       0,       0,         0,        0,          0.18,  0.77,   0.05,   0,      -1,  0,  0,  0,  0,       0;
       0,       0,         0,        0,          0,      0,      0,      0,       0, -1,  0,  0,  0,       0;
       0,       0,         0,        0,          0,      0,      0,      0,       0,  0, -1,  0,  0,       0;
       0,       0,         0,        0,          0,      0,      0,      0,       0,  0,  0, -1,  0,       0;
       0,      -1,         0,        0,          0,      0,      0,      0,       0,  0,  0,  0,  c(31),  0;
       0,       0,        -1,        0,          0,      0,      0,      0,       0,  0,  0,  0,  0,       c(31)]';

% --- Default inlet concentrations xi (14 x 1) ---------------------------
% Used as constant substrate composition for all feeding events.
xi = [7.42746759027778;
      0;
      0;
      1.03356401701299;
      129.110923498630;
      19.8236535343803;
      4.99887918689879;
      0.182468877424865;
      0.00960362512762449;
      3.02623121731264;
      0;
      0.0218481841639943;
      0;
      0];

% --- ODE and measurement function handles --------------------------------
% c, a, xi captured in closures; all helpers use the clean signatures:
%   odeFunc(x, u, theta)  ->  f (state derivative, n_states x 1)
%   measFunc(x, theta)    ->  g (6 outputs: gasflow, pCH4, pCO2, pH, SIN, Sac)
odeFunc  = @(x, u, theta) ADM1_R3_core_ode_sym_pi(x, u, xi, theta, c, a);
measFunc = @(x, theta)    ADM1_R3_core_mgl_sym_pi(x, theta, c);

% --- ODE solver options -------------------------------------------------
odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.5);

% --- Parameter vector: theta = [k_ch, k_pr, k_li, k_dec, k_m_ac, -------
%                                K_S_ac, K_I_nh3, DeltaS_ion, xi_SIN]  --
% Best abs-scale local estimate from ad_para_ident (Leander, IntBePro):
theta0 = [3.0035;                % k_ch      [1/d]
          0.0183;                % k_pr      [1/d]
          0.1842;                % k_li      [1/d]
          0.1446;                % k_dec     [1/d]
          1.3611;                % k_m_ac    [1/d]
          0.9558;                % K_S_ac    [g/L]
          0.0906;                % K_I_nh3   [g/L]
          0.3143 - 0.009125;     % DeltaS_ion [g/L]  (effective ion correction)
          1.0];                  % xi_SIN    [-]     (SIN inlet scaling)

thetaLB = [1e-3;  1e-3;  1e-3;  1e-4;  1e-2;  1e-3;  1e-3;  1e-4;  0.2];
thetaUB = [2e1;   2e1;   10;    5;     5;     1e1;   5;     1;     3  ];

% --- Output measurement noise standard deviations -----------------------
% One value per output; same sensor characteristics assumed for all windows.
% Units: [m^3/d, bar, bar, -, g/L, g/L]
sigmaY = [4e-4, 1.78e-2, 2.68e-2, 2e-2, 0.12, 5e-2];

% --- Meta struct for plot functions -------------------------------------
p.nParameters = numel(theta0);
p.names       = {'k_{ch}','k_{pr}','k_{li}','k_{dec}', ...
                  'k_{m,ac}','K_{S,ac}','K_{I,nh3}','DeltaS_{ion}','xi_{SIN}'};
p.units       = {'1/d','1/d','1/d','1/d','1/d','g/L','g/L','g/L','-'};
p.nOutputs    = 6;
p.outputNames = {'q_{gas}','p_{CH4}','p_{CO2}','pH','S_{IN}','S_{ac}'};
p.outputUnits = {'m^3/d','bar','bar','-','g/L','g/L'};

% --- Rough initial state and SS pre-simulation duration -----------------
x0_rough = xi;   % inlet concentrations give a reasonable biological operating point
t_ss     = 500;  % [d] steady-state pre-simulation duration for computeX0

%% -----------------------------------------------------------------------
%  3. LOAD & PREPARE DATASETS
% -----------------------------------------------------------------------

% Run split_data_pss.m first to generate these files from the raw MESS struct.
%
% load('../data/processed/data_init.mat',  'data_init');
% load('../data/processed/data_auto.mat',  'data_auto');
% load('../data/processed/data_cross.mat', 'data_cross');

% Unpack training (auto-validation) dataset:
%   tMeas{i}        -- (n_i x 1) measurement times per output [d]
%   yMeas{i}        -- (n_i x 1) measured values per output
%   t_feed_start    -- (n_events x 1) feed start times [d]
%   t_feed_end      -- (n_events x 1) feed end times [d]
%   u_feed_value    -- (n_events x 1) feed volume flow [m^3/d]
%   t0, tf          -- simulation window [d]
%
% tMeas        = data_auto.tMeas;
% yMeas        = data_auto.yMeas;
% t_feed_start = data_auto.t_feed_start;
% t_feed_end   = data_auto.t_feed_end;
% u_feed_value = data_auto.u_feed_value;
% t0           = data_auto.t0;
% tf           = data_auto.tf;
% n_out        = 6;

% Unpack cross-validation dataset (same fields, _CV suffix):
%
% tMeasCV        = data_cross.tMeas;
% yMeasCV        = data_cross.yMeas;
% t_feed_start_CV = data_cross.t_feed_start;
% t_feed_end_CV   = data_cross.t_feed_end;
% u_feed_value_CV = data_cross.u_feed_value;
% t0_CV           = data_cross.t0;
% tf_CV           = data_cross.tf;
% x0_CV          = [...];   % initial state for CV horizon (from state init)

% --- Measurement noise standard deviations (one per output) -------------
% sigmaY = [sigma_1, sigma_2, ..., sigma_nout];   % (1 x n_out)
% Same sensor noise assumed for both datasets; only sample counts differ.

% --- Build the long observation vector (training) -----------------------
% Entries are ordered by time, then by output index within each time point.
% At each measurement time only the outputs actually sampled there appear.
%
%   y_meas_long: (N_long x 1)  -- stacked measured values
%   t_meas_long: (N_long x 1)  -- corresponding time points
%   out_idx:     (N_long x 1)  -- output index (1..n_out) for each entry
%
% rows = [];
% for i = 1:n_out
%     ni   = numel(tMeas{i});
%     rows = [rows; tMeas{i}(:), i*ones(ni,1), yMeas{i}(:)];
% end
% rows        = sortrows(rows, [1, 2]);  % primary: time, secondary: output index
% t_meas_long = rows(:,1);
% out_idx     = rows(:,2);
% y_meas_long = rows(:,3);

% --- Scaling vector for the training long residual vector ---------------
% Each residual r_j is divided by sigma_i * sqrt(n_i), where i = out_idx(j)
% and n_i is the total number of samples for output i across the dataset.
%
% nSamples   = cellfun(@numel, tMeas);                               % (1 x n_out)
% scale_long = sigmaY(out_idx(:))' .* sqrt(nSamples(out_idx(:)))';  % (N_long x 1)

% --- Feeding event grid for training ------------------------------------
% t_events   = unique([t0; t_feed_start(:); t_feed_end(:); tf]);
% t_mid      = (t_events(1:end-1) + t_events(2:end)) / 2;
% u_segment  = zeros(numel(t_events)-1, 1);
% for i = 1:numel(t_feed_start)
%     active              = (t_mid >= t_feed_start(i)) & (t_mid <= t_feed_end(i));
%     u_segment(active)   = u_feed_value(i);
% end

% --- Build the long observation vector for cross-validation -------------
% Identical construction; uses CV measurement cell arrays.
%
% rowsCV        = [];
% for i = 1:n_out
%     ni        = numel(tMeasCV{i});
%     rowsCV    = [rowsCV; tMeasCV{i}(:), i*ones(ni,1), yMeasCV{i}(:)];
% end
% rowsCV           = sortrows(rowsCV, [1, 2]);
% t_meas_long_CV   = rowsCV(:,1);
% out_idx_CV       = rowsCV(:,2);
% y_meas_long_CV   = rowsCV(:,3);

% --- Scaling vector for the CV long residual vector ---------------------
% n_i is taken from the CV dataset, not the training dataset, so the
% normalisation reflects the actual CV sample density per output.
%
% nSamplesCV    = cellfun(@numel, tMeasCV);                                    % (1 x n_out)
% scale_long_CV = sigmaY(out_idx_CV(:))' .* sqrt(nSamplesCV(out_idx_CV(:)))'; % (N_long_CV x 1)

% --- Feeding event grid for cross-validation ----------------------------
% t_events_CV  = unique([t0_CV; t_feed_start_CV(:); t_feed_end_CV(:); tf_CV]);
% t_mid_CV     = (t_events_CV(1:end-1) + t_events_CV(2:end)) / 2;
% u_segment_CV = zeros(numel(t_events_CV)-1, 1);
% for i = 1:numel(t_feed_start_CV)
%     active                 = (t_mid_CV >= t_feed_start_CV(i)) & (t_mid_CV <= t_feed_end_CV(i));
%     u_segment_CV(active)   = u_feed_value_CV(i);
% end

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
%
% x0_rough = [...];   % (n_states x 1) rough initial state, e.g. from lit.
% t_ss     = 500;     % [d] pre-simulation duration for steady-state
%
% x0 = computeX0(theta0, data_init, t_ss, x0_rough, odeFunc, odeOpts);

%% -----------------------------------------------------------------------
%  5. PI #1 -- WEIGHTED LEAST SQUARES ON FULL PARAMETER VECTOR (fmincon)
% -----------------------------------------------------------------------

% Core simulation helper -- shared by costWLS and the FD sensitivity routine:
%
% See utils/simulateLong.m -- piecewise ODE integration then long-vector assembly.
% See utils/costWLS.m      -- thin wrapper: residuals + J = r'*r.
%
% objFun1 = @(theta) costWLS(theta, y_meas_long, t_meas_long, out_idx, ...
%               scale_long, t_events, u_segment, x0, odeFunc, measFunc, odeOpts);

% --- fmincon options ----------------------------------------------------
% opts1 = optimoptions('fmincon', ...
%     'Display',          'iter-detailed', ...
%     'Algorithm',        'interior-point', ...
%     'MaxFunctionEvaluations', 5000, ...
%     'OptimalityTolerance',    1e-8, ...
%     'StepTolerance',          1e-10);

% --- Run PI #1 ----------------------------------------------------------
% [thetaHat1, fval1, exitflag1, output1] = fmincon(objFun1, theta0, ...
%     [], [], [], [], thetaLB, thetaUB, [], opts1);

% thetaHat1: estimated parameter vector after PI #1

%% -----------------------------------------------------------------------
%  6. POST-PI #1 ANALYSIS
% -----------------------------------------------------------------------

% ---- 5a. Simulate with thetaHat1 and compute fit quality ---------------
% ySimLong1 has the same row ordering as y_meas_long (time-sorted).
%
% ySimLong1   = simulateLong(thetaHat1, t_meas_long, out_idx, t_events, ...
%                   u_segment, x0, odeFunc, measFunc, odeOpts);
% r_scaled1   = (ySimLong1 - y_meas_long) ./ scale_long;
% RMSE1_train = sqrt(mean(r_scaled1.^2));
% plotFit(t_meas_long, y_meas_long, ySimLong1, out_idx, p, 'PI #1 --training fit');

% ---- 5b. Scaled output sensitivity matrix via central finite differences -
%
% dydth1 has the same row ordering as the long observation vector:
%   rows: N_long entries, ordered by time then output index (same as y_meas_long)
%   cols: one per parameter
%   entry (j,k): d y_j/d theta_k * theta_k / scale_long(j)
%
% The relative parameter scaling (* theta_k) makes all columns dimensionless
% and comparable across parameters of different magnitudes and units.
% The / scale_long(j) row scaling is consistent with the cost function, so
% that FIM = dydth' * dydth directly gives the Fisher information in the
% same metric as the WLS cost.
%
% Step 1: compute raw (unscaled) FD sensitivities
%   dydth1_raw(j,k) = d y_j / d theta_k   [output units / parameter units]
%   rows follow the same ordering as y_meas_long (time-sorted, then out_idx)
%
% h_rel        = eps^(1/3);   % ~6e-6, optimal step for central differences
% N_long       = numel(t_meas_long);
% dydth1_raw   = zeros(N_long, p.nParameters);
%
% for k = 1:p.nParameters
%     delta_k       = h_rel * abs(thetaHat1(k));
%     if delta_k == 0;  delta_k = h_rel;  end   % guard for theta_k = 0
%
%     theta_fwd     = thetaHat1;  theta_fwd(k) = thetaHat1(k) + delta_k;
%     theta_bwd     = thetaHat1;  theta_bwd(k) = thetaHat1(k) - delta_k;
%
%     y_fwd = simulateLong(theta_fwd, t_meas_long, out_idx, t_events, ...
%                 u_segment, x0, odeFunc, measFunc, odeOpts);
%     y_bwd = simulateLong(theta_bwd, t_meas_long, out_idx, t_events, ...
%                 u_segment, x0, odeFunc, measFunc, odeOpts);
%
%     dydth1_raw(:,k) = (y_fwd - y_bwd) ./ (2*delta_k);
% end
%
% Step 2: scale rows and columns explicitly
%   Row scaling (output weight): divide row j by scale_long(j)
%     = sigma_{out_idx(j)} * sqrt(n_{out_idx(j)})
%     -- same denominator as the WLS cost, so the FIM is in the cost metric
%   Column scaling (parameter weight): multiply column k by thetaHat1(k)
%     -- a posteriori estimate from PI #1; makes columns dimensionless and
%        comparable across parameters of different magnitudes and units
%
% dydth1 = (dydth1_raw ./ scale_long) .* thetaHat1(:)';
%            ^-- (N_long x nP), broadcasts scale_long down each column
%                                ^-- broadcasts thetaHat1 across each row
%
% FIM1      = dydth1' * dydth1;          % Fisher Information Matrix
% C_theta1  = inv(FIM1);                 % Cramer-Rao lower bound
% stdTheta1 = sqrt(diag(C_theta1));      % std dev of parameter estimates
% plotUncertainty(thetaHat1, stdTheta1, p, 'PI #1 --parameter uncertainty');

% ---- 5c. Cross-validation with independent dataset ---------------------
% ySimLong1_CV = simulateLong(thetaHat1, t_meas_long_CV, out_idx_CV, ...
%                    t_events_CV, u_segment_CV, x0_CV, odeFunc, measFunc, odeOpts);
% r_scaled1_CV = (ySimLong1_CV - y_meas_long_CV) ./ scale_long_CV;
% RMSE1_cv     = sqrt(mean(r_scaled1_CV.^2));
% plotFit(t_meas_long_CV, y_meas_long_CV, ySimLong1_CV, out_idx_CV, p, ...
%         'PI #1 --cross-validation');

%% -----------------------------------------------------------------------
%  7. PARAMETER SUBSET SELECTION (PSS)
% -----------------------------------------------------------------------

% Uses the scaled sensitivity matrix from PI #1 as input.
% Applies SVD + QRP decomposition to identify the practically identifiable
% parameter subset (see subsetSelection.m -- Lopez et al., 2015).

% --- PSS thresholds (tune to your problem) ------------------------------
% kappa_max = 1e6;   % maximum allowable condition number
% gamma_max = 1e3;   % maximum allowable collinearity index

% --- Run PSS ------------------------------------------------------------
% Input is dydth1: the scaled sensitivity matrix from step 5b.
% Rows are scaled by 1/scale_long (output weight); columns by thetaHat1
% (a posteriori parameter estimate).  subsetSelection operates on this
% dimensionless matrix -- do NOT pass dydth1_raw here.
% [si, keep_idx, C_pp, pi_decomp, epsilon] = subsetSelection( ...
%     dydth1, kappa_max, gamma_max, p);

% keep_idx: indices of identifiable parameters
% The complement can be fixed at thetaHat1 or at literature values.

% --- Summarise PSS result -----------------------------------------------
% fprintf('\nPSS result: %d / %d parameters are identifiable.\n', ...
%         numel(keep_idx), p.nParameters);
% disp('Identifiable parameters:');
% disp(p.names(keep_idx));

% plotSingularValues(si, epsilon, p, 'PSS --singular value spectrum');
% plotVarianceDecomposition(pi_decomp, p, 'PSS --variance decomposition');

%% -----------------------------------------------------------------------
%  8. STATE INITIALISATION FOR PI #2
% -----------------------------------------------------------------------

% thetaHat1 may differ substantially from theta0, so x0 from section 4
% is no longer the best initial condition.  Repeat the same two-step
% procedure with thetaHat1 to obtain x0_2.
%
% x0_2 = computeX0(thetaHat1, data_init, t_ss, x0_rough, odeFunc, odeOpts);
%
% For cross-validation, data_cross immediately follows data_auto in time,
% so the terminal state of the auto simulation is used directly -- no
% separate steady-state pre-simulation is needed:
%
% x0_CV = xSol_auto_terminal;   % to be captured inside costWLS or stored
%                                % after the last simulateLong call in §6

%% -----------------------------------------------------------------------
%  9. PI #2 -- WEIGHTED LEAST SQUARES ON IDENTIFIABLE PARAMETER SUBSET
% -----------------------------------------------------------------------

% Non-identifiable parameters are fixed at thetaHat1 (their PI #1 values).
% The same long-vector cost function is reused; only the free variables change.

% --- Initialise fixed and free parameter vectors ------------------------
% thetaFixed  = thetaHat1;          % full vector; non-id entries stay fixed
% theta0_sub  = thetaHat1(keep_idx);
% thetaLB_sub = thetaLB(keep_idx);
% thetaUB_sub = thetaUB(keep_idx);

% --- Objective: embed the free sub-vector back into the full vector -----
% See utils/costWLS_sub.m
%
% objFun2 = @(theta_sub) costWLS_sub(theta_sub, keep_idx, thetaFixed, ...
%               y_meas_long, t_meas_long, out_idx, scale_long, ...
%               t_events, u_segment, x0, odeFunc, measFunc, odeOpts);

% --- Run PI #2 ----------------------------------------------------------
% [thetaSub2, fval2, exitflag2, output2] = fmincon(objFun2, theta0_sub, ...
%     [], [], [], [], thetaLB_sub, thetaUB_sub, [], opts1);

% --- Reconstruct the full parameter vector ------------------------------
% thetaHat2           = thetaFixed;
% thetaHat2(keep_idx) = thetaSub2;

%% -----------------------------------------------------------------------
%  10. POST-PI #2 ANALYSIS
% -----------------------------------------------------------------------

% ---- 8a. Fit quality ---------------------------------------------------
% ySimLong2   = simulateLong(thetaHat2, t_meas_long, out_idx, t_events, ...
%                   u_segment, x0, odeFunc, measFunc, odeOpts);
% r_scaled2   = (ySimLong2 - y_meas_long) ./ scale_long;
% RMSE2_train = sqrt(mean(r_scaled2.^2));
% plotFit(t_meas_long, y_meas_long, ySimLong2, out_idx, p, 'PI #2 --training fit');

% ---- 8b. Scaled output sensitivity matrix for the identifiable subset --
% Same FD procedure as 5b, but columns restricted to keep_idx.
% dydth2 has the same row structure as y_meas_long; it has numel(keep_idx)
% columns rather than p.nParameters.
%
% Step 1: raw FD sensitivities, columns restricted to keep_idx only
%
% dydth2_raw = zeros(N_long, numel(keep_idx));
%
% for ki = 1:numel(keep_idx)
%     k             = keep_idx(ki);
%     delta_k       = h_rel * abs(thetaHat2(k));
%     if delta_k == 0;  delta_k = h_rel;  end
%
%     theta_fwd     = thetaHat2;  theta_fwd(k) = thetaHat2(k) + delta_k;
%     theta_bwd     = thetaHat2;  theta_bwd(k) = thetaHat2(k) - delta_k;
%
%     y_fwd = simulateLong(theta_fwd, t_meas_long, out_idx, t_events, ...
%                 u_segment, x0, odeFunc, measFunc, odeOpts);
%     y_bwd = simulateLong(theta_bwd, t_meas_long, out_idx, t_events, ...
%                 u_segment, x0, odeFunc, measFunc, odeOpts);
%
%     dydth2_raw(:,ki) = (y_fwd - y_bwd) ./ (2*delta_k);
% end
%
% Step 2: apply output and parameter scaling
%   Row scaling: same output weights as the WLS cost (scale_long unchanged)
%   Column scaling: a posteriori estimates from PI #2 (thetaHat2(keep_idx))
%
% dydth2 = (dydth2_raw ./ scale_long) .* thetaHat2(keep_idx)';
%
% FIM2          = dydth2' * dydth2;
% C_theta2_sub  = inv(FIM2);
% stdTheta2_sub = sqrt(diag(C_theta2_sub));
% plotUncertainty(thetaHat2(keep_idx), stdTheta2_sub, p, ...
%                 'PI #2 --parameter uncertainty (identifiable subset)');

% ---- 8c. Cross-validation ----------------------------------------------
% ySimLong2_CV = simulateLong(thetaHat2, t_meas_long_CV, out_idx_CV, ...
%                    t_events_CV, u_segment_CV, x0_CV, odeFunc, measFunc, odeOpts);
% r_scaled2_CV = (ySimLong2_CV - y_meas_long_CV) ./ scale_long_CV;
% RMSE2_cv     = sqrt(mean(r_scaled2_CV.^2));
% plotFit(t_meas_long_CV, y_meas_long_CV, ySimLong2_CV, out_idx_CV, p, ...
%         'PI #2 --cross-validation');

%% -----------------------------------------------------------------------
%  11. COMPARISON: BEFORE vs. AFTER PSS
% -----------------------------------------------------------------------

% ---- 9a. Fit overlay (training + CV) -----------------------------------
% Show that ySimLong1 ~= ySimLong2 (fit quality preserved after PSS)
% plotFitComparison(t_meas_long, y_meas_long, ySimLong1, ySimLong2, out_idx, p, ...
%     'Fit comparison: PI #1 vs PI #2 (training)');
% plotFitComparison(t_meas_long_CV, y_meas_long_CV, ySimLong1_CV, ySimLong2_CV, ...
%     out_idx_CV, p, 'Fit comparison: PI #1 vs PI #2 (cross-validation)');

% ---- 9b. Confidence interval shrinkage ---------------------------------
% Compare stdTheta1 (all params, PI #1) vs stdTheta2_sub (subset, PI #2)
% Expected result: much smaller std devs after PSS because the ill-posed
%                  directions have been removed.
% plotCIComparison(thetaHat1, stdTheta1, thetaHat2, stdTheta2_sub, ...
%                  keep_idx, p, 'Parameter uncertainty: PI #1 vs PI #2');

% ---- 9c. Print summary table -------------------------------------------
% fprintf('\n=== RMSE Summary ===\n');
% fprintf('%-30s  Train  |  CV\n', 'Configuration');
% fprintf('%-30s  %.4f | %.4f\n', 'PI #1 (full set)',     RMSE1_train, RMSE1_cv);
% fprintf('%-30s  %.4f | %.4f\n', 'PI #2 (PSS subset)',   RMSE2_train, RMSE2_cv);
