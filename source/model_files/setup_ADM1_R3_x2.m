%% setup_ADM1_R3_x2.m
% Assemble all model-specific quantities for ADM1-R3-x2.
% Extends ADM1-R3 by promoting K_S_IN, k_La, and k_p from the fixed
% constant vector c into the tunable parameter vector theta.
%
% Author: Simon Hellmann. Created: 2026/05/24. Version: Matlab R2022b, Update 6
%
%% Output
%
%   c:          (19x1) time-invariant parameter vector
%   a:          (14x11) Petersen stoichiometry matrix
%   odeFunc:    @(x, u, xi, theta) state derivative  (wraps ADM1_R3_x2_core_ode)
%   measFunc:   @(x, theta) output vector            (wraps ADM1_R3_x2_core_output)
%   theta0:     (12x1) nominal tunable parameter vector
%   thetaLB:    (12x1) lower bounds for theta
%   thetaUB:    (12x1) upper bounds for theta
%   p:          metadata struct. Fields:
%               .nParameters    number of tunable parameters
%               .names          LaTeX parameter names  (1 x nP cell)
%               .units          parameter unit strings (1 x nP cell)
%               .nOutputs       number of model outputs
%               .outputNames    LaTeX output names     (1 x nO cell)
%               .outputUnits    output unit strings    (1 x nO cell)
%
%% Input
%
%   parameters_r3:  cell array of ADM1-R3 physico-chemical constants
%                   (from ADM1_parameters.mat, Weinrich 2017)
%   V_liq:          liquid volume                                [m^3]
%   V_gas:          gas headspace volume                         [m^3]
%
% Tunable parameters (theta):
%   th(1)   k_ch         hydrolysis rate, carbohydrates         [1/d]
%   th(2)   k_pr         hydrolysis rate, proteins              [1/d]
%   th(3)   k_li         hydrolysis rate, lipids                [1/d]
%   th(4)   k_dec        decay rate, all biomass                [1/d]
%   th(5)   k_m_ac       max. acetoclastic methanogenesis rate  [1/d]
%   th(6)   K_S_ac       half-saturation constant, acetate      [g/L]
%   th(7)   K_I_nh3      NH3 inhibition constant                [g/L]
%   th(8)   Delta_S_ion  effective ion correction (ion balance) [g/L]
%   th(9)   phi_IN       IN inlet scaling factor                [-]
%   th(10)  K_S_IN       N half-saturation, biomass growth      [g/L]
%   th(11)  k_La         liquid-gas mass transfer coefficient   [1/d]
%   th(12)  k_p          biogas extraction coefficient          [m^3/(bar*d)]
%

function [c, a, odeFunc, measFunc, theta0, thetaLB, thetaUB, p] = setup_ADM1_R3_x2( ...
    parameters_r3, V_liq, V_gas)

% --- Physico-chemical constants ------------------------------------------
p_atm = 1.01325;   % atmospheric pressure [bar]
T_K   = 273.15;    % standard temperature [K]
p_h2o = 0;         % water vapour partial pressure [bar] — 0: dry biogas assumption

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

% --- Assemble c (19 x 1) -------------------------------------------------
% K_S_IN, k_La, k_p are promoted to theta; c stores only the base scalar
% factors that the ODE needs to reconstruct k_La- and k_p-composite
% expressions from the current theta values at runtime.
% p_h2o and p_atm are stored explicitly so the gas ODE and q_gas output
% are correct for any p_h2o (including p_h2o > 0).
c = nan(19, 1);
c(1)  = 1 / V_liq;
c(2)  = n_ac;
c(3)  = 10^(-(3/2) * (pK_u_ac + pK_l_ac) / (pK_u_ac - pK_l_ac));
c(4)  = 4 * K_w;
c(5)  = K_H_ch4 * R_gas * T_K;   % Henry base const, CH4 (factor for k_La terms)
c(6)  = K_H_co2 * R_gas * T_K;   % Henry base const, CO2 (factor for k_La terms)
c(7)  = k_AB_ac;
c(8)  = k_AB_co2;
c(9)  = k_AB_IN;
c(10) = V_liq / V_gas;
c(11) = R_gas * T_K / M_ch4;     % RT_ch4 [bar*m^3/g]
c(12) = R_gas * T_K / M_co2;     % RT_co2 [bar*m^3/g]
c(13) = 1 / (V_gas * p_atm);     % base for k_p cubic/quadratic gas ODE terms
c(14) = k_AB_ac  * K_a_ac;
c(15) = k_AB_co2 * K_a_co2;
c(16) = k_AB_IN  * K_a_IN;
c(17) = 1 / p_atm;               % base for k_p in output q_gas
c(18) = p_h2o;                   % water vapour partial pressure [bar]
c(19) = p_atm;                   % atmospheric pressure [bar]

% --- Petersen stoichiometry matrix a (14 x 11) ---------------------------
a = [  0.6555,  0.081837,  0.2245,  -0.016932, -1,      0,      0,      0.11246, 0,  0,  0,  0,  0,        0;
       0.9947,  0.069636,  0.10291,  0.17456,   0,     -1,      0,      0.13486, 0,  0,  0,  0,  0,        0;
       1.7651,  0.19133,  -0.64716, -0.024406,  0,      0,     -1,      0.1621,  0,  0,  0,  0,  0,        0;
      -26.5447, 6.7367,   18.4808,  -0.15056,   0,      0,      0,      0,       1,  0,  0,  0,  0,        0;
       0,       0,         0,        0,          0.18,  0.77,   0.05,  -1,       0,  0,  0,  0,  0,        0;
       0,       0,         0,        0,          0.18,  0.77,   0.05,   0,      -1,  0,  0,  0,  0,        0;
       0,       0,         0,        0,          0,      0,      0,      0,      0, -1,  0,  0,  0,        0;
       0,       0,         0,        0,          0,      0,      0,      0,      0,  0, -1,  0,  0,        0;
       0,       0,         0,        0,          0,      0,      0,      0,      0,  0,  0, -1,  0,        0;
       0,      -1,         0,        0,          0,      0,      0,      0,      0,  0,  0,  0,  c(10),    0;
       0,       0,        -1,        0,          0,      0,      0,      0,      0,  0,  0,  0,  0,        c(10)]';

% --- ODE and measurement function handles --------------------------------
odeFunc  = @(x, u, xi, theta) ADM1_R3_x2_core_ode(x, u, xi, theta, c, a);
measFunc = @(x, theta)        ADM1_R3_x2_core_output(x, theta, c);

% --- Tunable parameter vector theta (12 x 1) -----------------------------
% th(1..9): identical to ADM1-R3.  th(10..12): promoted from parameters_r3.
theta0 = [parameters_r3{16,:};  % k_ch        [1/d]
          parameters_r3{21,:};  % k_pr        [1/d]
          parameters_r3{18,:};  % k_li        [1/d]
          parameters_r3{17,:};  % k_dec       [1/d]
          parameters_r3{19,:};  % k_m_ac      [1/d]
          parameters_r3{8,:};   % K_S_ac      [g/L]
          parameters_r3{4,:};   % K_I_nh3     [g/L]
          0.15 - 0.02;          % Delta_S_ion [g/L]
          1;                    % phi_IN      [-]
          parameters_r3{3,:};   % K_S_IN      [g/L]
          parameters_r3{15,:};  % k_La        [1/d]
          parameters_r3{20,:}]; % k_p         [m^3/(bar*d)]

thetaLB = [2e-2;  1e-2;  1e-2;  1e-3;  1e-1;  1e-3;  1e-2;  -1e-2;  1e-1;  1e-3;  1e1;  1e3];
thetaUB = [1e1;   1e1;   1;     1;     1;     2;     1e1;   1;      2;     1e2;   2e3;  1e7];

% --- Metadata struct -----------------------------------------------------
p.nParameters = numel(theta0);
p.names       = {'k_{ch}','k_{pr}','k_{li}','k_{dec}', ...
                  'k_{m,ac}','K_{S,ac}','K_{I,nh3}','\Delta S_{ion}','\phi_{IN}', ...
                  'K_{S,IN}','k_{La}','k_p'};
p.units       = {'1/d','1/d','1/d','1/d','1/d','g/L','g/L','mol/L','-','g/L','1/d','m^3/(bar*d)'};
p.nOutputs    = 6;
p.outputNames = {'q_{gas}','p_{CH4}','p_{CO2}','pH','S_{IN}','S_{ac}'};
p.outputUnits = {'m^3/d','bar','bar','-','g/L','g/L'};

end % fun
