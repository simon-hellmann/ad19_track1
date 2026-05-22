function [c, a, odeFunc, measFunc, theta0, thetaLB, thetaUB, p] = setup_ADM1_R3(parameters_r3, V_liq, V_gas)
% Assemble all model-specific quantities for ADM1-R3.
%
% Inputs:
%   parameters_r3  -- cell array of ADM1-R3 physico-chemical constants
%                     (from ADM1_parameters.mat, Weinrich 2017)
%   V_liq          -- liquid volume [m^3]
%   V_gas          -- gas headspace volume [m^3]
%
% Outputs:
%   c              -- (31x1) time-invariant parameter vector
%   a              -- (14x11) Petersen stoichiometry matrix
%   odeFunc        -- @(x,u,xi,theta) state derivative
%   measFunc       -- @(x,theta) output vector (6 outputs)
%   theta0         -- (9x1) nominal tunable parameter vector
%   thetaLB        -- (9x1) lower bounds for theta
%   thetaUB        -- (9x1) upper bounds for theta
%   p              -- metadata struct (names, units, counts)
%
% Tunable parameters (theta):
%   th(1) k_ch        [1/d]   -- hydrolysis rate, carbohydrates
%   th(2) k_pr        [1/d]   -- hydrolysis rate, proteins
%   th(3) k_li        [1/d]   -- hydrolysis rate, lipids
%   th(4) k_dec       [1/d]   -- decay rate (all biomass)
%   th(5) k_m_ac      [1/d]   -- max. acetoclastic methanogenesis rate
%   th(6) K_S_ac      [g/L]   -- half-saturation constant, acetate
%   th(7) K_I_nh3     [g/L]   -- NH3 inhibition constant
%   th(8) Delta_S_ion [g/L]   -- effective ion correction (ion balance)
%   th(9) phi_IN      [-]     -- IN inlet scaling factor

% --- Physico-chemical constants ------------------------------------------
p_atm    = 1.01325;   % atmospheric pressure [bar]
T_K      = 273.15;    % standard temperature [K]
p_h2o    = 0;         % dry biogas assumption [bar]

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

% --- Assemble c (31 x 1) -------------------------------------------------
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

% --- Petersen stoichiometry matrix a (14 x 11) ---------------------------
a = [  0.6555,  0.081837,  0.2245,  -0.016932, -1,      0,      0,      0.11246, 0,  0,  0,  0,  0,       0;
       0.9947,  0.069636,  0.10291,  0.17456,   0,     -1,      0,      0.13486, 0,  0,  0,  0,  0,       0;
       1.7651,  0.19133,  -0.64716, -0.024406,  0,      0,     -1,      0.1621,  0,  0,  0,  0,  0,       0;
      -26.5447, 6.7367,   18.4808,  -0.15056,   0,      0,      0,      0,       1,  0,  0,  0,  0,       0;
       0,       0,         0,        0,          0.18,  0.77,   0.05,  -1,       0,  0,  0,  0,  0,       0;
       0,       0,         0,        0,          0.18,  0.77,   0.05,   0,      -1,  0,  0,  0,  0,       0;
       0,       0,         0,        0,          0,      0,      0,      0,      0, -1,  0,  0,  0,       0;
       0,       0,         0,        0,          0,      0,      0,      0,      0,  0, -1,  0,  0,       0;
       0,       0,         0,        0,          0,      0,      0,      0,      0,  0,  0, -1,  0,       0;
       0,      -1,         0,        0,          0,      0,      0,      0,      0,  0,  0,  0,  c(31),   0;
       0,       0,        -1,        0,          0,      0,      0,      0,      0,  0,  0,  0,  0,       c(31)]';

% --- ODE and measurement function handles --------------------------------
odeFunc  = @(x, u, xi, theta) ADM1_R3_core_ode_sym_pi(x, u, xi, theta, c, a);
measFunc = @(x, theta)        ADM1_R3_core_mgl_sym_pi(x, theta, c);

% --- Tunable parameter vector theta (9 x 1) ------------------------------
% Nominal values from Weinrich (2017), except Delta_S_ion and phi_IN.
theta0 = [parameters_r3{16,:};  % k_ch        [1/d]
          parameters_r3{21,:};  % k_pr        [1/d]
          parameters_r3{18,:};  % k_li        [1/d]
          parameters_r3{17,:};  % k_dec       [1/d]
          parameters_r3{19,:};  % k_m_ac      [1/d]
          parameters_r3{8,:};   % K_S_ac      [g/L]
          parameters_r3{4,:};   % K_I_nh3     [g/L]
          0.15 - 0.02;          % Delta_S_ion [g/L]
          1];                   % phi_IN      [-]

thetaLB = [2e-2;  1e-2;  1e-2;  1e-3;  1e-1;  1e-3;  1e-2;  -1e-2;  1e-1];
thetaUB = [1e1;   1e1;   1;     1;     1;     2;     1e1;   1;      2  ];

% --- Metadata struct -----------------------------------------------------
p.nParameters = numel(theta0);
p.names       = {'k_{ch}','k_{pr}','k_{li}','k_{dec}', ...
                  'k_{m,ac}','K_{S,ac}','K_{I,nh3}','\Delta S_{ion}','\phi_{IN}'};
p.units       = {'1/d','1/d','1/d','1/d','1/d','g/L','g/L','mol/L','-'};
p.nOutputs    = 6;
p.outputNames = {'q_{gas}','p_{CH4}','p_{CO2}','pH','S_{IN}','S_{ac}'};
p.outputUnits = {'m^3/d','bar','bar','-','g/L','g/L'};

end
