%% validate_ADM1_R3_x1.m
% Cross-checks ADM1-R3-x1 against reference ADM1-R3 by evaluating both ODE
% and output functions at identical (x0, u, xi) conditions. Since x1 merely
% promotes K_S_IN, k_La, k_p from c into theta while leaving all equations
% structurally unchanged, both models must return numerically identical f and
% g vectors when theta_x1 carries the same nominal values.
%
% Author: Simon Hellmann. Created: 2026/05/24. Version: Matlab R2022b, Update 6
%
% no inputs and outputs (main file)
%

clc; clear; close all

[script_path, ~, ~] = fileparts(mfilename('fullpath'));
cd(script_path);

%% Setup and paths

addpath("../source/model_files");
addpath("../source/model_files/model_data");

tol_pass = 1e-8;   % absolute tolerance for numerical equivalence

%% Load parameters and assemble both models

load("ADM1_parameters.mat", 'parameters');
parameters_r3 = parameters.ADM1_R3;

V_liq = 0.012;   % [m^3] liquid volume
V_gas = 0.003;   % [m^3] gas headspace volume

[c_r3, ~, odeFunc_r3, measFunc_r3, theta0_r3, ~, ~, p_r3] = ...
    setup_ADM1_R3(parameters_r3, V_liq, V_gas);

[c_x1, ~, odeFunc_x1, measFunc_x1, theta0_x1, ~, ~, ~] = ...
    setup_ADM1_R3_x1(parameters_r3, V_liq, V_gas);

%% Define test conditions

% rough but physically representative initial state (from main_pss_workflow.m)
x0 = [0.049;   % S_ac      [g/L]
      0.012;   % S_ch4     [g/L]
      4.975;   % S_IC      [g/L]
      0.964;   % S_IN      [g/L]
      2.962;   % X_ch      [g/L]
      0.949;   % X_pr      [g/L]
      0.412;   % X_li      [g/L]
      1.926;   % X_bac     [g/L]
      0.522;   % X_ac      [g/L]
      0.049;   % S_ac-     [g/L]
      4.546;   % S_hco3-   [g/L]
      0.022;   % S_nh3     [g/L]
      0.358;   % S_ch4_gas [g/L]
      0.660];  % S_co2_gas [g/L]

u_test = 100/1000/1000;   % [g/d] -> [m^3/d] representative feed rate

xi_test = [0.002;   % S_ac  inlet [g/L]
           0;       % S_ch4 inlet [g/L]
           3.0;     % S_IC  inlet [g/L]
           1.5;     % S_IN  inlet [g/L]
           8.0;     % X_ch  inlet [g/L]
           4.0;     % X_pr  inlet [g/L]
           2.0;     % X_li  inlet [g/L]
           0;       % X_bac inlet [g/L]
           0];      % X_ac  inlet [g/L]

%% Evaluate ODE and compare f vectors

[f_r3, I_ac_r3, I_ph_r3, I_Nlim_r3, I_nh3_r3] = odeFunc_r3(x0, u_test, xi_test, theta0_r3);
[f_x1, I_ac_x1, I_ph_x1, I_Nlim_x1, I_nh3_x1] = odeFunc_x1(x0, u_test, xi_test, theta0_x1);

diff_f = f_x1 - f_r3;

state_names = {'S_ac','S_ch4','S_IC','S_IN','X_ch','X_pr','X_li','X_bac','X_ac', ...
               'S_ac-','S_hco3-','S_nh3','S_ch4_g','S_co2_g'};

fprintf("\n=== ODE: f_x1 vs f_r3 ===\n");
fprintf("%-10s  %13s  %13s  %13s\n", "State", "f_r3", "f_x1", "|diff|");
for k = 1:14
    fprintf("%-10s  %13.6e  %13.6e  %13.6e\n", ...
        state_names{k}, f_r3(k), f_x1(k), abs(diff_f(k)));
end

norm_diff_f = norm(abs(diff_f));
if norm_diff_f < tol_pass
    fprintf("ODE check:  PASS  (norm|diff| = %.2e < %.2e)\n", norm_diff_f, tol_pass);
else
    fprintf("ODE check:  FAIL  (norm|diff| = %.2e >= %.2e)\n", norm_diff_f, tol_pass);
end

%% Evaluate output and compare g vectors

g_r3 = measFunc_r3(x0, theta0_r3);
g_x1 = measFunc_x1(x0, theta0_x1);

diff_g = g_x1 - g_r3;

output_names = {'q_gas','p_CH4','p_CO2','pH','S_IN','S_ac'};

fprintf("\n=== Output: g_x1 vs g_r3 ===\n");
fprintf("%-8s  %13s  %13s  %13s\n", "Output", "g_r3", "g_x1", "|diff|");
for k = 1:p_r3.nOutputs
    fprintf("%-8s  %13.6e  %13.6e  %13.6e\n", ...
        output_names{k}, g_r3(k), g_x1(k), abs(diff_g(k)));
end

max_diff_g = max(abs(diff_g));
if max_diff_g < tol_pass
    fprintf("Output check:  PASS  (max|diff| = %.2e < %.2e)\n", max_diff_g, tol_pass);
else
    fprintf("Output check:  FAIL  (max|diff| = %.2e >= %.2e)\n", max_diff_g, tol_pass);
end

%% Inhibition factor comparison

fprintf("\n=== Inhibition factors ===\n");
fprintf("%-8s  %10s  %10s  %10s\n", "Factor", "R3", "x1", "|diff|");
fprintf("%-8s  %10.6f  %10.6f  %10.2e\n", "I_ph",   I_ph_r3,   I_ph_x1,   abs(I_ph_x1   - I_ph_r3));
fprintf("%-8s  %10.6f  %10.6f  %10.2e\n", "I_Nlim", I_Nlim_r3, I_Nlim_x1, abs(I_Nlim_x1 - I_Nlim_r3));
fprintf("%-8s  %10.6f  %10.6f  %10.2e\n", "I_nh3",  I_nh3_r3,  I_nh3_x1,  abs(I_nh3_x1  - I_nh3_r3));
fprintf("%-8s  %10.6f  %10.6f  %10.2e\n", "I_ac",   I_ac_r3,   I_ac_x1,   abs(I_ac_x1   - I_ac_r3));
