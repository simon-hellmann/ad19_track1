%% param_subset_selection.m
% Performs parameter subset selection based on scaled output sensitivity,
% singular value decomposition, QRP decomposition and epsilon threshold.
% Based on the measurement data provided, a subset of the total parameter
% vector is practically identifiable, and its complement is not.
%
%Diana C. Lopez C., Tilman Barz, Stefan Korkel, Gunter Wozny,
%Nonlinear ill-posed problem analysis in model-based parameter estimation and experimental design,
%Computers & Chemical Engineering,
%Volume 77,
%2015,
%Pages 24-42,
%ISSN 0098-1354,
%https://doi.org/10.1016/j.compchemeng.2015.03.002.
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Input
%
%   dy_dp_scaled:   double-scaled sensitivities (rows / scale_long, cols * thetaHat)
%                   size (n_data*n_outputs) x (n_params)
%   kappa_max:      maximum condition number (scalar)
%   gamma_max:      maximum collinearity index (scalar)
%   p:              parameter metadata struct with p.nParameters
%   thetaHat:       parameter estimates in physical units (n_params x 1);
%                   used to re-transform C_pp from dimensionless to physical units
%
%% Output
%
%   si:         singular values                                  (n_params x 1), descending
%   keep_idx:   indices of identifiable parameters, vector       (variable length)
%   C_pp:       parameter covariance in physical units           (n_params x n_params)
%   pi:         variance decomposition proportions               (n_params x n_params)
%               rows sum to 1; invariant to thetaHat scaling
%   epsilon:    identifiability threshold                        scalar
%   kappa:      condition number of sensitivity matrix
%   gamma:      collinearity index of sensitivity matrix
%

function [si, keep_idx, C_pp, pi, epsilon, kappa, gamma] = ...
    param_subset_selection(dy_dp_scaled, kappa_max, gamma_max, p, thetaHat)

% singular value decomposition:
[U,S,V] = svd(dy_dp_scaled,'econ');  %#ok<ASGLU>
si = diag(S);

% compute epsilon threshold and identify kept singular values:
kappa = si(1)/si(end);
gamma = 1/si(end);
epsilon_kappa = si(1)/kappa_max;
epsilon_gamma = 1/gamma_max;
epsilon  = max(epsilon_kappa, epsilon_gamma);   % eq. (13)
si_keep  = si(si > epsilon);                    %#ok<NASGU>

% re-order parameter vector via QRP:
[Q,R,P] = qr(dy_dp_scaled, 0);                 %#ok<ASGLU>
keep_idx = P(si > epsilon);

%% Postprocessing

% covariance and variance decomposition in scaled (dimensionless) space:
C_pp_scaled  = zeros(p.nParameters);
var_xi_theta = zeros(p.nParameters);

for col_k = 1:p.nParameters
    C_pp_scaled = C_pp_scaled + ...
        (V(:,col_k) * V(:,col_k).') / si(col_k)^2;   % eq. (17)
    for row_k = 1:p.nParameters
        var_xi_theta(row_k,col_k) = V(row_k,col_k)^2 / si(col_k)^2;  % eq. (18)
    end % for
end % for

% pi proportions: compute before re-transformation (invariant to scaling
% because numerator and denominator both scale by thetaHat(k)^2):
var_theta_scaled = diag(C_pp_scaled);
pi = var_xi_theta ./ var_theta_scaled;

% re-transform C_pp from dimensionless to physical parameter units,
% reversing the column scaling (.*thetaHat') applied before PSS.
D_theta = diag(thetaHat(:));
C_pp    = D_theta * C_pp_scaled * D_theta;

end % fun
