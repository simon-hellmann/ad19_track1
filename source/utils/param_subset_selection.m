%% subsetSelection.m
% Performs parameter subset selection based on scaled output sensitivity,
% singular value decomposiiton, QRP decomposition and epsilon threshold.
% Based on the measurement data provided, a subset of the total parameter
% vector is practically identifiable, and its complement is not. 
%
%Diana C. López C., Tilman Barz, Stefan Körkel, Günter Wozny,
%Nonlinear ill-posed problem analysis in model-based parameter estimation and experimental design,
%Computers & Chemical Engineering,
%Volume 77,
%2015,
%Pages 24-42,
%ISSN 0098-1354,
%https://doi.org/10.1016/j.compchemeng.2015.03.002.
%
%% Inputs
%
%   dy_do:      sensitivities (scaled), size (n_data*n_outputs)x(n_params)
%   kappa_max:  maximum condition number (scaler)
%   gamma_max:  maximum colinearity index (scaler)
%   p:          parameter structure, struct (variable)
%
%% Outputs
%
%   si:         singular values, size (n_params)x(n_params), descending
%   keep_idx:   indices of singular values to keep, vector (variable)
%   C:          parameter covariance, (n_params)x(n_params)
%   pi:         variance decomposition proportion, (n_params x 1)
%   epsilon:    epsilon-threshold, scalar

function [si, keep_idx, C_pp, pi, epsilon] = subsetSelection(dy_do, ...
    kappa_max, gamma_max, p)
    
    % singular value decomposition: 
    [U,S,V] = svd(dy_do,'econ');
    si = diag(S);

    % compute epsilon threshold and cut parameter vector: 
    kappa = si(1)/si(end);
    gamma = 1/si(end);
    epsilon = max(si(1)/kappa_max,1/gamma_max); % SH: eq. (13)
    si_keep = si((si>epsilon)); % all SV above the threshold
    
    % re-order parameter vector: 
    [Q,R,P] = qr(dy_do,0);
    keep_idx = P((si>epsilon)); 
    
    % initialize covariance matrices: 
    C_pp = zeros(p.nParameters); 
    var_xi_theta = zeros(p.nParameters);
    
    for col_k = 1:p.nParameters
        C_pp = C_pp + (V(:,col_k)*V(:,col_k).')/si(col_k)^2; % Cramer-Rao bound: C_pp = F^-1 = V*S^-2*V', eq. (17)
        for row_k = 1:p.nParameters
            var_xi_theta(row_k,col_k) = V(row_k,col_k)^2/si(col_k)^2; % eq. (18)
        end
    end
    
    var_theta = diag(C_pp);
    pi = var_xi_theta./var_theta; % SH: total covariance proportion per row
   
end