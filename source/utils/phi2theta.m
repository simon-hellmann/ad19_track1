%% phi2theta.m
% Convert a phi (log10-reparametrised) vector back to physical theta-space.
% Applies 10.^ only to the log-scaled entries; linear entries pass through unchanged.
% Works for the full parameter vector (log_mask = log_idx as integer indices) and
% for a subset vector (log_mask = keep_log_mask as a logical mask).
%
% Author: Simon Hellmann. Created: 2026/05/25. Version: Matlab R2022b, Update 6
%
%% Output
%
%   theta:      (n x 1) parameter vector in physical units
%
%% Input
%
%   phi:        (n x 1) parameter vector in phi-space
%   log_mask:   index vector or logical mask selecting the log10-scaled entries

function theta = phi2theta(phi, log_mask)

theta           = phi;
theta(log_mask) = 10.^(phi(log_mask));

end % fun
