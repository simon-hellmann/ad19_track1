%% plotVarianceDecomposition.m
% Stacked horizontal bar chart of the variance decomposition matrix.
% pi_decomp(i,k) is the proportion of parameter i's variance attributable
% to singular direction k.  Each row sums to 1.
%
% A parameter dominated by a single high-index (small) singular direction
% is poorly identifiable.
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — produces a figure)
%
%% Input
%
%   pi_decomp:  variance decomposition matrix   (n_params x n_params)
%   p:          meta struct with p.nParameters and p.names
%   titleStr:   figure title string
%

function plotVarianceDecomposition(pi_decomp, p, titleStr)

n = p.nParameters;

figure('Name',titleStr, 'NumberTitle','off');
barh(1:n, pi_decomp, 'stacked');
set(gca, 'YTick',1:n, 'YTickLabel',p.names, 'YDir','reverse');

xlabel('Cumulative variance proportion');
ylabel('Parameter');
title(titleStr, 'Interpreter','none');
xlim([0, 1]);
grid on; box on;

leg_labels = arrayfun(@(k) sprintf('SV %d', k), 1:n, 'UniformOutput',false);
legend(leg_labels, 'Location','eastoutside');

end % fun
