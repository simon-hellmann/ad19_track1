%% plotResidualBoxplot.m
% Boxplot of scaled residuals split by output channel.
%
% Author: Simon Hellmann. Created: 2026/05/23. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — produces a figure)
%
%% Input
%
%   r_scaled_cell:  (1 x n_out) cell — scaled residuals per output channel
%   p:              metadata struct. Fields used:
%                   .nOutputs     number of output channels
%                   .outputNames  LaTeX output names  (1 x nO cell)
%   title_str:      figure title string
%

function plotResidualBoxplot(r_scaled_cell, p, title_str)

% build concatenated residual vector and integer group index
all_r = [];
grp   = [];
for i = 1:p.nOutputs
    n_i   = numel(r_scaled_cell{i});
    all_r = [all_r; r_scaled_cell{i}(:)];
    grp   = [grp;   i * ones(n_i, 1)];
end

figure('Name',title_str, 'NumberTitle','off');
boxplot(all_r, grp, 'Labels',p.outputNames);
set(findobj(gca, 'type','text'), 'Interpreter','latex');
yline(0, 'k--');
ylabel('Scaled residual $\tilde{r}_k$ [-]', 'Interpreter','latex');
title(title_str, 'Interpreter','none');
grid on; box on;

end % fun
