%% plotCostChannels.m
% Bar plot of the WLS cost split by output channel.
%
% Author: Simon Hellmann. Created: 2026/05/23. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — produces a figure)
%
%% Input
%
%   J_ch:       (1 x n_out) WLS cost per output channel
%   p:          metadata struct. Fields used:
%               .nOutputs     number of output channels
%               .outputNames  LaTeX output names  (1 x nO cell)
%   title_str:  figure title string
%

function plotCostChannels(J_ch, p, title_str)

figure('Name',title_str, 'NumberTitle','off');

bar(J_ch);
set(gca, 'TickLabelInterpreter','latex');
xticks(1:p.nOutputs);
xticklabels(p.outputNames);
ylabel('WLS cost $J_k$ [-]', 'Interpreter','latex');
title(title_str, 'Interpreter','none');
grid on; box on;

end % fun
