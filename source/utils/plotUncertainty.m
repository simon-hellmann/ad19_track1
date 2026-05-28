%% plotUncertainty.m
% Horizontal bar chart of 1-sigma relative parameter uncertainty [%].
% Useful for assessing practical identifiability after each PI round.
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — produces a figure)
%
%% Input
%
%   thetaHat:   parameter estimates                                 (n_params x 1)
%   stdTheta:   parameter standard deviations from Cramer-Rao      (n_params x 1)
%   p:          meta struct with p.nParameters and p.names
%   titleStr:   figure title string
%

function plotUncertainty(thetaHat, stdTheta, p, titleStr)

relUnc = stdTheta(:) ./ abs(thetaHat(:)) * 100;  % [%]

figure('Name',titleStr, 'NumberTitle','off');
barh(1:p.nParameters, relUnc);
set(gca, 'YTick',1:p.nParameters, 'YTickLabel',p.names, 'YDir','reverse');
xlabel('1\sigma relative uncertainty [%]');
title(titleStr, 'Interpreter','none');
grid on; box on;
xlim([0, max(relUnc) * 1.15]);

end % fun
