function plotUncertainty(thetaHat, stdTheta, p, titleStr)
% Horizontal bar chart of 2-sigma relative parameter uncertainty [%].
% Useful for assessing practical identifiability after each PI round.
%
% Inputs:
%   thetaHat  -- parameter estimates (n_params x 1)
%   stdTheta  -- parameter standard deviations from Cramer-Rao (n_params x 1)
%   p         -- meta struct with p.nParameters and p.names
%   titleStr  -- figure title string

    relUnc = 2 * stdTheta(:) ./ abs(thetaHat(:)) * 100;   % [%]

    figure('Name', titleStr, 'NumberTitle', 'off');
    barh(1:p.nParameters, relUnc);
    set(gca, 'YTick', 1:p.nParameters, 'YTickLabel', p.names);
    xlabel('2\sigma relative uncertainty [%]');
    title(titleStr, 'Interpreter', 'none');
    grid on; box on;
    xlim([0, max(relUnc) * 1.15]);
end
