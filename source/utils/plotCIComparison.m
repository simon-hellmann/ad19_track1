function plotCIComparison(thetaHat1, stdTheta1, thetaHat2, stdTheta2_sub, ...
                          keep_idx, p, titleStr)
% Grouped horizontal bar chart comparing 2-sigma relative uncertainty
% before (PI #1, full set) and after (PI #2, PSS subset) for the
% identifiable parameters only.
%
% Expected result: much smaller bars for PI #2, confirming that removing
% ill-posed directions tightens confidence intervals substantially.
%
% Inputs:
%   thetaHat1      -- PI #1 estimates, full vector (n_params x 1)
%   stdTheta1      -- PI #1 std devs, full vector (n_params x 1)
%   thetaHat2      -- PI #2 estimates, full vector (n_params x 1)
%   stdTheta2_sub  -- PI #2 std devs, subset only (numel(keep_idx) x 1)
%   keep_idx       -- indices of identifiable parameters
%   p              -- meta struct with p.names
%   titleStr       -- figure title string

    relUnc1 = stdTheta1(keep_idx) ./ abs(thetaHat1(keep_idx)) * 100;
    relUnc2 = stdTheta2_sub(:)    ./ abs(thetaHat2(keep_idx)) * 100;
    n_keep  = numel(keep_idx);

    figure('Name', titleStr, 'NumberTitle', 'off');
    bh = barh(1:n_keep, [relUnc1(:), relUnc2(:)]);
    bh(1).FaceColor  = [0.2, 0.4, 0.8];
    bh(1).DisplayName = 'PI #1 (full set)';
    bh(2).FaceColor  = [0.8, 0.2, 0.2];
    bh(2).DisplayName = 'PI #2 (PSS subset)';

    set(gca, 'YTick', 1:n_keep, 'YTickLabel', p.names(keep_idx));
    xlabel('rel. std. deviation [%]');
    legend('Location', 'best');
    title(titleStr, 'Interpreter', 'none');
    grid on; box on;
    xlim([0, max([relUnc1; relUnc2]) * 1.15]);
end
