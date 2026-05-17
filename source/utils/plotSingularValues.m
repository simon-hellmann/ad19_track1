function plotSingularValues(si, epsilon, p, titleStr)
% Semilogy plot of the singular value spectrum with the epsilon threshold.
% Singular values above the threshold correspond to identifiable directions.
%
% Inputs:
%   si        -- singular values in descending order (n_params x 1)
%   epsilon   -- identifiability threshold scalar
%   p         -- meta struct with p.nParameters (used for x-axis ticks)
%   titleStr  -- figure title string

    n  = numel(si);
    idx = 1:n;

    figure('Name', titleStr, 'NumberTitle', 'off');
    semilogy(idx, si, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
             'DisplayName', 'Singular values');
    hold on;
    yline(epsilon, '--r', 'LineWidth', 1.5, ...
          'DisplayName', ['\epsilon = ', num2str(epsilon, '%.2e')]);

    n_keep = sum(si > epsilon);
    xline(n_keep + 0.5, ':k', 'LineWidth', 1, ...
          'DisplayName', sprintf('Cut (%d identifiable)', n_keep));

    xlabel('Singular value index (descending)');
    ylabel('Singular value magnitude');
    title(titleStr, 'Interpreter', 'none');
    legend('Location', 'northeast');
    grid on; box on;
    xticks(idx);
    xlim([0.5, n + 0.5]);
end
