function plotVarianceDecomposition(pi_decomp, p, titleStr)
% Stacked horizontal bar chart of the variance decomposition matrix.
% pi_decomp(i,k) is the proportion of parameter i's variance attributable
% to singular direction k.  Each row sums to 1.
%
% A parameter dominated by a single high-index (small) singular direction
% is poorly identifiable.  After PSS, the retained parameters should be
% predominantly explained by the leading singular directions.
%
% Inputs:
%   pi_decomp -- variance decomposition (n_params x n_params)
%   p         -- meta struct with p.nParameters and p.names
%   titleStr  -- figure title string

    n = p.nParameters;

    figure('Name', titleStr, 'NumberTitle', 'off');
    barh(1:n, pi_decomp, 'stacked');
    set(gca, 'YTick', 1:n, 'YTickLabel', p.names, 'YDir','reverse'); % flip order of y axis

    xlabel('Cumulative variance proportion');
    ylabel('Parameter');
    title(titleStr, 'Interpreter', 'none');
    xlim([0, 1]);
    grid on; box on;

    leg_labels = arrayfun(@(k) sprintf('SV %d', k), 1:n, 'UniformOutput', false);
    legend(leg_labels, 'Location', 'eastoutside');
end
