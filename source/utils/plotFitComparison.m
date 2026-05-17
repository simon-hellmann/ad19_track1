function plotFitComparison(t_meas_long, y_meas_long, y_sim_long1, y_sim_long2, ...
                           out_idx, p, titleStr)
% Overlay PI #1 and PI #2 simulations against measurements for each output.
% Shows that fit quality is preserved after PSS.
%
% Inputs:
%   t_meas_long   -- measurement times, long vector (N_long x 1)
%   y_meas_long   -- measured values, long vector (N_long x 1)
%   y_sim_long1   -- PI #1 simulated values (N_long x 1)
%   y_sim_long2   -- PI #2 simulated values (N_long x 1)
%   out_idx       -- output index for each entry (N_long x 1)
%   p             -- meta struct with p.nOutputs, p.outputNames, p.outputUnits
%   titleStr      -- figure title string

    figure('Name', titleStr, 'NumberTitle', 'off');

    for i = 1:p.nOutputs
        mask = (out_idx == i);
        if ~any(mask); continue; end

        subplot(p.nOutputs, 1, i);
        hold on; box on; grid on;

        plot(t_meas_long(mask), y_meas_long(mask), 'ko', ...
             'MarkerSize', 4, 'DisplayName', 'Measured');
        plot(t_meas_long(mask), y_sim_long1(mask), 'b-', ...
             'LineWidth', 1.5, 'DisplayName', 'PI #1 (full set)');
        plot(t_meas_long(mask), y_sim_long2(mask), 'r--', ...
             'LineWidth', 1.5, 'DisplayName', 'PI #2 (PSS subset)');

        ylabel([p.outputNames{i}, ' [', p.outputUnits{i}, ']']);
        if i == 1
            legend('Location', 'best');
        end
    end

    xlabel('Time [d]');
    sgtitle(titleStr, 'Interpreter', 'none');
end
