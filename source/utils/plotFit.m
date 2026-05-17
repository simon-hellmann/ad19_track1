function plotFit(t_meas_long, y_meas_long, y_sim_long, out_idx, p, titleStr)
% Plot measured vs simulated for each output channel.
% One subplot per output; measurements as markers, simulation as line.
%
% Inputs:
%   t_meas_long  -- measurement times, long vector (N_long x 1)
%   y_meas_long  -- measured values, long vector (N_long x 1)
%   y_sim_long   -- simulated values, long vector (N_long x 1)
%   out_idx      -- output index for each entry (N_long x 1)
%   p            -- parameter/output meta struct with fields:
%                     p.nOutputs, p.outputNames, p.outputUnits
%   titleStr     -- figure title string

    figure('Name', titleStr, 'NumberTitle', 'off');

    for i = 1:p.nOutputs
        mask = (out_idx == i);
        if ~any(mask); continue; end

        subplot(p.nOutputs, 1, i);
        hold on; box on; grid on;

        plot(t_meas_long(mask), y_meas_long(mask), 'o', ...
             'MarkerSize', 4, 'DisplayName', 'Measured');
        plot(t_meas_long(mask), y_sim_long(mask), '-', ...
             'LineWidth', 1.5, 'DisplayName', 'Simulated');

        ylabel([p.outputNames{i}, ' [', p.outputUnits{i}, ']']);
        if i == 1
            legend('Location', 'best');
        end
    end

    xlabel('Time [d]');
    sgtitle(titleStr, 'Interpreter', 'none');
end
