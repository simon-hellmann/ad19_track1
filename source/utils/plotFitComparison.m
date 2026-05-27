%% plotFitComparison.m
% Overlay PI #1 and PI #2 smooth simulations against measurements,
% with feed rate on top. Shows that fit quality is preserved after PSS.
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — produces a figure)
%
%% Input
%
%   t_meas_long:    (N_long x 1) measurement times, long vector        [d]
%   y_meas_long:    (N_long x 1) measured values, long vector
%   t_fine:         (n_fine x 1) equidistant fine time grid            [d]
%   y_fine1:        (n_fine x n_out) PI #1 simulated outputs on t_fine
%   y_fine2:        (n_fine x n_out) PI #2 simulated outputs on t_fine
%   out_idx:        (N_long x 1) output channel index per entry        [1..n_out]
%   p:              metadata struct. Fields used:
%                   .nOutputs     number of output channels
%                   .outputNames  LaTeX output names  (1 x nO cell)
%                   .outputUnits  output unit strings (1 x nO cell)
%   t_events:       (n_ev x 1) feed-event time grid                    [d]
%   u_segments:     (n_ev-1 x 1) feed rate per segment                 [m^3/d]
%   titleStr:       figure title string
%

function plotFitComparison(t_meas_long, y_meas_long, t_fine, y_fine1, y_fine2, ...
    out_idx, p, t_events, u_segments, titleStr)

n_rows = p.nOutputs + 1;

figure('Name',titleStr, 'NumberTitle','off');
plotFeed(subplot(n_rows, 1, 1), t_events, u_segments);

for i = 1:p.nOutputs
    mask = (out_idx == i);
    if ~any(mask); continue; end

    subplot(n_rows, 1, i + 1);
    hold on; box on; grid on;

    plot(t_meas_long(mask), y_meas_long(mask), 'ko', ...
         'MarkerSize',4, 'DisplayName','Measured');
    plot(t_fine, y_fine1(:,i), 'b-', ...
         'LineWidth',1.5, 'DisplayName','PI #1 (full set)');
    plot(t_fine, y_fine2(:,i), 'r--', ...
         'LineWidth',1.5, 'DisplayName','PI #2 (PSS subset)');

    ylabel([p.outputNames{i}, ' [', p.outputUnits{i}, ']']);
    xlim([t_events(1), t_events(end)]);
    if i == 1
        legend('Location','best');
    end
end

xlabel('Time [d]');
sgtitle(titleStr, 'Interpreter','none');

end % fun
