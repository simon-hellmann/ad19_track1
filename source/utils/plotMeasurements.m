%% plotMeasurements.m
% Plot raw measurement data for each output channel, with feed rate on top.
% One subplot per output channel plus one subplot for the feed rate.
%
% Author: Simon Hellmann. Created: 2026/05/23. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — produces a figure)
%
%% Input
%
%   t_meas_long:    (N_long x 1) measurement times, long vector        [d]
%   y_meas_long:    (N_long x 1) measured values, long vector
%   out_idx:        (N_long x 1) output channel index per entry        [1..n_out]
%   p:              metadata struct. Fields used:
%                   .nOutputs     number of output channels
%                   .outputNames  LaTeX output names  (1 x nO cell)
%                   .outputUnits  output unit strings (1 x nO cell)
%   t_events:       (n_ev x 1) feed-event time grid                    [d]
%   u_segments:     (n_ev-1 x 1) feed rate per segment                 [m^3/d]
%   title_str:      figure title string
%

function plotMeasurements(t_meas_long, y_meas_long, out_idx, p, ...
    t_events, u_segments, title_str)

n_rows = p.nOutputs + 1;

figure('Name',title_str, 'NumberTitle','off');
plotFeed(subplot(n_rows, 1, 1), t_events, u_segments);

for i = 1:p.nOutputs
    mask = out_idx == i;
    if ~any(mask); continue; end

    subplot(n_rows, 1, i + 1);
    hold on; box on; grid on;

    plot(t_meas_long(mask), y_meas_long(mask), 'ko', ...
         'MarkerSize',4, 'MarkerFaceColor','k');

    ylabel([p.outputNames{i}, ' [', p.outputUnits{i}, ']']);
    xlim([t_events(1), t_events(end)]);
end

xlabel('Time [d]');
sgtitle(title_str, 'Interpreter','none');

end % fun
