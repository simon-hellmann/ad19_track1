%% plotFeed.m
% Plot feed rate as a staircase into a given axes handle.
% Works as a reusable subplot helper or as a standalone figure.
%
% Author: Simon Hellmann. Created: 2026/05/23. Version: Matlab R2022b, Update 6
%
%% Output
%
%   (none — draws into ax or opens a new figure)
%
%% Input
%
%   ax:             axes handle to draw into; pass [] to open a new figure
%   t_events:       (n_ev x 1) feed-event time grid    [d]
%   u_segments:     (n_ev-1 x 1) feed rate per segment [m^3/d]
%

function plotFeed(ax, t_events, u_segments)

if isempty(ax)
    figure;
    ax = gca;
end

u_stairs = [u_segments; u_segments(end)];

axes(ax);
hold on; box on; grid on;
stairs(t_events, u_stairs, 'k-', 'LineWidth',1.5);
ylabel('u_{feed} [m^3/d]');
xlim([t_events(1), t_events(end)]);

end % fun
