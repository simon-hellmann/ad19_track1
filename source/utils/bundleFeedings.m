%% bundleFeedings.m
% Bundle consecutive feeding events separated by less than feeding_duration
% to prevent overlapping boluses in the ODE simulation.
%
% A new group starts whenever the gap between consecutive event start times
% reaches feeding_duration.  All events within a group are collapsed into
% one pulse: earliest start time, summed mass, mass-weighted average xi.
%
% Author: Simon Hellmann. Created: 2026/05/27. Version: Matlab R2022b, Update 6
%
%% Output
%
%   data_out:   PSS struct — identical to data_in with t_feed_start,
%               feed_mass, and xi_feed replaced by bundled arrays
%
%% Input
%
%   data_in:    PSS struct with fields
%               .t_feed_start   -- (n_ev x 1) feed event start times [d]
%               .feed_mass      -- (n_ev x 1) total mass per event   [kg]
%               .xi_feed        -- (n_ev x n_xi) inlet compositions
%               (all other fields are passed through unchanged)
%
%   feeding_duration:   bolus duration [d]; events with a gap smaller than
%                       this value are bundled (gap >= duration → no overlap)

function data_out = bundleFeedings(data_in, feeding_duration)

data_out = data_in;

n_ev = numel(data_in.t_feed_start);
if n_ev <= 1
    return
end

% Assign group IDs — new group when gap >= feeding_duration
group_id      = ones(n_ev, 1);
current_group = 1;
for i_ev = 2:n_ev
    if data_in.t_feed_start(i_ev) - data_in.t_feed_start(i_ev-1) >= feeding_duration
        current_group = current_group + 1;
    end
    group_id(i_ev) = current_group;
end % for
n_groups = current_group;

n_xi             = size(data_in.xi_feed, 2);
t_feed_start_out = nan(n_groups, 1);
feed_mass_out    = nan(n_groups, 1);
xi_feed_out      = nan(n_groups, n_xi);

for i_g = 1:n_groups
    mask_g       = group_id == i_g;
    mass_g       = data_in.feed_mass(mask_g);   % [kg]
    total_mass_g = sum(mass_g);

    t_feed_start_out(i_g) = min(data_in.t_feed_start(mask_g));
    feed_mass_out(i_g)    = total_mass_g;
    xi_feed_out(i_g, :)   = (mass_g' * data_in.xi_feed(mask_g, :)) / total_mass_g;
end % for

if n_groups < n_ev
    fprintf("  bundleFeedings: %d events → %d bundles (threshold %.1f min)\n", ...
        n_ev, n_groups, feeding_duration * 24*60);
end

data_out.t_feed_start = t_feed_start_out;
data_out.feed_mass    = feed_mass_out;
data_out.xi_feed      = xi_feed_out;

end % fun
