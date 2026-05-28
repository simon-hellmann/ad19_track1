%% mess2pssData.m
% Convert a cut MESS struct (output of cutMESS) into the PSS data format
% expected by main_pss_workflow.m.
%
% Output channels mapped to ADM1-R3-Core outputs:
%   i=1  gasflow [m^3/d]  i=2  p_CH4 [bar]  i=3  p_CO2 [bar]
%   i=4  pH [-]           i=5  S_IN [g/L]    i=6  S_ac [g/L]
%
% Author: Simon Hellmann. Created: 2026/05/17. Version: Matlab R2022b, Update 6
%
%% Output
%
%   data:           struct with fields:
%                   .tMeas{i}       relative measurement times per output  [d]
%                   .yMeas{i}       measured values per output
%                   .t_feed_start   feeding event start times              [d]
%                   .feed_mass      total substrate mass per event         [kg]
%                   .xi_feed        inlet concentration matrix (n_events x n_xi)
%                   .t0             0 (relative time origin)
%                   .tf             total window duration                  [d]
%                   .T_start        absolute start datetime
%                   .MESS           original MESS_cut for reference
%
%% Input
%
%   MESS_cut:       cut MESS struct from cutMESS (has .rel_time fields)
%

function data = mess2pssData(MESS_cut)

%% measurement cell arrays

data.tMeas = cell(6, 1);
data.yMeas = cell(6, 1);

data.tMeas{1} = MESS_cut.online_gas.rel_time_gasflow(:);
data.yMeas{1} = MESS_cut.online_gas.gasflow(:);

data.tMeas{2} = MESS_cut.online_gas.rel_time_gasflow(:);
data.yMeas{2} = MESS_cut.online_gas.methane(:);

data.tMeas{3} = MESS_cut.online_gas.rel_time_gasflow(:);
data.yMeas{3} = MESS_cut.online_gas.co2(:);

data.tMeas{4} = MESS_cut.online_liq.rel_time(:);
data.yMeas{4} = MESS_cut.online_liq.pH(:);

data.tMeas{5} = MESS_cut.offline_nh3.rel_time(:);
data.yMeas{5} = MESS_cut.offline_nh3.sin(:);

data.tMeas{6} = MESS_cut.offline_ac_eq.rel_time(:);
data.yMeas{6} = MESS_cut.offline_ac_eq.ac_eq(:);

%% feeding events

data.t_feed_start = MESS_cut.feed.rel_time(:);  % [d]
data.feed_mass    = MESS_cut.feed.mass(:);       % [kg]
data.xi_feed      = MESS_cut.feed.charac;        % (n_events x n_xi)

%% timing and metadata

data.t0      = 0;
data.tf      = MESS_cut.t_span;   % [d]
data.T_start = MESS_cut.T_start;
data.MESS    = MESS_cut;

end % fun
