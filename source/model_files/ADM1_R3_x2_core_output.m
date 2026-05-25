%% ADM1_R3_x2_core_output.m
% Output equation for ADM1-R3-x2. Same six outputs as ADM1-R3-Core;
% k_p (now th(12)) enters the q_gas expression.
%
% Author: Simon Hellmann. Created: 2026/05/24. Version: Matlab R2022b, Update 6
%
%% Output
%
%   g:          (6x1) output vector:
%               g(1)  q_gas   volumetric gas flow rate  [m^3/d]
%               g(2)  p_CH4   partial pressure, methane [bar]
%               g(3)  p_CO2   partial pressure, CO2     [bar]
%               g(4)  pH      pH value                  [-]
%               g(5)  S_IN    inorganic nitrogen        [g/L]
%               g(6)  S_ac    acetate                   [g/L]
%
%% Input
%
%   x:          (14x1) state vector                              [various]
%   th:         (12x1) tunable parameter vector (see setup_ADM1_R3_x2)
%   c:          (19x1) time-invariant parameters (see setup_ADM1_R3_x2)
%

function g = ADM1_R3_x2_core_output(x, th, c)

% ion balance:
Phi    = th(8) + (x(4) - x(12))/17 - x(11)/44 - x(10)/60;
% equivalent proton concentration:
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4));

% total biogas pressure including water vapour (c(18) = p_h2o, default 0):
p_tot = c(11)*x(13) + c(12)*x(14) + c(18);

% measurement equations
% q_gas = k_p/p_atm * p_tot * (p_tot - p_atm) = th(12)*c(17)*p_tot^2 - th(12)*p_tot
g = [th(12)*c(17)*p_tot^2 - th(12)*p_tot;  % q_gas [m^3/d]
     c(11)*x(13);                            % p_CH4 [bar]
     c(12)*x(14);                            % p_CO2 [bar]
     -log10(SHPlus);                         % pH    [-]
     x(4);                                   % S_IN  [g/L]
     x(1)];                                  % S_ac  [g/L]

end % fun
