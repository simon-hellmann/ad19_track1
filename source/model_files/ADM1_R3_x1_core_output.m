%% ADM1_R3_x1_core_output.m
% Output equation for ADM1-R3-x1. Identical to ADM1-R3-Core; K_S_IN does
% not appear in the output equations so no changes are required.
%
% Author: Simon Hellmann. Created: 2026/05/25. Version: Matlab R2022b, Update 6
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
%   th:         (10x1) tunable parameter vector (see setup_ADM1_R3_x1)
%   c:          (31x1) time-invariant parameters (see setup_ADM1_R3_x1)
%

function g = ADM1_R3_x1_core_output(x, th, c)

% ion balance:
Phi    = th(8) + (x(4) - x(12))/17 - x(11)/44 - x(10)/60;
% equivalent proton concentration:
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4));

% measurement equations (c indices unchanged from ADM1-R3-Core)
g = [c(13)*x(13)^2 + c(14)*x(13)*x(14) + c(15)*x(14)^2 + c(16)*x(13) + c(17)*x(14) + c(18); % q_gas
     c(19)*x(13);     % p_CH4
     c(20)*x(14);     % p_CO2
     -log10(SHPlus);  % pH
     x(4);            % S_IN
     x(1)];           % S_ac

end % fun
