%% ADM1_R3_core_output.m
% Output equation for ADM1-R3-Core. Returns the 6 measured outputs.
%
% Author: Simon Hellmann. Created: 2024/05/07. Version: Matlab R2022b, Update 6
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
%   th:         (9x1) tunable parameter vector  (see setup_ADM1_R3)
%   c:          (31x1) time-invariant parameters (see setup_ADM1_R3)
%

function g = ADM1_R3_core_output(x, th, c)

% ion balance:
Phi = th(8) + (x(4) - x(12))/17 - x(11)/44 - x(10)/60;
% equivalent proton concentration:
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4));

% measurement equations
g = [c(13)*x(13)^2 + c(14)*x(13)*x(14) + c(15)*x(14)^2 + c(16)*x(13) + c(17)*x(14) + c(18); % q_gas
     c(19)*x(13);     % p_CH4
     c(20)*x(14);     % p_CO2
     -log10(SHPlus);  % pH
     x(4);            % S_IN
     x(1)];           % S_ac

end % fun
