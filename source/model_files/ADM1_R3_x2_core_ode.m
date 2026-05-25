%% ADM1_R3_x2_core_ode.m
% ODE right-hand side for ADM1-R3-x2. Identical structure to ADM1-R3-Core,
% but K_S_IN (th(10)), k_La (th(11)), and k_p (th(12)) are tunable.
% k_La- and k_p-composite quantities are recomputed from theta at each call.
%
% Author: Simon Hellmann. Created: 2026/05/24. Version: Matlab R2022b, Update 6
%
%% Output
%
%   f:          (14x1) state derivative dx/dt
%   I_ac:       combined acetoclastic inhibition factor          [-]
%   I_ph:       pH inhibition sub-factor                         [-]
%   I_Nlim:     nitrogen limitation sub-factor                   [-]
%   I_nh3:      free ammonia inhibition sub-factor               [-]
%
%% Input
%
%   x:          (14x1) state vector                              [various]
%   u:          volumetric feed rate                             [m^3/d]
%   xi:         (n_xi x 1) inlet composition vector              [various]
%   th:         (12x1) tunable parameter vector (see setup_ADM1_R3_x2)
%   c:          (19x1) time-invariant parameters (see setup_ADM1_R3_x2)
%   a:          (14x11) Petersen stoichiometry  (see setup_ADM1_R3_x2)
%

function [f, I_ac, I_ph, I_Nlim, I_nh3] = ADM1_R3_x2_core_ode(x, u, xi, th, c, a)

% ion balance:
Phi    = th(8) + (x(4) - x(12))/17 - x(11)/44 - x(10)/60;
% equivalent proton concentration:
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4));

% inhibition factors (K_S_IN is now th(10)):
I_ph   = c(3)/(c(3) + SHPlus^(c(2)));
I_Nlim = x(4)/(x(4) + th(10));
I_nh3  = th(7)/(th(7) + x(12));
I_ac   = I_ph * I_Nlim * I_nh3;

% Compute k_La and k_p composites from tunable th(11) and th(12):
% kp_lin_coeff and kp_ph2o_contrib are general for any p_h2o (c(18));
% when c(18)=0 they reduce to kp_per_Vgas_patm*p_atm and 0 respectively.
kLa_Vliq_per_Vgas  = th(11) * c(10);                              % k_La * V_liq/V_gas
kp_per_Vgas_patm   = th(12) * c(13);                              % k_p / (V_gas * p_atm)
kp_lin_coeff       = kp_per_Vgas_patm * (c(19) - 2*c(18));        % coeff of linear-in-conc gas outflow terms
kp_ph2o_contrib    = kp_per_Vgas_patm * c(18) * (c(19) - c(18));  % water-vapour offset; zero for dry biogas

% dynamic equations
f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(5) + a(1,2)*th(2)*x(6) + a(1,3)*th(3)*x(7) + a(1,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(5) + a(2,2)*th(2)*x(6) + a(2,3)*th(3)*x(7) - th(11)*x(2) + th(11)*c(5)*x(13) + a(2,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;
     c(1)*(xi(3) - x(3))*u + a(3,1)*th(1)*x(5) + a(3,2)*th(2)*x(6) + a(3,3)*th(3)*x(7) - th(11)*x(3) + th(11)*x(11) + th(11)*c(6)*x(14) + a(3,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;
     c(1)*(xi(4)*th(9) - x(4))*u + a(4,1)*th(1)*x(5) + a(4,2)*th(2)*x(6) + a(4,3)*th(3)*x(7) + a(4,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;
     c(1)*(xi(5) - x(5))*u - th(1)*x(5) + a(5,5)*th(4)*x(8) + a(5,6)*th(4)*x(9);
     c(1)*(xi(6) - x(6))*u - th(2)*x(6) + a(6,5)*th(4)*x(8) + a(6,6)*th(4)*x(9);
     c(1)*(xi(7) - x(7))*u - th(3)*x(7) + a(7,5)*th(4)*x(8) + a(7,6)*th(4)*x(9);
     c(1)*(xi(8) - x(8))*u + a(8,1)*th(1)*x(5) + a(8,2)*th(2)*x(6) + a(8,3)*th(3)*x(7) - th(4)*x(8);
     c(1)*(xi(9) - x(9))*u + th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac - th(4)*x(9);
     c(14)*(x(1) - x(10)) - c(7)*x(10)*SHPlus;
     c(15)*(x(3) - x(11)) - c(8)*x(11)*SHPlus;
     c(16)*(x(4) - x(12)) - c(9)*x(12)*SHPlus;
     -kp_per_Vgas_patm*c(11)^2*x(13)^3 - 2*kp_per_Vgas_patm*c(11)*c(12)*x(13)^2*x(14) - kp_per_Vgas_patm*c(12)^2*x(13)*x(14)^2 + kp_lin_coeff*c(11)*x(13)^2 + kp_lin_coeff*c(12)*x(13)*x(14) + kLa_Vliq_per_Vgas*x(2) + (kp_ph2o_contrib - kLa_Vliq_per_Vgas*c(5))*x(13);
     -kp_per_Vgas_patm*c(12)^2*x(14)^3 - 2*kp_per_Vgas_patm*c(11)*c(12)*x(13)*x(14)^2 - kp_per_Vgas_patm*c(11)^2*x(13)^2*x(14) + kp_lin_coeff*c(12)*x(14)^2 + kp_lin_coeff*c(11)*x(13)*x(14) + kLa_Vliq_per_Vgas*x(3) - kLa_Vliq_per_Vgas*x(11) + (kp_ph2o_contrib - kLa_Vliq_per_Vgas*c(6))*x(14)];

end % fun
