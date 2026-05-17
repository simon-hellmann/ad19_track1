%% Version: Matlab R2022b, Update 6
% Author: Simon Hellmann
% first created: 2024-05-07

function [f, I_ac, I_ph,I_Nlim,I_nh3] = ADM1_R3_core_ode_sym_pi(x, u, xi, th, c, a)
% return symbolic expression (f) of the right-hand side of the ODE system 
% of the ADM1-R3-Core (no water, no ash, no 2. CH-fraction, no ions as 
% additional state) with constant values for u and xin 
% also return total inhibition function and its subordinate factors

% ion balance: 
Phi = th(8) + (x(4) - x(12))/17 - x(11)/44 - x(10)/60;
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 

% compute inhibition factors:
I_ph = c(3)/(c(3) + SHPlus^(c(2))); 
I_Nlim = x(4)/(x(4) + c(8)); 
I_nh3 = th(7)/(th(7) + x(12));
I_ac = I_ph * I_Nlim * I_nh3; % overall inhibition factor

% dynamic equations
f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(5) + a(1,2)*th(2)*x(6) + a(1,3)*th(3)*x(7) + a(1,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac; 
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(5) + a(2,2)*th(2)*x(6) + a(2,3)*th(3)*x(7) - c(5)*x(2) + c(6)*x(13) + a(2,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;
     c(1)*(xi(3) - x(3))*u + a(3,1)*th(1)*x(5) + a(3,2)*th(2)*x(6) + a(3,3)*th(3)*x(7) - c(5)*x(3) + c(5)*x(11) + c(7)*x(14) + a(3,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;
     c(1)*(xi(4)*th(9) - x(4))*u + a(4,1)*th(1)*x(5) + a(4,2)*th(2)*x(6) + a(4,3)*th(3)*x(7) + a(4,4)*th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac;    
     c(1)*(xi(5) - x(5))*u - th(1)*x(5) + a(5,5)*th(4)*x(8) + a(5,6)*th(4)*x(9);    
     c(1)*(xi(6) - x(6))*u - th(2)*x(6) + a(6,5)*th(4)*x(8) + a(6,6)*th(4)*x(9);
     c(1)*(xi(7) - x(7))*u - th(3)*x(7) + a(7,5)*th(4)*x(8) + a(7,6)*th(4)*x(9);
     c(1)*(xi(8) - x(8))*u + a(8,1)*th(1)*x(5) + a(8,2)*th(2)*x(6) + a(8,3)*th(3)*x(7) - th(4)*x(8);    
     c(1)*(xi(9) - x(9))*u + th(5)*x(1)*x(9)/(th(6) + x(1))*I_ac - th(4)*x(9);    
     c(28)*(x(1) - x(10)) - c(9)*x(10)*SHPlus;
     c(29)*(x(3) - x(11)) - c(10)*x(11)*SHPlus
     c(30)*(x(4) - x(12)) - c(11)*x(12)*SHPlus;
     c(21)*x(13)^3 + c(22)*x(13)^2*x(14) + c(23)*x(13)*x(14)^2 + c(24)*x(13)^2 + c(25)*x(13)*x(14) + c(12)*x(2) + c(26)*x(13);    
     c(23)*x(14)^3 + c(22)*x(13)*x(14)^2 + c(21)*x(13)^2*x(14) + c(25)*x(14)^2 + c(24)*x(13)*x(14) + c(12)*x(3) - c(12)*x(11) + c(27)*x(14)];    
end 