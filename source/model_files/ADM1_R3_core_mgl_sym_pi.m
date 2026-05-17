%% Version: Matlab R2022b, Update 6
% Author: Simon Hellmann
% first created: 2024-05-07

function g = ADM1_R3_core_mgl_sym_pi(x,th,c)
% delivers a symbolic expression (h) of measurement equation of the
% ADM1-R3-Core

% ion balance: 
Phi = th(8) + (x(4) - x(12))/17 - x(11)/44 - x(10)/60;
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 

% measurement equations
g = [c(13)*x(13)^2 + c(14)*x(13)*x(14) + c(15)*x(14)^2 + c(16)*x(13) + c(17)*x(14) + c(18); % volFlow
     c(19)*x(13);           % pch4
     c(20)*x(14);           % pco2
     -log10(SHPlus);        % pH
     x(4);                  % SIN
     x(1)];                 % Sac

end 