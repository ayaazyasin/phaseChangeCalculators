% Single iteration of Wasyner-TST loop
% Ayaaz Yasin - Nov 08, 2022
clear; %clc; 

% Fluid Properties
RV      = 1.55;                      % [kg/m^3], vapor density
RL      = 70.115;                    % [kg/m^3], liquid density
TV      = 21;                        % [K], vapor temp
TL      = 21.1;                      % [K], liquid temp
PV      = 121405;                    % [Pa], vapor pressure

% Constants
M       = 2.016e-3; 		         % [kg/mol], molar mass of H2
R       = 8.31446261815324;          % [J/K-mol], gas constant
hfg     = (451.98-6.4659)*1e3;       % [J/kg], latent heat of vaporization
VL      = 2.8751e-05;                % [m3/mol], molar volume of H2 at 21 K
sigma   = 0.0018029;                 % [N/m], surface tension gradient at 21 K
K       = 2;                         % [1/m], average curvature

% TST
[a,l]   = tst_alpha(RV,RL);          % calling TST-alpha function

% Wayner
alpha_coeff     = 2*a/(2-a);                                % [-]
S_coeff         = sqrt(M/(2*pi*R*TL));                      % [s/m]
Pc              = sigma*K;                                  % [Pa], capillary pressure
tempTerm        = (PV*M*hfg*(TL-TV))/(R*TV*TL);             % [Pa]
presTerm        = (VL*PV*Pc)/(R*TL);                        % [Pa]
mAreaFlux       = alpha_coeff*S_coeff*(tempTerm-presTerm);  % [kg/m2-s], mass flux

% Integration
r         = 5e-3;                   % [m], test cell radius
area      = pi*r^2;                 % [m^2], surface area assuming planar interface
mFlow     = mAreaFlux*area;         % [kg/s], total evaporative mass flow
mFlow_exp = 55e-9;                  % [kg/s], experimentally measured evaporative mass flow
err_ratio = mFlow/mFlow_exp;

% print statements
fprintf('<strong>Wayner calculation</strong>\n')
fprintf('alpha\t\t\t= %0.4f\t\t[-]\n',a)
fprintf('hfg\t\t\t\t= %0.0f\t\tJ/kg\n\n',hfg)
fprintf('alpha_coeff\t\t= %0.4f\t\t[-]\n',alpha_coeff)
fprintf('sqrt_coeff\t\t= %0.4f\t\ts/m\n',S_coeff)
fprintf('TL-TV\t\t\t= %0.4f\t\tK\n',TL-TV)
fprintf('temp_term\t\t= %0.4e\tPa\n',tempTerm)
fprintf('pres_term\t\t= %0.4e\tPa\n\n',presTerm)
fprintf('massAreaFlux\t= %0.4e\tkg/m^2.s\t<--\n',mAreaFlux)
fprintf('massFlow \t\t= %0.4e\tkg/s\n',mFlow)
fprintf('Error Ratio\t\t= %0.4e\tkg/s\n',err_ratio)
fprintf('---------------\n\n')