% Single iteration of Beta-form -TST loop
% Ayaaz Yasin - Nov 09, 2023
clear; %clc; 

% Fluid Properties
RV      = 1.55;                      % [kg/m^3], vapor density
RL      = 70.115;                    % [kg/m^3], liquid density
TV      = 21;                        % [K], vapor temp
TL      = 21.1;                      % [K], liquid temp
PV      = 121405;                    % [Pa], vapor pressure

% Constants
R       = 8.31446261815324;          % [J/K-mol], gas constant
hfg     = (451.98-6.4659)*1e3;       % [J/kg], latent heat of vaporization
sigma   = 0.0018029;                 % [N/m], surface tension gradient at 21 K
K       = 2;                         % [1/m], average curvature
kB      = 1.380649E-23; 		     % [kg.m2/s2.K], Boltzmann's constant
m       = 3.347447494734E-27; 		 % [kg] molecular mass of H2
Pc      = sigma*K;                   % [Pa], capillary pressure

% saturation vapor pressure
p1 = 144.5710;
p2 = -5.7056e+03;
p3 = 8.2028e+04;
p4 = -4.2358e+05;   
PVsat = p1.*TV.^3 + p2.*TV.^2 + p3.*TV + p4;

% TST
[a,l]   = tst_alpha(RV,RL);          % calling TST-alpha function

% beta-form mass flux equation
cR = sqrt(2*kB*TV/m);                                       % [m/s]
alpha_coeff     = 2*a/(2-a);                                % [-]
S_coeff         = PV/(sqrt(pi)*cR);                         % [s/m]
beta            = 1;                                        % [-]
tempRatio       = sqrt(TV/TL);                              % [-]

    term1  = PVsat/PV;
    term2a = 1 - (TV/TL);
    term2b = RV*hfg/PV;
    term3a = TV/TL;
    term3b = RV/RL;
    term3c = Pc/PV; % no disjoining pressure
presRatio       = term1 + term2a*term2b + term3a*term3b*term3c; % [-]
mAreaFlux       = alpha_coeff*S_coeff*(beta*presRatio*tempRatio-1);  % [kg/m2-s], mass flux

% Integration
r         = 5e-3;                   % [m], test cell radius
area      = pi*r^2;                 % [m^2], surface area assuming planar interface
mFlow     = mAreaFlux*area;         % [kg/s], total evaporative mass flow
mFlow_exp = 55e-9;                  % [kg/s], experimentally measured evaporative mass flow
err_ratio = mFlow/mFlow_exp;

% print statements
fprintf('<strong>Beta-form calculation</strong>\n')
fprintf('alpha\t\t\t= %0.4f\t\t[-]\n',a)
fprintf('hfg\t\t\t\t= %0.0f\t\tJ/kg\n',hfg)
fprintf('cR\t\t\t\t= %0.4e\tm/s\n\n',cR)
fprintf('alpha_coeff\t\t= %0.4f\t\t[-]\n',alpha_coeff)
fprintf('s_coeff\t\t\t= %0.4f\t\tkg/m^2.s\n',S_coeff)
fprintf('TL-TV\t\t\t= %0.4f\t\tK\n',TL-TV)
fprintf('beta\t\t\t= %0.4f\t\t[-]\n',beta)
fprintf('temp_ratio\t\t= %0.4e\tPa\n',tempRatio)
fprintf('pres_ratio\t\t= %0.4e\tPa\n',presRatio)
fprintf('\tPR_term1\t= %0.4e\t[-]\n',term1)
fprintf('\tPR_term2a\t= %0.4e\t[-]\n',term2a)
fprintf('\tPR_term2b\t= %0.4e\t[-]\n',term2b)
fprintf('\tPR_term3a\t= %0.4e\t[-]\n',term3a)
fprintf('\tPR_term3b\t= %0.4e\t[-]\n',term3b)
fprintf('\tPR_term3c\t= %0.4e\t[-]\n\n',term3c)
fprintf('massAreaFlux\t= %0.4e\tkg/m^2.s\t<--\n',mAreaFlux)
fprintf('massFlow \t\t= %0.4e\tkg/s\n',mFlow)
fprintf('Error Ratio\t\t= %0.4e\tkg/s\n',err_ratio)
fprintf('---------------\n\n')