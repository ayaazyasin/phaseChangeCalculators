% Single iteration of schrage-TST loop
% Ayaaz Yasin - Oct 31, 2022
clear; clc; 

% Fluid Properties
RV      = 1.55;                      % [kg/m^3], vapor density
RL      = 70.115;                    % [kg/m^3], liquid density
TV      = 21;                        % [K], vapor temp
TL      = 21.1;                      % [K], liquid temp
PV      = 121405;                    % [Pa], vapor pressure

% Vapor saturation temperature at vapor density
aa = 26.8923;
bb = 0.0100;
cc = -12.4620;
dd = -0.4054;
TVsat = aa*exp(bb*RV) + cc*exp(dd*RV);

% Liquid Saturation Pressure at current Liquid Temperature
p4 = -3.821302413486544e+05;
p3 =  7.673143119143602e+04;
p2 = -5.485344443202910e+03;
p1 =  1.415843632072226e+02;
PL =  p1*TL^3 + p2*TL^2 + p3*TL + p4;
% PL=PL/1000; % <----- uncomment to scale liquid pressure

% Constants
mBar = 3.347447494734E-27; 		% [kg] molecular mass of H2
kb = 1.380649E-23; 				% [m2 kg s-2 K-1] Boltzmann's constant

% TST
[a,l] = tst_alpha(RV,RL);

% Schrage
liq_term    = PL/sqrt(TL);
vap_term    = PV/sqrt(TV);
A           = 2*a/(2-a);
S           = sqrt(mBar/(2*pi*kb));
mAreaFlux   = A*S*(liq_term-vap_term);

% Integration
r         = 5e-3;                   % [m], test cell radius
area      = pi*r^2;                 % [m^2], surface area assuming planar interface
mFlow     = mAreaFlux*area;         % [kg/s], total evaporative mass flow
mFlow_exp = 55e-9;                  % [kg/s], experimentally measured evaporative mass flow
err_ratio = mFlow/mFlow_exp;

% print statements
fprintf('<strong>Schrage calculation</strong>\n')
fprintf('PL_sat\t\t\t= %0.0f\t\tPa\n',PL)
fprintf('alpha \t\t\t= %0.4f\t\t[-]\n',a)
fprintf('massAreaFlux \t= %0.4e\tkg/m^2.s\t<--\n',mAreaFlux)
fprintf('massFlow \t\t= %0.4e\tkg/s\n',mFlow)
fprintf('Error Ratio\t\t= %0.4e\tkg/s\n',err_ratio)
fprintf('---------------\n\n')
