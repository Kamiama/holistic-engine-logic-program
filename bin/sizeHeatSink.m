% Taken From
% https://www.sciencedirect.com/science/article/pii/S1359431122013096\
% Heister, Huzel, Bartz
% Assumptions: heat transfer to ambient environment is negligible, physical
% properties of the material(density, sp heat, thermal conductivity) are
% not affected by temperature.

% Chamber Parameters
At = 0.60132; % throat area (in2)
Dt = (sqrt(At/pi)); %throat diameter 
Pc = 300; % chamber pressure (psi)
cstar = 139.012;
R = 1; % throat radius of curvature

Tst = 3608.46; % stagnation temp (R)
mu = 9.3587*10^(-9); % viscosity
g = 32.17; % gravity, ft/s2
cpS = 0.12; % sp. heat capacity of 316 SS (BTU/lb-F) (0.5 J/g-C) 
cp = 0.55; %sp. heat capacity of gas (BTU/lb F)
M = 1; % mach number
Twg = Tst; % hot-gas-side local chamber-wall temp
Taw = Tst*0.9; % adiabatic wall temp of gas, = Tst*recovery factor (from 0.9-0.98)
gamma = 1.3;

A = At; % local cross sectional area
Pr = (4*gamma)/(9*gamma - 5); % Prandtl number, approx from Bartz
sigma = (((Twg/(2*Tst))*(1+((gamma-1)/2* M^2))+0.5)^-0.68) * (1+(gamma-1)/2 *M^2)^-0.12;
hg = ((0.026/(Dt^0.2)) * ((mu^0.2*cp)/Pr^0.6) * (Pc*g / cstar)^0.8 * (Dt/R)^0.1)*(At/A)^0.9 * sigma
q = hg*(Taw-Twg)


%{ 
% "Simplified 1D heat tranfer model" paper formulas - need to rewrite
T0 = 1;
delta = 1; % wall thickness [m]
lamda = 1; % thermal conduct coefficent of wall
dens = 1; %density of wall
mu1 = 1; % eigenvalue
dT = lamda/(dens*cp)*Twg; % 2.2 - needs to be reviewed/redone
Bi = hg*delta/lamda;
A1 = 2*(sin(mu1)/(mu1+sin(mu1)*cos(mu1)));
B1 = sin(mu1)/mu1;
q = (1-A1*(exp(-mu1^2*F0))*B1*dens*cp*delta*(T0-Tst)); 
%}