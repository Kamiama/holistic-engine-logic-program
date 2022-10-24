%% Engine Sizing Code CEA Wrapper
% Author: Kamon Blong (kblong@purdue.edu)
% First Created: 7/17/2022
% Last Updated: 

function [c_star, isp, exp_ratio, M, gamma, P, T, rho, mu, Pr, Mw, k, son, cp] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, sub, sup, file_name, throat, convert_to_imperial)

%{ 
Description: Sends engine performance and fuel properties to NASA CEA where
    a bunch of voodoo stoichiometry magic is performed to return relevant
    values for combustion properties of propellants.

Inputs:
- 

Outputs: 
- 
%}

% Change this variable to true to rerun CEA instead of using saved values
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;

%% CEA Inputs
% run setup
inp('type') = 'eq';              % sets the type of CEA calculation
inp('file_name') = file_name;    % input/output file name

% general parameters
inp('p') = P_c;                  % chamber pressure
inp('p_unit') = 'psi';           % chamber pressure units
if sub(1) && sup(1) ~= 0
    inp('sub') = sub;            % subsonic area ratios
    inp('sup') = sup;            % supersonic area ratios
elseif throat
    inp('sup') = 1;              % supersonic area ratios
elseif sub(1)
    inp('sub') = sub;            % subsonic area ratios
elseif sup(1)
    inp('sup') = sup;            % supersonic area ratios
else
    inp('pip') = P_c / P_e;      % pressure ratios  
end

% propellant inputs
inp('fuel') = fuel; % fuel name from thermo.inp
if fuel_weight(1) ~= 0
    inp('fuel_wt%') = fuel_weight;
end
inp('ox') = oxidizer; % ox name from thermo.inp
if fuel_temp(1) 
    inp('fuel_t') = fuel_temp;
end
if oxidizer_temp(1)
    inp('ox_t') = oxidizer_temp;
end
inp('o/f') = OF; % oxidizer / fuel ratio

%% Run CEA
if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% The output data structure, called 'data' in this case, is also a MATLAB
% map. 'data' contains a single entry for each of the CEA calculation types
% listed ('eq' and 'fr'). For instance, if only 'fr' is listed, then 'data'
% will only contain a single entry under data('fr').
data_eq = data('eq');

%% Extract Outputs
% Use keys(data_eq) or keys(data_fr) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.
isp = squeeze(data_eq('isp'));
c_star = squeeze(data_eq('cstar'));
exp_ratio = data_eq('ae/at');
M = squeeze(data_eq('mach'));
gamma = squeeze(data_eq('gammas'));
P = squeeze(data_eq('p'));
T = squeeze(data_eq('t'));
rho = squeeze(data_eq('rho'));
mu = squeeze(data_eq('visc'));
Pr = squeeze(data_eq('prandtl'));
Mw = squeeze(data_eq('m'));
k = squeeze(data_eq('k'));
son = squeeze(data_eq('son'));
cp = squeeze(data_eq('cp'));

oxidizer_temp = ('ox_t');
fuel_temp = ('fuel_t');

%% Convert & Correct Outputs
if convert_to_imperial
    %exp_ratio = 1 / M(end) * ((2 + (gamma(1) - 1 ) * M(end) ^ 2) / (gamma(1) + 1)) ^ ((gamma(1) + 1) / (2 * (gamma(1) - 1))); % calculate expansion ratio manually
    P = P.* 1.45038E-4; % convert [Pa] to [psi]
    isp = isp(end) / 9.8067; % normalize isp to seconds
    exp_ratio = exp_ratio(end); % select expansion ratio
    c_star = c_star(1) * 3.281; % convert [m/s] to [ft/s]
    T = T .* 1.8; % convert [K] to [R]
    rho = rho * 3.613E-5; % convert [kg/m^3] to [lbm/in^3]
    mu = mu * 1.450E-4; % convert [Pa-s] to [psi-s]
    Mw = Mw * 2.205; % convert [kg/kmol] to [lbm/kmol]
    k = k * 0.5782; % convert [W/m-K] to [Btu/hr-ft-R]
else
    isp = isp(end);
    c_star = c_star(1);
end