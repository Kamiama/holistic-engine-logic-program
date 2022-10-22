% CEA_ROCKET_EXAMPLE: Example file for the MATLAB CEA wrapper. For in-depth
% documentation read the headers of cea_rocket_run.m,
% cea_rocket_run_single.m, and cea_rocket_read.m

% Change this variable to true to rerun CEA instead of using saved values
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;
inp('type') = 'eq fr';              % Sets the type of CEA calculation
inp('p') = 1000;                    % Chamber pressure
inp('p_unit') = 'psi';              % Chamber pressure units
inp('pi/p') = 200;          % pressure ratio
inp('o/f') = 1.3;               % Mixture ratio
inp('sup') = 5:1:25;                     % Supersonic area ratios
% inp('pip') = [5];                 % Pressure ratios
inp('fuel') = 'NH3 (L)';            % Fuel name from thermo.inp
inp('fuel_t') = 239.72;                % Fuel inlet temperature
inp('ox') = 'O2 (L)';               % Ox name from thermo.inp
inp('ox_t') = 298;                  % Ox inlet temperature
inp('file_name') = 'HW6_X15.inp';   % Input/output file name
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
data_fr = data('fr');

% Use keys(data_eq) or keys(data_fr) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.
temperature = squeeze(data_eq('t'));

% Plots chamber temperature (hence the 1 in the first column) vs. O/F which
% corresponds to H2O2 concentration for a variety of supersonic area
% ratios.
percent_o2 = inp('o/f') ./ (inp('o/f') + 1) * 100;
plot(percent_o2, squeeze(temperature(1, :, :)));
ylabel('Chamber temperature (K)');
xlabel('%O_2');
leg = legend([{'stag'}, arrayfun(@num2str, inp('sup'), 'Uniform', false)], ...
    'Location', 'Best');
title('O_2 decomposition temperature for varying %O_2 with NH_3');
[leg,att] = legend('show');
title(leg, 'Area ratio')
leg.Title.Visible = 'on';