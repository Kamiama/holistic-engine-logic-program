%% NASA CEA Relation Plotter
% Author: Kamon Blong (kblong@purdue.edu)
% First Created: 4/10/2023
% Last Updated: 

%{ 
Description: Plots requested CEA outputs for given outputs - results are in
    displayed in metric
%}

clear;
clc;
close all;

CEA_input_name = 'plotcea';

%% Initialize Variables

% propellants
fuels = [{'CH3OH'}, {'C2H5OH'}, {'C3H8O,2propanol'}, {'Jet-A(L)'}];  % list of fuel strings for NASA CEA
oxidizer = {'O2(L)'}; % oxidizer formula for NASA CEA
fuel_weight = 0;      % fuel weights
fuel_temp = [273.15, 273.15, 273.15, 273.15];   % fuel temperature [K]
oxidizer_temp = 90.19;   % fuel temperature [K]

% performance parameters
P_c = 275;            % chamber pressure [psi]
P_e = 14.7;           % exit pressure [psi]

min_value_OF = .5;
max_value_OF = 3;
step_value_OF = .05;

OF_matrix = ones(1, round((max_value_OF - min_value_OF) / step_value_OF));
isp_matrix = ones(1, round((max_value_OF - min_value_OF) / step_value_OF));

% miscellaneous
legend_str = {'Data 1', 'Data 2'};

%% Iterate Over CEA Values & Plot
i = 1;
j = 1;

figure('Name', 'Fuel Trade Study')
hold("on")
legend(legend_str, 'location', 'Northwest');

% fuel / OF loop
for fuel = fuels
    OF = min_value_OF;
    [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp(j), oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 0, CEA_input_name);
    for OF = min_value_OF : step_value_OF : max_value_OF
        [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp(j), oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 0, CEA_input_name);
        OF_matrix(i) = min_value_OF + step_value_OF * (i - 1);
        isp_matrix(i) = isp(1) / 9.81;
        i = i + 1;
    end
    i = 1;
    plot(OF_matrix, isp_matrix, 'Linewidth', 4);
    propellant_string = num2cell(char(fuel + " / " + oxidizer), 2);
    legend_str{j} = [strjoin(propellant_string), j];
    legend(legend_str);

    [~, max_index] = max(isp_matrix);
    plot(OF_matrix(max_index), max(isp_matrix), 'k.', 'MarkerSize', 30, 'HandleVisibility', 'off')
  
    j = j + 1;
end

xlabel('Mixture Ratio')
ylabel('Ideal Isp (sec)')
set(gca, 'XLim', [0.75, 2.4], 'FontSize', 17)
set(gca, 'YLim', [205, 265])
grid on

%title("Propellant Performance vs Mixture Ratio: " + P_c + " psi P_c, " + P_e + " psi P_e")