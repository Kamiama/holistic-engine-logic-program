%% NASA CEA Relation Plotter
% Author: Kamon Blong (kblong@purdue.edu)
% First Created: 10/23/2022
% Last Updated: 

%{ 
Description: Plots requested CEA outputs for given outputs - results are in
    displayed in metric
%}

clear;
clc;

CEA_input_name = 'plotcea';

%% Assign Inputs
P_c = 275;            % chamber pressure [psi]
P_e = 18.5;           % exit pressure [psi]
OF = 1.3;             % OF ratio [N/A]
fuel = {'C3H8O,2propanol'};  % fuel formula for NASA CEA
fuel_temp = [293.15];   % fuel temperature [K]
fuel_weight = 0;      % fuel weights, does nothing now
oxidizer = {'O2(L)'}; % oxidizer formula for NASA CEA
oxidizer_temp = 0;   % fuel temperature [K]

iterating_value = "P_c";
min_value_OF = .2;
max_value_OF = 2;
step_value_OF = .05;

min_value_P_c = 50;
max_value_P_c = 350;
step_value_P_c = .5;


%% Iterate Over CEA Values
counter = 1;

if iterating_value == "OF"
    OF = min_value_OF;
    [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 0, CEA_input_name);
    for OF = min_value_OF : step_value_OF : max_value_OF
        [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 0, CEA_input_name);
        matrix(counter) = min_value_OF + step_value_OF * (counter - 1);
        T_matrix(counter) = T(1);
        isp_matrix(counter) = isp(1) / 9.81;
        counter = counter + 1;
    end

elseif iterating_value == "P_c"
    P_c = min_value_P_c;
    [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 0, CEA_input_name);
    for P_c = min_value_P_c : step_value_P_c : max_value_P_c
        [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 0, CEA_input_name);
        matrix(counter) = min_value_P_c + step_value_P_c * (counter - 1);
        T_matrix(counter) = T(1);
        isp_matrix(counter) = isp(1) / 9.81;
        counter = counter + 1;
    end
end

%% Plot Reults
figure('Name', 'PlotCEA')
hold("on")

yyaxis left
plot(matrix, isp_matrix, 'blue');
if iterating_value == "OF"
    xlabel('OF')
elseif iterating_value == "P_c"
    xlabel('P_c [psi]')
end
ylabel('Ideal Isp (sec)')
set(gca, 'Ycolor', 'k')

yyaxis right
plot(matrix, T_matrix, 'red');
if iterating_value == "OF"
    xlabel('OF')
elseif iterating_value == "P_c"
    xlabel('P_c [psi]')
end
ylabel('Combustion Temperature [K]')
set(gca, 'Ycolor', 'k')

grid on
legend('Ideal Isp', 'Combustion Temperature', 'Location', 'Northwest')

if iterating_value == "OF"
    title(iterating_value + " plot:     " + fuel + " / " + oxidizer + " @ " + P_c + " psi P_c, " + P_e + " psi P_e")
elseif iterating_value == "P_c"
    title(iterating_value + " plot:     " + fuel + " / " + oxidizer + " @ " + OF + " OF ratio, " + P_e + " psi P_e")
end

grid on