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

file_name = 'plotcea';

%% Assign Inputs
P_c = 300;            % chamber pressure [psi]
P_e = 5;           % exit pressure [psi]
fuel = {'C3H8O,1propanol'};  % fuel formula for NASA CEA
fuel_temp = 293.15;   % fuel temperature [K]
fuel_weight = 0;      % fuel weights, does nothing now
oxidizer = {'O2(L)'}; % oxidizer formula for NASA CEA
OF = 1.5;

iterating_value = "P_c";
min_value_OF = 1.1;
max_value_OF = 1.9;
step_value_OF = .02;

min_value_P_c = 100;
max_value_P_c = 1000;
step_value_P_c = 5;

%% Iterate Over CEA Values
counter = 1;

if iterating_value == "OF"
    OF = min_value_OF;
    [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, 0, OF, 0, 0, file_name, 0, 0);
    for OF = min_value_OF : step_value_OF : max_value_OF
        [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, 0, OF, 0, 0, file_name, 0, 0);
        matrix(counter) = min_value_OF + step_value_OF * (counter - 1);
        T_matrix(counter) = T(1);
        isp_matrix(counter) = isp(1) / 9.81;
        counter = counter + 1;
    end

elseif iterating_value == "P_c"
    P_c = min_value_OF;
    [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, 0, OF, 0, 0, file_name, 0, 0);
    for P_c = min_value_P_c : step_value_P_c : max_value_P_c
        [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, 0, OF, 0, 0, file_name, 0, 0);
        matrix(counter) = min_value_P_c + step_value_P_c * (counter - 1);
        T_matrix(counter) = T(1);
        isp_matrix(counter) = isp(1) / 9.81;
        counter = counter + 1;
    end
end

%% Plot Reults
figure(1)

hold("on")
subplot(1,2,1)
plot(matrix, T_matrix, 'black');
if iterating_value == "OF"
    title('OF vs Temperature')
    xlabel('OF')
elseif iterating_value == "P_c"
    title('P_c vs Temperature')
    xlabel('P_c [psi]')
end
ylabel('Temperature [K]')
grid on

subplot(1,2,2)
plot(matrix, isp_matrix, 'blue');
if iterating_value == "OF"
    title('OF vs isp')
    xlabel('OF')
elseif iterating_value == "P_c"
    title('P_c vs isp')
    xlabel('P_c [psi]')
end
ylabel('isp (sec)')

sgtitle(iterating_value + ": " + fuel + " / " + oxidizer)
grid on