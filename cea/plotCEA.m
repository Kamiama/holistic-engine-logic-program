%% NASA CEA Relation Plotter
% Author: Kamon Blong (kblong@purdue.edu)
% First Created: 10/23/2022
% Last Updated: 

%{ 
Description: Plots requested CEA outputs for given outputs - results are in
    displayed in metric
%}

file_name = 'plotcea';

%% Assign Inputs
P_c = 300;            % chamber pressure [psi]
P_e = 14.7;           % exit pressure [psi]
fuel = {'Jet-A(L)'};  % fuel formula for NASA CEA
fuel_weight = 0;      % fuel weights, does nothing now
oxidizer = {'O2(L)'}; % oxidizer formula for NASA CEA

iterating_value = "OF";
min_value = 1.5;
max_value = 3.5;
step_value = .05;

%% Iterate Over CEA Values
counter = 1;
for OF = min_value : step_value : max_value
    [c_star, isp, ~, ~, ~, P, T, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, 0, oxidizer, 0, OF, 0, 0, file_name, 0, 0);
    OF_matrix(counter) = min_value + step_value * (counter - 1);
    T_matrix(counter) = T(1);
    isp_matrix(counter) = isp(1);*
    counter = counter + 1;
end

%% Plot Reults
figure(1)
hold("on")
subplot(1,2,1)
plot(OF_matrix, T_matrix, 'black');
title('OF vs Temperature')
xlabel('OF')
ylabel('Temperature (K)')
grid on

subplot(1,2,2)
plot(OF_matrix, isp_matrix, 'blue');
title('OF vs isp')
xlabel('OF')
ylabel('isp (m/s)')
grid on