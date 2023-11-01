% Create a new figure with a specified name
close all
f = figure('Name', 'Temperature Plot');
allowable_temperature = 400;
known_mat = 460;
% Hold the plot to allow multiple data series to be plotted on the same axes
hold on;
set(gcf,'color','#eeeeee');
temps = readmatrix(pwd + "/Regen_workspace/heatflux_temps.xlsx");
T_7 = temps(1, 2:5001);
T_8 = temps(2, 2:5001);
T_9 = temps(3, 2:5001);
T_10 = temps(4, 2:5001);
T_11 = temps(5,2:5001);
T_12 = temps(6,2:5001);
T_13 = temps(7, 2:5001); 

% Set the font style for the plot
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Temperature plot with T_l removed for a more compact plot

p1 = plot(x_plot * 1000.*0.03937008, (T_7-273.15)*(9/5) +32, 'LineWidth', .8, 'DisplayName', '70%');
p2 = plot(x_plot * 1000.*0.03937008, (T_8-273.15)*(9/5) +32,'LineWidth', .8, 'DisplayName', '80%');
p3 = plot(x_plot * 1000.*0.03937008, (T_9-273.15)*(9/5) +32,  'LineWidth', .8, 'DisplayName', '90%');
% Draw a horizontal dashed line at the allowable temperature threshold
yline(allowable_temperature, '--r', 'Design temperature threshold (400°F)', 'LineWidth', 1);
yline(known_mat, '--r', 'Known material property threshold (460°F)', 'LineWidth', 1)
p4 = plot(x_plot * 1000.*0.03937008, (T_10-273.15)*(9/5) +32,'LineWidth', .8, 'DisplayName', '100%');
p5 = plot(x_plot * 1000.*0.03937008, (T_11-273.15)*(9/5) +32, 'LineWidth', .8, 'DisplayName', '110%');
p6 = plot(x_plot * 1000.*0.03937008, (T_12-273.15)*(9/5) +32, 'LineWidth', .8, 'DisplayName', '120%');
p7 = plot(x_plot * 1000.*0.03937008, (T_13-273.15)*(9/5) +32, 'LineWidth', .8, 'DisplayName', '130%');



% Set y-axis label and color
ylabel('Temperature [F]');
set(gca, 'Ycolor', 'k');
ylim([180 500])
%yticks([350 375 400 425 450 allowable_temperature 500])

% yyaxis right
% p3 = plot(x_plot .* 1000.*0.03937008, r_interpolated .* 1000.*0.03937008, 'black', 'LineStyle', '-',"DisplayName","chamber contour");
% ylabel('Radius [inches]')
% set(gca, 'Ycolor', 'k')
% ylim([.5 3.5])

% Enable grid lines for better visualization
grid on;

% Add labels and legend
xlabel('Location [inches]');
title('Temperature Variation: % Nominal Heat flux ');
legend([p7 p6 p5 p4 p3 p2 p1  ]); % Show legend with specified display names

% Customize legend position (optional)
legend('Location', 'Best');

% Adjust plot margins for better aesthetics (optional)
a = gca


% You can save the plot as an image file if needed (optional)
% saveas(gcf, 'temperature_plot.png');