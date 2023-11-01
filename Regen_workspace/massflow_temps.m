% Create a new figure with a specified name
close all
f = figure('Name', 'Temperature Plot');
allowable_temperature = 400;
known_mat = 460;
% Hold the plot to allow multiple data series to be plotted on the same axes
hold on;
set(gcf,'color','#eeeeee');
temps = readmatrix(pwd + "/Regen_workspace/mass_flow_data.xlsx");
T_7 = temps(1, 1:5000);
T_6 = temps(2, 1:5000);
T_5 = temps(3, 1:5000);
T_4 = temps(4, 1:5000);
T_3 = temps(5,1:5000);

% Set the font style for the plot
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Temperature plot with T_l removed for a more compact plot

p1 = plot(x_plot * 1000.*0.03937008, (T_7-273.15)*(9/5) +32, 'LineWidth', .8, 'DisplayName', '7 lb/s');
p2 = plot(x_plot * 1000.*0.03937008, (T_6-273.15)*(9/5) +32,'LineWidth', .8, 'DisplayName', '6 lb/s');
p3 = plot(x_plot * 1000.*0.03937008, (T_5-273.15)*(9/5) +32,  'LineWidth', .8, 'DisplayName', '5 lb/s');
% Draw a horizontal dashed line at the allowable temperature threshold
yline(allowable_temperature, '--r', 'Design temperature threshold (400°F)', 'LineWidth', 1);
yline(known_mat, '--r', 'Known material property threshold (460°F)', 'LineWidth', 1)
p4 = plot(x_plot * 1000.*0.03937008, (T_4-273.15)*(9/5) +32,'LineWidth', .8, 'DisplayName', '4 lb/s');
p5 = plot(x_plot * 1000.*0.03937008, (T_3-273.15)*(9/5) +32, 'LineWidth', .8, 'DisplayName', '3 lb/s');

% Set y-axis label and color
ylabel('Temperature [F]');
set(gca, 'Ycolor', 'k');
ylim([220 520])
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
title('Temperature Variation: Coolant Mass Flow Rate');
legend([p1 p2 p3 p4 p5 ]); % Show legend with specified display names

% Customize legend position (optional)
legend('Location', 'Best');

% Adjust plot margins for better aesthetics (optional)
a = gca


% You can save the plot as an image file if needed (optional)
% saveas(gcf, 'temperature_plot.png');