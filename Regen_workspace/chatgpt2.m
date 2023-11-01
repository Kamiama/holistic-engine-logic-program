% Create a new figure with a specified name
close all
f = figure('Name', 'Temperature Plot');
allowable_temperature = 400;
% Hold the plot to allow multiple data series to be plotted on the same axes
hold on;
set(gcf,'color','#eeeeee');
temps = readmatrix(pwd + "/Regen_workspace/tolerance_temps.xlsx");
T_basecase = temps(1, 2:5001);
T_tol1 = temps(2, 2:5001);
T_tol2 = temps(3, 2:5001);
T_wlbasecase = temps(4, 2:5001);
T_wltol1 = temps(5,2:5001);
T_wltol2 = temps(6,2:5001);
% Set the font style for the plot
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Temperature plot with T_l removed for a more compact plot

p1 = plot(x_plot * 1000.*0.03937008, (T_basecase-273.15)*(9/5) +32, 'r-', 'LineWidth', .8, 'DisplayName', 'T_wg nominal');
p2 = plot(x_plot * 1000.*0.03937008, (T_tol1-273.15)*(9/5) +32, 'r--', 'LineWidth', 1, 'DisplayName', 'T_wg average');
p3 = plot(x_plot * 1000.*0.03937008, (T_tol2-273.15)*(9/5) +32, 'r:', 'LineWidth', 1, 'DisplayName', 'T_wg double average');
% Draw a horizontal dashed line at the allowable temperature threshold
yline(allowable_temperature, '--r', 'Design temperature threshold (400°F)', 'LineWidth', 1.5);
p4 = plot(x_plot * 1000.*0.03937008, (T_wlbasecase-273.15)*(9/5) +32, 'm-', 'LineWidth', .8, 'DisplayName', 'T_wl nominal');
p5 = plot(x_plot * 1000.*0.03937008, (T_wltol1-273.15)*(9/5) +32, 'm--', 'LineWidth', .8, 'DisplayName', 'T_wl average');
p6 = plot(x_plot * 1000.*0.03937008, (T_wltol2-273.15)*(9/5) +32, 'm:', 'LineWidth', .8, 'DisplayName', 'T_wl double average');
% Set y-axis label and color
ylabel('Temperature [F]');
set(gca, 'Ycolor', 'k');
ylim([195 410])
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
title('Temperature Variation: Wall Thickness Tolerance');
legend([p1 p2 p3 p4 p5 p6]); % Show legend with specified display names

% Customize legend position (optional)
legend('Location', 'Best');

% Adjust plot margins for better aesthetics (optional)
a = gca


% You can save the plot as an image file if needed (optional)
% saveas(gcf, 'temperature_plot.png');