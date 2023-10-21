% Create a new figure with a specified name
figure('Name', 'Temperature Plot');

% Hold the plot to allow multiple data series to be plotted on the same axes
hold on;

% Set the font style for the plot
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

% Temperature plot with T_l removed for a more compact plot
yyaxis left; % Use dual y-axes
p1 = plot(x_plot * 1000.*0.03937008, T_wg, 'r-', 'LineWidth', 2.5, 'DisplayName', 'T_wg');
p2 = plot(x_plot * 1000.*0.03937008, T_wl, 'm-', 'LineWidth', 2, 'DisplayName', 'T_wl');

% Draw a horizontal dashed line at the allowable temperature threshold
yline(allowable_temperature, '--r', 'Max allowable wall temperature threshold (400°F)', 'LineWidth', 1.5);

% Set y-axis label and color
ylabel('Temperature [K]');
set(gca, 'Ycolor', 'k');
ylim([350 500])
yticks([350 375 400 425 450 allowable_temperature 500])

yyaxis right
p3 = plot(x_plot .* 1000.*0.03937008, r_interpolated .* 1000.*0.03937008, 'black', 'LineStyle', '-',"DisplayName","Chamber contour");
ylabel('Radius [inches]')
set(gca, 'Ycolor', 'k')
ylim([.5 3.5])

% Enable grid lines for better visualization
grid on;

% Add labels and legend
xlabel('Location [inches]');
title('Temperature Variation');
legend([p1 p2 p3]); % Show legend with specified display names

% Customize legend position (optional)
legend('Location', 'Best');

% Adjust plot margins for better aesthetics (optional)
a = gca


% You can save the plot as an image file if needed (optional)
% saveas(gcf, 'temperature_plot.png');
