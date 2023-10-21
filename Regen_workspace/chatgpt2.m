% Create a figure with a 16:9 aspect ratio and specified position
figure('Name', 'Channel Geometry', 'Position', [100, 100, 960, 540]);

% Plot channel dimensions with specified line styles and colors
hold on;
plot(x_plot * 1000, w_c_x * 1000, 'LineWidth', 2.5, 'DisplayName', 'Channel Width', 'Color', [0, 0.4470, 0.7410]);
plot(x_plot * 1000, h_c_x * 1000, 'LineWidth', 2, 'DisplayName', 'Channel Height', 'Color', [0.8500, 0.3250, 0.0980]);
plot(x_plot * 1000, t_w_x * 1000, 'LineWidth', 2, 'DisplayName', 'Wall Thickness', 'Color', [0.9290, 0.6940, 0.1250]);
plot(x_plot * 1000, rib_thickness * 1000, 'LineWidth', 2, 'DisplayName', 'Fin Thickness', 'Color', [0.3010, 0.7450, 0.9330]);

% Plot Chamber Contour using a different line style and color
yyaxis right;
plot(x_plot * 1000, r_interpolated * 1000, 'LineWidth', 2, 'DisplayName', 'Chamber Contour', 'Color', [0.4940, 0.1840, 0.5560], 'LineStyle', '--');

hold off;

% Customize the plot appearance
title("Channel Dimensions", 'FontSize', 16, 'FontWeight', 'bold');
xlabel("Location [mm]", 'FontSize', 14);
ylabel("Channel Dimensions [mm]", 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
yyaxis right;
ylabel('Chamber Contour [mm]', 'FontSize', 14);
set(gca, 'Ycolor', 'k');
axis equal;
legend('show', 'Location', 'northwest', 'FontSize', 12);
grid on;

% Set background color of the plot to white for a clean look
set(gcf, 'Color', 'w');