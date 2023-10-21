

% Create a figure with a 16:9 aspect ratio
figure('Name', 'Channel Geometry', 'Position', [100, 100, 960, 540]);
hold on
plot(x_plot * 1000, w_c_x * 1000, 'LineWidth', 2.5, 'DisplayName', 'Channel Width', 'Color', "b");
plot(x_plot * 1000, h_c_x * 1000, 'LineWidth', 2, 'DisplayName', 'Channel Height', 'Color', [0.4660, 0.6740, 0.1880]);
plot(x_plot * 1000, t_w_x * 1000, 'LineWidth', 2, 'DisplayName', 'Wall Thickness', 'Color', [0.9290, 0.6940, 0.1250]);
plot(x_plot * 1000, rib_thickness * 1000, 'LineWidth', 2, 'DisplayName', 'Fin Thickness', 'Color', [0.3010, 0.7450, 0.9330]);


yyaxis right;
plot(x_plot * 1000, r_interpolated * 1000, 'LineWidth', 2, 'DisplayName', 'Chamber Contour', 'Color', [0.4940, 0.1840, 0.5560]);

hold off;

% Customize the plot appearance
title("Channel Dimensions");
xlabel("Location [mm]");
ylabel("Channel Dimensions [mm]");
set(gca, 'FontName', 'Times New Roman');
yyaxis right;
ylabel('Chamber Contour [mm]');
set(gca, 'Ycolor', 'k');
axis equal;
legend('show', 'Location', 'northwest');
grid on;
