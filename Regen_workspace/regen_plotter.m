
% Water boiling point vector
watertable = readmatrix("water_boilingpoint.xlsx", "Range","A5:F98");
boiling_temp = zeros(size(P_l));
points_int = size(P_l);
points = points_int(2);
for i = 1: points
    boiling_temp(i) = interp1(watertable(:,2),watertable(:,6), P_l(i)* 1/6894.757 ,'linear','extrap');
end

% Plot heat transfer
mm2inch = 0.03937008;
close all;
f1 = figure('Name', 'Temperature Plot');
allowable_temperature = 400;
% Hold the plot to allow multiple data series to be plotted on the same axes
hold on;

% Set the font style for the plot
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
set(gcf,'color','#eeeeee');
%set(gca,'color','#eeeeee');

% Temperature plot with T_l removed for a more compact plot
yyaxis left; % Use dual y-axes
p1 = plot(x_plot * 1000.*0.03937008, (T_wg-273.15)*(9/5) +32, 'r-', 'LineWidth', 2.5, 'DisplayName', 'T_wg');
p2 = plot(x_plot * 1000.*0.03937008, (T_wl-273.15)*(9/5) +32, 'm-', 'LineWidth', 2, 'DisplayName', 'T_wl');

% Draw a horizontal dashed line at the allowable temperature threshold
yline(allowable_temperature, '--r', 'Design temperature threshold (400°F)', 'LineWidth', 1.5);

% Set y-axis label and color
ylabel('Temperature [F]');
set(gca, 'Ycolor', 'k');
ylim([150 450])
%yticks([350 375 400 425 450 allowable_temperature 500])

yyaxis right
p3 = plot(x_plot .* 1000.*0.03937008, r_interpolated .* 1000.*0.03937008, 'black', 'LineStyle', '-',"DisplayName","chamber contour");
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



f2 = figure('Name', 'Heat Transfer Plots');
subplot(2,2,[1,2])
hold on;
set(gca, 'FontName', 'Times New Roman')
% heat flux plot
%subplot(2,1,2)
yyaxis left
plot(x_plot .* 1000, qdot_g ./ 1000, 'red', 'LineStyle', '-');
ylabel('Heat Flux [kW/m^2]')
set(gca, 'Ycolor', 'k')
grid on

yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;

legend('Convective Heat Flux', 'Chamber Contour','Location','best')
title('Heat Flux Distribution')
xlabel('Location [mm]')

subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.*1000, h_g)
title("Gas Film Coeffcient [W/m^2-K]")
xlabel("Location [mm]");
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.*1000, h_l)
title("Liquid Film Coefficient [W/m^2-K]")
grid on


f3 = figure('Name','Water Flow');
subplot(2,2,[1,2])
hold on 
set(gca, 'FontName', 'Times New Roman')
set(gcf,'color','#eeeeee');
plot(x_plot.* 1000.*0.03937008, P_l * 1/6894.757, "LineWidth",1.5)
plot(x_plot.*1000.*0.03937008, P_g *1/6894.757, 'y',"LineWidth",1.5)
title("Coolant Pressure")
xlabel("Location [inches]")
ylabel("Pressure [PSI]")
yyaxis right
plot(x_plot .* 1000.*0.03937008, r_interpolated .* 1000.*0.03937008, 'black', 'LineStyle', '-');
ylabel('Radius [inches]');
set(gca, 'Ycolor', 'k');
ylim([.5 3.5]);
legend('Location', 'Best');
axis equal;
legend("Pressure curve","Engine pressure"," Chamber contour",'Location','best')
xlim([-7.5 2.5])
yticks([0 100 200 300 400 500])
grid on

subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000.*0.03937008, v.*3.28084,"LineWidth",1.2);
title("Coolant Velocity [ft/s]")
xlabel("Location [inches]")
yline(200,"--b","allowable velocity","LineWidth",1)
grid on
ylim([0 250])
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.*1000.*0.03937008, (T_l-273.15)*(9/5) +32, "LineWidth",1.2)
%plot(x_plot.*1000.*0.03937008, )
plot(x_plot.*1000.*0.03937008,boiling_temp, '--b', 'LineWidth', 1);
legend("Coolant temperature","Boiling point")
title("Coolant Temperature [F]");
xlabel("Location [inches]")
grid on




f4 = figure('Name','Channel Geometry');
hold on
grid on
set(gca, 'FontName', 'Times New Roman')
set(gcf,'color','#eeeeee');
plot(x_plot.* 1000.*0.03937008, w_c_x .*1000.*0.03937008, "LineWidth",1.5);
plot(x_plot.*1000.*0.03937008, h_c_x .*1000.*0.03937008, "LineWidth",1);
plot(x_plot.*1000.*0.03937008, t_w_x .* 1000.*0.03937008, "LineWidth",1);
plot(x_plot.*1000.*0.03937008, rib_thickness .*1000.*0.03937008, "LineWidth",1);
title("Channel Dimensions");
xlabel("Location [inches]");
ylabel("Channel Dimensions [inches]")
% ylim([0 3.5])
% yyaxis right
% p3 = plot(x_plot .* 1000.*0.03937008, r_interpolated .* 1000.*0.03937008, 'black', 'LineStyle', '-',"DisplayName","Chamber contour");
% ylabel('Radius [inches]')
% set(gca, 'Ycolor', 'k')
% ylim([.5 3.5])
legend('Channel width', 'Channel height', 'Wall thickness', 'Rib thickness');
legend('Location', 'Best');

    
f5 = figure('Name', 'Structural results (Stress)')
subplot(2,2,[1,2])
hold on
set(gca, 'FontName', 'Times New Roman');
set(gcf,'color','#eeeeee');
plot(x_plot.* 1000.*mm2inch, sigma_vl.*0.000001.*0.1450377377,'g',x_plot.* 1000.*mm2inch, sigma_vc*0.0000010.*.1450377377, "LineWidth",1.5)
plot(x_plot.* 1000.*mm2inch, yield.*.1450377377, "LineWidth", 1.5)
ylim([0 50])
title("Von Mises Stress: Yield Criterion")
xlabel("Location [inches]")
ylabel("[ksi]") 
yyaxis right
plot(x_plot .* 1000.*mm2inch, r_interpolated .* 1000.*mm2inch, 'black', 'LineStyle', '-');
ylabel('Radius [inches]')
ylim([.5 3.5])
set(gca, 'Ycolor', 'k')
axis equal;
legend('Von mises stress (Lands)','Von mises stress (channels)', 'Yield stress at temperature', 'Chamber contour','Location','best')
grid on
subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
set(gcf,'color','#eeeeee');
hold on
plot(x_plot.* 1000.*mm2inch, sigma_t*0.000001.*.1450377377,"b", "LineWidth",1.2);
plot(x_plot.* 1000.*mm2inch, sigma_tp*0.000001.*.1450377377,"m","LineWidth",1.2);
plot(x_plot.* 1000.*mm2inch, sigma_tt*0.000001.*.1450377377,"r","LineWidth",1.2);
legend("Total stress", "Pressure contribution", "Thermal contribution",'Location','best')
title("Tangential Stress [ksi]")
ylim([0 8])
xlabel("Location [inches]")
hold off
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
set(gcf,'color','#eeeeee')
plot(x_plot.* 1000.*mm2inch, sigma_lc*0.000001.*.1450377377,x_plot.* 1000.*mm2inch, sigma_ll*0.000001.*.1450377377, "LineWidth",1.2);
title("Longitudinal Stress [ksi]")
legend("At the channel", "At the lands",'Location','best');
xlabel("Location [inches]")
ylim([15 40])
grid on
% subplot(2,2,3)
% plot(x_plot.* 1000, sigmab*0.000001)
% title("Buckling Stress (MPA)")
% xlabel("Location [mm]")
f6 = figure('Name', 'Structural results (Strain)');
subplot(2,2,[1,2])
hold on
set(gca, 'FontName', 'Times New Roman')
set(gcf,'color','#eeeeee');
plot(x_plot.* 1000.*mm2inch, epsilon_vl*100,'g',x_plot.* 1000.*mm2inch, epsilon_vc*100,"LineWidth",1.5)
title("Von Mises Strain")
xlabel("Location [inches]")
ylabel("[%]")
ylim([0 .45])
yyaxis right
plot(x_plot .* 1000.*mm2inch, r_interpolated .* 1000.*mm2inch, 'black', 'LineStyle', '-');
ylabel('Radius [inches]')
set(gca, 'Ycolor', 'k')
ylim([.5 3.5])
axis equal;
legend('Von mises strain (lands)','Von mises strain (channels)', 'Chamber contour','Location','best')
grid on
subplot(2,2,3)
hold on 
set(gcf,'color','#eeeeee');
hold on
plot(x_plot.* 1000.*mm2inch, epsilon_t*100,"b","LineWidth",1.2);
plot(x_plot.* 1000.*mm2inch, epsilon_tp*100,"m","LineWidth",1.2);
plot(x_plot.* 1000.*mm2inch, epsilon_tt*100,"r","LineWidth",1.2);
legend("Total strain", "Pressure contribution", "Thermal contribution",'Location','best')
title("Tangential Strain [%]")
ylim([0 .08])
xlabel("Location [inches]")
hold off
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
set(gcf,'color','#eeeeee');
plot(x_plot.* 1000.*mm2inch, epsilon_lc*100,x_plot.* 1000.*mm2inch, epsilon_ll*100,"LineWidth",1.2);
title("Longitudinal Strain [%]")
legend("At the channel", "At the lands",'Location','best');
xlabel("Location [inches]")
ylim([.15 .40])
grid on

f7 = figure('Name', 'Fin results');
subplot(2,2,[1,2])
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, fin_q / 1000,'g')
title("Heat Flux")
xlabel("Location [mm]")
ylabel("[kW/m^2]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Fin Heat Flux', 'Chamber Contour','Location','best')
grid on
subplot(2,2,[3,4])
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, eta_fin);
plot(x_plot.* 1000, h_c_x./rib_thickness);
title("Fin Efficiency")
legend("Fin Efficiency", "Aspect Ratio",'Location','best');
xlabel("Location [mm]")
hold off
grid on
