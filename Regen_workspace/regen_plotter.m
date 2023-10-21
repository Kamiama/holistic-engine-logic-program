% Plot heat transfer
close all;
allowable_temperature = 477; % Kelvin

figure('Name', 'Temperature Plot');
hold on;
set(gca, 'FontName', 'Times New Roman')
% temperature plot
% subplot(2,1,1)
yyaxis left
p1 = plot(x_plot .* 1000, T_wg, 'red', 'LineStyle', '-',"LineWidth",2.5);
p2 = plot(x_plot .* 1000, T_wl, 'magenta', 'LineStyle', '-');
p3 = plot(x_plot .* 1000, T_l, 'blue', 'LineStyle', '-');
% plot(x_plot .*1000, allowable_temperature, 'red', 'LineStyle','-.')
%p5 = yline(allowable_temperature,'Allowable temperature threshold','red','LineStyle','-.', "LineWidth", 1.5);
yline(allowable_temperature,'-','Allowable temperature threshold (400 F)')
ylabel('Temperature [K]')
set(gca, 'Ycolor', 'k')
grid on

yyaxis right
%plot(x_contour .* 1000, r_contour .* 1000, 'black', 'LineStyle', '-');
p4 = plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;

legend([p1 p2 p3 p4], 'T_w_g', 'T_w_l', 'T_l', 'Chamber Contour', 'Location', 'southoutside', 'Orientation', 'horizontal','Location','best')
title('Temperature Distribution')
xlabel('Location [mm]')

figure('Name', 'Heat Transfer Plots');
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


figure('Name','Water Flow')
subplot(2,2,[1,2])
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, P_l * 1/6894.757)
plot(x_plot.*1000, P_g *1/6894.757, 'y')
title("Liquid Pressure Loss")
xlabel("Location [mm]")
ylabel("Pressure [PSI]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend("Pressure Curve","Gas Pressure"," Chamber Contour",'Location','best')
grid on

subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, v);
title("Coolant Velocity [m/s]")
xlabel("Location [mm]")
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.*1000, T_l)
title("Coolant Temperature [K]");
grid on




figure('Name','Channel Geometry');
subplot(2,2,[1,2]);
hold on
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, w_c_x .*1000);
plot(x_plot.*1000, h_c_x .*1000);
plot(x_plot.*1000, t_w_x .* 1000);
plot(x_plot.*1000, rib_thickness .*1000);
title("Channel Dimensions");
xlabel("Location [mm]");
ylabel("Channel Dimensions [mm]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Chamber Contour [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Channel Width', 'Channel Height', 'Wall Thickness', 'Fin Thickness', 'Chamber Contour','Location','northwest')
grid on

subplot(2,2,3);
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.*1000, h_c_x./rib_thickness);
title("Fin Aspect Ratio [mm]");
xlabel("Location [mm]");
grid on
subplot(2,2,4);
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, AR_channel);
title("Channel Aspect Ratio");
xlabel("Location [mm]");
grid on

    
figure('Name', 'Structural results (Stress)')
subplot(2,2,[1,2])
hold on
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, sigma_vl*0.000001,'g',x_plot.* 1000, sigma_vc*0.000001)
plot(x_plot.* 1000, yield)
title("Von Mises Stress")
xlabel("Location [mm]")
ylabel("[MPA]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Von Mises Stress (Lands)','Von Mises Stress (Channels)', 'Yield Stress at Temperature', 'Chamber Contour','Location','best')
grid on
subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
hold on
plot(x_plot.* 1000, sigma_t*0.000001,"b");
plot(x_plot.* 1000, sigma_tp*0.000001,"m");
plot(x_plot.* 1000, sigma_tt*0.000001,"r");
legend("Total Stress", "Pressure contribution", "Thermal Contribution",'Location','best')
title("Tangential Stress (MPA)")
hold off
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, sigma_lc*0.000001,x_plot.* 1000, sigma_ll*0.000001);
title("Longitudinal Stress (MPA)")
legend("At the channel", "At the lands",'Location','best');
xlabel("Location [mm]")
grid on
% subplot(2,2,3)
% plot(x_plot.* 1000, sigmab*0.000001)
% title("Buckling Stress (MPA)")
% xlabel("Location [mm]")
figure('Name', 'Structural results (Strain)')
subplot(2,2,[1,2])
hold on
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, epsilon_vl*100,'g',x_plot.* 1000, epsilon_vc*100)
title("Von Mises Strain [%]")
xlabel("Location [mm]")
ylabel("[%]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Von Mises Strain (Lands)','Von Mises Strain (Channels)', 'Chamber Contour','Location','best')
grid on
subplot(2,2,3)
hold on 
set(gca, 'FontName', 'Times New Roman')
hold on
plot(x_plot.* 1000, epsilon_t*100,"b");
plot(x_plot.* 1000, epsilon_tp*100,"m");
plot(x_plot.* 1000, epsilon_tt*100,"r");
legend("Total Strain", "Pressure contribution", "Thermal Contribution",'Location','best')
title("Tangential Strain (%)")
hold off
grid on
subplot(2,2,4)
hold on 
set(gca, 'FontName', 'Times New Roman')
plot(x_plot.* 1000, epsilon_lc*100,x_plot.* 1000, epsilon_ll*100);
title("Longitudinal Strain (%)")
legend("At the channel", "At the lands",'Location','best');
xlabel("Location [mm]")
grid on

figure('Name', 'Fin results')
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
