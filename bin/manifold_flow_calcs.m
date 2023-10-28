clear
clc
close all
u = convertUnits;

m_dot_array = [3, 4, 5, 6, 7];
for j = 1:1:size(m_dot_array,2)
    % Inlet environment definition
    m_dot = m_dot_array(j) * u.LB2KG /2; % kg/s
    T_l = 293.16; % K
    inlet_pressure = 500 * u.PSI2PA; % inlet pressure [PA]
    rho_l = py.CoolProp.CoolProp.PropsSI('D','T', T_l,'P', inlet_pressure, 'water');
    
    % Manifold Geometry Definition
    A_max = 0.07068583471 * u.IN2M^2;
    A_min = 0.015393804 * u.IN2M^2;  
    manifold_radius = 2; % [in]
    manifold_half_length = pi * manifold_radius * u.IN2M;
    
    % Manifold interpolation
    steps = 100;
    length = linspace(0, steps, steps+1);
    channel_mdot = m_dot - (m_dot - m_dot/steps) ./ steps .* (length); 
    cubic_area = (length + steps/2) .* (length - steps).^2 ./ (steps^2 * steps/2) * (A_max - A_min) + A_min;
    
    % Velocity calculations
    velo(j,:) = channel_mdot ./ (cubic_area * rho_l);
    
    % Pressure loss
    roughness_table = readmatrix(pwd + "/bin/surface_roughness.xlsx",'Range','A12:E16');
    e = roughness_table(1,2) .* 0.001; % Surface roughness (mm) [micrometer*conversion] @ 0 deg
    
    deltax = manifold_half_length / steps; 
    mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l, 'P', inlet_pressure, 'water'); % viscosity of bulk coolant [Pa-s]
    diam = sqrt(cubic_area / pi);
    ed = e./(diam.*1000); % 90 degrees
    
    Re_l = rho_l .* velo .* diam ./ mu_lb;
    
    pressure_drop(j, 1) = inlet_pressure;
    for i = 1:1:steps+1
        f(i) = moody(ed(i), Re_l(i)); % friction factor
        cf(i) = f(i)/4; % friction coefficient
        dP(i) = 2 * cf(i) * (deltax./diam(i)) * rho_l * velo(j, i)^2;
        pressure_drop(j, i+1) = pressure_drop(j, i) - dP(i);
    end
    
end

%% Figures
figure
subplot(2,2,[1 2])
hold on 
grid on 
set(gca, 'FontName', 'Times New Roman')
plot(length/100, linspace(A_max,A_max,steps+1) * u.M2IN^2, "red", 'LineStyle', "--", "Linewidth", 1)
plot(length/100, linspace(A_min,A_min,steps+1) * u.M2IN^2, "black", 'LineStyle', "--", "Linewidth", 1)
plot(length/100, cubic_area * u.M2IN^2, 'LineStyle', "-", "Linewidth", 2)
ylabel("Area (in^2)")
xlabel("Length")
legend("Max Area", "Min Area", "Manifold Contour", "Location", "west")
title("Water Manifold Contour")
ylim([0 A_max*u.M2IN^2+0.01])

subplot(2,2,3)
hold on 
grid on 
set(gca, 'FontName', 'Times New Roman')
yyaxis left
ylabel("Velocity (ft/s)")
plot(length/100, velo(end,:) * u.M2F, "Linewidth", 2)
yyaxis right
ylabel("Pressure (psi)")
plot(length/100, pressure_drop(end, 1:steps+1) * u.PA2PSI, "Linewidth", 2)
legend("Water Velocity", "Pressure Drop", "Location", "west")
xlabel("Length")
title("Water Manifold Pressure Loss (7 lb/s)")

subplot(2,2,4)
hold on 
grid on 
set(gca, 'FontName', 'Times New Roman')
plot(m_dot_array(1), (max(pressure_drop(1, 1:end))-min(pressure_drop(1, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(2), (max(pressure_drop(2, 1:end))-min(pressure_drop(2, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(3), (max(pressure_drop(3, 1:end))-min(pressure_drop(3, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(4), (max(pressure_drop(4, 1:end))-min(pressure_drop(4, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(5), (max(pressure_drop(5, 1:end))-min(pressure_drop(5, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
title("Pressure Loss with Mass Flow")
ylabel("Pressure (psi)")
xlabel("Mass Flow (lb/s)")

figure
subplot(2,2,[1 2])
hold on 
grid on 
set(gca, 'FontName', 'Times New Roman')
plot(m_dot_array(1), (max(pressure_drop(1, 1:end))-min(pressure_drop(1, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(2), (max(pressure_drop(2, 1:end))-min(pressure_drop(2, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(3), (max(pressure_drop(3, 1:end))-min(pressure_drop(3, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(4), (max(pressure_drop(4, 1:end))-min(pressure_drop(4, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
plot(m_dot_array(5), (max(pressure_drop(5, 1:end))-min(pressure_drop(5, 1:end))) * u.PA2PSI, "x", "Linewidth", 2)
title("Pressure Loss with Mass Flow")
ylabel("Pressure (psi)")
xlabel("Mass Flow (lb/s)")

subplot(2,2,[3 4])
hold on 
grid on 
set(gca, 'FontName', 'Times New Roman')
plot(m_dot_array(1), max(velo(1,:)) * u.M2F, "x", "Linewidth", 2)
plot(m_dot_array(2), max(velo(2,:)) * u.M2F, "x", "Linewidth", 2)
plot(m_dot_array(3), max(velo(3,:)) * u.M2F, "x", "Linewidth", 2)
plot(m_dot_array(4), max(velo(4,:)) * u.M2F, "x", "Linewidth", 2)
plot(m_dot_array(5), max(velo(5,:)) * u.M2F, "x", "Linewidth", 2)
title("Max Velocity with Mass Flow")
ylabel("Velocty (in/s)")
xlabel("Mass Flow (lb/s)")
