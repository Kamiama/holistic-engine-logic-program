% %% HELP Regenerative Cooling Sizing Code
% % Authors: Kamon Blong (kamon.blong@gmail.com), Jan Ayala, Andrew Radulovich, Alex Suppiah
% % First Created: 10/23/2022
% % Last Updated: 04/15/2023
% 
% function [] = sizeRegen(x_contour, r_contour, L_seg, R_t, nozzle_regen_pct, m_dot, P_c_max, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF)

   %{ 
    Description:
    This program sets fluid and material inputs and then iterates starting from a
    variable distance downstream of the throat and carries it through up to the injector face. 
    Channels have a constant width but have the option to change channel height to 
    decrease hydraulic diameter and increase fluid velocity (thereby increasing cooling
    effectiveness) in key locations. This method allows for optimal machining via CSJ
    with a slitter-saw CNC attatchment or via powder bed 3D printing.

    Program Methodology
    - take inputs from main sizing code
    - calculate optimal number of channels
    - prescribe channel inputs (pressure, temperature, mdot, ect)
    - run CEA and get bartz film coefficient axially along chamber
    - calculate minimum wall thickness at throat
    - begin iterating along channel length from nozzle end to the injector face
        - for each point, iterate through aspect ratio to find the maximum
        channel size possible without melting (below a certain prescribed size)
        - integrate heat entry and pressure loss for the given step size
        - repeat process all the way up the channel, checking for structural
        stability at each point - if there is a point where the channels overheat
        or the working fluid boils, lower the aspect ratio cap and repeat
    - run program again with new pressure drop value until the calculated exit 
    pressure equals the guessed exit pressure
    - once iterative sizing is complete, display axial wall thickness and
    channel depth on main engine contour plot
    
    Inputs:
    - x_contour: 
    - y_contour:
    - R_t: Throat Radius
    - nozzle_regen_pct: 
    - mdotf: Mass flow rate of fuel/coolant (lb/s)
    - P_c: Chamber Pressure
    - P_e: Exit Pressure
    - Oxidizer: 
    - Fuel:
    - OF_ratio: 
    - wall_material:
    
    Outputs: 
    - 
    
    Assumptions:
    - Steady state
    - No backside wall heat transfer
    - Equally distributed temperature inside channels
    - Wicking heat into fuel doesn't change bulk gas temperature

   %}

% bugfixing variables

clear;
clc; 
close all;

u = convertUnits;
contour = readmatrix('contour.xlsx');

% load("r_contour.mat");
% load("x_contour.mat");
r_contour = (contour(:,2) * u.IN2M)';
x_contour = (contour(:,1) * u.IN2M)';
R_t = .706; % [in]
P_c = 275; % [psi]
P_e = 18.5; % [psi]
fuel = 'C3H8O,2propanol';
fuel_weight = 0;
fuel_temp = 273.15; % [K]
oxidizer = 'O2(L)';
oxidizer_temp = 90.17; % [K]
OF = 1.3;

%L_seg = 0.0283; % [in]
coolant = 'Water';

%inlet_temperature = 293.16;
inlet_temperature = 303.862;
%inlet_pressure = 500 * u.PSI2PA; % [Pa]
inlet_pressure = 3322910; %[PA]


%% Initialize Variables

plots = 0;
coolantdirection = 0; % 1: Direction opposite of hot gas flow direction
                      % 0: Direction same as hot flow gas direction
% convert imperial units to metric

R_t = R_t * u.IN2M; % throat radius [m]
A_t = R_t ^ 2 * pi; % throat area [m^2]
P_c = P_c * u.PSI2PA; % chamber pressure [Pa]
P_e = P_e * u.PSI2PA; % exit pressure [Pa]

CEA_input_name = 'AAAAAA';


%% Parse variables

    % geometric properties
        A_local = []; % Local Cross Sectional Areas of Engine
        dlocal = [];  % Local Cross Sections of Channels
        characteristic_length = 10; % Longest length of Channels 
        %chamber_length = 7; %in
        D_t = R_t * 2; % diameter at nozzle throat [m]
        R_of_curve = 1.5 * D_t / 2; % [m]
        A_local = pi * (R_t) ^ 2; % local cross sectional areas of engine
    % temporary geometry initialization
        chamber_length = 0.0254 * 4.5427; % Chamber Length (m) [conversion * in]
        converging_length = 0.0254 * 1.8278; % Converging Length (m) [conversion * in]
        diverging_length = 0.0254 * 1.7191; % Diverging Length (m) [conversion * in]
        total_length = chamber_length + converging_length + diverging_length; % Total length (mm) [conversion * in]
        
        %[~, steps] = size(r_contour)
    % Discretize Length 
        steps = 2000; % Number of steps along chamber (Change resolution of simulation)
        deltax = (total_length/steps) % Change in distance per step [m]
        points = steps + 1; % Number of points along chamber
        

        x = 0:deltax:total_length; %Length Vector
        x_plot = (x - chamber_length - converging_length); %Length Vector adjusted so that 0 is at the throat (mm)
        
        %Initialize section length vectors
        x_chamber = [];
        x_converging = [];
        x_diverging = [];
        for i = x 
            if i <= chamber_length
            x_chamber = [x_chamber i];
            end 
            if (chamber_length < i) && (i <= chamber_length + converging_length)
            x_converging = [x_converging i];
            end 
            if i > (chamber_length + converging_length)
            x_diverging = [x_diverging i];
            end 
        end 

    % channel geometry: (1: chamber) (min: throat) (2: nozzle end)
        t_w = 0.001;        % [m]
        h_c = [.0035 .0014 .0035]; % [1 min 2] [m]    
        w_c = [.0055 0.0016 .004];     % [1 min 2] [m]
        %h_c = .0014;
        %w_c = .0016;
        %length = L_seg * u.IN2M;      % [m], arbitrary for now
        
        A = w_c .* h_c; % Channel Cross-sectional Area (m^2) [1 min 2]
        p_wet = 2*w_c + 2*h_c; % Wetted Perimeter of the pipe (m) [1 min 2]
        hydraulic_D = (4.*A)./p_wet; % hydraulic Diameter (m) [1 min 2]

        %Set channel width/height over channel length
        w_c_chamber = ones(1,size(x_chamber,2)).*w_c(1); % Channel width over Chamber length (constant)
        w_c_converging = ((w_c(2)-w_c(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
                 + ones(1,size(x_converging,2)).*w_c(1); % Channel width over converging length (linear interpolation)
        w_c_diverging = ((w_c(3)-w_c(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
                + ones(1,size(x_diverging,2)).*w_c(2);   % Channel width over diverging length (linear interpolation)
        w_c_x = [w_c_chamber w_c_converging w_c_diverging]; % Combine channel width vectors

        h_c_chamber = ones(1,size(x_chamber,2)).*h_c(1); % Channel height over Chamber length (constant_
        h_c_converging = ((h_c(2)-h_c(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
             + ones(1,size(x_converging,2)).*h_c(1);    % Channel height over converging length (linear interpolation)
        h_c_diverging = ((h_c(3)-h_c(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
             + ones(1,size(x_diverging,2)).*h_c(2); % Channel height over diverging length (linear interpolation)
        h_c_x = [h_c_chamber h_c_converging h_c_diverging]; % Combine channel height vectors

        A_x = (w_c_x .* h_c_x); % Channel area vector over channel length [m^2]
        p_wet_x = 2.*w_c_x + 2 .* h_c_x; % Wet perimeter over channel length [m]
        hydraulic_D_x = ((4.*(A_x))./p_wet_x); % Hydraulic Diameter over channel length [m]

        
       




   
   
    % working fluid properties
        

    % engine properties
        m_dot = 12 * u.LB2KG; % Coolant/fuel mass flow [kg/s], 1.2566

    % wall material properties
        k_w = 103; % Thermal Conductivity of Wall [W/m-K]
        E = 70E9; % in Pa
        CTE = 27E-6; % in 1/K
        nu = 0.3; 
        e = 24 * 0.001; % surface roughness (mm) [micrometer*conversion]

        qdot_tolerance = 0.0001;

%% THROAT CALCULATIONS FOR WALL THICKNESS
% determine optimal wall thickness based on structural calculations under throat conditions
% uses the incorrect assumption of nozzle inlet coolant properties

% Step 1: Prescribe initial properties
[c_star, ~, ~, M, gamma, P_g, T_g, ~, mu_g, Pr_g, ~, ~, ~, cp_g] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 1, 0, 0, CEA_input_name);

% Steps 2 & 3: Set channel inlet properties
num_channels = 40; % pi * (D_t + 0.8 * (D_t + 2 * t_w)) / (D_t + 2 * t_w); % number of Channels (EQ 6.30) (Change the coefficient later)
m_dot_CHANNEL = m_dot / num_channels; % mass flow of channel (EQ 6.31)


P_l = inlet_pressure; 
T_l = inlet_temperature; 

% Step 4: Take hot wall temperature guess and initialize loop

T_wg = 1000; % initial guess of wall side temperature [K]
T_wg_mn = 294.15; % minimum temperature bound
T_wg_mx = 2000; % maximum temperature bound

converged = 0; % wall temperature loop end condition
structurally_sound = 0; % wall thickness loop end condition
counter = 0; % counter for loop

% start loop with a minimum wall thickness value and increase wall thickness until structurally sound
%while ~(structurally_sound)

    % start loop to converge on wall temperature
    while ~(converged)
        % Step 5: Calculate gas film coefficient and gas-side convective heat flux
        sigma = (.5 * T_wg / T_g * (1 + (gamma - 1) / 2 * M ^ 2) + .5) ^ -.68 * (1 + (gamma - 1) / 2 * M ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
        h_g = (0.026 / D_t ^ 0.2) * (mu_g ^ 0.2 * cp_g / Pr_g ^ 0.6) * (P_g / c_star) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (A_t / A_t) ^ .9 * sigma; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
        r = Pr_g ^ (1 / 3); % recovery factor for a turbulent free boundary layer [N/A] - biased towards larger engines, very small engines should use Pr^.5 (Heister Table 6.2).
        T_r = T_g * (1 + (gamma - 1) / 2 * r * M ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 
        qdot_g = h_g * (T_r - T_wg); % gas convective heat flux [W/m^2] (Heister EQ 6.16).
    
        % Step 6: Calculate liquid wall temperature
        T_wl = T_wg - qdot_g * t_w / k_w; % liquid wall temperature calculated via conduction through wall [K] (Heister EQ 6.29).
    
        % Step 7: Calculate liquid film coefficient
        % run coolprop to get coolant properties
        mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l, 'P', P_l, coolant); % viscosity of bulk coolant [Pa-s]
        cp_l = py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l, 'P', P_l, coolant); % specific heat of coolant [J/kg-k] 
        k_l = py.CoolProp.CoolProp.PropsSI('L', 'T', T_l, 'P', P_l, coolant); % thermal conductivity of coolant [W/m-K]
            
        Re_l = (4 * m_dot_CHANNEL) / (pi * hydraulic_D(2) * mu_lb); % reynolds number for channel flow [] ALEX CITE SOURCE
        Pr_l = (cp_l * mu_lb) / k_l; % prantl number [] ALEX CITE SOURCE
        Nu_l = 0.023 * (Re_l ^ .8) * (Pr_l ^ .4) * (T_wl / T_l) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19).
        h_l = (Nu_l * k_l) / hydraulic_D(2); % liquid film coefficient [W/m^2-K] ALEX CITE SOURCE
        
        % Step 8: Calculate liquid-side convective heat flux
        qdot_l = h_l * (T_wl - T_l); % liquid convective heat flux [W/m^2] (Heister EQ 6.29).
    
        % Step 9: Check for convergence and continue loop / next step
        if abs(qdot_g - qdot_l) > qdot_tolerance % check for tolerance
    
            % convergence loop
            if qdot_g - qdot_l > 0
                T_wg_mn = T_wg;
            else 
                T_wg_mx = T_wg;
            end 
            T_wg = (T_wg_mx + T_wg_mn) / 2;
    
            counter = counter + 1;
        else 
            fprintf("Gas Side Wall Temp [K]: %0.2f\n", T_wg)
            
            converged = 1; % end loop
        end
    end
% 
%     % structural calculations for channel geometry...
%     % if structurally sound, end loop & apply saftey factor
%     % else, continue loop & increase wall thickness
% %end
% 
%% LOOP ALONG CHAMBER LENGTH
inlet_temperature = 293.16;
inlet_pressure = 500 * u.PSI2PA; % [Pa]


% Step 1: Prescribe initial properties

% Prescribe area ratios
r_interpolated = interp1(x_contour,r_contour,x_plot); % Linearly Interpolate r vector  
subsonic_area_ratios = (pi * r_interpolated(x_plot < 0) .^ 2) / A_t; % subsonic area ratios on discretized points
supersonic_area_ratios = (pi * r_interpolated(x_plot > 0) .^ 2) / A_t; %  supersonic area ratios on discretized points
A_ratio = [subsonic_area_ratios, supersonic_area_ratios];
for u = 1:length(subsonic_area_ratios)
    if subsonic_area_ratios(u) < 1
        subsonic_area_ratios(u) = 1;
    end
end


% axial coolant property matrices
P_l = zeros(1, points);
T_l = zeros(1, points);
rho_l = zeros(1, points);
v = zeros(1, points);

% axial cooling property matrices
qdot_l = zeros(1, points);
qdot_g = zeros(1, points);
T_wl = zeros(1, points);
T_wg = zeros(1, points);
h_g = zeros(1, points);
h_l = zeros(1, points);

% axial channel geometric property matrices
% A = ones(1, steps) .* A;


% axial combustion property matrices
c_star = zeros(1, points);
M = zeros(1, points);
gamma = zeros(1, points);
P_g = zeros(1, points);
T_g = zeros(1, points);
mu_g = zeros(1, points);
Pr_g = zeros(1, points);
cp_g = zeros(1, points);

i = 1;
for sub = subsonic_area_ratios
    [c_star(i), ~, ~, M(i), gamma(i), P_g(i), T_g(i), ~, mu_g(i), Pr_g(i), ~, ~, ~, cp_g(i)] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, sub, 0, 2, 0, 0, CEA_input_name);
    i = i + 1;
end
i = size(subsonic_area_ratios, 2) + 1;
for sup = supersonic_area_ratios
    [c_star(i), ~, ~, M(i), gamma(i), P_g(i), T_g(i), ~, mu_g(i), Pr_g(i), ~, ~, ~, cp_g(i)] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, sup, 2, 0, 0, CEA_input_name);
    i = i + 1;
end

% Steps 2 & 3: Set channel inlet properties
P_l(1) = inlet_pressure;
T_l(1) = inlet_temperature;

% Step 4: Take hot wall temperature guess and initialize loop
T_wg(1) = 1000; % initial guess of wall side temperature [K]

% perform cooling loop along the chamber
for i = 1:points % where i is the position along the chamber (1 = injector, end = nozzle)

    T_wg_mn = 280; % minimum temperature bound [K]
    T_wg_mx = 1500; % maximum temperature bound [K]
    
    converged = 0; % wall temperature loop end condition
    counter = 0; % counter for loop
    
    while ~(converged) && counter < 250
        % Step 5: Calculate gas film coefficient and gas-side convective heat flux
        sigma = (.5 * T_wg(i) / T_g(i) * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) + .5) ^ -.68 * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
        h_g(i) = (0.026 / D_t ^ 0.2) * (mu_g(i) ^ 0.2 * cp_g(i) / Pr_g(i) ^ 0.6) * (P_g(i) / c_star(i)) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / A_ratio(i)) ^ .9 * sigma; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
        r = Pr_g(i) ^ (1 / 3); % recovery factor for a turbulent free boundary layer [N/A] - biased towards larger engines, very small engines should use Pr^.5 (Heister Table 6.2).
        T_r = T_g(i) * (1 + (gamma(i) - 1) / 2 * r * M(i) ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 
        qdot_g(i) = h_g(i) * (T_r - T_wg(i)); % gas convective heat flux [W/m^2] (Heister EQ 6.16).
    
        % Step 6: Calculate liquid wall temperature
        T_wl(i) = T_wg(i) - qdot_g(i) * t_w / k_w; % liquid wall temperature calculated via conduction through wall [K] (Heister EQ 6.29).
    
        % Step 7: Calculate liquid film coefficient
        % run coolprop to get coolant properties
        mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l(i), 'P', P_l(i), coolant); % viscosity of bulk coolant [Pa-s]
        cp_l = py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l(i), 'P', P_l(i), coolant); % specific heat of coolant [J/kg-k] 
        k_l = py.CoolProp.CoolProp.PropsSI('L', 'T', T_l(i), 'P', P_l(i), coolant); % thermal conductivity of coolant [W/m-K]

        Re_l = (4 * m_dot_CHANNEL) / (pi * hydraulic_D_x(i) * mu_lb); % reynolds number for channel flow [] ALEX CITE SOURCE (I dont like this - Andrew)
        Pr_l = (cp_l * mu_lb) / k_l; % prantl number [] ALEX CITE SOURCE
        Nu_l = 0.023 * (Re_l ^ .8) * (Pr_l ^ .4) * (T_wl / T_l) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19).
        h_l = (Nu_l * k_l) / hydraulic_D_x(i); % liquid film coefficient [W/m^2-K] ALEX CITE SOURCE
    
        % Step 8: Calculate liquid-side convective heat flux
        qdot_l(i) = h_l * (T_wl(i) - T_l(i)); % liquid convective heat flux [W/m^2] (Heister EQ 6.29).

        % Step 9: Check for convergence and continue loop / next step
        if abs(qdot_g(i) - qdot_l(i)) > qdot_tolerance % check for tolerance
            
            % convergence loop
            if qdot_g(i) - qdot_l(i) > 0
                T_wg_mn = T_wg(i);
            else 
                T_wg_mx = T_wg(i);
            end 
            T_wg(i) = (3 * T_wg_mx + T_wg_mn) / 4;
    
            counter = counter + 1;
        else 
            if i < points
                % Step 10: End step & update fluid properties
                wall_area = w_c_x(i) * deltax;
                T_l(i+1) = T_l(i) + (1 / (m_dot_CHANNEL * cp_l)) * qdot_g(i) * wall_area; % new liquid temperature [K] (Heister EQ 6.39).
              
                % calculate skin friction factor
 
                rho_l(i) = py.CoolProp.CoolProp.PropsSI('D','T', T_l(i),'P', P_l(i),'Water');
                v(i) = m_dot_CHANNEL / rho_l(i) / A_x(i); % velocity at step [m/s]
                
                % Use moody diagram to find coefficient of friction
                ed = e/(hydraulic_D_x(i)*1000); 
                Re2 = (rho_l(i) * v(i)*(hydraulic_D_x(i)))/ mu_lb;
                f = moody(ed, Re2);
                cf = f/4;


                deltaP =   (2*cf*(deltax/(hydraulic_D_x(i))) * rho_l(i) *(v(i))^(2)); % Calculate change in pressure
                
                %deltaP =   (2*cf*(deltax/(hydraulic_D_x(i))) * rho
                %*(v(i))^(2)  + .5 * ((v_x(i)^2) -(v_x(i-1)^2))) % may need
                %to put this in different part of loop if we want to use
                %full bernoulli's


                %P_l(i+1) = P_l(i) - cf * deltax / hydraulic_D_x(i) * 2 * rho_l(i) * v(i)^2; % new liquid pressure (Heister EQ 6.36).
                % fprintf("Gas Side Wall Temp [K]: %0.2f\n", T_wg(i))
                %a = x_plot(i);
        
                % prepare for next step
                P_l(i+1) =   P_l(i) - deltaP; % Update pressure for next iteration
                T_wg(i+1) = T_wg(i);
            end
            
            converged = 1;
        end
    end
end

%% PLOT OUTPUTS

figure('Name', 'Temperature Plot');
hold on;

% temperature plot
% subplot(2,1,1)
yyaxis left
plot(x_plot .* 1000, T_wg, 'red', 'LineStyle', '-');
plot(x_plot .* 1000, T_wl, 'magenta', 'LineStyle', '-');
plot(x_plot .* 1000, T_l, 'blue', 'LineStyle', '-');
ylabel('Temperature [K]')
set(gca, 'Ycolor', 'k')
grid on

yyaxis right
%plot(x_contour .* 1000, r_contour .* 1000, 'black', 'LineStyle', '-');
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;

legend('T_w_g', 'T_w_l', 'T_l', 'Chamber Contour', 'Location', 'southoutside', 'Orientation', 'horizontal')
title('Temperature Distribution')
xlabel('Location [mm]')


figure(2)


subplot(2,2,1)
plot(x_plot.* 1000, A_x * 1000000);
title("Channel Area [mm^2]")
xlabel("Location [mm]")
subplot(2,2,2)
plot(x_plot.* 1000, v);
title("Coolant Velocity [m/s]")
xlabel("Location [mm]")
subplot(2,2,3)
plot(x_plot.* 1000, P_l * 1/6894.757)
title("Liquid Pressure [psi]")
xlabel("Location [mm]")

%figure('Name', 'Heat Flux Plot');
%hold on;
% % heat flux plot
% %subplot(2,1,2)
% yyaxis left
% plot(x_contour .* 1000, flip(qdot_g) ./ 1000, 'red', 'LineStyle', '-');
% ylabel('Heat Flux [kW/m^2]')
% set(gca, 'Ycolor', 'k')
% grid on
% 
% yyaxis right
% plot(x_contour .* 1000, r_contour .* 1000, 'black', 'LineStyle', '-');
% ylabel('Radius [mm]')
% set(gca, 'Ycolor', 'k')
% axis equal;
% 
% legend('Convective Heat Flux', 'Chamber Contour', 'Location', 'southoutside', 'Orientation', 'horizontal')
% title('Heat Flux Distribution')
% xlabel('Location [mm]')
    
%         %Step 12: Structural Analysis Checks
%         St = .5*(Pl- P_gas)*(width/wallthick)^2 + (E*a*gasheattransfer*wallthick)/(2*(1-v)kw); %Combined tangential stresses: Heister Eq 6.33 page 207
%         Sl = E*a*(Twl-fuel_temp); %Longtudinal thermal stress (Huzel and Huang, EQ 2-28, pg 92) The temperatures used here may not right for determining the delta T in this equation.
%         Sc = 4*Et*Ec*t_w/((((Et)^(1/2))*((Ec)^(1/2))^2)*(3*(1-v^2)*tube_radius)); %Critical Stress Buckling (Huzel and Huang, Eq 4-29, pg 93)

%% THERMAL FEA
if plots

    M = 150;
    N = 150;
    R1 = R_t; % inner radius 
    R2 = R_t + h_c * 4;  % outer radius
    nR = linspace(R1,R2,M);
    nT = linspace(-pi/num_channels, pi/num_channels + w_c / R_t, N);
    [R, T_g] = meshgrid(nR,nT) ;
    xg = R.*cos(T_g); 
    yg = R.*sin(T_g);
    xg = xg(:);
    yg = yg(:);
    
    % Define partial channel 
    M = 50;
    N = 50;
    R1 = R_t + t_w; % inner radius 
    R2 = R_t + t_w + h_c;  % outer radius
    x = linspace(R1,R2,M);
    y = linspace(-pi/num_channels - w_c / R_t / 2, -pi/num_channels + w_c / R_t / 2, N);
    [R, T_g] = meshgrid(x, y);
    x = R.*cos(T_g); 
    y = R.*sin(T_g);
    x = x(:);
    y = y(:);
    channel = alphaShape(x,y);
    in = inShape(channel,xg,yg);
    xg = xg(~in);
    yg = yg(~in);
    
    % Define full channel 
    M = 50;
    N = 50;
    R1 = R_t + t_w; % inner radius 
    R2 = R_t + t_w + h_c;  % outer radius
    x = linspace(R1,R2,M);
    y = linspace(pi/num_channels - w_c / R_t / 2, pi/num_channels + w_c / R_t / 2, N);
    [R, T_g] = meshgrid(x, y);
    x = R.*cos(T_g); 
    y = R.*sin(T_g);
    x = x(:);
    y = y(:);
    channel = alphaShape(x,y);
    in = inShape(channel,xg,yg);
    xg = xg(~in);
    yg = yg(~in);
    
    zg = ones(numel(xg),1);
    xg = repmat(xg,5,1);
    yg = repmat(yg,5,1);
    zg = zg*linspace(0,length,5);
    zg = zg(:);
    shp = alphaShape(xg,yg,zg);
    
    [elements,nodes] = boundaryFacets(shp);
    
    nodes = nodes';
    elements = elements';
    
    % Generate model
    model = createpde("thermal","steadystate");
    geometryFromMesh(model,nodes,elements);
    
    pdegplot(model,"FaceLabels","on","FaceAlpha",0.5)
    
    generateMesh(model,"Hmax",h_c/12);
    % figure
    % pdemesh(model)
    
    % Define material thermal properties
    thermalProperties(model,"ThermalConductivity",k_w);
    
    % Thermal boundary conditions
    thermalBC(model,"Face",10, ...
                     "ConvectionCoefficient",h_g, ...
                     "AmbientTemperature",T_r);
    thermalBC(model,"Face",[5 13 7 3 11 12 14], ...
                     "ConvectionCoefficient",h_l, ...
                     "AmbientTemperature",T_l(1));
    Rt = solve(model);
    
    figure
    pdeplot3D(model,"ColorMapData",Rt.Temperature)
    view([-90,90]);
    
    maxTempFEA = max(Rt.Temperature)
    
    %% Structural FEA
    model = createpde("structural","static-solid");
    geometryFromMesh(model,nodes,elements);
    generateMesh(model,"Hmax",1.84e-04);
    
    % Material properties
    structuralProperties(model,"YoungsModulus",E, ...
                                 "PoissonsRatio",nu, ...
                                 "CTE",CTE);
    model.ReferenceTemperature = 300 + 273.15; %in degrees K
    structuralBodyLoad(model,"Temperature",Rt);
    
    % Structural boundary conditions
    structuralBC(model,"Face",8,"Constraint","fixed");
    structuralBoundaryLoad(model,"Face",[5 13 7 3 11 12 14],"Pressure",P_l(end));
    structuralBoundaryLoad(model,"Face",10,"Pressure",P_c);
    
    % Solve structural
    Rts = solve(model);
    
    % Display results
    figure("units","normalized");
    hold on
    plot3(x(1),y(1),zg(end), "+", "LineWidth", 2,'Color','r')
    pdeplot3D(model,"ColorMapData",Rts.VonMisesStress, ...
                      "Deformation",Rts.Displacement, ...
                      "DeformationScaleFactor",2)
    view([-90,90]);
    caxis([1e6, 3e8])
    
    channelVonMises(1) = interpolateVonMisesStress(Rts,x(1),y(1),zg(end));
    channelVonMises(2) = interpolateVonMisesStress(Rts,x(1),y(end),zg(end));
    channelVonMises(3) = interpolateVonMisesStress(Rts,x(end),y(1),zg(end));
    channelVonMises(4) = interpolateVonMisesStress(Rts,xg(round(size(xg,1)-13200)),y(end),zg(end));
    maxVonMisesStressFEA = max(channelVonMises)

end
%%

%plot(x_standard, mu_g);
