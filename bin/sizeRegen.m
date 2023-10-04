% %% HELP Regenerative Cooling Sizing Code
% % Authors: Kamon Blong (kamon.blong@gmail.com), Jan Ayala, Andrew Radulovich, Alex Suppiah
% % First Created: 10/23/2022
% % Last Updated: 04/15/2023

   %{ 
    Description:
    This program calculates the heat flux and wall temperture across a
    regen engine. The user inputs engine definition parameters and channel
    inlet conditions. Steady state equilibrium equations are utilized to
    converge on a heatflux and temperature at each point. Calculations at
    the previous point are used as the initial conditions of the next as
    it moves down the engine. The program utilized CEA and CoolProp for
    combustion properties and coolant properties respectively.

    Program Methodology
    - define engine with user input
    - run CEA and get bartz film coefficient axially along chamber
    - begin iterating along channel length from nozzle end to the injector face
        - integrate heat entry and pressure loss for the given step size
        - repeat process all the way up the channel
    - repeat calculation at the throat using interpolated values from the
    previous step as the initial conditions
    - display heat transfer, temperatures, pressure drop, film coefficient,
    channel geometry, coolant velocity on the engine contour
    
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

%% INITIALIZATION
clear;
clc; 
close all;
u = convertUnits;
CEA_input_name = 'AAAAAA';

%% SIMULATION PARAMETERS (INPUTS)
plots = 0; %Do ansys or not ???? DUmb name
steps = 100; % Number of steps along chamber (Change resolution of simulation)
qdot_tolerance = 0.0001; % set heattransfer convergence tolerance


%% ENGINE DEFINITION (INPUTS)

% Engine Contour
contour = readmatrix('contour_100pts.xlsx'); % import engine contour
r_contour = (contour(:,2) * u.IN2M)'; % contour radius [m]
x_contour = (contour(:,1) * u.IN2M)'; % contour x-axis [m]
[R_t, t_local] = min(r_contour); % throat radius, throat location [m]
chamber_length = 0.0254 * 5.205; % chamber length (m) [conversion * in]
converging_length = 0.0254 * 1.8251; %  converging length (m) [conversion * in]
diverging_length = 0.0254 * 1.8557; % diverging length (m) [conversion * in]
total_length = chamber_length + converging_length + diverging_length; % total length (mm) 

% Propulsion Parameters
P_c = 250; % chamber pressure [psi] 
P_e = 17; % exit pressure [psi]
m_dot = 10 * u.LB2KG; % Coolant/fuel mass flow [kg/s], 1.2566
fuel = 'C3H8O,2propanol'; % fuel definition
oxidizer = 'O2(L)'; % oxidizer definition
fuel_weight = 0; % ???  
fuel_temp = 273.15; % [K]
oxidizer_temp = 90.17; % [K]
OF = 1.2; % oxidizer/fuel ratio

% material properties
k_w = 200; % thermal conductivity [W/m-K]
E = 76E9; % [Pa] 
CTE = 22.4E-6; % in 1/K 
nu = 0.3; % poissons ratio (guess)
%e = 24 * 0.001; % surface roughness (mm) [micrometer*conversion]
roughness_table = readmatrix(pwd + "/bin/surface_roughness.xlsx",'Range','A12:E16');
e = [roughness_table(2,2), roughness_table(5,2)] .* 0.001; %Surface roughness (mm) [micrometer*conversion] [45, 90]

% Cooling channel inlet initialization
coolant = 'Water'; %coolant definition
inlet_temperature = 293.16; % inlet temperature [K]
%inlet_temperature = 350; % inlet temperature [K]
inlet_pressure = 500 * u.PSI2PA; % inlet pressure [PA]
coolantdirection = 0; % 1: direction opposite of hot gas flow direction
                      % 0: direction same as hot flow gas direction
                        
% channel geometry: (1: chamber) (min: throat) (2: nozzle end)
%t_w = 0.0005; % inner wall thickness [m]
t_w = [.001 .0005];
inter_length = .010; % Length where wall thickness will interpolate between chamber and nozzle
h_c = [.0015 .0014 .0015]; % channel height [1 min 2] [m]    
w_c = [.0054 .0016 .003227];% channel width [1 min 2] [m]
num_channels = 46; 
t_w_c = .001778 ; % channel width at torch igniter
t_h_c = .003 ; % channel height at torch igniter
torch_loc = [2.2 3.2 3.7] .* 0.0254; % location of torch changing area [mm] [inch * conversion] [1 min 2]


%% Parse variables && initial calculations

% Convert imperial units to metric
A_t = (R_t ^ 2) * pi; % throat area [m^2]
P_c = P_c * u.PSI2PA; % chamber pressure [Pa]
P_e = P_e * u.PSI2PA; % exit pressure [Pa]

% Chamber Geometry 
D_t = R_t * 2; % diameter at nozzle throat [m]
R_of_curve = 1.5 * D_t / 2; % [m]
A_local = pi * (R_t) ^ 2; % local cross sectional areas of engine
        
% Discretize Chamber Length 
deltax = (total_length/steps); % change in distance per step [m]
points = steps + 1; % number of points along chamber
x = 0:deltax:total_length; % length vector
x_plot = (x - chamber_length - converging_length); % length vector adjusted so that 0 is at the throat (mm)

% Parse engine section length vectors
torch_conv_length = torch_loc(2) - torch_loc(1);
torch_div_length = torch_loc(3)-torch_loc(2);

x_chamber1 = []; % chamber length vector before igniter
x_torch_conv = []; % igniter convergence
x_torch_div = []; % igniter divergence
x_chamber2 = []; % chamber length vector after igniter
x_converging = [];% converging length vector
x_diverging = [];% diverging length vectir
for i = x 
    if i <= torch_loc(1)
        x_chamber1 = [x_chamber1 i];
    end
    if (torch_loc(1) < i) && (i <= torch_loc(2))
        x_torch_conv = [x_torch_conv i];
    end
    if (torch_loc(2) < i) && (i <= torch_loc(3))
        x_torch_div = [x_torch_div i];
    end
    if (torch_loc(3) < i) && (i <= chamber_length)
        x_chamber2 = [x_chamber2 i];
    end 
    if (chamber_length < i) && (i <= chamber_length + converging_length)
        x_converging = [x_converging i];
    end 
    if i > (chamber_length + converging_length)
        x_diverging = [x_diverging i];
    end 
end 

% parse channel geometry [1 min 2]
A = w_c .* h_c; % channel cross-sectional area (m^2) [1 min 2]
p_wet = 2*w_c + 2*h_c; % wetted perimeter of the pipe (m) [1 min 2]
hydraulic_D = (4.*A)./p_wet; % hydraulic diameter (m) [1 min 2]

% parse channel geometry over channel length
w_c_chamber1 = ones(1,size(x_chamber1,2)).*w_c(1); % channel width over chamber length (constant)
w_c_torch_conv = ((t_w_c-w_c(1))/(torch_conv_length)).*(x_torch_conv -x_torch_conv(1)) ... 
    + ones(1,size(x_torch_conv,2)).*w_c(1); % channel width over converging torch section
w_c_torch_div = ((w_c(1)-t_w_c)/(torch_div_length)).*(x_torch_div-x_torch_div(1))... 
    + ones(1,size(x_torch_div,2)).*t_w_c; % channel width over diverging torch section
w_c_chamber2 = ones(1,size(x_chamber2,2)).*w_c(1); % channel width over chamber length (constant)
w_c_converging = ((w_c(2)-w_c(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
         + ones(1,size(x_converging,2)).*w_c(1); % channel width over converging length (linear interpolation)
w_c_diverging = ((w_c(3)-w_c(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
        + ones(1,size(x_diverging,2)).*w_c(2);   % channel width over diverging length (linear interpolation)
w_c_x = [w_c_chamber1 w_c_torch_conv w_c_torch_div w_c_chamber2 w_c_converging w_c_diverging]; % combine channel width vectors

h_c_chamber1 = ones(1,size(x_chamber1,2)).*h_c(1); % channel height over chamber length (constant)
h_c_torch_conv = ((t_h_c-h_c(1))/(torch_conv_length)).*(x_torch_conv -x_torch_conv(1)) ... 
    + ones(1,size(x_torch_conv,2)).*h_c(1); % channel height over converging torch section
h_c_torch_div = ((h_c(1)-t_h_c)/(torch_div_length)).*(x_torch_div-x_torch_div(1))... 
    + ones(1,size(x_torch_div,2)).*t_h_c; % channel height over diverging torch section
h_c_chamber2 = ones(1,size(x_chamber2,2)).*h_c(1); % channel height over chamber length (constant
h_c_converging = ((h_c(2)-h_c(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
     + ones(1,size(x_converging,2)).*h_c(1);    % channel height over converging length (linear interpolation)
h_c_diverging = ((h_c(3)-h_c(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
     + ones(1,size(x_diverging,2)).*h_c(2); % channel height over diverging length (linear interpolation)
h_c_x = [h_c_chamber1 h_c_torch_conv h_c_torch_div h_c_chamber2 h_c_converging h_c_diverging]; % combine channel height vectors

A_x = (w_c_x .* h_c_x); % channel area vector over channel length [m^2]
p_wet_x = 2.*w_c_x + 2 .* h_c_x; % wet perimeter over channel length [m]
hydraulic_D_x = ((4.*(A_x))./p_wet_x); % bydraulic Diameter over channel length [m]

x_to_chamber2 = [x_chamber1 x_torch_conv x_torch_div];
x_to_converging = [x_to_chamber2 x_chamber2];
x_to_throat = [x_to_converging x_converging];
% calculate channel flow
m_dot_CHANNEL = m_dot / num_channels; % mass flow of channel (EQ 6.31)

% Wall thickness vector
x_inter_wall = [];
x_nozzle_wall = [];
for i = x
    if ((chamber_length) < i) && (i <= chamber_length + inter_length)
        x_inter_wall = [x_inter_wall i];
    end
    if(chamber_length + inter_length < i)
        x_nozzle_wall = [x_nozzle_wall i];
    end
end

t_w_chamber = t_w(1) * ones(size(x_to_converging));
t_w_inter = ((t_w(2)-t_w(1))/(inter_length)).*(x_inter_wall -x_inter_wall(1)) ... 
    + ones(size(inter_length)).*t_w(1); % channel width over converging torch section
t_w_nozzle = t_w(2) .* ones(size(x_nozzle_wall));
t_w_x = [t_w_chamber t_w_inter t_w_nozzle];

%% CHAMBER HEAT TRANSFER CALCULATIONS

% Step 1: Prescribe initial properties

% Prescribe area ratios
r_interpolated = interp1(x_contour,r_contour,x_plot,'linear','extrap'); % linearly interpolate radius vector  
subsonic_area_ratios = (pi * r_interpolated(x_plot < 0) .^ 2) / A_t; % subsonic area ratios on discretized points
supersonic_area_ratios = (pi * r_interpolated(x_plot > 0) .^ 2) / A_t; %  supersonic area ratios on discretized points
A_ratio = [subsonic_area_ratios, supersonic_area_ratios]; % area ratio vector [sub, sup]

% initialize property matrices
% axial coolant property matrices
P_l = zeros(1, points); % coolant pressure
T_l = zeros(1, points); % coolant temp
rho_l = zeros(1, points); % coolant density
v = zeros(1, points); % coolant ???

% axial cooling property matrices
qdot_l = zeros(1, points);  % liquid convective heat flux
qdot_g = zeros(1, points);  % gas convective heat flux
T_wl = zeros(1, points);    % liquid wall temperature
T_wg = zeros(1, points);    % gas wall temperature
h_g = zeros(1, points);     % gas film coefficient
sigma = zeros(1, points);   % film coefficient correction factor
h_l = zeros(1, points);     % liquid film coefficient

% axial combustion property matrices
c_star = zeros(1, points);  % characteristic velocity
M = zeros(1, points);       % mach number
gamma = zeros(1, points);   % ???
P_g = zeros(1, points);     % combustion pressure
T_g = zeros(1, points);     % combustion temperature
mu_g = zeros(1, points);    % combustion viscosity 
Pr_g = zeros(1, points);    % combustion prantl number
cp_g = zeros(1, points);    % combustion coefficient of pressure ???

% stress matricies
sigma_t = zeros(1,points); % tangential stress
sigma_tp = zeros(1,points); % tangential stress pressure
sigma_tt = zeros(1,points); % tangential stress temp
sigma_l = zeros(1,points); % longitudinal stress
sigmab = zeros(1,points); % buckling stress
sigma_v = zeros(1,points); % von mises stress

% call cea for all area ratios
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
        sigma(i) = (.5 * T_wg(i) / T_g(1) * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) + .5) ^ -.68 * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
        h_g(i) = (0.026 / D_t ^ 0.2) * (mu_g(i) ^ 0.2 * cp_g(i) / Pr_g(i) ^ 0.6) * (P_c / c_star(i)) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / A_ratio(i)) ^ .9 * sigma(i); % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
        r = Pr_g(i) ^ (1 / 3); % recovery factor for a turbulent free boundary layer [N/A] - biased towards larger engines, very small engines should use Pr^.5 (Heister Table 6.2).
        T_r = T_g(i) * (1 + (gamma(i) - 1) / 2 * r * M(i) ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 
        qdot_g(i) = h_g(i) * (T_r - T_wg(i)); % gas convective heat flux [W/m^2] (Heister EQ 6.16).
    
        % Step 6: Calculate liquid wall temperature
        T_wl(i) = T_wg(i) - qdot_g(i) * t_w_x(i) / k_w; % liquid wall temperature calculated via conduction through wall [K] (Heister EQ 6.29).
    
        % Step 7: Calculate liquid film coefficient
        % run coolprop to get coolant properties
        mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l(i), 'P', P_l(i), coolant); % viscosity of bulk coolant [Pa-s]
        cp_l = py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l(i), 'P', P_l(i), coolant); % specific heat of coolant [J/kg-k] 
        k_l = py.CoolProp.CoolProp.PropsSI('L', 'T', T_l(i), 'P', P_l(i), coolant); % thermal conductivity of coolant [W/m-K]
        rho_l(i) = py.CoolProp.CoolProp.PropsSI('D','T', T_l(i),'P', P_l(i),'Water'); % density of the coolant [???]
        v(i) = m_dot_CHANNEL / rho_l(i) / A_x(i); % velocity at step [m/s]
       
        Re_l = (rho_l(i) * v(i) * hydraulic_D_x(i)) / mu_lb; % reynolds number for channel flow [N/A] (Huzel and Huang , pg 90)
        Pr_l = (cp_l * mu_lb) / k_l; % prantl number [N/A] (Huzel and Huang, pg 90) 
        Nu_l = 0.023 * (Re_l ^ .8) * (Pr_l ^ .4) * (T_wl / T_l) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19). ****
        h_l(i) = (Nu_l * k_l) / hydraulic_D_x(i); % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)
    
        % Step 8: Calculate liquid-side convective heat flux
        qdot_l(i) = h_l(i) * (T_wl(i) - T_l(i)); % liquid convective heat flux [W/m^2] (Heister EQ 6.29).

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
                T_l(i+1) = T_l(i) + (1 / (m_dot_CHANNEL * cp_l)) * qdot_g(i) * wall_area; % new liquid temperature [K] (Heister EQ 6.39)
                
                % Use moody diagram to find coefficient of friction
                if (i < size(x_chamber1,2)) || (((size(x_to_chamber2,2)) <= i) && (i < size(x_to_converging,2)))
                    ed = e(2)/(hydraulic_D_x(i)*1000); %90 degrees
              
                else
                    ed = e(1)/(hydraulic_D_x(i)*1000); %45 degrees
                    
                end
                %ed = e/(hydraulic_D_x(i)*1000); % relative roughness
                f = moody(ed, Re_l); % friction factor
                cf = f/4; % friction coefficient

                if i > 1
                    deltaP =   (2*cf*(deltax/(hydraulic_D_x(i))) * rho_l(i) *(v(i))^(2)  + .5 * ((v(i)^2) -(v(i-1)^2))); % change in pressure (Bernoulli's equation)
                else
                    deltaP =   (2*cf*(deltax/(hydraulic_D_x(i))) * rho_l(i) *(v(i))^(2)); % change in pressure (Heister 6.36)
                end
                
                % calculate stesses 
                sigma_tp(i) = ( ((P_l(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) );
                sigma_tt(i) = (E*CTE*qdot_g(i)*t_w_x(i))/(2*(1-nu)*k_w);
                sigma_t(i) = ( ((P_l(i)-P_g(i))/2).*((w_c_x(i)./t_w_x(i)).^2) ) + (E*CTE*qdot_g(i)*t_w_x(i))/(2*(1-nu)*k_w); % tangential stress
                %sigma_l(i) = E*CTE*(T_wg(i)-T_wl(i)); % longitudinal stress (The temperatures here are wrong and I'm not sure this is applicable to rectagular channels
                sigma_l(i) = E*CTE*(T_wl(i)-300);
                %sigmab = ??? ; % buckling stress
                sigma_v(i) = sqrt(sigma_l(i)^2 + sigma_t(i)^2 - sigma_l(i)*sigma_t(i));
                % prepare for next step
                P_l(i+1) =   P_l(i) - deltaP; % Update pressure for next iteration
                T_wg(i+1) = T_wg(i);  % new gas wall temp guess based on current temp
            end  
            converged = 1;
        end
    end
end

%% THROAT HEAT TRANSFER CALCULATION
% determines heat transfer and temperature at the throat

% Step 1: Prescribe initial properties
[c_star_t, ~, ~, M_t, gamma_t, P_g_t, T_g_t, ~, mu_g_t, Pr_g_t, ~, ~, ~, cp_g_t] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 2, 0, 0, CEA_input_name); % ******

% Steps 2 & 3: Set channel inlet properties
P_l_t = (P_l(length(x_to_throat)) + P_l(length(x_to_throat) + 1) ) / 2; % coolant pressure at throat [Pa] (interpolated)
T_l_t = (T_l(length(x_to_throat)) + T_l(length(x_to_throat) + 1) ) / 2; % coolant tepemperature at throat [K] (interpolated)

% Step 4: Take hot wall temperature guess and initialize loop

T_wg_t = (T_wg(length(x_to_throat)) + T_wg(length(x_to_throat) + 1) ) / 2; % initial guess of wall side temperature at throat [K] (interpolated) 
T_wg_mn = 294.15; % minimum temperature bound
T_wg_mx = 2000; % maximum temperature bound

converged = 0; % wall temperature loop end condition
counter = 0; % counter for loop


% start loop to converge on wall temperature
while ~(converged)
    % Step 5: Calculate gas film coefficient and gas-side convective heat flux
    sigma = (.5 * T_wg_t / T_g(1) * (1 + (gamma_t - 1) / 2 * M_t ^ 2) + .5) ^ -.68 * (1 + (gamma_t - 1) / 2 * M_t ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
    h_g_t = (0.026 / D_t ^ 0.2) * (mu_g_t ^ 0.2 * cp_g_t / Pr_g_t ^ 0.6) * (P_c / c_star_t) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (A_t / A_t) ^ .9 * sigma; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
    r = Pr_g_t ^ (1 / 3); % recovery factor for a turbulent free boundary layer [N/A] - biased towards larger engines, very small engines should use Pr^.5 (Heister Table 6.2).
    T_r = T_g_t * (1 + (gamma_t - 1) / 2 * r * M_t ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 
    qdot_g_l = h_g_t * (T_r - T_wg_t); % gas convective heat flux [W/m^2] (Heister EQ 6.16).

    % Step 6: Calculate liquid wall temperature
    T_wl_t = T_wg_t - qdot_g_l * t_w(2) / k_w; % liquid wall temperature calculated via conduction through wall [K] (Heister EQ 6.29).

    % Step 7: Calculate liquid film coefficient
    % run coolprop to get coolant properties
    mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l_t, 'P', P_l_t, coolant); % viscosity of bulk coolant [Pa-s]
    cp_l = py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l_t, 'P', P_l_t, coolant); % specific heat of coolant [J/kg-k] 
    k_l = py.CoolProp.CoolProp.PropsSI('L', 'T', T_l_t, 'P', P_l_t, coolant); % thermal conductivity of coolant [W/m-K]
    rho_l_t = py.CoolProp.CoolProp.PropsSI('D','T', T_l_t,'P', P_l_t,'Water'); % density of the coolant [???]
    v_t = m_dot_CHANNEL / rho_l_t / A(2); % velocity at step [m/s]  

    Re_l = (rho_l_t * v_t * hydraulic_D(2)) / mu_lb; % reynolds number for channel flow [N/A] (Huzel and Huang , pg 90)
    Pr_l = (cp_l * mu_lb) / k_l; % prantl number [N/A] (Huzel and Huang, pg 90)
    Nu_l = 0.023 * (Re_l ^ .8) * (Pr_l ^ .4) * (T_wl_t / T_l_t) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19).
    h_l_t = (Nu_l * k_l) / hydraulic_D(2); % liquid film coefficient [W/m^2-K] (Heister EQ 6.19)
    
    % Step 8.5: Fin heat transfer, adiabatic tip, Biot << 0.1
    T_base = T_wl_t;                      % Temperature at fin base
    P_fin = 2 * deltax + 2 * w_c(2);    % Fin perimeter (step distance & channel width)
    A_c_fin = w_c(2) * deltax;          % Fin area at current step
    m_fin = sqrt(h_l_t * P_fin / (k_w * A_c_fin));                % Fin m
    M_fin = sqrt(h_l_t * P_fin * k_w * A_c_fin) * (T_base - T_l_t); % Fin M
    fin_q = M_fin * tanh(m_fin * h_c(2));               % Fin heat transfer rate
    eff_fin = tanh(m_fin * h_c(2)) / (M_fin * h_c(2));  % Fin efficienct
    biot_fin = h_l_t / k_w * (0.5 * (w_c(2) + h_c(2)));                      % Biot number to test assumptions

    % Step 8: Calculate liquid-side convective heat flux
    qdot_l = h_l_t * (T_wl_t - T_l_t) + fin_q; % liquid convective heat flux [W/m^2] (Heister EQ 6.29).

    % Step 9: Check for convergence and continue loop / next step
    if abs(qdot_g_l - qdot_l) > qdot_tolerance % check for tolerance

        % convergence loop
        if qdot_g_l - qdot_l > 0
            T_wg_mn = T_wg_t;
        else 
            T_wg_mx = T_wg_t;
        end 
        T_wg_t = (T_wg_mx + T_wg_mn) / 2;

        counter = counter + 1;
    else 
        fprintf("Gas Side Wall Temp [K]: %0.2f\n", T_wg_t)

        converged = 1; % end loop
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

figure('Name', 'Heat Transfer Plots');
subplot(2,2,[1,2])
hold on;
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

legend('Convective Heat Flux', 'Chamber Contour')
title('Heat Flux Distribution')
xlabel('Location [mm]')

subplot(2,2,3)
plot(x_plot.*1000, h_g)
title("Gas Film Coeffcient [W/m^2-K]")
xlabel("Location [mm]");
grid on
subplot(2,2,4)
plot(x_plot.*1000, h_l)
title("Liquid Film Coefficient [W/m^2-K]")
grid on



figure('Name','Water Flow')
subplot(2,2,[1,2])
plot(x_plot.* 1000, P_l * 1/6894.757)
title("Liquid Pressure Loss")
xlabel("Location [mm]")
ylabel("Pressure [PSI]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend("Pressure Curve"," Chamber Contour")
grid on

subplot(2,2,3)
plot(x_plot.* 1000, v);
title("Coolant Velocity [m/s]")
xlabel("Location [mm]")
grid on
subplot(2,2,4)
plot(x_plot.*1000, T_l)
title("Coolant Temperature [K]");
grid on




figure('Name','Channel Geometry');
subplot(2,2,[1,2]);
hold on
plot(x_plot.* 1000, w_c_x .*1000);
plot(x_plot.*1000, h_c_x .*1000);
title("Channel Dimensions");
xlabel("Location [mm]");
ylabel("Channel Dimensions [mm]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Chamber Contour [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Channel Width', 'Channel Height', 'Chamber Contour')
grid on

subplot(2,2,3);
plot(x_plot.*1000, t_w_x .* 1000);
title("Wall thickness [mm]");
xlabel("Location [mm]");
grid on
subplot(2,2,4);
plot(x_plot.* 1000, A_x .* 1000000);
title("Channel Area [mm^2]");
xlabel("Location [mm]");
grid on



    
figure('Name', 'Structural results')
subplot(2,2,[1,2])
plot(x_plot.* 1000, sigma_v*0.000001,'g')
title("Von Mises Stress")
xlabel("Location [mm]")
ylabel("[MPA]")
yyaxis right
plot(x_plot .* 1000, r_interpolated .* 1000, 'black', 'LineStyle', '-');
ylabel('Radius [mm]')
set(gca, 'Ycolor', 'k')
axis equal;
legend('Von Mises Stress','Chamber Contour')
grid on
subplot(2,2,3)
hold on
plot(x_plot.* 1000, sigma_t*0.000001,"b");
plot(x_plot.* 1000, sigma_tp*0.000001,"m");
plot(x_plot.* 1000, sigma_tt*0.000001,"r");
legend("Total Stress", "Pressure contribution", "Thermal Contribution")
title("Tangential Stress (MPA)")
hold off
grid on
subplot(2,2,4)
plot(x_plot.* 1000, sigma_l*0.000001);
title("Longitudinal Stress (MPA)")
xlabel("Location [mm]")
grid on
% subplot(2,2,3)
% plot(x_plot.* 1000, sigmab*0.000001)
% title("Buckling Stress (MPA)")
% xlabel("Location [mm]")






%         %Step 12: Structural Analysis Checks
%         St = .5*(Pl- P_gas)*(width/wallthick)^2 + (E*a*gasheattransfer*wallthick)/(2*(1-v)kw); %Combined tangential stresses: Heister Eq 6.33 page 207
%         Sl = E*a*(Twl-fuel_temp); %Longtudinal thermal stress (Huzel and Huang, EQ 2-28, pg 92) The temperatures used here may not right for determining the delta T in this equation.
%         Sc = 4*Et*Ec*t_w/((((Et)^(1/2))*((Ec)^(1/2))^2)*(3*(1-v^2)*tube_radius)); %Critical Stress Buckling (Huzel and Huang, Eq 4-29, pg 93)

%% THERMAL FEA
L_seg = 0.0283;
length = L_seg * u.IN2M;   
if plots

    M = 150;
    N = 150;
    R1 = R_t; % inner radius 
    R2 = R_t + h_c(2) * 4;  % outer radius
    nR = linspace(R1,R2,M);
    nT = linspace(-pi/num_channels, pi/num_channels + w_c(2) / R_t, N);
    [R, T_g] = meshgrid(nR,nT) ;
    xg = R.*cos(T_g); 
    yg = R.*sin(T_g);
    xg = xg(:);
    yg = yg(:);
    
    % Define partial channel 
    M = 50;
    N = 50;
    R1 = R_t + t_w; % inner radius 
    R2 = R_t + t_w + h_c(2);  % outer radius
    x = linspace(R1,R2,M);
    y = linspace(-pi/num_channels - w_c(2) / R_t / 2, -pi/num_channels + w_c(2) / R_t / 2, N);
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
    R2 = R_t + t_w + h_c(2);  % outer radius
    x = linspace(R1,R2,M);
    y = linspace(pi/num_channels - w_c(2) / R_t / 2, pi/num_channels + w_c(2) / R_t / 2, N);
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
    
    generateMesh(model,"Hmax",h_c(2)/12);
    % figure
    % pdemesh(model)
    
    % Define material thermal properties
    thermalProperties(model,"ThermalConductivity",k_w);
    
    % Thermal boundary conditions
    thermalBC(model,"Face",10, ...
                     "ConvectionCoefficient",h_g_t, ...
                     "AmbientTemperature",T_r);
    thermalBC(model,"Face",[5 13 7 3 11 12 14], ...
                     "ConvectionCoefficient",h_l_t, ...
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

