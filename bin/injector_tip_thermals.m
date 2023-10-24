clear;
clc;
u = convertUnits;

F_max = 550;
R_t = 0.7335 * u.IN2M;
A_t = pi * R_t ^ 2;
P_c = 250;
P_e = 17;
fuel = 'C3H8O,2propanol';
fuel_weight = 0;
fuel_temp = 293.15;
oxidizer = 'O2(L)';
oxidizer_temp = 90.17;
OF = 1.2;
CEA_input_name = 'test';

% Engine geometry definition
D_t = R_t * 2; 
R_of_curve = 1.5 * D_t / 2; % [m]
R_c = 1.875 * u.IN2M;
A_c = pi * R_c ^ 2; 

% Reference diameters
hydraulic_D = 0.5 * u.IN2M;   % LOx 

% Find mass flow per channel
channelMdot = 1.4837 * u.LB2KG; % ox

% Injector definition  
t_w = 0.001;  % [m]
dP_Pc = 0.25;
P_l = P_c / (1 - dP_Pc) * u.PSI2PA;    % [Pa]
T_l = 90;     % [K]

% Iteration initialization
T_wg = 1000; % initial guess of wall side temperature [K]
steps = 1;
distance = 0.01; 
delta_x = distance/steps; 

% wall material properties
k_w = 396; % Thermal Conductivity of Wall [W/m-K]

% CEA
P_c = P_c * u.PSI2PA;
[c_star, ~, ~, M, gamma, P_g, Tc_ns, ~, mu_g, Pr_g, ~, ~, ~, cp_g] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 1, 0, 0, CEA_input_name); 
% [c_star, ~, ~, M, gamma, P_g, T_g, ~, mu_g, Pr_g, ~, ~, ~, cp_g] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 2, 0, 0, CEA_input_name); 

qdot_tolerance = 0.0001;

for i = 1:steps
    converged = 0;
    T_wg_mn = 50; % minimum temperature bound
    T_wg_mx = 3000; % maximum temperature bound

    while ~(converged)
        % Step 5: Calculate gas film coefficient and gas-side convective heat flux
        sigma = (.5 * T_wg / Tc_ns * (1 + (gamma - 1) / 2 * M ^ 2) + .5) ^ -.68 * (1 + (gamma - 1) / 2 * M ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
        h_g = (0.026 / D_t ^ 0.2) * (mu_g ^ 0.2 * cp_g / Pr_g ^ 0.6) * (P_c / c_star) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (A_t / A_c) ^ .9 * sigma; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
        r = Pr_g ^ (1 / 3); % recovery factor for a turbulent free boundary layer [N/A] - biased towards larger engines, very small engines should use Pr^.5 (Heister Table 6.2).
        T_r = Tc_ns * (1 + (gamma - 1) / 2 * r * M ^ 2); % recovery temperature [K] - corrects for compressible boundry layers (Heister EQ 6.15). 
        qdot_g = h_g * (T_r - T_wg); % gas convective heat flux [W/m^2] (Heister EQ 6.16).
        
        % Step 6: Calculate liquid wall temperature
        T_wl = T_wg - qdot_g * t_w / k_w; % liquid wall temperature calculated via conduction through wall [K] (Heister EQ 6.29).
        
        % Step 7: Calculate liquid film coefficient
        % run coolprop to get coolant properties
        mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l, 'P', P_l, "O2"); % viscosity of bulk coolant [Pa-s]
        cp_l = py.CoolProp.CoolProp.PropsSI('C' , 'T', T_l, 'P', P_l, "O2"); % specific heat of coolant [J/kg-k] 
        k_l = py.CoolProp.CoolProp.PropsSI('L', 'T', T_l, 'P', P_l, "O2"); % thermal conductivity of coolant [W/m-K]
        
        Re_l = (4 * channelMdot) / (pi * hydraulic_D * mu_lb); % reynolds number for channel flow [] ALEX CITE SOURCE
        Pr_l = (cp_l * mu_lb) / k_l; % prantl number [] ALEX CITE SOURCE
        Nu_l = 0.023 * (Re_l ^ .8) * (Pr_l ^ .4) * (T_wl / T_l) ^ -.3; % nusselt number [N/A] - applicable for Re > 10,000, .7 < Pr < 160 (Heister EQ 6.19).
        h_l = (Nu_l * k_l) / hydraulic_D; % liquid film coefficient [W/m^2-K] ALEX CITE SOURCE
        
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
    
        else 
            wall_area = pi * hydraulic_D * delta_x; 
            T_l = T_l + (1 / (channelMdot * cp_l)) * qdot_g * wall_area; % new liquid temperature [K] (Heister EQ 6.39).

            rho_l = py.CoolProp.CoolProp.PropsSI('D','T', T_l,'P', P_l, 'O2');
            v = channelMdot / rho_l / (hydraulic_D ^ 2 / 4 * pi); % velocity at step [m/s]

            converged = 1; % end loop
        end
    end
end

%% Outputs
fprintf("Pintle Tip\n")
fprintf("Gas Side Wall Temp [K]: %0.2f\n", T_wg)
fprintf("Gas Heat Transfer Coefficient: %0.2f\n", h_g)
fprintf("Liquid Heat Transfer Coefficient: %0.2f\n", h_l)
fprintf("Gas Temp: %0.2f\n", Tc_ns)
fprintf("Ox Temp: %0.2f\n", T_l)



