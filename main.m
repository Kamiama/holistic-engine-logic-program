%% Holistic Engine Logic Program (HELP)
% Author: Kamon Blong (kamon.blong@gmail.com)
% First Created: 7/10/2022
% Last Updated: 
%{ 

Description: 

Assumptions: Currently assumes a stationary test engine with no regard for
changes in gravity or atmospheric pressure. Additionally assumes an
ideal rocket engine (perfect propellant mixing, gaseous flow, perfect gas
law, adiabatic, neglects friction and boundry layer effects, no
discontinuities, steady mDot, no transience, all exhaust contributes to
thrust, characteristics are uniform at any point along an axial slice of
the engine, and chemical equilibrium is constant axially throughout the engine 
(Sutton 46)). Without inputing a specific inlet temperature for
propellants, NASA CEA (and therefore this program) also assumes that
ordinary propellants are stored at room temperature and cryogens at their
boiling point.
%}

clear;
<<<<<<< Updated upstream
clc;
=======
%clc;
%close all;
>>>>>>> Stashed changes

%% Import Data
% import input properties
input_data = 'regen.xlsx';

%% Identify OS & Generate Paths 
% identify OS
if ismac || isunix
    function_path = append(pwd, '/bin');
    cea_path = append(pwd, '/cea');
    input_path = append(pwd, '/input');
    output_path = append(pwd, '/output');
elseif ispc
    function_path = append(pwd, '\bin');
    cea_path = append(pwd, '\cea');
    input_path = append(pwd, '\input');
    output_path = append(pwd, '\output');
end

% assign paths
main_path = cd;
addpath(function_path);
addpath(cea_path);
addpath(input_path)
addpath(output_path);
u = convertUnits;

%% Initialization
% initialize output document
run_description = input('Run Description: ', 's'); % run description that will be added to output file header
write_date = datestr(now, '_mm-dd-yy_HH:MM'); % date that will be added to output file header
file_name = strcat(run_description, write_date); % output file name

%% Program Run Definition
% parse input run properties
    sizeEngine = 1;
    sizeStructures = 1;
    enableREFPROP = 1;
    enableFigures = 1;
    enableDebug = 1;
    if enableREFPROP
        sizeFluids = 1;
        sizeInjector = 1;
        sizeCooling = 1;
    end

<<<<<<< Updated upstream
=======
    resolution_contour = 50;
    resolution_throttle = 7;

>>>>>>> Stashed changes
%% Parse Variables & Define Constants
% parse & convert input properties

    % performance properties
        F = xlsread(input_data,'Engine', 'C8');                    % engine thrust [lbf]
        P_c = xlsread(input_data, 'Engine', 'C9');                 % chamber pressure [psi]
        P_a = xlsread(input_data, 'Engine', 'C11');                % ambient pressure [psi]
        P_e = xlsread(input_data, 'Engine', 'C10');                % exit pressure [psi] (equals P_a for optimal expansion)
        eff_c_star = xlsread(input_data, 'Engine', 'C12') / 100;   % c* efficiency
        eff_c_f = xlsread(input_data, 'Engine', 'C13') / 100;      % cF efficiency
    
    % propellant properties
        [~, txt] = xlsread(input_data, 'Engine', 'H10');
        fuel1 = string(txt);                                       
        [~, txt] = xlsread(input_data, 'Engine', 'H12');
        fuel2 = string(txt);                                       
        fuel = [fuel1, fuel2];                                     % fuel formula for NASA CEA
        fuel_weight = xlsread(input_data,'Engine', 'H13');         % fuel weights
        [~, txt] = xlsread(input_data, 'Engine', 'H11');
        oxidizer = string(txt);                                    % oxidizer formula for NASA CEA
        fuel_temp = xlsread(input_data,'Engine', 'H17');           % inlet fuel temperature [K]
        oxidizer_temp = xlsread(input_data,'Engine', 'H18');       % inlet oxidizer temperature [K]
        OF = xlsread(input_data, 'Engine', 'H14');                 % O/F ratio
        time_burn = xlsread(input_data, 'Engine', 'C14');          % burn time [sec]

    % correlate NASA CEA propellant names with REFPROP names and assign temperature values
        if fuel(1) == "Jet-A(L)" || fuel(1) == "RP-1"
            fuel_REFPROP = 'RP1.mix';
            fuel_temp = 293.15;
        elseif fuel(1) == "C2H5OH(L)"
            fuel_REFPROP = 'ethanol';
            fuel_temp = 293.15;
        elseif fuel(1) == "CH4(L)"
            fuel_REFPROP = 'methane';
            fuel_temp = 111.6;
        elseif fuel(1) == "CH4"
            fuel_REFPROP = 'methane';
            fuel_temp = 293.15;
        elseif fuel(1) == "H2(L)"
            fuel_REFPROP = 'hydrogen';
            fuel_temp = 20.38;
        elseif fuel(1) == "H2"
            fuel_REFPROP = 'hydrogen';
            fuel_temp = 293.15;
        elseif fuel(1) == "C3H8O,1propanol"
            fuel_REFPROP = 'hydrogen';
            fuel_temp = 293.15;
        else
            if enableDebug
                fprintf("\nWarning: unrecognized propellant")
            end
        end
    
        if oxidizer(1) == "O2(L)"
            oxidizer_REFPROP = 'oxygen';
            oxidizer_temp = 90.17;
        elseif oxidizer(1) == "O2"
            oxidizer_REFPROP = 'oxygen';
            oxidizer_temp = 293.15;
        elseif oxidizer(1) == "N2O"
            oxidizer_REFPROP = 'oxygen';
            oxidizer_temp = 184.7;
        elseif oxidizer(1) == "H2O2(L)"
            oxidizer_REFPROP = 'h2o2';
            oxidizer_temp = 293.15;
        else
            if enableDebug
                fprintf("\nWarning: unrecognized propellant")
            end
        end
    
    % geometry properties

        % definition
        [~, txt] = xlsread(input_data, 'Engine', 'M10');      % geometry type
        geometry_type = string(txt);
        [~, txt] = xlsread(input_data, 'Engine', 'M11');      % sizing method
        sizing_method = string(txt);       
        [~, txt] = xlsread(input_data, 'Engine', 'M12');      % converging method
        converging_method =string(txt);   

        % general
        L_star = xlsread(input_data, 'Engine', 'M15');              % L*, characteristic combustion length [in]
        conv_angle = xlsread(input_data, 'Engine', 'M16');          % convergence 
        con_ratio = xlsread(input_data, 'Engine', 'M17');           % contraction ratio
        D_c = xlsread(input_data, 'Engine', 'M18');                 % chamber diameter [in]
        D_t = xlsread(input_data, 'Engine', 'M19');                 % throat diameter [in]

        % bell
        bell_pct = xlsread(input_data, 'Engine', 'M22') / 100;      % percent of bell nozzle

        % conical
        conical_half_angle = xlsread(input_data, 'Engine', 'M25');  % conical half angle [deg]
        L_crude_throat = xlsread(input_data, 'Engine', 'M26');      % length of straight throat section [in]
        bar_size = xlsread(input_data, 'Engine', 'M27'); % bar stock diameter [in]
            
    % injector properties

        % general parameters
        [~, txt] = xlsread(input_data, 'Injector', 'C8');
        injector_type = string(txt);                                % injector type
        d_P = 0; % pressure drop as percentage of chamber pressure [psi]
        spray_angle = 0; % spray angle [deg]

        % impinging jets
        if injector_type == "impinging"
            
        % pintle
        elseif injector_type == "pintle"

        % coaxial swirl
        elseif injector_type == "swirl"

        % coaxial shear
        elseif injector_type == "shear"

        end

    
    % cooling properties
        cooling_type = "regen";  % cooling method (regen or ablative)
        film_pct = 0;            % film cooling percent as fraction of fuel mass flow rate
        nozzle_regen_pct = 1;    % percent of nozzle to be regeneratively cooled

    % fluid properties
        vol_tank = 0;   % volume of propellant tanks [in^3]
        p_tank = 1000;  % tank pressure [psi]
    
    % define constants
        g = 32.174;     % gravitational constant [ft/sec^2]

%% Run NASA CEA
% sends P_c, P_e, and fuel characteristics from input.xlsx as inputs
% returns Isp, c*, and expansion ratio as outputs

% initialize CEA input/output files
CEA_name_prefix = strcat('CEA');
CEA_input_name = append(CEA_name_prefix, '.inp');
CEA_output_name = append(CEA_name_prefix, '.out');

% run CEA
[cea_c_star, cea_isp, exp_ratio, M, gamma, ~, T, ~, mu, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 0, 1);

%% Intermediate Calculations

% correct cea values & find effective exhaust velocity
isp = cea_isp * eff_c_star * eff_c_f;          % expected isp [sec]
c_star = cea_c_star * eff_c_star;              % expected c* [ft/s]
c = isp * g;                                   % effective exhaust velocity [ft/s]
    
% size m_dot normally
if sizing_method == "normal"
    m_dot = F / isp;                               % TOTAL mass flow rate [lbm/s]
    
% size m_dot based off of throat size
elseif sizing_method == "throat"
    A_t = (D_t / 2) ^ 2 * pi;                      % throat area [in^2]
    m_dot = A_t / c_star * P_c * g;                % TOTAL mass flow rate [lbm/s]
    F = m_dot * isp;                               % thrust [lbf]

% throw error for no selected sizing method
else
    error("No sizing method selected.")
end

% calculate mass flow rate of fuel, oxidizer, and film cooling
m_dot_FUEL = m_dot / (OF + 1);                 % FUEL mass flow rate [lbm/s]
m_dot_OX = m_dot_FUEL * OF;                    % OX mass flow rate [lbm/s]

m_dot_FILM = m_dot_FUEL * film_pct;            % FILM mass flow rate [lbm/s]
m_dot_INJ_FUEL = m_dot_FUEL * (1 - film_pct);  % regularly injected FUEL mass flow rate [lbm/s]

% output performance information
fprintf('\n-------- Performance Outputs --------\n')
disp(['                      Thrust (lbf): ' num2str(F)]);
disp(['            Chamber Pressure (psi): ' num2str(P_c)]);
disp(['               Exit Pressure (psi): ' num2str(P_e)]);
disp(['                              Fuel: ' fuel]);
disp(['                          Oxidizer: ' oxidizer]);
disp(['                         O/F Ratio: ' num2str(OF)]);
fprintf('\n')
disp(['                     CEA Isp (sec): ' num2str(cea_isp)]);
disp(['                Expected Isp (sec): ' num2str(isp)]);
disp([' Effective Exhaust Velocity (ft/s): ' num2str(c)]);
disp(['    Mass Flow Rate (total) (lbm/s): ' num2str(m_dot)]);
disp(['     Mass Flow Rate (fuel) (lbm/s): ' num2str(m_dot_FUEL)]);
disp([' Mass Flow Rate (oxidizer) (lbm/s): ' num2str(m_dot_OX)]);
disp(['         Adiabatic Flame Temps (R): ' num2str(T(1)) ' ' num2str(T(2))]);
disp(['                         Exit Mach: ' num2str(M(end))]);
disp(['                            Gammas: ' num2str(gamma(1)) ' ' num2str(gamma(2))]);
fprintf('-------- Performance Outputs --------\n')

%% Basic Geometry Calculations
A_t = m_dot * c_star / P_c / g;  % throat area [in^2]
A_e = A_t * exp_ratio;           % exit area [in^2]
D_t = 2 * sqrt(A_t / pi);        % throat diameter [in]

% size based off contraction ratio
if converging_method == "contraction"
    A_c = A_t * con_ratio;       % chamber area [in^2]

% size based off chamber diameter
elseif converging_method == "diameter"
    A_c = pi * (D_c / 2) ^ 2;    % chamber area [in^2]
    con_ratio = A_c / A_t;       % contraction ratio

% size based off a calculated contraction ratio
elseif converging_method == "auto"
    con_ratio = 8 * (2 * sqrt(A_t / pi)) ^ -.6 + 1.25;  % contraction ratio (I forget where I got this equation from, please don't kill me)
    A_c = A_t * con_ratio;                              % chamber area [in^2]

% size based off a calculated contraction ratio
elseif converging_method == "exit"
    A_c = A_e;                    % chamber area [in^2]
    con_ratio = exp_ratio;        % contraction ratio                              

% throw an error if no converging method is selected
else
    error("No converging method selected.")
end

D_e = 2 * sqrt(A_e / pi); % exit diameter [in]
D_c = 2 * sqrt(A_c / pi); % chamber diameter [in]
R_t = D_t / 2;            % throat radius [in]
R_e = D_e / 2;            % exit radius [in]

%% Generate Nozzle Contour
[x_contour, r_contour, L_c, L_total] = engineContour(geometry_type, bell_pct, R_t, exp_ratio, con_ratio, conv_angle, conical_half_angle, L_crude_throat, L_star, bar_size);

%% Cooling Calculations
if sizeCooling
    % factor in film cooling
    if film_pct
    %    [] = sizeFilm();
    end

    % size regenerative cooling
    if cooling_type == "regen"
        %[] = sizeRegen(x_contour, y_contour, R_t, nozzle_regen_pct, M, gamma, P, T, rho, mu, Pr, Mw, k, son, vel);

    % size ablative cooling
    elseif cooling_type == "ablative"
        %[] = sizeAblative(x_contour, y_contour);

    % analyze heat sink
    elseif cooling_type == "heat sink"
        %[] = sizeHeatSink(x_contour, y_contour);
    end 
end

%% Injector Sizing
if sizeInjector
    % size coaxial swirler
    if injector_type == "swirl"
        %[] = sizeCoaxSwirl();

    % size coaxial shear
    elseif injector_type == "shear"
        %[] = sizeCoaxShear();

    % size pintle
    elseif injector_type == "pintle"
        %[] = sizePintle();

    % size impinging jets
    elseif injector_type == "impinging"
        %[] = sizeImpingingJets();
    end

    %A_ox = m_dot_OX / (orificeDisCoef*sqrt(2 .* 32.2 .* densityOx .* dpOx*144)) * 144; % Lox total flow area [in^2]
    %A_fuel = m_dot_INJ_FUEL / (orificeDisCoef*sqrt(2 .* 32.2 .* densityFuel .* dpFuel*144))* 144; % Fuel total flow area [in^2]

    fprintf('\n-------- Injector Outputs --------\n')
    %disp(['             Fuel Flow Area (in^2): ' num2str(P_c)]);
    %disp(['         Oxidizer Flow Area (in^2): ' num2str(F)]);
    fprintf('-------- Injector Outputs --------\n')
end

%% Fluid Systems Sizing
% calculate required mass of propellants
if sizeFluids

    oxidizer_density = refpropm('D', 'T', oxidizer_temp, 'P', P_c * u.PSI2MPA, oxidizer_REFPROP) * u.KGM32LBIN3; % [lb/in^3]
    fuel_density = refpropm('D', 'T', fuel_temp, 'P', P_c * u.PSI2MPA, fuel_REFPROP) * u.KGM32LBIN3; % [lb/in^3]

    % size fluid system based off of tank size
    if vol_tank

        % oxidizer is limiting propellant
        if m_dot_OX / oxidizer_density > m_dot_FUEL / fuel_density 
            oxidizer_mass = oxidizer_density * vol_tank;
            oxidizer_vol = oxidizer_mass / oxidizer_density;
            time_burn = oxidizer_mass / m_dot_OX;
            fuel_mass = oxidizer_mass / OF;
            fuel_vol = fuel_mass / fuel_density;

        % fuel is limiting propellant
        else 
            fuel_mass = fuel_density * vol_tank;
            fuel_vol = fuel_mass / fuel_density;
            time_burn = fuel_mass / m_dot_FUEL;
            oxidizer_mass = fuel_mass * OF;
            oxidizer_vol = oxidizer_mass / oxidizer_density;
        end

    % size fluid system based off of burn time
    else
        oxidizer_mass = time_burn * m_dot_OX; 
        oxidizer_vol = oxidizer_mass / oxidizer_density;
        fuel_mass = time_burn * m_dot_FUEL;
        fuel_vol = fuel_mass / fuel_density;
    end
    impulse = time_burn * F; % total impulse [N-s]

    % output fluid system information
    fprintf('\n-------- Fluid System Outputs -------\n')
    disp(['                Mass Oxidizer (lb): ' num2str(oxidizer_mass)]);
    disp(['                    Mass Fuel (lb): ' num2str(fuel_mass)]);
    disp(['            Volume Oxidizer (in^3): ' num2str(oxidizer_vol)]);
    disp(['                Volume Fuel (in^3): ' num2str(fuel_vol)]);
    disp(['           Maximum Burn Time (sec): ' num2str(time_burn)]);
    disp(['               Total Impulse (sec): ' num2str(impulse)]);
    fprintf('-------- Fluid System Outputs -------\n')

end

%% Structural Sizing
if sizeStructures

    face_force = P_c * A_c; % injector face bolt force due to chamber pressure [lbf]

    % output fluid system information
    fprintf('\n--------- Structures Outputs --------\n')
    disp(['      Force on Injector Face (lbf): ' num2str(face_force)]);
    fprintf('--------- Structures Outputs --------\n')

end

%% Output Results

% figure(2)
% hold on
% 
% subplot(2,2,1)
% plot(x_contour, M, 'b');
% title('Mach')
% xlabel('Inches X')
% ylabel('Mach Number')
% grid on
% axis([x_contour(1) x_contour(end) .9*min(M) 1.1*max(M)])
% 
% subplot(2,2,2)
% plot(x_contour, P, 'g');
% title('Pressure')
% xlabel('Inches X')
% ylabel('Pressure [psi]')
% grid on
% axis([x_contour(1) x_contour(end) .9*min(P) 1.1*max(P)])
% 
% subplot(2,2,3)
% plot(x_contour, T, 'r');
% title('Temperature')
% xlabel('Inches X')
% ylabel('Temperature [R]')
% grid on
% axis([x_contour(1) x_contour(end) .9*min(T) 1.1*max(T)])
% 
% subplot(2,2,4)
% plot(x_contour, r_contour, 'black');
% title('Radius')
% xlabel('Inches X')
% ylabel('Inches Y')
% grid on
% axis([x_contour(1) x_contour(end) 0 x_contour(end)-x_contour(1)])

hold off