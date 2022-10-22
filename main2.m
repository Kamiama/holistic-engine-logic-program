% -------------------------------------
%%   Ultimate Engine Sizing Code v1.0
% Author: Kamon Blong (kblong@purdue.edu)
% First Created: 7/10/2022
% Last Updated: 
%{ 
Assumptions: Currently assumes a stationary test engine with no regard for
changes in gravity or atmospheric pressure. Additionally assumes an
ideal rocket engine (perfect propellant mixing, gaseous flow, perfect gas
law, adiabatic, neglects friction and boundry layer effects, no
discontinuities, steady mDot, no transience, all exhaust contributes to
thrust, characteristics are uniform at any point along an axial slice of
the engine, and chemical equilibrium is constant axially throughout the engine 
(Sutton 46)). Without tinputing a specific inlet temperature for
propellants, NASA CEA (and therefore this program) also assumes that
ordinary propellants are stored at room temperature and cryogens at their
boiling point.
%}

clear;
%clc;

%% Identify OS & Generate Paths 
% identify OS
if ismac || isunix
    function_path = append(pwd, '/bin');
    cea_path = append(pwd, '/cea');
    output_path = append(pwd, '/output');
elseif ispc
    function_path = append(pwd, '\bin');
    cea_path = append(pwd, '\cea');
    output_path = append(pwd, '\output');
end
% assign paths
main_path = cd;
addpath(function_path);
addpath(cea_path);
addpath(output_path);
u = convertUnits;

%% Initialization
% initialize output document
run_description = input('Run Description: ', 's'); % run description that will be added to output file header
write_date = datestr(now, '_mm-dd-yy_HH:MM'); % date that will be added to output file header
file_name = strcat(run_description, write_date); % output file name

%% Import Data
% import input properties
input_data = readmatrix('input.xlsx','NumHeaderLines',1);

%% Parse Variables & Define Constants
% parse & convert input properties

    % performance properties
    F = 2000; % engine thrust [lbf]
    P_c = 170; % chamber pressure [psi]
    P_a = 10.5; % ambient pressure [psi]
    P_e = P_a; % exit pressure [psi] (equals P_a for optimal expansion)
    eff_c_star = 0.9; % c* efficiency
    eff_c_f = 0.8; % cF efficiency
    
    % fuel properties
    fuel = {'C2H5OH(L)', 'C8H18(L),isooct'}; % fuel formula for NASA CEA
    fuel_weight = [98, 2]; % fuel weights, does nothing now
    fuel_density = 0.029046343; % fuel density [lb/in^3]
    oxidizer = 'O2(L)'; % oxidizer formula for NASA CEA
    oxidizer_density = 4.122E-2; % oxidizer density [lb/in^3]
    fuel_temp = 0; % inlet fuel temperature [K]
    oxidizer_temp = 0; % inlet oxidizer temperature [K]
    OF = 1.33; % O/F ratio
    time_burn = 2; % burn time [sec]
    
    % geometry properties
    geometry_type = "conical";
    bell_pct = .8; % percent of bell nozzle
    conical_half_angle = 15; % conical half angle [deg]
    L_crude_throat = .25; % length of straight throat section [in]

    theta_i = 20.75; % nozzle expansion (initial) angle [degrees]
    theta_e = 15.25; % nozzle exit angle [degrees] 
    % derived emperically from http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf,
    % matching graph to an equation is a possible improvement for the future

    D_t = 0;
    L_star = 75; % L*, characteristic combustion length [in]
    conv_angle = 30; % convergence 
    con_ratio = 0; % contraction ratio
    D_c = 7; % chamber diameter [in]
    
    % injector properties
    
    % cooling properties
    cooling_type = "regen"; % cooling method (regen or ablative)
    film_pct = .105; % film cooling percent as fraction of fuel mass flow rate
    nozzle_regen_pct = 1; % percent of nozzle to be regeneratively cooled

    % fluid properties
    vol_tank = 0; % volume of propellant tanks [in^3]
    p_tank = 1000; % tank pressure [psi]

% parse tooling and stock sizes
    bar_size = 0; % bar stock size [in]

% define constants
    g = 32.174; % gravitational constant [ft/sec^2]
    R_N2 = 0; % 

%% Run NASA CEA
% sends Pc, Pex, and fuel characteristics from input.xlsx as inputs
% returns Isp and c* as outputs

% initialize CEA input/output files
CEA_name_prefix = strcat('CEA');
CEA_input_name = append(CEA_name_prefix, '.inp');
CEA_output_name = append(CEA_name_prefix, '.out');

% run CEA
[cea_c_star, cea_isp, exp_ratio, M, gamma, P, T, rho, mu, Pr, Mw, k, son] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name);

%% Intermediate Calculations

% size normally
if D_t == 0
    % Correct values & find effective exhaust velocity
    isp = cea_isp * eff_c_star * eff_c_f;   % expected isp [sec]
    c_star = cea_c_star * eff_c_star; % expected c* [ft/s]
    c = isp * g; %+ A_e * (P_e - P_a); % effective exhaust velocity [ft/s]
    
    % Calculate mass flow rate of fuel, oxidizer, and film cooling
    m_dot = F / isp;                           % TOTAL mass flow rate [lbm/s]
    m_dot_FUEL = m_dot / (OF + 1);                 % FUEL mass flow rate [lbm/s]
    m_dot_OX = m_dot_FUEL * OF;                    % OX mass flow rate [lbm/s]
    
    m_dot_FILM = m_dot_FUEL * film_pct;            % FILM mass flow rate [lbm/s]
    m_dot_INJ_FUEL = m_dot_FUEL * (1 - film_pct);  % regularly injected FUEL mass flow rate [lbm/s]

    % run a second time for film cooling
    m_dot_FUEL_2 = m_dot_FUEL * (1 - film_pct);
    OF_2 = OF * m_dot_FUEL / m_dot_FUEL_2;
    m_dot_FUEL = m_dot_FUEL + m_dot_FILM;
    [cea_c_star, cea_isp, exp_ratio, M, gamma, P, T, rho, mu, Pr, Mw, k, son] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF_2, 0, 0, CEA_input_name);
    
    % Correct values & find effective exhaust velocity
    isp = cea_isp * eff_c_star * eff_c_f;   % expected isp [sec]
    c_star = cea_c_star * eff_c_star; % expected c* [ft/s]
    c = isp * g; %+ A_e * (P_e - P_a); % effective exhaust velocity [ft/s]

% size based off of throat size
else
    A_t = (D_t / 2) ^ 2 * pi;
    A_t = m_dot * c_star / P_c / g;
end

% Assign variables from CEA

% output performance information
fprintf('\n-------- Performance Outputs --------\n')
disp(['                      Thrust (lbf): ' num2str(F)]);
disp(['            Chamber Pressure (psi): ' num2str(P_c)]);
disp(['               Exit Pressure (psi): ' num2str(P_e)]);
disp(['                              Fuel: ' fuel]);
disp(['                          Oxidizer: ' oxidizer]);
disp(['                         O/F Ratio: ' num2str(OF_2)]);
fprintf('\n')
disp(['                     CEA Isp (sec): ' num2str(cea_isp)]);
disp(['                Expected Isp (sec): ' num2str(isp)]);
disp([' Effective Exhaust Velocity (ft/s): ' num2str(c)]);
disp(['    Mass Flow Rate (total) (lbm/s): ' num2str(m_dot)]);
disp(['     Mass Flow Rate (fuel) (lbm/s): ' num2str(m_dot_FUEL)]);
disp(['     Mass Flow Rate (film) (lbm/s): ' num2str(m_dot_FILM)]);
disp([' Mass Flow Rate (oxidizer) (lbm/s): ' num2str(m_dot_OX)]);
disp(['         Adiabatic Flame Temps (R): ' num2str(T(1)) ' ' num2str(T(2))]);
disp(['                         Exit Mach: ' num2str(M(end))]);
disp(['                            Gammas: ' num2str(gamma(1)) ' ' num2str(gamma(2))]);
fprintf('-------- Performance Outputs --------\n')

%% Preliminary Engine Geometry Calculations
A_t = m_dot * c_star / P_c / g; % throat area [in^2]
A_e = A_t * exp_ratio; % exit area [in^2]
D_t = 2 * sqrt(A_t / pi); % throat diameter [in]
if con_ratio
    A_c = A_t * con_ratio; % chamber area [in^2]
elseif D_c
    A_c = pi * (D_c / 2) ^ 2; % chamber area [in^2]
    con_ratio = A_c / A_t; % contraction ratio
else
    con_ratio = 8 * (2 * sqrt(A_t / pi)) ^ -.6 + 1.25; % contraction ratio
    A_c = A_t * con_ratio; % chamber area [in^2]
end

D_e = 2 * sqrt(A_e / pi); % exit diameter [in]
D_c = 2 * sqrt(A_c / pi); % chamber diameter [in]
R_t = D_t / 2; % throat radius [in]
R_e = D_e / 2; % exit radius [in]

%% Generate Nozzle Contour
[x_contour, r_contour, L_c, L_total] = engineContour(geometry_type, bell_pct, R_t, theta_i, theta_e, exp_ratio, con_ratio, conv_angle, conical_half_angle, L_crude_throat, L_star, bar_size);

%% Cooling Calculations

% continueScript; % ensure user wants to continue before running cooling section (it takes a little bit)

% % size respective cooling types
% if cooling_type == "regen"
%     subsonic_area_ratios = (pi * r_contour(x_contour < 0) .^ 2) / A_t;
%     supersonic_area_ratios = (pi * r_contour(x_contour > 0) .^ 2) / A_t;
%     [cea_c_star, cea_isp, exp_ratio, M, gamma, P, T, rho, mu, Pr, Mw, k, son] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, subsonic_area_ratios, supersonic_area_ratios, CEA_input_name);
%     vel = son .* M; % axial chamber velocity [ft/s]
%     %[] = regenCalculations(x_contour, y_contour, R_t, nozzle_regen_pct);
% elseif cooling_type == "ablative"
%     %[] = ablativeCalculations(x_contour, y_contour);
% end
% 
% % factor in film cooling
% if film_pct
% %    [] = filmCalculations();
% end

%% Injector Sizing

%% Fluid Systems Sizing
% calculate required mass of propellants
if vol_tank
    if m_dot_OX / oxidizer_density > m_dot_FUEL / fuel_density % oxidizer is limiting propellant
        oxidizer_mass = oxidizer_density * vol_tank;
        oxidizer_vol = oxidizer_mass / oxidizer_density;
        time_burn = oxidizer_mass / m_dot_OX;
        fuel_mass = oxidizer_mass / OF;
        fuel_vol = fuel_mass / fuel_density;
    else % fuel is limiting propellant
        fuel_mass = fuel_density * vol_tank;
        fuel_vol = fuel_mass / fuel_density;
        time_burn = fuel_mass / m_dot_FUEL;
        oxidizer_mass = fuel_mass * OF;
        oxidizer_vol = oxidizer_mass / oxidizer_density;
    end
else
    oxidizer_mass = time_burn * m_dot_OX; 
    oxidizer_vol = oxidizer_mass / oxidizer_density;
    fuel_mass = time_burn * m_dot_FUEL;
    fuel_vol = fuel_mass / fuel_density;
end
impulse = time_burn * F; % total impulse [N-s]

% calculate regulator flow rates
m_dot_N2 = 0;

% output fluid system information
fprintf('\n-------- Fluid System Outputs -------\n')
disp(['                Mass Oxidizer (lb): ' num2str(oxidizer_mass)]);
disp(['                    Mass Fuel (lb): ' num2str(fuel_mass)]);
disp(['            Volume Oxidizer (in^3): ' num2str(oxidizer_vol)]);
disp(['                Volume Fuel (in^3): ' num2str(fuel_vol)]);
disp(['           Maximum Burn Time (sec): ' num2str(time_burn)]);
disp(['               Total Impulse (sec): ' num2str(impulse)]);
fprintf('-------- Fluid System Outputs -------\n')

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