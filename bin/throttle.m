%% Engine Throttle Analysis
% Author: Kamon Blong (kamon.blong@gmail.com)
% First Created: 3/28/2023
% Last Updated: 4/10/2023

function [] = throttle(min_throttle_pct, F_max, eff_c_star, eff_c_f, P_c_max, P_e, P_a, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, CEA_input_name, resolution)

%{ 
Description: Takes throttle percent as an input and performs analysis to
    deduce engine performance over a given throttle range following the
    procedure outlined in chapter 4 of http://ltu.diva-portal.org/smash/get/diva2:1732226/FULLTEXT01.pdf

Assumptions:
    - Constant OF ratio across throttle range (possible improvement for
    the future)
    - Unchanging atmospheric conditions

Inputs:
    min_throttle_pct - minimum throttle percent as a decimal between 0 and 1
    F_max - maximum thrust [lbf]
    r_t - throat radius [in]
    P_c_max - maximum chamber pressure [psi]
    P_e - exit pressure [psi]
    fuel - fuel string for CEA [N/A]
    fuel_temp - fuel temperature [K]
    fuel_weight - fuel weights for fuel mixtures [N/A]
    oxidizer - oxidizer string for CEA [N/A]
    oxidizer_temp - oxidizer temperature [K]
    OF - oxidizer to fuel ratio [N/A]

Outputs: 
    - None
%}

% clear;
% clc;
% close all;

%% Identify OS & Generate Paths 
% add paths
current_dir = fileparts(mfilename('fullpath'));
main_dir = fileparts(current_dir);
cea_path = fullfile(main_dir, 'cea');
addpath(cea_path);

%% Initializion & Variable Definition

% % debugging values
% min_throttle_pct = .25;
% F_max = 550;
% P_c_max = 220;
% P_e = 18;
% P_a = 14.7;
% fuel = 'C3H8O,2propanol';
% fuel_weight = 0;
% fuel_temp = 293.15;
% oxidizer = 'O2(L)';
% oxidizer_temp = 90.17;
% OF = 1.3;
% CEA_input_name = 'test';
% eff_c_f = .9;
% eff_c_star = .9;

g = 32.174; % gravitational acceleration [ft/s^2]

P_sep = P_a * .4; % flow separation condition
debug = 0;

% initialize throttle range and desired variable matrices
throttle_pct = linspace(1, min_throttle_pct, resolution); % throttle percent matrix
Fs = F_max .* throttle_pct; % thrust matrix with correct linear correlation
Fs_temp = zeros(1, resolution); % temporary thrust matrix for iterative solving
P_cs = P_c_max .* throttle_pct; % initial chamber pressure matrix with initial assumption of linear correlation
P_es = zeros(1, resolution);
m_dots = zeros(1, resolution);
isps = zeros(1, resolution);
rel_isps = zeros(1, resolution);


%% Get Parameters at Full Throttle

% run CEA at full throttle conditions (F_max, P_c_max)
[cea_c_star, cea_isp, exp_ratio_cea, ~, gamma, ~, ~, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c_max, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 1, CEA_input_name);

% extract CEA values
c_star_max = cea_c_star * eff_c_star;
gamma_max = gamma(1); % chamber gamma

% solve for expansion ratio (Sutton EQ 3.25)
exp_ratio = 1/(((gamma_max+1)/2)^(1/(gamma_max-1)) * (P_e/P_c_max)^(1/gamma_max) * ((gamma_max+1)/(gamma_max-1)*(1-(P_e/P_c_max)^((gamma_max-1)/gamma_max)))^.5);

% solve for thrust coefficient (Sutton EQ 3.30)
c_f_max = (sqrt(((2*gamma_max^2)/(gamma_max-1)) * (2/(gamma_max+1))^((gamma_max+1)/(gamma_max-1)) * (1-(P_e/P_c_max)^((gamma_max-1)/gamma_max))) + (P_e-P_a)*exp_ratio/P_c_max) * eff_c_f;
%c_f_2 = (sqrt(((2*gamma_max^2)/(gamma_max-1)) * (2/(gamma_max+1))^((gamma_max+1)/(gamma_max-1)) * (1-(P_e/P_c_max)^((gamma_max-1)/gamma_max))) + (P_e-P_a)*exp_ratio_cea/P_c_max) * eff_c_f;

% basic calculations
A_t = F_max / c_f_max / P_c_max;
m_dot = A_t / c_star_max * P_c_max * g; % mass flow rate at full throttle [lb/s]
isp_max = c_star_max * c_f_max / g;

% isp_2 = cea_isp * eff_c_f * eff_c_star;
% m_dot_2 = F_max / isp_2;
% A_t_2 = m_dot_2 * c_star_max / P_c_max / g;

% isp_3 = c_star_max * c_f_2 / g;
% A_t_3 = F_max / c_f_2 / P_c_max;

% Initial Solution for Throttle Range

% initialize loop variables
index = 1;
%F = F_max;
tolerance = .02;
max_iter = 1000;

% iteratively loop through entire P_c matrix
for P_c = P_cs
    % initialize loop variables
    iter = 0;
    P_c_guess = P_c;
    P_c_mx = P_c_max + 50;
    P_c_mn = P_c;

    % loop
    while abs(Fs(index) - Fs_temp(index)) > tolerance && iter < max_iter
        % run CEA for given P_c_new
        [cea_c_star, ~, ~, ~, gamma, ~, ~, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 1, 1, CEA_input_name);
        
        % extract CEA values
        gamma = gamma(1);
        
        % solve for P_e numerically (Sutton EQ 3.25)
        func = @(P_e) 1/(((gamma+1)/2)^(1/(gamma-1)) * (P_e/P_c_guess)^(1/gamma) * ((gamma+1)/(gamma-1)*(1-(P_e/P_c_guess)^((gamma-1)/gamma)))^.5) - exp_ratio;
        x0 = P_a * 2 / 3; % initial guess for P_e
        if debug
            options = optimset('Display', 'iter'); % display solver progress
        else
            options = optimset('Display', 'off'); % display solver progress
        end
        P_e = fsolve(func, x0, options); % use the numerical solver fsolve to find the solution for P_e
    
        % solve for thrust coefficient (Sutton EQ 3.30)
        c_f = (sqrt(((2*gamma^2)/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1-(P_e/P_c_guess)^((gamma-1)/gamma))) + (P_e-P_a)*exp_ratio/P_c_guess) * eff_c_f;
  
        % calculate thrust
        Fs_temp(index) = A_t * c_f * P_c_guess;

        % update loop
        if Fs(index) > Fs_temp(index)
            P_c_mn = P_c_guess;
        else
            P_c_mx = P_c_guess;
        end
        P_c_guess = (P_c_mx + P_c_mn) / 2;
        iter = iter + 1;
    end

    % record variables to matrices
    Fs(index) = Fs_temp(index);
    m_dots(index) = A_t / (cea_c_star * eff_c_star) * P_c_guess * g;
    isps(index) = cea_c_star * c_f * eff_c_star / g;
    rel_isps(index) = isps(index) / isp_max;
    P_cs(index) = P_c_guess;
    P_es(index) = P_e;

    % increase matrix index
    index = index + 1;
end

% Calculate Lines of Best Fit

% thrust
coefficients = polyfit(throttle_pct, Fs, 1);
F_eq = sprintf('y = %.2f x + %.2f', coefficients(1), coefficients(2));

% chamber pressure
coefficients = polyfit(throttle_pct, P_cs, 1);
P_c_eq = sprintf('y = %.2f x + %.2f', coefficients(1), coefficients(2));

% exit pressure
coefficients = polyfit(throttle_pct, P_es, 1);
P_e_eq = sprintf('y = %.2f x + %.2f', coefficients(1), coefficients(2));

% mass flow rate
coefficients = polyfit(throttle_pct, m_dots, 1);
m_dots_eq = sprintf('y = %.2f x + %.2f', coefficients(1), coefficients(2));

%% Present Results

% text outputs
fprintf('\n----- Throttle Analysis Outputs -----\n')
disp(['              Minimum Throttle (%): ' num2str(min_throttle_pct)]);
disp(['              Minimum Thrust (lbf): ' num2str(min(Fs))]);
disp(['    Minimum Chamber Pressure (psi): ' num2str(min(P_cs))]);
disp(['       Minimum Exit Pressure (psi): ' num2str(min(P_es))]);
disp(['   Flow Separation Condition (psi): ' num2str(P_sep)]);
disp(['    Minimum Mass Flow Rate (lbm/s): ' num2str(min(m_dots))]);
disp(['                 Minimum Isp (sec): ' num2str(min(isps))]);
fprintf('\n           Thrust line of best fit: \n')
disp(['               ' F_eq]);
fprintf('\n Chamber pressure line of best fit: \n')
disp(['               ' P_c_eq]);
fprintf('\n    Exit pressure line of best fit: \n')
disp(['               ' P_e_eq]);
fprintf('\n   Mass flow rate line of best fit: \n')
disp(['               ' m_dots_eq]);
fprintf('\n')
fprintf('----- Throttle Analysis Outputs -----\n')

figure('Name', 'Throttle Analysis')
hold on

% chamber pressure plot
subplot(2,2,1)
plot(throttle_pct, P_cs, 'blue');
title('Chamber Pressure vs Throttle')
xlabel('Throttle %')
ylabel('Chamber Pressure [psi]')
axis([throttle_pct(end) throttle_pct(1) P_cs(end) P_cs(1)])
grid on

% mass flow rate plot
subplot(2,2,2)
plot(throttle_pct, m_dots, 'blue');
title('Mass Flow Rate vs Throttle')
xlabel('Throttle %')
ylabel('Mass Flow Rate [lb/s]')
axis([throttle_pct(end) throttle_pct(1) m_dots(end) m_dots(1)])
grid on

% exit pressure plot
subplot(2,2,3)
plot(throttle_pct, P_es, 'blue');
title('Exit Pressure vs Throttle')
xlabel('Throttle %')
ylabel('Exit Pressure [psi]')
if P_es(end) > P_sep
    axis([throttle_pct(end) throttle_pct(1) P_sep-1 P_es(1)])
else
    axis([throttle_pct(end) throttle_pct(1) P_es(end) P_es(1)])
end
% flow separation line
line([-100 100],[P_sep P_sep], 'Color', 'red', 'LineStyle', '--')
legend('', 'Flow Separation Condition', 'Location', 'Northwest')
grid on

% isp plot
subplot(2,2,4)
% relative isps
yyaxis left
plot(throttle_pct, rel_isps, 'blue');
title('Isp vs Throttle')
xlabel('Throttle %')
ylabel('Isp / Isp^n^o^m')
grid on
set(gca, 'XLim', [throttle_pct(end), throttle_pct(1)], 'YLim', [rel_isps(end), rel_isps(1)], 'Ycolor', 'k')
% actual isps
yyaxis right
isp_axis = round(linspace(min(flip(rel_isps .* isp_max)), max(flip(rel_isps .* isp_max)), 5) .* 10) ./ 10;
set(gca,'YTick', isp_axis, 'YLim', [min(isp_axis), max(isp_axis)], 'Ycolor', 'k')
ylabel('Isp [sec]')

sgtitle("Throttled Chamber Performance: " + F_max + "lbf, " + OF + " OF ratio, " + min_throttle_pct * 100 + "% min throttle")


