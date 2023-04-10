clear;
clc;

min_throttle_pct = .2;
F_max = 500;
A_t = 1.607;
P_c = 74.02;
P_e = 4.40;
P_a = 14.7;
fuel = 'C3H8O,2propanol';
fuel_weight = 0;
fuel_temp = 293.15;
oxidizer = 'O2(L)';
oxidizer_temp = 90.17;
OF = 1.3;
CEA_input_name = 'test';
eff_c_f = .9;
eff_c_star = .94;
g = 32.174; % gravitational acceleration [ft/s^2]
exp_ratio = 3.35;

[cea_c_star, ~, ~, ~, gamma, ~, ~, ~, ~, ~, ~, ~, ~, ~] = RunCEA(P_c, P_a, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 0, 1);
gamma = gamma(1);
c_f = (sqrt(((2*gamma^2)/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1-(P_e/P_c)^((gamma-1)/gamma))) + (P_e-P_a)*exp_ratio/P_c) * eff_c_f;
F = A_t * P_c * c_f;