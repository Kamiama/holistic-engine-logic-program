% temporary cooling script
function [] = coolingTemp(x_contour, r_contour, R_t, nozzle_regen_pct, m_dot, P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF)

%% Parse Variables
A_t = R_t ^ 2 * pi; % throat area [in^2]
D_t = R_t * 2; % throat diameter [in]
m_dot_FUEL = m_dot / (OF + 1); % fuel mass flow rate [lb/s]

%% Throat Conditions
% run CEA at throat
[c_star_t, ~, ~, M_t, gamma_t, P_t, T_t, rho_t, mu_t, Pr_gas_t, Mw_t, k_t, son_t, cp_t] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 1, 1);

% run bartz at throat
sigma = (.5 * T_t / T_t * (1 + (gamma_t - 1) / 2 * M_t ^ 2) + .5) ^ -.68 * (1 + (gamma_t - 1) / 2 * M_t ^ 2) ^ -.12; % film coefficient correction factor (Huzel & Huang 86).
h_g = (.026 / D_t ^ .2) * (mu_t ^ .2 * cp_t / Pr_gas_t ^ .6) * (P_c * g / c_star_t) ^ .8 * (D_t / radius_throat) ^ .1 * (A_t / A_t) ^ .9 * sigma; % film coefficient - bartz equation (Huzel & Huang 86).

num_channels = pi * (D_t + 0.8 * (D_t + 2 * t_w)) / (D_t + 2 * t_w); % number of Channels (EQ 6.30)

%sigma = (.5 * T_wg / T_c * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) + .5) ^ -.68 * (1 + (gamma(i) - 1) / 2 * M(i) ^ 2) ^ -.12; % film coefficient correction factor (Huzel & Huang 86).
%h_g = (.026 / D_t ^ .2) * (mu ^ .2 * cp / Pr ^ .6) * (P_c * g / c_star) ^ .8 * (D_t / radius_throat) ^ .1 * (A_t / A_t) ^ .9 * sigma; % film coefficient - bartz equation (Huzel & Huang 86).
            