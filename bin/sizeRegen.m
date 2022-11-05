%% HELP Regenerative Cooling Sizing Code
% Author: Kamon Blong (kamon.blong@gmail.com)
% First Created: 10/23/2022
% Last Updated: 10/23/2022

function [] = sizeRegen(x_contour, y_contour, R_t, nozzle_regen_pct)

%{ 
Description: Currently sizes a single-pass regen system that injects fuel at a 
    variable distance downstream of the throat and carries it through up to the injector face. 
    Channels have a constant width but have the option to change channel height to 
    decrease hydraulic diameter and increase fluid velocity (thereby increasing cooling
    effectiveness) in key locations. This method allows for optimal machining via CSJ
    with a slitter-saw CNC attatchment or via powder bed 3D printing.

Inputs:
- 

Outputs: 
- 

Assumptions:
- Steady state
- No backside wall heat transfer
- Equally distributed temperature inside channels
- Wicking heat into fuel doesn't change bulk gas temperature
%}

%% Parse variables

% define constants
g = 9.81; % gravitational accelaration [m/s]

%% Parse material properties

%% Calculate minimum wall thickness
% determine optimal wall thickness based on throat conditions (because this is the most extreme location in the engine)

% initial CEA run to find throat conditions (OUTPUTS IN METRIC)
[~, ~, ~, M, gamma, P, T, rho, mu, Pr, Mw, k, son, cp] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 1, 0);

% perform initial heat transfer and structural calculations to find wall thickness
r = Pr ^ (1 / 3); % recovery factor - biased towards larger engines, very small engines should use Pr^.5 (Heister 196).
T_r = T * (1 + r * (gamma - 1) / 2 * M ^ 2) / (1 * (gamma - 1) / 2 * M ^ 2); % recovery temp (adiabatic wall temp) - corrects for compressible boundry layers (Huzel & Huang 85).


sigma = (.5 * T_wg / T_c * (1 + (gamma - 1) / 2 * M ^ 2) + .5) ^ -.68 * (1 + (gamma - 1) / 2 * M ^ 2) ^ -.12; % film coefficient correction factor (Huzel & Huang 86).
h_g = (.026 / D_t ^ .2) * (mu ^ .2 * cp / Pr ^ .6) * (P_c * g / c_star) ^ .8 * (D_t / radius_throat) ^ .1 * (A_t / A_t) ^ .9 * sigma; % film coefficient - bartz equation (Huzel & Huang 86).
q_dot = 

%% Set inlet properties
% necessary initial flow properties: m_dot, pressure, temperature (all other fluid properties can be derived from pressure and temperature
% necessary geometric properties: channel shape (I will kill you if it's anything but a square), channel diameter, wall thickness (calculated via structural limitations)
% necessary thermal properties: maximum hot wall temperature for your wall material

max_T_wg = 120; % maximum hot wall temperature [K]

% calculate number of channels and mass flow rate through each
num_channels = pi * (D_t + .8 * (D_ch + 2 * t_w)) / (D_ch + 2 * t_w); % number of channels (Heister 206).
m_dot_CHANNEL = m_dot / num_channels;


subsonic_area_ratios = (pi * r_contour(x_contour < 0) .^ 2) / A_t;
supersonic_area_ratios = (pi * r_contour(x_contour > 0) .^ 2) / A_t;
[~, ~, ~, M, gamma, P, T, rho, mu, Pr, Mw, k, son] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 0);

%% Begin iterative sizing
% start at manifold and step up along chamber length

% engine_pos = final_pos
% for engine_pos ~= final_pos
    % make T_wg guess (gas-side wall temperature)

    % calculate h_g (film coefficient) from q_dot = h_g * (T_r - T_wg_guess) - T_r is the "gas temperature" axially along the chamber from CEA
%% Regen Sizing (Heister Procedure)
%%Initialization
stepnum = 100;
Twgi = 1000; %Initial Guess of Gas side Temperature
gasheattransfer = 1000;
liqheatransfer = 0;

%Geometric
wallthick = 100; %Wall Thickness
d = 5; %?????? d in num channel equation
engine_diameter = 10; %????? In num channel equation
Alocal = []; %Local Cross Sectional Areas of Engine
A_t = 10;  %Area at throat
dlocal = [];  %Local Cross Sections of Channels

Dstar = %Diameter at Nozzle Throat
characteristic_length = 10 %Longest length of Channels 


%Coolant properties
viscosity = 10; %Viscosity
Cp = 100; %Specific heat at constant pressure
kc = 100; %Thermal Conductivity of Coolant
Tl = 293;

%Engine Properties
Pi = 20; %Inflow Pressure
Ti = 200; %Inflow Temp
Tgas = []; %Gas temperature from 1-D 
mdotf = 30; %Coolant/fuel mass flow
chamber_pressure = 100; %Chamber Pressure
Cstar = 1000; %Characteristic Velocity
throat_radius_curvature = 10; %Throat radius Curvature

kw = 100; %Thermal Conductivity of Wall


%Step 2: Calculate Number of channels and channel mass flow rate
numchannels = pi * (engine_diameter + 0.8 * (Dstar + 2 * wallthick)) / (Dstar + 2 * wallthick); %Number of Channels (EQ 6.30) (Change the coefficient later)
mdotchan = mdotf / numchannels; %Mass flow of channel (EQ 6.31)
%Step 3: Begin Stepping Down tube/channel
for i = [1:stepsize]
    [~, ~, ~, M, gamma, P, T, rho, mu, Pr, Mw, k, son, cp] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 1, 0);

    Pr_gas = densityfs * viscosityfs * D;
    Pr_liquid = density(Tl(i));
    while liqheattransfer < gasheattransfer * upperbound || liqheattransfer > gasheattransfer * lowerbound
        %Step 5: Calculate Gas film coefficient and heat transfer 
        hg = 0.026 / Dstar^0.2 * (viscosity^0.2 * Cp / prandtl^0.6) * ...
        (chamber_pressure * g / Cstar)^0.8 * (Dstar / throat_radius_curvature)^0.1 * (Astar/Alocal[i])^0.9 * (meandens / freedens)^0.8 * (meanvisc / freevisc)^0.2; %Gas Film Coefficent Bartz 
        gasheattransfer = heatcoef * (Tgas(i) - Twgi);  %Gas Heat Transfer (EQ 6.16)

        %Step 6: Calculate Liquid Wall Temperature from conduction 
        Twl(i) = gaswalltemp - gasheattrans * wallthick / thermconduct; %Liquid Wall Temp (EQ 6.29)

        %Step 7: Compute Liquid Film Coefficient
        hl = (a * re^m * pr^n * (visfree / viscwall)^b) * k / L;    %Liquid Film Coefficient (EQ 6.19) 

        %Step 8: Compute Heat Flux. Compare it to step 5
        liqheatransfer = hl * (Tl(i) - Twl(i)); %Liquid Heat Transfer (EQ  6.29)

        %Step 4: Guess Gass Wall Temperature

        twg(i) = Tgas(i) - liqheatransfer / hg; %Guess Gass Wall temp using liquid heat transfer value

        %Step 9: Run loop until step 5/8 get same value for heat flux
    end
    Tl(i + 1) = T(i) + 1 / (mdotchan * Cp);

