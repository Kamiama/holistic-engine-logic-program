%% HELP Regenerative Cooling Sizing Code
% Authors: Kamon Blong (kamon.blong@gmail.com), Jan Ayala, Andrew Radulovich
% First Created: 10/23/2022
% Last Updated: 11/14/2022

function [] = sizeRegen(x_contour, y_contour, R_t, nozzle_regen_pct, mdotf, P_c, P_e, Oxidizer, Fuel, OF_ratio, wall_material)
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

%% Parse variables

    % geometric properties
        t_w = 100; % wall thickness [in]
        D_c = 10; % 
        A_local = []; %Local Cross Sectional Areas of Engine
        A_t =  m_dot * c_star / P_c / g;    % throat area [in^2]
        dlocal = [];  %Local Cross Sections of Channels
        characteristic_length = 10; %Longest length of Channels 
        chamber_length =;
        A_t =  m_dot * c_star / P_c / g;    % throat area [in^2]
        D_t = R_t * 2; % diameter at nozzle throat [in]
        A_local = pi() * (R_t)^2; % local cross sectional areas of engine
    
    
    % working fluid properties
        liq_mu = 10; %Viscosity
        Cp = 100; %Specific heat at constant pressure
        k_fluid = 100; %Thermal Conductivity of Coolant
        Tl = 293;
    
    % engine properties
        Pi = 20; %Inflow Pressure
        Ti = 200; %Inflow Temp
        Tgas = []; %Gas temperature from 1-D 
        mdotf = 30; %Coolant/fuel mass flow
        chamber_pressure = 100; %Chamber Pressure
        Cstar = 1000; %Characteristic Velocity
        throat_radius_curvature = 10; %Throat radius Curvature

    % wall material properties
        kw = 100; %Thermal Conductivity of Wall
        E = 10; %Youngs modulus of elaticisty of the tube wall material
        a = 10; %Thermal expansion coefficient of the tube wall material
        v = 24; %Poisson ration of the tube wall material

    % define constants
        g = 9.81; % gravitational accelaration [m/s]

%% Calculate minimum wall thickness
% determine optimal wall thickness based on throat conditions (because this is the most extreme location in the engine)

    %% CEA Inputs
    % initial CEA run to find throat conditions
    [~, ~, ~, M, gamma, P, T, rho, mu, Pr_gas, Mw, k, son, cp] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 1, 1);
    
    % perform initial heat transfer and structural calculations to find wall thickness
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

%% Begin iterative sizing
% start at manifold and step up along chamber length

% engine_pos = final_pos
% for engine_pos ~= final_pos
    % make T_wg guess (gas-side wall temperature)

    % calculate h_g (film coefficient) from q_dot = h_g * (T_r - T_wg_guess) - T_r is the "gas temperature" axially along the chamber from CEA

    %% Initialization
    stepnum = 100;
    T_wg_i = 1000; %Initial Guess of Gas side Temperature
    gasheattransfer = 1000;
    liqheatransfer = 0;
    
    
    
    
    
    
    %Step 2: Calculate Number of channels and channel mass flow rate
    numchannels = pi * (D_t + 0.8 * (D_t + 2 * t_w)) / (D_t + 2 * t_w); %Number of Channels (EQ 6.30) (Change the coefficient later)
    mdotchan = mdotf / numchannels; %Mass flow of channel (EQ 6.31)
    %Step 3: Begin Stepping Down tube/channel
    for i = [1:stepsize]
        [~, ~, ~, ~, gamma, P_gas, T_gas, density, mu_gas, Pr_gas, Mw, k, son, cp] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, CEA_input_name, 1, 0);
        r = Pr_gas ^ (1 / 3); % recovery factor - biased towards larger engines, very small engines should use Pr^.5 (Heister 196).
        T_r = T_gas * (1 + r * (gamma - 1) / 2 * M ^ 2) / (1 * (gamma - 1) / 2 * M ^ 2); % recovery temp (adiabatic wall temp) - corrects for compressible boundry layers (Huzel & Huang 85).
        liq_mu = exp(3.402 + 0.0132 * Pl(i) + (957.3 + 3.090 * liq_pressure -0.0542 *liq_pressure ^2) / (Tl - 57.35)); % J. Chem. Eng. Data 2022, 67, 9, 2242–2256
        Pr_liquid = liq_mu * Cp / kc;
        while liqheattransfer < gasheattransfer * upperbound || liqheattransfer > gasheattransfer * lowerbound
            %Step 5: Calculate Gas film coefficient and heat transfer 
            sigma = (.5 * T_wg / T_c * (1 + (gamma - 1) / 2 * M ^ 2) + .5) ^ -.68 * (1 + (gamma - 1) / 2 * M ^ 2) ^ -.12; % film coefficient correction factor (Huzel & Huang 86).
            h_g = (.026 / D_t ^ .2) * (mu ^ .2 * cp / Pr ^ .6) * (P_c * g / c_star) ^ .8 * (D_t / radius_throat) ^ .1 * (A_t / A_t) ^ .9 * sigma; % film coefficient - bartz equation (Huzel & Huang 86).
            gasheattransfer = h_g * (T_r - Twgi);  %Gas Heat Transfer (EQ 6.16)
    
            %Step 6: Calculate Liquid Wall Temperature from conduction 
            Twl(i) = gaswalltemp - gasheattrans * t_w / thermconduct; %Liquid Wall Temp (EQ 6.29)
    
            %Step 7: Compute Liquid Film Coefficient
            hl = (a * Re_liquid^m * Pr_liquid^n * (visfree / viscwall)^b) * k / L;    %Liquid Film Coefficient (EQ 6.19) 
    
            %Step 8: Compute Heat Flux. Compare it to step 5
            liqheatransfer = hl * (Twl(i) - Tl(i)); %Liquid Heat Transfer (EQ  6.29)
    
            %Step 4: Guess Gass Wall Temperature
    
            twg(i) = Tgas(i) - liqheatransfer / hg; %Guess Gass Wall temp using liquid heat transfer value
    
            %Step 9: Run loop until step 5/8 get same value for heat flux
        end
        %Step 10: Obtain new liquid temperature/pressure
        Tl(i + 1) = T(i) + 1 / (mdotchan * Cp)* gasheattransfer * chamber_length/stepnum * D_c/numchannels;  %Heister Eq. 6.39 (pg 212) wall_length/stepnum term needs to be fixed % IS IT D_C OR D_T
        Pl(i + 1) = Pl(i) - cf(i) * (chamber_length / (stepnum *chan_diam)) * 2 * liq_density * liq_velocity^2; %Channel Height???
    
        %Step 12: Structural Analysis Checks
        St = .5*(Pl- P_gas)*(width/wallthick)^2 + (E*a*gasheattransfer*wallthick)/(2*(1-v)kw); %Combined tangential stresses: Heister Eq 6.33 page 207
        Sl = E*a*(Twl-Tl); %Longtudinal thermal stress (Huzel and Huang, EQ 2-28, pg 92) The temperatures used here may not right for determining the delta T in this equation.
        Sc = 4*Et*Ec*t_w/((((Et)^(1/2))*((Ec)^(1/2))^2)*(3*(1-v^2)*tube_radius)); %Critical Stress Buckling (Huzel and Huang, Eq 4-29, pg 93)




end
