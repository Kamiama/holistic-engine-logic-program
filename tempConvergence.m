%CEA Inputs
fuel = 'Jet-A(L)';
oxidizer = 'O2(L)';
oxidizer_temp = 90.19; %K
P_c = 300; %psi
P_e = 101325;
OF = 2;
g = 9.81; %ft/s^2
fuel_pressure = 2000000; %liquid Pressure Pa
fuel_temp = 293;
temp_wall_gas = 2000;
D_t = 1; %diameter of the throat
[~, ~, ~, ~, gamma, P_gas, temp_chamber, density, mu_gas, Pr_gas, Mw, k, M, Cp] ...
        = RunCEA(P_c, P_e, fuel, 0, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 'test', 1, 0);
%Heat transfer Inputs
mass_flow = 0.453592; %kg/s
Area_throat = 0.001122578;%m^2
liquidthermconduct = 0.145; %W/m-K
metalthermconduct = 16.1; %W/m-K
liq_mu = exp(3.402 + 0.0132 * fuel_pressure + (957.3 + 3.090 * fuel_pressure -0.0542 * fuel_pressure ^2) / (fuel_temp - 57.35)); % J. Chem. Eng. Data 2022, 67, 9, 2242–2256
Pr_liquid = liq_mu * Cp / liquidthermconduct;
c_star = P_c * Area_throat / mass_flow; %Characteristic Velocity
t_w = .2; %wall thickness
a = .023;
m = .8;
n = .4;
liqheattransfer = 0;
gasheattransfer = 1;
throat_radius = 0.018796; %m
radius_throat_curve = 1.5 * throat_radius;
Area_ratio = 5.09;
Pr_liquid = 1.7;
Re_liquid = 15000;
channel_width = 0.0042333164; %m
while liqheattransfer ~= gasheattransfer
    %Step 5: Calculate Gas film coefficient and heat transfer 
    sigma = (.5 .* temp_wall_gas / temp_chamber .* (1 + (gamma - 1) / 2 .* M .^ 2) + .5) .^ -.68 .* (1 + (gamma - 1) / 2 .* M .^ 2) .^ -.12; % film coefficient correction factor (Huzel & Huang 86).
    h_g = (.026 / D_t .^ .2) .* (mu_gas .^ .2 .* Cp / Pr_gas .^ .6) .* (P_c .* g / c_star) .^ .8 .* (D_t / radius_throat_curve) .^ .1 .* Area_ratio .^ .9 .* sigma; % film coefficient - bartz equation (Huzel & Huang 86).
    gasheattransfer = h_g .* (temp_chamber - temp_wall_gas);  %Gas Heat Transfer (EQ 6.16)

    %Step 6: Calculate Liquid Wall Temperature from conduction 
    Twl = temp_wall_gas - gasheattransfer .* t_w / metalthermconduct; %Liquid Wall Temp (EQ 6.29)

    %Step 7: Compute Liquid Film Coefficient
    hl = (.023 * Re_liquid^.8 * Pr_liquid^n * (Twl / fuel_temp)^(-.3)) * liquidthermconduct / channel_width;    %Liquid Film Coefficient (EQ 6.19) 

    %Step 8: Compute Heat Flux. Compare it to step 5
    liqheatransfer = hl .* (Twl - fuel_temp); %Liquid Heat Transfer (EQ  6.29)

    %Step 4: Guess Gass Wall Temperature

    temp_wall_gas = temp_chamber - liqheatransfer ./ h_g; %Guess Gass Wall temp using liquid heat transfer value

    %Step 9: Run loop until step 5/8 get same value for heat flux
    fprintf('%.2f\n', temp_wall_gas)
end