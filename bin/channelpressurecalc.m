% %Pressure Drop over Regen Channel

%% ____________________________________
%% Inputs 
close;
clear;
clc;
%Sim Inputs
iterations = 100000;

%Chamber Parameters
chan_num = 40;
length = 25.4 * 6.6175; %Total length (mm) [conversion * in]
chamber_length = 25.4 * 3.2967; %Chamber Length (mm) [conversion * in]
converging_length = 25.4 * 1.8271; %Converging Length (mm) [conversion * in]
diverging_length = 25.4 * 1.4937; %Diverging Length (mm) [conversion * in]
k = 110; %Thermal Conductivity (W/(m K))
mdot = 2.7393; %Total Propellant Massflowrate (lbm/s)
%cf = .01; %Coefficient of Friction
e = 24 * 0.001; %Surface roughness (mm) [micrometer*conversion]
%Channel Geometery
t = 1; %Inner wall thickness (mm)
Rib_height = [3.5 1.4 3.5]; %Rib Height (mm) [1 min 2]
chan_width = [5.5 1.6 4]; % Width of Channel (mm) [1 min 2]

%Inlet Parameters and coolant definition
P_inlet = 500; %Absolute inlet pressure (psia)
T = 273.15; %Inlet Temperature (K)
Rmdot = 4.38; %Relative Mass flow rate (lb/s)
rho = 999.85; %Coolant Density kg/m^3
mew = 0.001741; %Dynamic Viscosity of water at 274 K (Ns/m^2)
coolantdirection = 1; % 1: Direction opposite of hot gas flow direction
                      % 0: Direction same as hot flow gas direction

%% _____________________________________
%% Initialization
deltax = (length/iterations) ; %m
cmdot =  (Rmdot * mdot)/((2.2046)*(chan_num)); %Channel Mass Flow Rate (kg/s)

A = chan_width .* Rib_height; %Channel Cross-sectional Area (mm^2) [1 min 2]
p_wet = 2*chan_width + 2*Rib_height; %Wetted Perimeter of the pipe (mm)
dh = (4.*A)./p_wet; %hydraulic Diameter

x = linspace(0,length,iterations); %Length Vector

x_chamber = [];
x_converging = [];
x_diverging = [];


for i = x 
    if i <= chamber_length
        x_chamber = [x_chamber i];
    end 
    if (chamber_length < i) && (i <= chamber_length + converging_length)
        x_converging = [x_converging i];
    end 
    if i > (chamber_length + converging_length)
        x_diverging = [x_diverging i];
    end 
end 

chan_width_chamber = ones(1,size(x_chamber,2)).*chan_width(1); %Rib height over Chamber length 
chan_width_converging = ((chan_width(2)-chan_width(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
    + ones(1,size(x_converging,2)).*chan_width(1);
chan_width_diverging = ((chan_width(3)-chan_width(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
    + ones(1,size(x_diverging,2)).*chan_width(2);
chan_width_x = [chan_width_chamber chan_width_converging chan_width_diverging];

Rib_height_chamber = ones(1,size(x_chamber,2)).*Rib_height(1); %Rib height over Chamber length 
Rib_height_converging = ((Rib_height(2)-Rib_height(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
    + ones(1,size(x_converging,2)).*Rib_height(1);
Rib_height_diverging = ((Rib_height(3)-Rib_height(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
    + ones(1,size(x_diverging,2)).*Rib_height(2);
Rib_height_x = [Rib_height_chamber Rib_height_converging Rib_height_diverging];

Ax = (chan_width_x .* Rib_height_x); %mm^2
p_wet_x = 2.*chan_width_x + 2 .* Rib_height_x; %mm
dh_x = ((4.*(Ax))./p_wet_x); %mm

v_x = (cmdot./((Ax .* (1*10^(-6))) .* rho)) ; %m/s  

%% Pump Table
pump_table_import = readmatrix(pwd + "\Pump1.xlsx");

pump_table(:,1) = pump_table_import(:,3);
pump_table(:,2) = pump_table_import(:,1);


%% _____________________________
%% Find Pressure
P_x = zeros(1,size(x,2));

if coolantdirection == 1
    P_x(size(x,2)) = P_inlet;
    for i = 1:size(x,2)-1
        j = size(x,2)-i;  
        ed = e/dh_x(j+1);
        Re = (rho * v_x(j+1)*(dh_x(j+1) * .001))/ mew;
        f = moody(ed, Re);
        cf = f/4;
        deltaP =   (2*cf*(deltax/(dh_x(j+1))) * rho *(v_x(j+1))^(2)  ) * (1/6894.7573);
        P_x(j) =   P_x(j+1) - deltaP;
    end
else
    P_x(1) = P_inlet;
    for j = 2:iterations
        ed = e/dh_x(j-1);
        Re = (rho * v_x(j-1)*(dh_x(j-1) * .001))/ mew;
        f = moody(ed, Re);
        cf = f/4;
        deltaP =   (2*cf*(deltax/(dh_x(j-1))) * rho *(v_x(j-1))^(2)  ) * (1/6894.7573);
        P_x(j) =   P_x(j-1) - deltaP;
    end 
end

%% Iterate over pump curve 
% deltaP_total = P_x(1) - P_x(size(P_x,2));
% pmdot_gal = interp1(pump_table(:,1), pump_table(:,2), deltaP_total);
% pmdot = .06 * pmdot_gal;
% cmdot_new = pmdot / 20;
% if cmdot_new == cmdot 
%     i = true; 
% else 
%     i = false;
% end
% while i == false 
%     cmdot = cmdot_new;
%     v_x = (cmdot./((Ax .* (1*10^(-6))) .* rho)) ; %m/s  
%     P_x(1) = P_inlet;
%     for k = 2:iterations
%         deltaP =   (2*cf*(deltax/(dh_x(k-1))) * rho *(v_x(k-1))^(2)  ) * (1/6894.7573);
%         P_x(k) =   P_x(k-1) - deltaP;
%     end 
%     deltaP_total = P_x(1) - P_x(size(P_x,2));
%     cmdot_new = .06* interp1(pump_table(:,1), pump_table(:,2), deltaP_total)
% %     if cmdot_new == cmdot 
% %         i = true; 
% %     else 
% %         i = false;
% %     end
% i = true
% end 
% 





%% Plot
figure;
subplot(2,2,1);
plot(x, P_x);
title("Pressure Drop along Channels");
xlabel("distance along chamber (mm)");
ylabel("pressure (psia)");
xlim([0,x(size(x,2))]);
grid on

subplot(2,2,2);
plot(x,v_x);
title("Velocity in Channels");
xlabel("distance along chamber (mm)");
ylabel("velocity (m/s)");
xlim([0,x(size(x,2))]);
grid on

subplot(2,2,3);
plot(x,Ax);
title("Cross-sectional Area of Channels");
xlabel("distance along chamber (mm)");
ylabel("Area (mm^2)");
xlim([0,x(size(x,2))]);
grid on

subplot(2,2,4);
plot(x,dh_x);
title("Hydraulic Diameter of Channels");
xlabel("distance along chamber (mm)");
ylabel("Diameter (mm)");
xlim([0,x(size(x,2))]);
grid on
