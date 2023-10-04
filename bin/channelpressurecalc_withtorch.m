% %Pressure Drop over Regen Channel with torch igniter

%% ____________________________________
%% Inputs 
close;
clear;
clc;
%Sim Inputs
iterations = 100000;

%Chamber Parameters
chan_num = 46;
length = 25.4 * 8.9583; %Total length (mm) [conversion * in]
chamber_length = 25.4 * 5.205; %Chamber Length (mm) [conversion * in]
converging_length = 25.4 * 1.8251; %Converging Length (mm) [conversion * in]
diverging_length = 25.4 * 1.8557; %Diverging Length (mm) [conversion * in]
k = 163; %Thermal Conductivity (W/(m K))
mdot = 2.7202; %Total Propellant Massflowrate (lbm/s)
%cf = .01; %Coefficient of Friction
roughness_table = readmatrix(pwd + "/bin/surface_roughness.xlsx",'Range','A12:E16');
e = [roughness_table(2,2), roughness_table(5,2)] .* 0.001; %Surface roughness (mm) [micrometer*conversion] [45, 90]
%Channel Geometery
t = .5; %Inner wall thickness (mm)
Rib_height = [1.5 1.4 1.5]; %Rib Height (mm) [1 min 2]
chan_width = [5.4 1.6 3.227]; % Width of Channel (mm) [1 min 2]
%torch_width = 1.778;
torch_width = 1.778;
torch_height = 3;
pt_width = 5;
pt_height = 3.5;
torch_loc = [2.2 3.2 3.7] .* 25.4;
torch_conv_length = torch_loc(2) - torch_loc(1);
torch_div_length = torch_loc(3)-torch_loc(2);

%Inlet Parameters and coolant definition
P_inlet = 500; %Absolute inlet pressure (psia)
T = 273.15; %Inlet Temperature (K)
Rmdot = 3.6762; %Relative Mass flow rate (lb/s)
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

x_chamber1 = [];
x_torch_conv = [];
x_torch_div = [];
x_chamber2 = [];
x_converging = [];
x_diverging = [];


for i = x 
    if i <= torch_loc(1)
        x_chamber1 = [x_chamber1 i];
    end
    if (torch_loc(1) < i) && (i <= torch_loc(2))
        x_torch_conv = [x_torch_conv i];
    end
    if (torch_loc(2) < i) && (i <= torch_loc(3))
        x_torch_div = [x_torch_div i];
    end
    if (torch_loc(3) < i) && (i <= chamber_length)
        x_chamber2 = [x_chamber2 i];
    end 
    if (chamber_length < i) && (i <= chamber_length + converging_length)
        x_converging = [x_converging i];
    end 
    if i > (chamber_length + converging_length)
        x_diverging = [x_diverging i];
    end 
end 

chan_width_chamber1 = ones(1,size(x_chamber1,2)).*chan_width(1); %Rib height over Chamber length 
chan_width_torch_conv = ((torch_width-chan_width(1))/(torch_conv_length)).*(x_torch_conv -x_torch_conv(1)) ... 
    + ones(1,size(x_torch_conv,2)).*chan_width(1);
chan_width_torch_div = ((chan_width(1)-torch_width)/(torch_div_length)).*(x_torch_div-x_torch_div(1))... 
    + ones(1,size(x_torch_div,2)).*torch_width;
chan_width_chamber2 = ones(1,size(x_chamber2,2)).*chan_width(1); %Rib height over Chamber length 
chan_width_converging = ((chan_width(2)-chan_width(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
    + ones(1,size(x_converging,2)).*chan_width(1);
chan_width_diverging = ((chan_width(3)-chan_width(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
    + ones(1,size(x_diverging,2)).*chan_width(2);
chan_width_x = [chan_width_chamber1 chan_width_torch_conv chan_width_torch_div chan_width_chamber2 chan_width_converging chan_width_diverging];

Rib_height_chamber1 = ones(1,size(x_chamber1,2)).*Rib_height(1); %Rib height over Chamber length 
Rib_height_torch_conv = ((torch_height-Rib_height(1))/(torch_conv_length)).*(x_torch_conv -x_torch_conv(1)) ... 
    + ones(1,size(x_torch_conv,2)).*Rib_height(1);
Rib_height_torch_div = ((Rib_height(1)-torch_height)/(torch_div_length)).*(x_torch_div-x_torch_div(1))... 
    + ones(1,size(x_torch_div,2)).*torch_height;
Rib_height_chamber2 = ones(1,size(x_chamber2,2)).*Rib_height(1); %Rib height over Chamber length 
Rib_height_converging = ((Rib_height(2)-Rib_height(1))/(converging_length)).*(x_converging -x_converging(1)) ... 
    + ones(1,size(x_converging,2)).*Rib_height(1);
Rib_height_diverging = ((Rib_height(3)-Rib_height(2))/(diverging_length)).*(x_diverging-x_diverging(1))... 
    + ones(1,size(x_diverging,2)).*Rib_height(2);
Rib_height_x = [Rib_height_chamber1 Rib_height_torch_conv Rib_height_torch_div Rib_height_chamber2 Rib_height_converging Rib_height_diverging];

Ax = (chan_width_x .* Rib_height_x); %mm^2
p_wet_x = 2.*chan_width_x + 2 .* Rib_height_x; %mm
dh_x = ((4.*(Ax))./p_wet_x); %mm

v_x = (cmdot./((Ax .* (1*10^(-6))) .* rho)) ; %m/s  

%% Pump Table
% pump_table_import = readmatrix(pwd + "\Pump1.xlsx");
% 
% pump_table(:,1) = pump_table_import(:,3);
% pump_table(:,2) = pump_table_import(:,1);


%% _____________________________
%% Find Pressure
x_to_chamber2 = [x_chamber1 x_torch_conv x_torch_div];
x_to_converging = [x_to_chamber2 x_chamber2];

P_x = zeros(1,size(x,2));

if coolantdirection == 1
    P_x(size(x,2)) = P_inlet;
    for i = 1:size(x,2)-1
        j = size(x,2)-i; 
        if (i < size(x_chamber1,2)) || (((size(x_to_chamber2,2)) <= i) && (i < size(x_to_converging,2)))
            ed = e(2)/dh_x(j+1); %90 degrees
            ninety = i
        else
            ed = e(1)/dh_x(j+1); %45 degrees
            fourty_five = i
        end
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
