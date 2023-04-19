%
clc;
clear;
%% make contour:
load('SSME_contour2.mat');

[throat_r,throat_rINDEX] = min(r_contour);
throat_a = throat_r^2*pi;
stop = length(r_contour);

i = throat_rINDEX;
test = 0;
while test<1
    i = i+1;
    ratio = r_contour(i)^2*pi / throat_a;
    if ratio >= 4.48
        test = 1;
    end
end

rCC = r_contour(1:i);
rN = r_contour(i:stop);
xCC = x_contour(1:i);
xN = x_contour(i:stop);




% Split for CEA (@ throat)
radius_contour_index = r_contour(1:throat_rINDEX);
radius_nozzle_index = r_contour(throat_rINDEX:stop);
x_contour_index = x_contour(1:throat_rINDEX);
x_nozzle_index = x_contour(throat_rINDEX:stop);

% Output contour
fprintf('Distance of nozzle split in meters:\n')
disp(xN(1))


figure(1)
plot(xCC,rCC,'r')
hold on
plot(xN,rN,'b')
legend('combustion chamber','nozzle')
xlabel('distance from throat (m)')
ylabel('radius from centerline (m)')
title('nozzle contour')
hold off


%% Initial Conditions:
mcc_mdot_eng = 29; %lb/s
mcc_mdot = mcc_mdot_eng * 0.453592; %kg/s
n_mdot_eng = 47; %lb/s
n_mdot = n_mdot_eng * 0.453592; %kg/s
mcc_AR = 1.5;
n_AR = 1;
mcc_N = 430;
n_N = 1080;


%% CEA cc

sub = (radius_contour_index.^2*pi)/(throat_a);
sup = (radius_nozzle_index.^2*pi)/(throat_a);

% % CEA_ROCKET_EXAMPLE: Example file for the MATLAB CEA wrapper. For in-depth
% % documentation read the headers of cea_rocket_run.m,
% % cea_rocket_run_single.m, and cea_rocket_read.m
% 
% % Change this variable to true to rerun CEA instead of using saved values
% CEA_RUN = true;
% CEA_SAVE_FILE = 'cea.mat';
% 
% % The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% % hash in C) as input. The dictionary uses MATLAB character arrays as the
% % keys, and the value data type varies by which key is used. Details of
% % each key are listed in cea_rocket_run.m
% % For example: inp('key') = value.
% inp = containers.Map;
% inp('type') = 'eq fr';              % Sets the type of CEA calculation
% inp('p') = 2994;                % Chamber pressure
% inp('p_unit') = 'psi';              % Chamber pressure units
% inp('o/f') = 6.03;               % Mixture ratio
% inp('sub') = sub;               % Subsonic area ratios
% % inp('sup') = sup;               % Supersonic ratios
% inp('fuel') = 'H2';             % Fuel name from thermo.inp
% inp('fuel_t') = 800;                % Fuel inlet temperature
% inp('ox') = 'O2(L)';              % Ox name from thermo.inp
% inp('ox_t') = 90;                  % Ox inlet temperature
% inp('file_name') = 'aae539HW3P1_readme1.inp';    % Input/output file name
% if CEA_RUN
%     data = cea_rocket_run(inp);     % Call the CEA MATLAB code
%     save(CEA_SAVE_FILE, 'data');
% else
%     load(CEA_SAVE_FILE);
% end
CEA_input_name = 'AAAAAA';
i = stop;
for subs = sub
    [c_star(i), ~, ~, mach(i), gamma(i), p(i), t(i), rho(i), mu(i), pr(i), mw(i), k(i), son(i), cp_g(i)] = RunCEA(2994, 0, 'H2', 0, 800, 'O2(L)', 90, 6.03, sub, 0, 1, 1, 0, CEA_input_name);
    i = i - 1;
end
for sups = sup
    [c_star(i), ~, ~, mach(i), gamma(i), p(i), t(i), rho(i), mu(i), pr(i), mw(i), k(i), son(i), cp_g(i)] = RunCEA(2994, 0, 'H2', 0, 800, 'O2(L)', 90, 6.03, 0, sup, 1, 1, 0, CEA_input_name);
    i = i - 1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT
% 
% data_eq = data('eq');
% data_fr = data('fr');
% 
% 
% gamma = [gamma;gamma1];
% p = squeeze(data_eq('p')) / 10^6;
% p1 = [p(1);p(1)];
% p = [p;p1];
% mach = squeeze(data_eq('mach'));
% mach1 = [mach(1);mach(1)];
% mach = [mach;mach1];
% t = squeeze(data_eq('t'));
% t1 = [t(1);t(1)];
% t = [t;t1];
% pr = squeeze(data_eq('prandtl'));
% pr1 = [pr(1);pr(1)];
% pr = [pr;pr1];
% k = squeeze(data_eq('k'));
% k1 = [k(1);k(1)];
% k = [k;k1];
% rho = squeeze(data_eq('rho'));
% rho1 = [rho(1);rho(1)];
% rho = [rho;rho1];
% mu = squeeze(data_eq('visc'));
% mu1 = [mu(1);mu(1)];
% mu = [mu;mu1]; 
% son = squeeze(data_eq('son'));
% % son1 = zeros(1,19);
% son1 = [son(1);son(1)];
% % son1(:) = son(1);
% son = [son;son1];



% %% cea N
% %  CEA_ROCKET_EXAMPLE: Example file for the MATLAB CEA wrapper. For in-depth
% % documentation read the headers of cea_rocket_run.m,
% % cea_rocket_run_single.m, and cea_rocket_read.m
% 
% % Change this variable to true to rerun CEA instead of using saved values
% CEA_RUN = true;
% CEA_SAVE_FILE = 'cea.mat';
% 
% % The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% % hash in C) as input. The dictionary uses MATLAB character arrays as the
% % keys, and the value data type varies by which key is used. Details of
% % each key are listed in cea_rocket_run.m
% % For example: inp('key') = value.
% inp = containers.Map;
% inp('type') = 'eq fr';              % Sets the type of CEA calculation
% inp('p') = 2994;                % Chamber pressure
% inp('p_unit') = 'psi';              % Chamber pressure units
% inp('o/f') = 6.03;               % Mixture ratio
% % inp('sub') = sub;               % Subsonic area ratios
% inp('sup') = sup;               % Supersonic ratios
% inp('fuel') = 'H2';             % Fuel name from thermo.inp
% inp('fuel_t') = 800;                % Fuel inlet temperature
% inp('ox') = 'O2(L)';              % Ox name from thermo.inp
% inp('ox_t') = 90;                  % Ox inlet temperature
% inp('file_name') = 'aae539HW3P1_readme2.inp';    % Input/output file name
% if CEA_RUN
%     data = cea_rocket_run(inp);     % Call the CEA MATLAB code
%     save(CEA_SAVE_FILE, 'data');
% else
%     load(CEA_SAVE_FILE);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT
% 
% data_eq = data('eq');
% data_fr = data('fr');
% 
% stop = length(x_nozzle_index);
% gamma2 = squeeze(data_eq('gammas'));
% gamma1 = gamma2(2:stop+1);
% gamma = [gamma;gamma1];
% p2 = squeeze(data_eq('p')) / 10^6;
% p1 = p2(2:stop+1);
% p = [p;p1];
% mach2 = squeeze(data_eq('mach'));
% mach1 = mach2(2:stop+1);
% mach = [mach;mach1];
% t2 = squeeze(data_eq('t'));
% t1 = t2(2:stop+1);
% t = [t;t1];
% pr2 = squeeze(data_eq('prandtl'));
% pr1 = pr2(2:stop+1);
% pr = [pr;pr1];
% k2 = squeeze(data_eq('k'));
% k1 = k2(2:stop+1);
% k = [k;k1];
% rho2 = squeeze(data_eq('rho'));
% rho1 = rho2(2:stop+1);
% rho = [rho;rho1];
% mu2 = squeeze(data_eq('visc'));
% mu1 = mu2(2:stop+1);
% mu = [mu;mu1]./2; %%%%% for testing
% son2 = squeeze(data_eq('son'));
% son1 = son2(2:stop+1);
% son = [son;son1];

%}%{
%% Loop start
geometryNUM = 2; %set to number of geometries input

mu0 = mu(1);
t0 = t(1);
mach(19) = 1;
mach(20) = 1;
v = son.*mach;
pl = zeros(1,20);
tl = zeros(1,20);
pl(34) = 5956 * 6894.76; %units = pa 
tl(34) = 51.4833; %units = k 
pl(205) = 5956 * 6894.76; %units = pa 
tl(205) = 51.4833; %units = k 
Tw_guess = 400;

%
%% N parameters
N = n_N;
AR = 1.8; %%% fix MCC
mdot = n_mdot / N;
D = rN.*2;
s = (pi() .* D ./ N);
d0 = s;
th = 0.25.*d0; %m
di = d0 - 2.*th;
W = AR .* di;
A = di .* W + pi.*(di./2).^2;
AN =A;
stop = length(p);
Km = 24; %%%%%%%%%%%%%%%%%%%%%
topI = 35;
per = 1.*W + pi.*di;
Dl = 4.*A ./ (per);
DlN = Dl;
x_contour_index = xN;    
iedit = length(p) - length(rN);
f = .005;
TW_guess = 400;
% eta = 1;
a =.05;
b=2;
eta = 1 - exp(-a.*(rN./Dl).^2)/b;

etaN=eta;
%% N loop
ii = 1;
for i =stop:-1:topI+2
    check = 0;
    

       
    while check==0
       Tam(ii) = .5*(t(i) + Tw_guess);
       Rhoam(ii) = t(i) / Tam(ii)*rho(i);  %ideal gas law
       w = log(mu0/mu(i))/log(t0/t(i));
       Muam(ii) = mu0*(Tam(ii)/t0)^w;
       re(ii) = rho(i)*v(i)*D(i-iedit) / mu(i);
       re(ii) = re(ii)/2;
       Nug(ii) = eta(i-iedit)*(.026*re(ii)^0.8*pr(i).^0.4.*(Rhoam(ii)/rho(i))^0.8.*(Muam(ii)/mu0)^0.2);
       hg(ii) = Nug(ii)*k(i)/D(i-iedit);
       r = pr(i)^(1/3);
       Tr(ii) = t(i)*(1+(gamma(i)-1)/2 * mach(i)^2*r);
       qdot1(ii) = hg(ii)*(Tr(ii)-Tw_guess);
       Twl(ii) = Tw_guess - qdot1(ii)*th(i-iedit)/Km;
       prl(ii) = py.CoolProp.CoolProp.PropsSI('PRANDTL','P',pl(i),'T',tl(i),'Hydrogen');
       rhol(ii) = py.CoolProp.CoolProp.PropsSI('D','P',pl(i),'T',tl(i),'Hydrogen');
       mul(ii) = py.CoolProp.CoolProp.PropsSI('V','P',pl(i),'T',tl(i),'Hydrogen');
       cp(ii) = py.CoolProp.CoolProp.PropsSI('C','P',pl(i),'T',tl(i),'Hydrogen');
       Kl(ii) = py.CoolProp.CoolProp.PropsSI('L','P',pl(i),'T',tl(i),'Hydrogen');
       Vl(ii) = mdot/(rhol(ii)*A(i-iedit));
       rel(ii) = rhol(ii)*Vl(ii)*Dl(i-iedit) / mul(ii);
       Nuc(ii) = 0.023*rel(ii)^0.8*prl(ii)^0.4*(Twl(ii)/tl(i))^-0.3;
       hl(ii) = Nuc(ii)*Kl(ii) / Dl(i-iedit);
       qdot2(ii) = hl(ii)*(Twl(ii) - tl(i));

       if abs(qdot1(ii) - qdot2(ii))/qdot1(ii) <= 0.01
           check = 1;
       elseif qdot1(ii) - qdot2(ii) > 0
           Tw_guess = Tw_guess+.1;
           check = 0;
       elseif qdot2(ii) - qdot1(ii) > 0
           Tw_guess = Tw_guess-.1;
           check = 0;
       else
           disp('error in nozzle loop')
       end
    end
    TwN(ii) = Tw_guess;
    Tw_guess = TwN(ii);
    if i > 37
        tl(i-1) = tl(i) + qdot2(ii)*s(i-iedit)*(x_contour_index(i-iedit)-x_contour_index(i-1-iedit))/ (mdot*cp(ii));
        f(ii) = (.0014+.125/rel(ii)^0.32)/12;
        pl(i-1) = pl(i) - f(ii)*(x_contour_index(i-iedit)-x_contour_index(i-1-iedit))/(Dl(i-iedit)-Dl(i-1-iedit)) * 2*rhol(ii)*Vl(ii)^2;
        ii = ii +1;
    end

end

TwN = flip(TwN);

% plot(x, TwN)


%}
%
ii = 185;

%% MCC parameters
th = .0008; %m
N = mcc_N;
mdot = mcc_mdot / mcc_N;
D = rCC.*2;
s = (pi().*D./N)/2;
AR = 1.5;
w = AR .*s;
A = s.*w;
stop = length(rCC);
Km = 360; %%%%%%%%%%%%%%%%%%%%%
topI = 1;
Dl = 4.*A ./ (2.*s+2.*w);
x_contour_index = xN;
TW_guess = 1200;
eta =1;
a =.05;
b=2;
eta = 1 - exp(-a.*(rCC./Dl).^2)/b;



%% MCC loop 1
c=1;
for i =stop:-1:21
    check = 0;
    

       
    while check==0
       Tam(ii) = .5*(t(i) + Tw_guess);
       Rhoam(ii) = t(i) / Tam(ii)*rho(i);  %ideal gas law
       if mu0 == mu(i)
           w = .75;

       elseif t0 == t(i)
           w = .75;
       else
           w = log(mu0/mu(i))/log(t0/t(i));
       end
       
       Muam(ii) = mu0*(Tam(ii)/t0)^w;
       if Muam(ii) < 0.000001
          Muam(ii) = 3.53-05;
       end
       re(ii) = rho(i)*v(i)*D(i ) / mu(i);
       re(ii) = re(ii)/2;
       Nug(ii) = eta(i)*(.026*re(ii)^0.8*pr(i).^0.4.*(Rhoam(ii)/rho(i))^0.8.*(Muam(ii)/mu0)^0.2);
       hg(ii) = Nug(ii)*k(i)/D(i );
       r = pr(i)^(1/3);
       Tr(ii) = t(i)*(1+(gamma(i)-1)/2 * mach(i)^2*r);
       qdot1(ii) = hg(ii)*(Tr(ii)-Tw_guess);
       Twl(ii) = Tw_guess - qdot1(ii)*th/Km;
       prl(ii) = py.CoolProp.CoolProp.PropsSI('PRANDTL','P',pl(i),'T',tl(i),'Hydrogen');
       rhol(ii) = py.CoolProp.CoolProp.PropsSI('D','P',pl(i),'T',tl(i),'Hydrogen');
       mul(ii) = py.CoolProp.CoolProp.PropsSI('V','P',pl(i),'T',tl(i),'Hydrogen');
       cp(ii) = py.CoolProp.CoolProp.PropsSI('C','P',pl(i),'T',tl(i),'Hydrogen');
       Kl(ii) = py.CoolProp.CoolProp.PropsSI('L','P',pl(i),'T',tl(i),'Hydrogen');
       Vl(ii) = mdot/(rhol(ii)*A(i ));
       rel(ii) = rhol(ii)*Vl(ii)*Dl(i ) / mul(ii);
       Nuc(ii) = 0.023*rel(ii)^0.8*prl(ii)^0.4*(Twl(ii)/tl(i))^-0.3;
       hl(ii) = Nuc(ii)*Kl(ii) / Dl(i );
       qdot2(ii) = hl(ii)*(Twl(ii) - tl(i));

       if abs(qdot1(ii) - qdot2(ii))/qdot1(ii) <= 0.01
           check = 1;
       elseif qdot1(ii) - qdot2(ii) > 0
           Tw_guess = Tw_guess+.05;
           check = 0;
       elseif qdot2(ii) - qdot1(ii) > 0
           Tw_guess = Tw_guess-.05;
           check = 0;
       else
           disp('error in MCC loop')
       end
%        disp(i)
    end
    Tw1(c) = Tw_guess;
%     Tw_guess = Tw(ii);
    if i > 1
        tl(i-1) = tl(i) + qdot2(ii)*s(i )*(x_contour_index(i )-x_contour_index(i-1 ))/ (mdot*cp(ii));
        f(ii) = (.0014+.125/rel(ii)^0.32)/18;
        pl(i-1) = pl(i) - f(ii)*(x_contour_index(i )-x_contour_index(i-1 ))/(Dl(i )-Dl(i-1 )) * 2*rhol(ii)*Vl(ii)^2;
        ii = ii +1;
    end
    c = c+1;
end

TwMCC1 = flip(Tw1);
%} 
%

%% MCC loop 2
ii = ii +3;
tl(17) = tl(i-1);
pl(17) = pl(i-1);
c=1;
for i =17:-1:topI+1
    
    check = 0;
    

       
    while check==0
       Tam(ii) = .5*(t(i) + Tw_guess);
       Rhoam(ii) = t(i) / Tam(ii)*rho(i);  %ideal gas law
       if mu0 == mu(i)
           w = .75;

       elseif t0 == t(i)
           w = .75;
       else
           w = log(mu0/mu(i))/log(t0/t(i));
       end
       
       Muam(ii) = mu0*(Tam(ii)/t0)^w;
       if Muam(ii) < 0.000001
          Muam(ii) = 3.53-05;
       end
       re(ii) = rho(i)*v(i)*D(i ) / mu(i);
       re(ii) = re(ii)/2;
       Nug(ii) = eta(i)*(.026*re(ii)^0.8*pr(i).^0.4.*(Rhoam(ii)/rho(i))^0.8.*(Muam(ii)/mu0)^0.2);
       hg(ii) = Nug(ii)*k(i)/D(i );
       r = pr(i)^(1/3);
       Tr(ii) = t(i)*(1+(gamma(i)-1)/2 * mach(i)^2*r);
       qdot1(ii) = hg(ii)*(Tr(ii)-Tw_guess);
       Twl(ii) = Tw_guess - qdot1(ii)*th/Km;
       prl(ii) = py.CoolProp.CoolProp.PropsSI('PRANDTL','P',pl(i),'T',tl(i),'Hydrogen');
       rhol(ii) = py.CoolProp.CoolProp.PropsSI('D','P',pl(i),'T',tl(i),'Hydrogen');
       mul(ii) = py.CoolProp.CoolProp.PropsSI('V','P',pl(i),'T',tl(i),'Hydrogen');
       cp(ii) = py.CoolProp.CoolProp.PropsSI('C','P',pl(i),'T',tl(i),'Hydrogen');
       Kl(ii) = py.CoolProp.CoolProp.PropsSI('L','P',pl(i),'T',tl(i),'Hydrogen');
       Vl(ii) = mdot/(rhol(ii)*A(i ));
       rel(ii) = rhol(ii)*Vl(ii)*Dl(i ) / mul(ii);
       Nuc(ii) = 0.023*rel(ii)^0.8*prl(ii)^0.4*(Twl(ii)/tl(i))^-0.3;
       hl(ii) = Nuc(ii)*Kl(ii) / Dl(i );
       qdot2(ii) = hl(ii)*(Twl(ii) - tl(i));

       if abs(qdot1(ii) - qdot2(ii))/qdot1(ii) <= 0.01
           check = 1;
       elseif qdot1(ii) - qdot2(ii) > 0
           Tw_guess = Tw_guess+.05;
           check = 0;
       elseif qdot2(ii) - qdot1(ii) > 0
           Tw_guess = Tw_guess-.05;
           check = 0;
       else
           disp('error in MCC loop')
       end
%        disp(i)
    end
    Tw2(c) = Tw_guess;
%     Tw_guess = Tw(ii);
    if i > 1
        tl(i-1) = tl(i) + qdot2(ii)*s(i )*(x_contour_index(i )-x_contour_index(i-1 ))/ (mdot*cp(ii));
        f(ii) = (.0014+.125/rel(ii)^0.32)/20;
        pl(i-1) = pl(i) - f(ii)*(x_contour_index(i )-x_contour_index(i-1 ))/(Dl(i-1 )-Dl(i )) * 2*rhol(ii)*Vl(ii)^2;
        ii = ii +1;
    end
    c = c+1;
end
%}

TwMCC2 = flip(Tw2);

TwMCC = [TwMCC2,TwMCC1];
XMCC = [xCC(2:17),xCC(21:34)];
Twl_Mcc = flip([Twl(185:198),Twl(202:217)]);
Twl_N = flip(Twl(1:169));
Tl_Mcc = ([tl(2:17),tl(21:34)]);
Tl_N = (tl(37:205));
plMcc = ([pl(2:17),pl(21:34)]);
plN = (pl(37:205));

figure(2)
plot(XMCC,TwMCC)
hold on
plot(x_contour_index,TwN)
plot(XMCC,Twl_Mcc)
plot(x_contour_index,Twl_N)
plot(XMCC,Tl_Mcc)
plot(x_contour_index,Tl_N)
plot(xCC,t(1:202))
xline(xN(1),'-','split')
ylabel('Temperature (k)')
xlabel('Nozzle distance (m)')
legend('Wall Gas temp (MCC)','Wall Gas temp (N)','Wall liquid temp (MCC)','Wall Liquid temp (N)','Liquid temp (MCC)','Liquid temp (N)')
hold off



AMcc = ([A(2:17),A(21:34)]);
figure(4)
plot(XMCC,AMcc)
hold on
plot(x_contour_index,AN)
xline(xN(1),'-','split')
ylabel('area (m^2)')
xlabel('Nozzle distance (m)')
title('cross-sectional area vs nozzle')
hold off

DlMcc = ([Dl(2:17),Dl(21:34)]);
figure(5)
plot(XMCC,DlMcc)
hold on
plot(x_contour_index,DlN)
xline(xN(1),'-','split')
ylabel('diameter (m)')
xlabel('Nozzle distance (m)')
title('hydrolic diameter vs nozzle')
hold off


rhoMcc = flip([rhol(185:198),rhol(202:217)]);
rhoN = flip(rhol(1:169));
figure(6)
plot(XMCC,rhoMcc)
hold on
plot(x_contour_index,rhoN)
xline(xN(1),'-','split')
ylabel('density (kg/m^2)')
xlabel('Nozzle distance (m)')
title('density vs nozzle')
hold off



fMcc = flip(f(185:198));
fN = flip(f(1:169));
figure(7)
plot(XMCC(1:14),fMcc)
hold on
plot(x_contour_index,fN)
xline(xN(1),'-','split')
ylabel('friction factor')
xlabel('Nozzle distance (m)')
title('friction factor vs nozzle')
hold off


vMcc = flip([Vl(185:198),Vl(202:217)]);
vN = flip(Vl(1:169));
figure(8)
plot(XMCC,vMcc)
hold on
plot(x_contour_index,vN)
xline(xN(1),'-','split')
ylabel('Coolant Velocity (m/s)')
xlabel('Nozzle distance (m)')
title('coolant velocity vs nozzle')
hold off


machMcc = ([mach(2:17);mach(21:34)]);
machN = (mach(37:205));
figure(9)
plot(XMCC,machMcc)
hold on
plot(x_contour_index,machN)
xline(xN(1),'-','split')
ylabel('mach number ()')
xlabel('Nozzle distance (m)')
title('mach vs nozzle')
hold off


qdotMcc = flip([qdot1(185:198),qdot1(202:217)]);
qdotN = flip(qdot1(1:169));
figure(10)
plot(XMCC,machMcc)
hold on
plot(x_contour_index,machN)
xline(xN(1),'-','split')
ylabel('qdot (Heat Flux)')
xlabel('Nozzle distance (m)')
title('heat flux vs nozzle')
hold off


hgMcc = flip([hg(185:198),hg(202:217)]);
hgN = flip(hg(1:169));
hlMcc = flip([hl(185:198),hl(202:217)]);
hlN = flip(hl(1:169));
figure(11)
plot(XMCC,hlMcc,'b')
hold on
plot(x_contour_index,hlN,'b')
plot(XMCC,hgMcc,'r')
plot(x_contour_index,hgN,'r')
xline(xN(1),'-','split')
ylabel('heat transfer coefficient (H)')
xlabel('Nozzle distance (m)')
legend('gas side','','coolant side')
title('h vs nozzle')
hold off

etaMcc = ([eta(2:17),eta(21:34)]);
% etaN = (eta(37:205));
figure(12)
plot(XMCC,etaMcc,'b')
hold on
plot(x_contour_index,etaN,'b')
xline(xN(1),'-','split')
ylabel('eta')
xlabel('Nozzle distance (m)')
legend('gas side','','coolant side')
title('eta vs nozzle')
hold off

