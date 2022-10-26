%% Unit Conversions Function
% Author: Liam Schenk
% First Created: 7/10/2022
% Last Updated: 10/26/2022 by Kamon Blong

function u = convertUnits()

    % LENGTH CONVERSIONS:
    % F: feet [ft]
    % M: meters [m]
    % IN: inches [in]
    % MI: miles [mi]
    % KM: kilometers [km]
    u.F2M  = 0.3048;            % Meters per foot
    u.M2F  = 1 / u.F2M;         % Feet per meter
    u.F2IN = 12.0;              % Inches per foot
    u.IN2F = 1 / u.F2IN;        % Feet per inches
    u.M2IN  = u.M2F * u.F2IN;   % Inches per meter
    u.IN2M  = 1 / u.M2IN;       % Meters per inch
    u.MI2F = 5280;              % Feet per mile
    u.F2MI = 1 / u.MI2F;        % Miles per foot
    u.M2KM  = 1 / 1000;         % Kilometers per meter
    u.KM2M  = 1 / u.M2KM;       % Meters per kilometer
    u.M2MI  = u.M2F * u.F2MI;   % Miles per meter
    u.MI2M  = 1 / u.M2MI;       % Meters per mile
    u.F2KM = u.F2M * u.M2KM;    % Feet per kilometer
    u.KM2F = 1 / u.F2KM;        % Kilometers per foot

    % FORCE CONVERSIONS:
    % LBF: pound-force [lbf]
    % N: Newton [N]
    u.LBF2N = 4.4482216152605;  % Newtons per pound-force
    u.N2LBF = 1 / u.LBF2N;      % Pound-force per Newton

    % PRESSURE CONVERSIONS:
    % PA: Pascal [Pa]
    % PSI: pound-force per square inch [psi]
    % ATM: atmospheres [atm]
    % B: bars [Bar]
    u.PA2PSI = u.N2LBF / (u.M2IN^2); % PSI per Pascal
    u.PSI2PA = 1 / u.PA2PSI;         % Pascal per PSI
    u.ATM2B = 1.01325;               % Bars per atmosphere
    u.B2ATM = 1 / u.ATM2B;           % Atmospheres per Bar
    u.ATM2PA = 101325;               % Pascals per atmosphere
    u.PA2ATM = 1 / u.ATM2PA;         % Atmospheres per Pascal
    u.ATM2PSI = 14.69595;            % PSI per atmosphere
    u.PSI2ATM = 1 / u.ATM2PSI;       % Atmospheres per PSI
    u.B2PA = 100000;                 % Pascal per Bar
    u.PA2B = 1 / u.B2PA;             % Bars per Pascal
    u.B2PSI = 14.50377;              % PSI per Bar
    u.PSI2B = 1 / u.B2PSI;           % Bar per PSI
    u.MPA2PSI = 145.0377;            % Megapascals per PSI
    u.PSI2MPA = 1 / u.MPA2PSI;       % PSI per Megapascals

    % MASS CONVERSIONS:
    % LB: pound [lbm]
    % KG: kilogram [kg]
    % S: slugs [slug]
    % G: gram [g]
    u.LB2KG = 0.45359237;       % Kilograms per pound
    u.KG2LB = 1 / u.LB2KG;      % Pounds per kilogram
    u.KG2S = 0.0685218;         % Slugs per kilogram
    u.S2KG = 1 / u.KG2S;        % Kilograms per slug
    u.S2LB = 32.17405;          % Pounds per slug
    u.LB2S = 1 / u.S2LB;        % Slugs per pound
    u.KG2G = 1000;              % Grams per kilogram
    u.G2KG = 1 / u.KG2G;        % Kilograms per gram

    % TEMPERATURE CONVERSIONS:
    % K: Kelvin [K]
    % R: degrees Rankine [degR]
    u.K2R = 9/5;                % Degrees Rankine per Kelvin
    u.R2K = 1 / u.K2R;          % Kelvin per Degree Rankine

    % ANGLE CONVERSIONS:
    % RAD: radian [Rad]
    % DEG: degree [Deg]
    u.RAD2DEG = 180 / pi;       % Degrees per radian
    u.DEG2RAD = pi / 180;       % Radians per degree

    % VOLUME CONVERSIONS:
    % CF: cubic feet [ft^3]
    % CI: cubic inches [in^3]
    % CM: cubic meters [m^3]
    u.CF2CI = 1728;             % Cubic inches per cubic feet
    u.CI2CF = 1 / u.CF2CI;      % Cubic feet per cubic inches
    u.CM2CI = 61024;            % Cubic inches per cubic meter
    u.CI2CM = 1 / u.CM2CI;      % Cubic meters per cubic inch
    u.CM2CF = 35.315;           % Cubic feet per cubic meter
    u.CF2CM = 1 / u.CM2CF;      % Cubic meters per cubic feet

    % AREA CONVERSIONS:
    % SM: square meter [m^2]
    % SF: square foot [ft^2]
    % SI: square inch [in^2]
    u.SM2SF = 10.76391;         % Square feet per square meter
    u.SF2SM = 1 / u.SM2SF;      % Square meter per square foot
    u.SF2SI = 144;              % Square inches per square foot
    u.SI2SF = 1 / u.SF2SI;      % Square feet per square inch
    
    % ENERGY CONVERSIONS:
    % J: joules [J]
    % KJ: kilojoules [kJ]
    % WH: watt hour [Wh]
    u.J2KJ = 0.001;             % Kilojoules per joule
    u.KJ2J = 1 / u.J2KJ;        % Joules per kilojoule
    u.KJ2WH = 0.2777778;        % Watt hours per kilojoule
    u.WH2KJ = 1 / u.KJ2WH;      % Kilojoules per watt hour
    u.J2WH = 0.0002777778;      % Watt hour per joule
    u.WH2J = 1 / u.J2WH;        % Joules per watt hour

    % VELOCITY CONVERSIONS:
    % MPS: meters per second [m/s]
    % FPS: feet per second [ft/s]
    % MPH: miles per hour [mi/h]
    % KPH: kilometers per hour [km/h]
    u.MPH2KPH = 1.609344;       % Kilometers per hour per miph
    u.KPH2MPH = 1 / u.MPH2KPH;  % Miles per hour per kmph
    u.MPS2FPS = 3.281;          % Feet per second per mps
    u.FPS2MPS = 1 / u.MPS2FPS;  % Meters per second per fps
    u.MPH2FPS = 1.467;          % Feet per second per miph
    u.FPS2MPH = 1 / u.MPH2FPS;  % Miles per hour per fps

    % DENSITY CONVERSIONS:
    u.KGM32LBIN3 = 1 / 27680;
    u.LBIN32KGM3 = 1 / u.KGM32LBIN3;

end