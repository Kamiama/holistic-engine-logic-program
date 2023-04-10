%% HELP Engine Contour Generator
% Author: Kamon Blong (kamon.blong@gmail.com)
% Collaborators: Andrej Damjanov, Chris Bermack, Jackson Ferry-Zamora
% First Created: 7/17/2022
% Last Updated: 10/23/2022

function [x_total, r_total, L_c, L_total] = engineContour(geometry_type, bell_pct, R_t, exp_ratio, con_ratio, conv_angle, conical_half_angle, L_crude_throat, L_star, bar_size)

%{ 
Description: 

Inputs:
- 

Outputs: 
- 
%}

% interpolate initial and exit rao angles
[theta_i, theta_e] = raoAngleInterpolation(exp_ratio); 

% set conical half angle
if geometry_type == "conical"
    theta_i = conical_half_angle;
end

% convert angle inputs to radians
theta_i = theta_i * pi / 180;          % nozzle expansion (initial) angle [rad]
theta_e = theta_e * pi / 180;          % nozzle exit angle [rad] 
conv_angle = conv_angle * pi /180;     % convergence angle [rad]

% set default converging/diverging fillets (do not change unless you know what you're doing)
R_converging_fillet = 1.5;             % radius modifier for converging fillet [in]
R_diverging_fillet = .382;             % radius modifier for diverging fillet [in]
R_chamber_fillet = .75;                 % radius modifier for chamber fillet [in]

% perform intermediate geometry calculations
A_t = (R_t) .^ 2 .* pi;                % throat area [in^2]
A_c = A_t * con_ratio;                 % chamber area [in^2]
R_c = sqrt (A_c / pi);                 % chamber radius [in]
R_e = sqrt(exp_ratio) * R_t;           % exit radius [in]
A_e = (R_e).^2.*pi;                    % exit area [in^2]
L_nozzle = (R_e - R_t) / tan(pi / 12); % nozzle length [in] (Sutton 77).

%% Start generating chamber and throat geometry, working axially backwards from the end of the diverging section

% generate curve for conical and bell geometry
if geometry_type == "conical" || geometry_type == "bell"
    % diverging section - generate a arc .382 times the throat radius (Sutton 77).
    nozzle_pos = linspace(-pi / 2, theta_i - pi / 2); % arc angle is between -pi/2 and convergance angle
    x_diverging = R_diverging_fillet * R_t * cos(nozzle_pos);
    y_diverging = R_diverging_fillet * R_t * sin(nozzle_pos) + R_diverging_fillet * R_t + R_t;
    
    % throat converging section - generate a arc 1.5 times greater than throat radius (Sutton 77).
    nozzle_pos = linspace(-conv_angle - pi/2, -pi / 2); % arc angle is between convergence angle and -pi/2
    x_converging_throat = R_converging_fillet * R_t * cos(nozzle_pos);
    y_converging_throat = R_converging_fillet * R_t * sin(nozzle_pos) + R_converging_fillet * R_t + R_t;
    
    % chamber converging section - generate an arc to fillet between the chamber wall and throat converging section
    % this part isn't very important for performance, any decisions here should be made based off of size constraints
    nozzle_pos = linspace(pi / 2, pi / 2 + conv_angle); % arc angle is between convergence angle and -pi/2
    y_converging_chamber = R_chamber_fillet * R_t * sin(nozzle_pos) - R_chamber_fillet * R_t + R_c;
    x_converging_chamber = -(R_chamber_fillet * R_t * cos(nozzle_pos)) - (min(y_converging_chamber) - max(y_converging_throat)) / tan(conv_angle) + min(x_converging_throat) + R_chamber_fillet * R_t * cos(pi / 2 + conv_angle);
    
    % linear converging section - generate a line between the end points of the chamber and throat converging sections
    x_converging_linear = [max(x_converging_chamber), min(x_converging_throat)];
    y_converging_linear = [min(y_converging_chamber), max(y_converging_throat)];
    
    % combine throat contour into one geometry
    x_throat = [x_converging_chamber, x_converging_linear, x_converging_throat];
    y_throat = [y_converging_chamber, y_converging_linear, y_converging_throat];
    
    % calculate chamber volume & length
    vol_chamber = L_star * A_t;                                                                    % volume of entire chamber [in^3]
    L_converging = -min(x_throat);                                                                 % length of converging section [in]

    vol_converging_chamber = pi * trapz(fliplr((double(x_converging_chamber)) .^ 2), fliplr((double(y_converging_chamber)) .^ 2));
    vol_converging_linear = pi / 3 * (max(x_converging_linear) - min(x_converging_linear)) * (R_c ^ 2 + R_c*R_t + R_t ^ 2);
    vol_converging_throat = pi * trapz(fliplr((double(x_converging_throat)) .^ 2), fliplr((double(y_converging_throat)) .^ 2));
    
    vol_converging = vol_converging_chamber + vol_converging_linear + vol_converging_throat;       % volume of converging section [in^3]
    vol_cylindrical = vol_chamber - vol_converging;                                                % volume of cylindrical section [in^3]
    L_c = vol_cylindrical / (pi * R_c ^ 2);                                                        % chamber length [in]
    
    % straight chamber section - generate a line between injector face and chamber converging arc
    x_chamber = [min(x_throat) - L_c, min(x_throat)];
    y_chamber = [R_c, R_c];

% generate curve for crude conical geometry
elseif geometry_type == "crude"
    % straight throat section
    x_throat_linear = [0, L_crude_throat];
    y_throat_linear = [R_t, R_t];
    
    % linear converging section - generate a line between the end points of the chamber and throat converging sections
    x_converging_linear = [-(R_c - R_t) / tan(conv_angle), 0];
    y_converging_linear = [R_c, R_t];
    
    % consolidate converging and straight throat section to fit other combination format
    x_throat = [x_converging_linear, x_throat_linear];
    y_throat = [y_converging_linear, y_throat_linear];

    % calculate chamber volume & length
    vol_chamber = L_star * A_t; % volume of entire chamber [in^3]
    L_converging = -min(x_converging_linear); % length of converging section [in]
    vol_converging = pi / 3 * -min(x_converging_linear) * (R_c ^ 2 + R_c * R_t + R_t ^ 2); % volume of converging section [in^3]
    vol_cylindrical = vol_chamber - vol_converging; % volume of cylindrical section [in^3]
    L_c = vol_cylindrical / (pi * R_c ^ 2); % chamber length [in]
    
    % straight chamber section - generate a line between injector face and chamber converging arc
    x_chamber = [min(x_converging_linear) - L_c, min(x_converging_linear)];
    y_chamber = [R_c, R_c];
end

%% Generate nozzle contour

% generate curve for conical geometry
if geometry_type == "conical"
    % diverging nozzle section
    x_nozzle = [max(x_diverging), L_nozzle + max(x_diverging)];
    y_nozzle = [max(y_diverging), R_e];

% generate curve for bell geometry
elseif geometry_type == "bell"
    L_nozzle = bell_pct * L_nozzle; % truncate nozzle length for bell contour
    
    % bell section - generate a bezier curve (http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf)
        
    % calculate N, E, and Q coordinates for bezier curve equation
    x_N = 0.382 * R_t * cos(theta_i - pi / 2);
    y_N = 0.382 * R_t * sin(theta_i - pi / 2) + 0.382 * R_t + R_t;
    
    x_E = bell_pct * (((sqrt(exp_ratio) - 1) * R_t) / tan(pi / 12));
    y_E = R_e;
    
    m1 = tan(theta_i);
    m2 = tan(theta_e);
    c1 = y_N - m1 * x_N;
    c2 = y_E - m2 * x_E;
    x_Q = (c2 - c1) / (m1 - m2);
    y_Q = (m1 * c2 - m2 * c1) / (m1 - m2);

    % solve bezier curve for nozzle contour
    nozzle_pos = linspace(0, 1); % nozzle section made up of 100 points
    x_nozzle = (1 - nozzle_pos) .^ 2 .* x_N + 2 .* (1 - nozzle_pos) .* nozzle_pos .* x_Q + nozzle_pos .^ 2 .* x_E;
    y_nozzle = (1 - nozzle_pos) .^ 2 .* y_N + 2 .* (1 - nozzle_pos) .* nozzle_pos .* y_Q + nozzle_pos .^ 2 .* y_E;

% generate curve for crude conical geometry
elseif geometry_type == "crude"
    % diverging nozzle section
    x_nozzle = [L_crude_throat, L_nozzle + L_crude_throat];
    y_nozzle = [R_t, R_e];
end

%% Combine nozzle and chamber contour and output results

% combine all parts
x_total = [x_chamber, x_throat, x_nozzle];
r_total = [y_chamber, y_throat, y_nozzle];

% remove repeated indices
[unused, unique_index] = unique(x_total);
unique_index_sorted = sortrows(unique_index);
x_total = x_total(unique_index_sorted);
r_total = r_total(unique_index_sorted);

% send to .txt file for solidworks export
dlmwrite('nozzle_contour.txt', [x_total.', r_total.', zeros(length(x_total), 1)], 'delimiter', '\t', 'precision', 5);

% find total length of engine from injector plate to end of nozzle
L_total = x_total(end) - x_total(1);

% print important outputs
fprintf('\n---------- Contour Outputs ----------\n')
disp(['                     Geometry Type: ' num2str(geometry_type)]);
disp(['               Chamber radius (in): ' num2str(R_c)]);
disp(['                Throat radius (in): ' num2str(R_t)]);
disp(['                  Exit radius (in): ' num2str(R_e)]);
disp(['                 Total Length (in): ' num2str(L_total)]);
disp(['               Chamber Length (in): ' num2str(L_c)]);
disp(['            Converging Length (in): ' num2str(L_converging)]);
disp(['                Nozzle Length (in): ' num2str(L_nozzle)]);
disp(['               Chamber area (in^2): ' num2str(A_c)]);
disp(['                Throat area (in^2): ' num2str(A_t)]);
disp(['                  Exit area (in^2): ' num2str(A_e)]);
disp(['             Chamber Volume (in^3): ' num2str(vol_chamber)]);
disp(['          Converging Volume (in^3): ' num2str(vol_converging)]);
disp(['                           L* (in): ' num2str(L_star)]);
disp(['                   Expansion Ratio: ' num2str(exp_ratio)]);
disp(['                 Contraction Ratio: ' num2str(con_ratio)]);
disp(['           Convergence Angle (deg): ' num2str(conv_angle*180/pi)]);
disp(['                     Theta i (deg): ' num2str(theta_i*180/pi)]);
disp(['                     Theta e (deg): ' num2str(theta_e*180/pi)]);
fprintf('---------- Contour Outputs ----------\n')

%% Display Results
figure('Name', 'Engine Contour Plot')
hold on
grid on

% plot contour
plot(x_total, r_total, 'b');
plot(x_total, -r_total, 'b');

% plot centerline
line([-100 100],[0 0], 'Color', 'black', 'LineStyle', '--')

% plot round bar stock
if bar_size ~= 0
    line([min(x_total) max(x_total)],[bar_size/2 bar_size/2], 'Color', 'blue')
    line([min(x_total) max(x_total)],[-bar_size/2 -bar_size/2], 'Color', 'blue')
end

% labels!!!!!!
title('Engine Contour Plot')
xlabel('Inches X')
ylabel('Inches Y')
axis([x_total(1)-L_total*.1 x_total(end)+L_total*.1 -L_total*.6 L_total*.6])
hold off