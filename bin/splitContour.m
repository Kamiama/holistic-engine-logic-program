%% HELP Contour Splitter
% Author: Kamon Blong (kamon.blong@gmail.com)
% First Created: 04/14/2023
% Last Updated: 04/14/2023

%{ 

Description: Takes in engine geometry with a variable distance between
    points and converts it to a contour with equal-length segments.

%}

function [x_contour_new, r_contour_new, L_seg] = splitContour(x_contour, r_contour, resolution, debug)

%% Initialize Variables

num_points = resolution; % number of equally spaced points

%% Calculations

% calculate the cumulative distance for each point
cumulative_distance = zeros(1, length(x_contour));
for i = 2:length(x_contour)
    cumulative_distance(i) = cumulative_distance(i-1) + sqrt((x_contour(i)-x_contour(i-1))^2 + (r_contour(i)-r_contour(i-1))^2);
end

% Create an interpolation function for x and r points based on cumulative distance
x_interp = interp1(cumulative_distance, x_contour, 'linear', 'pp');
r_interp = interp1(cumulative_distance, r_contour, 'linear', 'pp');

% calculate the new cumulative distances for the equally spaced points
L_seg = cumulative_distance(end) / (num_points - 1)
new_cumulative_distance = 0:L_seg:cumulative_distance(end);

% interpolate the new x and r points using the new cumulative distances
x_contour_new = ppval(x_interp, new_cumulative_distance);
r_contour_new = ppval(r_interp, new_cumulative_distance);

%% Display Outputs

if debug
    % print the equal spacing distance
    disp(['Equal spacing distance: ', num2str(L_seg)]);
    
    % plot the original and new points
    figure;
    plot(x_contour, r_contour, 'bo-', 'LineWidth', 1);
    hold on;
    plot(x_contour_new, r_contour_new, 'ro-', 'LineWidth', 1);
    legend('Original Points', 'Equally Spaced Points');
    xlabel('Inches X');
    ylabel('Inches Y');
    axis equal;
    grid on;
    title('Engine Contour with Equally Spaced Points');
end