%% Geometric Mapping - MZI Project Validation (Updated for Taper)
% This script traces the waveguide centers and widths to verify the geometry.

%% 1. Input Parameters (Matching current BPM script)
R = 5e-6;                   % Curvature radius [m]
d = 1.9735e-6;              % Maximum center-to-center distance [m]
core_width = 1.413e-6;      % Initial Waveguide width [m]
dz = 0.01e-6;               % High resolution for geometry validation

% Calculated curved section longitudinal length (h)
h = sqrt((d*R)/2 - d^2/16);

% Segment longitudinal lengths
Lz_seg1 = 3e-6;
Lz_seg4 = 16e-6 - 4*h - 6e-6;
Lz_seg7 = 3e-6;  

%% 2. Array Initialization
Z_total = [];
X_right = [];
W_total = []; % New array to track dynamic width
current_Z_offset = 0;

%% 3. Segment Logic (Tracing the waveguide center AND width)
% Seg 1: Initial Straight Section (3 um)
z = 0:dz:Lz_seg1;
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, zeros(size(z))];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 2: Opening Curve 1 (Concave Up)
z = 0:dz:h;
x = R * (1 - cos(asin(z / R)));
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, x];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 3: Opening Curve 2 (Concave Down)
z = 0:dz:h;
x = (d/2 - 2*R) + R * (1 - cos(asin((z - h) / R) - pi));
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, x];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 4: Straight Offset Section (100 um)
z = 0:dz:Lz_seg4;
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, ones(size(z)) * d/2];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 5: Return Curve 1 (Concave Up)
z = 0:dz:h;
x = d/2 - R * (1 - cos(asin(z / R)));
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, x];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 6: Return Curve 2 (Concave Down)
z = 0:dz:h;
x = R * (1 - cos(asin((z - h) / R)));
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, x];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 7: Final Straight Section (3 um)
z = 0:dz:Lz_seg7;
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, zeros(size(z))];
W_total = [W_total, ones(size(z)) * core_width];

%% 4. Plotting and Visualization
f = figure('Color','w', 'Name', 'MZI Geometry Validation', 'Position', [100, 100, 1000, 400]);   
ax = gca;
ax.Color = 'w';                
ax.XColor = 'k';               
ax.YColor = 'k';
set(gcf, 'InvertHardcopy', 'off');   

% Draw the waveguide core areas (using dynamic W_total)
fill([Z_total, fliplr(Z_total)]*1e6, [(X_right-W_total/2), fliplr(X_right+W_total/2)]*1e6, ...
     [0.8 0.8 1], 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Right Arm Core'); hold on;
     
fill([Z_total, fliplr(Z_total)]*1e6, [(-X_right-W_total/2), fliplr(-X_right+W_total/2)]*1e6, ...
     [1 0.8 0.8], 'EdgeColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Left Arm Core');

% Draw centerlines
plot(Z_total*1e6, X_right*1e6, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(Z_total*1e6, -X_right*1e6, 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Formatting
axis equal;
grid on;
xlabel('Longitudinal Position Z (\mu m)', 'Color', 'k');
ylabel('Transversal Position X (\mu m)', 'Color', 'k');
title('MZI Device Mask Layout - Tapered Design (Top View)', 'Color', 'k');

% Legend
lgd = legend('show', 'Location', 'northeastoutside');
lgd.Box = 'off';
lgd.Color = 'none';
lgd.TextColor = 'k';

% Safety Check: Display lengths in command window
fprintf('--- Geometry Validation ---\n');
fprintf('Calculated h length (per curve): %.4f um\n', h*1e6);
fprintf('Total Curved Length (4xh): %.4f um\n', 4*h*1e6);
fprintf('Total Device Length: %.2f micrometers\n', Z_total(end)*1e6);