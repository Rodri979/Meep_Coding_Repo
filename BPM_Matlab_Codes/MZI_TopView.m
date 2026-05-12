%% OPTIMIZED MZI - Geometric Mapping & Validation
% This script traces the waveguide centers and widths to verify the 
% 8-segment optimized geometry, including the taper.

%% 1. Input Parameters (Matching current BPM script)
R = 41.5e-6;                % [m] OPTIMIZED Curvature radius
d = 1.9735e-6;              % [m] Maximum center-to-center distance
core_width = 1.1413e-6;     % [m] Initial Waveguide width
taper_width_out = 0.3e-6;   % [m] Final Taper width
dz = 0.01e-6;               % [m] High resolution for geometry validation

% Calculated curved section longitudinal length (h)
h = sqrt((d*R)/2 - d^2/16);

% Segment longitudinal lengths
Lz_seg1 = 3e-6;
Lz_seg4 = 1.337e-6;         % OPTIMIZED intermediate straight segment
Lz_seg7 = 20e-6;            % Taper length
Lz_seg8 = 3e-6;             % Final measurement straight segment

%% 2. Array Initialization
Z_total = [];
X_right = [];
W_total = []; % Tracks dynamic width for the taper
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

% Seg 4: Straight Offset Section (Arms separated by d)
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

% Seg 6: Return Curve 2 (Concave Down - Combiner closes)
z = 0:dz:h;
x = R * (1 - cos(asin((z - h) / R)));
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, x];
W_total = [W_total, ones(size(z)) * core_width];
current_Z_offset = Z_total(end);

% Seg 7: Taper / Funil (Linear width reduction)
z = 0:dz:Lz_seg7;
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, zeros(size(z))]; % Merged at center
% Linear reduction equation matching the BPM simulation:
w_taper = core_width - (core_width - taper_width_out) * (z / Lz_seg7);
W_total = [W_total, w_taper];
current_Z_offset = Z_total(end);

% Seg 8: Final Straight Measurement Section
z = 0:dz:Lz_seg8;
Z_total = [Z_total, z + current_Z_offset];
X_right = [X_right, zeros(size(z))];
W_total = [W_total, ones(size(z)) * taper_width_out];

%% 4. Plotting and Visualization
f = figure('Color','w', 'Name', 'Optimized MZI Geometry', 'Position', [100, 100, 1200, 400]);   
ax = gca;
ax.Color = 'w';                
ax.XColor = 'k';               
ax.YColor = 'k';
set(gcf, 'InvertHardcopy', 'off');   

% Draw the waveguide core areas
fill([Z_total, fliplr(Z_total)]*1e6, [(X_right-W_total/2), fliplr(X_right+W_total/2)]*1e6, ...
     [0.2 0.6 1], 'EdgeColor', 'b', 'FaceAlpha', 0.4, 'DisplayName', 'Right Arm Core'); hold on;
     
fill([Z_total, fliplr(Z_total)]*1e6, [(-X_right-W_total/2), fliplr(-X_right+W_total/2)]*1e6, ...
     [1 0.4 0.4], 'EdgeColor', 'r', 'FaceAlpha', 0.4, 'DisplayName', 'Left Arm Core');

% Draw centerlines
plot(Z_total*1e6, X_right*1e6, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(Z_total*1e6, -X_right*1e6, 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Formatting
axis equal;
grid on;
xlabel('Longitudinal Position Z (\mu m)', 'Color', 'k', 'FontWeight', 'bold');
ylabel('Transversal Position X (\mu m)', 'Color', 'k', 'FontWeight', 'bold');
title('Optimized MZI Device Mask Layout (Top View)', 'Color', 'k', 'FontSize', 12);

% Legend
lgd = legend('show', 'Location', 'best');
lgd.Box = 'off';
lgd.Color = 'none';
lgd.TextColor = 'k';

% Safety Check: Display lengths in command window
fprintf('\n--- Optimized Geometry Validation ---\n');
fprintf('Calculated h length (per curve): %.4f um\n', h*1e6);
fprintf('Total Curved Length (4xh): %.4f um\n', 4*h*1e6);
fprintf('Segment 4 Length: %.4f um\n', Lz_seg4*1e6);
fprintf('Taper Length: %.4f um\n', Lz_seg7*1e6);
fprintf('Total Device Length: %.2f micrometers\n', Z_total(end)*1e6);