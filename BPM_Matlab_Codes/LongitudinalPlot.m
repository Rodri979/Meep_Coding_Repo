%% INTERACTIVE LONGITUDINAL VIEWER (X-Z PLANE)
%{
- This is a longitudinal viewer for the electric field X-Z planes. In its
user interface, you can select three types of plots:
    - real part of electric field (Re(E));
    - absolute intensity (|E|^2);
    - phase;
- You can also set the Y value of the plane to plot.

INSTRUCTIONS:
- This script must be executed AFTER the FD_BPM(P) solver finishes.
- You need to set P.storeE3D to true before the previous simulation
starts. If not set, run the simulation again.
- The 'P.updates' attribute sets how many points of the e-field are stored
in the 3D maps along the Z-direction and influences directly in the plot. 
Change it to adjust the resolution.
- Do NOT run 'clear' before running this script, as it uses existant 
workspace variables. If ran, run the simulation again.
%}

% Change the name of the model below if not P (BPMmatlab.model object 
% as initialized before). 
BPM_model = P;

% 1. Concatenate data from all segments into 3D blocks
E_total = cat(3, BPM_model.E3D{:});
n_total = cat(3, BPM_model.n3D{:});

% Convert axis vectors to micrometers for readability
z_coords = BPM_model.z * 1e6; 
x_coords = BPM_model.x * 1e6; 
y_coords = BPM_model.y * 1e6; 

% 2. Find the central Y index (usually near Y = 0)
[~, current_y_idx] = min(abs(y_coords)); 

% 3. Create Colormaps
% Default colormaps for Intensity and Phase
cmap_intensity = jet(256);
cmap_phase = hsv(256);

% Custom Colormap: Blue (-) -> Gray (0) -> Red (+) for the Real Part
c_blue = [0 0 1];
c_gray = [0.8 0.8 0.8];
c_red  = [1 0 0];
n_colors = 128;

% Interpolate between colors
cmap_blue_gray = [linspace(c_blue(1), c_gray(1), n_colors)', ...
                  linspace(c_blue(2), c_gray(2), n_colors)', ...
                  linspace(c_blue(3), c_gray(3), n_colors)'];
cmap_gray_red  = [linspace(c_gray(1), c_red(1), n_colors)', ...
                  linspace(c_gray(2), c_red(2), n_colors)', ...
                  linspace(c_gray(3), c_red(3), n_colors)'];
              
% Combine into a single divergent colormap
cmap_real = [cmap_blue_gray; cmap_gray_red(2:end, :)];

% Struct to easily pass colormaps to the callback function
cmaps = struct('intensity', cmap_intensity, 'real', cmap_real, 'phase', cmap_phase);

% 4. Create the GUI Figure
fig = figure('Name', 'MZI Longitudinal Mapping', 'Position', [100, 100, 1100, 750], 'Color', 'w');

% -- Subplot 1: Refractive Index (Fixed at the top) --
ax1 = subplot(2,1,1);
h_n_img = imagesc(z_coords, x_coords, squeeze(n_total(:, current_y_idx, :)));
axis xy tight; colormap(ax1, 'bone'); colorbar;
title(ax1, ['Refractive Index Profile @ Y = ', num2str(y_coords(current_y_idx), '%.2f'), ' \mu m']);
ylabel('X (\mu m)');

% -- Subplot 2: Electric Field (Dynamic based on dropdown) --
ax2 = subplot(2,1,2);
% Starts showing Intensity by default
h_E_img = imagesc(z_coords, x_coords, squeeze(abs(E_total(:, current_y_idx, :)).^2));
axis xy tight; colormap(ax2, cmap_intensity); colorbar;
title(ax2, ['Intensity (|E|^2) @ Y = ', num2str(y_coords(current_y_idx), '%.2f'), ' \mu m']);
xlabel('Z (\mu m)'); ylabel('X (\mu m)');

% 5. User Interface (UI) Controls
% Slider for the Y-axis position (Height)
sld = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', length(y_coords), ...
    'Value', current_y_idx, ...
    'SliderStep', [1/length(y_coords) 0.1], ...
    'Position', [450 15 200 20]);

uicontrol('Style','text', 'Position',[450 38 200 20], 'String','Adjust Y Slice (Height)', ...
    'BackgroundColor', 'w', 'FontWeight', 'bold');

% Dropdown menu to choose the field visualization type
drop = uicontrol('Style', 'popupmenu', ...
    'String', {'Intensity (|E|^2)', 'Real Part (Re(E))', 'Phase (\angle E)'}, ...
    'Position', [700 15 200 20], ...
    'BackgroundColor', 'w');

uicontrol('Style','text', 'Position',[700 38 200 20], 'String','Field Visualization Type', ...
    'BackgroundColor', 'w', 'FontWeight', 'bold');

% Connect callbacks to UI components
sld.Callback = @(src, event) updatePlots(sld, drop, h_n_img, h_E_img, ax1, ax2, n_total, E_total, y_coords, cmaps);
drop.Callback = @(src, event) updatePlots(sld, drop, h_n_img, h_E_img, ax1, ax2, n_total, E_total, y_coords, cmaps);

% 6. Update Function (Callback)
function updatePlots(slider, dropdown, h_n, h_E, ax1, ax2, n_dat, E_dat, y_vals, cmaps)
    % Get current slider and dropdown values
    idx = round(slider.Value);
    plotType = dropdown.Value;
    
    % Update the Refractive Index plot (always the same metric)
    h_n.CData = squeeze(n_dat(:, idx, :));
    title(ax1, ['Refractive Index Profile @ Y = ', num2str(y_vals(idx), '%.2f'), ' \mu m']);
    
    % Extract the current X-Z slice of the complex E-field
    slice_E = squeeze(E_dat(:, idx, :));
    
    % Update the bottom plot based on the Dropdown selection
    switch plotType
        case 1 % Intensity
            h_E.CData = abs(slice_E).^2;
            colormap(ax2, cmaps.intensity);
            clim(ax2, 'auto'); % Let MATLAB define the color limits
            title(ax2, ['Intensity (|E|^2) @ Y = ', num2str(y_vals(idx), '%.2f'), ' \mu m']);
            
        case 2 % Real Part
            h_E.CData = real(slice_E);
            colormap(ax2, cmaps.real);
            % To ensure gray is exactly at zero, limits (CLim) must be symmetric
            max_val = max(abs(real(slice_E(:)))); 
            if max_val == 0; max_val = 1; end % Prevent division by zero error
            clim(ax2, [-max_val max_val]);
            title(ax2, ['Real Part (Re(E)) @ Y = ', num2str(y_vals(idx), '%.2f'), ' \mu m']);
            
        case 3 % Phase
            h_E.CData = angle(slice_E);
            colormap(ax2, cmaps.phase);
            clim(ax2, [-pi pi]); % Phase always ranges from -pi to pi
            title(ax2, ['Phase (\angle E) @ Y = ', num2str(y_vals(idx), '%.2f'), ' \mu m']);
    end
end