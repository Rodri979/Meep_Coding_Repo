%% Material Optical Properties Plotter (Interactive)
% This script reads optical constants from TXT files and allows the user 
% to dynamically switch between Permittivity (e1, e2) and Refractive Index (n, k).

% File paths
file_ablated = 'material_data/Ablated.txt';
file_bare = 'material_data/Bare.txt';

% Call the plotting function for each material
plotMaterialData(file_ablated, 'Ablated Material');
plotMaterialData(file_bare, 'Bare Material');

%% --- HELPER FUNCTION ---
function plotMaterialData(filepath, matName)
    % 1. Read data robustly (Compatible with R2019b and older)
    rawData = importdata(filepath);
    
    if isstruct(rawData)
        data = rawData.data; % Extract only the numeric matrix
    else
        data = rawData;
    end
    
    % 2. Extract Raw Permittivity Columns
    wl   = data(:, 1); % Wavelength (nm)
    e1_o = data(:, 2); % e1, Ordinary
    e1_e = data(:, 3); % e1, Extra-ordinary
    e2_o = data(:, 4); % e2, Ordinary
    e2_e = data(:, 5); % e2, Extra-ordinary

    % 3. Calculate Refractive Index (n, k)
    % Ordinary components
    emag_o = sqrt(e1_o.^2 + e2_o.^2);
    n_o = sqrt((emag_o + e1_o) / 2);
    k_o = sqrt((emag_o - e1_o) / 2);
    
    % Extra-ordinary components
    emag_e = sqrt(e1_e.^2 + e2_e.^2);
    n_e = sqrt((emag_e + e1_e) / 2);
    k_e = sqrt((emag_e - e1_e) / 2);

    % 4. Create the Interactive Figure
    figure('Name', matName, 'Color', 'w', 'Position', [100, 100, 800, 700]);

    % -- TOP PLOT: Ordinary Properties --
    ax1 = subplot(2, 1, 1);
    h1_real = plot(wl, e1_o, 'ob', 'MarkerSize', 5, 'MarkerFaceColor', 'b'); hold on;
    h1_imag = plot(wl, e2_o, 'sr', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    grid on;
    title(ax1, [matName, ' - Ordinary Permittivity'], Color='k');
    ylabel(ax1, 'Permittivity (\epsilon)', Color='k');
    legend(ax1, '\epsilon_1 (Real Part)', '\epsilon_2 (Imaginary Part)', 'Location', 'best', Color='k');
    ax1.FontSize = 11;
    ax1.Color = [0.85 0.85 0.85];
    ax1.XColor = 'k';
    ax1.YColor = 'k';

    % -- BOTTOM PLOT: Extra-ordinary Properties --
    ax2 = subplot(2, 1, 2);
    h2_real = plot(wl, e1_e, 'ob', 'MarkerSize', 5, 'MarkerFaceColor', 'b'); hold on;
    h2_imag = plot(wl, e2_e, 'sr', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    grid on;
    title(ax2, [matName, ' - Extra-ordinary Permittivity'], Color='k');
    xlabel(ax2, 'Wavelength (nm)', Color='k');
    ylabel(ax2, 'Permittivity (\epsilon)', Color='k');
    legend(ax2, '\epsilon_1 (Real Part)', '\epsilon_2 (Imaginary Part)', 'Location', 'best', Color='k');
    ax2.FontSize = 11;
    ax2.Color = [0.85 0.85 0.85];
    ax2.XColor = 'k';
    ax2.YColor = 'k';

    % 5. UI Controls: Dropdown Menu
    % Placed at the top of the window
    drop = uicontrol('Style', 'popupmenu', ...
        'String', {'Permittivity (\epsilon_1, \epsilon_2)', 'Refractive Index (n, k)'}, ...
        'Position', [300, 660, 200, 25], ...
        'BackgroundColor', 'w', 'FontSize', 10);
        
    uicontrol('Style','text', 'Position',[300, 685, 200, 15], ...
        'String','Select Optical Property', 'BackgroundColor', 'w', 'FontWeight', 'bold');

    % Connect the callback to update the plots when the dropdown changes
    drop.Callback = @(src, event) updatePlots(drop, ax1, ax2, ...
        h1_real, h1_imag, h2_real, h2_imag, ...
        e1_o, e2_o, n_o, k_o, ...
        e1_e, e2_e, n_e, k_e, matName);
end

%% --- CALLBACK FUNCTION ---
function updatePlots(dropdown, ax1, ax2, h1_real, h1_imag, h2_real, h2_imag, ...
                     e1_o, e2_o, n_o, k_o, e1_e, e2_e, n_e, k_e, matName)
                 
    plotChoice = dropdown.Value;
    
    if plotChoice == 1
        % --- Switch to Permittivity ---
        h1_real.YData = e1_o; h1_imag.YData = e2_o;
        h2_real.YData = e1_e; h2_imag.YData = e2_e;
        
        % Update Labels and Titles
        ylabel(ax1, 'Permittivity (\epsilon)', Color='k');
        ylabel(ax2, 'Permittivity (\epsilon)', Color='k');
        title(ax1, [matName, ' - Ordinary Permittivity'], Color='k');
        title(ax2, [matName, ' - Extra-ordinary Permittivity'], Color='k');
        legend(ax1, '\epsilon_1 (Real Part)', '\epsilon_2 (Imaginary Part)', 'Location', 'best', Color='k');
        legend(ax2, '\epsilon_1 (Real Part)', '\epsilon_2 (Imaginary Part)', 'Location', 'best', Color='k');
        
    elseif plotChoice == 2
        % --- Switch to Refractive Index ---
        h1_real.YData = n_o;  h1_imag.YData = k_o;
        h2_real.YData = n_e;  h2_imag.YData = k_e;
        
        % Update Labels and Titles
        ylabel(ax1, 'Refractive Index');
        ylabel(ax2, 'Refractive Index');
        title(ax1, [matName, ' - Ordinary Refractive Index']);
        title(ax2, [matName, ' - Extra-ordinary Refractive Index']);
        legend(ax1, 'n (Real Part)', 'k (Imaginary Part)', 'Location', 'best', Color='k');
        legend(ax2, 'n (Real Part)', 'k (Imaginary Part)', 'Location', 'best', Color='k');
    end
    
    % Autoscale the axes to fit the new data perfectly
    axis(ax1, 'auto');
    axis(ax2, 'auto');
end