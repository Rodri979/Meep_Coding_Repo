%% OPTIMIZED MACH-ZEHNDER INTERFEROMETER (MZI)
%{
This is a BPM-Matlab implementation of a MZI after the optimizations
described on the Report. The ablated section is still 1um high. A new
segment was added: the taper.

    1: initial 3um-long straight section with 1.1413um of width;

    2 and 3: curved sections with the same setup as the original MZI.
    OPTIMIZATION: arcs now have a radius of 41.5um.

    4: intermediate straight segment, consisting of two straight waveguides
    that represent the MZI's arms. It is 1.337 um long. One of them may
    have different refractive index (deltaN = 0.8384) to induce destructive
    interference;

    5 and 6: combiner section. Same arc radius as sections 2 and 3.

    7: OPTIMIZATION: taper section to force interaction between the fields. 
    Input width is 1.1413um and output is 0.3um.

    8: final 3um-long straight section with 1.1413um width and 1um height.
%}
P = BPMmatlab.model;

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true;    % If false, BPM-Matlab will leave one processor unused. 
% Useful for doing other work on the PC while simulations are running.

P.useGPU = false;       % (Default: false) Use CUDA acceleration for NVIDIA GPUs.
% Available on Windows only.

%% Visualization parameters
P.updates = 1;      % Number of times to update plot. Must be at least 1, 
% showing the final state.
P.plotEmax = 0.5;   % Max of color scale in the intensity plot, relative to
% the peak of initial intensity

% Uncomment line below and line 'finalizeVideo(P)' to save a video of 
% the simulation after it's complete.

% P.saveVideo = true;

% Uncomment line below to save a 3D map of the electric field, needed to
% use the Longitudinal Viewer.

% P.storeE3D = true;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 8e-6;       % [m] x side length of main area
P.Ly_main = 8e-6;       % [m] y side length of main area
P.Nx_main = 400;        % x resolution of main area
P.Ny_main = 400;        % y resolution of main area
P.padfactor = 1.5;      % How much absorbing padding to add on the sides 
                        % of the main area (1 means no padding, 2 means the 
                        % absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.1e-6;   % [m] z step size to aim for
P.alpha = 3e14;         % [1/m^3] "Absorption coefficient" per squared unit 
                        % length distance out from edge of main area

%% Resolution (px/m)
res = P.Nx_main/P.Lx_main;

%% Problem definition
wl = 2242;
P.lambda = wl * 1e-9; % [m] Wavelength

P.n_background = 1; % [] (may be complex) Background refractive index 
                    % (in this case, the air)

%% Splitter and combiner geometric parameters

% OPTIMIZATION 1: new curved sections' radius
R = 41.5e-6;           % [m] circular section radius
d = 1.9735e-6;      % [m] maximum distance between the two arms centers 

curved_sections_length = sqrt((d*R)/2 - d^2/16);

%% Refractive Indices specifications
core_height = 1 * 1e-6;
core_width = 1.1413 * 1e-6;

bare_e_array = readPermittivityTable('material_data/Bare.txt', wl);
ablated_e_array = readPermittivityTable('material_data/Ablated.txt', wl);
fprintf("Ablated e1 = %.6f , e2 = %.6f\n", ablated_e_array(1), ablated_e_array(2));
fprintf("Bare e1 = %.6f , e2 = %.6f\n", bare_e_array(1), bare_e_array(2));

ablated_n_tilde = permittivityToRefrIdx(ablated_e_array(1), ablated_e_array(2), 0);
bare_n_tilde = permittivityToRefrIdx(bare_e_array(1), bare_e_array(2), 1);

P.n_0 = real(ablated_n_tilde); % [] reference refractive index

sio2_n = sio2RefrIdx(P.lambda);


%% Core Indices for Power Calculation
center = P.Nx_main * P.padfactor / 2;
x_core_start = ceil(center - res*core_width/2);
x_core_end   = ceil(center + res*core_width/2);
y_core_start = ceil(center - res*core_height/2);
y_core_end   = ceil(center + res*core_height/2);

%% Segment 1
P.Lz = 3e-6; % [m] z propagation distances for this segment
P = initializeRIfromFunction(P,@calcRIsegs1and8,{core_height, core_width, ...
                             ablated_n_tilde, bare_n_tilde, sio2_n});
% Finding modes
P = findModes(P, 1, 'plotModes', false);

% Initial electric field is the fundamental mode
P.E = P.modes(1);

% E-field power decay analysis - INITIAL STATE
field_intensity = abs(P.E.field).^2;
in_total_power_sum = sum(field_intensity(:)); 

% Note: on E.field, X and Y indices are flipped! Insert column (horizontal)
% boundaries first, then row (vertical) boundaries.
core_field_intensity = abs(P.E.field(x_core_start:x_core_end, y_core_start:y_core_end)).^2;
in_core_power_sum = sum(core_field_intensity(:));

%% Segment 1
P.figTitle = 'Segment 1';

P = FD_BPM(P);

%% Segment 2 
P.figTitle = 'Segment 2';

% [m] z propagation distances for this segment
P.Lz = curved_sections_length;

% Number of slices in z-direction for this segment.
Nz_n = ceil(P.Lz / P.dz_target);

P = initializeRIfromFunction(P,@calcRIseg2,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n},...
                             Nz_n);
P = FD_BPM(P);

%% Segment 3
P.figTitle = 'Segment 3';

% [m] z propagation distances for this segment
P.Lz = curved_sections_length;

P = initializeRIfromFunction(P,@calcRIseg3,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n, d}, ...
                             Nz_n);
P = FD_BPM(P);

%% Segment 4
P.figTitle = 'Segment 4';

% [m] z propagation distances for this segment
P.Lz = 1.337e-6;

% deltaN = 0.8384 -> destructive interference. Set it up here and comment
% 'ablated_n_offset = 0;' line:

% ablated_n_offset = 0.8384;

ablated_n_offset = 0;

P = initializeRIfromFunction(P, @calcRIseg4, {core_height, core_width,...
                                 ablated_n_tilde, bare_n_tilde, sio2_n, d, ablated_n_offset});
P = FD_BPM(P);

%% Segment 5
P.figTitle = 'Segment 5';

% [m] z propagation distances for this segment
P.Lz = curved_sections_length;

P = initializeRIfromFunction(P,@calcRIseg5,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n, d}, ...
                             Nz_n);
P = FD_BPM(P);

%% Segment 6
P.figTitle = 'Segment 6';

% [m] z propagation distances for this segment
P.Lz = curved_sections_length;

P = initializeRIfromFunction(P,@calcRIseg6,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n, d}, ...
                             Nz_n);
P = FD_BPM(P);

%% Segment 7 (OPTIMIZATION: TAPER)
P.figTitle = 'Segment 7 (Taper)';

% [m] z propagation distances for this segment
P.Lz = 20e-6; 

% Number of slices in z-direction for this segment.
Nz_taper = ceil(P.Lz / P.dz_target);

taper_width_out = 0.3e-6;

P = initializeRIfromFunction(P, @calcRItaper, {core_height, core_width, ...
                             taper_width_out, ablated_n_tilde, bare_n_tilde, ...
                             sio2_n, P.Lz}, Nz_taper);
P = FD_BPM(P);

%% Segment 8
P.figTitle = 'Segment 8';

% [m] z propagation distances for this segment
P.Lz = 3e-6; 
P = initializeRIfromFunction(P, @calcRIsegs1and8, {core_height, taper_width_out, ...
                             ablated_n_tilde, bare_n_tilde, sio2_n});
P = FD_BPM(P);

% finalizeVideo(P);

x_core_start = ceil(center - res*taper_width_out/2);
x_core_end   = ceil(center + res*taper_width_out/2);

% E-field power decay analysis - FINAL STATE
field_intensity = abs(P.E.field).^2;
out_total_power_sum = sum(field_intensity(:));

core_field_intensity = abs(P.E.field(x_core_start:x_core_end, y_core_start:y_core_end)).^2;
out_core_power_sum = sum(core_field_intensity(:));

fprintf('--- E-field Power Numerical Analysis ---\n');

% (Metric no. 1) P_core_in/P_total_in
guided_in_total_in = in_core_power_sum/in_total_power_sum;
fprintf('Guided_in / Total_in ');
fprintf('(must be close to 1 for guided modes): ');
fprintf('%.8f (%.2f%%)\n', guided_in_total_in, guided_in_total_in*100)

% (Metric no. 2) P_total_out/P_total_in
total_out_total_in = P.powers(end);
fprintf('Remaining total power (solver): %.8f (%.2f%%)\n', ...
    total_out_total_in, total_out_total_in*100);

% (Metric no. 3) P_core_out/P_core_in
guided_out_guided_in = out_core_power_sum/in_core_power_sum;
fprintf('Guided_out / Guided_in: %.8f (%.2f%%)\n', guided_out_guided_in, ...
    guided_out_guided_in*100);

% (Metric no. 4) P_core_out/P_total_out
guided_out_total_out = out_core_power_sum/out_total_power_sum;
fprintf('Guided_out / Total_out: %.8f (%.2f%%)\n\n', guided_out_total_out, ...
    guided_out_total_out*100);

%% Read tables
function e = readPermittivityTable(path, x)
T = readtable(path);

T.Properties.VariableNames = {'wl', 'e1_ord', 'e1_ext', 'e2_ord', 'e2_ext'};
A = T(:,{'wl', 'e1_ext','e2_ext'});
B = A(A.wl == x, :);
e = [B.e1_ext B.e2_ext];
e = [e(1,1), e(1,2)];
end

%% USER DEFINED RI FUNCTIONS
%% RI function for segments 1 and 8
function n = calcRIsegs1and8(X,Y,n_background,nParameters)
% n may be complex

% core limits
core_height = nParameters{1};
core_width = nParameters{2};

% ablated section
ablated_n_tilde = nParameters{3};
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
n(X > -(core_width/2) & X < (core_width/2) & Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;

% metal bare
bare_n_tilde = nParameters{4};
n((X < -(core_width/2) | X > core_width/2) & Y > -(core_height/2) & Y < core_height/2) = bare_n_tilde;

% SiO2 substract
n(Y < -(core_height/2)) = nParameters{5};
end

%% RI function for segment 2
function n = calcRIseg2(X,Y,Z,n_background,nParameters)
% Start by setting all pixels to n_background 
n = n_background*ones(size(X)); 

% Define the parameters
core_height = nParameters{1};
core_width = nParameters{2};
ablated_n_tilde = nParameters{3};
bare_n_tilde = nParameters{4};
R = nParameters{5};
sio2_n = nParameters{6};

% Bare and SiO2 sections
n(Y > -(core_height/2) & Y < core_height/2) = bare_n_tilde;
n(Y < -(core_height/2)) = sio2_n;

% Function to define central X positions of ablated sections
R_central_x = R * (1 - cos(asin(Z / R)));
L_central_x = -R_central_x;

% Right arm
n(X > (R_central_x - core_width/2) & X < (R_central_x + core_width/2) & ...
Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;

% Left arm
n(X > (L_central_x - core_width/2) & X < (L_central_x + core_width/2) & ...
Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;
end

%% RI function for segment 3
function n = calcRIseg3(X,Y,Z,n_background,nParameters)
% Start by setting all pixels to n_background 
n = n_background*ones(size(X)); 

% Define the parameters
core_height = nParameters{1};
core_width = nParameters{2};
ablated_n_tilde = nParameters{3};
bare_n_tilde = nParameters{4};
R = nParameters{5};
sio2_n = nParameters{6};
d = nParameters{7};

% Auxiliar parameter h (Z-shift to start curve at X = d/4)
h = sqrt((d*R)/2 - d^2/16);

% Bare and SiO2 sections
n(Y > -(core_height/2) & Y < core_height/2) = bare_n_tilde;
n(Y < -(core_height/2)) = sio2_n;

% Function to define central X positions of ablated sections
R_central_x = (d/2 - 2*R) + R * (1 - cos(asin((Z - h) / R) - pi));
L_central_x = -R_central_x;

% Right arm
n(X > (R_central_x - core_width/2) & X < (R_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;

% Left arm
n(X > (L_central_x - core_width/2) & X < (L_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;
end

%% RI function for segment 4
function n = calcRIseg4(X,Y,n_background,nParameters)
% Start by setting all pixels to n_background 
n = n_background*ones(size(X)); 

% Define the parameters
core_height = nParameters{1};
core_width = nParameters{2};
ablated_n_tilde = nParameters{3};
bare_n_tilde = nParameters{4};
sio2_n = nParameters{5};
d = nParameters{6};
n_offset = nParameters{7};

% Bare and SiO2 sections
n(Y > -(core_height/2) & Y < core_height/2) = bare_n_tilde;
n(Y < -(core_height/2)) = sio2_n;

R_central_x = d/2;
L_central_x = -R_central_x;
% Right arm
n(X > (R_central_x - core_width/2) & X < (R_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde + n_offset;

% Left arm
n(X > (L_central_x - core_width/2) & X < (L_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;
end

%% RI function for segment 5
function n = calcRIseg5(X,Y,Z,n_background,nParameters)
% Start by setting all pixels to n_background 
n = n_background*ones(size(X)); 

% Define the parameters
core_height = nParameters{1};
core_width = nParameters{2};
ablated_n_tilde = nParameters{3};
bare_n_tilde = nParameters{4};
R = nParameters{5};
sio2_n = nParameters{6};
d = nParameters{7};

% Bare and SiO2 sections
n(Y > -(core_height/2) & Y < core_height/2) = bare_n_tilde;
n(Y < -(core_height/2)) = sio2_n;

% Function to define central X positions of ablated sections
R_central_x = d/2 - R * (1 - cos(asin(Z / R)));
L_central_x = -R_central_x;

% Right arm
n(X > (R_central_x - core_width/2) & X < (R_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;

% Left arm
n(X > (L_central_x - core_width/2) & X < (L_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;
end

%% RI function for segment 6
function n = calcRIseg6(X,Y,Z,n_background,nParameters)
% Start by setting all pixels to n_background 
n = n_background*ones(size(X)); 

% Define the parameters
core_height = nParameters{1};
core_width = nParameters{2};
ablated_n_tilde = nParameters{3};
bare_n_tilde = nParameters{4};
R = nParameters{5};
sio2_n = nParameters{6};
d = nParameters{7};

% Auxiliar parameter h (Z-shift to start curve at X = d/4)
h = sqrt((d*R)/2 - d^2/16);

% Bare and SiO2 sections
n(Y > -(core_height/2) & Y < core_height/2) = bare_n_tilde;
n(Y < -(core_height/2)) = sio2_n;

% Function to define central X positions of ablated sections
R_central_x = R * (1 - cos(asin(Z - h) / R));
L_central_x = -R_central_x;

% Right arm
n(X > (R_central_x - core_width/2) & X < (R_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;

% Left arm
n(X > (L_central_x - core_width/2) & X < (L_central_x + core_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n_tilde;
end

%% RI function for Segment 7 (Taper)
function n = calcRItaper(X, Y, Z, n_background, nParameters)
% Start by setting all pixels to n_background 
n = n_background * ones(size(X)); 

% Define the parameters
core_height = nParameters{1};
width_in = nParameters{2}; % Initial width (1.413 um)
width_out = nParameters{3}; % Final width (0.3 um)
ablated_n = nParameters{4};
bare_n = nParameters{5};
sio2_n = nParameters{6};
Lz_total = nParameters{7}; % Funnel total width

% Calculate the current width (Linear Taper)
current_width = width_in - (width_in - width_out) * (Z / Lz_total);

% Bare and SiO2 sections
n(Y > -(core_height/2) & Y < core_height/2) = bare_n;
n(Y < -(core_height/2)) = sio2_n;

% Tapered Core
n(X > -(current_width/2) & X < (current_width/2) & ...
  Y > -(core_height/2) & Y < core_height/2) = ablated_n;
end

function complex_n = permittivityToRefrIdx(e1,e2,absorption)
%{
Given the real and imaginary electric permittivity values of a material, 
calculates its complex refractive index.

Parameters:
    - e1: real electric permittivity coefficient;
    - e2: imaginary electric permittivity coefficient;
    - absorption: 0 (if no absorption considered) or 1 (if absorption
    considered)

Returns:
    - the complex refractive index, in 1x1 complex double form (if
    absorption = 1) or 1x1 double (if absorption = 0)
%}
e_mag = sqrt(e1^2 + e2^2);
n = sqrt((e_mag + e1) / 2);
k = sqrt((e_mag - e1) / 2);
complex_n = n + absorption*1i*k;
end

function n = sio2RefrIdx(x)
%{
Calculates the refractive index of extraordinary SiO2 using the Sellmeier's
dispersion formula (Radhakrishnan, 1951), for a given wavelength.

Parameters:
    - x: wavelength [m]

Returns:
    - n: refractive index of SiO2.
%}
x = x/(1e-6);
n = sqrt(1+0.665721./(1-(0.060./x).^2)+0.503511./(1-(0.106./x).^2)+0.214792./(1-(0.119./x).^2)+0.539173./(1-(8.792./x).^2)+1.807613./(1-(19.70./x).^2));
end



