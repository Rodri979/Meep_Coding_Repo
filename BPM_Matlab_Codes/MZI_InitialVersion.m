P = BPMmatlab.model;

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 1;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 0.5; % Max of color scale in the intensity plot, relative to the peak of initial intensity
% P.saveVideo = true;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 8e-6;        % [m] x side length of main area
P.Ly_main = 8e-6;        % [m] y side length of main area
P.Nx_main = 400;          % x resolution of main area
P.Ny_main = 400;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Resolution (px/m)
res = P.Nx_main/P.Lx_main;

%% Problem definition
wl = 2242;
P.lambda = wl * 1e-9; % [m] Wavelength

P.n_background = 1; % [] (may be complex) Background refractive index (in this case, the air)

R = 5e-6; % [m] circular section radius
d = 1.9735e-6; % [m] maximum distance between the two arms centers 
curved_sections_length = sqrt((d*R)/2 - d^2/16);

%% Waveguide Refractive Indices specifications
core_height = 1 * 1e-6;
core_width = 1.1413 * 1e-6;

%% Core Indices for Power Calculation
center = P.Nx_main * P.padfactor / 2;
x_core_start = ceil(center - res*core_width/2);
x_core_end   = ceil(center + res*core_width/2);
y_core_start = ceil(center - res*core_height/2);
y_core_end   = ceil(center + res*core_height/2);

bare_e_array = readPermittivityTable('material_data/Bare.txt', wl);
ablated_e_array = readPermittivityTable('material_data/Ablated.txt', wl);
fprintf("Ablated e1 = %.6f , e2 = %.6f\n", ablated_e_array(1), ablated_e_array(2));
fprintf("Bare e1 = %.6f , e2 = %.6f\n", bare_e_array(1), bare_e_array(2));

ablated_n_tilde = permittivityToRefrIdx(ablated_e_array(1), ablated_e_array(2), 0);
bare_n_tilde = permittivityToRefrIdx(bare_e_array(1), bare_e_array(2), 1);

P.n_0 = real(ablated_n_tilde); % [] reference refractive index

sio2_n = sio2RefrIdx(P.lambda);

%% Segment 1
P.Lz = 3e-6; % [m] z propagation distances for this segment
P = initializeRIfromFunction(P,@calcRIsegs1and7,{core_height, core_width, ...
                             ablated_n_tilde, bare_n_tilde, sio2_n});
% Finding modes
P = findModes(P, 1, 'plotModes', false);

% Initial electric field is the fundamental mode
P.E = P.modes(1);

% E-field power decay analysis - INITIAL STATE
field_intensity = abs(P.E.field).^2;
in_total_power_sum = sum(field_intensity(:)); 

% Note: on E.field, X and Y indices are flipped! 
core_field_intensity = abs(P.E.field(x_core_start:x_core_end, y_core_start:y_core_end)).^2;
in_core_power_sum = sum(core_field_intensity(:));

P.figTitle = 'Segment 1';

P = FD_BPM(P);

%% Segment 2 
P.figTitle = 'Segment 2';
P.Lz = curved_sections_length;
Nz_n = ceil(P.Lz / P.dz_target);

P = initializeRIfromFunction(P,@calcRIseg2,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n},...
                             Nz_n);
P = FD_BPM(P);

%% Segment 3
P.figTitle = 'Segment 3';
P.Lz = curved_sections_length;

P = initializeRIfromFunction(P,@calcRIseg3,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n, d}, ...
                             Nz_n);
P = FD_BPM(P);


%% Segment 4
P.figTitle = 'Segment 4';
seg_4_len = 16e-6 - 4*curved_sections_length - 6e-6;
P.Lz = seg_4_len;
fprintf("Segment 4 length = %.3f um\n", seg_4_len*1e6);

% ablated_n_offset = 0.8384;
ablated_n_offset = 0;

P = initializeRIfromFunction(P, @calcRIseg4, {core_height, core_width,...
                                 ablated_n_tilde, bare_n_tilde, sio2_n, d, ablated_n_offset});
P = FD_BPM(P);

%% Segment 5
P.figTitle = 'Segment 5';
P.Lz = curved_sections_length;

P = initializeRIfromFunction(P,@calcRIseg5,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n, d}, ...
                             Nz_n);
P = FD_BPM(P);

%% Segment 6
P.figTitle = 'Segment 6';
P.Lz = curved_sections_length;

P = initializeRIfromFunction(P,@calcRIseg6,{core_height, core_width,...
                             ablated_n_tilde, bare_n_tilde, R, sio2_n, d}, ...
                             Nz_n);
P = FD_BPM(P);

elecFieldCrossSec_cpx = P.E.field(265:335,300);
electricFieldCrossSec = real(elecFieldCrossSec_cpx);

%% Segment 7
P.figTitle = 'Segment 7';

P.Lz = 3e-6; % [m] z propagation distances for this segment
P = initializeRIfromFunction(P,@calcRIsegs1and7,{core_height, core_width, ...
                             ablated_n_tilde, bare_n_tilde, sio2_n});
P = FD_BPM(P);
% 
% % finalizeVideo(P);
% 
% field_intensity = abs(P.E.field).^2;
% out_total_power_sum = sum(field_intensity(:));
% 
% core_field_intensity = abs(P.E.field(x_core_start:x_core_end, y_core_start:y_core_end)).^2;
% out_core_power_sum = sum(core_field_intensity(:));
% 
fprintf('--- E-field Power Numerical Analysis ---\n');
fprintf('Remaining total power (solver): %.8f\n', P.powers(end));
% 
% guided_in_total_in = in_core_power_sum/in_total_power_sum;
% fprintf('Guided_in / Total_in (must be close to 1): %.8f (%.2f%%)\n', guided_in_total_in, guided_in_total_in*100);
% 
% guided_initial = out_core_power_sum/in_total_power_sum;
% fprintf('Guided_out / Total_in: %.8f (%.2f%%)\n', guided_initial, guided_initial*100);
% 
% guided_final = out_core_power_sum/out_total_power_sum;
% fprintf('Guided_out / Total_out: %.8f (%.2f%%)\n', guided_final, guided_final*100);
% 
% guided_final_2 = out_core_power_sum/in_core_power_sum;
% fprintf('Guided_out / Guided_in: %.8f (%.2f%%)\n\n', guided_final_2, guided_final_2*100);

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
%% RI function for segments 1 and 7
function n = calcRIsegs1and7(X,Y,n_background,nParameters)
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

function complex_n = permittivityToRefrIdx(e1,e2,absorption)
e_mag = sqrt(e1^2 + e2^2);
n = sqrt((e_mag + e1) / 2);
k = sqrt((e_mag - e1) / 2);
complex_n = n + absorption*1i*k;
end

function n = sio2RefrIdx(x) 
x = x/(1e-6);
n = sqrt(1+0.665721./(1-(0.060./x).^2)+0.503511./(1-(0.106./x).^2)+0.214792./(1-(0.119./x).^2)+0.539173./(1-(8.792./x).^2)+1.807613./(1-(19.70./x).^2));
end

