function Plot_3d_w_black_lines_gif(eps_h5, eps_code, ez_h5, ez_code, slice, slice_num, max_field, min_field, resolution, title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function, much like Plot_3d_w_black_lines.m, overlays a plot of a
%   field profile over a plot of geometry. This function, however, is only
%   used in the gif_plotter function and is not meant to be used on its own.
%   The difference is that this plot sets a preset scale for the field colors
%   so that a gif mantains a constant field scaling throughout and does not
%   change between pngs. The max amplitude for a field in the scaling is 
%   input in max_field and the min_field for a scaling (could just be set to 0
%   for a energy density graph) is the minimum field value for scaling. In 
%   gif_plotter.m it cycles through each field profile, finds the absolute min
%   and max field values and then these are input into this function to set a
%   common scaling.
%
%   eps_h5 -> 3d h5 datafile for the dielectric strucutre
%
%   ez_h5 -> 3d h5 datafile for efield 
%
%   slice -> 'x', 'y', or 'z'. Determines the plane that the slice will be a
%   part of (if 'x' is chosen then the yz plane will be used, 'y' 
%   then the xz plane)
% 
%   slice_num -> The number that the chosen slice dimension will be set to,
%   for example if slice_num = 75 and slice = x then the plot will be output
%   at x = 75
%
%   max_field -> Maximum field amplitude for scaling purposes (in gif_plotter.m
%   this is found through a loop and input to this function)
%
%   min_field -> Minimum field amplutude (0 for energy density plots). In 
%   gif_plotter.m this is found through a loop and input to this function.
%
%   resolution -> The resolution that the simulation was run in meep
%
%   title -> Title of plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_data = h5read(eps_h5, eps_code);
ez_data = h5read(ez_h5, ez_code);
w_b_color = [1 1 1; 0 0 0]; % Defining a colormap with black for high values
b_w_color = [0 0 0; 1 1 1]; % Defining a colormap with white for high values
data_size = size(eps_data);

if(length(data_size) == 3)
    switch slice
        case 'z'
            eps_slice = squeeze(eps_data(slice_num, :, :));
            ez_slice = squeeze(ez_data(slice_num, :, :));
        case 'y'
            eps_slice = squeeze(eps_data(:, slice_num, :));
            ez_slice = squeeze(ez_data(:, slice_num, :));
        case 'x'
            eps_slice = squeeze(eps_data(:, :, slice_num));
            ez_slice = squeeze(ez_data(:, :, slice_num));
        otherwise
            warning('%s is an unexpected slice dimension please choose x, y, or z', slice)
            return
    end
elseif((length(data_size) == 2) && (data_size(1) ~= 1))
    eps_slice = eps_data;
    ez_slice = ez_data;
end

figure;

%plot first data 
ax1 = axes; 
%im = imagesc(ax1,flipud(floor(eps_slice/max(eps_slice(:)))));
im = imagesc(ax1,flipud(eps_slice/max(eps_slice(:))));
im.AlphaData = 0.6; % change this value to change the background image transparency 
axis square; 
hold all; 
   
%plot second data 
ax2 = axes; 
im1 = imagesc(ax2,flipud(ez_slice)); 
im1.AlphaData = 0.75; % change this value to change the foreground image transparency 
axis square; 
    
%link axes 
linkaxes([ax1,ax2]) 
    
%%Hide the top axes 
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = []; 
    
colormap(ax1,flipud(gray)) 
colormap(ax2,'jet')
caxis([min_field, max_field])
set([ax1,ax2],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]); 
    
sgtitle(title);
    
    
ylabel(cb2, '|E|', 'Position', [.5 .0], 'Rotation', 0, 'FontSize', 16);

min_2 = -1 * round((size(eps_slice,2) - 1)/(2 * resolution),1);
max_2 = -1 * min_2;

min_1 = -1 * round((size(eps_slice,1) - 1)/(2 * resolution), 1);
max_1 = -1 * min_1;

set(ax1, 'XTick', linspace(1,size(eps_slice,2),5), 'XTickLabel', linspace(min_2,max_2,5), 'FontSize', 12);
set(ax1, 'YTick', linspace(1,size(eps_slice,1),5), 'YTickLabel', linspace(max_1, min_1, 5), 'FontSize', 12);

switch slice
    case 'z'
        xlabel(ax1, 'X (\mum)', 'FontSize', 14);
        ylabel(ax1, 'Y (\mum)', 'FontSize', 14);
    case 'y'
        xlabel(ax1, 'X (\mum)', 'FontSize', 14);
        ylabel(ax1, 'Z (\mum)', 'FontSize', 14);
    case 'x'
        xlabel(ax1, 'Y (\mum)', 'FontSize', 14);
        ylabel(ax1, 'Z (\mum)', 'FontSize', 14);
end
    
drawnow
hold off;

end



