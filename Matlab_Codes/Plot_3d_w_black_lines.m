function Plot_3d_w_black_lines(eps_h5, eps_code, ez_h5, ez_code, slice, slice_num, title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eps_h5 -> 3d h5 datafile for the dielectric strucutre
%
% ez_h5 -> 3d h5 datafile for efield 
%
% slice -> 'x', 'y', or 'z'. Determines the plane that the slice will be a
%   part of (if 'x' is chosen then the yz plane will be used, 'y' 
%   then the xz plane)
% 
% slice_num -> The number that the chosen slice dimension will be set to,
% for example if slice_num = 75 and slice = x then the plot will be output
% at x = 75
%
% title -> Title of plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_data = h5read(eps_h5, eps_code);
ez_data = h5read(ez_h5, ez_code);
w_b_color = [1 1 1; 0 0 0]; % Defining a colormap with black for high values
b_w_color = [0 0 0; 1 1 1]; % Defining a colormap with white for high values

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
im1 = imagesc(ax2,flipud(abs(ez_slice))); 
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
    
set([ax1,ax2],'Position',[.17 .11 .685 .815]); 
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]); 
    
sgtitle(title);
    
    
ylabel(cb2, '|E|', 'Position', [.5 .0], 'Rotation', 0, 'FontSize', 16);
    
set(ax1, 'XTick', [1 round(size(eps_slice,2)/2) size(eps_slice,2)], 'XTickLabel', {'-5' '0' '5'}, 'FontSize', 12);
set(ax1, 'YTick', 1:round(size(eps_slice,1)/4):size(eps_slice,2), 'YTickLabel', linspace(-2, 2, 5), 'FontSize', 12);

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



