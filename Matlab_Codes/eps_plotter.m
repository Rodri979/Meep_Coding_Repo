function eps_plotter(eps_data, slice, slice_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% eps_data -> eps h5 datafile. Input the name of your eps h5 file in quotes
% ex: "straight_waveguide_recreate-eps-000000.00.h5"
%
% slice -> 'x', 'y', or 'z' slice of data you want to take. This combined with
% slice num is like chosing a plane that goes through a specified point.
% For example if slice = 'x' and slice_num = 37 then you are selecting the
% yz plane which cuts through x = 37. FYI the x,y,z values here do not
% exactly correlate with xyz values in meep so play around with slice until
% you get the visualization you are looking for.
%
% slice_num -> The value to which the "slice" variable will be assigned.
% e.g. if slice = 'y' and slice_num = 53 then the plane that satisfies 
% y = 53 will be shown
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_data = h5read(eps_data, "/eps");

switch slice
    case 'z'
        eps_slice = squeeze(eps_data(slice_num, :, :));
    case 'y'
        eps_slice = squeeze(eps_data(:, slice_num, :));
    case 'x'
        eps_slice = squeeze(eps_data(:, :, slice_num));
    otherwise
        warning('%s is an unexpected slice dimension please choose x, y, or z', slice)
        return
end

% Visualize the 2D slice
figure;
imagesc(flipud(eps_slice));
colorbar; 
colormap(flipud(gray));
title("Eps 2d Slice at " + slice + " = " + string(slice_num));
axis square;
switch slice
    case 'z'
        xlabel('x');
        ylabel('y');
    case 'y'
        xlabel('x');
        ylabel('z');
    case 'x'
        xlabel('y');
        ylabel('z');
end
end