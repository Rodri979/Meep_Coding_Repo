function eps_plotter(eps_data, slice, slice_num, resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% eps_data -> eps h5 datafile. Input the name of your eps h5 file in quotes
% ex: "straight_waveguide_recreate-eps-000000.00.h5"
%
% slice -> 'x', 'y', or 'z' slice of data you want to take. This combined with
% slice num is like chosing a plane that goes through a specified point.
% For example if slice = 'x' and slice_num = 37 then you are selecting the
% yz plane which cuts through x = 37.
%
% slice_num -> The value to which the "slice" variable will be assigned.
% e.g. if slice = 'y' and slice_num = 53 then the plane that satisfies 
% y = 53 will be shown
%
% resolution -> The value of the resolution which the simulation was run in
% meep. Used to calculate the size (in um) of the dimensions. If the refrence
% dimension is not 1 um, you can change the hardcoded um dimension in the 
% bottom switch case for now.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_data = h5read(eps_data, "/eps");
data_size = size(eps_data)
if(length(data_size) == 3)
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
elseif((length(data_size) == 2) && (data_size(1) ~= 1))
    eps_slice = eps_data
end


% Visualize the 2D slice
figure;
imagesc(flipud(eps_slice));
colorbar; 
colormap(flipud(gray));
title("Eps 2d Slice at " + slice + " = " + string(slice_num));
axis square;
min_2 = -1 * round((size(eps_slice,2) - 1)/(2 * resolution),1);
max_2 = -1 * min_2;

min_1 = -1 * round((size(eps_slice,1) - 1)/(2 * resolution), 1);
max_1 = -1 * min_1;

xticks(linspace(1,size(eps_slice,2),5))
xticklabels(linspace(min_2,max_2,5))

yticks(linspace(1,size(eps_slice,1),5))
yticklabels(linspace(max_1, min_1, 5))

switch slice
    case 'z'
        xlabel('X (\mum)');
        ylabel('Y (\mum)');
    case 'y'
        xlabel('X (\mum)');
        ylabel('Z (\mum)');
    case 'x'
        xlabel('Y (\mum)');
        ylabel('Z (\mum)');
end

end
