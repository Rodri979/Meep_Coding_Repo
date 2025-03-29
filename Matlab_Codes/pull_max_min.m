function [maximum, minimum] = pull_max_min(ez_h5, ez_code, slice, slice_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function pulls the max field profile value and min field profile
%   value in an h5 profile and output them. This function is used in
%   gif_plotter.m to find the overall max and min field profile values in 
%   a group of h5 files (through a loop).
%
%   ez_h5 -> 3d h5 datafile for efield
%
%   ez_code -> either "/denergy" for a dpwr h5, "/ez" for and ez file,
%   "/ex" for an ex file. h5read needs this code to read in the data
%
%   slice -> 'x', 'y', or 'z'. Determines the plane that the slice will be a
%   part of (if 'x' is chosen then the yz plane will be used, 'y' 
%   then the xz plane)
% 
%   slice_num -> The number that the chosen slice dimension will be set to,
%   for example if slice_num = 75 and slice = x then the plot will be output
%   at x = 75
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ez_data = h5read(ez_h5, ez_code);

switch slice
    case 'z'
        ez_slice = squeeze(ez_data(slice_num, :, :));
    case 'y'
        ez_slice = squeeze(ez_data(:, slice_num, :));
    case 'x'
        ez_slice = squeeze(ez_data(:, :, slice_num));
    otherwise
        warning('%s is an unexpected slice dimension please choose x, y, or z', slice)
        return
end

maximum = max(ez_slice, [], "all");
minimum = min(ez_slice, [], "all");

end
