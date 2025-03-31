function spit_fire_plotter(directory, mode_data, eps_file, slice, resolution, slice_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   spit_fire_plotter plots all modes in a directory
%   
%   directory -> input of form "directory" (MUST BE IN DOUBLE QUOTES) which
%   has more directories inside named kz_# with # given by the k values in 
%   mode_data
%   
%   mode_data -> Data of the format output by mode_finder (this is the same
%   data used to run the mode images)
%
%   slice -> 'x', 'y', or 'z'. Determines the plane that the slice will be a
%   part of (if 'x' is chosen then the yz plane will be used, 'y' 
%   then the xz plane). Slice to be input into Plot_3d_w_black_lines.m
% 
%   slice_num -> The number that the chosen slice dimension will be set to,
%   for example if slice_num = 75 and slice = x then the plot will be output
%   at x = 75. Slice_num to be put into Plot_3d_w_black_lines.m
%
%   resolution -> The value of the resolution which the simulation was run in
%   meep. Used to calculate the size (in um) of the dimensions. If the refrence
%   dimension is not 1 um, you can change the hardcoded um dimension in the 
%   bottom switch case in Plot_3d_w_black_lines for now.
%
%   eps_file -> The epsilon structure h5 file with which the field will be
%   overlayed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_struct = dir(directory);
pic_names = [" " " " " "];
pic_index = 1;

for dir_num = linspace(3,length(dir_struct),(length(dir_struct)-2))
    k_name = dir_struct(dir_num).name;
    k_dir = dir(directory + "/" + k_name);
    for index = 5:10:length(k_dir)
        pic_names(pic_index, 3) = k_dir(index).name;
        pic_names(pic_index, 2) = round(mode_data(pic_index, 2), 3);
        pic_names(pic_index, 1) = mode_data(pic_index, 1);
        pic_index = pic_index + 1;
    end
end
   disp(pic_names)
for graph_index = linspace(1, length(pic_names(:,1)), length(pic_names(:,1)))
    Plot_3d_w_black_lines(eps_file, "/eps", "./" + directory + "/kz_" + pic_names(graph_index, 1) + "/" + pic_names(graph_index, 3), "/denergy", slice, slice_num, resolution, "|E| at z=25 for k=" + pic_names(graph_index, 1) + " and w=" + pic_names(graph_index, 2))
end

end
