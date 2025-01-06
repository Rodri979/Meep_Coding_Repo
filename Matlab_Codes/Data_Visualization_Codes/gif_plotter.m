function gif_plotter(directory, eps_file, slice, slice_num, resolution, title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   gif_plotter plots all modes in a directory. IT MUST ONLY CONTATION THE 
%   H5 denergy FILES OUTPUT FROM A MEEP SIMULATION. It is meant to plot all
%   files which are automatically saved (in the same directory as the h5 
%   files) as pngs and named (thus if a mistake was made and it needs to
%   run again delete all pngs in the directory and start over). The named
%   png files can then be converted into a gif using the convert *.png
%   name.gif function in bash.
%
%   THIS IS PART OF A PACKAGE OF CODE. IT NEEDS pull_max_min.m AND
%   Plot_3d_w_black_lines_gif.m to run 
%
%   **** NOTE
%   Plot_3d_w_black_lines_gif.m IS NOT the same as Plot_3d_w_black_lines.m
%
%   *** YOU MUST EDIT THE PATH IN LINES 56, 59, and 69 which read 
%   "/scratch/bell/rodri979/meep_files/matlab/" TO MATCH THE PATH 
%   DIRECTORY WHERE THE FOLDER CONTIANING ONLY MODE h5 FILES IS 
%   STORED ON YOUR DEVICE
%   
%   directory -> input of form "directory" (MUST BE IN DOUBLE QUOTES). This
%   is the name of the folder which has in it ONLY the h5 denergy files 
%   output from a meep simulation
%
%   slice -> 'x', 'y', or 'z'. Determines the plane that the slice will be a
%   part of (if 'x' is chosen then the yz plane will be used, 'y' 
%   then the xz plane). Slice to be input into Plot_3d_w_black_lines.m
% 
%   slice_num -> The number that the chosen slice dimension will be set to,
%   for example if slice_num = 75 and slice = x then the plot will be output
%   at x = 75. Slice_num to be put into Plot_3d_w_black_lines.m
%
%   resolution -> The resolution that the simulation was run in meep
%
%   title -> The title on all plots. It will be the title seen on the gif.
%
%   eps_file -> The epsilon structure h5 file with which the field will be
%   overlayed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_struct = dir(directory);
pic_names = [" " " " " "];
pic_index = 1;
%max_field = -1;
%min_field = 100000000000000;

for dir_num = linspace(3,length(dir_struct),(length(dir_struct)-2))
    k_name = dir_struct(dir_num).name;
    k_dir = dir(directory + "/" + k_name);
    for index = 1:1:length(k_dir)
        pic_names(pic_index) = k_dir(index).name;
        pic_index = pic_index + 1;
    end
end

%disp(pic_names)

[max_field, min_field] = pull_max_min("/scratch/bell/rodri979/meep_files/matlab/" + directory + "/" + pic_names(1), "/denergy", slice, slice_num);
   
for graph_index = 2:1:length(pic_names)
    [temp_max, temp_min] = pull_max_min("/scratch/bell/rodri979/meep_files/matlab/" + directory + "/" + pic_names(graph_index), "/denergy", slice, slice_num);
    if temp_max > max_field
        max_field = temp_max;
    end
    if temp_min < min_field
        min_field = temp_min;
    end
end

for graph_index = 1:1:length(pic_names)
    Plot_3d_w_black_lines_gif(eps_file, "/eps", "/scratch/bell/rodri979/meep_files/matlab/" + directory + "/" + pic_names(graph_index), "/denergy", slice, slice_num, max_field, min_field, resolution, title)
    fprintf(num2str(graph_index) + ": " + pic_names(graph_index))
    fprintf("\n")
    png_name = directory + "/" + slice + "_" + num2str(slice_num) + "_" + strrep(pic_names(graph_index),"h5","png");
    saveas(figure(graph_index),png_name, 'png')
end

end

