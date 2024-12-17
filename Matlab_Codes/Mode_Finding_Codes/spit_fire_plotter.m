function spit_fire_plotter(directory, mode_data, eps_file, slice, slice_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   spit_fire_plotter plots all modes in a directory
%
%   *** YOU MUST EDIT THE PATH IN LINE 33 which reads "/scratch/bell/rodri979/
%   meep_files/matlab/" TO MATCH THE DIRECTORY WHERE h5 DATA GENERATED FROM 
%   pic_run.sh IS ON YOUR SPECIFIC COMPUTER***
%   
%   directory -> input of form "directory" (MUST BE IN DOUBLE QUOTES) which
%   has more directories inside named kz_# with # given by the k values in 
%   mode_data
%   
%   mode_data -> Data of the format output by mode_finder (this is the same
%   data used to run the mode images)
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
    Plot_3d_w_black_lines(eps_file, "/eps", "/scratch/bell/rodri979/meep_files/matlab/" + directory + "/kz_" + pic_names(graph_index, 1) + "/" + pic_names(graph_index, 3), "/denergy", slice, slice_num, "|E| at z=25 for k=" + pic_names(graph_index, 1) + " and w=" + pic_names(graph_index, 2))
end

end
