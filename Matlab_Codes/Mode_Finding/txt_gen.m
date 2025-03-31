function txt_gen(data, name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   txt_gen creates a txt file with k values on the left column and
%   corresponding frequency values on the right column
%   
%   data -> Output data of the form output by load_harminv_dat
%   name -> Name of file in the form 'name.txt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
used_data = [data(:,1) data(:,2)];
writematrix(used_data, name, 'Delimiter', ' ')
end