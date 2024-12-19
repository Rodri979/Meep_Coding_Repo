# Meep_Codeing_Repo
A repository for the meep files run for photonics research.

This repository is split into two sections. One section contains scheme code and code used specifically in Purdue's Bell Cluster to run a meep scheme control file. Examples of code in Meep_Codes contain code used for finding resosant modes (SiO2_boundary.ctl) and bash scripts designed for use on Purdue's bell cluster which can submit an arbitrary amount of jobs at once (pic run). The scheme codes themselves are skeleton programs which can be modified to fit any particular geometry, boundary layers, source etc...

Matlab_Codes contain matlab programs that are used in conjunction with the scheme meep programs. Some of these include a function which sifts through meep "harminv" resonant mode data and find modes which fit a particular catagory  (i.e. Q > 300, amplitude > 0.01 etc...) [mode_finder.m], functions which ouput a txt file of all the particular interesting modes to be used by pic_run [txt_gen.m]. This also contatins mode visualization functions that can plot 2d slices of 3d geometry (eps_plotter.m) and mode profiles (Plot_3d_w_black_lines.m).
