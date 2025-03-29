# Meep_Codeing_Repo
A repository for the meep files run for photonics research.

This repository is split into two sections. One section contains scheme code and code used specifically in Purdue's Bell Cluster to run a meep scheme control file. Examples of code in Meep_Codes contain code used for finding resosant modes (SiO2_boundary.ctl) and bash scripts designed for use on Purdue's bell cluster which can submit an arbitrary amount of jobs at once (pic run). The scheme codes themselves are skeleton programs which can be modified to fit any particular geometry, boundary layers, source etc...

Matlab_Codes contain matlab programs that are used in conjunction with the scheme meep programs. Some of these include a function which sifts through meep "harminv" resonant mode data and find modes which fit a particular catagory  (i.e. Q > 300, amplitude > 0.01 etc...) [mode_finder.m], functions which ouput a txt file of all the particular interesting modes to be used by pic_run [txt_gen.m]. This also contatins mode visualization functions that can plot 2d slices of 3d geometry (eps_plotter.m) and mode profiles (Plot_3d_w_black_lines.m).

The new structure of this repo makes it harder to determine the use of all functions but easier to use the repo in practice below is a structure that lists the functions of certain codes:

# Meep_Codes
## Mode Finding Codes
- pic_run.sh -> Runs many mode outputing programs (in either scheme or python which output the energy desity of modes) in parallel and outputs corresponding h5 files at a desired directrory (typically inside the matlab directory) in structured folders named kz_#. The modes are previously found by Harminv and list output by txt_gen.m (params.txt) is input to pic_run.sh to tell it what modes to run. The matlab program spit_fire_plotter.m plots a common 2D slice of all the modes in these folders and they can be sifted manually to find modes of interest.
- submit_job.sh -> Run this to run pic_run.sh, outputs the length of params.txt if it is needed.
## Gif Plotting Codes
- straight_wagveguide_recreate_rough -> Can create a gif starting when the source is turned off and continuming for num_periods with num_pictures output per num_periods
## Flux Codes
- straight_waveguide_recreate_rough -> Can output the flux leaving the index waveguide until the field inside decays by a certain amount.

# Matlab Codes
