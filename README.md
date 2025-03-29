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
## Data Visualization Codes
- Plot_3d_w_black_lines.m -> Outputs a heat map of a 2D slice from a 3D h5 file. Overlays field data on top of eps (structure) data. Must have both a field profile file (or an energy density file) and an eps file to run. The eps_code will always be "/eps", and the ez_code for an energy density file is "/denergy".
- Plot_3d_w_blak_lines_gif.m -> Again outputs a heat map of a 2D slice from a 3D h5 file. The only difference from Plot_3d_w_black_lines.m is that this is used for plotting gifs (used inside of gif_plotter). It takes a minimum field value and maximum field value which manually sets the scale for the heat map, used to keep a common scale in gifs.
- eps_plotter.m -> Plots a 2D slice of a 3D h5 eps file to visualize structures. The higher the eps the darker it will appear in the plot.
- gif_plotter -> Plots a gif of a 2D slice in a set of 3D h5 files. The field or energy density files must all be inside of a folder within the matlab directory, the eps file MUST NOT be in this folder but inside the matlab directory itself. One inputs a directory and it plots and saves the 2D heat map as a png for all h5 files in the directory. Once this has been done one can use imagemagick to convert the outputted pngs into a gif.
- pull_max_min.m -> This file is used in gif_plotter.m to take the maximum and minimum field values of a slice. Then gif_plotter.m loops through all slices in a directory and finds the absolute maximum and minimum values. It uses these to set a common scale for all outputted pngs to make a gif.
## Dispersion Codes
- Dispersion_Models.m -> Contains the parameters for several example dispersion models that can be inplemented in meep, and displays their dispersion plots to be compared to experiment. This can be added to as more dispersion models are obtained.
- View_Dispersion.m -> Takes the values that meep will use to define a lorentzian dispersion model and will output a plot of a given range that can be used to compare to experiment. If this plot matches experiment then these values are good to use for a meep lorentzian dispersion model.
- fitNlorentzian -> This uses stocasitc gradient descent to approximate any given expiremental data to a linear combination of lorentzians. It inputs these values to View_Dispersion.m for comparison and then outputs the parameter values to be used in meep simulation.
