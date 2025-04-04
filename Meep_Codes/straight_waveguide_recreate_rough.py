import meep as mp
import argparse
import numpy as np

def main(args):
    resolution = 50   # pixels/um
    n_bto = 2.278 # BTO index
    n_sio2 = 1.44 # SiO2 index
    lip_width = args.lip_width_um   # Setting lip_width (um)
    lip_height = args.lip_height_um # Setting lip height (um)
    bto_height = args.bto_height_um # Setting bto height (um)
    sio2_height = args.sio2_height_um # Setting SiO2 height (um)
    output_gif = args.output_gif # True or false value if a gif is to be output
    sio2_offset = (-2-0.5)*0.5 # Value to bring sio2 down to make the structure in the center of simulation zone
    freq = args.freq # Setting center frequency
    nfreq = args.nfreq # Number of frequencies to calcualte flux for
    df = args.df # Frequnecy spread of source
    num_pics = args.num_pics # Number of h5 files output per period
    num_periods = args.num_periods # Number of periods of h5 files outputed
    kz = args.kz # kz value for boundary conditions
    avg_num_bumps = args.avg_num_bumps # mean of normal distribution of bumps, each face will have random number
    std_num_bumps = args.std_num_bumps # standard deviatation
    avg_r_bto_air = args.avg_r_bto_air # Mean radius of bump on a bto to air interface
    std_r_bto_air = args.std_r_bto_air # Std of bump on a bto-air interface
    avg_r_bto_sio2 = args.avg_r_bto_sio2 # Mean radius of bump on a bto-sio2 interface
    std_r_bto_sio2 = args.std_r_bto_sio2 # Std radius of bump on bto_sio2 interface
    sx = 10 # size of cell in x direction (perpendicular to wvg.)
    sy = 4 # size of cell in y direction (perpendicular to wvg.)
    sz = 1 # size of cell in z direction (parallel to wvg.)
    pad = 0.1234 # padding between last hole and PML edge
    dpml = 1 # PML thickness
    epsilon = args.epsilon # Defining epsilon value away from faces for flux planes (+ is bigger then waveguide - is smaller)

    cell = mp.Vector3(sx,sy,sz)

    sio2_block = mp.Block(size=mp.Vector3(mp.inf,sio2_height,mp.inf),center=mp.Vector3(0,sio2_offset,0),material=mp.Medium(index=n_sio2))
    bto_block = mp.Block(size=mp.Vector3(mp.inf,bto_height,mp.inf),center=mp.Vector3(0,sio2_offset+0.5*sio2_height+0.5*bto_height,0),material=mp.Medium(index=n_bto))
    lip_block = mp.Block(size=mp.Vector3(lip_width,lip_height,mp.inf),center=mp.Vector3(0,sio2_offset+0.5*sio2_height+bto_height+0.5*lip_height,0),material=mp.Medium(index=n_bto))
    geometry = [sio2_block,bto_block,lip_block] 

    top_area = lip_width*sz # Getting area of top face to normalize bump density of each face
    lip_side_area = lip_height*sz # Getting area of side of lip to normailize bump density
    area_ratio = lip_side_area/top_area

# BUMPS ON TOP
    top_num_bumps = int(round(np.random.normal(avg_num_bumps,std_num_bumps))) # defining number of bumps on top of structure
    top_x_positions = np.random.uniform(-lip_width/2,lip_width/2,top_num_bumps) # i random distribution of x-coordinates for the bumps in the lip_height and bto_height
    top_z_positions = np.random.uniform(-sz*0.5,sz*0.5,top_num_bumps) # random distribution of z-coordinates for on the top of the waveguide
    top_y_pos = sio2_offset + sio2_height*0.5 + bto_height + lip_height # Placing top bumps on top plane
    top_r = np.random.normal(avg_r_bto_air,std_r_bto_air,top_num_bumps) # defining the radius of the bumps as a normal random variable

# BUMPS ON LEFT SIDE
    left_num_bumps= int(round(np.random.normal(avg_num_bumps*area_ratio,std_num_bumps*area_ratio))) # defining number of bumps on left of structure
    left_y_positions = np.random.uniform(sio2_offset + sio2_height*0.5 + bto_height,sio2_offset + sio2_height*0.5 + bto_height+ lip_height,left_num_bumps) # random distribution of x-coordinates for the bumps in the lip_height and bto_height
    left_z_positions = np.random.uniform(-sz*0.5,sz*0.5,left_num_bumps) # random distribution of z-coordinates for on the top of the waveguide
    left_x_pos = -lip_width*0.5 # Placing top bumps on top plane
    left_r = np.random.normal(avg_r_bto_air,std_r_bto_air,left_num_bumps) # defining the radius of the bumps as a normal random variable

# BUMPS ON RIGHT SIDE
    right_num_bumps = int(round(np.random.normal(avg_num_bumps*area_ratio/top_,std_num_bumps))) # defining number of bumps on right of structure
    right_y_positions = np.random.uniform(sio2_offset + sio2_height*0.5 + bto_height,sio2_offset + sio2_height*0.5 + bto_height + lip_height,right_num_bumps) # i random distribution of x-coordinates for the bumps in the lip_height and bto_height
    right_z_positions = np.random.uniform(-sz*0.5,sz*0.5,right_num_bumps) # random distribution of z-coordinates for on the top of the waveguide
    right_x_pos = lip_width*0.5 # Placing top bumps on top plane
    right_r = np.random.normal(avg_r_bto_air,std_r_bto_air,right_num_bumps) # defining the radius of the bumps as a normal random variable

# BUMPS ON BOTTOM
    bot_num_bumps = int(round(np.random.normal(avg_num_bumps,std_num_bumps))) # defining number of bumps on bottom of structure
    bot_x_positions = np.random.uniform(-lip_width/2,lip_width/2,bot_num_bumps) # i random distribution of x-coordinates for the bumps in the lip_height and bto_height
    bot_z_positions = np.random.uniform(-sz*0.5,sz*0.5,bot_num_bumps) # random distribution of z-coordinates for on the top of the waveguide
    bot_y_pos = sio2_offset + sio2_height*0.5 # Placing top bumps on top plane
    bot_r = np.random.normal(avg_r_bto_sio2,std_r_bto_sio2,bot_num_bumps) # defining the radius of the bumps as a normal random variable

    # SYNTAX FOR APPENDING STRUCTURES
    # ADDING TOP BUMPS
    for i in range(top_num_bumps):
            coin_flip = np.random.randint(2) # Flipping a coing to see if bump cuts into or grows out of structure
            if coin_flip == 0:
                geometry.append(mp.Sphere(radius=top_r[i], center=mp.Vector3(top_x_positions[i],top_y_pos,top_z_positions[i]),material=mp.Medium(index=n_bto))) # Sphere of BTO
            else:
                geometry.append(mp.Sphere(radius=top_r[i], center=mp.Vector3(top_x_positions[i],top_y_pos,top_z_positions[i]),material=mp.Medium(index=1))) # Sphere of air

    # ADDING LEFT BUMPS
    for i in range(left_num_bumps):
            coin_flip = np.random.randint(2) # Flipping a coing to see if bump cuts into or grows out of structure
            if coin_flip == 0:
                geometry.append(mp.Sphere(radius=left_r[i], center=mp.Vector3(left_x_pos,left_y_positions[i],left_z_positions[i]),material=mp.Medium(index=n_bto))) # Sphere of BTO
            else:
                geometry.append(mp.Sphere(radius=left_r[i], center=mp.Vector3(left_x_pos,left_y_positions[i],left_z_positions[i]),material=mp.Medium(index=1))) # Sphere of air

    # ADDING RIGHT BUMPS
    for i in range(right_num_bumps):
            coin_flip = np.random.randint(2) # Flipping a coing to see if bump cuts into or grows out of structure
            if coin_flip == 0:
                geometry.append(mp.Sphere(radius=right_r[i], center=mp.Vector3(right_x_pos,right_y_positions[i],right_z_positions[i]),material=mp.Medium(index=n_bto))) # Sphere of BTO
            else:
                geometry.append(mp.Sphere(radius=right_r[i], center=mp.Vector3(right_x_pos,right_y_positions[i],right_z_positions[i]),material=mp.Medium(index=1))) # Sphere of air

    # ADDING BOTTOM BUMPS
    for i in range(bot_num_bumps):
            coin_flip = np.random.randint(2) # Flipping a coing to see if bump cuts into or grows out of structure
            if coin_flip == 0:
                geometry.append(mp.Sphere(radius=bot_r[i], center=mp.Vector3(bot_x_positions[i],bot_y_pos,bot_z_positions[i]),material=mp.Medium(index=n_bto))) # Sphere of BTO
            else:
                geometry.append(mp.Sphere(radius=bot_r[i], center=mp.Vector3(bot_x_positions[i],bot_y_pos,bot_z_positions[i]),material=mp.Medium(index=n_sio2))) # Sphere of sio2


    pml_layers = [mp.PML(thickness=dpml, direction=mp.X),mp.PML(thickness=dpml, direction=mp.Y)] # Setting pml

    src = [mp.Source(mp.GaussianSource(freq, fwidth=df),
                     component=mp.Hz,
                     center=mp.Vector3(0,sio2_offset+0.5*sio2_height+0.5*(bto_height+lip_height),-0.5*sz+pad),
                     size=mp.Vector3(lip_width*0.7,0.7*(lip_height+bto_height),0))] 

    sim = mp.Simulation(cell_size=cell,
                    geometry=geometry,
                    boundary_layers=pml_layers,
                    sources=src,
                    resolution=resolution)

    sim.k_point = mp.Vector3(0,0,kz)

    left_flux = mp.FluxRegion(center=mp.Vector3(-0.5*lip_width-epsilon,sio2_offset+0.5*(sio2_height+bto_height+lip_height),0),
                         size=mp.Vector3(0,lip_height + bto_height + 2*epsilon, sz), weight=-1)
    right_flux = mp.FluxRegion(center=mp.Vector3(0.5*lip_width+epsilon,sio2_offset+0.5*(sio2_height+bto_height+lip_height),0)                         ,size=mp.Vector3(0,lip_height + bto_height + 2*epsilon, sz), weight=1)
    top_flux = mp.FluxRegion(center=mp.Vector3(0,sio2_offset+0.5*sio2_height+bto_height+lip_height+epsilon,0)                                         ,size=mp.Vector3(lip_width+2*epsilon, 0, sz), weight=1)
    bottom_flux = mp.FluxRegion(center=mp.Vector3(0,sio2_offset+0.5*sio2_height-epsilon,0)                                                            ,size=mp.Vector3(lip_width+2*epsilon, 0, sz), weight=-1)

    # Adding flux
    flux = sim.add_flux(freq, df, nfreq, left_flux, right_flux, top_flux, bottom_flux)
    
    if(output_gif):
        sim.run(mp.at_beginning(mp.output_epsilon),
                mp.after_sources(mp.at_every(1/(freq*num_pics), mp.output_dpwr)),until_after_sources=num_periods/freq)
    else:
        sim.run(mp.at_beginning(mp.output_epsilon),
                until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(0,sio2_offset+0.5*(sio2_height+bto_height+lip_height), -0.5*sz+1.5*pad), 1e-2))

        sim.display_fluxes(flux)  # print out the flux spectrum

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-freq', type=float, default=0.65, help='pulse center frequency')
    parser.add_argument('-df', type=float, default=0.5, help='pulse frequency width')
    parser.add_argument('-lip_width_um', type=float, default=1.5, help='with of top part of structure')
    parser.add_argument('-lip_height_um', type=float, default=0.175, help='height of top part of structure')
    parser.add_argument('-bto_height_um', type=float, default=0.325, help='intitial BTO layer thickness')
    parser.add_argument('-sio2_height_um', type=float, default=2, help='SiO2 layer thickness')
    parser.add_argument('-output_gif', type=bool, default=False, help='Boolean to output a field profile gif or not')
    parser.add_argument('-kz', type=float, default=1.25, help='k in z direction')
    parser.add_argument('-num_pics', type=float, default=10, help='Gif resolution - number of pics per period')
    parser.add_argument('-nfreq', type=int, default=300, help='Number of frequencies to output in flux calc')
    parser.add_argument('-num_periods', type=float, default=5, help='Number of periods to record for gif')
    parser.add_argument('-epsilon', type=float, default=0.1, help='Length of epsilon away from faces for flux planes')
    parser.add_argument('-avg_num_bumps', type=float, default=50, help='Average of normal dist for num bumps top/bot face')
    parser.add_argument('-std_num_bumps', type=float, default=5, help='Std for normal dist of num bumps on each face')
    parser.add_argument('-avg_r_bto_air', type=float, default=0.1, help='Mean radius of bump for BTO-air interface')
    parser.add_argument('-std_r_bto_air', type=float, default=0.02, help='Std for bumps at BTO-air interface')
    parser.add_argument('-avg_r_bto_sio2', type=float, default=0.1, help='Mean radius of bump for BTO-sio2 interface')
    parser.add_argument('-std_r_bto_sio2', type=float, default=0.02, help='Std for bumps at BTO-sio2 interface')
    args = parser.parse_args()
    main(args)

# (define-param BTO_eps_0 1)
# (define-param BTO-frq1 (/ 1 0.223))
# (define-param BTO-gam1 0)
# (define-param BTO-sig1 4.187)

