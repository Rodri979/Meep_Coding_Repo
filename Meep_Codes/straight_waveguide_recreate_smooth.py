import meep as mp
import argparse

def main(args):
    resolution = 50   # pixels/um

    n_bto = 2.278 # BTO index
    n_sio2 = 1.44 # SiO2 index
    lip_width = args.lip_width_um   # Setting lip_width (um)
    lip_height = args.lip_height_um # Setting lip height (um)
    bto_height = args.bto_height_um # Setting bto height (um)
    sio2_height = args.sio2_height_um # Setting SiO2 height (um)
    output_gif = args.output_gif # True or false value if a gif is to be output
    freq = args.freq # Setting center frequency
    nfreq = args.nfreq # Number of frequencies to calcualte flux for
    df = args.df # Frequnecy spread of source
    num_pics = args.num_pics # Number of h5 files output per period
    num_periods = args.num_periods # Number of periods of h5 files outputed
    kz = args.kz # kz value for boundary conditions
    
    sx = 10 # size of cell in x direction (perpendicular to wvg.)
    sy = 4 # size of cell in y direction (perpendicular to wvg.)
    sz = 1 # size of cell in z direction (parallel to wvg.)
    pad = 0.1234 # padding between last hole and PML edge
    dpml = 1 # PML thickness
    sio2_offset = (sio2_height-sy)*0.5 # Value to bring sio2 down to make the structure in the center of simulation zone
    epsilon = args.epsilon # Defining epsilon value away from faces for flux planes (+ is bigger then waveguide - is smaller)

    cell = mp.Vector3(sx,sy,sz)

    sio2_block = mp.Block(size=mp.Vector3(mp.inf,sio2_height,mp.inf),center=mp.Vector3(0,sio2_offset,0),material=mp.Medium(index=n_sio2))
    bto_block = mp.Block(size=mp.Vector3(mp.inf,bto_height,mp.inf),center=mp.Vector3(0,sio2_offset+0.5*sio2_height+0.5*bto_height,0),material=mp.Medium(index=n_bto))
    lip_block = mp.Block(size=mp.Vector3(lip_width,lip_height,mp.inf),center=mp.Vector3(0,sio2_offset+0.5*sio2_height+bto_height+0.5*lip_height,0),material=mp.Medium(index=n_bto))
    geometry = [sio2_block,bto_block,lip_block]

    # SYNTAX FOR APPENDING STRUCTURES
    # for i in range(N):
    #        geometry.append(mp.Cylinder(r, center=mp.Vector3(d/2+i)))
    #        geometry.append(mp.Cylinder(r, center=mp.Vector3(-(d/2+i))))

    pml_layers = [mp.PML(thickness=dpml, direction=mp.X),mp.PML(thickness=dpml, direction=mp.Y)] # Setting pml

    src = [mp.Source(mp.GaussianSource(freq, fwidth=df),
                     component=mp.Hz,
                     center=mp.Vector3(0,sio2_offset+0.5*sio2_height+0.5*(bto_height+lip_height),-0.5*sz+pad),
                     size=mp.Vector3(lip_width*0.8,0.8*(lip_height+bto_height),0))] 

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
    parser.add_argument('-sio2_height_um', type=float, default=1.5, help='SiO2 layer thickness')
    parser.add_argument('-output_gif', type=bool, default=False, help='Boolean to output a field profile gif or not')
    parser.add_argument('-kz', type=float, default=1.25, help='k in z direction')
    parser.add_argument('-num_pics', type=float, default=10, help='Gif resolution - number of pics per period')
    parser.add_argument('-nfreq', type=int, default=300, help='Number of frequencies to output in flux calc')
    parser.add_argument('-num_periods', type=float, default=5, help='Number of periods to record for gif')
    parser.add_argument('-epsilon', type=float, default=0.1, help='Length of epsilon away from faces for flux planes')
    args = parser.parse_args()
    main(args)

# (define-param BTO_eps_0 1)
# (define-param BTO-frq1 (/ 1 0.223))
# (define-param BTO-gam1 0)
# (define-param BTO-sig1 4.187)

