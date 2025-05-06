# Change to the directory from which you originally submitted this job.
cd $SLURM_SUBMIT_DIR

salloc -A pbermel --time=5-00:00:00 --ntasks=10

conda activate pmp

mpirun -np=10 python straight_waveguide_recreate_rough.py -df=0.1 -kz=1.35 -avg_r_bto_air=0.05 -std_r_bto_air=0.01 -avg_r_bto_sio2=0.05 -std_r_bto_sio2=0.01 | tee bump_rad_0.05_50_bumps_rough_fluxs_kz1.35.out

conda deactivate

echo "Done"
