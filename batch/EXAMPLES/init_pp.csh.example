#!/bin/csh
cd ../source_p3m
make clean
make -j8 -f Make_NOMHD_PP
cd ../utils/cic_power
rm -f cic_power
mpif77 -fpp -g -O3 -xT -DBINARY -DPPINT cic_power.f90 -o cic_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -fpp -g -O3 -xT -DBINARY -DDEBUG cic_init_power.f90 -o cic_init_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
cd ../dist_init
rm -f dist_init
mpif77 -fpp -g -O3 -xT -DBINARY -dyncom "rvar" dist_init_dm.f90 -o dist_init  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
cd ../halo_merge
rm -f halo_merge
mpif77 -fpp -DBINARY indexedsort.f90 halo_merge.f90 -o halo_merge
cd ../recompose
rm -f recompose
mpif77 -fpp -DBINARY -DPPINT recompose.f90 -o recompose
cd ../PSvsSim
rm -f PSvsSim
ifort -fpp PS_vs_simul_new.f90 deltac.f growth.f sigma_cobe_CMBfast.f spline.f splint.f -o PSvsSim
cd ../pgm_proj
rm -f pgm_proj
ifort -fpp -DBINARY pgm_proj.f90 -o pgm_proj
cd ../dm_slice_sample
rm -f slice_sample
ifort -fpp -DBINARY -DPPINT slice_sample.f90 -o slice_sample
cd ..
ifort mem_usage.f90 -o mem_usage
exit 0
