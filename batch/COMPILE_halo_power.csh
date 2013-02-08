# Load proper modules

#module clear
#module unload fftw
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw #/2.1.5-intel10

# Could add -DKAISER or -DNGP

cd ../utils/cic_power

rm -f halo_power

#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT -mt_mpi cic_power.f90 -o cic_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DPPINT -DNGP halo_power.f90 -o halo_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DPPINT -DNGP cic_power_halo_corr.f90 -o cic_power_halo  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT -DNGP -mt_mpi halo_power.f90 -o halo_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DDEBUG  -mt_mpi cic_init_power.f90 -o cic_init_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DDEBUG -DNGP -mt_mpi cic_init_power.f90 -o ngp_init_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
cd ../../batch/

echo "Ready to run" 
