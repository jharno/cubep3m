# Load proper modules

#module clear
#module unload fftw
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw #/2.1.5-intel10

# Could add -DKAISER or -DNGP

cd ../utils/cic_power

rm -f cic_pow
rm -f ngp_pow
rm -f cic_init_power
rm -f ngp_init_power

#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT cic_power.f90 -o cic_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DPPINT -DNGP cic_power.f90 -o ngp_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DPPINT -DNGP cic_power.f90 -o ngp_power  -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DDEBUG cic_init_power.f90 -o cic_init_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DDEBUG -DNGP cic_init_power.f90 -o ngp_init_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DDEBUG -DNGP cic_init_power.f90 -o ngp_init_power  -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
cd ../../batch/

echo "Ready to run" 
