# Load proper modules

#module clear
#module unload fftw
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw #/2.1.5-intel10

# Could add -DKAISER or -DNGP

cd ../utils/cic_power

rm -f cic_pow_mhd
rm -f ngp_pow_mhd
#rm -f cic_init_power
#rm -f ngp_init_power



#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT cic_power_mhd.f90 -o ngp_power_mhd  -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DPPINT -DNGP cic_power_mhd.f90 -o ngp_power_mhd  -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DDEBUG cic_init_power_mhd.f90 -o cic_init_power_mhd   -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -DDEBUG -DNGP cic_init_power_mhd.f90 -o ngp_init_power_mhd  -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl




cd ../../batch/

echo "Ready to run" 
