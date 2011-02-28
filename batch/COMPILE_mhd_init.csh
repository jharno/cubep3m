#module unload fftw
#module purge
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw/2.1.5-intel-intelmpi3 #/2.1.5-intel10
#module load intelmpi

#module load extras

cd ../utils/mhd_init
rm -f mhd_init
mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -mt_mpi mhd_init.f90 -o mhd_init  -lm -ldl 
cd ../../batch
