#module unload fftw
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw #/2.1.5-intel10

#module load extras

cd ../utils/dist_init
rm -f dist_init
mpif77 -shared-intel -fpp -CB -fpe0 -g -O3 -xhost -DBINARY -mt_mpi dist_init_dm.f90 -o dist_init  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl #-DMY_TRANSFER
#mpif77 -shared-intel -fpp -CB -fpe0 -g -O0 -xhost -DBINARY -mt_mpi dist_init_dm.f90 -o dist_init -L/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib -I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include  -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl 


#mpxlf_r -fpp -q64 -O3 -qhot -qarch=pwr6 -qtune=pwr6 -bdatapsize:64k -bstacksize:64k -DBINARY -DPPINT dist_init_dm.f90 -o dist_init -L$SCINET_EXTRAS_LIB -I$SCINET_EXTRAS_INC
cd ../../batch
