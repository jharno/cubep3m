cd ../utils/cic_power

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include

rm -f cic_power_p3dfft
rm -f ngp_power_p3dfft
rm -f cic_power_dmnu
rm -f ngp_power_dmnu

# To run with FFTW2:
#mpif77 -shared-intel -fpp -g -O3 -xhost -DSLAB -DBINARY -DPPINT -DNGP   -mt_mpi cic_power_p3dfft.f90 -o ngp_power_p3dfft  -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

# To run with p3dfft
mpif90 -shared-intel -fpp -g -O3 -DNGP -DBINARY -mt_mpi cic_power_p3dfft.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_power_p3dfft -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
#mpif90 -shared-intel -fpp -g -O3 -DNGP -mt_mpi cic_power_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_power_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_power_p3dfft.csh"

