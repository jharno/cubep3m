cd ../utils/cic_power

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

rm -f cic_power_p3dfft_dana
rm -f ngp_power_p3dfft_dana

mpif90 -mcmodel=medium -shared-intel -fpp -g -O3 -DNGP -DCROSS -Dwrite_den -mt_mpi cic_power_p3dfft_dana_2.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_power_p3dfft_dana -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_power_p3dfft_dana.csh"

