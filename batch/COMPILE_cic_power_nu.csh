cd ../utils/cic_power

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

rm -f cic_power_dmnu
rm -f ngp_power_dmnu

mpif90 -shared-intel -openmp -mcmodel=medium -fpp -g -O3 -DNGP -DLOGBIN -DGROUPS -mt_mpi cic_power_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_power_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_power_p3dfft.csh"

