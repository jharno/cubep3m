cd ../utils/cic_velpower

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

rm -f cic_veldivg
rm -f ngp_veldivg

mpif90 -shared-intel -fpp -g -O3 -openmp -mcmodel=medium -DNGP -DLOGBIN -DCURL -mt_mpi indexedsort.f90 cic_veldivg_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_veldivg -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_veldivg.csh" 

