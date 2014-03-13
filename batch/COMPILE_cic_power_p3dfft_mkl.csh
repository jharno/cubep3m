cd ../utils/cic_power

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include
MKL_FFTW_INC=/home/p/pen/haoran/mkl_fftw/include/fftw
MKL_FFTW_LIB=/home/p/pen/haoran/mkl_fftw/lib/intel64/

rm -f cic_power_p3dfft
rm -f ngp_power_p3dfft
rm -f cic_power_dmnu
rm -f ngp_power_dmnu

mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -mt_mpi cic_power_p3dfft.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_power_p3dfft -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 
#mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -mt_mpi cic_power_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_power_dmnu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64

cd ../../batch/

echo "Sourced COMPILE_cic_power_p3dfft.csh"

