cd ../utils/dist_init_ng

# Load proper modules
module purge
module load intel intelmpi 

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include
MKL_FFTW_INC=/home/p/pen/haoran/mkl_fftw/include/fftw
MKL_FFTW_LIB=/home/p/pen/haoran/mkl_fftw/lib/intel64/

rm -f dist_init_p3dfft

mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -Dwrite_den -DBINARY -mt_mpi dist_init_ng_p3dfft.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_p3dfft -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 

cd ../../batch

echo "Sourced dist_init_p3dfft.csh"

