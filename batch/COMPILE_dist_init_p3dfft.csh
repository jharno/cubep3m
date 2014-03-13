cd ../utils/dist_init_ng

# Load proper modules
#module purge
#module load intel intelmpi fftw/3.3.0-intel-impi

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

rm -f dist_init_p3dfft

mpif77 -shared-intel -fpp -g -O3 -xhost -Dwrite_den -DBINARY -mt_mpi dist_init_ng_p3dfft.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o dist_init_p3dfft -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f  

cd ../../batch

echo "Sourced dist_init_p3dfft.csh"

