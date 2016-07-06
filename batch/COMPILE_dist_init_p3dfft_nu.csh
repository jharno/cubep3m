cd ../utils/dist_init

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

rm -f dist_init_dmnu_dm
rm -f dist_init_dmnu_nu

mpif90 -mcmodel=medium -shared-intel -fpp -g -O3 -xhost -mt_mpi -openmp dist_init_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o dist_init_dmnu_dm -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
mpif90 -mcmodel=medium -shared-intel -fpp -g -O3 -xhost -mt_mpi -openmp -DNEUTRINOS -DVELTRANSFER dist_init_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o dist_init_dmnu_nu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch

echo "Sourced dist_init_p3dfft.csh"

