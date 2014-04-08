cd ../utils/cic_velpower

# Library links
P3DFFT_LIB=/home/p/pen/emberson/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/home/p/pen/emberson/lib/p3dfft_2.5.1/include

rm -f cic_velpower
rm -f ngp_velpower
rm -f cic_velpower_dmnu
rm -f ngp_velpower_dmnu

# NOTE: ON THE BGQ NEED TO CHANGE -openmp TO -qsmp=omp

mpif90 -shared-intel -fpp -g -O3 -openmp -DNGP -DDIAG -mt_mpi indexedsort.f90 cic_velpower.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_velpower -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
#mpif90 -shared-intel -fpp -g -O3 -DNGP -DGAUSSIAN_SMOOTH -DMOMENTUM -mt_mpi cic_velpower_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_velpower_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
mpif90 -shared-intel -fpp -g -O3 -openmp -DNGP -mt_mpi indexedsort.f90 cic_velpower_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o ngp_velpower_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_velpower_p3dfft.csh"

