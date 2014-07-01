cd ../utils/dist_init

# Library links
P3DFFT_LIB=/vol-th/home/bnu_ztj_haoran/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/vol-th/home/bnu_ztj_haoran/lib/p3dfft_2.5.1/include
MKL_FFTW_INC=/vol-th/intel/composer_xe_2013_sp1.1.106/mkl/include/fftw
MKL_FFTW_LIB=/vol-th/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64/

rm -f dist_init_dmnu_dm
rm -f dist_init_dmnu_nu

echo "compile dist_init_dmnu_dm ..."
mpif90 -shared-intel -fpp -g -O3 -xhost -openmp -i_dynamic -mkl -mcmodel=medium -DDEBUG_LOW dist_init_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_dmnu_dm -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 
echo "done"

echo "compile dist_init_dmnu_nu ..."
mpif90 -shared-intel -fpp -g -O3 -xhost -openmp -i_dynamic -mkl -mcmodel=medium -DDEBUG_LOW -DNEUTRINOS -DVELTRANSFER dist_init_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_dmnu_nu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 
echo "done"

cd ../../batch
