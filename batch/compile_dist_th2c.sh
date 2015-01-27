cd ../utils/dist_init

# Library links
P3DFFT_LIB=/HOME/bnu_ztj_1/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/HOME/bnu_ztj_1/lib/p3dfft_2.5.1/include
MKL_FFTW_INC=/opt/intel/composer_xe_2013_sp1.2.144/mkl/include/fftw
MKL_FFTW_LIB=/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/

rm -f dist_init_dmnu_dm
rm -f dist_init_dmnu_nu

echo
echo "compiling ../utils/dist_init/dist_init_dmnu_dm"
#mpif90 -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium dist_init_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_dmnu_dm -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 
mpif90 -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium  dist_init_dmnu.f90 -DZIP -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_dmnu_dm -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 

echo
echo "done"
echo

echo "compiling ../utils/dist_init/dist_init_dmnu_nu"
#mpif90 -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNEUTRINOS -DVELTRANSFER dist_init_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_dmnu_nu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64
mpif90 -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium         -DNEUTRINOS -DVELTRANSFER dist_init_dmnu.f90 -DZIP -I$P3DFFT_INC -I$MKL_FFTW_INC -o dist_init_dmnu_nu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64 
echo
echo "done"
echo
# the largest array that can be set as a private variable in a threaded region is 1000 MB
#export OMP_STACKSIZE="1000M"
#echo "export OMP_STACKSIZE='1000M' done"
echo

cd ../../batch
