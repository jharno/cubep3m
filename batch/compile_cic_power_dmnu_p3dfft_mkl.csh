cd ../utils/cic_power

P3DFFT_LIB=/vol-th/home/bnu_ztj_haoran/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/vol-th/home/bnu_ztj_haoran/lib/p3dfft_2.5.1/include
MKL_FFTW_INC=/vol-th/intel/composer_xe_2013_sp1.1.106/mkl/include/fftw
MKL_FFTW_LIB=/vol-th/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64

rm -f ngp_power_dmnu

echo "compile ngp_power_dmnu ..."
mpif90 -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -DGROUPS -Dwrite_den cic_power_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_power_dmnu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64

echo "done"

cd ../../batch/
