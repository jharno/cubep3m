cd ../utils/cic_power

# Library links
P3DFFT_LIB=/HOME/bnu_ztj_1/lib/p3dfft_2.5.1/lib
P3DFFT_INC=/HOME/bnu_ztj_1/lib/p3dfft_2.5.1/include
MKL_FFTW_INC=/opt/intel/composer_xe_2013_sp1.2.144/mkl/include/fftw
MKL_FFTW_LIB=/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/


rm -f ngp_power_dmnu

echo "compile ngp_power_dmnu ..."
mpif90 -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DZIP -DZIPDM -DNGP -DGROUPS -DLOGBIN -DSPLITHALOS -Dwrite_den ../cic_velpower/indexedsort.f90 cic_power_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_power_dmnu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64

echo "done"

cd ../../batch/
