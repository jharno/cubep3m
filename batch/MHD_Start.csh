let nodes='1'

module load intel
module load lam/lam-7.1.4-intel-11.1.072
module load fftw/2.1.5-intel-11.1.072-lam 


lamboot

source COMPILE_dist_init.csh
mpirun -np $nodes ./dist_init
source COMPILE_mhd_init.csh
mpirun -np $nodes ./mhd_init
#cd ../source_threads
source COMPILE_cubep3m_mhd.csh
time mpirun -np $nodes ./cubep3m
#cd ../batch
source MHD_ANALYZE.csh

lamhalt
