let nodes='8'

source COMPILE_dist_init.csh
mpirun -np $nodes ./dist_init
source COMPILE_mhd_init.csh
mpirun -np $nodes ./mhd_init
cd ../source_threads
make clean
make cubep3m -j8 -f makefile_mhd_pp
time mpirun -np $nodes ./cubep3m
cd ../batch
