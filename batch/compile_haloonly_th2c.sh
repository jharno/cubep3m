cd ../source_threads
make clean

make -f Makefile_haloonly NEUTRINOS=0
cd ../batch
