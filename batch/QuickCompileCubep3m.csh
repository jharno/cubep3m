cd ../source_threads/
make -j8 cubepm -f Make_P5M_THREADS
cd ../batch/

module unload intel
module unload lam
module unload fftw

module load intel
module load lam
module load fftw

lamboot
