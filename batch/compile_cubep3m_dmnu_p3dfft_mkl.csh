##!/bin/csh
## p3dfft, use MKL, offload enabled

cd ../source_threads

make clean;
make -j

cd ../batch
export OFFLOAD_REPORT=3
echo "Sourced COMPILE_cubep3m_dmnu_p3dfft_mkl.csh"
#export LD_DYNAMIC_WEAK=1
