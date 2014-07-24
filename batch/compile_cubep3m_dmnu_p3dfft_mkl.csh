##!/bin/csh
## p3dfft, use MKL, offload enabled

cd ../source_threads

make clean;
make -j NEUTRINOS=0
echo
echo "compiled dark matter executable ../source_threads/cubep3m_dm"
echo
make clean;
make -j
echo
echo "compiled dark matter plus neutrino executable ../source_threads/cubep3m_nu"
echo

cd ../batch

export OFFLOAD_REPORT=3

echo
echo "Sourced COMPILE_cubep3m_dmnu_p3dfft_mkl.csh"
echo

#export LD_DYNAMIC_WEAK=1
