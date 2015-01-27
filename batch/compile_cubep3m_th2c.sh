##!/bin/csh
## p3dfft, use MKL, offload enabled

cd ../source_threads

make clean;
make -f Makefile_th2c  -j NEUTRINOS=0
echo
echo "compiled dark matter executable ../source_threads/cubep3m_dm"
echo
make clean;
make -f Makefile_th2c  -j
echo
echo "compiled dark matter plus neutrino executable ../source_threads/cubep3m_nu"
echo

cd ../batch

echo

#export LD_DYNAMIC_WEAK=1
