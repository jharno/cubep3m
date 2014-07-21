#!/bin/bash

cd ../source_threads_nu

make clean
make -f Makefile_GPC_nu NEUTRINOS=0
make clean
make -f Makefile_GPC_nu

cd ../batch

echo "Compiled cubep3m dark matter and neutrino executables"

