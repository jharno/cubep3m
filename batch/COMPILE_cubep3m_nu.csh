#!/bin/bash

cd ../source_threads_nu

make clean
make -f Makefile_GPC_nu NEUTRINOS=0
echo
echo "COMPILED DARK MATTER EXECUTABLE"
echo
make clean
make -f Makefile_GPC_nu
echo
echo "COMPILED DARK MATTER PLUS NEUTRINO EXECUTABLE"
echo

cd ../batch

