#!/bin/csh

cd ../source_threads
make clean

cd ../utils/cic_power
rm -f cic_power
rm -f cic_init_power
rm -f mem_usage

cd ../dist_init
rm -f dist_init

cd ../halo_merge
rm -f halo_merge

cd ../recompose
rm -f recompose

cd ../PSvsSim
rm -f PSvsSim

cd ../pgm_proj
rm -f pgm_proj

cd ../dm_slice_sample
rm -f slice_sample

cd ../../batch

exit 0
