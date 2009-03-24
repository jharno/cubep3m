# Load proper modules

#module clear
module unload fftw
module load intel #/intel-10.0b.017                                     
module load lam #/lam-7.1.3-intel
module load fftw #/2.1.5-intel10

module list

source ./init_pp_threads.csh
#source ./init_pp_threads_minimal.csh

#cd ../batch/

echo "Ready to run" 
