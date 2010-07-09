# Load proper modules

#module clear
#module unload fftw
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw #/2.1.5-intel10

cd ../utils/halo_merge
#parallel:
#mpif77 -shared-intel -fpp -DBINARY indexedsort.f90 halo_merge.f90 -o halo_merge 
#serial:
ifort -shared-intel -fpp -DPID_FLAG -CB indexedsort.f90 serial_halo_merge.f90 -o s_halo_merge    
cd ../../batch/

echo "Ready to run" 
