#!/bin/csh 
#PBS -l nodes=1:ppn=8
#PBS -q workq 
#PBS -r n
##PBS -l other=scratch-2week
#PBS -l walltime=10:00:00
#PBS -N Compile1000Binaries

#cd $PBS_O_WORKDIR

# Make sure parameters starts at the right output path  

#module unload fftw
#module load fftw



#set Run = '1'
#set MaxRun = '30'
#set Version = '1'
set MinVersion = '1'
set MaxVersion = '5'


sed '3 s/Run-1/Run-'"$MinVersion"'/' ../parameters -i  # Comment to run on same seeds
sed '7 s/Run-1/Run-'"$MinVersion"'/' ../parameters -i 
sed '11 s/Run-1/Run-'"$MinVersion"'/' ../parameters -i 

#sed '3 s/LBG1/LBG'"$MinVersion"'/' ../parameters -i  # Comment to run on same seeds
#sed '7 s/LBG1/LBG'"$MinVersion"'/' ../parameters -i 
#sed '11 s/LBG1/LBG'"$MinVersion"'/' ../parameters -i 


#while($Run <= $MaxRun)

#    @ NextRun = ($Run + 1)

    #set Version = '1'
    set Version = $MinVersion    

    while($Version <= $MaxVersion)
	
	@ NextVersion = ($Version + 1)
	if($NextVersion > $MaxVersion)then
	    @ NextVersion = ($NextVersion - $MaxVersion)
	endif

	set suffix = $Version''

	echo cubep3m_$suffix

	#source COMPILE.csh	
	source COMPILE_cubep3m.csh
	source COMPILE_cic_power.csh
	source COMPILE_dist_init.csh
	#source COMPILE_pgm_proj.csh
	#source COMPILE_halo_power.csh

	#cd ../utils/halo_merge
	#rm -f halo_merge
	#mpif77 -shared-intel -fpp -DBINARY indexedsort.f90 halo_merge.f90 -o halo_merge

	mv ../source_threads/cubep3m ../source_threads/cubep3m_$suffix 
	mv ../utils/dist_init/dist_init ../utils/dist_init/dist_init_$suffix
	#mv ../utils/halo_merge/s_halo_merge ../utils/halo_merge/s_halo_merge_$suffix	
	#mv ../utils/cic_power/cic_power ../utils/cic_power/cic_power_$suffix
        mv ../utils/cic_power/ngp_power ../utils/cic_power/ngp_power_$suffix
	#mv ../utils/cic_power/cic_init_power ../utils/cic_power/cic_init_power_$suffix
	#mv ../utils/cic_power/ngp_init_power ../utils/cic_power/ngp_init_power_$suffix
	#mv ../utils/cic_power/halo_power ../utils/cic_power/halo_power_$suffix
	#cd ../../batch

	#grep pp_range ../source_threads/cubepm.par
	#sed '88 s/pp_range = '"$Version"'/pp_range = '"$NextVersion"'/' ../source_threads/cubepm.par -i 


	#grep V ../parameters
        grep Run ../parameters  
	sed '3 s/Run-'"$Version"'/Run-'"$NextVersion"'/' ../parameters -i  
	sed '7 s/Run-'"$Version"'/Run-'"$NextVersion"'/' ../parameters -i 
	sed '11 s/Run-'"$Version"'/Run-'"$NextVersion"'/' ../parameters -i 
	#sed '3 s/LBG'"$Version"'/LBG'"$NextVersion"'/' ../parameters -i  
	#sed '7 s/LBG'"$Version"'/LBG'"$NextVersion"'/' ../parameters -i 
	#sed '11 s/LBG'"$Version"'/LBG'"$NextVersion"'/' ../parameters -i 

	@ Version = ($Version + 1)

    end

#    sed '3 s/RUN-'"$Run"'/RUN-'"$NextRun"'/' ../parameters -i  # Comment to run on same seeds
#    sed '7 s/RUN-'"$Run"'/RUN-'"$NextRun"'/' ../parameters -i 
#    sed '11 s/RUN-'"$Run"'/RUN-'"$NextRun"'/' ../parameters -i 


#    @ Run = ($Run + 1)

#end

@ NextVersion = (1)    
#@ NextRun = (1)    

sed '3 s/Run-'"$Version"'/Run-'"$NextVersion"'/' ../parameters -i  # Comment to run on same seeds
sed '7 s/Run-'"$Version"'/Run-'"$NextVersion"'/' ../parameters -i 
sed '11 s/Run-'"$Version"'/Run-'"$NextVersion"'/' ../parameters -i 

#sed '3 s/LBG'"$Version"'/LBG'"$NextVersion"'/' ../parameters -i  # Comment to run on same seeds
#sed '7 s/LBG'"$Version"'/LBG'"$NextVersion"'/' ../parameters -i 
#sed '11 s/LBG'"$Version"'/LBG'"$NextVersion"'/' ../parameters -i 
