#!/bin/bash
      #Unset all MPI Affinities
export       MV2_USE_AFFINITY=0
export    MV2_ENABLE_AFFINITY=0
export    VIADEV_USE_AFFINITY=0
export VIADEV_ENABLE_AFFINITY=0
                  #TasksPerNode
TPN=`echo $PE | sed 's/way//'`
[ ! $TPN ] && echo TPN NOT defined!
[ ! $TPN ] && exit 1
        #Get which ever MPI's rank is set     
[ "x$PMI_RANK" != "x" ] && RANK=$PMI_RANK
[ "x$MPI_RANK" != "x" ] && RANK=$MPI_RANK
[ "x$OMPI_MCA_ns_nds_vpid" != "x" ] \
            && RANK=$OMPI_MCA_ns_nds_vpid

socket=$(( $RANK % $TPN  ))

numactl -C $socket -i all ./a.out
