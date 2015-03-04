use omp_lib
include 'mpif.h'
#include "cubepm.fh"
call mpi_initialize
call variable_initialize
call particle_initialize
call link_list
call particle_pass
call link_list
call halofind
call mpi_barrier(mpi_comm_world,ierr)
if (rank==0) print*, 'halofind done'
call mpi_finalize(ierr)
end
