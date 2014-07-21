!! calculate the maximum dt based on the coarse mesh force
  subroutine coarse_max_dt 
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    integer(kind=4) :: i,j,k
    real(kind=4) :: force,max_force

    call system_clock(count=count_i)

    max_force=0.0

    !$omp parallel do default(shared) &
    !$omp private(i,j,k,force) &
    !$omp reduction(MAX:max_force)
    do k=1,nc_node_dim
      do j=1,nc_node_dim
        do i=1,nc_node_dim
          force=sqrt(force_c(1,i,j,k)**2+force_c(2,i,j,k)**2 + &
                force_c(3,i,j,k)**2)
          if (force.gt.max_force) max_force=force
        enddo
      enddo
    enddo
    !$omp end parallel do
    if (rank == 0) write(*,*)  'Max coarse force= ', max_force

    call mpi_reduce(max_force,dt_c_acc,1,mpi_real,mpi_max,0, &
                    mpi_comm_world,ierr)
    dt_c_acc=sqrt(real(mesh_scale,kind=4)/(dt_c_acc*a_mid*G))
    call mpi_bcast(dt_c_acc,1,mpi_real,0,mpi_comm_world,ierr)

    if (rank == 0) write(*,*) 'maximum dt from coarse grid=',dt_c_acc
   
    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('c max dt',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'coarse max dt finished',real(count_f-count_i)/real(count_r)
#endif

  end subroutine coarse_max_dt 

