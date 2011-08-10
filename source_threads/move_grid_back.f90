!! move grid back
  subroutine move_grid_back
    implicit none
#ifdef DISP_MESH 
    include 'mpif.h'
#endif
    include 'cubepm.fh'

    integer(4) :: i

#ifdef DISP_MESH 
 
    if (rank==0) then      
       print *,'substracting  shake offset of :',shake_offset
    endif

    call mpi_bcast(shake_offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call system_clock(count=count_i)
    
    !$omp parallel do default(shared) private(i)
    do i=1,np_local
      xv(1:3,i)=xv(1:3,i)-shake_offset(:)
    enddo
    !$omp end parallel do
    shake_offset = 0.0

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('move_grid_back',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'move grid finished',real(count_f-count_i)/real(count_r)
#endif

#else
    if(rank==0) write(*,*)  '*** Could not move back, no off set to start with! ***'
#endif

  end subroutine move_grid_back
