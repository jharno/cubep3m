!! update particle positions
  subroutine update_position
    implicit none
#ifdef DISP_MESH 
    include 'mpif.h'
#endif
    include 'cubepm.fh'

    integer(4) :: i,j
#ifdef DISP_MESH 
    real(4), dimension(3) :: offset

    if (rank==0) then
      call random_number(offset)
      offset=(offset-0.5)*mesh_scale
      shake_offset=shake_offset+offset
      print *,'current shake offset:',shake_offset
    endif
    if (pair_infall_no_shake.and.pair_infall .or. pp_test) offset=0.0
    call mpi_bcast(offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

    call system_clock(count=count_i)

    !$omp parallel do default(shared) private(i)
    do i=1,np_local
#ifdef DISP_MESH 
      xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)+offset(:)
#else
      xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)
#endif
    enddo
    !$omp end parallel do

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('pos updt',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'position update finished',real(count_f-count_i)/real(count_r)
#endif

  end subroutine update_position
