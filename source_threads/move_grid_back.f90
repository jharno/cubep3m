!! move grid back
  subroutine move_grid_back
    implicit none
#ifdef DISP_MESH 
    include 'mpif.h'
#endif
    include 'cubepm.fh'

!#ifdef READ_SEED
!    character(len=max_path) :: seedfile
!    integer(4) :: seedsize
!    integer(4), allocatable, dimension(:) :: iseed
!#endif

    integer(4) :: i!,j
!#ifdef DISP_MESH 
!    real(4), dimension(3) :: current_offset

    if (rank==0) then
       
!This will always use the same random number at each time step. 
!It surely introduces a bias, but is good for testing code. 

!#ifdef READ_SEED
!       
!       call random_seed
!       call random_seed(size=seedsize)
!       allocate(iseed(seedsize))
!
!       seedfile = ic_path//'seed0.init'  !or any other seed files available
       !open(11,file=seedfile)
       !write(*,*) 'opened ',seedfile
       !do i = 1,seedsize 
         ! read(11,*) j,iseed(i)
       !enddo
       !close(11)

       !call random_seed(put=iseed(1:seedsize))

!#endif

!       call random_number(offset)
!       offset=(offset-0.5)*mesh_scale
!       shake_offset=shake_offset+offset
       print *,'substracting  shake offset of :',shake_offset
    endif
!    if (pair_infall_no_shake.and.pair_infall .or. pp_test) offset=0.0
!    call mpi_bcast(offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!#endif

    call mpi_bcast(shake_offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call system_clock(count=count_i)
    
    !$omp parallel do default(shared) private(i)
    do i=1,np_local
!#ifdef DISP_MESH 
!      xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)+offset(:)
      xv(1:3,i)=xv(1:3,i)-shake_offset(:)
!#else
!      xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)
!#endif
    enddo
    !$omp end parallel do
    shake_offset = 0.0

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('move_grid_back',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'move grid finished',real(count_f-count_i)/real(count_r)
#endif

  end subroutine move_grid_back
