!! update particle positions
  subroutine update_position
    use omp_lib
#ifdef FFTMKL
   use MKL_DFTI
#endif
implicit none
#ifdef DISP_MESH 
    include 'mpif.h'
#endif
#    include "cubepm.fh"

#ifdef READ_SEED
    character(len=max_path) :: seedfile
    integer(4) :: seedsize
    integer(4), allocatable, dimension(:) :: iseed!,old
#endif
#ifdef READOFFSET
  logical :: openoffset=.true.

  integer(4) :: ipass
#endif
!    integer k_seed, clock
!    integer, dimension(8) :: values

    integer(4) :: i,j
#ifdef DISP_MESH 
    real(4), dimension(3) :: offset

if (rank==0) then

       
# ifdef READ_SEED

       !This will always use the same random number at each time step. 
       !It surely introduces a bias, but is good for testing code.        
       call random_seed
       call random_seed(size=seedsize)

       allocate(iseed(seedsize))
       !allocate(old(seedsize))


       seedfile = ic_path//'node0/seed0.init'  !or any other seed files available
       open(11,file=seedfile)
       write(*,*) 'opened ',seedfile
       do i = 1,seedsize 
          read(11,*) j,iseed(i)
       enddo
       close(11)

       call random_seed(put=iseed(1:seedsize))
# else

       !call date_and_time(values=values)
       !if(rank==0) write(*,*) values(7:8), iseed      
       !call random_seed(put=values(7:8)) 

# endif


# ifndef READOFFSET
  call random_number(offset)
# else
  if (openoffset) then
    open(91,file='offsets.dat',status='old',access='stream')
    openoffset=.false.
    print*,'opened offset file: offsets.dat'
    ipass=0
    do ipass=1,nts+cur_checkpoint+10
      read(91) offset
      print*, 'past timesteps',ipass,'offset=', offset
    enddo
  endif
  read(91,end=999) offset; goto 999
  998 print*,'offsets.dat too short, need more random numbers.'; call mpi_abort(mpi_comm_world,ierr,ierr)
  999 print*,'update_position: get offset(1:3) =',offset
# endif

# ifdef NEUTRINOS
     offset=(offset-0.5)*mesh_scale ! no shake offset
# else
     offset=(offset-0.5)*mesh_scale*4.0  - shake_offset
# endif
     shake_offset=shake_offset+offset
     print*,'current shake offset:',shake_offset

endif ! rank==0

if (pair_infall_no_shake.and.pair_infall .or. pp_test) offset=0.0

call mpi_bcast(offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(shake_offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)

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
