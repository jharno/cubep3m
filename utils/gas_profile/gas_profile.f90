!! gas_profile.f90 : Hugh Merz Sept 9th, 2005
!! Compile with: mpif77 -fpp -g -w -O3 -axN -DBINARY gas_profile.f90 -o gas_profile 

program gas_profile
  implicit none
  include 'mpif.h'

! frequently changed parameters are found in this header file:
  include '../../parameters'

  logical, parameter :: ic_ps=.false. !.true.

  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints'

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! internals
  integer, parameter :: max_checkpoints=100
  real, dimension(max_checkpoints) :: z_checkpoint
  integer num_checkpoints, cur_checkpoint

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim

  integer(4) ::  mpi_comm_cart, cart_rank, rank, ierr, cart_coords(3)

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.14159

  !! arrays
  real, dimension(5,nc_node_dim,nc_node_dim,nc_node_dim) :: u 
  real, dimension(5,nc) :: u_profile,u_profileg

  !! Common block
  common u,u_profile,u_profileg

!!---start main--------------------------------------------------------------!!

  call mpi_initialize
  if (rank == 0) call writeparams
  if (ic_ps) then
    call initvar
    call read_gas_ic 
    call profile
    if (rank == 0) call writeprofile
  else 
    call read_checkpoint_list
    do cur_checkpoint=1,num_checkpoints
      call initvar
      call read_gas
      call profile
      if (rank == 0) call writeprofile
    enddo
  endif
  call mpi_finalize(ierr)

contains

!!---------------------------------------------------------------------------!!
  subroutine read_gas_ic
    implicit none

      character(len=80) :: fn1,fn2,fn3
      integer :: fstat,ierr
   
      if (rank==0) print *,'reading in gas initial conditions'
      fn1=ic_path//'mhd_ic'
      write(fn2,'(i10,".dat")'),rank
      fn2=adjustl(fn2)
      fn3=fn1(1:len_trim(fn1))//fn2(1:len_trim(fn2))
#ifdef BINARY
      open(186,file=fn3,form='binary',status='old',iostat=fstat)
#else
      open(186,file=fn3,form='unformatted',status='old',iostat=fstat)
#endif
      if (fstat /= 0) then
        print *,'error opening gas initial conditions'
        print *,'rank',rank,'file:',fn3
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(186) u
      close(186)

  end subroutine read_gas_ic

!!---------------------------------------------------------------------------!!
  subroutine profile 
    implicit none
    integer :: i,j,k,ig,jg,kg

    real time1,time2
    call cpu_time(time1)

    do k=1,nc_node_dim
       kg=k+nc_node_dim*cart_coords(1)
       do j=1,nc_node_dim
          jg=j+nc_node_dim*cart_coords(2)
          if (jg /= kg) cycle
  !        print *,kg,jg,cart_rank
          do i=1,nc_node_dim
             ig=i+nc_node_dim*cart_coords(3)
             if (ig==jg) then
!               print *,ig,cart_rank,u(1,i,j,k)
               u_profile(:,ig)=u(:,i,j,k)
             endif
          enddo
       enddo
    enddo

    call mpi_reduce(u_profile,u_profileg,5*nc,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

    if (rank==0) u_profile=u_profileg 
 
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called profile'
    return

  end subroutine profile 
!!---------------------------------------------------------------------------!!

  subroutine mpi_initialize
    implicit none
    
    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder
  
!! set up global mpi communicator

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (nodes_returned /= nodes ) then
      write(*,*) 'dist_init compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'dist_init nodes=',nodes 
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    if (mod(nc,nodes) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into slabs'
      write(*,*) 'nc=',nc,'nodes=',nodes,'mod(nc,nodes) != 0'
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    if (rank==0) then
      write(*,*) 'dist_init running on',nodes,'nodes'
      write(*,*) 'using cubic distribution:',nodes_dim,'nodes per dimension'
      write(*,*) nc,'cells in mesh'
    endif

!! create cartesian communicator based on cubic decomposition

    dims(:) = nodes_dim
    periodic(:) = .true.
    reorder = .false.
    ndim = 3

    call mpi_cart_create(mpi_comm_world, ndim,dims, periodic, &
                       reorder, mpi_comm_cart, ierr)
    call mpi_comm_rank(mpi_comm_cart, cart_rank, ierr)
    call mpi_cart_coords(mpi_comm_cart, cart_rank, ndim,  &
                         cart_coords, ierr)

    print *,rank,cart_rank,cart_coords

  end subroutine mpi_initialize

!!---------------------------------------------------------------------------!!

  subroutine read_checkpoint_list
!! read in list of checkpoints to calculate spectra for
    implicit none

    integer :: i,fstat

    if (rank == 0) then
      open(11,file=checkpoints,status='old',iostat=fstat)
      if (fstat /= 0) then
        print *,'error opening checkpoint list file'
        print *,'rank',rank,'file:',checkpoints
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      do num_checkpoints=1,max_checkpoints
        read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
      enddo
  41  num_checkpoints=num_checkpoints-1
  51  close(11)
      print *,'checkpoints to recompose:'
      do i=1,num_checkpoints
        write(*,'(f5.1)') z_checkpoint(i)
      enddo
    endif

    call mpi_bcast(num_checkpoints,1,mpi_integer,0,mpi_comm_world,ierr)

  end subroutine read_checkpoint_list

!!---------------------------------------------------------------------------!!

  subroutine read_gas
    implicit none
    
    real z_write
    integer j,fstat
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name

      type comm_wld
        integer :: g !global number of zones in one direction
        integer :: r !global index of index 0 in local array
        integer :: m,n !start and end of local section without buffers
        integer :: l !dimension of array inclusive buffers
        integer, dimension(4) :: requests !communication handles
      end type comm_wld


    !! these are unnecessary headers from the checkpoint
    real(4) :: cur_t 
    integer(4) :: cur_iter

    type(comm_wld) ::  nxx,nyy,nzz 

!! generate checkpoint names on each node
    if (rank==0) then
      z_write = z_checkpoint(cur_checkpoint)
      print *,'calculating spectrum for z=',z_write
    endif

    call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

    write(z_string,'(f7.3)') z_write
    z_string=adjustl(z_string)

    write(rank_string,'(i4)') rank
    rank_string=adjustl(rank_string)

    check_name=output_path//z_string(1:len_trim(z_string))//'mhd'// &
               rank_string(1:len_trim(rank_string))//'.dat'

!! open checkpoint    
#ifdef BINARY
    open(unit=21,file=check_name,status='old',iostat=fstat,form='binary')
#else
    open(unit=21,file=check_name,status='old',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint'
      write(*,*) 'rank',rank,'file:',check_name
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! read in checkpoint header data
     read(21) cur_iter,cur_t,nxx,nyy,nzz
!! read in gas array   
    read(21) u 
    close(21)
 
  end subroutine read_gas

!-------------------------------------------------------------------!

  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nodes   ', nodes
    write(*,*) 'nc      ', nc
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams

!!------------------------------------------------------------------!!

  subroutine writeprofile
    implicit none
    integer      :: k
    character*80 :: fn
    character*7  :: z_write
    real time1,time2
    call cpu_time(time1)

    !! Output gas profile 

    if (ic_ps) then
      fn=output_path//'ic_gas_profile.dat'
    else
      write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
      z_write=adjustl(z_write)
      fn=output_path//z_write(1:len_trim(z_write))//'gas_profile.dat'
    endif
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=1,nc
       write(11,*) u_profile(:,k)
    enddo
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write gas profile'
    return
  end subroutine writeprofile

!!------------------------------------------------------------------!!

  subroutine initvar
    implicit none
    integer :: k

    real time1,time2
    call cpu_time(time1)

    do k=1,nc_node_dim
       u(:,:,:,k)=0
    enddo
    do k=1,nc
       u_profile(:,k)=0.0
    enddo
    
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

end program gas_profile
