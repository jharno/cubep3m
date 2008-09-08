!! gas_power.f90 : Hugh Merz Sept 9th, 2005
!! Compile with: mpif77 -fpp -g -w -O3 -axN -DBINARY gas_power.f90 -o gas_power  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

program gas_power 
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
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! parallelization variables
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
  integer(4), dimension(3) :: slab_coord
  integer(4) :: slab_rank, mpi_comm_cart, rank, ierr

  !#ifdef IA64  !!! REMEMBER TO CHANGE ON 64 BIT SYSTEMS
  !#  integer(8) :: plan, iplan
  !#else
    integer(4) :: plan, iplan
  !#endif

  logical :: firstfftw

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.14159

  !! Power spectrum arrays
  real, dimension(2,nc) :: pkg

  !! Fourier transform arrays
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
  real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work

  !! arrays
  real, dimension(5,nc_node_dim,nc_node_dim,nc_node_dim) :: u 

  !! Equivalence arrays to save memory
  equivalence (slab_work,recv_cube) 

  !! Common block
  common u,slab_work,cube,slab,pkg

!!---start main--------------------------------------------------------------!!

  call mpi_initialize
  if (rank == 0) call writeparams
  firstfftw=.true.  ! initialize fftw so that it generates the plans
  if (ic_ps) then
    call initvar
    call read_gas_ic 
    call overdensity
    if (rank == 0) call writepowerspectra
  else 
    call read_checkpoint_list
    do cur_checkpoint=1,num_checkpoints
      call initvar
      call read_gas
      call overdensity
      if (rank == 0) call writepowerspectra
    enddo
  endif
  call cp_fftw(0)
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
  subroutine overdensity
    implicit none
    integer :: i,j,k
    real    :: d,dmin,dmax,sum_gas,sum_gas_local,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart

    real time1,time2
    call cpu_time(time1)

    sum_gas_local=sum(u(1,:,:,:)) 
    call mpi_reduce(sum_gas_local,sum_gas,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) print *,'gas total mass=',sum_gas

    !! Convert gas density field to delta field
    cube=u(1,:,:,:)*omega_m/omega_b-1.0

    dmin=0
    dmax=0
    dsum=0
    dvar=0

    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim
             d=cube(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo

    call mpi_reduce(dsum,dsumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dvar,dvart,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dmin,dmint,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
    call mpi_reduce(dmax,dmaxt,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

    if (rank==0) then
      dsum=dsumt/nc**3
      dvar=sqrt(dvart/nc**3)
      write(*,*)
      write(*,*) 'gas min    ',dmint
      write(*,*) 'gas max    ',dmaxt
      write(*,*) 'Delta sum ',real(dsum)
      write(*,*) 'Delta var ',real(dvar)
      write(*,*)
    endif
 
    !! Forward FFT gas delta field
    call cp_fftw(1)

    !! Compute dm power spectrum
    call powerspectrum(slab,pkg)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called overdensity'
    return

  end subroutine overdensity
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

!! calculate coordinates within slab for cube processes

    slab_coord(3) = rank / nodes_slab
    slab_rank = rank - slab_coord(3) * nodes_slab
    slab_coord(2) = slab_rank / nodes_dim
    slab_coord(1) = slab_rank - slab_coord(2) * nodes_dim
    do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_neighbor(i,j) = i + j * nodes_dim + slab_coord(3) &
                           * nodes_slab
      enddo
    enddo

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

!!---------------------------------------------------------------------------!!

  subroutine pack_slab
!! pack cubic data into slab decomposition for fftw transform
    implicit none

    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag
      
    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status
        
    num_elements = nc_node_dim * nc_node_dim * nc_slab
                       
!! swap data           
        
    do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_slice = i + j * nodes_dim
        tag = rank**2
        rtag= slab_neighbor(i,j)**2
        call mpi_isend(cube(1,1,slab_slice*nc_slab + 1), num_elements, &
                       mpi_real, slab_neighbor(i,j), tag, mpi_comm_world, &
                       requests(slab_slice+1),ierr)
        call mpi_irecv(recv_cube(1,1,1,slab_slice), &
                       num_elements, mpi_real, slab_neighbor(i,j),rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
      enddo
    enddo
    
    call mpi_waitall(2*nodes_dim**2, requests, wait_status, ierr)

!! place data in the slab

    do j = 0, nodes_dim - 1
      j0 = j * nc_node_dim + 1
      j1 = (j + 1) * nc_node_dim
      do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        slab_slice = i + j * nodes_dim
        slab(i0:i1,j0:j1,:) = recv_cube(:,:,:,slab_slice)
      enddo
    enddo
      
  end subroutine pack_slab
    
!-------------------------------------------------------------------!

  subroutine unpack_slab
!! unpack slab data into cubic decomposition following fftw transform
    implicit none
      
    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag
    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status
      
!! place data in the recv_cube buffer
      
    do j = 0, nodes_dim - 1
      j0 = j * nc_node_dim + 1
      j1 = (j + 1) * nc_node_dim
      do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        slab_slice = i + j * nodes_dim
        recv_cube(:,:,:,slab_slice) = slab(i0:i1,j0:j1,:)
      enddo
    enddo

    num_elements = nc_node_dim * nc_node_dim * nc_slab

!! swap data

   do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_slice = i + j * nodes_dim
        tag  = rank**2
        rtag = slab_neighbor(i,j)**2
        call mpi_isend(recv_cube(1,1,1,slab_slice), num_elements, &
                       mpi_real, slab_neighbor(i,j), tag, mpi_comm_world, &
                       requests(slab_slice+1),ierr)
        call mpi_irecv(cube(1,1,slab_slice * nc_slab +1), &
                       num_elements, mpi_real, slab_neighbor(i,j), rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
      enddo
    enddo

    call mpi_waitall(2*nodes_dim**2,requests, wait_status, ierr)

  end subroutine unpack_slab

!-------------------------------------------------------------------!

  subroutine cp_fftw(command)
!! calculate fftw transform
!! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
    implicit none
!#ifdef IA64
!    include '/scratch/merz/fftw/include/fftw_f77.i'
!#else
    include '/home/merz/include/fftw_f77.i'
!#endif

    integer(4), parameter :: order=FFTW_NORMAL_ORDER ! FFTW_TRANSPOSED_ORDER

    integer(4) :: i
    integer(4) :: command

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'starting fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

! initialize plan variables for fftw

    if (firstfftw) then
      call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nc, &
            nc,nc, FFTW_REAL_TO_COMPLEX, FFTW_MEASURE)
      call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nc, &
            nc,nc, FFTW_COMPLEX_TO_REAL, FFTW_MEASURE)
#ifdef DEBUG_LOW
      print *,'finished initialization of fftw',rank
#endif
      firstfftw=.false.
    endif

! giver

    if (command /= 0) then

!! call pack routine if we are going forward

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'starting pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif
      if (command > 0) call pack_slab

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'finished forward slab pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    if (command > 0) then
      call rfftwnd_f77_mpi(plan,1,slab,slab_work,1,order)
    else
      call rfftwnd_f77_mpi(iplan,1,slab,slab_work,1,order)
      slab=slab/real(nc*nc*nc)
    endif

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'finished fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

!! unpack the slab data

      if (command < 0) call unpack_slab

    else

! if command = 0 we delete the plans

      call rfftwnd_f77_mpi_destroy_plan(iplan)
      call rfftwnd_f77_mpi_destroy_plan(plan)
    endif

  end subroutine cp_fftw

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

  subroutine writepowerspectra
    implicit none
    integer      :: k
    real         :: kr
    character*80 :: fn
    character*7  :: z_write
    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is dm p(k)
    !! 3rd is standard deviation

    if (ic_ps) then
      fn=output_path//'ic_gasps.dat'
    else
      write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
      z_write=adjustl(z_write)
      fn=output_path//z_write(1:len_trim(z_write))//'gasps.dat'
    endif
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       write(11,*) kr,pkg(:,k)
    enddo
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write power spectra'
    return
  end subroutine writepowerspectra

!!------------------------------------------------------------------!!

  subroutine powerspectrum(delta,pk)
    implicit none
    real, dimension(2,nc)       :: pk
    real, dimension(nc+2,nc,nc_slab) :: delta

    integer :: i,j,k,kg
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow
    real, dimension(3,nc,nc_slab) :: pkt
    real, dimension(3,nc) :: pktsum

    real time1,time2
    call cpu_time(time1)

    pkt=0.0
    pktsum=0.0
    !! Compute power spectrum
    do k=1,nc_slab
       kg=k+nc_slab*rank
       if (kg .lt. hc+2) then
          kz=kg-1
       else
          kz=kg-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          do i=1,nc+2,2
             kx=(i-1)/2
             kr=sqrt(kx**2+ky**2+kz**2)
             if (kr .ne. 0) then
                k1=ceiling(kr)
                k2=k1+1
                w1=k1-kr
                w2=1-w1
                pow=sum((delta(i:i+1,j,k)/ncr**3)**2)
                pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                pkt(3,k1,k)=pkt(3,k1,k)+w1
                pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                pkt(3,k2,k)=pkt(3,k2,k)+w2
             endif
          enddo
       enddo
    enddo

    !! Merge power spectrum from threads
    do k=2,nc_slab
       pkt(:,:,1)=pkt(:,:,1)+pkt(:,:,k)
    enddo

    !! Reduce to rank 0
    call mpi_reduce(pkt(:,:,1),pktsum,3*nc,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

    !! Divide by weights
    !! pk(1,k) stores pk(k)
    !! pk(2,k) stores standard deviation
    if (rank == 0) then
      do k=1,nc
        if (pktsum(3,k) .eq. 0) then
          pk(:,k)=0
        else
          pk(:,k)=pktsum(1:2,k)/pktsum(3,k)
          pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))
          pk(1:2,k)=4*pi*(k-1)**3*pk(1:2,k)
       endif
      enddo
    endif

    call mpi_bcast(pk,2*nc,mpi_real,0,mpi_comm_world,ierr)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called power spectrum'
    return
  end subroutine powerspectrum

!!------------------------------------------------------------------!!

  subroutine initvar
    implicit none
    integer :: k

    real time1,time2
    call cpu_time(time1)

    do k=1,nc_slab
       slab_work(:,:,k)=0
    enddo
    do k=1,nc_node_dim
       u(:,:,:,k)=0
    enddo
    do k=1,nc_node_dim
       cube(:,:,k)=0
    enddo
    do k=1,nc_slab
       slab(:,:,k)=0
    enddo
    do k=1,nc
       pkg(:,k)=0
    enddo
    
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

end program gas_power 
