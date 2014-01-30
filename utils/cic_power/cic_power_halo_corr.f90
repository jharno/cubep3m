!! cic_power.f90 Parallelized: Hugh Merz Jun 15, 2005
!! modified by Vincent Desjacques to calculate halo PS: May 1st, 2008
!! Compile with: mpif77 -fpp -g -w -O3 -mcmodel=medium -shared-intel -DDEBUG_LOW cic_power.f90 -o cic_power -L/share/home/00506/merzh/lib/fftw-2.1.5/lib -I/share/home/00506/merzh/lib/fftw-2.1.5/include -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

!!old: mpif77 -fpp -g -w -O3 -axN cic_power.f90 -o cic_power  -L/home/merz/lib/fftw-2.1.5_intel8/lib -I/home/merz/lib/fftw-2.1.5_intel8/include -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

program cic_power_halo
  implicit none
  include 'mpif.h'

! frequently changed parameters are found in this header file:
  include '../../parameters'
!  include '../parameterfiles/ics.param'

  logical, parameter :: correct_kernel=.false.

  character(len=*), parameter :: halofinds=cubepm_root//'/input/checkpoints_ppext2'
!  character(len=*), parameter :: halofinds='../parameterfiles/checkpoints_halos'
!  character(len=*), parameter :: halofinds=cubepm_root//'/input/halofinds'

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc
  real, parameter    :: npr=np

  !! internals
  integer, parameter :: max_checkpoints=100
  real, dimension(max_checkpoints) :: z_checkpoint
  integer num_checkpoints, cur_checkpoint
  integer, parameter :: num_massbin=4
  integer            :: cur_massbin, ntotal_massbin(num_massbin), nlocal_massbin(num_massbin)
  real, parameter    :: min_mass=160. ! 20 particles in fine grid unit
!  real, parameter    :: dlmass=0.5
  real, parameter    :: dlmass=1.0

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: np_node_dim = np/nodes_dim

!  integer(4), parameter :: max_np = density_buffer * ( ((nf_tile-2*nf_buf)*tiles_node_dim/2)**3 + &
!                                  (8*nf_buf**3 + 6*nf_buf*(((nf_tile-2*nf_buf)*tiles_node_dim)**2) + &
!                                  12*(nf_buf**2)*((nf_tile-2*nf_buf)*tiles_node_dim))/8.0 )
!  integer(4), parameter :: np_buffer=int(2./3.*max_np)

  integer(4), parameter :: max_np = np_node_dim*np_node_dim*(np_node_dim/20)
  integer(4), parameter :: np_buffer=int(2./3.*max_np)

  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! parallelization variables
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local

  integer(8) :: plan, iplan

  logical :: firstfftw,do_poisson

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter haloes arrays
  integer(4), parameter :: nvh=7  ! 3*pos+1*mass+3*vel
  real, dimension(7,max_np) :: xvp
  real, dimension(4,np_buffer) :: xp_buf  ! 3*pos+1*mass
  real, dimension(4*np_buffer) :: send_buf, recv_buf
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: den 
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: den_buf 

  !! Power spectrum arrays
  real, dimension(3,nc) :: pkhh
  real, dimension(3,nc) :: poisson
#ifdef PLPLOT
  real*8, dimension(3,nc) :: pkplot
#endif

  !! Fourier transform arrays
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
  real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work

  !! Equivalence arrays to save memory
  equivalence (den,slab_work,recv_cube,xp_buf) 
!  equivalence (slab,xvp,cube)  !! merz --  not sure if xvp is larger than slab?????
!  equivalence (xvp,slab,cube) 
  equivalence (slab,cube)

  !! Common block
#ifdef PLPLOT
!  common xvp,send_buf,slab_work,den_buf,den,cube,slab,xp_buf,recv_buf,pkhh,pkplot
  common xvp,send_buf,den_buf,den,recv_buf,pkhh,pkplot,poisson,do_poisson
#else
!  common xvp,send_buf,slab_work,den_buf,den,cube,slab,xp_buf,recv_buf,pkhh
  common xvp,send_buf,den_buf,den,recv_buf,pkhh,poisson
#endif

!!---start main--------------------------------------------------------------!!

  call mpi_initialize
  if (rank == 0) call writeparams
  firstfftw=.true.  ! initialize fftw so that it generates the plans
  call read_checkpoint_list
  do cur_checkpoint=1,num_checkpoints
    call initvar
    if (rank == 0)print*,'finished initvar'
    call read_haloes
    if (rank == 0)print*,'finished read_halos'
    call pass_haloes
    if (rank == 0)print*,'finished pass_halos'
    do_poisson=.false.
    do cur_massbin=1,num_massbin
      call darkmatterhalo
      if (rank == 0) call writepowerspectra
    enddo
    do_poisson=.true.
    do cur_massbin=1,num_massbin
      call PoissonNoise
      if (rank == 0) call writepowerspectra
    enddo
  enddo
  call cp_fftw(0)
  call mpi_finalize(ierr)

contains

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
      write(*,*) 'cic_power_halo compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'cic_power_halo nodes=',nodes 
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
      write(*,*) 'cic_power_halo running on',nodes,'nodes'
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

! cart_neighbor(1) -> down (negative z)
! cart_neighbor(2) -> up (positive z)
! cart_neighbor(3) -> back (negative y)
! cart_neighbor(4) -> front (positive y)
! cart_neighbor(5) -> left (negative x)
! cart_neighbor(6) -> right (positive x)

    do i = 0, ndim-1
      call mpi_cart_shift(mpi_comm_cart, i, 1, cart_neighbor(2*(i+1)-1), &
                          cart_neighbor(2*(i+1)), ierr)
    enddo

#ifdef DEBUG_LOW
  do i=0,nodes-1
    if (i==rank) write(*,'(8i4)') rank,cart_rank,cart_neighbor
    call mpi_barrier(mpi_comm_world,ierr)
  enddo
#endif

  end subroutine mpi_initialize

!!---------------------------------------------------------------------------!!

  subroutine read_checkpoint_list
!! read in list of checkpoints to calculate spectra for
    implicit none

    integer :: i,fstat

    if (rank == 0) then
!      open(11,file=checkpoints,status='old',iostat=fstat)
      open(11,file=halofinds,status='old',iostat=fstat)
      if (fstat /= 0) then
        print *,'error opening checkpoint list file'
!        print *,'rank',rank,'file:',checkpoints
        print *,'rank',rank,'file:',halofinds
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

  subroutine read_haloes
    implicit none
    
    real z_write,np_total,a
    integer(4) :: j,fstat,np_buf,cube_rank
    integer(4),dimension(3) :: cube_coord
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=255) :: check_name
    real, dimension(28) :: halo_input_buffer
    !real, dimension(17) :: halo_input_buffer

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

!    check_name=output_path//'halo_data/'//z_string(1:len_trim(z_string))//'halo'// &
    check_name=output_path//z_string(1:len_trim(z_string))//'halo'// &
               rank_string(1:len_trim(rank_string))//'.dat'

!! cube coordinates
    cube_coord(3)=rank/nodes_dim/nodes_dim
    cube_rank=rank-cube_coord(3)*nodes_dim*nodes_dim
    cube_coord(2)=cube_rank/nodes_dim
    cube_coord(1)=cube_rank-cube_coord(2)*nodes_dim

!! read number of haloes   
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
    read(21) np_buf
    np_buf=0
    do while(.true.)
      read(21,end=103) halo_input_buffer
      np_buf=np_buf+1
      !! position in local grid
!      xvp(1,np_buf)=halo_input_buffer(1)-cube_coord(1)*nc_node_dim
!      xvp(2,np_buf)=halo_input_buffer(2)-cube_coord(2)*nc_node_dim
!      xvp(3,np_buf)=halo_input_buffer(3)-cube_coord(3)*nc_node_dim
      xvp(1:3,np_buf)=modulo(halo_input_buffer(1:3),real(nc_node_dim))
      !! mass in fine grid units
      xvp(4,np_buf)=halo_input_buffer(17)    !28 coluns
      !xvp(4,np_buf)=halo_input_buffer(15)   !17 columns
      !! velocity
      xvp(5:7,np_buf)=halo_input_buffer(7:9)
    enddo
103 close(21)
    np_local=np_buf

!! tally up total number of haloes
    call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
                         mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'number of haloes =', int(np_total,8)

#ifdef KAISER

    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...
   
    a = 1.0/(1.0 + z_write)
 
    xvp(3,1:np_local)=xvp(3,1:np_local) + xvp(7,1:np_local)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))  

    call pass_haloes

    if(rank==0) then
       write(*,*) '**********************'
       write(*,*) 'Included Kaiser Effect'
       write(*,*) 'Omega_m =', omega_m, 'a =', a
       !write(*,*) '1/H(z) =', 1.5*sqrt(omegam/cubepm_a)
       write(*,*) '1/H(z) =', 1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
       write(*,*) '**********************'
    endif
#endif

  end subroutine read_haloes

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
       
#ifdef DEBUG_LOW
  print *,'rank',rank,'starting swap'
#endif
 
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

#ifdef DEBUG_LOW
  print *,'rank',rank,'starting wait'
#endif
    
    call mpi_waitall(2*nodes_dim**2, requests, wait_status, ierr)

#ifdef DEBUG_LOW
  print *,'rank',rank,'finished wait'
#endif

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
      
#ifdef DEBUG_LOW
  print *,'rank',rank,'finished wait'
#endif

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

    call mpi_barrier(mpi_comm_world,ierr)

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
    include 'fftw_f77.i'

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
            nc,nc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
      call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nc, &
            nc,nc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)
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
    write(*,*) 'np      ', np
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
#ifdef PLPLOT
    integer :: kp
#endif
    real         :: kr
    character*512 :: fn
    character*7  :: z_write
    character*5  :: lmo_write,lmu_write
    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is dm p(k)
    !! 3rd is standard deviation

    write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
    z_write=adjustl(z_write)
    write(lmo_write,'(f5.3)') alog10(min_mass)+(cur_massbin-1.)*dlmass
    lmo_write=adjustl(lmo_write)
    write(lmu_write,'(f5.3)') alog10(min_mass)+cur_massbin*dlmass
    lmu_write=adjustl(lmu_write)


!---- Specify the filename here:

!    fn=output_path//z_write(1:len_trim(z_write))//'cicps_halo_'//lmo_write(1:len_trim(lmo_write))//'_'lmu_write(1:len_trim(lmu_write))//'.dat'
#ifdef KAISER
    if(do_poisson) then
       fn=output_path//z_write(1:len_trim(z_write))//'ngpps_halo_'//lmo_write(1:len_trim(lmo_write))//'_'//lmu_write(1:len_trim(lmu_write))//'-RSD-poisson.dat'
    else
       fn=output_path//z_write(1:len_trim(z_write))//'ngpps_halo_'//lmo_write(1:len_trim(lmo_write))//'_'//lmu_write(1:len_trim(lmu_write))//'-RSD.dat'
    endif
#else
    if(do_poisson) then
       fn=output_path//z_write(1:len_trim(z_write))//'ngpps_halo_'//lmo_write(1:len_trim(lmo_write))//'_'//lmu_write(1:len_trim(lmu_write))//'-poisson.dat'
    else
       fn=output_path//z_write(1:len_trim(z_write))//'ngpps_halo_'//lmo_write(1:len_trim(lmo_write))//'_'//lmu_write(1:len_trim(lmu_write))//'.dat'
    endif
#endif

    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)

    if(do_poisson) then

    do k=2,hc+1
       kr=2*pi*(k-1)/box
#ifdef NGP
       !write(11,*) kr,pkhh(:,k-1),poisson(:,k-1)
       write(11,*) pkhh(3,k-1),poisson(1:2,k-1)
#else
       write(11,*) kr,poisson(:,k)
#endif

#ifdef PLPLOT
       kp=k-1
       pkplot(1,kp)=real(kr,kind=8)
       pkplot(2:3,kp)=real(poisson(:,k),kind=8)
#endif
    enddo

    else

    do k=2,hc+1
       kr=2*pi*(k-1)/box
#ifdef NGP
       !write(11,*) kr,pkhh(:,k-1),poisson(:,k-1)
       write(11,*) pkhh(3,k-1),pkhh(1:2,k-1)
#else
       write(11,*) kr,pkhh(:,k)
#endif

#ifdef PLPLOT
       kp=k-1
       pkplot(1,kp)=real(kr,kind=8)
       pkplot(2:3,kp)=real(pkhh(:,k),kind=8)
#endif
    enddo

    endif




    close(11)

#ifdef PLPLOT
    kp=3
    call plot_power(kp,hc,pkplot(:,:hc),fn(1:len_trim(fn)-4))
#endif

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write power spectra'
    return
  end subroutine writepowerspectra

!!------------------------------------------------------------------!!

  subroutine darkmatterhalo
    implicit none
    integer :: i,j,k
    integer :: i1,j1,k1
    real    :: d,dmin,dmax,sum_halo,sum_halo_local,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart,vfactor
    real, dimension(3) :: dis

    real time1,time2
    call cpu_time(time1)

    !! Initialize density field to be zero
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo

    !! Assign masses to grid to compute halo power spectrum
    call cicmass
    write(*,*) 'nhalos', cur_massbin, nlocal_massbin(cur_massbin), rank

    call mpi_barrier(mpi_comm_world,ierr)

    call mpi_reduce(nlocal_massbin(cur_massbin),ntotal_massbin(cur_massbin),1,mpi_int,mpi_sum,0,mpi_comm_world,ierr) 
    call mpi_bcast(ntotal_massbin(cur_massbin),1,mpi_int,0,mpi_comm_world,ierr)
    if(rank==0) write(*,*) '*** Total number of halos in the bin ***', cur_massbin, ntotal_massbin(cur_massbin)

    !! have to accumulate buffer density 
    call mesh_buffer
    cube=den(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)

    sum_halo_local=sum(cube) 
    call mpi_reduce(sum_halo_local,sum_halo,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) vfactor=dble(sum_halo)/dble(nc)**3
    if (rank == 0) print*,rank,'sum_halo=',sum_halo
    if (rank == 0) print*,rank,'nc=',nc
    if (rank == 0) print*,rank,'vfactor=',vfactor
    call mpi_bcast(vfactor,1,mpi_double_precision,0,mpi_comm_world,ierr)
    call mpi_bcast(sum_halo,1,mpi_real,0,mpi_comm_world,ierr)

    if (sum_halo > 0) then

      !! Convert halo density to delta field
      dmin=0
      dmax=0
      dsum=0
      dvar=0

      do k=1,nc_node_dim
         do j=1,nc_node_dim
            do i=1,nc_node_dim
               cube(i,j,k)=cube(i,j,k)/vfactor-1.0
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
        dsum=dsumt/real(nc)**3
        dvar=sqrt(dvart/real(nc)**3)
        write(*,*)
        write(*,*) 'Min density    ',dmint
        write(*,*) 'Max density    ',dmaxt
        write(*,*) 'Delta sum      ',real(dsum,8)
        write(*,*) 'Delta var      ',real(dvar,8)
        write(*,*)
      endif
 
      !! Forward FFT halo delta field
      call cp_fftw(1)
      if(rank==0) write(*,*) 'Done FFT'

      !! Compute halo power spectrum
      call powerspectrum(slab,pkhh)

      call cpu_time(time2)
      time2=(time2-time1)
      if (rank == 0) write(*,"(f8.2,a)") time2,'  Called cic halo'

   endif   

   return
  end subroutine darkmatterhalo

!------------------------------------------------------------!

  subroutine PoissonNoise
    implicit none
    integer :: i,j,k, fstat
    integer :: i1,j1,k1
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,z_write, sum_halo,sum_halo_local
    real*8  :: dsum,dvar,dsumt,dvart, vfactor
    real, dimension(3) :: dis
    real time1,time2

    call cpu_time(time1)

    write(*,*) 'Calculate Poisson Noise'

    !! Initialized density field to be zero
    !! could do OMP loop here
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo

    ! Share evenly the halos amongst nodes
    np_local = ntotal_massbin(cur_massbin)/nodes

    ! Assign Halos to random positions between 0 and nc_node_dim
    call random_number(xvp(1:3,:))
    xvp(1:3,:) = xvp(1:3,:)*nc_node_dim

    !! Assign masses to grid to compute dm power spectrum
    call cicmass

    !! have to accumulate buffer density 
    call mesh_buffer
    cube=den(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)


    sum_halo_local=sum(cube) 
    call mpi_reduce(sum_halo_local,sum_halo,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) vfactor=dble(sum_halo)/dble(nc)**3
    if (rank == 0) print*,rank,'sum_halo=',sum_halo
    if (rank == 0) print*,rank,'nc=',nc
    if (rank == 0) print*,rank,'vfactor=',vfactor
    call mpi_bcast(vfactor,1,mpi_double_precision,0,mpi_comm_world,ierr)
    call mpi_bcast(sum_halo,1,mpi_real,0,mpi_comm_world,ierr)

    if (sum_halo > 0) then

      sum_dm_local=sum(cube) 
      call mpi_reduce(sum_dm_local,sum_dm,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
      if (rank == 0) print *,'DM total mass=',sum_dm

      !! Convert dm density field to delta field
      dmin=0
      dmax=0
      dsum=0
      dvar=0

      do k=1,nc_node_dim
         do j=1,nc_node_dim
            do i=1,nc_node_dim
               cube(i,j,k)=cube(i,j,k)/vfactor-1.0
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
        dsum=dsumt/real(nc)**3
        dvar=sqrt(dvart/real(nc)**3)
        write(*,*)
        write(*,*) 'Poisson DM min    ',dmint
        write(*,*) 'Poisson DM max    ',dmaxt
        write(*,*) 'Poisson Delta sum ',real(dsum,8)
        write(*,*) 'Poisson Delta var ',real(dvar,8)
        write(*,*)
      endif
 
      !! Forward FFT dm delta field
      call cp_fftw(1)

      !! Compute dm power spectrum
      call powerspectrum(slab,poisson)

    endif

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine PoissonNoise

!!------------------------------------------------------------------!!
  subroutine pass_haloes
    implicit none

    integer i,pp,np_buf,np_exit,np_final,npo,npi
    real x(4),lb,ub
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03

    lb=0.0
    ub=real(nc_node_dim)

    np_buf=0
    pp=1
    do
      if (pp > np_local) exit
      x=xvp(:4,pp)
      if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
          x(3) < lb .or. x(3) >= ub ) then
!        write (*,*) 'PARTICLE OUT',xv(:,pp)
        np_buf=np_buf+1
        if (np_buf > np_buffer) then
          print *,rank,'np_buffer =',np_buffer,'exceeded - np_buf =',np_buf
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif 
        xp_buf(:,np_buf)=xvp(:4,pp)
        xvp(:,pp)=xvp(:,np_local)
        np_local=np_local-1
        cycle 
      endif
      pp=pp+1
    enddo
 
    call mpi_reduce(np_buf,np_exit,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr) 

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'np_exit=',np_buf
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    if (rank == 0) print *,'total exiting haloes =',np_exit

! pass +x

    tag=11 
    npo=0
    pp=1
    do 
      if (pp > np_buf) exit
      if (xp_buf(1,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*4+1:npo*4)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'np_out=',npo
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(6), &
                              tag,cart_neighbor(5),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*4,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*4,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*4+1:pp*4)
      xp_buf(1,np_buf+pp)=max(xp_buf(1,np_buf+pp)-ub,lb)
    enddo

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'x+ np_local=',np_local
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:4,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'x+ np_exit=',np_buf,np_local
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

! pass -x

    tag=12
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(1,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*4+1:npo*4)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle 
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(5), &
                              tag,cart_neighbor(6),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*4,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*4,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*4+1:pp*4)
      xp_buf(1,np_buf+pp)=min(xp_buf(1,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:4,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
  
    np_buf=np_buf+npi

! pass +y

    tag=13 
    npo=0
    pp=1
    do 
      if (pp > np_buf) exit
      if (xp_buf(2,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*4+1:npo*4)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle 
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(4), &
                              tag,cart_neighbor(3),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*4,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*4,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*4+1:pp*4)
      xp_buf(2,np_buf+pp)=max(xp_buf(2,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:4,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

! pass -y

    tag=14
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(2,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*4+1:npo*4)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(3), &
                              tag,cart_neighbor(4),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*4,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*4,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*4+1:pp*4)
      xp_buf(2,np_buf+pp)=min(xp_buf(2,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:4,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
  
    np_buf=np_buf+npi

! pass +z

    tag=15 
    npo=0
    pp=1
    do 
      if (pp > np_buf) exit
      if (xp_buf(4,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*4+1:npo*4)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle 
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(2), &
                              tag,cart_neighbor(1),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*4,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*4,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*4+1:pp*4)
      xp_buf(3,np_buf+pp)=max(xp_buf(3,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:4,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

! pass -z

    tag=16
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(3,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*4+1:npo*4)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(1), &
                              tag,cart_neighbor(2),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*4,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*4,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*4+1:pp*4)
      xp_buf(3,np_buf+pp)=min(xp_buf(3,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:4,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
  
    np_buf=np_buf+npi

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'particles left in buffer=',np_buf
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    call mpi_reduce(np_buf,np_exit,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) print *,'total buffered haloes =',np_exit

    call mpi_reduce(np_local,np_final,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
      print *,'total haloes =',np_final
    endif
 
!!  Check for particles out of bounds

    do i=1,np_local
      if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
          xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
          xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
        print *,'halo out of bounds',rank,i,xvp(:3,i),nc_node_dim
      endif
    enddo

  end subroutine pass_haloes

!------------------------------------------------------------!

  subroutine cicmass
    implicit none

    integer :: i,i1,i2,j1,j2,k1,k2
    real    :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2,vf,v(3),lmass,min_lmass,l_lmass,u_lmass

    min_lmass=alog10(min_mass)
    l_lmass=min_lmass+(cur_massbin-1.)*dlmass
    u_lmass=min_lmass+cur_massbin*dlmass

    if (rank == 0) print*,'halo mass range = ',l_lmass,u_lmass

    print*,rank,'starting cic interpolation for',np_local,' haloes'

    if(np_local>=1)then

    do i=1,np_local

       lmass=alog10(xvp(4,i))
 
!!       print*,'inside loop',lmass,np_local

       if ((lmass.ge.l_lmass.and.lmass.lt.u_lmass) .or. ntotal_massbin(cur_massbin) .gt. 0) then 

!!          print*,'inside loop: if',lmass,l_lmass, u_lmass	
          nlocal_massbin(cur_massbin) = nlocal_massbin(cur_massbin) + 1

          x=xvp(1,i)-0.5
          y=xvp(2,i)-0.5
          z=xvp(3,i)-0.5

          i1=floor(x)+1
          i2=i1+1
          dx1=i1-x
          dx2=1-dx1
          j1=floor(y)+1
          j2=j1+1
          dy1=j1-y
          dy2=1-dy1
          k1=floor(z)+1
          k2=k1+1
          dz1=k1-z
          dz2=1-dz1

          if (i1 < 0 .or. i2 > nc_node_dim+1 .or. j1 < 0 .or. &
              j2 > nc_node_dim+1 .or. k1 < 0 .or. k2 > nc_node_dim+1) then 
             print *,'halo out of bounds',i1,i2,j1,j2,k1,k2,nc_node_dim
          endif 

          den(i1,j1,k1)=den(i1,j1,k1)+dx1*dy1*dz1
          den(i2,j1,k1)=den(i2,j1,k1)+dx2*dy1*dz1
          den(i1,j2,k1)=den(i1,j2,k1)+dx1*dy2*dz1
          den(i2,j2,k1)=den(i2,j2,k1)+dx2*dy2*dz1
          den(i1,j1,k2)=den(i1,j1,k2)+dx1*dy1*dz2
          den(i2,j1,k2)=den(i2,j1,k2)+dx2*dy1*dz2
          den(i1,j2,k2)=den(i1,j2,k2)+dx1*dy2*dz2
          den(i2,j2,k2)=den(i2,j2,k2)+dx2*dy2*dz2

       endif

    enddo

    end if

    print*,rank,'done cic interpolation'

    return
  end subroutine cicmass

!!--------------------------------------------------------------!!

  subroutine powerspectrum(delta,pk)
    implicit none
    real, dimension(3,nc)       :: pk
    real, dimension(nc+2,nc,nc_slab) :: delta

    real, parameter :: pi4=3.14159
    integer :: i,j,k,kg
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow, x,y,z,sync_x, sync_y, sync_z, kernel
    real, dimension(3,nc,nc_slab) :: pkt
    real, dimension(3,nc) :: pktsum
    real, dimension(nc) :: kcen, kcount
    real    :: kavg

    real time1,time2
    call cpu_time(time1)

    pkt=0.0
    pktsum=0.0

    kcen(:)   = 0.
    kcount(:) = 0.


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
             if(kx.eq.0 .and. ky <=0 .and. kz <=0)cycle;
             if(kx.eq.0 .and. ky >0 .and. kz <0)cycle;
             if (kr .ne. 0) then
                k1=ceiling(kr)
                k2=k1+1

                x = pi*real(kx)/ncr
                y = pi*real(ky)/ncr
                z = pi*real(kz)/ncr

                if(x==0) then
                   sync_x = 1
                else
                   sync_x = sin(x)/x
                endif
                if(y==0) then
                   sync_y = 1
                else
                   sync_y = sin(y)/y
                endif
                if(z==0) then
                   sync_z = 1
                else
                   sync_z = sin(z)/z
                endif

                kernel = sync_x*sync_y*sync_z


#ifdef NGP
                w1=1
                w2=0
#else
                w1=k1-kr
                w2=1-w1
#endif
                pow=sum((delta(i:i+1,j,k)/real(ncr)**3)**2)/kernel**4
                pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                pkt(3,k1,k)=pkt(3,k1,k)+w1
                pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                pkt(3,k2,k)=pkt(3,k2,k)+w2

                kcen(k1) = kcen(k1) + w1 * kr
                kcen(k2) = kcen(k2) + w2 * kr

                kcount(k1) = kcount(k1) + w1
                kcount(k2) = kcount(k2) + w2


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
        ! TESTING
        !write(18,*) k,nc,pktsum(1,k),pktsum(2,k),pktsum(3,k)
        if (pktsum(3,k) .eq. 0) then
          pk(1:2,k)=0.
        else
          pk(1:2,k)=pktsum(1:2,k)/pktsum(3,k)
          pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))

          kavg = kcen(k) / kcount(k)
          pk(3,k) = 2. * pi * kavg / box

          !write(19,*) pk(3,k),nc,pk(1,k),pk(2,k)
! TESTING
!          pk(1:2,k)=4.*pi*real(k-1)**3*pk(1:2,k)
#ifdef NGP
          pk(1:2,k)=4.*pi4*kavg**3*pk(1:2,k)
#else
	  pk(1:2,k)=4.*pi4*real(k-1)**3*pk(1:2,k)
#endif
          !write(20,*) pk(3,k),nc,pk(1,k),pk(2,k)
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

    do k=1,max_np
       xvp(:,k)=0
    enddo
    do k=1,nc_slab
       slab_work(:,:,k)=0
    enddo
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo
    do k=1,nc_node_dim
       cube(:,:,k)=0
    enddo
    do k=1,nc_slab
       slab(:,:,k)=0
    enddo
    do k=1,np_buffer
       xp_buf(:,k)=0
    enddo
    do k=1,3*np_buffer
       recv_buf(k)=0
    enddo
    do k=1,nc
       pkhh(:,k)=0
    enddo

    nlocal_massbin(:)=0   
    ntotal_massbin(:)=0   
 
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

!!------------------------------------------------------------------!!

subroutine mesh_buffer
!! mesh_buffer -- buffer cubic decomposition mesh
  implicit none

  integer(4) :: buffer_size
  integer(4) :: tag
  integer(4) :: status(MPI_STATUS_SIZE)

    buffer_size = (nc_node_dim + 2)**2

  tag=64

!! send to node in -x

    den_buf(:,:)=den(0,:,:)
    call mpi_sendrecv_replace(den_buf,buffer_size,mpi_real, &
                              cart_neighbor(5),tag,cart_neighbor(6), &
                              tag,mpi_comm_cart,status,ierr)
    den(nc_node_dim,:,:)=den(nc_node_dim,:,:)+den_buf(:,:)

!! send to node in +x
   
      den_buf(:,:)=den(nc_node_dim+1,:,:)
      call mpi_sendrecv_replace(den_buf,buffer_size,mpi_real, &
                              cart_neighbor(6),tag,cart_neighbor(5), &
                              tag,mpi_comm_cart,status,ierr)
      den(1,:,:)=den(1,:,:)+den_buf(:,:)

!! send to node in -y

      den_buf(:,:)=den(:,0,:)
      call mpi_sendrecv_replace(den_buf,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)
      den(:,nc_node_dim,:)=den(:,nc_node_dim,:)+den_buf(:,:)

!! send to node in +y

      den_buf(:,:)=den(:,nc_node_dim+1,:)
      call mpi_sendrecv_replace(den_buf,buffer_size,mpi_real, &
                              cart_neighbor(4),tag,cart_neighbor(3), &
                              tag,mpi_comm_cart,status,ierr)
      den(:,1,:)=den(:,1,:)+den_buf(:,:)

!! send to node in -z
    
      den_buf(:,:)=den(:,:,0)
      call mpi_sendrecv_replace(den_buf,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)
      den(:,:,nc_node_dim)=den(:,:,nc_node_dim)+den_buf(:,:)

!! send to node in +z

      den_buf(:,:)=den(:,:,nc_node_dim+1)
      call mpi_sendrecv_replace(den_buf,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &
                              tag,mpi_comm_cart,status,ierr)

      den(:,:,1)=den(:,:,1)+den_buf(:,:)

  end subroutine mesh_buffer

end program cic_power_halo 
