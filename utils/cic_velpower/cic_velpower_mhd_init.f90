!! cic_velpower_kernel_corr.f90 Parallelized: JD Emberson Jan 19, 2012
!! Compile with: mpif77 -fpp -g -w -O3 -axN cic_velpower_kernel_corr.f90 -o cic_velpower  -L/home/merz/lib/fftw-2.1.5_intel8/lib -I/home/merz/lib/fftw-2.1.5_intel8/include -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

program cic_velpower_mhd 

  implicit none

  include 'mpif.h'

  !! frequently changed parameters are found in this header file:
  include '../../parameters'

  logical, parameter :: correct_kernel=.false.

  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints'

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc / 2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc
  real, parameter    :: npr=np

  !! internals
  integer, parameter :: max_checkpoints = 100
  real, dimension(max_checkpoints) :: z_checkpoint
  integer num_checkpoints, cur_checkpoint

  !! have velocity power spectra for each x, y, z
  integer cur_dimension

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc / nodes_dim
  integer(4), parameter :: np_node_dim = np / nodes_dim
  integer(4), parameter :: np_buffer = 4 * np_node_dim**3
  integer(4), parameter :: max_np = np_node_dim**3 + np_buffer
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! parallelization variables
  integer(4), dimension(0:nodes_dim-1, 0:nodes_dim-1) :: slab_neighbor
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local

  integer(8) :: plan, iplan

  logical :: firstfftw

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter arrays
  real, dimension(6, max_np) :: xvp
  real, dimension(6, np_buffer) :: xp_buf
  real, dimension(6 * np_buffer) :: send_buf, recv_buf

  !! Power spectrum arrays
  real, dimension(3, 2, nc) :: pkcurldm
  real, dimension(3, 2, nc) :: pkmomdim
  real, dimension(2, nc) :: pkdivdm
  real, dimension(2, nc) :: pkmassdm
  real, dimension(2, nc) :: pkdm
#ifdef PLPLOT
  real*8, dimension(3, nc) :: pkplot
#endif

  !! Internal energy array
  real, dimension(5, nc_node_dim, nc_node_dim, nc_node_dim) :: u

  !! Fourier transform arrays
  real, dimension(nc_node_dim, nc_node_dim, nc_node_dim) :: cube
  real, dimension(nc_node_dim, nc_node_dim, nc_slab, 0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2, nc, nc_slab) :: slab, slab_work

  !! Density arrays
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: den
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: den_buf

  !! Array containing (x, y, z) components of the momentum density field
  real, dimension(3, 0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momden
  real, dimension(nc_node_dim**2) :: momden_send_buff
  real, dimension(nc_node_dim**2) :: momden_recv_buff

  !! Array containing (x, y, z) components of the curl of the momentum density field
  real, dimension(3, 0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momcurl

  !! Array containing divergence of the momentum density field
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momdiv

  !! Array containing the matter density field
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: massden

  !! Equivalence arrays to save memory
  equivalence (slab_work,recv_cube) 

  !! Common block
#ifdef PLPLOT
  common xvp, send_buf, slab_work, den_buf, den, cube, slab, xp_buf, recv_buf, pkcurldm, pkdivdm, pkdm, pkplot
#else
  common xvp, send_buf, slab_work, den_buf, den, cube, slab, xp_buf, recv_buf, pkcurldm, pkdivdm, pkdm
#endif

  type comm_wld
    integer :: g                        ! Global number of zones in one direction
    integer :: r                        ! Global index of index 0 in local array
    integer :: m, n                     ! Start and end of local section without buffers
    integer :: l                        ! Dimension of array inclusive buffers
    integer, dimension(4) :: requests   ! Communication handles
  end type comm_wld

  type(comm_wld) :: nx,ny,nz

! -------------------------------------------------------------------------------------------------------
! MAIN
! -------------------------------------------------------------------------------------------------------

  call mpi_initialize

  if (rank == 0) call writeparams

  firstfftw = .true.  ! initialize fftw so that it generates the plans

  call initvar
  call momentum_density
  call pass_momdensity 

  !
  ! Compute curl power spectrum
  !

  call momentum_curl

  do cur_dimension = 1, 3 !! Each curl component 

    call gasmatter(0)

  enddo

  if (rank == 0) call writepowerspectra(0)

  !
  ! Compute divergence power spectrum
  !

  call momentum_divergence
  call gasmatter(1)

  if (rank == 0) call writepowerspectra(1)

  !
  ! Compute cross-power spectrum between divergence and density field
  !

  do cur_dimension = 1, 3

      call gasmatter(2)

  enddo

  if (rank == 0) call writepowerspectra(2)

  call cp_fftw(0)
  call mpi_finalize(ierr)

! -------------------------------------------------------------------------------------------------------
! SUBROUTINES
! -------------------------------------------------------------------------------------------------------

contains

! -------------------------------------------------------------------------------------------------------

subroutine mpi_initialize

    implicit none
     
    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder
  
    !! Set up global mpi communicator

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)

    call mpi_comm_size(mpi_comm_world, nodes_returned, ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)
    
    if (nodes_returned /= nodes ) then
      write(*,*) 'cic_pow compiled for a different number of nodes'
      write(*,*) 'mpirun nodes = ', nodes_returned, ' cic_pow nodes = ',nodes 
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif
    
    if (mod(nc, nodes) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into slabs'
      write(*,*) 'nc = ', nc, ' nodes = ', nodes, ' mod(nc, nodes) != 0'
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)

    if (rank == 0) then
      write(*,*) 'cic_pow running on ', nodes, ' nodes'
      write(*,*) 'using cubic distribution: ', nodes_dim, ' nodes per dimension'
      write(*,*) nc, ' cells in mesh'
    endif

    !! Calculate coordinates within slab for cube processes

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

    !! Create cartesian communicator based on cubic decomposition

    dims(:) = nodes_dim
    periodic(:) = .true.
    reorder = .false.
    ndim = 3

    call mpi_cart_create(mpi_comm_world, ndim, dims, periodic, &
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
    do i = 0, nodes-1
        if (i == rank) write(*, '(8i4)') rank, cart_rank, cart_neighbor
        call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif

end subroutine mpi_initialize

! -------------------------------------------------------------------------------------------------------

subroutine read_checkpoint_list
    !
    ! Read in the list of redshift checkpoints for which to calculate spectra for 
    !

    implicit none

    integer :: i, fstat

    if (rank == 0) then
      
        open(11, file=checkpoints, status='old', iostat=fstat)
   
        !! Check for opening error 
        if (fstat /= 0) then
            print *,'ERROR: Cannot open checkpoint list file'
            print *,'rank ', rank, ' file: ', checkpoints
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
      
        !! Read the redshifts
        do num_checkpoints = 1, max_checkpoints 
            read(unit=11, err=51, end=41, fmt='(f20.10)') z_checkpoint(num_checkpoints)
        enddo

        !! Tabulate total number of checkpoints
   41  num_checkpoints = num_checkpoints - 1
   51  close(11)
      
        !! Print to screen
        print *, 'Checkpoints to recompose:'
        do i = 1, num_checkpoints
            write(*, '(f5.1)') z_checkpoint(i)
        enddo
    
    endif

    call mpi_bcast(num_checkpoints, 1, mpi_integer, 0, mpi_comm_world, ierr)

end subroutine read_checkpoint_list

! -------------------------------------------------------------------------------------------------------

subroutine initvar
    !
    ! Initialize data arrays
    !

    implicit none

    integer :: k

    real time1,time2
    call cpu_time(time1)

    !! Internal energy array
    do k = 1, nc_node_dim
        u(:, :, :, k) = 0.
    enddo

    !! Momentum density and curl arrays
    do k = 0, nc_node_dim + 1
        momden(:, :, :, k) = 0.
        momcurl(:, :, :, k) = 0.
        momdiv(:, :, k) = 0.
        massden(:, :, k) = 0.
    enddo

    !! Fourier transform arrays
    do k = 1, nc_slab
       slab_work(:, :, k)=0
    enddo

    do k = 1, nc_node_dim
       cube(:, :, k) = 0
    enddo

    do k = 1, nc_slab
       slab(:, :, k) = 0
    enddo

    do k = 1, np_buffer
       xp_buf(:, k) = 0
    enddo

    do k = 1, 3 * np_buffer
       recv_buf(k) = 0
    enddo

    !! Power spectrum arrays
    do k = 1, nc
        pkcurldm(:, :, k) = 0.
        pkdivdm(:, k) = 0.
        pkdm(:, k) = 0.
    enddo

    call cpu_time(time2)
    time2=(time2-time1)

    if (rank == 0) write(*, "(f8.2,a)") time2, '  Called init var'

    return

end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine momentum_density 
    !
    ! Read internal energy of mhd particles and assign to momentum density grid
    !

    implicit none
    
    real z_write, np_total
    integer l, k, j, i, fstat
    character(len=4) :: rank_string
    character(len=100) :: check_name

    write(rank_string,'(i4)') rank
    rank_string=adjustl(rank_string)

    check_name=output_path//'mhd_ic'// &
               rank_string(1:len_trim(rank_string))//'.dat'

    !
    ! Open and read the file    
    !

#ifdef BINARY
    open(unit=21, file=check_name, status='old', iostat=fstat, form='binary')
#else
    open(unit=21, file=check_name, status='old', iostat=fstat, form='unformatted')
#endif

    !! Check for opening error
    if (fstat /= 0) then
      write(*,*) 'ERROR: Cannot open checkpoint position file'
      write(*,*) 'rank', rank, ' file: ',check_name
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    !! Read in the checkpoint body
    read(21) u
    
    close(21)
    
    !
    ! Assign to the momentum density field
    !

    do j = 1, 3

        momden(j, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) = -u(1+j, :, :, :)

    enddo 

end subroutine momentum_density

! -------------------------------------------------------------------------------------------------------

subroutine mass_density
    !
    ! Read internal energy of mhd particles and assign to momentum density grid
    !

    implicit none

    real z_write, np_total
    integer l, k, j, i, fstat
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name

    !! These are unnecessary headers from the checkpoint
    integer(4) :: cur_iter, cur_t

    !! Generate checkpoint names on each node
    if (rank==0) then
      z_write = z_checkpoint(cur_checkpoint)
      print *,'Calculating spectrum for z = ',z_write
    endif

    call mpi_bcast(z_write, 1, mpi_real, 0, mpi_comm_world, ierr)

    !! Determine the file name
    write(z_string,'(f7.3)') z_write
    z_string=adjustl(z_string)

    write(rank_string,'(i4)') rank
    rank_string=adjustl(rank_string)

    check_name = output_path//z_string(1:len_trim(z_string))//'mhd'// &
               rank_string(1:len_trim(rank_string))//'.dat'

    !
    ! Open and read the file    
    !

#ifdef BINARY
    open(unit=21, file=check_name, status='old', iostat=fstat, form='binary')
#else
    open(unit=21, file=check_name, status='old', iostat=fstat, form='unformatted')
#endif

    !! Check for opening error
    if (fstat /= 0) then
      write(*,*) 'ERROR: Cannot open checkpoint position file'
      write(*,*) 'rank', rank, ' file: ',check_name
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    !! Read in checkpoint header data
    read(21) cur_iter, cur_t, nx, ny, nz

    !! Read in the checkpoint body
    read(21) u

    close(21)

    !
    ! Assign to the mass density field
    !

    massden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) = u(1, :, :, :)

end subroutine mass_density

! -------------------------------------------------------------------------------------------------------

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
    
! -------------------------------------------------------------------------------------------------------

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

! -------------------------------------------------------------------------------------------------------

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

! -------------------------------------------------------------------------------------------------------

  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nodes   ', nodes
    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams

! -------------------------------------------------------------------------------------------------------

subroutine writepowerspectra(command)
    !
    ! Writes the dimensionless power spectrum for the curl/divergence of the momentum density field
    !    

    implicit none

    integer      :: i, j, k
#ifdef PLPLOT
    integer :: kp
#endif
    real         :: kr
    character*180 :: fn
    character*7  :: prefix
    real    :: vsim2phys
    integer(4) :: command ! 0 for curl, 1 for divergence, 2 for divergence-matter cross 

    real time1,time2
    call cpu_time(time1)

    vsim2phys = 300. * sqrt(omega_m) * box * (1. + z_i) / 2. / nc

    !
    ! Determine name of output file
    !

#ifdef NGP 
    if (command == 0) then
        prefix = 'ngpcvps'
    else if (command == 1) then
        prefix = 'ngpdvps'
    else if (command == 2) then
        prefix = 'ngptvps'
    endif
#else
    if (command == 0) then
        prefix = 'ciccvps'
    else if (command == 1) then
        prefix = 'cicdvps'
    else if (command == 2) then
        prefix = 'cictvps'
    endif
#endif

#ifdef KAISER
   fn=output_path//prefix//'-RSD_mhd_init.dat' 
#else
   fn=output_path//prefix//'_mhd_init.dat' 
#endif

    !
    ! Asign data to be written
    !

    do i = 1, nc
        pkdm(:, i) = 0.
    enddo

    if (command == 0) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkcurldm(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkcurldm(j, 2, i)
            enddo
        enddo

    else if (command == 1) then

        do i = 1, nc
            pkdm(1, i) = pkdivdm(1, i)
            pkdm(2, i) = pkdivdm(2, i)
        enddo

    else if (command == 2) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkmomdim(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkmomdim(j, 2, i)
            enddo
        enddo

    endif

    !
    ! Convert to physical units in km/s
    !

    do i = 1, nc

        pkdm(1, i) = vsim2phys**2 * pkdm(1, i)
        pkdm(2, i) = vsim2phys**2 * pkdm(2, i)

    enddo

    !
    ! Write to output file with column ordering [k, p(k), sigma(k)]
    !

    write(*,*) 'Writing ', fn
    open(11, file=fn, recl=500)

    do k = 2, hc + 1

       kr = 2 * pi * (k-1) / box
#ifdef NGP
       write(11,*) kr, pkdm(:,k-1)
#else
       write(11,*) kr, pkdm(:,k)
#endif
#ifdef PLPLOT
       kp = k-1
       pkplot(1, kp)=real(kr, kind=8)
       pkplot(2:3, kp)=real(pkdm(:, k), kind=8)
#endif
    enddo
    close(11)

#ifdef PLPLOT
    kp = 3
    call plot_power(kp, hc, pkplot(:, :hc), fn(1:len_trim(fn)-4))
#endif

    call cpu_time(time2)
    time2=time2-time1

    write(*,"(f8.2,a)") time2,'  Called writepowerspectra'

    return

end subroutine writepowerspectra

! -------------------------------------------------------------------------------------------------------

subroutine gasmatter(command)

    implicit none

    integer :: i, j, k
    integer :: i1, j1, k1
    real    :: d, dmin, dmax, sum_dm, sum_dm_local, dmint, dmaxt
    real*8  :: dsum, dvar, dsumt, dvart
    real, dimension(3) :: dis

    integer(4) :: command ! 0 for curl, 1 for divergence, 2 for matter 

    real time1,time2
    call cpu_time(time1)

    !
    ! Initialize density array to zero 
    !

    do k = 0, nc_node_dim + 1
        den(:, :, k) = 0.
    enddo

    !
    ! Assign data to density grid
    !

    if (command == 0) then !! Fill with given curl component 

        den(:, :, :) = momcurl(cur_dimension, :, :, :)

    else if (command == 1) then !! Fill with divergence 

        den(:, :, :) = momdiv(:, :, :)

    else if (command == 2) then !! Fill with matter density field

        den(:, :, :) = momden(cur_dimension, :, :, :) 

    endif

    !
    ! Perpare for Fourier transform
    !

    !! Have to accumulate buffer density 
    call mesh_buffer
    cube = den(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)

    sum_dm_local = sum(cube)
    call mpi_reduce(sum_dm_local, sum_dm, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
    if (rank == 0) print *,'Gas total momentum density =', sum_dm

    !! Convert gas density field to delta field
    dmin = 0
    dmax = 0
    dsum = 0
    dvar = 0

    do k = 1, nc_node_dim
       do j = 1, nc_node_dim
          do i = 1, nc_node_dim
             !cube(i, j, k) = cube(i, j, k) - 1.
             d = cube(i, j, k)
             dsum = dsum + d
             dvar = dvar + d*d
             dmin = min(dmin, d)
             dmax = max(dmax, d)
          enddo
       enddo
    enddo

    call mpi_reduce(dsum, dsumt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(dvar, dvart, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(dmin, dmint, 1, mpi_real, mpi_min, 0, mpi_comm_world, ierr)
    call mpi_reduce(dmax, dmaxt, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)

    if (rank == 0) then

      dsum = dsumt / nc**3
      dvar = sqrt(dvart / nc**3)
      write(*,*)
      write(*,*) 'Gas min   ', dmint
      write(*,*) 'Gas max   ', dmaxt
      write(*,*) 'Delta sum ', real(dsum)
      write(*,*) 'Delta var ', real(dvar)
      write(*,*)

    endif

    ! 
    ! Forward FFT gas delta field
    !    

    call cp_fftw(1)

    !
    ! Compute power spectrum
    !

    if (command == 0) then

        call powerspectrum(slab, pkcurldm(cur_dimension, :, :))

    else if (command == 1) then

        call powerspectrum(slab, pkdivdm)

    else if (command == 2) then

        call powerspectrum(slab, pkmomdim(cur_dimension, :, :)) 

    endif

    call cpu_time(time2)
    time2=(time2-time1)

    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called gasmatter'

    return

end subroutine gasmatter

! -------------------------------------------------------------------------------------------------------

subroutine pass_momdensity
    !
    ! Fill momdensity buffer with data from adjacent nodes.
    !

    implicit none

    integer :: k, j, i, ind
    integer :: num2send, npi
    integer, dimension(mpi_status_size) :: status, sstatus, rstatus
    integer :: tag, srequest, rrequest, sierr, rierr, ierr

    !
    ! Pass each x, y, z component of the momentum density separately
    !

    num2send = nc_node_dim**2

    do i = 1, 3

        !
        ! Pass +x
        ! 

        tag = 111

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden_send_buff(ind) = momden(i, nc_node_dim, j, k)
                ind = ind + 1
            enddo
        enddo

        npi = num2send

        call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(6), &
                                  tag, cart_neighbor(5), tag, mpi_comm_world, &
                                  status, ierr)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden(i, 0, j, k) = momden_recv_buff(ind)
                ind = ind + 1
            enddo
        enddo

        !
        ! Pass -x
        ! 

        tag = 112

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden_send_buff(ind) = momden(i, 1, j, k)
                ind = ind + 1
            enddo
        enddo

        npi = num2send

        call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(5), &
                                  tag, cart_neighbor(6), tag, mpi_comm_world, &
                                  status, ierr)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden(i, nc_node_dim+1, j, k) = momden_recv_buff(ind)
                ind = ind + 1
            enddo
        enddo

        !
        ! Pass +y
        ! 

        tag = 113

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden_send_buff(ind) = momden(i, j, nc_node_dim, k)
                ind = ind + 1
            enddo
        enddo

        npi = num2send

        call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(4), &
                                  tag, cart_neighbor(3), tag, mpi_comm_world, &
                                  status, ierr)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden(i, j, 0, k) = momden_recv_buff(ind)
                ind = ind + 1
            enddo
        enddo

        !
        ! Pass -y
        ! 

        tag = 114

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden_send_buff(ind) = momden(i, j, 1, k)
                ind = ind + 1
            enddo
        enddo

        npi = num2send

        call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(3), &
                                  tag, cart_neighbor(4), tag, mpi_comm_world, &
                                  status, ierr)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden(i, j, nc_node_dim+1, k) = momden_recv_buff(ind)
                ind = ind + 1
            enddo
        enddo

        !
        ! Pass +z
        ! 

        tag = 115

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden_send_buff(ind) = momden(i, j, k, nc_node_dim)
                ind = ind + 1
            enddo
        enddo

        npi = num2send

        call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(2), &
                                  tag, cart_neighbor(1), tag, mpi_comm_world, &
                                  status, ierr)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden(i, j, k, 0) = momden_recv_buff(ind)
                ind = ind + 1
            enddo
        enddo

        !
        ! Pass -z
        ! 

        tag = 116

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden_send_buff(ind) = momden(i, j, k, 1)
                ind = ind + 1
            enddo
        enddo

        npi = num2send

        call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(1), &
                                  tag, cart_neighbor(2), tag, mpi_comm_world, &
                                  status, ierr)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        ind = 1
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                momden(i, j, k, nc_node_dim+1) = momden_recv_buff(ind)
                ind = ind + 1
            enddo
        enddo

    enddo

    return

end subroutine pass_momdensity

! -------------------------------------------------------------------------------------------------------

subroutine momentum_curl
    !
    ! Calculates the curl of the momentum density field
    !

    implicit none
    integer :: k, j, i
    real :: dx

    !! Spacing between grid points
    dx = 2. 

    !! Compute derivatives using dy/dx(a) = y(a+1) - y(a-1) / x(a+1) - x(a-1)

    do i = 1, nc_node_dim
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim

                !! x component: dp_z/dy - dp_y/dz
                momcurl(1, i, j, k) = (momden(3, i, j+1, k) - momden(3, i, j-1, k) - &
                                          momden(2, i, j, k+1) + momden(2, i, j, k-1)) / dx

                !! y component: dp_x/dz - dp_z/dx
                momcurl(2, i, j, k) = (momden(1, i, j, k+1) - momden(1, i, j, k-1) - &
                                          momden(3, i+1, j, k) + momden(3, i-1, j, k)) / dx

                !! z component: dp_y/dx - dp_x/dy
                momcurl(3, i, j, k) = (momden(2, i+1, j, k) - momden(2, i-1, j, k) - &
                                          momden(1, i, j+1, k) + momden(1, i, j-1, k)) / dx

            enddo
        enddo
    enddo

end subroutine momentum_curl

! -------------------------------------------------------------------------------------------------------

subroutine momentum_divergence
    !
    ! Calculates the divergence of the momentum density field
    !

    implicit none
    integer :: k, j, i
    real :: dx

    !! Spacing between grid points
    dx = 2. 

    !! Compute derivatives using dy/dx(a) = y(a+1) - y(a-1) / x(a+1) - x(a-1)

    do i = 1, nc_node_dim
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim

                momdiv(i, j, k) = (momden(1, i+1, j, k) - momden(1, i-1, j, k) + &
                                    momden(2, i, j+1, k) - momden(2, i, j-1, k) + &
                                    momden(3, i, j, k+1) - momden(3, i, j, k-1)) / dx

            enddo
        enddo
    enddo

end subroutine momentum_divergence

! -------------------------------------------------------------------------------------------------------

  subroutine powerspectrum(delta,pk)
    implicit none
    real, dimension(2,nc)       :: pk
    real, dimension(nc+2,nc,nc_slab) :: delta

    integer :: i,j,k,kg
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow, x,y,z,sync_x, sync_y, sync_z,kernel
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
                w1=k1-kr
                w2=1-w1
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
                !write(*,*) i,j,k,kx,ky,kz,x,y,z,kernel
                !pause
#ifdef NGP
                w1=1
                w2=0
#endif                
                pow=sum((delta(i:i+1,j,k)/ncr**3)**2)/kernel
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
            if (pktsum(3,k) .eq. 0) then
                pk(:,k)=0
            else
                pk(:,k)=pktsum(1:2,k)/pktsum(3,k)
                pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))

                kavg = kcen(k) / kcount(k)

#ifdef NGP
                pk(1:2,k)=4*pi*(kavg)**3*pk(1:2,k)
#else
                pk(1:2,k)=4*pi*(kavg-1.)**3*pk(1:2,k)
#endif
            endif
        enddo
    endif

    call mpi_bcast(pk,2*nc,mpi_real,0,mpi_comm_world,ierr)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called power spectrum'
    return
  end subroutine powerspectrum

! -------------------------------------------------------------------------------------------------------

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

end program cic_velpower_mhd
 
