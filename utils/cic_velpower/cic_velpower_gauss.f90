! ------------------------------------------------------------------------------------------------------- 
! cic_velpower.f90 ... Last edit: February 11, 2013 by JD Emberson
! -------------------------------------------------------------------------------------------------------

program cic_velpower 

  implicit none

  include 'mpif.h'
  include '../../parameters'

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
  real, parameter :: pi = 3.14159
  real, parameter :: cell_smooth = 1.0 
  real, parameter :: R_smooth = box / nc * cell_smooth 

  !! Dark matter arrays
  real, dimension(6, max_np) :: xvp
  real, dimension(6, np_buffer) :: xp_buf
  real, dimension(6 * np_buffer) :: send_buf, recv_buf

  !! Power spectrum arrays
  real, dimension(3, 3, nc) :: pkcurldm
  real, dimension(3, 3, nc) :: pkmomdim
  real, dimension(3, nc) :: pkdivdm
  real, dimension(3, nc) :: pkdm

  !! Fourier transform arrays
  real, dimension(nc_node_dim, nc_node_dim, nc_node_dim) :: cube
  real, dimension(nc_node_dim, nc_node_dim, nc_slab, 0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2, nc, nc_slab) :: slab, slab_work

  !! Array containing (x, y, z) components of the momentum density field
  real, dimension(3, 0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momden
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: momden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: momden_recv_buff

  !! Array containing (x, y, z) components of the curl of the momentum density field
  real, dimension(3, nc_node_dim, nc_node_dim, nc_node_dim) :: momcurl

  !! Array containing divergence of the momentum density field
  real, dimension(nc_node_dim, nc_node_dim, nc_node_dim) :: momdiv

  !! Array containing the matter density field
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: massden
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: massden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: massden_recv_buff

  !! Equivalence arrays to save memory
  equivalence (slab_work,recv_cube) 

  !! Common block
  common xvp, send_buf, slab_work, cube, slab, xp_buf, recv_buf, pkcurldm, pkmomdim, pkdivdm, pkdm

! -------------------------------------------------------------------------------------------------------
! MAIN
! -------------------------------------------------------------------------------------------------------

  call mpi_initialize

  if (rank == 0) call writeparams

  firstfftw = .true.  ! initialize fftw so that it generates the plans

  call read_checkpoint_list

  do cur_checkpoint = 1, num_checkpoints

    call initvar
    call read_particles
    call pass_particles
    call momentum_density
    call buffer_momdensity
    call pass_momdensity

    !
    ! Compute curl power spectrum
    !

    call momentum_curl

    do cur_dimension = 1, 3 !! Each curl component 

      call darkmatter(0)
    
    enddo

    if (rank == 0) call writepowerspectra(0)

    !
    ! Compute divergence power spectrum
    !

    call momentum_divergence
    call darkmatter(1)

    if (rank == 0) call writepowerspectra(1)

    !
    ! Compute power spectrum of total velocity field
    !

    !! Convert momentum field into velocity field 
    call mass_density
    call buffer_massdensity
    call momentum2velocity

    do cur_dimension = 1, 3

        call darkmatter(2)

    enddo

    if (rank == 0) call writepowerspectra(2)

#ifdef write_vel
    call writevelocityfield
#endif

  enddo

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

    !! Particle positions and velocities
    do k = 1, max_np
       xvp(:, k) = 0.
    enddo

    !! Momentum and matter density arrays
    do k = 0, nc_node_dim + 1
        momden(:, :, :, k) = 0.
        massden(:, :, k) = 0.
    enddo

    !! Divergence and curl arrays
    do k = 1, nc_node_dim
        momcurl(:, :, :, k) = 0.
        momdiv(:, :, k) = 0.
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
        pkmomdim(:, :, k) = 0.
        pkdivdm(:, k) = 0.
        pkdm(:, k) = 0.
    enddo

    return

end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine read_particles
    !
    ! Read x, y, z positions and velocities and store in xvp
    !

    implicit none
    
    real z_write, np_total
    integer j, fstat
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name

    !! These are unnecessary headers from the checkpoint
    real(4) :: a, t, tau, dt_f_acc, dt_c_acc, dt_pp_acc, mass_p
    integer(4) :: nts, sim_checkpoint, sim_projection, sim_halofind

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

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

    check_name = output_path//z_string(1:len_trim(z_string))//'xv'// &
               rank_string(1:len_trim(rank_string))//'.dat'

    !! Open the file    
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
#ifdef PPINT
    read(21) np_local, a, t, tau, nts, dt_f_acc, dt_pp_acc, dt_c_acc, sim_checkpoint, &
               sim_projection, sim_halofind, mass_p
#else
    read(21) np_local, a, t, tau, nts, dt_f_acc, dt_c_acc, sim_checkpoint, &
               sim_projection, sim_halofind, mass_p
#endif

    !! Check for memory problems
    if (np_local > max_np) then
      write(*,*) 'ERROR: Too many particles to store in memory!'
      write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    !! Tally up total number of particles
    call mpi_reduce(real(np_local, kind=4), np_total, 1, mpi_real, &
                         mpi_sum, 0, mpi_comm_world, ierr)
    
    if (rank == 0) write(*,*) 'Total number of particles = ', int(np_total,4)

    !! Read positions and velocities
    read(21) xvp(:, :np_local)

    close(21)
 
#ifdef KAISER

    !! Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !! Converting seconds into simulation time units
    !! cancels the H0...
    
    !! xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) + xv(6,ip+1:ip+nploc(i))*1.5*sqrt(omegam/cubepm_a)
    !! xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) + xv(6,ip+1:ip+nploc(i))*1.5/sqrt(cubepm_a*(1+cubepm_a*omegak/omegam + omegav/omegam*cubepm_a**3))
    
    xvp(3, :) = xvp(3, :) + xvp(6, :) * 1.5 / sqrt(a * (1 + a * (1 - omega_m - omega_l) / omega_m + omega_l / omega_m * a**3))  

    call pass_particles

    if(rank == 0) then
       write(*,*) '**********************'
       write(*,*) 'Included Kaiser Effect'
       write(*,*) 'Omega_m = ', omega_m, ' a = ', a
       !write(*,*) '1/H(z) =', 1.5*sqrt(omegam/cubepm_a)
       write(*,*) '1 / H(z) = ', 1.5 / sqrt(a * (1 + a * (1 - omega_m - omega_l) / omega_m + omega_l / omega_m * a**3))
       write(*,*) '**********************'
    endif
#endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished read_particles ... elapsed time = ", time2-time1

    return
    
end subroutine read_particles

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

    write(*,*) 'nodes   ', nodes
    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*)

    return
  end subroutine writeparams

! -------------------------------------------------------------------------------------------------------

subroutine writepowerspectra(command)
    !
    ! Writes the dimensionless power spectrum for the curl/divergence of the momentum density field
    !    

    implicit none
    
    integer      :: i, j, k
    character*180 :: fn
    character*7  :: prefix
    character*7  :: z_write
    real    :: vsim2phys, zcur
    integer(4) :: command

    !! Determine conversion factor for sim velocity to physical
    zcur      = z_checkpoint(cur_checkpoint)
    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

    if (rank == 0) write(*, *) "---> vsim2phys = ", vsim2phys

    !
    ! Determine name of output file
    !

    write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
    z_write=adjustl(z_write)
    
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
   fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD.dat' 
#else
   fn=output_path//z_write(1:len_trim(z_write))//prefix//'.dat' 
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
            pkdm(3, i) = pkcurldm(1, 3, i)
        enddo

    else if (command == 1) then

        do i = 1, nc
            pkdm(1, i) = pkdivdm(1, i)
            pkdm(2, i) = pkdivdm(2, i)
            pkdm(3, i) = pkdivdm(3, i)
        enddo

    else if (command == 2) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkmomdim(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkmomdim(j, 2, i)
            enddo
            pkdm(3, i) = pkmomdim(1, 3, i)
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

#ifdef NGP
       write(11,*) pkdm(3,k-1), pkdm(1:2,k-1)
#else
       write(11,*) pkdm(3,k), pkdm(1:2,k)
#endif

    enddo
    close(11)

    return

end subroutine writepowerspectra

! -------------------------------------------------------------------------------------------------------

subroutine darkmatter(command)

    implicit none

    integer :: i, j, k
    integer :: i1, j1, k1
    real    :: d, dmin, dmax, sum_dm, sum_dm_local, dmint, dmaxt
    real*8  :: dsum, dvar, dsumt, dvart
    real, dimension(3) :: dis

    integer(4) :: command ! 0 for curl, 1 for divergence, 2 for momentum

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize FFT array to zero 
    !

    do k = 1, nc_node_dim
        cube(:, :, k) = 0.
    enddo

    !
    ! Assign data to density grid
    !

    if (command == 0) then !! Fill with given curl component 

        cube(:, :, :) = momcurl(cur_dimension, :, :, :) 

    else if (command == 1) then !! Fill with divergence 

        cube(:, :, :) = momdiv(:, :, :) 

    else if (command == 2) then !! Fill with momentum density field

        cube(:, :, :) = momden(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)

    endif

    !
    ! Calculate some statistics
    !

    sum_dm_local = sum(cube) 
    call mpi_reduce(sum_dm_local, sum_dm, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
    if (rank == 0) print *, "CUBE total sum = ", sum_dm, " command = ", command

    dmin = 0
    dmax = 0
    dsum = 0
    dvar = 0

    do k = 1, nc_node_dim
       do j = 1, nc_node_dim
          do i = 1, nc_node_dim
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
      write(*,*) 'Darkmatter command ', command
      write(*,*) 'Cube min    ', dmint
      write(*,*) 'Cube max    ', dmaxt
      write(*,*) 'Cube sum ', real(dsum)
      write(*,*) 'Cube var ', real(dvar)
      write(*,*)

    endif

    ! 
    ! Forward FFT dm delta field
    !    

    call cp_fftw(1)

    !
    ! Compute power spectrum
    !
    
    if (command == 0) then
 
        call powerspectrum(slab, pkcurldm(cur_dimension, :, :), command)

    else if (command == 1) then

        call powerspectrum(slab, pkdivdm, command)

    else if (command == 2) then

        call powerspectrum(slab, pkmomdim(cur_dimension, :, :), command)

    endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished darkmatter ... elapsed time = ", time2-time1

    return

end subroutine darkmatter

! -------------------------------------------------------------------------------------------------------

subroutine pass_particles
    !
    ! Pass particles inside buffer space to their appropriate nodes.
    !
    
    implicit none

    integer i,pp,np_buf,np_exit,np_final,npo,npi
    real x(6),lb,ub
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Identify particles within the buffer
    !

    lb = 0.
    ub = real(nc_node_dim)

    np_buf = 0
    pp = 1

    do
    
        if (pp > np_local) exit

        !! Read its position  
        x = xvp(:, pp)
        
        !! See if it lies within the buffer
        if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
            x(3) < lb .or. x(3) >= ub ) then
       
            !write (*,*) 'PARTICLE OUT', xvp(:, pp)
        
            !! Make sure we aren't exceeding the maximum
            np_buf = np_buf + 1
        
            if (np_buf > np_buffer) then
                print *, rank, 'np_buffer =', np_buffer, 'exceeded - np_buf =', np_buf
                call mpi_abort(mpi_comm_world, ierr, ierr)
            endif 

            xp_buf(:, np_buf) = xvp(:, pp)
            xvp(:, pp)        = xvp(:, np_local)
            np_local          = np_local - 1
        
            cycle 
      
        endif
      
        pp = pp + 1
    
    enddo
 
    call mpi_reduce(np_buf, np_exit, 1, mpi_integer, mpi_sum, 0, &
                    mpi_comm_world, ierr) 

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'np_exit = ', np_buf
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    if (rank == 0) print *,'Total exiting particles = ',np_exit

    !
    ! Pass +x
    !

    !! Find particles that need to be passed

    tag = 11 
    npo = 0
    pp  = 1
    do 
      if (pp > np_buf) exit
      if (xp_buf(1, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'np_out=', npo
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(6), &
                              tag, cart_neighbor(5), tag, mpi_comm_world, &
                              status, ierr) 

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf + pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(1, np_buf + pp) = max(xp_buf(1, np_buf+pp) - ub, lb)
    enddo

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'x+ np_local=', np_local
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    pp = 1

    do 
      if (pp > npi) exit 
      x = xp_buf(:, np_buf + pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle
      endif
      pp = pp + 1
    enddo
   
    np_buf = np_buf + npi

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'x+ np_exit=', np_buf, np_local
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    !
    ! Pass -x
    !

    tag = 12
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(1, pp) < lb) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle 
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(5), &
                              tag, cart_neighbor(6), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(1, np_buf+pp) = min(xp_buf(1,np_buf+pp) + ub, ub-eps)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle 
      endif
      pp = pp + 1
    enddo
  
    np_buf = np_buf + npi

    !
    ! Pass +y
    !

    tag = 13 
    npo = 0
    pp  = 1
    do 
      if (pp > np_buf) exit
      if (xp_buf(2, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle 
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(4), &
                              tag, cart_neighbor(3), tag, mpi_comm_world, &
                              status, ierr) 

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(2, np_buf+pp) = max(xp_buf(2, np_buf+pp)-ub, lb)
    enddo

    pp = 1
    do 
      if (pp > npi) exit 
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi-1
        cycle 
      endif
      pp = pp + 1
    enddo
   
    np_buf = np_buf + npi

    !
    ! Pass -y
    !

    tag = 14
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(2,pp) < lb) then
        npo = npo+1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(3), &
                              tag, cart_neighbor(4), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(2, np_buf+pp) = min(xp_buf(2, np_buf+pp)+ub, ub-eps)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local+1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp = pp + 1
    enddo
  
    np_buf = np_buf + npi

    !
    ! Pass +z
    !

    tag = 15 
    npo = 0
    pp  = 1
    do 
      if (pp > np_buf) exit
      if (xp_buf(3, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle 
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(2), &
                              tag,cart_neighbor(1),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*6+1:pp*6)
      xp_buf(3,np_buf+pp)=max(xp_buf(3,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

    !
    ! Pass -z
    !

    tag=16
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(3,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*6+1:npo*6)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*6+1:pp*6)
      xp_buf(3,np_buf+pp)=min(xp_buf(3,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:,np_local)=x
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

    if (rank == 0) print *,'total buffered particles =',np_exit

    call mpi_reduce(np_local,np_final,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
      print *,'total particles =',np_final
      if (np_final /= np**3) then
        print *,'ERROR: total number of particles incorrect after passing'
      endif
    endif
 
!!  Check for particles out of bounds

    do i=1,np_local
      if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
          xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
          xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
        print *,'particle out of bounds',rank,i,xvp(:3,i),nc_node_dim
      endif
    enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished pass_particles ... elapsed time = ", time2-time1

    return

end subroutine pass_particles

! -------------------------------------------------------------------------------------------------------

subroutine mass_density
    !
    ! Bin particles in position space to generate mass density field 
    ! 

    implicit none

    real, parameter :: mp = (ncr / np)**3

    integer :: i, j, i1, i2, j1, j2, k1, k2
    real    :: x, y, z, dx1, dx2, dy1, dy2, dz1, dz2

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    do i = 1, np_local

        !! Read particle position
        x = xvp(1, i) - 0.5
        y = xvp(2, i) - 0.5
        z = xvp(3, i) - 0.5

        !! Determine particle grid location
        i1 = floor(x) + 1
        i2 = i1 + 1
        dx1 = i1 - x
        dx2 = 1 - dx1
        j1 = floor(y) + 1
        j2 = j1 + 1
        dy1 = j1 - y
        dy2 = 1 - dy1
        k1 = floor(z) + 1
        k2 = k1 + 1
        dz1 = k1 - z
        dz2 = 1 - dz1

        if (i1 < 0 .or. i2 > nc_node_dim+1 .or. j1 < 0 .or. &
            j2 > nc_node_dim+1 .or. k1 < 0 .or. k2 > nc_node_dim+1) then
                print *,'WARNING: Particle out of bounds', i1, i2, j1, j2, k1, k2, nc_node_dim
        endif

       dz1 = mp * dz1
       dz2 = mp * dz2

       massden(i1, j1, k1) = massden(i1, j1, k1) + dx1 * dy1 * dz1
       massden(i2, j1, k1) = massden(i2, j1, k1) + dx2 * dy1 * dz1
       massden(i1, j2, k1) = massden(i1, j2, k1) + dx1 * dy2 * dz1
       massden(i2, j2, k1) = massden(i2, j2, k1) + dx2 * dy2 * dz1
       massden(i1, j1, k2) = massden(i1, j1, k2) + dx1 * dy1 * dz2
       massden(i2, j1, k2) = massden(i2, j1, k2) + dx2 * dy1 * dz2
       massden(i1, j2, k2) = massden(i1, j2, k2) + dx1 * dy2 * dz2
       massden(i2, j2, k2) = massden(i2, j2, k2) + dx2 * dy2 * dz2

    enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished mass_density ... elapsed time = ", time2-time1

    return

end subroutine mass_density

! -------------------------------------------------------------------------------------------------------

subroutine momentum_density
    !
    ! Bin particles in position space to generate the 3D momentum density field
    ! 

    implicit none

    real, parameter :: mp = (ncr / np)**3

    integer :: i, j, i1, i2, j1, j2, k1, k2
    real    :: x, y, z, dx1, dx2, dy1, dy2, dz1, dz2
    real    :: dv1, dv2, v(3)
    real    :: vsim2phys, zcur
    real    :: vrmsx, vrmsy, vrmsz

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Determine RMS velocity in each dimension
    !

    vrmsx = 0.
    vrmsy = 0.
    vrmsz = 0.

    do i = 1, np_local

        !! Read particle position
        x = xvp(1, i) - 0.5 
        y = xvp(2, i) - 0.5
        z = xvp(3, i) - 0.5

        !! Determine particle grid location
        i1 = floor(x) + 1
        i2 = i1 + 1
        dx1 = i1 - x
        dx2 = 1 - dx1
        j1 = floor(y) + 1
        j2 = j1 + 1
        dy1 = j1 - y
        dy2 = 1 - dy1
        k1 = floor(z) + 1
        k2 = k1 + 1
        dz1 = k1 - z
        dz2 = 1 - dz1

        if (i1 < 0 .or. i2 > nc_node_dim+1 .or. j1 < 0 .or. &
            j2 > nc_node_dim+1 .or. k1 < 0 .or. k2 > nc_node_dim+1) then
                print *,'WARNING: Particle out of bounds', i1, i2, j1, j2, k1, k2, nc_node_dim
        endif

        !! Read particle velocities 
        v(1) = xvp(4, i)
        v(2) = xvp(5, i)
        v(3) = xvp(6, i)

        vrmsx = vrmsx + v(1)**2
        vrmsy = vrmsy + v(2)**2
        vrmsz = vrmsz + v(3)**2

        !! Momentum density field in each dimension
        do j = 1, 3

            dv1 = dz1 * mp * v(j)
            dv2 = dz2 * mp * v(j)

            momden(j, i1, j1, k1) = momden(j, i1, j1, k1) + dx1 * dy1 * dv1
            momden(j, i2, j1, k1) = momden(j, i2, j1, k1) + dx2 * dy1 * dv1
            momden(j, i1, j2, k1) = momden(j, i1, j2, k1) + dx1 * dy2 * dv1
            momden(j, i2, j2, k1) = momden(j, i2, j2, k1) + dx2 * dy2 * dv1
            momden(j, i1, j1, k2) = momden(j, i1, j1, k2) + dx1 * dy1 * dv2
            momden(j, i2, j1, k2) = momden(j, i2, j1, k2) + dx2 * dy1 * dv2
            momden(j, i1, j2, k2) = momden(j, i1, j2, k2) + dx1 * dy2 * dv2
            momden(j, i2, j2, k2) = momden(j, i2, j2, k2) + dx2 * dy2 * dv2

        enddo

    enddo

    !
    ! Determine root mean square velocities
    !

    ! First get conversion factor to physical velocity in km/s
    zcur      = z_checkpoint(cur_checkpoint)
    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc !! Takes h = 1

    vrmsx = vsim2phys * sqrt(vrmsx / np_local)
    vrmsy = vsim2phys * sqrt(vrmsy / np_local)
    vrmsz = vsim2phys * sqrt(vrmsz / np_local)

    write(*, *) "RMS physical velocities (rank, z, vx, vy, vz): ", rank, zcur, vrmsx, vrmsy, vrmsz

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished momentum_density ... elapsed time = ", time2-time1

    return

end subroutine momentum_density

! -------------------------------------------------------------------------------------------------------

subroutine momentum2velocity
    !
    ! Converts momentum density field v*(1+delta) to velocity field by dividing
    ! by density field (1+delta). Both fields are first smoothed with a Gaussian
    ! filter.
    !

    implicit none

    integer :: m, i, j, k
    integer(8) :: ecount_local, ecount

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! First smooth the density field
    !

    cube(:, :, :) = massden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) 
    call Gaussian_filter 
    massden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) = cube(:, :, :)     

    !
    ! Now smooth each component of the momentum field
    !

    do m = 1, 3

        cube(:, :, :) = momden(m, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
        call Gaussian_filter
        momden(m, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) = cube(:, :, :)

    enddo

    !
    ! Divide the latter by the former to get the velocity field
    !

    ecount_local = 0

    do i = 1, nc_node_dim
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                if (massden(i, j, k) .ne. 0.) then
                    do m = 1, 3
                        momden(m, i, j, k) = momden(m, i, j, k) / &
                                             massden(i, j, k)
                    enddo
                else
                    ecount_local = ecount_local + 1
                endif
            enddo
        enddo
    enddo

    call mpi_reduce(ecount_local, ecount, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    time2 = mpi_wtime(ierr)

    if (rank == 0) write(*, *) "Empty cells in the velocity field = ", ecount
    if (rank == 0) write(*, *) "Finished momentum2velocity ... elapsed time = ", time2-time1

    return

end subroutine momentum2velocity

! -------------------------------------------------------------------------------------------------------

subroutine Gaussian_filter
    !
    ! Smooths the provided real space field (in cube) with a Gaussian filter in
    ! Fourier space (via slab) and transforms back to real space (through cube).
    !

    implicit none

    integer :: i, j, k, kg
    real :: fr, kr, kx, ky, kz

    !! Forward transform cube to slab
    call cp_fftw(1)

    !! Convolve Fourier field with a Gaussian filter

    do k = 1, nc_slab
        kg = k + nc_slab*rank
        if (kg .lt. hc+2) then
            kz = kg - 1
        else
            kz = kg - 1 - nc
        endif
        do j = 1, nc
            if (j .lt. hc+2) then
                ky = j - 1
            else
                ky = j - 1 - nc
            endif
            do i = 1, nc+2, 2
                kx = (i-1)/2
                kr = sqrt(kx**2 + ky**2 + kz**2)
                fr = 2. * pi * kr / box

                slab(i:i+1,j,k) = slab(i:i+1,j,k) * exp(-0.5*(fr*R_smooth)**2)
         
            enddo
        enddo
    enddo

    !! Backward transform slab to cube
    call cp_fftw(-1)

    return

end subroutine Gaussian_filter

! -------------------------------------------------------------------------------------------------------

subroutine buffer_momdensity 
    !
    ! Accumulate buffer from adjacent nodes into physical volume.
    !

    implicit none

    integer :: i
    integer, dimension(mpi_status_size) :: status, sstatus, rstatus
    integer :: tag, srequest, rrequest, sierr, rierr, ierr
    integer, parameter :: num2send = (nc_node_dim + 2)**2

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    do i = 1, 3

        !
        ! Pass +x
        ! 

        tag = 111

        momden_send_buff(:, :) = momden(i, nc_node_dim+1, :, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, 1, :, :) = momden(i, 1, :, :) + momden_recv_buff(:, :)

        !
        ! Pass -x
        ! 

        tag = 112

        momden_send_buff(:, :) = momden(i, 0, :, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, nc_node_dim, :, :) = momden(i, nc_node_dim, :, :) + momden_recv_buff(:, :)

        !
        ! Pass +y
        ! 

        tag = 113

        momden_send_buff(:, :) = momden(i, :, nc_node_dim+1, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, 1, :) = momden(i, :, 1, :) + momden_recv_buff(:, :)

        !
        ! Pass -y
        ! 

        tag = 114

        momden_send_buff(:, :) = momden(i, :, 0, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, nc_node_dim, :) = momden(i, :, nc_node_dim, :) + momden_recv_buff(:, :)

        !
        ! Pass +z
        ! 

        tag = 115

        momden_send_buff(:, :) = momden(i, :, :, nc_node_dim+1)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, :, 1) = momden(i, :, :, 1) + momden_recv_buff(:, :)

        !
        ! Pass -z
        ! 

        tag = 116

        momden_send_buff(:, :) = momden(i, :, :, 0)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, :, nc_node_dim) = momden(i, :, :, nc_node_dim) + momden_recv_buff(:, :)

    enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished buffer_momdensity ... elapsed time = ", time2-time1

    return

end subroutine buffer_momdensity

! -------------------------------------------------------------------------------------------------------

subroutine buffer_massdensity
    !
    ! Accumulate buffer from adjacent nodes into physical volume.
    !

    implicit none

    integer, dimension(mpi_status_size) :: status, sstatus, rstatus
    integer :: tag, srequest, rrequest, sierr, rierr, ierr
    integer, parameter :: num2send = (nc_node_dim + 2)**2

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Pass +x
    ! 

    tag = 111

    massden_send_buff(:, :) = massden(nc_node_dim+1, :, :)

    call mpi_isend(massden_send_buff, num2send, mpi_real, cart_neighbor(6), &
               tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(massden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
               tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    massden(1, :, :) = massden(1, :, :) + massden_recv_buff(:, :)

    !
    ! Pass -x
    ! 

    tag = 112

    massden_send_buff(:, :) = massden(0, :, :)

    call mpi_isend(massden_send_buff, num2send, mpi_real, cart_neighbor(5), &
               tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(massden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
               tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    massden(nc_node_dim, :, :) = massden(nc_node_dim, :, :) + massden_recv_buff(:, :)

    !
    ! Pass +y
    ! 

    tag = 113

    massden_send_buff(:, :) = massden(:, nc_node_dim+1, :)

    call mpi_isend(massden_send_buff, num2send, mpi_real, cart_neighbor(4), &
               tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(massden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
               tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    massden(:, 1, :) = massden(:, 1, :) + massden_recv_buff(:, :)

    !
    ! Pass -y
    ! 

    tag = 114

    massden_send_buff(:, :) = massden(:, 0, :)

    call mpi_isend(massden_send_buff, num2send, mpi_real, cart_neighbor(3), &
               tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(massden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
               tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    massden(:, nc_node_dim, :) = massden(:, nc_node_dim, :) + massden_recv_buff(:, :)

    !
    ! Pass +z
    ! 

    tag = 115

    massden_send_buff(:, :) = massden(:, :, nc_node_dim+1)

    call mpi_isend(massden_send_buff, num2send, mpi_real, cart_neighbor(2), &
               tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(massden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
               tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    massden(:, :, 1) = massden(:, :, 1) + massden_recv_buff(:, :)

    !
    ! Pass -z
    ! 

    tag = 116

    massden_send_buff(:, :) = massden(:, :, 0)

    call mpi_isend(massden_send_buff, num2send, mpi_real, cart_neighbor(1), &
               tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(massden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
               tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    massden(:, :, nc_node_dim) = massden(:, :, nc_node_dim) + massden_recv_buff(:, :)

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished buffer_massdensity ... elapsed time = ", time2-time1

    return

end subroutine buffer_massdensity

! -------------------------------------------------------------------------------------------------------

subroutine pass_momdensity
    !
    ! Pass physical boundaries to adjacent nodes for finite differencing later
    !

    implicit none

    integer :: i
    integer, dimension(mpi_status_size) :: status, sstatus, rstatus
    integer :: tag, srequest, rrequest, sierr, rierr, ierr
    integer, parameter :: num2send = (nc_node_dim + 2)**2

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    do i = 1, 3

        !
        ! Pass +x
        ! 

        tag = 111

        momden_send_buff(:, :) = momden(i, nc_node_dim, :, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, 0, :, :) = momden_recv_buff(:, :)

        !
        ! Pass -x
        ! 

        tag = 112

        momden_send_buff(:, :) = momden(i, 1, :, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, nc_node_dim+1, :, :) = momden_recv_buff(:, :)

        !
        ! Pass +y
        ! 

        tag = 113

        momden_send_buff(:, :) = momden(i, :, nc_node_dim, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, 0, :) = momden_recv_buff(:, :)

        !
        ! Pass -y
        ! 

        tag = 114

        momden_send_buff(:, :) = momden(i, :, 1, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, nc_node_dim+1, :) = momden_recv_buff(:, :)

        !
        ! Pass +z
        ! 

        tag = 115

        momden_send_buff(:, :) = momden(i, :, :, nc_node_dim)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, :, 0) = momden_recv_buff(:, :)

        !
        ! Pass -z
        ! 

        tag = 116

        momden_send_buff(:, :) = momden(i, :, :, 1)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(i, :, :, nc_node_dim+1) = momden_recv_buff(:, :)

    enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished pass_momdensity ... elapsed time = ", time2-time1

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

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Compute derivatives using dy/dx(a) = y(a+1) - y(a-1) / x(a+1) - x(a-1)
    ! 

    !! Spacing between grid points
    dx = 2. * box / nc

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

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished momentum_curl ... elapsed time = ", time2-time1

    return

end subroutine momentum_curl

! -------------------------------------------------------------------------------------------------------

subroutine momentum_divergence
    !
    ! Calculates the divergence of the momentum density field
    !

    implicit none
    integer :: k, j, i
    real :: dx

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Compute derivatives using dy/dx(a) = y(a+1) - y(a-1) / x(a+1) - x(a-1)
    !

    !! Spacing between grid points
    dx = 2. * box / nc

    do i = 1, nc_node_dim
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim

                momdiv(i, j, k) = (momden(1, i+1, j, k) - momden(1, i-1, j, k) + &
                                    momden(2, i, j+1, k) - momden(2, i, j-1, k) + &
                                    momden(3, i, j, k+1) - momden(3, i, j, k-1)) / dx

            enddo
        enddo
    enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished momentum_divergence ... elapsed time = ", time2-time1

    return

end subroutine momentum_divergence

! -------------------------------------------------------------------------------------------------------

  subroutine powerspectrum(delta, pk, command)
    implicit none
    real, dimension(3, nc)       :: pk
    real, dimension(nc+2, nc, nc_slab) :: delta
    integer :: command

    integer :: i, j, k, kg
    integer :: k1, k2
    real    :: kr, kx, ky, kz, w1, w2, pow, x, y, z, sync_x, sync_y, sync_z, kernel
    real, dimension(3, nc, nc_slab) :: pkt
    real, dimension(3, nc) :: pktsum
    real, dimension(nc) :: kcen, kcount
    real    :: kavg

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    pkt = 0.0
    pktsum = 0.0

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
#ifdef NGP
                w1=1
                w2=0
#endif                
                pow=sum((delta(i:i+1,j,k)/ncr**3)**2)/kernel**4
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

    !! NOTE: Binning the Fourier transform of the curl/divergence of the
    !! momentum field introduces a factor of k^2 over what would be obtained
    !! from separting the parallel and perpendicular components of the
    !! transformed field in Fourier space. We must therefore divide by k^2 when
    !! when constructing the power spectrum in this way.

    if (rank == 0) then
        do k=1,nc
            if (pktsum(3,k) .eq. 0) then
                pk(:,k)=0
            else
                pk(1:2,k)=pktsum(1:2,k)/pktsum(3,k)
                pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))

                kavg = kcen(k) / kcount(k)
                pk(3,k) = 2. * pi * kavg / box

#ifdef NGP
                pk(1:2,k)=4*pi*(kavg)**3*pk(1:2,k)
#else
                pk(1:2,k)=4*pi*(kavg-1.)**3*pk(1:2,k)
#endif

                !! Divide by k^2 for divergence and curl components 
                if (command == 0 .or. command == 1) then

                    pk(1:2, k) = pk(1:2, k) / pk(3, k)**2

                endif

            endif
        enddo
    endif

    call mpi_bcast(pk,3*nc,mpi_real,0,mpi_comm_world,ierr)

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished powerspectrum ... elapsed time = ", time2-time1

    return

  end subroutine powerspectrum

! -------------------------------------------------------------------------------------------------------

subroutine writevelocityfield

    implicit none

    integer :: m, i, j, k, fstat
    character(len=180) :: fn
    character(len=7)   :: z_write
    character(len=4)   :: rank_string
    character(len=1)   :: dim_string
    real :: vsim2phys, zcur

    !
    ! Determine conversion to proper velocity [km/s]
    !

    zcur      = z_checkpoint(cur_checkpoint)
    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

    if (rank == 0) write(*,*) "zcur = ", zcur, "vsim2phys = ", vsim2phys

    !
    ! Checkpoint and rank strings
    !

    write(z_write, '(f7.3)') z_checkpoint(cur_checkpoint)
    z_write = adjustl(z_write)

    write(rank_string, '(i4)') rank
    rank_string = adjustl(rank_string)

    !
    ! Write out velocity field for each dimension
    !

    do m = 1, 3

        if (m == 1) dim_string = "x"
        if (m == 2) dim_string = "y"
        if (m == 3) dim_string = "z"

        fn = output_path//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//&
             rank_string(1:len_trim(rank_string))//".bin"

        open(unit=11, file=fn, status="replace", iostat=fstat, form="binary")

        do k = 1, nc_node_dim
            do j = 1, nc_node_dim

                write(11) momden(m, 1:nc_node_dim, j, k) * vsim2phys

            enddo
        enddo

        close(11)

    enddo

    return

end subroutine writevelocityfield

! -------------------------------------------------------------------------------------------------------

end program cic_velpower 
