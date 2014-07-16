!!
!! cic_velpower_dmnu.f90
!!
!! Program to compute the velocity power spectra of both dark matter and neutrinos as well as their 
!! cross velocity power spectra. This is to be used in conjunction with particle checkpoint files
!! produced by a -DNEUTRINOS cubep3m simulation.
!!
!! * Using FFTW on the SciNet GPC compile with:
!!   mpif90 -shared-intel -fpp -g -O3 -DNGP -DGAUSSIAN_SMOOTH -mt_mpi cic_velpower_dmnu.f90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_velpower_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Using MKL on the SciNet GPC compile with:
!!   mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -DGAUSSIAN_SMOOTH -mt_mpi 
!!        cic_velpower_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_velpower_dmnu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64
!!
!! * Using FFTW on the SciNet BGQ compile with:
!!   mpif90 -q64 -O3 -qhot -qarch=qp -qtune=qp -WF,-DNGP,-DGAUSSIAN_SMOOTH cic_velpower_dmnu.F90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_velpower_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Optional flags:
!!   -DNGP: Uses NGP interpolation for binning of the power spectrum. 
!!   -DSLAB: Alternatively run with FFTW slab decomposition instead of P3DFFT pencil decomposition.
!!   -DMOMENTUM: Uses the approach of computing the velocity field by dividing the momentum density
!!               field by the density field. If this flag is not used then the velocity field is computed
!!               by taking the average velocity of the nearest particles to each grid centre.
!!   -DGAUSSIAN_SMOOTH: Smooths the density field with a Gaussian filter to avoid division by empty 
!!                      cells when converting from momentum density to velocity. 
!!   -Dwrite_vel: Writes gridded velocity fields to binary files. 
!!   -DKAISER: Adjusts for redshift space distortions.
!!   -DDEBUG: Output useful debugging information.

program cic_velpower 
#ifndef MOMENTUM
  use omp_lib
#endif

  implicit none

  include 'mpif.h'
  include '../../parameters'

  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints'

#ifndef MOMENTUM
  !! Threading
  integer(4), parameter :: nt = 8

  !! Number of nearest particles from grid centre to determine average velocity
  integer, parameter :: N_closest_nu = 30
  integer, parameter :: N_closest_dm = 4
  integer, parameter :: N_closest_auto = 1
#endif

#ifdef GAUSSIAN_SMOOTH
  !! Gaussian smoothing parameters (to go from momentum density field to velocity field)
  real, parameter :: cell_smooth = 1.0
  real, parameter :: R_smooth = box / nc * cell_smooth
#endif

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc
  real, parameter    :: npr=np

#ifdef LOGBIN
  integer, parameter :: numbins = hc / 48 
#endif

  !! internals
  integer, parameter :: max_checkpoints=100
  real, dimension(max_checkpoints) :: z_checkpoint
  integer num_checkpoints, cur_checkpoint

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: np_node_dim = np/nodes_dim
  integer(4), parameter :: max_np = density_buffer * ( ((nf_tile-2*nf_buf)*tiles_node_dim/2)**3 + &
                                  (8*nf_buf**3 + 6*nf_buf*(((nf_tile-2*nf_buf)*tiles_node_dim)**2) + &
                                  12*(nf_buf**2)*((nf_tile-2*nf_buf)*tiles_node_dim))/8.0 )
  integer(4), parameter :: np_buffer=int(2./3.*max_np)
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! For storage of dark matter particles (usually small since r_n_1_3 generally > 1)
  integer(4), parameter :: max_np_dm = max_np / ratio_nudm_dim**3

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local,np_local_dm
  integer(8) :: plan, iplan
  logical :: firstfftw

  !! have velocity power spectra for each x, y, z
  integer cur_dimension

  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(6,max_np_dm) :: xvp_dm
  real, dimension(6,np_buffer) :: xp_buf
  real, dimension(6*np_buffer) :: send_buf, recv_buf

  !! Power spectrum arrays
  real, dimension(3, 3, nc) :: pkmomdim
  real, dimension(3, 3, nc) :: pkmomdim_dm
  real, dimension(3, 3, nc) :: pkmomdim_dmnu
  real, dimension(3, 3, nc) :: pkmomdim_rel
  real, dimension(3, nc) :: pkdm

  !! For pencils decomposition:
#ifdef SLAB
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
#else
  integer(4), parameter   :: nodes_pen = nodes_dim
  integer(4), parameter   :: nc_pen = nc_node_dim / nodes_dim
  integer(4), parameter   :: dim_y = nodes_dim
  integer(4), parameter   :: dim_z = nodes_dim**2
  integer(4) :: pen_dims(2), istart(3), iend(3), isize(3), fstart(3), fend(3), fsize(3), mypadd
  integer(4), dimension(0:nodes_dim-1) :: pen_neighbor_to
  integer(4), dimension(0:nodes_dim-1) :: pen_neighbor_fm
#endif

  !! Fourier transform arrays
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
#ifdef SLAB
  real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2,nc,nc_slab) :: slab, slab2, slab3, slab_work
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab, slab2, slab3
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1) :: recv_cube
#endif

#ifdef MOMENTUM
  !! Array containing the matter density field
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: massden
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: massden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: massden_recv_buff
#else
  !! Parameters for linked list
  integer(4), parameter :: nfine_buf = 8
  integer(4), parameter :: mesh_scale = 2
  integer(4), parameter :: nc_buf = nfine_buf / mesh_scale
  integer(4), parameter :: nm_node_dim = nc_node_dim / mesh_scale
  integer(4), parameter :: hoc_nc_l = 1 - nc_buf
  integer(4), parameter :: hoc_nc_h = nm_node_dim + nc_buf
  integer(4), parameter :: hoc_pass_depth = 2*nc_buf
  real(4), parameter    :: rnf_buf = real(nfine_buf)
  integer(4), parameter :: num_ngbhs = (2*nc_buf+1)**3
  integer(4), parameter :: max_npart_search = int(num_ngbhs*real(max_np)/real(nc_node_dim)**2)

  !! Linked list and head-of-chain arrays
  integer(4), dimension(max_np) :: ll
  integer(4), dimension(max_np_dm) :: ll_dm
  integer(4) :: hoc(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h)
  integer(4) :: hoc_dm(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h)
  integer(4) :: ipos(nt, max_npart_search)
  real(4) :: rpos(nt, 2, max_npart_search)
  integer(4) :: cell_search(num_ngbhs, 3)
  logical :: done_shell(num_ngbhs)

  !! Array storing what group each particle belongs to (for autocorrelation break each species 
  !! into two groups and then cross-correlate to remove shot noise)
  integer(1) :: GID(max_np), GID_dm(max_np_dm)
  integer(1) :: GID_buf(np_buffer), GID_send_buf(np_buffer), GID_recv_buf(np_buffer)

  integer(1), parameter :: g0 = 0
  integer(1), parameter :: g1 = 1
#endif

  !! Only work with one component of the velocity field at a time 
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momden
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: momden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: momden_recv_buff

  !! For threading 
  integer :: kt
#ifdef SLAB
  integer, parameter :: kt_stop = nc_slab
#else
  integer, parameter :: kt_stop = nc_pen+2
#endif

  !! Equivalence arrays to save memory
  !! Must make sure that ONLY the largest array in each equivalence statement is on the common block. 
  !! Fill in check_equivalence.py with simulation parameters and run if unsure.
  equivalence(recv_cube, xp_buf) !! May sometimes need to be changed (see check_equivalence.py)
  equivalence(momden, send_buf) !! May sometimes need to be changed (see check_equivalence.py)
#ifdef SLAB
  equivalence(slab_work, recv_buf) !! May sometimes need to be changed (see check_equivalence.py)
#endif
  equivalence(slab, cube) !! Slab will always be larger than cube
#ifdef MOMENTUM
  equivalence(massden_send_buff, momden_send_buff) !! These are the same size
  equivalence(massden_recv_buff, momden_recv_buff) !! These are the same size
#endif

  !! Common block
#ifdef SLAB
#ifdef MOMENTUM
  common xvp, xvp_dm, massden, recv_cube, momden, slab_work, slab, slab2, slab3, massden_send_buff, massden_recv_buff, pkmomdim, pkmomdim_dm, pkmomdim_dmnu, pkmomdim_rel, pkdm
#else
  common xvp, xvp_dm, hoc, hoc_dm, ll, ll_dm, ipos, rpos, cell_search, done_shell, recv_cube, momden, slab_work, slab, slab2, slab3, momden_send_buff, momden_recv_buff, pkmomdim, pkmomdim_dm, pkmomdim_dmnu, pkmomdim_rel, pkdm
  common /gvar/ GID, GID_dm, GID_buf, GID_send_buf, GID_recv_buf 
#endif
#else
#ifdef MOMENTUM
  common xvp, xvp_dm, massden, recv_cube, momden, recv_buf, slab, slab2, slab3, massden_send_buff, massden_recv_buff, pkmomdim, pkmomdim_dm, pkmomdim_dmnu, pkmomdim_rel, pkdm
#else
  common xvp, xvp_dm, hoc, hoc_dm, ll, ll_dm, ipos, rpos, cell_search, done_shell, recv_cube, momden, recv_buf, slab, slab2, slab3, momden_send_buff, momden_recv_buff, pkmomdim, pkmomdim_dm, pkmomdim_dmnu, pkmomdim_rel, pkdm
  common /gvar/ GID, GID_dm, GID_buf, GID_send_buf, GID_recv_buf
#endif
#endif

! -------------------------------------------------------------------------------------------------------
! MAIN
! -------------------------------------------------------------------------------------------------------

  call mpi_initialize
#ifndef MOMENTUM
  call omp_set_num_threads(nt)
  call initialize_random_number
#endif

  if (rank == 0) call writeparams

  firstfftw = .true.  ! initialize fftw so that it generates the plans

  call read_checkpoint_list

  do cur_checkpoint = 1, num_checkpoints

    call initvar

    !! Read and pass neutrino particles
    call read_particles(0)
    call pass_particles(0)
#ifndef MOMENTUM
    call assign_groups(0)
    call link_list(0)
    call buffer_particles(0)
#endif

    !! Read and pass dark matter particles (stored in xvp_dm)
    call read_particles(1)
    call pass_particles(1)
#ifndef MOMENTUM
    call assign_groups(1)
    call link_list(1)
    call buffer_particles(1)
#endif

#ifndef MOMENTUM
    !
    ! Compute power spectra for both neutrinos and dark matter by separting each
    ! species into two groups and computing the cross-spectrum of the groups.
    ! This is done to remove shot noise which does not correlate between groups.
    ! At the same time also compute the auto-power of the relative velocity field.
    ! This also has to be done with the cross of two separate groups to remove noise.
    !

    do cur_dimension = 1, 3 !! Each x, y, z dimension

        !! --------------------------------------------------------------------------------
        !! NEUTRINOS
        !! --------------------------------------------------------------------------------

        !! Determine cur_dimension component of neutrino velocity for group 0
        call velocity_density(0, g0, N_closest_auto)
        call buffer_momdensity

        !! Fourier transform velocity field for neutrino group 0 and store in slab2
        call darkmatter
        !$omp parallel do num_threads(nt) default(shared) private(kt)
        do kt = 1, kt_stop
            slab2(:,:,kt) = slab(:,:,kt)
        enddo
        !$omp end parallel do 

        !! Determine cur_dimension component of neutrino velocity for group 1
        call velocity_density(0, g1, N_closest_auto)
        call buffer_momdensity

        !! Fourier transform velocity field for neutrino group 1 and compute cross-correlation with group 0
        call darkmatter
        call powerspectrum(slab, slab2, pkmomdim(cur_dimension,:,:))

        !! --------------------------------------------------------------------------------
        !! DARK MATTER
        !! --------------------------------------------------------------------------------

        !! Determine cur_dimension component of dark matter velocity for group 0 
        call velocity_density(1, g0, N_closest_auto)
        call buffer_momdensity

        !! Fourier transform velocity field for dark matter group 0 
        call darkmatter
    
        !! slab3 now stores the relative velocity difference for group 0 (dark matter - neutrino)
        !$omp parallel do num_threads(nt) default(shared) private(kt)
        do kt = 1, kt_stop
            slab3(:,:,kt) = slab(:,:,kt) - slab2(:,:,kt)
        enddo
        !$omp end parallel do 

        !! Now store dark matter group 0 in slab2 (overwrites neutrino group 0)
        !$omp parallel do num_threads(nt) default(shared) private(kt)
        do kt = 1, kt_stop
            slab2(:,:,kt) = slab(:,:,kt)
        enddo
        !$omp end parallel do 

        !! Determine cur_dimension component of dark matter velocity for group 1
        call velocity_density(1, g1, N_closest_auto)
        call buffer_momdensity

        !! Fourier transform velocity field for group 1 and compute cross-correlation with group 0
        call darkmatter
        call powerspectrum(slab, slab2, pkmomdim_dm(cur_dimension,:,:))

        !! --------------------------------------------------------------------------------
        !! DARK MATTER - NEUTRINO RELATIVE VELOCITY
        !! --------------------------------------------------------------------------------

        !! Store dark matter group 1 in slab2
        !$omp parallel do num_threads(nt) default(shared) private(kt)
        do kt = 1, kt_stop
            slab2(:,:,kt) = slab(:,:,kt)
        enddo
        !$omp end parallel do 
        
        !! Recompute cur_dimension component of neutrino velocity for group 1
        call velocity_density(0, g1, N_closest_auto)
        call buffer_momdensity

        !! Fourier transform velocity field for neutrino group 1
        call darkmatter 

        !! Now slab2 will store the relative velocity difference for group 1 (dark matter - neutrino)
        !$omp parallel do num_threads(nt) default(shared) private(kt)
        do kt = 1, kt_stop
            slab2(:,:,kt) = slab2(:,:,kt) - slab(:,:,kt)
        enddo
        !$omp end parallel do

        !! Compute relative velocity power spectra
        call powerspectrum(slab2, slab3, pkmomdim_rel(cur_dimension,:,:))

    enddo

    !! Put all particles into group 0 for neutrino-dark matter cross-power
    call clear_groups

#endif

    do cur_dimension = 1, 3 !! Each x, y, z dimension

#ifdef MOMENTUM
        !! Compute momentum density field for this dimension
        call momentum_density(0)
        call buffer_momdensity

        !! Convert momentum field into velocity field
        call mass_density(0)
        call buffer_massdensity
        call momentum2velocity
#else
        !! Determine cur_dimension component of velocity for all particles (all in group 0 now)
        call velocity_density(0, g0, N_closest_nu)
        call buffer_momdensity
#endif

#ifdef write_vel
        call writevelocityfield(0)
#endif

        !! Fourier transform velocity field and compute power spectrum
        call darkmatter
#ifdef MOMENTUM
        call powerspectrum(slab, slab, pkmomdim(cur_dimension,:,:))
#endif

        !! Store neutrino Fourier field in slab 2 for cross spectra later
        !$omp parallel do num_threads(nt) default(shared) private(kt)
        do kt = 1, kt_stop
            slab2(:,:,kt) = slab(:,:,kt)
        enddo
        !$omp end parallel do 

#ifdef MOMENTUM
        !! Compute dark matter momentum density field for this dimension
        call momentum_density(1)
        call buffer_momdensity

        !! Convert dark matter momentum field into velocity field
        call mass_density(0)
        call buffer_massdensity
        call momentum2velocity
#else
        !! Determine cur_dimension component of velocity for all particles (all in group 0 now)
        call velocity_density(1, g0, N_closest_dm)
        call buffer_momdensity
#endif

#ifdef write_vel
        call writevelocityfield(1)
#endif

        !! Fourier transform velocity field and compute power spectrum
        call darkmatter
#ifdef MOMENTUM
        call powerspectrum(slab, slab, pkmomdim_dm(cur_dimension,:,:))
#endif

        !! Compute cross power spectra
        call powerspectrum(slab, slab2, pkmomdim_dmnu(cur_dimension,:,:))        

    enddo

    !! Write out the power spectra
    if (rank == 0) then
        call writepowerspectra(0)
        call writepowerspectra(1)
        call writepowerspectra(2)
        call writepowerspectra(3)
    endif

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
#ifdef SLAB
      do i = 0, nodes_dim - 1
        slab_neighbor(i,j) = i + j * nodes_dim + slab_coord(3) &
                           * nodes_slab
      enddo
#else
        pen_neighbor_to(j) = nodes_slab*slab_coord(3) + slab_coord(2) + j*nodes_dim
        pen_neighbor_fm(j) = nodes_slab*slab_coord(3) + j + nodes_dim*slab_coord(1)
#endif
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

    integer :: k, j, i
    integer :: ind, ibuf

    !! Particle positions and velocities
    do k = 1, max_np
       xvp(:, k) = 0.
    enddo

    !! Momentum and matter density arrays
    do k = 0, nc_node_dim + 1
#ifdef MOMENTUM
        massden(:, :, k) = 0.
        massden_send_buff(:, k) = 0.
        massden_recv_buff(:, k) = 0.
#else
        momden_send_buff(:, k) = 0.
        momden_recv_buff(:, k) = 0.
#endif
        momden(:, :, k) = 0.
    enddo

    !! Fourier transform arrays
#ifdef SLAB
    do k = 1, nc_slab
       slab_work(:, :, k)=0
    enddo
#endif

    do k = 1, nc_node_dim
       cube(:, :, k) = 0
    enddo

#ifdef SLAB
    do k = 1, nc_slab
#else
    do k = 1, nc_pen+2
#endif
       slab(:, :, k) = 0
    enddo

    do k = 1, np_buffer
       xp_buf(:, k) = 0
    enddo

    do k = 1, 6*np_buffer
        send_buf(k) = 0.
        recv_buf(k) = 0.
    enddo

    do k = 1, 3 * np_buffer
       recv_buf(k) = 0
    enddo

#ifdef SLAB
    do k = 0, nodes_slab-1
#else
    do k = 0, nodes_pen-1
#endif
        recv_cube(:, :, :, k) = 0.
    enddo

    !! Power spectrum arrays
    do k = 1, nc
        pkmomdim(:, :, k) = 0.
        pkmomdim_dm(:, :, k) = 0.
        pkmomdim_dmnu(:, :, k) = 0.
        pkmomdim_rel(:, :, k) = 0.
        pkdm(:, k) = 0.
    enddo

#ifndef MOMENTUM
    !! Array to store the order to search the relative index of neighbouring
    !cells
    ind = 1
    do ibuf = 0, nc_buf
        do k = -ibuf, ibuf
            do j = -ibuf, ibuf
                do i = -ibuf, ibuf
                    if (.not. (abs(i) < ibuf .and. abs(j) < ibuf .and. abs(k) < ibuf)) then
                        cell_search(ind, 1) = i
                        cell_search(ind, 2) = j
                        cell_search(ind, 3) = k
                        done_shell(ind) = .false.
                        ind = ind + 1
                    endif
                enddo
            enddo
        enddo
        done_shell(ind-1) = .true.
    enddo
#endif

    return

end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine read_particles(command)
    !
    ! Read x, y, z positions and velocities and store in xvp
    !

    implicit none
    
    real z_write, np_total
    integer j, fstat
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name
    integer(4) :: command

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

    if (command == 0) then
        if(z_write .eq. z_i) then
           check_name=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        else
           check_name=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        endif
    else
        if(z_write .eq. z_i) then
           check_name=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        else
           check_name=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        endif
    endif

    !! Open the file    

    open(unit=21,file=check_name,status="old",iostat=fstat,access="stream")

    !! Check for opening error
    if (fstat /= 0) then
      write(*,*) 'ERROR: Cannot open checkpoint position file'
      write(*,*) 'rank', rank, ' file: ',check_name
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    !! Read in checkpoint header data
    if (command == 0) then
        read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
                   sim_projection,sim_halofind,mass_p
    else
        read(21) np_local_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
                   sim_projection,sim_halofind,mass_p
    endif

    !! Check for memory problems
    if (command == 0) then
        if (np_local > max_np) then
          write(*,*) 'ERROR: Too many particles to store in memory!'
          write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
          call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm > max_np_dm) then
          write(*,*) 'ERROR: Too many particles to store in memory!'
          write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
          call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    !! Tally up total number of particles
    if (command == 0) then
        call mpi_reduce(real(np_local, kind=4), np_total, 1, mpi_real, &
                         mpi_sum, 0, mpi_comm_world, ierr)
    else
        call mpi_reduce(real(np_local_dm, kind=4), np_total, 1, mpi_real, &
                         mpi_sum, 0, mpi_comm_world, ierr)

    endif    

    if (rank == 0) write(*,*) 'Total number of particles = ', int(np_total,8)

    if (command == 0) then
        do j=1, np_local
            read(21) xvp(:,j)
        enddo
    else
        do j=1, np_local_dm
            read(21) xvp_dm(:,j)
        enddo
    endif

    close(21)
 
#ifdef KAISER

    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...

    if (command == 0) then
        xvp(3,:)=xvp(3,:) + xvp(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
    else
        xvp_dm(3,:)=xvp_dm(3,:) + xvp_dm(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
    endif

    call pass_particles(command)

    if(rank==0) then
       write(*,*) '**********************'
       write(*,*) 'Included Kaiser Effect'
       write(*,*) 'Omega_m =', omega_m, 'a =', a
       !write(*,*) '1/H(z) =', 1.5*sqrt(omegam/cubepm_a)
       write(*,*) '1/H(z) =', 1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
       write(*,*) '**********************'
    endif
#endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished read_particles ... elapsed time = ", time2-time1

    return
    
end subroutine read_particles

! -------------------------------------------------------------------------------------------------------

#ifdef SLAB
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
#else
subroutine pack_pencils
    !
    ! Pack cubic data into pencils for p3dfft transform.
    !

    implicit none

    integer(4) :: i,j,k,i0,i1,k1
    integer(4) :: pen_slice,tag,rtag
    integer(4) :: num_elements !possibly should be double so as not to 
                               !overflow for large runs, but integer*8 
                               !might be unsupported by MPI

    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    !
    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
    !

    breakup = 1
    num_elements = nc_node_dim * nc_node_dim * nc_pen
    passGB = 4. * num_elements / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements / breakup

    !
    ! Send the data from cube to recv_cube
    !

    do k = 1, breakup
        nc_pen_break = nc_pen/breakup*(k-1)
        do j = 0, nodes_dim - 1
            pen_slice = j
            tag  = rank**2
            rtag = pen_neighbor_fm(j)**2
            call mpi_isend(cube(1,1, pen_slice*nc_pen + nc_pen_break + 1), num_elements, &
                           mpi_real, pen_neighbor_to(j), tag, mpi_comm_world, &
                           requests(pen_slice+1),ierr)
            call mpi_irecv(recv_cube(1,1,1+nc_pen_break,pen_slice), &
                           num_elements, mpi_real, pen_neighbor_fm(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim, requests, wait_status, ierr)

    enddo

    !
    ! Place this data into the pencils (stored in the slab array)
    !

    do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        pen_slice = i

        do k = 1, nc_pen
            do j = 1, nc_node_dim
                slab(i0:i1,j,k) = recv_cube(:,j,k,pen_slice)
            enddo
        enddo

    enddo

end subroutine pack_pencils
#endif
    
! -------------------------------------------------------------------------------------------------------

#ifdef SLAB
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
#else
subroutine unpack_pencils
    !
    ! Unpack data from the pencils back into the cubic decompisition following
    ! p3dfft transform.
    !

    implicit none

    integer(4) :: i,j,k,i0,i1,k1
    integer(4) :: pen_slice,num_elements,tag,rtag
    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    !
    ! Place data in the recv_cube buffer
    !

    do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        pen_slice = i
        do k = 1, nc_pen
            do j = 1, nc_node_dim
                recv_cube(:, j, k, pen_slice) = slab(i0:i1, j, k)
            enddo
        enddo
    enddo

    !
    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
    !

    breakup = 1
    num_elements = nc_node_dim * nc_node_dim * nc_pen
    passGB = 4. * num_elements / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements / breakup

    !
    ! Put this data back into cube
    !

    do k = 1, breakup
        nc_pen_break = nc_pen/breakup*(k-1)
        do j = 0, nodes_dim - 1
            pen_slice = j
            tag  = rank**2
            rtag = pen_neighbor_to(j)**2
            call mpi_isend(recv_cube(1,1,1+nc_pen_break,pen_slice), num_elements, &
                           mpi_real, pen_neighbor_fm(j), tag, mpi_comm_world, &
                           requests(pen_slice+1),ierr)
            call mpi_irecv(cube(1,1,pen_slice*nc_pen + nc_pen_break + 1), &
                           num_elements, mpi_real, pen_neighbor_to(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim,requests, wait_status, ierr)

    enddo

end subroutine unpack_pencils
#endif

! -------------------------------------------------------------------------------------------------------

subroutine cp_fftw(command)
    !
    ! Calculate fftw transform using P3DFFT.
    ! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
    !

#ifndef SLAB
    use p3dfft
#endif
    implicit none
#ifdef SLAB
    include 'fftw_f77.i'
    integer(4), parameter :: order=FFTW_NORMAL_ORDER ! FFTW_TRANSPOSED_ORDER
#endif

    integer(4) :: i
    integer(4) :: command

#ifdef DEBUG
    do i=0,nodes-1
      if (rank == i) print *,'starting fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

! initialize plan variables for fftw

    if (firstfftw) then
#ifdef SLAB
      call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nc, &
            nc,nc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
      call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nc, &
            nc,nc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)
#else
        pen_dims = (/dim_y,dim_z/)
        call p3dfft_setup(pen_dims, nc, nc, nc, .true.)
        call p3dfft_get_dims(istart, iend, isize, 1, mypadd)
        call p3dfft_get_dims(fstart, fend, fsize, 2)
#endif
#ifdef DEBUG
      print *,'finished initialization of fftw',rank
#endif
      firstfftw=.false.
    endif

! giver

    if (command /= 0) then

!! call pack routine if we are going forward

#ifdef DEBUG
    do i=0,nodes-1
      if (rank == i) print *,'starting pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

#ifdef SLAB
      if (command > 0) call pack_slab
#else
      if (command > 0) call pack_pencils
#endif

#ifdef DEBUG
    do i=0,nodes-1
      if (rank == i) print *,'finished forward slab pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    if (command > 0) then
#ifdef SLAB
        call rfftwnd_f77_mpi(plan,1,slab,slab_work,1,order)
#else
        call ftran_r2c(slab, slab, "fft")
#endif
    else
#ifdef SLAB
        call rfftwnd_f77_mpi(iplan,1,slab,slab_work,1,order)
#else
        call btran_c2r(slab, slab, "tff")
#endif
        slab=slab/(real(nc)*real(nc)*real(nc))
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (rank == i) print *,'finished fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

!! unpack the slab data

#ifdef SLAB
      if (command < 0) call unpack_slab
#else
      if (command < 0) call unpack_pencils
#endif

    else

! if command = 0 we delete the plans

#ifdef SLAB
        call rfftwnd_f77_mpi_destroy_plan(iplan)
        call rfftwnd_f77_mpi_destroy_plan(plan)
#else
        call p3dfft_clean
#endif

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
    character*6  :: prefix
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
    prefix = 'ngpvps'
#else
    prefix = 'cicvps'
#endif

#ifdef KAISER
   if (command == 0) then !! Neutrino power spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD_nu.dat' 
   else if (command == 1) then !! Dark matter power spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD.dat'
   else if (command == 2) then !! Dark matter-neutrino cross spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD_dmnu.dat'
   else if (command == 3) then !! Dark matter-neutrino relative spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD_rel.dat'
   endif
#else
   if (command == 0) then !! Neutrino power spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'_nu.dat'
   else if (command == 1) then !! Dark matter power spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'.dat'
   else if (command == 2) then !! Dark matter-neutrino cross spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'_dmnu.dat'
   else if (command == 3) then !! Dark matter-neutrino relative spectra
       fn=output_path//z_write(1:len_trim(z_write))//prefix//'_rel.dat'
   endif
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
                pkdm(1, i) = pkdm(1, i) + pkmomdim(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkmomdim(j, 2, i)
            enddo
            pkdm(3, i) = pkmomdim(1, 3, i)
        enddo

    else if (command == 1) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkmomdim_dm(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkmomdim_dm(j, 2, i)
            enddo
            pkdm(3, i) = pkmomdim_dm(1, 3, i)
        enddo

    else if (command == 2) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkmomdim_dmnu(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkmomdim_dmnu(j, 2, i)
            enddo
            pkdm(3, i) = pkmomdim_dmnu(1, 3, i)
        enddo

    else if (command == 3) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkmomdim_rel(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkmomdim_rel(j, 2, i)
            enddo
            pkdm(3, i) = pkmomdim_rel(1, 3, i)
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

#ifdef LOGBIN
    do k = 2, numbins + 1
#else
    do k = 2, hc + 1
#endif

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


subroutine darkmatter

    implicit none

    integer :: i, j, k
    integer :: i1, j1, k1
    real    :: d, dmin, dmax, sum_dm, sum_dm_local, dmint, dmaxt
    real*8  :: dsum, dvar, dsumt, dvart
    real, dimension(3) :: dis

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

    cube(:, :, :) = momden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)

    !
    ! Calculate some statistics
    !

    sum_dm_local = sum(cube) 
    call mpi_reduce(sum_dm_local, sum_dm, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
    if (rank == 0) print *, "CUBE total sum = ", sum_dm

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
    
    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished darkmatter ... elapsed time = ", time2-time1

    return

end subroutine darkmatter

! -------------------------------------------------------------------------------------------------------

subroutine pass_particles(command)
    !
    ! Pass particles inside buffer space to their appropriate nodes.
    !
    
    implicit none

    integer i,pp,np_buf,np_exit,np_final,npo,npi
    real x(6),lb,ub
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03
    integer(4) :: command

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

        if (command == 0) then
            if (pp > np_local) exit
        else
            if (pp > np_local_dm) exit
        endif

        !! Read its position  
        if (command == 0) then
            x = xvp(:, pp)
        else
            x = xvp_dm(:, pp)
        endif
        
        !! See if it lies within the buffer
        if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
            x(3) < lb .or. x(3) >= ub ) then
       
            !! Make sure we aren't exceeding the maximum
            np_buf = np_buf + 1
        
            if (np_buf > np_buffer) then
                print *, rank, 'np_buffer =', np_buffer, 'exceeded - np_buf =', np_buf
                call mpi_abort(mpi_comm_world, ierr, ierr)
            endif 

            if (command == 0) then
                xp_buf(:, np_buf) = xvp(:, pp)
                xvp(:, pp)        = xvp(:, np_local)
                np_local          = np_local - 1
            else
                xp_buf(:, np_buf) = xvp_dm(:, pp)
                xvp_dm(:, pp)     = xvp_dm(:, np_local_dm)
                np_local_dm       = np_local_dm - 1
            endif       
 
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
      if (rank == i .and. command == 0) print *, rank, 'x+ np_local=', np_local
      if (rank == i .and. command == 1) print *, rank, 'x+ np_local_dm=', np_local_dm
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    pp = 1

    do 
      if (pp > npi) exit 
      x = xp_buf(:, np_buf + pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        if (command == 0) then
            np_local = np_local + 1
            xvp(:, np_local) = x
        else
            np_local_dm = np_local_dm + 1
            xvp_dm(:, np_local_dm) = x
        endif
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle
      endif
      pp = pp + 1
    enddo
   
    np_buf = np_buf + npi

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i .and. command == 0) print *, rank, 'x+ np_exit=', np_buf, np_local
      if (rank == i .and. command == 1) print *, rank, 'x+ np_exit=', np_buf, np_local_dm
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
        if (command == 0) then
            np_local = np_local + 1
            xvp(:, np_local) = x
        else
            np_local_dm = np_local_dm + 1
            xvp_dm(:, np_local_dm) = x
        endif
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
        if (command == 0) then
            np_local = np_local + 1
            xvp(:, np_local) = x
        else
            np_local_dm = np_local_dm + 1
            xvp_dm(:, np_local_dm) = x
        endif
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
        if (command == 0) then
            np_local = np_local+1
            xvp(:, np_local) = x
        else
            np_local_dm = np_local_dm+1
            xvp_dm(:, np_local_dm) = x
        endif
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
        if (command == 0) then
            np_local=np_local+1
            xvp(:,np_local)=x
        else
            np_local_dm=np_local_dm+1
            xvp_dm(:,np_local_dm)=x
        endif
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
        if (command == 0) then
            np_local=np_local+1
            xvp(:,np_local)=x
        else
            np_local_dm=np_local_dm+1
            xvp_dm(:,np_local_dm)=x
        endif
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

    if (command == 0) then
        call mpi_reduce(np_local,np_final,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)
    else
        call mpi_reduce(np_local_dm,np_final,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)
    endif

    if (rank == 0) then
      print *,'total particles =',np_final
      if (np_final /= np**3) then
        print *,'ERROR: total number of particles incorrect after passing'
      endif
    endif
 
!!  Check for particles out of bounds

    if (command == 0) then
        do i=1,np_local
          if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
              xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
              xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
            print *,'particle out of bounds',rank,i,xvp(:3,i),nc_node_dim
          endif
        enddo
    else
        do i=1,np_local_dm
          if (xvp_dm(1,i) < 0 .or. xvp_dm(1,i) >= nc_node_dim .or. &
              xvp_dm(2,i) < 0 .or. xvp_dm(2,i) >= nc_node_dim .or. &
              xvp_dm(3,i) < 0 .or. xvp_dm(3,i) >= nc_node_dim) then
            print *,'particle out of bounds',rank,i,xvp_dm(:3,i),nc_node_dim
          endif
        enddo
    endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished pass_particles ... elapsed time = ", time2-time1

    return

end subroutine pass_particles

! -------------------------------------------------------------------------------------------------------

#ifndef MOMENTUM
subroutine link_list(command)
    !
    ! Makes linked list of particles in each coarse mesh cell.
    !

    implicit none

    integer(4) :: i, j, k, pp
    integer(4) :: command
    real(8) :: sec1, sec2

    sec1 = mpi_wtime(ierr)

    !
    ! Initialize arrays to zero
    !

    if (command == 0) then
        do k = 1, max_np
            ll(k) = 0
        enddo
        do k = hoc_nc_l, hoc_nc_h
            hoc(:, :, k) = 0
        enddo
    else
        do k = 1, max_np_dm
            ll_dm(k) = 0
        enddo
        do k = hoc_nc_l, hoc_nc_h
            hoc_dm(:, :, k) = 0
        enddo
    endif

    !
    ! Now construct linked list
    !

    pp = 1
    if (command == 0) then
        do
            if (pp > np_local) exit
            i = floor(xvp(1, pp)/mesh_scale) + 1
            j = floor(xvp(2, pp)/mesh_scale) + 1
            k = floor(xvp(3, pp)/mesh_scale) + 1
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write(*,*) "WARNING: LinkList Particle Deleted: ", xvp(1:3,pp)
                xvp(:,pp) = xvp(:,np_local)
                np_local = np_local - 1
                cycle
            else
                ll(pp) = hoc(i, j, k)
                hoc(i, j, k) = pp
            endif
            pp = pp + 1
        enddo
    else
        do
            if (pp > np_local_dm) exit
            i = floor(xvp_dm(1, pp)/mesh_scale) + 1
            j = floor(xvp_dm(2, pp)/mesh_scale) + 1
            k = floor(xvp_dm(3, pp)/mesh_scale) + 1
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write(*,*) "WARNING: LinkList Particle Deleted: ", xvp_dm(1:3,pp)
                xvp_dm(:,pp) = xvp_dm(:,np_local_dm)
                np_local_dm = np_local_dm - 1
                cycle
            else
                ll_dm(pp) = hoc_dm(i, j, k)
                hoc_dm(i, j, k) = pp
            endif
            pp = pp + 1
        enddo
    endif

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished link_list ... elapsed time = ", sec2 - sec1

end subroutine link_list

! -------------------------------------------------------------------------------------------------------

subroutine buffer_particles(command)
    !
    ! Exchange coarse mesh buffers for finding halos near node edges. Add these
    ! particles to the linked list.
    !

    implicit none

    integer :: pp, i, j, k
    integer :: np_buf, nppx, npmx, nppy, npmy, nppz, npmz
    integer, dimension(mpi_status_size) :: status, sstatus, rstatus
    integer :: tag, srequest, rrequest, sierr, rierr
    real(4), parameter :: eps = 1.0e-03
    integer(4) :: np_local_start, np_buffer_sent_local
    integer(8) :: np_buffer_sent
    real(8) :: sec1, sec2
    integer(4) :: command

    sec1 = mpi_wtime(ierr)

    if (command == 0) then
        np_local_start = np_local
    else
        np_local_start = np_local_dm
    endif

    !
    ! Pass +x
    ! 

    tag = 11
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_h-hoc_pass_depth, hoc_nc_h
                if (command == 0) then
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp(1, pp) >= nc_node_dim-rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                            GID_send_buf(np_buf) = GID(pp)
                        endif
                        pp = ll(pp)
                    enddo
                else
                    pp = hoc_dm(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp_dm(1, pp) >= nc_node_dim-rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
                            GID_send_buf(np_buf) = GID_dm(pp)
                        endif
                        pp = ll_dm(pp)
                    enddo
                endif
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    nppx = np_buf

    call mpi_sendrecv_replace(nppx, 1, mpi_integer, cart_neighbor(6), tag, &
                              cart_neighbor(5), tag, mpi_comm_world, status, ierr)

    if (command == 0) then
        if (np_local+nppx > max_np) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (nppx): ", rank, np_local, nppx, max_np
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm+nppx > max_np_dm) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (nppx): ", rank, np_local_dm, nppx, max_np_dm
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, nppx*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    call mpi_isend(GID_send_buf, np_buf, mpi_integer1, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(GID_recv_buf, nppx, mpi_integer1, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    if (command == 0) then
        do i = 1, nppx
            xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            xvp(1, np_local+i) = max(xvp(1,np_local+i)-nc_node_dim, -rnf_buf)
            GID(np_local+i) = GID_recv_buf(i)
        enddo
    else
        do i = 1, nppx
            xvp_dm(:, np_local_dm+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            xvp_dm(1, np_local_dm+i) = max(xvp_dm(1,np_local_dm+i)-nc_node_dim, -rnf_buf)
            GID_dm(np_local_dm+i) = GID_recv_buf(i)
        enddo
    endif

    if (command == 0) then
        np_local = np_local + nppx
    else
        np_local_dm = np_local_dm + nppx
    endif

    !
    ! Pass -x
    !

    tag = 12
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_l, hoc_nc_l+hoc_pass_depth
                if (command == 0) then
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp(1, pp) < rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                            GID_send_buf(np_buf) = GID(pp)
                        endif
                        pp = ll(pp)
                    enddo
                else
                    pp = hoc_dm(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp_dm(1, pp) < rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
                            GID_send_buf(np_buf) = GID_dm(pp)
                        endif
                        pp = ll_dm(pp)
                    enddo
                endif
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    npmx = np_buf

    call mpi_sendrecv_replace(npmx, 1, mpi_integer, cart_neighbor(5), tag, &
                              cart_neighbor(6), tag, mpi_comm_world, status, ierr)

    if (command == 0) then
        if (np_local+npmx > max_np) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (npmx): ", rank, np_local, npmx, max_np
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm+npmx > max_np_dm) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (npmx): ", rank, np_local_dm, npmx, max_np_dm
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npmx*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    call mpi_isend(GID_send_buf, np_buf, mpi_integer1, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(GID_recv_buf, npmx, mpi_integer1, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    if (command == 0) then
        do i = 1, npmx
            xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            GID(np_local+i) = GID_recv_buf(i)
            if (abs(xvp(1, np_local+i)) .lt. eps) then
                if(xvp(1, np_local+i) < 0.) then
                    xvp(1, np_local+i) = -eps
                else
                    xvp(1, np_local+i) = eps
                endif
            endif
            xvp(1, np_local+i) = min(xvp(1, np_local+i)+real(nc_node_dim,4), &
                                     nc_node_dim+rnf_buf-eps)
        enddo
    else
        do i = 1, npmx
            xvp_dm(:, np_local_dm+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            GID_dm(np_local_dm+i) = GID_recv_buf(i)
            if (abs(xvp_dm(1, np_local_dm+i)) .lt. eps) then
                if(xvp_dm(1, np_local_dm+i) < 0.) then
                    xvp_dm(1, np_local_dm+i) = -eps
                else
                    xvp_dm(1, np_local_dm+i) = eps
                endif
            endif
            xvp_dm(1, np_local_dm+i) = min(xvp_dm(1, np_local_dm+i)+real(nc_node_dim,4), &
                                     nc_node_dim+rnf_buf-eps)
        enddo
    endif

    if (command == 0) then
        np_local = np_local + npmx
    else
        np_local_dm = np_local_dm + npmx
    endif

    !
    ! Add additional particles to the linked list
    !

    if (command == 0) then
        pp = np_local-npmx-nppx + 1
        do
            if (pp > np_local) exit
            i = floor(xvp(1, pp)/mesh_scale) + 1
            j = floor(xvp(2, pp)/mesh_scale) + 1
            k = floor(xvp(3, pp)/mesh_scale) + 1
#ifdef DIAG
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write (*, *) 'BUFFER PARTICLE DELETED', xvp(1:3, pp)
                xvp(:,pp) = xvp(:,np_local)
                np_local = np_local - 1
                cycle
            endif
#endif
            ll(pp) = hoc(i, j, k)
            hoc(i, j, k) = pp
            pp = pp + 1
        enddo
    else
        pp = np_local_dm-npmx-nppx + 1
        do
            if (pp > np_local_dm) exit
            i = floor(xvp_dm(1, pp)/mesh_scale) + 1
            j = floor(xvp_dm(2, pp)/mesh_scale) + 1
            k = floor(xvp_dm(3, pp)/mesh_scale) + 1
#ifdef DIAG
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write (*, *) 'BUFFER PARTICLE DELETED', xvp_dm(1:3, pp)
                xvp_dm(:,pp) = xvp_dm(:,np_local_dm)
                np_local_dm = np_local_dm - 1
                cycle
            endif
#endif
            ll_dm(pp) = hoc_dm(i, j, k)
            hoc_dm(i, j, k) = pp
            pp = pp + 1
        enddo
    endif

    !
    ! Pass +y
    ! 

    tag = 13
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_h-hoc_pass_depth, hoc_nc_h
            do i = hoc_nc_l, hoc_nc_h
                if (command == 0) then
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp(2, pp) >= nc_node_dim-rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                            GID_send_buf(np_buf) = GID(pp)
                        endif
                        pp = ll(pp)
                    enddo
                else
                    pp = hoc_dm(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp_dm(2, pp) >= nc_node_dim-rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
                            GID_send_buf(np_buf) = GID_dm(pp)
                        endif
                        pp = ll_dm(pp)
                    enddo
                endif
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    nppy = np_buf

    call mpi_sendrecv_replace(nppy, 1, mpi_integer, cart_neighbor(4), tag, &
                              cart_neighbor(3), tag, mpi_comm_world, status, ierr)

    if (command == 0) then
        if (np_local+nppy > max_np) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (nppy): ", rank, np_local, nppy, max_np
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm+nppy > max_np_dm) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (nppy): ", rank, np_local_dm, nppy, max_np_dm
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, nppy*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    call mpi_isend(GID_send_buf, np_buf, mpi_integer1, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(GID_recv_buf, nppy, mpi_integer1, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    if (command == 0) then
        do i = 1, nppy
            xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            xvp(2, np_local+i) = max(xvp(2,np_local+i)-nc_node_dim, -rnf_buf)
            GID(np_local+i) = GID_recv_buf(i)
        enddo
    else
        do i = 1, nppy
            xvp_dm(:, np_local_dm+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            xvp_dm(2, np_local_dm+i) = max(xvp_dm(2,np_local_dm+i)-nc_node_dim, -rnf_buf)
            GID_dm(np_local_dm+i) = GID_recv_buf(i)
        enddo
    endif

    if (command == 0) then 
        np_local = np_local + nppy
    else
        np_local_dm = np_local_dm + nppy
    endif

    !
    ! Pass -y
    !

    tag = 14
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_l+hoc_pass_depth
            do i = hoc_nc_l, hoc_nc_h
                if (command == 0) then 
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp(2, pp) < rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                            GID_send_buf(np_buf) = GID(pp)
                        endif
                        pp = ll(pp)
                    enddo
                else
                    pp = hoc_dm(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp_dm(2, pp) < rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
                            GID_send_buf(np_buf) = GID_dm(pp)
                        endif
                        pp = ll_dm(pp)
                    enddo
                endif
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    npmy = np_buf

    call mpi_sendrecv_replace(npmy, 1, mpi_integer, cart_neighbor(3), tag, &
                              cart_neighbor(4), tag, mpi_comm_world, status, ierr)

    if (command == 0) then 
        if (np_local+npmy > max_np) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (npmy): ", rank, np_local, npmy, max_np
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm+npmy > max_np_dm) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (npmy): ", rank, np_local_dm, npmy, max_np_dm
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npmy*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    call mpi_isend(GID_send_buf, np_buf, mpi_integer1, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(GID_recv_buf, npmy, mpi_integer1, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    if (command == 0) then
        do i = 1, npmy
            xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            GID(np_local+i) = GID_recv_buf(i)
            if (abs(xvp(2, np_local+i)) .lt. eps) then
                if(xvp(2, np_local+i) < 0.) then
                    xvp(2, np_local+i) = -eps
                else
                    xvp(2, np_local+i) = eps
                endif
            endif
            xvp(2, np_local+i) = min(xvp(2,np_local+i)+real(nc_node_dim,4), &
                                     nc_node_dim+rnf_buf-eps)
        enddo
    else
        do i = 1, npmy
            xvp_dm(:, np_local_dm+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            GID_dm(np_local_dm+i) = GID_recv_buf(i)
            if (abs(xvp_dm(2, np_local_dm+i)) .lt. eps) then
                if(xvp_dm(2, np_local_dm+i) < 0.) then
                    xvp_dm(2, np_local_dm+i) = -eps
                else
                    xvp_dm(2, np_local_dm+i) = eps
                endif
            endif
            xvp_dm(2, np_local_dm+i) = min(xvp_dm(2,np_local_dm+i)+real(nc_node_dim,4), &
                                     nc_node_dim+rnf_buf-eps)
        enddo
    endif

    if (command == 0) then
        np_local = np_local + npmy
    else
        np_local_dm = np_local_dm + npmy
    endif

    !
    ! Add additional particles to the linked list 
    !

    if (command == 0) then
        pp = np_local-npmy-nppy + 1
        do
            if (pp > np_local) exit
            i = floor(xvp(1, pp)/mesh_scale) + 1
            j = floor(xvp(2, pp)/mesh_scale) + 1
            k = floor(xvp(3, pp)/mesh_scale) + 1
#ifdef DIAG
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write (*, *) 'BUFFER PARTICLE DELETED', xvp(1:3, pp)
                xvp(:,pp) = xvp(:,np_local)
                np_local = np_local - 1
                cycle
            endif
#endif
            ll(pp) = hoc(i, j, k)
            hoc(i, j, k) = pp
            pp = pp + 1
        enddo
    else
        pp = np_local_dm-npmy-nppy + 1
        do
            if (pp > np_local_dm) exit
            i = floor(xvp_dm(1, pp)/mesh_scale) + 1
            j = floor(xvp_dm(2, pp)/mesh_scale) + 1
            k = floor(xvp_dm(3, pp)/mesh_scale) + 1
#ifdef DIAG
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write (*, *) 'BUFFER PARTICLE DELETED', xvp_dm(1:3, pp)
                xvp_dm(:,pp) = xvp_dm(:,np_local_dm)
                np_local_dm = np_local_dm - 1
                cycle
            endif
#endif
            ll_dm(pp) = hoc_dm(i, j, k)
            hoc_dm(i, j, k) = pp
            pp = pp + 1
        enddo
    endif

    !
    ! Pass +z
    ! 

    tag = 15
    np_buf = 0

    do k = hoc_nc_h-hoc_pass_depth, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_l, hoc_nc_h
                if (command == 0) then
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp(3, pp) >= nc_node_dim-rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                            GID_send_buf(np_buf) = GID(pp)
                        endif
                        pp = ll(pp)
                    enddo
                else
                    pp = hoc_dm(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp_dm(3, pp) >= nc_node_dim-rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
                            GID_send_buf(np_buf) = GID_dm(pp)
                        endif
                        pp = ll_dm(pp)
                    enddo
                endif
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    nppz = np_buf

    call mpi_sendrecv_replace(nppz, 1, mpi_integer, cart_neighbor(2), tag, &
                              cart_neighbor(1), tag, mpi_comm_world, status, ierr)

    if (command == 0) then
        if (np_local+nppz > max_np) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (nppz):", rank, np_local, nppz, max_np
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm+nppz > max_np_dm) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (nppz):", rank, np_local_dm, nppz, max_np_dm
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, nppz*6, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    call mpi_isend(GID_send_buf, np_buf, mpi_integer1, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(GID_recv_buf, nppz, mpi_integer1, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    if (command == 0) then
        do i = 1, nppz
            xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            xvp(3, np_local+i) = max(xvp(3,np_local+i)-nc_node_dim, -rnf_buf)
            GID(np_local+i) = GID_recv_buf(i)
        enddo
    else
        do i = 1, nppz
            xvp_dm(:, np_local_dm+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            xvp_dm(3, np_local_dm+i) = max(xvp_dm(3,np_local_dm+i)-nc_node_dim, -rnf_buf)
            GID_dm(np_local_dm+i) = GID_recv_buf(i)
        enddo
    endif

    if (command == 0) then
        np_local = np_local + nppz
    else
        np_local_dm = np_local_dm + nppz
    endif

    !
    ! Pass -z
    !

    tag = 16
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_l+hoc_pass_depth
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_l, hoc_nc_h
                if (command == 0) then
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp(3, pp) < rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                            GID_send_buf(np_buf) = GID(pp)
                        endif
                        pp = ll(pp)
                    enddo
                else
                    pp = hoc_dm(i, j, k)
                    do
                        if (pp == 0) exit
                        if (xvp_dm(3, pp) < rnf_buf) then
                            np_buf = np_buf + 1
                            send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
                            GID_send_buf(np_buf) = GID_dm(pp)
                        endif
                        pp = ll_dm(pp)
                    enddo
                endif
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    npmz = np_buf

    call mpi_sendrecv_replace(npmz, 1, mpi_integer, cart_neighbor(1), tag, &
                              cart_neighbor(2), tag, mpi_comm_world, status, ierr)

    if (command == 0) then
        if (np_local+npmz > max_np) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (npmz): ", rank, np_local, npmz, max_np
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_dm+npmz > max_np_dm) then
            write(*,*) "ERROR: Not enough space to receive buffer particles (npmz): ", rank, np_local_dm, npmz, max_np_dm
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npmz*6, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    call mpi_isend(GID_send_buf, np_buf, mpi_integer1, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(GID_recv_buf, npmz, mpi_integer1, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    if (command == 0) then
        do i = 1, npmz
            xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            GID(np_local+i) = GID_recv_buf(i)
            if (abs(xvp(3, np_local+i)) .lt. eps) then
                if(xvp(3, np_local+i) < 0.) then
                    xvp(3, np_local+i) = -eps
                else
                    xvp(3, np_local+i) = eps
                endif
            endif
            xvp(3,np_local+i) = min(xvp(3,np_local+i)+real(nc_node_dim,4), &
                                     nc_node_dim+rnf_buf-eps)
        enddo
    else
        do i = 1, npmz
            xvp_dm(:, np_local_dm+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
            GID_dm(np_local_dm+i) = GID_recv_buf(i)
            if (abs(xvp_dm(3, np_local_dm+i)) .lt. eps) then
                if(xvp_dm(3, np_local_dm+i) < 0.) then
                    xvp_dm(3, np_local_dm+i) = -eps
                else
                    xvp_dm(3, np_local_dm+i) = eps
                endif
            endif
            xvp_dm(3,np_local_dm+i) = min(xvp_dm(3,np_local_dm+i)+real(nc_node_dim,4), &
                                     nc_node_dim+rnf_buf-eps)
        enddo
    endif

    if (command == 0) then
        np_local = np_local + npmz
    else
        np_local_dm = np_local_dm + npmz
    endif

    !
    ! Add additional particles to the linked list 
    !

    if (command == 0) then
        pp = np_local-npmz-nppz + 1
        do
            if (pp > np_local) exit
            i = floor(xvp(1, pp)/mesh_scale) + 1
            j = floor(xvp(2, pp)/mesh_scale) + 1
            k = floor(xvp(3, pp)/mesh_scale) + 1
#ifdef DIAG
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write (*, *) 'BUFFER PARTICLE DELETED', xvp(1:3, pp)
                xvp(:,pp) = xvp(:,np_local)
                np_local = np_local - 1
                cycle
            endif
#endif
            ll(pp) = hoc(i, j, k)
            hoc(i, j, k) = pp
            pp = pp + 1
        enddo
    else
        pp = np_local_dm-npmz-nppz + 1
        do
            if (pp > np_local_dm) exit
            i = floor(xvp_dm(1, pp)/mesh_scale) + 1
            j = floor(xvp_dm(2, pp)/mesh_scale) + 1
            k = floor(xvp_dm(3, pp)/mesh_scale) + 1
#ifdef DIAG
            if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
                j < hoc_nc_l .or. j > hoc_nc_h .or. &
                k < hoc_nc_l .or. k > hoc_nc_h) then
                write (*, *) 'BUFFER PARTICLE DELETED', xvp_dm(1:3, pp)
                xvp_dm(:,pp) = xvp_dm(:,np_local_dm)
                np_local_dm = np_local_dm - 1
                cycle
            endif
#endif
            ll_dm(pp) = hoc_dm(i, j, k)
            hoc_dm(i, j, k) = pp
            pp = pp + 1
        enddo
    endif

    !
    ! Collect some statistics
    !

    if (command == 0) then
        np_buffer_sent_local = np_local - np_local_start
    else
        np_buffer_sent_local = np_local_dm - np_local_start
    endif
    call mpi_reduce(int(np_buffer_sent_local, kind=8), np_buffer_sent, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    if (rank == 0) write(*,*) "Total Buffer Particles: ", np_buffer_sent

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished buffer_particles ... elapsed time = ", sec2 - sec1

end subroutine buffer_particles

! -------------------------------------------------------------------------------------------------------

subroutine initialize_random_number
    !
    ! Initialize random number generator
    !

    implicit none

    integer :: i, j, k
    real(4) :: x
    integer :: seedsize, clock
    integer, dimension(8) :: values
    integer(4), allocatable, dimension(:) :: iseed
    integer(4), allocatable, dimension(:) :: iseed_all

    call random_seed()
    call random_seed(size=seedsize)

    allocate(iseed(seedsize))
    allocate(iseed_all(seedsize*nodes))

    call system_clock(count=clock)
    iseed = clock + 37 * (/ (i - 1, i = 1, 2) /)

    call date_and_time(values=values)
    if(rank==0) write(*,*) values(7:8), iseed

    call random_seed(put=values(7:8))
    call random_seed(put=iseed(1:seedsize))

    if (rank == 0) then
        write(*,*) 'Generating seeds'
        do j = 1, nodes
            do k = 1, seedsize
                call random_number(x)
                iseed_all((j-1)*seedsize+k) = int(x*huge(0))
            enddo
        enddo
    endif
    call mpi_scatter(iseed_all, seedsize, mpi_integer, iseed, seedsize, mpi_integer, 0, mpi_comm_world, ierr)

    call random_seed(put=iseed(1:seedsize))

    return

end subroutine initialize_random_number

! -------------------------------------------------------------------------------------------------------

subroutine assign_groups(command)
    !
    ! To remove shot noise in autocorrelation we randomly divide all particles
    ! into one of two groups (controlled by their value of GID). Then we
    ! cross-correlate the two groups together.
    !

    implicit none

    integer(4) :: command
    integer :: k, np_tot, np0, np1
    integer(8) :: np1_tot, np2_tot, npa_tot
    real(4) :: r 
    integer(1) :: g

    if (command == 0) then
        np_tot = np_local
    else 
        np_tot = np_local_dm
    endif

    np0 = 0
    np1 = 0

    do k = 1, np_tot

        call random_number(r)
        g = int(r+0.5)

        if (g == 0) then
            np0 = np0 + 1
        else if (g == 1) then
            np1 = np1 + 1
        endif

        if (command == 0) then
            GID(k) = g
        else 
            GID_dm(k) = g
        endif

    enddo

    call mpi_reduce(int(np_tot,kind=8), npa_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(int(np0,kind=8), np1_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(int(np1,kind=8), np2_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    if (rank == 0) then
        write(*,*) "Groups assigned: ", np1_tot, np2_tot, npa_tot
    endif

    return

end subroutine assign_groups

! -------------------------------------------------------------------------------------------------------

subroutine clear_groups
    !
    ! Assign all particles to group 0
    !

    implicit none

    GID(:) = 0
    GID_dm(:) = 0

    return

end subroutine clear_groups

! -------------------------------------------------------------------------------------------------------

subroutine velocity_density(command, glook, nfind)
    !
    ! Determine velocity field by taking the average velocity of the closest
    ! particles to each grid centre.
    ! 

    implicit none

    integer :: ind, i, j, k, pp, thread
    integer :: ic, jc, kc
    integer :: npart
    real, dimension(3) :: dr, rc
    real :: vx
    integer(4) :: command, nfind
    integer(1) :: glook
    logical :: converged
    integer(8) :: num_notconverged, num2

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize to zeros
    !

    do i = 0, nc_node_dim + 1
        momden(:, :, i) = 0.
    enddo

    !
    ! For each cell find the closest nfind particles and average their velocities
    !

    num_notconverged = 0

    !$omp  parallel num_threads(nt) default(shared) private(i, j, k, thread, rc, npart, ic, jc, kc, ind, pp, dr, vx, converged) reduction(+:num_notconverged)
    thread = 1
    thread = omp_get_thread_num() + 1
    !$omp do
    do k = 1, nc_node_dim
        do j = 1, nc_node_dim
            do i = 1, nc_node_dim

                !! Centre of cell
                rc(1) = i - 0.5
                rc(2) = j - 0.5
                rc(3) = k - 0.5

                !! Number of total particles found
                npart = 0

                !! Determine whether or not this procedure converged
                converged = .false.

                if (command == 0) then
                    do ind = 1, num_ngbhs

                        !! Cell indices to search within
                        ic = int((i-1)/mesh_scale) + 1 + cell_search(ind, 1)
                        jc = int((j-1)/mesh_scale) + 1 + cell_search(ind, 2)
                        kc = int((k-1)/mesh_scale) + 1 + cell_search(ind, 3)

                        !! Loop over particles in this cell
                        pp = hoc(ic, jc, kc)
                        do
                            if (pp == 0) exit
                            if (GID(pp) == glook) then
                                npart = npart + 1
                                if (npart > max_npart_search) then
                                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_search
                                    call mpi_abort(mpi_comm_world, ierr, ierr)
                                endif
                                dr(:) = xvp(:3,pp) - rc(:)
                                rpos(thread, 1, npart) = xvp(3+cur_dimension,pp)
                                rpos(thread, 2, npart) = dr(1)**2 + dr(2)**2 + dr(3)**2
                            endif
                            pp = ll(pp)
                        enddo !! pp

                        !! Exit if we have found at least nfind particles and we have completed an entire shell.
                        if (npart >= nfind .and. done_shell(ind)) then
                            converged = .true.
                            exit
                        endif

                    enddo !! ind
                else
                    do ind = 1, num_ngbhs

                        !! Cell indices to search within
                        ic = int((i-1)/mesh_scale) + 1 + cell_search(ind, 1)
                        jc = int((j-1)/mesh_scale) + 1 + cell_search(ind, 2)
                        kc = int((k-1)/mesh_scale) + 1 + cell_search(ind, 3)

                        !! Loop over particles in this cell
                        pp = hoc_dm(ic, jc, kc)
                        do
                            if (pp == 0) exit
                            if (GID_dm(pp) == glook) then
                                npart = npart + 1
                                if (npart > max_npart_search) then
                                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_search
                                    call mpi_abort(mpi_comm_world, ierr, ierr)
                                endif
                                dr(:) = xvp_dm(:3,pp) - rc(:)
                                rpos(thread, 1, npart) = xvp_dm(3+cur_dimension,pp)
                                rpos(thread, 2, npart) = dr(1)**2 + dr(2)**2 + dr(3)**2
                            endif
                            pp = ll_dm(pp)
                        enddo !! pp

                        !! Exit if we have found at least nfind particles and we have completed an entire shell.
                        if (npart >= nfind .and. done_shell(ind)) then
                            converged = .true.
                            exit
                        endif 

                    enddo !! ind
                endif

                !! Sort particles from closest to furthest from the cell centre 
                ipos(thread, :npart) =  (/ (ic, ic=1, npart) /)
                call indexedsort(npart, rpos(thread,2,:npart), ipos(thread,:npart))

                !! Get the average of the closest particles 
                vx = 0.
                do ind = 1, nfind 
                    vx = vx + rpos(thread, 1, ipos(thread, ind))
                enddo
                momden(i, j, k) = vx / nfind 

                !! Count number of cells that did not converge 
                if (converged .eqv. .false.) then
                    num_notconverged = num_notconverged + 1
                endif

            enddo ! i
        enddo ! j 
    enddo ! k
    !$omp end do
    !$omp end parallel

    call mpi_reduce(num_notconverged, num2, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    time2 = mpi_wtime(ierr)
    if (rank == 0) then 
        write(*, *) "Finished velocity_density ... elapsed time = ", time2-time1
        if (num2 /= 0) then
            write(*,*) "WARNING: num_notconverged = ", num2, " ... consider increasing nfine_buf !!"
            write(*,*) "         command, glook, nfind = ", command, glook, nfind 
        endif
    endif

    return

end subroutine velocity_density

#endif

! -------------------------------------------------------------------------------------------------------

#ifdef MOMENTUM
subroutine mass_density(command)
    !
    ! Bin particles in position space to generate mass density field 
    ! 

    implicit none

    real :: mp
    integer :: i, j, i1, i2, j1, j2, k1, k2
    real    :: x, y, z, dx1, dx2, dy1, dy2, dz1, dz2

    integer(4) :: command, iend

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize grid to 0
    !

    do k1 = 0, nc_node_dim+1
        massden(:,:,k1) = 0.
    enddo

    if (command == 0) then 
        mp = (ncr/np)**3
    else !! Dark matter are further reduced by factor of r_n_1_3 
        mp = (ncr/(np/ratio_nudm_dim))**3
    endif

    if (command == 0) then
        iend = np_local
    else
        iend = np_local_dm
    endif

    do i = 1, iend

        !! Read particle position
        if (command == 0) then
            x = xvp(1, i) - 0.5
            y = xvp(2, i) - 0.5
            z = xvp(3, i) - 0.5
        else
            x = xvp_dm(1, i) - 0.5
            y = xvp_dm(2, i) - 0.5
            z = xvp_dm(3, i) - 0.5
        endif

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

subroutine momentum_density(command)
    !
    ! Bin particles in position space to generate the 3D momentum density field
    ! 

    implicit none

    real :: mp
    integer :: i, j, i1, i2, j1, j2, k1, k2
    real    :: x, y, z, dx1, dx2, dy1, dy2, dz1, dz2
    real    :: dv1, dv2, v(3)
    real    :: vsim2phys, zcur
    real    :: vrmsx, vrmsy, vrmsz

    integer(4) :: command, iend

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize grid to 0
    !

    do k1 = 0, nc_node_dim+1
        momden(:,:,k1) = 0.
    enddo

    if (command == 0) then 
        mp = (ncr/np)**3
    else !! Dark matter are further reduced by factor of r_n_1_3 
        mp = (ncr/(np/ratio_nudm_dim))**3
    endif

    !
    ! Determine RMS velocity in each dimension
    !

    vrmsx = 0.
    vrmsy = 0.
    vrmsz = 0.

    if (command == 0) then
        iend = np_local
    else
        iend = np_local_dm
    endif

    do i = 1, iend

        !! Read particle position
        if (command == 0) then
            x = xvp(1, i) - 0.5 
            y = xvp(2, i) - 0.5
            z = xvp(3, i) - 0.5
            v(1) = xvp(4, i)
            v(2) = xvp(5, i)
            v(3) = xvp(6, i)
        else
            x = xvp_dm(1, i) - 0.5
            y = xvp_dm(2, i) - 0.5
            z = xvp_dm(3, i) - 0.5
            v(1) = xvp_dm(4, i)
            v(2) = xvp_dm(5, i)
            v(3) = xvp_dm(6, i)
        endif

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

        vrmsx = vrmsx + v(1)**2
        vrmsy = vrmsy + v(2)**2
        vrmsz = vrmsz + v(3)**2

        dv1 = dz1 * mp * v(cur_dimension)
        dv2 = dz2 * mp * v(cur_dimension)

        momden(i1, j1, k1) = momden(i1, j1, k1) + dx1 * dy1 * dv1
        momden(i2, j1, k1) = momden(i2, j1, k1) + dx2 * dy1 * dv1
        momden(i1, j2, k1) = momden(i1, j2, k1) + dx1 * dy2 * dv1
        momden(i2, j2, k1) = momden(i2, j2, k1) + dx2 * dy2 * dv1
        momden(i1, j1, k2) = momden(i1, j1, k2) + dx1 * dy1 * dv2
        momden(i2, j1, k2) = momden(i2, j1, k2) + dx2 * dy1 * dv2
        momden(i1, j2, k2) = momden(i1, j2, k2) + dx1 * dy2 * dv2
        momden(i2, j2, k2) = momden(i2, j2, k2) + dx2 * dy2 * dv2

    enddo

    !
    ! Determine root mean square velocities
    !

    ! First get conversion factor to physical velocity in km/s
    zcur      = z_checkpoint(cur_checkpoint)
    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc !! Takes h = 1

    vrmsx = vsim2phys * sqrt(vrmsx / iend)
    vrmsy = vsim2phys * sqrt(vrmsy / iend)
    vrmsz = vsim2phys * sqrt(vrmsz / iend)

    if (rank == 0 .and. cur_dimension == 0) write(*, *) "RMS physical velocities (rank, z, vx, vy, vz): ", rank, zcur, vrmsx, vrmsy, vrmsz

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished momentum_density ... elapsed time = ", time2-time1

    return

end subroutine momentum_density

! -------------------------------------------------------------------------------------------------------

subroutine momentum2velocity
    !
    ! Converts momentum density field v*(1+delta) to velocity field by dividing
    ! by density field (1+delta)
    !

    implicit none

    integer :: m, i, j, k
    integer(8) :: ecount_local, ecount

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

#ifdef GAUSSIAN_SMOOTH
    !
    ! First smooth the density field
    !

    cube(:, :, :) = massden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
    call Gaussian_filter
    massden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) = cube(:, :, :)

    !
    ! Now smooth each component of the momentum field
    !

    cube(:, :, :) = momden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
    call Gaussian_filter
    momden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) = cube(:, :, :)

#endif

    !
    ! Divide the latter by the former to get the velocity field
    !

    ecount_local = 0

    do i = 1, nc_node_dim
        do j = 1, nc_node_dim
            do k = 1, nc_node_dim
                if (massden(i, j, k) .ne. 0.) then
                    momden(i, j, k) = momden(i, j, k) / massden(i, j, k)
                else
                    ecount_local = ecount_local + 1
                endif
            enddo
        enddo
    enddo

    call mpi_reduce(ecount_local, ecount, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    time2 = mpi_wtime(ierr)

    if (rank == 0) write(*, *) "Empty cells in the density field = ", ecount
    if (rank == 0) write(*, *) "Finished momentum2velocity ... elapsed time = ", time2-time1

    return

end subroutine momentum2velocity

! -------------------------------------------------------------------------------------------------------

#ifdef GAUSSIAN_SMOOTH
subroutine Gaussian_filter
    !
    ! Smooths the provided real space field (in cube) with a Gaussian filter in
    ! Fourier space (via slab) and transforms back to real space (through cube).
    !

    implicit none

    integer :: i, j, k, kg, ig, mg, jg
    real :: fr, kr, kx, ky, kz
#ifndef SLAB
    integer :: ind, dx, dxy
#endif

    !! Forward transform cube to slab
    call cp_fftw(1)

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

    !! Convolve Fourier field with a Gaussian filter

#ifdef SLAB
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
#else
    do k = 1, nc_pen+mypadd
        do j = 1, nc_node_dim
            do i = 1, nc, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc
                endif
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                kx = (ig-1)/2
#endif
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
#endif
#endif

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

        !
        ! Pass +x
        ! 

        tag = 111

        momden_send_buff(:, :) = momden(nc_node_dim+1, :, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(1, :, :) = momden(1, :, :) + momden_recv_buff(:, :)

        !
        ! Pass -x
        ! 

        tag = 112

        momden_send_buff(:, :) = momden(0, :, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(nc_node_dim, :, :) = momden(nc_node_dim, :, :) + momden_recv_buff(:, :)

        !
        ! Pass +y
        ! 

        tag = 113

        momden_send_buff(:, :) = momden(:, nc_node_dim+1, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(:, 1, :) = momden(:, 1, :) + momden_recv_buff(:, :)

        !
        ! Pass -y
        ! 

        tag = 114

        momden_send_buff(:, :) = momden(:, 0, :)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(:, nc_node_dim, :) = momden(:, nc_node_dim, :) + momden_recv_buff(:, :)

        !
        ! Pass +z
        ! 

        tag = 115

        momden_send_buff(:, :) = momden(:, :, nc_node_dim+1)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(:, :, 1) = momden(:, :, 1) + momden_recv_buff(:, :)

        !
        ! Pass -z
        ! 

        tag = 116

        momden_send_buff(:, :) = momden(:, :, 0)

        call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        momden(:, :, nc_node_dim) = momden(:, :, nc_node_dim) + momden_recv_buff(:, :)

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished buffer_momdensity ... elapsed time = ", time2-time1

    return

end subroutine buffer_momdensity

! -------------------------------------------------------------------------------------------------------

#ifdef MOMENTUM
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
#endif

! -------------------------------------------------------------------------------------------------------

  subroutine powerspectrum(delta, delta2, pk)
    implicit none
    real, dimension(3, nc)       :: pk
#ifdef SLAB
    real, dimension(nc+2,nc,nc_slab) :: delta, delta2
#else
    real, dimension(nc, nc_node_dim, nc_pen+2) :: delta, delta2
#endif

    integer :: i, j, k, kg, ig, mg, jg
    integer :: k1, k2
    real    :: kr, kx, ky, kz, w1, w2, pow, x, y, z, sync_x, sync_y, sync_z, kernel
#ifdef SLAB
    real, dimension(3,nc,nc_slab) :: pkt
#else
    real, dimension(3,nc,nc_pen+2) :: pkt
#endif
    real, dimension(3, nc) :: pktsum

    real(8), dimension(nc) :: kcen, kcount
    real(8), dimension(nc) :: kcensum, kcountsum
    real    :: kavg
    integer :: ind, dx, dxy

#ifdef LOGBIN
    real :: k1r
#endif

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    pkt = 0.0
    pktsum = 0.0

    kcen(:)   = 0.
    kcount(:) = 0.
    kcensum(:) = 0.
    kcountsum(:) = 0.

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

#ifdef SLAB
    !! Compute power spectrum
    !COULD OMP DO PARALLEL THIS LOOP?
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
#else
    !! Compute power spectrum
    !! Cannot thread because of ind index
    do k = 1, nc_pen+mypadd
        do j = 1, nc_node_dim
            do i = 1, nc, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc
                endif
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                kx = (ig-1)/2
#endif
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if(kx.eq.0 .and. ky <=0 .and. kz <=0)cycle;
                if(kx.eq.0 .and. ky >0 .and. kz <0)cycle;
#ifndef LOGBIN
                if (kr .ne. 0) then
                    k1=ceiling(kr)
                    w1=k1-kr
#else
                if (kr > 1. .and. kr <= hcr) then
                    k1r = log10(kr)/log10(hcr)*numbins
                    k1  = ceiling(k1r)
                    w1 = k1-k1r
#endif
                    k2=k1+1
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
                    pow=sum((delta(i:i+1,j,k)*delta2(i:i+1,j,k)/real(ncr)**6))/kernel**4
                    pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                    pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                    pkt(3,k1,k)=pkt(3,k1,k)+w1
                    pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                    pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                    pkt(3,k2,k)=pkt(3,k2,k)+w2

#ifdef LOGBIN
                    kcen(k1) = kcen(k1) + w1 * log10(kr)
                    kcen(k2) = kcen(k2) + w2 * log10(kr)
#else
                    kcen(k1) = kcen(k1) + w1 * kr
                    kcen(k2) = kcen(k2) + w2 * kr
#endif

                    kcount(k1) = kcount(k1) + w1
                    kcount(k2) = kcount(k2) + w2

                endif
            enddo
        enddo
    enddo

    !! Merge power spectrum from threads
#ifdef SLAB
    do k=2,nc_slab
#else
    do k=2,nc_pen+2
#endif
       pkt(:,:,1)=pkt(:,:,1)+pkt(:,:,k)
    enddo

    !! Reduce to rank 0
    call mpi_reduce(pkt(:,:,1),pktsum,3*nc,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(kcen, kcensum, nc, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(kcount, kcountsum, nc, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)

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

#ifdef LOGBIN
                kavg = real(10**(kcensum(k) / kcountsum(k)), kind=4)
#else
                kavg = real(kcensum(k) / kcountsum(k), kind=4)
#endif
                pk(3,k) = 2. * pi * kavg / box

#ifdef NGP
                pk(1:2,k)=4*pi*(kavg)**3*pk(1:2,k)
#else
                pk(1:2,k)=4*pi*(kavg-1.)**3*pk(1:2,k)
#endif

            endif
        enddo
    endif

    call mpi_bcast(pk,3*nc,mpi_real,0,mpi_comm_world,ierr)

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished powerspectrum ... elapsed time = ", time2-time1

    return

  end subroutine powerspectrum

! -------------------------------------------------------------------------------------------------------

subroutine writevelocityfield(command)

    implicit none

    integer :: m, i, j, k, fstat
    character(len=180) :: fn
    character(len=7)   :: z_write
    character(len=4)   :: rank_string
    character(len=1)   :: dim_string
    real :: vsim2phys, zcur

    integer :: command

    !
    ! Determine conversion to proper velocity [km/s]
    !

    if (rank == 0)  zcur = z_checkpoint(cur_checkpoint)
    call mpi_bcast(zcur, 1, mpi_real, 0, mpi_comm_world, ierr)

    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

    if (rank == 0) write(*,*) "zcur = ", zcur, "vsim2phys = ", vsim2phys

    !
    ! Checkpoint and rank strings
    !

    write(z_write, '(f7.3)') zcur
    z_write = adjustl(z_write)

    write(rank_string, '(i4)') rank
    rank_string = adjustl(rank_string)

    !
    ! Write out velocity field for each dimension
    !

        m = cur_dimension 

        if (m == 1) dim_string = "x"
        if (m == 2) dim_string = "y"
        if (m == 3) dim_string = "z"

        if (command == 0) then
            fn = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_write(1:len_trim(z_write))//&
                 "vel"//dim_string//&
                 rank_string(1:len_trim(rank_string))//"_nu.bin"
        else
            fn = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_write(1:len_trim(z_write))//&
                 "vel"//dim_string//&
                 rank_string(1:len_trim(rank_string))//".bin"
        endif

        open(unit=11, file=fn, status="replace", iostat=fstat, access="stream") 

        do k = 1, nc_node_dim
            do j = 1, nc_node_dim

                write(11) momden(1:nc_node_dim, j, k) * vsim2phys

            enddo
        enddo

        close(11)

    return

end subroutine writevelocityfield

! -------------------------------------------------------------------------------------------------------

end program cic_velpower 
