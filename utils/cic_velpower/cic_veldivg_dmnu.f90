!!
!! cic_crossvel_dmnu.f90
!!
!! Program to compute the curl and divergence components of both the dark matter and neutrino
!! velocity fields as well as for the relative velocity field between the two species.
!!
!! * Using FFTW on the SciNet GPC compile with:
!!   mpif90 -shared-intel -fpp -g -O3 -DNGP -mt_mpi cic_crossvel_dmnu.f90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_crossvel_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Using MKL on the SciNet GPC compile with:
!!   mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -mt_mpi 
!!        cic_crossvel_dmnu.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_crossvel_dmnu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64
!!
!! * Using FFTW on the SciNet BGQ compile with:
!!   mpif90 -q64 -O3 -qhot -qarch=qp -qtune=qp -WF,-DNGP cic_crossvel_dmnu.F90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_crossvel_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Optional flags:
!!   -DNGP: Uses NGP interpolation for binning of the power spectrum. 
!!   -DSLAB: Alternatively run with FFTW slab decomposition instead of P3DFFT pencil decomposition.
!!   -DCOARSE_HACK: Enables a hack to coarsen the grid for which the particles are interpolated to.
!!                  Must be used in conjunction to a change made in the parameters file (see below).
!!   -DKAISER: Adjusts for redshift space distortions.
!!   -DDEBUG: Output useful debugging information.
!!   -DCURL: Enable the computation of the curl component of absolute fields (not relative fields). 
!!           This serves mainly as a consistency check since curl = total - divergence
!!   -Dwrite_vel: Writes velocity fields to binary files (be careful, this may be a lot!)

program cic_crossvel

  use omp_lib

  implicit none

  include 'mpif.h'
#ifndef COARSE_HACK
  include '../../parameters'
#else
  !! parameters_hack should have the following line changed from:
  !! integer(4),   parameter :: nc = (nf_tile-2*nf_buf)*tiles_node_dim*nodes_dim 
  !! to something like:
  !! integer, parameter :: coarsen_factor = 4
  !! integer(4),   parameter :: nc = (nf_tile-2*nf_buf)*tiles_node_dim*nodes_dim / coarsen_factor
  include '../../parameters_hack'
#endif

  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints_nu'
  logical, parameter :: turn_off_halos = .false.  !! Can set this True if halo portion is slow
    
  !! Threading
  integer(4), parameter :: nt = 8

  !! Number of nearest particles from grid centre to determine average velocity
  integer, parameter :: N_closest_nu = 1
  integer, parameter :: N_closest_dm = 1
  integer, parameter :: N_closest_h  = 1
  integer, parameter :: max_N_closest = max(N_closest_nu, N_closest_dm, N_closest_h) + 1 

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
#ifndef COARSE_HACK
  integer, parameter :: np=hc
#else
  integer, parameter :: np=hc*coarsen_factor
#endif
  real, parameter    :: npr=np

#ifdef LOGBIN
  integer, parameter :: numbins = 32 
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
  integer(4), parameter :: max_np_h = 1000000

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local,np_local_dm,np_local_h,np_add,np_add_dm,np_add_h
  integer(8) :: plan, iplan, np_h
  logical :: firstfftw

  !! have velocity power spectra for each x, y, z
  integer cur_dimension

  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(6,max_np_dm) :: xvp_dm
  real, dimension(6,max_np_h)  :: xvmp_h
  real, dimension(6,np_buffer) :: xp_buf
  real, dimension(6*np_buffer) :: send_buf, recv_buf

  !! Power spectrum arrays
  real, dimension(3, 3, nc, 6) :: pkvel, pkcurl
  real, dimension(3, nc, 6) :: pkdivg
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
  real, dimension(nc+2,nc,nc_slab) :: slab, slab2, slab_work
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab, slab2
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1) :: recv_cube
#endif

  !! Two velocity fields with buffers
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: velden, velden2
#ifdef CURL
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: velden3
#endif
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: velden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: velden_recv_buff

  !! Curl and divergence arrays (work with one component of curl at a time)
  real, dimension(nc_node_dim, nc_node_dim, nc_node_dim) :: veldivg

  !! Parameters for linked list
  integer(4), parameter :: nfine_buf = 16
  integer(4), parameter :: mesh_scale = 4
  integer(4), parameter :: nc_buf = nfine_buf / mesh_scale
  integer(4), parameter :: nm_node_dim = nc_node_dim / mesh_scale
  integer(4), parameter :: hoc_nc_l = 1 - nc_buf
  integer(4), parameter :: hoc_nc_h = nm_node_dim + nc_buf
  integer(4), parameter :: hoc_pass_depth = 2*nc_buf
  real(4), parameter    :: rnf_buf = real(nfine_buf)
  integer(4), parameter :: num_ngbhs = (2*nc_buf+1)**3
  integer(4), parameter :: nfine_buf_h = 128
  integer(4), parameter :: mesh_scale_h = 8
  integer(4), parameter :: nc_buf_h = nfine_buf_h / mesh_scale_h
  integer(4), parameter :: nm_node_dim_h = nc_node_dim / mesh_scale_h
  integer(4), parameter :: hoc_nc_l_h = 1 - nc_buf_h
  integer(4), parameter :: hoc_nc_h_h = nm_node_dim_h + nc_buf_h
  integer(4), parameter :: hoc_pass_depth_h = 2*nc_buf_h
  real(4), parameter    :: rnf_buf_h = real(nfine_buf_h)
  integer(4), parameter :: num_ngbhs_h = (2*nc_buf_h+1)**3

  integer(4) :: hoc(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h)
  integer(4) :: hoc_dm(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h)
  integer(4) :: hoc_h(hoc_nc_l_h:hoc_nc_h_h, hoc_nc_l_h:hoc_nc_h_h, hoc_nc_l_h:hoc_nc_h_h)
  real(4)    :: rpos(nt, 2, max_N_closest)
  real(4)    :: rpos_h(nt, 2, max_N_closest)
  integer(4) :: cell_search(3, num_ngbhs)
  real(4)    :: cell_search_r(num_ngbhs)
  integer(4) :: cell_search_h(3, num_ngbhs_h)
  real(4)    :: cell_search_r_h(num_ngbhs_h)

  integer(4) :: np_groups(0:1), np_groups_dm(0:1), np_groups_h(0:1)
  integer(4) :: np_buf_groups(0:1), np_buf_groups_dm(0:1), np_buf_groups_h(0:1)
  integer(4) :: np_groups_tot(0:1), np_groups_tot_dm(0:1), np_groups_tot_h(0:1)

  integer(1), parameter :: g0 = 0
  integer(1), parameter :: g1 = 1

#ifdef CURL
  integer(4), dimension(3, 2) :: curlcom 
#endif

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
  equivalence(recv_cube, xp_buf) !! May sometimes need to be changed (see check_equivalence_divgdmnu.py)
  equivalence(velden, send_buf) !! May sometimes need to be changed (see check_equivalence_divgdmnu.py)
  equivalence(velden2, recv_buf) !! May sometimes need to be changed (see check_equivalence_divgdmnu.py)
  equivalence(hoc, hoc_dm, hoc_h)
#ifdef CURL
  equivalence(velden3, slab, cube) !! May sometimes need to be changed (see check_equivalence_divgdmnu.py) 
#else
  equivalence(slab, cube) !! Slab will always be larger than cube
#endif

  !! Common block
#ifdef SLAB
  common xvp, xvp_dm, xvmp_h, recv_cube, velden, velden2, hoc, veldivg, slab2, velden_send_buff, velden_recv_buff, slab_work
#else
  common xvp, xvp_dm, xvmp_h, recv_cube, velden, velden2, hoc, veldivg, slab2, velden_send_buff, velden_recv_buff
#endif
  common /pvar/ pkdm, pkdivg, pkcurl, pkvel 
#ifdef CURL
  common /cvar/ velden3
#else
  common /cvar/ slab
#endif

! -------------------------------------------------------------------------------------------------------
! MAIN
! -------------------------------------------------------------------------------------------------------

  call mpi_initialize
  call omp_set_num_threads(nt)
  call initialize_random_number

  if (rank == 0) call writeparams

  firstfftw = .true.  ! initialize fftw so that it generates the plans

  call read_checkpoint_list

  do cur_checkpoint = 1, num_checkpoints

    ! --------------------------------------------------------------------------------
    ! Initialize some variables
    ! --------------------------------------------------------------------------------

    call initvar
    call init_cell_search(0)
    call init_cell_search(2)

    ! --------------------------------------------------------------------------------
    ! Read, pass, and sort neutrino particles into two separate groups
    ! --------------------------------------------------------------------------------

    if (rank == 0) write(*,*) "Reading, passing, and sorting neutrino particles ... "
    call read_particles(0)
    call pass_particles(0)
    call order_xvp_groups(0)
    call buffer_particles_groups(0, g0)
    call buffer_particles_groups(0, g1)
    call order_xvp_ll(0, g0)
    call order_xvp_ll(0, g1)
    if (rank == 0) write(*,*)
       
    ! --------------------------------------------------------------------------------
    ! Read, pass, and sort dark matter particles into two separate groups
    ! --------------------------------------------------------------------------------

    if (rank == 0) write(*,*) "Reading, passing, and sorting dark matter particles ... "
    call read_particles(1)
    call pass_particles(1)
    call order_xvp_groups(1)
    call buffer_particles_groups(1, g0)
    call buffer_particles_groups(1, g1)
    call order_xvp_ll(1, g0)
    call order_xvp_ll(1, g1)
    if (rank == 0) write(*,*)

    ! --------------------------------------------------------------------------------
    ! Read, pass, and sort halo particles into two separate groups
    ! --------------------------------------------------------------------------------

    if (.not. turn_off_halos) then
        if (rank == 0) write(*,*) "Reading, passing, and sorting halo particles ... "
        call read_particles(2)
        call pass_particles(2)
        call order_xvp_groups(2)
        call buffer_particles_groups(2, g0)
        call buffer_particles_groups(2, g1)
        call order_xvp_ll(2, g0)
        call order_xvp_ll(2, g1)
        if (rank == 0) write(*,*)
    endif

    ! --------------------------------------------------------------------------------
    ! Compute total velocity power spectra for each of the three different species 
    ! (0 for neutrinos, 1 for dark matter, 2 for halos) as well as the total relative
    ! velocity power spectra. In each case, each species is divided into two random
    ! groups of particles in order to remove shot noise, which is uncorrelated between 
    ! the two groups.
    ! --------------------------------------------------------------------------------
    
    do cur_dimension = 1, 3 !! Each x, y, z dimension
    
        ! --------------------------------------------------------------------------------                                                            
        ! NEUTRINOS                                                                                                                                   
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing neutrino total velocity for dim = ", cur_dimension
        call velocity_density(cur_dimension, 0, g0, N_closest_nu)
#ifdef write_vel
        call writevelocityfield(0, g0)
#endif
        call darkmatter(0)
        call swap_slab12(0)
        call velocity_density(cur_dimension, 0, g1, N_closest_nu)
#ifdef write_vel
        call writevelocityfield(0, g1)
#endif
        call darkmatter(0)
        call powerspectrum(slab, slab2, pkvel(cur_dimension,:,:,1), 0)
        if (rank == 0) write(*,*)
        
        ! --------------------------------------------------------------------------------                                                            
        ! DARK MATTER
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing dark matter total velocity for dim = ", cur_dimension
        call velocity_density(cur_dimension, 1, g0, N_closest_dm)
#ifdef write_vel
        call writevelocityfield(1, g0)
#endif
        call darkmatter(0)
        call swap_slab12(0)
        call velocity_density(cur_dimension, 1, g1, N_closest_dm)
#ifdef write_vel
        call writevelocityfield(1, g1)
#endif
        call darkmatter(0)
        call powerspectrum(slab, slab2, pkvel(cur_dimension,:,:,2), 0)
        if (rank == 0) write(*,*)

        ! --------------------------------------------------------------------------------                                                            
        ! HALOS
        ! --------------------------------------------------------------------------------

        if (.not. turn_off_halos) then
            if (rank == 0) write(*,*) "Computing halo total velocity for dim = ", cur_dimension
            call velocity_density(cur_dimension, 2, g0, N_closest_h)
#ifdef write_vel
            call writevelocityfield(2, g0)
#endif
            call darkmatter(0)
            call swap_slab12(0)
            call velocity_density(cur_dimension, 2, g1, N_closest_h)
#ifdef write_vel
            call writevelocityfield(2, g1)
#endif
            call darkmatter(0)
            call powerspectrum(slab, slab2, pkvel(cur_dimension,:,:,3), 0)
            if (rank == 0) write(*,*)
        endif

        ! --------------------------------------------------------------------------------                                                            
        ! NEUTRINO-DARK MATTER RELATIVE FIELD
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing nu-dm total velocity for dim = ", cur_dimension
        call velocity_density(cur_dimension, 0, g0, N_closest_nu)
        call swap_velden12(0)
        call velocity_density(cur_dimension, 1, g0, N_closest_dm)
        call relative_velocity
#ifdef write_vel
        call writevelocityfield(3, g0)
#endif
        call darkmatter(0)
        call swap_slab12(0)
        call velocity_density(cur_dimension, 0, g1, N_closest_nu)
        call swap_velden12(0)
        call velocity_density(cur_dimension, 1, g1, N_closest_dm)
        call relative_velocity
#ifdef write_vel
        call writevelocityfield(3, g1)
#endif
        call darkmatter(0)
        call powerspectrum(slab, slab2, pkvel(cur_dimension,:,:,4), 0)
        if (rank == 0) write(*,*)

        if (.not. turn_off_halos) then
            ! --------------------------------------------------------------------------------                                                            
            ! NEUTRINO-HALO RELATIVE FIELD
            ! --------------------------------------------------------------------------------

            if (rank == 0) write(*,*) "Computing nu-halo total velocity for dim = ", cur_dimension
            call velocity_density(cur_dimension, 0, g0, N_closest_nu)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g0, N_closest_h)
            call relative_velocity
#ifdef write_vel
            call writevelocityfield(4, g0)
#endif
            call darkmatter(0)
            call swap_slab12(0)
            call velocity_density(cur_dimension, 0, g1, N_closest_nu)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g1, N_closest_h)
            call relative_velocity
#ifdef write_vel
            call writevelocityfield(4, g1)
#endif
            call darkmatter(0)
            call powerspectrum(slab, slab2, pkvel(cur_dimension,:,:,5), 0)
            if (rank == 0) write(*,*)

            ! --------------------------------------------------------------------------------                                                            
            ! DARK MATTER-HALO RELATIVE FIELD
            ! --------------------------------------------------------------------------------

            if (rank == 0) write(*,*) "Computing dm-halo total velocity for dim = ", cur_dimension
            call velocity_density(cur_dimension, 1, g0, N_closest_dm)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g0, N_closest_h)
            call relative_velocity
#ifdef write_vel
            call writevelocityfield(5, g0)
#endif
            call darkmatter(0)
            call swap_slab12(0)
            call velocity_density(cur_dimension, 1, g1, N_closest_dm)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g1, N_closest_h)
            call relative_velocity
#ifdef write_vel
            call writevelocityfield(5, g1)
#endif
            call darkmatter(0)
            call powerspectrum(slab, slab2, pkvel(cur_dimension,:,:,6), 0)
            if (rank == 0) write(*,*)
        endif

    enddo

    ! --------------------------------------------------------------------------------
    ! Write the total velocity power spectra for each field
    ! --------------------------------------------------------------------------------

    if (rank == 0) then
        if (.not. turn_off_halos) then
            do kt = 1, 6 
                call writepowerspectra(0, kt)
            enddo
        else
            call writepowerspectra(0, 1)
            call writepowerspectra(0, 2)
            call writepowerspectra(0, 4)
        endif
    endif

    ! --------------------------------------------------------------------------------
    ! Compute the power spectrum of the divergence component of velocity for each of 
    ! the three fields as well as the various relative velocity fields. Again we use
    ! the random group method to remove shot noise.
    ! --------------------------------------------------------------------------------

    ! --------------------------------------------------------------------------------                                                            
    ! NEUTRINOS                                                                                                                                   
    ! --------------------------------------------------------------------------------

    if (rank == 0) write(*,*) "Computing neutrino velocity divergence" 
    do cur_dimension = 1, 3 
        call velocity_density(cur_dimension, 0, g0, N_closest_nu)
        call pass_veldensity
        call velocity_divergence
    enddo
    call darkmatter(1)
    call swap_slab12(0)
    do cur_dimension = 1, 3
        call velocity_density(cur_dimension, 0, g1, N_closest_nu)
        call pass_veldensity
        call velocity_divergence
    enddo
    call darkmatter(1)
    call powerspectrum(slab, slab2, pkdivg(:,:,1), 1) 
    if (rank == 0) call writepowerspectra(1, 1)
    if (rank == 0) write(*,*)

    ! --------------------------------------------------------------------------------                                                            
    ! DARK MATTER
    ! --------------------------------------------------------------------------------

    if (rank == 0) write(*,*) "Computing dark matter velocity divergence"
    do cur_dimension = 1, 3
        call velocity_density(cur_dimension, 1, g0, N_closest_dm)
        call pass_veldensity
        call velocity_divergence
    enddo
    call darkmatter(1)
    call swap_slab12(0)
    do cur_dimension = 1, 3
        call velocity_density(cur_dimension, 1, g1, N_closest_dm)
        call pass_veldensity
        call velocity_divergence
    enddo
    call darkmatter(1)
    call powerspectrum(slab, slab2, pkdivg(:,:,2), 1)
    if (rank == 0) call writepowerspectra(1, 2)
    if (rank == 0) write(*,*)

    ! --------------------------------------------------------------------------------                                                            
    ! HALOS
    ! --------------------------------------------------------------------------------

    if (.not. turn_off_halos) then
        if (rank == 0) write(*,*) "Computing halo velocity divergence"
        do cur_dimension = 1, 3
            call velocity_density(cur_dimension, 2, g0, N_closest_h)
            call pass_veldensity
            call velocity_divergence
        enddo
        call darkmatter(1)
        call swap_slab12(0)
        do cur_dimension = 1, 3
            call velocity_density(cur_dimension, 2, g1, N_closest_h)
            call pass_veldensity
            call velocity_divergence
        enddo
        call darkmatter(1)
        call powerspectrum(slab, slab2, pkdivg(:,:,3), 1)
        if (rank == 0) call writepowerspectra(1, 3)
        if (rank == 0) write(*,*)
    endif

    ! --------------------------------------------------------------------------------                                                            
    ! NEUTRINO-DARK MATTER RELATIVE FIELD
    ! --------------------------------------------------------------------------------

    if (rank == 0) write(*,*) "Computing nu-dm velocity divergence"
    do cur_dimension = 1, 3
        call velocity_density(cur_dimension, 0, g0, N_closest_nu)
        call swap_velden12(0)
        call velocity_density(cur_dimension, 1, g0, N_closest_dm)
        call relative_velocity 
        call pass_veldensity
        call velocity_divergence
    enddo
    call darkmatter(1)
    call swap_slab12(0)
    do cur_dimension = 1, 3
        call velocity_density(cur_dimension, 0, g1, N_closest_nu)
        call swap_velden12(0)
        call velocity_density(cur_dimension, 1, g1, N_closest_dm)
        call relative_velocity
        call pass_veldensity
        call velocity_divergence
    enddo
    call darkmatter(1)
    call powerspectrum(slab, slab2, pkdivg(:,:,4), 1)
    if (rank == 0) call writepowerspectra(1, 4)
    if (rank == 0) write(*,*)

    if (.not. turn_off_halos) then
        ! --------------------------------------------------------------------------------                                                            
        ! NEUTRINO-HALO RELATIVE FIELD
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing nu-halo velocity divergence"
        do cur_dimension = 1, 3
            call velocity_density(cur_dimension, 0, g0, N_closest_nu)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g0, N_closest_h)
            call relative_velocity
            call pass_veldensity
            call velocity_divergence
        enddo
        call darkmatter(1)
        call swap_slab12(0)
        do cur_dimension = 1, 3
            call velocity_density(cur_dimension, 0, g1, N_closest_nu)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g1, N_closest_h)
            call relative_velocity
            call pass_veldensity
            call velocity_divergence
        enddo
        call darkmatter(1)
        call powerspectrum(slab, slab2, pkdivg(:,:,5), 1)
        if (rank == 0) call writepowerspectra(1, 5)
        if (rank == 0) write(*,*)

        ! --------------------------------------------------------------------------------                                                            
        ! DARK MATTER-HALO RELATIVE FIELD
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing dm-halo velocity divergence"
        do cur_dimension = 1, 3
            call velocity_density(cur_dimension, 1, g0, N_closest_dm)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g0, N_closest_h)
            call relative_velocity
            call pass_veldensity
            call velocity_divergence
        enddo
        call darkmatter(1)
        call swap_slab12(0)
        do cur_dimension = 1, 3
            call velocity_density(cur_dimension, 1, g1, N_closest_dm)
            call swap_velden12(0)
            call velocity_density(cur_dimension, 2, g1, N_closest_h)
            call relative_velocity
            call pass_veldensity
            call velocity_divergence
        enddo
        call darkmatter(1)
        call powerspectrum(slab, slab2, pkdivg(:,:,6), 1)
        if (rank == 0) call writepowerspectra(1, 6)
        if (rank == 0) write(*,*)
    endif

#ifdef CURL
    ! --------------------------------------------------------------------------------
    ! Compute the power spectrum of the curl component of velocity for each of the 
    ! three fields as well as the various relative velocity fields. Again we use the
    ! random group method to remove shot noise.
    ! --------------------------------------------------------------------------------

    do cur_dimension = 1, 3

        ! --------------------------------------------------------------------------------                                                            
        ! NEUTRINOS                                                                                                                                   
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing neutrino curl velocity for dim = ", cur_dimension
        call velocity_density(curlcom(cur_dimension,1), 0, g0, N_closest_nu)
        call pass_veldensity
        call swap_velden13(0) 
        call velocity_density(curlcom(cur_dimension,2), 0, g0, N_closest_nu)
        call pass_veldensity
        call velocity_curl
        call darkmatter(1)
        call swap_slab12(0)
        call velocity_density(curlcom(cur_dimension,1), 0, g1, N_closest_nu)
        call pass_veldensity
        call swap_velden13(0)
        call velocity_density(curlcom(cur_dimension,2), 0, g1, N_closest_nu)
        call pass_veldensity
        call velocity_curl
        call darkmatter(1)
        call powerspectrum(slab, slab2, pkcurl(cur_dimension,:,:,1), 1)
        if (rank == 0) write(*,*) 

        ! --------------------------------------------------------------------------------                                                            
        ! DARK MATTER
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing dark matter curl velocity for dim = ", cur_dimension
        call velocity_density(curlcom(cur_dimension,1), 1, g0, N_closest_dm)
        call pass_veldensity
        call swap_velden13(0)
        call velocity_density(curlcom(cur_dimension,2), 1, g0, N_closest_dm)
        call pass_veldensity
        call velocity_curl
        call darkmatter(1)
        call swap_slab12(0)
        call velocity_density(curlcom(cur_dimension,1), 1, g1, N_closest_dm)
        call pass_veldensity
        call swap_velden13(0)
        call velocity_density(curlcom(cur_dimension,2), 1, g1, N_closest_dm)
        call pass_veldensity
        call velocity_curl
        call darkmatter(1)
        call powerspectrum(slab, slab2, pkcurl(cur_dimension,:,:,2), 1)
        if (rank == 0) write(*,*)

        ! --------------------------------------------------------------------------------                                                            
        ! HALOS
        ! --------------------------------------------------------------------------------

        if (.not. turn_off_halos) then
            if (rank == 0) write(*,*) "Computing halo curl velocity for dim = ", cur_dimension
            call velocity_density(curlcom(cur_dimension,1), 2, g0, N_closest_h)
            call pass_veldensity
            call swap_velden13(0)
            call velocity_density(curlcom(cur_dimension,2), 2, g0, N_closest_h)
            call pass_veldensity
            call velocity_curl
            call darkmatter(1)
            call swap_slab12(0)
            call velocity_density(curlcom(cur_dimension,1), 2, g1, N_closest_h)
            call pass_veldensity
            call swap_velden13(0)
            call velocity_density(curlcom(cur_dimension,2), 2, g1, N_closest_h)
            call pass_veldensity
            call velocity_curl
            call darkmatter(1)
            call powerspectrum(slab, slab2, pkcurl(cur_dimension,:,:,3), 1)
            if (rank == 0) write(*,*)
        endif

        ! --------------------------------------------------------------------------------                                                            
        ! NEUTRINO-DARK MATTER RELATIVE FIELD
        ! --------------------------------------------------------------------------------

        if (rank == 0) write(*,*) "Computing nu-dm curl velocity for dim = ", cur_dimension
        call velocity_density(curlcom(cur_dimension,1), 0, g0, N_closest_nu)
        call pass_veldensity
        call swap_velden12(0)
        call velocity_density(curlcom(cur_dimension,1), 1, g0, N_closest_dm)
        call pass_veldensity
        call relative_velocity
        call swap_velden13(0)
        call velocity_density(curlcom(cur_dimension,2), 0, g0, N_closest_nu)
        call pass_veldensity
        call swap_velden12(0)
        call velocity_density(curlcom(cur_dimension,2), 1, g0, N_closest_dm)
        call pass_veldensity
        call relative_velocity
        call velocity_curl
        call darkmatter(1)
        call swap_slab12(0)
        call velocity_density(curlcom(cur_dimension,1), 0, g1, N_closest_nu)
        call pass_veldensity
        call swap_velden12(0)
        call velocity_density(curlcom(cur_dimension,1), 1, g1, N_closest_dm)
        call pass_veldensity
        call relative_velocity
        call swap_velden13(0)
        call velocity_density(curlcom(cur_dimension,2), 0, g1, N_closest_nu)
        call pass_veldensity
        call swap_velden12(0)
        call velocity_density(curlcom(cur_dimension,2), 1, g1, N_closest_dm)
        call pass_veldensity
        call relative_velocity
        call velocity_curl
        call darkmatter(1)
        call powerspectrum(slab, slab2, pkcurl(cur_dimension,:,:,4), 1)
        if (rank == 0) write(*,*)

        if (.not. turn_off_halos) then
            ! --------------------------------------------------------------------------------                                                            
            ! NEUTRINO-HALO RELATIVE FIELD
            ! --------------------------------------------------------------------------------

            if (rank == 0) write(*,*) "Computing nu-halo curl velocity for dim = ", cur_dimension
            call velocity_density(curlcom(cur_dimension,1), 0, g0, N_closest_nu)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,1), 2, g0, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call swap_velden13(0)
            call velocity_density(curlcom(cur_dimension,2), 0, g0, N_closest_nu)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,2), 2, g0, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call velocity_curl
            call darkmatter(1)
            call swap_slab12(0)
            call velocity_density(curlcom(cur_dimension,1), 0, g1, N_closest_nu)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,1), 2, g1, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call swap_velden13(0)
            call velocity_density(curlcom(cur_dimension,2), 0, g1, N_closest_nu)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,2), 2, g1, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call velocity_curl
            call darkmatter(1)
            call powerspectrum(slab, slab2, pkcurl(cur_dimension,:,:,5), 1)
            if (rank == 0) write(*,*)

            ! --------------------------------------------------------------------------------                                                            
            ! DARK MATTER-HALO RELATIVE FIELD
            ! --------------------------------------------------------------------------------

            if (rank == 0) write(*,*) "Computing dm-halo curl velocity for dim = ", cur_dimension
            call velocity_density(curlcom(cur_dimension,1), 1, g0, N_closest_dm)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,1), 2, g0, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call swap_velden13(0)
            call velocity_density(curlcom(cur_dimension,2), 1, g0, N_closest_dm)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,2), 2, g0, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call velocity_curl
            call darkmatter(1)
            call swap_slab12(0)
            call velocity_density(curlcom(cur_dimension,1), 1, g1, N_closest_dm)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,1), 2, g1, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call swap_velden13(0)
            call velocity_density(curlcom(cur_dimension,2), 1, g1, N_closest_dm)
            call pass_veldensity
            call swap_velden12(0)
            call velocity_density(curlcom(cur_dimension,2), 2, g1, N_closest_h)
            call pass_veldensity
            call relative_velocity
            call velocity_curl
            call darkmatter(1)
            call powerspectrum(slab, slab2, pkcurl(cur_dimension,:,:,6), 1)
            if (rank == 0) write(*,*)
        endif

    enddo

    ! --------------------------------------------------------------------------------
    ! Write the total velocity power spectra for each field
    ! --------------------------------------------------------------------------------

    if (rank == 0) then 
        if (.not. turn_off_halos) then
            do kt = 1, 6
                call writepowerspectra(2, kt)
            enddo
        else
            call writepowerspectra(2, 1)
            call writepowerspectra(2, 2)
            call writepowerspectra(2, 4)
        endif
    endif
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
    do k = 1, max_np_dm
        xvp_dm(:, k) = 0.
    enddo
    do k = 1, max_np_h
        xvmp_h(:, k) = 0.
    enddo

    !! Velocity density arrays
    do k = 0, nc_node_dim + 1
        velden(:, :, k) = 0.
        velden2(:, :, k) = 0.
#ifdef CURL
        velden3(:, :, k) = 0.
#endif
        velden_send_buff(:, k) = 0.
        velden_recv_buff(:, k) = 0.
    enddo
    do k = 1, nc_node_dim
        veldivg(:, :, k) = 0.
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
       slab(:, :, k) = 0.
       slab2(:, :, k) = 0.
    enddo

    !! Buffer arrays
    do k = 1, np_buffer
       xp_buf(:, k) = 0
    enddo
    do k = 1, 6*np_buffer
        send_buf(k) = 0.
        recv_buf(k) = 0.
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
        pkvel(:, :, :, k) = 0.
        pkdivg(:, :, k) = 0.
        pkcurl(:, :, :, k) = 0.
        pkdm(:, k) = 0.
    enddo

#ifdef CURL
    curlcom(1,1) = 3
    curlcom(1,2) = 2
    curlcom(2,1) = 1
    curlcom(2,2) = 3
    curlcom(3,1) = 2
    curlcom(3,2) = 1
#endif

    return

end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine init_cell_search(command)
  !
  ! Initialize cell_search (for command == 0/1) or cell_search_h (for command == 2)
  ! and also cell_search_r or cell_search_r_h
  !
  ! cell_search contains the index of the neighbour cells to one cell
  ! these indexes are relative to this cell : (0,0,0) represents this cell
  ! they are ordered by their minimal distance to the cell (stored in cell_search_r), closest first
  !

  implicit none
  integer command

  integer :: k, j, i
  integer :: k1,j1,i1,kc,jc,ic
  real    :: dmin,dmax
  integer :: ind, ibuf
  integer, parameter :: max_num_ngbhs = max(num_ngbhs, num_ngbhs_h)

  integer(4) :: cell_search_work(3, max_num_ngbhs)
  integer(4) :: sorted_indexes(num_ngbhs_h)

  if (command == 0 .or. command == 1) then
     ind = 1
     do ibuf = 0, nc_buf
        do k = -ibuf, ibuf
           do j = -ibuf, ibuf
              do i = -ibuf, ibuf
                 if (.not. (abs(i) < ibuf .and. abs(j) < ibuf .and. abs(k) < ibuf)) then
                    cell_search_work(1, ind) = i
                    cell_search_work(2, ind) = j
                    cell_search_work(3, ind) = k
                    dmin = sqrt(real(i)**2+ real(j)**2+real(k)**2)
                    do k1 = -1,1
                       do j1 = -1,1
                          do i1 = -1,1
                             do kc = -1,1
                                do jc = -1,1
                                   do ic = -1,1
                                      dmin = min(dmin, sqrt(real(i+0.5*(i1-ic))**2 + real(j+0.5*(j1-jc))**2 + real(k+0.5*(k1-kc))**2))
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                    cell_search_r(ind) = dmin * mesh_scale
                    sorted_indexes(ind) = ind
                    ind = ind + 1
                 endif
              enddo
           enddo
        enddo
     enddo

     cell_search_r(0) = 0.
     call indexedsort(num_ngbhs, cell_search_r(:num_ngbhs), sorted_indexes(:num_ngbhs))

     do ind = 1, num_ngbhs
        cell_search(:,ind) = cell_search_work(:,sorted_indexes(ind))
     enddo

  else if (command == 2) then

     ind = 1
     do ibuf = 0, nc_buf_h
        do k = -ibuf, ibuf
           do j = -ibuf, ibuf
              do i = -ibuf, ibuf
                 if (.not. (abs(i) < ibuf .and. abs(j) < ibuf .and. abs(k) < ibuf)) then
                    cell_search_work(1, ind) = i
                    cell_search_work(2, ind) = j
                    cell_search_work(3, ind) = k
                    dmin = sqrt(real(i)**2+ real(j)**2+real(k)**2)
                    do k1 = -1,1
                       do j1 = -1,1
                          do i1 = -1,1
                             do kc = -1,1
                                do jc = -1,1
                                   do ic = -1,1
                                      dmin = min(dmin, sqrt(real(i+0.5*(i1-ic))**2 + real(j+0.5*(j1-jc))**2 + real(k+0.5*(k1-kc))**2))
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                    cell_search_r_h(ind) = dmin * mesh_scale_h
                    sorted_indexes(ind) = ind
                    ind = ind + 1
                 endif
              enddo
           enddo
        enddo
     enddo

     cell_search_r_h(0) = 0.0
     call indexedsort(num_ngbhs_h, cell_search_r_h(:num_ngbhs_h), sorted_indexes(:num_ngbhs))

     do ind = 1, num_ngbhs_h
        cell_search_h(:,ind) = cell_search_work(:,sorted_indexes(ind))
    enddo
 endif

end subroutine init_cell_search

! -------------------------------------------------------------------------------------------------------

subroutine read_particles(command)
    !
    ! Read x, y, z positions and velocities and store in xvp
    !

    implicit none
    
    real z_write
    integer(8) :: np_total
    integer j, fstat
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name
    integer(4) :: command

    !! These are unnecessary headers from the checkpoint
    real(4) :: a, t, tau, dt_f_acc, dt_c_acc, dt_pp_acc, mass_p
    integer(4) :: nts, sim_checkpoint, sim_projection, sim_halofind

    !! unnecessary data for halos
    real(4)                :: garbage1
    real(4), dimension(3)  :: garbage3
    real(4), dimension(4)  :: garbage4
    real(4), dimension(15) :: garbage15

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
           check_name=ic_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        else
           check_name=output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        endif
    else if (command == 1) then
        if(z_write .eq. z_i) then
           check_name=ic_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        else
           check_name=output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        endif
    else if (command == 2) then
        check_name=output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_string(1:len_trim(z_string))//'halo'//&
                   rank_string(1:len_trim(rank_string))//'.dat'
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
    else if (command == 1) then
        read(21) np_local_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
                   sim_projection,sim_halofind,mass_p
    else
        read(21) np_local_h,t, tau
    endif

    !! Check for memory problems
    if (command == 0) then
        if (np_local > max_np) then
          write(*,*) 'ERROR: Too many particles to store in memory!'
          write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
          call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else if (command == 1) then
        if (np_local_dm > max_np_dm) then
          write(*,*) 'ERROR: Too many particles to store in memory!'
          write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
          call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    else
        if (np_local_h > max_np_h) then
          write(*,*) 'ERROR: Too many halos to store in memory!'
          write(*,*) 'rank', rank, 'np_local_h', np_local_h, 'max_np_h', max_np_h
          call mpi_abort(mpi_comm_world, ierr, ierr)
        endif
    endif

    !! Tally up total number of particles
    if (command == 0) then
        call mpi_reduce(int(np_local, kind=8), np_total, 1, mpi_integer8, &
                         mpi_sum, 0, mpi_comm_world, ierr)
    else if (command == 1) then
        call mpi_reduce(int(np_local_dm, kind=8), np_total, 1, mpi_integer8, &
                         mpi_sum, 0, mpi_comm_world, ierr)
    else
        call mpi_reduce(int(np_local_h, kind=8), np_total, 1, mpi_integer8, &
                        mpi_sum, 0, mpi_comm_world, ierr)
        call mpi_allreduce(int(np_local_h,kind=8), np_h, 1, mpi_integer8, mpi_sum, &
                           mpi_comm_world, ierr)
    endif

    if (rank == 0) write(*,*) 'Total number of particles = ', np_total

    if (command == 0) then
        do j=1, np_local
            read(21) xvp(:,j)
        enddo
    else if (command == 1) then
        do j=1, np_local_dm
            read(21) xvp_dm(:,j)
        enddo
    else
        !! Read halo global coordinates and velocities
        do j=1, np_local_h
            read(21) xvmp_h(1:3, j)
            read(21) garbage4 !! m (x2) + r (x2) 
            read(21) garbage3 !! xbar (x3)
            read(21) xvmp_h(4:6, j)
            read(21) garbage3 !! angular momentum
            read(21) garbage3 !! var in vel
            read(21) garbage3 !! var in pos
            read(21) garbage3 !! moment of inertia
            read(21) garbage3 !! moment of inertia
            read(21) garbage3 !! xbar nu
            read(21) garbage3 !! vbar nu
            read(21) garbage1 !! nbr nu
        enddo
        !! Convert global halo coordinates to local node coordinates
        do j=1 , np_local_h
            xvmp_h(1:3, j) = xvmp_h(1:3, j) - slab_coord(:)*nc_node_dim
        enddo
    endif

    close(21)

#ifdef COARSE_HACK
    if (command == 0) then
        do j=1, np_local
            xvp(1:6,j) = xvp(1:6,j)/coarsen_factor
        enddo
    else if (command == 1) then
        do j=1, np_local_dm
            xvp_dm(1:6,j) = xvp_dm(1:6,j)/coarsen_factor
        enddo
    else
        do j=1, np_local_h
            xvmp_h(1:6,j) = xvmp_h(1:6,j)/coarsen_factor
        enddo
    endif
#endif
 
#ifdef KAISER

    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...

    if (command == 0) then
        xvp(3,:)=xvp(3,:) + xvp(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
    else if (command == 1) then
        xvp_dm(3,:)=xvp_dm(3,:) + xvp_dm(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
    else
     xvmp_h(3,:) = xvmp_h(3,:) + xvmp_h(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
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
    integer(8) :: num_elements_i8
    integer(4) :: num_elements
    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    !
    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
    !

    breakup = 1
    num_elements_i8 = int(nc_node_dim,kind=8) * nc_node_dim * nc_pen
    passGB = 4. * num_elements_i8 / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements_i8 / breakup

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
    integer(8) :: num_elements_i8
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
    num_elements_i8 = int(nc_node_dim,kind=8) * nc_node_dim * nc_pen
    passGB = 4. * num_elements_i8 / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements_i8 / breakup

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

subroutine writepowerspectra(command, index)
    !
    ! Writes the dimensionless power spectrum for the curl/divergence of the momentum density field
    !    

    implicit none
    
    integer      :: i, j, k
    character*180 :: fn
    character*3  :: prefix
    character*7  :: z_write
    real    :: vsim2phys, zcur
    integer(4) :: command, index

    !! Determine conversion factor for sim velocity to physical
    zcur      = z_checkpoint(cur_checkpoint)
    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

    if (rank == 0) write(*, *) "---> vsim2phys = ", vsim2phys

    !
    ! Determine name of output file
    !

    !! Determine suffix
    fn = '.dat'
#ifdef KAISER
    fn = '-RSD'//fn
#endif

    !! Determine which field we are dealing with
    if (index == 1) then
        fn = '_nunu'//fn
    else if (index == 2) then
        fn = '_dmdm'//fn
    else if (index == 3) then
        fn = '_haha'//fn
    else if (index == 4) then
        fn = '_nudm'//fn
    else if (index == 5) then
        fn = '_nuha'//fn
    else if (index == 6) then
        fn = '_dmha'//fn
    endif

    !! Determine whether we are dealing with total velocity or divergence
    if (command == 0) then 
        fn = '_totvel'//fn
    else if (command == 1) then 
        fn = '_divvel'//fn
#ifdef CURL
    else if (command == 2) then
        fn = '_crlvel'//fn
#endif
    endif

    write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
    z_write=adjustl(z_write)
    
#ifdef NGP 
    prefix = 'ngp'
#else
    prefix = 'cic'
#endif

    fn = output_path//z_write(1:len_trim(z_write))//prefix//fn

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
                pkdm(1, i) = pkdm(1, i) + pkvel(j, 1, i, index) 
                pkdm(2, i) = pkdm(2, i) + pkvel(j, 2, i, index)
            enddo
            pkdm(3, i) = pkvel(1, 3, i, index)
        enddo

    else if (command == 1) then

        do i = 1, nc
            pkdm(1, i) = pkdivg(1, i, index)
            pkdm(2, i) = pkdivg(2, i, index)
            pkdm(3, i) = pkdivg(3, i, index)
        enddo

#ifdef CURL
    else if (command == 2) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkcurl(j, 1, i, index)
                pkdm(2, i) = pkdm(2, i) + pkcurl(j, 2, i, index)
            enddo
            pkdm(3, i) = pkcurl(1, 3, i, index)
        enddo
#endif

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


subroutine darkmatter(command)

    implicit none

    integer :: i, j, k, command
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

    if (command == 0) then !! Take velocity field
        cube(:, :, :) = velden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
    else !! Take divergence field
        cube(:, :, :) = veldivg(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
    endif

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
      write(*,*) 'Cube min    ', dmint
      write(*,*) 'Cube max    ', dmaxt
      write(*,*) 'Cube sum ', real(dsum)
      write(*,*) 'Cube var ', real(dvar)
    endif

    ! 
    ! Forward FFT dm delta field
    !    

    call cp_fftw(1)

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

    integer i,pp,np_buf,np_exit,npo,npi
    integer(8) :: np_final
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
        else if (command == 1) then
            if (pp > np_local_dm) exit
        else
            if (pp > np_local_h) exit
        endif

        !! Read its position  
        if (command == 0) then
            x = xvp(:, pp)
        else if (command == 1) then
            x = xvp_dm(:, pp)
        else
            x = xvmp_h(:, pp)
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
            else if (command == 1) then
                xp_buf(:, np_buf) = xvp_dm(:, pp)
                xvp_dm(:, pp)     = xvp_dm(:, np_local_dm)
                np_local_dm       = np_local_dm - 1
            else
                xp_buf(:, np_buf) = xvmp_h(:, pp)
                xvmp_h(:, pp)     = xvmp_h(:, np_local_h)
                np_local_h       = np_local_h - 1
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
        else if (command == 1) then
            np_local_dm = np_local_dm + 1
            xvp_dm(:, np_local_dm) = x
        else
            np_local_h = np_local_h + 1
            xvmp_h(:, np_local_h) = x
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
        else if (command == 1) then
            np_local_dm = np_local_dm + 1
            xvp_dm(:, np_local_dm) = x
        else
            np_local_h = np_local_h + 1
            xvmp_h(:, np_local_h) = x
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
        else if (command == 1) then
            np_local_dm = np_local_dm + 1
            xvp_dm(:, np_local_dm) = x
        else
            np_local_h = np_local_h + 1
            xvmp_h(:, np_local_h) = x
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
        else if (command == 1) then
            np_local_dm = np_local_dm+1
            xvp_dm(:, np_local_dm) = x
        else
            np_local_h = np_local_h+1
            xvmp_h(:, np_local_h) = x
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
        else if (command == 1) then
            np_local_dm=np_local_dm+1
            xvp_dm(:,np_local_dm)=x
        else
            np_local_h=np_local_h+1
            xvmp_h(:,np_local_h)=x
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
        else if (command == 1) then
            np_local_dm=np_local_dm+1
            xvp_dm(:,np_local_dm)=x
        else
            np_local_h=np_local_h+1
            xvmp_h(:,np_local_h)=x
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

  if (command == 0) then
     call mpi_reduce(int(np_local,kind=8),np_final,1,mpi_integer8,mpi_sum,0, &
          mpi_comm_world,ierr)
  else if (command == 1) then
     call mpi_reduce(int(np_local_dm,kind=8),np_final,1,mpi_integer8,mpi_sum,0, &
          mpi_comm_world,ierr)
  else
     call mpi_reduce(int(np_local_h,kind=8),np_final,1,mpi_integer8,mpi_sum,0, &
          mpi_comm_world,ierr)
  endif

  if (rank == 0) then
     if (command == 0) then
        if (np_final /= int(np,kind=8)**3) then
           print *,'ERROR: total number of neutrinos incorrect after passing'
        endif
     else if (command == 1) then
        if (np_final /= (int(np,kind=8)/ratio_nudm_dim)**3) then
           print *,'ERROR: total number of dark matter particles incorrect after passing'
        endif
     else
        if (np_final /= np_h) then
           print *,'ERROR: total number of halos incorrect after passing'
        endif
     endif
  endif
 
!!  Check for particles out of bounds

  if (command == 0) then
     do i=1,np_local
        if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
             xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
             xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
           print *,'neutrino out of bounds',rank,i,xvp(:3,i),nc_node_dim
        endif
     enddo
  else if (command == 1) then
     do i=1,np_local_dm
        if (xvp_dm(1,i) < 0 .or. xvp_dm(1,i) >= nc_node_dim .or. &
             xvp_dm(2,i) < 0 .or. xvp_dm(2,i) >= nc_node_dim .or. &
             xvp_dm(3,i) < 0 .or. xvp_dm(3,i) >= nc_node_dim) then
           print *,'dark matter particle out of bounds',rank,i,xvp_dm(:3,i),nc_node_dim
        endif
     enddo
  else
     do i=1,np_local_h
        if (xvmp_h(1,i) < 0 .or. xvmp_h(1,i) >= nc_node_dim .or. &
             xvmp_h(2,i) < 0 .or. xvmp_h(2,i) >= nc_node_dim .or. &
             xvmp_h(3,i) < 0 .or. xvmp_h(3,i) >= nc_node_dim) then
           print *,'halo out of bounds',rank,i,xvmp_h(:3,i),nc_node_dim
        endif
     enddo
  endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) 'total buffered particles =', np_exit
    if (rank == 0) write(*,*) 'total particles          =', np_final
    if (rank == 0) write(*, *) "Finished pass_particles ... elapsed time = ", time2-time1

    return

end subroutine pass_particles

! -------------------------------------------------------------------------------------------------------

subroutine buffer_particles_groups(command, glook)
  !
  ! Exchange coarse mesh buffers for finding halos near node edges.
  ! Only deals one group at a time
  ! 

  !! It would be possible to exchange only the nearest particles ...

  implicit none

  integer :: pp,pp_down,pp_up
  integer :: np_buf, nppx, npmx, nppy, npmy, nppz, npmz
  integer, dimension(mpi_status_size) :: status, sstatus, rstatus
  integer :: tag, srequest, rrequest, sierr, rierr
  real(4), parameter :: eps = 1.0e-03
  integer(4) :: np_buffer_sent_local
  integer(8) :: np_buffer_sent
  integer(4) :: command
  integer(1) :: glook

  real(8) :: time1, time2
  time1 = mpi_wtime(ierr)

  ! --------------------------------------------------------
  ! Pass +x
  ! --------------------------------------------------------
  tag = 11
  np_buf = 0

  !! Determine bouding indices in particle array of local particles on the node.
  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1)
     endif
  endif

  !! Search for particles that need to be passed to the adjacent node.
  if (command == 0) then
     do pp = pp_down, pp_up
        if (xvp(1, pp) >= nc_node_dim-rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
        endif
     enddo
  else if (command == 1) then
     do pp = pp_down,pp_up
        if (xvp_dm(1, pp) >= nc_node_dim-rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
        endif
     enddo
  else if (command == 2) then
     do pp = pp_down, pp_up
        if (xvmp_h(1, pp) >= nc_node_dim-rnf_buf_h) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvmp_h(:, pp)
        endif
     enddo
  endif

  if (6*np_buf > np_buffer) then
     write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  nppx = np_buf

  call mpi_sendrecv_replace(nppx, 1, mpi_integer, cart_neighbor(6), tag, &
       cart_neighbor(5), tag, mpi_comm_world, status, ierr)

  if (command == 0) then
     if (np_local+nppx > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer nu particles (nppx): ", rank, np_local, nppx, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm+nppx > max_np_dm) then
        write(*,*) "ERROR: Not enough space to receive buffer dm particles (nppx): ", rank, np_local_dm, nppx, max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else
     if (np_local_h+nppx > max_np_h) then
        write(*,*) "ERROR: Not enough space to receive buffer halo (nppx): ", rank, np_local_h, nppx, max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(6), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(recv_buf, nppx*6, mpi_real, cart_neighbor(5), &
       tag, mpi_comm_world, rrequest, rierr)

  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)

  !! Start filling in buffer particles at pp_down+1
  if (glook == g0) then
     if (command == 0) then
        pp_down = np_groups(g0)
     else if (command == 1) then
        pp_down = np_groups_dm(g0)
     else if (command == 2) then
        pp_down = np_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_down = max_np-np_groups(g1)
     else if (command == 1) then
        pp_down = max_np_dm-np_groups_dm(g1)
     else if (command == 2) then
        pp_down = max_np_h-np_groups_h(g1)
     endif
     pp_down = pp_down - nppx
  endif

  if (command == 0) then
     do pp = 1, nppx
        xvp(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvp(1, pp_down+pp) = max(xvp(1,pp_down+pp)-nc_node_dim, -rnf_buf)
     enddo
  else if (command == 1) then
     do pp = 1, nppx
        xvp_dm(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvp_dm(1, pp_down+pp) = max(xvp_dm(1,pp_down+pp)-nc_node_dim, -rnf_buf)
     enddo
  else if (command == 2) then
     do pp = 1, nppx
        xvmp_h(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvmp_h(1, pp_down+pp) = max(xvmp_h(1,pp_down+pp)-nc_node_dim, -rnf_buf_h)
     enddo
  endif

  if (command == 0) then
     np_local = np_local + nppx
     np_buf_groups(glook) = nppx
  else if (command == 1) then
     np_local_dm = np_local_dm + nppx
     np_buf_groups_dm(glook) = nppx
  else if (command == 2) then
     np_local_h = np_local_h + nppx
     np_buf_groups_h(glook) = nppx
  endif

  ! --------------------------------------------------------
  ! Pass -x
  ! --------------------------------------------------------

  tag = 12
  np_buf = 0

  !! Determine bouding indices in particle array of local particles on the node.
  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1)
     endif
  endif

  !! Search for particles that need to be passed to the adjacent node.
  if (command == 0) then
     do pp = pp_down, pp_up
        if (xvp(1, pp) < rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
        endif
     enddo
  else if (command == 1) then
     do pp = pp_down,pp_up
        if (xvp_dm(1, pp) < rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
        endif
     enddo
  else if (command == 2) then
     do pp = pp_down,pp_up
        if (xvmp_h(1, pp) < rnf_buf_h) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvmp_h(:, pp)
        endif
     enddo
  endif

  if (6*np_buf > np_buffer) then
     write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  npmx = np_buf

  call mpi_sendrecv_replace(npmx, 1, mpi_integer, cart_neighbor(5), tag, &
       cart_neighbor(6), tag, mpi_comm_world, status, ierr)

  if (command == 0) then
     if (np_local+npmx > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer nu particles (npmx): ", rank, np_local, npmx, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm+npmx > max_np_dm) then
        write(*,*) "ERROR: Not enough space to receive buffer dm particles (npmx): ", rank, np_local_dm, npmx, max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 2) then
     if (np_local_h+npmx > max_np_h) then
        write(*,*) "ERROR: Not enough space to receive buffer halos (npmx): ", rank, np_local_h, npmx, max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(5), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(recv_buf, npmx*6, mpi_real, cart_neighbor(6), &
       tag, mpi_comm_world, rrequest, rierr)

  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)

  !! Start filling in buffer particles at pp_down+1
  if (glook == g0) then
     if (command == 0) then
        pp_down = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_down = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_down = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_down = max_np-np_groups(g1)-np_buf_groups(g1)
     else if (command == 1) then
        pp_down = max_np_dm-np_groups_dm(g1)-np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_down = max_np_h-np_groups_h(g1)-np_buf_groups_h(g1)
     endif
     pp_down = pp_down - npmx
  endif

  if (command == 0) then
     do pp = 1, npmx
        xvp(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvp(1, pp_down+pp)) .lt. eps) then
           if(xvp(1, pp_down+pp) < 0.) then
              xvp(1, pp_down+pp) = -eps
           else
              xvp(1, pp_down+pp) = eps
           endif
        endif
        xvp(1, pp_down+pp) = min(xvp(1, pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf-eps)
     enddo
  else if (command == 1) then
     do pp = 1, npmx
        xvp_dm(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvp_dm(1, pp_down+pp)) .lt. eps) then
           if(xvp_dm(1, pp_down+pp) < 0.) then
              xvp_dm(1, pp_down+pp) = -eps
           else
              xvp_dm(1, pp_down+pp) = eps
           endif
        endif
        xvp_dm(1, pp_down+pp) = min(xvp_dm(1, pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf-eps)
     enddo
  else if (command == 2) then
     do pp = 1, npmx
        xvmp_h(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvmp_h(1, pp_down+pp)) .lt. eps) then
           if(xvmp_h(1, pp_down+pp) < 0.) then
              xvmp_h(1, pp_down+pp) = -eps
           else
              xvmp_h(1, pp_down+pp) = eps
           endif
        endif
        xvmp_h(1, pp_down+pp) = min(xvmp_h(1, pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf_h-eps)
     enddo
  endif

  if (command == 0) then
     np_local = np_local + npmx
     np_buf_groups(glook) = np_buf_groups(glook) + npmx
  else if (command == 1) then
     np_local_dm = np_local_dm + npmx
     np_buf_groups_dm(glook) = np_buf_groups_dm(glook) + npmx
  else if (command == 2) then
     np_local_h = np_local_h + npmx
     np_buf_groups_h(glook) = np_buf_groups_h(glook) + npmx
  endif

  ! --------------------------------------------------------
  ! Pass +y
  ! --------------------------------------------------------

  tag = 13
  np_buf = 0

  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
  endif

  if (command == 0) then
     do pp = pp_down, pp_up
        if (xvp(2, pp) >= nc_node_dim-rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
        endif
     enddo
  else if (command == 1) then
     do pp = pp_down, pp_up
        if (xvp_dm(2, pp) >= nc_node_dim-rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
        endif
     enddo
  else if (command == 2) then
     do pp = pp_down, pp_up
        if (xvmp_h(2, pp) >= nc_node_dim-rnf_buf_h) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvmp_h(:, pp)
        endif
     enddo
  endif

  if (6*np_buf > np_buffer) then
     write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  nppy = np_buf

  call mpi_sendrecv_replace(nppy, 1, mpi_integer, cart_neighbor(4), tag, &
       cart_neighbor(3), tag, mpi_comm_world, status, ierr)

  if (command == 0) then
     if (np_local+nppy > max_np) then
        write(*,*) "ERROR: Not enough space to receive nu buffer particles (nppy): ", rank, np_local, nppy, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm+nppy > max_np_dm) then
        write(*,*) "ERROR: Not enough space to receive dm buffer particles (nppy): ", rank, np_local_dm, nppy, max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 2) then
     if (np_local_h+nppy > max_np_h) then
        write(*,*) "ERROR: Not enough space to receive buffer halos (nppy): ", rank, np_local_h, nppy, max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(4), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(recv_buf, nppy*6, mpi_real, cart_neighbor(3), &
       tag, mpi_comm_world, rrequest, rierr)

  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)

  if (glook == g0) then
     if (command == 0) then
        pp_down = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_down = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_down = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_down = max_np-np_groups(g1)-np_buf_groups(g1)
     else if (command == 1) then
        pp_down = max_np_dm-np_groups_dm(g1)-np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_down = max_np_h-np_groups_h(g1)-np_buf_groups_h(g1)
     endif
     pp_down = pp_down - nppy
  endif

  if (command == 0) then
     do pp = 1, nppy
        xvp(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvp(2, pp_down+pp) = max(xvp(2,pp_down+pp)-nc_node_dim, -rnf_buf)
     enddo
  else if (command == 1) then
     do pp = 1, nppy
        xvp_dm(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvp_dm(2, pp_down+pp) = max(xvp_dm(2,pp_down+pp)-nc_node_dim, -rnf_buf)
     enddo
  else
     do pp = 1, nppy
        xvmp_h(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvmp_h(2, pp_down+pp) = max(xvmp_h(2,pp_down+pp)-nc_node_dim, -rnf_buf_h)
     enddo
  endif

  if (command == 0) then
     np_local = np_local + nppy
     np_buf_groups(glook) = np_buf_groups(glook) + nppy
  else if (command == 1) then
     np_local_dm = np_local_dm + nppy
     np_buf_groups_dm(glook) = np_buf_groups_dm(glook) + nppy
  else if (command == 2) then
     np_local_h = np_local_h + nppy
     np_buf_groups_h(glook) = np_buf_groups_h(glook) + nppy
  endif

  ! --------------------------------------------------------  
  ! Pass -y
  ! --------------------------------------------------------  

  tag = 14
  np_buf = 0

  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
     pp_up = pp_up - nppy
  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
     pp_down = pp_down + nppy
  endif

  if (command == 0) then
     do pp = pp_down, pp_up
        if (xvp(2, pp) < rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
        endif
     enddo
  else if (command == 1) then
     do pp = pp_down, pp_up
        if (xvp_dm(2, pp) < rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
        endif
     enddo
  else if (command == 2) then
     do pp = pp_down, pp_up
        if (xvmp_h(2, pp) < rnf_buf_h) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvmp_h(:, pp)
        endif
     enddo
  endif

  if (6*np_buf > np_buffer) then
     write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  npmy = np_buf

  call mpi_sendrecv_replace(npmy, 1, mpi_integer, cart_neighbor(3), tag, &
       cart_neighbor(4), tag, mpi_comm_world, status, ierr)

  if (command == 0) then
     if (np_local+npmy > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer nu particles (npmy): ", rank, np_local, npmy, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm+npmy > max_np_dm) then
        write(*,*) "ERROR: Not enough space to receive buffer dm particles (npmy): ", rank, np_local_dm, npmy, max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 2) then
     if (np_local_h+npmy > max_np_h) then
        write(*,*) "ERROR: Not enough space to receive buffer halos (npmy): ", rank, np_local_h, npmy, max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(3), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(recv_buf, npmy*6, mpi_real, cart_neighbor(4), &
       tag, mpi_comm_world, rrequest, rierr)

  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)

  if (glook == g0) then
     if (command == 0) then
        pp_down = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_down = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_down = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_down = max_np-np_groups(g1)-np_buf_groups(g1)
     else if (command == 1) then
        pp_down = max_np_dm-np_groups_dm(g1)-np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_down = max_np_h-np_groups_h(g1)-np_buf_groups_h(g1)
     endif
     pp_down = pp_down - npmy
  endif

  if (command == 0) then
     do pp = 1, npmy
        xvp(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvp(2, pp_down+pp)) .lt. eps) then
           if(xvp(2, pp_down+pp) < 0.) then
              xvp(2, pp_down+pp) = -eps
           else
              xvp(2, pp_down+pp) = eps
           endif
        endif
        xvp(2, pp_down+pp) = min(xvp(2,pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf-eps)
     enddo
  else if (command == 1) then
     do pp = 1, npmy
        xvp_dm(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvp_dm(2, pp_down+pp)) .lt. eps) then
           if(xvp_dm(2, pp_down+pp) < 0.) then
              xvp_dm(2, pp_down+pp) = -eps
           else
              xvp_dm(2, pp_down+pp) = eps
           endif
        endif
        xvp_dm(2, pp_down+pp) = min(xvp_dm(2,pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf-eps)
     enddo
  else if (command == 2) then
     do pp = 1, npmy
        xvmp_h(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvmp_h(2, pp_down+pp)) .lt. eps) then
           if(xvmp_h(2, pp_down+pp) < 0.) then
              xvmp_h(2, pp_down+pp) = -eps
           else
              xvmp_h(2, pp_down+pp) = eps
           endif
        endif
        xvmp_h(2, pp_down+pp) = min(xvmp_h(2,pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf_h-eps)
     enddo
  endif

  if (command == 0) then
     np_local = np_local + npmy
     np_buf_groups(glook) = np_buf_groups(glook) + npmy
  else if (command == 1) then
     np_local_dm = np_local_dm + npmy
     np_buf_groups_dm(glook) = np_buf_groups_dm(glook) + npmy
  else
     np_local_h = np_local_h + npmy
     np_buf_groups_h(glook) = np_buf_groups_h(glook) + npmy
  endif

  ! --------------------------------------------------------    
  ! Pass +z
  ! --------------------------------------------------------    

  tag = 15
  np_buf = 0

  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
  endif

  if (command == 0) then
     do pp = pp_down, pp_up
        if (xvp(3, pp) >= nc_node_dim-rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
        endif
     enddo
  else if (command == 1) then
     do pp = pp_down, pp_up
        if (xvp_dm(3, pp) >= nc_node_dim-rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
        endif
     enddo
  else if (command == 2) then
     do pp = pp_down, pp_up
        if (xvmp_h(3, pp) >= nc_node_dim-rnf_buf_h) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvmp_h(:, pp)
        endif
     enddo
  endif

  if (6*np_buf > np_buffer) then
     write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  nppz = np_buf

  call mpi_sendrecv_replace(nppz, 1, mpi_integer, cart_neighbor(2), tag, &
       cart_neighbor(1), tag, mpi_comm_world, status, ierr)

  if (command == 0) then
     if (np_local+nppz > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer nu particles (nppz):", rank, np_local, nppz, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm+nppz > max_np_dm) then
        write(*,*) "ERROR: Not enough space to receive buffer dm particles (nppz):", rank, np_local_dm, nppz, max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 2) then
     if (np_local_h+nppz > max_np_h) then
        write(*,*) "ERROR: Not enough space to receive buffer halos (nppz):", rank, np_local_h, nppz, max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(2), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(recv_buf, nppz*6, mpi_real, cart_neighbor(1), &
       tag, mpi_comm_world, rrequest, rierr)

  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)

  if (glook == g0) then
     if (command == 0) then
        pp_down = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_down = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_down = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_down = max_np-np_groups(g1)-np_buf_groups(g1)
     else if (command == 1) then
        pp_down = max_np_dm-np_groups_dm(g1)-np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_down = max_np_h-np_groups_h(g1)-np_buf_groups_h(g1)
     endif
     pp_down = pp_down - nppz
  endif

  if (command == 0) then
     do pp = 1, nppz
        xvp(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvp(3, pp_down+pp) = max(xvp(3,pp_down+pp)-nc_node_dim, -rnf_buf)
     enddo
  else if (command == 1) then
     do pp = 1, nppz
        xvp_dm(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvp_dm(3, pp_down+pp) = max(xvp_dm(3,pp_down+pp)-nc_node_dim, -rnf_buf)
     enddo
  else if (command == 2) then
     do pp = 1, nppz
        xvmp_h(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        xvmp_h(3, pp_down+pp) = max(xvmp_h(3,pp_down+pp)-nc_node_dim, -rnf_buf_h)
     enddo
  endif

  if (command == 0) then
     np_local = np_local + nppz
     np_buf_groups(glook) = np_buf_groups(glook) + nppz
  else if (command == 1) then
     np_local_dm = np_local_dm + nppz
     np_buf_groups_dm(glook) = np_buf_groups_dm(glook) + nppz
  else if (command == 2) then
     np_local_h = np_local_h + nppz
     np_buf_groups_h(glook) = np_buf_groups_h(glook) + nppz
  endif

  ! --------------------------------------------------------    
  ! Pass -z
  ! --------------------------------------------------------    

  tag = 16
  np_buf = 0

  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
     pp_up = pp_up - nppz
  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
     pp_down = pp_down + nppz
  endif

  if (command == 0) then
     do pp = pp_down, pp_up
        if (xvp(3, pp) < rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
        endif
     enddo
  else if (command == 1) then
     do pp = pp_down, pp_up
        if (xvp_dm(3, pp) < rnf_buf) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvp_dm(:, pp)
        endif
     enddo
  else if (command == 2) then
     do pp = pp_down, pp_up
        if (xvmp_h(3, pp) < rnf_buf_h) then
           np_buf = np_buf + 1
           send_buf((np_buf-1)*6+1:np_buf*6) = xvmp_h(:, pp)
        endif
     enddo
  endif

  if (6*np_buf > np_buffer) then
     write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  npmz = np_buf

  call mpi_sendrecv_replace(npmz, 1, mpi_integer, cart_neighbor(1), tag, &
       cart_neighbor(2), tag, mpi_comm_world, status, ierr)

  if (command == 0) then
     if (np_local+npmz > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer nu particles (npmz): ", rank, np_local, npmz, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm+npmz > max_np_dm) then
        write(*,*) "ERROR: Not enough space to receive buffer dm particles (npmz): ", rank, np_local_dm, npmz, max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else
     if (np_local_h+npmz > max_np_h) then
        write(*,*) "ERROR: Not enough space to receive halos (npmz): ", rank, np_local_h, npmz, max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(1), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(recv_buf, npmz*6, mpi_real, cart_neighbor(2), &
       tag, mpi_comm_world, rrequest, rierr)

  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)

  if (glook == g0) then
     if (command == 0) then
        pp_down = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_down = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_down = np_groups_h(g0) + np_buf_groups_h(g0)
     endif
  else if (glook == g1) then
     if (command == 0) then
        pp_down = max_np-np_groups(g1)-np_buf_groups(g1)
     else if (command == 1) then
        pp_down = max_np_dm-np_groups_dm(g1)-np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_down = max_np_h-np_groups_h(g1)-np_buf_groups_h(g1)
     endif
     pp_down = pp_down - npmz
  endif

  if (command == 0) then
     do pp = 1, npmz
        xvp(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvp(3, pp_down+pp)) .lt. eps) then
           if(xvp(3, pp_down+pp) < 0.) then
              xvp(3, pp_down+pp) = -eps
           else
              xvp(3, pp_down+pp) = eps
           endif
        endif
        xvp(3,pp_down+pp) = min(xvp(3,pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf-eps)
     enddo
  else if (command == 1) then
     do pp = 1, npmz
        xvp_dm(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvp_dm(3, pp_down+pp)) .lt. eps) then
           if(xvp_dm(3, pp_down+pp) < 0.) then
              xvp_dm(3, pp_down+pp) = -eps
           else
              xvp_dm(3, pp_down+pp) = eps
           endif
        endif
        xvp_dm(3,pp_down+pp) = min(xvp_dm(3,pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf-eps)
     enddo
  else if (command == 2) then
     do pp = 1, npmz
        xvmp_h(:, pp_down+pp) = recv_buf((pp-1)*6+1:(pp-1)*6+6)
        if (abs(xvmp_h(3, pp_down+pp)) .lt. eps) then
           if(xvmp_h(3, pp_down+pp) < 0.) then
              xvmp_h(3, pp_down+pp) = -eps
           else
              xvmp_h(3, pp_down+pp) = eps
           endif
        endif
        xvmp_h(3,pp_down+pp) = min(xvmp_h(3,pp_down+pp)+real(nc_node_dim,4), &
             nc_node_dim+rnf_buf_h-eps)
     enddo
  endif

  if (command == 0) then
     np_local = np_local + npmz
     np_buf_groups(glook) = np_buf_groups(glook) + npmz
  else if (command == 1) then
     np_local_dm = np_local_dm + npmz
     np_buf_groups_dm(glook) = np_buf_groups_dm(glook) + npmz
  else
     np_local_h = np_local_h + npmz
     np_buf_groups_h(glook) = np_buf_groups_h(glook) + npmz
  endif

  ! --------------------------------------------------------    
  ! Calculate some statistics
  ! --------------------------------------------------------    

  if (command == 0) then
     np_buffer_sent_local = np_buf_groups(glook)
  else if (command == 1) then
     np_buffer_sent_local = np_buf_groups_dm(glook)
  else
     np_buffer_sent_local = np_buf_groups_h(glook)
  endif
  call mpi_reduce(int(np_buffer_sent_local, kind=8), np_buffer_sent, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

  time2 = mpi_wtime(ierr)
  if (rank == 0) then
     write(*,*) 'Total Buffer Particles : ', np_buffer_sent
     write(*,*) 'Finished buffer_particles_groups ... time elapsed = ', time2-time1
  endif

end subroutine buffer_particles_groups

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

subroutine velocity_density(cdim, command, glook, nfind)
    !
    ! Determine velocity field by taking the average velocity of the closest particles to each grid centre.
    ! 

    implicit none

    integer :: ind, i, j, k, pp, thread, ni
    integer :: ic, jc, kc, pp_up, pp_down, pp_down0
    integer :: npart
    real, dimension(3) :: dr, rc
    real :: vx, dmax, xvtemp(6)
    integer(4) :: cdim, command, nfind
    integer(1) :: glook
    logical :: converged, pb
    integer(8) :: num_notconverged, num2, sameR, num1
    real*8 :: rsum, rsumt
    real    :: v, r2, vswap, r2swap
    integer :: m,mm

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

  if (rank == 0) write(*,*) "Starting velocity_density c, g, n = ", command, glook, nfind

  ! --------------------------------------------------------------------------------------
  ! Initialize to zeros + make hoc
  ! --------------------------------------------------------------------------------------

  call make_hoc(command, glook)

  do i = 0, nc_node_dim + 1
     velden(:, :, i) = 0.
  enddo

  if (glook == g0) then
     pp_down0 = 1
  else if (glook == g1) then
     if (command == 0) then
        pp_down0 = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_down0 = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_down0 = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
  endif

  ! --------------------------------------------------------------------------------------
  ! For each cell find the closest nfind particles and average their velocities
  ! --------------------------------------------------------------------------------------

  num_notconverged = 0
  sameR = 0
  rsum = 0.

  !$omp  parallel num_threads(nt) default(shared) private(i, j, k, thread, rc, npart, ic, jc, kc, ind, pp, dr, vx, converged, pb, &
  !$omp&      dmax, v, r2, vswap, r2swap, m, mm, ni, xvtemp, pp_up, pp_down) reduction(+:num_notconverged) reduction(+:sameR) reduction(+:rsum)
  thread = 1
  thread = omp_get_thread_num() + 1

  if (command == 0) then !! NEUTRINOS

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
              pb = .false.

              !! If dmax > 0, it will represent the max distance to search. Gets set once nfind particles are found
              dmax = -1.0 

              do ind = 1, num_ngbhs

                 !! If we are too far, we can already exit the loop
                 if (dmax >= 0 .and. cell_search_r(ind) >= dmax) exit 

                 !! Cell indices to search within
                 ic = int((i-1)/mesh_scale) + 1 + cell_search(1, ind)
                 jc = int((j-1)/mesh_scale) + 1 + cell_search(2, ind)
                 kc = int((k-1)/mesh_scale) + 1 + cell_search(3, ind)

                 !! Start and finish particle indices within this cell
                 if ( ic == hoc_nc_l .and. jc == hoc_nc_l .and. kc == hoc_nc_l ) then
                    pp_down = pp_down0
                 else if (ic == hoc_nc_l .and. jc == hoc_nc_l) then
                    pp_down = hoc(hoc_nc_h, hoc_nc_h, kc-1) + 1
                 else if (ic == hoc_nc_l) then
                    pp_down = hoc(hoc_nc_h, jc-1, kc) + 1
                 else
                    pp_down = hoc(ic-1, jc, kc) + 1
                 endif
                    pp_up = hoc(ic, jc, kc)

                 !! Loop over particles in this coarse mesh cell
                 do pp = pp_down, pp_up

                    xvtemp(:) = xvp(:,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvtemp(3+cdim)
                    npart = npart + 1

                    if (npart == 1) then !! If it's the first particle we found, just add it to the list 
                       rpos(thread, 1, 1) = v
                       rpos(thread, 2, 1) = r2
                    else if ((npart - 1) < nfind) then !! Else just append it to its right place
                       m = 1
                       do while (m < npart)
                          if (r2 >= rpos(thread, 2, m)) exit
                          m = m + 1
                       enddo
                       do while (m < npart)
                          vswap = rpos(thread, 1, m)
                          r2swap = rpos(thread, 2, m)
                          rpos(thread, 1, m) = v
                          rpos(thread, 2, m) = r2
                          v = vswap
                          r2 = r2swap
                          m = m + 1
                       enddo
                       rpos(thread, 1, npart) = v
                       rpos(thread, 2, npart) = r2
                    else if (r2 < rpos(thread, 2, 1)) then !! Else push it to its right place and get rid of the furtest particle
                       m = 2
                       do while (m <= nfind)
                          if (r2 >= rpos(thread, 2, m)) exit
                          m = m+1
                       enddo
                       if (m == 2) then
                          pb = .false.
                       endif
                       do mm = 1, m-1
                          vswap = rpos(thread, 1, m-mm)
                          r2swap = rpos(thread, 2, m-mm)
                          rpos(thread, 1, m-mm) = v
                          rpos(thread, 2, m-mm) = r2
                          v = vswap
                          r2 = r2swap
                       enddo
                    else if (r2 == rpos(thread, 2, 1)) then
                       pb = .true.
                    endif
                 enddo !! pp

                 if (npart >= nfind) then
                    dmax = sqrt(rpos(thread, 2, 1))
                    converged = .true.
                 endif

              enddo !! ind

              !! Average velocity and distance of these particles
              v = 0.
              r2 = 0.
              if (npart > 0) then
                 ni = min(nfind, npart)
                 do ind = 1, ni
                    v = v + rpos(thread, 1, ind)
                 enddo
                 velden(i, j, k) = v / ni 
                 do ind = 1, ni 
                    r2 = r2 + sqrt(rpos(thread, 2, ind))
                 enddo
                 rsum = rsum + r2 / ni 
              else
                 rsum = rsum + nfine_buf*sqrt(3.0)
              endif

              !! Tally number of cells for which the closest nfind particles were not found
              if (.not. converged ) then
                 num_notconverged = num_notconverged + 1
              endif
              
              !! Tally the number of cells for which the first particle not included has the same distance as the last particle included
              if (pb) then
                 sameR = sameR + 1
              endif

           enddo !! i
        enddo !! j
     enddo !! k
     !$omp end do  

  else if (command == 1) then !! DARK MATTER

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
              pb = .false.

              !! If dmax > 0, it will represent the max distance to search. Gets set once nfind particles are found
              dmax = -1.0

              do ind = 1, num_ngbhs

                 !! If we are too far, we can already exit the loop
                 if (dmax >= 0 .and. cell_search_r(ind) >= dmax) exit

                 !! Cell indices to search within
                 ic = int((i-1)/mesh_scale) + 1 + cell_search(1, ind)
                 jc = int((j-1)/mesh_scale) + 1 + cell_search(2, ind)
                 kc = int((k-1)/mesh_scale) + 1 + cell_search(3, ind)

                 !! Start and finish particle indices within this cell
                 if ( ic == hoc_nc_l .and. jc == hoc_nc_l .and. kc == hoc_nc_l ) then
                    pp_down = pp_down0
                 else if (ic == hoc_nc_l .and. jc == hoc_nc_l) then
                    pp_down = hoc_dm(hoc_nc_h, hoc_nc_h, kc-1) + 1
                 else if (ic == hoc_nc_l) then
                    pp_down = hoc_dm(hoc_nc_h, jc-1, kc) + 1
                 else
                    pp_down = hoc_dm(ic-1, jc, kc) + 1
                 endif
                    pp_up = hoc_dm(ic, jc, kc)

                 !! Loop over particles in this coarse mesh cell
                 do pp = pp_down, pp_up

                    xvtemp(:) = xvp_dm(:,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvtemp(3+cdim)
                    npart = npart + 1

                    if (npart == 1) then !! If it's the first particle we found, just add it to the list 
                       rpos(thread, 1, 1) = v
                       rpos(thread, 2, 1) = r2
                    else if ((npart - 1) < nfind) then !! Else just append it to its right place
                       m = 1
                       do while (m < npart)
                          if (r2 >= rpos(thread, 2, m)) exit
                          m = m + 1
                       enddo
                       do while (m < npart)
                          vswap = rpos(thread, 1, m)
                          r2swap = rpos(thread, 2, m)
                          rpos(thread, 1, m) = v
                          rpos(thread, 2, m) = r2
                          v = vswap
                          r2 = r2swap
                          m = m + 1
                       enddo
                       rpos(thread, 1, npart) = v
                       rpos(thread, 2, npart) = r2
                    else if (r2 < rpos(thread, 2, 1)) then !! Else push it to its right place and get rid of the furtest particle
                       m = 2
                       do while (m <= nfind)
                          if (r2 >= rpos(thread, 2, m)) exit
                          m = m+1
                       enddo
                       if (m == 2) then
                          pb = .false.
                       endif
                       do mm = 1, m-1
                          vswap = rpos(thread, 1, m-mm)
                          r2swap = rpos(thread, 2, m-mm)
                          rpos(thread, 1, m-mm) = v
                          rpos(thread, 2, m-mm) = r2
                          v = vswap
                          r2 = r2swap
                       enddo
                    else if (r2 == rpos(thread, 2, 1)) then
                       pb = .true.
                    endif
                 enddo !! pp

                 if (npart >= nfind) then
                    dmax = sqrt(rpos(thread, 2, 1))
                    converged = .true.
                 endif

              enddo !! ind

              !! Average velocity and distance of these particles
              v = 0.
              r2 = 0.
              if (npart > 0) then
                 ni = min(nfind, npart)
                 do ind = 1, ni
                    v = v + rpos(thread, 1, ind)
                 enddo
                 velden(i, j, k) = v / ni
                 do ind = 1, ni 
                    r2 = r2 + sqrt(rpos(thread, 2, ind))
                 enddo
                 rsum = rsum + r2 / ni
              else
                 rsum = rsum + nfine_buf*sqrt(3.0)
              endif

              !! Tally number of cells for which the closest nfind particles were not found
              if (.not. converged ) then
                 num_notconverged = num_notconverged + 1
              endif

              !! Tally the number of cells for which the first particle not included has the same distance as the last particle included
              if (pb) then
                 sameR = sameR + 1
              endif

           enddo !! i
        enddo !! j
     enddo !! k
     !$omp end do 

  else if (command == 2) then !! HALOS

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
              pb = .false.

              !! If dmax > 0, it will represent the max distance to search. Gets set once nfind particles are found
              dmax = -1.0

              do ind = 1, num_ngbhs_h

                 !! If we are too far, we can already exit the loop
                 if (dmax >= 0 .and. cell_search_r_h(ind) >= dmax) exit

                 !! Cell indices to search within
                 ic = int((i-1)/mesh_scale_h) + 1 + cell_search_h(1, ind)
                 jc = int((j-1)/mesh_scale_h) + 1 + cell_search_h(2, ind)
                 kc = int((k-1)/mesh_scale_h) + 1 + cell_search_h(3, ind)

                 !! Start and finish particle indices within this cell
                 if ( ic == hoc_nc_l_h .and. jc == hoc_nc_l_h .and. kc == hoc_nc_l_h ) then
                    pp_down = pp_down0
                 else if (ic == hoc_nc_l_h .and. jc == hoc_nc_l_h) then
                    pp_down = hoc_h(hoc_nc_h_h, hoc_nc_h_h, kc-1) + 1
                 else if (ic == hoc_nc_l_h) then
                    pp_down = hoc_h(hoc_nc_h_h, jc-1, kc) + 1
                 else
                    pp_down = hoc_h(ic-1, jc, kc) + 1
                 endif
                    pp_up = hoc_h(ic,jc,kc)

                 !! Loop over particles in this coarse mesh cell
                 do pp = pp_down, pp_up

                    xvtemp(:) = xvmp_h(:,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvtemp(3+cdim)
                    npart = npart + 1

                    if (npart == 1) then !! If it's the first particle we found, just add it to the list 
                       rpos(thread, 1, 1) = v
                       rpos(thread, 2, 1) = r2
                    else if ((npart - 1) < nfind) then !! Else just append it to its right place
                       m = 1
                       do while (m < npart)
                          if (r2 >= rpos(thread, 2, m)) exit
                          m = m + 1
                       enddo
                       do while (m < npart)
                          vswap = rpos(thread, 1, m)
                          r2swap = rpos(thread, 2, m)
                          rpos(thread, 1, m) = v
                          rpos(thread, 2, m) = r2
                          v = vswap
                          r2 = r2swap
                          m = m + 1
                       enddo
                       rpos(thread, 1, npart) = v
                       rpos(thread, 2, npart) = r2
                    else if (r2 < rpos(thread, 2, 1)) then !! Else push it to its right place and get rid of the furtest particle
                       m = 2
                       do while (m <= nfind)
                          if (r2 >= rpos(thread, 2, m)) exit
                          m = m+1
                       enddo
                       if (m == 2) then
                          pb = .false.
                       endif
                       do mm = 1, m-1
                          vswap = rpos(thread, 1, m-mm)
                          r2swap = rpos(thread, 2, m-mm)
                          rpos(thread, 1, m-mm) = v
                          rpos(thread, 2, m-mm) = r2
                          v = vswap
                          r2 = r2swap
                       enddo
                    else if (r2 == rpos(thread, 2, 1)) then
                       pb = .true.
                    endif
                 enddo !! pp

                 if (npart >= nfind) then
                    dmax = sqrt(rpos(thread, 2, 1))
                    converged = .true.
                 endif

              enddo !! ind

              !! Average velocity and distance of these particles
              v = 0.
              r2 = 0.
              if (npart > 0) then
                 ni = min(nfind, npart)
                 do ind = 1, ni
                    v = v + rpos(thread, 1, ind)
                 enddo
                 velden(i, j, k) = v / ni
                 do ind = 1, ni 
                    r2 = r2 + sqrt(rpos(thread, 2, ind))
                 enddo
                 rsum = rsum + r2 / ni 
              else
                 rsum = rsum + nfine_buf*sqrt(3.0)
              endif

              !! Tally number of cells for which the closest nfind particles were not found
              if (.not. converged ) then
                 num_notconverged = num_notconverged + 1
              endif

              !! Tally the number of cells for which the first particle not included has the same distance as the last particle included
              if (pb) then
                 sameR = sameR + 1
              endif

           enddo !! i
        enddo !! j
     enddo !! k
     !$omp end do  

  endif !! command
  !$omp end parallel

  !! Tally up statistics across nodes
  call mpi_reduce(num_notconverged, num2, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(sameR, num1, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(rsum, rsumt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  rsum = rsumt / nc**3

  time2 = mpi_wtime(ierr)
  if (rank == 0) then
     if (num2 /= 0) then
        write(*,*) "WARNING: num_notconverged = ", num2, " ... consider increasing nfine_buf !!"
        write(*,*) "nfind = ", nfind
     endif
     write(*,*) "rmean = ", rsum
     if (num1 /= 0) then
        write(*,*) 'WARNING: num sameR = ', num1
     endif
     write(*,*) "Finished velocity_density ... elapsed time = ", time2-time1
  endif

  return

end subroutine velocity_density

! -------------------------------------------------------------------------------------------------------

subroutine pass_veldensity
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

        velden_send_buff(:, :) = velden(nc_node_dim, :, :)

        call mpi_isend(velden_send_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(velden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        velden(0, :, :) = velden_recv_buff(:, :)

        !
        ! Pass -x
        ! 

        tag = 112

        velden_send_buff(:, :) = velden(1, :, :)

        call mpi_isend(velden_send_buff, num2send, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(velden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        velden(nc_node_dim+1, :, :) = velden_recv_buff(:, :)

        !
        ! Pass +y
        ! 

        tag = 113

        velden_send_buff(:, :) = velden(:, nc_node_dim, :)

        call mpi_isend(velden_send_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(velden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        velden(:, 0, :) = velden_recv_buff(:, :)

        !
        ! Pass -y
        ! 

        tag = 114

        velden_send_buff(:, :) = velden(:, 1, :)

        call mpi_isend(velden_send_buff, num2send, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(velden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        velden(:, nc_node_dim+1, :) = velden_recv_buff(:, :)

        !
        ! Pass +z
        ! 

        tag = 115

        velden_send_buff(:, :) = velden(:, :, nc_node_dim)

        call mpi_isend(velden_send_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(velden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        velden(:, :, 0) = velden_recv_buff(:, :)

        !
        ! Pass -z
        ! 

        tag = 116

        velden_send_buff(:, :) = velden(:, :, 1)

        call mpi_isend(velden_send_buff, num2send, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
        call mpi_irecv(velden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)
        call mpi_wait(srequest, sstatus, sierr)
        call mpi_wait(rrequest, rstatus, rierr)

        velden(:, :, nc_node_dim+1) = velden_recv_buff(:, :)

    enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished pass_veldensity ... elapsed time = ", time2-time1

    return

end subroutine pass_veldensity

! -------------------------------------------------------------------------------------------------------

subroutine velocity_divergence
    !
    ! Calculates the divergence of the velocity field for one component at a
    ! time. If we are starting with the x component (cur_dimension == 1) then
    ! first initialize the velocity divergence array to 0. Otherwise, we add
    ! on to the array. Afterwards the field is Fourier transformed and the 
    ! power spectrum is computed.
    !

    implicit none
    integer :: k, j, i
    real :: dx

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize to zeros if we are doing the x component
    !

    if (cur_dimension == 1) then
        do k = 1, nc_node_dim
            veldivg(:, :, k) = 0.
        enddo
    endif

    !
    ! Compute derivatives using dy/dx(a) = y(a+1) - y(a-1) / x(a+1) - x(a-1)
    !

    !! Spacing between grid points
    dx = 2. * box / nc

    if (cur_dimension == 1) then
        do i = 1, nc_node_dim
            do j = 1, nc_node_dim
                do k = 1, nc_node_dim
                    veldivg(i, j, k) = (velden(i+1, j, k) - velden(i-1, j, k)) / dx
                enddo
            enddo
        enddo
    else if (cur_dimension == 2) then
        do i = 1, nc_node_dim
            do j = 1, nc_node_dim
                do k = 1, nc_node_dim
                    veldivg(i, j, k) = veldivg(i, j, k) + (velden(i, j+1, k) - velden(i, j-1, k)) / dx
                enddo
            enddo
        enddo
    else if (cur_dimension == 3) then
        do i = 1, nc_node_dim
            do j = 1, nc_node_dim
                do k = 1, nc_node_dim
                    veldivg(i, j, k) = veldivg(i, j, k) + (velden(i, j, k+1) - velden(i, j, k-1)) / dx
                enddo
            enddo
        enddo
    endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished velocity_divergence ... elapsed time = ", time2-time1

    return

end subroutine velocity_divergence

! -------------------------------------------------------------------------------------------------------
#ifdef CURL
subroutine velocity_curl
    !
    ! Calculates the cur_dimension componenent of the curl of the velocity field
    ! with first curlcom dimension stored in velden3 and second curlcom
    ! dimension stored in velden. 

    implicit none
    integer :: k, j, i
    real :: dx

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize to zeros using veldivg array to store curl componenet 
    !
    
    do k = 1, nc_node_dim
        veldivg(:, :, k) = 0.
    enddo

    !
    ! Compute derivatives using dy/dx(a) = y(a+1) - y(a-1) / x(a+1) - x(a-1)
    ! 

    !! Spacing between grid points
    dx = 2. * box / nc

    if (cur_dimension == 1) then
        do i = 1, nc_node_dim
            do j = 1, nc_node_dim
                do k = 1, nc_node_dim
                    veldivg(i, j, k) = ( (velden3(i,j+1,k)-velden3(i,j-1,k)) - (velden(i,j,k+1)-velden(i,j,k-1)) ) / dx
                enddo
            enddo
        enddo
    else if (cur_dimension == 2) then
        do i = 1, nc_node_dim
            do j = 1, nc_node_dim
                do k = 1, nc_node_dim
                    veldivg(i, j, k) = ( (velden3(i,j,k+1)-velden3(i,j,k-1)) - (velden(i+1,j,k)-velden(i-1, j, k)) ) / dx
                enddo
            enddo
        enddo
    else if (cur_dimension == 3) then
        do i = 1, nc_node_dim
            do j = 1, nc_node_dim
                do k = 1, nc_node_dim
                    veldivg(i, j, k) = ( (velden3(i+1,j,k)-velden3(i-1,j,k)) - (velden(i,j+1,k)-velden(i,j-1,k)) ) / dx
                enddo
            enddo
        enddo
    endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished velocity_curl ... elapsed time = ", time2-time1

    return

end subroutine velocity_curl
#endif
! -------------------------------------------------------------------------------------------------------

  subroutine powerspectrum(delta, delta2, pk, command)
    implicit none
    real, dimension(3, nc)       :: pk
#ifdef SLAB
    real, dimension(nc+2,nc,nc_slab) :: delta, delta2
#else
    real, dimension(nc, nc_node_dim, nc_pen+2) :: delta, delta2
#endif
    integer :: command

    integer :: i, j, k, kg, ig, mg, jg
    integer :: k1, k2
    real    :: kr, kx, ky, kz, w1, w2, pow, x, y, z, sync_x, sync_y, sync_z, kernel
    real    :: kxx, kyy, kzz, krr
#ifdef SLAB
    real, dimension(3,nc,nc_slab) :: pkt
#else
    real, dimension(3,nc,nc_pen+2) :: pkt
#endif
    real, dimension(3, nc) :: pktsum

    real(8), dimension(nc) :: kcen, kcount
    real(8), dimension(nc) :: kcen2
    real(8), dimension(nc) :: kcensum, kcountsum
    real(8), dimension(nc) :: kcensum2
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
    kcen2(:)   = 0.
    kcensum2(:) = 0.

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

                    kxx = sin(2.*pi*kx/ncr)
                    kyy = sin(2.*pi*ky/ncr)
                    kzz = sin(2.*pi*kz/ncr)
                    krr = sqrt(kxx**2+kyy**2+kzz**2)

#ifdef LOGBIN
                    kcen(k1) = kcen(k1) + w1 * log10(kr)
                    kcen(k2) = kcen(k2) + w2 * log10(kr)
                    kcen2(k1) = kcen2(k1) + w1 * log10(krr)
                    kcen2(k2) = kcen2(k2) + w2 * log10(krr)
#else
                    kcen(k1) = kcen(k1) + w1 * kr
                    kcen(k2) = kcen(k2) + w2 * kr
                    kcen2(k1) = kcen2(k1) + w1 * krr
                    kcen2(k2) = kcen2(k2) + w2 * krr
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
    call mpi_reduce(kcen2, kcensum2, nc, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
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

                !! Divide by k^2 for divergence 
                if (command == 1) then

#ifdef LOGBIN
                    kavg = real(10**(kcensum2(k) / kcountsum(k)), kind=4)
#else
                    kavg = real(kcensum2(k) / kcountsum(k), kind=4)
#endif
                    kavg = kavg / box * ncr
                    pk(1:2, k) = pk(1:2, k) / kavg**2 

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

#ifdef write_vel
subroutine writevelocityfield(command, glook)

    implicit none

    integer :: m, i, j, k, fstat
    character(len=180) :: fn
    character(len=7)   :: z_write
    character(len=4)   :: rank_string
    character(len=1)   :: dim_string
    character(len=2)   :: g_string
    real :: vsim2phys, zcur

    integer :: command
    integer(1) :: glook

    !
    ! Determine conversion to proper velocity [km/s]
    !

    if (rank == 0)  zcur = z_checkpoint(cur_checkpoint)
    call mpi_bcast(zcur, 1, mpi_real, 0, mpi_comm_world, ierr)

    vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

    !
    ! Checkpoint and rank strings
    !

    write(z_write, '(f7.3)') zcur
    z_write = adjustl(z_write)

    write(rank_string, '(i4)') rank
    rank_string = adjustl(rank_string)

    !
    ! Write out velocity field for given dimension
    !

    m = cur_dimension 

    if (m == 1) dim_string = "x"
    if (m == 2) dim_string = "y"
    if (m == 3) dim_string = "z"
    if (glook == g0) g_string = "g0"
    if (glook == g1) g_string = "g1"

    if (command == 0) then
        fn = output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//rank_string(1:len_trim(rank_string))//"_nu"//g_string//".bin"
    else if (command == 1) then
        fn = output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//rank_string(1:len_trim(rank_string))//g_string//".bin"
    else if (command == 2) then
        fn = output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//rank_string(1:len_trim(rank_string))//"_halo"//g_string//".bin"
    else if (command == 3) then
        fn = output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//rank_string(1:len_trim(rank_string))//"_nu-dm"//g_string//".bin"
    else if (command == 4) then
        fn = output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//rank_string(1:len_trim(rank_string))//"_nu-ha"//g_string//".bin"
    else if (command == 5) then
        fn = output_path//'/node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//&
             "vel"//dim_string//rank_string(1:len_trim(rank_string))//"_dm-ha"//g_string//".bin"
    endif

    open(unit=11, file=fn, status="replace", iostat=fstat, access="stream") 

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim

            write(11) velden(1:nc_node_dim, j, k) * vsim2phys

        enddo
    enddo

    close(11)

    return

end subroutine writevelocityfield
#endif

! -------------------------------------------------------------------------------------------------------

subroutine order_xvp_groups(command)
  !
  ! Split the xvp array in 2 :
  ! Group 0 at the beginning of the array
  ! Group 1 at the end of the array
  !

  implicit none

  integer(4) :: command
  integer :: np0, np1, pp
  integer(8) :: np1_tot, np2_tot
  real(4) :: r
  integer(1) :: g
  real :: xvtemp(6)

  real(8) time1, time2
  time1 = mpi_wtime(ierr)

  ! --------------------------------------------------------    
  ! Assign groups + put group1 particles at the end of the xvp array
  ! --------------------------------------------------------    

  np0 = 0
  np1 = 0

  pp = 1

  if (command == 0) then
     do
        if (np0 + np1 >= np_local) exit
        call random_number(r)
        g = int(r+0.5)

        if (g == g0) then
           np0 = np0 + 1
           pp = pp + 1
        else if (g == g1) then
           xvtemp(:) = xvp(:,np_local - np1)
           xvp(:,max_np-np1) = xvp(:,pp)
           xvp(:,pp) = xvtemp(:)
           np1 = np1 + 1
        endif
     enddo
  else if (command == 1) then
     do
        if (np0 + np1 >= np_local_dm) exit
        call random_number(r)
        g = int(r+0.5)

        if (g == g0) then
           np0 = np0 + 1
           pp = pp + 1
        else if (g == g1) then
           xvtemp(:) = xvp_dm(:,np_local_dm - np1)
           xvp_dm(:,max_np_dm-np1) = xvp_dm(:,pp)
           xvp_dm(:,pp) = xvtemp(:)
           np1 = np1 + 1
        endif
     enddo
  else if (command == 2) then
     do
        if (np0 + np1 >= np_local_h) exit
        call random_number(r)
        g = int(r+0.5)

        if (g == g0) then
           np0 = np0 + 1
           pp = pp + 1
        else if (g == g1) then
           xvtemp(:) = xvmp_h(:,np_local_h - np1)
           xvmp_h(:,max_np_h-np1) = xvmp_h(:,pp)
           xvmp_h(:,pp) = xvtemp(:)
           np1 = np1 + 1
        endif
     enddo
  endif

  ! --------------------------------------------------------    
  ! Assign np_groups in consequence
  ! --------------------------------------------------------    

  if (command == 0) then
     np_groups(0) = np0
     np_groups(1) = np1
  else if (command == 1) then
     np_groups_dm(0) = np0
     np_groups_dm(1) = np1
  else if (command == 2) then
     np_groups_h(0) = np0
     np_groups_h(1) = np1
  endif

  call mpi_allreduce(int(np0,kind=8), np1_tot, 1, mpi_integer8, mpi_sum, mpi_comm_world, ierr)
  call mpi_allreduce(int(np1,kind=8), np2_tot, 1, mpi_integer8, mpi_sum, mpi_comm_world, ierr)

  if (command == 0) then
     np_groups_tot(0) = np1_tot
     np_groups_tot(1) = np2_tot
  else if (command == 1) then
     np_groups_tot_dm(0) = np1_tot
     np_groups_tot_dm(1) = np2_tot
  else if (command == 2) then
     np_groups_tot_h(0) = np1_tot
     np_groups_tot_h(1) = np2_tot
  endif

  time2 = mpi_wtime(ierr)
  if (rank == 0) then
     write(*,*) 'np_groups_tot(0) = ', np1_tot
     write(*,*) 'np_groups_tot(1) = ', np2_tot
     write(*,*) 'Total = ', np1_tot+np2_tot
     write(*,*) "Finished order_xvp_groups ... elapsed time = ", time2-time1 
  endif

end subroutine order_xvp_groups

! -------------------------------------------------------------------------------------------------------

subroutine order_xvp_ll(command, glook)
  !
  ! Ordre the xvp array acoording to the particles's coarse grid cell appartenance
  ! Treat both groups separately
  !

  implicit none

  integer    :: k,j,i,k2,j2,i2,pp,pp2,base,command,pp_up,pp_down
  integer(8) :: index, index2
  integer(1) :: glook
  real       :: xvtemp(6), xvswap(6)
  logical    :: found,equal

  real(8) :: time1, time2 
  time1 = mpi_wtime(ierr)

  ! --------------------------------------------------------
  ! Make appropriate hoc array 
  ! --------------------------------------------------------

  call make_hoc(command, glook)

  ! --------------------------------------------------------
  ! Determine particle indices we need to look at
  ! --------------------------------------------------------

  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0) + np_buf_groups_h(g0)
     endif

  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
  endif

  ! --------------------------------------------------------
  ! Sort the xvp_array
  ! --------------------------------------------------------
    
  if (command == 0) then

     !! Start at the beginning of xvp_array, looking at the first coarse cell
     pp = pp_down
     k2 = hoc_nc_l
     j2 = hoc_nc_l
     i2 = hoc_nc_l
     base = hoc_nc_h - hoc_nc_l + 1
     index2 = 0
     xvtemp(:) = xvp(:,pp)

     do
        !! if we are looking at the last coarse cell, it means that everything is already sorted
        if (k2 == hoc_nc_h .and. j2 == hoc_nc_h .and. i2 == hoc_nc_h) exit

        !! else we want to fill the coarse cell we are looking at
        do
           !! if the coarse cell is filled, jump to next coarse cell
           if (pp > hoc(i2, j2, k2)) exit

           !! coarse cell the current particle should be
           i = floor(xvtemp(1)/mesh_scale) + 1
           j = floor(xvtemp(2)/mesh_scale) + 1
           k = floor(xvtemp(3)/mesh_scale) + 1

           index = (i-hoc_nc_l) + (j-hoc_nc_l) * base + (k-hoc_nc_l)*base*base
           !! if the particle is already sorted
           !! jump to next particle
           if ( index <= index2 ) then
              pp = pp + 1
              xvtemp(:) = xvp(:,pp)
           !! else we put it in its right cell
           else
              !! determines its position
              pp2 = hoc(i,j,k)
              !! swap with the particle at its position
              xvswap(:) = xvp(:,pp2)
              xvp(:,pp) = xvswap(:)
              xvp(:,pp2) = xvtemp(:)
              !! new particle we are looking at
              xvtemp(:) = xvswap(:)
              !! update the position for the cell
              hoc(i,j,k) = hoc(i,j,k) - 1
           endif

        enddo

        !! now we have pp = hoc(i2,j2,k2) + 1,
        !! which means that the cell is completed,
        !! we just have to jump to the next cell
        if ( j2 == hoc_nc_h .and. i2 == hoc_nc_h ) then
           k2 = k2 + 1
           j2 = hoc_nc_l
           i2 = hoc_nc_l
        else if (i2 == hoc_nc_h ) then
           j2 = j2 + 1
           i2 = hoc_nc_l
        else
           i2 = i2 + 1
        endif
        index2 = (i2-hoc_nc_l) + (j2-hoc_nc_l) * base + (k2-hoc_nc_l)*base*base

     enddo

  else if (command == 1) then

     pp = pp_down
     k2 = hoc_nc_l
     j2 = hoc_nc_l
     i2 = hoc_nc_l
     base = hoc_nc_h - hoc_nc_l + 1
     index2 = 0
     xvtemp(:) = xvp_dm(:,pp)

     do
        if (k2 == hoc_nc_h .and. j2 == hoc_nc_h .and. i2 == hoc_nc_h) exit

        do
           if (pp > hoc_dm(i2, j2, k2)) exit

           i = floor(xvtemp(1)/mesh_scale) + 1
           j = floor(xvtemp(2)/mesh_scale) + 1
           k = floor(xvtemp(3)/mesh_scale) + 1

           index = (i-hoc_nc_l) + (j-hoc_nc_l) * base + (k-hoc_nc_l)*base*base
           if ( index <= index2 ) then
              pp = pp + 1
              xvtemp(:) = xvp_dm(:,pp)
           else
              pp2 = hoc_dm(i,j,k)
              xvswap(:) = xvp_dm(:,pp2)
              xvp_dm(:,pp) = xvswap(:)
              xvp_dm(:,pp2) = xvtemp(:)
              xvtemp(:) = xvswap(:)
              hoc_dm(i,j,k) = hoc_dm(i,j,k) - 1
           endif

        enddo
        if ( j2 == hoc_nc_h .and. i2 == hoc_nc_h ) then
           k2 = k2 + 1
           j2 = hoc_nc_l
           i2 = hoc_nc_l
        else if (i2 == hoc_nc_h ) then
           j2 = j2 + 1
           i2 = hoc_nc_l
        else
           i2 = i2 + 1
        endif
        index2 = (i2-hoc_nc_l) + (j2-hoc_nc_l) * base + (k2-hoc_nc_l)*base*base
     enddo

  else if (command == 2) then

     pp = pp_down
     k2 = hoc_nc_l_h
     j2 = hoc_nc_l_h
     i2 = hoc_nc_l_h
     base = hoc_nc_h_h - hoc_nc_l_h + 1
     index2 = 0
     xvtemp(:) = xvmp_h(:,pp)

     do
        if (k2 == hoc_nc_h_h .and. j2 == hoc_nc_h_h .and. i2 == hoc_nc_h_h) exit

        do
           if (pp > hoc_h(i2, j2, k2)) exit

           i = floor(xvtemp(1)/mesh_scale_h) + 1
           j = floor(xvtemp(2)/mesh_scale_h) + 1
           k = floor(xvtemp(3)/mesh_scale_h) + 1

           index = (i-hoc_nc_l_h) + (j-hoc_nc_l_h) * base + (k-hoc_nc_l_h)*base*base
           if ( index <= index2 ) then
              pp = pp + 1
              xvtemp(:) = xvmp_h(:,pp)
           else
              pp2 = hoc_h(i,j,k)
              xvswap(:) = xvmp_h(:,pp2)
              xvmp_h(:,pp) = xvswap(:)
              xvmp_h(:,pp2) = xvtemp(:)
              xvtemp(:) = xvswap(:)
              hoc_h(i,j,k) = hoc_h(i,j,k) - 1
           endif

        enddo
        if ( j2 == hoc_nc_h_h .and. i2 == hoc_nc_h_h ) then
           k2 = k2 + 1
           j2 = hoc_nc_l_h
           i2 = hoc_nc_l_h
        else if (i2 == hoc_nc_h_h ) then
           j2 = j2 + 1
           i2 = hoc_nc_l_h
        else
           i2 = i2 + 1
        endif
        index2 = (i2-hoc_nc_l_h) + (j2-hoc_nc_l_h) * base + (k2-hoc_nc_l_h)*base*base
     enddo

  endif

  time2 = mpi_wtime(ierr)
  if (rank == 0) then
     write(*,*) 'Number of particles in order_xvp_ll: ', (pp_up-pp_down+1)
     write(*,*) "Finished order_xvp_ll ... time elapsed = ", time2-time1
  endif

end subroutine order_xvp_ll

! -------------------------------------------------------------------------------------------------------

subroutine make_hoc(command, glook)
  !
  ! Reconstructs hoc just before velocity density (this way hoc can be equivalenced)
  !

  implicit none
  integer      :: command,pp,pp_down,pp_up,k,j,i
  real         :: xvtemp(6)
  integer(1)   :: glook

  ! --------------------------------------------------------
  ! Determine xvp spot we have to look at
  ! --------------------------------------------------------

  if (glook == g0) then
     pp_down = 1
     if (command == 0) then
        pp_up = np_groups(g0) + np_buf_groups(g0)
     else if (command == 1) then
        pp_up = np_groups_dm(g0) + np_buf_groups_dm(g0)
     else if (command == 2) then
        pp_up = np_groups_h(g0) + np_buf_groups_h(g0)
     endif

  else if (glook == g1) then
     if (command == 0) then
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif
  endif

  ! --------------------------------------------------------
  ! Actually constructs hoc
  ! --------------------------------------------------------

  if (command == 0) then !! NEUTRINOS

     !! Clear hoc
     do k = hoc_nc_l, hoc_nc_h
        hoc(:,:,k) = 0
     enddo

     !! Count particles in each coarse cell                                                                                                     
     do pp = pp_down,pp_up
        xvtemp(1:3) = xvp(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale) + 1
        j = floor(xvtemp(2)/mesh_scale) + 1
        k = floor(xvtemp(3)/mesh_scale) + 1
        hoc(i,j,k) = hoc(i,j,k) + 1
     enddo

     !! cum sum hoc --> hoc will point to the end of memory spot reserved for each coarse cell                                                
     pp = pp_down - 1
     do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
           do i = hoc_nc_l, hoc_nc_h
              pp = hoc(i,j,k) + pp
              hoc(i,j,k) = pp
           enddo
        enddo
     enddo

  else if (command == 1) then !! DARK MATTER

     !! Clear hoc
     do k = hoc_nc_l, hoc_nc_h
        hoc_dm(:,:,k) = 0
     enddo

     !! Count particles in each coarse cell
     do pp = pp_down,pp_up
        xvtemp(1:3) = xvp_dm(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale) + 1
        j = floor(xvtemp(2)/mesh_scale) + 1
        k = floor(xvtemp(3)/mesh_scale) + 1
        hoc_dm(i,j,k) = hoc_dm(i,j,k) + 1
     enddo

     !! cum sum hoc --> hoc will point to the end of memory spot reserved for each coarse cell
     pp = pp_down - 1
     do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
           do i = hoc_nc_l, hoc_nc_h
              pp = hoc_dm(i,j,k) + pp
              hoc_dm(i,j,k) = pp
           enddo
        enddo
     enddo

  else if (command == 2) then !! HALOS

     !! Clear hoc
     do k = hoc_nc_l_h, hoc_nc_h_h
        hoc_h(:,:,k) = 0
     enddo

     !! Count particles in each coarse cell
     do pp = pp_down,pp_up
        xvtemp(1:3) = xvmp_h(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale_h) + 1
        j = floor(xvtemp(2)/mesh_scale_h) + 1
        k = floor(xvtemp(3)/mesh_scale_h) + 1
        hoc_h(i,j,k) = hoc_h(i,j,k) + 1
     enddo

     !! cum sum hoc --> hoc will point to the end of memory spot reserved for each coarse cell
     pp = pp_down - 1
     do k = hoc_nc_l_h, hoc_nc_h_h
        do j = hoc_nc_l_h, hoc_nc_h_h
           do i = hoc_nc_l_h, hoc_nc_h_h
              pp = hoc_h(i,j,k) + pp
              hoc_h(i,j,k) = pp
           enddo
        enddo
     enddo

  endif

  return

end subroutine make_hoc

! -------------------------------------------------------------------------------------------------------

subroutine swap_slab12(command)
    !
    ! Swaps the data between slab and slab2 depending on command.
    !

    implicit none
    
    integer :: command

    if (command == 0) then !! Place slab into slab2

        !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
        do kt = 1, kt_stop
            slab2(:,:,kt) = slab(:,:,kt)
        enddo
        !$omp end parallel do

    else if (command == 1) then !! Place slab2 into slab

        !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
        do kt = 1, kt_stop
            slab(:,:,kt) = slab2(:,:,kt)
        enddo
        !$omp end parallel do

    endif
end subroutine swap_slab12

! -------------------------------------------------------------------------------------------------------

subroutine swap_velden12(command)
    !
    ! Swaps data between velden and velden2 depending on command.
    !

    implicit none

    integer :: command

    if (command == 0) then !! Place velden into velden2

        !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
        do kt = 0, nc_node_dim+1
            velden2(:,:,kt) = velden(:,:,kt)
        enddo
        !$omp end parallel do

    else if (command == 1) then !! Place velden2 into velden

        !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
        do kt = 0, nc_node_dim+1
            velden(:,:,kt) = velden2(:,:,kt)
        enddo
        !$omp end parallel do

    endif

end subroutine swap_velden12

! -------------------------------------------------------------------------------------------------------
#ifdef CURL
subroutine swap_velden13(command)
    !
    ! Swaps data between velden and velden3 depending on command.
    !

    implicit none

    integer :: command

    if (command == 0) then !! Place velden into velden3

        !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
        do kt = 0, nc_node_dim+1
            velden3(:,:,kt) = velden(:,:,kt)
        enddo
        !$omp end parallel do

    else if (command == 1) then !! Place velden3 into velden

        !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
        do kt = 0, nc_node_dim+1
            velden(:,:,kt) = velden3(:,:,kt)
        enddo
        !$omp end parallel do

    endif

end subroutine swap_velden13
#endif
! -------------------------------------------------------------------------------------------------------

subroutine relative_velocity
    !
    ! Computes the difference velden-velden2 and stores the result in velden.
    !

    implicit none

    !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
    do kt = 0, nc_node_dim+1
        velden(:,:,kt) = velden(:,:,kt) - velden2(:,:,kt)
    enddo
    !$omp end parallel do

end subroutine relative_velocity

! -------------------------------------------------------------------------------------------------------

end program cic_crossvel

