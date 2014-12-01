!!
!! cic_velpower_dmnu.f90
!!
!! Program to compute velocity power spectra and linear velocity power spectra
!! for dark matter, neutrinos and dark matter halos.
!!  
!! This is to be used in conjunction with particle checkpoint files
!! produced by a -DNEUTRINOS cubep3m simulation.
!!
!! * Using FFTW on the SciNet GPC compile with:
!!   mpif90 -shared-intel -fpp -g -O3 -openmp -mcmodel=medium -mt_mpi indexedsort.f90 cic_velpower_dmnu.f90 -I$SCINET_FFTW_INC 
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
!!           => THIS DOES NOT WORK
!!   -DKAISER: Adjusts for redshift space distortions.
!!   -DDEBUG: Output useful debugging information.
!!   -DDEBUG_LOW: idem
!!   -DLOGBIN
!!   -DNPLINLOGMAX: Uses the algotithm in indexedsort.f90 to search for the nearest particle
!!   -DNP_FAST: Uses an array to speed up the algorithm searching for the nearest particle
!!   -DCATCH_SEED: Catches the random number generator seed from file
!!   -DWRITE_SEED: Writes the seed to a file
!!   -DRADIUS: save the distance to the nearest particle instead of the velocity
!!             useful to make a plot... do not forget to activate write_vel
!!   -Dwrite_vel: Writes gridded velocity fields to binary files. 
!!   -Dwrite_den: Writes gridded density fields to binary files. 
!!   -DUCCFD: Uses Coarse Cells For Density (see subroutine densityfield)
!!   -Dcompute_denPS: Computes density power spectra
 
program cic_velpower_halo
  use omp_lib

  implicit none

  include 'mpif.h'
  include '../../parameters'

  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints'

  !! Threading
  integer(4), parameter :: nt = 8

  !! Number of nearest particles from grid centre to determine average velocity
  integer, parameter :: N_closest_nu = 1
  integer, parameter :: N_closest_dm = 1 
  integer, parameter :: N_closest_h = 1
  integer, parameter :: N_closest_auto_nu = 1
  integer, parameter :: N_closest_auto_dm = 1
  integer, parameter :: N_closest_auto_h = 1
  integer, parameter :: max_N_closest = 1

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc
  real, parameter    :: npr=np

#ifdef LOGBIN
  integer, parameter :: numbins = 16 !hc / 48 
#endif

  !! checkpoints parameters
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

  !! For storage of dark matter particles and halos (usually small since ratio_nudm_dim generally > 1)
  integer(4), parameter :: max_np_dm = max_np / ratio_nudm_dim**3
  integer(4), parameter :: max_np_h = 1000000

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local,np_local_dm,np_local_h,np_h
  integer(8) :: plan, iplan
  logical :: firstfftw

  !! have velocity power spectra for each x, y, z
  integer cur_dimension

  !! Other parameters
  real, parameter :: pi=3.14159
  
  !! Particles arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(6,max_np_dm) :: xvp_dm
  real, dimension(6,max_np_h)  :: xvmp_h
  real, dimension(6,np_buffer) :: xp_buf
  real, dimension(6*np_buffer) :: send_buf, recv_buf

  !! Power spectrum arrays
  real, dimension(3, 3, nc, 0:34) :: pk_all 
  real, dimension(3, nc) :: pkdm
  logical, dimension(0:34) :: computePS

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

  !! Parameters for the coarse grid
  integer(4), parameter :: nfine_buf = 16
  integer(4), parameter :: mesh_scale = 4
  integer(4), parameter :: nc_buf = nfine_buf / mesh_scale
  integer(4), parameter :: nm_node_dim = nc_node_dim / mesh_scale
  integer(4), parameter :: hoc_nc_l = 1 - nc_buf
  integer(4), parameter :: hoc_nc_h = nm_node_dim + nc_buf
  integer(4), parameter :: hoc_pass_depth = 2*nc_buf
  real(4), parameter    :: rnf_buf = real(nfine_buf)
  integer(4), parameter :: num_ngbhs = (2*nc_buf+1)**3

  integer(4), parameter :: nfine_buf_h = 256
  integer(4), parameter :: mesh_scale_h = 32
  integer(4), parameter :: nc_buf_h = nfine_buf_h / mesh_scale_h
  integer(4), parameter :: nm_node_dim_h = nc_node_dim / mesh_scale_h
  integer(4), parameter :: hoc_nc_l_h = 1 - nc_buf_h
  integer(4), parameter :: hoc_nc_h_h = nm_node_dim_h + nc_buf_h
  integer(4), parameter :: hoc_pass_depth_h = 2*nc_buf_h
  real(4), parameter    :: rnf_buf_h = real(nfine_buf_h)
  integer(4), parameter :: num_ngbhs_h = (2*nc_buf_h+1)**3

  !! To order cell_search
  integer(4) :: cell_search_work(3, num_ngbhs)
  integer(4) :: cell_search_work_h(3, num_ngbhs_h)
#ifdef NP_LINLOGMAX
  integer(4) :: cell_search_r_max_work(num_ngbhs)
  integer(4) :: cell_search_r_max_work_h(num_ngbhs_h)
#endif
  integer(4) :: sorted_indexes(num_ngbhs)
  integer(4) :: sorted_indexes_h(num_ngbhs_h)

  !! hoc array (see remake_hoc)
  integer(4) :: hoc(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h)
  integer(4) :: hoc_dm(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h)
  integer(4) :: hoc_h(hoc_nc_l_h:hoc_nc_h_h, hoc_nc_l_h:hoc_nc_h_h, hoc_nc_l_h:hoc_nc_h_h)

  !! Array to search the nearest particle
#ifdef NP_FAST
  integer(4), parameter :: max_npart_cell_search = int(real(max_np)/real(nm_node_dim)**2)
  integer(4), parameter :: max_npart_cell_search_h = int(real(max_np_h)/real(nm_node_dim)**2)
  real(4) :: xvp_cell(6, max_npart_cell_search, nt)
  real(4) :: xvp_cell_h(6, max_npart_cell_search_h, nt)
#endif

#ifdef NP_LINLOGMAX
  integer(4), parameter :: max_npart_search = int(num_ngbhs*real(max_np)/real(nc_node_dim)**2)
  integer(4), parameter :: max_npart_search_h = int(num_ngbhs_h*real(max_np_h)/real(nc_node_dim)**2)

  integer(4) :: ipos(max_npart_search, nt)
  real(4)    :: rpos(2, max_npart_search, nt)
  integer(4) :: ipos_h(max_npart_search_h, nt)
  real(4)    :: rpos_h(2, max_npart_search_h, nt)

  real(4)    :: cell_search_r_max(num_ngbhs)
  real(4)    :: cell_search_r_max_h(num_ngbhs_h)
#else
  real(4)    :: rpos(2, max_N_closest, nt)
  real(4)    :: rpos_h(2, max_N_closest, nt)
#endif
  integer(4) :: cell_search(3, num_ngbhs)
  real(4)    :: cell_search_r(num_ngbhs)

  integer(4) :: cell_search_h(3, num_ngbhs_h)
  real(4)    :: cell_search_r_h(num_ngbhs_h)

  !! Number of particles in each group / in the buffer arrea
  integer(4) :: np_groups(0:1), np_groups_dm(0:1), np_groups_h(0:1)
  integer(4) :: np_buf_groups(0:1), np_buf_groups_dm(0:1), np_buf_groups_h(0:1)
  integer(4) :: np_groups_tot(0:1), np_groups_tot_dm(0:1), np_groups_tot_h(0:1)

  integer(1), parameter :: g0 = 0
  integer(1), parameter :: g1 = 1

  !! Random number generator parameters
#ifdef CATCH_SEED
  logical, parameter :: generate_seed = .false. !! determine if the program should catch the seed from a file
#else
  logical, parameter :: generate_seed = .true. !! determine if the program should catch the seed from a file
#endif

#ifdef WRITE_SEED
  logical, parameter :: write_seed = .true. !! determine if the program should write the seed into a file                                   
#else
  logical, parameter :: write_seed = .false. !! determine if the program should write the seed into a file                                   
#endif

  integer(4), allocatable, dimension(:) :: randomseed
  integer :: randomseedsize

  !! Array storing the velocity fields
  !! (only work with one component of the velocity field at a time) 
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

  !! For timing
  real(8) time1, time2, timeStart, timeCheckpoint

!-------------------------------------------------------------------------------------------------------
! Equivalence statement. Place here the code returned by the python script. 
!-------------------------------------------------------------------------------------------------------

equivalence (xp_buf, slab, cube, hoc, hoc_dm, hoc_h)
equivalence (send_buf, slab3)
equivalence (recv_buf, momden, recv_cube)
equivalence (cell_search_work_h, cell_search_work)
equivalence (sorted_indexes_h, sorted_indexes)

common xvp, recv_buf, send_buf, xp_buf, slab2, xvp_dm, xvmp_h, momden_recv_buff, momden_send_buff, pk_all, cell_search_h, cell_search_work_h, cell_search_r_h, sorted_indexes_h, cell_search, cell_search_r, rpos_h, rpos

!-------------------------------------------------------------------------------------------------------
! MAIN
!-------------------------------------------------------------------------------------------------------

  call mpi_initialize
  call omp_set_num_threads(nt)
  call initialize_random_number

  timeStart = mpi_wtime(ierr)

  if (rank == 0) call writeparams

  firstfftw = .true.  ! initialize fftw so that it generates the plans

  call read_checkpoint_list

  do cur_checkpoint = 1, num_checkpoints

     if (rank == 0) then
        timeCheckpoint = mpi_wtime(ierr)
        write(*,*) '****************************************************'
        write(*,*) 'TIME = ', timeCheckpoint - timeStart
        write(*,*) '****************************************************'
        write(*,*) 'STARTING CHECKPOINT : ',cur_checkpoint
        write(*,*) 'z = ', z_checkpoint(cur_checkpoint)
        write(*,*) '****************************************************'
     endif

     !--------------------------------------------------------------------------------                                                            
     ! INPUT
     !--------------------------------------------------------------------------------

     computePS(0)  = .true.  ! nu
     computePS(1)  = .true.  ! dm
     computePS(2)  = .true.  ! dm x nu
     computePS(3)  = .true.  ! rel
     computePS(4)  = .true.  ! halo
     computePS(5)  = .true.  ! nulin grad
     computePS(6)  = .true.  ! dmlin grad
     computePS(7)  = .true.  ! halolin grad
     computePS(8)  = .true.  ! nuh grad
     computePS(9)  = .true.  ! nudm grad
     computePS(10) = .true.  ! dm x halo
     computePS(11) = .true.  ! halo x halo lin
     computePS(12) = .true.  ! dm x halo lin
     computePS(13) = .true.  ! dm x dm lin
     computePS(14) = .true.  ! nu x nuh
     computePS(15) = .true.  ! nu x nudm
     computePS(16) = .true.  ! relh grad
     computePS(17) = .true.  ! reldm grad
     computePS(18) = .true.  ! rel x relh
     computePS(19) = .true.  ! rel x reldm
     computePS(20) = .true.  ! rel2
     computePS(21) = .true.  ! rel2 x relh
     computePS(22) = .true.  ! halo x nu
     computePS(23) = .true.  ! nulin vtf
     computePS(24) = .true.  ! dmlin vtf
     computePS(25) = .true.  ! dmlin x nulin vtf
     computePS(26) = .true.  ! halolin vtf
     computePS(27) = .true.  ! dmlin x halolin vtf
     computePS(28) = .true.  ! nudm vtf
     computePS(29) = .true.  ! reldm vtf
     computePS(30) = .true.  ! nuh vtf
     computePS(31) = .true.  ! relh vtf 
     computePS(32) = .true.  ! halolin x nu
     computePS(33) = .true.  ! dm x nuh
     computePS(34) = .true.  ! halolin x nuh grad

     !--------------------------------------------------------------------------------                                                            
     ! Initialize variables and read particles files                                                                                      
     !--------------------------------------------------------------------------------
     
     call initvar
     call init_cell_search(0)
     call init_cell_search(2)

     if (rank == 0) then
        timeCheckpoint = mpi_wtime(ierr)
        write(*,*) '****************************************************'
        write(*,*) 'TIME = ', timeCheckpoint - timeStart
        write(*,*) '****************************************************'
        write(*,*) 'STARTING READING PARTICLES ...'
        write(*,*) '****************************************************'
     endif

     call read_particles_files

     !-------------------------------------------------------------------                                                                            
     ! Auto power spectra loop                                                                                                          
     !-------------------------------------------------------------------
     ! Compute auto-spectra by dividing the particles into two groups and
     ! computing the cross spectrum of the two groups.
     ! This is done to remove shot noise which does not correlate 
     ! between the two groups.
     !
     ! Compute one dimension at a time.
     !-------------------------------------------------------------------                                                                           
     
     if (rank == 0) then
        timeCheckpoint = mpi_wtime(ierr)
        write(*,*) '****************************************************'
        write(*,*) 'TIME = ', timeCheckpoint - timeStart
        write(*,*) '****************************************************'
        write(*,*) 'ENTERING AUTO-POWER SPECTRA LOOP ...'
        write(*,*) '****************************************************'
     endif

     call emptyPSArray
     do cur_dimension = 1, 3
        call auto_power_loop(cur_dimension)
     enddo

     !-------------------------------------------------------------------                                                                            
     ! VTF Auto power spectra loop                                                                                                          
     !-------------------------------------------------------------------                                                                           
     ! Compute linear velocity spectra without the gradient method
     ! (this does not produce a velocity field)
     !-------------------------------------------------------------------

     if (rank == 0) then
        timeCheckpoint = mpi_wtime(ierr)
        write(*,*) '****************************************************'
        write(*,*) 'TIME = ', timeCheckpoint - timeStart
        write(*,*) '****************************************************'
        write(*,*) 'ENTERING VTF AUTO-POWER SPECTRA LOOP ...'
        write(*,*) '****************************************************'
     endif

     call emptyPSArray
     call auto_power_vtf

     !-------------------------------------------------------------------                                                                            
     ! Put all particles into group 0 for cross-power                                                                         
     !-------------------------------------------------------------------                                                                            

     call clear_groups
     call order_xvp_ll(0, g0, .false.)
     call order_xvp_ll(1, g0, .false.)
     call order_xvp_ll(2, g0, .false.)

#ifdef compute_denPS
     call emptyPSArray
     call cross_power_density
#endif

     !-------------------------------------------------------------------                                                                            
     ! Cross power spectra loop                                                                                                          
     !-------------------------------------------------------------------                                                                           
     
     if (rank == 0) then
        timeCheckpoint = mpi_wtime(ierr)
        write(*,*) '****************************************************'
        write(*,*) 'TIME = ', timeCheckpoint - timeStart
        write(*,*) '****************************************************'
        write(*,*) 'ENTERING CROSS-POWER SPECTRA LOOP ...'
        write(*,*) '****************************************************'
     endif

     call emptyPSArray
     do cur_dimension = 1, 3
        call cross_power_loop(cur_dimension)
     enddo
  enddo

  if (rank == 0) then
     timeCheckpoint = mpi_wtime(ierr)
     write(*,*) '****************************************************'
     write(*,*) 'TIME = ', timeCheckpoint - timeStart
     write(*,*) '****************************************************'
  endif

  call cp_fftw(0)
  call mpi_finalize(ierr)

!--------------------------------------------------------------------------------------------------------
! SUBROUTINES
!--------------------------------------------------------------------------------------------------------

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

  integer :: k
  
  !! Particle positions and velocities
  do k = 1, max_np
     xvp(:, k) = 0.0
  enddo
  do k = 1, max_np_dm
     xvp_dm(:,k) = 0.0
  enddo
  do k = 1, max_np_h
     xvmp_h(:,k) = 0.0
  enddo
  
  !! Momentum and matter density arrays
  do k = 0, nc_node_dim + 1
     momden_send_buff(:, k) = 0.
     momden_recv_buff(:, k) = 0.
     momden(:, :, k) = 0.
  enddo
  
  !! Fourier transform arrays
#ifdef SLAB
  do k = 1, nc_slab
     slab_work(:, :, k)=0.0
  enddo
#endif
  
  do k = 1, nc_node_dim
     cube(:, :, k) = 0.0
  enddo
  
  do k = 1, kt_stop
     slab(:, :, k) = 0.0
     slab2(:, :, k) = 0.0
     slab3(:, :, k) = 0.0
  enddo
  
  do k = 1, np_buffer
     xp_buf(:, k) = 0.0
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
     pkdm(:, k) = 0.
  enddo
  do k = 0, 34
     pk_all(:, :, :, k) = 0.0
  enddo
  
  !! Arrays to search the nearest particle
  do k = 1, num_ngbhs
     cell_search_work(:,k) = 0
#ifdef NP_LINLIGMAX
     cell_search_r_max_work(k) = 0
#endif
     sorted_indexes(k) = 0
  enddo

  do k = 1, num_ngbhs_h
     cell_search_work_h(:,k) = 0
#ifdef NP_LINLIGMAX
     cell_search_r_max_work_h(k) = 0
#endif
     sorted_indexes_h(k) = 0
  enddo

  do k = hoc_nc_l, hoc_nc_h
     hoc(:,:,k) = 0
     hoc_dm(:,:,k) = 0
  enddo

  do k = hoc_nc_l_h, hoc_nc_h_h
     hoc_h(:,:,k) = 0
  enddo

#ifdef NP_FAST

  do k = 1, nt
     xvp_cell(:,:,k) = 0.0
     xvp_cell_h(:,:,k) = 0.0
  enddo

#endif

#ifdef NP_LINLOGMAX

  do k = 1, nt
     ipos(:,k) = 0
     rpos(:,:,k) = 0.0
     ipos_h(:,k) = 0
     rpos_h(:,:,k) = 0.0
  enddo

  do k = 1, num_ngbhs
     cell_search_r_max(k) = 0.0
  enddo
  do k = 1, num_ngbhs_h
     cell_search_r_max_h(k) = 0.0
  enddo

#else

  do k = 1, nt
     rpos(:,:,k) = 0.0
     rpos_h(:,:,k) = 0.0
  enddo

#endif

  do k = 1, num_ngbhs
     cell_search(:,k) = 0
     cell_search_r(k) = 0.0
  enddo

  do k = 1, num_ngbhs_h
     cell_search_h(:,k) = 0
     cell_search_r_h(k) = 0.0
  enddo

  np_groups(:) = 0
  np_groups_dm(:) = 0
  np_groups_h(:) = 0
  np_buf_groups(:) = 0
  np_buf_groups_dm(:) = 0
  np_buf_groups_h(:) = 0
  np_groups_tot(:) = 0
  np_groups_tot_dm(:) = 0
  np_groups_tot_h(:) = 0

  return
end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine emptyPSArray
  !
  ! Reinitialize power spectra arrays
  !
  implicit none
  integer k

  do k = 1, nc
     pkdm(:, k) = 0.
  enddo
  do k = 0, 34
     pk_all(:, :, :, k) = 0.0
  enddo

end subroutine emptyPSArray

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
  integer :: k1,j1,i1, kc,jc,ic
  real    :: dmin,dmax
  integer :: ind, ibuf

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
#ifdef NP_LINLOGMAX
                    dmax = dmin
#endif
                    do k1 = -1,1
                       do j1 = -1,1
                          do i1 = -1,1
                             do kc = -1, 1
                                do jc = -1, 1
                                   do ic = -1, 1
                                      dmin = min(dmin, sqrt(real(i+0.5*(i1-ic))**2 + real(j+0.5*(j1-jc))**2 + real(k+0.5*(k1-kc))**2))
#ifdef NP_LINLOGMAX
                                      dmax = max(dmax, sqrt(real(i+0.5*(i1-ic))**2 + real(j+0.5*(j1-jc))**2 + real(k+0.5*(k1-kc))**2))
#endif
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                    cell_search_r(ind) = dmin * mesh_scale
#ifdef NP_LINLOGMAX
                    cell_search_r_max_work(ind) = dmax * mesh_scale
#endif
                    sorted_indexes(ind) = ind
                    ind = ind + 1
                 endif
              enddo
           enddo
        enddo
     enddo

     cell_search_r(0) = 0.0
     call indexedsort(num_ngbhs, cell_search_r(:num_ngbhs), sorted_indexes(:num_ngbhs))

     do ind = 1, num_ngbhs
        cell_search(:,ind) = cell_search_work(:,sorted_indexes(ind))
#ifdef NP_LINLOGMAX
        cell_search_r_max(ind) = cell_search_r_max_work(sorted_indexes(ind))
#endif
     enddo

  else if (command == 2) then

     ind = 1
     do ibuf = 0, nc_buf_h
        do k = -ibuf, ibuf
           do j = -ibuf, ibuf
              do i = -ibuf, ibuf
                 if (.not. (abs(i) < ibuf .and. abs(j) < ibuf .and. abs(k) < ibuf)) then
                    cell_search_work_h(1, ind) = i
                    cell_search_work_h(2, ind) = j
                    cell_search_work_h(3, ind) = k
                    dmin = sqrt(real(i)**2+ real(j)**2+real(k)**2)
#ifdef NP_LINLOGMAX
                    dmax = dmin
#endif
                    do k1 = -1,1
                       do j1 = -1,1
                          do i1 = -1,1
                             do kc = -1,1
                                do jc = -1,1
                                   do ic = -1,1
                                      dmin = min(dmin, sqrt(real(i+0.5*(i1-ic))**2 + real(j+0.5*(j1-jc))**2 + real(k+0.5*(k1-kc))**2))
#ifdef NP_LINLOGMAX
                                      dmax = max(dmax, sqrt(real(i+0.5*(i1-ic))**2 + real(j+0.5*(j1-jc))**2 + real(k+0.5*(k1-kc))**2))
#endif
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                    cell_search_r_h(ind) = dmin * mesh_scale_h
#ifdef NP_LINLOGMAX
                    cell_search_r_max_work_h(ind) = dmax * mesh_scale_h
#endif
                    sorted_indexes_h(ind) = ind
                    ind = ind + 1
                 endif
              enddo
           enddo
        enddo
     enddo

     cell_search_r_h(0) = 0.0
     call indexedsort(num_ngbhs_h, cell_search_r_h(:num_ngbhs_h), sorted_indexes_h(:num_ngbhs))

     do ind = 1, num_ngbhs_h
        cell_search_h(:,ind) = cell_search_work_h(:,sorted_indexes_h(ind))
#ifdef NP_LINLOGMAX
        cell_search_r_max_h(ind) = cell_search_r_max_work_h(sorted_indexes_h(ind))
#endif
    enddo
 endif

end subroutine init_cell_search

! -------------------------------------------------------------------------------------------------------

subroutine read_particles(command)
  !
  ! Read xv file and store position and velocity in the xvp arrays
  ! 

  implicit none
    
  real z_write, np_total
  integer j, fstat
  character(len=7) :: z_string
  character(len=4) :: rank_string
  character(len=100) :: check_name
  integer(4) :: command

  !! These are unnecessary headers from the checkpoint
  real(4) :: a, t, tau, dt_f_acc, dt_c_acc, dt_pp_acc, mass_p, z_current
  integer(4) :: nts, sim_checkpoint, sim_projection, sim_halofind
  
  !! unnecessary data for halos
  real(4)                :: garbage1
  real(4), dimension(3)  :: garbage3
  real(4), dimension(4)  :: garbage4
  real(4), dimension(15) :: garbage15

  time1 = mpi_wtime(ierr)

  z_current = z_checkpoint(cur_checkpoint)

  if (rank == 0) then
     z_write = z_checkpoint(cur_checkpoint)
     write(*,*) 'Starting read_particles ...'
  endif

  call mpi_bcast(z_write, 1, mpi_real, 0, mpi_comm_world, ierr)

  !-----------------------------------------------------------
  ! Determine file name
  !-----------------------------------------------------------

  write(z_string,'(f7.3)') z_write
  z_string=adjustl(z_string)
  
  write(rank_string,'(i4)') rank
  rank_string=adjustl(rank_string)

  check_name = '.dat'

  if (command == 0) then
     check_name = '_nu'//check_name
  endif

  check_name = rank_string(1:len_trim(rank_string))//check_name
  
  if (command == 0 .or. command == 1) then
     check_name = 'xv'//check_name
  else if (command == 2) then
     check_name = 'halo'//check_name
  endif

  check_name = 'node'//rank_string(1:len_trim(rank_string))//'/'//z_string(1:len_trim(z_string))//check_name
  
  if (command == 0 .and. z_write == z_i_nu) then
     check_name = ic_path//check_name
  else if (command == 1 .and. z_write == z_i) then
     check_name = ic_path//check_name
  else
     check_name = output_path//check_name
  endif

  !-----------------------------------------------------------
  ! Open file
  !-----------------------------------------------------------
  
  open(unit=21,file=check_name,status="old",iostat=fstat,access="stream")

  if (fstat /= 0) then
     write(*,*) 'ERROR: Cannot open particle file'
     write(*,*) 'rank', rank, ' file: ',check_name
     call mpi_abort(mpi_comm_world, ierr, ierr)
  endif

  !-----------------------------------------------------------
  ! Read header data
  !-----------------------------------------------------------

  if (command == 0) then
     read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
          sim_projection,sim_halofind,mass_p
  else if (command == 1) then
     read(21) np_local_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
          sim_projection,sim_halofind,mass_p
  else
     read(21) np_local_h,t, tau
  endif

  !-----------------------------------------------------------
  ! Check for memory problems
  !-----------------------------------------------------------

  if (command == 0) then
     if (np_local > max_np) then
        write(*,*) 'ERROR: Too many neutrinos to store in memory!'
        write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else if (command == 1) then
     if (np_local_dm > max_np_dm) then
        write(*,*) 'ERROR: Too many dm particles to store in memory!'
        write(*,*) 'rank', rank, 'np_local_dm', np_local_dm, 'max_np_dm', max_np_dm
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  else
     if (np_local_h > max_np_h) then
        write(*,*) 'ERROR: Too many halos to store in memory!'
        write(*,*) 'rank', rank, 'np_local_h', np_local_h, 'max_np_h', max_np_h
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif
  endif

  !-----------------------------------------------------------
  ! Tally up total number of particles
  !-----------------------------------------------------------

  if (command == 0) then
     call mpi_reduce(real(np_local, kind=4), np_total, 1, mpi_real, &
          mpi_sum, 0, mpi_comm_world, ierr)
  else if (command == 1) then
     call mpi_reduce(real(np_local_dm, kind=4), np_total, 1, mpi_real, &
          mpi_sum, 0, mpi_comm_world, ierr)
  else
     call mpi_reduce(real(np_local_h, kind=4), np_total, 1, mpi_real, &
          mpi_sum, 0, mpi_comm_world, ierr)
     np_h = np_local_h
     call mpi_allreduce(np_local_h,np_h,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  endif

  if (command == 0) then
     do j=1, np_local
        read(21) xvp(:,j)
     enddo
  else if (command == 1) then
     do j=1, np_local_dm
        read(21) xvp_dm(:,j)
     enddo
  else
     do j=1, np_local_h
        read(21) xvmp_h(1:3,j)
        read(21) garbage4 !! m (x2) + r (x2) 
        read(21) garbage3 !! xbar (x3)
        read(21) xvmp_h(4:6,j)
        read(21) garbage3 !! angular momentum
        read(21) garbage3 !! var in vel
        read(21) garbage3 !! var in pos
        read(21) garbage3 !! moment of inertia
        read(21) garbage3 !! moment of inertia
        read(21) garbage3 !! xbar nu
        read(21) garbage3 !! vbar nu
        read(21) garbage1 !! nbr nu
     enddo
  endif

  close(21)

  !-----------------------------------------------------------
  ! Get halo local coordinates
  !-----------------------------------------------------------

  if (command == 2) then
     do j=1, np_local_h
        xvmp_h(1:3,j) = xvmp_h(1:3, j) - slab_coord(:)*nc_node_dim
     enddo
  endif
 
#ifdef KAISER
  !-----------------------------------------------------------
  ! Include Kaiser effect
  !-----------------------------------------------------------
  ! Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
  ! Converting seconds into simulation time units
  ! cancels the H0...
  !-----------------------------------------------------------

  if (command == 0) then
     xvp(3,:) = xvp(3,:) + xvp(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
  else if (command == 1) then
     xvp_dm(3,:) = xvp_dm(3,:) + xvp_dm(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))
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
  if (rank == 0) then
     write(*,*) 'Finished read_particles ... '
     if (command == 0) write(*,*) 'particles read : neutrinos'
     if (command == 1) write(*,*) 'particles read : dark matter'
     if (command == 2) write(*,*) 'particles read : halos'
     write(*,*) 'Total number of particles = ', int(np_total,8)
     write(*,*) 'Time elapsed = ', time2-time1
     write(*,*) '--------------------------'
  endif
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

  write(*,*) '****************************************************'
  write(*,*) 'cic_velpower_halo running ...'
  write(*,*) '****************************************************'
  write(*,*) 'nodes   :', nodes
  write(*,*) 'nc      :', nc
  write(*,*) 'np      :', np
  write(*,*) 'box     :',box
  write(*,*) '****************************************************'

  return
end subroutine writeparams

! -------------------------------------------------------------------------------------------------------

subroutine writepowerspectra(command, convert)
  !
  ! Writes the dimensionless power spectrum for the curl/divergence of the momentum density field
  !    

  implicit none
  
  integer      :: i, j, k
  character*180 :: fn
  character*6  :: prefix
  character*7  :: z_write
  real    :: vsim2phys, zcur
  integer(4) :: command, convert

  !-----------------------------------------------------------
  ! Determine conversion factor for sim velocity to physical
  !-----------------------------------------------------------

  zcur      = z_checkpoint(cur_checkpoint)
  vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

  !-----------------------------------------------------------
  ! Determine name of output file
  !-----------------------------------------------------------
  
  write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
  z_write=adjustl(z_write)
    
#ifdef NGP 
  prefix = 'ngpvps'
#else
  prefix = 'cicvps'
#endif

  if (command == 0) then
     fn = '_nu.dat'
  else if (command == 1) then
     fn = '_dm.dat'
  else if (command == 2) then
     fn = '_dmnu.dat'
  else if (command == 3) then
     fn = '_rel.dat'
  else if (command == 4) then
     fn = '_halo.dat'
  else if (command == 5) then
     fn = '_nulin_grad.dat'
  else if (command == 6) then
     fn = '_dmlin_grad.dat'
  else if (command == 7) then
     fn = '_halolin_grad.dat'
  else if (command == 8) then
     fn = '_nuhalo_grad.dat'
  else if (command == 9) then
     fn = '_nudm_grad.dat'
  else if (command == 10) then
     fn = '_dmhalo.dat'
  else if (command == 11) then
     fn = '_halohalolin.dat'
  else if (command == 12) then
     fn = '_dmhalolin.dat'
  else if (command == 13) then
     fn = '_dmdmlin.dat'
  else if (command == 14) then
     fn = '_nunuh.dat'
  else if (command == 15) then
     fn = '_nunudm.dat'
  else if (command == 16) then
     fn = '_relh_grad.dat'
  else if (command == 17) then
     fn = '_reldm_grad.dat'
  else if (command == 18) then
     fn = '_relrelh.dat'
  else if (command == 19) then
     fn = '_relreldm.dat'
  else if (command == 20) then
     fn = '_rel2.dat'
  else if (command == 21) then
     fn = '_rel2relh.dat'
  else if (command == 22) then
     fn = '_halonu.dat'
  else if (command == 23) then
     fn = '_nulin_vtf.dat'
  else if (command == 24) then
     fn = '_dmlin_vtf.dat'
  else if (command == 25) then
     fn = '_dmlinnulin_vtf.dat'
  else if (command == 26) then
     fn = '_halolin_vtf.dat'
  else if (command == 27) then
     fn = '_dmlinhalolin_vtf.dat'
  else if (command == 28) then
     fn = '_nudm_vtf.dat'
  else if (command == 29) then
     fn = '_reldm_vtf.dat'
  else if (command == 30) then
     fn = '_nuhalo_vtf.dat'
  else if (command == 31) then
     fn = '_relh_vtf.dat'
  else if (command == 32) then
     fn = '_halolinnu.dat'
  else if (command == 33) then
     fn = '_dmnuhalo.dat'
  else if (command == 34) then
     fn = '_halolinnuh_grad.dat'
  else
     fn = '.dat'
  endif

#ifdef KAISER
  fn = '-RSD_'//fn
#endif
  fn = output_path//z_write(1:len_trim(z_write))//prefix//fn

  !-----------------------------------------------------------
  ! Asign data to be written
  !-----------------------------------------------------------

  do i = 1, nc
     pkdm(:, i) = 0.
  enddo

  do i = 1, nc
     do j = 1, 3
        pkdm(1, i) = pkdm(1, i) + pk_all(j, 1, i, command)
        pkdm(2, i) = pkdm(2, i) + pk_all(j, 2, i, command)
     enddo
     pkdm(3, i) = pk_all(1, 3, i, command)
  enddo

  !-----------------------------------------------------------
  ! Convert to physical units in km/s if needed
  !-----------------------------------------------------------
 
  if (convert == 1) then 
     do i = 1, nc
        pkdm(1, i) = vsim2phys**2 * pkdm(1, i)
        pkdm(2, i) = vsim2phys**2 * pkdm(2, i) 
     enddo
  endif

  !-----------------------------------------------------------
  ! Write to output file with column ordering [k, p(k), sigma(k)]
  !-----------------------------------------------------------

  write(*,*) 'Writing power spectrum to file ...'
  write(*,*) 'file :', fn
  write(*,*) '--------------------------'

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

subroutine writedenpowerspectra(command)
  !
  ! Writes the dimensionless power spectrum
  !    

  implicit none
  
  integer      :: i, j, k
  character*180 :: fn
  character*5  :: prefix
  character*7  :: z_write
  real    :: zcur
  integer(4) :: command

  !-----------------------------------------------------------
  ! Determine current redshift
  !-----------------------------------------------------------

  zcur = z_checkpoint(cur_checkpoint)
  
  !-----------------------------------------------------------
  ! Determine name of output file
  !-----------------------------------------------------------
  
  write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
  z_write=adjustl(z_write)
    
#ifdef NGP 
  prefix = 'ngpps'
#else
  prefix = 'cicps'
#endif

  if (command == 0) then
     fn = '_nu.dat'
  else if (command == 1) then
     fn = '_dm.dat'
  else if (command == 2) then
     fn = '_dmnu.dat'
  else if (command == 3) then
     fn = '_rel.dat'
  else if (command == 4) then
     fn = '_halo.dat'
  else if (command == 5) then
     fn = '_nu.dat'
  else if (command == 6) then
     fn = '_dm.dat'
  else if (command == 7) then
     fn = '_halo.dat'
  else if (command == 8) then
     fn = '_nuhalo_grad.dat'
  else if (command == 9) then
     fn = '_nudm_grad.dat'
  else if (command == 10) then
     fn = '_dmhalo.dat'
  else if (command == 11) then
     fn = '_halohalolin.dat'
  else if (command == 12) then
     fn = '_dmhalolin.dat'
  else if (command == 13) then
     fn = '_dmdmlin.dat'
  else if (command == 14) then
     fn = '_nunuh.dat'
  else if (command == 15) then
     fn = '_nunudm.dat'
  else if (command == 16) then
     fn = '_relh_grad.dat'
  else if (command == 17) then
     fn = '_reldm_grad.dat'
  else if (command == 18) then
     fn = '_relrelh.dat'
  else if (command == 19) then
     fn = '_relreldm.dat'
  else if (command == 20) then
     fn = '_rel2.dat'
  else if (command == 21) then
     fn = '_rel2relh.dat'
  else if (command == 22) then
     fn = '_halonu.dat'
  else if (command == 23) then
     fn = '_nulin_vtf.dat'
  else if (command == 24) then
     fn = '_dmlin_vtf.dat'
  else if (command == 25) then
     fn = '_dmlinnulin_vtf.dat'
  else if (command == 26) then
     fn = '_halolin_vtf.dat'
  else if (command == 27) then
     fn = '_dmlinhalolin_vtf.dat'
  else if (command == 28) then
     fn = '_nudm_vtf.dat'
  else if (command == 29) then
     fn = '_reldm_vtf.dat'
  else if (command == 30) then
     fn = '_nuhalo_vtf.dat'
  else if (command == 31) then
     fn = '_relh_vtf.dat'
  else if (command == 32) then
     fn = '_halolinnu.dat'
  else if (command == 33) then
     fn = '_dmnuhalo.dat'
  else if (command == 34) then
     fn = '_halolinnuh_grad.dat'
  else
     fn = '.dat'
  endif

#ifdef KAISER
  fn = '-RSD_'//fn
#endif
  fn = output_path//z_write(1:len_trim(z_write))//prefix//fn

  !-----------------------------------------------------------
  ! Asign data to be written
  !-----------------------------------------------------------

  do i = 1, nc
     pkdm(:, i) = 0.
  enddo

  do i = 1, nc
     do j = 1, 3
        pkdm(1, i) = pkdm(1, i) + pk_all(j, 1, i, command)
        pkdm(2, i) = pkdm(2, i) + pk_all(j, 2, i, command)
     enddo
     pkdm(3, i) = pk_all(1, 3, i, command)
  enddo

  !-----------------------------------------------------------
  ! Write to output file with column ordering [k, p(k), sigma(k)]
  !-----------------------------------------------------------

  write(*,*) 'Writing density power spectrum to file ...'
  write(*,*) 'file :', fn
  write(*,*) '--------------------------'

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

end subroutine writedenpowerspectra

! -------------------------------------------------------------------------------------------------------

subroutine darkmatter

  implicit none
  
  integer :: i, j, k
  real    :: d, dmin, dmax, sum_dm, sum_dm_local, dmint, dmaxt
  real*8  :: dsum, dvar, dsumt, dvart

  time1 = mpi_wtime(ierr)  
  if (rank == 0) write(*, *) 'Starting darkmatter ...'
  
  !--------------------------------------------------------
  ! Initialize FFT array to zero 
  !--------------------------------------------------------

  do k = 1, nc_node_dim
     cube(:, :, k) = 0.
  enddo

  !--------------------------------------------------------
  ! Assign data to density grid
  !--------------------------------------------------------
  
  cube(:, :, :) = momden(1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
  
  !--------------------------------------------------------
  ! Calculate some statistics
  !--------------------------------------------------------
  
  sum_dm_local = sum(cube) 
  call mpi_reduce(sum_dm_local, sum_dm, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
  if (rank == 0) write(*,*) "CUBE total sum = ", sum_dm
  
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
     write(*,*) 'Cube min ', dmint
     write(*,*) 'Cube max ', dmaxt
     write(*,*) 'Cube sum ', real(dsum)
     write(*,*) 'Cube var ', real(dvar)
  endif

  !--------------------------------------------------------
  ! Forward FFT 
  !--------------------------------------------------------
  
  call cp_fftw(1)
  
  time2 = mpi_wtime(ierr)

  if (rank == 0) then
     write(*,*) 'Finished darkmatter ...'
     write(*,*) "Elapsed time = ", time2 - time1
     write(*,*) '--------------------------'
  endif

  return
end subroutine darkmatter

! -------------------------------------------------------------------------------------------------------

subroutine pass_particles(command)
  !
  ! Pass particles inside buffer space to their appropriate nodes.
  !

  ! Note : it would be possible to remove the xp_buf array by doing as in subroutine buffer_particles_groups
    
  implicit none

  integer i,pp,np_buf,np_exit,np_final,npo,npi
  real x(6),lb,ub
  integer, dimension(mpi_status_size) :: status,sstatus,rstatus
  integer :: tag,srequest,rrequest,sierr,rierr
  real(4), parameter :: eps = 1.0e-03
  integer(4) :: command

  time1 = mpi_wtime(ierr)

  if (rank == 0) write(*,*) 'Starting pass_particles ...'

  !--------------------------------------------------------
  ! Identify particles within the buffer
  !--------------------------------------------------------

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

     !! Read its position and velocity  
     if (command == 0) then
        x = xvp(:, pp)
     else if (command == 1) then
        x = xvp_dm(:, pp)
     else
        x = xvmp_h(1:6,pp)
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

  if (rank == 0) write(*,*) 'Total exiting particles = ',np_exit

  !--------------------------------------------------------
  ! Pass +x
  !--------------------------------------------------------

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

  !--------------------------------------------------------
  ! Pass -x
  !--------------------------------------------------------

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

  !--------------------------------------------------------
  ! Pass +y
  !--------------------------------------------------------

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

  !--------------------------------------------------------
  ! Pass -y
  !--------------------------------------------------------

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
  
  !--------------------------------------------------------
  ! Pass +z
  !--------------------------------------------------------

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

  !--------------------------------------------------------
  ! Pass -z
  !--------------------------------------------------------

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
     call mpi_reduce(np_local,np_final,1,mpi_integer,mpi_sum,0, &
          mpi_comm_world,ierr)
  else if (command == 1) then
     call mpi_reduce(np_local_dm,np_final,1,mpi_integer,mpi_sum,0, &
          mpi_comm_world,ierr)
  else
     call mpi_reduce(np_local_h,np_final,1,mpi_integer,mpi_sum,0, &
          mpi_comm_world,ierr)
  endif
  
  if (rank == 0) then
     if (command == 0) then
        if (np_final /= np**3) then
           print *,'ERROR: total number of neutrinos incorrect after passing'
        endif
     else if (command == 1) then 
        if (np_final /= (np/ratio_nudm_dim)**3) then
           print *,'ERROR: total number of dark matter particles incorrect after passing'
        endif
     else
        if (np_final /= np_h) then
           print *,'ERROR: total number of halos incorrect after passing'
        endif
     endif
  endif
 
  !--------------------------------------------------------
  ! Check for particles out of bounds
  !--------------------------------------------------------

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
  if (rank == 0) then
     write(*,*) 'Finished pass_particles ...'
     write(*,*) 'time elapsed             = ', time2-time1
     write(*,*) 'total buffered particles =',np_exit
     write(*,*) 'total particles          =',np_final
     write(*,*) '--------------------------'
  endif
  return
  
end subroutine pass_particles

! -------------------------------------------------------------------------------------------------------

subroutine order_xvp_ll(command, glook, test_order)
  !
  ! Ordre the xvp array acoording to the particles's coarse grid cell appartenance
  ! Treat both groups separately
  !

  implicit none

  integer    :: k,j,i,k2,j2,i2,pp,pp2,base,command,pp_up,pp_down
  integer(8) :: index, index2
  integer(1) :: glook
  real       :: xvtemp(6), xvswap(6)
  logical    :: test_order

  real(8)    :: sec1, sec2, sec3

  sec1 = mpi_wtime(ierr)

  if (rank == 0) then
     write(*,*) 'Starting order_xvp_ll ...'
  endif

  !--------------------------------------------------------
  ! Determine extent of the xvp array we have to work on
  !--------------------------------------------------------

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

     !--------------------------------------------------------
     ! Constructs hoc
     !--------------------------------------------------------

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

     !--------------------------------------------------------
     ! Actually sort the xvp_array
     !--------------------------------------------------------

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

     !--------------------------------------------------------
     ! Test that the xvp array has the right order
     !--------------------------------------------------------

     !! Now the xvp array should be ordered, we want to test it

     if (test_order) then
        if (rank == 0) then
           write(*,*) 'Checking sorting is right ...'
        endif

        index2 = 0
        do pp = pp_down, pp_up
           xvtemp(1:3) = xvp(1:3,pp)
           
           i = floor(xvtemp(1)/mesh_scale) + 1
           j = floor(xvtemp(2)/mesh_scale) + 1
           k = floor(xvtemp(3)/mesh_scale) + 1
           
           index = (i-hoc_nc_l) + (j-hoc_nc_l) * base + (k-hoc_nc_l)*base*base
           
           if (index >= index2) then
              index2 = index
           else
              write(*,*) 'WARNING : xvp array not correctly sorted !!! CHECK ALGORITHM',pp, index, index2
           endif
        enddo
     endif
  
  else if (command == 1) then

     !--------------------------------------------------------
     ! Do the same for dark matter
     !--------------------------------------------------------

     do k = hoc_nc_l, hoc_nc_h
        hoc_dm(:,:,k) = 0
     enddo

     do pp = pp_down,pp_up
        xvtemp(1:3) = xvp_dm(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale) + 1
        j = floor(xvtemp(2)/mesh_scale) + 1
        k = floor(xvtemp(3)/mesh_scale) + 1
        hoc_dm(i,j,k) = hoc_dm(i,j,k) + 1
     enddo
     
     pp = pp_down - 1
     do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
           do i = hoc_nc_l, hoc_nc_h
              pp = hoc_dm(i,j,k) + pp
              hoc_dm(i,j,k) = pp
           enddo
        enddo
     enddo

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

     if (test_order) then
        if (rank == 0) then
           write(*,*) 'Checking sorting is right ...'
        endif

        index2 = 0
        do pp = pp_down, pp_up
           xvtemp(1:3) = xvp_dm(1:3,pp)
           
           i = floor(xvtemp(1)/mesh_scale) + 1
           j = floor(xvtemp(2)/mesh_scale) + 1
           k = floor(xvtemp(3)/mesh_scale) + 1
           
           index = (i-hoc_nc_l) + (j-hoc_nc_l) * base + (k-hoc_nc_l)*base*base
           
           if (index >= index2) then
              index2 = index
           else
              write(*,*) 'WARNING : xvp array not correctly sorted !!! CHECK ALGORITHM',pp, index, index2
           endif
        enddo
     endif

  else if (command == 2) then

     !--------------------------------------------------------
     ! Do the same for halos
     !--------------------------------------------------------

     do k = hoc_nc_l_h, hoc_nc_h_h
        hoc_h(:,:,k) = 0
     enddo

     do pp = pp_down,pp_up
        xvtemp(1:3) = xvmp_h(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale_h) + 1
        j = floor(xvtemp(2)/mesh_scale_h) + 1
        k = floor(xvtemp(3)/mesh_scale_h) + 1
        hoc_h(i,j,k) = hoc_h(i,j,k) + 1
     enddo
     
     pp = pp_down - 1
     do k = hoc_nc_l_h, hoc_nc_h_h
        do j = hoc_nc_l_h, hoc_nc_h_h
           do i = hoc_nc_l_h, hoc_nc_h_h
              pp = hoc_h(i,j,k) + pp
              hoc_h(i,j,k) = pp
           enddo
        enddo
     enddo

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

     if (test_order) then
        if (rank == 0) then
           write(*,*) 'Checking sorting is right ...'
        endif

        index2 = 0
        do pp = pp_down, pp_up
           xvtemp(1:3) = xvmp_h(1:3,pp)
           
           i = floor(xvtemp(1)/mesh_scale_h) + 1
           j = floor(xvtemp(2)/mesh_scale_h) + 1
           k = floor(xvtemp(3)/mesh_scale_h) + 1
           
           index = (i-hoc_nc_l_h) + (j-hoc_nc_l_h) * base + (k-hoc_nc_l_h)*base*base
           
           if (index >= index2) then
              index2 = index
           else
              write(*,*) 'WARNING : xvp array not correctly sorted !!! CHECK ALGORITHM',pp, index, index2
           endif
        enddo
     endif

  endif

  sec2 = mpi_wtime(ierr)
  if (rank == 0) then
     write(*,*) 'Finished order_xvp_ll'
     write(*,*) 'Group processed     : ', glook
     write(*,*) 'Number of particles : ', (pp_up-pp_down+1)
     write(*,*) "Elapsed time        = ", sec2 - sec1
     write(*,*) '--------------------------'
  endif

end subroutine order_xvp_ll

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

  time1 = mpi_wtime(ierr)

  !--------------------------------------------------------    
  ! Assign groups + put group1 particles at the end of the xvp array
  !--------------------------------------------------------    

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
           xvp(:,max_np-np1) = xvp(:,pp)
           xvp(:,pp) = xvp(:,np_local - np1)
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
           xvp_dm(:,max_np_dm-np1) = xvp_dm(:,pp)
           xvp_dm(:,pp) = xvp_dm(:,np_local_dm - np1)
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
           xvmp_h(:,max_np_h-np1) = xvmp_h(:,pp)
           xvmp_h(:,pp) = xvmp_h(:,np_local_h - np1)
           np1 = np1 + 1
        endif
     enddo
  endif

  !--------------------------------------------------------    
  ! Assign np_groups in consequence
  !--------------------------------------------------------    

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
     write(*,*) 'Finished order_xvp_groups ...'
     write(*,*) 'np_groups_tot(0) = ', np1_tot
     write(*,*) 'np_groups_tot(1) = ', np2_tot
     write(*,*) 'Total = ', np1_tot+np2_tot
     write(*,*) 'Elapsed time = ', time2 - time1
     write(*,*) '--------------------------'
  endif

end subroutine order_xvp_groups

! -------------------------------------------------------------------------------------------------------

subroutine remake_hoc(command, glook)
  !
  ! Reconstructs hoc just before velocity density
  ! That way, hoc can be equivalenced
  !

  implicit none
  integer      :: command,pp,pp_down,pp_up,k,j,i
  real         :: xvtemp(6)
  integer(1)   :: glook

  !--------------------------------------------------------
  ! Determine xvp spot we have to look at
  !--------------------------------------------------------

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

  !--------------------------------------------------------
  ! Actually constructs hoc
  !--------------------------------------------------------

  if (command == 0) then
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
     
  else if (command == 1) then
     do k = hoc_nc_l, hoc_nc_h
        hoc_dm(:,:,k) = 0
     enddo

     do pp = pp_down,pp_up
        xvtemp(1:3) = xvp_dm(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale) + 1
        j = floor(xvtemp(2)/mesh_scale) + 1
        k = floor(xvtemp(3)/mesh_scale) + 1
        hoc_dm(i,j,k) = hoc_dm(i,j,k) + 1
     enddo

     pp = pp_down - 1
     do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
           do i = hoc_nc_l, hoc_nc_h
              pp = hoc_dm(i,j,k) + pp
              hoc_dm(i,j,k) = pp
           enddo
        enddo
     enddo

  else if (command == 2) then
     do k = hoc_nc_l_h, hoc_nc_h_h
        hoc_h(:,:,k) = 0
     enddo

     do pp = pp_down,pp_up
        xvtemp(1:3) = xvmp_h(1:3,pp)
        i = floor(xvtemp(1)/mesh_scale_h) + 1
        j = floor(xvtemp(2)/mesh_scale_h) + 1
        k = floor(xvtemp(3)/mesh_scale_h) + 1
        hoc_h(i,j,k) = hoc_h(i,j,k) + 1
     enddo

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

end subroutine remake_hoc

! -------------------------------------------------------------------------------------------------------

subroutine buffer_particles_groups(command,glook)
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

  time1 = mpi_wtime(ierr)
  
  if (rank == 0) then 
     write(*,*) 'Starting buffer_particles_groups ...'
  endif
  
  !--------------------------------------------------------
  ! Pass +x
  !--------------------------------------------------------

  tag = 11
  np_buf = 0
  
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

  !--------------------------------------------------------
  ! Pass -x
  !--------------------------------------------------------
  
  tag = 12
  np_buf = 0
  
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

  !--------------------------------------------------------
  ! Pass +y
  !--------------------------------------------------------

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

  !--------------------------------------------------------  
  ! Pass -y
  !--------------------------------------------------------  
  
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

  !--------------------------------------------------------    
  ! Pass +z
  !--------------------------------------------------------    

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

  !--------------------------------------------------------    
  ! Pass -z
  !--------------------------------------------------------    
  
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

  !--------------------------------------------------------    
  ! Calculate some statistics
  !--------------------------------------------------------    
  
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
     write(*,*) 'Finished buffer_particles_groups ...'
     write(*,*) 'Group processed        : ',glook
     write(*,*) 'Time elapsed           : ',time2 - time1
     write(*,*) 'Total Buffer Particles : ', np_buffer_sent
     write(*,*) '--------------------------'
  endif

end subroutine buffer_particles_groups

! -------------------------------------------------------------------------------------------------------

subroutine initialize_random_number
  !
  ! Initialize random number generator
  !

  implicit none

  integer :: i, j, k, fstat
  real(4) :: x
  integer :: seedsize, clock
  integer, dimension(8) :: values
  integer(4), allocatable, dimension(:) :: iseed
  integer(4), allocatable, dimension(:) :: iseed_all
  character(len=4)   :: rank_string
  character(len=180) :: fn

  if (rank == 0) write(*,*) 'Starting initialize_random_number ...'
    
  write(rank_string, '(i4)') rank
  rank_string = adjustl(rank_string)

  if (generate_seed) then

     call random_seed()
     call random_seed(size=seedsize)

     allocate(iseed(seedsize))
     allocate(randomseed(seedsize))
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

     randomseedsize = seedsize
     randomseed(1:seedsize) = iseed(1:seedsize)

     call random_seed(put=iseed(1:seedsize))
  else
     fn = output_path//'node'//rank_string(1:len_trim(rank_string))//'/seed'//rank_string(1:len_trim(rank_string))

     open(unit=21,file=fn,status="old",iostat=fstat,access="stream")

     if (fstat /= 0) then
        write(*,*) 'ERROR: Cannot open seed file'
        write(*,*) 'rank', rank, ' file: ',fn
        call mpi_abort(mpi_comm_world, ierr, ierr)
     endif

     if (rank == 0) write(*,*) 'Catching seeds from file'

     read(21) randomseedsize
     allocate(randomseed(randomseedsize))
     do k = 1, randomseedsize
        read(21) randomseed(k)
     enddo
     close(21)

     call random_seed(put=randomseed(1:randomseedsize))
  endif

  if (write_seed) then

     if (rank == 0) write(*,*) 'Writing seeds to file'
     fn = output_path//'node'//rank_string(1:len_trim(rank_string))//'/seed'//rank_string(1:len_trim(rank_string))
     open(unit=11, file=fn, status="replace", iostat=fstat, access="stream")
     write(11) randomseedsize
     do k = 1, randomseedsize
        write(11) randomseed(k)
     enddo
     close(11)
  endif

  if (rank == 0) then
     write(*,*) 'Finished initialize_random_number ... '
     write(*,*) 'seedsize : ', randomseedsize
     write(*,*) '--------------------------'
  endif

  return
end subroutine initialize_random_number

! -------------------------------------------------------------------------------------------------------

subroutine reinitialize_random_number
  !                                                                                                                                                 
  ! Re initialize the random number generator using the seed                                                                                        
  !                                                                                                                                            

  call random_seed(put=randomseed(1:randomseedsize))

  if (rank == 0) then
     write(*,*) 'Finished reinitialize_random_number ... '
     write(*,*) 'seedsize : ', randomseedsize
     write(*,*) '--------------------------'
  endif

end subroutine reinitialize_random_number

! -------------------------------------------------------------------------------------------------------

subroutine clear_groups
  !
  ! Assign all particles to group 0
  ! 
  
  implicit none
  integer :: pp_start, pp_up, pp_down, pp, command

  time1 = mpi_wtime(ierr)  

  !--------------------------------------------------------    
  ! Reput all particles at the beginnig of the xvp array
  !--------------------------------------------------------    

  do command = 0,2

     if (command == 0) then
        pp_start = np_groups(g0) + np_buf_groups(g0) + 1
        pp_up = max_np
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
     else if (command == 1) then
        pp_start = np_groups_dm(g0) + np_buf_groups_dm(g0) + 1
        pp_up = max_np_dm
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
     else if (command == 2) then
        pp_start = np_groups_h(g0) + np_buf_groups_h(g0) + 1
        pp_up = max_np_h
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
     endif

     if (command == 0) then 
        do pp = 0, (pp_up-pp_down)
           xvp(:,pp_start+pp) = xvp(:,pp_down+pp) 
        enddo
     else if (command == 1) then
        do pp = 0, (pp_up-pp_down)
           xvp_dm(:,pp_start+pp) = xvp_dm(:,pp_down+pp) 
        enddo
     else if (command == 2) then
        do pp = 0, (pp_up-pp_down)
           xvmp_h(:,pp_start+pp) = xvmp_h(:,pp_down+pp) 
        enddo
     endif

  enddo
    
  !--------------------------------------------------------    
  ! Update np_groups ...
  !--------------------------------------------------------    

  np_groups(g0) = np_groups(g0) + np_groups(g1)
  np_groups(g1) = 0
  
  np_groups_dm(g0) = np_groups_dm(g0) + np_groups_dm(g1)
  np_groups_dm(g1) = 0
  
  np_groups_h(g0) = np_groups_h(g0) + np_groups_h(g1)
  np_groups_h(g1) = 0
  
  np_groups_tot(g0) = np_groups_tot(g0) + np_groups_tot(g1)
  np_groups_tot(g1) = 0
  
  np_groups_tot_dm(g0) = np_groups_tot_dm(g0) + np_groups_tot_dm(g1)
  np_groups_tot_dm(g1) = 0
  
  np_groups_tot_h(0) = np_groups_tot_h(g0) + np_groups_tot_h(g1)
  np_groups_tot_h(g1) = 0

  np_buf_groups(g0) = np_buf_groups(g0) + np_buf_groups(g1)
  np_buf_groups(g1) = 0

  np_buf_groups_dm(g0) = np_buf_groups_dm(g0) + np_buf_groups_dm(g1)
  np_buf_groups_dm(g1) = 0

  np_buf_groups_h(g0) = np_buf_groups_h(g0) + np_buf_groups_h(g1)
  np_buf_groups_h(g1) = 0

  time2 = mpi_wtime(ierr)  
  if (rank == 0) then
     write(*,*) 'Finished clear_groups ...'     
     write(*,*) 'Time elapsed       : ',time2 - time1
     write(*,*) 'Total neutrinos    : ',np_groups_tot(0)
     write(*,*) 'Total dm particles : ',np_groups_tot_dm(0)
     write(*,*) 'Total halos        : ',np_groups_tot_h(0)
     write(*,*) '--------------------------'
  endif

  return  
end subroutine clear_groups

! -------------------------------------------------------------------------------------------------------

subroutine velocity_density_V2(command, glook, nfind, dir)
  !
  ! Determine velocity field by taking the average velocity of the closest
  ! particles to each grid centre.
  ! 
  
  implicit none
  
  integer :: ind, i, j, k, pp, thread, dir
  integer :: ic, jc, kc, pp_up, pp_down, pp_down0
  integer :: npart
  real, dimension(3) :: dr, rc
  real :: dmax,xvtemp(6)
  integer(4) :: command, nfind
  integer(1) :: glook
  logical :: converged, pb
  integer(8) :: num_notconverged, num2, sameR, num1, corner, num3
  real*8 :: rsum, rsumt
  real    :: v, r2, vswap, r2swap
  integer :: m,mm
  
  time1 = mpi_wtime(ierr)
  
  !-----------------------------------------------------------
  ! Initialize to zeros + remake hoc
  !-----------------------------------------------------------
  
  call remake_hoc(command, glook)

  do i = 0, nc_node_dim + 1
     momden(:, :, i) = 0.
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


  !-----------------------------------------------------------
  ! For each cell find the closest nfind particles and average their velocities
  !-----------------------------------------------------------

  num_notconverged = 0
  sameR = 0
  corner = 0
  rsum = 0.0

  !$omp  parallel num_threads(nt) default(shared) private(i, j, k, thread, rc, npart, ic, jc, kc, ind, pp, dr, converged, pb, dmax, v, r2, vswap, r2swap, m, mm, xvtemp, pp_up, pp_down) reduction(+:num_notconverged) reduction(+:sameR) reduction(+:rsum) reduction(+:corner)

  thread = 1
  thread = omp_get_thread_num() + 1

  if (command == 0) then  
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

              !! if dmax > 0, it will represent the max distance to search
              dmax = -1.0 !! OR  dmax = cell_search_r(num_ngbhs) + 1.0 OR dmax = 2*nfine_buf + 1.0 --> do not do that with NP_LINLOGMAX

              do ind = 1, num_ngbhs

                 !! if we are too far, we can already exit the loop
                 if (dmax >= 0 .and. cell_search_r(ind) >= dmax) exit !! OR JUST if(cell_search_r(ind) >= dmax) ...

                 !! Cell indices to search within
                 ic = int((i-1)/mesh_scale) + 1 + cell_search(1, ind)
                 jc = int((j-1)/mesh_scale) + 1 + cell_search(2, ind)
                 kc = int((k-1)/mesh_scale) + 1 + cell_search(3, ind)

                 if ( ic == hoc_nc_l .and. jc == hoc_nc_l .and. kc == hoc_nc_l ) then
                    pp_down = pp_down0
                    pp_up = hoc(ic,jc,kc)
                 else if (ic == hoc_nc_l .and. jc == hoc_nc_l) then
                    pp_down = hoc(hoc_nc_h, hoc_nc_h, kc-1) + 1
                    pp_up = hoc(hoc_nc_l, hoc_nc_l, kc)
                 else if (ic == hoc_nc_l) then
                    pp_down = hoc(hoc_nc_h, jc-1, kc) + 1
                    pp_up = hoc(hoc_nc_l, jc, kc)
                 else
                    pp_down = hoc(ic-1, jc, kc) + 1
                    pp_up = hoc(ic, jc, kc)
                 endif

#ifdef NP_LINLOGMAX

#ifdef NP_FAST
                 if ( (pp_up - pp_down + 1) > max_npart_cell_search) then
                    write(*,*) "ERROR: npart > max_npart_cell_search. Consider increasing max_npart_cell_search !!", max_npart_cell_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 xvp_cell(:,1:(pp_up - pp_down + 1),thread) = xvp(:,pp_down:pp_up)

                 do pp = 1, (pp_up - pp_down + 1)
                    dr(1:3) = xvp_cell(1:3, pp, thread) - rc(1:3)
                    xvp_cell(2+dir, pp, thread) = dr(1)**2 + dr(2)**2 + dr(3)**2
                 enddo

                 if (npart + pp_up - pp_down + 1 > max_npart_search) then
                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 do pp = 1, (pp_up - pp_down + 1)
                    npart = npart + 1
                    rpos(1:2, npart, thread) = xvp_cell(2+dir:3+dir, pp, thread)
                 enddo
#else
                 if (npart + pp_up - pp_down + 1 > max_npart_search) then
                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 do pp = pp_down, pp_up
                    xvtemp(1:6) = xvp(1:6,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    xvtemp(2+dir) = dr(1)**2 + dr(2)**2 + dr(3)**2
                    npart = npart + 1
                    rpos(1:2, npart, thread) = xvtemp(2+dir:3+dir)
                 enddo
#endif
#else

#ifdef NP_FAST
                 if ( (pp_up - pp_down + 1) > max_npart_cell_search) then
                    write(*,*) "ERROR: npart > max_npart_cell_search. Consider increasing max_npart_cell_search !!", max_npart_cell_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 !! Copy the particles to a temporary array
                 xvp_cell(:,1:(pp_up - pp_down + 1),thread) = xvp(:,pp_down:pp_up)

                 !! Compute distance to the center of the cell
                 !! and insert the nearest particles into the table
                 do pp = 1,(pp_up - pp_down + 1)
                    dr(:) = xvp_cell(1:3, pp, thread) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvp_cell(3+dir, pp, thread)
#else
                 do pp = pp_down, pp_up
                    xvtemp(1:6) = xvp(1:6,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvtemp(3+dir)
#endif
                    npart = npart + 1

                    !! if it's the first particle we found, just add it to the list                                                                    
                    if (npart == 1) then
                       rpos(1, 1, thread) = v
                       rpos(2, 1, thread) = r2
                       
                       !! else if we haven't found at least nfind particles, just append it to its right place                              
                    else if ((npart - 1) < nfind) then
                       m = 1
                       do while (m < npart)
                          if (r2 > rpos(2, m, thread)) then
                             exit
                          else if (r2 == rpos(2, m, thread)) then
                             exit
                          endif
                          m = m+1
                       enddo
                       
                       do while (m < npart)
                          vswap = rpos(1, m, thread)
                          r2swap = rpos(2, m, thread)
                          rpos(1, m, thread) = v
                          rpos(2, m, thread) = r2
                          v = vswap
                          r2 = r2swap
                          m = m+1
                       enddo
                       rpos(1, npart, thread) = v
                       rpos(2, npart, thread) = r2
                       
                       !! else push it to its right place and get rid of the furtest particle                                     
                    else if (r2 < rpos(2, 1, thread)) then
                       m = 2
                       do while (m <= nfind)
                          if (r2 > rpos(2, m, thread)) then
                             exit
                          else if (r2 == rpos(2, m, thread)) then
                             exit
                          endif
                          m = m+1
                       enddo

                       if (m == 2) then
                          pb = .false.
                       else if (rpos(2, 1, thread) == rpos(2, 2, thread)) then
                          pb = .true.
                       else if (pb) then
                          pb = .false.
                       endif

                       do mm = 1, m-1
                          vswap = rpos(1, m-mm, thread)
                          r2swap = rpos(2, m-mm, thread)
                          rpos(1, m-mm, thread) = v
                          rpos(2, m-mm, thread) = r2
                          v = vswap
                          r2 = r2swap
                       enddo
                    else if (r2 == rpos(2, 1, thread)) then
                       pb = .true.
                    endif
                 enddo !! pp
#endif

#ifdef NP_LINLOGMAX
                 if (npart >= nfind .and. dmax < 0) then
                    dmax = cell_search_r_max(ind)
                    converged = .true.
                 endif
#else
                 if (npart >= nfind) then
                    dmax = sqrt(rpos(2, 1, thread))
                    converged = .true.
                 endif
#endif                 
              enddo !! ind

#ifdef NP_LINLOGMAX
              ipos(:npart, thread) =  (/ (ic, ic=1, npart) /)
              call indexedsort(npart, rpos(2,:npart, thread), ipos(:npart, thread))
              v = 0.0
              r2 = 0.0
              if (npart > 0) then
                 do ind = 1, min(nfind,npart)
                    v = v + rpos(1, ipos(ind, thread), thread)
                 enddo
                 r2 = sqrt(rpos(2, ipos(min(nfind, npart), thread), thread))
#ifdef RADIUS
                 momden(i, j, k) = r2
#else
                 momden(i, j, k) = v / real(min(nfind,npart))
#endif            
                 if (npart > nfind) then
                    if ( rpos(2, ipos(nfind+1, thread), thread) == rpos(2, ipos(nfind, thread), thread) ) then
                       pb = .true.
                    endif
                 endif
              else
                 r2 = rnf_buf*sqrt(3.0)
#ifdef RADIUS
                 momden(i, j, k) = r2
#endif
              endif
#else
              v = 0.0
              r2 = 0.0
              if (npart > 0) then
                 do ind = 1, min(nfind, npart)
                    v = v + rpos(1, ind, thread)
                 enddo
                 r2 = sqrt(rpos(2, 1, thread))
#ifdef RADIUS
                 momden(i, j, k) = r2
#else
                 momden(i, j, k) = v / real(min(nfind, npart))
#endif            
              else
                 r2 = rnf_buf*sqrt(3.0)
#ifdef RADIUS
                 momden(i, j, k) = r2
#endif
              endif
#endif
              rsum = rsum + r2
              if (.not. converged ) then
                 num_notconverged = num_notconverged + 1
              endif

              if (r2 >= rnf_buf) then
                 corner = corner + 1
              endif

              if (pb) then
                 sameR = sameR + 1
              endif

           enddo !! i
        enddo !! j
     enddo !! k
     !$omp end do  

  else if (command == 1) then  
     !$omp do
     do k = 1, nc_node_dim
        do j = 1, nc_node_dim
           do i = 1, nc_node_dim

              rc(1) = i - 0.5
              rc(2) = j - 0.5
              rc(3) = k - 0.5

              npart = 0
              converged = .false.
              pb = .false.
              dmax = -1.0

              do ind = 1, num_ngbhs
                 if (dmax >= 0 .and. cell_search_r(ind) >= dmax) exit
                 ic = int((i-1)/mesh_scale) + 1 + cell_search(1, ind)
                 jc = int((j-1)/mesh_scale) + 1 + cell_search(2, ind)
                 kc = int((k-1)/mesh_scale) + 1 + cell_search(3, ind)

                 if ( ic == hoc_nc_l .and. jc == hoc_nc_l .and. kc == hoc_nc_l ) then
                    pp_down = pp_down0
                    pp_up = hoc_dm(ic,jc,kc)
                 else if (ic == hoc_nc_l .and. jc == hoc_nc_l) then
                    pp_down = hoc_dm(hoc_nc_h, hoc_nc_h, kc-1) + 1
                    pp_up = hoc_dm(hoc_nc_l, hoc_nc_l, kc)
                 else if (ic == hoc_nc_l) then
                    pp_down = hoc_dm(hoc_nc_h, jc-1, kc) + 1
                    pp_up = hoc_dm(hoc_nc_l, jc, kc)
                 else
                    pp_down = hoc_dm(ic-1, jc, kc) + 1
                    pp_up = hoc_dm(ic, jc, kc)
                 endif

#ifdef NP_LINLOGMAX
#ifdef NP_FAST
                 if ( (pp_up - pp_down + 1) > max_npart_cell_search) then
                    write(*,*) "ERROR: npart > max_npart_cell_search. Consider increasing max_npart_cell_search !!", max_npart_cell_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 xvp_cell(:,1:(pp_up - pp_down + 1),thread) = xvp_dm(:,pp_down:pp_up)

                 do pp = 1, (pp_up - pp_down + 1)
                    dr(1:3) = xvp_cell(1:3, pp, thread) - rc(1:3)
                    xvp_cell(2+dir, pp, thread) = dr(1)**2 + dr(2)**2 + dr(3)**2
                 enddo

                 if (npart + pp_up - pp_down + 1 > max_npart_search) then
                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 do pp = 1, (pp_up - pp_down + 1)
                    npart = npart + 1
                    rpos(1:2, npart, thread) = xvp_cell(2+dir:3+dir, pp, thread)
                 enddo
#else
                 if (npart + pp_up - pp_down + 1 > max_npart_search) then
                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 do pp = pp_down, pp_up
                    xvtemp(1:6) = xvp_dm(1:6,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    xvtemp(2+dir) = dr(1)**2 + dr(2)**2 + dr(3)**2
                    npart = npart + 1
                    rpos(1:2, npart, thread) = xvtemp(2+dir:3+dir)
                 enddo
#endif
#else

#ifdef NP_FAST
                 if ( (pp_up - pp_down + 1) > max_npart_cell_search) then
                    write(*,*) "ERROR: npart > max_npart_search. Consider increasing max_npart_search !!", max_npart_cell_search
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 xvp_cell(:,1:(pp_up - pp_down + 1),thread) = xvp_dm(:,pp_down:pp_up)

                 do pp = 1,(pp_up - pp_down + 1)
                    dr(:) = xvp_cell(1:3,pp,thread) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvp_cell(3+dir,pp,thread)
#else
                 do pp = pp_down, pp_up
                    xvtemp(:) = xvp_dm(:,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvtemp(3+dir)
#endif                    
                    npart = npart + 1

                    if (npart == 1) then
                       rpos(1, 1, thread) = v
                       rpos(2, 1, thread) = r2
                       
                    else if ((npart - 1) < nfind) then
                       m = 1
                       do while (m < npart)
                          if (r2 >= rpos(2, m, thread)) exit
                          m = m+1
                       enddo
                       
                       do while (m < npart)
                          vswap = rpos(1, m, thread)
                          r2swap = rpos(2, m, thread)
                          rpos(1, m, thread) = v
                          rpos(2, m, thread) = r2
                          v = vswap
                          r2 = r2swap
                          m = m+1
                       enddo
                       rpos(1, npart, thread) = v
                       rpos(2, npart, thread) = r2
                       
                    else if (r2 < rpos(2, 1, thread)) then
                       m = 2
                       do while (m <= nfind)
                          if (r2 >= rpos(2, m, thread)) exit
                          m = m+1
                       enddo

                       if (m == 2) then
                          pb = .false.
                       else if (rpos(2, 1, thread) == rpos(2, 2, thread)) then
                          pb = .true.
                       else if (pb) then
                          pb = .false.
                       endif

                       do mm = 1, m-1
                          vswap = rpos(1, m-mm, thread)
                          r2swap = rpos(2, m-mm, thread)
                          rpos(1, m-mm, thread) = v
                          rpos(2, m-mm, thread) = r2
                          v = vswap
                          r2 = r2swap
                       enddo

                    else if (r2 == rpos(2, 1, thread)) then
                       pb = .true.
                    endif
                 enddo !! pp
#endif
#ifdef NP_LINLOGMAX
                 if (npart >= nfind .and. dmax < 0) then
                    dmax = cell_search_r_max(ind)
                    converged = .true.
                 endif
#else
                 if (npart >= nfind) then
                    dmax = sqrt(rpos(2, 1, thread))
                    converged = .true.
                 endif
#endif                 
              enddo !! ind

#ifdef NP_LINLOGMAX
              ipos(:npart, thread) =  (/ (ic, ic=1, npart) /)
              call indexedsort(npart, rpos(2,:npart, thread), ipos(:npart, thread))
              v = 0.0
              r2 = 0.0
              if (npart > 0) then
                 do ind = 1, min(nfind,npart)
                    v = v + rpos(1, ipos(ind, thread), thread)
                 enddo
                 r2 = sqrt(rpos(2, ipos(min(nfind, npart), thread), thread))
#ifdef RADIUS
                 momden(i, j, k) = r2
#else
                 momden(i, j, k) = v / real(min(nfind, npart))
#endif
                 if (npart > nfind) then
                    if ( rpos(2, ipos(nfind+1, thread), thread) == rpos(2, ipos(nfind, thread), thread) ) then
                       pb = .true.
                    endif
                 endif
              else
                 r2 = rnf_buf*sqrt(3.0)
#ifdef RADIUS
                 momden(i, j, k) = r2
#endif
              endif
#else
              v = 0.0
              r2 = 0.0
              if (npart > 0) then
                 do ind = 1, min(nfind, npart)
                    v = v + rpos(1, ind, thread)
                 enddo
                 r2 = sqrt(rpos(2, 1, thread))
#ifdef RADIUS
                 momden(i, j, k) = r2
#else
                 momden(i, j, k) = v / real(min(nfind, npart))
#endif
              else
                 r2 = rnf_buf*sqrt(3.0)
#ifdef RADIUS
                 momden(i, j, k) = r2
#endif
              endif
#endif
              rsum = rsum + r2
              if (.not. converged ) then
                 num_notconverged = num_notconverged + 1
              endif

              if (r2 >= rnf_buf) then
                 corner = corner + 1
              endif

              if (pb) then
                 sameR = sameR + 1
              endif

           enddo !! i
        enddo !! j
     enddo !! k
     !$omp end do  

  else if (command == 2) then  
     !$omp do
     do k = 1, nc_node_dim
        do j = 1, nc_node_dim
           do i = 1, nc_node_dim

              rc(1) = i - 0.5
              rc(2) = j - 0.5
              rc(3) = k - 0.5

              npart = 0
              converged = .false.
              pb = .false.
              dmax = -1.0

              do ind = 1, num_ngbhs_h
                 if (dmax >= 0 .and. cell_search_r_h(ind) >= dmax) exit
                 ic = int((i-1)/mesh_scale_h) + 1 + cell_search_h(1, ind)
                 jc = int((j-1)/mesh_scale_h) + 1 + cell_search_h(2, ind)
                 kc = int((k-1)/mesh_scale_h) + 1 + cell_search_h(3, ind)

                 if ( ic == hoc_nc_l_h .and. jc == hoc_nc_l_h .and. kc == hoc_nc_l_h ) then
                    pp_down = pp_down0
                    pp_up = hoc_h(ic,jc,kc)
                 else if (ic == hoc_nc_l_h .and. jc == hoc_nc_l_h) then
                    pp_down = hoc_h(hoc_nc_h_h, hoc_nc_h_h, kc-1) + 1
                    pp_up = hoc_h(hoc_nc_l_h, hoc_nc_l_h, kc)
                 else if (ic == hoc_nc_l_h) then
                    pp_down = hoc_h(hoc_nc_h_h, jc-1, kc) + 1
                    pp_up = hoc_h(hoc_nc_l_h, jc, kc)
                 else
                    pp_down = hoc_h(ic-1, jc, kc) + 1
                    pp_up = hoc_h(ic, jc, kc)
                 endif

#ifdef NP_LINLOGMAX
#ifdef NP_FAST
                 if ( (pp_up - pp_down + 1) > max_npart_cell_search_h) then
                    write(*,*) "ERROR: npart > max_npart_cell_search_h. Consider increasing max_npart_cell_search_h !!", max_npart_cell_search_h
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 xvp_cell_h(:,1:(pp_up - pp_down + 1),thread) = xvmp_h(:,pp_down:pp_up)

                 do pp = 1, (pp_up - pp_down + 1)
                    dr(1:3) = xvp_cell_h(1:3, pp, thread) - rc(1:3)
                    xvp_cell_h(2+dir, pp, thread) = dr(1)**2 + dr(2)**2 + dr(3)**2
                 enddo

                 if (npart + pp_up - pp_down + 1 > max_npart_search_h) then
                    write(*,*) "ERROR: npart > max_npart_search_h. Consider increasing max_npart_search_h !!", max_npart_search_h
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 do pp = 1, (pp_up - pp_down + 1)
                    npart = npart + 1
                    rpos_h(1:2, npart, thread) = xvp_cell_h(2+dir:3+dir, pp, thread)
                 enddo
#else
                 if (npart + pp_up - pp_down + 1 > max_npart_search_h) then
                    write(*,*) "ERROR: npart > max_npart_search_h. Consider increasing max_npart_search_h !!", max_npart_search_h
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 do pp = pp_down, pp_up
                    xvtemp(1:6) = xvmp_h(1:6,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    xvtemp(2+dir) = dr(1)**2 + dr(2)**2 + dr(3)**2
                    npart = npart + 1
                    rpos_h(1:2, npart, thread) = xvtemp(2+dir:3+dir)
                 enddo
#endif
#else

#ifdef NP_FAST
                 if ( (pp_up - pp_down + 1) > max_npart_cell_search_h) then
                    write(*,*) "ERROR: npart > max_npart_search_h. Consider increasing max_npart_search !!", max_npart_cell_search_h
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                 endif

                 xvp_cell_h(:,1:(pp_up - pp_down + 1), thread) = xvmp_h(:,pp_down:pp_up)

                 do pp = 1,(pp_up - pp_down + 1)
                    dr(:) = xvp_cell_h(1:3, pp, thread) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvp_cell_h(3+dir, pp, thread)
#else
                 do pp = pp_down, pp_up
                    xvtemp(:) = xvmp_h(:,pp)
                    dr(:) = xvtemp(1:3) - rc(:)
                    r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                    v = xvtemp(3+dir)
#endif                    
                    npart = npart + 1

                    if (npart == 1) then
                       rpos_h(1, 1, thread) = v
                       rpos_h(2, 1, thread) = r2
                       
                    else if ((npart - 1) < nfind) then
                       m = 1
                       do while (m < npart)
                          if (r2 >= rpos(2, m, thread)) exit
                          m = m+1
                       enddo
                       
                       do while (m < npart)
                          vswap = rpos_h(1, m, thread)
                          r2swap = rpos_h(2, m, thread)
                          rpos_h(1, m, thread) = v
                          rpos_h(2, m, thread) = r2
                          v = vswap
                          r2 = r2swap
                          m = m+1
                       enddo
                       rpos_h(1, npart, thread) = v
                       rpos_h(2, npart, thread) = r2
                       
                    else if (r2 < rpos_h(2, 1, thread)) then
                       m = 2
                       do while (m <= nfind)
                          if (r2 >= rpos_h(2, m, thread)) exit
                          m = m+1
                       enddo

                       if (m == 2) then
                          pb = .false.
                       else if (rpos_h(2, 1, thread) == rpos_h(2, 2, thread)) then
                          pb = .true.
                       else if (pb) then
                          pb = .false.
                       endif

                       do mm = 1, m-1
                          vswap = rpos_h(1, m-mm, thread)
                          r2swap = rpos_h(2, m-mm, thread)
                          rpos_h(1, m-mm, thread) = v
                          rpos_h(2, m-mm, thread) = r2
                          v = vswap
                          r2 = r2swap
                       enddo

                    else if (r2 == rpos_h(2, 1, thread)) then
                       pb = .true.
                    endif
                 enddo !! pp
#endif
#ifdef NP_LINLOGMAX
                 if (npart >= nfind .and. dmax < 0) then
                    dmax = cell_search_r_max_h(ind)
                    converged = .true.
                 endif
#else
                 if (npart >= nfind) then
                    dmax = sqrt(rpos_h(2, 1, thread))
                    converged = .true.
                 endif
#endif                 
              enddo !! ind

#ifdef NP_LINLOGMAX
              ipos_h(:npart, thread) =  (/ (ic, ic=1, npart) /)
              call indexedsort(npart, rpos_h(2,:npart, thread), ipos_h(:npart, thread))
              v = 0.0
              r2 = 0.0
              if (npart > 0) then
                 do ind = 1, min(nfind,npart)
                    v = v + rpos_h(1, ipos_h(ind, thread), thread)
                 enddo
                 r2 = sqrt(rpos_h(2, ipos_h(min(nfind, npart), thread), thread))
#ifdef RADIUS
                 momden(i, j, k) = r2
#else
                 momden(i, j, k) = v / real(min(nfind,npart))
#endif
                 if (npart > nfind) then
                    if ( rpos_h(2, ipos_h(nfind+1, thread), thread) == rpos_h(2, ipos_h(nfind, thread), thread) ) then
                       pb = .true.
                    endif
                 endif
              else
                 r2 = rnf_buf_h*sqrt(3.0)
#ifdef RADIUS
                 momden(i, j, k) = r2
#endif
              endif
#else
              v = 0.0
              r2 = 0.0
              if (npart > 0) then
                 do ind = 1, min(nfind, npart)
                    v = v + rpos_h(1, ind, thread)
                 enddo
                 r2 = sqrt(rpos_h(2, 1, thread))
#ifdef RADIUS
                 momden(i,j,k) = r2
#else
                 momden(i,j,k) = v / real(min(nfind, npart))
#endif
              else
                 r2 = rnf_buf_h*sqrt(3.0)
#ifdef RADIUS
                 momden(i,j,k) = r2
#endif
              endif
#endif
              rsum = rsum + r2
              if (.not. converged ) then
                 num_notconverged = num_notconverged + 1
              endif

              if (r2 >= rnf_buf_h) then
                 corner = corner + 1
              endif

              if (pb) then
                 sameR = sameR + 1
              endif

           enddo !! i
        enddo !! j
     enddo !! k
     !$omp end do  

  endif !! command
  !$omp end parallel

  call mpi_reduce(num_notconverged, num2, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(sameR, num1, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(corner, num3, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

  call mpi_reduce(rsum, rsumt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  rsum = rsumt / nc**3

  time2 = mpi_wtime(ierr)

  if (rank == 0) then
     write(*,*) 'Finished velocity_density ...'

     if (command == 0) write(*,*) 'command = neutrinos'
     if (command == 1) write(*,*) 'command = dark matter'
     if (command == 2) write(*,*) 'command = halos'

     if (glook == g0) write(*,*) 'group = g0'
     if (glook == g1) write(*,*) 'group = g1'

     write(*,*) 'r_max mean = ', rsum

     if (num2 /= 0) then
        write(*,*) "WARNING: num_notconverged = ", num2, " ... consider increasing nfine_buf !!"
        write(*,*) "nfind = ", nfind 
     endif

     if (num1 /= 0) then
        write(*,*) 'WARNING: num sameR = ', num1
     endif

     if (num3 /= 0) then
        write(*,*) 'WARNING: num in the corner = ', num3
     endif

     write(*,*) 'Time elapsed : ', time2 - time1
     write(*,*) '--------------------------'
  endif

  return
end subroutine velocity_density_V2

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

  time1 = mpi_wtime(ierr)

  !-----------------------------------------------------------
  ! Pass +x
  !-----------------------------------------------------------
  
  tag = 111
  
  momden_send_buff(:, :) = momden(nc_node_dim+1, :, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(6), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(1, :, :) = momden(1, :, :) + momden_recv_buff(:, :)
  
  !-----------------------------------------------------------
  ! Pass -x
  !-----------------------------------------------------------
  
  tag = 112
  
  momden_send_buff(:, :) = momden(0, :, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(5), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(nc_node_dim, :, :) = momden(nc_node_dim, :, :) + momden_recv_buff(:, :)
  
  !-----------------------------------------------------------
  ! Pass +y
  !-----------------------------------------------------------
  
  tag = 113
  
  momden_send_buff(:, :) = momden(:, nc_node_dim+1, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(4), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, 1, :) = momden(:, 1, :) + momden_recv_buff(:, :)
  
  !-----------------------------------------------------------
  ! Pass -y
  !-----------------------------------------------------------
  
  tag = 114
  
  momden_send_buff(:, :) = momden(:, 0, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(3), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, nc_node_dim, :) = momden(:, nc_node_dim, :) + momden_recv_buff(:, :)
  
  !-----------------------------------------------------------
  ! Pass +z
  !-----------------------------------------------------------
  
  tag = 115
  
  momden_send_buff(:, :) = momden(:, :, nc_node_dim+1)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(2), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, :, 1) = momden(:, :, 1) + momden_recv_buff(:, :)
  
  !-----------------------------------------------------------
  ! Pass -z
  !-----------------------------------------------------------
  
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
  return
end subroutine buffer_momdensity

! -------------------------------------------------------------------------------------------------------

subroutine writebuffer_momden
  !
  ! Accumulate buffer from physical volume to adjacent nodes.
  !

  implicit none
  
  integer :: i
  integer, dimension(mpi_status_size) :: status, sstatus, rstatus
  integer :: tag, srequest, rrequest, sierr, rierr, ierr
  integer, parameter :: num2send = (nc_node_dim + 2)**2

  time1 = mpi_wtime(ierr)

  !------------------------------------------------------------
  ! Pass +x
  !------------------------------------------------------------
  
  tag = 111
  
  momden_send_buff(:, :) = momden(nc_node_dim, :, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(6), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(5), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(0, :, :) = momden_recv_buff(:, :)
  
  !------------------------------------------------------------
  ! Pass -x
  !------------------------------------------------------------
  
  tag = 112
  
  momden_send_buff(:, :) = momden(1, :, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(5), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(6), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(nc_node_dim+1, :, :) = momden_recv_buff(:, :)
  
  !------------------------------------------------------------
  ! Pass +y
  !------------------------------------------------------------
  
  tag = 113
  
  momden_send_buff(:, :) = momden(:, nc_node_dim, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(4), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(3), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, 0, :) = momden_recv_buff(:, :)
  
  !------------------------------------------------------------
  ! Pass -y
  !------------------------------------------------------------
  
  tag = 114
  
  momden_send_buff(:, :) = momden(:, 1, :)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(3), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(4), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, nc_node_dim+1, :) = momden_recv_buff(:, :)
  
  !------------------------------------------------------------
  ! Pass +z
  !------------------------------------------------------------
  
  tag = 115
  
  momden_send_buff(:, :) = momden(:, :, nc_node_dim)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(2), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(1), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, :, 0) = momden_recv_buff(:, :)
  
  !------------------------------------------------------------
  ! Pass -z
  !------------------------------------------------------------
  
  tag = 116
  
  momden_send_buff(:, :) = momden(:, :, 1)
  
  call mpi_isend(momden_send_buff, num2send, mpi_real, cart_neighbor(1), &
       tag, mpi_comm_world, srequest, sierr)
  call mpi_irecv(momden_recv_buff, num2send, mpi_real, cart_neighbor(2), &
       tag, mpi_comm_world, rrequest, rierr)
  call mpi_wait(srequest, sstatus, sierr)
  call mpi_wait(rrequest, rstatus, rierr)
  
  momden(:, :, nc_node_dim+1) = momden_recv_buff(:, :)
  
  time2 = mpi_wtime(ierr)
  return
end subroutine writebuffer_momden

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

    time1 = mpi_wtime(ierr)
    
    if (rank == 0) write(*,*) 'Starting powerspectrum ... '

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
    if (rank == 0) then
       write(*,*) 'Finished powerspectrum ... '
       write(*,*) 'Elapsed time = ', time2 - time1
       write(*,*) '--------------------------'
    endif
    return

  end subroutine powerspectrum

! -------------------------------------------------------------------------------------------------------

subroutine writevelocityfield(command, component, cubeOrMomden)
  !
  ! Writes the velocity field for the given component
  !

  implicit none

  integer :: j, k, fstat
  character(len=180) :: fn
  character(len=7)   :: z_write
  character(len=4)   :: rank_string
  real :: vsim2phys, zcur
  integer :: command, component, cubeOrMomden

  !------------------------------------------------------------
  ! Determine conversion to proper velocity [km/s]
  !------------------------------------------------------------

  if (rank == 0)  zcur = z_checkpoint(cur_checkpoint)
  call mpi_bcast(zcur, 1, mpi_real, 0, mpi_comm_world, ierr)

  vsim2phys = 300. * sqrt(omega_m) * box * (1. + zcur) / 2. / nc

  !------------------------------------------------------------
  ! Checkpoint and rank strings
  !------------------------------------------------------------

  write(z_write, '(f7.3)') zcur
  z_write = adjustl(z_write)
  
  write(rank_string, '(i4)') rank
  rank_string = adjustl(rank_string)

  !------------------------------------------------------------
  ! Write out velocity field for each dimension
  !------------------------------------------------------------

  if (command == 0) then
     fn = '_nu.bin'
  else if (command == 1) then
     fn = '_dm.bin'
  else if (command == 2) then
     fn = '_halo.bin'
  else if (command == 3) then
     fn = '_dmlin.bin'
  else if (command == 4) then
     fn = '_halolin.bin'
  else if (command == 5) then
     fn = '_nudm.bin'
  else if (command == 6) then
     fn = '_nuh.bin'
  else
     fn = '.bin'
  endif

  fn = rank_string(1:len_trim(rank_string))//fn

  if (component == 1) then
     fn = 'x'//fn
  else if (component == 2) then
     fn = 'y'//fn
  else if (component == 3) then
     fn = 'z'//fn
  endif

  fn = output_path//'node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//'vel'//fn

  open(unit=11, file=fn, status="replace", iostat=fstat, access="stream") 

  if (cubeOrMomden == 0) then
     do k = 1, nc_node_dim
        do j = 1, nc_node_dim
#ifndef RADIUS      
           write(11) momden(1:nc_node_dim, j, k) * vsim2phys
#else
           write(11) momden(1:nc_node_dim, j, k)
#endif
        enddo
     enddo

  else if (cubeOrMomden == 1) then
     do k = 1, nc_node_dim
        do j = 1, nc_node_dim      
           write(11) cube(1:nc_node_dim, j, k) * vsim2phys
        enddo
     enddo
  endif

  close(11)

  if (rank == 0) then
     write(*,*) 'Writing velocity field to file ...'
     write(*,*) 'file : ', fn
     write(*,*) '--------------------------'
  endif

  return
end subroutine writevelocityfield

! -------------------------------------------------------------------------------------------------------

subroutine writedensityfield(command, cubeOrMomden)
  !
  ! Writes the density field 
  !

  implicit none

  integer :: j, k, fstat
  character(len=180) :: fn
  character(len=7)   :: z_write
  character(len=4)   :: rank_string
  real :: zcur
  integer :: command, cubeOrMomden

  !------------------------------------------------------------
  ! Determine redshift
  !------------------------------------------------------------

  if (rank == 0)  zcur = z_checkpoint(cur_checkpoint)
  call mpi_bcast(zcur, 1, mpi_real, 0, mpi_comm_world, ierr)

  !------------------------------------------------------------
  ! Checkpoint and rank strings
  !------------------------------------------------------------

  write(z_write, '(f7.3)') zcur
  z_write = adjustl(z_write)
  
  write(rank_string, '(i4)') rank
  rank_string = adjustl(rank_string)

  !------------------------------------------------------------
  ! Write out density field for each dimension
  !------------------------------------------------------------

  if (command == 0) then
     fn = '_nu.bin'
  else if (command == 1) then
     fn = '_dm.bin'
  else if (command == 2) then
     fn = '_halo.bin'
  else if (command == 3) then
     fn = '_dmlin.bin'
  else if (command == 4) then
     fn = '_halolin.bin'
  else if (command == 5) then
     fn = '_nudm.bin'
  else if (command == 6) then
     fn = '_nuh.bin'
  else
     fn = '.bin'
  endif

  fn = rank_string(1:len_trim(rank_string))//fn

  fn = output_path//'node'//rank_string(1:len_trim(rank_string))//'/'//z_write(1:len_trim(z_write))//'den'//fn

  open(unit=11, file=fn, status="replace", iostat=fstat, access="stream") 

  if (cubeOrMomden == 0) then
     do k = 1, nc_node_dim
        do j = 1, nc_node_dim
           write(11) momden(1:nc_node_dim, j, k)
        enddo
     enddo

  else if (cubeOrMomden == 1) then
     do k = 1, nc_node_dim
        do j = 1, nc_node_dim      
           write(11) cube(1:nc_node_dim, j, k)
        enddo
     enddo
  endif

  close(11)

  if (rank == 0) then
     write(*,*) 'Writing density field to file ...'
     write(*,*) 'file : ', fn
     write(*,*) '--------------------------'
  endif

  return
end subroutine writedensityfield

! -------------------------------------------------------------------------------------------------------

#ifdef UCCFD
subroutine densityfield(command, glook)
#else
subroutine densityfield_V2(command, glook)
#endif
  !
  ! Compute the density field for the given species and given group
  !

  implicit none
  real       :: mp,xvtemp(3),dx1(3),dx2(3)
  integer    :: ic,jc,kc,pp,pp_up,pp_down,pp_down0,index1(3),index2(3),command, thread
  integer(1) :: glook

  time1 = mpi_wtime(ierr)

  if (command == 0) then
     mp = ncr**3 / real(np_groups_tot(glook))
  else if (command == 1) then
     mp = ncr**3 / real(np_groups_tot_dm(glook))
  else if (command == 2) then
     mp = ncr**3 / real(np_groups_tot_h(glook))
  endif
  
  !------------------------------------------------------------
  ! Initialize to zero + remake hoc
  !------------------------------------------------------------

  call remake_hoc(command, glook)

  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                 
  do kt = 0, nc_node_dim+1
     momden(:,:,kt) = 0.0
  enddo
  !$omp end parallel do

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

  !------------------------------------------------------------
  ! Assign mass to grid cell
  !------------------------------------------------------------

  thread = 1

  if (command == 0) then
     
     ! In principle, this could be omp-threaded
     ! since fine mesh cells are included in coarse cells.
    
     do kc = 0, nm_node_dim+1 
        !! we here assume nfine_buf > 0, with this, you avoid to have to call buffer_momdensity after density_field
        !! if, in subroutine mesh_buffer, we ONLY send the nearest particles, this part would have to become do kc = 1, nc_node_dim
        !! and we would have to call a buffer_momdensity at the end 

        do jc = 0, nm_node_dim+1
           do ic = 0, nm_node_dim+1
              
              if ( ic == hoc_nc_l .and. jc == hoc_nc_l .and. kc == hoc_nc_l ) then
                 pp_down = pp_down0
                 pp_up = hoc(ic,jc,kc)
              else if (ic == hoc_nc_l .and. jc == hoc_nc_l) then
                 pp_down = hoc(hoc_nc_h, hoc_nc_h, kc-1) + 1
                 pp_up = hoc(hoc_nc_l, hoc_nc_l, kc)
              else if (ic == hoc_nc_l) then
                 pp_down = hoc(hoc_nc_h, jc-1, kc) + 1
                 pp_up = hoc(hoc_nc_l, jc, kc)
              else
                 pp_down = hoc(ic-1, jc, kc) + 1
                 pp_up = hoc(ic, jc, kc)
              endif

#ifdef NP_FAST
              if ( (pp_up - pp_down + 1) > max_npart_cell_search) then
                 write(*,*) "ERROR: npart > max_npart_cell_search. Consider increasing max_npart_cell_search !!", max_npart_cell_search
                 call mpi_abort(mpi_comm_world, ierr, ierr)
              endif

              !! Copy the particles to a temporary array
              xvp_cell(1:3,1:(pp_up - pp_down + 1),thread) = xvp(1:3,pp_down:pp_up)

              do pp = 1,(pp_up - pp_down + 1)
                 index1(1:3) = 1 + floor(xvp_cell(1:3, pp, thread) - 0.5)
#else
              do pp = pp_down, pp_up
                 xvtemp(1:3) = xvp(1:3,pp) - 0.5
                 index1(1:3) = 1 + floor(xvtemp(1:3))
#endif
                 index2(1:3) = 1 + index1(1:3)
                 dx1(1:3) = real(index1(1:3)) - xvtemp(1:3)
                 dx2(1:3) = 1.0 - dx1(1:3)

                 if ( index1(1) >= 0 .and. index2(1) <= nc_node_dim + 1 .and. &
                      index1(2) >= 0 .and. index2(2) <= nc_node_dim + 1 .and. &
                      index1(3) >= 0 .and. index2(3) <= nc_node_dim + 1 ) then
                    
                    momden(index1(1),index1(2),index1(3)) = momden(index1(1),index1(2),index1(3)) +mp*dx1(1)*dx1(2)*dx1(3)
                    momden(index1(1),index1(2),index2(3)) = momden(index1(1),index1(2),index2(3)) +mp*dx1(1)*dx1(2)*dx2(3)
                    momden(index1(1),index2(2),index1(3)) = momden(index1(1),index2(2),index1(3)) +mp*dx1(1)*dx2(2)*dx1(3)
                    momden(index1(1),index2(2),index2(3)) = momden(index1(1),index2(2),index2(3)) +mp*dx1(1)*dx2(2)*dx2(3)

                    momden(index2(1),index1(2),index1(3)) = momden(index2(1),index1(2),index1(3)) +mp*dx2(1)*dx1(2)*dx1(3)
                    momden(index2(1),index1(2),index2(3)) = momden(index2(1),index1(2),index2(3)) +mp*dx2(1)*dx1(2)*dx2(3)
                    momden(index2(1),index2(2),index1(3)) = momden(index2(1),index2(2),index1(3)) +mp*dx2(1)*dx2(2)*dx1(3)
                    momden(index2(1),index2(2),index2(3)) = momden(index2(1),index2(2),index2(3)) +mp*dx2(1)*dx2(2)*dx2(3)

                 endif
              enddo !!pp
           enddo !! i
        enddo !! j
     enddo !! k

  else if (command == 1) then
     
     do kc = 0, nm_node_dim+1
        do jc = 0, nm_node_dim+1
           do ic = 0, nm_node_dim+1
              
              if ( ic == hoc_nc_l .and. jc == hoc_nc_l .and. kc == hoc_nc_l ) then
                 pp_down = pp_down0
                 pp_up = hoc_dm(ic,jc,kc)
              else if (ic == hoc_nc_l .and. jc == hoc_nc_l) then
                 pp_down = hoc_dm(hoc_nc_h, hoc_nc_h, kc-1) + 1
                 pp_up = hoc_dm(hoc_nc_l, hoc_nc_l, kc)
              else if (ic == hoc_nc_l) then
                 pp_down = hoc_dm(hoc_nc_h, jc-1, kc) + 1
                 pp_up = hoc_dm(hoc_nc_l, jc, kc)
              else
                 pp_down = hoc_dm(ic-1, jc, kc) + 1
                 pp_up = hoc_dm(ic, jc, kc)
              endif

#ifdef NP_FAST
              if ( (pp_up - pp_down + 1) > max_npart_cell_search) then
                 write(*,*) "ERROR: npart > max_npart_cell_search. Consider increasing max_npart_cell_search !!", max_npart_cell_search
                 call mpi_abort(mpi_comm_world, ierr, ierr)
              endif

              xvp_cell(1:3,1:(pp_up - pp_down + 1),thread) = xvp_dm(1:3,pp_down:pp_up)

              do pp = 1,(pp_up - pp_down + 1)
                 index1(1:3) = 1 + floor(xvp_cell(1:3, pp, thread) - 0.5)
#else
              do pp = pp_down, pp_up
                 xvtemp(1:3) = xvp_dm(1:3,pp) - 0.5
                 index1(1:3) = 1 + floor(xvtemp(1:3))
#endif
                 index2(1:3) = 1 + index1(1:3)
                 dx1(1:3) = real(index1(1:3)) - xvtemp(1:3)
                 dx2(1:3) = 1.0 - dx1(1:3)

                 if ( index1(1) >= 0 .and. index2(1) <= nc_node_dim + 1 .and. &
                      index1(2) >= 0 .and. index2(2) <= nc_node_dim + 1 .and. &
                      index1(3) >= 0 .and. index2(3) <= nc_node_dim + 1 ) then
                    
                    momden(index1(1),index1(2),index1(3)) = momden(index1(1),index1(2),index1(3)) +mp*dx1(1)*dx1(2)*dx1(3)
                    momden(index1(1),index1(2),index2(3)) = momden(index1(1),index1(2),index2(3)) +mp*dx1(1)*dx1(2)*dx2(3)
                    momden(index1(1),index2(2),index1(3)) = momden(index1(1),index2(2),index1(3)) +mp*dx1(1)*dx2(2)*dx1(3)
                    momden(index1(1),index2(2),index2(3)) = momden(index1(1),index2(2),index2(3)) +mp*dx1(1)*dx2(2)*dx2(3)

                    momden(index2(1),index1(2),index1(3)) = momden(index2(1),index1(2),index1(3)) +mp*dx2(1)*dx1(2)*dx1(3)
                    momden(index2(1),index1(2),index2(3)) = momden(index2(1),index1(2),index2(3)) +mp*dx2(1)*dx1(2)*dx2(3)
                    momden(index2(1),index2(2),index1(3)) = momden(index2(1),index2(2),index1(3)) +mp*dx2(1)*dx2(2)*dx1(3)
                    momden(index2(1),index2(2),index2(3)) = momden(index2(1),index2(2),index2(3)) +mp*dx2(1)*dx2(2)*dx2(3)

                 endif
              enddo !!pp
           enddo !! i
        enddo !! j
     enddo !! k

  else if (command == 2) then
     
     do kc = 0, nm_node_dim_h+1
        do jc = 0, nm_node_dim_h+1
           do ic = 0, nm_node_dim_h+1
              
              if ( ic == hoc_nc_l_h .and. jc == hoc_nc_l_h .and. kc == hoc_nc_l_h ) then
                 pp_down = pp_down0
                 pp_up = hoc_h(ic,jc,kc)
              else if (ic == hoc_nc_l_h .and. jc == hoc_nc_l_h) then
                 pp_down = hoc_h(hoc_nc_h_h, hoc_nc_h_h, kc-1) + 1
                 pp_up = hoc_h(hoc_nc_l_h, hoc_nc_l_h, kc)
              else if (ic == hoc_nc_l_h) then
                 pp_down = hoc_h(hoc_nc_h_h, jc-1, kc) + 1
                 pp_up = hoc_h(hoc_nc_l_h, jc, kc)
              else
                 pp_down = hoc_h(ic-1, jc, kc) + 1
                 pp_up = hoc_h(ic, jc, kc)
              endif

#ifdef NP_FAST
              if ( (pp_up - pp_down + 1) > max_npart_cell_search_h) then
                 write(*,*) "ERROR: npart > max_npart_search_h. Consider increasing max_npart_search !!", max_npart_cell_search_h
                 call mpi_abort(mpi_comm_world, ierr, ierr)
              endif

              xvp_cell_h(:,1:(pp_up - pp_down + 1), thread) = xvmp_h(:,pp_down:pp_up)

              do pp = 1,(pp_up - pp_down + 1)
                 index1(1:3) = 1 + floor(xvp_cell_h(1:3, pp, thread) - 0.5)
#else
              do pp = pp_down, pp_up
                 xvtemp(1:3) = xvmp_h(1:3,pp) - 0.5
                 index1(1:3) = 1 + floor(xvtemp(1:3))
#endif
                 index2(1:3) = 1 + index1(1:3)
                 dx1(1:3) = real(index1(1:3)) - xvtemp(1:3)
                 dx2(1:3) = 1.0 - dx1(1:3)

                 if ( index1(1) >= 0 .and. index2(1) <= nc_node_dim + 1 .and. &
                      index1(2) >= 0 .and. index2(2) <= nc_node_dim + 1 .and. &
                      index1(3) >= 0 .and. index2(3) <= nc_node_dim + 1 ) then
                    
                    momden(index1(1),index1(2),index1(3)) = momden(index1(1),index1(2),index1(3)) +mp*dx1(1)*dx1(2)*dx1(3)
                    momden(index1(1),index1(2),index2(3)) = momden(index1(1),index1(2),index2(3)) +mp*dx1(1)*dx1(2)*dx2(3)
                    momden(index1(1),index2(2),index1(3)) = momden(index1(1),index2(2),index1(3)) +mp*dx1(1)*dx2(2)*dx1(3)
                    momden(index1(1),index2(2),index2(3)) = momden(index1(1),index2(2),index2(3)) +mp*dx1(1)*dx2(2)*dx2(3)

                    momden(index2(1),index1(2),index1(3)) = momden(index2(1),index1(2),index1(3)) +mp*dx2(1)*dx1(2)*dx1(3)
                    momden(index2(1),index1(2),index2(3)) = momden(index2(1),index1(2),index2(3)) +mp*dx2(1)*dx1(2)*dx2(3)
                    momden(index2(1),index2(2),index1(3)) = momden(index2(1),index2(2),index1(3)) +mp*dx2(1)*dx2(2)*dx1(3)
                    momden(index2(1),index2(2),index2(3)) = momden(index2(1),index2(2),index2(3)) +mp*dx2(1)*dx2(2)*dx2(3)

                 endif
              enddo !!pp
           enddo !! i
        enddo !! j
     enddo !! k

  endif !! command


  time2 = mpi_wtime(ierr)

  if (rank == 0) then
     write(*,*) 'Finished densityfield ...'
     write(*,*) 'Time elapsed : ', time2 - time1
     write(*,*) '--------------------------'
  endif
  return

#ifdef UCCFD
end subroutine densityfield
#else
end subroutine densityfield_V2
#endif

! -------------------------------------------------------------------------------------------------------
#ifdef UCCFD
subroutine densityfield_V2(command, glook)
#else
subroutine densityfield(command, glook)
#endif
  !
  ! Alternative method that does not use the coarse grid
  !

  implicit none
  real       :: mp,xvtemp(3),dx1(3),dx2(3)
  integer    :: pp,pp_up,pp_down,index1(3),index2(3),command
  integer(1) :: glook

  time1 = mpi_wtime(ierr)

  if (command == 0) then
     mp = ncr**3 / real(np_groups_tot(glook))
  else if (command == 1) then
     mp = ncr**3 / real(np_groups_tot_dm(glook))
  else if (command == 2) then
     mp = ncr**3 / real(np_groups_tot_h(glook))
  endif
  
  !------------------------------------------------------------
  ! Initialize to zero 
  !------------------------------------------------------------

  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                 
  do kt = 0, nc_node_dim+1
     momden(:,:,kt) = 0.0
  enddo
  !$omp end parallel do

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
        pp_down = max_np + 1 - np_groups(g1) - np_buf_groups(g1)
        pp_up = max_np
     else if (command == 1) then
        pp_down = max_np_dm + 1 - np_groups_dm(g1) - np_buf_groups_dm(g1)
        pp_up = max_np_dm
     else if (command == 2) then
        pp_down = max_np_h + 1 - np_groups_h(g1) - np_buf_groups_h(g1)
        pp_up = max_np_h
     endif
  endif

  !------------------------------------------------------------
  ! Assign mass to grid cell
  !------------------------------------------------------------

  if (command == 0) then
     
     do pp = pp_down, pp_up
        xvtemp(1:3) = xvp(1:3,pp) - 0.5
        index1(1:3) = 1 + floor(xvtemp(1:3))
        index2(1:3) = 1 + index1(1:3)
        dx1(1:3) = real(index1(1:3)) - xvtemp(1:3)
        dx2(1:3) = 1.0 - dx1(1:3)
        
        if ( index1(1) >= 0 .and. index2(1) <= nc_node_dim + 1 .and. &
             index1(2) >= 0 .and. index2(2) <= nc_node_dim + 1 .and. &
             index1(3) >= 0 .and. index2(3) <= nc_node_dim + 1 ) then
           
           momden(index1(1),index1(2),index1(3)) = momden(index1(1),index1(2),index1(3)) +mp*dx1(1)*dx1(2)*dx1(3)
           momden(index1(1),index1(2),index2(3)) = momden(index1(1),index1(2),index2(3)) +mp*dx1(1)*dx1(2)*dx2(3)
           momden(index1(1),index2(2),index1(3)) = momden(index1(1),index2(2),index1(3)) +mp*dx1(1)*dx2(2)*dx1(3)
           momden(index1(1),index2(2),index2(3)) = momden(index1(1),index2(2),index2(3)) +mp*dx1(1)*dx2(2)*dx2(3)
           
           momden(index2(1),index1(2),index1(3)) = momden(index2(1),index1(2),index1(3)) +mp*dx2(1)*dx1(2)*dx1(3)
           momden(index2(1),index1(2),index2(3)) = momden(index2(1),index1(2),index2(3)) +mp*dx2(1)*dx1(2)*dx2(3)
           momden(index2(1),index2(2),index1(3)) = momden(index2(1),index2(2),index1(3)) +mp*dx2(1)*dx2(2)*dx1(3)
           momden(index2(1),index2(2),index2(3)) = momden(index2(1),index2(2),index2(3)) +mp*dx2(1)*dx2(2)*dx2(3)

        endif
     enddo !!pp
     
  else if (command == 1) then
     
     do pp = pp_down, pp_up
        xvtemp(1:3) = xvp_dm(1:3,pp) - 0.5
        index1(1:3) = 1 + floor(xvtemp(1:3))
        index2(1:3) = 1 + index1(1:3)
        dx1(1:3) = real(index1(1:3)) - xvtemp(1:3)
        dx2(1:3) = 1.0 - dx1(1:3)
        
        if ( index1(1) >= 0 .and. index2(1) <= nc_node_dim + 1 .and. &
             index1(2) >= 0 .and. index2(2) <= nc_node_dim + 1 .and. &
             index1(3) >= 0 .and. index2(3) <= nc_node_dim + 1 ) then
           
           momden(index1(1),index1(2),index1(3)) = momden(index1(1),index1(2),index1(3)) +mp*dx1(1)*dx1(2)*dx1(3)
           momden(index1(1),index1(2),index2(3)) = momden(index1(1),index1(2),index2(3)) +mp*dx1(1)*dx1(2)*dx2(3)
           momden(index1(1),index2(2),index1(3)) = momden(index1(1),index2(2),index1(3)) +mp*dx1(1)*dx2(2)*dx1(3)
           momden(index1(1),index2(2),index2(3)) = momden(index1(1),index2(2),index2(3)) +mp*dx1(1)*dx2(2)*dx2(3)
           
           momden(index2(1),index1(2),index1(3)) = momden(index2(1),index1(2),index1(3)) +mp*dx2(1)*dx1(2)*dx1(3)
           momden(index2(1),index1(2),index2(3)) = momden(index2(1),index1(2),index2(3)) +mp*dx2(1)*dx1(2)*dx2(3)
           momden(index2(1),index2(2),index1(3)) = momden(index2(1),index2(2),index1(3)) +mp*dx2(1)*dx2(2)*dx1(3)
           momden(index2(1),index2(2),index2(3)) = momden(index2(1),index2(2),index2(3)) +mp*dx2(1)*dx2(2)*dx2(3)

        endif
     enddo !!pp

  else if (command == 2) then
     
     do pp = pp_down, pp_up
        xvtemp(1:3) = xvmp_h(1:3,pp) - 0.5
        index1(1:3) = 1 + floor(xvtemp(1:3))
        index2(1:3) = 1 + index1(1:3)
        dx1(1:3) = real(index1(1:3)) - xvtemp(1:3)
        dx2(1:3) = 1.0 - dx1(1:3)
        
        if ( index1(1) >= 0 .and. index2(1) <= nc_node_dim + 1 .and. &
             index1(2) >= 0 .and. index2(2) <= nc_node_dim + 1 .and. &
             index1(3) >= 0 .and. index2(3) <= nc_node_dim + 1 ) then
           
           momden(index1(1),index1(2),index1(3)) = momden(index1(1),index1(2),index1(3)) +mp*dx1(1)*dx1(2)*dx1(3)
           momden(index1(1),index1(2),index2(3)) = momden(index1(1),index1(2),index2(3)) +mp*dx1(1)*dx1(2)*dx2(3)
           momden(index1(1),index2(2),index1(3)) = momden(index1(1),index2(2),index1(3)) +mp*dx1(1)*dx2(2)*dx1(3)
           momden(index1(1),index2(2),index2(3)) = momden(index1(1),index2(2),index2(3)) +mp*dx1(1)*dx2(2)*dx2(3)
           
           momden(index2(1),index1(2),index1(3)) = momden(index2(1),index1(2),index1(3)) +mp*dx2(1)*dx1(2)*dx1(3)
           momden(index2(1),index1(2),index2(3)) = momden(index2(1),index1(2),index2(3)) +mp*dx2(1)*dx1(2)*dx2(3)
           momden(index2(1),index2(2),index1(3)) = momden(index2(1),index2(2),index1(3)) +mp*dx2(1)*dx2(2)*dx1(3)
           momden(index2(1),index2(2),index2(3)) = momden(index2(1),index2(2),index2(3)) +mp*dx2(1)*dx2(2)*dx2(3)

        endif
     enddo !!pp
     
  endif !! command

  time2 = mpi_wtime(ierr)

  if (rank == 0) then
     write(*,*) 'Finished densityfield_V2 ...'
     write(*,*) 'Time elapsed : ', time2 - time1
     write(*,*) '--------------------------'
  endif
  return

#ifdef UCCFD
end subroutine densityfield_V2
#else
end subroutine densityfield
#endif

! -------------------------------------------------------------------------------------------------------

subroutine densitydarkmatter
  !
  ! Converts density field to density contrast field and apply Fourier transform
  !

  implicit none
  integer :: k,j,i
  real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,z_write
  real*8  :: dsum,dvar,dsumt,dvart

  time1 = mpi_wtime(ierr)

  !------------------------------------------------------------
  ! Assign data
  !------------------------------------------------------------

  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                
  do kt = 1, nc_node_dim
     cube(:,:,kt) = momden(1:nc_node_dim,1:nc_node_dim,kt)
  enddo
  !$omp end parallel do                                                                                                                          
  
  sum_dm_local=sum(cube)
  call mpi_reduce(sum_dm_local,sum_dm,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) print *,'Total mass=',sum_dm

  !------------------------------------------------------------
  ! Convert density field to delta field
  !------------------------------------------------------------

  dmin=0
  dmax=0
  dsum=0
  dvar=0
  
  do k=1,nc_node_dim
     do j=1,nc_node_dim
        do i=1,nc_node_dim
           cube(i,j,k)=cube(i,j,k)-1.0
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
  
  !------------------------------------------------------------    
  ! Forward FFT delta field
  !------------------------------------------------------------
  
  call cp_fftw(1)
  
  time2 = mpi_wtime(ierr)

  if (rank == 0) then

     dsum=dsumt/real(nc)**3
     dvar=sqrt(dvart/real(nc)**3)

     write(*,*) 'Finished densitydarkmatter ...'
     write(*,*) 'Cube min     :',dmint
     write(*,*) 'Cube max     :',dmaxt
     write(*,*) 'Delta sum    :',real(dsum,8)
     write(*,*) 'Delta var    :',real(dvar,8)
     write(*,*) 'Time elapsed : ', time2 - time1
     write(*,*) '--------------------------'
  endif
  return
  
end subroutine densitydarkmatter

!----------------------------------------------------------------------------------

subroutine veltransfer(command)
  !
  ! Apply a transfer function in Fourier space
  ! command determines the transfer function to apply
  ! (see below for more details)
  !

  implicit none
  integer                           :: command

  !! IT WOULD MAKE SENSE TO DEAL WITH THESE FILES IN A BETTER WAY...
  character(*),     parameter       :: fntf0 = 'nu_transfer_out_z0.dat'
  character(*),     parameter       :: fntf1 = 'nu_transfer_out_z1.dat'
  character(*),     parameter       :: fntf10 = 'nu_transfer_out_z10.dat' 
  character*30                      :: fntf
  
  character(*),     parameter       :: vfntf0 = 'nu_veltransfer_out_z0.dat'
  character(*),     parameter       :: vfntf1 = 'nu_veltransfer_out_z1.dat'
  character(*),     parameter       :: vfntf10 = 'nu_veltransfer_out_z10.dat' 
  character*30                      :: vfntf

  integer,          parameter       :: nk = 705, nucol = 6, dmcol = 2
  real, dimension(7,nk)             :: tf, vtf

  integer                           :: k, j, i, l, kg, jg, ig, mg, ioerr, ii, thread
  real(4)                           :: kz, ky, kx, kr, interpHTF_nu,interpMTF, kphys, interpVTF, interpVTF_nu
  integer                           :: dx, dxy, ind
  
  character*180                     :: fn
  character(len=6)                  :: rank_s
  real                              :: z_current, w1, w2, vsim2phys
  
  time1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) 'Starting veltransfer ...'

  !------------------------------------------------------------
  ! Assign file names         
  !------------------------------------------------------------

  z_current = z_checkpoint(cur_checkpoint)
  vsim2phys = 300. * sqrt(omega_m) * box * (1. + z_current) / 2. / nc

  if (z_current == 10.0) then
     fntf = fntf10
     vfntf = vfntf10
  else if (z_current == 1.0) then
     fntf = fntf1
     vfntf = vfntf1
  else
     fntf = fntf0
     vfntf = vfntf0
  endif

  !------------------------------------------------------------
  ! Read density transfer function         
  !------------------------------------------------------------
  
  if (rank ==0) then
     write(*,*) 'Reading ', fntf
     open(11, file = fntf)
     do k = 1,nk
        read(11,*) tf(1,k), tf(2,k), tf(3,k), tf(4,k), tf(5,k), tf(6,k), tf(7,k)
     end do
     close(11)
  end if

  call mpi_bcast(tf, 7*nk, mpi_real, 0, mpi_comm_world, ierr)

  !------------------------------------------------------------
  ! Read velocity transfer function                                                                                                                      
  !------------------------------------------------------------

  if (rank == 0) then
     write(*,*) 'Reading ',vfntf
     open(11,file=vfntf)
     do k=1,nk
        read(11,*) vtf(1,k),vtf(2,k),vtf(3,k),vtf(4,k),vtf(5,k),vtf(6,k),vtf(7,k)
     end do
     close(11)

     !------------------------------------------------------------
     ! Put vtf in simulation's unit
     !------------------------------------------------------------

     do k=1,nk
        do j = 2,7
           vtf(j,k) = vtf(j,k) / vsim2phys
        enddo
     enddo

  endif

  call mpi_bcast(vtf, 7*nk, mpi_real, 0, mpi_comm_world, ierr)
    
  !------------------------------------------------------------
  ! Multiply slab by the appropriate transfer function                                                                                              
  !------------------------------------------------------------

#ifndef SLAB
  dx  = fsize(1)
  dxy = dx * fsize(2)
  ind = 0
#endif
#ifdef SLAB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                 
!! DON'T USE THIS
!! IT DOESN'T WORK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k = 1, nc_slab
     kg=k+nc_slab*rank
     if (kg .lt. hc+2) then
        kz=kg-1
     else
        kz=kg-1-nc
     endif
     kz=2*sin(pi*kz/ncr)
     do j = 1, nc
        if (j .lt. hc+2) then
           ky=j-1
        else
           ky=j-1-nc
        endif
        ky=2*sin(pi*ky/ncr)
        do i = 1, nc+2, 2
           kx=(i-1)/2
           kx=2*sin(pi*kx/ncr)
#else
  !$omp parallel default(shared) private(k, j, i, kg, mg, jg, ig, ind, kz, ky, kx, kphys, kr, interpMTF, interpHTF_nu, interpVTF, interpVTF_nu, w1, w2, ii) 
  !!!omp parallel default(private) shared(nc_pen, mypadd, nc_node_dim, nc, dxy, dx, hc, box, pi, ncr, tf, nk, slab, slab2)
  !$omp do schedule(dynamic) 
  do k = 1, nc_pen+mypadd
     ind = (k-1)*nc_node_dim*nc/2
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
           
           kphys = 2*pi*sqrt(kx**2+ky**2+kz**2)/box
           
           kx = 2*sin(pi*kx/ncr)
           ky = 2*sin(pi*ky/ncr)
           kz = 2*sin(pi*kz/ncr)
           
#endif
           kr = sqrt(kx**2+ky**2+kz**2)

           !------------------------------------------------------------
           ! Find appropriate velocity transfer function
           !------------------------------------------------------------
           
           if (kphys <= vtf(1,1)) then
              interpVTF = vtf(dmcol, 1)
              interpVTF_nu = vtf(nucol, 1)
           else if (kphys >= vtf(1, nk)) then
              interpVTF = vtf(dmcol, nk)
              interpVTF_nu = vtf(nucol, nk)
           else
              do ii = 1, nk-1
                 if (kphys <= vtf(1,ii+1)) then
                    w1 = vtf(1,ii+1) - kphys
                    w2 = kphys - vtf(1,ii)
                    interpVTF = (w1*vtf(dmcol,ii)+w2*vtf(dmcol,ii+1))/(vtf(1,ii+1)-vtf(1,ii))
                    interpVTF_nu = (w1*vtf(nucol,ii)+w2*vtf(nucol,ii+1))/(vtf(1,ii+1)-vtf(1,ii))
                    exit
                 endif
              enddo
           endif
           
           !------------------------------------------------------------
           ! Find appropriate density transfer function
           !------------------------------------------------------------

           if (kphys <= tf(1,1)) then
              interpMTF = tf(dmcol, 1)
              interpHTF_nu = tf(nucol, 1)
           else if (kphys >= tf(1, nk)) then
              interpMTF = tf(dmcol, nk)
              interpHTF_nu = tf(nucol, nk)
           else
              do ii = 1, nk-1
                 if (kphys <= tf(1,ii+1)) then
                    w1 = tf(1,ii+1) - kphys
                    w2 = kphys - tf(1,ii)
                    interpMTF = (w1*tf(dmcol,ii)+w2*tf(dmcol,ii+1))/(tf(1,ii+1)-tf(1,ii))
                    interpHTF_nu = (w1*tf(nucol,ii)+w2*tf(nucol,ii+1))/(tf(1,ii+1)-tf(1,ii))
                    exit
                 endif
              enddo
           endif
           
           if (kr == 0) then
              kr = 0.000001
           endif
           
           !------------------------------------------------------------
           ! Apply transfer function
           !------------------------------------------------------------

           if (command == 0) then !! dm density => dm velocity potential                                                                       
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpVTF/interpMTF / kr
              
           else if (command == 1) then !! dm velocity potential => nu velocity potential                                                         
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpMTF/interpHTF_nu * interpVTF_nu/interpVTF
              
           else if (command == 2) then !! dm density => nu velocity potential                                                                     
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpVTF_nu/interpMTF / kr
              
           else if (command == 3) then !! dm density => dm velocity
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpVTF/interpMTF
              
           else if (command == 4) then !! dm density => nu velocity
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpVTF_nu/interpMTF
              
           else if (command == 5) then !! dm velocity => nu velocity
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpVTF_nu/interpVTF
              
           else if (command == 6) then !! nu density => nu velocity
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * interpVTF_nu/interpHTF_nu
              
           else if (command == 7) then !! nu velocity => (dm-nu) relative velocity
              slab(i:i+1, j, k) = (-1) * slab(i:i+1, j, k) * (interpVTF - interpVTF_nu)/interpVTF_nu
              
           end if
        enddo !! i
     enddo !! j
  enddo !! k
#ifndef SLAB
  !$omp end do                                                                                                                                          
  !$omp end parallel                                                                                                                                
#endif

  time2 = mpi_wtime(ierr)
  
  if (rank == 0) then
     write(*,*) 'Finished veltransfer ...'
     write(*,*) 'command      : ', command
     write(*,*) 'Time elapsed : ', time2 - time1
     write(*,*) '--------------------------'
  endif

end subroutine veltransfer

! -------------------------------------------------------------------------------------------------------                                               

subroutine numerical_gradient(dir, command)
  !
  ! Computes numerical gradient in direction dir
  !

  implicit none
  integer ::  dir, command
  integer :: kt, k, j, i 

  ! slab contains dm or nu velocity potential (after veltransfer)                                                                                 
  time2 = mpi_wtime(ierr)
  
  !------------------------------------------------------------
  ! Fourier transform backward slab to cube
  !------------------------------------------------------------                                                        
  call cp_fftw(-1)

  !------------------------------------------------------------
  ! Stores cube in momden to compute gradient                                                      
  !------------------------------------------------------------

  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
  do kt = 0, nc_node_dim+1
     momden(:,:,kt) = 0.0
  enddo
  !$omp end parallel do

  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
  do kt = 1, nc_node_dim
     momden(1:nc_node_dim,1:nc_node_dim,kt) = cube(:,:,kt)
  enddo
  !$omp end parallel do

  call writebuffer_momden

  !------------------------------------------------------------
  ! Numerical gradient
  !------------------------------------------------------------

  if (dir == 1) then
     !------------------------------------------------------------
     ! Along X
     !------------------------------------------------------------
     !$omp parallel do num_threads(nt) default(shared) private(k,j,i)                                                                                
     do k = 1,nc_node_dim
        do j = 1,nc_node_dim
           do i = 1,nc_node_dim
              cube(i,j,k) = (momden(i+1,j,k)-momden(i-1,j,k))/2.
           enddo
        enddo
     enddo
     !$omp end parallel do                                                                                                                          
  else if (dir == 2) then
     !------------------------------------------------------------
     ! Along Y
     !------------------------------------------------------------
     !$omp parallel do num_threads(nt) default(shared) private(k,j,i)                                                                               
     do k = 1,nc_node_dim
        do j = 1,nc_node_dim
           do i = 1,nc_node_dim
              cube(i,j,k) = (momden(i,j+1,k)-momden(i,j-1,k))/2.
           enddo
        enddo
     enddo
     !$omp end parallel do                                                                                                                          
  else if (dir == 3) then
     !------------------------------------------------------------
     ! Along Z, in 2 steps for threading
     !------------------------------------------------------------
     !$omp parallel do num_threads(nt) default(shared) private(k,j,i)                                                                                 
     do k = 1,nc_node_dim
        do j = 1,nc_node_dim
           do i = 1,nc_node_dim
              cube(i,j,k) = momden(i,j,k+1)/2.
           enddo
        enddo
     enddo
     !$omp end parallel do                                                                                                                           
     !$omp parallel do num_threads(nt) default(shared) private(k,j,i)                                                                                   
     do k = 1,nc_node_dim
        do j = 1,nc_node_dim
           do i = 1,nc_node_dim
              cube(i,j,k) = cube(i,j,k) - momden(i,j,k-1)/2.
           enddo
        enddo
     enddo
     !$omp end parallel do                                                                                                                           
  endif

  !------------------------------------------------------------
  ! Writes velocity field to file if needed
  !------------------------------------------------------------
#ifdef write_vel
  if (command >= 0) then
     call writevelocityfield(command, dir, 1)
  endif
#endif

  !------------------------------------------------------------
  ! Fourier transform cube to slab
  !------------------------------------------------------------
  call cp_fftw(1)
  
  time2 = mpi_wtime(ierr)
  
  if (rank == 0) then
     write(*,*) 'Finished Numerical_gradient ...'
     write(*,*) 'Time elapsed : ', time2 - time1
     write(*,*) '--------------------------'
  endif

end subroutine numerical_gradient

! -------------------------------------------------------------------------------------------------------

subroutine read_particles_files
  !
  ! Read all particles files
  !

  implicit none

  !--------------------------------------------------------------------------------                                                            
  ! Neutrinos                                                      
  !--------------------------------------------------------------------------------
     
  call read_particles(0)
  call pass_particles(0)
  call order_xvp_groups(0)
  call buffer_particles_groups(0,g0)
  call buffer_particles_groups(0,g1)
  call order_xvp_ll(0,g0,.false.)
  call order_xvp_ll(0,g1,.false.)

  !--------------------------------------------------------------------------------                                                            
  ! Dark matter                                                      
  !--------------------------------------------------------------------------------
     
  call read_particles(1)
  call pass_particles(1)
  call order_xvp_groups(1)
  call buffer_particles_groups(1,g0)
  call buffer_particles_groups(1,g1)
  call order_xvp_ll(1,g0,.false.)
  call order_xvp_ll(1,g1,.false.)

  !--------------------------------------------------------------------------------                                                            
  ! Halos                                                      
  !--------------------------------------------------------------------------------
     
  call read_particles(2)
  call pass_particles(2)
  call order_xvp_groups(2)
  call buffer_particles_groups(2,g0)
  call buffer_particles_groups(2,g1)
  call order_xvp_ll(2,g0,.false.)
  call order_xvp_ll(2,g1,.false.)
  
  
end subroutine read_particles_files

! -------------------------------------------------------------------------------------------------------

subroutine auto_power_loop(component)
  !
  ! Auto-power loop
  !

  implicit none
  integer component

  if (rank == 0) then
     timeCheckpoint = mpi_wtime(ierr)
     write(*,*) '****************************************************'
     write(*,*) 'TIME = ', timeCheckpoint - timeStart
     write(*,*) '****************************************************'
     write(*,*) 'AUTO-POWER FOR COMPONENT : ', component
     write(*,*) '****************************************************'
  endif
  
  !--------------------------------------------------------------------------------                                                            
  ! NEUTRINOS                                                                                                                                   
  !--------------------------------------------------------------------------------
  ! g0 nu --> slab --> slab2
  ! g1 nu --> slab
  !--------------------------------------------------------------------------------

  if (computePS(0) .or. computePS(3)) then
     call velocity_density_V2(0,g0,N_closest_auto_nu,component)
     call darkmatter

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                            
  endif

  if (computePS(0)) then
     call velocity_density_V2(0,g1,N_closest_auto_nu,component)
     call darkmatter

  !--------------------------------------------------------------------------------
  ! g1 nu x g0 nu
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab2, pk_all(component,:,:,0))
     if (rank == 0) call writepowerspectra(0,1)
  endif

  !--------------------------------------------------------------------------------                                                            
  ! DARK MATTER                                                                                                                                   
  !--------------------------------------------------------------------------------
  ! g0 dm --> slab
  ! g0 dm - nu --> slab3
  ! g0 dm --> slab2
  ! g1 dm --> slab
  !--------------------------------------------------------------------------------
  
  if (computePS(1) .or. computePS(3)) then
     call velocity_density_V2(1,g0,N_closest_auto_dm,component)
     call darkmatter
  endif

  if (computePS(3)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                  
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt) - slab2(:,:,kt)
     enddo
     !$omp end parallel do
  endif

  if (computePS(1)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                     
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                             
  endif

  if (computePS(1) .or. computePS(3)) then
     call velocity_density_V2(1,g1,N_closest_auto_dm,component)
     call darkmatter
  endif

  !--------------------------------------------------------------------------------
  ! g1 dm x g0 dm
  !--------------------------------------------------------------------------------
  if (computePS(1)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,1))
     if (rank == 0) call writepowerspectra(1,1)
  endif

  !--------------------------------------------------------------------------------                                                                
  ! DARK MATTER - NEUTRINO RELATIVE VELOCITY                                                                                                       
  !--------------------------------------------------------------------------------
  ! g1 dm --> slab2
  ! g1 nu --> slab
  ! g1 dm - nu --> slab2
  !--------------------------------------------------------------------------------
       
  if (computePS(3)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                    
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                          
  endif

  if (computePS(3) .or. computePS(20)) then
     call velocity_density_V2(0, g1, N_closest_auto_nu, component)
     call darkmatter
  endif

  if (computePS(3)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                  
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab2(:,:,kt) - slab(:,:,kt)
     enddo
     !$omp end parallel do

  !--------------------------------------------------------------------------------
  ! g1 dm - nu x g0 dm - nu
  !--------------------------------------------------------------------------------
     call powerspectrum(slab2, slab3, pk_all(component,:,:,3))
     if (rank == 0) call writepowerspectra(3,1)
  endif

  !--------------------------------------------------------------------------------
  ! g1 nu --> slab3
  !--------------------------------------------------------------------------------

  if (computePS(20)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do 
  endif

  !--------------------------------------------------------------------------------                                                            
  ! HALOS                                                                                                                                   
  !--------------------------------------------------------------------------------
  ! g0 h --> slab --> slab2
  ! g1 h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(4) .or. computePS(20)) then
     call velocity_density_V2(2, g0, N_closest_auto_h,component)
     call darkmatter

     !$omp parallel do num_threads(nt) default(shared) private(kt)
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do 

     call velocity_density_V2(2, g1, N_closest_auto_h,component)
     call darkmatter
  endif

  !--------------------------------------------------------------------------------
  ! g1 h x g0 h
  !--------------------------------------------------------------------------------
  if (computePS(4)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,4))
     if (rank == 0) call writepowerspectra(4,1)
  endif

  !--------------------------------------------------------------------------------
  ! HALO - NEUTRINOS RELATIVE VELOCITY                                                                                                          
  !--------------------------------------------------------------------------------
  ! g1 h - nu --> slab3
  ! g0 nu --> slab
  ! g0 h - nu --> slab2
  !--------------------------------------------------------------------------------

  if (computePS(20)) then

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                       
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt) - slab3(:,:,kt)
     enddo
     !$omp end parallel do

     call velocity_density_V2(0, g0, N_closest_auto_nu, component)
     call darkmatter

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                       
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab2(:,:,kt) - slab(:,:,kt)
     enddo
     !$omp end parallel do

  !--------------------------------------------------------------------------------
  ! g0 h - nu x g1 h - nu
  !--------------------------------------------------------------------------------
     call powerspectrum(slab2, slab3, pk_all(component,:,:,20))
     if (rank == 0) call writepowerspectra(20,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^DM                                                                   
  !--------------------------------------------------------------------------------
  ! g0 dm^dm --> slab --> slab3
  ! g1 dm^dm --> slab
  !--------------------------------------------------------------------------------

#ifdef compute_denPS
  if (component == 1) then
     call densityfield(1,g0)
     call densitydarkmatter
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        
     call densityfield(1,g1)
     call densitydarkmatter
     call powerspectrum(slab, slab3, pk_all(component,:,:,6))
     call writedenpowerspectra(6)
  endif
#endif


  if (computePS(6) .or. computePS(17)) then
     call densityfield(1, g0)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        
  endif

  if (computePS(6)) then
     call densityfield(1, g1)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)

  !--------------------------------------------------------------------------------
  ! g1 dm^dm x g0 dm^dm
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(component,:,:,6))
     if (rank == 0) call writepowerspectra(6,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^DM
  !--------------------------------------------------------------------------------
  ! g1 nu^dm --> slab --> slab2
  ! g0 nu^dm --> slab
  ! g0 dm^dm - nu^dm --> slab3
  !--------------------------------------------------------------------------------

  if (computePS(9) .or. computePS(17)) then
     call densityfield(1, g1)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, -1)

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do            

     call densityfield(1, g0)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, -1)
  endif
  
  if (computePS(17)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab3(:,:,kt) - slab(:,:,kt)
     enddo
     !$omp end parallel do
  endif
      
  !--------------------------------------------------------------------------------
  ! g0 nu^dm x g1 nu^dm
  !--------------------------------------------------------------------------------
  if (computePS(9)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,9))
     if (rank == 0) call writepowerspectra(9,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^DM - NU^DM RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! g1 dm^dm --> slab
  ! g1 dm^dm - nu^dm --> slab
  !--------------------------------------------------------------------------------

  if (computePS(17)) then
     call densityfield(1, g1)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab(:,:,kt) = slab(:,:,kt) - slab2(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                   

  !--------------------------------------------------------------------------------
  ! g1 dm^dm - nu^dm x g0 dm^dm - nu^dm
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(component,:,:,17))
     if (rank == 0) call writepowerspectra(17,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^H
  !--------------------------------------------------------------------------------
  ! g0 dm^h --> slab --> slab3
  ! g1 dm^h --> slab
  !--------------------------------------------------------------------------------

#ifdef compute_denPS
  if (component == 1) then
     call densityfield(2,g0)
     call densitydarkmatter
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        
     call densityfield(2,g1)
     call densitydarkmatter
     call powerspectrum(slab, slab3, pk_all(component,:,:,7))
     call writedenpowerspectra(7)
  endif
#endif

  if (computePS(7) .or. computePS(16)) then
     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)
  
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do           
  endif

  if (computePS(7)) then
     call densityfield(2, g1)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)
     
  !--------------------------------------------------------------------------------
  ! g1 dm^h x g0 dm^h
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(component,:,:,7))
     if (rank == 0) call writepowerspectra(7,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^H
  !--------------------------------------------------------------------------------
  ! g1 nu^h --> slab --> slab2
  ! g0 nu^h --> slab
  ! g0 dm^h - nu^h --> slab3
  !--------------------------------------------------------------------------------

  if (computePS(8) .or. computePS(16)) then
     call densityfield(2, g1)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, -1)
  
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do

     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, -1)
  endif
     
  if (computePS(16)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab3(:,:,kt) - slab(:,:,kt)
     enddo
     !$omp end parallel do 
  endif

  !--------------------------------------------------------------------------------
  ! g0 nu^h x g1 nu^h
  !--------------------------------------------------------------------------------
  if (computePS(8)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,8))
     if (rank == 0) call writepowerspectra(8,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^H - NU^H RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! g1 dm^h --> slab
  ! g1 dm^h - nu^h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(16)) then
     call densityfield(2, g1)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab(:,:,kt) = slab(:,:,kt) - slab2(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                   

  !--------------------------------------------------------------------------------
  ! g1 dm^h - nu^h x g0 dm^h - nu^h
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(component,:,:,16))
     if (rank == 0) call writepowerspectra(16,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^NU
  !--------------------------------------------------------------------------------
  ! g0 nu^nu --> slab --> slab3
  ! g1 nu^nu --> slab
  !--------------------------------------------------------------------------------

#ifdef compute_denPS
  if (component == 1) then
     call densityfield(0,g0)
     call densitydarkmatter
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        
     call densityfield(0,g1)
     call densitydarkmatter
     call powerspectrum(slab, slab3, pk_all(component,:,:,5))
     call writedenpowerspectra(5)
  endif
#endif


  if (computePS(5)) then
     call densityfield(1, g0)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        
     
     call densityfield(1, g1)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)

  !--------------------------------------------------------------------------------
  ! g1 nu^nu x g0 nu^nu
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(component,:,:,5))
     if (rank == 0) call writepowerspectra(5,1)
  endif

end subroutine auto_power_loop

! -------------------------------------------------------------------------------------------------------

subroutine auto_power_vtf
  !
  ! Auto-power, linear, without gradient 
  !
  implicit none

  if (rank == 0) then
     timeCheckpoint = mpi_wtime(ierr)
     write(*,*) '****************************************************'
     write(*,*) 'TIME = ', timeCheckpoint - timeStart
     write(*,*) '****************************************************'
     write(*,*) 'VTF AUTO-POWER'
     write(*,*) '****************************************************'
  endif

  !--------------------------------------------------------------------------------
  ! DM^DM
  !--------------------------------------------------------------------------------
  ! g0 dm^dm --> slab --> slab3
  ! g1 dm^dm --> slab
  !--------------------------------------------------------------------------------

  if (computePS(24) .or. computePS(28) .or. computePS(29)) then
     call densityfield(1, g0)
     call densitydarkmatter
     call veltransfer(3)
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        

     call densityfield(1, g1)
     call densitydarkmatter
     call veltransfer(3)
  endif

  !--------------------------------------------------------------------------------
  ! g1 dm^dm x g0 dm^dm
  !--------------------------------------------------------------------------------
  if (computePS(24)) then
     call powerspectrum(slab, slab3, pk_all(1,:,:,24))
     if (rank == 0) call writepowerspectra(24,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^DM
  !--------------------------------------------------------------------------------
  ! g1 nu^dm --> slab --> slab2
  ! g0 nu^dm --> slab
  !--------------------------------------------------------------------------------

  if (computePS(28) .or. computePS(29)) then
     call veltransfer(5)
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab(:,:,kt) = slab3(:,:,kt)
     enddo
     !$omp end parallel do
     call veltransfer(5)
  endif

  !--------------------------------------------------------------------------------
  ! g0 nu^dm x g1 nu^dm
  !--------------------------------------------------------------------------------
  if (computePS(28)) then
     call powerspectrum(slab, slab2, pk_all(1,:,:,28))
     if (rank == 0) call writepowerspectra(28,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^DM - NU^DM RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! g0 dm^dm - nu^dm --> slab --> slab3
  ! g1 dm^dm - nu^dm --> slab
  !--------------------------------------------------------------------------------

  if (computePS(29)) then
     call veltransfer(7)
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab(:,:,kt) = slab2(:,:,kt)
     enddo
     !$omp end parallel do
     call veltransfer(7)

  !--------------------------------------------------------------------------------
  ! g1 dm^dm - nu^dm x g0 dm^dm - nu^dm
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(1,:,:,29))
     if (rank == 0) call writepowerspectra(29,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^H
  !--------------------------------------------------------------------------------
  ! g0 dm^h --> slab --> slab3
  ! g1 dm^h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(26) .or. computePS(30) .or. computePS(31)) then
     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(3)
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        

     call densityfield(2, g1)
     call densitydarkmatter
     call veltransfer(3)
  endif

  !--------------------------------------------------------------------------------
  ! g1 dm^h x g0 dm^h
  !--------------------------------------------------------------------------------
  if (computePS(26)) then
     call powerspectrum(slab, slab3, pk_all(1,:,:,26))
     if (rank == 0) call writepowerspectra(26,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^H
  !--------------------------------------------------------------------------------
  ! g1 nu^h --> slab --> slab2
  ! g0 nu^h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(30) .or. computePS(31)) then
     call veltransfer(5)
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab(:,:,kt) = slab3(:,:,kt)
     enddo
     !$omp end parallel do
     call veltransfer(5)
  endif

  !--------------------------------------------------------------------------------
  ! g0 nu^h x g1 nu^h
  !--------------------------------------------------------------------------------
  if (computePS(30)) then
     call powerspectrum(slab, slab2, pk_all(1,:,:,30))
     if (rank == 0) call writepowerspectra(30,1)
  endif
 
  !--------------------------------------------------------------------------------
  ! DM^H - NU^H RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! g0 dm^h - nu^h --> slab --> slab3
  ! g1 dm^h - nu^h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(31)) then
     call veltransfer(7)
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab(:,:,kt) = slab2(:,:,kt)
     enddo
     !$omp end parallel do
     call veltransfer(7)
     
  !--------------------------------------------------------------------------------
  ! g1 dm^h - nu^h x g0 dm^h - nu^h
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(1,:,:,31))
     if (rank == 0) call writepowerspectra(31,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^NU
  !--------------------------------------------------------------------------------
  ! g0 nu^nu --> slab --> slab3
  ! g1 nu^nu --> slab
  !--------------------------------------------------------------------------------

  if (computePS(23)) then
     call densityfield(0, g0)
     call densitydarkmatter
     call veltransfer(6)
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                        
     
     call densityfield(0, g1)
     call densitydarkmatter
     call veltransfer(6)
     
  !--------------------------------------------------------------------------------
  ! g1 nu^nu x g0 nu^nu
  !--------------------------------------------------------------------------------
     call powerspectrum(slab, slab3, pk_all(1,:,:,23))
     if (rank == 0) call writepowerspectra(23,1)
  endif
  
end subroutine auto_power_vtf

! -------------------------------------------------------------------------------------------------------

subroutine cross_power_loop(component)
  !
  ! Cross-power loop
  !
  implicit none
  integer component

  if (rank == 0) then
     timeCheckpoint = mpi_wtime(ierr)
     write(*,*) '****************************************************'
     write(*,*) 'TIME = ', timeCheckpoint - timeStart
     write(*,*) '****************************************************'
     write(*,*) 'CROSS-POWER LOOP FOR COMPONENT : ', component
     write(*,*) '****************************************************'
  endif

  !--------------------------------------------------------------------------------
  ! NEUTRINOS
  !--------------------------------------------------------------------------------
  ! nu --> slab --> slab2
  !--------------------------------------------------------------------------------

  if (computePS(2) .or. computePS(22) .or. computePS(32) .or. computePS(14) &
       .or. computePS(15) .or. computePS(19) .or. computePS(18)) then 

     call velocity_density_V2(0,g0,N_closest_nu,component)
#ifdef write_vel
     call writevelocityfield(0, component, 0)
#endif
     call darkmatter

     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                           
  endif

  !--------------------------------------------------------------------------------
  ! DARK MATTER
  !--------------------------------------------------------------------------------
  ! dm --> slab --> slab3
  !--------------------------------------------------------------------------------

  if (computePS(2) .or. computePS(10) .or. computePS(12) .or. computePS(33) &
       .or. computePS(13) .or. computePS(19) .or. computePS(18)) then

     call velocity_density_V2(1,g0,N_closest_dm,component)
#ifdef write_vel
     call writevelocityfield(1, component, 0)
#endif
     call darkmatter
  
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                      
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do                                                                                                                 
  endif

  !--------------------------------------------------------------------------------
  ! HALOS
  !--------------------------------------------------------------------------------
  ! h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(10) .or. computePS(22)) then

     call velocity_density_V2(2,g0,N_closest_h,component)
#ifdef write_vel
     call writevelocityfield(2, component, 0)
#endif
     call darkmatter
  endif

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! dm x nu
  ! dm x halo
  ! nu x halo
  !--------------------------------------------------------------------------------

  if (computePS(2)) then
     call powerspectrum(slab3, slab2, pk_all(component,:,:, 2))
     if (rank == 0) call writepowerspectra(2,1)
  endif
  if (computePS(10)) then
     call powerspectrum(slab3, slab, pk_all(component,:,:, 10))
     if (rank == 0) call writepowerspectra(10,1)
  endif
  if (computePS(22)) then
     call powerspectrum(slab2, slab, pk_all(component,:,:, 22))
     if (rank == 0) call writepowerspectra(22,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^H
  !--------------------------------------------------------------------------------
  ! dm^h --> slab
  !--------------------------------------------------------------------------------
  
  if (computePS(12) .or. computePS(32)) then
     call densityfield(2, g0)
#ifdef write_den
     if (component == 1) then
        call writedensityfield(2,0)
     endif
#endif
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, 4)
  endif

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! dm^h x dm
  ! dm^h x nu
  !--------------------------------------------------------------------------------

  if (computePS(12)) then
     call powerspectrum(slab, slab3, pk_all(component,:,:,12))
     if (rank == 0) call writepowerspectra(12,1)
  endif
  if (computePS(32)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,32))
     if (rank == 0) call writepowerspectra(32,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^H
  !--------------------------------------------------------------------------------
  ! nu^h --> slab
  !--------------------------------------------------------------------------------

  if (computePS(14) .or. computePS(33)) then
     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, 6)
  endif

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! nu^h x dm
  ! nu^h x nu
  !--------------------------------------------------------------------------------

  if (computePS(33)) then
     call powerspectrum(slab, slab3, pk_all(component,:,:,33))
     if (rank == 0) call writepowerspectra(33,1)
  endif
  if (computePS(14)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,14))
     if (rank == 0) call writepowerspectra(14,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^DM
  !--------------------------------------------------------------------------------
  ! dm^dm --> slab
  !--------------------------------------------------------------------------------

  if (computePS(13)) then
     call densityfield(1, g0)
#ifdef write_den
     if (component == 1) then
        call writedensityfield(1,0)
     endif
#endif
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, 3)

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! dm^dm x dm
  !--------------------------------------------------------------------------------

     call powerspectrum(slab, slab3, pk_all(component,:,:,13))
     if (rank == 0) call writepowerspectra(13,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^DM
  !--------------------------------------------------------------------------------
  ! nu^dm --> slab
  !--------------------------------------------------------------------------------

  if (computePS(15) .or. computePS(19)) then
     call densityfield(1, g0)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, 5)

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! nu^dm x nu
  !--------------------------------------------------------------------------------

     call powerspectrum(slab, slab2, pk_all(component,:,:,15))
     if (rank == 0) call writepowerspectra(15,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM - NU RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! dm - nu --> slab3
  !--------------------------------------------------------------------------------

  if (computePS(18) .or. computePS(19)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                            
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab3(:,:,kt) - slab2(:,:,kt)
     enddo
     !$omp end parallel do
  endif

  !--------------------------------------------------------------------------------
  ! DM^DM - NU^DM RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! nu^dm --> slab2
  ! dm^dm --> slab
  ! dm^dm - nu^dm --> slab2
  !--------------------------------------------------------------------------------

  if (computePS(19)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                            
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do
     
     call densityfield(1, g0)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt) - slab2(:,:,kt)
     enddo
     !$omp end parallel do

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! (dm^dm - nu^dm) x (dm - nu)
  !--------------------------------------------------------------------------------

     call powerspectrum(slab2, slab3, pk_all(component,:,:,19))
     if (rank == 0) call writepowerspectra(19,1)
  endif

  !--------------------------------------------------------------------------------
  ! NU^H
  !--------------------------------------------------------------------------------
  ! nu^h --> slab --> slab2
  !--------------------------------------------------------------------------------

  if (computePS(18) .or. computePS(21) .or. computePS(34)) then
     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(2)
     call numerical_gradient(component, -1)

     !$omp parallel do num_threads(nt) default(shared) private(kt)
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do

  !--------------------------------------------------------------------------------
  ! DM^H
  !--------------------------------------------------------------------------------
  ! dm^h --> slab
  !--------------------------------------------------------------------------------

     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)
  endif

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! dm^h x nu^h
  !--------------------------------------------------------------------------------

  if (computePS(34)) then
     call powerspectrum(slab, slab2, pk_all(component,:,:,34))
     if (rank == 0) call writepowerspectra(34,1)
  endif

  !--------------------------------------------------------------------------------
  ! DM^H - NU^H RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! dm^h - nu^h --> slab2
  !--------------------------------------------------------------------------------

  if (computePS(18) .or. computePS(21)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                             
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt) - slab2(:,:,kt)
     enddo
     !$omp end parallel do
  endif

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! (dm^h - nu^h) x (dm - nu)
  !--------------------------------------------------------------------------------

  if (computePS(18)) then
     call powerspectrum(slab2, slab3, pk_all(component,:,:,18))
     if (rank == 0) call writepowerspectra(18,1)
  endif

  !--------------------------------------------------------------------------------
  ! HALO - NU RELATIVE VELOCITY
  !--------------------------------------------------------------------------------
  ! nu --> slab --> slab3
  ! h --> slab
  ! h - nu --> slab3
  !--------------------------------------------------------------------------------

  if (computePS(21)) then
     call velocity_density_V2(0,g0,N_closest_nu,component)
     call darkmatter
     
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                             
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt)
     enddo
     !$omp end parallel do
  endif

  if (computePS(21) .or. computePS(11)) then
     call velocity_density_V2(2,g0,N_closest_h,component)
     call darkmatter
  endif

  if (computePS(21)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                     
     do kt = 1, kt_stop
        slab3(:,:,kt) = slab(:,:,kt) - slab3(:,:,kt) 
     enddo
     !$omp end parallel do

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! (dm^h - nu^h) x (h - nu)
  ! h --> slab2
  !--------------------------------------------------------------------------------

     call powerspectrum(slab2, slab3, pk_all(component,:,:,21))
     if (rank == 0) call writepowerspectra(21,1)
  endif

  if (computePS(11)) then
     !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                     
     do kt = 1, kt_stop
        slab2(:,:,kt) = slab(:,:,kt) 
     enddo
     !$omp end parallel do

  !--------------------------------------------------------------------------------
  ! DM^H
  !--------------------------------------------------------------------------------
  ! dm^h --> slab
  !--------------------------------------------------------------------------------

     call densityfield(2, g0)
     call densitydarkmatter
     call veltransfer(0)
     call numerical_gradient(component, -1)

  !--------------------------------------------------------------------------------
  ! Cross PS
  !--------------------------------------------------------------------------------
  ! h x dm^h
  !--------------------------------------------------------------------------------

     call powerspectrum(slab2, slab, pk_all(component,:,:,11))
     if (rank == 0) call writepowerspectra(11,1)
  endif

end subroutine cross_power_loop

! -------------------------------------------------------------------------------------------------------

subroutine cross_power_density
  !
  ! Loop for xpower spectra for density
  !
  implicit none

  if (rank == 0) then
     timeCheckpoint = mpi_wtime(ierr)
     write(*,*) '****************************************************'
     write(*,*) 'TIME = ', timeCheckpoint - timeStart
     write(*,*) '****************************************************'
     write(*,*) 'CROSS-POWER LOOP FOR DENSITY'
     write(*,*) '****************************************************'
  endif

  !--------------------------------------------------------------------------------
  ! Neutrinos
  !--------------------------------------------------------------------------------

  call densityfield(0,g0)
  call densitydarkmatter
  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                          
  do kt = 1, kt_stop
     slab3(:,:,kt) = slab(:,:,kt)
  enddo
  !$omp end parallel do                                                                                                                                  

  !--------------------------------------------------------------------------------
  ! Dark matter
  !--------------------------------------------------------------------------------

  call densityfield(1,g0)
  call densitydarkmatter
  !$omp parallel do num_threads(nt) default(shared) private(kt)                                                                                          
  do kt = 1, kt_stop
     slab2(:,:,kt) = slab(:,:,kt)
  enddo
  !$omp end parallel do                                                                                                                                  

  !--------------------------------------------------------------------------------
  ! Halos
  !--------------------------------------------------------------------------------

  call densityfield(2,g0)
  call densitydarkmatter


  !--------------------------------------------------------------------------------
  ! Cross power spectra
  !--------------------------------------------------------------------------------
  ! dm x nu
  ! dm x halo
  !--------------------------------------------------------------------------------

  call powerspectrum(slab2, slab3, pk_all(1,:,:,2))
  call writedenpowerspectra(2)

  call powerspectrum(slab2, slab, pk_all(1,:,:,10))
  call writedenpowerspectra(10)


end subroutine cross_power_density

end program cic_velpower_halo
