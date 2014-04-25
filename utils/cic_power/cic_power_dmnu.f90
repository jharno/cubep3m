!!
!! cic_power_dmnu.f90
!!
!! Program to compute the power spectra of both dark matter and neutrinos as well as their
!! cross power spectra. This is to be used in conjunction with the particle checkpoint files
!! produced by a -DNEUTRINOS cubep3m simulation.
!!  
!! * Using FFTW on the SciNet GPC compile with:
!!   mpif90 -shared-intel -fpp -g -O3 -DNGP -mt_mpi cic_power_dmnu.f90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_power_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Using MKL on the SciNet GPC compile with:
!!   mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -mt_mpi cic_power_dmnu.f90 
!!        -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_power_dmnu -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64
!!
!! * Using FFTW on the SciNet BGQ compile with:
!!   mpif90 -q64 -O3 -qhot -qarch=qp -qtune=qp -WF,-DNGP cic_power_dmnu.F90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_power_dmnu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Optional flags:
!!   -DNGP: Uses NGP interpolation for binning of the power spectrum. 
!!   -DSLAB: Alternatively run with FFTW slab decomposition instead of P3DFFT pencil decomposition.
!!   -DLOGBIN: Use this to compute power spectrum with logarithmically spaced bins. Can change number of bins below.
!!   -DCOARSE_HACK: Enables a hack to coarsen the grid for which the particles are interpolated to.
!!                  Must be used in conjunction to a change made in the parameters file (see below).
!!   -Dwrite_den: Writes gridded density fields to binary files. 
!!   -Dwrite_poisson: Writes the gridded Poisson field to binary files.
!!   -DGROUPS: Instead of using PoissinNoise subroutine, remove noise by splitting each population into two groups.
!!   -DKAISER: Adjusts for redshift space distortions.
!!   -DPLPLOT: Plots power spectra at end.
!!   -DDEBUG: Output useful debugging information.

program cic_power_dmnu

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

  logical, parameter :: correct_kernel=.false.
  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints'

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

  integer(4), parameter :: max_np_dm = max_np / ratio_nudm_dim**3
  integer(4), parameter :: np_buffer=int(2./3.*max_np)

  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local, np_local_dm
#ifdef GROUPS
  real(4) :: np_groups(0:1), np_groups_dm(0:1)
#endif
  integer(8) :: plan, iplan
  logical :: firstfftw

  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(6,max_np_dm) :: xvp_dm
  real, dimension(3,np_buffer) :: xp_buf
  real, dimension(3*np_buffer) :: send_buf, recv_buf
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: den 
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: den_buf 

  !! Power spectrum arrays
  real, dimension(3,nc) :: pkdm
  real, dimension(3,nc) :: poisson_dm, poisson_nu
#ifdef PLPLOT
  real*8, dimension(3,nc) :: pkplot
#endif

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
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work, slab2
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab, slab2
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1) :: recv_cube
#endif

#ifdef GROUPS
  !! Array storing what group each particle belongs to (for autocorrelation break each species 
  !! into two groups and then cross-correlate to remove shot noise)
  integer(1) :: GID(max_np)
#endif

  !! Equivalence arrays to save memory
#ifdef SLAB
  equivalence (den,slab_work,recv_cube,xp_buf)
#else
  equivalence (den,recv_cube,xp_buf) 
#endif
#ifndef GROUPS
  equivalence (xvp,slab,cube)
#else
  equivalence (slab,cube)
#endif

  !! Common block
#ifdef PLPLOT
 common xvp,xvp_dm,send_buf,den_buf,den,recv_buf,pkdm,poisson_dm,poisson_nu,pkplot,slab2
#else
 common xvp,xvp_dm,send_buf,den_buf,den,recv_buf,pkdm,poisson_dm,poisson_nu,slab2
#endif
#ifdef GROUPS
  common /gvar/ slab, GID
#endif

!!---start main--------------------------------------------------------------!!

    call mpi_initialize
#ifdef GROUPS
  call initialize_random_number
#endif

    if (rank == 0) call writeparams
    firstfftw=.true.  ! initialize fftw so that it generates the plans
  
    call read_checkpoint_list
  
    do cur_checkpoint = 1, num_checkpoints
        
        call initvar

#ifndef GROUPS 
        !
        ! Compute Poisson noise first so don't overwrite slab's later
        !

        call PoissonNoise(0) 
        call powerspectrum(slab,slab,poisson_nu)
        call PoissonNoise(1) 
        call powerspectrum(slab,slab,poisson_dm)
#endif

        !
        ! Compute neutrino power spectrum
        !

        call read_particles(0)
        call pass_particles(0)
#ifdef GROUPS
        !! Randomly assign neutrinos to one of two groups 
        call assign_groups(0)
        !! Fourier transform density field for group 0 and store in slab2
        call darkmatter(0, 0, 0)
        slab2(:,:,:) = slab(:,:,:)
        !! Fourier transform density field for group 1 and compute correlation with group 0
        call darkmatter(0, 1, 0)
        call powerspectrum(slab,slab2,pkdm)
#else
        call darkmatter(0)
        call powerspectrum(slab,slab,pkdm)
#endif
        if (rank == 0) call writepowerspectra(0)
    
#ifndef GROUPS
        !! Store dark matter Fourier field in slab2 for cross spectra later
        slab2(:,:,:) = slab(:,:,:) 
#endif    

        !
        ! Compute dark matter power spectrum
        !
 
        call read_particles(1)
        call pass_particles(1)
#ifdef GROUPS
        !! Randomly assign dark matter to one of two groups 
        call assign_groups(1)
        !! Fourier transform density field for group 0 and store in slab2
        call darkmatter(1, 0, 0)
        slab2(:,:,:) = slab(:,:,:)
        !! Fourier transform density field for group 1 and compute correlation with group 0
        call darkmatter(1, 1, 0)
        call powerspectrum(slab,slab2,pkdm)
#else
        call darkmatter(1)
        call powerspectrum(slab,slab,pkdm)
#endif
        if (rank == 0) call writepowerspectra(1)

        !
        ! Compute cross power spectra
        !

#ifdef GROUPS
        !! Put all particles into group 0
        call clear_groups
        !! Fourier transform neutrino density field and store in slab2
        call darkmatter(0, 0, 1)
        slab2(:,:,:) = slab(:,:,:)
        !! Fourier transform dark matter density field and cross with neutrinos
        call darkmatter(1, 0, 1) 
#endif
        call powerspectrum(slab,slab2,pkdm)
        if (rank == 0) call writepowerspectra(2)

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
      write(*,*) 'cic_pow compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'cic_pow nodes=',nodes 
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
      write(*,*) 'cic_pow running on',nodes,'nodes'
      write(*,*) 'using cubic distribution:',nodes_dim,'nodes per dimension'
      write(*,*) nc,'cells in mesh'
    endif

!! calculate coordinates within slab for cube processes

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

#ifdef DEBUG
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

  subroutine read_particles(command)
    implicit none
    
    real z_write,np_total
    integer i,j,fstat, blocksize, nplow, nphigh, num_writes
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name
    integer :: command

    !! these are unnecessary headers from the checkpoint
    real(4) :: a,t,tau,dt_f_acc,dt_c_acc,dt_pp_acc,mass_p
    integer(4) :: nts,sim_checkpoint,sim_projection,sim_halofind

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

    if (command == 0) then
        if(z_write .eq. z_i) then
           check_name=ic_path//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        else
           check_name=output_path//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        endif
    else
        if(z_write .eq. z_i) then
           check_name=ic_path//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        else
           check_name=output_path//z_string(1:len_trim(z_string))//'xv'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        endif
    endif

!! open checkpoint    

    open(unit=21,file=check_name,status="old",iostat=fstat,access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint'
      write(*,*) 'rank',rank,'file:',check_name
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! read in checkpoint header data

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

#ifdef COARSE_HACK
    if (command == 0) then
        do j=1, np_local
            xvp(1:3,j) = xvp(1:3,j)/coarsen_factor
        enddo
    else
        do j=1, np_local_dm
            xvp_dm(1:3,j) = xvp_dm(1:3,j)/coarsen_factor
        enddo
    endif
#endif

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

  end subroutine read_particles

!!---------------------------------------------------------------------------!!

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

!-------------------------------------------------------------------!

#ifdef SLAB
subroutine unpack_slab
!! unpack slab data into cubic decomposition following fftw transform
    implicit none

    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag
    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(mpi_status_size,2*nodes_dim**2) :: wait_status

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

!-------------------------------------------------------------------!

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

!-------------------------------------------------------------------!

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

!!------------------------------------------------------------------!!

  subroutine writepowerspectra(command)
    implicit none
    integer      :: k
#ifdef PLPLOT
    integer :: kp
#endif
    real         :: kr
    character*180 :: fn
    character*5  :: prefix
    character*7  :: z_write
    integer :: command
    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is dm d2(k)
    !! 3rd is standard deviation

    write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
    z_write=adjustl(z_write)
    
#ifdef NGP 
    prefix='ngpps'
#else
    prefix='cicps'
#endif

#ifdef KAISER
    if (command == 1) then !! Dark matter power spectra
        fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD.dat' 
    else if (command == 0) then !! Neutrino power spectra
        fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD_nu.dat'
    else !! Dark matter-neutrino cross spectra
        fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD_dmnu.dat'
    endif
#else
    if (command == 1) then !! Dark matter power spectra
        fn=output_path//z_write(1:len_trim(z_write))//prefix//'.dat'
    else if (command == 0) then !! Neutrino power spectra
        fn=output_path//z_write(1:len_trim(z_write))//prefix//'_nu.dat'
    else !! Dark matter-neutrino cross spectra
        fn=output_path//z_write(1:len_trim(z_write))//prefix//'_dmnu.dat'
    endif
#endif

    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
#ifdef LOGBIN
    do k=2,numbins+1
#else
    do k=2,hc+1
#endif
       kr=2*pi*(k-1)/box
#ifdef NGP
#ifdef GROUPS
         write(11,*) pkdm(3,k-1),pkdm(1:2,k-1)
#else
       if (command == 1) then
         write(11,*) pkdm(3,k-1),pkdm(1:2,k-1),poisson_dm(1:2,k-1)
       else if (command == 0) then
         write(11,*) pkdm(3,k-1),pkdm(1:2,k-1),poisson_nu(1:2,k-1)
       else !! No poission noise with cross spectra
         write(11,*) pkdm(3,k-1),pkdm(1:2,k-1)
       endif
#endif
#else
#ifdef GROUPS
         write(11,*) pkdm(3,k),pkdm(1:2,k)
#else
       if (command == 1) then
         write(11,*) pkdm(3,k),pkdm(1:2,k),poisson_dm(1:2,k)
       else if (command == 0) then
         write(11,*) pkdm(3,k),pkdm(1:2,k),poisson_nu(1:2,k)
       else !! No poission noise with cross spectra
         write(11,*) pkdm(3,k),pkdm(1:2,k)
       endif
#endif
#endif
#ifdef PLPLOT
       kp=k-1
       pkplot(1,kp)=real(kr,kind=8)
       pkplot(2:3,kp)=real(pkdm(:,k),kind=8)
#endif
    enddo
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

#ifdef GROUPS
subroutine darkmatter(command, glook, canwrite)
#else
subroutine darkmatter(command)
#endif

    implicit none
    integer :: i,j,k, fstat
    integer :: i1,j1,k1
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,z_write
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name
    integer :: command
#ifdef GROUPS
    integer(4) :: glook, canwrite
#endif

    real time1,time2
    call cpu_time(time1)

    !! Initialized density field to be zero
    !! could do OMP loop here
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo

    !! Assign masses to grid to compute dm power spectrum
#ifdef GROUPS
    call cicmass(command, glook)
#else
    call cicmass(command)
#endif

    !! have to accumulate buffer density 
    call mesh_buffer
    cube=den(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)

#ifdef write_den
#ifdef GROUPS
    if (canwrite == 1) then
#endif
    !! generate checkpoint names on each node
        if (rank==0) then
           z_write = z_checkpoint(cur_checkpoint)
           print *,'Wrinting density to file for z = ',z_write
        endif

        call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

        write(z_string,'(f7.3)') z_write
        z_string=adjustl(z_string)
        write(rank_string,'(i4)') rank
        rank_string=adjustl(rank_string)

#ifdef KAISER
        if (command == 1) then !! Dark matter
        check_name=output_path//z_string(1:len_trim(z_string))//'den'// &
                   rank_string(1:len_trim(rank_string))//'-rsd.dat'
        else !! Neutrinos
        check_name=output_path//z_string(1:len_trim(z_string))//'den'// &
                   rank_string(1:len_trim(rank_string))//'-rsd_nu.dat'
        endif
#else 
        if (command == 1) then !! Dark matter
        check_name=output_path//z_string(1:len_trim(z_string))//'den'// &
                   rank_string(1:len_trim(rank_string))//'.dat'
        else !! Neutrinos
        check_name=output_path//z_string(1:len_trim(z_string))//'den'// &
                   rank_string(1:len_trim(rank_string))//'_nu.dat'
        endif
#endif

    !! open and write density file   
        open(unit=21,file=check_name,status="replace",iostat=fstat,access="stream")

        if (fstat /= 0) then
          write(*,*) 'error opening density file'
          write(*,*) 'rank',rank,'file:',check_name
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        write(21) cube
#ifdef GROUPS
    endif
#endif
#endif

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

    if (rank==0) then
      dsum=dsumt/real(nc)**3
      dvar=sqrt(dvart/real(nc)**3)
      write(*,*)
      if (command == 1) then
          write(*,*) 'DM min    ',dmint
          write(*,*) 'DM max    ',dmaxt
          write(*,*) 'Delta sum ',real(dsum,8)
          write(*,*) 'Delta var ',real(dvar,8)
      else
          write(*,*) 'NU min    ',dmint
          write(*,*) 'NU max    ',dmaxt
          write(*,*) 'Delta sum ',real(dsum,8)
          write(*,*) 'Delta var ',real(dvar,8)
      endif
      write(*,*)
    endif
 
    !! Forward FFT dm delta field
    call cp_fftw(1)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine darkmatter

!!------------------------------------------------------------------!!

#ifndef GROUPS
  subroutine PoissonNoise(command)
    implicit none
    integer :: i,j,k, fstat
    integer :: i1,j1,k1
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,z_write
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name
    integer :: command
    real time1,time2

    call cpu_time(time1)

    !! Initialized density field to be zero
    !! could do OMP loop here
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo

    ! Randomize positions across all nodes, hence np_local should be equal
    if (command == 0) then 
        np_local = np_node_dim**3 
    else !! Dark matter are further reduced by factor of ratio_nudm_dim 
        np_local_dm = (np_node_dim/ratio_nudm_dim)**3
    endif

    ! Assign particles to random positions between 0 and nc_node_dim
    if (command == 0) then
        call random_number(xvp(1:3,:np_local))
        xvp(1:3,:np_local) = xvp(1:3,:np_local)*nc_node_dim
    else
        call random_number(xvp_dm(1:3,:np_local_dm))
        xvp_dm(1:3,:np_local_dm) = xvp_dm(1:3,:np_local_dm)*nc_node_dim
    endif

    !! Assign masses to grid to compute dm power spectrum
    call cicmass(command)

    !! have to accumulate buffer density 
    call mesh_buffer
    cube=den(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)

#ifdef write_poisson
!! generate checkpoint names on each node
    if (rank==0) then
       z_write = z_checkpoint(cur_checkpoint)
       print *,'Wrinting density to file for z = ',z_write
    endif

    call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

    write(z_string,'(f7.3)') z_write
    z_string=adjustl(z_string)
    write(rank_string,'(i4)') rank
    rank_string=adjustl(rank_string)

#ifdef KAISER
    if (command == 1) then !! Dark matter
    check_name=output_path//z_string(1:len_trim(z_string))//'den-poisson'// &
               rank_string(1:len_trim(rank_string))//'-rsd.dat'
    else !! Neutrinos
    check_name=output_path//z_string(1:len_trim(z_string))//'den-poisson'// &
               rank_string(1:len_trim(rank_string))//'-rsd_nu.dat'
    endif
#else 
    if (command == 1) then !! Dark matter
    check_name=output_path//z_string(1:len_trim(z_string))//'den-poisson'// &
               rank_string(1:len_trim(rank_string))//'.dat'
    else !! Neutrinos
    check_name=output_path//z_string(1:len_trim(z_string))//'den-poisson'// &
               rank_string(1:len_trim(rank_string))//'_nu.dat'
    endif
#endif

!! open and write density file   
    open(unit=21,file=check_name,status="replace",iostat=fstat,access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening density file'
      write(*,*) 'rank',rank,'file:',check_name
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    write(21) cube
#endif

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

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return

  end subroutine PoissonNoise
#endif

!!------------------------------------------------------------------!!

subroutine pass_particles(command)
    !
    ! Pass particles inside buffer space to their appropriate nodes.
    !

    implicit none

    integer i,pp,np_buf,np_exit,np_final,npo,npi
    real x(3),lb,ub
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
            x = xvp(1:3, pp)
        else
            x = xvp_dm(1:3, pp)
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
                xp_buf(:, np_buf) = xvp(1:3, pp)
                xvp(:, pp)        = xvp(:, np_local)
                np_local          = np_local - 1
            else
                xp_buf(:, np_buf) = xvp_dm(1:3, pp)
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
      if (rank == i .and. command == 1) print *, rank, 'x+ np_local_dm=',
np_local_dm
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
            xvp(1:3, np_local) = x
        else
            np_local_dm = np_local_dm + 1
            xvp_dm(1:3, np_local_dm) = x
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
            xvp(1:3, np_local) = x
        else
            np_local_dm = np_local_dm + 1
            xvp_dm(1:3, np_local_dm) = x
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
            xvp(1:3, np_local) = x
        else
            np_local_dm = np_local_dm + 1
            xvp_dm(1:3, np_local_dm) = x
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
            xvp(1:3, np_local) = x
        else
            np_local_dm = np_local_dm+1
            xvp_dm(1:3, np_local_dm) = x
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
            xvp(1:3,np_local)=x
        else
            np_local_dm=np_local_dm+1
            xvp_dm(1:3,np_local_dm)=x
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
            xvp(1:3,np_local)=x
        else
            np_local_dm=np_local_dm+1
            xvp_dm(1:3,np_local_dm)=x
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

!------------------------------------------------------------!

#ifdef GROUPS
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

!------------------------------------------------------------!

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

        GID(k) = g

    enddo

    call mpi_allreduce(int(np_tot,kind=8), npa_tot, 1, mpi_integer8, mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(int(np0,kind=8), np1_tot, 1, mpi_integer8, mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(int(np1,kind=8), np2_tot, 1, mpi_integer8, mpi_sum, mpi_comm_world, ierr)

    if (command == 0) then
        np_groups(0) = real(np1_tot)
        np_groups(1) = real(np2_tot)
    else
        np_groups_dm(0) = real(np1_tot)
        np_groups_dm(1) = real(np2_tot)
    endif

    if (rank == 0) then
        write(*,*) "Groups assigned: ", np1_tot, np2_tot, npa_tot
    endif

    return

end subroutine assign_groups

!------------------------------------------------------------!

subroutine clear_groups
    !
    ! Assign all particles to group 0
    !

    implicit none

    GID(:) = 0

    np_groups(0) = real(np)**3 
    np_groups_dm(0) = (real(np)/ratio_nudm_dim)**3 
    np_groups(1) = 0.
    np_groups_dm(1) = 0.

    return

end subroutine clear_groups
#endif

!------------------------------------------------------------!

#ifdef GROUPS
  subroutine cicmass(command, glook)
#else
  subroutine cicmass(command)
#endif
    implicit none
    real :: mp
    integer :: i,i1,i2,j1,j2,k1,k2,np_total
    real    :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2,vf,v(3)
    integer :: command
#ifdef GROUPS
    integer(4) :: glook
#endif

    if (command == 1) then !! Dark matter are further reduced by factor of ratio_nudm_dim 
#ifdef GROUPS
        mp = ncr**3 / np_groups_dm(glook)
#else
        mp = (ncr/(np/ratio_nudm_dim))**3
#endif
        np_total = np_local_dm
    else
#ifdef GROUPS
        mp = ncr**3 / np_groups(glook)
#else
        mp = (ncr/np)**3
#endif
        np_total = np_local
    endif 

    do i=1,np_total 
#ifdef GROUPS
       if (GID(i) /= glook) cycle
#endif

       if (command == 0) then
         x=xvp(1,i)-0.5
         y=xvp(2,i)-0.5
         z=xvp(3,i)-0.5
       else
         x=xvp_dm(1,i)-0.5
         y=xvp_dm(2,i)-0.5
         z=xvp_dm(3,i)-0.5
       endif

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
         print *,'particle out of bounds',i1,i2,j1,j2,k1,k2,nc_node_dim
       endif 

       dz1=mp*dz1
       dz2=mp*dz2
       den(i1,j1,k1)=den(i1,j1,k1)+dx1*dy1*dz1
       den(i2,j1,k1)=den(i2,j1,k1)+dx2*dy1*dz1
       den(i1,j2,k1)=den(i1,j2,k1)+dx1*dy2*dz1
       den(i2,j2,k1)=den(i2,j2,k1)+dx2*dy2*dz1
       den(i1,j1,k2)=den(i1,j1,k2)+dx1*dy1*dz2
       den(i2,j1,k2)=den(i2,j1,k2)+dx2*dy1*dz2
       den(i1,j2,k2)=den(i1,j2,k2)+dx1*dy2*dz2
       den(i2,j2,k2)=den(i2,j2,k2)+dx2*dy2*dz2
    enddo

    return
  end subroutine cicmass

!!--------------------------------------------------------------!!

  subroutine powerspectrum(delta,delta2,pk)
    implicit none
    real, dimension(3,nc)       :: pk
#ifdef SLAB
    real, dimension(nc+2,nc,nc_slab) :: delta
    real, dimension(nc+2,nc,nc_slab) :: delta2
#else
    real, dimension(nc, nc_node_dim, nc_pen+2) :: delta
    real, dimension(nc, nc_node_dim, nc_pen+2) :: delta2
#endif

    integer :: i,j,k,kg,ig,mg,jg
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow, x,y,z,sync_x, sync_y,sync_z,kernel
#ifdef SLAB
    real, dimension(3,nc,nc_slab) :: pkt
#else
    real, dimension(3,nc,nc_pen+2) :: pkt
#endif
    real, dimension(3,nc) :: pktsum

    real(8), dimension(nc) :: kcen, kcount
    real(8), dimension(nc) :: kcensum, kcountsum
    real    :: kavg
    integer :: ind, dx, dxy

#ifdef LOGBIN
    real :: k1r
#endif

    real time1,time2
    call cpu_time(time1)

    pkt=0.0
    pktsum=0.0

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
          pk(1:2,k)=4.*pi*(kavg)**3*pk(1:2,k)
#else
          pk(1:2,k)=4.*pi*(kavg-1.)**3*pk(1:2,k)
#endif

       endif
      enddo
    endif

    call mpi_bcast(pk,3*nc,mpi_real,0,mpi_comm_world,ierr)

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
#ifdef SLAB
    do k=1,nc_slab
       slab_work(:,:,k)=0
    enddo
#endif
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo
    do k=1,nc_node_dim
       cube(:,:,k)=0
    enddo
#ifdef SLAB
    do k=1,nc_slab
#else
    do k=1,nc_pen+2
#endif
       slab(:,:,k)=0
       slab2(:,:,k)=0
    enddo
    do k=1,np_buffer
       xp_buf(:,k)=0
    enddo
    do k=1,3*np_buffer
       recv_buf(k)=0
    enddo
    do k=1,nc
       pkdm(:,k)=0
       poisson_dm(:,k)=0
       poisson_nu(:,k)=0
    enddo    
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

end program cic_power_dmnu
