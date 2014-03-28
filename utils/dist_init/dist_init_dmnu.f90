!! init.f90 written by Hy Trac
!! dist_init.f90 Parallelized: Hugh Merz Jun 4, 2005
!! Updated Jun 2 2006 -- minimal memory usage, DM only !!
!!-----------------------------------------------------!!
!! Compile with: mpif77 -fpp -g -w -O3 -axN dist_init.f90 -o dist_init  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
!! On Ranger compile with:
!! mpif90 -Wl,-rpath, -L$TACC_MKL_LIB -I$TACC_MKL_INC -fpp -g -O2 -xW -shared-intel -DBINARY dist_init_dm.f90 -o dist_init -L$TACC_MKL_LIB -lmkl_em64t -openmp


program dist_init 

  implicit none

  include 'mpif.h'
  include '../../parameters'

  integer,parameter  :: nt=1
#ifdef NEUTRINOS
  logical, parameter :: generate_seeds=.false.
#else
  logical, parameter :: generate_seeds=.true.
#endif
  logical, parameter :: correct_kernel=.true.

  real, parameter :: ns = 0.96 
  real, parameter :: s8 = 0.779
  real, parameter :: omegal=omega_l 
  real, parameter :: omegam=1.0-omegal 

  real, parameter :: redshift=z_i 
  real, parameter :: scalefactor=1/(1+redshift)

#ifdef WDM
  real, parameter :: m_wdm = 1. ! keV
  real, parameter :: mu = 1.12
  real, parameter :: alpha = 0.049 * m_wdm**-1.11 * (omegam/0.25)**0.11 * (H0/100./0.7)**1.22
#endif

#ifdef BROKEN_POWER
  real, parameter :: dns = -0.05
  real, parameter :: k0  = 0.06
#endif

#ifdef NEUTRINOS
  logical, parameter :: nu_random = .false.
  logical, parameter :: nu_relfd = .true.
  real(4), parameter :: Vphys2sim = 1.0/(300. * sqrt(omega_m) * box * (1. + z_i) / 2. / nc)!(180.8892437/mass_neutrino)/(box*300.0*(omega_m)**0.5/2.0/nc)
  integer, parameter :: nv=10000
  character(*), parameter :: cdfTable = 'CDFTable.txt'
  real(4), dimension(2,nv) :: cdf    !Col1 is v, col2 is cdf
#endif

  !! NEUTRINO sims may use this:
  integer, parameter      :: nk=708
  character(*), parameter :: fntf = 'nu_transfer_out_z10.dat'

  !! Transfer function file
!  integer, parameter      :: nk=922
!  character(*), parameter :: fntf='camb_jddefault_transfer_z0.dat'

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
#ifdef NEUTRINOS
  integer, parameter :: np=hc 
#else
  integer, parameter :: np=hc/ratio_nudm_dim
#endif
  real, parameter    :: npr=np

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: np_node_dim = np/nodes_dim
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local,wc_counter, count_i,count_f,count_r
  logical :: firstfftw
  integer(8) :: plan, iplan

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.141592654

  !! Power spectrum arrays
!  real, dimension(5,nk) :: tf   !cmbfast
  real, dimension(7,nk) :: tf    !CAMB
  real, dimension(2,nc) :: pkm,pkn

  !! For pencils decomposition:
#ifdef SLAB
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
  integer(4), parameter :: nc_slab = nc / nodes
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
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work
  real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab, slab_work
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1)      :: recv_cube
#endif
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: phi
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: phi_buf

  !! Timing variables
  real(8) :: sec1, sec2

  !! Equivalence arrays to save memory
  equivalence (phi,slab_work,recv_cube) 
  equivalence (slab,cube)

  !! Common block
  common /rvar/ tf, pkm, pkn, phi_buf
  common / equiv1 / phi
  common / equiv2 / slab

  call mpi_initialize

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STARTING INITIAL CONDITIONS: ", sec1

  !$ call omp_set_num_threads(nt)
  if (rank == 0) call writeparams

  call initvar
  call transferfnc
  call noisemap
  call deltafield
  call potentialfield
  call dm

  if (rank == 0) call writepowerspectra

  call di_fftw(0)

  sec2 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STOPPING INITIAL CONDITIONS: ", sec2
  if (rank == 0) write(*,*) "ELAPSED TIME: ", sec2-sec1

  call mpi_finalize(ierr)

contains

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
#ifdef SLAB
    if (mod(nc,nodes) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into slabs'
      write(*,*) 'nc=',nc,'nodes=',nodes,'mod(nc,nodes) != 0'
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#else
    if (mod(nc,nodes_dim**2) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into pencils'
      write(*,*) 'nc=',nc, 'nodes_dim**2=',nodes_dim**2,'mod(nc,nodes_dim**2)=', mod(nc,nodes_dim**2)
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#endif
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

#ifdef DEBUG_LOW
  do i=0,nodes-1
    if (i==rank) write(*,'(8i5)') rank,cart_rank,cart_neighbor
    call mpi_barrier(mpi_comm_world,ierr)
  enddo
#endif

  end subroutine mpi_initialize

!-------------------------------------------------------------------!

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

    num_elements = nc_node_dim * nc_node_dim * nc_pen

    !
    ! Send the data from cube to recv_cube
    !

    do j = 0, nodes_dim - 1
        pen_slice = j
        tag  = rank**2
        rtag = pen_neighbor_fm(j)**2
        call mpi_isend(cube(1,1, pen_slice*nc_pen + 1), num_elements, &
                       mpi_real, pen_neighbor_to(j), tag, mpi_comm_world, &
                       requests(pen_slice+1),ierr)
        call mpi_irecv(recv_cube(1,1,1,pen_slice), &
                       num_elements, mpi_real, pen_neighbor_fm(j),rtag, &
                       mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                       ierr)
    enddo

    call mpi_waitall(2*nodes_dim, requests, wait_status, ierr)

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

    num_elements = nc_node_dim * nc_node_dim * nc_pen

    !
    ! Put this data back into cube
    !

    do j = 0, nodes_dim - 1
        pen_slice = j
        tag  = rank**2
        rtag = pen_neighbor_to(j)**2
        call mpi_isend(recv_cube(1,1,1,pen_slice), num_elements, &
                       mpi_real, pen_neighbor_fm(j), tag, mpi_comm_world, &
                       requests(pen_slice+1),ierr)
        call mpi_irecv(cube(1,1,pen_slice*nc_pen + 1), &
                       num_elements, mpi_real, pen_neighbor_to(j),rtag, &
                       mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                       ierr)
    enddo

    call mpi_waitall(2*nodes_dim,requests, wait_status, ierr)

end subroutine unpack_pencils
#endif

!-------------------------------------------------------------------!

subroutine di_fftw(command)
    !
    ! Calculate fftw transform using P3DFFT or FFTW slab decomposition (ifdef SLAB)
    ! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
    !

#ifndef SLAB
    use p3dfft
#endif
    implicit none
#ifdef SLAB
    include 'fftw_f77.i'
    integer(4), parameter :: order=FFTW_NORMAL_ORDER
#endif

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
#ifdef SLAB
      call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nc, &
            nc,nc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE) !_MEASURE)
      call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nc, &
            nc,nc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE) !_MEASURE)
#else
        pen_dims = (/dim_y,dim_z/)
        call p3dfft_setup(pen_dims, nc, nc, nc, .true.)
        call p3dfft_get_dims(istart, iend, isize, 1, mypadd)
        call p3dfft_get_dims(fstart, fend, fsize, 2)
#endif
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

#ifdef SLAB
      if (command > 0) call pack_slab
#else
      if (command > 0) call pack_pencils 
#endif

#ifdef DEBUG_LOW
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

#ifdef DEBUG_LOW
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

end subroutine di_fftw

!-------------------------------------------------------------------!

  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nodes   ', nodes
    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*)
    write(*,*) 'n_s      ',ns
    write(*,*) 'sigma_8  ',s8
    write(*,*) 'omega_m  ',omegam
    write(*,*) 'omega_l  ',omegal
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*) 'redshift ',redshift
    write(*,*)

#ifdef MY_TRANSFER
    write(*,*) 'power index',power_index
#endif

#ifdef WDM
    write(*,*) "m_wdm = ", m_wdm
    write(*,*) "mu    = ", mu
    write(*,*) "alpha = ", alpha
#endif

#ifdef BROKEN_POWER
    write(*,*) "dns = ", dns
    write(*,*) "k0  = ", k0
#endif

#ifdef NEUTRINOS
    write(*,*) 'ratio_nudm_dim ', ratio_nudm_dim 
    write(*,*) 'Vphys2sim ',Vphys2sim
    write(*,*) 'nu_random ', nu_random
    write(*,*) 'nu_relfd ', nu_relfd    
#endif

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams

!-------------------------------------------------------------------!

  subroutine writepowerspectra
    implicit none
    integer      :: k
    real         :: kr
    character*180 :: fn

    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is matter p(k)
    !! 3rd is standard deviation
    !! 6th is noise p(k)
    !! 7th is noise standard deviation
#ifdef NEUTRINOS
    fn=output_path//'pk_nu.init'
#else
    fn=output_path//'pk.init'
#endif
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       write(11,*) kr,pkm(:,k),pkn(:,k)
    enddo
    close(11)

    !! Output cmbfast power spectrum
#ifdef NEUTRINOS
    fn=output_path//'pk0_nu.init'
#else
    fn=output_path//'pk0.init'
#endif
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,2*nc+1
       kr=2*pi*(k-1)/(2*box)
       write(11,*) kr,power(kr,1,2),power(kr,1,3)
    enddo
    close(11)
    
    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write power spectra'
    return
  end subroutine writepowerspectra

!!------------------------------------------------------------------!!

  subroutine transferfnc
    implicit none
    integer :: i,k
    real    :: kr,kmax,dummy
    real*8  :: v8

    real time1,time2
    call cpu_time(time1)

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (i==rank) print *,'rank:',rank,'starting transferfnc'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo 
#endif 

    if (rank == 0) then
      !! Get transfer function from CMBFAST or CAMB
      write(*,*) 'Reading ',fntf
      open(11,file=fntf)
!      read(11,*) tf
      do k=1,nk
         read(11,*) tf(1,k),tf(2,k),tf(3,k),tf(4,k),tf(5,k),tf(6,k),tf(7,k)
      end do
      close(11)

      !! Compute \Delta^2
      do k=1,nk
         kr     =tf(1,k)
         tf(2,k)=kr**(3+ns)*tf(2,k)**2/(2*pi**2)
         tf(3,k)=kr**(3+ns)*tf(3,k)**2/(2*pi**2)
#ifdef NEUTRINOS
         tf(6,k)=kr**(3+ns)*tf(6,k)**2/(2*pi**2)
#endif

#ifdef MY_TRANSFER
         tf(2:3,k) = kr**(3-power_index)/(2*pi**2) 
         write(*,*)'Over-written by power law'
#endif

#ifdef WDM
        tf(2:3,k) = tf(2:3,k)*wdm_transfer(kr)**2
#endif
#ifdef BROKEN_POWER
       tf(2:3,k) = tf(2:3,k)*broken_power_transfer(kr)**2
#endif

      enddo

      !! Compute dk
      tf(4,1)=tf(1,2)/2
      do k=2,nk-1
         tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
      enddo
      tf(4,nk)=tf(1,nk)-tf(1,nk-1)

      !! Compute variance in 8 h^-1 Mpc spheres
      v8=0
      kmax=2*pi*sqrt(3.)*hc/box
      do k=1,nk
         if (tf(1,k) .gt. kmax) exit
         v8=v8+tf(2,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
      enddo

      !! Normalize to \sigma_8
      !! Include growth factor
#ifdef NEUTRINOS
      !! Set dm transfer function = neutrino one
      !! This step makes sure that the neutrinos have the correct
      !! sigma8 normalization
      tf(2,:) = tf(6,:)*(s8**2/v8)*Dgrow(scalefactor)**2
#else
      tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*Dgrow(scalefactor)**2
#endif
    endif

    call mpi_barrier(mpi_comm_world,ierr)

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (i==rank) print *,rank,'finished mpi_barrier',ierr
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_bcast(tf,7*nk,mpi_real,0,mpi_comm_world,ierr) 
   
#ifdef DEBUG_LOW
    print *,rank,'finished mpi_bcast',ierr
#endif

    !! tf(1,i) stores k
    !! tf(2,i) stores \Delta^2_m
    !! tf(3,i) stores \Delta^2_b
    !! tf(4,i) stores dk

    call cpu_time(time2)
    time2=time2-time1
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called transfer fnc'
    return
  end subroutine transferfnc

!------------------------------------------------------------------!

  subroutine noisemap
    implicit none
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s
    integer      :: i,j,k,seedsize,clock
    integer, dimension(8) :: values

#ifdef POWER_CHECK
    integer :: ig,jg,kg
#endif
    real         :: x,x1,x2
    character*180 :: fn
    integer(4),allocatable,dimension(:) :: iseed
    integer(4),allocatable,dimension(:) :: iseed_all

    real time1,time2
    call cpu_time(time1)

    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)

#ifdef DEBUG_LOW
    print *,'rank:',rank,'starting noisemap'
#endif 

    call random_seed() 
#ifdef IBM
    call random_seed(generator=2)
#endif
    call random_seed(size=seedsize)

    allocate(iseed(seedsize))
    allocate(iseed_all(seedsize*nodes))
 
    call system_clock(count=clock)
    iseed = clock + 37 * (/ (i - 1, i = 1, 2) /)
 
    call date_and_time(values=values)
    if(rank==0) write(*,*) values(7:8), iseed

    call random_seed(put=values(7:8)) 
    call random_seed(put=iseed(1:seedsize))

    if (generate_seeds) then
#ifdef NODE_DEPENDANT
       if (rank == 0) then
         write(*,*) 'Generating seeds'
         do k=1,seedsize
            do j=0,rank
              call random_number(x)
            enddo
            iseed(k)=int(x*huge(0)+time1)+rank
         enddo
       endif
#else
       if (rank==0) then
         write(*,*) 'Generating seeds'
         do j=1,nodes
           do k=1,seedsize
             call random_number(x)
             iseed_all((j-1)*seedsize+k)=int(x*huge(0))
           enddo
         enddo
       endif
       call mpi_scatter(iseed_all,seedsize,mpi_integer,iseed,seedsize,mpi_integer,0,mpi_comm_world,ierr)
#endif
    else
       fn=output_path//'seed'//rank_s(1:len_trim(rank_s))//'.init'
       print *, 'rank',rank,'Reading ',fn(1:len_trim(fn))
       open(11,file=fn)
       do k=1,seedsize
          read(11,*) i,iseed(k)
       enddo
       close(11)
    endif

    !! Generate random reals between 0 and 1
    call random_seed(put=iseed(1:seedsize))
    call random_number(cube)

    fn=output_path//'seed'//rank_s(1:len_trim(rank_s))//'.init'
    print *, 'rank',rank,'Writing ',fn(1:len_trim(fn))
    open(11,file=fn)
    do k=1,seedsize
       write(11,*) k,iseed(k)
    enddo
    close(11)
   
    deallocate(iseed)
    deallocate(iseed_all)

    if (rank==0) then
       write(*,*) 'starting generation of Gaussian random numbers'    
    end if

    !! Generate random Gaussian numbers
    !$omp parallel do default(shared) private(i,j,k,x1,x2)
    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim,2
             x1=2*pi*cube(i,j,k)
             if (cube(i+1,j,k)<=0.0) cube(i+1,j,k)=0.0001
             x2=sqrt(-2*log(cube(i+1,j,k)))
             cube(i,j,k)=x2*cos(x1)
             cube(i+1,j,k)=x2*sin(x1)
          enddo
       enddo
    enddo
    !$omp end parallel do

    if (rank==0) then
       write(*,*) 'finished generation of Gaussian random numbers'
    end if

#ifdef POWER_CHECK
    do k=1,nc_node_dim
      kg=k+cart_coords(1)*nc_node_dim
      do j=1,nc_node_dim
        jg=j+cart_coords(2)*nc_node_dim
        do i=1,nc_node_dim
          ig=i+cart_coords(3)*nc_node_dim
          cube(i,j,k)=mod(ig,2)+mod(jg,4)+mod(kg,8)
        enddo
      enddo
    enddo
#endif

    if (rank==0) then
       write(*,*) 'starting forward FFT'
    end if

    !! Forward FFT white noise field
    call di_fftw(1)
    !! noise is now in slab array

    if (rank==0) then
       write(*,*) 'finished forward FFT'
    end if

    !! Generate noise spectrum
    call powerspectrum(pkn)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called noisemap'
    return
  end subroutine noisemap

!-----------------------------------------------------------------------!

  subroutine deltafield
    implicit none
    integer :: i,j,k,kg,mg,jg,ig,fstat
    real    :: kr,kx,ky,kz
    real    :: powb,powm
    real    :: d,dmin,dmax,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart
    character(len=4) :: rank_string
    character(len=100) :: check_name

    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    !! Interpolate \Delta^2 to grid
    !! Determine Fourier modes \delta(k)

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

#ifdef SLAB
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(kr,kx,ky,kz,powb,powm,kg)
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
             if (kr .eq. 0) then
                slab(i:i+1,j,k)=0.
             else
                powm=power(2*pi*kr/box,1,2)/(4*pi*kr**3)
                slab(i:i+1,j,k)=sqrt(powm*ncr**3)*slab(i:i+1,j,k)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
#else
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
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr .eq. 0) then
                    slab(i:i+1,j,k)=0.
                else
                    powm=power(2*pi*kr/box,1,2)/(4*pi*kr**3)
                    slab(i:i+1,j,k)=sqrt(powm*ncr**3)*slab(i:i+1,j,k)
                endif
            enddo
        enddo
    enddo
#endif

    !! Calculate matter delta power spectrum

    call powerspectrum(pkm)

    !! Calculate matter delta field statistics

    call di_fftw(-1)

    dmin=0
    dmax=0
    dsum=0
    dvar=0

    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) &
    !$omp& reduction(max:dmax) &
    !$omp& reduction(+:dsum,dvar)
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
    !$omp end parallel do

    call mpi_reduce(dsum,dsumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dvar,dvart,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dmin,dmint,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
    call mpi_reduce(dmax,dmaxt,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

    if (rank==0) then
      dsum=dsumt/real(nc)**3
      dvar=sqrt(dvart/real(ncr)**3)
      write(*,*)
      write(*,*) 'Delta min    ',dmint
      write(*,*) 'Delta max    ',dmaxt
      write(*,*) 'Delta sum ',real(dsum)
      write(*,*) 'Delta var ',real(dvar)
      write(*,*)
    endif

#ifdef write_den 
    if (rank == 0) then
        print *,'Writing density contrast to file'
    endif

    write(rank_string,'(i4)') rank
    rank_string = adjustl(rank_string)

#ifdef NEUTRINOS
    check_name = output_path//'initdeltafield'// &
        rank_string(1:len_trim(rank_string))//'_nu.bin'
#else    
    check_name = output_path//'initdeltafield'// &
        rank_string(1:len_trim(rank_string))//'.bin'
#endif

    open(unit=21, file=check_name, status='replace', iostat=fstat, access='stream')

    if (fstat /= 0) then
        write(*,*) 'error opening density file'
        write(*,*) 'rank', rank, 'file:', check_name
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim

            write(21) cube(:, j, k)

        enddo
    enddo

#endif

    call di_fftw(1)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called delta field'
    return
  end subroutine deltafield

!!------------------------------------------------------------------!!

  subroutine potentialfield
    implicit none

    integer :: i,j,k,ioerr
    integer :: im,ip,ig,jm,jp,jg,km,kp,kg,mg
    real    :: r,x,y,z
    real    :: kr,ksq,kx,ky,kz
    real    :: phi8,phi8tot
    character*180 :: fn
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s

    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    if (correct_kernel) then
      !! write delta to disk so we can build kernel
      write(*,*) 'Caching Delta on disk'
      write(rank_s,'(i6)') rank
      rank_s=adjustl(rank_s)
      fn=scratch_path//'delta'//rank_s(1:len_trim(rank_s))

      open(11,file=fn,status='replace',iostat=ioerr, access='stream')

      if (ioerr /= 0) then
        print *,'error opening Delta cache file:',fn
        stop
      endif
#ifdef SLAB
      do k=1,nc_slab
#else
      do k=1,nc_pen+2
#endif
        write(11) slab(:,:,k)
      enddo
      close(11)
    else
      slab_work=slab
    endif

    !! Construct uncorrected potential kernel in Fourier space

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

#ifdef SLAB
    !$omp parallel do default(shared) private(i,j,k,ip,kg) &
    !$omp& private(ksq,kx,ky,kz)
    do k=1,nc_slab
       kg=k+nc_slab*rank
       if (kg .lt. hc+2) then
          kz=kg-1
       else
          kz=kg-1-nc
       endif
       kz=2*sin(pi*kz/ncr)
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          ky=2*sin(pi*ky/ncr)
          do i=1,nc+2,2
             ip=i+1
             kx=(i-1)/2
             kx=2*sin(pi*kx/ncr)
             ksq=kx**2+ky**2+kz**2
             if (ksq .eq. 0) then
                slab(i:ip,j,k)=0
             else
                slab(i,j,k)=-4*pi/ksq
                slab(ip,j,k)=0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
#else
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
                kz=2*sin(pi*kz/ncr)
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                ky=2*sin(pi*ky/ncr)
                kx = (ig-1)/2
                kx=2*sin(pi*kx/ncr)
                ksq=kx**2+ky**2+kz**2
                ip=i+1
                if (ksq .eq. 0) then
                    slab(i:ip,j,k)=0
                else
                    slab(i,j,k)=-4*pi/ksq
                    slab(ip,j,k)=0
                endif
            enddo
        enddo
    enddo
#endif

    if (correct_kernel) then

       !! Inverse FFT potential kernel
       call di_fftw(-1) 

       phi8=0.0

       if (nc_node_dim < 9) then
         print *,'warning: mesh too small in potential kernel correction'
         call mpi_abort(mpi_comm_world,ierr,ierr)
       endif

       if (cart_coords(1) == 0 .and. cart_coords(2) == 0 .and. &
           cart_coords(3) == 0) phi8=cube(9,1,1)+cube(1,9,1)+cube(1,1,9)

       if (cart_coords(3) == nodes_dim-1 .and. cart_coords(2) == 0 .and. &
           cart_coords(1) == 0) phi8=phi8+cube(nc_node_dim-7,1,1)

       if (cart_coords(3) == 0 .and. cart_coords(2) == nodes_dim-1 .and. &
           cart_coords(1) == 0) phi8=phi8+cube(1,nc_node_dim-7,1)

       if (cart_coords(3) == 0 .and. cart_coords(2) == 0 .and. &
           cart_coords(1) == nodes_dim -1) phi8=phi8+cube(1,1,nc_node_dim-7)

       call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       if (rank == 0) phi8=phi8tot/6.0      
       call mpi_bcast(phi8,1,mpi_real,0,mpi_comm_world,ierr)

 
       !! Construct Ewald potential kernel in real space
       !$omp parallel do default(shared) private(i,j,k) &
       !$omp& private(r,x,y,z,ig,jg,kg)
       do k=1,nc_node_dim
          kg=k+nc_node_dim*cart_coords(1)
          if (kg .lt. hc+2) then
             z=kg-1
          else
             z=kg-1-nc
          endif
          do j=1,nc_node_dim
             jg=j+nc_node_dim*cart_coords(2)
             if (jg .lt. hc+2) then
                y=jg-1
             else
                y=jg-1-nc
             endif
             do i=1,nc_node_dim
                ig=i+nc_node_dim*cart_coords(3)
                if (ig .lt. hc+2) then
                   x=ig-1
                else
                   x=ig-1-nc
                endif
                r=sqrt(x**2+y**2+z**2)
                if (r .gt. 8) then
                   cube(i,j,k)=cube(i,j,k)-(phi8+1/8.)
                else
                   if (r .eq. 0) then
                      cube(i,j,k)=-2.5
                   else
                      cube(i,j,k)=-1/r
                   endif
                endif
             enddo
          enddo
       enddo
       !$omp end parallel do

#ifdef DEBUG_KC
!       print *,'rank',rank,'phi8=',phi8
!       if (rank == 0) then
!         do i=1,16
!           print *,rank,cube(i,i,i)
!         enddo
!       endif
       print *,'rank',rank,'sum kern',sum(cube)
#endif

       !! Forward FFT potential kernel
       call di_fftw(1)


#ifdef DEBUG_KC
      phi8=0.0
      phi8tot=0.0
#ifdef SLAB
      do k=1,nc_slab
        do j=1,nc
          do i=1,nc+2
#else
      do k=1,nc_pen+2
        do j=1,nc_node_dim
          do i=1,nc
#endif
            phi8=phi8+slab(i,j,k)
          enddo
        enddo
      enddo
      call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
      if (rank==0) print *,'total slab=',phi8tot
      phi8=0.0
      phi8tot=0.0
#ifdef SLAB
      do k=1,nc_slab
        do j=1,nc
          do i=1,nc+2
#else
      do k=1,nc_pen+2
        do j=1,nc_node_dim
          do i=1,nc
#endif
            phi8=phi8+delta(i,j,k)
          enddo
        enddo
      enddo
      call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
      if (rank==0) print *,'total delta=',phi8tot
#endif

      open(11,file=fn,status='old',iostat=ioerr,access='stream')

      if (ioerr /= 0) then
        print *,'error opening Delta cache file:',fn
        stop
      endif
#ifdef SLAB
      do k=1,nc_slab
#else
      do k=1,nc_pen+2
#endif
        read(11) slab_work(:,:,k)
      enddo
      close(11)

    endif

    !! Complex multiply density field with potential kernel
    !$omp parallel do default(shared) private(i,j,k)
#ifdef SLAB
    do k=1,nc_slab
       do j=1,nc
          do i=1,nc+2,2
#else
    do k=1,nc_pen+2
       do j=1,nc_node_dim
          do i=1,nc,2
#endif
             slab(i:i+1,j,k)=slab_work(i:i+1,j,k)*slab(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    !! Inverse FFT potential field
    call di_fftw(-1)

    !! put cube in phi
    phi=0.0
    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube
    call mesh_buffer

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called potential field'
    return
  end subroutine potentialfield

!!------------------------------------------------------------------!!
    !! Dark matter data
    !! xvp(1:3) = x,y,z in grid units from 0 to nc
    !! xvp(4:6) = vx,vy,vz in km/s
  subroutine dm
    implicit none
    integer :: i,j,k,ioerr 
    integer :: i1,j1,k1,lb,ub
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,vf,xvp(6)
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis,x
    real*8, dimension(3) :: xav
    character*180 :: fn
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s

#ifdef NEUTRINOS
    integer :: l
    real(4) :: rnum1,rnum2,rnum3,rnum4
#endif

    real time1,time2
    call cpu_time(time1)
 

    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)
#ifdef NEUTRINOS
    fn=scratch_path//'xv'//rank_s(1:len_trim(rank_s))//'_nu.ic'
#else
    fn=scratch_path//'xv'//rank_s(1:len_trim(rank_s))//'.ic'
#endif

    open(11,file=fn,status='replace',iostat=ioerr,access='stream')

    if (ioerr /= 0) then
      print *,'error opening:',fn
      stop
    endif

#ifdef NEUTRINOS
    !! Read cdfTable
    if (rank == 0) then
        write(*,*) 'Reading ',cdfTable
        open(12, file=cdfTable)
        do k = 1, nv
        read(12,*) cdf(1,k), cdf(2,k)
        enddo
        close(12)
    endif
    call mpi_bcast(cdf, 2*nv, mpi_real, 0, mpi_comm_world, ierr)
#endif

    np_local=np_node_dim**3
    write(11) np_local

    !! Displace particles
    !! Finite-difference potential to get displacement field
    vf=vfactor(scalefactor)
    do k=1,np_node_dim
       k1=(nc/np)*(k-1)+1
       do j=1,np_node_dim
          j1=(nc/np)*(j-1)+1
          do i=1,np_node_dim

             i1=(nc/np)*(i-1)+1
             dis(1)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)
             dis(2)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)
             dis(3)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)

             xvp(1)=dis(1)+(i1-0.5)
             xvp(2)=dis(2)+(j1-0.5)
             xvp(3)=dis(3)+(k1-0.5)

#ifdef NEUTRINOS
             if (nu_random) then !! Uses Gaussian velocity distribution

                !uniform random #
                call random_number(rnum1)
                call random_number(rnum2)

                !convert to gaussian 
                rnum3=2*pi*rnum1
                rnum4=sqrt(-2*log(rnum2))

                rnum1=rnum4*cos(rnum3)
                rnum2=rnum4*sin(rnum3)

                !Convert to real gaussian using dispersion
                !scale factor cancels out in velocity dispersion and conversion factor
                !sigma_nu = 181km/s / (mnu * a)
                rnum1 = rnum1*Vphys2sim*(180.8892437/mass_neutrino/scalefactor/3.0**0.5) !Divide by sqrt3  
                rnum2 = rnum2*Vphys2sim*(180.8892437/mass_neutrino/scalefactor/3.0**0.5) !for each individual component

                xvp(4)=dis(1)*vf + rnum1
                xvp(5)=dis(2)*vf + rnum2

                call random_number(rnum3)
                call random_number(rnum4)

                rnum3=2*pi*rnum3
                rnum4=sqrt(-2*log(rnum4))

                rnum3=rnum4*cos(rnum3)
                rnum3 = rnum3*Vphys2sim*(180.8892437/mass_neutrino/scalefactor/3.0**0.5)

                xvp(6)=dis(3)*vf + rnum3

              else if (nu_relfd) then !! Use Fermi-Dirac velocity distribution

                call random_number(rnum1) !|vel|
                call random_number(rnum2) !costheta
                call random_number(rnum3) !phi

                !Convert to fermi-dirac
                do l=1,nv
                    if (cdf(2,l).GT.rnum1) then
                        rnum4 = cdf(2,l)-rnum1
                        if ( (l.NE.1) .AND. (rnum1-cdf(2,l-1) .LT. rnum4) ) then
                            rnum1 = cdf(1,l-1)
                        else
                            rnum1 = cdf(1,l)
                        endif
                        exit
                    endif
                enddo

                !convert to fermi-dirac with units rnum1<-rnum1*c*(kT/m)
                rnum1 = rnum1 *(50.2476/mass_neutrino/scalefactor) !*ckT/m = 50.25

                rnum2 = rnum2*2.0-1.0 !convert to range 1,-1
                rnum3 = rnum3*2.0*pi  !convert to range 0,2pi

                xvp(4)=dis(1)*vf + rnum1*(1.0-rnum2**2)**0.5*cos(rnum3)*Vphys2sim
                xvp(5)=dis(2)*vf + rnum1*(1.0-rnum2**2)**0.5*sin(rnum3)*Vphys2sim
                xvp(6)=dis(3)*vf + rnum1*rnum2*Vphys2sim

              else !! Uses linear theory velocity

                xvp(4)=dis(1)*vf
                xvp(5)=dis(2)*vf
                xvp(6)=dis(3)*vf

              endif
#else
             xvp(4)=dis(1)*vf
             xvp(5)=dis(2)*vf
             xvp(6)=dis(3)*vf
#endif

             write(11) xvp

          enddo
       enddo
    enddo

    close(11)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine dm

!!--------------------------------------------------------------!!

  subroutine powerspectrum(pk)
    implicit none
    real, dimension(2,nc)       :: pk

    integer :: i,j,k,kg,mg,jg,ig
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow
#ifdef SLAB
    real, dimension(3,nc,nc_slab) :: pkt
#else
    real, dimension(3,nc,nc_pen+2) :: pkt
#endif
    real, dimension(3,nc) :: pktsum
    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

    pkt=0.0
    pktsum=0.0

    !! Compute power spectrum
#ifdef SLAB
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(kr,kx,ky,kz,kg,k1,k2,w1,w2,pow)
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
                pow=sum((slab(i:i+1,j,k)/ncr**3)**2)
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
    !$omp end parallel do
#else
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
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr .ne. 0) then
                    k1=ceiling(kr)
                    k2=k1+1
                    w1=k1-kr
                    w2=1-w1
                    pow=sum((slab(i:i+1,j,k)/ncr**3)**2)
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
#endif

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
          pk(1:2,k)=4.*pi*(real(k)-1.)**3*pk(1:2,k)
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

#ifdef DEBUG_LOW
    do k=0,nodes-1
      if (rank==k) print *,'rank:',rank,'starting initvar'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    !! Initialize large arrays in parallel
    !$omp parallel default(shared) private(k)
    !$omp do
    do k=0,nc_node_dim+1
       phi(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
#ifdef SLAB
    do k=1,nc_slab
#else
    do k=1,nc_pen+2
#endif
       slab(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nk
       tf(:,k)=0
    enddo
    !$omp end do 
    !$omp do
    do k=1,nc
       pkm(:,k)=0
    enddo
    !$omp end do 
    !$omp do
    do k=1,nc
       pkn(:,k)=0
    enddo
    !$omp end do 
    !$omp end parallel
    
    !! Initialize fftw so that it generates the plans!
    firstfftw=.true.

#ifdef DEBUG_LOW
    do k=0,nodes-1
      if (k==rank) print *,rank,'finished initvar'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

!!------------------------------------------------------------------!!

subroutine mesh_buffer
  implicit none

  integer(4) :: buffer_size
  integer(4) :: tag
  integer(4) :: status(MPI_STATUS_SIZE)

#ifdef DEBUG_BUFFER_MESH
    write(*,*) 'calling mesh_buffer on rank ', rank
#endif

  buffer_size = (nc_node_dim + 2)**2

  if (rank==0) write(*,*) 'buffer_size =', buffer_size


  tag=64

!! send to node in -x

      phi_buf(:,:)=phi(1,:,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(5),tag,cart_neighbor(6), &
                              tag,mpi_comm_cart,status,ierr)

      phi(nc_node_dim+1,:,:)=phi(nc_node_dim+1,:,:)+phi_buf(:,:)

#ifdef DEBUG_BUFFER_MESH
    write(*,*) 'finished -x exchanged on node ', rank
#endif


!! send to node in +x
   
      phi_buf(:,:)=phi(nc_node_dim,:,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(6),tag,cart_neighbor(5), &
                              tag,mpi_comm_cart,status,ierr)

      phi(0,:,:)=phi(0,:,:)+phi_buf(:,:)

!! send to node in -y

      phi_buf(:,:)=phi(:,1,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,nc_node_dim+1,:)=phi(:,nc_node_dim+1,:)+phi_buf(:,:)

!! send to node in +y

      phi_buf(:,:)=phi(:,nc_node_dim,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(4),tag,cart_neighbor(3), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,0,:)=phi(:,0,:)+phi_buf(:,:)

!! send to node in -z
    
      phi_buf(:,:)=phi(:,:,1)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,:,nc_node_dim+1)=phi(:,:,nc_node_dim+1)+phi_buf(:,:)

!! send to node in +z

      phi_buf(:,:)=phi(:,:,nc_node_dim)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,:,0)=phi(:,:,0)+phi_buf(:,:)

  end subroutine mesh_buffer

!!------------------------------------------------------------------!!

  function power(kr,ix,iy)
   implicit none
    real    :: kr
    integer :: ix,iy

    integer :: i,i1,i2
    real    :: x,y,x1,x2,y1,y2
    real    :: power

    i1=1
    i2=nk
    do while (i2-i1 .gt. 1)
       i=(i1+i2)/2
       if (kr .gt. tf(ix,i)) then
          i1=i
       else
          i2=i
       endif
    enddo

    x1=log(tf(ix,i1))
    y1=log(tf(iy,i1))
    x2=log(tf(ix,i2))
    y2=log(tf(iy,i2))
    x =log(kr)
    y =y1+(y2-y1)*(x-x1)/(x2-x1)
    power=exp(y)
    
    return
  end function power

!!------------------------------------------------------------------!!

  function Dgrow(a)
    implicit none
    real, parameter :: om=omegam
    real, parameter :: ol=omegal
    real a
    real Dgrow

    real g,ga,hsq,oma,ola

    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    Dgrow=a*ga/g
    return
  end function Dgrow

!!------------------------------------------------------------------!!

  function vfactor(a)
    implicit none
    real :: a

    real :: H,km,lm
    real :: vfactor

    lm=omegal/omegam
    km=(1-omegam-omegal)/omegam
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H

    return
  end function vfactor

!!------------------------------------------------------------------!!

  function tophat(x)
    implicit none
    real :: x,tophat

    if (x .ne. 0) then
       tophat=3*(sin(x)-cos(x)*x)/x**3
    else
       tophat=1
    endif

    return
  end function tophat

!!------------------------------------------------------------------!!

#ifdef WDM
function wdm_transfer(x)
    !
    ! Returns WDM transfer function using equations (5) and (6) of Schneider, Smith, and Reed (2013)
    !

    implicit none

    real :: x, wdm_transfer

    wdm_transfer = (1. + (alpha*x)**(2.*mu))**(-5./mu)

    return

end function wdm_transfer
#endif

!!------------------------------------------------------------------!!

#ifdef BROKEN_POWER
function broken_power_transfer(x)
    !
    ! Applies a broken power law at a small scale pivot (k0) as an alternative to WDM.
    !

    implicit none

    real :: x, broken_power_transfer

    if (x < k0) then
        broken_power_transfer = 1.
    else
        broken_power_transfer = (x/k0)**(dns/2.)    
    endif

    return

end function broken_power_transfer
#endif

end program dist_init 
