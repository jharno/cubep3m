!!
!! cic_power_p3dfft.f90
!!
!! Program to compute the dark matter power spectra from cubep3m particle checkpoints.
!!  
!! * Using FFTW on the SciNet GPC compile with:
!!   mpif90 -shared-intel -fpp -g -O3 -DNGP -mt_mpi cic_power.f90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_power -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Using MKL on the SciNet GPC compile with:
!!   mpiifort -shared-intel -fpp -g -O3 -xhost -i_dynamic -mkl -mcmodel=medium -DNGP -mt_mpi 
!!        cic_power.f90 -I$P3DFFT_INC -I$MKL_FFTW_INC -o ngp_power -L$P3DFFT_LIB -L$MKL_FFTW_LIB -lp3dfft -lmkl_intel_lp64
!!
!! * Using FFTW on the SciNet BGQ compile with:
!!   mpif90 -q64 -O3 -qhot -qarch=qp -qtune=qp -WF,-DNGP cic_power.F90 -I$SCINET_FFTW_INC 
!!        -I$P3DFFT_INC -o ngp_power -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
!!
!! * Optional flags:
!!   -DNGP: Uses NGP interpolation for binning of the power spectrum. 
!!   -DSLAB: Alternatively run with FFTW slab decomposition instead of P3DFFT pencil decomposition.
!!   -Dwrite_den: Writes gridded density fields to binary files. 
!!   -Dwrite_poisson: Writes the gridded Poisson field to binary files.
!!   -DMOVE_GRID_BACK: Undos the random grid displacements by reading offsetfile. 
!!   -DKAISER: Adjusts for redshift space distortions.
!!   -DPLPLOT: Plots power spectra at end.
!!   -DDEBUG: Output useful debugging information.

program cic_power 

  implicit none

  include 'mpif.h'
  include '../../parameters'

  logical, parameter :: correct_kernel=.false.
  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints_high'

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

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local
  integer(8) :: plan, iplan
  logical :: firstfftw

  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(3,np_buffer) :: xp_buf
  real, dimension(3*np_buffer) :: send_buf, recv_buf
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: den 
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: den_buf 

  !! Power spectrum arrays
  real, dimension(3,nc) :: pkdm
  real, dimension(3,nc) :: poisson
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
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1) :: recv_cube
#endif

  !! parameters for moving the displaced mesh back to its unperturbed state
  character(len=*), parameter :: offsetfile=cubepm_root//"shake_offsets.dat"
  real, dimension(max_checkpoints, 3) :: shake_offsets

  !! Equivalence arrays to save memory
#ifdef SLAB
  equivalence (den,slab_work,recv_cube,xp_buf)
#else
  equivalence (den,recv_cube,xp_buf) 
#endif
  equivalence (xvp,slab,cube)

  !! Common block
#ifdef PLPLOT
 common xvp,send_buf,den_buf,den,recv_buf,pkdm,pkplot
#else
  common xvp,send_buf,den_buf,den,recv_buf,pkdm,poisson
#endif

!!---start main--------------------------------------------------------------!!

  call mpi_initialize
  if (rank == 0) call writeparams
  firstfftw=.true.  ! initialize fftw so that it generates the plans
  call read_checkpoint_list
#ifdef MOVE_GRID_BACK
  call read_shake_offsets_list
#endif
  do cur_checkpoint=1,num_checkpoints
    call initvar
    call read_particles
    call pass_particles
    call darkmatter
    call PoissonNoise
    if (rank == 0) call writepowerspectra
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

subroutine read_shake_offsets_list
    !
    ! Read in list of cumulative shake offsets for each checkpoint in order to
    ! move grid back to unperturbed state.
    !

    implicit none

    integer :: k, j, i, fstat

    open(11, file=offsetfile, status='old', iostat=fstat)

    do i = 1, num_checkpoints

        read(11,*) shake_offsets(i, :)

    enddo

    close(11)

    !! Make the shake offset within the range [-nc, +nc]

    do i = 1, num_checkpoints
        do j = 1, 3
            if (shake_offsets(i, j) < 0.) then
                do k = 1, 100
                    if (shake_offsets(i, j) < -nc) then
                        shake_offsets(i, j) = shake_offsets(i, j) + nc
                    else
                        exit
                    endif
                enddo
            else
                do k = 1, 100
                    if (shake_offsets(i, j) > nc) then
                        shake_offsets(i, j) = shake_offsets(i, j) - nc
                    else
                        exit
                    endif
                enddo
            endif
        enddo
    enddo

    return

end subroutine read_shake_offsets_list

!!---------------------------------------------------------------------------!!

  subroutine read_particles
    implicit none
    
    real z_write,np_total
    integer i,j,fstat, blocksize, nplow, nphigh, num_writes
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name

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

    if(z_write .eq. z_i) then
       check_name=ic_path//'xv'//rank_string(1:len_trim(rank_string))//'.ic'
    else
       check_name=output_path//z_string(1:len_trim(z_string))//'xv'// &
               rank_string(1:len_trim(rank_string))//'.dat'
    endif

!! open checkpoint    

    open(unit=21,file=check_name,status="old",iostat=fstat,access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint'
      write(*,*) 'rank',rank,'file:',check_name
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! read in checkpoint header data
    if(z_write .eq. z_i)then
       read(21) np_local
       a = 1.0/(1.0 + z_i)
    else
    read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
               sim_projection,sim_halofind,mass_p
    endif

    if (np_local > max_np) then
      write(*,*) 'too many particles to store'
      write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! tally up total number of particles
    call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
                         mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'number of particles =', int(np_total,8)

    do j=1, np_local
        read(21) xvp(:,j)
    enddo

    close(21)

#ifdef MOVE_GRID_BACK
    if (rank == 0) write(*,*) "Moving grid back with shake offset = ", shake_offsets(cur_checkpoint, :)
    do j = 1, np_local
        xvp(1, j) = xvp(1, j) - shake_offsets(cur_checkpoint, 1)
        xvp(2, j) = xvp(2, j) - shake_offsets(cur_checkpoint, 2)
        xvp(3, j) = xvp(3, j) - shake_offsets(cur_checkpoint, 3)
    enddo
#endif
 
#ifdef KAISER

    !Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !Converting seconds into simulation time units
    !cancels the H0...
    
    xvp(3,:)=xvp(3,:) + xvp(6,:)*1.5/sqrt(a*(1+a*(1-omega_m-omega_l)/omega_m + omega_l/omega_m*a**3))  

    call pass_particles

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
        recv_cube(:,:,:, pen_slice) = slab(i0:i1,:,:)
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

  subroutine writepowerspectra
    implicit none
    integer      :: k
#ifdef PLPLOT
    integer :: kp
#endif
    real         :: kr
    character*180 :: fn
    character*5  :: prefix
    character*7  :: z_write
    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is dm d2(k)
    !! 3rd is standard deviation
    !! 4th is Poisson d2(k)
    !! 5th is standard deviation on Poisson d2(k)

    write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
    z_write=adjustl(z_write)
    
#ifdef NGP 
    prefix='ngpps'
#else
    prefix='cicps'
#endif

#ifdef KAISER
    fn=output_path//z_write(1:len_trim(z_write))//prefix//'-RSD.dat' 
#else
    fn=output_path//z_write(1:len_trim(z_write))//prefix//'.dat' 
#endif

    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
#ifdef NGP
       write(11,*) pkdm(3,k-1),pkdm(1:2,k-1) ,poisson(1:2,k-1)
#else
       write(11,*) pkdm(3,k),pkdm(1:2,k),poisson(1:2,k)
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

  subroutine darkmatter
    implicit none
    integer :: i,j,k, fstat
    integer :: i1,j1,k1
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,z_write
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name


    real time1,time2
    call cpu_time(time1)

    !! Initialized density field to be zero
    !! could do OMP loop here
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo

    !! Assign masses to grid to compute dm power spectrum
    call cicmass

    !! have to accumulate buffer density 
    call mesh_buffer
    cube=den(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)

#ifdef write_den
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
    check_name=output_path//z_string(1:len_trim(z_string))//'den'// &
               rank_string(1:len_trim(rank_string))//'-rsd.dat'
#else 
    check_name=output_path//z_string(1:len_trim(z_string))//'den'// &
               rank_string(1:len_trim(rank_string))//'.dat'
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
      write(*,*) 'DM min    ',dmint
      write(*,*) 'DM max    ',dmaxt
      write(*,*) 'Delta sum ',real(dsum,8)
      write(*,*) 'Delta var ',real(dvar,8)
      write(*,*)
    endif
 
    !! Forward FFT dm delta field
    call cp_fftw(1)

    !! Compute dm power spectrum
    call powerspectrum(slab,pkdm)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine darkmatter

!!------------------------------------------------------------------!!

  subroutine PoissonNoise
    implicit none
    integer :: i,j,k, fstat
    integer :: i1,j1,k1
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,z_write
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name
    real time1,time2

    call cpu_time(time1)

    !! Initialized density field to be zero
    !! could do OMP loop here
    do k=0,nc_node_dim+1
       den(:,:,k)=0
    enddo

    ! Randomize positions across all nodes, hence np_local should be equal
    np_local = np_node_dim**3 

    ! Assign particles to random positions between 0 and nc_node_dim
    call random_number(xvp(1:3,:np_local))
    xvp(1:3,:np_local) = xvp(1:3,:np_local)*nc_node_dim

    !! Assign masses to grid to compute dm power spectrum
    call cicmass

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
    check_name=output_path//z_string(1:len_trim(z_string))//'den-poisson'// &
               rank_string(1:len_trim(rank_string))//'-rsd.dat'
#else 
    check_name=output_path//z_string(1:len_trim(z_string))//'den-poisson'// &
               rank_string(1:len_trim(rank_string))//'.dat'
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

    !! Compute dm power spectrum
    call powerspectrum(slab,poisson)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine PoissonNoise
!------------------------------------------------------------!

  subroutine pass_particles
    implicit none

    integer i,pp,np_buf,np_exit,npo,npi
    integer*8 np_final
    real x(3),lb,ub
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03

    lb=0.0
    ub=real(nc_node_dim)

    np_buf=0
    pp=1
    do
      if (pp > np_local) exit
      x=xvp(:3,pp)
      if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
          x(3) < lb .or. x(3) >= ub ) then
!        write (*,*) 'PARTICLE OUT',xv(:,pp)
        np_buf=np_buf+1
        if (np_buf > np_buffer) then
          print *,rank,'np_buffer =',np_buffer,'exceeded - np_buf =',np_buf
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif 
        xp_buf(:,np_buf)=xvp(:3,pp)
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

    if (rank == 0) print *,'total exiting particles =',np_exit

! pass +x

    tag=11 
    npo=0
    pp=1
    do 
      if (pp > np_buf) exit
      if (xp_buf(1,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
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
        xvp(:3,np_local)=x
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
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(1,np_buf+pp)=min(xp_buf(1,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
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
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(2,np_buf+pp)=max(xp_buf(2,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
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
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(2,np_buf+pp)=min(xp_buf(2,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
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
      if (xp_buf(3,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(3,np_buf+pp)=max(xp_buf(3,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
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
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
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

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(3,np_buf+pp)=min(xp_buf(3,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
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
      print *,'total particles =',int(np_final,8)
      if (np_final /= (int(np,8))**3) then
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

  end subroutine pass_particles

!------------------------------------------------------------!

  subroutine cicmass
    implicit none
    real, parameter :: mp=(ncr/np)**3

    integer :: i,i1,i2,j1,j2,k1,k2
    real    :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2,vf,v(3)

    do i=1,np_local
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

  subroutine powerspectrum(delta,pk)
    implicit none
    real, dimension(3,nc)       :: pk
#ifdef SLAB
    real, dimension(nc+2,nc,nc_slab) :: delta
#else
    real, dimension(nc, nc_node_dim, nc_pen+2) :: delta
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

    real, dimension(nc) :: kcen, kcount
    real    :: kavg
    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    pkt=0.0
    pktsum=0.0

    kcen(:)   = 0.
    kcount(:) = 0.

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
          pk(1:2,k)=pktsum(1:2,k)/pktsum(3,k)
          pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))

          kavg = kcen(k) / kcount(k)
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
    enddo
    do k=1,np_buffer
       xp_buf(:,k)=0
    enddo
    do k=1,3*np_buffer
       recv_buf(k)=0
    enddo
    do k=1,nc
       pkdm(:,k)=0
    enddo    
    do k=1,nc
       poisson(:,k)=0
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

end program cic_power 
