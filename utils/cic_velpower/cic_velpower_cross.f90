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

  !! Power spectrum arrays
  real, dimension(3, 3, nc) :: pkdimdm0, pkdimdm1, pkdimdm2
  real, dimension(3, nc) :: pkdm

  !! Fourier transform arrays
  real, dimension(nc_node_dim, nc_node_dim, nc_node_dim) :: cube1
  real, dimension(nc_node_dim, nc_node_dim, nc_node_dim) :: cube2
  real, dimension(nc_node_dim, nc_node_dim, nc_slab, 0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2, nc, nc_slab) :: slab1, slab_work
  real, dimension(nc+2, nc, nc_slab) :: slab2

  !! Array containing (x, y, z) components of the momentum density field
  real, dimension(3, 0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momden1
  real, dimension(3, 0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: momden2
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: momden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: momden_recv_buff

  !! Array containing mass density field
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: massden1
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1, 0:nc_node_dim+1) :: massden2
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: massden_send_buff
  real, dimension(0:nc_node_dim+1, 0:nc_node_dim+1) :: massden_recv_buff

  !! Equivalence arrays to save memory
  equivalence (slab_work, recv_cube) 

  !! Common block
  common slab_work, cube1, cube2, slab1, slab2, pkdimdm0, pkdimdm1, pkdimdm2, pkdm

! -------------------------------------------------------------------------------------------------------
! MAIN
! -------------------------------------------------------------------------------------------------------

  call mpi_initialize

  if (rank == 0) call writeparams

  firstfftw = .true.  ! initialize fftw so that it generates the plans

  call read_checkpoint_list

  do cur_checkpoint = 1, num_checkpoints

    call initvar

    !
    ! Compute power spetrum of first velocity field
    !

    call read_velocity_field(1)

    do cur_dimension = 1, 3

        call darkmatter(1)

    enddo

    if (rank == 0) call writepowerspectra(1)

    !
    ! Compute power spetrum of second velocity field
    !

    call read_velocity_field(2)

    do cur_dimension = 1, 3

        call darkmatter(2)

    enddo

    if (rank == 0) call writepowerspectra(2)

    !
    ! Compute cross-power spectrum
    !

    do cur_dimension = 1, 3

        call darkmatter(0)
  
    enddo

    if (rank == 0) call writepowerspectra(0)

  enddo

  call cp_fftw(0, 0)
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

    !! Momentum and matter density arrays
    do k = 0, nc_node_dim + 1
        momden1(:, :, :, k) = 0.
        massden1(:, :, k) = 0.
        momden2(:, :, :, k) = 0.
        massden2(:, :, k) = 0.
    enddo

    !! Fourier transform arrays
    do k = 1, nc_slab
       slab_work(:, :, k)=0
    enddo

    !! Power spectrum arrays
    do k = 1, nc
        pkdimdm0(:, :, k) = 0.
        pkdm(:, k) = 0.
        pkdimdm1(:, :, k) = 0.
        pkdimdm2(:, :, k) = 0.
    enddo

    return

end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine read_velocity_field(command)
    !
    ! Read binary file containing velocity field for the given rank.
    !

    implicit none

    integer :: m, i, j, k, fstat
    character(len=180) :: fn
    character(len=7)   :: z_write
    character(len=4)   :: rank_string
    character(len=1)   :: dim_string

    integer :: command ! should be 1 or 2

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

        if (command == 1) then

            fn = output_path//z_write(1:len_trim(z_write))//&
                 "vel"//dim_string//&
                 rank_string(1:len_trim(rank_string))//".bin"

            open(unit=11, file=fn, status="old", iostat=fstat, form="binary")

            do k = 1, nc_node_dim
                do j = 1, nc_node_dim

                    read(11) momden1(m, 1:nc_node_dim, j, k)

                enddo
            enddo

            close(11)

        else if (command == 2) then

            fn = output_path//z_write(1:len_trim(z_write))//&
                 "vel"//dim_string//&
                 rank_string(1:len_trim(rank_string))//".bin"

            open(unit=11, file=fn, status="old", iostat=fstat, form="binary")

            do k = 1, nc_node_dim
                do j = 1, nc_node_dim

                    read(11) momden2(m, 1:nc_node_dim, j, k) 

                enddo
            enddo

            close(11)

        endif

    enddo

    return

end subroutine read_velocity_field

! -------------------------------------------------------------------------------------------------------

  subroutine pack_slab(command)
!! pack cubic data into slab decomposition for fftw transform
    implicit none

    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag
      
    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status
        
    integer :: command ! should be 1 or 2

    num_elements = nc_node_dim * nc_node_dim * nc_slab
                       
!! swap data           
        
    do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_slice = i + j * nodes_dim
        tag = rank**2
        rtag= slab_neighbor(i,j)**2

        if (command == 1) then
            call mpi_isend(cube1(1,1,slab_slice*nc_slab + 1), num_elements, &
                       mpi_real, slab_neighbor(i,j), tag, mpi_comm_world, &
                       requests(slab_slice+1),ierr)

        else if (command == 2) then
            call mpi_isend(cube2(1,1,slab_slice*nc_slab + 1), num_elements, &
                       mpi_real, slab_neighbor(i,j), tag, mpi_comm_world, &
                       requests(slab_slice+1),ierr)
        endif

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

        if (command == 1) then
            slab1(i0:i1,j0:j1,:) = recv_cube(:,:,:,slab_slice)
        else if (command == 2) then
            slab2(i0:i1,j0:j1,:) = recv_cube(:,:,:,slab_slice)
        endif

      enddo
    enddo
      
  end subroutine pack_slab
    
! -------------------------------------------------------------------------------------------------------

  subroutine unpack_slab(command)
!! unpack slab data into cubic decomposition following fftw transform
    implicit none
      
    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag
    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status
      
    integer :: command

!! place data in the recv_cube buffer
      
    do j = 0, nodes_dim - 1
      j0 = j * nc_node_dim + 1
      j1 = (j + 1) * nc_node_dim
      do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        slab_slice = i + j * nodes_dim

        if (command == 1) then
            recv_cube(:,:,:,slab_slice) = slab1(i0:i1,j0:j1,:)
        else if (command == 2) then
            recv_cube(:,:,:,slab_slice) = slab2(i0:i1,j0:j1,:)
        endif

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
    
        if (command == 1) then
            call mpi_irecv(cube1(1,1,slab_slice * nc_slab +1), &
                       num_elements, mpi_real, slab_neighbor(i,j), rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
        else if (command == 2) then
            call mpi_irecv(cube2(1,1,slab_slice * nc_slab +1), &
                       num_elements, mpi_real, slab_neighbor(i,j), rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
        endif
      enddo
    enddo

    call mpi_waitall(2*nodes_dim**2,requests, wait_status, ierr)

  end subroutine unpack_slab

! -------------------------------------------------------------------------------------------------------

  subroutine cp_fftw(command, command2)
!! calculate fftw transform
!! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
    implicit none
    include 'fftw_f77.i'

    integer(4), parameter :: order=FFTW_NORMAL_ORDER ! FFTW_TRANSPOSED_ORDER

    integer(4) :: i
    integer(4) :: command, command2

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
      if (command > 0) call pack_slab(command2)

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'finished forward slab pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    if (command2 == 1) then
        if (command > 0) then
          call rfftwnd_f77_mpi(plan,1,slab1,slab_work,1,order)
        else
          call rfftwnd_f77_mpi(iplan,1,slab1,slab_work,1,order)
          slab1=slab1/real(nc*nc*nc)
        endif
    else if (command2 == 2) then
        if (command > 0) then
          call rfftwnd_f77_mpi(plan,1,slab2,slab_work,1,order)
        else
          call rfftwnd_f77_mpi(iplan,1,slab2,slab_work,1,order)
          slab2=slab2/real(nc*nc*nc)
        endif
    endif

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'finished fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

!! unpack the slab data

      if (command < 0) call unpack_slab(command2)

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
    character*8  :: prefix
    character*7  :: z_write
    integer(4) :: command

    !
    ! Determine name of output file
    !

    write(z_write,'(f7.3)') z_checkpoint(cur_checkpoint)
    z_write=adjustl(z_write)
    
#ifdef NGP 
    if (command == 0) then
        prefix = 'ngpvps12'
    else if (command == 1) then
        prefix = 'ngpvps11'
    else if (command == 2) then
        prefix = 'ngpvps22'
    endif
#else
    if (command == 0) then
        prefix = 'cicvps12'
    else if (command == 1) then
        prefix = 'cicvps11'
    else if (command == 2) then
        prefix = 'cicvps22'
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
                pkdm(1, i) = pkdm(1, i) + pkdimdm0(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkdimdm0(j, 2, i)
            enddo
            pkdm(3, i) = pkdimdm0(1, 3, i)
        enddo

    else if (command == 1) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkdimdm1(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkdimdm1(j, 2, i)
            enddo
            pkdm(3, i) = pkdimdm1(1, 3, i)
        enddo

    else if (command == 2) then

        !! Sum over all three dimensions 
        do i = 1, nc
            do j = 1, 3
                pkdm(1, i) = pkdm(1, i) + pkdimdm2(j, 1, i)
                pkdm(2, i) = pkdm(2, i) + pkdimdm2(j, 2, i)
            enddo
            pkdm(3, i) = pkdimdm2(1, 3, i)
        enddo

    endif

    !
    ! Write to output file with column ordering [k, p(k), sigma(k)]
    !

    !! Will be in whatever velocity units the binary files were in

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
    real    :: d1, ddmin1, ddmax1, sum_dm1, sum_dm_local1, dmint1, dmaxt1
    real*8  :: dsum1, dvar1, dsumt1, dvart1
    real    :: d2, ddmin2, ddmax2, sum_dm2, sum_dm_local2, dmint2, dmaxt2
    real*8  :: dsum2, dvar2, dsumt2, dvart2
    real, dimension(3) :: dis

    integer(4) :: command ! should be 0, 1, or 2 

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    !
    ! Initialize FFT arrays to zero 
    !

    do k = 1, nc_node_dim
        cube1(:, :, k) = 0.
        cube2(:, :, k) = 0.
    enddo

    do k = 1, nc_slab
       slab1(:, :, k) = 0.
       slab2(:, :, k) = 0.
    enddo

    !
    ! Assign data to be transformed
    !

    if (command == 0) then !! Computing cross-power spectrum 

        cube1(:, :, :) = momden1(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim) 
        cube2(:, :, :) = momden2(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)

    else if (command == 1) then !! Computing power spectrum of first velocity field 

        cube1(:, :, :) = momden1(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
        cube2(:, :, :) = momden1(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)

    else if (command == 2) then !! Computing power spectrum of second velocity field 

        cube1(:, :, :) = momden2(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)
        cube2(:, :, :) = momden2(cur_dimension, 1:nc_node_dim, 1:nc_node_dim, 1:nc_node_dim)

    endif

    !
    ! Calculate some statistics
    !

    sum_dm_local1 = sum(cube1) 
    call mpi_reduce(sum_dm_local1, sum_dm1, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
    if (rank == 0) print *, "CUBE1 total sum = ", sum_dm1, " command = ", command

    sum_dm_local2 = sum(cube2)
    call mpi_reduce(sum_dm_local2, sum_dm2, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
    if (rank == 0) print *, "CUBE2 total sum = ", sum_dm2, " command = ", command

    ddmin1 = 0
    ddmax1 = 0
    dsum1 = 0
    dvar1 = 0

    ddmin2 = 0
    ddmax2 = 0
    dsum2 = 0
    dvar2 = 0

    do k = 1, nc_node_dim
       do j = 1, nc_node_dim
          do i = 1, nc_node_dim

             d1 = cube1(i, j, k)
             dsum1 = dsum1 + d1
             dvar1 = dvar1 + d1*d1
             ddmin1 = min(ddmin1, d1)
             ddmax1 = max(ddmax1, d1)

             d2 = cube2(i, j, k)
             dsum2 = dsum2 + d2
             dvar2 = dvar2 + d2*d2
             ddmin2 = min(ddmin2, d2)
             ddmax2 = max(ddmax2, d2)

          enddo
       enddo
    enddo

    call mpi_reduce(dsum1, dsumt1, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(dvar1, dvart1, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(ddmin1, dmint1, 1, mpi_real, mpi_min, 0, mpi_comm_world, ierr)
    call mpi_reduce(ddmax1, dmaxt1, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)

    call mpi_reduce(dsum2, dsumt2, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(dvar2, dvart2, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(ddmin2, dmint2, 1, mpi_real, mpi_min, 0, mpi_comm_world, ierr)
    call mpi_reduce(ddmax2, dmaxt2, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)

    if (rank == 0) then

      dsum1 = dsumt1 / nc**3
      dvar1 = sqrt(dvart1 / nc**3)
      dsum2 = dsumt2 / nc**3
      dvar2 = sqrt(dvart2 / nc**3)

      write(*,*)
      write(*,*) 'Darkmatter command ', command
      write(*,*) 'Cube1 min    ', dmint1
      write(*,*) 'Cube1 max    ', dmaxt1
      write(*,*) 'Cube1 sum ', real(dsum1)
      write(*,*) 'Cube1 var ', real(dvar1)
      write(*,*) 'Cube2 min    ', dmint2
      write(*,*) 'Cube2 max    ', dmaxt2
      write(*,*) 'Cube2 sum ', real(dsum2)
      write(*,*) 'Cube2 var ', real(dvar2)
      write(*,*)

    endif

    ! 
    ! Forward FFT dm delta field
    !    

    call cp_fftw(1, 1) ! packs cube1 into slab1
    call cp_fftw(1, 2) ! packs cube2 into slab2

    !
    ! Compute power spectrum
    !
   
    if (command == 0) then

        call powerspectrum(slab1, slab2, pkdimdm0(cur_dimension, :, :))
 
    else if (command == 1) then
    
        call powerspectrum(slab1, slab2, pkdimdm1(cur_dimension, :, :))

    else if (command == 2) then

        call powerspectrum(slab1, slab2, pkdimdm2(cur_dimension, :, :))

    endif

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished darkmatter ... elapsed time = ", time2-time1

    return

end subroutine darkmatter

! -------------------------------------------------------------------------------------------------------

  subroutine powerspectrum(delta1, delta2, pk)
    implicit none
    real, dimension(3, nc)       :: pk
    real, dimension(nc+2, nc, nc_slab) :: delta1, delta2
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
                pow=sum((delta1(i:i+1,j,k)/ncr**3)*(delta2(i:i+1,j,k)/ncr**3))/kernel**4
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

            endif
        enddo
    endif

    call mpi_bcast(pk,3*nc,mpi_real,0,mpi_comm_world,ierr)

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*, *) "Finished powerspectrum ... elapsed time = ", time2-time1

    return

  end subroutine powerspectrum

! -------------------------------------------------------------------------------------------------------

end program cic_velpower 
