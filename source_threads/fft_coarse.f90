!*******************************************************
! fftw3ds.f90 - cubic decomposition subroutines for fftw

  subroutine pack_slab
!! pack cubic data into slab decomposition for fftw transform
    use omp_lib
    implicit none
    include 'mpif.h'
    include 'cubepm.fh'

    integer(4) :: i,j,k,i0,j0,i1,j1,k1
!    integer(4) :: slab_slice,num_elements,tag,rtag
    integer(4) :: slab_slice,tag,rtag
    integer(4) :: num_elements !possibly should be double so as not to 
                               !overflow for large runs, but integer*8 
                               !might be unsupported by MPI

    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status

    num_elements = nc_node_dim * nc_node_dim * nc_slab

!! swap data

    do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_slice = i + j * nodes_dim
        tag = rank**2
        rtag= slab_neighbor(i,j)**2
        call mpi_isend(rho_c(1,1,slab_slice*nc_slab + 1), num_elements, &
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

  subroutine unpack_slab 
!! unpack slab data into cubic decomposition following fftw transform
    use omp_lib
    implicit none
    include 'mpif.h'
    include 'cubepm.fh'
    
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
        call mpi_irecv(rho_c(1,1,slab_slice * nc_slab +1), &
                       num_elements, mpi_real, slab_neighbor(i,j), rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
      enddo
    enddo

    call mpi_waitall(2*nodes_dim**2,requests, wait_status, ierr)

  end subroutine unpack_slab

  subroutine cubepm_fftw(command)
!! calculate fftw transform
!! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
    use omp_lib
    use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'
    include 'fftw3-mpi.f03'
    include 'cubepm.fh'
    !include 'fftw_f77.i'

    !integer(4), parameter :: order=FFTW_NORMAL_ORDER ! FFTW_TRANSPOSED_ORDER

    integer(4) :: i
    integer(4) :: command 

    type(C_PTR),save :: plan, iplan

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'starting fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

! initialize plan variables for fftw

    if (firstfftw) then

      call fftwf_mpi_init
      !call fftwf_mpi_plan_dft_r2c_3d(plan, nc_dim, nc_dim, nc_dim, slab, slab,&
      !      mpi_comm_world, FFTW_ESTIMATE)


      plan =  fftwf_mpi_plan_dft_r2c_3d(nc_dim, nc_dim, nc_dim, slab, slab_cmplx, mpi_comm_world, FFTW_ESTIMATE)
      iplan =  fftwf_mpi_plan_dft_c2r_3d( nc_dim, nc_dim, nc_dim, slab_cmplx, slab, mpi_comm_world, FFTW_ESTIMATE)
      !plan =  fftwf_mpi_plan_dft_r2c_3d(nc_dim, nc_dim, nc_dim, slab, slab, mpi_comm_world, FFTW_ESTIMATE)
      !iplan =  fftwf_mpi_plan_dft_c2r_3d( nc_dim, nc_dim, nc_dim, slab, slab, mpi_comm_world, FFTW_ESTIMATE)



      !call fftwf_mpi_plan_dft_c2r_3d(iplan, nc_dim, nc_dim, nc_dim, slab, slab,&
      !      mpi_comm_world, FFTW_ESTIMATE)
      !call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nc_dim, &
      !      nc_dim,nc_dim, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE) !_MEASURE)
      !call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nc_dim, &
      !      nc_dim,nc_dim, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE) !_MEASURE)
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
      !call rfftwnd_f77_mpi(plan,1,slab,slab_work,1,order)
      call sfftw_execute_dft_r2c(plan, slab, slab_cmplx)
      !slab =  sfftw_execute_dft_r2c(plan, slab)
    else
      !call rfftwnd_f77_mpi(iplan,1,slab, slab_work,1,order)
      call sfftw_execute_dft_c2r(iplan, slab_cmplx, slab)
      !slab =  sfftw_execute_dft_c2r(iplan, slab)
      slab=slab/(real(nc_dim)*real(nc_dim)*real(nc_dim))
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

      !call rfftwnd_f77_mpi_destroy_plan(iplan)
      !call rfftwnd_f77_mpi_destroy_plan(plan)
      call sfftw_destroy_plan(plan)
      call sfftw_destroy_plan(iplan)
      call fftwf_mpi_cleanup
    endif
  
  end subroutine cubepm_fftw
