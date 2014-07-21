!*******************************************************
! p3dfft_coarse.f90 - pensil decomposition subroutines 
! for fft on the coarse mesh. 
! Written by Joachim Harnois-Deraps, 2013/02
! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
!*******************************************************

  subroutine cubepm_fftw(command)

    use p3dfft

    implicit none
    include 'mpif.h'
#    include "cubepm.fh"

    integer(4) :: i
    integer(4) :: command 

!-------------------------------------
! initialize pencil variables for p3dfft

    if (firstfftw) then

       ! call mpi_dims_create(nodes_dim, DECOMP, dims, ierr) 
       ! This previous call looks for optimal decomposition, not necessarily optimal for us.
       ! We went for simplicity and hard coded dims in cubepm.par. 
      
       pen_dims = (/dim_y,dim_z/)
       call p3dfft_setup(pen_dims, nc_dim, nc_dim, nc_dim, .true.)
       call p3dfft_get_dims(istart, iend, isize, 1, mypadd)
       call p3dfft_get_dims(fstart, fend, fsize, 2)

       ! Now each node knows about the actual coordinates of the pencil it is associated with.
       call mpi_barrier(mpi_comm_world, ierr)

       firstfftw=.false.

    endif

!-----------------------

    if (command /= 0) then

      !! Call pack routine if we are going forward
      if (command > 0) then 
        call pack_pencils
      endif

      !! Do the FFT
      if (command > 0) then 
        call ftran_r2c(slab, slab, "fft")
      else 
        call btran_c2r(slab, slab, "tff")
        slab=slab/(real(nc_dim)*real(nc_dim)*real(nc_dim))
      endif

     !! Unpack the pencil data 
      if (command < 0) then 
         call unpack_pencils
      endif
    
    else !! If command == 0 we delete the plans
        call p3dfft_clean
    endif
  
  end subroutine cubepm_fftw

!------------------------
  subroutine pack_pencils
    !
    ! Pack cubic data into pencils for p3dfft transform.
    !

    use omp_lib
    implicit none
    include 'mpif.h'
#    include "cubepm.fh"

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
            call mpi_isend(rho_c(1,1, pen_slice*nc_pen + nc_pen_break + 1), num_elements, &
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
!------------------
  subroutine unpack_pencils
    !
    ! Unpack data from the pencils back into the cubic decompisition following
    ! p3dfft transform.
    !

    use omp_lib
    implicit none
    include 'mpif.h'
#    include "cubepm.fh"

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
            call mpi_irecv(rho_c(1,1,pen_slice*nc_pen + nc_pen_break + 1), &
                           num_elements, mpi_real, pen_neighbor_to(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim,requests, wait_status, ierr)

    enddo

  end subroutine unpack_pencils

