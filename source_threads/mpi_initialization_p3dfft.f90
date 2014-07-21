!! initialize mpi parameters
  subroutine mpi_initialize

    implicit none
    include 'mpif.h'
#    include "cubepm.fh"

    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder

!! set up global mpi communicator

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (nodes_returned /= nodes ) then
      write(*,*) 'cubepm compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'cubepm nodes=',nodes 
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    if (mod(nc_dim,nodes_dim**2) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into pencils'
      write(*,*) 'nc=',nf_physical_dim, 'nc coarse=',nc_dim,'nodes_dim**2=',nodes_dim**2,'mod(nc,nodes_dim**2)=', mod(nc_dim,nodes_dim**2)
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#ifdef DIAG
    if (rank==0) then
      write(*,*) 'cubepm running on',nodes,'nodes with',cores,'cores per node'
      write(*,*) 'using cubic distribution:',nodes_dim,'nodes per dimension'
      write(*,*) nf_physical_dim,'cells in mesh'
    endif
#endif

!! calculate coordinates within slab for cube processes

    slab_coord(3) = rank / nodes_slab
    slab_rank = rank - slab_coord(3) * nodes_slab
    slab_coord(2) = slab_rank / nodes_dim
    slab_coord(1) = slab_rank - slab_coord(2) * nodes_dim

    do j = 0, nodes_dim - 1
        pen_neighbor_to(j) = nodes_slab*slab_coord(3) + slab_coord(2) + j*nodes_dim
        pen_neighbor_fm(j) = nodes_slab*slab_coord(3) + j + nodes_dim*slab_coord(1)
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
    if (i==rank) write(*,'(8i4)') rank,cart_rank,cart_neighbor
    call mpi_barrier(mpi_comm_world,ierr)
  enddo 
#endif

  end subroutine mpi_initialize
