#include "./preprocessor"
module mMPI
	use Parameters
	use Variables
	implicit none

	include 'mpif.h'

	integer, dimension(6) :: cart_neighbor
	integer, dimension(3) :: slab_coord, cart_coords
	integer :: slab_rank, mpi_comm_cart, cart_rank, rank, rank_io, ierr
	character(len=8) :: rank_s

	real(8) :: time_start, time_end

	public :: rank
	public :: rank_s
	public :: start_mpi
	public :: end_mpi
	public :: mpi_error_stop

	public :: slab_coord

contains

	subroutine start_mpi
		implicit none
		integer :: j,nodes_returned, ndim
		integer, dimension(3) :: dims
		logical :: reorder
		logical, dimension(3) :: periodic

#if VERBOSITY > 1
		write(*,*) "[Mod - mMPI] Beginning subroutine start_mpi"
#endif

		call mpi_init(ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_init')

		call mpi_comm_size(mpi_comm_world, nodes_returned, ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_comm_size')

		if (nodes_returned /= nodes) then
			write(*,*) nodes_returned, nodes
			call mpi_error_stop('Error in subroutine &
			&start_mpi: unexpected number of nodes returned')
		end if

		call mpi_comm_rank(mpi_comm_world,rank,ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_comm_rank')

		slab_coord(3) = rank / nodes_slab
		slab_rank = rank - slab_coord(3) * nodes_slab
		slab_coord(2) = slab_rank / nodes_dim
		slab_coord(1) = slab_rank - slab_coord(2) * nodes_dim

		dims(:) = nodes_dim
		periodic(:) = .true.
		reorder = .false.
		ndim = 3

		call mpi_cart_create(mpi_comm_world, ndim, dims, periodic, reorder, mpi_comm_cart, ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_cart_create')

		call mpi_comm_rank(mpi_comm_cart, cart_rank, ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_comm_rank')

		call mpi_cart_coords(mpi_comm_cart, cart_rank, ndim, cart_coords, ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_cart_coords')

if (nodes_dim /= Nglobal_dim) then
  rank_io=Nglobal_dim**2*cart_coords(1)+Nglobal_dim*cart_coords(2)+cart_coords(3)
  write(rank_s,'(i8)') rank_io
else
  write(rank_s,'(i8)') rank
endif

		do j=0, ndim-1
			call mpi_cart_shift(mpi_comm_cart, j, 1, cart_neighbor(2*(j+1)-1), &
				& cart_neighbor(2*(j+1)), ierr)
			if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_cart_shift')
		end do

		time_start = mpi_wtime(ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_wtime')

		call omp_set_num_threads(Nomp)

#if VERBOSITY > 1
		write(*,*) "[Mod - mMPI] Finishing subroutine start_mpi: rank=",rank
#if (DEBUG)
			call mpi_barrier(mpi_comm_world,ierr)
			if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&start_mpi at call to mpi_barrier in DEBUG')
#endif
#endif
	end subroutine start_mpi

	subroutine get_time(time)
		implicit none
		real(8) :: time
		time_end = mpi_wtime(ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&get_time at call to mpi_wtime')
		time = time_end-time_start
	end subroutine get_time

        subroutine one_write(str)
                implicit none
                character(len=*), intent(in) :: str
                if (rank.eq.0) write(*,*) str
        end subroutine one_write

	subroutine end_mpi
		implicit none
		time_end = mpi_wtime(ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&end_mpi at call to mpi_wtime')
#if VERBOSITY > 0
		if (rank==0) write(*,*) 'Total time taken=',time_end-time_start
#endif
		call mpi_finalize(ierr)
		if (ierr /= mpi_success) call mpi_error_stop('Error in subroutine &
			&end_mpi at call to mpi_finalize')
	end subroutine end_mpi

	subroutine mpi_error_stop(expl)
		implicit none
		character(len=*), intent(in) :: expl
		write(*,*) '[Mod - mMPI]'
		write(*,*) '-->'//expl
		call mpi_abort(mpi_comm_world, ierr, ierr)
	end subroutine mpi_error_stop

end module mMPI
