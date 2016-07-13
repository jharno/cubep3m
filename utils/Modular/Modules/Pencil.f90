module Pencil
	use Parameters
	use Variables
	use mMPI
	use p3dfft
	implicit none
	public

	integer :: plan, iplan
	logical :: firstfftw = .true.
	logical :: fourier_space =.false.

	integer, parameter :: nodes_pen = nodes_dim
	integer, parameter :: nc_pen = nc_node_dim / nodes_dim
	integer, parameter :: dim_y = nodes_dim
	integer, parameter :: dim_z = nodes_dim**2
	integer :: pen_dims(2), istart(3), iend(3), isize(3), fstart(3), fend(3), fsize(3), mypadd
	integer, dimension(0:nodes_dim-1) :: pen_neighbor_to
	integer, dimension(0:nodes_dim-1) :: pen_neighbor_fm


	real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
	real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1) :: recv_cube

	real, dimension(nc, nc_node_dim, nc_pen+2) :: slab,slab2,slab3
	equivalence(cube, slab)

	public :: nc_pen
	public :: cube
	public :: slab!, slab2
	public :: setup_pencil
	public :: cp_fftw

contains

	subroutine setup_pencil
		implicit none
		integer :: j

		if (mod(nc,nodes_dim**2) /= 0) call pencil_error_stop('Error in subroutine &
			&setup_pencil: mesh does not evenly decompose')

		fourier_space = .false.

		do j=0, nodes_dim -1
			pen_neighbor_to(j) = nodes_slab*slab_coord(3) + slab_coord(2) + j*nodes_dim
			pen_neighbor_fm(j) = nodes_slab*slab_coord(3) + j + nodes_dim*slab_coord(1)
		end do

		!$omp workshare
			cube=0.0
		!$omp end workshare

		if (firstfftw) then
			pen_dims = (/dim_y, dim_z/)
			call p3dfft_setup(pen_dims,nc,nc,nc,.true.)
			call p3dfft_get_dims(istart,iend,isize,1,mypadd)
			call p3dfft_get_dims(fstart,fend,fsize,2)
			firstfftw = .false.
		end if

	end subroutine setup_pencil

	subroutine reset_pencil
		implicit none
		fourier_space=.false.
	end subroutine reset_pencil

 subroutine force_pencil
   implicit none
   fourier_space=.true.
 end subroutine force_pencil

	subroutine pack_pencils
		!Copied from JDs code, blame him

	    ! Pack cubic data into pencils for p3dfft transform.
	    implicit none

	    integer(4) :: i,j,k,i0,i1,k1
	    integer(4) :: pen_slice,tag,rtag
	    integer(4) :: num_elements
	    !possibly should be double so as not to overflow for large runs, but integer*8 might be
	    !unsupported by MPI

	    integer(4), dimension(2*nodes_dim) :: requests
	    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
	    integer(4) nc_pen_break, breakup
	    real(4) :: passGB

	    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
	    breakup = 1
	    num_elements = nc_node_dim * nc_node_dim * nc_pen
	    passGB = 4. * num_elements / 1024.**3
	    if (passGB > 1.) then
	        breakup = 2**ceiling(log(passGB)/log(2.))
	    end if
	    num_elements = num_elements / breakup

	    ! Send the data from cube to recv_cube
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
	        end do
	        call mpi_waitall(2*nodes_dim, requests, wait_status, ierr)
	    end do

	    ! Place this data into the pencils (stored in the slab array)
	    !$omp parallel do default(none) shared(slab, recv_cube) private(i,j,k,i0,i1,pen_slice)
	    do i = 0, nodes_dim - 1
	        i0 = i * nc_node_dim + 1
	        i1 = (i + 1) * nc_node_dim
	        pen_slice = i
	        do k = 1, nc_pen
	            do j = 1, nc_node_dim
	                slab(i0:i1,j,k) = recv_cube(:,j,k,pen_slice)
	            end do
	        end do
	    end do
	    !$omp end parallel do

	end subroutine pack_pencils

	subroutine unpack_pencils
	    !Copied from JDs code, blame him

	    ! Unpack data from the pencils back into the cubic decompisition following
	    ! p3dfft transform.
	    implicit none

	    integer(4) :: i,j,k,i0,i1,k1
	    integer(4) :: pen_slice,num_elements,tag,rtag
	    integer(4), dimension(2*nodes_dim) :: requests
	    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
	    integer(4) nc_pen_break, breakup
	    real(4) :: passGB

	    ! Place data in the recv_cube buffer
	    !$omp parallel do default(none) shared(recv_cube, slab) private(i,j,k,i0,i1,pen_slice)
	    do i = 0, nodes_dim - 1
	        i0 = i * nc_node_dim + 1
	        i1 = (i + 1) * nc_node_dim
	        pen_slice = i
	        do k = 1, nc_pen
	            do j = 1, nc_node_dim
	                recv_cube(:, j, k, pen_slice) = slab(i0:i1, j, k)
	            end do
	        end do
	    end do

	    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
	    breakup = 1
	    num_elements = nc_node_dim * nc_node_dim * nc_pen
	    passGB = 4. * num_elements / 1024.**3
	    if (passGB > 1.) then
	        breakup = 2**ceiling(log(passGB)/log(2.))
	    end if
	    num_elements = num_elements / breakup

	    ! Put this data back into cube
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
	        end do
	        call mpi_waitall(2*nodes_dim,requests, wait_status, ierr)
	    end do

	end subroutine unpack_pencils

	subroutine cp_fftw(direction)
		implicit none
		integer :: i
		integer :: direction

		if (direction/=0) then
			if (direction>0) then
				if (fourier_space .eqv. .true.) call pencil_error_stop('Attempting Forward FFT &
				&while already in Fourier space')
				call pack_pencils
				call ftran_r2c(slab,slab,"fft")
				fourier_space = .true.
			else
				if (fourier_space .eqv. .false.) call pencil_error_stop('Attempting Backward FFT &
				&while already in real space')
				call btran_c2r(slab,slab,"tff")
				slab = slab/(real(nc)**3)
				call unpack_pencils
				fourier_space = .false.
			end if
		else
			call p3dfft_clean
		end if

	end subroutine cp_fftw

	subroutine pencil_error_stop(expl)
		implicit none
		character(len=*) :: expl
		write(*,*) '[Mod - Pencil]'
		write(*,*) '-->'//expl
		call mpi_abort(mpi_comm_world, ierr, ierr)
	end subroutine pencil_error_stop

end module Pencil
