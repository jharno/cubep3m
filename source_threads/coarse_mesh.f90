!! coarse mesh velocity update
  subroutine coarse_mesh
   implicit none
   include 'mpif.h'
#   include "cubepm.fh"

#ifdef MHD
   integer :: nerr,nerrl
   real :: cmax, cmaxl
#endif

#ifdef DIAG
   integer(4) :: i,j,k
   real(8) :: sumrhoc,sumrhof
#endif

#ifdef DEBUG_RHOC
   integer(4) :: p_in_cell,pp
#endif

#ifdef DEBUG
   if (rank == 0) write(*,*) 'starting coarse mesh'
#endif

   call coarse_mass

#ifdef DIAG
    sumrhof=0.0
    sumrhoc=0.0
    do k=1,nc_node_dim
      do j=1,nc_node_dim
        do i=1,nc_node_dim
          sumrhof=sumrhof+real(rho_c(i,j,k),kind=8)
        enddo
      enddo
    enddo
    call mpi_reduce(sumrhof,sumrhoc,1,mpi_double_precision, &
                    mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'sum of rho_c=',sumrhoc
    if (rank == 0) write(*,*) 'rank', rank, 'min/max rho_c =', minval(rho_c), maxval(rho_c)
#endif

#ifdef DEBUG
   if (rank == 0) write(*,*) 'finished coarse mass'
#endif

#ifdef DEBUG_RHOC
   do k=1,nc_node_dim
     do j=1,nc_node_dim
       do i=1,nc_node_dim
         if (rho_c(i,j,k) /= 0.0) write(*,*) rank,'rho_c',&
           i+nc_node_dim*cart_coords(3),j+nc_node_dim*cart_coords(2), &
           k+nc_node_dim*cart_coords(1), rho_c(i,j,k)
       enddo
     enddo
   enddo
!   write(*,*) 'finished density check'
   do k=0,nc_node_dim+1
     do j=0,nc_node_dim+1
       do i=0,nc_node_dim+1
         p_in_cell=0
         pp=hoc(i,j,k)
         do
           if(pp==0) exit
           p_in_cell=p_in_cell+1
           pp=ll(pp)
         enddo
         if (p_in_cell /= 0) write(*,*) rank,'p_in_cell',&
           i+nc_node_dim*cart_coords(3),j+nc_node_dim*cart_coords(2), &
           k+nc_node_dim*cart_coords(1), p_in_cell 
       enddo
     enddo
   enddo
#endif

   call coarse_force

#ifdef DEBUG
   if (rank == 0) print *, 'finished coarse force'
#endif

   call coarse_force_buffer

#ifdef DEBUG
   if (rank == 0) print *, 'finished coarse force buffer'
#endif

   call coarse_max_dt

#ifdef DEBUG
   if (rank == 0) print *, 'finished coarse max dt'
#endif

#ifdef MHD
   call coarse_velocity(cmax,nerr)
   cmaxl=cmax
   nerrl=nerr
   call mpi_reduce(cmaxl,cmax,1,mpi_real,mpi_max,0,mpi_comm_cart,ierr)
   call mpi_reduce(nerrl,nerr,1,mpi_integer,mpi_sum,0,mpi_comm_cart,ierr)
   if (rank==0) print *,'after coarse velocity fluid stats',cmax/freeze,dt*cmax,nerr
#else
   if (coarse_vel_update) call coarse_velocity
#endif

#ifdef DEBUG
   if (rank == 0) print *, 'finished coarse velocity'
#endif

  end subroutine coarse_mesh
