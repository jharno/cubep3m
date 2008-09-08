subroutine densitystatistics
  implicit none
 
  include 'mpif.h'
  include 'dist_init.fh'

  integer, parameter :: kpt=nc/nt

  integer i,j,k
  real d,dmin,dmax
  real*8 dsum,dvar

  real(4) :: red4
  real(8) :: red8

  real time1,time2
  call cpu_time(time1)

  !! Inverse Fourier transform to get real space dm density field in cube
  call fftw(-1) 

  dmin=0
  dmax=0
  dsum=0
  dvar=0
  !$omp parallel do default(shared) private(i,j,k,d) &
  !$omp& reduction(min:dmin) reduction(max:dmax) reduction(+:dsum,dvar)
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
 
  !! Now we need to reduce these values to the master node

  if (rank == 0) write(*,*)
  call mpi_reduce(dmin,red4,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'DM min   ',red4
  call mpi_reduce(dmax,red4,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'DM max   ',red4
  call mpi_reduce(dsum,red8,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'Deltasum ',real(red8)
  call mpi_reduce(dvar,red8,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
  red8=sqrt(red8/nc**3)
  if (rank == 0) then
    write(*,*) 'Deltavar ',real(red8)
    write(*,*)
  endif

  !! Now for the gas

  !! copy init back into slab
  slab=init

  !! transform back to get real space gas density field in cube
  call fftw(-1)

  dmin=0
  dmax=0
  dsum=0
  dvar=0
  !$omp parallel do default(shared) private(i,j,k,d) &
  !$omp& reduction(min:dmin) reduction(max:dmax) reduction(+:dsum,dvar)
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

  !! Now we need to reduce these values to the master node

  if (rank == 0) write(*,*)
  call mpi_reduce(dmin,red4,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)  
  if (rank == 0) write(*,*) 'Gas min  ',red4
  call mpi_reduce(dmax,red4,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)  
  if (rank == 0) write(*,*) 'Gas max  ',red4
  call mpi_reduce(dsum,red8,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'Deltasum ',real(red8)
  call mpi_reduce(dvar,red8,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
  red8=sqrt(red8/nc**3)
  if (rank == 0) then
    write(*,*) 'Deltavar ',real(red8)
    write(*,*)
  endif

  call cpu_time(time2)
  time2=(time2-time1)/nt
  if (rank == 0) write(*,"(f8.2,a)") time2,'  Called density statistics'

end subroutine densitystatistics
