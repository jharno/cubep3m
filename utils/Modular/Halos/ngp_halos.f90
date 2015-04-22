#include './preprocessor'
program ngp_halos
  use Parameters
  use Variables
  use mMPI
  use HaloReader
  use FieldIO
  implicit none

  real, dimension(Ncells, Ncells, Ncells) :: grid

  real, dimension(:,:), allocatable :: xvh
  real, dimension(:), allocatable :: mhg,mhm
  character(len=*), parameter :: jdir = '/scratch/p/pen/emberson/cubep3m/timing/cubep3m_movie_th2' 
  character(len=1000) :: fn

  integer, parameter :: group = 2 !Number between 1 and numgp
  integer, parameter :: numgp = 2
  character(len=2) :: group_str
  real, parameter :: maxCDF = (1.0*group)/(numgp)
  real, parameter :: minCDF = (1.0*(group-1.0))/(numgp)

  integer, parameter :: ncdf = 1000
  real, dimension(2,ncdf) :: cdf
  character(len=*), parameter :: fcdf = '/scratch/p/pen/dinman/test_th2/halo_cdf/cdf_dmnu_gpc.dat'
  
  integer :: h,nhl,nhg,nhm,nhr,nind,stat

  integer :: group_halos = 0
  integer :: non_group_halos = 0
  integer :: skipped_halos = 0
  integer :: nx, ny, nz

  integer :: locl,glob
  real :: mp

  call start_mpi

  !Open halo file
  fn = jdir//'/node'//trim(adjustl(rank_s))//'/0.000halo'//trim(adjustl(rank_s))//'.dat'
  open(unit=11,file=fn,status='old',iostat=stat,access='stream')
  read(11) nhl
  close (11)

  !Compute global number of halos
  call mpi_allreduce(nhl,nhg,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  if (rank==0) write(*,*) 'Total # of halos: ',nhg

  !Compute maximum number of halos
  call mpi_allreduce(nhl,nhm,1,mpi_integer,mpi_max,mpi_comm_world,ierr)
  if (rank==0) write(*,*) 'Max # of halos per node: ',nhm

  !Allocate array
  allocate(xvh(35,nhl))
  
  !Read halos
  call read_halo_file(trim(adjustl(fn)),xvh,nhl)

  !Check halos correct size
  locl = int(minval(xvh(4,:)/64.0))
  glob=0
  call mpi_allreduce(locl,glob,1,mpi_integer,mpi_min,mpi_comm_world,ierr)
  if (rank==0) write(*,*) 'Min halo mass', glob, minval(xvh(4,:)/64)

  !Get total particles in halos
  locl = int(sum(xvh(4,:)/64.0))
  call mpi_allreduce(locl,glob,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  if (rank==0) write(*,*) 'Total number of dm particles in halos (virial): ',glob
  locl = int(sum(xvh(5,:)/64.0))
  call mpi_allreduce(locl,glob,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  if (rank==0) write(*,*) 'Total number of dm particles in halos (odc): ',glob
  if (rank==0) write(*,*) 'Total number of dm particles: ',(nc)**3

  locl = int(sum(xvh(35,:)))
  call mpi_allreduce(locl,glob,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  if (rank==0) write(*,*) 'Total number of nu particles near halos: ',glob
  if (rank==0) write(*,*) 'Total number of nu particles: ',(2*nc)**3

!  !Convert halo coordinates to local coordinates and masses to log10 masses
!  !$omp parallel do default(none) shared(nhl,xvh) private(h)
  do h=1,nhl
!     xvh(1:3,h) = local_coordinates(xvh(1:3,h)/4.0)
     xvh(1:3,h) = xvh(1:3,h) - slab_coord(:)*nc_node_dim*4
     xvh(1:3,h) = xvh(1:3,h)/4.0
     xvh(8:10,h) = local_coordinates(xvh(8:10,h)/4.0)
     xvh(29:31,h) = local_coordinates(xvh(29:31,h)/4.0)
     xvh(4,h) = log10(xvh(4,h))
     xvh(5,h) = log10(xvh(5,h))
  end do
!  !$omp end parallel do

  !Read halo CDF files
  open(unit=11,file=fcdf,iostat=stat,status='old')
  do h=1,ncdf
     read(11,*) cdf(:,h)
  end do
  close(11)

  !Initialize grid
  grid = 0.0

  if (rank==0) write(*,*) 'Assigning groups in:',minCDF,maxCDF

  !Add halos in group to grid
  do h=1, nhl

     nx = 1+floor(xvh(1,h))
     ny = 1+floor(xvh(2,h))
     nz = 1+floor(xvh(3,h))

     !Only use halos on this node

     if ( nx .lt. 1 .or. nx .gt. Ncells ) cycle
     if ( ny .lt. 1 .or. ny .gt. Ncells ) cycle
     if ( nz .lt. 1 .or. nz .gt. Ncells ) cycle     

     !Check to see if in group
     nind = nindex(xvh(5,h),cdf(1,:))

     if (cdf(2,nind) .gt. minCDF .and. cdf(2,nind) .le. maxCDF) then
        !In group
        grid(nx,ny,nz) = grid(nx,ny,nz) + 1.0
        group_halos = group_halos+1
     else 
        !Not in group
        non_group_halos = non_group_halos+1
     end if

  end do
  
  !Some statistics
  skipped_halos = nhg
  call mpi_allreduce(group_halos,glob,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  if (rank==0) write(*,*) "# of halos in this group: ", glob
  skipped_halos = skipped_halos-glob
  call mpi_allreduce(non_group_halos,glob,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  if (rank==0) write(*,*) "# of halos not in this group: ", glob
  skipped_halos = skipped_halos-glob
  if (rank==0) write(*,*) "# of halos not on correct node: ", skipped_halos
  if (rank==0) write(*,*) 'Expected # of halos in group: ', 1.0*(nhg-skipped_halos)*(maxCDF-minCDF)
  
  !Compute density contrast
  nhg = sum(grid)
  call mpi_allreduce(nhg,glob,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
  mp = (1.0*Ncells)**3*Nmpi/glob

  !$omp workshare
  grid = mp*grid - 1.0
  !$omp end workshare

  !Write out delta
  write(group_str,'(I2)') group
  fn = dir//'fields/ha/den/0.000halo'//trim(adjustl(rank_s))//'_g'//trim(adjustl(group_str))//'.dat'
  call write_field3(grid,fn)

  deallocate(xvh)
  call end_mpi

contains
  
  !Function to find which index in the cdf table a given mass corresponds too
  function nindex(m,col) result(index)
    implicit none
    real, intent(in) :: m
    real, dimension(:) :: col
    integer :: i, index
    real :: dm

    do i=1,size(col)
       dm = m - col(i)
       if (dm .lt. 0) then
          index = i
          if ( i.gt. 1 ) then
             if ( m - col(i-1) .lt. abs(dm) ) index = i-1
          end if
          return
       end if
    end do
    index = size(col)

  end function nindex
  
end program ngp_halos
