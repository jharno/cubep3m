#include "./preprocessor"

!Program to determine ngp_density from zip2 and zip3 files
!Output: den and den_nu files containing the density contrast

program ngp_density
  use Parameters
  use Variables
  use mMPI
  use FieldIO

  implicit none

  character(len=100) :: fn2, fn3
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube

  call start_mpi
  if (rank==0) write(*,*) 'Starting program ngp_density'

  !DM density
  if (rank==0) write(*,*) 'Computing dark matter density contrast'
  !!Read
  !$omp workshare
    cube=0.0
  !$omp end workshare
  fn2 = dir//'dmnu/dm/zip2/'//redshift//'zip2_'//trim(adjustl(rank_s))//'.dat'
  fn3 = dir//'dmnu/dm/zip3/'//redshift//'zip3_'//trim(adjustl(rank_s))//'.dat'
  call read_zip23(cube,trim(adjustl(fn2)),trim(adjustl(fn3)))

  !!Normalize
  call normalize_den(cube)

  !!Write
  fn2 = dir//'fields/dm/den/'//redshift//'den'//trim(adjustl(rank_s))//'.dat'
  call write_field3(cube,trim(adjustl(fn2)))

  !NU density
  if (rank==0) write(*,*) 'Computing dark matter density contrast'
  !!Read
  !$omp workshare
    cube=0.0
  !$omp end workshare
  fn2 = dir//'dmnu/nu/zip2/'//redshift//'zip2_'//trim(adjustl(rank_s))//'_nu.dat'
  fn3 = dir//'dmnu/nu/zip3/'//redshift//'zip3_'//trim(adjustl(rank_s))//'_nu.dat'
  call read_zip23(cube,trim(adjustl(fn2)),trim(adjustl(fn3)))

  !!Normalize
  call normalize_den(cube)

  !!Write
  fn2 = dir//'fields/nu/den/'//redshift//'den'//trim(adjustl(rank_s))//'_nu.dat'
  call write_field3(cube,trim(adjustl(fn2)))

  if (rank==0) write(*,*) 'Finished program ngp_density'
  call end_mpi

contains

  subroutine read_zip23(grid,fileA,fileB)
    implicit none
    character(len=*), intent(in) :: fileA,fileB
    real, dimension(:,:,:), intent(out) :: grid
    integer :: i,j,k,n,stat
    integer(1) :: i1(4)
    integer(4) :: i4

    equivalence(i1,i4)

    if (rank==0) write(*,*) 'Entering subroutine read_zip23'

    !Size of array
    n=size(grid,dim=1)
    if (n/=Ncells) write(*,*) 'SEVERE ERROR IN read_zip23: n/=Ncells',n,Ncells

    !Open files
    open(unit=11,file=trim(adjustl(fileA)),status='old',iostat=stat,access='stream')
    open(unit=12,file=trim(adjustl(fileB)),status='old',iostat=stat,access='stream')

    do k=1,n
      do j=1,n
        do i=1,n
          i1=0
          i4=0
#ifdef BGQ
          read(11) i1(4)
#else
          read(11) i1(1)
#endif
          if(i4.eq.255) read(12) i4
          grid(i,j,k) = real(i4)
        end do
      end do
    end do

    !Check to see that you finished reading zip2 and zip3 files
    read(11, end=96) i1
    write(*,*) "Error with file "//trim(adjustl(fileA))
    call mpi_abort(mpi_comm_world,ierr,ierr)
    96 close(11)

    read(12, end=97) i1
    write(*,*) "Error with file "//trim(adjustl(fileB))
    call mpi_abort(mpi_comm_world,ierr,ierr)
    97 close(12)

    if (rank==0) write(*,*) 'Finished subroutine read_zip23'

  end subroutine read_zip23

  subroutine normalize_den(grid)
    implicit none
    real, dimension(:,:,:), intent(inout) :: grid
    real :: locmpi, glompi, mp, ntot

    if (rank==0) write(*,*) 'Entering subroutine normalize_den'

    !Diagnostics
    locmpi = maxval(grid)
    call mpi_allreduce(locmpi,glompi,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    if(ierr/=mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (rank.eq.0) write(*,*) '-->Max value of grid = ', glompi

    locmpi = minval(grid)
    call mpi_allreduce(locmpi,glompi,1,mpi_real,mpi_min,mpi_comm_world,ierr)
    if(ierr/=mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (rank.eq.0) write(*,*) '-->Min value of cube = ', glompi

    locmpi = sum(grid)
    call mpi_allreduce(locmpi,glompi,1,mpi_real,mpi_sum,mpi_comm_world,ierr)
    if (ierr/=mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)
    ntot = glompi

    !Compute mass
    mp = (1.0*Ncells)**3*Nmpi/Ntot
    if (rank.eq.0) then
      write(*,*) '-->Total number of particles = ', ntot
      write(*,*) '-->Number of cells = ', (1.0*Ncells)**3*Nmpi
      write(*,*) '-->Particle mass = ', mp
    end if

    !Compute delta
    !$omp workshare
      grid = mp*grid - 1.0
    !$omp end workshare

    if (rank==0) write(*,*) 'Finished subroutine normalize_den'

  end subroutine normalize_den

end program ngp_density
