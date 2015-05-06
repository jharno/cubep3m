#include "./preprocessor"
module Zip
  use Parameters
  use Variables
  use mMPI
  implicit none
  
contains

  !Read in zip23 files and convert them to delta
  subroutine ngp_delta_zip23(grid, file2, file3)
    implicit none
    
    !Argument variables
    character(len=*), intent(in) :: file2, file3
    real, dimension(:,:,:), intent(out) :: grid

    !Zip 23 variables
    integer :: i,j,k,n,stat
    integer(1) :: i1(4)
    integer(4) :: i4

    !MPI variables
    real :: locmpi, glompi, mp, ntot
    real(8) :: lr8, gr8

    equivalence(i1,i4)

#if (VERBOSITY>0)
    if (rank==0) write(*,*) 'Entering subroutine ngp_delta_zip23'
#endif
    !Size of array                                                                                                                                                                           
    n=size(grid,dim=1)
    if (n/=Ncells) write(*,*) 'SEVERE ERROR IN module Zip in subroutine ngp_delta_zip23: n/=Ncells',n,Ncells

    !Open files                                                                                                                                                                              
    open(unit=11,file=trim(adjustl(file2)),status='old',iostat=stat,access='stream')
    if (stat/=0) then
       write(*,*) "Error in module Zip in subroutine ngp_delta_zip23 opening file "//trim(adjustl(file2))
       call mpi_abort(mpi_comm_world,ierr,ierr)
    end if
    open(unit=12,file=trim(adjustl(file3)),status='old',iostat=stat,access='stream')
    if (stat/=0) then
       write(*,*) "Error in module Zip in subroutine ngp_delta_zip23 opening file "//trim(adjustl(file3))
       call mpi_abort(mpi_comm_world,ierr,ierr)
    end if

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
    write(*,*) "Error with file "//trim(adjustl(file2))
    call mpi_abort(mpi_comm_world,ierr,ierr)
    96 close(11)

    read(12, end=97) i1
    write(*,*) "Error with file "//trim(adjustl(file3))
    call mpi_abort(mpi_comm_world,ierr,ierr)
    97 close(12)

    !Normalize density field
    !Diagnostics                                                                                                                                                                             
    locmpi = maxval(grid)
    call mpi_allreduce(locmpi,glompi,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    if(ierr/=mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
#if (VERBOSITY>0)
    if (rank.eq.0) write(*,*) '-->Max value of grid = ', glompi
#endif
    locmpi = minval(grid)
    call mpi_allreduce(locmpi,glompi,1,mpi_real,mpi_min,mpi_comm_world,ierr)
    if(ierr/=mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
#if (VERBOSITY>0)
    if (rank.eq.0) write(*,*) '-->Min value of cube = ', glompi
#endif

!    locmpi = sum(grid)
!    call mpi_allreduce(locmpi,glompi,1,mpi_real,mpi_sum,mpi_comm_world,ierr)
!    if (ierr/=mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)
!    ntot = glompi

    lr8 = sum(grid*1.d0)
    call mpi_allreduce(lr8,gr8,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    if (ierr/=mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)
    ntot = real(gr8,kind=4)

    !Compute mass                                                                                                                                                                            
    mp = (1.0*Ncells)**3*Nmpi/Ntot
#if (VERBOSITY>0)
    if (rank.eq.0) then
      write(*,*) '-->Total number of particles = ', ntot
      write(*,*) '-->Number of cells = ', (1.0*Ncells)**3*Nmpi
      write(*,*) '-->Particle mass = ', mp
    end if
#endif

    !Compute delta                                                                                                                                                                           
    !$omp workshare
      grid = mp*grid - 1.0
    !$omp end workshare

#if (VERBOSITY>0)
    if (rank==0) write(*,*) 'Finished subroutine ngp_delta_zip23'
#endif
  end subroutine ngp_delta_zip23

end module Zip
