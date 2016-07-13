module FieldIO
  use Parameters
  use Variables
  use mMPI
  implicit none

contains

  subroutine read_field3(grid,fileI)
    implicit none
    real, dimension(:,:,:), intent(out) :: grid
    character(len=*), intent(in) :: fileI
    integer :: stat
    !if (rank==0) write(*,*) 'Entering subroutine read_field3'
    open(unit=21,file=trim(adjustl(fileI)),status='old',iostat=stat,access='stream')
    if (stat/=0) then
      write(*,*) 'ERROR in module FieldIO in subroutine read_filed3 opening file: '//trim(adjustl(fileI))
      call mpi_abort(mpi_comm_world,ierr,ierr)
    end if
    read(21) grid
    close(21)
    !if (rank==0) write(*,*) 'Finished subroutine read_field3'
  end subroutine read_field3

  subroutine write_field3(grid,fileO)
    implicit none
    character(len=*), intent(in) :: fileO
    real, dimension(:,:,:), intent(in) :: grid
    integer :: stat
    !if (rank==0) write(*,*) 'Entering subroutine write_den'
    open(unit=21,file=trim(adjustl(fileO)),status='replace',iostat=stat,access='stream')
    if (stat/=0) then
      write(*,*) 'ERROR in FieldIO in subroutine write_den opening file: '//trim(adjustl(fileO))
      call mpi_abort(mpi_comm_world,ierr,ierr)
    end if
    write(21) grid
    close(21)
    !if (rank==0) write(*,*) 'Finished subroutine write_den'
  end subroutine write_field3

end module FieldIO
