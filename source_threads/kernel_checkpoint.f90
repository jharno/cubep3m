! write kernel to disk
  subroutine kernel_checkpoint(mode)
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    logical :: mode
    character (len=max_path) :: ofile
    character (len=7) :: z_s  
    integer(kind=4) :: i,j,k,tag,fstat
    real, dimension(:,:,:,:), allocatable :: full_kern
    integer, dimension(mpi_status_size) :: status

!! Create checkpoint file name

    if (rank==0) then

      ofile=output_path//'coarse_kernel.dat'

!! Open checkpoint

      open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening coarse kernel file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      allocate(full_kern(3,nc_dim/2+1,nc_dim,nc_dim))

    endif

!real(4), dimension(3,nc_dim/2+1,nc_dim,nc_slab)       :: kern_c

    if (mode) then
      do i=1,nodes-1
        tag=i
        if (rank==0) then
          call mpi_recv(full_kern(1,1,1,i*nc_slab+1),3*(nc_dim/2+1)*nc_dim*nc_slab,mpi_real,i,tag,mpi_comm_world,status,ierr)
        elseif (rank == i) then
          call mpi_send(kern_c,3*(nc_dim/2+1)*nc_dim*nc_slab,mpi_real,0,tag,mpi_comm_world,ierr)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        if (rank==0)  print *,'finished rank',i
      enddo
      if (rank==0) then
        write(12) full_kern
        close(12)
        deallocate(full_kern)
      endif
    else
      read(12) full_kern 
      close(12)

      !! comparison code 

    endif

    write(*,*) 'Finished coarse kernel checkpoint:',rank

    stop

  end subroutine kernel_checkpoint

