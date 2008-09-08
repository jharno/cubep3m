implicit none
include 'mpif.h'
character(len=7) rank_s
character(len=80) ofile
integer :: fstat,rank,ierr,nn,np,npt
call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
call mpi_comm_size(mpi_comm_world,nn,ierr)
      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)

      ofile='/scratch/merz/cubepm/xv'//rank_s(1:len_trim(rank_s))//'.ic'
!      print *,'opening particle list:',ofile(1:len_trim(ofile))
      open(unit=20,file=ofile,form='binary',iostat=fstat,status='old')
      if (fstat /= 0) then
        write(*,*) 'error opening initial conditions'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(20) np
close(20)
call mpi_reduce(np,npt,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
if (rank==0) print *,npt
call mpi_finalize(ierr)
end    
