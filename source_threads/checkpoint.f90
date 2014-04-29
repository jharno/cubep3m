! write checkpoints to disk
  subroutine checkpoint

#ifdef MHD
    use mpi_tvd_mhd
#endif

    implicit none

    include 'mpif.h'
#ifdef PPINT
    include 'cubep3m.fh'
#else
    include 'cubepm.fh'
#endif

    character (len=max_path) :: ofile,ofile2
    character (len=4) :: rank_s
    character (len=7) :: z_s  

    integer(kind=4) :: i,j,fstat,blocksize,num_writes,nplow,nphigh
    real(kind=4) :: z_write

!! label files with the same z as in the checkpoints file

    if (rank == 0) z_write=z_checkpoint(cur_checkpoint)
    call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

!! most linux systems choke when writing more than 2GB of data
!! in one write statement, so break up into blocks < 2GB 

!    blocksize=(2047*1024*1024)/24
! reduced to 32MB chunks because of intel compiler
    blocksize=(32*1024*1024)/24
    num_writes=np_local/blocksize+1

!! Create checkpoint file name

    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)

    write(z_s,'(f7.3)') z_write
    z_s=adjustl(z_s)

    ofile=output_path//z_s(1:len_trim(z_s))//'xv'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#ifdef PID_FLAG
    ofile2=output_path//z_s(1:len_trim(z_s))//'PID'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#endif

!! Open checkpoint

#ifdef BINARY
    open (unit=12,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=12,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! Increment checkpoint counter so restart works on next checkpoint

    cur_checkpoint=cur_checkpoint+1

!! This is the file header

#ifdef PPINT
    write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p
#else
    write(12) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p
#endif
!! Particle list
!    do i=0,nodes-1
!      if (rank==i) print *,rank,num_writes,np_local
!      call mpi_barrier(mpi_comm_world,ierr)
!    enddo
! the intel compiler puts arrays on the stack to write to disk
! this causes memory problems
!    write(12) xv(:,:np_local)
    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
!!      print *,rank,nplow,nphigh,np_local
      do j=nplow,nphigh
        !write(12) xv(:,j)
        write(12) xv(1:3,j) - shake_offset, xv(4:6,j)
      enddo
    enddo

    close(12)

!! Open PID file
#ifdef PID_FLAG

#ifdef BINARY
    open (unit=15,file=ofile2,status='replace',iostat=fstat,form='binary')
#else
    open (unit=15,file=ofile2,status='replace',iostat=fstat,form='unformatted')
#endif

    if (fstat /= 0) then
      write(*,*) 'error opening PID file for write'
      write(*,*) 'rank',rank,'file:',ofile2
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! This is the file header

#ifdef PPINT
    write(15) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p
#else
    write(15) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p
#endif
!! Particle list
!    do i=0,nodes-1
!      if (rank==i) print *,rank,num_writes,np_local
!      call mpi_barrier(mpi_comm_world,ierr)
!    enddo
! the intel compiler puts arrays on the stack to write to disk
! this causes memory problems
!    write(12) xv(:,:np_local)
    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
!!      print *,rank,nplow,nphigh,np_local
      do j=nplow,nphigh
        write(15) PID(j)
      enddo
    enddo

    close(15)

#endif
    

#ifdef MHD
!! Write gas checkpoint
    call mpi_tvd_mhd_state_output(output_path,nts,t,z_s)
#endif

    write(*,*) 'Finished checkpoint:',rank

    checkpoint_step=.false.

  end subroutine checkpoint
