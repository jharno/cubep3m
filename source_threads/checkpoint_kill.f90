! write checkpoints to disk
  subroutine checkpoint_kill

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
#ifdef NEUTRINOS
    integer(4) :: np_dm, np_nu, ind_check1, ind_check2
    character (len=max_path) :: ofile_nu, ofile2_nu
#endif

!! label files with the same z as in the checkpoints file

    if (rank == 0) z_write = 1./a - 1. 
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

    ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xvres'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#ifndef NEUTRINOS
#ifdef PID_FLAG
    ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PIDres'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#endif
#endif
#ifdef NEUTRINOS
    ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xvres'// &
         rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif

!! Open checkpoint
    open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! Increment checkpoint counter so restart works on next checkpoint

!!    cur_checkpoint=cur_checkpoint+1

!! This is the file header

#ifdef NEUTRINOS

    !! Open neutrino checkpoint file
    open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile_nu
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    !! Determine how many dark matter and neutrino particles this rank has
    np_dm = 0
    np_nu = 0

    do i = 1, np_local
        if (PID(i) == 1) then
            np_dm = np_dm + 1
        else
            np_nu = np_nu + 1
        endif
    enddo

    if (rank == 0) write(*,*) "checkpoint np_dm, np_nu = ", np_dm, np_nu


    write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p
    write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p

#else

    write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p

#endif

#ifdef NEUTRINOS

    ind_check1 = 0
    ind_check2 = 0

    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
      do j=nplow,nphigh
        if (PID(j) == 1) then
#ifdef DISP_MESH
            write(12) xv(1:3,j) - shake_offset
            write(12) xv(4:6,j)
#else
            write(12) xv(:,j)
#endif
            ind_check1 = ind_check1 + 1
        else
#ifdef DISP_MESH
            write(22) xv(1:3,j) - shake_offset
            write(22) xv(4:6,j)
#else
            write(22) xv(:,j)
#endif
            ind_check2 = ind_check2 + 1
        endif
      enddo
    enddo

    close(12)
    close(22)

    !! Consistency check
    if (ind_check1 .ne. np_dm .or. ind_check2 .ne. np_nu) then
        write(*,*) "Dark Matter checkpoint error: ind_checks ", ind_check1, np_dm, ind_check2, np_nu
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#else

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
#ifdef DISP_MESH
        write(12) xv(1:3,j) - shake_offset
        write(12) xv(4:6,j)
#else
        write(12) xv(:,j)
#endif
      enddo
    enddo

    close(12)

#endif

!! Open PID file
#ifndef NEUTRINOS
#ifdef PID_FLAG

    open(unit=15, file=ofile2, status="replace", iostat=fstat, access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening PID file for write'
      write(*,*) 'rank',rank,'file:',ofile2
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! This is the file header

    write(15) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p

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
#endif

#ifdef MHD
!! Write gas checkpoint
    call mpi_tvd_mhd_state_output(output_path,nts,t,z_s)
#endif

    write(*,*) 'Finished checkpoint_kill:',rank

    !checkpoint_step=.false.

  end subroutine checkpoint_kill

subroutine read_remaining_time
    !
    ! Subroutine that reads the first line of killtime_path which is to contain
    ! the amount of time (in seconds) remaining in the job. Will checkpoint with
    ! kill_remaining (parameter in cubepm.par) seconds remaining. 
    !

    implicit none

    integer(4) :: time_lefti
    real(4) :: time_left
    real(4) :: time_left_default = 48.*3600. !! Maximum walltime of GPC is taken if problems exist

    include 'mpif.h'
#ifdef PPINT
    include 'cubep3m.fh'
#else
    include 'cubepm.fh'
#endif

    if (rank == 0) then !! Only master rank will read this file

        open(31, file=killtime_path, status="old", iostat=ierr)
        if (ierr /= 0) then !! If file does not exist then assume default time remaining 
            time_left = time_left_default
            write(*,*) "WARNING: Could not open file ", killtime_path
            write(*,*) "Assuming the time left (in hours) is ", time_left/3600.
        else !! Read the file
            read(unit=31,fmt='(I20)') time_lefti
            close(31)
            time_left = 1.*time_lefti
        endif 

        !! Kill the job in this amount of time from now.
        kill_time = time_left - kill_remaining
        write(*,*) "Killing job ", kill_time/3600., " hours from now"

        !! Consistency check
        if (kill_time <= 0.) then
            kill_time = time_left_default
            write(*,*) "WARNING: kill_time did not make sense. Taking kill_time (in hours): ", kill_time/3600. 
            write(*,*) "You should check these: ", time_left, time_left_default, kill_remaining
        endif
    
    endif

    !! Broadcast to everyone
    call mpi_bcast(kill_time, 1, mpi_real, 0, mpi_comm_world, ierr)

end subroutine

