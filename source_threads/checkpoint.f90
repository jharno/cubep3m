! write checkpoints to disk
  subroutine checkpoint

#ifdef MHD
    use mpi_tvd_mhd
#endif

    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    character (len=max_path) :: ofile,ofile2
    character (len=4) :: rank_s
    character (len=7) :: z_s  

    integer(kind=4) :: i,j,fstat,blocksize,num_writes,nplow,nphigh
    integer(kind=4) :: cur_proj,cur_halo
    real(kind=4) :: z_write
#ifdef NEUTRINOS
    integer(4) :: np_dm, np_nu, ind_check1, ind_check2
    character (len=max_path) :: ofile_nu, ofile2_nu
#endif

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

    ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#ifdef NEUTRINOS
    ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
         rank_s(1:len_trim(rank_s))//'_nu.dat'
#ifdef NUPID
    ofile2_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PID'// &
         rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif
#else
#ifdef PID_FLAG
    ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PID'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#endif
#endif

!! Open checkpoint
    open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! Increment checkpoint counter so restart works on next checkpoint

    cur_checkpoint=cur_checkpoint+1
    cur_proj = cur_projection
    cur_halo = cur_halofind
    if (projection_step) cur_proj = cur_proj + 1
    if (halofind_step) cur_halo = cur_halo + 1

!! This is the file header

#ifdef NEUTRINOS

    !! Open neutrino checkpoint file
    open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream") 
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile_nu
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#ifdef NUPID
    !! Open neutrino PID file
    open(unit=23, file=ofile2_nu, status="replace", iostat=fstat, access="stream")
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile2_nu
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#endif

    !! Determine how many dark matter and neutrino particles this rank has
    np_dm = 0
    np_nu = 0

    do i = 1, np_local
#ifdef NUPID
        if (PID(i) == 0) then
#else
        if (PID(i) == 1) then
#endif
            np_dm = np_dm + 1
        else
            np_nu = np_nu + 1
        endif
    enddo

    if (rank == 0) write(*,*) "checkpoint np_dm, np_nu = ", np_dm, np_nu

    write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_proj,cur_halo,mass_p
    write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_proj,cur_halo,mass_p
#ifdef NUPID
    write(23) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_proj,cur_halo,mass_p
#endif

#else

    write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_proj,cur_halo,mass_p

#endif


#ifndef SLOWXV

#  ifdef DISP_MESH                                                                                                           
   xv(1:3,:) = xv(1:3,:) + spread(shake_offset,dim=2,ncopies=np_local) 
#  endif

#  ifdef NEUTRINOS
  ind_check1 = 0
  ind_check2 = 0
do i=1,num_writes
  nplow=(i-1)*blocksize+1
  nphigh=min(i*blocksize,np_local)
  write(12) pack(xv(:,nplow:nphigh),spread(PID(nplow:nphigh)==1,dim=1,ncopies=6))
  ind_check1 = ind_check1 + count(PID(nplow:nphigh)==1)
  write(22) pack(xv(:,nplow:nphigh),spread(PID(nplow:nphigh)/=1,dim=1,ncopies=6))
  ind_check2 = ind_check2 + count(PID(nplow:nphigh)/=1)
#  ifdef NUPID
  write(23) PID(j)
#  endif
enddo
close(12); close(22)
#  ifdef NUPID
close(23)
#  endif 
if (ind_check1 .ne. np_dm .or. ind_check2 .ne. np_nu) then
    write(*,*) "Dark Matter checkpoint error: ind_checks ", ind_check1, np_dm, ind_check2, np_nu
    call mpi_abort(mpi_comm_world,ierr,ierr)
endif

#  else
! NEUTRINOS

do i=1,num_writes
  nplow=(i-1)*blocksize+1
  nphigh=min(i*blocksize,np_local)
  write(12) xv(:,nplow:nphigh)
enddo
close(12)
#  endif
! NEUTRINOS

#  ifdef DISP_MESH
       xv(1:3,:) = xv(1:3,:) - spread(shake_offset,dim=2,ncopies=np_local)
#  endif

#else
! SLOWXV
#ifdef NEUTRINOS

    ind_check1 = 0
    ind_check2 = 0
    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
      do j=nplow,nphigh

#ifdef NUPID   !!! DARK MATTER
        if (PID(j) == 0) then
#else
        if (PID(j) == 1) then
#endif

#ifdef DISP_MESH
            write(12) xv(1:3,j) - shake_offset
            write(12) xv(4:6,j)
#else
            write(12) xv(:,j)
#endif
            ind_check1 = ind_check1 + 1
        else  !!! NEUTRINOS
#ifdef DISP_MESH
            write(22) xv(1:3,j) - shake_offset
            write(22) xv(4:6,j)
#else
            write(22) xv(:,j)
#endif
#  ifdef NUPID
            write(23) PID(j)
#  endif
            ind_check2 = ind_check2 + 1
        endif
      enddo
    enddo

    close(12)
    close(22)
#  ifdef NUPID
    close(23)
#  endif

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
# ifdef DISP_MESH
        write(12) xv(1:3,j) - shake_offset
        write(12) xv(4:6,j)
# else
        write(12) xv(:,j)
# endif
      enddo
    enddo

    close(12)

#endif
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
              cur_proj,cur_halo,mass_p

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

if (rank==0) then
   print*, 'current steps recorded in xv file:'
   print*, 'cur_checkpoint =', cur_checkpoint
   print*, 'cur_projection =', cur_proj
   print*, 'cur_halofind   =', cur_halo
endif

    write(*,*) 'Finished checkpoint:',rank

    checkpoint_step=.false.

  end subroutine checkpoint
