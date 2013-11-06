!! initialize particle list
  subroutine particle_initialize
    use omp_lib
    implicit none

    include 'mpif.h'
    include 'cubepm.fh'

    real(4) :: rnum,z_write
    integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh
    integer*8 :: np_total,npl8
    character(len=max_path) :: ofile
    character(len=4) :: rank_s
    character(len=7) :: z_s
#ifdef CHECK_IP
    real(8) :: xva(6)
#endif

    fstat=0

    np_local=(nf_physical_node_dim/2)**3

!! set flag in cubepm.par to select desired initial conditions
 
    if (random_ic) then

      do i=1,np_local
        do j=1,3
          call random_number(rnum)
          if (rnum >= 0.9999) then
            rnum=0.995
          endif
          xv(j,i)=rnum*nf_physical_node_dim
        enddo
        xv(4:6,i)=0.0
      enddo

    elseif (grid_ic) then

      pp=0
      do k=0,nf_physical_node_dim-1,2
        do j=0,nf_physical_node_dim-1,2
          do i=0,nf_physical_node_dim-1,2
            pp=pp+1
            xv(1,pp) = i + 0.5
            xv(2,pp) = j + 0.5
            xv(3,pp) = k + 0.5
            xv(4:6,pp) = 0.0
          enddo
        enddo
      enddo
      if (pp .ne. np_local) then
        write(*,*) 'pp ne np_local',pp,np_local,'rank',rank
        stop
      endif
      write(*,*) 'rank',rank,'finished grid_ic'

    elseif (pairwise_ic.or.pair_infall) then

      if (rank == 0) then
        np_local=2
        xv(:3,1:2) = 1.0
        xv(4:,1:2) = 0.0
      else
        np_local=0
      endif


    elseif (pp_test) then

      if (rank==0) then
        np_local=4
        xv(1:3,1)=0.3+nf_physical_node_dim/2

        xv(1:3,2)=0.2+nf_physical_node_dim/2

        xv(1,3)=0.2+nf_physical_node_dim/2
        xv(2:3,3)=0.3+nf_physical_node_dim/2

        xv(1,4)=0.3+nf_physical_node_dim/2
        xv(2:3,4)=0.2+nf_physical_node_dim/2

        xv(4:,1:4)=0.0
      else
        np_local=0
      endif

    elseif (restart_ic) then

      if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
      call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

      write(z_s,'(f7.3)') z_write
      z_s=adjustl(z_s)

      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)

      ofile=output_path//z_s(1:len_trim(z_s))//'xv'// &
            rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
      open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
      open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

#ifdef PPINT
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#else
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#endif
      if (rank == 0) print *,'restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

!     blocksize=(2047*1024*1024)/24
!     reduced to 32MB chunks because of intel compiler
      blocksize=(32*1024*1024)/24
      num_writes=np_local/blocksize+1


      !read(21) xv(:,:np_local)
      do i = 1,num_writes
         nplow=(i-1)*blocksize+1
         nphigh=min(i*blocksize,np_local)
!!       print *,rank,nplow,nphigh,np_local
         do j=nplow,nphigh
            read(21) xv(:,j)
         enddo
      enddo
 
      close(21)

#ifdef PID_FLAG

      ofile=output_path//z_s(1:len_trim(z_s))//'PID'// &
            rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
      open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
      open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

#ifdef PPINT
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#else
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#endif

      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      !read(21) PID(:np_local)
      do i = 1,num_writes
         nplow=(i-1)*blocksize+1
         nphigh=min(i*blocksize,np_local)
!!       print *,rank,nplow,nphigh,np_local
         do j=nplow,nphigh
            read(21) PID(j)
         enddo
      enddo
      close(21)

#endif


#ifdef CHECKPOINT_KILL

    elseif (restart_kill) then

      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)

      ofile=output_path//reskill_prefix//'xvres'// &
            rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
      open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
      open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

#ifdef PPINT
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#else
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#endif

      if (rank == 0) print *,'restarting simulation from z=', reskill_prefix
      !if (rank == 0) print *,'restarting simulation from z=',z_checkpoint(cur_checkpoint-1)
      if (rank == 0) print *,'current checkpoint, proj and halo entries are:', cur_checkpoint, &
               cur_projection,cur_halofind
      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

!     blocksize=(2047*1024*1024)/24
!     reduced to 32MB chunks because of intel compiler
      blocksize=(32*1024*1024)/24
      num_writes=np_local/blocksize+1

      !read(21) xv(:,:np_local)
      do i = 1,num_writes
         nplow=(i-1)*blocksize+1
         nphigh=min(i*blocksize,np_local)
!!       print *,rank,nplow,nphigh,np_local
         do j=nplow,nphigh
            read(21) xv(:,j)
         enddo
      enddo
      close(21)

#ifdef PID_FLAG

      ofile=output_path//reskill_prefix//'PIDres'// &
            rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
      open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
      open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

#ifdef PPINT
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#else
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p
#endif

      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      !read(21) PID(:np_local)
      do i = 1,num_writes
         nplow=(i-1)*blocksize+1
         nphigh=min(i*blocksize,np_local)
!!       print *,rank,nplow,nphigh,np_local
         do j=nplow,nphigh
            read(21) PID(j)
         enddo
      enddo
      close(21)

#endif

#endif

    elseif (shake_test_ic) then
      np_local=1
      xv(:,1)=(/0.0,0.0,0.0,0.,0.,0./) 
    else

!! Read in initial conditions

      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)

      ofile=ic_path//'xv'//rank_s(1:len_trim(rank_s))//'.ic'
      print *,'opening particle list:',ofile(1:len_trim(ofile))
#ifdef BINARY
      open(unit=20,file=ofile,form='binary',iostat=fstat,status='old')
#else
      open(unit=20,file=ofile,form='unformatted',iostat=fstat,status='old')
#endif
      if (fstat /= 0) then
         write(*,*) 'error opening initial conditions'
         write(*,*) 'rank',rank,'file:',ofile
         call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(20) np_local
      !      np_local=100**3
      if (np_local > max_np) then
         write(*,*) 'too many particles to store'
         write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
         call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
#ifndef BINARY
      ! this was for the old initial condition generator 
       !     read(20) xv(:3,:np_local)
            do i=1,np_local
              read(20) xv(:,i)
            enddo
#else
      read(20) xv(:,:np_local)
#endif
      close(20)

#ifdef PID_FLAG
      write(*,*) 'np_local before delete', np_local, 'rank =', rank
      !call delete_particles
      !write(*,*) 'np_local after delete', np_local, 'rank =', rank


      do i=1,np_local
         ! Assign a uniqe ID to particles in physical volumes.        
         ! This assumes that every node starts with the same np_local
         !PID(i) = int(i + rank*np_local,kind=8)
         PID(i) = int(i,kind=8) + int(rank*np_local,kind=8)
      enddo

#ifdef DEBUG_PID_INIT
      write(*,*) 'np_local=',np_local
      write(*,*) 'xv',xv(1:3,1:10), 'PID', PID(1:10) !
      write(*,*) 'xv',xv(1:3,np_local), 'PID', PID(np_local) ! should be non-zero
      write(*,*) 'xv',xv(1:3,np_local+1), 'PID', PID(np_local+1) ! should be zero
#endif

    if(rank==0)write(*,*) 'PID initialized' 

    ! Write initial PIDs to a file:
    fstat = 0

#ifdef BINARY
      open(unit=21,file=ic_path//'PID'//rank_s(1:len_trim(rank_s))//'.ic',form='binary',iostat=fstat)
#else
      open(unit=21,file=ic_path//'PID'//rank_s(1:len_trim(rank_s))//'.ic',form='unformatted',iostat=fstat)
#endif
      if (fstat /= 0) then
         write(*,*) 'error writing initial PID'
         write(*,*) 'rank',rank,'file:',ic_path//'PID'//rank_s(1:len_trim(rank_s))//'.ic'
         call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      write(21) np_local
      write(21) PID(:np_local)
      close(21)
      !stop

#endif
    
   endif

!! calculate total number of particles and particle mass

    npl8=int(np_local,kind=8)
    call mpi_reduce(npl8,np_total,1,MPI_INTEGER8, &
                           mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'number of particles =', np_total
    call mpi_bcast(np_total,1,MPI_INTEGER8,0,mpi_comm_world,ierr)

    if (.not.restart_ic) then
       if (pairwise_ic) then
          mass_p=10000.0
       elseif (pair_infall) then
          mass_p=pair_infall_mass
       elseif (pp_test) then
          mass_p=10000.0/4.
       else
!#ifdef MHD
!          mass_p = (real(nf_physical_dim)**3 / real(np_total))!*(1-omega_b/omega_m)
!#else
          mass_p = real(nf_physical_dim)**3 / real(np_total)
!#endif
       endif
    endif

    if (rank == 0) write(*,*) 'particle mass=', mass_p
    if (rank == 0) write(*,*) 'total dark matter mass =', mass_p * np_total

!! This is to scale the initial conditions if we are doing testing with another data-set
#ifdef SCALED_IC
    do i=1,np_local
      xv(:,i)=xv(:,i)/4.0
    enddo
#endif

! this is to test if we can reconstruct particle ID from drifting vs non-drifting populations
#ifdef X_DRIFT
    do i=1,np_local
      xv(4,i)=xv(4,i) + 10
    enddo
#endif

!! check to make sure no particles are out of bounds
!! this is not really necessary, as link_list can handle out of bounds
#ifdef CHECK_IP 
    if (.not.restart_ic) then
      xva=0.0
      do i=1,np_local
        xva=xva+xv(:,i)
        do j=1,3
          if (xv(j,i) < 0.0 .or. xv(j,i) >= real(nf_physical_node_dim,4)) then
            write(*,*) 'particle out of bounds'
            write(*,*) xv(:,i)
            call mpi_abort(mpi_comm_world,ierr,ierr)
          endif
        enddo
      enddo
      xva=xva/real(np_local)
      print *,rank,xva 
    endif
#endif

  end subroutine particle_initialize
