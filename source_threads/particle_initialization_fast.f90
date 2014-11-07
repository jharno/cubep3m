!! initialize particle list
  subroutine particle_initialize
    use omp_lib
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    real(4) :: rnum,z_write,dummy, v_r2i
    integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh
    integer*8 :: np_total,npl8
    character(len=max_path) :: ofile1,ofile2,ofile3
    character(len=4) :: rank_s
    character(len=7) :: z_s, z_s2
    integer(4) :: np_nu
#ifdef NEUTRINOS
    integer(4) :: np_dm
    integer(8) :: np_total_nu

    !! Check that the PID flag is also defined
#ifndef PID_FLAG
    write(*,*) "ERROR: Using Neutrinos but PID is not enabled !!"
    call mpi_abort(mpi_comm_world,ierr,ierr)
#endif
#endif

fstat=0
np_local=(nf_physical_node_dim/2)**3

if (random_ic) then
elseif (restart_ic) then ! for dm+nu

  if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)

  write(rank_s,'(i4)') rank
  rank_s=adjustl(rank_s)

  ofile1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
  ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoc'//rank_s(1:len_trim(rank_s))//'.dat'
  ofile3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoce'//rank_s(1:len_trim(rank_s))//'.dat'


  open(unit=21, file=ofile1, status="old", iostat=fstat, access="stream")
  open(unit=22, file=ofile2, status="old", iostat=fstat, access="stream")
  open(unit=23, file=ofile3, status="old", iostat=fstat, access="stream")

  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i
  read(22) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i
  read(23) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i

do k=1,nc_node_dim
do j=1,nc_node_dim
do i=1,nc_node_dim
  !! get number of particles in the coarse grid
  read(22) rho_c_dm_int1
  if (rho_c_dm_int1<127) then
    rho_c_dm_int4=rho_c_dm_int1+128
  else
    read(23) rho_c_dm_int4
  endif

  do l=1,rho_c_dm_int4
    np_local=np_local+1
    xv(1:3,l)=
    xv(4:6,l)=
  enddo
enddo ! i
enddo ! j
enddo ! k 

!     reduced to 32MB chunks because of intel compiler
      blocksize=(32*1024*1024)/24
      num_writes=np_local/blocksize+1

      !read(21) xv(:,:np_local)
      do i = 1,num_writes
         nplow=(i-1)*blocksize+1
         nphigh=min(i*blocksize,np_local)
         do j=nplow,nphigh
            read(21) xv(:,j)
         enddo
      enddo
      close(21)

#ifdef NEUTRINOS
        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'
        if (rank==0) print*, 'opening neutrino checkpoint file:',ofile

        open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")

        !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
        read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
      
        if (rank == 0) print *,'neutrinos restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

!       reduced to 32MB chunks because of intel compiler
        blocksize=(32*1024*1024)/24
        num_writes=np_nu/blocksize+1

        do i = 1,num_writes
            nplow=(i-1)*blocksize+1 + np_local
            nphigh=min(i*blocksize,np_nu) + np_local
            do j=nplow,nphigh
                read(21) xv(:,j)
            enddo
        enddo
        close(21)
#endif

#ifdef CHECKPOINT_KILL

    elseif (restart_kill) then

      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)

      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'xvres'// &
            rank_s(1:len_trim(rank_s))//'.dat'

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")

      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p

      if (rank == 0) print *,'restarting simulation from z=', reskill_prefix
      !if (rank == 0) print *,'restarting simulation from z=',z_checkpoint(cur_checkpoint-1)
      if (rank == 0) print *,'current checkpoint, proj and halo entries are:', cur_checkpoint, &
               cur_projection,cur_halofind

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

#ifdef NEUTRINOS
        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'xvres'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'

        open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")

        !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
        read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

        if (rank == 0) print *,'neutrinos restarting simulation from z=',reskill_prefix

!       reduced to 32MB chunks because of intel compiler
        blocksize=(32*1024*1024)/24
        num_writes=np_nu/blocksize+1

        do i = 1,num_writes
            nplow=(i-1)*blocksize+1 + np_local
            nphigh=min(i*blocksize,np_nu) + np_local
            do j=nplow,nphigh
                read(21) xv(:,j)
            enddo
        enddo

        close(21)


#else
#ifdef PID_FLAG

      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'PIDres'// &
            rank_s(1:len_trim(rank_s))//'.dat'

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")

      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p


!     reduced to 32MB chunks because of intel compiler
      blocksize=(32*1024*1024)/24
      num_writes=np_local/blocksize+1

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

! ---------------------------------------------------------------------------------------------
! Done read in checkpoint kill
! ---------------------------------------------------------------------------------------------

#endif

    else

! ---------------------------------------------------------------------------------------------
! Read in ICs
! ---------------------------------------------------------------------------------------------

      write(z_s,'(f7.3)') z_i
      z_s=adjustl(z_s)
#ifdef NUPID
      write(z_s2,'(f7.3)') z_i_nu
      z_s2=adjustl(z_s2)
#endif

      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)

      ofile=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
      if (rank==0) print *,'opening particle list:',ofile(1:len_trim(ofile))

      open(unit=20, file=ofile, status="old", iostat=fstat, access="stream")
      read(20) np_local,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
      read(20) xv(:,:np_local)
      close(20)

#ifdef PID_FLAG
        write(*,*) 'np_local before delete', np_local, 'rank =', rank
        !call delete_particles
        !write(*,*) 'np_local after delete', np_local, 'rank =', rank

        do i=1,np_local
            ! Assign a uniqe ID to particles in physical volumes.        
            ! This assumes that every node starts with the same np_local
            !PID(i) = int(i + rank*np_local,kind=8)
            PID(i) = int(i,kind=8) + int(rank*int(np_local,kind=8),kind=8)
        enddo

        if(rank==0)write(*,*) 'PID initialized' 

        !
        ! Write initial ICs to a file
        !

        fstat = 0

        open(unit=21, file=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'PID'//rank_s(1:len_trim(rank_s))//'.ic', iostat=fstat, access="stream")
        write(21) np_local
        write(21) PID(:np_local)
        close(21)

#endif

! ---------------------------------------------------------------------------------------------
! Done read in ICs
! ---------------------------------------------------------------------------------------------

   endif

!! calculate total number of particles and particle mass

    npl8=int(np_local,kind=8)
    call mpi_reduce(npl8,np_total,1,MPI_INTEGER8, &
                           mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'number of particles =', np_total
    call mpi_bcast(np_total,1,MPI_INTEGER8,0,mpi_comm_world,ierr)

          mass_p = real(nf_physical_dim)**3 / real(np_total)
    endif

    if (rank == 0) write(*,*) 'particle mass=', mass_p
    if (rank == 0) write(*,*) 'total dark matter mass =', mass_p * np_total

#ifdef NEUTRINOS

    !
    ! Print some stats to screen
    !

    npl8=int(np_nu,kind=8)
    call mpi_reduce(npl8,np_total_nu,1,MPI_INTEGER8, mpi_sum,0,mpi_comm_world,ierr)
    !! Append np_nu to np_local (must be done after PIDs and after mass calculation)
    np_dm = np_local
    np_local = np_local + np_nu
    if (rank == 0) then
        write(*,*) "np_dm = ", np_dm
        write(*,*) "np_nu = ", np_nu
        write(*,*) "np_local = ", np_local
        write(*,*) "np_dm_total = ", np_total
        write(*,*) "np_nu_total = ", np_total_nu
    endif
#endif

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

  end subroutine particle_initialize
