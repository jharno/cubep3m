!! initialize particle list
  subroutine particle_initialize
    use omp_lib
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    real(4) :: rnum,z_write,dummy
    integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh
    integer*8 :: np_total,npl8
    character(len=max_path) :: ofile
    character(len=6) :: rank_s
    character(len=7) :: z_s, z_s2
    integer(4) :: np_nu
#ifdef CHECK_IP
    real(8) :: xva(6)
#endif
#ifdef NEUTRINOS
    integer(4) :: np_dm
    integer(8) :: np_total_nu

    !! Check that the PID flag is also defined
#ifndef PID_FLAG
    write(*,*) "ERROR: Using Neutrinos but PID is not enabled !!"
    call mpi_abort(mpi_comm_world,ierr,ierr)
#endif
#endif
#if defined(ZIP) || defined(ZIPDM)
    character(len=max_path) :: f_zip0,f_zip1,f_zip2,f_zip3
    integer :: fstat0, fstat1, fstat2, fstat3, l
    real(4) :: v_r2i
    integer(4) :: np_uzip
    integer(1) :: xi1(4,3), rhoc_i1(4), test_i1
    integer(4) :: xi4(3), rr_i4
    integer(2) :: vi2(3)
    equivalence(xi1,xi4)
    equivalence(rr_i4,rhoc_i1)
#endif
print*,'particle_initialize'
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

! ---------------------------------------------------------------------------------------------
! Read in checkpoint
! ---------------------------------------------------------------------------------------------

      if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
      call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
      write(z_s,'(f7.3)') z_write
      z_s=adjustl(z_s)

      write(rank_s,'(i6)') rank
      rank_s=adjustl(rank_s)

      !
      ! Read dark matter particles
      !

#ifndef ZIPDM

      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
            rank_s(1:len_trim(rank_s))//'.dat'
      if (rank==0) print*, 'opening dark matter checkpoint file:',ofile

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p

      if (rank == 0) print *,'restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

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

#ifdef ZIP
      np_uzip = np_local !! Need to do if ZIP is used for neutrinos 
#endif

#else

    f_zip0 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip0_'//rank_s(1:len_trim(rank_s))//'.dat'
    f_zip1 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'.dat'
    f_zip2 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'.dat'
    f_zip3 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'.dat'
    open(10, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
    open(11, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
    open(12, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
    open(13, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')

    if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
        write(*,*) 'error opening dm zip checkpoint'
        write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    !! Read header (which incluedes velocity conversion factor and shake_offset)

    read(10) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
    read(11) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

    if (rank == 0) print *,'restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

    if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    !! Read zipped positions and velocities

    np_uzip = 0

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim
            do i = 1, nc_node_dim
                rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#ifndef BGQ
                read(12) rhoc_i1(1) !! Get number of particles in the coarse grid
#else
                read(12) rhoc_i1(4)
#endif
                if (rr_i4==255) read(13) rr_i4
                do l = 1, rr_i4
                    xi1=0; xi4=0
                    np_uzip = np_uzip + 1
#ifndef BGQ
                    read(10) xi1(1,:)
#else
                    read(10) xi1(4,:)
#endif
                    read(11) vi2
                    xv(1:3, np_uzip) = mesh_scale * ( (xi4+0.5)/256. + (/i,j,k/) - 1 )
                    xv(4:6, np_uzip) = vi2 / v_r2i
                enddo
            enddo
        enddo
    enddo

    !! Consistency checks
    read(10, end=701) test_i1
    print*, 'ERROR: rank', rank, ': file not ending:', f_zip0
    call mpi_abort(mpi_comm_world, ierr, ierr)
701 close(10)
    read(11, end=801) test_i1
    print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
    call mpi_abort(mpi_comm_world, ierr, ierr)
801 close(11) ; close(12) ; close(13)

    if (np_uzip /= np_local) then
        write(*,*) "Something wrong with reading dark matter zipped files: ", rank, np_local, np_uzip
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#endif

      !
      ! Read neutrino particles 
      !

#ifdef NEUTRINOS

#ifndef ZIP

        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'
        if (rank==0) print*, 'opening neutrino checkpoint file:',ofile

        open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error opening checkpoint'
            write(*,*) 'rank',rank,'file:',ofile
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
        read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
      
        if (rank == 0) print *,'neutrinos restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

        if (np_local+np_nu > max_np) then
            write(*,*) 'too many particles to store'
            write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

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

    !! Open neutrino files for reading

    f_zip0 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip0_'//rank_s(1:len_trim(rank_s))//'_nu.dat' 
    f_zip1 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
    f_zip2 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
    f_zip3 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'_nu.dat'

    open(20, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
    open(21, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
    open(22, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
    open(23, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')
    if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
        write(*,*) 'error opening dm zip checkpoint'
        write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    !! Read header (which incluedes velocity conversion factor and shake_offset). Only store the local number
    !! of particles and velocity conversion factor. All other info already read from dark matter checkpoint file.

    read(20) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,v_r2i,dummy,dummy,dummy
    read(21) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

    if (rank == 0) print *,'neutrinos restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

    if (np_local+np_nu > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'np_nu',np_nu,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    !! Read zipped positions and velocities
    do k = 1, nc_node_dim
        do j = 1, nc_node_dim
            do i = 1, nc_node_dim
                rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#ifndef BGQ
                read(22) rhoc_i1(1) !! Get number of particles in the coarse grid
#else
                read(22) rhoc_i1(4)
#endif
                if (rr_i4==255) read(23) rr_i4
                do l = 1, rr_i4
                    xi1=0; xi4=0
                    np_uzip = np_uzip + 1
#ifndef BGQ
                    read(20) xi1(1,:)
#else
                    read(20) xi1(4,:)
#endif
                    read(21) vi2
                    xv(1:3, np_uzip) = mesh_scale * ( (xi4+0.5)/256. + (/i,j,k/) - 1 )
                    xv(4:6, np_uzip) = vi2 / v_r2i
                enddo
            enddo
        enddo
    enddo
    !! Consistency checks
    read(20, end=702) test_i1
    print*, 'ERROR: rank', rank, ': file not ending:', f_zip0
    call mpi_abort(mpi_comm_world, ierr, ierr)
702 close(20) 
    read(21, end=802) test_i1
    print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
    call mpi_abort(mpi_comm_world, ierr, ierr)
802 close(21) ; close(22) ; close(23)

    if (np_uzip /= np_local+np_nu) then
        write(*,*) "Something wrong with reading neutrino zipped files: ", rank, np_local, np_nu, np_uzip
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#endif

#ifdef NUPID
      !! DM PIDs are 0 
      do j = 1, np_local
        PID(j) = 0
      enddo

      !! Read neutrino PIDs
      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PID'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      read(21) np_dm,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

      if (np_dm /= np_nu) then
        write(*,*) "ERROR WITH NEUTRINO PID FILE: ", rank, np_dm, np_nu, ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      blocksize=(32*1024*1024)/24
      num_writes=np_nu/blocksize+1

      do i = 1,num_writes
         nplow=(i-1)*blocksize+1 + np_local
         nphigh=min(i*blocksize,np_nu) + np_local
         do j=nplow,nphigh
            read(21) PID(j)
         enddo
      enddo
      close(21)
#endif

#else

#ifdef PID_FLAG

      !
      ! Read dark matter IDs (if not a neutrino simulation) 
      !

      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PID'// &
            rank_s(1:len_trim(rank_s))//'.dat'

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p

      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

!     reduced to 32MB chunks because of intel compiler
      blocksize=(32*1024*1024)/24
      num_writes=np_local/blocksize+1

      do i = 1,num_writes
         nplow=(i-1)*blocksize+1
         nphigh=min(i*blocksize,np_local)
         do j=nplow,nphigh
            read(21) PID(j)
         enddo
      enddo
      close(21)

#endif

#endif

! ---------------------------------------------------------------------------------------------
! Done read in checkpoint
! ---------------------------------------------------------------------------------------------

#ifdef CHECKPOINT_KILL

! ---------------------------------------------------------------------------------------------
! Read in checkpoint kill
! ---------------------------------------------------------------------------------------------

    elseif (restart_kill) then

      write(rank_s,'(i6)') rank
      rank_s=adjustl(rank_s)

      !
      ! Read dark matter particles
      !

#ifdef ZIPDM
     if (restart_kill_zip .eqv. .false.) then !! Happens if checkpoint_kill occured in particle_pass
#endif

          ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'xvres'// &
                rank_s(1:len_trim(rank_s))//'.dat'

          open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
          if (fstat /= 0) then
            write(*,*) 'error opening checkpoint'
            write(*,*) 'rank',rank,'file:',ofile
            call mpi_abort(mpi_comm_world,ierr,ierr)
          endif

          read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
                   cur_projection,cur_halofind,mass_p

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

#ifdef ZIP
          np_uzip = np_local !! Need to do if ZIP is used for neutrinos 
#endif

#ifdef ZIPDM
    else

        f_zip0 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres0_'//rank_s(1:len_trim(rank_s))//'.dat'
        f_zip1 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres1_'//rank_s(1:len_trim(rank_s))//'.dat'
        f_zip2 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres2_'//rank_s(1:len_trim(rank_s))//'.dat'
        f_zip3 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres3_'//rank_s(1:len_trim(rank_s))//'.dat'

        open(10, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
        open(11, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
        open(12, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
        open(13, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')
        if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
            write(*,*) 'error opening dm zip checkpoint'
            write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        !! Read header (which incluedes velocity conversion factor and shake_offset)

        read(10) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
        read(11) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

        if (rank == 0) print *,'restarting simulation from z=',reskill_prefix
        if (rank == 0) print *,'current checkpoint, proj and halo entries are:', cur_checkpoint, cur_projection,cur_halofind

        if (np_local > max_np) then
            write(*,*) 'too many particles to store'
            write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        !! Read zipped positions and velocities
        np_uzip = 0
        do k = 1, nc_node_dim
            do j = 1, nc_node_dim
                do i = 1, nc_node_dim
                    rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#ifndef BGQ
                    read(12) rhoc_i1(1) !! Get number of particles in the coarse grid
#else
                    read(12) rhoc_i1(4)
#endif
                    if (rr_i4==255) read(13) rr_i4
                    do l = 1, rr_i4
                        xi1=0; xi4=0
                        np_uzip = np_uzip + 1
#ifndef BGQ
                        read(10) xi1(1,:)
#else
                        read(10) xi1(4,:)
#endif
                        read(11) vi2
                        xv(1:3, np_uzip) = mesh_scale * ( (xi4+0.5)/256. + (/i,j,k/) - 1 )
                        xv(4:6, np_uzip) = vi2 / v_r2i
                    enddo
                enddo
            enddo
        enddo

        !! Consistency checks
        read(10, end=703) test_i1
        print*, 'ERROR: rank', rank, ': file not ending:', f_zip0
        call mpi_abort(mpi_comm_world, ierr, ierr)
703     close(10)
        read(11, end=803) test_i1
        print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
        call mpi_abort(mpi_comm_world, ierr, ierr)
803     close(11) ; close(12) ; close(13)

        if (np_uzip /= np_local) then
            write(*,*) "Something wrong with reading dark matter zipped files: ", rank, np_local, np_uzip
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

    endif
#endif

#ifdef NEUTRINOS

      !
      ! Read neutrino particles 
      !

#ifdef ZIP
     if (restart_kill_zip .eqv. .false.) then !! Happens if checkpoint_kill occured in particle_pass
#endif

        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'xvres'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'

        open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error opening checkpoint'
            write(*,*) 'rank',rank,'file:',ofile
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
        read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

        if (rank == 0) print *,'neutrinos restarting simulation from z=',reskill_prefix

        if (np_local+np_nu > max_np) then
            write(*,*) 'too many particles to store'
            write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

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

#ifdef ZIP
    else

        !! Open neutrino files for reading

        f_zip0 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres0_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        f_zip1 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres1_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        f_zip2 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres2_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        f_zip3 = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'zipres3_'//rank_s(1:len_trim(rank_s))//'_nu.dat'

        open(20, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
        open(21, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
        open(22, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
        open(23, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')
        if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
            write(*,*) 'error opening dm zip checkpoint'
            write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        !! Read header (which incluedes velocity conversion factor and shake_offset). Only store the local number
        !! of particles and velocity conversion factor. All other info already read from dark matter checkpoint file.

        read(20) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,v_r2i,dummy,dummy,dummy
        read(21) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

        if (rank == 0) print *,'neutrinos restarting simulation from z=',reskill_prefix

        if (np_local+np_nu > max_np) then
            write(*,*) 'too many particles to store'
            write(*,*) 'rank',rank,'np_local',np_local,'np_nu',np_nu,'max_np',max_np
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        !! Read zipped positions and velocities

        do k = 1, nc_node_dim
            do j = 1, nc_node_dim
                do i = 1, nc_node_dim
                    rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#ifndef BGQ
                    read(22) rhoc_i1(1) !! Get number of particles in the coarse grid
#else
                    read(22) rhoc_i1(4)
#endif
                    if (rr_i4==255) read(23) rr_i4
                    do l = 1, rr_i4
                        xi1=0; xi4=0
                        np_uzip = np_uzip + 1
#ifndef BGQ
                        read(20) xi1(1,:)
#else
                        read(20) xi1(4,:)
#endif
                        read(21) vi2
                        xv(1:3, np_uzip) = mesh_scale * ( (xi4+0.5)/256. + (/i,j,k/) - 1 )
                        xv(4:6, np_uzip) = vi2 / v_r2i
                    enddo
                enddo
            enddo
        enddo

        !! Consistency checks
        read(20, end=704) test_i1
        print*, 'ERROR: rank', rank, ': file not ending:', f_zip0
        call mpi_abort(mpi_comm_world, ierr, ierr)
704     close(20); 
        read(21, end=804) test_i1
        print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
        call mpi_abort(mpi_comm_world, ierr, ierr)
804     close(21) ; close(22) ; close(23)

        if (np_uzip /= np_local+np_nu) then
            write(*,*) "Something wrong with reading neutrino zipped files: ", rank, np_local, np_nu, np_uzip
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

    endif
#endif

#ifdef NUPID
      !! DM PIDs are 0 
      do j = 1, np_local
        PID(j) = 0
      enddo

      !! Read neutrino PIDs
      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'PIDres'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      read(21) np_dm,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

      if (np_dm /= np_nu) then
        write(*,*) "ERROR WITH NEUTRINO PID FILE: ", rank, np_dm, np_nu, ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      blocksize=(32*1024*1024)/24
      num_writes=np_nu/blocksize+1

      do i = 1,num_writes
         nplow=(i-1)*blocksize+1 + np_local
         nphigh=min(i*blocksize,np_nu) + np_local
         do j=nplow,nphigh
            read(21) PID(j)
         enddo
      enddo
      close(21)
#endif

#else

#ifdef PID_FLAG

      !
      ! Read dark matter IDs (if not a neutrino simulation) 
      !

      ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'PIDres'// &
            rank_s(1:len_trim(rank_s))//'.dat'

      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
               cur_projection,cur_halofind,mass_p

      if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

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

    elseif (shake_test_ic) then
      np_local=1
      xv(:,1)=(/0.0,0.0,0.0,0.,0.,0./) 
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

      write(rank_s,'(i6)') rank
      rank_s=adjustl(rank_s)

      !
      ! Read dark matter particles
      !

#ifndef ZIPDM

      ofile=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
      if (rank==0) print *,'opening particle list:',ofile(1:len_trim(ofile))

      open(unit=20, file=ofile, status="old", iostat=fstat, access="stream")
      if (fstat /= 0) then
         write(*,*) 'error opening initial conditions'
         write(*,*) 'rank',rank,'file:',ofile
         call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(20) np_local,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
      if (np_local > max_np) then
         write(*,*) 'too many particles to store'
         write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
         call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(20) xv(:,:np_local)
      close(20)

#else

    f_zip0 = ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip0_'//rank_s(1:len_trim(rank_s))//'.dat' 
    f_zip1 = ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'.dat'
    f_zip2 = ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'.dat'
    f_zip3 = ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'.dat'
    open(10, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
    open(11, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
    open(12, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
    open(13, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')
    if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
        write(*,*) 'error opening dm zip checkpoint'
        write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    !! Read header (which incluedes velocity conversion factor and shake_offset)
    read(10) np_local,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,v_r2i,dummy,dummy,dummy
    read(11) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
    if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    !! Read zipped positions and velocities !! fixed for big endian
    np_uzip = 0
    do k = 1, nc_node_dim
        do j = 1, nc_node_dim
            do i = 1, nc_node_dim
                rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#ifndef BGQ
                read(12) rhoc_i1(1) !! Get number of particles in the coarse grid
#else
                read(12) rhoc_i1(4) !! big endian
#endif
!print*,'rhoc_i1, rr_i4 =',rhoc_i1,rr_i4
                if (rr_i4==255) read(13) rr_i4
                do l = 1, rr_i4
                    xi1=0; xi4=0
                    np_uzip = np_uzip + 1
#ifndef BGQ
                    read(10) xi1(1,:)
#else
                    read(10) xi1(4,:) !! big endian
#endif
!print*,'xi1, xi4 =',xi1,xi4
                    read(11) vi2
                    xv(1:3, np_uzip) = mesh_scale * ( (xi4+0.5)/256. + (/i,j,k/) - 1 )
                    xv(4:6, np_uzip) = vi2 / v_r2i
                enddo
            enddo
        enddo
    enddo

    !! Consistency checks

    read(10, end=705) test_i1
    print*, 'ERROR: rank', rank, ': file not ending:', f_zip0
    call mpi_abort(mpi_comm_world, ierr, ierr)
705 close(10)
    read(11, end=805) test_i1
    print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
    call mpi_abort(mpi_comm_world, ierr, ierr)
805 close(11) ; close(12) ; close(13)

    if (np_uzip /= np_local) then
        write(*,*) "Something wrong with reading dark matter zipped files: ", rank, np_local, np_uzip
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#endif

#ifdef NUPID
        !! Write PIDs to a binary file so that when neutrinos restart at z_i_nu they can read in the PIDs.
        np_nu = np_local * ratio_nudm_dim**3
        dummy = 0
        ofile=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s2(1:len_trim(z_s2))//'PID'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        open(unit=21, file=ofile, iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error writing initial PID', ofile
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        write(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
        do i = 1, np_nu
            write(21) int(i,kind=8) + int(rank*int(np_nu,kind=8),kind=8)
        enddo
        close(21)
#endif

#ifdef PID_FLAG
        if (rank == 0) write(*,*) 'np_local before delete', np_local, 'rank =', rank
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
        if (fstat /= 0) then
            write(*,*) 'error writing initial PID'
            write(*,*) 'rank',rank,'file:',ic_path//'PID'//rank_s(1:len_trim(rank_s))//'.ic'
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
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

    if (.not.restart_ic) then
       if (pairwise_ic) then
          mass_p=10000.0
       elseif (pair_infall) then
          mass_p=pair_infall_mass
       elseif (pp_test) then
          mass_p=10000.0/4.
       else
          mass_p = real(nf_physical_dim)**3 / real(np_total)
       endif
    endif

    if (rank == 0) write(*,*) 'particle mass=', mass_p
    if (rank == 0) write(*,*) 'total dark matter mass =', mass_p * np_total

#ifdef NEUTRINOS
#ifndef NUPID
    !
    ! Assign 1 byte integer PIDs
    !

    do i = 1, np_local
        PID(i) = 1 !! All dark matter will have a PID of 1
    enddo
    do i = np_local+1, np_local+np_nu
        PID(i) = 2 !! All neutrinos will have a PID of 2 
    enddo
#endif

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
