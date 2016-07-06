! write checkpoints to disk
#if defined(ZIP) || defined(ZIPDM)
  subroutine checkpoint_kill(dozip)
#else
  subroutine checkpoint_kill
#endif

#ifdef MHD
    use mpi_tvd_mhd
#endif

    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

#if defined(ZIP) || defined(ZIPDM)
    logical :: dozip
#endif

    character (len=max_path) :: ofile,ofile2
    character (len=6) :: rank_s
    character (len=7) :: z_s  

    integer(kind=4) :: i,j,fstat,blocksize,num_writes,nplow,nphigh
    real(kind=4) :: z_write
#ifdef NEUTRINOS
    integer(4) :: np_dm, np_nu, ind_check1, ind_check2
    character (len=max_path) :: ofile_nu, ofile2_nu
#ifdef NUPID
    integer, parameter :: pdm = 0
#else
    integer, parameter :: pdm = 1
#endif
#endif
#if defined(ZIP) || defined(ZIPDM)
    integer :: l, k
    integer :: ind, ind_nu
    integer :: fstat0, fstat1, fstat2, fstat3
    character (len=max_path) :: fdm_zip0,fdm_zip1,fdm_zip2,fdm_zip3
    character (len=max_path) :: fnu_zip0,fnu_zip1,fnu_zip2,fnu_zip3
    real(4) :: v_resolution = 32767.499
    real(4) :: vmax, vmax_local, v_r2i, v_r2i_nu
    integer(4) :: rhoc_dm_i4, rhoc_nu_i4
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

    !
    ! Checkpoint file names
    !

    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)

    write(z_s,'(f7.3)') z_write
    z_s=adjustl(z_s)

    !! Dark matter file
#ifdef ZIPDM
    if (dozip) then
        fdm_zip0=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip0res_'//rank_s(1:len_trim(rank_s))//'.dat'
        fdm_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1res_'//rank_s(1:len_trim(rank_s))//'.dat'
        fdm_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2res_'//rank_s(1:len_trim(rank_s))//'.dat'
        fdm_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3res_'//rank_s(1:len_trim(rank_s))//'.dat'
    else
#endif
        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xvres'//rank_s(1:len_trim(rank_s))//'.dat'
#ifdef ZIPDM
    endif
#endif

#ifdef NEUTRINOS
    !! Neutrino file
#ifdef ZIP
    if (dozip) then
        fnu_zip0=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip0res_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        fnu_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1res_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        fnu_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2res_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
        fnu_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3res_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
    else
#endif
    ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xvres'//rank_s(1:len_trim(rank_s))//'_nu.dat'
#ifdef ZIP
    endif
#endif
    !! Neutrino PID file
#ifdef NUPID
    ofile2_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PIDres'// &
         rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif
#else
    !! Dark matter PID file (if not a neutrino sim)
#ifdef PID_FLAG
    ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PIDres'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#endif
#endif

#if defined(ZIP) || defined(ZIPDM)
    !
    ! Determine velocity conversion factor
    !
#ifdef NEUTRINOS
    vmax_local = maxval(abs(xv(4:6,1:np_local)),mask=spread(PID(1:np_local)==pdm,dim=1,ncopies=3))
    call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
    v_r2i = v_resolution/vmax
    vmax_local = maxval(abs(xv(4:6,1:np_local)),mask=spread(PID(1:np_local)/=pdm,dim=1,ncopies=3))
    call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
    v_r2i_nu = v_resolution/vmax
#else
    vmax_local = maxval(abs(xv(4:6,1:np_local)))
    call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
    v_r2i = v_resolution/vmax
#endif
#endif

    !
    ! Open checkpoint files
    !

    !! Dark matter 
#ifdef ZIPDM
    if (dozip) then
        open(unit=10, file=fdm_zip0, status="replace", iostat=fstat0, access="stream", buffered='yes')
        open(unit=11, file=fdm_zip1, status="replace", iostat=fstat1, access="stream", buffered='yes')
        open(unit=12, file=fdm_zip2, status="replace", iostat=fstat2, access="stream", buffered='yes')
        open(unit=13, file=fdm_zip3, status="replace", iostat=fstat3, access="stream", buffered='yes')
        if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
          write(*,*) 'error opening zipped checkpoint file for write'
          write(*,*) 'rank',rank,'file:',fdm_zip0,fdm_zip1,fdm_zip2,fdm_zip3
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
    else
#endif
        open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
          write(*,*) 'error opening checkpoint file for write'
          write(*,*) 'rank',rank,'file:',ofile
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
#ifdef ZIPDM
    endif
#endif

#ifdef NEUTRINOS

    !! Neutrinos 
#ifdef ZIP
    if (dozip) then
        open(unit=20, file=fnu_zip0, status="replace", iostat=fstat0, access="stream", buffered='yes')
        open(unit=21, file=fnu_zip1, status="replace", iostat=fstat1, access="stream", buffered='yes')
        open(unit=22, file=fnu_zip2, status="replace", iostat=fstat2, access="stream", buffered='yes')
        open(unit=23, file=fnu_zip3, status="replace", iostat=fstat3, access="stream", buffered='yes')
        if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
          write(*,*) 'error opening zipped checkpoint file for write'
          write(*,*) 'rank',rank,'file:',fnu_zip1,fnu_zip2,fnu_zip3
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
    else
#endif
        open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
          write(*,*) 'error opening checkpoint file for write'
          write(*,*) 'rank',rank,'file:',ofile_nu
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
#ifdef ZIP
    endif
#endif

#ifdef NUPID
    !! Neutrino PID file 
    open(unit=25, file=ofile2_nu, status="replace", iostat=fstat, access="stream")
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile2_nu
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#endif

    !
    ! Write file headers
    !

    !! Determine how many dark matter and neutrino particles this rank has
    np_dm = count(PID(1:np_local) == pdm)
    np_nu = np_local - np_dm
    if (rank == 0) write(*,*) "checkpoint np_dm, np_nu = ", np_dm, np_nu

#ifdef ZIPDM
    if (dozip) then
        write(10) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
        write(11) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
    else
#endif
        write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p
#ifdef ZIPDM
    endif
#endif

#ifdef ZIP
    if (dozip) then
        write(20) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i_nu,shake_offset
        write(21) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i_nu,shake_offset
    else
#endif
        write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p
#ifdef ZIP
    endif
#endif

#ifdef NUPID
    write(25) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p
#endif

#else

#ifdef ZIPDM
    if (dozip) then
        write(10) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
        write(11) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
    else
#endif
        write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p
#ifdef ZIPDM
    endif
#endif
#endif

#ifdef NEUTRINOS

    !
    ! Write data for neutrino simulation
    !

    ind_check1 = 0
    ind_check2 = 0

#if defined(ZIP) || defined(ZIPDM)
    if (dozip) then

        !! Zipped format

        do k = 1, nc_node_dim
            do j = 1, nc_node_dim
                !$omp parallel default(shared) private(i,l,ind,ind_nu)
                !$omp do schedule(static)
                do i = 1, nc_node_dim
                    ind = 1 ; ind_nu = 1
                    rhoc_i4(i) = 0 ; rhoc_i4_nu(i) = 0
                    l = hoc(i,j,k)
                    do while (l > 0)
                        if (PID(l) == pdm) then
                            rhoc_i4(i) = rhoc_i4(i) + 1
#ifdef ZIPDM
                            pos_i1(:, ind, i) = int(mod(xv(1:3,l)/mesh_scale,1.)*256, kind=1)
                            vel_i2(:, ind, i) = nint(xv(4:6,l)*v_r2i, kind=2)
#else
                            pos_i1(:, ind, i) = xv(1:6,l)
#endif
                            ind = ind + 1
                        else
                            rhoc_i4_nu(i) = rhoc_i4_nu(i) + 1
#ifdef ZIP
                            pos_i1_nu(:, ind_nu, i) = int(mod(xv(1:3,l)/mesh_scale,1.)*256, kind=1)
                            vel_i2_nu(:, ind_nu, i) = nint(xv(4:6,l)*v_r2i_nu, kind=2)
#else
                            pos_i1_nu(:, ind_nu, i) = xv(1:6,l)
#endif
#ifdef NUPID
                            pid_i8_nu(ind_nu, i) = PID(l)
#endif
                            ind_nu = ind_nu + 1
                        endif
                        l = ll(l)
                    enddo !! l
                    if (ind > max_cell_np+1 .or. ind_nu > max_cell_np+1) then
                        write(*,*) "ERROR: max_cell_np too small !!", rank, ind, ind_nu, max_cell_np
                        call mpi_abort(mpi_comm_world, ierr, ierr)
                    endif
                enddo !! i
                !$omp end do
                !$omp end parallel
                ind_check1 = ind_check1 + sum(rhoc_i4)
                ind_check2 = ind_check2 + sum(rhoc_i4_nu)
#ifdef ZIPDM
                write(12) int(min(rhoc_i4, 255), kind=1)
                write(13) pack(rhoc_i4, rhoc_i4>254)
#endif
#ifdef ZIP
                write(22) int(min(rhoc_i4_nu, 255), kind=1)
                write(23) pack(rhoc_i4_nu, rhoc_i4_nu>254)
#endif
                do i = 1, nc_node_dim
#ifdef ZIPDM
                    write(10) pos_i1(:, 1:rhoc_i4(i), i)
                    write(11) vel_i2(:, 1:rhoc_i4(i), i)
#else
                    write(12) pos_i1(:, 1:rhoc_i4(i), i)
#endif
#ifdef ZIP
                    write(20) pos_i1_nu(:, 1:rhoc_i4_nu(i), i)
                    write(21) vel_i2_nu(:, 1:rhoc_i4_nu(i), i)
#else
                    write(22) pos_i1_nu(:, 1:rhoc_i4_nu(i), i)
#endif
#ifdef NUPID
                    write(25) pid_i8_nu(1:rhoc_i4_nu(i), i)
#endif
                enddo
            enddo !! j
        enddo !! k

        close(12)
#ifdef ZIPDM
        close(10); close(11); close(13)
#endif
        close(22)
#ifdef ZIP
        close(20); close(21); close(23)
#endif

    else
#endif

    !! Regular format

        do i=1,num_writes
          nplow=(i-1)*blocksize+1
          nphigh=min(i*blocksize,np_local)
          do j=nplow,nphigh
#ifdef NUPID
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
            else
#ifdef DISP_MESH
                write(22) xv(1:3,j) - shake_offset
                write(22) xv(4:6,j)
#else
                write(22) xv(:,j)
#endif
#ifdef NUPID
                write(25) PID(j)
#endif
                ind_check2 = ind_check2 + 1
            endif
          enddo
        enddo

        close(12)
        close(22)

#ifdef ZIP
    endif
#endif

#ifdef NUPID
    !! Close neutrino PID file
    close(25)
#endif

    !! Consistency check
    if (ind_check1 .ne. np_dm .or. ind_check2 .ne. np_nu) then
        write(*,*) "Dark Matter checkpoint error: ind_checks ", ind_check1, np_dm, ind_check2, np_nu
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#else

    !
    ! Write data for non-neutrino simulation
    !

#ifdef ZIPDM
    if (dozip) then

        !! Zipped format

#ifdef PID_FLAG
    !! Open PID file and write header
    open(unit=15, file=ofile2, status="replace", iostat=fstat, access="stream", buffered="yes")
    if (fstat /= 0) then
      write(*,*) 'error opening PID file for write'
      write(*,*) 'rank',rank,'file:',ofile2
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    write(15) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p
#endif

        do k = 1, nc_node_dim
            do j = 1, nc_node_dim
                !$omp parallel default(shared) private(i,l,ind)
                !$omp do schedule(static)
                do i = 1, nc_node_dim
                    ind = 1
                    rhoc_i4(i) = 0
                    l = hoc(i,j,k)
                    do while (l > 0)
                        rhoc_i4(i) = rhoc_i4(i) + 1
                        pos_i1(:, ind, i) = int(mod(xv(1:3,l)/mesh_scale,1.)*256, kind=1)
                        vel_i2(:, ind, i) = nint(xv(4:6,l)*v_r2i, kind=2)
#ifdef PID_FLAG
                        pid_i8(ind, i) = PID(l)
#endif
                        ind = ind + 1
                        l = ll(l)
                    enddo !! l
                    if (ind > max_cell_np+1) then
                        write(*,*) "ERROR: max_cell_np too small !!", rank, ind, max_cell_np
                        call mpi_abort(mpi_comm_world, ierr, ierr)
                    endif
                enddo !! i
                !$omp end do
                !$omp end parallel
                write(12) int(min(rhoc_i4, 255), kind=1)
                write(13) pack(rhoc_i4, rhoc_i4>254)
                do i = 1, nc_node_dim
                    write(10) pos_i1(:, 1:rhoc_i4(i), i)
                    write(11) vel_i2(:, 1:rhoc_i4(i), i)
#ifdef PID_FLAG
                    write(15) pid_i8(1:rhoc_i4(i), i)
#endif
                enddo
            enddo !! j
        enddo !! k

        close(10); close(11); close(12); close(13)
#ifdef PID_FLAG
        close(15)
#endif

    else
#endif

#ifdef DISP_MESH
    do j = 1, np_local
        xv(1:3,j) = xv(1:3,j) - shake_offset
    enddo
#endif

        do i=1,num_writes
          nplow=(i-1)*blocksize+1
          nphigh=min(i*blocksize,np_local)
          write(12) xv(:,nplow:nphigh)
        enddo

        close(12)

#ifdef DISP_MESH
    do j = 1, np_local
        xv(1:3,j) = xv(1:3,j) + shake_offset
    enddo
#endif

#ifdef ZIPDM
    endif
#endif

#endif

    !
    ! Write dark matter IDs (if not a neutrino simulation)
    !

#ifndef NEUTRINOS
#ifdef PID_FLAG
#ifdef ZIPDM
    if (.not. dozip) then
#endif
    open(unit=15, file=ofile2, status="replace", iostat=fstat, access="stream")

    if (fstat /= 0) then
      write(*,*) 'error opening PID file for write'
      write(*,*) 'rank',rank,'file:',ofile2
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    !! This is the file header

    write(15) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
              cur_projection,cur_halofind,mass_p

    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
      write(15) PID(nplow:nphigh)
    enddo

    close(15)
#ifdef ZIPDM
    endif
#endif
#endif    
#endif

    !
    ! Write gas checkpoint
    !

#ifdef MHD
    call mpi_tvd_mhd_state_output(output_path,nts,t,z_s)
#endif

    write(*,*) 'Finished checkpoint_kill:',rank

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
#    include "cubepm.fh"

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

