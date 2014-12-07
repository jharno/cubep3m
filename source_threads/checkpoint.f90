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
#ifdef NUPID
    integer, parameter :: pdm = 0
#else
    integer, parameter :: pdm = 1
#endif
#endif
#ifdef ZIP
    integer :: l, k
    integer :: fstat1, fstat2, fstat3
    character (len=max_path) :: fdm_zip1,fdm_zip2,fdm_zip3
    character (len=max_path) :: fnu_zip1,fnu_zip2,fnu_zip3
    integer(4) :: v_resolution = 16384
    real(4) :: vmax, vmax_local, v_r2i
    integer(4) :: rhoc_dm_i4, rhoc_nu_i4
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

    !
    ! Checkpoint file names
    !

    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)

    write(z_s,'(f7.3)') z_write
    z_s=adjustl(z_s)

    !! Dark matter file
#ifndef ZIP
    ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#else
    fdm_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'.dat'
    fdm_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'.dat'
    fdm_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'.dat'
#endif

#ifdef NEUTRINOS
    !! Neutrino file
#ifndef ZIP
    ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
         rank_s(1:len_trim(rank_s))//'_nu.dat'
#else
    fnu_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
    fnu_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
    fnu_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif
    !! Neutrino PID file
#ifdef NUPID
    ofile2_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PID'// &
         rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif
#else
    !! Dark matter PID file (if not a neutrino sim)
#ifdef PID_FLAG
    ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'PID'// &
         rank_s(1:len_trim(rank_s))//'.dat'
#endif
#endif

#ifdef ZIP
    !
    ! Determine velocity conversion factor
    !

    vmax_local = maxval(abs(xv(4:6,1:np_local)))
    call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
    v_r2i = v_resolution/vmax
#endif

    !
    ! Increment checkpoint counter so restart works on next checkpoint
    !
    
    cur_checkpoint=cur_checkpoint+1
    cur_proj = cur_projection
    cur_halo = cur_halofind
    if (projection_step) cur_proj = cur_proj + 1
    if (halofind_step) cur_halo = cur_halo + 1

    !
    ! Open checkpoint files
    !

    !! Dark matter 
#ifndef ZIP
    open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#else
    open(unit=11, file=fdm_zip1, status="replace", iostat=fstat1, access="stream", buffered='yes')
    open(unit=12, file=fdm_zip2, status="replace", iostat=fstat2, access="stream", buffered='yes')
    open(unit=13, file=fdm_zip3, status="replace", iostat=fstat3, access="stream", buffered='yes')
    if (fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
      write(*,*) 'error opening zipped checkpoint file for write'
      write(*,*) 'rank',rank,'file:',fdm_zip1,fdm_zip2,fdm_zip3
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#endif

#ifdef NEUTRINOS

    !! Neutrinos 
#ifndef ZIP
    open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream") 
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint file for write'
      write(*,*) 'rank',rank,'file:',ofile_nu
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#else
    open(unit=21, file=fnu_zip1, status="replace", iostat=fstat1, access="stream", buffered='yes')
    open(unit=22, file=fnu_zip2, status="replace", iostat=fstat2, access="stream", buffered='yes')
    open(unit=23, file=fnu_zip3, status="replace", iostat=fstat3, access="stream", buffered='yes')
    if (fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
      write(*,*) 'error opening zipped checkpoint file for write'
      write(*,*) 'rank',rank,'file:',fnu_zip1,fnu_zip2,fnu_zip3
      call mpi_abort(mpi_comm_world,ierr,ierr)
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

#ifndef ZIP
    write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
    write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
#else
    write(11) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p,v_r2i,shake_offset
    write(21) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p,v_r2i,shake_offset
#endif

#ifdef NUPID
    write(25) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
#endif

#else
#ifndef ZIP
    write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
#else
    write(11) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p,v_r2i,shake_offset
#endif
#endif

#ifdef NEUTRINOS

    !
    ! Write data for neutrino simulation
    !

    ind_check1 = 0
    ind_check2 = 0

#ifdef ZIP

    !! Zipped format

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim
            do i = 1, nc_node_dim

                rhoc_dm_i4 = 0 ; rhoc_nu_i4 = 0

                l = hoc(i,j,k)
                do while (l > 0)
                    if (PID(l) == pdm) then
                        rhoc_dm_i4 = rhoc_dm_i4 + 1
                        write(11) int(mod(xv(1:3,l)/mesh_scale,1.)*256, kind=1)
                        write(11) int(xv(4:6,l)*v_r2i, kind=2)
                        ind_check1 = ind_check1 + 1
                    else
                        rhoc_nu_i4 = rhoc_nu_i4 + 1
                        write(21) int(mod(xv(1:3,l)/mesh_scale,1.)*256, kind=1)
                        write(21) int(xv(4:6,l)*v_r2i, kind=2)
#ifdef NUPID
                        write(25) PID(l)
#endif
                        ind_check2 = ind_check2 + 1
                    endif
                    l = ll(l)
                enddo !! l

                if (rhoc_dm_i4 < 255) then
                    write(12) int(rhoc_dm_i4, kind=1) ! write density in int1
                else
                    write(12) int(255, kind=1)
                    write(13) rhoc_dm_i4
                endif
                if (rhoc_nu_i4 < 255) then
                    write(22) int(rhoc_nu_i4, kind=1)
                else
                    write(22) int(255, kind=1)
                    write(23) rhoc_nu_i4
                endif

            enddo !! i
        enddo !! j
    enddo !! k

    close(11); close(12); close(13)
    close(21); close(22); close(23)

#else

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

#ifdef ZIP

    !! Zipped format

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim
            do i = 1, nc_node_dim

                rhoc_dm_i4 = 0

                l = hoc(i,j,k)
                do while (l > 0)
                    rhoc_dm_i4 = rhoc_dm_i4 + 1
                    write(11) int(mod(xv(1:3,l)/mesh_scale,1.)*256, kind=1)
                    write(11) int(xv(4:6,l)*v_r2i, kind=2)
                    l = ll(l)
                enddo !! l

                if (rhoc_dm_i4 < 255) then
                    write(12) int(rhoc_dm_i4, kind=1) ! write density in int1
                else
                    write(12) int(255, kind=1)
                    write(13) rhoc_dm_i4
                endif

            enddo !! i
        enddo !! j
    enddo !! k

    close(11); close(12); close(13)

#else

    !! Regular format

    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
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

#endif

    !
    ! Write dark matter IDs (if not a neutrino simulation)
    !

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

    do i=1,num_writes
      nplow=(i-1)*blocksize+1
      nphigh=min(i*blocksize,np_local)
      do j=nplow,nphigh
        write(15) PID(j)
      enddo
    enddo

    close(15)
#endif
#endif    

    !
    ! Write gas checkpoint
    !

#ifdef MHD
    call mpi_tvd_mhd_state_output(output_path,nts,t,z_s)
#endif

    !
    ! Print info to screen
    !

    if (rank==0) then
        print*, 'current steps recorded in xv file:'
        print*, 'cur_checkpoint =', cur_checkpoint
        print*, 'cur_projection =', cur_proj
        print*, 'cur_halofind   =', cur_halo
    endif

    write(*,*) 'Finished checkpoint:',rank

    checkpoint_step=.false.

  end subroutine checkpoint
