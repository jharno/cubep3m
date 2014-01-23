!! initialize variables
  subroutine variable_initialize
    implicit none

    include 'mpif.h'
    include 'cubepm.fh'

    integer(4) :: i,fstat

    ierr=-1
    fstat=0

!! initial scalefactor, see cubepm.par

    a=a_i
    tau=-3.0/sqrt(a_i)

!! set these variables prohibitively high such that they don't 
!! bugger up the first half-step update 
!! (1st sweep starts with a half-step position update)

    dt_f_acc=1000.0
    dt_c_acc=1000.0
#ifdef PPINT
    dt_pp_acc=1000.0
#ifdef PP_EXT
    dt_pp_ext_acc=1000.0
#endif
#endif

!! zero everything else 

    dt_old=0.0
    dt=0.0
    da=0.0
    nts=0
    t=0.0
    np_local=0
    mass_p=0.0
    if (pair_infall) then
      cur_sep=pair_infall_sep
    else
      cur_sep=min_sep 
    endif

    rho_f=0.0
    kern_f=0.0
    force_f=0.0
    kern_c=0.0
    slab=0.0
    send_buf=0.0
    recv_buf=0.0
    xv=0.0
    ll=0
    hoc=0
    np_buf=0
    final_step=.false.
    shake_offset=0.0
#ifdef PID_FLAG 
    PID=0
    send_buf_PID=0.0
    recv_buf_PID=0.0
#endif
    
    if (rank == 0) then

!! read in when to store projections

      open(11,file=projections, status='old', iostat=fstat)
      if (fstat /= 0) then
        write(*,*) 'error opening projection list file'
        write(*,*) 'rank',rank,'file:',projections
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif 
      do num_projections=1,max_input
        read(unit=11,err=61,end=71,fmt='(f20.10)') z_projection(num_projections)
        print*, 'projections ',z_projection(num_projections)
      enddo
    71  num_projections=num_projections-1
    61  close(11)

      if (num_projections.eq.max_input) then
        write(*,*) 'too many projections to store > ',max_input
        call mpi_abort(mpi_comm_world,i,ierr)
      endif

      print*,num_projections,' projections to do'
      if (num_projections.eq.1) then
        write(*,*) 'problem reading projections '
        call mpi_abort(mpi_comm_world,i,ierr)
      endif


      do i=1,num_projections
        a_projection(i)=1.0/(1.0+z_projection(i))
      enddo

      if (.not. pairwise_ic) then
        if (num_projections > 0) then
          write(*,*) 'Projections performed at:'
          write(*,*) 'z        a'
          do i=1,num_projections
            write(*,'(f8.4,2x,f6.4)') z_projection(i),a_projection(i)
          enddo
        else
          a_projection(1)=100.0
          write(*,*) 'no projections to be stored'
        endif
      endif

!! read in when to store checkpoints

      open(11,file=checkpoints,status='old',iostat=fstat)
      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint list file'
        write(*,*) 'rank',rank,'file:',checkpoints
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif 
      do num_checkpoints=1,max_input
        read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
        print*,'checkoints', z_checkpoint(num_checkpoints)
      enddo
    41  num_checkpoints=num_checkpoints-1
    51  close(11)

      print*,num_checkpoints,' checkpoints to do'
      if (num_checkpoints.eq.1) then
        write(*,*) 'problem reading checkpoints '
        call mpi_abort(mpi_comm_world,i,ierr)
      endif

!! temporary rig

!      z_checkpoint=z_projection
!      num_checkpoints=num_projections

      if (num_checkpoints.eq.max_input) then
        write(*,*) 'too many checkpoints to store > ',max_input
        call mpi_abort(mpi_comm_world,i,ierr)
      endif

      if (z_checkpoint(1).gt.z_i) then
        write(*,*) 'z_initial less than first checkpoint, exiting',z_i,z_checkpoint(1)
        call mpi_abort(mpi_comm_world,i,ierr)
      endif

      do i=1,num_checkpoints
        a_checkpoint(i)=1.0/(1.0+z_checkpoint(i))
      enddo

      if ( pairwise_ic ) then
        write(*,*) '------------------------------'
        write(*,*) 'Conducting Pairwise Force Test'
        write(*,*) '------------------------------'
      else
        write(*,*) 'Starting simulation at:'
        write(*,*) 'z      a'
        write(*,'(f8.4,2x,f6.4)') z_i,a_i
        if (num_checkpoints > 0) then
          write(*,*) 'Checkpointing performed at:'
          do i=1,num_checkpoints
            write(*,'(f8.4,2x,f6.4)') z_checkpoint(i),a_checkpoint(i)
          enddo
        else
          write(*,*) 'no checkpoints to be stored'
        endif
      endif

!! read in when to store halo catalogs

      open(11,file=halofinds, status='old', iostat=fstat)
      if (fstat /= 0) then
        write(*,*) 'error opening halo catalog list file'
        write(*,*) 'rank',rank,'file:',halofinds
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      do num_halofinds=1,max_input
        read(unit=11,err=81,end=91,fmt='(f20.10)') z_halofind(num_halofinds)
        print*,'check halofinds ', z_halofind(num_halofinds)
      enddo
    91  num_halofinds=num_halofinds-1
    81  close(11)

      if (num_halofinds.eq.max_input) then
        write(*,*) 'too many halo catalogs to store > ',max_input
        call mpi_abort(mpi_comm_world,i,ierr)
      endif

        print*,num_halofinds,' halofinds to do'
        if (num_halofinds.eq.1) then
        write(*,*) 'problem reading halofinds '
        call mpi_abort(mpi_comm_world,i,ierr)
      endif


      do i=1,num_halofinds
        a_halofind(i)=1.0/(1.0+z_halofind(i))
      enddo

      if (.not. pairwise_ic) then
        if (num_halofinds > 0) then
          write(*,*) 'Halo catalogs generated at:'
          write(*,*) 'z        a'
          do i=1,num_halofinds
            write(*,'(f7.4,2x,f6.4)') z_halofind(i),a_halofind(i)
          enddo
        else
          a_halofind(1)=100.0
          write(*,*) 'no halo catalogs to be stored'
        endif
      endif
   endif

    cur_halofind=1
    cur_projection=1
    cur_checkpoint=1

! Initialize first time through fftw flag

    firstfftw=.true.
    firstfftw2=.true.

    ! NEW TRICK BY JHD TO REMOVE MEMORY RACING CONDITION ON THREADED PLAN
    ! CREATION. 
    do i = 1,cores 
      call cubepm_fftw2('o',i)
    enddo

! Initialize halo finding arrays

    call initialize_halofind

  end subroutine variable_initialize
