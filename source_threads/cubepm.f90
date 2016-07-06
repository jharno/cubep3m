!!  cubep3m - cubical decomposition 2-level particle mesh algorithm with particle-particle interactions
!! Hugh Merz :: merz@cita.utoronto.ca :: 2006 11 02 
program cubep3m
  use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif

#ifdef MHD
  use mpi_tvd_mhd
#endif

  implicit none
  include 'mpif.h'
#  include "cubepm.fh"

#ifdef MHD
  integer(4) :: nc_to_mhd(3)
#endif

#ifdef WRITELOG
  integer(4) :: fstat, np_max
#endif

  real(4) :: t_elapsed
  external t_elapsed

  real(8) :: sec1, sec2, sec01, sec02, seconds2wait
  real(8) :: sec1a, sec2a
  logical(kind=4) :: i_continue

#ifdef CHECKPOINT_KILL
  logical :: kill_step, kill_step_done
  kill_step_done = .false.
#endif

  call mpi_initialize

  if (rank == 0) call datestamp

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STARTING CUBEP3M: ", sec1

#ifdef CHECKPOINT_KILL
  call read_remaining_time
#endif

  call t_start(wc_counter)

  call memory_usage

  call variable_initialize

  if (rank == 0) write(*,*) 'finished variable init',t_elapsed(wc_counter)

#ifdef NESTED_OMP
    call omp_set_num_threads(cores*nested_threads)
    call omp_set_nested(.true.)
#else
    call omp_set_num_threads(cores)
#endif

if (rank == 0) write(*,*) 'finished omp call',t_elapsed(wc_counter)

  call coarse_kernel
if (rank == 0) write(*,*) 'finished coarse kernel init',t_elapsed(wc_counter)

  call fine_kernel
if (rank == 0) write(*,*) 'finished kernel init',t_elapsed(wc_counter)

#ifdef KERN_DUMP
  call kernel_checkpoint(.true.)
#endif

#ifdef MEMORY_WAIT
if (rank==0) then
  seconds2wait=0
  print*, 'before particle_initialize, wait for',seconds2wait,' seconds'
  sec01=mpi_wtime(ierr)
  i_continue=.false.
  do while (i_continue == .false.)
    sec02 = mpi_wtime(ierr) - sec01
    if (sec02 > seconds2wait) i_continue=.true.
  enddo
endif
call mpi_barrier(mpi_comm_world,ierr)
#endif

  call particle_initialize

  if (rank == 0) write(*,*) 'finished initializing particles',t_elapsed(wc_counter)

#ifdef MHD
  nc_to_mhd=nf_physical_node_dim
  call mpi_tvd_mhd_init(nc_to_mhd,mpi_comm_cart,cart_rank, &
                        cart_coords,cart_neighbor,nodes_dim)
  if (rank == 0) write(*,*) 'finished initializing mhd',t_elapsed(wc_counter)

  if (grid_ic) then
    u=0.0
    u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)=1.0
    u(5,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)=0.0001
    print *,rank,sum(u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
    b=0.0
  elseif (random_ic) then
    u=0.0
    u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)=1.0
    u(5,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)=0.0001
    print *,rank,sum(u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
    b=0.0
  else
    call mpi_tvd_mhd_ic(ic_path)
  endif

  if (rank == 0) write(*,*) 'finished reading mhd initial conditions',t_elapsed(wc_counter)

  call comm_bufferupdate(u,b,nx,ny,nz)
  call transposef(u,b,nx%l,ny%l,nz%l)
  call comm_bufferupdate(u,b,ny,nz,nx)
  call transposef(u,b,ny%l,nz%l,nx%l)
  call comm_bufferupdate(u,b,nz,nx,ny)
  call transposef(u,b,nz%l,nx%l,ny%l)

  if (rank == 0) write(*,*) 'finished updating buffers with initial conditions',t_elapsed(wc_counter)
#endif

 !if(rank ==0) write(*,*)'Calling init_projection.f90'
 call link_list
 call init_projection


  if (rank == 0) write(*,*) 'starting main loop'

#ifdef WRITELOG
  if (rank==0) then
    fstat=0
    open(unit=76,file=logfile,status='replace',iostat=fstat,form='formatted')
  endif
#endif

  do 
    call timestep
    sec1a = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "TIMESTEP_TIME [hrs] = ", (sec1a - sec1) / 3600.


#ifdef WRITELOG
    if(rank==0) then
       write(unit=76,fmt='(i6,2x)',    advance='no' )  nts
       write(unit=76,fmt='(f10.6,2x)', advance='no' )  1.0/a-1.0
       write(unit=76,fmt='(f10.6,2x)', advance='no' )  (sec1a-sec1)/3600.
       write(unit=76,fmt='(f10.6)',    advance='yes')  min_den_buf
    endif
#endif

#ifdef MHD

! Note that the dimensions of u and b don't match the
! declaration in this routine on return from transpose.
! u and b simply act as pointers and the dimensions are
! correctly interpreted in sweep.

! first pass gas update

! Forward Sweep
    call sweep(forward,u,b,nx,ny,nz,dt_gas)  !x
    call transposef(u,b,nx%l,ny%l,nz%l)
    call sweep(forward,u,b,ny,nz,nx,dt_gas)  !y
    call transposef(u,b,ny%l,nz%l,nx%l)
    call sweep(forward,u,b,nz,nx,ny,dt_gas)  !z

    if (rank == 0) write(*,*) 'finished forward gas sweep',t_elapsed(wc_counter)

! Backward Sweep
    call sweep(backward,u,b,nz,nx,ny,dt_gas) !z
    call transposeb(u,b,nz%l,nx%l,ny%l)
    call sweep(backward,u,b,ny,nz,nx,dt_gas) !y
    call transposeb(u,b,ny%l,nz%l,nx%l)
    call sweep(backward,u,b,nx,ny,nz,dt_gas) !x

    if (rank == 0) write(*,*) 'finished backward gas sweep',t_elapsed(wc_counter)

#ifdef DEBUG_DEN_BUF
    print *,rank,sum(u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
#endif

    call gas_density_buffer

    if (rank == 0) write(*,*) 'finished gas density buffer',t_elapsed(wc_counter)

#endif

    call particle_mesh
    if (rank == 0) write(*,*) 'finished particle mesh',t_elapsed(wc_counter)

#ifdef MHD
! second pass gas update
! Forward Sweep
    call sweep(forward,u,b,nx,ny,nz,dt_gas)  !x
    call transposef(u,b,nx%l,ny%l,nz%l)
    call sweep(forward,u,b,ny,nz,nx,dt_gas)  !y
    call transposef(u,b,ny%l,nz%l,nx%l)
    call sweep(forward,u,b,nz,nx,ny,dt_gas)  !z

    if (rank == 0) write(*,*) 'finished forward gas sweep',t_elapsed(wc_counter)

! Backward Sweep
    call sweep(backward,u,b,nz,nx,ny,dt_gas) !z
    call transposeb(u,b,nz%l,nx%l,ny%l)
    call sweep(backward,u,b,ny,nz,nx,dt_gas) !y
    call transposeb(u,b,ny%l,nz%l,nx%l)
    call sweep(backward,u,b,nx,ny,nz,dt_gas) !x

    if (rank == 0) write(*,*) 'finished backward gas sweep',t_elapsed(wc_counter)
#endif

#ifdef ZOOMCHECK
    call zoomcheckpoint
#endif
#ifdef BOXDUMP
    call boxcheckpoint
#endif

#ifdef CHECKPOINT_KILL
    !! Determine if it is time to write a checkpoint before being killed
    kill_step = .false.
    sec1a = mpi_wtime(ierr)
    if (rank == 0) then
        if ((sec1a - sec1) .ge. kill_time) kill_step = .true.
    endif
    call mpi_bcast(kill_step, 1, mpi_logical, 0, mpi_comm_world, ierr)

    if (checkpoint_step.or.projection_step.or.halofind_step.or.kill_step) then
#else
    if (checkpoint_step.or.projection_step.or.halofind_step) then
#endif

!! advance the particles to the end of the current step.

      dt_old = 0.0
      call update_position

#if defined(ZIP) || defined(ZIPDM)
    force_grid_back = .true.
    call move_grid_back
    force_grid_back = .false.
#endif
    call link_list

#ifdef CHECKPOINT_KILL
      if (kill_step .eqv. .true. .and. kill_step_done .eqv. .false.) then
        sec1a = mpi_wtime(ierr)
        if (rank == 0) write(*,*) "STARTING CHECKPOINT_KILL: ", sec1a
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.true.)
#else
        call checkpoint_kill
#endif
        sec2a = mpi_wtime(ierr)
        if (rank == 0) write(*,*) "STOPPING CHECKPOINT_KILL: ", sec2a
        if (rank == 0) write(*,*) "ELAPSED CHECKPOINT_KILL TIME: ", sec2a-sec1a
        kill_step_done = .true. ! Don't want to keep doing this
      endif
#endif

      if (checkpoint_step) then
        sec1a = mpi_wtime(ierr)
        if (rank == 0) write(*,*) "STARTING CHECKPOINT: ", sec1a
        call checkpoint
        if (rank == 0) write(*,*) 'finished checkpoint',t_elapsed(wc_counter)
        sec2a = mpi_wtime(ierr)
        if (rank == 0) write(*,*) "STOPPING CHECKPOINT: ", sec2a
        if (rank == 0) write(*,*) "ELAPSED CHECKPOINT TIME: ", sec2a-sec1a
      endif

      if (projection_step.or.halofind_step) then
        call particle_pass
        if (halofind_step) then
          sec1a = mpi_wtime(ierr)
          if (rank == 0) write(*,*) "STARTING HALOFIND: ", sec1a
          call halofind
          if (rank == 0) write(*,*) 'finished halofind',t_elapsed(wc_counter)
          sec2a = mpi_wtime(ierr)
          if (rank == 0) write(*,*) "STOPPING HALOFIND: ", sec2a
          if (rank == 0) write(*,*) "ELAPSED HALOFIND TIME: ", sec2a-sec1a
        endif

        if (projection_step) then
          call projection
          if (rank == 0) write(*,*) 'finished projection',t_elapsed(wc_counter)
        endif

        if(superposition_test)then
           !fine_clumping=0.0
           !call link_list
           !call particle_pass
           !call halofind
           if(rank==0)write(*,*) 'Calling report force'
           call report_force
           if(rank==0)then
              write(*,*) 'Called report force'       
              write(*,*) '*** Ending simulation here ***'
           endif
              !stop
              !call  mpi_finalize(ierr)
              !exit           
           stop
        endif

        !! Clean up ghost particles
        call delete_particles

      endif

      dt = 0.0

    endif

    if (nts == max_nts .or. final_step .or. a .gt. 1.0) exit

  enddo


#ifdef WRITELOG
    if(rank==0) close(76)
#endif


#ifdef TIMING
  if (rank==0) then
    print *,'cubep3m finished:' 
    call datestamp
  endif
#endif


#ifdef MHD
  call mpi_tvd_mhd_finalize
#endif

  call cubepm_fftw(0)
  do ierr=1,cores 
    call cubepm_fftw2('q',ierr)
  enddo

  sec2 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STOPPING CUBEP3M: ", sec2
  if (rank == 0) write(*,*) "ELAPSED CUBEP3M TIME: ", sec2-sec1

  call mpi_finalize(ierr)

  if (rank == 0) call datestamp

contains
  
  subroutine memory_usage
!! calculate and display common block memory usage
    implicit none

    real(4) :: memory_used


    memory_used = real((nc_dim+2)*(nc_dim)*(nc_slab)) &    !slab
       + 3.0*real((nc_dim/2+1)*nc_dim*nc_slab) &  !kern_c
       + real((nf_tile+2)*nf_tile**2) & !rho_f
       + real(2*nf_tile*(nf_tile/2+1)*nf_tile) & !rho_ft
       + real(3*(nf_tile-2*nf_buf+2)**3) &             !force_f
       + real(nf_tile*(nf_tile/2+1)*nf_tile) & !kern_f
       + real(max_buf) & !send_buf
       + real(max_buf) & !recv_buf
       + 6.0*real(max_np) & !xv
       + real(max_np) &	 !ll
       + real((hoc_nc_h-hoc_nc_l)**3) !hoc


#ifdef PID_FLAG
!!  should also add PIDsend/PIDrecv buffers, but they're small 
    memory_used = memory_used + 2.0*real(max_np) !PID
#endif



#ifdef MHD
!! should also add send/recv mhd buffers, but they're small 
    memory_used = memory_used + 8*real(nf_physical_node_dim+6)**3 !u+b
#endif

    if (rank==0) then
      write(*,*) 'majority of memory used / node :'
      write(*,*) memory_used*4.0,'bytes'
      write(*,*) memory_used*4.0/1024/1024/1024,'GB' 
    endif

  end subroutine memory_usage

end program cubep3m
