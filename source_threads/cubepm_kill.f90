!!  cubep3m - cubical decomposition 2-level particle mesh algorithm with particle-particle interactions
!! Hugh Merz :: merz@cita.utoronto.ca :: 2006 11 02 
!! Joachim Harnois-Deraps :: jharno@cita.utoronto.ca :: 2013 01 10
!! Now has particle ID, extended pp force available, and MHD has been awaken
!! from the grave, and is now compatible with the threaded version.
program cubep3m
  use omp_lib

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

#ifdef MHD
  integer(4) :: nc_to_mhd(3),i,j,k
  real(4)    :: total_baryon, local_baryon
#endif

  real(4) :: t_elapsed
  external t_elapsed
  real(4) :: z_write
  character(len=7) :: z_s

#ifdef CHECKPOINT_KILL
  logical :: kill_step, kill_step_done
  ! Time in seconds to write checkpoint before being killed:
  !real, parameter :: kill_time = 5.*60               ! 5 minutes
  real, parameter :: kill_time = (48. - 2.) * 3600. ! 46 hours
  real(8) :: sec1, sec2
#endif

  call datestamp 

  call mpi_initialize

  call t_start(wc_counter)

  call memory_usage

  call variable_initialize

  
#ifdef CHECKPOINT_KILL
  sec1 = mpi_wtime(ierr)
  kill_step_done = .false.
#endif

  if (rank == 0) write(*,*) 'finished variable init',t_elapsed(wc_counter)

  !$ call omp_set_num_threads(cores)

if (rank == 0) write(*,*) 'finished omp call',t_elapsed(wc_counter)

  call coarse_kernel
if (rank == 0) write(*,*) 'finished coarse kernel init',t_elapsed(wc_counter)

  call fine_kernel
if (rank == 0) write(*,*) 'finished kernel init',t_elapsed(wc_counter)

#ifdef KERN_DUMP
  call kernel_checkpoint(.true.)
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
  elseif (restart_ic)then
    if (rank==0) z_write = z_checkpoint(restart_checkpoint)
    call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

    write(z_s,'(f7.3)') z_write
    z_s=adjustl(z_s)
    call mpi_tvd_mhd_restart(output_path, z_s)
  else
    call mpi_tvd_mhd_ic(ic_path)
  endif

  if (rank == 0) write(*,*) 'finished reading mhd initial conditions',t_elapsed(wc_counter)
  local_baryon = sum(u(1,:,:,:))
  call mpi_reduce(local_baryon,total_baryon,1,MPI_REAL, mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'total baryon mass = ', total_baryon
  call mpi_bcast(total_baryon,1,MPI_REAL,0,mpi_comm_world,ierr)

#ifdef MHD_DEBUG
  write(*,*) rank, 'mean local mass = ' , sum(u(1,:,:,:))/nf_physical_node_dim**3
  write(*,*) rank, 'min/max den = ', minval( u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)), maxval(u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
  write(*,*) rank, 'mean local den*vx = ' , sum(u(2,:,:,:))/nf_physical_node_dim**3
  write(*,*) rank, 'min/max den*vx = ', minval(u(2,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)), maxval(u(2,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)) 
  write(*,*) rank, 'mean local den*vy = ' , sum(u(3,:,:,:))/nf_physical_node_dim**3
  write(*,*) rank, 'min/max den*vy = ', minval( u(3,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)), maxval( u(3,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
  write(*,*) rank, 'mean local den*vz = ' , sum(u(4,:,:,:))/nf_physical_node_dim**3
  write(*,*) rank, 'min/max den*vz = ', minval( u(4,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)), maxval( u(4,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
  write(*,*) rank, 'mean local E = ' , sum(u(5,:,:,:))/nf_physical_node_dim**3
  write(*,*) rank, 'min/max E = ', minval( u(5,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)), maxval(u(5,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
  do k=nz%m,nz%n
    do j=ny%m,ny%n
      do i=nx%m,nx%n
         if(abs(u(3,i,j,k)) > 20) write(*,*) 'velocity too high at ' ,rank, i,j,k,u(3,i,j,k)
  
      enddo
    enddo
  enddo
#endif

  call comm_bufferupdate(u,b,nx,ny,nz)
  call transposef(u,b,nx%l,ny%l,nz%l)
  call comm_bufferupdate(u,b,ny,nz,nx)
  call transposef(u,b,ny%l,nz%l,nx%l)
  call comm_bufferupdate(u,b,nz,nx,ny)
  call transposef(u,b,nz%l,nx%l,ny%l)

  if (rank == 0) write(*,*) 'finished updating buffers with initial conditions',t_elapsed(wc_counter)
#endif

#ifdef init_proj_only
  if(rank ==0) write(*,*)'Calling init_projection.f90'
    call particle_pass
    call link_list
    call init_projection
    call cubepm_fftw(0)
  do ierr=1,cores
    call cubepm_fftw2('q',ierr)
  enddo
  call mpi_finalize(ierr)
  call datestamp
  if(rank==0) write(*,*)'Quiting after projecting'
#endif

#ifdef proj_only_restart
  if (rank == 0) write(*,*) 'Entering finer projection then quit'
  call link_list
  call particle_pass
  call projection
  call mpi_finalize(ierr)
  call datestamp
#endif

  if (rank == 0) write(*,*) 'starting main loop'
  do 
    call timestep

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

#ifdef CHECKPOINT_KILL
    !! Determine if it is time to write a checkpoint before being killed
    kill_step = .false.
    sec2 = mpi_wtime(ierr)
    if (rank == 0) then
        if ((sec2 - sec1) .ge. kill_time) kill_step = .true.
    endif
    call mpi_bcast(kill_step, 1, mpi_logical, 0, mpi_comm_world, ierr)

    if (checkpoint_step.or.projection_step.or.halofind_step.or.kill_step) then
#else
    if (checkpoint_step.or.projection_step.or.halofind_step) then
#endif

!! advance the particles to the end of the current step.

      dt_old = 0.0
      call update_position

#ifdef MOVE_GRID_BACK
      call move_grid_back
#endif

#ifdef CHECKPOINT_KILL
      if (kill_step .and. .not. kill_step_done) then

        call checkpoint_kill
        if (rank == 0) write(*,*) 'finished checkpoint_kill',t_elapsed(wc_counter)

        kill_step_done = .true. ! Don't want to keep doing this

      endif
#endif

      !dt = 0.0

      if (checkpoint_step) then
        call checkpoint
        if (rank == 0) write(*,*) 'finished checkpoint',t_elapsed(wc_counter)
      endif

      if (projection_step.or.halofind_step) then

!! Update ll and pass particles

        call link_list
        call particle_pass

        fine_clumping=0.0
        if (halofind_step) then
          call halofind
          if (rank == 0) write(*,*) 'finished halofind',t_elapsed(wc_counter)
        endif

        call mpi_barrier(mpi_comm_world,ierr)

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

    if (nts == max_nts .or. final_step .or. (a .gt. 1.0)) exit
  enddo

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
  call mpi_finalize(ierr)

  call datestamp

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
