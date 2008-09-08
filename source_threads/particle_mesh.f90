!! main particle mesh subroutine
  subroutine particle_mesh
    implicit none
    include 'mpif.h'
    include 'cubepm.fh'

    integer(4) :: i,j,k
    integer(4), dimension(3) :: tile
    integer(4) :: thread
    real(4) :: f_force_max_node
    real(4) :: pp_force_max_node
    integer(4) :: omp_get_thread_num
    external omp_get_thread_num
 
#ifdef MHD
    integer(4) :: nerrl,nerr
    real(4) :: cmaxl,cmax
#endif

#ifdef DIAG
    real(8) :: sumrhof
#endif

!! start of particle mesh.  All particles are within (1:nc_node_dim]

    if (pairwise_ic) then
      call set_pair
    elseif (pair_infall) then
      if (nts==1) then
        call set_pair
      endif
      call update_position
    else
      call update_position
    endif

!! particles must not have advanced past hoc_nc_l:hoc_nc_h

    call link_list

    call particle_pass

#ifdef MHD
    nerr=0
    cmax=1e-5
#endif

    !$omp parallel default(shared) &
    !$omp private(i,j,tile,thread)
    thread=1
    !$ thread=omp_get_thread_num()+1
    f_mesh_mass(thread)=0.0
    f_force_max(thread)=0.0
#ifdef PPINT
    pp_force_max(thread)=0.0
#endif
    !$omp do 
    do i=1,tiles_node
      tile(3) = (i-1) / (tiles_node_dim * tiles_node_dim)
      j = i - tile(3) * tiles_node_dim * tiles_node_dim
      tile(2) = (j-1) /  tiles_node_dim
      j = j - tile(2) * tiles_node_dim
      tile(1) = j - 1
#ifdef MHD
      call fine_mesh(tile,cmax,nerr,thread)
#else
      call fine_mesh(tile,thread)
#endif
    enddo
    !$omp end do
    !$omp end parallel

#ifdef MHD
    cmaxl=cmax
    nerrl=nerr
    call mpi_reduce(cmaxl,cmax,1,mpi_real,mpi_max,0,mpi_comm_cart,ierr)
    call mpi_reduce(nerrl,nerr,1,mpi_integer,mpi_sum,0,mpi_comm_cart,ierr)

    if (rank == 0) then
      print *,'fluid stats',cmax/freeze,dt*cmax,nerr
    endif
#endif

!! calculate maximum dt from fine mesh force

    f_force_max_node=maxval(f_force_max)

    call mpi_reduce(f_force_max_node,dt_f_acc,1,mpi_real,mpi_max,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
! I"m pretty sure this is incorrect!!! (no wonder this always seemed wrong :-/ )
!      dt_f_acc=1.0/sqrt(min(0.0001,dt_f_acc)*a_mid*G)
      dt_f_acc=1.0/sqrt(max(0.0001,dt_f_acc)*a_mid*G)
      write(*,*) 'maximum timestep from fine force=',dt_f_acc
    endif

    call mpi_bcast(dt_f_acc,1,mpi_real,0,mpi_comm_world,ierr)

#ifdef PPINT

!! calculate maximum dt from particle-particle force
    
    pp_force_max_node=maxval(pp_force_max)

    call mpi_reduce(pp_force_max_node,dt_pp_acc,1,mpi_real,mpi_max,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
      dt_pp_acc=sqrt(dt_pp_scale*rsoft)/max(sqrt(dt_pp_acc*a_mid*G),1e-3)
      write(*,*) 'maximum timestep from pp force=',dt_pp_acc
    endif
   
    call mpi_bcast(dt_pp_acc,1,mpi_real,0,mpi_comm_world,ierr)

    if (pp_test) then
      do i=1,np_local
        print *,i,xv(:,i)
      enddo 
    endif

#endif

!! calculate mass of fine mesh

#ifdef DIAG
    call mpi_reduce(sum(f_mesh_mass),sumrhof,1,mpi_double_precision, &
                    mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'sum of rho_f=',sumrhof
#endif

#ifdef DEBUG_VEL
    write(*,*) rank,xv(:,1:np_local)
#endif

    call coarse_mesh

!! delete all particles outside (1:nc_node_dim]

    call delete_particles

    if (pairwise_ic.or.pair_infall) then
      call report_pair
    endif

  end subroutine particle_mesh
