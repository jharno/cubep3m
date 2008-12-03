!! main particle mesh subroutine
  subroutine particle_mesh
    implicit none
    include 'mpif.h'
    include 'cubepm.fh'

    integer(4) :: i,j,k,dd,cur_tile
    integer(4), dimension(3) :: tile
    integer(4) :: thread
    real(4) :: f_force_max_node
    real(4) :: pp_force_max_node
#ifdef DOPENMP
    integer(4) :: omp_get_thread_num
    external omp_get_thread_num
#endif
! these are for fine mesh
    integer(4) :: pp,ii,im,i3
    integer(4), dimension(3) :: cic_l, cic_h 
! these are for fine ngp mass
    integer(4), dimension(3) :: i1
    real(4),    dimension(3) :: x, offset
! these are for fine velocity
    real(4), dimension(3) :: dx1, dx2
    real(4) :: dVc
    integer(4), dimension(3) :: i2
    integer(4) :: jm,km,ip,jp,kp
    real(4) :: force_mag
#ifdef PPINT
    integer pp1,pp2,ipl(mesh_scale,mesh_scale,mesh_scale)
    real sep(3), force_pp(3), rmag, pp_force_mag, v_init(3) ! pp_force_accum(3)
#endif
 
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

!!!!  IS FFTW2 THREADSAFE???  NEED TO TEST !!!!

    !$omp parallel default(shared) &
    !$omp private(cur_tile,i,j,k,tile,thread,pp,ii,im,i3,cic_l,cic_h,i1,x, &
    !$              offset,dx1,dx2,dVc,i2,jm,km,ip,jp,kp,force_mag,pp1,pp2,ipl,& 
    !$              sep,force_pp,rmag,pp_force_mag,v_init)
!!$    !$omp private(cur_tile,i,j,k,tile,thread,pp,ii,im,i3,cic_l,cic_h,i1,x,offset) &
!!$    !$omp private(dx1,dx2,dVc,i2,jm,km,ip,jp,kp,force_mag,pp1,pp2,ipl,sep,force_pp)&
!!$    !$omp private(rmag,pp_force_mag,v_init)  !for tests on the laptop only!
    thread=1
#ifdef DOPENMP
    thread = omp_get_thread_num() + 1
#endif    
    f_mesh_mass(thread)=0.0
    f_force_max(thread)=0.0
#ifdef PPINT
    pp_force_max(thread)=0.0
#endif
    !$omp do 
    do cur_tile=1,tiles_node
      tile(3) = (cur_tile-1) / (tiles_node_dim * tiles_node_dim)
      j = cur_tile - tile(3) * tiles_node_dim * tiles_node_dim
      tile(2) = (j-1) /  tiles_node_dim
      j = j - tile(2) * tiles_node_dim
      tile(1) = j - 1
!#ifdef MHD
!      call fine_mesh(tile,cmax,nerr,thread)
!#else
!      call fine_mesh(tile,thread)
!#endif

!! normalize fine mesh density
      rho_f(:,:,:,thread)= 0.0
!! calculate coarse mesh offsets
#ifdef NGP
      cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
      cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1
#else
      cic_l(:) = nc_tile_dim * tile(:) + 1 - nc_buf
      cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf
#endif
!! calculate fine mesh density for tile
      do k = cic_l(3), cic_h(3)
        do j = cic_l(2), cic_h(2)
          do i = cic_l(1), cic_h(1)
            pp=hoc(i,j,k)
#ifdef NGP
!            call fine_ngp_mass(pp,tile,thread)

            offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf
! removed the half-cell offset so that fine mesh cells will line up with coarse mesh cells
!    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf - 0.5

            do
              if (pp == 0) exit
              x(:) = xv(1:3,pp) + offset(:)
              i1(:) = floor(x(:)) + 1
              rho_f(i1(1),i1(2),i1(3),thread) = rho_f(i1(1),i1(2),i1(3),thread)+mass_p
              pp = ll(pp)
            enddo

#else
            if (i == cic_l(1) .or. i == cic_h(1) .or. &
                j == cic_l(2) .or. j == cic_h(2) .or. &
                k == cic_l(3) .or. k == cic_h(3)) then
              call fine_cic_mass_boundry(pp,tile,thread)
            else
              call fine_cic_mass(pp,tile,thread)
            endif
#endif
          enddo
        enddo
      enddo
! sum total fine mesh mass
#ifdef DIAG 
      do k=1+nf_buf,nf_tile-nf_buf
        do j=1+nf_buf,nf_tile-nf_buf
          do i=1+nf_buf,nf_tile-nf_buf
            f_mesh_mass(thread)=f_mesh_mass(thread)+real(rho_f(i,j,k,thread),kind=8)
          enddo
        enddo
      enddo
#endif
!! transform and calculate fine mesh force 
      call cubepm_fftw2('f',thread)
#ifdef DEBUG
      print *,'rank',rank,'finished first fft'
#endif
      cmplx_rho_f(:,:,:,thread)=rho_f(:,:,:,thread)

      do i3=1,3
        do k = 1, nf_tile
          do j = 1, nf_tile
            do i = 1, nf_tile/2+1
              ii=2*i
              im=ii-1
              rho_f(im,j,k,thread)=-cmplx_rho_f(ii,j,k,thread)*kern_f(i3,i,j,k)
              rho_f(ii,j,k,thread)=cmplx_rho_f(im,j,k,thread)*kern_f(i3,i,j,k)
            enddo
          enddo
        enddo

#ifdef DEBUG
        print *,'rank',rank,'thread',thread,'finished convolve'
#endif
        call cubepm_fftw2('b',thread)
#ifdef DEBUG
        print *,'rank',rank,'thread',thread,'finished second fft'
#endif

        force_f(i3,:,:,:,thread) = rho_f(nf_buf-1:nf_tile-nf_buf+1,nf_buf-1:nf_tile-nf_buf+1, &
                               nf_buf-1:nf_tile-nf_buf+1,thread)
      enddo
! fine velocity
   
!! calculate max fine mesh force   
      if (.not.pair_infall) then
        force_mag=0.0
        do k=nf_buf-1,nf_tile-nf_buf+1
          do j=nf_buf-1,nf_tile-nf_buf+1
            do i=nf_buf-1,nf_tile-nf_buf+1
              force_mag=sqrt(force_f(1,i,j,k,thread)**2 &
                  +force_f(2,i,j,k,thread)**2+force_f(3,i,j,k,thread)**2)
              if (force_mag > f_force_max(thread)) f_force_max(thread)=force_mag
            enddo
          enddo
        enddo
      endif

!! update dark matter velocity

      offset(:) = real(nf_buf) - tile(:) * nf_physical_tile_dim

!      print *,'thread',thread,'offset',offset,'tile',tile

! removed the half cell offset so that fine mesh cells will line up with coarse mesh cells
!    offset(:) = - 0.5 + real(nf_buf) - tile(:) * nf_physical_tile_dim

      do k = tile(3) * nc_tile_dim + 1, (tile(3) + 1) * nc_tile_dim
        do j = tile(2) * nc_tile_dim + 1, (tile(2) + 1) * nc_tile_dim
          do i = tile(1) * nc_tile_dim + 1, (tile(1) + 1) * nc_tile_dim
            pp = hoc(i,j,k)
#ifdef PPINT
#ifdef DEBUG
            if (pp /= 0) print *,pp,i,j,k
#endif
            ipl=0
#endif
            do
              if (pp == 0) exit
              x(:) = xv(1:3,pp) + offset(:)
              i1(:) = floor(x(:)) + 1

! arghh!!

!              do dd=1,3
!                if (i1(dd) < nf_buf-1 .or. i1(dd) > nf_tile-nf_buf+1) then
!                  print *,'out',thread,i1,x,pp,i,j,k,hoc(i,j,k),xv(:,pp)            
!                  stop
!                endif
!              enddo

#ifdef NGP
              if (pp_test) print *,'before ngp',pp,xv(:,pp)
#ifdef DEBUG
              print *,'force',i1,force_f(:,i1(1),i1(2),i1(3),thread)
#endif
              if (ngp_fmesh_force) xv(4:6,pp)=xv(4:6,pp)+force_f(1:3,i1(1),i1(2),i1(3),thread) * &
                         a_mid * G * dt
              if (pair_infall) then
                force_mag=sqrt(force_f(1,i1(1),i1(2),i1(3),thread)**2+force_f(2,i1(1),i1(2),i1(3),thread)**2+ &
                               force_f(3,i1(1),i1(2),i1(3),thread)**2)
                if (force_mag > f_force_max(thread)) f_force_max(thread)=force_mag
              endif
              if (pp_test) print *,'before pp',pp,xv(:,pp)
#ifdef PPINT
              do im=1,3
                i2(im)=mod(i1(im)-1,mesh_scale)+1
              enddo
              ipl(i2(1),i2(2),i2(3))=ipl(i2(1),i2(2),i2(3))+1
              if (ipl(i2(1),i2(2),i2(3))>max_llf) then
                print *,'exceeded max_llf',max_llf,i1,i2,ipl
                stop
              endif
              llf(ipl(i2(1),i2(2),i2(3)),i2(1),i2(2),i2(3),thread)=pp
#endif
#else
              i2(:) = i1(:) + 1
              dx1(:) = i1(:) - x(:)
              dx2(:) = 1.0 - dx1(:)

              dVc = a_mid * G * dt * dx1(1) * dx1(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i1(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx1(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i1(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx1(1) * dx2(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i2(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx2(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i2(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx1(1) * dx1(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i1(2),i2(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx1(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i1(2),i2(3),thread) * dVc
              dVc = a_mid * G * dt * dx1(1) * dx2(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i2(2),i2(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx2(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i2(2),i2(3),thread) * dVc
#endif
              pp = ll(pp)
            enddo
#ifdef PPINT
            do km=1,mesh_scale
              do jm=1,mesh_scale
                do im=1,mesh_scale
                  if (pp_test .and. ipl(im,jm,km) /= 0) print *,'ipl',im,jm,km,ipl(im,jm,km)
#ifdef DEBUG
                  if ( ipl(im,jm,km) > 1) print *,'ipl',rank,i,j,k,im,jm,km,ipl(im,jm,km)
#endif
                  pp_force_accum(:,:ipl(im,jm,km),thread)=0.0
                  do ip=1,ipl(im,jm,km)-1
                    pp1=llf(ip,im,jm,km,thread)
                    do jp=ip+1,ipl(im,jm,km)
                      pp2=llf(jp,im,jm,km,thread)
                      sep=xv(:3,pp1)-xv(:3,pp2)
                      rmag=sqrt(sep(1)*sep(1)+sep(2)*sep(2)+sep(3)*sep(3))
                      if (rmag>rsoft) then
                        force_pp=mass_p*(sep/(rmag*pp_bias)**3)  !mass_p divides out below
                        pp_force_accum(:,ip,thread)=pp_force_accum(:,ip,thread)-force_pp
                        pp_force_accum(:,jp,thread)=pp_force_accum(:,jp,thread)+force_pp
                        if (pp_force_flag) then
                          xv(4:,pp1)=xv(4:,pp1)-force_pp*a_mid*G*dt
                          xv(4:,pp2)=xv(4:,pp2)+force_pp*a_mid*G*dt
                        endif
                      endif
                    enddo
                  enddo
                  do ip=1,ipl(im,jm,km)
                    pp_force_mag=sqrt(pp_force_accum(1,ip,thread)**2+pp_force_accum(2,ip,thread)**2+pp_force_accum(3,ip,thread)**2)
                    if (pp_force_mag>pp_force_max(thread)) pp_force_max(thread)=pp_force_mag
                  enddo
                enddo
              enddo
            enddo
#endif
          enddo
        enddo
      enddo
! end fine velocity

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
