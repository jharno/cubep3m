!! update velocity with fine mesh contribution for particles in tile
#ifdef MHD
  subroutine fine_velocity(tile,cmax,nerr,thread)
    use mpi_tvd_mhd
!    use omp_lib
#else
  subroutine fine_velocity(tile,thread)
!    use omp_lib
#endif
    implicit none

    include 'cubepm.fh'

    integer(4), dimension(3) :: tile
    integer(4) :: thread
    real(4), dimension(3) :: offset, x, dx1, dx2
    real(4) :: dVc
    integer(4), dimension(3) :: i1, i2
    integer(4) :: i,j,k,im,jm,km,ip,jp,kp,pp
    real(4) :: force_mag

#ifdef PPINT
    integer pp1,pp2,ipl(mesh_scale,mesh_scale,mesh_scale)
    real sep(3), force_pp(3), rmag, pp_force_mag, v_init(3) ! pp_force_accum(3)
#endif

#ifdef MHD
    real, parameter :: gg=gamma*(gamma-1)
    integer nerr,q,kff,ku,jff,ju,iff,iu,u_offset(3)
    real cmax,c(3),v(3),cs,dv(3)
    real gaz(5),acc(3)
    real, parameter ::T_CMB = 2.725 ! in K
    real, parameter :: k_B = 1.38065E-23 ! in J/K
    real, parameter  :: h = 0.701
    real, parameter :: mu = 1.22 ! reduced mass
    real, parameter :: mproton = 1.6726E-27 ! in kg
    real :: Nprime, Ephys2sim, Econst
    real    :: E_thermal
    real    :: E_kinetic
    real    :: z

    z = 1./a - 1.
    E_thermal = 0.

#endif

#ifdef FF_BOUND_DEBUG
    if (rank==0) then
     open(51,file='fforce.dat',form='formatted',status='replace')
     do i=nf_buf-1,nf_tile-nf_buf+1
       write(51,*) force_f(:,i,i,i,thread)
     enddo
     close(51)
    endif
#endif

!    if (.not.pair_infall) then
!      force_mag=0.0
!      do k=nf_buf-1,nf_tile-nf_buf+1
!        do j=nf_buf-1,nf_tile-nf_buf+1
!          do i=nf_buf-1,nf_tile-nf_buf+1
!            force_mag=sqrt(force_f(1,i,j,k,thread)**2 &
!                +force_f(2,i,j,k,thread)**2+force_f(3,i,j,k,thread)**2)
!            if (force_mag > f_force_max(thread)) f_force_max(thread)=force_mag
!          enddo
!        enddo
!      enddo
!    endif

#ifdef MHD
!! update gas velocity

!! need to update all gas cells local to the current tile
!! gas in local cells / tile depends on the tile position within the local u array:

    u_offset(1)=nx%m-1+tile(1)*nf_physical_tile_dim
    u_offset(2)=ny%m-1+tile(2)*nf_physical_tile_dim
    u_offset(3)=nz%m-1+tile(3)*nf_physical_tile_dim

    do k=1,nf_physical_tile_dim
      kff=k+nf_buf
      ku=k+u_offset(3)
      do j=1,nf_physical_tile_dim
        jff=j+nf_buf
        ju=j+u_offset(2)
        do i=1,nf_physical_tile_dim
          iff=i+nf_buf
          iu=i+u_offset(1)
!! NEED TO SORT OUT TIMESTEP dt IN BELOW (should be 2 * dt!)
!          acc= a_mid * G * 2.0 * dt * force_f(:,iff,jff,kff,thread)
          acc= a_mid * G * dt * force_f(:,iff,jff,kff,thread)
          gaz=u(:,iu,ju,ku)
          v=gaz(2:4)/gaz(1)
          cs=sqrt(abs(gg*(gaz(5)/gaz(1)-sum(v**2)/2)))
          c=cfactor*(abs(v+acc)+cs)
          cmax=max(cmax,maxval(c))
          do q=1,3
            if (c(q) .lt. 0.9/dt) then
              dv(q)=acc(q)
            else
              nerr=nerr+1
              dv(q)=acc(q)-sign(c(q)-0.9/dt,acc(q))
            endif
          enddo
#ifndef NO_MHD_GRAV_FINE
          u(5,iu,ju,ku)=gaz(5)+sum((gaz(2:4)+gaz(1)*dv/2)*dv) 
          u(2:4,iu,ju,ku)=gaz(2:4)+gaz(1)*dv
#endif

#ifdef CMB_coupling
    
          ! Temperature coupling to CMB for z > z_dec = 150 
          ! T_CMB(z) evolves as (1+z)*T_CMB 

          if (z > 150.) then
      
            !! Nprime is the number of physical particles represented by each sim particle
            !! Ephys2sim converts physical energy units (Joules) to simulation units
            !! Econst stores the remaing numerical factors from Nprime and Ephys2sim

            Econst = (4. / 9.) * 1.e-10
            Nprime = omega_b * box**3 / mu / mproton / real(nc)**3
            Ephys2sim = a**2 * real(nc)**5 / omega_m**2 / box**5
            E_thermal = u(1,iu,ju,ku) * Econst * Nprime * k_B * T_CMB * (1. + z) * Ephys2sim

            v=u(2:4,iu,ju,ku)/u(1,iu,ju,ku)
            E_kinetic=u(1,iu,ju,ku)*sum(v**2)/2

            u(5,iu,ju,ku) = E_kinetic + E_thermal


            ! Else it is no longer coupled to the CMB. 
            ! Setting E_thermal to zero would mean zero coupling AND zero
            ! temperature in the baryons, which we don't usually want.             

          endif 
#endif
        enddo
      enddo
    enddo
            ! using the line below to check the order of Kinetic/Thermal Energy
            ! print*,'CMB coupling,E_k=',E_kinetic,'E_thermal=',E_thermal
#endif
#ifdef DM_VEL_OLD
!! update dark matter velocity

    offset(:) = real(nf_buf) - tile(:) * nf_physical_tile_dim
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
!            print *,pp,i,j,k
#ifdef NGP
            if (pp_test) print *,'before ngp',pp,xv(:,pp)
#ifdef DEBUG
            print *,'force',i1,force_f(:,i1(1),i1(2),i1(3),thread)
#endif
            if (ngp_fmesh_force) xv(4:6,pp)=xv(4:6,pp)+force_f(1:3,i1(1),i1(2),i1(3),thread) * &
                       a_mid * G * dt
!            xv(4:6,pp)=xv(4:6,pp)+force_f(1:3,i1(1),i1(2),i1(3)) * a_mid * G * dt
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
        !    vel_state(:,ipl(i2(1),i2(2),i2(3)),i2(1),i2(2),i2(3))=xv(4:,pp)
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
! pp_force_mag, pp_force_max
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
                  !if (ip==1) v_init=xv(4:,pp1) 
                  do jp=ip+1,ipl(im,jm,km)
                    pp2=llf(jp,im,jm,km,thread)
                    sep=xv(:3,pp1)-xv(:3,pp2)
                    rmag=sqrt(sep(1)*sep(1)+sep(2)*sep(2)+sep(3)*sep(3))
                  !  print *,'rmag',rmag,'sep',sep
                    if (rmag>rsoft) then
                  !    force_pp=a_mid*G*mass_p*(sep/rmag**3)
!                      force_pp=mass_p*(((sep+sep_bias)/rmag**3)+pp_bias)  !mass_p divides out below
                      force_pp=mass_p*(sep/(rmag*pp_bias)**3)  !mass_p divides out below
                      pp_force_accum(:,ip,thread)=pp_force_accum(:,ip,thread)-force_pp
                      pp_force_accum(:,jp,thread)=pp_force_accum(:,jp,thread)+force_pp 
                   !!   if (ip==1) then
                   !!     pp_force_accum=pp_force_accum+force_pp
                   !!   endif
                 !    print *,'force_pp',force_pp
                      if (pp_force_flag) then
                        xv(4:,pp1)=xv(4:,pp1)-force_pp*a_mid*G*dt
                        xv(4:,pp2)=xv(4:,pp2)+force_pp*a_mid*G*dt
                      endif
                    endif                    
                  enddo
       !!           if (ip==1) then
              !      pp_force_mag=sqrt( (xv(4,pp1)-v_init(1))**2 + (xv(5,pp1)-v_init(2))**2 + &
              !                         (xv(6,pp1)-v_init(3))**2 ) / (a_mid*G*dt) 
        !!            pp_force_mag=sqrt(pp_force_accum(1)**2+pp_force_accum(2)**2+pp_force_accum(3)**2)
         !!         endif
                enddo 
                do ip=1,ipl(im,jm,km)
                  pp_force_mag=sqrt(pp_force_accum(1,ip,thread)**2+pp_force_accum(2,ip,thread)**2+pp_force_accum(3,ip,thread)**2)
                  if (pp_force_mag>pp_force_max(thread)) pp_force_max(thread)=pp_force_mag
                enddo
         !       do ip=1,ipl(im,jm,km)
         !         pp1=llf(ip,im,jm,km)
         !         pp_force_mag=sqrt( (xv(4,pp1)-vel_state(1,ip,im,jm,km))**2 + (xv(5,pp1)-vel_state(2,ip,im,jm,km))**2 + &
         !                              (xv(6,pp1)-vel_state(3,ip,im,jm,km))**2 ) / (a_mid*G*dt)
         !!       if (pp_force_mag>pp_force_max) pp_force_max=pp_force_mag
         !       enddo
              enddo
            enddo
          enddo
#endif
        enddo
      enddo
    enddo
! end of #ifdef DM_VEL_OLD:
#endif
end subroutine fine_velocity
