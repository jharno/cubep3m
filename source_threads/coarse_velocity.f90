!! update coarse mesh velocity
#ifdef MHD
  subroutine coarse_velocity(cmax,nerr)
  use mpi_tvd_mhd
  use omp_lib
#else
  subroutine coarse_velocity
    use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif
#endif
    implicit none

#    include "cubepm.fh"

    integer(4) :: i,j,k,pp
    integer(4), dimension(3) :: i1,i2
    real(4), dimension(3) :: x,dx1,dx2,dV

#ifdef MHD
    real(4), parameter :: gg=gamma*(gamma-1)
    integer(4) :: iu,ju,ku,ifc,jfc,kfc,nerr,q
    real(4), dimension(3) :: v,c,acc
    real(4), dimension(5) :: gaz
    real(4) :: cs,cmax
#endif

    call system_clock(count=count_i)

#ifdef MHD
    nerr=0
    cmax=1e-5
    do k=1,nf_physical_node_dim
      kfc=((k-1)/mesh_scale)+1
      ku=nz%m+k-1
      do j=1,nf_physical_node_dim
        jfc=((j-1)/mesh_scale)+1
        ju=ny%m+j-1
        do i=1,nf_physical_node_dim
          ifc=((i-1)/mesh_scale)+1
          iu=nx%m+i-1
!! Not sure about this... 2*dt ??
!          acc= a_mid * G * 2.0 * dt * force_c(:,ifc,jfc,kfc)

!------------------------------------------------
#ifdef MHD_CIC
          !if (rank==0) write(*,*) 'Applying coarse force with CIC on baryons'

          !
          ! Grid offset of -0.25 resolves baryon wiggles 
          !


          x(1) = (1.0/real(mesh_scale)) * i - 0.25 !0.5
          x(2) = (1.0/real(mesh_scale)) * j - 0.25 !0.5
          x(3) = (1.0/real(mesh_scale)) * k - 0.25 !0.5
          i1(:) = floor(x(:)) + 1
          i2(:) = i1(:) + 1
          dx1(:) = i1(:) - x(:)
          dx2(:) = 1.0 - dx1(:)

          gaz=u(:,iu,ju,ku)
          v=gaz(2:4)/gaz(1)
          cs=sqrt(abs(gg*(gaz(5)/gaz(1)-sum(v**2)/2)))

          acc = 0.0
          acc = acc + a_mid * G * dt * dx1(1) * dx1(2) * dx1(3) * force_c(:,ifc-1,jfc-1,kfc-1)
          acc = acc + a_mid * G * dt * dx2(1) * dx1(2) * dx1(3) * force_c(:,ifc,jfc-1,kfc-1)
          acc = acc + a_mid * G * dt * dx1(1) * dx2(2) * dx1(3) * force_c(:,ifc-1,jfc,kfc-1)
          acc = acc + a_mid * G * dt * dx2(1) * dx2(2) * dx1(3) * force_c(:,ifc,jfc,kfc-1)
          acc = acc + a_mid * G * dt * dx1(1) * dx1(2) * dx2(3) * force_c(:,ifc-1,jfc-1,kfc)
          acc = acc + a_mid * G * dt * dx2(1) * dx1(2) * dx2(3) * force_c(:,ifc,jfc-1,kfc)
          acc = acc + a_mid * G * dt * dx1(1) * dx2(2) * dx2(3) * force_c(:,ifc-1,jfc,kfc)
          acc = acc + a_mid * G * dt * dx2(1) * dx2(2) * dx2(3) * force_c(:,ifc,jfc,kfc)
          
          c=cfactor*(abs(v+acc)+cs)
          cmax=max(cmax,maxval(c))
          do q=1,3
            if (c(q) .lt. 0.9/dt) then
              dV(q)=acc(q)
            else
              nerr=nerr+1
              dV(q)=acc(q)-sign(c(q)-0.9/dt,acc(q))
            endif
          enddo
#ifndef NO_MHD_GRAV_COARSE
          u(5,iu,ju,ku)=gaz(5)+sum((gaz(2:4)+gaz(1)*dV/2)*dV)
          u(2:4,iu,ju,ku)=gaz(2:4)+gaz(1)*dV
#endif

#ifdef DEBUG_MHD_CIC
          write(*,*) iu,ju,ku,ifc,jfc,kfc,i1,i2,dx1,dx2,x,acc,dV
#endif

#else
 !---------------------
 ! NGP:
          acc= a_mid * G * dt * force_c(:,ifc,jfc,kfc) 
          gaz=u(:,iu,ju,ku)
          v=gaz(2:4)/gaz(1)
          cs=sqrt(abs(gg*(gaz(5)/gaz(1)-sum(v**2)/2)))
          c=cfactor*(abs(v+acc)+cs)
          cmax=max(cmax,maxval(c))
          do q=1,3
            if (c(q) .lt. 0.9/dt) then
              dV(q)=acc(q)
            else
              nerr=nerr+1
              dV(q)=acc(q)-sign(c(q)-0.9/dt,acc(q))
            endif
          enddo
#ifndef NO_MHD_GRAV_COARSE
          u(5,iu,ju,ku)=gaz(5)+sum((gaz(2:4)+gaz(1)*dV/2)*dV)
          u(2:4,iu,ju,ku)=gaz(2:4)+gaz(1)*dV
#endif
!---------------------
#endif
        enddo
      enddo
    enddo

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('cm g vel',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'coarse gas velocity finished',real(count_f - count_i)/real(count_r)
#endif
    call system_clock(count=count_i)

#endif

    !$omp parallel do default(shared) private(i,j,k,pp,x,i1,i2,dx1,dx2,dV)
    do k=1,nc_node_dim
      do j=1,nc_node_dim
        do i=1,nc_node_dim
          pp=hoc(i,j,k)
          do; if (pp == 0) exit
            x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5
            i1(:) = floor(x(:)) + 1
            i2(:) = i1(:) + 1
#ifdef COARSE_NGP
            dx1(:) = 0.0
            dx2(:) = 1.0
#else
            dx1(:) = i1(:) - x(:)
            dx2(:) = 1.0 - dx1(:)
#endif
            dV = a_mid * G * dt * dx1(1) * dx1(2) * dx1(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i1(2),i1(3)) * dV
            dV = a_mid * G * dt * dx2(1) * dx1(2) * dx1(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i1(2),i1(3)) * dV
            dV = a_mid * G * dt * dx1(1) * dx2(2) * dx1(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i2(2),i1(3)) * dV
            dV = a_mid * G * dt * dx2(1) * dx2(2) * dx1(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i2(2),i1(3)) * dV
            dV = a_mid * G * dt * dx1(1) * dx1(2) * dx2(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i1(2),i2(3)) * dV
            dV = a_mid * G * dt * dx2(1) * dx1(2) * dx2(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i1(2),i2(3)) * dV
            dV = a_mid * G * dt * dx1(1) * dx2(2) * dx2(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i2(2),i2(3)) * dV
            dV = a_mid * G * dt * dx2(1) * dx2(2) * dx2(3)
            xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i2(2),i2(3)) * dV
#ifdef DEBUG_VEL
            write(*,*) 'pp,aGdt,i1(3),i2(3),dx1(3),dx2(3),force_c(1,8)'
            write(*,*) pp,a_mid*G*dt,i1(:),i2(:),dx1(:),dx2(:), &
                       force_c(1,i1(1):i2(1),i1(2):i2(2),i1(3):i2(3))
#endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('cm   vel',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'coarse dark matter velocity finished' &
                            ,real(count_f - count_i)/real(count_r)
#endif

  end subroutine coarse_velocity
