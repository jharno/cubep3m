!! add mass to coarse mesh density
  subroutine coarse_cic_mass(pp)
    use omp_lib
    implicit none

#ifdef PPINT
    include 'cubep3m.fh'
#else
    include 'cubepm.fh'
#endif

    integer(4) :: pp
    integer(4), dimension(3) :: i1,i2
    real(4), dimension(3) :: x,dx1,dx2

    do
      if (pp == 0) exit
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
#ifdef MHD
      dx1(1) = mass_p * dx1(1)*(1.0-omega_b/omega_m)
      dx2(1) = mass_p * dx2(1)*(1.0-omega_b/omega_m)
#else
      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)
#endif

#ifdef DEBUG_CCIC
!      if (i1(1) == 1 .and. i1(2) == 1 .and. &
!          i1(3) == 1 .or. i2(1) == 1 .and. &
!          i2(2) == 1 .and. i2(3) == 1) then
        write(*,*) 'internal',pp,x(:),i1(:),i2(:),dx1(:),dx2(:)
!      endif
#endif


      rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) &
                               + dx1(1) * dx1(2) * dx1(3)

      rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) &
                               + dx2(1) * dx1(2) * dx1(3)

      rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) &
                               + dx1(1) * dx2(2) * dx1(3)

      rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) &
                               + dx2(1) * dx2(2) * dx1(3)

      rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) &
                               + dx1(1) * dx1(2) * dx2(3)

      rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) &
                               + dx2(1) * dx1(2) * dx2(3)

      rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) &
                               + dx1(1) * dx2(2) * dx2(3)

      rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) &
                               + dx2(1) * dx2(2) * dx2(3)

      pp = ll(pp)
    enddo

  end subroutine coarse_cic_mass
