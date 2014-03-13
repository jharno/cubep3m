!! add mass to fine mesh density within tile
  subroutine fine_cic_mass(pp,tile,thread)
    implicit none

    include 'cubepm.fh'

    integer(4) :: pp,thread
    integer(4), dimension(3) :: tile

    integer(4), dimension(3) :: i1, i2
    real(4), dimension(3) :: x, offset, dx1, dx2

    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf !- 0.5 

#ifndef NEUTRINOS

    do
      if (pp == 0) exit
      x(:) = xv(1:3,pp) + offset(:)
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      rho_f(i1(1),i1(2),i1(3),thread) = rho_f(i1(1),i1(2),i1(3),thread) &
                                       + dx1(1) * dx1(2) * dx1(3)
      rho_f(i2(1),i1(2),i1(3),thread) = rho_f(i2(1),i1(2),i1(3),thread) &
                                       + dx2(1) * dx1(2) * dx1(3)
      rho_f(i1(1),i2(2),i1(3),thread) = rho_f(i1(1),i2(2),i1(3),thread) &
                                       + dx1(1) * dx2(2) * dx1(3)
      rho_f(i2(1),i2(2),i1(3),thread) = rho_f(i2(1),i2(2),i1(3),thread) &
                                       + dx2(1) * dx2(2) * dx1(3)
      rho_f(i1(1),i1(2),i2(3),thread) = rho_f(i1(1),i1(2),i2(3),thread) &
                                       + dx1(1) * dx1(2) * dx2(3)
      rho_f(i2(1),i1(2),i2(3),thread) = rho_f(i2(1),i1(2),i2(3),thread) &
                                       + dx2(1) * dx1(2) * dx2(3)
      rho_f(i1(1),i2(2),i2(3),thread) = rho_f(i1(1),i2(2),i2(3),thread) &
                                       + dx1(1) * dx2(2) * dx2(3)
      rho_f(i2(1),i2(2),i2(3),thread) = rho_f(i2(1),i2(2),i2(3),thread) &
                                       + dx2(1) * dx2(2) * dx2(3)
      pp = ll(pp)
    enddo

#else

  if (.not. doing_halofind) then

    do
      if (pp == 0) exit
      x(:) = xv(1:3,pp) + offset(:)
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      if (PID(pp) > np_dm_total) then !! This is a neutrino
          dx1(1) = mass_p * dx1(1) * mpfac_nt
          dx2(1) = mass_p * dx2(1) * mpfac_nt
      else !! This is a dark matter particle
          dx1(1) = mass_p * dx1(1) * mpfac_dm
          dx2(1) = mass_p * dx2(1) * mpfac_dm
      endif

      rho_f(i1(1),i1(2),i1(3),thread) = rho_f(i1(1),i1(2),i1(3),thread) &
                                       + dx1(1) * dx1(2) * dx1(3)
      rho_f(i2(1),i1(2),i1(3),thread) = rho_f(i2(1),i1(2),i1(3),thread) &
                                       + dx2(1) * dx1(2) * dx1(3)
      rho_f(i1(1),i2(2),i1(3),thread) = rho_f(i1(1),i2(2),i1(3),thread) &
                                       + dx1(1) * dx2(2) * dx1(3)
      rho_f(i2(1),i2(2),i1(3),thread) = rho_f(i2(1),i2(2),i1(3),thread) &
                                       + dx2(1) * dx2(2) * dx1(3)
      rho_f(i1(1),i1(2),i2(3),thread) = rho_f(i1(1),i1(2),i2(3),thread) &
                                       + dx1(1) * dx1(2) * dx2(3)
      rho_f(i2(1),i1(2),i2(3),thread) = rho_f(i2(1),i1(2),i2(3),thread) &
                                       + dx2(1) * dx1(2) * dx2(3)
      rho_f(i1(1),i2(2),i2(3),thread) = rho_f(i1(1),i2(2),i2(3),thread) &
                                       + dx1(1) * dx2(2) * dx2(3)
      rho_f(i2(1),i2(2),i2(3),thread) = rho_f(i2(1),i2(2),i2(3),thread) &
                                       + dx2(1) * dx2(2) * dx2(3)
      pp = ll(pp)
    enddo

  else !! When doing halofind we don't want to include neutrinos in the search for density peaks

    do
      if (pp == 0) exit
      x(:) = xv(1:3,pp) + offset(:)
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      if (PID(pp) <= np_dm_total) then !! This is a dark matter particle 
          dx1(1) = mass_p * dx1(1)
          dx2(1) = mass_p * dx2(1) 

          rho_f(i1(1),i1(2),i1(3),thread) = rho_f(i1(1),i1(2),i1(3),thread) &
                                       + dx1(1) * dx1(2) * dx1(3)
          rho_f(i2(1),i1(2),i1(3),thread) = rho_f(i2(1),i1(2),i1(3),thread) &
                                       + dx2(1) * dx1(2) * dx1(3)
          rho_f(i1(1),i2(2),i1(3),thread) = rho_f(i1(1),i2(2),i1(3),thread) &
                                       + dx1(1) * dx2(2) * dx1(3)
          rho_f(i2(1),i2(2),i1(3),thread) = rho_f(i2(1),i2(2),i1(3),thread) &
                                       + dx2(1) * dx2(2) * dx1(3)
          rho_f(i1(1),i1(2),i2(3),thread) = rho_f(i1(1),i1(2),i2(3),thread) &
                                       + dx1(1) * dx1(2) * dx2(3)
          rho_f(i2(1),i1(2),i2(3),thread) = rho_f(i2(1),i1(2),i2(3),thread) &
                                       + dx2(1) * dx1(2) * dx2(3)
          rho_f(i1(1),i2(2),i2(3),thread) = rho_f(i1(1),i2(2),i2(3),thread) &
                                       + dx1(1) * dx2(2) * dx2(3)
          rho_f(i2(1),i2(2),i2(3),thread) = rho_f(i2(1),i2(2),i2(3),thread) &
                                       + dx2(1) * dx2(2) * dx2(3)

      endif

      pp = ll(pp)
    enddo

  endif

#endif

  end subroutine fine_cic_mass
