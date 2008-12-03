!! add mass to fine mesh density on boundry of tile
  subroutine fine_cic_mass_boundry(pp,tile,thread)
    implicit none

    include 'cubepm.fh'

    integer(4) :: pp,thread
    integer(4), dimension(3) :: tile
    integer(4), dimension(3) :: i1, i2
    real(4), dimension(3) :: x, offset, dx1, dx2

    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf !- 0.5 

    do
      if (pp == 0) exit
      x(:) = xv(1:3,pp) + offset(:)
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      if (i1(3) >= 1) then
        if (i1(2) >= 1) then
          if (i1(1) >= 1) rho_f(i1(1),i1(2),i1(3),thread) = &
            rho_f(i1(1),i1(2),i1(3),thread) + dx1(1) * dx1(2) * dx1(3)
          if (i2(1) <= nf_tile) rho_f(i2(1),i1(2),i1(3),thread) = &
            rho_f(i2(1),i1(2),i1(3),thread) + dx2(1) * dx1(2) * dx1(3)
        endif
        if (i2(2) <= nf_tile) then
          if (i1(1) >= 1) rho_f(i1(1),i2(2),i1(3),thread) = &
            rho_f(i1(1),i2(2),i1(3),thread) + dx1(1) * dx2(2) * dx1(3)
          if (i2(1) <= nf_tile) rho_f(i2(1),i2(2),i1(3),thread) = &
            rho_f(i2(1),i2(2),i1(3),thread) + dx2(1) * dx2(2) * dx1(3)
        endif
      endif

      if (i2(3) <= nf_tile) then
        if (i1(2) >= 1) then
          if (i1(1) >= 1) rho_f(i1(1),i1(2),i2(3),thread) = &
            rho_f(i1(1),i1(2),i2(3),thread) + dx1(1) * dx1(2) * dx2(3)
          if (i2(1) <= nf_tile) rho_f(i2(1),i1(2),i2(3),thread) = &
            rho_f(i2(1),i1(2),i2(3),thread) + dx2(1) * dx1(2) * dx2(3)
        endif
        if (i2(2) <= nf_tile) then
          if (i1(1) >= 1) rho_f(i1(1),i2(2),i2(3),thread) = &
            rho_f(i1(1),i2(2),i2(3),thread) + dx1(1) * dx2(2) * dx2(3)
          if (i2(1) <= nf_tile) rho_f(i2(1),i2(2),i2(3),thread) = &
              rho_f(i2(1),i2(2),i2(3),thread) + dx2(1) * dx2(2) * dx2(3)
        endif
      endif

      pp = ll(pp)
    enddo
  end subroutine fine_cic_mass_boundry
