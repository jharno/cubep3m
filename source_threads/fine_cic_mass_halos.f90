!! add mass to fine mesh density within tile
  subroutine fine_cic_mass_halos(pp,tile)
    implicit none

    include 'cubep3m.fh'

    integer(4) :: pp,thread
    integer(4), dimension(3) :: tile

    integer(4), dimension(3) :: i1, i2
    real(4), dimension(3) :: x, offset, dx1, dx2

!offset of tile on global grid
    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf !- 0.5 

    do
      if (pp == 0) exit
      x(:) = xv(1:3,pp)*finer_halo_grid + offset(:)
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      rho_f_halos(i1(1),i1(2),i1(3)) = rho_f_halos(i1(1),i1(2),i1(3)) &
                                       + dx1(1) * dx1(2) * dx1(3)
      rho_f_halos(i2(1),i1(2),i1(3)) = rho_f_halos(i2(1),i1(2),i1(3)) &
                                       + dx2(1) * dx1(2) * dx1(3)
      rho_f_halos(i1(1),i2(2),i1(3)) = rho_f_halos(i1(1),i2(2),i1(3)) &
                                       + dx1(1) * dx2(2) * dx1(3)
      rho_f_halos(i2(1),i2(2),i1(3)) = rho_f_halos(i2(1),i2(2),i1(3)) &
                                       + dx2(1) * dx2(2) * dx1(3)
      rho_f_halos(i1(1),i1(2),i2(3)) = rho_f_halos(i1(1),i1(2),i2(3)) &
                                       + dx1(1) * dx1(2) * dx2(3)
      rho_f_halos(i2(1),i1(2),i2(3)) = rho_f_halos(i2(1),i1(2),i2(3)) &
                                       + dx2(1) * dx1(2) * dx2(3)
      rho_f_halos(i1(1),i2(2),i2(3)) = rho_f_halos(i1(1),i2(2),i2(3)) &
                                       + dx1(1) * dx2(2) * dx2(3)
      rho_f_halos(i2(1),i2(2),i2(3)) = rho_f_halos(i2(1),i2(2),i2(3)) &
                                       + dx2(1) * dx2(2) * dx2(3)
      pp = ll(pp)
    enddo

  end subroutine fine_cic_mass_halos
