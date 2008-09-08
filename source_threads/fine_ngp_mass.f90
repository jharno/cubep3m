!! add mass to fine mesh density within tile using nearest gridpoint scheme
  subroutine fine_ngp_mass(pp,tile,thread)
    implicit none

    include 'cubepm.fh'

    integer(4)               :: pp,thread
    integer(4), dimension(3) :: tile,i1
    real(4),    dimension(3) :: x, offset

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

  end subroutine fine_ngp_mass
