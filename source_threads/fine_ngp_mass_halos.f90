!! add mass to fine mesh density within tile using nearest gridpoint scheme
  subroutine fine_ngp_mass_halos(pp,tile)
    implicit none

    include 'cubep3m.fh'

    integer(4)               :: pp !,thread
    integer(4), dimension(3) :: tile,i1
    real(4),    dimension(3) :: x, offset

!offset of tile on local patch grid
    offset(:)= - tile(:) * nf_physical_tile_dim_halos + nf_buf_halos 
!    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf 

!    if (rank == 0)print*,'inside fine_ngp_mass_halos'

    do
      if (pp == 0) exit
!      print*,'check old values',floor(xv(1:3,pp) - tile(:) * nf_physical_tile_dim + nf_buf) + 1
      x(:) = xv(1:3,pp)*finer_halo_grid + offset(:)
      i1(:) = floor(x(:)) + 1
      if(i1(1)==0.or.i1(2)==0.or.i1(3)==0)print*,'check ngp assignment',finer_halo_grid,xv(1:3,pp),offset(:)
!      print*,'fine grid assignment check',i1,thread,xv(1:3,pp)
      rho_f_halos(i1(1),i1(2),i1(3)) = rho_f_halos(i1(1),i1(2),i1(3))+mass_p
      pp = ll(pp)
    enddo

  end subroutine fine_ngp_mass_halos
