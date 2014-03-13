!! write density projections to disk
  subroutine projection_half_finer 
    use omp_lib
#ifdef FFTMKL
   use MKL_DFTI
#endif
    implicit none

    include 'mpif.h'
    include 'cubepm.fh'

    character (len=max_path) :: ofile
    character (len=7) :: z_s

    integer(4), parameter :: finer_factor = proj_finer_factor!4
    !real, dimension(nc*finer_factor, nc*finer_factor):: rho_pxy_finer,rho_pxz_finer, rho_pyz_finer, rp_buf_finer

    integer(4) :: i,j,fstat
    integer(4), dimension(3) :: tile
    real(8) :: rho_node, rho_tot

!! Initialize projection variables
    call mpi_barrier(mpi_comm_world,ierr)
    !if(rank==0) write(*,*) ' *** Inside finer projection ***'
   

    rho_node=0.0
    rho_tot=0.0
    rho_pxy_finer=0.0
    rho_pxz_finer=0.0
    rho_pyz_finer=0.0

!! Construct density projections for each node by looping of (regular) fine mesh
! tiles 
    
    do i=1,tiles_node
      tile(3) = (i-1) / (tiles_node_dim * tiles_node_dim)
      j = i - tile(3) * tiles_node_dim * tiles_node_dim
      tile(2) = (j-1) /  tiles_node_dim
      j = j - tile(2) * tiles_node_dim
      tile(1) = j - 1
      call build_projection_half_finer(tile,rho_node, finer_factor)
      !call build_projection_half_finer(tile,rho_node, finer_factor, rho_pxy_finer,rho_pxz_finer, rho_pyz_finer)
      !write(*,*) 'Done tile ', tile, rank
    enddo

    call mpi_reduce(rho_node,rho_tot,1,mpi_double_precision, &
                    mpi_sum,0,mpi_comm_world,ierr)
    if (rank==0) write(*,*) 'total projected mass=',rho_tot

!! accumulate on master node
    !write(*,*) 'Accumulating on the first node'

    rp_buf_finer=0.0
    call mpi_reduce(rho_pxy_finer,rp_buf_finer,nf_physical_dim*nf_physical_dim*finer_factor*finer_factor, &
                    mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) rho_pxy_finer = rp_buf_finer
   
    rp_buf_finer=0.0
    call mpi_reduce(rho_pxz_finer,rp_buf_finer,nf_physical_dim*nf_physical_dim*finer_factor*finer_factor, &
                    mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) rho_pxz_finer = rp_buf_finer

    rp_buf_finer=0.0
    call mpi_reduce(rho_pyz_finer,rp_buf_finer,nf_physical_dim*nf_physical_dim*finer_factor*finer_factor, &
                    mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) rho_pyz_finer = rp_buf_finer

!! Create projection files 

    if (rank == 0) then

      write(z_s,'(f7.3)') z_projection(cur_projection)
      z_s=adjustl(z_s)

      ofile=output_path//z_s(1:len_trim(z_s))//'proj_half_finer_xy.dat'
#ifdef BINARY
      open (unit=12,file=ofile,status='replace',iostat=fstat,form='binary')
#else
      open (unit=12,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
      if (fstat /= 0) then
        write(*,*) 'error opening projection file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      ofile=output_path//z_s(1:len_trim(z_s))//'proj_half_finer_xz.dat'
#ifdef BINARY
      open (unit=13,file=ofile,status='replace',iostat=fstat,form='binary')
#else
      open (unit=13,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
      if (fstat /= 0) then
        write(*,*) 'error opening projection file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      ofile=output_path//z_s(1:len_trim(z_s))//'proj_half_finer_yz.dat'
#ifdef BINARY
      open (unit=14,file=ofile,status='replace',iostat=fstat,form='binary')
#else
      open (unit=14,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
      if (fstat /= 0) then
        write(*,*) 'error opening projection file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

!! This is the file header

      write(12) a
      write(13) a
      write(14) a

!! Write out data

      write(12) rho_pxy_finer
      write(13) rho_pxz_finer
      write(14) rho_pyz_finer

      close(12)
      close(13)
      close(14)

    endif

   ! write(*,*) 'Finished half finer projection:',rank

!! Increment projection counter already done in regular projection routine

    !cur_projection=cur_projection+1
    !projection_step=.false.

  end subroutine projection_half_finer

  subroutine build_projection_half_finer(tile,rho_node, finer_factor)!, rho_pxy_finer,rho_pxz_finer, rho_pyz_finer)
  !subroutine build_projection_half_finer(tile,rho_node, finer_factor, rho_pxy_finer,rho_pxz_finer, rho_pyz_finer)
    use omp_lib
#ifdef FFTMKL
   use MKL_DFTI
#endif
    implicit none
    include 'cubepm.fh'

    integer(4) :: os_x, os_y, os_z
    integer(4) :: i,j,k,pp, finer_factor
    integer(4), dimension(3) :: cic_l,cic_h,tile
    real(8) :: rho_node
    !real, dimension(finer_factor*nf_tile, finer_factor*nf_tile, finer_factor*nf_tile):: rho_f_finer  

    !real, dimension(nc*finer_factor, nc*finer_factor):: rho_pxy_finer,rho_pxz_finer, rho_pyz_finer
    !thread=1
 
    rho_f_finer(:,:,:)=0.0

    !write(*,*) 'Calling finer build projection for tile ', tile, rank
!! we only need to do 1 extra coarse cell outside of 
!! physical fine mesh to populate non-ghost density

    cic_l(:) = nc_tile_dim * tile(:) 
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + 1

!! calculate sub-fine mesh density for the current tile

    do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
        do i = cic_l(1), cic_h(1)
          pp=hoc(i,j,k)
          !call sub_fine_cic_mass(rho_f_finer,finer_factor, pp,tile)
          call sub_fine_cic_mass(finer_factor, pp,tile)
        enddo
      enddo
    enddo
    !write(*,*) 'Got the CIC mass', tile, rank 
!! calculate offsets for tile in global coords, in sub-fine grid units

    os_x=(tile(1)*nf_physical_tile_dim + &
              cart_coords(3)*nf_physical_node_dim)*finer_factor
    os_y=(tile(2)*nf_physical_tile_dim + &
              cart_coords(2)*nf_physical_node_dim)*finer_factor
    os_z=(tile(3)*nf_physical_tile_dim + &
              cart_coords(1)*nf_physical_node_dim)*finer_factor
    !write(*,*) 'Offsets', os_x, os_y, os_z, tile, rank 
!! add mass to projections

    do k=1,nf_physical_tile_dim*finer_factor
      do j=1,nf_physical_tile_dim*finer_factor
        do i=1,nf_physical_tile_dim*finer_factor
!#ifdef REDUCE_PROJ
!! bad hack to reduce the volume.  Will be nc/nodes_dim thick now.
! Need to set this manually to the number of nodes I want to include in the collapse.
! i.e. for a 4**3 node run, half the box would be cart_coords == 0 and 1.
 
          if (cart_coords(1)==0 .or. cart_coords(1) == 1) then
            rho_pxy_finer(os_x+i,os_y+j)=rho_pxy_finer(os_x+i,os_y+j)+ &
                                   rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k)
          endif
          if (cart_coords(2)==0 .or. cart_coords(2) == 1) then
            rho_pxz_finer(os_x+i,os_z+k)=rho_pxz_finer(os_x+i,os_z+k)+ &
                                   rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k)
          endif
          if (cart_coords(3)==0 .or. cart_coords(3) == 1) then
            rho_pyz_finer(os_y+j,os_z+k)=rho_pyz_finer(os_y+j,os_z+k)+ &
                                   rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k)
          endif
!#else
!          rho_pxy(os_x+i,os_y+j)=rho_pxy(os_x+i,os_y+j)+ &
!                                 rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k,1)
!          rho_pxz(os_x+i,os_z+k)=rho_pxz(os_x+i,os_z+k)+ &
!                                 rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k,1)
!          rho_pyz(os_y+j,os_z+k)=rho_pyz(os_y+j,os_z+k)+ &
!                                 rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k,1)          
!#endif          
          rho_node=rho_node+real(rho_f_finer(nf_buf*finer_factor+i,nf_buf*finer_factor+j,nf_buf*finer_factor+k),8)
        enddo
      enddo
    enddo

  end subroutine build_projection_half_finer

!----------

!! add mass to fine mesh density within tile
  subroutine sub_fine_cic_mass(finer_factor, pp,tile)
  !subroutine sub_fine_cic_mass(rho_f_finer, finer_factor, pp,tile)
    use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif
    implicit none

    include 'cubepm.fh'

    integer(4) :: pp, finer_factor!
    integer(4), dimension(3) :: tile

    integer(4), dimension(3) :: i1, i2
    real(4), dimension(3) :: x, offset, dx1, dx2
    !real, dimension(finer_factor*nf_tile, finer_factor*nf_tile,finer_factor*nf_tile):: rho_f_finer


    ! Bring particles from local node coordidinates to local tile coordinates:
    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf !- 0.5 

    do
      if (pp == 0) exit
      x(:) = (xv(1:3,pp) + offset(:))*finer_factor
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      rho_f_finer(i1(1),i1(2),i1(3)) = rho_f_finer(i1(1),i1(2),i1(3)) &
                                       + dx1(1) * dx1(2) * dx1(3)
      rho_f_finer(i2(1),i1(2),i1(3)) = rho_f_finer(i2(1),i1(2),i1(3)) &
                                       + dx2(1) * dx1(2) * dx1(3)
      rho_f_finer(i1(1),i2(2),i1(3)) = rho_f_finer(i1(1),i2(2),i1(3)) &
                                       + dx1(1) * dx2(2) * dx1(3)
      rho_f_finer(i2(1),i2(2),i1(3)) = rho_f_finer(i2(1),i2(2),i1(3)) &
                                       + dx2(1) * dx2(2) * dx1(3)
      rho_f_finer(i1(1),i1(2),i2(3)) = rho_f_finer(i1(1),i1(2),i2(3)) &
                                       + dx1(1) * dx1(2) * dx2(3)
      rho_f_finer(i2(1),i1(2),i2(3)) = rho_f_finer(i2(1),i1(2),i2(3)) &
                                       + dx2(1) * dx1(2) * dx2(3)
      rho_f_finer(i1(1),i2(2),i2(3)) = rho_f_finer(i1(1),i2(2),i2(3)) &
                                       + dx1(1) * dx2(2) * dx2(3)
      rho_f_finer(i2(1),i2(2),i2(3)) = rho_f_finer(i2(1),i2(2),i2(3)) &
                                       + dx2(1) * dx2(2) * dx2(3)
      pp = ll(pp)
    enddo

  end subroutine sub_fine_cic_mass

