!! write density projections to disk
  subroutine init_projection 
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    character (len=max_path) :: ofile
    character (len=7) :: z_s

    integer(4) :: i,j,fstat
    integer(4), dimension(3) :: tile
    real(8) :: rho_node, rho_tot

#ifdef NEUTRINOS
    integer(1) :: ispec
#endif

!! Initialize projection variables

#ifndef NOPROJ

#ifdef NEUTRINOS
    !! Loop over particle species
    do ispec = 1, 2
#endif

   if(rank == 0) write(*,*) 'Inside init_projection.f90'

    rho_node=0.0
    rho_tot=0.0
    rho_pxy=0.0
    rho_pxz=0.0
    rho_pyz=0.0

!! Construct density projections for each node
    
    do i=1,tiles_node
      tile(3) = (i-1) / (tiles_node_dim * tiles_node_dim)
      j = i - tile(3) * tiles_node_dim * tiles_node_dim
      tile(2) = (j-1) /  tiles_node_dim
      j = j - tile(2) * tiles_node_dim
      tile(1) = j - 1
#ifdef NEUTRINOS
      call build_init_projection(tile,rho_node,ispec)
#else
      call build_init_projection(tile,rho_node)
#endif
    enddo

    call mpi_reduce(rho_node,rho_tot,1,mpi_double_precision, &
                    mpi_sum,0,mpi_comm_world,ierr)
    if (rank==0) write(*,*) 'total projected mass=',rho_tot

!! accumulate on master node
 
    rp_buf=0.0
    call mpi_reduce(rho_pxy,rp_buf,nf_physical_dim*nf_physical_dim, &
                    mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) rho_pxy = rp_buf
   
    rp_buf=0.0
    call mpi_reduce(rho_pxz,rp_buf,nf_physical_dim*nf_physical_dim, &
                    mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) rho_pxz = rp_buf

    rp_buf=0.0
    call mpi_reduce(rho_pyz,rp_buf,nf_physical_dim*nf_physical_dim, &
                    mpi_real,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) rho_pyz = rp_buf

!! Create projection files 

    if (rank == 0) then

      write(z_s,'(f7.3)') z_i
      z_s=adjustl(z_s)

#ifdef NEUTRINOS
      if (ispec == 1) then !! Dark matter
        ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_xy.dat'
      else !! Neutrinos
        ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_xy_nu.dat'
      endif
#else
      ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_xy.dat'
#endif

      open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening projection file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

#ifdef NEUTRINOS
      if (ispec == 1) then !! Dark matter
        ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_xz.dat'
      else !! Neutrinos
        ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_xz_nu.dat'
      endif
#else
      ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_xz.dat'
#endif

      open(unit=13, file=ofile, status="replace", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening projection file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

#ifdef NEUTRINOS
      if (ispec == 1) then !! Dark matter
        ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_yz.dat'
      else !! Neutrinos
        ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_yz_nu.dat'
      endif
#else
      ofile=output_path//z_s(1:len_trim(z_s))//'init_proj_yz.dat'
#endif

      open(unit=14, file=ofile, status="replace", iostat=fstat, access="stream")
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

      write(12) rho_pxy
      write(13) rho_pxz
      write(14) rho_pyz

      close(12)
      close(13)
      close(14)

    endif

#ifdef NEUTRINOS
    enddo !! ispec
#endif

    if (rank == 0) write(*,*) 'Finished projection:',rank

#endif

!! Increment projection counter 

    !cur_projection=cur_projection+1

    projection_step=.false.

    !write(*,*)'wrote initial projection, STOPPING HERE FOR NOW'
    !stop
  end subroutine init_projection 

#ifndef NOPROJ
#ifdef NEUTRINOS
  subroutine build_init_projection(tile,rho_node,ispec)
#else
  subroutine build_init_projection(tile,rho_node)
#endif
    implicit none
#    include "cubepm.fh"

    integer(4) :: os_x, os_y, os_z
    integer(4) :: i,j,k,pp,thread
    integer(4), dimension(3) :: cic_l,cic_h,tile
    real(8) :: rho_node
#ifdef NEUTRINOS
    integer(1) :: ispec
#endif
  
    thread=1
 
    rho_f(:,:,:,thread)=0.0

!! we only need to do 1 extra coarse cell outside of 
!! physical fine mesh to populate non-ghost density

    cic_l(:) = nc_tile_dim * tile(:) 
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + 1

!! calculate fine mesh density for tile

    do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
        do i = cic_l(1), cic_h(1)
          pp=hoc(i,j,k)
#ifdef NEUTRINOS
          call fine_cic_mass(pp,tile,thread,ispec)
#else
          call fine_cic_mass(pp,tile,thread)
#endif
        enddo
      enddo
    enddo

!! calculate offsets for tile in global coords

    os_x=tile(1)*nf_physical_tile_dim + &
              cart_coords(3)*nf_physical_node_dim
    os_y=tile(2)*nf_physical_tile_dim + &
              cart_coords(2)*nf_physical_node_dim    
    os_z=tile(3)*nf_physical_tile_dim + &
              cart_coords(1)*nf_physical_node_dim

!! add mass to projections

    do k=1,nf_physical_tile_dim
      do j=1,nf_physical_tile_dim
        do i=1,nf_physical_tile_dim
!! bad hack to reduce the volume.  Will be nc/nodes_dim thick now.
          if (cart_coords(1)==0) then
            rho_pxy(os_x+i,os_y+j)=rho_pxy(os_x+i,os_y+j)+ &
                                   rho_f(nf_buf+i,nf_buf+j,nf_buf+k,1)
          endif
          if (cart_coords(2)==0) then
            rho_pxz(os_x+i,os_z+k)=rho_pxz(os_x+i,os_z+k)+ &
                                   rho_f(nf_buf+i,nf_buf+j,nf_buf+k,1)
          endif
          if (cart_coords(3)==0) then
            rho_pyz(os_y+j,os_z+k)=rho_pyz(os_y+j,os_z+k)+ &
                                   rho_f(nf_buf+i,nf_buf+j,nf_buf+k,1)
          endif
          rho_node=rho_node+real(rho_f(nf_buf+i,nf_buf+j,nf_buf+k,1),8)
        enddo
      enddo
    enddo

  end subroutine build_init_projection 
#endif
