!! add gas density to fine mesh tile 
subroutine fine_gas_mass(tile, thread)
    use mpi_tvd_mhd
!    use omp_lib
    implicit none

    include 'cubepm.fh'

    integer(4), dimension(3) :: tile
    integer(4) :: rlb,rli,rui,rub,rn,rbt,thread
    integer(4), dimension(3) :: igl,igh,lgl,lgh,ugl,ugh,bil,bih

#ifdef DEBUG_GAS_BUFFER
    integer(4) :: i,j,k
#endif

   ! :: gas density buffer arrays

!  6 cell faces
real(4), dimension(nf_buf,nf_physical_node_dim,nf_physical_node_dim) :: gb_xp,gb_xm
real(4), dimension(nf_physical_node_dim,nf_buf,nf_physical_node_dim) :: gb_yp,gb_ym
real(4), dimension(nf_physical_node_dim,nf_physical_node_dim,nf_buf) :: gb_zp,gb_zm
!  12 cell edges
real(4), dimension(nf_buf,nf_buf,nf_physical_node_dim) :: gb_xpyp,gb_xpym,gb_xmyp,gb_xmym
real(4), dimension(nf_buf,nf_physical_node_dim,nf_buf) :: gb_xpzp,gb_xpzm,gb_xmzp,gb_xmzm
real(4), dimension(nf_physical_node_dim,nf_buf,nf_buf) :: gb_ypzp,gb_ypzm,gb_ymzp,gb_ymzm
!  8 cell corners
real(4), dimension(nf_buf,nf_buf,nf_buf) :: gb_xmymzm,gb_xmymzp,gb_xmypzm, &
         gb_xmypzp,gb_xpymzm,gb_xpymzp,gb_xpypzm,gb_xpypzp

common /gasbuffer/ gb_xp,gb_xm,gb_yp,gb_ym,gb_zp,gb_zm,gb_xpyp,gb_xpym,gb_xmyp,gb_xmym, &
                gb_xpzp,gb_xpzm,gb_xmzp,gb_xmzm,gb_ypzp,gb_ypzm,gb_ymzp,gb_ymzm, &
                gb_xmymzm,gb_xmymzp,gb_xmypzm,gb_xmypzp,gb_xpymzm,gb_xpymzp,gb_xpypzm, &
                gb_xpypzp

!      print *,'rank',rank,'tile',tile,'starting init'

!! offsets for rho_f boundaries
    rlb=nf_buf              !! rho_f upper lower buffer cell
    rli=rlb+1               !! rho_f lower internal cell
    rui=nf_tile-nf_buf      !! rho_f upper internal cell 
    rub=rui+1               !! rho_f lower upper buffer cell
    rn=nf_tile              !! rho_f max cell

!! if there is only one fine mesh tile / node then the density allocation
!! is straightforward, but if there are more than one the allocation
!! depends on the tile location within the local portion of the mesh

    if (tiles_node_dim == 1) then
  !! internal
      rho_f(rli:rui,rli:rui,rli:rui,thread)=u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
  !! xmymzm
      rho_f(:rlb,:rlb,:rlb,thread)=gb_xmymzm
  !! xmymzp
      rho_f(:rlb,:rlb,rub:,thread)=gb_xmymzp
  !! xmypzm
      rho_f(:rlb,rub:,:rlb,thread)=gb_xmypzm
  !! xmypzp 
      rho_f(:rlb,rub:,rub:,thread)=gb_xmypzp
  !! xpymzm
      rho_f(rub:rn,:rlb,:rlb,thread)=gb_xpymzm
  !! xpymzp
      rho_f(rub:rn,:rlb,rub:,thread)=gb_xpymzp
  !! xpypzm
      rho_f(rub:rn,rub:,:rlb,thread)=gb_xpypzm
  !! xpypzp 
      rho_f(rub:rn,rub:,rub:,thread)=gb_xpypzp
  !! xmym
      rho_f(:rlb,:rlb,rli:rui,thread)=gb_xmym
  !! xmyp
      rho_f(:rlb,rub:,rli:rui,thread)=gb_xmyp
  !! xmzm
      rho_f(:rlb,rli:rui,:rlb,thread)=gb_xmzm
  !! xmzp
      rho_f(:rlb,rli:rui,rub:,thread)=gb_xmzp
  !! xpym
      rho_f(rub:rn,:rlb,rli:rui,thread)=gb_xpym
  !! xpyp
      rho_f(rub:rn,rub:,rli:rui,thread)=gb_xpyp
  !! xpzm
      rho_f(rub:rn,rli:rui,:rlb,thread)=gb_xpzm
  !! xpzp
      rho_f(rub:rn,rli:rui,rub:,thread)=gb_xpzp
  !! ymzm
      rho_f(rli:rui,:rlb,:rlb,thread)=gb_ymzm
#ifdef DEBUG_GAS_BUFFER
      if (rank==0) then
      do k=1,rlb
        do j=1,rlb
          do i=rli,rui
            print *,'gb_ymzm',i,j,k,rho_f(i,j,k,thread)
          enddo
        enddo
      enddo
      endif
#endif  
  !! ymzp
      rho_f(rli:rui,:rlb,rub:,thread)=gb_ymzp
  !! ypzm
      rho_f(rli:rui,rub:,:rlb,thread)=gb_ypzm
  !! ypzp 
      rho_f(rli:rui,rub:,rub:,thread)=gb_ypzp
  !! xm
      rho_f(:rlb,rli:rui,rli:rui,thread)=gb_xm
  !! xp
      rho_f(rub:rn,rli:rui,rli:rui,thread)=gb_xp
  !! ym
      rho_f(rli:rui,:rlb,rli:rui,thread)=gb_ym
  !! yp
      rho_f(rli:rui,rub:,rli:rui,thread)=gb_yp
  !! zm
      rho_f(rli:rui,rli:rui,:rlb,thread)=gb_zm
  !! zp
      rho_f(rli:rui,rli:rui,rub:,thread)=gb_zp

    else

  !! offset for reading in the upper portion of any buffer
      rbt=nf_physical_node_dim+nf_buf+1-nf_physical_tile_dim

  !! offsets for reading internally out of a buffer
      bil(:)=tile(:)*nf_physical_tile_dim-nf_buf+1
      bih(:)=tile(:)*nf_physical_tile_dim-nf_buf+nf_tile

  !! offsets for gas array

  !! internal tile
      igl(1)=nx%m-nf_buf+tile(1)*nf_physical_tile_dim
      igh(1)=nx%m+nf_buf+(tile(1)+1)*nf_physical_tile_dim-1
      igl(2)=ny%m-nf_buf+tile(2)*nf_physical_tile_dim
      igh(2)=ny%m+nf_buf+(tile(2)+1)*nf_physical_tile_dim-1
      igl(3)=nz%m-nf_buf+tile(3)*nf_physical_tile_dim
      igh(3)=nz%m+nf_buf+(tile(3)+1)*nf_physical_tile_dim-1

  !! lower tile
      lgl(1)=nx%m
      lgh(1)=nx%m+nf_buf+nf_physical_tile_dim-1
      lgl(2)=ny%m
      lgh(2)=ny%m+nf_buf+nf_physical_tile_dim-1
      lgl(3)=nz%m
      lgh(3)=nz%m+nf_buf+nf_physical_tile_dim-1

  !! upper tile
      ugl(1)=nx%n-nf_physical_tile_dim-nf_buf+1
      ugh(1)=nx%n
      ugl(2)=ny%n-nf_physical_tile_dim-nf_buf+1
      ugh(2)=ny%n
      ugl(3)=nz%n-nf_physical_tile_dim-nf_buf+1
      ugh(3)=nz%n

!      print *,'rank',rank,'tile',tile,'finished init'

  !! fill rho_f with appropriate data from buffers and gas density 
  !! depending on which tile we are working on. 
      if (tile(1)==0) then
        if (tile(2)==0) then
          if (tile(3)==0) then
  !internal
            rho_f(rli:rn,rli:,rli:,thread)= u(1,lgl(1):lgh(1),lgl(2):lgh(2),lgl(3):lgh(3))
  !xmymzm
            rho_f(:rlb,:rlb,:rlb,thread)=gb_xmymzm
  !xmym
            rho_f(:rlb,:rlb,rli:,thread)=gb_xmym(:,:,:rui)
  !xmzm
            rho_f(:rlb,rli:,:rlb,thread)=gb_xmzm(:,:rui,:)
  !ymzm
            rho_f(rli:rn,:rlb,:rlb,thread)=gb_ymzm(:rui,:,:)
  !xm
            rho_f(:rlb,rli:,rli:,thread)=gb_xm(:,:rui,:rui)
  !ym
            rho_f(rli:rn,:rlb,rli:,thread)=gb_ym(:rui,:,:rui)
  !zm
            rho_f(rli:rn,rli:,:rlb,thread)=gb_zm(:rui,:rui,:)

          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(rli:rn,rli:,:rui,thread)= u(1,lgl(1):lgh(1),lgl(2):lgh(2),ugl(3):ugh(3))        
  !xmymzp
            rho_f(:rlb,:rlb,rub:,thread)=gb_xmymzp
  !xmym
            rho_f(:rlb,:rlb,:rui,thread)=gb_xmym(:,:,rbt:)
  !xmzp
            rho_f(:rlb,rli:,rub:,thread)=gb_xmzp(:,:rui,:)
  !ymzp
            rho_f(rli:rn,:rlb,rub:,thread)=gb_ymzp(:rui,:,:)
  !xm
            rho_f(:rlb,rli:,:rui,thread)=gb_xm(:,:rui,rbt:)
  !ym
            rho_f(rli:rn,:rlb,:rui,thread)=gb_ym(:rui,:,rbt:)
  !zp
            rho_f(rli:rn,rli:,rub:,thread)=gb_zm(:rui,:rui,:)
          else
  !internal
            rho_f(rli:rn,rli:,:,thread)= u(1,lgl(1):lgh(1),lgl(2):lgh(2),igl(3):igh(3)) 
  !xmym
            rho_f(:rlb,:rlb,:,thread)= gb_xmym(:,:,bil(3):bih(3))
  !xm
            rho_f(:rlb,rli:,:,thread)= gb_xm(:,:rui,bil(3):bih(3))
  !ym
            rho_f(rli:rn,:rlb,:,thread)= gb_ym(:rui,:,bil(3):bih(3))
          endif
        elseif (tile(2)==tiles_node_dim-1) then
          if (tile(3)==0) then
  !internal
            rho_f(rli:rn,:rui,rli:,thread)= u(1,lgl(1):lgh(1),ugl(2):ugh(2),lgl(3):lgh(3))
  !xmypzm
            rho_f(:rlb,rub:,:rlb,thread)=gb_xmypzm
  !xmyp
            rho_f(:rlb,rub:,rli:,thread)=gb_xmyp(:,:,:rui)
  !xmzm
            rho_f(:rlb,:rui,:rlb,thread)=gb_xmzm(:,rbt:,:)
  !ypzm
            rho_f(rli:rn,rub:,:rlb,thread)=gb_ypzm(:rui,:,:)
  !xm
            rho_f(:rlb,:rui,rli:,thread)=gb_xm(:,rbt:,:rui)
  !yp
            rho_f(rli:rn,rub:,rli:,thread)=gb_yp(:rui,:,:rui)
  !zm
            rho_f(rli:rn,:rui,:rlb,thread)=gb_zm(:rui,rbt:,:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(rli:rn,:rui,:rui,thread)= u(1,lgl(1):lgh(1),ugl(2):ugh(2),ugl(3):ugh(3))
  !xmypzp
            rho_f(:rlb,rub:,rub:,thread)=gb_xmypzp
  !xmyp
            rho_f(:rlb,rub:,:rui,thread)=gb_xmyp(:,:,rbt:) 
  !xmzp
            rho_f(:rlb,:rui,rub:,thread)=gb_xmzp(:,rbt:,:)
  !ypzp
            rho_f(rli:rn,rub:,rub:,thread)=gb_ypzp(:rui,:,:)
  !xm
            rho_f(:rlb,:rui,:rui,thread)=gb_xm(:,rbt:,rbt:)
  !yp
            rho_f(rli:rn,rub:,:rui,thread)=gb_yp(:rui,:,rbt:)
  !zp
            rho_f(rli:rn,:rui,rub:,thread)=gb_zp(:rui,rbt:,:)
          else
  !internal 
            rho_f(rli:rn,:rui,:,thread)= u(1,lgl(1):lgh(1),ugl(2):ugh(2),igl(3):igh(3))
  !xmyp
            rho_f(:rlb,rub:,:,thread)= gb_xmyp(:,:,bil(3):bih(3))
  !xm
            rho_f(:rlb,:rui,:,thread)= gb_xm(:,rbt:,bil(3):bih(3))
  !yp
            rho_f(rli:rn,:rlb,:,thread)= gb_yp(:rui,:,bil(3):bih(3))
          endif
        else
          if (tile(3)==0) then
  !internal
            rho_f(rli:rn,:,rli:,thread)= u(1,lgl(1):lgh(1),igl(2):igh(2),lgl(3):lgh(3))
  !xmzm
            rho_f(:rlb,:,:rlb,thread)= gb_xmzm(:,bil(2):bih(2),:)
  !xm
            rho_f(:rlb,:,rli:,thread)= gb_xm(:,bil(2):bih(2),:rui)
  !zm
            rho_f(rli:rn,:,:rlb,thread)= gb_zm(:rui,bil(2):bih(2),:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(rli:rn,:,:rui,thread)= u(1,lgl(1):lgh(1),igl(2):igh(2),ugl(3):ugh(3))
  !xmzp
            rho_f(:rlb,:,rub:,thread)= gb_xmzp(:,bil(2):bih(2),:)
  !xm
            rho_f(:rlb,:,:rui,thread)= gb_xm(:,bil(2):bih(2),rbt:)
  !zp
            rho_f(rli:rn,:,rub:,thread)= gb_zp(:rui,bil(2):bih(2),:)
          else
  !internal
            rho_f(rli:rn,:,:,thread)= u(1,lgl(1):lgh(1),igl(2):igh(2),igl(3):igh(3))
  !xm
            rho_f(:rlb,:,:,thread)= gb_xm(:,bil(2):bih(2),bil(3):bih(3))
          endif
        endif
      elseif (tile(1)==tiles_node_dim-1) then
        if (tile(2)==0) then
          if (tile(3)==0) then
  !internal
            rho_f(:rui,rli:,rli:,thread)= u(1,ugl(1):ugh(1),lgl(2):lgh(2),lgl(3):lgh(3))
  !xpymzm
            rho_f(rub:rn,:rlb,:rlb,thread)= gb_xpymzm 
  !xpym
            rho_f(rub:rn,:rlb,rli:,thread)= gb_xpym(:,:,:rui)
  !xpzm
            rho_f(rub:rn,rli:,:rlb,thread)= gb_xpzm(:,:rui,:)
  !ymzm
            rho_f(:rui,:rlb,:rlb,thread)= gb_ymzm(rbt:,:,:)
  !xp
            rho_f(rub:rn,rli:,rli:,thread)= gb_xp(:,:rui,:rui) 
  !ym
            rho_f(:rui,:rlb,rli:,thread)= gb_ym(rbt:,:,:rui)
  !zm
            rho_f(:rui,rli:,:rlb,thread)= gb_zm(rbt:,:rui,:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(:rui,rli:,:rui,thread)= u(1,ugl(1):ugh(1),lgl(2):lgh(2),ugl(3):ugh(3))
  !xpymzp
            rho_f(rub:rn,:rlb,rub:,thread)= gb_xpymzp
  !xpym
            rho_f(rub:rn,:rlb,:rui,thread)= gb_xpym(:,:,rbt:) 
  !xpzp
            rho_f(rub:rn,rli:,rub:,thread)= gb_xpzp(:,:rui,:)
  !ymzp
            rho_f(:rui,:rlb,rub:,thread)= gb_ymzp(rbt:,:,:)
  !xp
            rho_f(rub:rn,rli:,:rui,thread)= gb_xp(:,:rui,rbt:)
  !ym
            rho_f(:rui,:rlb,:rui,thread)= gb_ym(rbt:,:,rbt:)
  !zp
            rho_f(:rui,rli:,rub:,thread)= gb_zm(rbt:,:rui,:)
          else
  !internal
            rho_f(:rui,rli:,:,thread)= u(1,ugl(1):ugh(1),lgl(2):lgh(2),igl(3):igh(3))
  !xpym
            rho_f(rub:rn,:rlb,:,thread)= gb_xpym(:,:,bil(3):bih(3))
  !xp
            rho_f(rub:rn,rli:,:,thread)= gb_xp(:,:rui,bil(3):bih(3))
  !ym
            rho_f(:rui,:rlb,:,thread)= gb_ym(rbt:,:,bil(3):bih(3))
          endif
        elseif (tile(2)==tiles_node_dim-1) then
          if (tile(3)==0) then
  !internal
            rho_f(:rui,:rui,rli:,thread)= u(1,ugl(1):ugh(1),ugl(2):ugh(2),lgl(3):lgh(3))
  !xpypzm
            rho_f(rub:rn,rub:,:rlb,thread)= gb_xpypzm
  !xpyp
            rho_f(rub:rn,rub:,rli:,thread)= gb_xpyp(:,:,:rui)
  !xpzm
            rho_f(rub:rn,:rui,:rlb,thread)= gb_xpzm(:,rbt:,:)
  !ypzm
            rho_f(:rui,rub:,:rlb,thread)= gb_ypzm(rbt:,:,:)
  !xp
            rho_f(rub:rn,:rui,rli:,thread)= gb_xp(:,rbt:,:rui)
  !yp
            rho_f(:rui,rub:,rli:,thread)= gb_yp(rbt:,:,:rui)
  !zm
            rho_f(:rui,:rui,:rlb,thread)= gb_zm(rbt:,rbt:,:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(:rui,:rui,:rui,thread)= u(1,ugl(1):ugh(1),ugl(2):ugh(2),ugl(3):ugh(3))
  !xpypzp
            rho_f(rub:rn,rub:,rub:,thread)= gb_xpypzp
  !xpyp
            rho_f(rub:rn,rub:,:rui,thread)= gb_xpyp(:,:,rbt:)
  !xpzp
            rho_f(rub:rn,:rui,rub:,thread)= gb_xpzp(:,rbt:,:)
  !ypzp
            rho_f(:rui,rub:,rub:,thread)= gb_ypzp(rbt:,:,:)
  !xp
            rho_f(rub:rn,:rui,:rui,thread)= gb_xp(:,rbt:,rbt:)
  !yp
            rho_f(:rui,rub:,:rui,thread)= gb_yp(rbt:,:,rbt:)
  !zp
            rho_f(:rui,:rui,rub:,thread)= gb_zp(rbt:,rbt:,:)
          else
  !internal
            rho_f(:rui,:rui,:,thread)= u(1,ugl(1):ugh(1),ugl(2):ugh(2),igl(3):igh(3))
  !xpyp
            rho_f(rub:rn,rub:,:,thread)= gb_xpyp(:,:,bil(3):bih(3))
  !xp 
            rho_f(rub:rn,:rui,:,thread)= gb_xp(:,rbt:,bil(3):bih(3))
  !yp
            rho_f(:rui,rub:,:,thread)= gb_yp(rbt:,:,bil(3):bih(3))
          endif
        else
          if (tile(3)==0) then
  !internal
            rho_f(:rui,:,rli:,thread)= u(1,ugl(1):ugh(1),igl(2):igh(2),lgl(3):lgh(3))
  !xpzm
            rho_f(rub:rn,:,:rlb,thread)= gb_xpzm(:,bil(2):bih(2),:)
  !xp
            rho_f(rub:rn,:,rli:,thread)= gb_xp(:,bil(2):bih(2),:rui)
  !zm
            rho_f(:rui,:,:rlb,thread)= gb_zm(rbt:,bil(2):bih(2),:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(:rui,:,:rui,thread)= u(1,ugl(1):ugh(1),igl(2):igh(2),ugl(3):ugh(3))
  !xpzp
            rho_f(rub:rn,:,rub:,thread)= gb_xpzp(:,bil(2):bih(2),:)
  !xp
            rho_f(rub:rn,:,:rui,thread)= gb_xp(:,bil(2):bih(2),rbt:)
  !zp
            rho_f(:rui,:,rub:,thread)= gb_zp(rbt:,bil(2):bih(2),:)
          else
  !internal
            rho_f(:rui,:,:,thread)= u(1,ugl(1):ugh(1),igl(2):igh(2),igl(3):igh(3))
  !xp
            rho_f(rub:rn,:,:,thread)= gb_xp(:,bil(2):bih(2),bil(3):bih(3))
          endif
        endif
      else
        if (tile(2)==0) then
          if (tile(3)==0) then
  !internal
            rho_f(:rn,rli:,rli:,thread)= u(1,igl(1):igh(1),lgl(2):lgh(2),lgl(3):lgh(3))
  !ymzm
            rho_f(:rn,:rlb,:rlb,thread)= gb_ymzm(bil(1):bih(1),:,:)
  !ym
            rho_f(:rn,:rlb,rli:,thread)= gb_ym(bil(1):bih(1),:,:rui)
  !zm
            rho_f(:rn,rli:,:rlb,thread)= gb_zm(bil(1):bih(1),:rui,:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(:rn,rli:,:rui,thread)= u(1,igl(1):igh(1),lgl(2):lgh(2),ugl(3):ugh(3))
  !ymzp
            rho_f(:rn,:rlb,rub:,thread)= gb_ymzp(bil(1):bih(1),:,:)
  !ym
            rho_f(:rn,:rlb,:rui,thread)= gb_ym(bil(1):bih(1),:,rbt:)
  !zp
            rho_f(:rn,rli:,rub:,thread)= gb_zp(bil(1):bih(1),:rui,:)
          else
  !internal
            rho_f(:rn,rli:,:,thread)= u(1,igl(1):igh(1),lgl(2):lgh(2),igl(3):igh(3))
  !ym
            rho_f(:rn,:rlb,:,thread)= gb_ym(bil(1):bih(1),:,bil(3):bih(3))
          endif
        elseif (tile(2)==tiles_node_dim-1) then
          if (tile(3)==0) then
  !internal
            rho_f(:rn,:rui,rli:,thread)= u(1,igl(1):igh(1),ugl(2):ugh(2),lgl(3):lgh(3))
  !ypzm
            rho_f(:rn,rub:,:rlb,thread)= gb_ypzm(bil(1):bih(1),:,:)
  !yp
            rho_f(:rn,rub:,rli:,thread)= gb_yp(bil(1):bih(1),:,:rui)
  !zm
            rho_f(:rn,:rui,:rlb,thread)= gb_zm(bil(1):bih(1),rbt:,:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(:rn,:rui,:rui,thread)= u(1,igl(1):igh(1),ugl(2):ugh(2),ugl(3):ugh(3))
  !ypzp
            rho_f(:rn,rub:,rub:,thread)= gb_ypzp(bil(1):bih(1),:,:)
  !yp
            rho_f(:rn,rub:,:rui,thread)= gb_yp(bil(1):bih(1),:,rbt:)
  !zp
            rho_f(:rn,:rui,rub:,thread)= gb_zp(bil(1):bih(2),rbt:,:)
          else
  !internal
            rho_f(:rn,:rui,:,thread)= u(1,igl(1):igh(1),ugl(2):ugh(2),igl(3):igh(3))
 !yp
            rho_f(:rn,rub:,:,thread)= gb_yp(bil(1):bih(1),:,bil(3):bih(3))
          endif
        else
          if (tile(3)==0) then
  !internal
            rho_f(:rn,:,rli:,thread)= u(1,igl(1):igh(1),igl(2):igh(2),lgl(3):lgh(3))
  !zm
            rho_f(:rn,:,:rlb,thread)= gb_zm(bil(1):bih(1),bil(2):bih(2),:)
          elseif(tile(3)==tiles_node_dim-1) then
  !internal
            rho_f(:rn,:,:rui,thread)= u(1,igl(1):igh(1),igl(2):igh(2),ugl(3):ugh(3))
  !zp
            rho_f(:rn,:,rub:,thread)= gb_zp(bil(1):bih(1),bil(2):bih(2),:)
          else
  !interior
            rho_f(:rn,:,:,thread)=u(1,igl(1):igh(1),igl(2):igh(2),igl(3):igh(3))
          endif
        endif
      endif
  
  endif

  end subroutine fine_gas_mass
