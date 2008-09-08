subroutine gas_density_buffer
use mpi_tvd_mhd
implicit none
include 'cubepm.fh'
include 'mpif.h'
integer,parameter :: gnf=nf_physical_node_dim

integer :: i,buffer_size,tag
integer :: status(mpi_status_size)

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
!! Arbitrary TAG

   tag=41

!! General idea is to populate buffers with local data, and then pass the data
!! to the disired process.  

   gb_xm=u(1,nx%n-nf_buf+1:nx%n,       ny%m:ny%n,nz%m:nz%n)
!   gb_xp=u(1,         nx%m:nx%m+nf_buf,ny%m:ny%n,nz%m:nz%n)
   gb_xp=u(1,         nx%m:nx%m+nf_buf-1,ny%m:ny%n,nz%m:nz%n)

   gb_ym=u(1,nx%m:nx%n,ny%n-nf_buf+1:ny%n,       nz%m:nz%n)
!   gb_yp=u(1,nx%m:nx%n,         ny%m:ny%m+nf_buf,nz%m:nz%n)
   gb_yp=u(1,nx%m:nx%n,         ny%m:ny%m+nf_buf-1,nz%m:nz%n)

   gb_zm=u(1,nx%m:nx%n,ny%m:ny%n,nz%n-nf_buf+1:nz%n)
!   gb_zp=u(1,nx%m:nx%n,ny%m:ny%n,         nz%m:nz%m+nf_buf)
   gb_zp=u(1,nx%m:nx%n,ny%m:ny%n,         nz%m:nz%m+nf_buf-1)

#ifdef DEBUG_DEN_BUF
  if (rank==0) then
   print *,'initial u(1)',sum(u(1,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n))
   print *,sum(u(1,nx%m:nx%n,ny%n-nf_buf+1:ny%n,       nz%m:nz%n))
   print *,sum(u(1,nx%m:nx%n,         ny%m:ny%m+nf_buf,nz%m:nz%n))
   print *,'initial buffers',sum(gb_xm),sum(gb_xp),sum(gb_ym),sum(gb_yp),sum(gb_zm),sum(gb_zp)
  endif
  !stop
#endif

   buffer_size=gnf*gnf*nf_buf

!send to node in -x
   call mpi_sendrecv_replace(gb_xp,buffer_size,mpi_real, &
                             cart_neighbor(5),tag,cart_neighbor(6), &
                             tag,mpi_comm_cart,status,ierr)

!now need to send edges from gb_xp slab to face-centered neighbors
   buffer_size=gnf*nf_buf*nf_buf
   gb_xpyp=gb_xp(:,:nf_buf,:)
   gb_xpym=gb_xp(:,gnf-nf_buf+1:,:)
   gb_xpzp=gb_xp(:,:,:nf_buf)
   gb_xpzm=gb_xp(:,:,gnf-nf_buf+1:) 

   call mpi_sendrecv_replace(gb_xpyp,buffer_size,mpi_real, &
                             cart_neighbor(3),tag,cart_neighbor(4), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xpym,buffer_size,mpi_real, &
                             cart_neighbor(4),tag,cart_neighbor(3), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xpzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xpzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)

!now need to send corners from gb_xpyp and gb_xpym edges to diagonals 
   buffer_size=nf_buf*nf_buf*nf_buf
   gb_xpypzp=gb_xpyp(:,:,:nf_buf)
   gb_xpypzm=gb_xpyp(:,:,gnf-nf_buf+1:)
   gb_xpymzp=gb_xpym(:,:,:nf_buf)
   gb_xpymzm=gb_xpym(:,:,gnf-nf_buf+1:)

   call mpi_sendrecv_replace(gb_xpypzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xpypzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xpymzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xpymzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)

!send to node in +x
   buffer_size=gnf*gnf*nf_buf
   call mpi_sendrecv_replace(gb_xm,buffer_size,mpi_real, &
                             cart_neighbor(6),tag,cart_neighbor(5), &
                             tag,mpi_comm_cart,status,ierr)

!now need to send edges from gb_xm slab to face-centered neighbors
   buffer_size=gnf*nf_buf*nf_buf
   gb_xmyp=gb_xm(:,:nf_buf,:)
   gb_xmym=gb_xm(:,gnf-nf_buf+1:,:)
   gb_xmzp=gb_xm(:,:,:nf_buf)
   gb_xmzm=gb_xm(:,:,gnf-nf_buf+1:)

   call mpi_sendrecv_replace(gb_xmyp,buffer_size,mpi_real, &
                             cart_neighbor(3),tag,cart_neighbor(4), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xmym,buffer_size,mpi_real, &
                             cart_neighbor(4),tag,cart_neighbor(3), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xmzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xmzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)

!now need to send corners from gb_xmyp and gb_xmym edges to diagonals
   buffer_size=nf_buf*nf_buf*nf_buf
   gb_xmypzp=gb_xmyp(:,:,:nf_buf)
   gb_xmypzm=gb_xmyp(:,:,gnf-nf_buf+1:)
   gb_xmymzp=gb_xmym(:,:,:nf_buf)
   gb_xmymzm=gb_xmym(:,:,gnf-nf_buf+1:)
   
   call mpi_sendrecv_replace(gb_xmypzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xmypzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xmymzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_xmymzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)

!send to node in -y
!    buffer_size=gnf*nf_buf*nf_buf
    buffer_size=gnf*gnf*nf_buf
    call mpi_sendrecv_replace(gb_yp,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)

!now need to send edges from gb_yp slab to face-centered neighbors
   buffer_size=gnf*nf_buf*nf_buf
   gb_ypzp=gb_yp(:,:,:nf_buf)
   gb_ypzm=gb_yp(:,:,gnf-nf_buf+1:)

   call mpi_sendrecv_replace(gb_ypzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_ypzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)

!send to node in +y
   buffer_size=gnf*gnf*nf_buf
   call mpi_sendrecv_replace(gb_ym,buffer_size,mpi_real, &
                             cart_neighbor(4),tag,cart_neighbor(3), &
                             tag,mpi_comm_cart,status,ierr)

!now need to send edges from gb_ym slab to face-centered neighbors
   buffer_size=gnf*nf_buf*nf_buf
   gb_ymzp=gb_ym(:,:,:nf_buf)
   gb_ymzm=gb_ym(:,:,gnf-nf_buf+1:)

   call mpi_sendrecv_replace(gb_ymzp,buffer_size,mpi_real, &
                             cart_neighbor(1),tag,cart_neighbor(2), &
                             tag,mpi_comm_cart,status,ierr)
   call mpi_sendrecv_replace(gb_ymzm,buffer_size,mpi_real, &
                             cart_neighbor(2),tag,cart_neighbor(1), &
                             tag,mpi_comm_cart,status,ierr)

!send to node in -z
   buffer_size=gnf*gnf*nf_buf
   call mpi_sendrecv_replace(gb_zp,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)
!send to node in +z
   call mpi_sendrecv_replace(gb_zm,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &                   
                              tag,mpi_comm_cart,status,ierr)                             

#ifdef DEBUG_GAS_BUFFER
print *,cart_rank,'finished passing'
call mpi_barrier(mpi_comm_cart,ierr)
if (cart_rank == 0) then
  print *,'me',u(1,1,1,1)
  print *,'xm',gb_xm(1,1,1)
  print *,'xp',gb_xp(1,1,1)
  print *,'ym',gb_ym(1,1,1)
  print *,'yp',gb_yp(1,1,1)
  print *,'zm',gb_zm(1,1,1)
  print *,'zp',gb_zp(1,1,1)
  print *,'xmym',gb_xmym(1,1,1)
  print *,'xmyp',gb_xmyp(1,1,1)
  print *,'xmzm',gb_xmzm(1,1,1)
  print *,'xmzp',gb_xmzp(1,1,1)
  print *,'xpym',gb_xpym(1,1,1)
  print *,'xpyp',gb_xpyp(1,1,1)
  print *,'xpzm',gb_xpzm(1,1,1)
  print *,'xpzp',gb_xpzp(1,1,1)
  print *,'ymzm',gb_ymzm(1,1,1)
  print *,'ymzp',gb_ymzp(1,1,1)
  print *,'ypzm',gb_ypzm(1,1,1)
  print *,'ypzp',gb_ypzp(1,1,1)
  print *,'xmymzm',gb_xmymzm(1,1,1)
  print *,'xmymzp',gb_xmymzp(1,1,1)
  print *,'xmypzm',gb_xmypzm(1,1,1)
  print *,'xmypzp',gb_xmypzp(1,1,1)
  print *,'xpymzm',gb_xpymzm(1,1,1)
  print *,'xpymzp',gb_xpymzp(1,1,1)
  print *,'xpypzm',gb_xpypzm(1,1,1)
  print *,'xpypzp',gb_xpypzp(1,1,1)
endif
#endif

end subroutine gas_density_buffer
