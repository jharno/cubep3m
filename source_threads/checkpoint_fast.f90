! write checkpoints to disk, using integers, 9 bytes per particle, instead of 24.                                                                                                            
subroutine checkpoint_fast
implicit none
include 'mpif.h'
#include "cubepm.fh"

character (len=max_path) :: fdm_zip1,fdm_zip2,fdm_zip3
character (len=4) :: rank_s
character (len=7) :: z_s  
integer(kind=4) :: i,j,k,l,fstat,blocksize,num_writes,nplow,nphigh
integer(kind=4) :: cur_proj,cur_halo
real(kind=4) :: z_write
integer(4) :: v_resolution=16384
integer(4) :: rhoc_dm_i4, rhoc_nu_i4, np_dm, np_nu
real(4) :: vmax, vmax_local, v_r2i
#ifdef NEUTRINOS
    integer(4) :: ind_check1, ind_check2
    character (len=max_path) :: fnu_zip1,fnu_zip2,fnu_zip3
#endif

if (rank == 0) z_write=z_checkpoint(cur_checkpoint)
print*,'checkpoint_fast =', cur_checkpoint, z_write
call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
vmax_local = maxval(abs(xv(4:6,1:np_local)))
call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
v_r2i = v_resolution/vmax

write(rank_s,'(i4)') rank; rank_s=adjustl(rank_s)
write(z_s,'(f7.3)') z_write; z_s=adjustl(z_s)

fdm_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'.dat'
fdm_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'.dat'
fdm_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'.dat'
open(unit=11, file=fdm_zip1, status="replace", iostat=fstat, access="stream", buffered='yes')
open(unit=12, file=fdm_zip2, status="replace", iostat=fstat, access="stream", buffered='yes')
open(unit=13, file=fdm_zip3, status="replace", iostat=fstat, access="stream", buffered='yes')
np_dm=np_local
#ifdef NEUTRINOS
  fnu_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
  fnu_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
  fnu_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
  open(unit=21, file=fnu_zip1, status="replace", iostat=fstat, access="stream", buffered='yes')
  open(unit=22, file=fnu_zip2, status="replace", iostat=fstat, access="stream", buffered='yes')
  open(unit=23, file=fnu_zip3, status="replace", iostat=fstat, access="stream", buffered='yes')
  np_dm = count(PID(1:np_local)==1); np_nu = count(PID(1:np_local)==2)
#endif

cur_checkpoint=cur_checkpoint+1
cur_proj = cur_projection
cur_halo = cur_halofind
if (projection_step) cur_proj = cur_proj + 1
if (halofind_step) cur_halo = cur_halo + 1

write(11) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p,v_r2i,shake_offset
print*,'checkpoint_fast: headers =',np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p,v_r2i,shake_offset

#ifdef NEUTRINOS
  shake_offset=0
  write(21) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p, v_r2i,shake_offset
  print*,'checkpoint_fast: headers_nu =',np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p,v_r2i,shake_offset
#endif
if (rank == 0) then
  print*, "checkpoint_fast: np_dm, np_nu =", np_dm, np_nu
  print*, "checkpoint_fast: vmax =", vmax
  print*, "v_r2i =",v_r2i
endif

print*, 'will do checkpoint_fast:'
print*,'min max x  =', minval(xv(1,1:np_local)), maxval(xv(1,1:np_local))
print*,'min max y  =', minval(xv(2,1:np_local)), maxval(xv(2,1:np_local))
print*,'min max z  =', minval(xv(3,1:np_local)), maxval(xv(3,1:np_local))
print*,'min max vx =', minval(xv(4,1:np_local)), maxval(xv(4,1:np_local))
print*,'min max vy =', minval(xv(5,1:np_local)), maxval(xv(5,1:np_local))
print*,'min max vz =', minval(xv(6,1:np_local)), maxval(xv(6,1:np_local))

!!!!!!!!!!!!!!!! Loop over physical coarse grids !!!!!!!!!!!!!!!!
do k=1,nc_node_dim
do j=1,nc_node_dim
do i=1,nc_node_dim
  rhoc_dm_i4=0; rhoc_nu_i4=0
  l=hoc(i,j,k) ! set tail particle
  do while (l>0) ! iterate over all particles
#   ifdef NEUTRINOS
      if (PID(l)==1) then ! dm
        rhoc_dm_i4=rhoc_dm_i4+1 ! increment of density
        write(11) int( mod( xv(1:3,l)/mesh_scale, 1. ) * 256 ,kind=1) ! write x
        write(11) int( xv(4:6,l) * v_r2i ,kind=2) ! write v
      else ! nu
        rhoc_nu_i4=rhoc_nu_i4+1
        write(21) int( mod( xv(1:3,l)/mesh_scale, 1. ) * 256 ,kind=1)
        write(21) int( xv(4:6,l) * v_r2i ,kind=2)
      endif ! PID
#   else
      rhoc_dm_i4=rhoc_dm_i4+1 ! dm only
      write(11) int( mod( xv(1:3,l)/mesh_scale, 1. ) * 256 ,kind=1)
      write(11) int( xv(4:6,l) * v_r2i ,kind=2)
#   endif
    l=ll(l) ! prev particle
  enddo ! while

  if (rhoc_dm_i4<255) then
    write(12) int(rhoc_dm_i4,kind=1) ! write density in int1
  else
    write(12) int(255,kind=1)
    write(13) rhoc_dm_i4 ! [255,]
  endif
# ifdef NEUTRINOS
    if (rhoc_nu_i4<255) then
      write(22) int(rhoc_nu_i4,kind=1)
    else
      write(22) int(255,kind=1)
      write(23) rhoc_nu_i4
    endif
# endif
enddo
enddo
enddo ! k
!!!!!!!!!!!!!!!! Loop over physical coarse grids !!!!!!!!!!!!!!!!

close(11); close(12); close(13)
#ifdef NEUTRINOS
  close(21); close(22); close(23)
#endif

#ifdef MHD
!! Write gas checkpoint
    call mpi_tvd_mhd_state_output(output_path,nts,t,z_s)
#endif

if (rank==0) then
   print*, 'fast checkpoint done:'
   print*, 'cur_checkpoint =', cur_checkpoint
   print*, 'cur_projection =', cur_proj
   print*, 'cur_halofind   =', cur_halo
endif

    write(*,*) 'Finished fast checkpoint:',rank

    checkpoint_step=.false.

  end subroutine checkpoint_fast
