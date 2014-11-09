! write checkpoints to disk, using integers, 9 bytes per particle, instead of 24.                                                                                                            
subroutine checkpoint
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
integer(4) :: rhoc_dm_i4
real(4) :: vmax, vmax_local, v_r2i

#ifdef NEUTRINOS
    integer(4) :: np_dm, np_nu, ind_check1, ind_check2, rhoc_nu_i4
    character (len=max_path) :: fnu_zip1,fnu_zip2,fnu_zip3
#endif

if (rank == 0) z_write=z_checkpoint(cur_checkpoint)
call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

vmax_local = maxval(abs(xv(4:6,:)))
call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)

v_r2i = v_resolution/vmax

!! most linux systems choke when writing more than 2GB of data
!! in one write statement, so break up into blocks < 2GB 
!    blocksize=(2047*1024*1024)/24
! reduced to 32MB chunks because of intel compiler
blocksize=(32*1024*1024)/24
num_writes=np_local/blocksize+1

write(rank_s,'(i4)') rank
rank_s=adjustl(rank_s)
write(z_s,'(f7.3)') z_write
z_s=adjustl(z_s)

fdm_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'.dat'
fdm_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'.dat'
fdm_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'.dat'
#ifdef NEUTRINOS
fnu_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
fnu_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
fnu_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif

open(unit=12, file=fdm_zip1, status="replace", iostat=fstat, access="stream")
open(unit=13, file=fdm_zip2, status="replace", iostat=fstat, access="stream")
open(unit=14, file=fdm_zip3, status="replace", iostat=fstat, access="stream")
#ifdef NEUTRINOS
open(unit=22, file=fnu_zip1, status="replace", iostat=fstat, access="stream")
open(unit=23, file=fnu_zip2, status="replace", iostat=fstat, access="stream")
open(unit=24, file=fnu_zip3, status="replace", iostat=fstat, access="stream")
#endif

cur_checkpoint=cur_checkpoint+1
cur_proj = cur_projection
cur_halo = cur_halofind
if (projection_step) cur_proj = cur_proj + 1
if (halofind_step) cur_halo = cur_halo + 1


#ifdef NEUTRINOS
!! Determine how many dark matter and neutrino particles this rank has
np_dm = count(PID==1)
np_nu = count(PID==0)
if (rank == 0) write(*,*) "checkpoint np_dm, np_nu = ", np_dm, np_nu
write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p, v_r2i
write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p, v_r2i
#else
write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p, v_r2i
#endif

do k=1,nc_node_dim
do j=1,nc_node_dim
do i=1,nc_node_dim

  rhoc_dm_i4=0
# ifdef NEUTRINOS
    rhoc_nu_i4=0
# endif
  l=hoc(i,j,k) ! set tail particle

  do while (l>0) ! iterate over all particles
#   ifdef NEUTRINOS
      if (PID(l)==1) then ! dm
        rhoc_dm_i4=rhoc_dm_i4+1 ! increment of density
        write(12) int( fraction( xv(1:3,l)/mesh_scale ) * 256 ,kind=1) ! write x
        write(12) int( xv(4:6,l) * v_r2i ,kind=2) ! write v
      else ! nu
        rhoc_nu_i4=rhoc_nu_i4+1
        write(22) int( fraction( xv(1:3,l)/mesh_scale ) * 256 ,kind=1)
        write(22) int( xv(4:6,l) * v_r2i ,kind=2)
      endif
#   else
      rhoc_dm_i4=rhoc_dm_i4+1 ! dm only
      write(12) int( fraction( xv(1:3,l)/mesh_scale ) * 256 ,kind=1)
      write(12) int( xv(4:6,l) * v_r2i ,kind=2)
#   endif
    l=ll(l) ! prev particle
  enddo ! while

!! write dm coarse density field to files
  if (rhoc_dm_i4<255) then
    write(13) int(rhoc_dm_i4,kind=1) ! write density in int1
  else
    write(13) int(255,kind=1)
    write(14) rhoc_dm_i4 ! [255,]
  endif

# ifdef NEUTRINOS
!! write also neutrino density field
  if (rhoc_nu_i4<255) then
    write(23) int(rhoc_nu_i4,kind=1)
  else
    write(23) int(255,kind=1)
    write(24) rhoc_nu_i4
  endif
#endif

enddo
enddo
enddo

close(12); close(13); close(14)
#ifdef NEUTRINOS
close(22); close(23); close(24)
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

  end subroutine checkpoint
