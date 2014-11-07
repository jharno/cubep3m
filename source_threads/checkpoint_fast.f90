! write checkpoints to disk, using integers, 9 bytes per particle, instead of 24.                                                                                                            
subroutine checkpoint
implicit none
include 'mpif.h'
#include "cubepm.fh"

character (len=max_path) :: ofile_dm,ofile_rho_c_dm,ofile_rho_c_dm_e
character (len=4) :: rank_s
character (len=7) :: z_s  
integer(kind=4) :: i,j,k,l,fstat,blocksize,num_writes,nplow,nphigh
integer(kind=4) :: cur_proj,cur_halo
real(kind=4) :: z_write
integer(1) :: rho_c_dm_int1, rho_c_nu_int1, x_int1(3)
integer(4) :: rho_c_dm_int4, rho_c_nu_int4
integer(1), parameter :: nmax_int1=127
real(4) :: vmax, vmax_local, v_r2i
integer(2) :: v_res_int2=16384, v_int2(3)

#ifdef NEUTRINOS
    integer(4) :: np_dm, np_nu, ind_check1, ind_check2
    character (len=max_path) :: ofile_nu,ofile_rho_c_nu,ofile_rho_c_nu_e
#endif

if (rank == 0) z_write=z_checkpoint(cur_checkpoint)
call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

vmax_local = maxval(abs(xv(4:6,:)))
call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)

v_r2i = real(v_res_int2)/vmax

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

ofile_dm=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
ofile_rho_c_dm=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoc'//rank_s(1:len_trim(rank_s))//'.dat'
ofile_rho_c_dm_e=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoce'//rank_s(1:len_trim(rank_s))//'.dat'
#ifdef NEUTRINOS
ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'_nu.dat'
ofile_rho_c_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoc'//rank_s(1:len_trim(rank_s))//'_nu.dat'
ofile_rho_c_nu_e=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoce'//rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif

open(unit=12, file=ofile_dm, status="replace", iostat=fstat, access="stream")
open(unit=13, file=ofile_rho_c_dm, status="replace", iostat=fstat, access="stream")
open(unit=14, file=ofile_rho_c_dm_e, status="replace", iostat=fstat, access="stream")
#ifdef NEUTRINOS
open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream")
open(unit=23, file=ofile_rho_c_nu, status="replace", iostat=fstat, access="stream")
open(unit=24, file=ofile_rho_c_nu_e, status="replace", iostat=fstat, access="stream")
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

  rho_c_dm_int1=-127; rho_c_nu_int1=-127 ! [-127,128] ! reset coarse density
  rho_c_dm_int4=0; rho_c_nu_int4=0
  l=hoc(i,j,k) ! set tail particle

  do while (l>0) ! iterate over all particles
#   ifdef NEUTRINOS
      x_int1(:) = int(fraction( xv(1:3,l)/mesh_scale ) * 256 - 128, kind=1) ! [-128,127]
      v_int2(:) = int(xv(4:6,l) * v_r2i, kind=2) ! [-16384,16384], note that int2:= [-32768,32767]
      if PID(l)=1 then ! dm
        rho_c_dm_int4=rho_c_dm_int4+1 ! increment of density
        write(12) x_int1(:) ! write x
        write(12) v_int2(:) ! write v
      else ! nu
        rho_c_nu_int4=rho_c_nu_int4+1
        write(22) x_int1(:)
        write(22) v_int2(:)
      endif
#   else
      rho_c_dm_int4=rho_c_dm_int4+1 ! dm only
      x_int1(:) = int(fraction( xv(1:3,l)/mesh_scale ) * 256 - 128, kind=1)
      v_int2(:) = int(xv(4:6,l) * v_r2i, kind=2)
      write(12) x_int1(:)
      write(12) v_int2(:)
#   endif
    l=ll(l) ! prev particle
  enddo ! while

!! write dm coarse density field to files
  if (rho_c_dm_int4<255) then
    rho_c_dm_int1=rho_c_dm_int4-128 ! [-128,126]
    write(13) rho_c_dm_int1 ! write density in int1
  else
    write(13) nmax_int1 ! 127
    write(14) rho_c_dm_int4 ! [255,+infty]
  endif

# ifdef NEUTRINOS
!! write also neutrino density field
  if (rho_c_nu_int4<255) then
    rho_c_nu_int1=rho_c_nu_int4-128
    write(23) rho_c_nu_int1
  else
    write(23) nmax_int1
    write(24) rho_c_nu_int4
  endif
#endif

enddo ! i
enddo ! j
enddo ! k

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
