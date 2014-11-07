!! initialize particle list
  subroutine particle_initialize
    use omp_lib
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    real(4) :: rnum,z_write,dummy, v_r2i
    integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh, l
    integer(8) :: np_total, npl8, np_uzip
    character(len=max_path) :: ofile1,ofile2,ofile3
    character(len=4) :: rank_s
    character(len=7) :: z_s, z_s2
integer(1) :: xi1,yi1,zi1, rho_c_dm_int1
integer(4) :: xi4,yi4,zi4, rho_c_dm_int4
integer(2) :: vi2(3)
    
#ifdef NEUTRINOS
  integer(4) :: np_dm, np_nu, np_uzip
  integer(8) :: np_total_nu
#endif

equivalence(xi1,xi4)
equivalence(yi1,yi4)
equivalence(zi1,zi4)
equivalence(rho_c_dm_int1,rho_c_dm_int4)
#ifdef NEUTRINOS
  equivalence(rho_c_nu_int1,rho_c_nu_int4)
#endif

if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
write(z_s,'(f7.3)') z_write
z_s=adjustl(z_s)
write(rank_s,'(i4)') rank
rank_s=adjustl(rank_s)

ofile1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoc'//rank_s(1:len_trim(rank_s))//'.dat'
ofile3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoce'//rank_s(1:len_trim(rank_s))//'.dat'
open(unit=21, file=ofile1, status="old", iostat=fstat, access="stream")
open(unit=22, file=ofile2, status="old", iostat=fstat, access="stream")
open(unit=23, file=ofile3, status="old", iostat=fstat, access="stream")
read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i

np_uzip=0
! unzip dm particles
do k=1,nc_node_dim
do j=1,nc_node_dim
do i=1,nc_node_dim
  read(22) rho_c_dm_int1 ! get number of particles in the coarse grid
  if (rho_c_dm_int4>=255) read(23) rho_c_dm_int4
  do l=1,rho_c_dm_int4
    np_uzip=np_uzip+1
    read(21) xi1,yi1,zi1, vi2
    xv(1,np_uzip) = mesh_scale * ( xi4/256. + i )
    xv(2,np_uzip) = mesh_scale * ( yi4/256. + j )
    xv(3,np_uzip) = mesh_scale * ( zi4/256. + k )
    xv(4:6,np_uzip) = vi2 / v_r2i
    np_uzip=np_uzip+1
  enddo
enddo
enddo
enddo

print*,'particle_init: np_local=', np_local
print*,'particle_init: np_uzip=', np_uzip

close(21); close(22); close(23)

#ifdef NEUTRINOS
! unzip nu particles
np_dm=np_local
PID(1:np_dm)=1
ofile1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'_nu.dat'
ofile2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoc'//rank_s(1:len_trim(rank_s))//'_nu.dat'
ofile3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'rhoce'//rank_s(1:len_trim(rank_s))//'_nu.dat'
open(unit=21, file=ofile1, status="old", iostat=fstat, access="stream")
open(unit=22, file=ofile2, status="old", iostat=fstat, access="stream")
open(unit=23, file=ofile3, status="old", iostat=fstat, access="stream")
read(21) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i

do k=1,nc_node_dim
do j=1,nc_node_dim
do i=1,nc_node_dim
  read(22) rho_c_dm_int1 ! get number of particles in the coarse grid
  if (rho_c_dm_int4>=255) read(23) rho_c_dm_int4
  do l=1,rho_c_dm_int4
    np_uzip=np_uzip+1
    read(21) xi1,yi1,zi1, vi2
    xv(1,np_uzip) = mesh_scale * ( xi4/256. + i )
    xv(2,np_uzip) = mesh_scale * ( yi4/256. + j )
    xv(3,np_uzip) = mesh_scale * ( zi4/256. + k )
    xv(4:6,np_uzip) = vi2 / v_r2i
    np_uzip=np_uzip+1
  enddo
enddo
enddo
enddo
np_local=np_dm+np_nu
PID(np_dm+1:np_local)=2
print*,'particle_init: np_nu=', np_nu
print*,'particle_init: np_local=', np_local
print*,'particle_init: np_uzip=', np_uzip

close(21); close(22); close(23)
#endif

end subroutine particle_initialize
