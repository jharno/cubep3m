!! initialize particle list
  subroutine particle_initialize_fast
    use omp_lib
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

real(4) :: rnum,z_write,dummy, v_r2i
integer(4) :: i,j,k,l,pp,fstat,blocksize,num_writes,nplow,nphigh
integer(8) :: np_total, npl8
character(len=max_path) :: f_zip1,f_zip2,f_zip3
character(len=4) :: rank_s
character(len=7) :: z_s, z_s2
integer(1) :: xi1(4,3), rhoc_i1(4), test_i1
integer(4) :: xi4(3), rhoc_i4=0, np_uzip=0, np_dm=0, np_nu=0
integer(2) :: vi2(3)
equivalence(xi1,xi4)
equivalence(rhoc_i4,rhoc_i1)

np_local=0; xv=0

if (restart_ic) then
  if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
  write(z_s,'(f7.3)') z_write; z_s=adjustl(z_s)
else
  write(z_s,'(f7.3)') z_i
endif
write(rank_s,'(i4)') rank; rank_s=adjustl(rank_s)
! output_path = ic_path
f_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'.dat'
f_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'.dat'
f_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'.dat'
open(11, file=f_zip1, status="old", iostat=fstat, access="stream", buffered='yes')
open(12, file=f_zip2, status="old", iostat=fstat, access="stream", buffered='yes')
open(13, file=f_zip3, status="old", iostat=fstat, access="stream", buffered='yes')
print*,'particle_initialize_fast: opening file:',f_zip1
read(11) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i,shake_offset
print*,'particle_initialize_fast: headers =',np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
!!a=0.999
cur_checkpoint=1
cur_projection=1
cur_halofind=1
do k=1,nc_node_dim
do j=1,nc_node_dim
do i=1,nc_node_dim
  rhoc_i4=0; xi4=0 ! clean up, very imortant.
  read(12) rhoc_i1(1) ! get number of particles in the coarse grid
  if (rhoc_i4==255) read(13) rhoc_i4
  do l=1,rhoc_i4
    np_uzip=np_uzip+1
    read(11) xi1(1,:), vi2
    xv(1:3,np_uzip) = mesh_scale * ( xi4/256. + (/i,j,k/) - 1 )
    xv(4:6,np_uzip) = vi2 / v_r2i
  enddo ! l
enddo
enddo
enddo! k
read(11,end=701) test_i1
print*,'ERROR: rank',rank,': file not ending:',f_zip1
call mpi_abort(mpi_comm_world,ierr,ierr)
701 close(11);close(12);close(13)

#ifdef NEUTRINOS
  f_zip1=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip1_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
  f_zip2=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip2_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
  f_zip3=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'zip3_'//rank_s(1:len_trim(rank_s))//'_nu.dat'
  open(21, file=f_zip1, status="old", iostat=fstat, access="stream", buffered='yes')
  open(22, file=f_zip2, status="old", iostat=fstat, access="stream", buffered='yes')
  open(23, file=f_zip3, status="old", iostat=fstat, access="stream", buffered='yes')
  read(21) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p, v_r2i,shake_offset
  print*,'particle_initialize_fast: headers_nu =',np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset

  do k=1,nc_node_dim
  do j=1,nc_node_dim
  do i=1,nc_node_dim
    rhoc_i4=0; xi4=0 ! clean up, very imortant.
    read(22) rhoc_i1(1)
    if (rhoc_i4==255) read(23) rhoc_i4
    do l=1,rhoc_i4
      np_uzip=np_uzip+1
      read(21) xi1(1,:), vi2
      xv(1:3,np_uzip) = mesh_scale * ( xi4/256. + (/i,j,k/) - 1 )
      xv(4:6,np_uzip) = vi2 / v_r2i
    enddo
  enddo
  enddo
  enddo
  PID(1:np_dm)=1; PID(np_dm+1:np_dm+np_nu)=2
  read(21,end=702) test_i1
  print*,'ERROR: rank',rank,': file not ending:',f_zip1
  call mpi_abort(mpi_comm_world,ierr,ierr)
  702 close(21); close(22); close(23)
#endif

np_local=np_dm+np_nu
if (rank==0) then
  print*,'particle_initialize_fast done:'
  print*,'np_local, np_uzip  =', np_local, np_uzip
  print*,'np_dm, np_nu =', np_dm, np_nu
endif

print*, 'done particle_initialize_fast:'
print*,'min max x  =', minval(xv(1,1:np_local)), maxval(xv(1,1:np_local))
print*,'min max y  =', minval(xv(2,1:np_local)), maxval(xv(2,1:np_local))
print*,'min max z  =', minval(xv(3,1:np_local)), maxval(xv(3,1:np_local))
print*,'min max vx =', minval(xv(4,1:np_local)), maxval(xv(4,1:np_local))
print*,'min max vy =', minval(xv(5,1:np_local)), maxval(xv(5,1:np_local))
print*,'min max vz =', minval(xv(6,1:np_local)), maxval(xv(6,1:np_local))                                                                          
pause
end subroutine particle_initialize_fast
