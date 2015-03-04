subroutine particle_initialize
use omp_lib
implicit none
include 'mpif.h'
#include "cubepm.fh"

real(4) :: rnum,z_write,dummy
integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh
integer*8 :: np_total,npl8
character(len=max_path) :: ofile
character(len=6) :: rank_s
character(len=7) :: z_s, z_s2
integer(4) :: np_nu
#ifdef CHECK_IP
  real(8) :: xva(6)
#endif
#ifdef NEUTRINOS
    integer(4) :: np_dm
    integer(8) :: np_total_nu
  !! Check that the PID flag is also defined
# ifndef PID_FLAG
      write(*,*) "ERROR: Using Neutrinos but PID is not enabled !!"
      call mpi_abort(mpi_comm_world,ierr,ierr)
# endif
#endif
#if defined(ZIP) || defined(ZIPDM)
  character(len=max_path) :: f_zip0,f_zip1,f_zip2,f_zip3
  integer :: fstat0, fstat1, fstat2, fstat3, l
  real(4) :: v_r2i,v_r2i_nu
  integer(4) :: np_uzip
  integer(1) :: xi1(4,3), rhoc_i1(4), test_i1
  integer(4) :: xi4(3), rr_i4
  integer(2) :: vi2(3)
  equivalence(xi1,xi4)
  equivalence(rr_i4,rhoc_i1)
#endif
  integer(4) :: hostnm
  character(len=100) :: myhost
  real(8) :: t1, t2

if (rank==0) print*,'particle_initialize'
fstat=0
np_local=(nf_physical_node_dim/2)**3

!! set flag in cubepm.par to select desired initial conditions 
if (random_ic) then !! use random numbers as particle x and v=0
      xv=0
      call random_number(xv(:3,:np_local))
      xv=xv*nf_physical_node_dim
elseif (grid_ic) then !! put particles at the middle of the grid and v=0
      pp=0
      do k=0,nf_physical_node_dim-1,2
        do j=0,nf_physical_node_dim-1,2
          do i=0,nf_physical_node_dim-1,2
            pp=pp+1
            xv(1,pp) = i + 0.5
            xv(2,pp) = j + 0.5
            xv(3,pp) = k + 0.5
            xv(4:6,pp) = 0.0
          enddo
        enddo
      enddo
      if (pp .ne. np_local) then
        write(*,*) 'pp ne np_local',pp,np_local,'rank',rank
        stop
      endif
      write(*,*) 'rank',rank,'finished grid_ic'
elseif (pairwise_ic.or.pair_infall) then !! pairwise
      if (rank == 0) then
        np_local=2
        xv(:3,1:2) = 1.0
        xv(4:,1:2) = 0.0
      else
        np_local=0
      endif
elseif (pp_test) then !! pp test
      if (rank==0) then
        np_local=4
        xv(1:3,1)=0.3+nf_physical_node_dim/2
        xv(1:3,2)=0.2+nf_physical_node_dim/2
        xv(1,3)=0.2+nf_physical_node_dim/2
        xv(2:3,3)=0.3+nf_physical_node_dim/2
        xv(1,4)=0.3+nf_physical_node_dim/2
        xv(2:3,4)=0.2+nf_physical_node_dim/2
        xv(4:,1:4)=0.0
      else
        np_local=0
      endif
elseif (shake_test_ic) then !! shake test
        np_local=1; xv(:,1)=0
else !!================ physical IC ================!!
# ifdef SUBV
    write(rank_s,'(i6)') rank_global
# else
    write(rank_s,'(i6)') rank
# endif
  rank_s=adjustl(rank_s)
  if (restart_ic) then
    if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
    call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
    write(z_s,'(f7.3)') z_write
    if (rank==0) print*, 'restarting simulation from z = ',z_s
  elseif (restart_kill) then
    z_s=reskill_prefix
    if (rank==0) print*, 'rescuing simulation from z = ',z_s
  else !! normally
    write(z_s,'(f7.3)') z_i
    if (rank==0) print*, 'starting simulation normally from z = ',z_s
  endif
  
  t1 = mpi_wtime(ierr)

  z_s=adjustl(z_s)

# ifdef ZIP
    !! Reading ZIPed dm
    f_zip0=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip0_'//trim(rank_s)//'.dat'
    f_zip1=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip1_'//trim(rank_s)//'.dat'
    f_zip2=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip2_'//trim(rank_s)//'.dat'
    f_zip3=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip3_'//trim(rank_s)//'.dat'
    open(10, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
    open(11, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
    open(12, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
    open(13, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')
    if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
      write(*,*) 'error opening dm zip checkpoint'
      write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
      call system('hostname')
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    read(10) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset
    read(11) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
    !print*, 'headers=',np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset; pause
    if (np_local > max_np) then
      write(*,*) 'too many particles to store'
      write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    np_uzip=0 !! reset
    do k = 1, nc_node_dim
    do j = 1, nc_node_dim
    do i = 1, nc_node_dim
      rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#     ifdef BGQ
        read(12) rhoc_i1(4)
#     else
        read(12) rhoc_i1(1)
#     endif
      if (rr_i4==255) read(13) rr_i4
      do l = 1, rr_i4
        xi1=0; xi4=0 !! Clean up
        np_uzip = np_uzip + 1
#       ifdef BGQ
          read(10) xi1(4,:)
#       else
          read(10) xi1(1,:)
#       endif
        read(11) vi2
        xv(1:3, np_uzip) = mesh_scale*((xi4+0.5)/256.+ (/i,j,k/)-1)
        xv(4:6, np_uzip) = vi2/v_r2i
      enddo
    enddo
    enddo
    enddo 
    read(10, end=701) test_i1 !! Consistency checks
      print*, 'ERROR: rank', rank, ': file not ending:', f_zip0                                                                            
      call mpi_abort(mpi_comm_world, ierr, ierr)
701 close(10)
    read(11, end=801) test_i1
      print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
      call mpi_abort(mpi_comm_world, ierr, ierr)
801 close(11);close(12);close(13)
    if (np_uzip /= np_local) then
      write(*,*) "Inconsistency! rank, np_local, np_uzip=",rank, np_local, np_uzip
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#   ifdef NEUTRINOS
      !! Reading ZIPed nu
      f_zip0=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip0_'//trim(rank_s)//'_nu.dat'
      f_zip1=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip1_'//trim(rank_s)//'_nu.dat'
      f_zip2=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip2_'//trim(rank_s)//'_nu.dat'
      f_zip3=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip3_'//trim(rank_s)//'_nu.dat'
      open(10, file=f_zip0, status="old", iostat=fstat0, access="stream", buffered='yes')
      open(11, file=f_zip1, status="old", iostat=fstat1, access="stream", buffered='yes')
      open(12, file=f_zip2, status="old", iostat=fstat2, access="stream", buffered='yes')
      open(13, file=f_zip3, status="old", iostat=fstat3, access="stream", buffered='yes')
      if (fstat0 /= 0 .or. fstat1 /= 0 .or. fstat2 /= 0 .or. fstat3 /= 0) then
        write(*,*) 'error opening dm zip checkpoint'
        write(*,*) 'rank',rank,'files:',f_zip0,f_zip1,f_zip2,f_zip3
        call system('hostname')
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(10) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,cur_checkpoint,cur_projection,cur_halofind,dummy,v_r2i_nu,dummy,dummy,dummy
      read(11) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
      if (np_local+np_nu > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'np_nu',np_nu,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      do k = 1, nc_node_dim
      do j = 1, nc_node_dim
      do i = 1, nc_node_dim
        rhoc_i1=0; rr_i4=0 !! Clean up, very imortant.
#       ifdef BGQ
          read(12) rhoc_i1(4)
#       else
          read(12) rhoc_i1(1)
#       endif
        if (rr_i4==255) read(13) rr_i4
        do l = 1, rr_i4
          xi1=0; xi4=0
          np_uzip = np_uzip + 1
#         ifdef BGQ
            read(10) xi1(4,:)
#         else
            read(10) xi1(1,:)
#         endif
          read(11) vi2
          xv(1:3, np_uzip) = mesh_scale * ( (xi4+0.5)/256.+(/i,j,k/)-1)
          xv(4:6, np_uzip) = vi2 / v_r2i_nu
        enddo
      enddo
      enddo
      enddo
      read(10, end=702) test_i1 !! Consistency checks
        print*, 'ERROR: rank', rank, ': file not ending:', f_zip0
        call mpi_abort(mpi_comm_world, ierr, ierr)
702   close(10)
      read(11, end=802) test_i1
        print*, 'ERROR: rank', rank, ': file not ending:', f_zip1
        call mpi_abort(mpi_comm_world, ierr, ierr)
802   close(11); close(12); close(13)
      if (np_uzip /= np_local+np_nu) then
        write(*,*) "Inconsistency! rank, np_local, np_uzip=", rank, np_local, np_nu, np_uzip
        call mpi_abort(mpi_comm_world,ierr,ierr)                                                                                         
      endif
#   endif !! NEUTRINOS

# else !! ndef ZIP

    !! read in ICs by traditional xv's
    ofile=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'xv'//trim(rank_s)//'.dat'
    if (rank==0) print*, 'opening dark matter checkpoint file:',ofile
    open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint'
      write(*,*) 'rank',rank,'file:',ofile
      call system('hostname')
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p
    if (np_local > max_np) then
      write(*,*) 'too many particles to store'
      write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    blocksize=(32*1024*1024)/24 !! reduced to 32MB chunks for intel compiler
    num_writes=np_local/blocksize+1
    do i = 1,num_writes
      nplow=(i-1)*blocksize+1; nphigh=min(i*blocksize,np_local)
      read(21) xv(:,nplow:nphigh)
    enddo
    close(21)
#   ifdef NEUTRINOS
      ofile=output_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'xv'//trim(rank_s)//'_nu.dat'
      if (rank==0) print*, 'opening neutrino checkpoint file:',ofile
      open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call system('hostname')
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(21) np_nu,a,t,tau,dummy,dummy,dummy,dummy,cur_checkpoint,cur_projection,cur_halofind,dummy
      if (np_local+np_nu > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      blocksize=(32*1024*1024)/24
      num_writes=np_nu/blocksize+1
      do i = 1,num_writes
        nplow=(i-1)*blocksize+1 + np_local; nphigh=min(i*blocksize,np_nu) + np_local
        read(21) xv(:,nplow:nphigh)
      enddo
      close(21)
#   endif !! NEUTRINOS

# endif !! ZIP
t2 = mpi_wtime(ierr)

ierr = hostnm(myhost)
#ifdef SUBV
  print*, 'rank:',rank,' rank_global:',rank_global, trim(myhost), (t2-t1)/60.,'min'
#else
  print*, 'done reading particles by rank:',rank, trim(myhost), (t2-t1)/60.,' min'
#endif

call mpi_barrier(mpi_comm_world,ierr)
endif

npl8=int(np_local,kind=8) !! calculate total number of particles and particle mass
call mpi_reduce(npl8,np_total,1,MPI_INTEGER8,mpi_sum,0,mpi_comm_world,ierr)
call mpi_bcast(np_total,1,MPI_INTEGER8,0,mpi_comm_world,ierr)
if (rank == 0) write(*,*) 'number of particles =', np_total
if (.not.restart_ic) then
  if (pairwise_ic) then
    mass_p=10000.0
  elseif (pair_infall) then
    mass_p=pair_infall_mass
  elseif (pp_test) then
    mass_p=10000.0/4.
  else
    mass_p = real(nf_physical_dim)**3 / real(np_total)
  endif
endif

    if (rank == 0) write(*,*) 'particle mass=', mass_p
    if (rank == 0) write(*,*) 'total dark matter mass =', mass_p * np_total

#ifdef NEUTRINOS
# ifndef NUPID
    PID(:np_local)=1
    PID(np_local+1:np_local+np_nu)=2
# endif
  npl8=int(np_nu,kind=8)
  call mpi_reduce(npl8,np_total_nu,1,MPI_INTEGER8, mpi_sum,0,mpi_comm_world,ierr)
  !! Append np_nu to np_local (must be done after PIDs and after mass calculation)
  np_dm = np_local; np_local = np_local + np_nu
  if (rank == 0) then
    write(*,*) "np_dm = ", np_dm
    write(*,*) "np_nu = ", np_nu
    write(*,*) "np_local = ", np_local
    write(*,*) "np_dm_total = ", np_total
    write(*,*) "np_nu_total = ", np_total_nu
  endif
#endif !! NEUTRINOS

!! This is to scale the initial conditions if we are doing testing with another data-set
#ifdef SCALED_IC
    do i=1,np_local
      xv(:,i)=xv(:,i)/4.0
    enddo
#endif
! this is to test if we can reconstruct particle ID from drifting vs non-drifting populations
#ifdef X_DRIFT
    do i=1,np_local
      xv(4,i)=xv(4,i) + 10
    enddo
#endif

!! check to make sure no particles are out of bounds
!! this is not really necessary, as link_list can handle out of bounds
#ifdef CHECK_IP 
    if (.not.restart_ic) then
      xva=0.0
      do i=1,np_local
        xva=xva+xv(:,i)
        do j=1,3
          if (xv(j,i) < 0.0 .or. xv(j,i) >= real(nf_physical_node_dim,4)) then
            write(*,*) 'particle out of bounds'
            write(*,*) xv(:,i)
            call mpi_abort(mpi_comm_world,ierr,ierr)
          endif
        enddo
      enddo
      xva=xva/real(np_local)
      print *,rank,xva 
    endif
#endif

if (rank==0) then !! print some information on the screen
  print*,'particle_initialize done'
  print*,'headers=',np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,' /',v_r2i_nu,shake_offset
  print*,'Scale factor a =',a
  print*,'t =',t
  print*,'Current timestep nts =',nts
  print*,'dt information',dt_f_acc,dt_pp_acc,dt_c_acc
  print*,'Current ck,proj,halo =',cur_checkpoint,cur_projection,cur_halofind
  print*,'Particle mass mass_p =',mass_p
  print*,'v_r2i. v_r2i_nu =',v_r2i,v_r2i_nu
  print*,'Current shake_offset =',shake_offset
endif

end subroutine particle_initialize
