!! calculate coarse mesh density
  subroutine coarse_mass

#ifdef MHD
   use mpi_tvd_mhd
#endif

   implicit none

    include 'mpif.h'
#   include "cubepm.fh"

   integer(4) :: i,j,k,k0,pp,ii,jj,kk
   integer(4), dimension(3) :: i1,i2
   real(4), dimension(3) :: x,dx1,dx2

#ifdef COARSEPROJ
    character (len=max_path) :: ofile,oofile
    character (len=6) :: step_string
    character (len=6) :: rank_s
    character (len=200) :: cmd
    integer :: fstat
    real(8) :: cpsum_local, cpsum
    integer, parameter :: ijstart = -1
    integer, parameter :: ijstop  = nc_node_dim + 2
#else
    integer, parameter :: ijstart = 0
    integer, parameter :: ijstop  = nc_node_dim + 1 
#endif
    integer(4) :: hostnm
    character(len=100) :: myhost
   call system_clock(count=count_i)
   ierr = hostnm(myhost)
   rho_c= 0.0 !- mass_p * (mesh_scale / 2)**3  

#ifdef MHD
   !add gas mass

#ifdef MHD_CIC_MASS
    do k=4,nf_physical_node_dim-4
      !kfc=((k-1)/mesh_scale)+1
      !ku=nz%m+k-1
      do j=4,nf_physical_node_dim-4
       ! jfc=((j-1)/mesh_scale)+1
        !ju=ny%m+j-1
        do i=4,nf_physical_node_dim-4
        !  ifc=((i-1)/mesh_scale)+1
        !  iu=nx%m+i-1
          x(1) = (1.0/real(mesh_scale)) * i - 0.5
          x(2) = (1.0/real(mesh_scale)) * j - 0.5
          x(3) = (1.0/real(mesh_scale)) * k - 0.5
          i1(:) = floor(x(:)) + 1
          i2(:) = i1(:) + 1
          dx1(:) = i1(:) - x(:)
          dx2(:) = 1.0 - dx1(:)
          rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) &
                               + dx1(1) * dx1(2) * dx1(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) &
                               + dx2(1) * dx1(2) * dx1(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) &
                               + dx1(1) * dx2(2) * dx1(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) &
                               + dx2(1) * dx2(2) * dx1(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) &
                               + dx1(1) * dx1(2) * dx2(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) &
                               + dx2(1) * dx1(2) * dx2(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) &
                               + dx1(1) * dx2(2) * dx2(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) &
                               + dx2(1) * dx2(2) * dx2(3) * u(1,nx%m+(i-1),ny%m+(j-1),nz%m+(k-1))
          enddo
        enddo
      enddo

#else

   do k=1,nc_node_dim
     do j=1,nc_node_dim
       do i=1,nc_node_dim
         rho_c(i,j,k)=sum( u(1,nx%m+(i-1)*mesh_scale : nx%m+i*mesh_scale-1, &
                               ny%m+(j-1)*mesh_scale : ny%m+j*mesh_scale-1, &
                               nz%m+(k-1)*mesh_scale : nz%m+k*mesh_scale-1))
       enddo
     enddo
   enddo
#endif

   rho_c = rho_c*(omega_b/omega_m)

#endif

#ifdef COARSEPROJ
    !! Determine if want to do projection here
    doCoarseProj = .false.
    if (mod(nts,writeCoarseProjEverySteps) == 0 .and. a >= writeCoarseProjAboveA) then 
        doCoarseProj = .true.
        !! Initialize grids
        crhoproj = 0. ; crhoprojsum = 0.        
#ifdef NEUTRINOS
        crhoproj_nu = 0. ; crhoprojsum_nu = 0.
#endif
        !! Completely remove z offset while binning x, y to nearest coarse mesh cell
        call mpi_bcast(shake_offset, 3, mpi_real, 0, mpi_comm_world, ierr)
        soffcproj(1) = (shake_offset(1) - mesh_scale*int(shake_offset(1)/mesh_scale))/mesh_scale 
        soffcproj(2) = (shake_offset(2) - mesh_scale*int(shake_offset(2)/mesh_scale))/mesh_scale
        soffcproj(3) = shake_offset(3)/mesh_scale

        if (rank == 2) write(*,*) "shake_offset = ", shake_offset
        if (rank == 2) write(*,*) "soffcproj = ", soffcproj

    endif
#endif 

    do k0 = 0, mesh_scale-1 
        !$omp parallel do schedule(dynamic) default(shared) private(i,j,k,pp)
        do k = k0, nc_node_dim + 1, mesh_scale 
            do j = ijstart, ijstop 
                do i = ijstart, ijstop 
                    pp=hoc(i,j,k)
                    if (i <= 1 .or. i >= nc_node_dim .or. &
                        j <= 1 .or. j >= nc_node_dim .or. &
                        k <= 1 .or. k >= nc_node_dim) then
                      call coarse_cic_mass_boundry(pp)
                    else
                        call coarse_cic_mass(pp)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    enddo

#ifdef COARSEPROJ
    if (doCoarseProj) then
        !! Sum up grid as consistency check
        cpsum_local = sum(crhoproj(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim))
        call mpi_reduce(cpsum_local, cpsum, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr) 
        if (rank == 0) write(*,*) "Dark matter coarsened projection sum = ", cpsum
# ifdef NEUTRINOS 
        cpsum_local = sum(crhoproj_nu(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim))
        call mpi_reduce(cpsum_local, cpsum, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
        if (rank == 0) write(*,*) "Neutrino coarsened projection sum = ", cpsum
# endif
        !! Write data using only those on the top xy slab of the box
        if (cart_coords(1) == 0) then
            !! Rank and time step strings
            write(rank_s,'(i6)') rank
            rank_s = adjustl(rank_s)
            write(step_string, "(i6.6)") nts
            !! Sum up grids into projections
            do k = nc_coarseproj_start, nc_coarseproj_stop 
                do j = 1, nc_node_dim
                    do i = 1, nc_node_dim
                        crhoprojsum(i, j) = crhoprojsum(i,j) + crhoproj(i,j,k) 
# ifdef NEUTRINOS
                        crhoprojsum_nu(i,j) = crhoprojsum_nu(i,j) + crhoproj_nu(i,j,k)
# endif
                    enddo
                enddo
            enddo
            !! Write dark matter projection
            oofile=shm_path//'cproj'//trim(rank_s)//'_'//step_string//'.dat'
            ofile=output_path//'coarseproj/node'//trim(rank_s)//'/cproj'//trim(rank_s)//'_'//step_string//'.dat'
#           ifdef SHMCP
              open(unit=12, file=oofile,status="replace",iostat=fstat,access="stream")
#           else
              open(unit=12, file=ofile,status="replace",iostat=fstat,access="stream")
#           endif
            if (fstat /= 0) then
              write(*,*) 'error opening coarse projection file for write'
#             ifdef SHMCP
                write(*,*) trim(myhost),' rank',rank,'file:',oofile
#             else
                write(*,*) trim(myhost),' rank',rank,'file:',ofile
#             endif
              call mpi_abort(mpi_comm_world,ierr,ierr)
            endif
            write(12) shake_offset
            write(12) soffcproj
            write(12) crhoprojsum 
            close(12)
# ifdef NEUTRINOS
            !! Write neutrino projection
            oofile=shm_path//'cproj'//trim(rank_s)//'_'//step_string//'_nu.dat'
            ofile=output_path//'coarseproj/node'//trim(rank_s)//'/cproj'//trim(rank_s)//'_'//step_string//'_nu.dat'
#           ifdef SHMCP
              open(unit=13,file=oofile,status="replace",iostat=fstat,access="stream")
#           else
              open(unit=13,file=ofile,status="replace",iostat=fstat,access="stream")
#           endif
            if (fstat /= 0) then
              write(*,*) 'error opening coarse projection file for write'
#             ifdef SHMCP
                write(*,*) trim(myhost),' rank',rank,'file:',oofile
#             else
                write(*,*) trim(myhost),' rank',rank,'file:',ofile
#             endif
              call mpi_abort(mpi_comm_world,ierr,ierr)
            endif
            write(13) shake_offset
            write(13) soffcproj
            write(13) crhoprojsum_nu
            close(13)
# endif
        endif !! cart_coord
        if (rank == 0) write(*,*) "Coarse Projection written"
# ifdef SHMCP
    cmd='mv '//trim(shm_path)//'*cproj* '//trim(output_path)//'coarseproj/node'//trim(rank_s)//'/ &'
    call system(cmd)
# endif
        call mpi_barrier(mpi_comm_world, ierr)
    endif !! doCoarseProj
#endif

#ifdef DEBUG_VEL
   do k=1,nc_node_dim
     do j=1,nc_node_dim
       do i=1,nc_node_dim
         if (rho_c(i,j,k) /= 0.0 ) write(*,*) 'rhoc',rank,i,j,k,rho_c(i,j,k)
       enddo
     enddo
   enddo
#endif

  call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
  call mpi_time_analyze('cm  mass',real(count_f-count_i)/real(count_r),rank,nodes)
#else
  if (rank==0) write(*,*) 'coarse mass finished',real(count_f-count_i)/real(count_r)
#endif


  end subroutine coarse_mass
