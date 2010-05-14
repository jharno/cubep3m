!! calculate coarse mesh density
  subroutine coarse_mass

#ifdef MHD
   use mpi_tvd_mhd
#endif

   implicit none

#ifdef PPINT
   include 'cubep3m.fh'
#else
   include 'cubepm.fh'
#endif

   integer(4) :: i,j,k,pp

   call system_clock(count=count_i)

   rho_c= 0.0 !- mass_p * (mesh_scale / 2)**3  

#ifdef MHD
   !add gas mass

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

! unfortunately even unrolling this didn't prevent shared accesses to rho_c
! unroll it more.  do every third. Still didn't work
   !!$omp parallel do default(shared) private(i,j,k,pp)
   do k = 0, nc_node_dim + 1 !,2
     do j = 0, nc_node_dim + 1
       do i = 0, nc_node_dim + 1
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
   !!$omp end parallel do

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
