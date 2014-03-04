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

   integer(4) :: i,j,k,k0,pp,ii,jj,kk
   integer(4), dimension(3) :: i1,i2
   real(4), dimension(3) :: x,dx1,dx2


   call system_clock(count=count_i)

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

    do k0 = 0, mesh_scale-1 
        !$omp parallel do schedule(dynamic) default(shared) private(i,j,k,pp)
        do k = k0, nc_node_dim + 1, mesh_scale 
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
        !$omp end parallel do
    enddo

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
