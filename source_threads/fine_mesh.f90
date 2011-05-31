!! calculate fine mesh density, as well as fine
!! mesh velocity update #ifdef MHD
#ifdef MHD
  subroutine fine_mesh(tile,cmax,nerr,thread)
    use mpi_tvd_mhd
#else
  subroutine fine_mesh(tile,thread)
#endif
    implicit none

    include 'cubepm.fh'
    integer(4) :: thread
    integer(4) :: pp, i,j,k,ii,im,i2
    integer(4), dimension(3) :: cic_l, cic_h, tile

#ifdef MHD
    real(4) :: cmax
    integer(4) :: nerr
#endif

   if (thread==1) call system_clock(count=count_i)

!! normalize fine mesh density
    rho_f(:,:,:,thread)= 0.0

#ifdef MHD
!    print *,rank,'gas mass init rho_f:',sum(rho_f)
    call fine_gas_mass(tile,thread)
!    print *,rank,'gas mass rho_f:',sum(rho_f)
    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('fmg mass',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0 .and. thread==1) write(*,*) 'fine mesh gas mass finished',real(count_f-count_i)/real(count_r)
#endif
    call system_clock(count=count_i)
#endif

#ifdef NGP
    cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1
#else
    cic_l(:) = nc_tile_dim * tile(:) + 1 - nc_buf
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf
#endif

!! calculate fine mesh density for tile
    do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
        do i = cic_l(1), cic_h(1)
          pp=hoc(i,j,k)
#ifdef NGP
          call fine_ngp_mass(pp,tile,thread)
#else
          if (i == cic_l(1) .or. i == cic_h(1) .or. &
              j == cic_l(2) .or. j == cic_h(2) .or. &
              k == cic_l(3) .or. k == cic_h(3)) then
            call fine_cic_mass_boundry(pp,tile,thread)
          else
            call fine_cic_mass(pp,tile,thread)
          endif
#endif
        enddo
      enddo
    enddo

    
    if (thread==1) call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
! this need to be fixed for threaded version 
!    call mpi_time_analyze('fm  mass',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0 .and. thread==1) write(*,*) 'fine mesh mass finished',real(count_f-count_i)/real(count_r)
#endif


#ifdef DIAG 
    do k=1+nf_buf,nf_tile-nf_buf
      do j=1+nf_buf,nf_tile-nf_buf
        do i=1+nf_buf,nf_tile-nf_buf
          f_mesh_mass(thread)=f_mesh_mass(thread)+real(rho_f(i,j,k,thread),kind=8)
        enddo
      enddo
    enddo
#endif

#ifdef RHO_FBOUND_DEBUG
   if (rank==0) then
     open(52,file='ffmass.dat',form='formatted',status='replace')
     do i=1,nf_tile
       write(52,*) rho_f(i,i,i)
     enddo 
     close(52) 
   endif
#endif
#ifdef BOUND_DEBUG
  if (rank==0) then
    do k=1,nf_tile
       do j=1,nf_tile
         do i=1,nf_tile
           if (rho_f(i,j,k)>1.0) then
             print *,i,j,k,rho_f(i,j,k)
           endif
         enddo
       enddo
     enddo
   endif
#endif
#ifdef BOUND_DEBUG
  if (rank==0) then
    do k=1,nf_tile
       do j=1,nf_tile
         do i=1,nf_tile
           if (rho_f(i,j,k)/=1.0) then
              print *,i,j,k,rho_f(i,j,k)
              stop
           endif
         enddo
       enddo
     enddo
   endif
#endif

   if (thread == 1) call system_clock(count=count_i)

!! transform and calculate fine mesh force 
    !$omp critical 
    call cubepm_fftw2('f',thread)
    !$omp end critical 

#ifdef DEBUG
    print *,'rank',rank,'finished first fft'
#endif
    cmplx_rho_f(:,:,:,thread)=rho_f(:,:,:,thread)

    do i2=1,3
      do k = 1, nf_tile
        do j = 1, nf_tile
          do i = 1, nf_tile/2+1
            ii=2*i
            im=ii-1
            rho_f(im,j,k,thread)=-cmplx_rho_f(ii,j,k,thread)*kern_f(i2,i,j,k)
            rho_f(ii,j,k,thread)=cmplx_rho_f(im,j,k,thread)*kern_f(i2,i,j,k)
          enddo
        enddo
      enddo

#ifdef DEBUG
    print *,'rank',rank,'thread',thread,'finished convolve'
#endif

      !$omp critical 
      call cubepm_fftw2('b',thread)
      !$omp end critical 

#ifdef DEBUG
    print *,'rank',rank,'thread',thread,'finished second fft'
#endif

      force_f(i2,:,:,:,thread) = rho_f(nf_buf-1:nf_tile-nf_buf+1,nf_buf-1:nf_tile-nf_buf+1, &
                               nf_buf-1:nf_tile-nf_buf+1,thread)
    enddo

    if (thread == 1) call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
! see above
!    call mpi_time_analyze('fm force',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0 .and. thread==1) write(*,*) 'finished fine mesh force',real(count_f-count_i)/real(count_r)
#endif

    if (thread == 1) call system_clock(count=count_i)

#ifdef MHD
    call fine_velocity(tile,cmax,nerr)
#else
    call fine_velocity(tile,thread)
#endif
    
    if (thread == 1) call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
! see above
!    call mpi_time_analyze('fm   vel',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0 .and. thread==1) write(*,*) 'finished fine mesh velocity',real(count_f-count_i)/real(count_r)
#endif


  end subroutine fine_mesh
