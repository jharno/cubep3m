!! main particle mesh subroutine
  subroutine particle_mesh

    implicit none

    include 'mpif.h'
    include 'cubepm.fh'

    integer(4) :: i, j, k, cur_tile, thread
    integer(4), dimension(3) :: tile
    real(4), dimension(3) :: offset
    integer(4) :: k0, i3
    integer(4), dimension(3) :: cic_l, cic_h 
    integer(4), dimension(3, cores) :: cic_fine_l, cic_fine_h
    real(4) :: force_mag
    real(4) :: f_force_max_node
    real(4) :: pp_force_max_node
    real(4) :: pp_ext_force_max_node
    integer(4) :: omp_get_thread_num
    external omp_get_thread_num
#ifdef NESTED_OMP
    real(4), dimension(nested_threads) :: pp_ext_sum
#else
    real(4), dimension(1) :: pp_ext_sum
#endif
#ifdef MHD
    integer(4) :: nerrl, nerr
    real(4) :: cmaxl, cmax
#endif
#ifdef DIAG
    real(8) :: sumrhof
#endif

    !! Variables that appear in the nested parts (give them all _n suffix) 
    integer(4) :: i_n, j_n, k_n, thread_n, n_pairs_n 
    real(4), dimension(3) :: x_n, offset_n, dx1_n, dx2_n
    integer(4) :: pp_n, pp1_n, pp2_n
    integer(4), dimension(3) :: i1_n, i2_n
    integer(4) :: ii_n, jj_n, kk_n, im_n, jm_n, km_n, ip_n, jp_n, kp_n
    integer(4) :: ip_min_n, ip_max_n, jp_min_n, jp_max_n, kp_min_n, kp_max_n
    integer(4) :: ipl_n(mesh_scale, mesh_scale, mesh_scale)
    real(4), dimension(3) :: sep_n, force_pp_n
    real(4) :: pp_force_mag_n, force_mag_n, rmag_n, dVc_n, pp_ext_sum_k_n

!! start of particle mesh.  All particles are within (1:nc_node_dim]

    if (pairwise_ic) then
      call set_pair
    elseif (pair_infall) then
      if (nts==1) then
        call set_pair
      endif
      call update_position
    else
      call update_position
    endif

!! particles must not have advanced past hoc_nc_l:hoc_nc_h

    call link_list

    call particle_pass

#ifdef MHD
    nerr=0
    cmax=1e-5
#endif

#ifdef NESTED_OMP
    !$omp  parallel num_threads(cores) default(shared) &
    !$omp& private(cur_tile, i, j, k, k0, tile, thread, i3, cic_l, cic_h, offset, force_mag, pp_ext_sum)
#else
    !$omp  parallel num_threads(cores) default(shared) &
    !$omp& private(cur_tile, i, j, k, k0, tile, thread, i3, cic_l, cic_h, offset, force_mag, pp_ext_sum, &
    !$omp&         thread_n, pp_n, pp1_n, pp2_n, ii_n, im_n, jm_n, km_n, x_n, dx1_n, dx2_n, dVc_n, i2_n, i1_n, offset_n, &
    !$omp&         ip_n, jp_n, kp_n, ipl_n, sep_n, force_pp_n, rmag_n, force_mag_n, pp_force_mag_n, n_pairs_n, &
    !$omp&         pp_ext_sum_k_n, i_n, j_n, k_n, ip_min_n, ip_max_n, jp_min_n, jp_max_n, kp_min_n, kp_max_n, jj_n, kk_n)
#endif
    thread = 1
    thread = omp_get_thread_num() + 1
    f_mesh_mass(thread)=0.0
    f_force_max(thread, :)=0.0
#ifdef PPINT
    pp_force_max(thread, :)=0.0
#endif
    !$omp do schedule(dynamic) 
    do cur_tile=1,tiles_node
      tile(3) = (cur_tile-1) / (tiles_node_dim * tiles_node_dim)
      j = cur_tile - tile(3) * tiles_node_dim * tiles_node_dim
      tile(2) = (j-1) /  tiles_node_dim
      j = j - tile(2) * tiles_node_dim
      tile(1) = j - 1

!#ifdef MHD
!      call fine_mesh(tile,cmax,nerr,thread)
!#else
!      call fine_mesh(tile,thread)
!#endif

!! normalize fine mesh density
     rho_f(:,:,:,thread)= 0.0

#ifdef MHD
#ifdef DEBUG_PP_MESH
!    print *,rank,'gas mass init rho_f:',sum(rho_f)
    print *,'Entering fine_gas_mass', cur_tile, thread, rank
#endif

    call fine_gas_mass(tile,thread) 
    rho_f(:,:,:,thread) = rho_f(:,:,:,thread)*(omega_b/omega_m) 

#ifdef DEBUG_PP_MESH
    print *,cur_tile ,thread,rank,'gas mass rho_f:',sum(rho_f)
#endif

#endif


!! calculate coarse mesh offsets
#ifdef NGP
      cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
      cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1
#else
      cic_l(:) = nc_tile_dim * tile(:) + 1 - nc_buf
      cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf
#endif
!! calculate fine mesh density for tile
#ifdef NESTED_OMP
    do k0 = 0, mesh_scale-1 
        !$omp parallel num_threads(nested_threads) default(shared) private(i_n, j_n, k_n, pp_n, offset_n, x_n, i1_n) 
        !$omp do schedule(dynamic)
        do k_n = cic_l(3) + k0, cic_h(3), mesh_scale
#else
        do k_n = cic_l(3), cic_h(3)
#endif
            do j_n = cic_l(2), cic_h(2)
                do i_n = cic_l(1), cic_h(1)
                    pp_n = hoc(i_n, j_n, k_n)
#ifdef NGP
!                   call fine_ngp_mass(pp_n,tile,thread)

                    offset_n(:)= - tile(:) * nf_physical_tile_dim + nf_buf
! removed the half-cell offset so that fine mesh cells will line up with coarse mesh cells
!    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf - 0.5

                    do; if (pp_n == 0) exit
              
                        x_n = xv(1:3, pp_n) + offset_n(:)
                        i1_n(:) = floor(x_n(:)) + 1

#ifdef MHD
                        rho_f(i1_n(1),i1_n(2),i1_n(3),thread) = rho_f(i1_n(1),i1_n(2),i1_n(3),thread)+mass_p*(1.0-omega_b/omega_m)
#else
                        rho_f(i1_n(1),i1_n(2),i1_n(3),thread) = rho_f(i1_n(1),i1_n(2),i1_n(3),thread)+mass_p
#endif
                        pp_n = ll(pp_n)

                    enddo !! pp_n

#else
                    if (i_n == cic_l(1) .or. i_n == cic_h(1) .or. &
                        j_n == cic_l(2) .or. j_n == cic_h(2) .or. &
                        k_n == cic_l(3) .or. k_n == cic_h(3) ) then
                            call fine_cic_mass_boundry(pp_n, tile, thread)
                    else
                            call fine_cic_mass(pp_n, tile, thread)
                    endif
#endif
                enddo !! i_n
            enddo !! j_n
        enddo !! k_n
#ifdef NESTED_OMP
        !$omp end do
        !$omp end parallel
    enddo !! k0
#endif

! sum total fine mesh mass
#ifdef DIAG 
      do k=1+nf_buf,nf_tile-nf_buf
        do j=1+nf_buf,nf_tile-nf_buf
          do i=1+nf_buf,nf_tile-nf_buf
            f_mesh_mass(thread)=f_mesh_mass(thread)+real(rho_f(i,j,k,thread),kind=8)
          enddo
        enddo
      enddo
#endif
!! transform and calculate fine mesh force 
      call cubepm_fftw2('f',thread)

#ifdef DEBUG_PP_MESH
      print *,'(tile,thread,rank) = (',cur_tile,thread,rank,') finished first fft'
#endif
      cmplx_rho_f(:,:,:,thread)=rho_f(:,:,:,thread)

      do i3 = 1, 3
#ifdef NESTED_OMP
        !$omp parallel num_threads(nested_threads) default(shared) private(k_n, j_n, i_n, ii_n, im_n)
        !$omp do
#endif
        do k_n = 1, nf_tile
          do j_n = 1, nf_tile
            do i_n = 1, nf_tile/2+1
              ii_n = 2*i_n
              im_n = ii_n - 1
              rho_f(im_n, j_n, k_n, thread) = -cmplx_rho_f(ii_n, j_n, k_n, thread) * kern_f(i3, i_n, j_n, k_n)
              rho_f(ii_n, j_n, k_n, thread) =  cmplx_rho_f(im_n, j_n, k_n, thread) * kern_f(i3, i_n, j_n, k_n)
            enddo
          enddo
        enddo
#ifdef NESTED_OMP
        !$omp end do 
        !$omp end parallel
#endif

#ifdef DEBUG_PP_MESH
        print *,'(tile,thread, rank) = (',cur_tile, thread, rank,') finished convolve'
#endif
        call cubepm_fftw2('b',thread)
#ifdef DEBUG_PP_MESH
        print *,'(tile,thread,rank) = *',cur_tile,thread,rank,') finished second fft'
#endif

#ifdef NESTED_OMP
        !$omp parallel num_threads(nested_threads) default(shared) private(i_n, j_n, k_n)
        !$omp do
        do k_n = nf_buf-1, nf_tile-nf_buf+1
            do j_n = nf_buf-1, nf_tile-nf_buf+1
                do i_n = nf_buf-1, nf_tile-nf_buf+1
                    force_f(i3, i_n, j_n, k_n, thread) = rho_f(i_n, j_n, k_n, thread)
                enddo
            enddo
        enddo
        !$omp end do 
        !$omp end parallel
#else
        force_f(i3,:,:,:,thread) = rho_f(nf_buf-1:nf_tile-nf_buf+1,nf_buf-1:nf_tile-nf_buf+1, &
                               nf_buf-1:nf_tile-nf_buf+1,thread)
#endif

      enddo
! fine velocity
   
!! calculate max fine mesh force   
      if (.not.pair_infall) then
        force_mag=0.0
        do k=nf_buf-1,nf_tile-nf_buf+1
          do j=nf_buf-1,nf_tile-nf_buf+1
            do i=nf_buf-1,nf_tile-nf_buf+1
              force_mag=sqrt(force_f(1,i,j,k,thread)**2 &
                  +force_f(2,i,j,k,thread)**2+force_f(3,i,j,k,thread)**2)
              if (force_mag > maxval(f_force_max(thread, :))) f_force_max(thread, 1)=force_mag
            enddo
          enddo
        enddo
      endif

!! update dark matter velocity

      offset(:) = real(nf_buf) - tile(:) * nf_physical_tile_dim

!      print *,'thread',thread,'offset',offset,'tile',tile

! removed the half cell offset so that fine mesh cells will line up with coarse mesh cells
!    offset(:) = - 0.5 + real(nf_buf) - tile(:) * nf_physical_tile_dim

#ifdef NESTED_OMP
    do k0 = 0, mesh_scale-1 
      !$omp  parallel num_threads(nested_threads) default(shared) &
      !$omp& private(i_n, j_n, k_n, pp_n, x_n, i1_n, im_n, jm_n, km_n, i2_n, force_mag_n, dx1_n, dx2_n, dVc_n, &
      !$omp&         ipl_n, ip_n, jp_n, kp_n, pp1_n, pp2_n, sep_n, rmag_n, force_pp_n, pp_force_mag_n, thread_n)
      thread_n = 1
      thread_n = omp_get_thread_num() + 1
      !$omp do
      do k_n = k0 + tile(3) * nc_tile_dim + 1, (tile(3) + 1) * nc_tile_dim, mesh_scale 
#else
      thread_n = 1
      do k_n = tile(3) * nc_tile_dim + 1, (tile(3) + 1) * nc_tile_dim
#endif
        do j_n = tile(2) * nc_tile_dim + 1, (tile(2) + 1) * nc_tile_dim
          do i_n = tile(1) * nc_tile_dim + 1, (tile(1) + 1) * nc_tile_dim
            pp_n = hoc(i_n, j_n, k_n)
#ifdef PPINT
#ifdef DEBUG_PP_MESH_INTENSE
            if (pp_n /= 0) print *, pp_n, i_n, j_n, k_n
#endif
            ipl_n = 0
#endif
            do; if (pp_n == 0) exit
              x_n = xv(1:3, pp_n) + offset(:)
              i1_n(:) = floor(x_n(:)) + 1

!********************************************************
#ifdef PPINT

              do im_n = 1, 3
                i2_n(im_n) = mod(i1_n(im_n)-1,mesh_scale) + 1
              enddo
              ipl_n(i2_n(1),i2_n(2),i2_n(3)) = ipl_n(i2_n(1),i2_n(2),i2_n(3))+1
              if (ipl_n(i2_n(1),i2_n(2),i2_n(3))>max_llf) then
                print *,'exceeded max_llf',max_llf,i1_n,i2_n,ipl_n
                stop
              endif
              llf(ipl_n(i2_n(1),i2_n(2),i2_n(3)),i2_n(1),i2_n(2),i2_n(3),thread,thread_n)=pp_n
#endif
!********************************************************

#ifdef NGP
              if (pp_test) print *,'before ngp',pp_n,xv(:,pp_n)
#ifdef DEBUG_PP_MESH_INTENSE
              print *,'force',i1_n,force_f(:,i1_n(1),i1_n(2),i1_n(3),thread)
#endif
              if (ngp_fmesh_force) xv(4:6,pp_n)=xv(4:6,pp_n)+force_f(1:3,i1_n(1),i1_n(2),i1_n(3),thread) * &
                         a_mid * G * dt
              if (pair_infall) then
                force_mag_n=sqrt(force_f(1,i1_n(1),i1_n(2),i1_n(3),thread)**2+force_f(2,i1_n(1),i1_n(2),i1_n(3),thread)**2+ &
                               force_f(3,i1_n(1),i1_n(2),i1_n(3),thread)**2)
                if (force_mag_n > f_force_max(thread, thread_n)) f_force_max(thread, thread_n)=force_mag_n
              endif
              if (pp_test) print *,'before pp',pp_n,xv(:,pp_n)

#else
              i2_n(:)  = i1_n(:) + 1
              dx1_n(:) = i1_n(:) - x_n(:)
              dx2_n(:) = 1.0 - dx1_n(:)

              dVc_n = a_mid * G * dt * dx1_n(1) * dx1_n(2) * dx1_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i1_n(1),i1_n(2),i1_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx2_n(1) * dx1_n(2) * dx1_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i2_n(1),i1_n(2),i1_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx1_n(1) * dx2_n(2) * dx1_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i1_n(1),i2_n(2),i1_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx2_n(1) * dx2_n(2) * dx1_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i2_n(1),i2_n(2),i1_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx1_n(1) * dx1_n(2) * dx2_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i1_n(1),i1_n(2),i2_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx2_n(1) * dx1_n(2) * dx2_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i2_n(1),i1_n(2),i2_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx1_n(1) * dx2_n(2) * dx2_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i1_n(1),i2_n(2),i2_n(3),thread) * dVc_n
              dVc_n = a_mid * G * dt * dx2_n(1) * dx2_n(2) * dx2_n(3)
              xv(4:6,pp_n) = xv(4:6,pp_n) &
                         + force_f(1:3,i2_n(1),i2_n(2),i2_n(3),thread) * dVc_n
#endif
              pp_n = ll(pp_n)
            enddo

!***********
#ifdef PPINT
!***********
            do km_n = 1, mesh_scale
              do jm_n = 1, mesh_scale
                do im_n = 1, mesh_scale
                  if (pp_test .and. ipl_n(im_n, jm_n, km_n) /= 0) print *,'ipl_n', im_n, jm_n, km_n, ipl_n(im_n, jm_n, km_n)
#ifdef DEBUG_PP_MESH
                  if ( ipl_n(im_n, jm_n, km_n) > 1) print *,'ipl', rank, i_n, j_n, k_n, im_n, jm_n, km_n, ipl_n(im_n, jm_n, km_n)
#endif
                  pp_force_accum(:, :ipl_n(im_n,jm_n,km_n), thread, thread_n) = 0.
                  do ip_n = 1, ipl_n(im_n,jm_n,km_n) - 1
                    pp1_n = llf(ip_n, im_n, jm_n, km_n, thread, thread_n)
                    do jp_n = ip_n+1, ipl_n(im_n,jm_n,km_n)
                      pp2_n = llf(jp_n, im_n, jm_n, km_n, thread, thread_n)
                      sep_n = xv(:3,pp1_n) - xv(:3,pp2_n)                      
                      rmag_n = sqrt(sep_n(1)*sep_n(1) + sep_n(2)*sep_n(2) + sep_n(3)*sep_n(3))
                      if (rmag_n > rsoft) then
#ifdef MHD
                        force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)*(1.0 - omega_b/omega_m)
#else          
                        force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)  !mass_p divides out below
#endif
                        pp_force_accum(:, ip_n, thread, thread_n) = pp_force_accum(:, ip_n, thread, thread_n) - force_pp_n
                        pp_force_accum(:, jp_n, thread, thread_n) = pp_force_accum(:, jp_n, thread, thread_n) + force_pp_n
                        if (pp_force_flag) then
                          xv(4:,pp1_n) = xv(4:,pp1_n) - force_pp_n*a_mid*G*dt
                          xv(4:,pp2_n) = xv(4:,pp2_n) + force_pp_n*a_mid*G*dt
                        endif
                      endif
                    enddo
                  enddo
                  do ip_n=1, ipl_n(im_n,jm_n,km_n)
                    pp_force_mag_n = sqrt(pp_force_accum(1,ip_n,thread,thread_n)**2 + pp_force_accum(2,ip_n,thread,thread_n)**2 + pp_force_accum(3,ip_n,thread,thread_n)**2)
                    if (pp_force_mag_n > pp_force_max(thread, thread_n)) pp_force_max(thread, thread_n)=pp_force_mag_n
                  enddo
                enddo
              enddo
            enddo
#endif
!****************
! end of pp force
!****************
          enddo
        enddo
      enddo
#ifdef NESTED_OMP 
      !$omp end do
      !$omp end parallel
    enddo
#endif 

! end fine velocity on dm

#ifdef MHD
    !write(*,*)  'Calling fine_velocity for MHD'
    !! NOTE: THIS SUBROUTINE WILL NEED TO BE MODIFIED FOR NESTED_THREADS 
    call fine_velocity(tile,cmax,nerr,thread)
    write(*,*)  'Called fine_velocity for MHD on tile', cur_tile
#endif


#ifdef PP_EXT

!*****************
!Extended pp force
!*****************

!**********************************
! 1-Create link list on fine mesh *
!**********************************

#ifdef DEBUG_PP_EXT
      write(*,*) 'Creating link_list on fine fine mesh for :'
      write(*,*) 'tile',tile, 'thread', thread
#endif
 
      hoc_fine(:,:,:,thread)=0
 
      ! Get physical coordinate of boundaries of tiles in fine mesh units
      cic_fine_l(:,thread) = tile(:)*nf_physical_tile_dim + 1
      cic_fine_h(:,thread) = (tile(:) + 1)*nf_physical_tile_dim
      
      ! Include pp_range
      cic_fine_l(:,thread) = cic_fine_l(:,thread) - pp_range 
      cic_fine_h(:,thread) = cic_fine_h(:,thread) + pp_range 
      
#ifdef DEBUG_PP_EXT
      write(*,*) 'cic_fine_l =', cic_fine_l(:,thread) 
      write(*,*) 'cic_fine_h =', cic_fine_h(:,thread) 
#endif
            
#ifdef NESTED_OMP
    !$omp parallel num_threads(nested_threads) default(shared) private(i_n, j_n, k_n, ii_n, jj_n, kk_n, pp_n)
    !$omp do schedule(dynamic)
#endif
    do k_n = cic_l(3), cic_h(3)
        do j_n = cic_l(2), cic_h(2)
            do i_n = cic_l(1), cic_h(1)

                pp_n = hoc(i_n, j_n, k_n)

                do; if (pp_n == 0) exit
                    ii_n = floor(xv(1,pp_n))+1
                    jj_n = floor(xv(2,pp_n))+1
                    kk_n = floor(xv(3,pp_n))+1

                    if (ii_n >= cic_fine_l(1,thread) .and. ii_n <= cic_fine_h(1,thread) .and. &
                        jj_n >= cic_fine_l(2,thread) .and. jj_n <= cic_fine_h(2,thread) .and. &
                        kk_n >= cic_fine_l(3,thread) .and. kk_n <= cic_fine_h(3,thread)) then

                        ll_fine(pp_n,thread) = hoc_fine(ii_n-cic_fine_l(1,thread)+1, jj_n-cic_fine_l(2,thread)+1, kk_n-cic_fine_l(3,thread)+1,thread)
                        hoc_fine(ii_n-cic_fine_l(1,thread)+1, jj_n-cic_fine_l(2,thread)+1, kk_n-cic_fine_l(3,thread)+1,thread)=pp_n

                    endif

                    pp_n = ll(pp_n)
                enddo

            enddo
        enddo
    enddo 
#ifdef NESTED_OMP
    !$omp end do
    !$omp end parallel
#endif

#ifdef DEBUG_LINK_LIST_FINE
     
      write(*,*) 'hoc_fine(1,1,1,',thread,') =',hoc_fine(1,1,1,thread) 
      write(*,*) 'hoc_fine(2,1,1,',thread,') =',hoc_fine(2,1,1,thread) 
      write(*,*) 'hoc_fine(3,1,1,',thread,') =',hoc_fine(3,1,1,thread) 
      write(*,*) 'hoc_fine(1,2,1,',thread,') =',hoc_fine(1,2,1,thread) 
      write(*,*) 'hoc_fine(2,2,1,',thread,') =',hoc_fine(2,2,1,thread) 
      write(*,*) 'hoc_fine(3,2,1,',thread,') =',hoc_fine(3,2,1,thread) 
      write(*,*) 'hoc_fine(1,3,1,',thread,') =',hoc_fine(1,3,1,thread) 
      write(*,*) 'hoc_fine(2,3,1,',thread,') =',hoc_fine(2,3,1,thread) 
      write(*,*) 'hoc_fine(3,3,1,',thread,') =',hoc_fine(3,3,1,thread) 
      write(*,*) 'hoc_fine(1,1,2,',thread,') =',hoc_fine(1,1,2,thread) 
      write(*,*) 'hoc_fine(2,1,2,',thread,') =',hoc_fine(2,1,2,thread) 
      write(*,*) 'hoc_fine(3,1,2,',thread,') =',hoc_fine(3,1,2,thread) 
      write(*,*) 'hoc_fine(1,2,2,',thread,') =',hoc_fine(1,2,2,thread) 
      write(*,*) 'hoc_fine(2,2,2,',thread,') =',hoc_fine(2,2,2,thread) 
      write(*,*) 'hoc_fine(3,2,2,',thread,') =',hoc_fine(3,2,2,thread) 
      write(*,*) 'hoc_fine(1,3,2,',thread,') =',hoc_fine(1,3,2,thread) 
      write(*,*) 'hoc_fine(2,3,2,',thread,') =',hoc_fine(2,3,2,thread) 
      write(*,*) 'hoc_fine(3,3,2,',thread,') =',hoc_fine(3,3,2,thread) 
      write(*,*) 'hoc_fine(1,1,3,',thread,') =',hoc_fine(1,1,3,thread) 
      write(*,*) 'hoc_fine(2,1,3,',thread,') =',hoc_fine(2,1,3,thread) 
      write(*,*) 'hoc_fine(3,1,3,',thread,') =',hoc_fine(3,1,3,thread) 
      write(*,*) 'hoc_fine(1,2,3,',thread,') =',hoc_fine(1,2,3,thread) 
      write(*,*) 'hoc_fine(2,2,3,',thread,') =',hoc_fine(2,2,3,thread) 
      write(*,*) 'hoc_fine(3,2,3,',thread,') =',hoc_fine(3,2,3,thread) 
      write(*,*) 'hoc_fine(1,3,3,',thread,') =',hoc_fine(1,3,3,thread) 
      write(*,*) 'hoc_fine(2,3,3,',thread,') =',hoc_fine(2,3,3,thread) 
      write(*,*) 'hoc_fine(3,3,3,',thread,') =',hoc_fine(3,3,3,thread) 
      
            
      ! loop over fine cells in the tile, including the pp_range
      do k=1,nf_physical_tile_dim+2*pp_range
         do j=1,nf_physical_tile_dim+2*pp_range
            do i=1,nf_physical_tile_dim+2*pp_range
               
               write(*,*) 'hoc_fine(',i,j,k,thread,') = ',hoc_fine(i,j,k,thread) 
               pause
               
            enddo
         enddo
      enddo
      
#endif

!*******************************************************
! 2-Loop over pp in present and neighbouring fine cells*
!*******************************************************

#ifdef NESTED_OMP
    !$omp parallel num_threads(nested_threads) default(shared) private(k_n)
    !$omp do
#endif
    do k_n = 1, max_np 
      pp_ext_force_accum(:,k_n,thread) = 0 
    enddo
#ifdef NESTED_OMP
    !$omp end do
    !$omp end parallel
#endif

    if(pp_range .ne. 0)then

#ifdef NESTED_OMP
        do k0 = 0, pp_range
            !$omp  parallel num_threads(nested_threads) default(shared) & 
            !$omp& private(i_n, j_n, k_n, pp1_n, kp_min_n, kp_max_n, kp_n, jp_min_n, jp_max_n, &
            !$omp&         jp_n, ip_min_n, ip_max_n, ip_n, pp2_n, n_pairs_n, sep_n, rmag_n, force_pp_n) 
            !$omp do
            do k_n=1+k0,nf_physical_tile_dim+pp_range,pp_range+1 ! We never loop towards smaller z
#else
            do k_n=1,nf_physical_tile_dim+pp_range ! We never loop towards smaller z
#endif
                do j_n=1,nf_physical_tile_dim+2*pp_range
                    do i_n=1,nf_physical_tile_dim+2*pp_range
               
                        pp1_n = hoc_fine(i_n, j_n, k_n, thread) 
                        if(pp1_n == 0) cycle 
 
                        kp_min_n = k_n
                        kp_max_n = k_n + pp_range
                  
                        do kp_n = kp_min_n, kp_max_n

                            if(kp_n == k_n) then
                                jp_min_n = j_n
                            else
                                jp_min_n = j_n - pp_range
                                if(jp_min_n <= 0) jp_min_n = 1
                            endif
                            jp_max_n = j_n + pp_range
                            if(jp_max_n > nf_physical_tile_dim+2*pp_range) jp_max_n = nf_physical_tile_dim+2*pp_range

                            do jp_n = jp_min_n, jp_max_n 

                                if((kp_n == k_n .and. jp_n == j_n))then
                                    ip_min_n = i_n+1
                                else
                                    ip_min_n = i_n - pp_range
                                    if(ip_min_n <= 0) ip_min_n = 1
                                endif
                                ip_max_n = i_n + pp_range
                                if(ip_max_n > nf_physical_tile_dim+2*pp_range) ip_max_n = nf_physical_tile_dim+2*pp_range
                        
                                do ip_n = ip_min_n, ip_max_n 
                           
                                    pp2_n = hoc_fine(ip_n, jp_n, kp_n, thread) 
                                    if(pp2_n == 0) cycle                                                     

#ifdef DEBUG_PP_EXT
                                    write(*,*) 'Found a pair of cells with sep :',ip_n-i_n,jp_n-j_n,kp_n-k_n,'on thread', thread
                                    write(*,*) 'at position (i,j,k)=',i_n,j_n,k_n
                                    write(*,*) 'on tile :',tile
                                    write(*,*) 'with hoc_fine:', pp1_n, pp2_n
                                    write(*,*) 'xv(hoc1)=',xv(:3,pp1_n)
                                    write(*,*) 'xv(hoc2)=',xv(:3,pp2_n)
                                    write(*,*) 'll_fine=',ll_fine(1:20,thread)                           
#endif
                                    n_pairs_n = 0
                        
                                    do; if(pp1_n == 0)exit
                           
                                        do ; if(pp2_n == 0)exit

                                            n_pairs_n = n_pairs_n + 1
                           
                                            !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                                            sep_n = xv(:3,pp1_n) - xv(:3,pp2_n)
                                            rmag_n = sqrt(sep_n(1)*sep_n(1) + sep_n(2)*sep_n(2) + sep_n(3)*sep_n(3))

                                            if (rmag_n > rsoft) then

                                                if(rmag_n>real(nf_cutoff)+sqrt(3.0))then
                                                    force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)
                                                else
                                                    force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)*(1 - (7.0/4.0)*(rmag_n*pp_bias/(nf_cutoff))**3 + &
                                                                 (3.0/4.0)*(rmag_n*pp_bias/(nf_cutoff))**5)  !mass_p divides out below
                                                endif
#ifdef MHD
                                                force_pp_n = force_pp_n*(1.0 - omega_b/omega_m)                
#endif

                                                !force_pp = force_pp - mass_p*( -7*rmag/(4*nf_cutoff**3) + 3*rmag**3/(4*nf_cutoff**5))
                                                !force_pp = force_pp + sep*mass_p*(7/(4*nf_cutoff**3) - 3*rmag**2/(4*nf_cutoff**5))
                                                pp_ext_force_accum(:, pp1_n, thread) = pp_ext_force_accum(:, pp1_n, thread) - force_pp_n
                                                pp_ext_force_accum(:, pp2_n, thread) = pp_ext_force_accum(:, pp2_n, thread) + force_pp_n
                              
                                                if (pp_ext_force_flag) then
                                 
                                                    ! Update only particles in physical space
                                                    if((pp_range<i_n).and.(i_n<=nf_physical_tile_dim+pp_range) .and.&
                                                       (pp_range<j_n).and.(j_n<=nf_physical_tile_dim+pp_range).and.&
                                                       (pp_range<k_n).and.(k_n<=nf_physical_tile_dim+pp_range)) then 
                                    
                                                        !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                                                        xv(4:,pp1_n) = xv(4:,pp1_n) - force_pp_n*a_mid*G*dt
                                                    endif
                                 
                                                    if((pp_range<ip_n).and.(ip_n<=nf_physical_tile_dim+pp_range) .and.&
                                                       (pp_range<jp_n).and.(jp_n<=nf_physical_tile_dim+pp_range).and.&
                                                       (pp_range<kp_n).and.(kp_n<=nf_physical_tile_dim+pp_range)) then 

                                                        !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                                                        xv(4:,pp2_n)=xv(4:,pp2_n) + force_pp_n*a_mid*G*dt
                                                    endif
                             
                                                endif !! pp_ext_force_flag
                           
                                            endif !! (rmag_n > rsoft)
                           
                                            !HERE, I SHOULD TRY TO AVOID READING THE ll_fine VARIABLE TO INCREASE PROCESSOR SPEED
                                            pp2_n = ll_fine(pp2_n, thread)                           
                           
                                        enddo !! pp2_n

                                        pp2_n = hoc_fine(ip_n, jp_n, kp_n, thread)

                                        !HERE, I SHOULD TRY TO AVOID READING THE ll_fine VARIABLE TO INCREASE PROCESSOR SPEED
                                        pp1_n = ll_fine(pp1_n, thread)
                        
                                    enddo !! pp1_n
#ifdef DEBUG_PP_EXT
                                    write(*,*) 'n_pairs in that cell couple =', n_pairs_n
#endif
                                    ! Restore pp1_n value for the next pp2_n iteration
                                    pp1_n = hoc_fine(i_n, j_n, k_n, thread)
                        
                                enddo !! ip_n                    
                            enddo !! jp_n
                        enddo !! kp_n                  
                    enddo !! i_n
                enddo !! j_n
            enddo !! k_n
#ifdef NESTED_OMP
            !$omp end do
            !$omp end parallel
        enddo !! k0
#endif

     endif !! (pp_range .ne. 0)

     pp_ext_sum(:) = 0.
#ifdef NESTED_OMP
     !$omp parallel num_threads(nested_threads) default(shared) private(k_n, pp_ext_sum_k_n, thread_n)
     thread_n = 1
     thread_n = omp_get_thread_num() + 1
     !$omp do
#else
     thread_n = 1
#endif
     do k_n = 1, max_np
        pp_ext_sum_k_n = pp_ext_force_accum(1,k_n,thread)**2 + pp_ext_force_accum(2,k_n,thread)**2 + pp_ext_force_accum(3,k_n,thread)**2
        if (pp_ext_sum_k_n > pp_ext_sum(thread_n)) pp_ext_sum(thread_n) = pp_ext_sum_k_n
     enddo
#ifdef NESTED_OMP
     !$omp end do
     !$omp end parallel
#endif

     pp_ext_force_max(thread) = sqrt(maxval(pp_ext_sum(:)))
     
#ifdef DEBUG_PP_EXT
      write(*,*) 'pp_ext_force_max(',thread,') =', pp_ext_force_max(thread)
#endif
      
#endif
      
   enddo
   !$omp end do
   !$omp end parallel
   
#ifdef MHD
    cmaxl=cmax
    nerrl=nerr
    call mpi_reduce(cmaxl,cmax,1,mpi_real,mpi_max,0,mpi_comm_cart,ierr)
    call mpi_reduce(nerrl,nerr,1,mpi_integer,mpi_sum,0,mpi_comm_cart,ierr)

    if (rank == 0) then
      print *,'fluid stats',cmax/freeze,dt*cmax,nerr
    endif
#endif

!! calculate maximum dt from fine mesh force

    f_force_max_node=maxval(f_force_max)

    call mpi_reduce(f_force_max_node,dt_f_acc,1,mpi_real,mpi_max,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
! I"m pretty sure this is incorrect!!! (no wonder this always seemed wrong :-/ )
!      dt_f_acc=1.0/sqrt(min(0.0001,dt_f_acc)*a_mid*G)
      dt_f_acc=1.0/sqrt(max(0.0001,dt_f_acc)*a_mid*G)
      write(*,*) 'maximum timestep from fine force=',dt_f_acc
    endif

    call mpi_bcast(dt_f_acc,1,mpi_real,0,mpi_comm_world,ierr)

#ifdef PPINT

!! calculate maximum dt from particle-particle force
    
    pp_force_max_node=maxval(pp_force_max)

    call mpi_reduce(pp_force_max_node,dt_pp_acc,1,mpi_real,mpi_max,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
      dt_pp_acc=sqrt(dt_pp_scale*rsoft)/max(sqrt(dt_pp_acc*a_mid*G),1e-3)
      write(*,*) 'maximum timestep from pp force=',dt_pp_acc
    endif
   

    call mpi_bcast(dt_pp_acc,1,mpi_real,0,mpi_comm_world,ierr)

    if (pp_test) then
      do i=1,np_local
        print *,i,xv(:,i)
      enddo 
    endif

#endif

#ifdef PP_EXT

    pp_ext_force_max_node=maxval(pp_ext_force_max)

    call mpi_reduce(pp_ext_force_max_node,dt_pp_ext_acc,1,mpi_real,mpi_max,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
      !write(*,*) dt_pp_scale, rsoft,dt_pp_ext_acc,a_mid,G, max(sqrt(dt_pp_ext_acc*a_mid*G),1e-3)
      dt_pp_ext_acc=sqrt(dt_pp_scale*rsoft)/max(sqrt(dt_pp_ext_acc*a_mid*G),1e-3)
      write(*,*) 'maximum timestep from pp ext force=',dt_pp_ext_acc
    endif
   
    call mpi_bcast(dt_pp_ext_acc,1,mpi_real,0,mpi_comm_world,ierr)

#endif

!! calculate mass of fine mesh

#ifdef DIAG
    call mpi_reduce(sum(f_mesh_mass),sumrhof,1,mpi_double_precision, &
                    mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'sum of rho_f=',sumrhof
#endif

#ifdef DEBUG_VEL
    write(*,*) rank,xv(:,1:np_local)
#endif

    call coarse_mesh

! undo the random shift
#ifdef MOVE_GRID_BACK
    call move_grid_back
#endif    

!! delete all particles outside (1:nc_node_dim]

    call delete_particles

    if (pairwise_ic.or.pair_infall) then
      call report_pair
    endif

  end subroutine particle_mesh
