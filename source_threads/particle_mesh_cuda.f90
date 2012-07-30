!! main particle mesh subroutine
  subroutine particle_mesh
    implicit none
    include 'mpif.h'
    include 'cubepm.fh'

    integer(4) :: i,j,k,dd,cur_tile, np_tile, np_tile_buf, n_pairs
    integer(4), dimension(3) :: tile
    integer(4) :: thread
    real(4) :: f_force_max_node
    real(4) :: pp_force_max_node
    real(4) :: pp_ext_force_max_node
    integer(4) :: omp_get_thread_num
    external omp_get_thread_num

! these are for fine mesh
    integer(4) :: pp,ii,im,i3
    integer(4), dimension(3) :: cic_l, cic_h 
    integer(4), dimension(3,cores) :: cic_fine_l, cic_fine_h

! these are for fine ngp mass
    integer(4), dimension(3) :: i1
    real(4),    dimension(3) :: x, offset
    real(8),    dimension(3) :: tmpx
! these are for fine velocity
    real(4), dimension(3) :: dx1, dx2
    real(4) :: dVc
    integer(4), dimension(3) :: i2
    integer(4) :: jm,km,ip,jp,kp, ip_min, jp_min, kp_min,ip_max, jp_max, kp_max
    real(4) :: force_mag
!#ifdef PPINT
    integer pp1,pp2,ipl(mesh_scale,mesh_scale,mesh_scale)
    real sep(3), force_pp(3), rmag, pp_force_mag, v_init(3) ! pp_force_accum(3)
!#endif
 
#ifdef pp_ext_on_GPU
    real(4), dimension(3,10000,cores):: x1_gpu, x2_gpu, f1
    integer(4), dimension(10000,cores):: pp1_ll 
#endif
    integer(4) n1,n2

#ifdef MHD
    integer(4) :: nerrl,nerr
    real(4) :: cmaxl,cmax
#endif

#ifdef DIAG
    real(8) :: sumrhof
#endif

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

!!!!  IS FFTW2 THREADSAFE???  NEED TO TEST !!!!

    !$omp parallel default(shared) &
    !$omp& private(cur_tile,i,j,k,tile,thread,pp,ii,im,i3,cic_l,cic_h,i1,x,tmpx, &
    !$omp&              offset,dx1,dx2,dVc,i2,jm,km,ip,jp,kp,force_mag,pp1,pp2,ipl,& 
    !$omp&              sep,force_pp,rmag,pp_force_mag,v_init, np_tile, &
    !$omp&              ip_min, jp_min, kp_min,ip_max, jp_max, kp_max, n_pairs,n1)
    thread=1
    thread = omp_get_thread_num() + 1
    f_mesh_mass(thread)=0.0
    f_force_max(thread)=0.0
#ifdef PPINT
    pp_force_max(thread)=0.0
#endif
    !$omp do 
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
      do k = cic_l(3), cic_h(3)
        do j = cic_l(2), cic_h(2)
          do i = cic_l(1), cic_h(1)
            pp=hoc(i,j,k)
#ifdef NGP
!            call fine_ngp_mass(pp,tile,thread)

            offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf
! removed the half-cell offset so that fine mesh cells will line up with coarse mesh cells
!    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf - 0.5

            do; if (pp == 0) exit
              !x(:) = xv(1:3,pp) + offset(:)
              tmpx(:) = real(xv(1:3,pp), kind=8) + real(offset(:), kind=8)
              x = real(tmpx,kind=4)
              i1(:) = floor(tmpx(:)) + 1
              !i1(:) = floor(x(:)) + 1

#ifdef MHD
              rho_f(i1(1),i1(2),i1(3),thread) = rho_f(i1(1),i1(2),i1(3),thread)+mass_p*(1.0-omega_b/omega_m)
#else
              rho_f(i1(1),i1(2),i1(3),thread) = rho_f(i1(1),i1(2),i1(3),thread)+mass_p
#endif
              pp = ll(pp)
            enddo

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

      do i3=1,3
        do k = 1, nf_tile
          do j = 1, nf_tile
            do i = 1, nf_tile/2+1
              ii=2*i
              im=ii-1
              rho_f(im,j,k,thread)=-cmplx_rho_f(ii,j,k,thread)*kern_f(i3,i,j,k)
              rho_f(ii,j,k,thread)=cmplx_rho_f(im,j,k,thread)*kern_f(i3,i,j,k)
            enddo
          enddo
        enddo

#ifdef DEBUG_PP_MESH
        print *,'(tile,thread, rank) = (',cur_tile, thread, rank,') finished convolve'
#endif
        call cubepm_fftw2('b',thread)
#ifdef DEBUG_PP_MESH
        print *,'(tile,thread,rank) = *',cur_tile,thread,rank,') finished second fft'
#endif

        force_f(i3,:,:,:,thread) = rho_f(nf_buf-1:nf_tile-nf_buf+1,nf_buf-1:nf_tile-nf_buf+1, &
                               nf_buf-1:nf_tile-nf_buf+1,thread)
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
              if (force_mag > f_force_max(thread)) f_force_max(thread)=force_mag
            enddo
          enddo
        enddo
      endif

!! update dark matter velocity

      offset(:) = real(nf_buf) - tile(:) * nf_physical_tile_dim

!      print *,'thread',thread,'offset',offset,'tile',tile

! removed the half cell offset so that fine mesh cells will line up with coarse mesh cells
!    offset(:) = - 0.5 + real(nf_buf) - tile(:) * nf_physical_tile_dim

      do k = tile(3) * nc_tile_dim + 1, (tile(3) + 1) * nc_tile_dim
        do j = tile(2) * nc_tile_dim + 1, (tile(2) + 1) * nc_tile_dim
          do i = tile(1) * nc_tile_dim + 1, (tile(1) + 1) * nc_tile_dim
            pp = hoc(i,j,k)
#ifdef PPINT
#ifdef DEBUG_PP_MESH_INTENSE
            if (pp /= 0) print *,pp,i,j,k
#endif
            ipl=0
#endif
            do; if (pp == 0) exit
              tmpx(:) = real(xv(1:3,pp), kind=8) + real(offset(:), kind=8)
              x = real(tmpx,kind=4)
              i1(:) = floor(tmpx(:)) + 1
              !x(:) = xv(1:3,pp) + offset(:)
              !i1(:) = floor(x(:)) + 1

! arghh!!

!              do dd=1,3
!                if (i1(dd) < nf_buf-1 .or. i1(dd) > nf_tile-nf_buf+1) then
!                  print *,'out',thread,i1,x,pp,i,j,k,hoc(i,j,k),xv(:,pp)            
!                  stop
!                endif
!              enddo

#ifdef NGP
              if (pp_test) print *,'before ngp',pp,xv(:,pp)
#ifdef DEBUG_PP_MESH_INTENSE
              print *,'force',i1,force_f(:,i1(1),i1(2),i1(3),thread)
#endif
              if (ngp_fmesh_force) xv(4:6,pp)=xv(4:6,pp)+force_f(1:3,i1(1),i1(2),i1(3),thread) * &
                         a_mid * G * dt
              if (pair_infall) then
                force_mag=sqrt(force_f(1,i1(1),i1(2),i1(3),thread)**2+force_f(2,i1(1),i1(2),i1(3),thread)**2+ &
                               force_f(3,i1(1),i1(2),i1(3),thread)**2)
                if (force_mag > f_force_max(thread)) f_force_max(thread)=force_mag
              endif
              if (pp_test) print *,'before pp',pp,xv(:,pp)
!********************************************************
#ifdef PPINT

              do im=1,3
                i2(im)=mod(i1(im)-1,mesh_scale)+1
              enddo
              ipl(i2(1),i2(2),i2(3))=ipl(i2(1),i2(2),i2(3))+1
              if (ipl(i2(1),i2(2),i2(3))>max_llf) then
                print *,'exceeded max_llf',max_llf,i1,i2,ipl
                stop
              endif
              llf(ipl(i2(1),i2(2),i2(3)),i2(1),i2(2),i2(3),thread)=pp
#endif
!********************************************************

#else
              i2(:) = i1(:) + 1
              dx1(:) = i1(:) - x(:)
              dx2(:) = 1.0 - dx1(:)

              dVc = a_mid * G * dt * dx1(1) * dx1(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i1(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx1(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i1(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx1(1) * dx2(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i2(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx2(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i2(2),i1(3),thread) * dVc
              dVc = a_mid * G * dt * dx1(1) * dx1(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i1(2),i2(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx1(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i1(2),i2(3),thread) * dVc
              dVc = a_mid * G * dt * dx1(1) * dx2(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i1(1),i2(2),i2(3),thread) * dVc
              dVc = a_mid * G * dt * dx2(1) * dx2(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) &
                         + force_f(1:3,i2(1),i2(2),i2(3),thread) * dVc
#endif
              pp = ll(pp)
            enddo

!***********
#ifdef PPINT
!***********
            do km=1,mesh_scale
              do jm=1,mesh_scale
                do im=1,mesh_scale
                  if (pp_test .and. ipl(im,jm,km) /= 0) print *,'ipl',im,jm,km,ipl(im,jm,km)
#ifdef DEBUG_PP_MESH
                  if ( ipl(im,jm,km) > 1) print *,'ipl',rank,i,j,k,im,jm,km,ipl(im,jm,km)
#endif
                  pp_force_accum(:,:ipl(im,jm,km),thread)=0.0
                  do ip=1,ipl(im,jm,km)-1
                    pp1=llf(ip,im,jm,km,thread)
                    do jp=ip+1,ipl(im,jm,km)
                      pp2=llf(jp,im,jm,km,thread)
                      sep=xv(:3,pp1)-xv(:3,pp2)                      
                      !write(*,*) 'xv1=', xv(:,pp1)
                      !write(*,*) 'xv2=', xv(:,pp2)                      
                      rmag=sqrt(sep(1)*sep(1)+sep(2)*sep(2)+sep(3)*sep(3))
                      if (rmag>rsoft) then
#ifdef MHD
              force_pp=mass_p*(sep/(rmag*pp_bias)**3)*(1.0 - omega_b/omega_m)
#else          
              force_pp=mass_p*(sep/(rmag*pp_bias)**3)  !mass_p divides out below
#endif
                        pp_force_accum(:,ip,thread)=pp_force_accum(:,ip,thread)-force_pp
                        pp_force_accum(:,jp,thread)=pp_force_accum(:,jp,thread)+force_pp
                        if (pp_force_flag) then
                          xv(4:,pp1)=xv(4:,pp1)-force_pp*a_mid*G*dt
                          xv(4:,pp2)=xv(4:,pp2)+force_pp*a_mid*G*dt
                        endif
                      endif
                    enddo
                  enddo
                  do ip=1,ipl(im,jm,km)
                    pp_force_mag=sqrt(pp_force_accum(1,ip,thread)**2+pp_force_accum(2,ip,thread)**2+pp_force_accum(3,ip,thread)**2)
                    if (pp_force_mag>pp_force_max(thread)) pp_force_max(thread)=pp_force_mag
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
! end fine velocity on dm

#ifdef MHD
    !write(*,*)  'Calling fine_velocity for MHD'
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
      np_tile=0
 
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
      
      
      pp=1
      do; if (pp > np_local) exit
         i=floor(xv(1,pp))+1
         j=floor(xv(2,pp))+1
         k=floor(xv(3,pp))+1       
         
         if (i < cic_fine_l(1,thread) .or. i > cic_fine_h(1,thread) .or. &
              j < cic_fine_l(2,thread) .or. j > cic_fine_h(2,thread) .or. &
              k < cic_fine_l(3,thread) .or. k > cic_fine_h(3,thread)) then

            !write (*,*) 'PARTICLE NOT IN CURRENT TILE RANGE',xv(:,pp) 
            np_tile_buf=np_tile_buf+1
            !         goto 94
            !cycle
         else
            ll_fine(pp,thread) = hoc_fine(i-cic_fine_l(1,thread)+1, j-cic_fine_l(2,thread)+1, k-cic_fine_l(3,thread)+1,thread)
            hoc_fine(i-cic_fine_l(1,thread)+1, j-cic_fine_l(2,thread)+1, k-cic_fine_l(3,thread)+1,thread)=pp
            !ll_fine(pp,thread) = hoc_fine(i-cic_fine_l(1,thread)+1, j-cic_fine_l(2,thread)+1, k-cic_fine_l(3,thread)+1,thread)
            !hoc_fine(i-cic_fine_l(1,thread)+1, j-cic_fine_l(2,thread)+1, k-cic_fine_l(3,thread)+1)=pp
            
            !write(*,*)'hoc_fine(',i,j,k,thread,')= ',hoc_fine(i-cic_fine_l(1,thread)+1, j-cic_fine_l(2,thread)+1, k-cic_fine_l(3,thread)+1,thread)
            !write(*,*) 'pp =', pp
            !pause
            
            np_tile=np_tile+1
            
         endif
         pp=pp+1
      enddo
      
#ifdef DEBUG_LINK_LIST_FINE

      write(*,*) 'np_tile =', np_tile
      
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
#ifdef pp_ext_on_GPU 

      write(*,*)'trying gpu on cubep3m'
      
      pp_ext_force_accum(:,:,thread) = 0.0
      n_pairs = 0

      if(pp_range .ne. 0)then

!        write(*,*)'Entering first loop'
         do k=1,nf_physical_tile_dim+2*pp_range ! We do loop towards smaller z in GPU
            do j=1,nf_physical_tile_dim+2*pp_range
               do i=1,nf_physical_tile_dim+2*pp_range
                  

                  !write(*,*) 'make list of particles from group 1 in '
                  !write(*,*) 'i,j,k = ',i,j,k
                  pp1 = hoc_fine(i,j,k,thread)
                  n1 = 0
                  !pp1_ll(:,thread)  = 0
                  !write(*,*) 'pp1 = ',pp1
                  if(pp1 == 0) cycle
                  do; if(pp1==0)exit
                     
                     !write(*,*) 'Got a particle! ', pp1,xv(1,pp1),xv(2,pp1),xv(3,pp1)
                     n1 = n1+1 
                     x1_gpu(1,n1,thread) = xv(1,pp1)
                     x1_gpu(2,n1,thread) = xv(2,pp1)
                     x1_gpu(3,n1,thread) = xv(3,pp1)

                     ! Create a list of the particle position in the 
                     ! xv array to find it in the end
                     pp1_ll(n1,thread) = pp1

                     pp1 = ll_fine(pp1,thread)
                  enddo
                  !write(*,*) 'Got x1_gpu,n1 =' , n1                                    
                 
 
                  ! loop over groups 2, and for each,
                  ! first make a list, then calculate the force  on x1,
                  ! and update that force  
                  ! for gpu code, I only store the resulting force on f1,
                  ! so I can't use Newton's third law  and get the force on group 2.
                  ! That one must be recalculated. Hence the extended force loops 
                  ! over all fine cells around, with no worries about double counting.

                  n2=0

                  kp_min = k-pp_range
                  if(kp_min <=0) kp_min = 1
                  kp_max = k+pp_range
                  if(kp_max > nf_physical_tile_dim+2*pp_range) kp_max = nf_physical_tile_dim+2*pp_range
                  do kp = kp_min,kp_max                    
                     jp_min = j-pp_range
                     if(jp_min <=0) jp_min = 1
                     jp_max = j+pp_range
                     if(jp_max > nf_physical_tile_dim+2*pp_range) jp_max = nf_physical_tile_dim+2*pp_range
                     do jp = jp_min,jp_max
                        ip_min = i-pp_range
                        if(ip_min <=0) ip_min = 1
                        ip_max = i+pp_range
                        if(ip_max > nf_physical_tile_dim+2*pp_range) ip_max = nf_physical_tile_dim+2*pp_range
                        do ip=ip_min,ip_max
 
                           if(i.eq.ip .and. j.eq.jp .and. k.eq.kp)  cycle
 
                           pp2 = hoc_fine(ip,jp,kp,thread)
                           !pp1 = ll_fine(pp2,thread)
                           !n2 = 0
                           if(pp2 == 0) cycle
                           do; if(pp2==0)exit
                              n2 = n2+1
                              x2_gpu(1,n2,thread) = xv(1,pp2)
                              x2_gpu(2,n2,thread) = xv(2,pp2)
                              x2_gpu(3,n2,thread) = xv(3,pp2)
                              pp2 = ll_fine(pp2,thread) 
                           enddo
                        enddo
                     enddo
                  enddo
                  !write(*,*) 'Got x2_gpu, n2 = ' ,  n2
                  if (n2==0)cycle
                                             
                  ! use GPU packakage here:
                  f1(:,:,thread) = 0.0;
                  !write(*,*) 'Calling CUDA'
                  call pp_force_c(n1,x1_gpu(1,1:n1,thread),x1_gpu(2,1:n1,thread),x1_gpu(3,1:n1,thread),f1(1,1:n1,thread),f1(2,1:n1,thread),f1(3,1:n1,thread),n2,x2_gpu(1,1:n2,thread),x2_gpu(2,1:n2,thread),x2_gpu(3,1:n2,thread))
                  !write(*,*) 'pp_gpu :' ,x1_gpu(1,1,thread),x1_gpu(2,1,thread),x1_gpu(3,1,thread),f1(1,1,thread),f1(2,1,thread),f1(3,1,thread),x2_gpu(1,1,thread),x2_gpu(2,1,thread),x2_gpu(3,1,thread)

                  ! The output in the f1 vector is \vec{r}/r**3.
                  ! So (f dot f) = r**2/r**6 = r**(-4) 
                  !write(*,*) 'Entering ii loop'

                  do ii = 1,n1
                     if(f1(1,ii,thread)==0.0 .and. f1(2,ii,thread)==0 .and. f1(3,ii,thread)==0) then
                        write(*,*)'Not sure why I am here, but should skip!'
                        cycle
                     endif
                     rmag = (f1(1,ii,thread)*f1(1,ii,thread) + f1(2,ii,thread)*f1(2,ii,thread) + f1(3,ii,thread)*f1(3,ii,thread))**(-1.0/4.0)
                     !write(*,*) 'dr/r**3 = ',thread , f1(:,ii,thread), 'rmag = ', rmag                      

                     if (rmag>rsoft) then
                        if(rmag>real(nf_cutoff)+sqrt(3.0))then
                           force_pp=mass_p*f1(:,ii,thread)*(pp_bias)**(-3)
                        else
                           force_pp=mass_p*f1(:,ii,thread)*(pp_bias)**(-3)*(1 - (7.0/4.0)*(rmag*pp_bias/(nf_cutoff))**3 + &
                              (3.0/4.0)*(rmag*pp_bias/(nf_cutoff))**5)  !mass_p divides out below
                        endif
                        !write(*,*)'force_pp =',  force_pp

                        pp_ext_force_accum(:,pp1_ll(ii,thread),thread)=pp_ext_force_accum(:,pp1_ll(ii,thread),thread)+force_pp
                        !pp_ext_force_accum(:,pp2,thread)=pp_ext_force_accum(:,pp2,thread)+force_pp
                        if (pp_ext_force_flag) then
  
                           ! Update only particles in physical space
                           if((pp_range<i).and.(i<=nf_physical_tile_dim+pp_range) .and.&
                              (pp_range<j).and.(j<=nf_physical_tile_dim+pp_range) .and.&
                              (pp_range<k).and.(k<=nf_physical_tile_dim+pp_range)) then

                              !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                              xv(4:,pp1_ll(ii,thread))=xv(4:,pp1_ll(ii,thread))+force_pp*a_mid*G*dt
                           endif
  
                           !if((pp_range<ip).and.(ip<=nf_physical_tile_dim+pp_range) .and.&
                           !   (pp_range<jp).and.(jp<=nf_physical_tile_dim+pp_range).and.&
                           !   (pp_range<kp).and.(kp<=nf_physical_tile_dim+pp_range)) then
                           !
                           !   !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                           !   xv(4:,pp2)=xv(4:,pp2)+force_pp*a_mid*G*dt
                           !endif
                        endif
                     endif
                  enddo
                  !write(*,*)'Done reassigning loop  on  cell',  i,j,k
                  
               enddo
            enddo
         enddo
      endif
      write(*,*) 'Done extended pp loop'

      pp_ext_force_max(thread) = maxval(sqrt(pp_ext_force_accum(1,:,thread)**2 + pp_ext_force_accum(2,:,thread)**2 + pp_ext_force_accum(3,:,thread)**2))
      !pp_ext_force_max(thread) = 0.0001


#else
      pp_ext_force_accum(:,:,thread) = 0 
      n_pairs = 0

      if(pp_range .ne. 0)then

         do k=1,nf_physical_tile_dim+pp_range ! We never loop towards smaller z
            do j=1,nf_physical_tile_dim+2*pp_range
               do i=1,nf_physical_tile_dim+2*pp_range
               
                  pp1 = hoc_fine(i,j,k,thread) 
                  if(pp1 == 0) cycle 
 
                  kp_min = k
                  kp_max = k+pp_range
                  do kp = kp_min,kp_max
                     if(kp==k)then
                        jp_min = j
                     else
                        jp_min = j-pp_range
                        if(jp_min <=0) jp_min = 1
                     endif
                     jp_max = j+pp_range
                     if(jp_max > nf_physical_tile_dim+2*pp_range) jp_max = nf_physical_tile_dim+2*pp_range
                     do jp = jp_min,jp_max 
                        if((kp==k .and. jp==j))then
                           ip_min = i+1
                        else
                           ip_min = i-pp_range
                           if(ip_min <=0) ip_min = 1
                        endif
                        ip_max = i+pp_range
                        if(ip_max > nf_physical_tile_dim+2*pp_range) ip_max = nf_physical_tile_dim+2*pp_range
                        do ip=ip_min,ip_max 
                           
                           !write(*,*) '(i,j,k)= ', i,j,k
                           !write(*,*) '(ip,jp,kp)  =',  ip,jp,kp
                           !pause
                           
                           pp2 = hoc_fine(ip,jp,kp,thread) 
                           if(pp2 == 0) cycle                                                     
                              
 
#ifdef DEBUG_PP_EXT
                           write(*,*) 'Found a pair of cells with sep :',ip-i,jp-j,kp-k,'on thread', thread
                           write(*,*) 'at position (i,j,k)=',i,j,k
                           write(*,*) 'on tile :',tile
                           write(*,*) 'with hoc_fine:', pp1, pp2
                           write(*,*) 'xv(hoc1)=',xv(:3,pp1)
                           write(*,*) 'xv(hoc2)=',xv(:3,pp2)
                           write(*,*) 'll_fine=',ll_fine(1:20,thread)
                           !pause
                           
#endif
                           n_pairs = 0
                             
                           do; if(pp1==0)exit
                           do ; if(pp2==0)exit

                           n_pairs = n_pairs+1
                           
                           !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                           sep = xv(:3,pp1) - xv(:3,pp2)
                           rmag=sqrt(sep(1)*sep(1)+sep(2)*sep(2)+sep(3)*sep(3))

                           !write(*,*) 'n_pairs =',n_pairs
                           !write(*,*) 'sep =',sep
                           !write(*,*) 'rmag =',rmag

                           if (rmag>rsoft) then
                              if(rmag>real(nf_cutoff)+sqrt(3.0))then
                                 force_pp=mass_p*(sep/(rmag*pp_bias)**3)
                              else
                                 force_pp=mass_p*(sep/(rmag*pp_bias)**3)*(1 - (7.0/4.0)*(rmag*pp_bias/(nf_cutoff))**3 + &
                                      (3.0/4.0)*(rmag*pp_bias/(nf_cutoff))**5)  !mass_p divides out below
                              endif
#ifdef MHD
                              force_pp = force_pp*(1.0 - omega_b/omega_m)                
#endif

                              !force_pp = force_pp - mass_p*( -7*rmag/(4*nf_cutoff**3) + 3*rmag**3/(4*nf_cutoff**5))
                              !force_pp = force_pp + sep*mass_p*(7/(4*nf_cutoff**3) - 3*rmag**2/(4*nf_cutoff**5))
                              pp_ext_force_accum(:,pp1,thread)=pp_ext_force_accum(:,pp1,thread)-force_pp
                              pp_ext_force_accum(:,pp2,thread)=pp_ext_force_accum(:,pp2,thread)+force_pp
                              if (pp_ext_force_flag) then
                                 
                                 ! Update only particles in physical space
                                 if((pp_range<i).and.(i<=nf_physical_tile_dim+pp_range) .and.&
                                      (pp_range<j).and.(j<=nf_physical_tile_dim+pp_range).and.&
                                      (pp_range<k).and.(k<=nf_physical_tile_dim+pp_range)) then 
                                    
                                    !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                                    xv(4:,pp1)=xv(4:,pp1)-force_pp*a_mid*G*dt
                                 endif
                                 
                                 if((pp_range<ip).and.(ip<=nf_physical_tile_dim+pp_range) .and.&
                                      (pp_range<jp).and.(jp<=nf_physical_tile_dim+pp_range).and.&
                                      (pp_range<kp).and.(kp<=nf_physical_tile_dim+pp_range)) then 

                                    !HERE, I SHOULD TRY TO AVOID READING THE xv VARIABLE TO INCREASE PROCESSOR SPEED
                                    xv(4:,pp2)=xv(4:,pp2)+force_pp*a_mid*G*dt
                                 endif
                              endif
                           endif
                           
                                !HERE, I SHOULD TRY TO AVOID READING THE ll_fine VARIABLE TO INCREASE PROCESSOR SPEED
				pp2=ll_fine(pp2,thread)                           
                           	enddo

				pp2 = hoc_fine(ip,jp,kp,thread)

                                !HERE, I SHOULD TRY TO AVOID READING THE ll_fine VARIABLE TO INCREASE PROCESSOR SPEED
				pp1=ll_fine(pp1,thread)
				enddo
#ifdef DEBUG_PP_EXT
                                write(*,*) 'n_pairs in that cell couple =',n_pairs
#endif
				! Restore pp1 value for the next pp2 iteration
				pp1 = hoc_fine(i,j,k,thread)
                        enddo
                     enddo
                  enddo
                  
               enddo
            enddo
         enddo
      endif
      
      pp_ext_force_max(thread) = maxval(sqrt(pp_ext_force_accum(1,:,thread)**2 + pp_ext_force_accum(2,:,thread)**2 + pp_ext_force_accum(3,:,thread)**2))

#endif 
!  endif pp_ext_on_GPU

#ifdef DEBUG_PP_EXT
      write(*,*) 'pp_ext_force_max(',thread,') =', pp_ext_force_max(thread)!,'n_pairs =', n_pairs
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
