!! Halofinding subroutine
 
  subroutine halofind 
    implicit none

    include 'mpif.h'
    include 'cubep3m.fh'

    real(4) :: z_write
    integer(4) :: i,j,k,fstat
    integer(4), dimension(3) :: tile

    character (len=max_path) :: ofile
    character (len=7) :: z_s
    character (len=3) :: t_s
    character (len=5) :: r_s

    integer(4) :: hi,pp,mass_cell
    integer(4), dimension(3,2) :: search_limit
    integer(8) :: imass
    real(4) :: r,radius_calc,v_disp,clump_factor
    real(4), dimension(3) :: x_mean,v_mean,v2_mean,l,dx,offset
    real(8) :: cfmassl,cfmassl2,cftmassl,cftmassl2

    offset(1)=cart_coords(3)*nf_physical_node_dim
    offset(2)=cart_coords(2)*nf_physical_node_dim
    offset(3)=cart_coords(1)*nf_physical_node_dim

    nhalo=0
    cfmass=0.0
    cfmass2=0.0
    cftmass=0.0
    cftmass2=0.0

!! find halos in the local volume

    do i=1,tiles_node
      tile(3) = (i-1) / (tiles_node_dim * tiles_node_dim)
      j = i - tile(3) * tiles_node_dim * tiles_node_dim
      tile(2) = (j-1) /  tiles_node_dim
      j = j - tile(2) * tiles_node_dim
      tile(1) = j - 1
      call find_halos(tile)
    enddo

!! open halo file

    z_write=z_halofind(cur_halofind)
    call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
    write(z_s,'(f7.3)') z_write 
    z_s=adjustl(z_s)
    
    write(r_s,'(i5)') rank
    r_s=adjustl(r_s)

    ofile=output_path//z_s(1:len_trim(z_s))//'halo'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=12,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=12,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening halo catalog for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif

!! calculate clumping factor on fine mesh 

    cfmassl=0.0 
    cfmassl2=0.0 
    call mpi_reduce(cfmass,cfmassl,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(cfmass2,cfmassl2,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    cftmassl=0.0 
    cftmassl2=0.0 
    call mpi_reduce(cftmass,cftmassl,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(cftmass2,cftmassl2,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    if (rank==0) then
      cfmass=(cfmassl2*nf_physical_dim**3)/cfmassl**2
      cftmass=(cftmassl2*nf_physical_dim**3)/cftmassl**2
      ofile=output_path//"fine_structure.dat"
      open(124,file=ofile,access='append',form='formatted')
      write(124,*) cftmass,cfmass,z_write
      close(124)
    endif

!! Write out fine clumping

    ofile=output_path//z_s(1:len_trim(z_s))//'fc'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening fine clumping file for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif
    write(72) fine_clumping 
    close(72)

!! construct coarse density and velocity field

   rho_c=0.0 !Have to fill it with zeros first?
   velocity_field=0.0

    do k = 0, nc_node_dim + 1
      do j = 0, nc_node_dim + 1
        do i = 0, nc_node_dim + 1
          pp=hoc(i,j,k)
          if (i <= 1 .or. i >= nc_node_dim .or. &
              j <= 1 .or. j >= nc_node_dim .or. &
              k <= 1 .or. k >= nc_node_dim) then
            call coarse_cic_mass_vel_boundry(pp)
          else
            call coarse_cic_mass_vel(pp)
          endif
        enddo
      enddo
    enddo
    
!! write coarse density
    ofile=output_path//z_s(1:len_trim(z_s))//'rho_c'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening coarse density file for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif
    write(72) rho_c
    close(72)

#ifdef COARSEST
!! Write out density field on the coarsened grid (for RT)

    ! redistribute data to coarsest grid
    do i=1,n_coarsest
       do j=1,n_coarsest
          do k=1,n_coarsest
             density_field_coarsest(i,j,k)=&
                  sum(rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                  coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                  coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))
          end do
       end do
    end do

!! write coarsest density
    ofile=output_path//z_s(1:len_trim(z_s))//'rho_coarsest'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening coarse density file for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif
    write(72) density_field_coarsest
    close(72)

!! Write out fine clumping on the coarsened grid (for RT)

    ! redistribute data to coarsest grid
    do i=1,n_coarsest
       do j=1,n_coarsest
          do k=1,n_coarsest
             fine_clumping_coarsest(i,j,k)=&
                  sum(fine_clumping(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                  coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                  coarsest_scale*k-coarsest_scale+1:coarsest_scale*k)&
                  *rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                  coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                  coarsest_scale*k-coarsest_scale+1:coarsest_scale*k)**2)/coarsest_scale**3
             fine_clumping_coarsest(i,j,k)=fine_clumping_coarsest(i,j,k)/&
                  (density_field_coarsest(i,j,k)/coarsest_scale**3)**2

!clumping3(i,j,k)=sum(1.0d0*clumping2(2*i-1:2*i,2*j-1:2*j,2*k-1:2*k)*rho2(2*i-1:2*i,2*j-1:2*j,2*k-1:2*k)**2)/8.
!              clumping3(i,j,k)=clumping3(i,j,k)/(rho3(i,j,k)/8.)**2
          end do
       end do
    end do


    ofile=output_path//z_s(1:len_trim(z_s))//'fc_coarsest'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening coarsened fine clumping file for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif
    write(72) fine_clumping_coarsest
    close(72)
#endif


#ifdef HALO_VEL_FIELD

!!  write velocity field 
    ofile=output_path//z_s(1:len_trim(z_s))//'vel'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening velocity field file for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif
    write(72) velocity_field 
    close(72)

#ifdef COARSEST
    ! mean bulk velocity: we need to divide the values in the file (which are sum over all 
    ! particles) by the number of particles in the cell 
    ! NOTE: assumes (as currently set) that coarsened_velocity_scale = mesh_scale (i.e. that 
    ! velocity_field (each direction) and rho_c have the same dimensions) 

    do i=1,nc_node_dim
       do j=1,nc_node_dim
          do k=1,nc_node_dim
             !note: the 8 is to convert grid masses to particle masses
             velocity_field(1,i,j,k) = velocity_field(1,i,j,k)/(rho_c(i,j,k)/8.)
             velocity_field(2,i,j,k) = velocity_field(2,i,j,k)/(rho_c(i,j,k)/8.)
             velocity_field(3,i,j,k) = velocity_field(3,i,j,k)/(rho_c(i,j,k)/8.)
          end do
       end do
    end do
    

    !! Write out velocity field on the coarsened grid (for RT)
    
    ! redistribute data to coarsest grid
    do i=1,n_coarsest
       do j=1,n_coarsest
          do k=1,n_coarsest
             if(sum(rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                  coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                  coarsest_scale*k-coarsest_scale+1:coarsest_scale*k)).gt.0.0d0)then
                velocity_field_coarsest(1,i,j,k)= &
                     sum(velocity_field(1,coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k)&
                     *rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))/&
                     sum(rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))!mean bulk velocity per coarsest cell
                velocity_field_coarsest(2,i,j,k)= &
                     sum(velocity_field(2,coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k)&
                     *rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))/&
                     sum(rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))!mean bulk velocity per coarsest cell
                velocity_field_coarsest(3,i,j,k)= &
                     sum(velocity_field(3,coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k)&
                     *rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))/&
                     sum(rho_c(coarsest_scale*i-coarsest_scale+1:coarsest_scale*i,&
                     coarsest_scale*j-coarsest_scale+1:coarsest_scale*j,&
                     coarsest_scale*k-coarsest_scale+1:coarsest_scale*k))!mean bulk velocity per coarsest cell
              else
                 velocity_field_coarsest(1,i,j,k)=0.0d0
                 velocity_field_coarsest(2,i,j,k)=0.0d0
                 velocity_field_coarsest(3,i,j,k)=0.0d0
              end if
           end do
        end do
     end do

!!  write coarsest velocity field
    ofile=output_path//z_s(1:len_trim(z_s))//'vel_coarsest'//  &
          r_s(1:len_trim(r_s))//'.dat'
#ifdef BINARY
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='binary')
#else
    open (unit=72,file=ofile,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening velocity field file for write'
      write(*,*) 'rank',rank,'file:',ofile
      stop
    endif
    write(72) velocity_field_coarsest
    close(72)

#endif
#endif

!! Do halo analysis, dump to file

    write(12) nhalo

    do hi=1,nhalo     
       radius_calc=(halo_mass(hi)/halo_odc/(4.0*pi/3.0))**(1.0/3.0)
       search_limit(:,1)=int(1.0/real(mesh_scale)*halo_pos(:,hi)-1.0 &
            /real(mesh_scale)*radius_calc-1.0)
       search_limit(:,2)=int(1.0/real(mesh_scale)*halo_pos(:,hi)+1.0 &
            /real(mesh_scale)*radius_calc+1.0)
       imass=0
       x_mean=0.0
       v_mean=0.0
       v2_mean=0.0
       l=0.0
      do k=search_limit(3,1),search_limit(3,2)
        if (k < hoc_nc_l .or. k > hoc_nc_h) then
          print *,'halo analysis out of bounds in z dimension',k
          cycle
        endif
        do j=search_limit(2,1),search_limit(2,2)
          if (j < hoc_nc_l .or. j > hoc_nc_h) then
            print *,'halo analysis out of bounds in y dimension',j
            cycle
          endif
          do i=search_limit(1,1),search_limit(1,2)
            if (i < hoc_nc_l .or. i > hoc_nc_h) then
              print *,'halo analysis out of bounds in x dimension',i
              cycle
            endif
            pp=hoc(i,j,k)
            mass_cell=0
            do
              if (pp==0) exit
              dx=halo_pos(:,hi)-xv(:3,pp)
              r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
              if(r<radius_calc)then
                 mass_cell=mass_cell+1
                 imass=imass+1
                 x_mean=x_mean+xv(:3,pp)
                 v_mean=v_mean+xv(4:,pp)
                 v2_mean=v2_mean+xv(4:,pp)*xv(4:,pp)
                 l(1)=l(1)+(dx(3)*xv(5,pp)-dx(2)*xv(6,pp))
                 l(2)=l(2)+(dx(1)*xv(6,pp)-dx(3)*xv(4,pp))
                 l(3)=l(3)+(dx(2)*xv(4,pp)-dx(1)*xv(5,pp))
              end if
              pp=ll(pp)
           enddo
!            halo_particle_count ????
          enddo
        enddo
      enddo

         halo_pos(:,hi)=halo_pos(:,hi)+offset 
         x_mean=x_mean/real(imass)+offset
         v_mean=v_mean/real(imass)
         v2_mean=v2_mean/real(imass)
         l=l/real(imass)
         v_disp=sqrt(v2_mean(1)+v2_mean(2)+v2_mean(3))
         !!      write stuff to file or store in common arrays for Ifront 
         !!      to store in common block arrays, each array should be
         !!      max_maxima in size
         if (halo_write .and. imass>0 .and. halo_mass(hi)>160) &
              write(12) halo_pos(:,hi),x_mean,v_mean,l,v_disp,radius_calc,&
              halo_mass(hi),imass*mass_p,halo_mass1(hi)
         !! these plus nhalo should be passed to C^2Ray for each iteration
         halo_x_mean(:,hi)=x_mean
         halo_v_mean(:,hi)=v_mean
         halo_l(:,hi)=l
         halo_v_disp(hi)=v_disp
         halo_radius_calc(hi)=radius_calc
         halo_imass(hi)=imass*mass_p
    enddo         

!! close halo catalog

    close(12)

    write(*,*) 'Finished halofind:',rank

!! Increment halofind counter 

    cur_halofind=cur_halofind+1

    halofind_step =.false.

  end subroutine halofind 

!-----------------------------------------------------------------------------!

  subroutine find_halos(tile)
    implicit none
    include 'cubep3m.fh'

!!  Halo finding variables

    integer(4), dimension(3) :: offset
    integer(4), dimension(3) :: cic_l,cic_h,tile

    integer(4) :: i,j,k,pp,ix,iy,iz,ix0,iy0,iz0,k1,j1,i1
    integer(4) :: ic,iloc,ilcic,jlcic,klcic,thread

    real(4)    :: denmax,r,amass,amtot,fcf,fcf2

    real(4),dimension(3) :: x,fx,x0

    real(4) :: para_inter 
    external para_inter 

!! calculate offsets for tile in local coords

    offset=tile*nf_physical_tile_dim-nf_buf
!    offset_fine=tile*nf_physical_tile_dim_halos-nf_buf_halos
 
!! initialize density 

    thread=1

    rho_f(:,:,:,thread)=0.0
    rho_f_halos(:,:,:)=0.0

!! limits for mass assignment.  Ignore outmost buffer
!! cells (4 fine cells).

    cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1

!! calculate fine mesh density for tile

    do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
        do i = cic_l(1), cic_h(1)
          pp=hoc(i,j,k)
#ifdef NGPH    !NGP/CIC used only for halo finding
          call fine_ngp_mass(pp,tile,thread)
#else
          call fine_cic_mass(pp,tile,thread)
#endif
        enddo
      enddo
    enddo

    do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
        do i = cic_l(1), cic_h(1)
          pp=hoc(i,j,k)
#ifdef NGPH    !NGP/CIC used only for halo finding
          call fine_ngp_mass_halos(pp,tile)
#else
          call fine_cic_mass_halos(pp,tile,thread)
#endif
        enddo
      enddo
    enddo



!    print*,'check fine densities',sum(rho_f),sum(rho_f_halos)
!    do i=1,20
!       j=10
!       k=10
!       if (rank == 0) &
!            print*, 'check values',i,j,k,rho_f(i,j,k,thread), &
!            sum(rho_f_halos(finer_halo_grid*(i-1)+1:finer_halo_grid*i,&
!            finer_halo_grid*(j-1)+1:finer_halo_grid*j,&
!            finer_halo_grid*(k-1)+1:finer_halo_grid*k))
!    end do

!! calculate clumping factor before halo extraction
!! Find density maxima

    ic=0
    do k=1+nf_buf,nf_buf+nf_physical_tile_dim
      do j=1+nf_buf,nf_buf+nf_physical_tile_dim
        do i=1+nf_buf,nf_buf+nf_physical_tile_dim
          cftmass=cftmass+rho_f(i,j,k,thread)
          cftmass2=cftmass2+rho_f(i,j,k,thread)**2
          denmax=maxval(rho_f(i-1:i+1,j-1:j+1,k-1:k+1,thread))
          if (denmax == rho_f(i,j,k,thread) .and. denmax > den_peak_cutoff) then
            if (ic > max_maxima -2) then
              write(*,*) 'too many halos'
              exit
            endif
            ic=ic+1

            ipeak(:,ic)=(/i,j,k/)
            den_peak(ic)=denmax

!parabolic interpolation for finding the peak position/halo center
            if (para_inter_hc) then
              x(1)=real(i-1)-0.5
              x(2)=real(i)-0.5
              x(3)=real(i+1)-0.5
              fx(1)=rho_f(i-1,j,k,thread)
              fx(2)=rho_f(i,j,k,thread)
              fx(3)=rho_f(i+1,j,k,thread)
              peak_pos(1,ic) = para_inter(x,fx)
              x(1)=real(j-1)-0.5
              x(2)=real(j)-0.5
              x(3)=real(j+1)-0.5
              fx(1)=rho_f(i,j-1,k,thread)
              fx(2)=rho_f(i,j,k,thread)
              fx(3)=rho_f(i,j+1,k,thread)
              peak_pos(2,ic) = para_inter(x,fx)
              x(1)=real(k-1)-0.5
              x(2)=real(k)-0.5
              x(3)=real(k+1)-0.5
              fx(1)=rho_f(i,j,k-1,thread)
              fx(2)=rho_f(i,j,k,thread)
              fx(3)=rho_f(i,j,k+1,thread)
              peak_pos(3,ic) = para_inter(x,fx)
            else
              peak_pos(1,ic) = real(i)-0.5
              peak_pos(2,ic) = real(j)-0.5
              peak_pos(3,ic) = real(k)-0.5
            endif

         endif
        enddo
      enddo
    enddo

!! sort density maxima 

    isortpeak(:ic)=(/ (i,i=1,ic) /)
    call indexedsort(ic,den_peak(:),isortpeak(:))
    ipeak(:,:ic)=ipeak(:,isortpeak(:ic))
    peak_pos(:,:ic)=peak_pos(:,isortpeak(:ic))

! on finer grid

    peak_pos_fine(:,:ic)=finer_halo_grid*peak_pos(:,:ic)
    ipeak_fine(:,:ic)=floor(peak_pos_fine(:,:ic))+1

!    do j=1,20
!       if(rank == 0) print*,'check peak positions',ic,peak_pos(:,j),&
!            peak_pos_fine(:,j),ipeak(:,j),ipeak_fine(:,j)
!    end do
    
!! find mass in each halo

    do iloc=ic,1,-1
      ix0=ipeak_fine(1,iloc)
      iy0=ipeak_fine(2,iloc)
      iz0=ipeak_fine(3,iloc)
      amtot=0
      do i=1,irtot
        ix=ix0+idist(1,i)
        if (ix < 5 .or. ix > nf_tile_halos-4) cycle
        iy=iy0+idist(2,i)
        if (iy < 5 .or. iy > nf_tile_halos-4) cycle
        iz=iz0+idist(3,i)
        if (iz < 5 .or. iz > nf_tile_halos-4) cycle
        amass=rho_f_halos(ix,iy,iz)
        rho_f_halos(ix,iy,iz)=0.0
!        amass=rho_f(ix,iy,iz,thread)
!        rho_f(ix,iy,iz,thread)=0.0
        amtot=amtot+amass
        if (complete_shell.and.rdist(i)==rdist(i+1)) cycle 
        if (i > 18 .and.(amtot/(real(i)))*finer_halo_grid**3 < halo_odc) then
           actual_odc=(amtot/(real(i)))*finer_halo_grid**3!remember what the actual overdensity inside the found halo is
           if (amtot >= min_halo_particles*mass_p) then
              nhalo=nhalo+1
              if (nhalo > max_maxima) then
                 print *,'too many halos to store',rank,nhalo,max_maxima
                 stop
              endif
           endif
           exit
        endif
     enddo

      if (i>irtot) then
        i=irtot
        write(*,*) 'ran out of irtot'
      endif

      if (nhalo > 0) then
        halo_pos(:,nhalo)=peak_pos(:,iloc)+offset !(/os_x,os_y,os_z/)
        halo_mass1(nhalo)=amtot
!        write(*,*) 'check mass rescaling 0',overdens(5),mass_rescaling(5)
!        write(*,*) 'check mass rescaling 1',amtot,halo_nondim_mass,actual_odc
        do iii=1,9999
           if((overdens(iii)-actual_odc)*(overdens(iii+1)-actual_odc)<0)then
              halo_mass(nhalo)=amtot*(halo_nondim_mass/mass_rescaling(iii+1))
!              write(*,*) 'check mass rescaling 2',iii,amtot,halo_mass(nhalo),mass_rescaling(iii+1)
!              write(*,*) 'check mass rescaling 3',actual_odc,overdens(iii),overdens(iii+1)
              exit
           end if
        end do
      endif

!halo_peak(nhalo)=den_peak(iloc)
!halo_radius(nhalo)=rdist(i)

! We want to store these variables rather than write them to disk
!
!      write(12,'(6f20.10)') halo_pos(1,iloc)+os_x,halo_pos(2,iloc)+os_y, &
!                                    halo_pos(3,iloc)+os_z,amtot,den_peak(iloc), &
!                                    rdist(i)
    enddo

!! clumping factor calculation  

    do k=1+nf_buf,nf_buf+nf_physical_tile_dim
      do j=1+nf_buf,nf_buf+nf_physical_tile_dim
        do i=1+nf_buf,nf_buf+nf_physical_tile_dim
          cfmass=cfmass+rho_f(i,j,k,thread)
          cfmass2=cfmass2+rho_f(i,j,k,thread)**2
        enddo
      enddo
    enddo

!! fine clumping factor

    do k=1+nf_buf,nf_buf+nf_physical_tile_dim,fine_clumping_scale
      iz=(k-nf_buf+tile(3)*nf_physical_tile_dim)/fine_clumping_scale+1
      do j=1+nf_buf,nf_buf+nf_physical_tile_dim,fine_clumping_scale
        iy=(j-nf_buf+tile(2)*nf_physical_tile_dim)/fine_clumping_scale+1
        do i=1+nf_buf,nf_buf+nf_physical_tile_dim,fine_clumping_scale
          ix=(i-nf_buf+tile(1)*nf_physical_tile_dim)/fine_clumping_scale+1
          fcf=0.0
          fcf2=0.0
          do k1=0,fine_clumping_scale-1
            do j1=0,fine_clumping_scale-1
              do i1=0,fine_clumping_scale-1
                fcf=fcf+rho_f(i+i1,j+j1,k+k1,thread)
                fcf2=fcf2+(rho_f(i+i1,j+j1,k+k1,thread))**2
              enddo
            enddo
          enddo
          fine_clumping(ix,iy,iz)=(fcf2*real(fine_clumping_scale)**3)/fcf**2
        enddo
      enddo
    enddo

  end subroutine find_halos 

!-----------------------------------------------------------------------------!
!! parabolic interpolation function

  real(4) function para_inter(x,fx)
  real(4), dimension(3) :: x,fx

  para_inter = x(2)-0.5*((x(2)-x(1))**2*(fx(2)-fx(3)) &
             -(x(2)-x(3))**2*(fx(2)-fx(1)))/((x(2)-x(1)) &
             *(fx(2)-fx(3))-(x(2)-x(3))*(fx(2)-fx(1)))

  end function para_inter

!-----------------------------------------------------------------------------!

!! Initialize halo finding arrays

  subroutine initialize_halofind

    implicit none
    include 'cubep3m.fh'
    include 'mpif.h'

    integer(4) :: ii, i, j, k, fstat
    real(4) :: r

! Loop through a box of length 2*nc_halo_max
! if cell is within sphere of radius = box length / 2
! include distince in rdist at entry ii
! ordered bottom left to top right

    ii=0
    do i=-nc_halo_max,nc_halo_max
      do j=-nc_halo_max,nc_halo_max
        do k=-nc_halo_max,nc_halo_max
          r=sqrt(real(i)**2+real(j)**2+real(k)**2)
          if (r>nc_halo_max) cycle
          ii=ii+1
          if (ii>nlist) then
            write(*,*) 'ii exceeded ',nlist
            pause
          endif
          idist(:,ii)=(/i,j,k/)
          rdist(ii)=r
        enddo
      enddo
    enddo
    irtot=ii

! sorts the rdist array from lowest to highest radial position
! from center of sphere saves rdist array position values in idist

    isortdist(:ii)=(/ (i,i=1,ii) /)
    call indexedsort(ii,rdist,isortdist)
    idist(:,:ii)=idist(:,isortdist(:ii))


    if (rank == 0) then
       
       ! read in table of dimensionless mass vs. overdensity for TIS+1/r^2 halo model
       ! for rescaling the halo mass in case of overshooting in overdensity
       open(11,file='table_M_Delta.dat', status='old', iostat=fstat)
       if (fstat /= 0) then
          write(*,*) 'error opening mass-overdensity table file'
          write(*,*) 'rank',rank
          call mpi_abort(mpi_comm_world,ierr,ierr)
       endif
       do iii=1,10000
          read(11,*) mass_rescaling(iii),overdens(iii)
!          write(*,*) 'table check',mass_rescaling(iii),overdens(iii)
       enddo
81     close(11)
    end if
    
    call mpi_bcast(mass_rescaling,10000,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(overdens,10000,mpi_real,0,mpi_comm_world,ierr)
    
  end subroutine initialize_halofind

!-----------------------------------------------------------------------------!

!! construct coarse mesh density and velocity field 
  subroutine coarse_cic_mass_vel(pp)
    implicit none

    include 'cubep3m.fh'

    integer(4), parameter :: CV=fine_clumping_scale/mesh_scale
    integer(4) :: pp
    integer(4), dimension(3) :: i1,i2
    real(4), dimension(3) :: x,dx1,dx2,vx1,vx2

    do
      if (pp == 0) exit
      x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5 !maybe this -0.5 is the problem??
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1.0 - dx1(:)

#ifdef HALO_VEL_FIELD
      vx1=xv(4:6,pp)*dx1(1)
      vx2=xv(4:6,pp)*dx2(1)
#endif

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) &
                               + dx1(1) * dx1(2) * dx1(3)
      rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) &
                               + dx2(1) * dx1(2) * dx1(3)
      rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) &
                               + dx1(1) * dx2(2) * dx1(3)
      rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) &
                               + dx2(1) * dx2(2) * dx1(3)
      rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) &
                               + dx1(1) * dx1(2) * dx2(3)
      rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) &
                               + dx2(1) * dx1(2) * dx2(3)
      rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) &
                               + dx1(1) * dx2(2) * dx2(3)
      rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) &
                               + dx2(1) * dx2(2) * dx2(3)

#ifdef HALO_VEL_FIELD
      velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
         velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx1*dx1(2)*dx1(3)
      velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
         velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx2*dx1(2)*dx1(3)
      velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
         velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx1*dx2(2)*dx1(3)
      velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
         velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx2*dx2(2)*dx1(3)
      velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
         velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx1*dx1(2)*dx2(3)
      velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
         velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx2*dx1(2)*dx2(3)
      velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
         velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx1*dx2(2)*dx2(3)
      velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
         velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx2*dx2(2)*dx2(3)
#endif

      pp = ll(pp)
    enddo

  end subroutine coarse_cic_mass_vel

!-----------------------------------------------------------------------------!

!! construct coarse mesh density and velocity field along nodal boundry
  subroutine coarse_cic_mass_vel_boundry(pp)
    implicit none

    include 'cubep3m.fh'

    integer(4), parameter :: CV=fine_clumping_scale/mesh_scale
    integer(4) :: pp
    integer(4), dimension(3) :: i1,i2
    real(4), dimension(3) :: x,dx1,dx2, vx1,vx2

    do
      if (pp == 0) exit
      x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5 !maybe this -0.5 is the problem??
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1.0 - dx1(:)

#ifdef HALO_VEL_FIELD
      vx1=xv(4:6,pp)*dx1(1)
      vx2=xv(4:6,pp)*dx2(1)
#endif

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      if (i1(3) >= 1 .and. i1(3) <= nc_node_dim) then
        if (i1(2) >= 1 .and. i1(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
             rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) + &
                  dx1(1) * dx1(2) * dx1(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
                  velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx1*dx1(2)*dx1(3)
#endif
          end if
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
             rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) + &
                  dx2(1) * dx1(2) * dx1(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
                  velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx2*dx1(2)*dx1(3)
#endif
          endif
       end if
       if (i2(2) >= 1 .and. i2(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
             rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) + &
                  dx1(1) * dx2(2) * dx1(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
                  velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx1*dx2(2)*dx1(3)
#endif
          end if
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
             rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) + &
                  dx2(1) * dx2(2) * dx1(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1) = &
                  velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i1(3)-1)/CV+1)+vx2*dx2(2)*dx1(3)
#endif
          end if
       endif
    endif

    if (i2(3) >= 1 .and. i2(3) <= nc_node_dim) then
       if (i1(2) >= 1 .and. i1(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
             rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) + &
                  dx1(1) * dx1(2) * dx2(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
                  velocity_field(:,(i1(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx1*dx1(2)*dx2(3)
#endif
          end if
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
             rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) + &
                  dx2(1) * dx1(2) * dx2(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
                  velocity_field(:,(i2(1)-1)/CV+1,(i1(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx2*dx1(2)*dx2(3)
#endif
          endif
       end if
       if (i2(2) >= 1 .and. i2(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
             rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) + &
                  dx1(1) * dx2(2) * dx2(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
                  velocity_field(:,(i1(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx1*dx2(2)*dx2(3)
#endif
          end if
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
             rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) + &
                  dx2(1) * dx2(2) * dx2(3)
#ifdef HALO_VEL_FIELD
             velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1) = &
                  velocity_field(:,(i2(1)-1)/CV+1,(i2(2)-1)/CV+1,(i2(3)-1)/CV+1)+vx2*dx2(2)*dx2(3)
#endif
          end if
       endif
    endif
    
    pp = ll(pp)
 enddo
  
end subroutine coarse_cic_mass_vel_boundry
