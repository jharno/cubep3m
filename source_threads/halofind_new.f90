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

! now gone!!

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

!!$!! Do halo analysis, dump to file
!!$
!!$    write(12) nhalo
!!$
!!$    do hi=1,nhalo     
!!$       radius_calc=(halo_mass(hi)/halo_odc/(4.0*pi/3.0))**(1.0/3.0)
!!$       search_limit(:,1)=int(1.0/real(mesh_scale)*halo_pos(:,hi)-1.0 &
!!$            /real(mesh_scale)*radius_calc-1.0)
!!$       search_limit(:,2)=int(1.0/real(mesh_scale)*halo_pos(:,hi)+1.0 &
!!$            /real(mesh_scale)*radius_calc+1.0)
!!$       imass=0
!!$       x_mean=0.0
!!$       v_mean=0.0
!!$       v2_mean=0.0
!!$       l=0.0
!!$      do k=search_limit(3,1),search_limit(3,2)
!!$        if (k < hoc_nc_l .or. k > hoc_nc_h) then
!!$          print *,'halo analysis out of bounds in z dimension',k
!!$          cycle
!!$        endif
!!$        do j=search_limit(2,1),search_limit(2,2)
!!$          if (j < hoc_nc_l .or. j > hoc_nc_h) then
!!$            print *,'halo analysis out of bounds in y dimension',j
!!$            cycle
!!$          endif
!!$          do i=search_limit(1,1),search_limit(1,2)
!!$            if (i < hoc_nc_l .or. i > hoc_nc_h) then
!!$              print *,'halo analysis out of bounds in x dimension',i
!!$              cycle
!!$            endif
!!$            pp=hoc(i,j,k)
!!$            mass_cell=0
!!$            do
!!$              if (pp==0) exit
!!$              dx=halo_pos(:,hi)-xv(:3,pp)
!!$              r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
!!$              if(r<radius_calc)then
!!$                 mass_cell=mass_cell+1
!!$                 imass=imass+1
!!$                 x_mean=x_mean+xv(:3,pp)
!!$                 v_mean=v_mean+xv(4:,pp)
!!$                 v2_mean=v2_mean+xv(4:,pp)*xv(4:,pp)
!!$                 l(1)=l(1)+(dx(3)*xv(5,pp)-dx(2)*xv(6,pp))
!!$                 l(2)=l(2)+(dx(1)*xv(6,pp)-dx(3)*xv(4,pp))
!!$                 l(3)=l(3)+(dx(2)*xv(4,pp)-dx(1)*xv(5,pp))
!!$              end if
!!$              pp=ll(pp)
!!$           enddo
!!$!            halo_particle_count ????
!!$          enddo
!!$        enddo
!!$      enddo
!!$
!!$         halo_pos(:,hi)=halo_pos(:,hi)+offset 
!!$         x_mean=x_mean/real(imass)+offset
!!$         v_mean=v_mean/real(imass)
!!$         v2_mean=v2_mean/real(imass)
!!$         l=l/real(imass)
!!$         v_disp=sqrt(v2_mean(1)+v2_mean(2)+v2_mean(3))
!!$         !!      write stuff to file or store in common arrays for Ifront 
!!$         !!      to store in common block arrays, each array should be
!!$         !!      max_maxima in size
!!$         if (halo_write .and. imass>0 .and. halo_mass(hi)>160) &
!!$              write(12) halo_pos(:,hi),x_mean,v_mean,l,v_disp,radius_calc,&
!!$              halo_mass(hi),imass*mass_p,halo_mass1(hi)
!!$         !! these plus nhalo should be passed to C^2Ray for each iteration
!!$         halo_x_mean(:,hi)=x_mean
!!$         halo_v_mean(:,hi)=v_mean
!!$         halo_l(:,hi)=l
!!$         halo_v_disp(hi)=v_disp
!!$         halo_radius_calc(hi)=radius_calc
!!$         halo_imass(hi)=imass*mass_p
!!$    enddo         
!!$
!!$!! close halo catalog
!!$
!!$    close(12)

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

    integer(4) :: i,j,k,pp,ix,iy,iz,ix0,iy0,iz0,k1,j1,i1,pp1
    integer(4) :: ic,iloc,ilcic,jlcic,klcic,thread

    real(4)    :: denmax,r,amass,amtot,fcf,fcf2

    real(4),dimension(3) :: x,fx,x0,dx

    real(4) :: para_inter 
    external para_inter 

!! calculate offsets for tile in local coords

    offset=tile*nf_physical_tile_dim-nf_buf
!    offset_fine=tile*nf_physical_tile_dim_halos-nf_buf_halos
 
!! initialize density 

    thread=1

!! limits for particle search.  
!! Ignore outmost buffer cells (4 fine cells). (is this still necessary?)

    cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1

!#ifdef SO_NEW    ! if defined use the new Spherical Overdensity halo 
                 ! finding algorithm, based on Tinker et al. 2008

!! calculate density associated with each particle in the
!! local volume
 
!    is_in_halo=.false. !initially no particles are yet in halos
    dens=1 !initially no particles are yet in halos, set density to 1
           !dens=0 will indicate that this particle is inside a halo 
           !dens>1 means it is a local peak

    ic=0 !will contain (local for this cuboid) number of peaks

    do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
        do i = cic_l(1), cic_h(1)
           !start with the first particle in the chain for this coarse grid cell
           pp=hoc(i,j,k) !head-of-chain for the linked list of particles inside
                         !the (i,j,k) coarse grid cell
           do
              if (pp==0) exit
              if(dens(pp)==1)then !particle does not yet belong to a halo 
                 flag=0
                 !go over all particles inside this and neighboring coarse grid cells
                 !
                 !Note: some work might be saved here by a more complex code (by looking
                 !for particles only in the cells that have overlap with search volume
                 do k1 = k-1, k+1
                    do j1 = j-1, j+1
                       do i1 = i-1, i+1
                          pp1=hoc(i1,j1,k1) !first particle in linked list
                          do
                             if (pp1==0) exit
                             !xv(1:3,pp1)=particle coordinates
                             dx=xv(1:3,pp1)-xv(1:3,pp)
                             r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2) !distance between the two particles
                             if(r<r_soft_peak)then
                                if(flag==0)then !peak is not yet counted
                                   if (ic > max_maxima -2) then
                                      write(*,*) 'too many halos'
                                      exit
                                   endif
                                   ic=ic+1
                                   flag=1 !makes sure we do not count this peak twice
                                end if
                                dens(pp)=dens(pp)+1
                                dens(pp1)=0 !mark particle, so it is not considered as a halo center
                             end if
                             pp1 = ll(pp1)!go to next particle
                          end do
                       end do
                    end do
                 end do
                 ipeak(:,ic)=(/i,j,k/)!which cell is the peak in?
                 den_peak(ic)=dens(pp) !what is the value of the peak?
                 peak_pointer(ic)=pp !this will tell us which particle marks the peak
                 !
                 pp = ll(pp)!go to next particle
              else
                 pp = ll(pp) !just go to next particle
              end if
           end do
        enddo
     enddo
  enddo

!how all particles with dens>1 are local peaks, locations are still rough

!! sort density maxima 

    isortpeak(:ic)=(/ (i,i=1,ic) /)
    call indexedsort(ic,den_peak(:),isortpeak(:))
    ipeak(:,:ic)=ipeak(:,isortpeak(:ic))
    peak_pointer(1:ic)=peak_pointer(isortpeak(1:ic))
    peak_pos(:,:ic)=peak_pos(:,isortpeak(:ic))

!#endif 

!! find mass in each halo

    nhalo=0 !halo number counter

    do iloc=ic,1,-1 !loop over peaks in decreasing density order
      ix0=ipeak(1,iloc)
      iy0=ipeak(2,iloc)
      iz0=ipeak(3,iloc)
      amtot=0

      ! WE ASSUME NO HALOS HAVE A RADIUS < nf_buf
      !(i.e. 24 cells) - usually a fine assumption,
      !but could be violated for largest halos in small boxes at late times

      do i=1,irtot !loop over cells at increasing distance from the one 
                   !containing the halo center and add particles into radial mass bins
        ix=ix0+idist(1,i)
        if (ix < 5 .or. ix > nf_tile_halos-4) cycle
        iy=iy0+idist(2,i)
        if (iy < 5 .or. iy > nf_tile_halos-4) cycle
        iz=iz0+idist(3,i)
        if (iz < 5 .or. iz > nf_tile_halos-4) cycle
        !we are in cell (ix,iy,iz) going outwards

        pp=hoc(ix,iy,iz)
        do   !loop over all particles in this cell
           if (pp==0) exit
           dx=xv(1:3,peak_pointer(iloc))-xv(1:3,pp)           
           r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2) !distance between the two particles

           mass_bins(floor(r*real(max_bins)/real(nf_buf))) = mass_bins(floor(r*real(max_bins)/real(nf_buf)))+1

           pp = ll(pp)
        end do
     enddo

     if (i>irtot) then
        i=irtot
        write(*,*) 'ran out of irtot'
     endif
     
     ! Check the overdensity criterion
     
     do ii=2,max_bins
        !here binned mass becomes cumulative and once divided by volume - overdensity
        mass_bins(ii)=(mass_bins(ii)+mass_bins(ii-1))
        if((8.*mass_bins(ii-1)/(4.*pi/3.*rad_bins(ii)**3)-halo_odc)* &
             (8.*mass_bins(ii)/(4.*pi/3.*rad_bins(ii)**3)-halo_odc)<0.)then
           !the overdensity threshold is reached. 
           !The 8 is to convert the mass from particles to grid masses.
           nhalo=nhalo+1
           halo_pos(:,nhalo)=peak_pos(:,iloc)+offset !(/os_x,os_y,os_z/)
           halo_mass(nhalo)=mass_bins(ii)
           halo_radius(nhalo)=rad_bins(ii)
           exit
        end if
     end do
     
     !refine the position of the halo center

     rad_current = halo_radius(nhalo)/3. !start from 1/3 of the halo radius

22   halo_x_mean(1:3,nhalo)=0.0
     num_particles=0
     diff=0.0d0

     !first converge on center-of-mass within R_Delta/3
     do while(diff > rsoft) !iterate to convergence. IS THIS CONDITION TOO STRINGENT?
        
        center_of_mass=0.0d0
        
        do i=1,irtot !loop over cells at increasing distance from the one 
           !containing the halo center and add particles into radial mass bins
           ix=ix0+idist(1,i)
           if (ix < 5 .or. ix > nf_tile_halos-4) cycle
           iy=iy0+idist(2,i)
           if (iy < 5 .or. iy > nf_tile_halos-4) cycle
           iz=iz0+idist(3,i)
           if (iz < 5 .or. iz > nf_tile_halos-4) cycle
           !we are in cell (ix,iy,iz) going outwards
           
           pp=hoc(ix,iy,iz)
           do   !loop over all particles in this cell
              if (pp==0) exit
              dx=xv(1:3,peak_pointer(iloc))-xv(1:3,pp)           
              r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2) !distance between the two particles
              
              if(r<rad_current)then
                 center_of_mass = center_of_mass + xv(:,peak_pointer(iloc))
                 num_particles=num_particles+1
              end do
              pp = ll(pp)
           end do
        enddo
        
        !center of mass position for particles inside R_Delta/3
        halo_x_mean(:,nhalo) = center_of_mass(:)/num_particles
        
        diff = sqrt((halo_x_mean(1,nhalo)-ix0)**2+ &
             (halo_x_mean(2,nhalo)-iy0)**2+ &
             (halo_x_mean(3,nhalo)-iz0)**2)
        
        ix0 = halo_x_mean(1,nhalo)
        iy0 = halo_x_mean(2,nhalo) 
        iz0 = halo_x_mean(3,nhalo) 
     end do

     !the center of mass is now converged 

     !Now refine halo center position

     if(num_particles>20 .and. rad_current > R_Delta/15.)then !reduce the radius and find center-of-mass again
        rad_current=0.99*rad_current
        go to 22
     end if

     ! Now, re-grow the sphere around the newly-found, converged center of mass

           do i=1,irtot !loop over cells at increasing distance from the one 
                   !containing the halo center and add particles into radial mass bins
        ix=ix0+idist(1,i)
        if (ix < 5 .or. ix > nf_tile_halos-4) cycle
        iy=iy0+idist(2,i)
        if (iy < 5 .or. iy > nf_tile_halos-4) cycle
        iz=iz0+idist(3,i)
        if (iz < 5 .or. iz > nf_tile_halos-4) cycle
        !we are in cell (ix,iy,iz) going outwards

        pp=hoc(ix,iy,iz)
        do   !loop over all particles in this cell
           if (pp==0) exit
           dx=xv(1:3,peak_pointer(iloc))-xv(1:3,pp)           
           r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2) !distance between the two particles

           mass_bins(floor(r*real(max_bins)/real(nf_buf))) = mass_bins(floor(r*real(max_bins)/real(nf_buf)))+1

           pp = ll(pp)
        end do
     enddo

     if (i>irtot) then
        i=irtot
        write(*,*) 'ran out of irtot'
     endif
     
     ! Check the overdensity criterion
     
     do ii=2,max_bins
        !here binned mass becomes cumulative and once divided by volume - overdensity
        mass_bins(ii)=(mass_bins(ii)+mass_bins(ii-1))
        if((8.*mass_bins(ii-1)/(4.*pi/3.*rad_bins(ii)**3)-halo_odc)* &
             (8.*mass_bins(ii)/(4.*pi/3.*rad_bins(ii)**3)-halo_odc)<0.)then
           !the overdensity threshold is reached. 
           !The 8 is to convert the mass from particles to grid masses.
           halo_mass(nhalo)=mass_bins(ii)
           halo_radius(nhalo)=rad_bins(ii)
           exit
        end if
     end do

! HERE WE SHOULD CHECK FOR A HALO CENTER BEING INSIDE ANOTHER HALO'S RADIUS




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

!! Do halo analysis, dump to file. Moved here from main halofind code. FIX

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


!halo_peak(nhalo)=den_peak(iloc)
!halo_radius(nhalo)=rdist(i)

! We want to store these variables rather than write them to disk
!
!      write(12,'(6f20.10)') halo_pos(1,iloc)+os_x,halo_pos(2,iloc)+os_y, &
!                                    halo_pos(3,iloc)+os_z,amtot,den_peak(iloc), &
!                                    rdist(i)
    enddo

!! clumping factor calculation  

! now gone!!

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

    do i=1,max_bins !set radial bins for mass binning around each density peak
       rad_bins(i)=real(i)*real(nf_buf)/real(max_bins) !WE ASSUME NO HALOS HAVE A RADIUS < nf_buf
    end do


    
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
