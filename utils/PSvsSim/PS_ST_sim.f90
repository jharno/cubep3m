!ifort   PS_ST_sim.f90 deltac.f90 growth.f90 sigma_cobe_CMBfast.f90 spline.f90 splint.f90 -o PSvsSim
program PS_check
  !reads in halo catalogues and compares to the predictions of the
  !Press-Schechter, Sheth-Tormen and Tinker approximations
  implicit none
  include '../../parameters'
  !
  !  integer(4), dimension(4) :: nhalo
  integer(4) :: nhalo, num_z
  real(4), dimension(3) :: halo_pos ! halo position (cells)
  real(4), dimension(3) :: x_mean ! centre of mass position 
  real(4), dimension(3) :: v_mean ! velocity
  real(4), dimension(3) :: l      ! angular momentum
  real(4), dimension(3) :: v_disp      ! velocity dispersion
  real(4) :: radius_calc ! radius
  real(4) :: halo_mass,halo_mass_uncorrected ! mass calculated on the grid (in grid masses)
  real(4) :: imass       ! mass calculated from the particles (in grid masses)
  real(4), dimension(3) :: var      ! some kind of xyz variance of the halo shape
  real(4), dimension(6) :: I_ij     ! some kind of mass matrix of the halo shape
!!$  !real(4), dimension(MAX_ODC) :: odc=(/178.0,130.0,100.0,30.0/)
  !
  integer :: i,j,k,nn(3),ll,nbf, nbo, ii, lob
!  integer(kind=4) :: odc
  real(4) :: correction
  real(kind=4) :: dummy!,halocatalogue(8,max_halos),ratio
  !
  real*8 :: M, mass0, m_grid, mass_halos, m_min, m_max, m_box, m_grid0
  !  real*8,parameter :: omegabh2=0.02156,omega0=0.27,lambda0=0.73,h=0.7,an=1.0d0,tcmb=2.73
  real*8,parameter :: h=0.701,omegabh2=omega_b*h**2,omega0=omega_m,lambda0=omega_l,an=0.96d0,tcmb=2.73,s8=0.817 !wmap5
  !real*8,parameter :: omegabh2=0.0223,omega0=0.24,lambda0=0.76,h=0.73,an=0.95d0,tcmb=2.73,s8=0.74 !wmap3
  !real*8,parameter :: omegabh2=0.02156,omega0=0.27,lambda0=0.73,h=0.7,an=0.96d0,tcmb=2.73,s8=0.8 !wmap3+
  !  real*8,parameter :: omegabh2=0.02156,omega0=0.3,lambda0=0.7,h=0.7,an=1.0d0,tcmb=2.73,s8=0.9  !wmap1
  !real*8 :: boxsize=150 !Mpc/h 
  real*8,parameter :: boxsize=box !100 !64 !100 ! by T.T. Lu !11.4 !Mpc/h
  !  real*8 :: n_box=2048 !cells/side
  !real*8 :: n_box=6144 !cells/side
  real*8,parameter :: n_box=1024 !1024 !128 !3456 !128 ! by T.T. Lu !cells/side
  
  real*8 :: omegab, h0, deltac, growth, mass8, sigma8, length_scale
  real*8 :: rho_crit_0, rho_crit_z, rho_crit_0_MpcM_sun
  real*8 :: rho_bar, deltac_z, z, nu, nu_old,a, m_ave
  real*8,parameter :: M_sun = 1.989d33, cm_to_Mpc=3.086d24
  real*8,parameter :: pi=3.1415927, grav_const=6.672d-8 !cm^3/g/s^2
  integer,parameter :: dim = 50, dim_ps=1000
  real*8 :: n_fine(dim_ps+1), mms_fine(0:dim_ps+1), slope(dim_ps+1),sigma(dim_ps+1)
  real*8 :: n(dim+1), n_simul(dim+1), mms(dim+1)
  !  real*8 :: f(dim+1),cl(dim+1), coll_frac, integr, M0,rho_bg_z
  integer,parameter :: max_size = 240
  !  real*8 :: kw(max_size), Deltak1(max_size), Deltak(max_size)
  real*8 :: mass(max_size), sigma1(max_size), sigma_cobe(max_size)
  real*8 :: table_flat(1000,2),table_open(1000,2)
  integer(kind=4),parameter :: MSL = 512
  character(MSL) :: ifile1,ifile2,ifile3,ifile4,nfile1,nfile2,nfile3,nfile4,zstring
  character(MSL) :: ifile5,ifile6,ifile7,ifile8,mfile1,mfile2,mfile3,mfile4
  !ST parameters
  real*8,parameter ::  alphapar=0.707, ppar=0.3,  Apar=0.322
  real*8 :: n_ST_fine(dim_ps+1), n_ST(dim+1) 
  !Tinker2008 parameters, assuming \Delta=200 halos make it to the catalogs
  real*8,parameter ::  A0_T=0.186, a_0_T=1.47,  b0_T=2.57, c_T=1.19
  real*8,parameter ::  alpha_T = 0.4845! log_alpha = -[0.75/(log(200/75))]**1.2 = -0.7247
  real*8 :: Az_T, a_z_T, bz_T 
  real*8 :: n_Tinker_fine(dim_ps+1), n_Tinker(dim+1) 

  integer :: halo_type
  real(8) :: coll_frac(3) !collapsed fractions in halos: low-mass (1),high-mass(2),minihalos (3)
  real*8,parameter :: large_halo=1.e9
  real*8,parameter :: minihalo=1.e8
  real*4 reading_buffer(28)
  !
  common /interp/ table_flat, nbo, nbf
  !
  ! Cosmological parameters.
  !
  h0 = h*100.
  omegab = omegabh2/h**2
  rho_crit_0=1.88d-29*h**2!g/cm^3
  rho_crit_0_MpcM_sun=rho_crit_0*cm_to_Mpc**3/M_sun
  rho_bar = rho_crit_0_MpcM_sun*omega0
  print*,'rho_bar=',rho_bar
  !
  M_box = rho_bar*(boxsize/h)**3 ! In M_solar
  M_grid = M_box/real(n_box)**3
  !  M_grid0 =1.341e5*(boxsize/35)**3!M_solar 
  print*,'M_grid=',m_grid, 'M_box=', M_box
    !stop
  
  ! Read interpolation table
  open(unit=1,file='deltac_flat.dat',status='old')
  nbf=0
  do i=1,1000
     read(1,*,end=6) dummy, table_flat(i,1), table_flat(i,2)
     nbf=nbf+1
  enddo
6 close(unit=1)
  write(*,*) 'read deltac_flat.dat'


  open(20,file='ms.z0_wmap3plus') !contains mass scale vs. sigma
  do i=1,max_size
     read(20,*) mass(i),sigma1(i)
     mass(i)=mass(i)/h!since mass is in 1/h units in sigma file from cmbfast
     !print*, i,mass(i),sigma1(i)
  end do
  !  stop
  close(20)

  write(*,*) 'read ms.z0_wmap3plus'

  call spline(mass,sigma1,max_size,1d31,1d31,sigma_cobe)
  !
  ! Check the normalization of sigma
  !
  ! Normalization length scale in Mpc.
  !
  length_scale=8/h
  mass8 = 4.*pi*omega0*h**2/3*(length_scale/1.533e-04)**3
  call splint(mass,sigma1,sigma_cobe,max_size,mass8,sigma8)

  sigma1=sigma1*s8/sigma8

  print*,'check, sigma8=',sigma8,'mass8=',mass8,'length=',length_scale
  
  !  m_min=800.*m_grid
    m_min=160.*m_grid
 ! m_min=20.*m_grid !! by T.T. Lu
  m_max=1e15
    !print*, m_min,m_max
  !  pause

  write(*,*)' Loop over mass scales'
  do i=1,dim+1
     mms(i)=m_min*10.**(6.*(real(i-1))/dim)!mass bins for plotting
  end do
  do i=1,dim_ps+1
     !mass bins for PS integrals over above mass bins
     mms_fine(i)=m_min*10.**(6.*(real(i-1))/dim_ps)
     call splint(mass,sigma1,sigma_cobe,max_size,mms_fine(i),sigma(i))
  end do
  !stop
  !
  print*,'which redshift to do?'
  !read(*,*) z
  open(1,file='redshifts1.dat')
!  open(1,file='redshifts.dat')
  open(2,file='z_halonum_fcoll.dat')
  read(1,*) num_z
  
  do j=1,num_z
  !do j=num_z-1,num_z

     !write(*,*)'Testing on the last z slice only'
     !
     coll_frac=0.0d0
     !
     read(1,*) z
     !
     print*,'deltac',deltac(omega0,lambda0,z)
     print*,'growth',growth(omega0,lambda0,z),'z=',z
     deltac_z = deltac(omega0,lambda0,z)*growth(omega0,lambda0,z)
     
     !compute dlog(1/sigma) / dlogM
     do i=2,dim_ps+1
        slope(i) = -(log(sigma(i))-log(sigma(i-1)))/(log(mms_fine(i))-log(mms_fine(i-1)))
     end do

     if(z.ge.10.)then
        write(zstring,'(f6.3)') z
     else
        write(zstring,'(f5.3)') z
     endif

     !-----------
     !Write down fine grid predictions to file:
     ifile1 = output_path//trim(zstring)//"theory_fine.dat"
     nfile1 = output_path//trim(zstring)//"theory_fine_mult.dat"

     open(unit=50, file=ifile1)
     open(unit=51, file=nfile1)
     write(50,*)z     

     !----------------
     ! Compute Press-Schechter (PS), Sheth-Tormen (ST) and Tinker (T) on a fine grid.
     ! Here, nu is the collapse fraction at the given redshift, 
     ! n_fine is dn(M)/dM in the PS formalism, and n_ST_fine is the ST
     ! counterpart. We use the center of each bin in the calculation, so mms_fine(ii-1) -> m_ave

     do ii=2,dim_ps+1
        nu = deltac_z/sigma(ii)
        m_ave=(mms_fine(ii)+mms_fine(ii-1))/2.
        n_fine(ii)=sqrt(2./pi)*rho_bar/m_ave**2*nu*exp(-nu**2/2.)*slope(ii)!m_sun^-1 Mpc^-3
        ! n_fine(ii)=sqrt(2./pi)*rho_bar/mms_fine(ii-1)**2*nu*exp(-nu**2/2.)*slope(ii)!m_sun^-1 Mpc^-3
        n_ST_fine(ii)=Apar*(1+(sqrt(alphapar)*nu)**(-2.*ppar))*sqrt(2./pi)*rho_bar/m_ave**2* &
             nu*sqrt(alphapar)*exp(-alphapar*nu**2/2.)*slope(ii)!m_sun^-1 Mpc^-3

        Az_T = A0_T*(1+z)**(-0.14)
        a_z_T = a_0_T*(1+z)**(-0.06)
        bz_T = b0_T*(1+z)**(-alpha_T)
        n_Tinker_fine(ii) = Az_T*( (sigma(ii)/bz_T)**(-a_z_T) + 1)*exp(-c_T/sigma(ii)**2)*rho_bar/m_ave**2*slope(ii)

        write(50,*) real(m_ave,4),real(n_fine(ii),4) ,real(n_ST_fine(ii),4),real(n_Tinker_fine(ii),4)
        write(51,*) real(m_ave,4),real(n_fine(ii)*m_ave**2/rho_bar,4), real(n_ST_fine(ii)*m_ave**2/rho_bar,4),real(n_Tinker_fine(ii)*m_ave**2/rho_bar,4)
     end do
     close(50)
     close(51)
     print*,'done fine grid'
     !----------------------
     !stop
     !----------------------

     !do number of halos integrals over coarse mass bins     
     do i=1,dim
        n(i) = 0.0d0
        n_ST(i) = 0.0d0
        n_Tinker(i) = 0.0d0
     end do

     ! n is the PS dn(M_c)/dM_c on the coarse mass bin dM_c = mms(ii) - mms(ii-1)
     ! Similarly for ST and Tinker predictions. To get these numbers, we first calculate 
     ! dn(M_c) = mean (over fine bins i inside coarse bin ii) of dn(M).
     ! So we multiply n_fine(i) by dM (mms_fine(i) - mms_fine(i-1) to get dN(M(i)).
     ! We then get the value at the center of each bin (n_fine(i) + n(fine(i-1))/2, 
     ! and sum over all contributions in the coarse bin ii. This gives dN(M_c).
     ! Finally, we divide by dM_c to get dn(M_c)/dM_c.

     do ii=2,dim+1
        do i=2,dim_ps       
           if(mms_fine(i)>mms(ii-1).and.mms_fine(i)<=mms(ii))then
              n(ii) = n(ii)+(n_fine(i)+n_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
              n_ST(ii) = n_ST(ii)+(n_ST_fine(i)+n_ST_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
              n_Tinker(ii) = n_Tinker(ii)+(n_Tinker_fine(i)+n_Tinker_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
           end if
        end do
        n(ii)=n(ii)/(mms(ii)-mms(ii-1))
        n_ST(ii)=n_ST(ii)/(mms(ii)-mms(ii-1))
        n_Tinker(ii)=n_Tinker(ii)/(mms(ii)-mms(ii-1))
     end do
     !stop
     !------------------------------------------------------------------  
     
     if(z.ge.10.)then
        write(zstring,'(f6.3)') z 
     else
        write(zstring,'(f5.3)') z 
     endif

     ifile1 = output_path//trim(zstring)//"halo.dat"
     nfile1 = output_path//trim(zstring)//"PS_ST_sim.dat"
     mfile1 = output_path//trim(zstring)//"mult_ST_sim.dat"
 
     
     print*,trim(ifile1)
     print*,trim(nfile1),trim(mfile1)
     
     open(unit=31,file=ifile1,form='formatted')
     open(unit=101,file=nfile1, form='formatted')
     open(unit=102,file=mfile1, form='formatted')
     
     !odc=1
     !do odc=1,4
        do i=1,dim+1
           n_simul(i) = 0.0d0
           nn=0
        end do
     !end do

     ! Not sure what that is...
     !mass correction factor: set all at 130 overdensity
     correction=1. 
!!$     correction(1)=0.893 !178
!!$     correction(2)=1.!130
!!$     correction(3)=1.103 !100
!!$     correction(4)=1.792 !30



     !-------------------------
     !  Start Reading Halo File

     !  read(21) nhalo
     !  print*,'number of halos=', nhalo
     do !ll=1,nhalo
        ! read(21,end=133) halo_pos(:), x_mean(:,:), v_mean(:,:), l(:,:), &
        !                v_disp(:), radius_calc(:), halo_mass(:), imass(:)
        ! read(21,'(6f20.10)',end=133) halo_pos(:), halo_mass, radius_calc, v_disp !
        ! read(unit=31,end=133,fmt=*) &
       !do

        read(31,'(28f20.4)',end=133) reading_buffer
        !read(31, end=133) reading_buffer
        halo_pos(:) = reading_buffer(1:3)
        x_mean(:) = reading_buffer(4:6)
        v_mean(:) = reading_buffer(7:9)
        l(:)= reading_buffer(10:12)
        halo_mass = reading_buffer(17)
        imass = reading_buffer(18)
 
        !write(*,*) 'halo_mass = ', halo_mass ,'sim units'
        !write(*,*) 'halo_mass = ', halo_mass*M_grid/correction ,'M_sun'
        ! print*,'check',halo_pos(:), x_mean(:), v_mean(:), l(:), &
        !      v_disp, radius_calc, halo_mass, imass        
        !stop

        !if(halo_mass > 160. .and. imass >0)then
        if(halo_mass > 20. .and. imass >0)then ! by T.T. Lu
!           nn=nn+1
           mass0=halo_mass*M_grid/correction
           !        mass0=halo_mass_uncorrected*M_grid/correction
           !        mass0=imass*M_grid/correction
           !mass_halos = mass_halos+mass0
           !     print*,'mass check',l,mass_halos,mass
           if(mass0.gt.m_max) m_max=mass0
           if(mass0.lt.m_min.and. mass0 .gt.0.) m_min=mass0        
           do ii=2,dim+1
              if(mass0>mms(ii-1).and.mass0<=mms(ii)) n_simul(ii) = n_simul(ii)+1
           end do

           !----------------------
           ! Calculate the collapsed fractions
           if ( halo_mass*M_grid > large_halo ) then 
              halo_type=2
              ! these are the larger, never suppressed halos
           elseif ( halo_mass*M_grid < large_halo .and. halo_mass*M_grid > minihalo) then
              halo_type=1
                 ! these are the smaller, Pop. II/III, suppressed halos
                 ! include only well-resolved halos, with >20 particles
           else !this is a minihalo
              halo_type=3
           endif
           !write(*,*) 'halo_type = ', halo_type
           
           coll_frac(halo_type)=coll_frac(halo_type)+halo_mass*M_grid           
           nn(halo_type)=nn(halo_type)+1
        endif
     end do
     100  print*,'N=',nn
133  continue

     coll_frac=coll_frac/M_box
     print*,'coll. fracs =',coll_frac
     
101  print*,'N=',nn
     print*,'min,max',m_min,m_max
     !     m_min=1e8
     !     m_max=1e13

     !-------------
     ! Write to file
     do ii=2,dim+1
        m_ave=(mms(ii)+mms(ii-1))/2.
        ! write(101,*) mms(ii-1),n(ii),(n_simul(ii,odc)/(boxsize/h)**3/(mms(ii)-mms(ii-1)),odc=1,4)
        ! write(102,*) mms(ii-1),n(ii)*mms(ii-1)**2/rho_bar,(n_simul(ii,odc)/(boxsize/h)**3/(mms(ii)-mms(ii-1))*mms(ii-1)**2/rho_bar,odc=1,4)
        write(101,'(5E15.4)') real(mms(ii-1)),real(n(ii)),real(n_ST(ii)),n_simul(ii)/(boxsize/h)**3/(mms(ii)-mms(ii-1)),real(n_Tinker(ii))
        write(102,'(5E15.4)') mms(ii-1),real(n(ii)*mms(ii-1)**2/rho_bar),real(n_ST(ii)*mms(ii-1)**2/rho_bar),&
             real(n_simul(ii)/(boxsize/h)**3/(mms(ii)-mms(ii-1))*mms(ii-1)**2/rho_bar), real(n_Tinker(ii)*mms(ii-1)**2/rho_bar)
     end do
     !--------------

     write(2,*) z,real(nn),real(coll_frac)
     close(unit=21)
     close(unit=101)
     close(unit=102)
  end do
10   format(4(e14.5))
     close(unit=2) ! by T.T.Lu
     close(unit=1) ! by T.T.Lu
     stop
     
end program PS_check
