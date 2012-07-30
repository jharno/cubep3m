program PS_check
  !reads in halo catalogues and compares to the predictions of the
  !Press-Schechter and Sheth-Tormen approximations
  implicit none
  !
  !  integer(4), dimension(4) :: nhalo
  integer(4) :: nhalo, num_z
  real(4), dimension(3) :: halo_pos ! halo position (cells)
  real(4), dimension(3) :: x_mean ! centre of mass position 
  real(4), dimension(3) :: v_mean ! velocity
  real(4), dimension(3) :: l      ! angular momentum
  real(4) :: v_disp      ! velocity dispersion
  real(4) :: radius_calc ! radius
  real(4) :: halo_mass,halo_mass_uncorrected ! mass calculated on the grid (in grid masses)
  real(4) :: imass       ! mass calculated from the particles (in grid masses)
!!$  !real(4), dimension(MAX_ODC) :: odc=(/178.0,130.0,100.0,30.0/)
  !
  integer :: i,j,k,nn(3),ll,nbf, nbo, ii, lob
!  integer(kind=4) :: odc
  real(4) :: correction
  real(kind=4) :: dummy!,halocatalogue(8,max_halos),ratio
  !
  real*8 :: M, mass0, m_grid, mass_halos, m_min, m_max, m_box, m_grid0
  !  real*8,parameter :: omegabh2=0.02156,omega0=0.27,lambda0=0.73,h=0.7,an=1.0d0,tcmb=2.73
  !  real*8,parameter :: omegabh2=0.0223,omega0=0.24,lambda0=0.76,h=0.73,an=0.95d0,tcmb=2.73,s8=0.74 !wmap3
  real*8,parameter :: omegabh2=0.02156,omega0=0.27,lambda0=0.73,h=0.7,an=0.96d0,tcmb=2.726,s8=0.8 !wmap3+
  !  real*8,parameter :: omegabh2=0.02156,omega0=0.3,lambda0=0.7,h=0.7,an=1.0d0,tcmb=2.73,s8=0.9  !wmap1
  real*8 :: boxsize=425. !Mpc h^-1
  !  real*8 :: n_box=2048 !cells/side
  real*8 :: n_box=10976 !cells/side
  !
  real*8 :: omegab, h0, deltac, growth, mass8, sigma8, length_scale
  real*8 :: rho_crit_0, rho_crit_z, rho_crit_0_MpcM_sun
  real*8 :: rho_bar, deltac_z, z, nu, nu_old,a, m_ave
  real*8,parameter :: M_sun = 1.989d33, cm_to_Mpc=3.086d24
  real*8,parameter :: pi=3.1415927, grav_const=6.672d-8 !cm^3/g/s^2
  integer,parameter :: dim = 50, dim_ps=1000
  real*8 :: n_fine(dim_ps+1), mms_fine(0:dim_ps+1), slope(dim_ps+1),sigma(dim_ps)
  real*8 :: n(dim+1), n_simul(dim+1), mms(dim+1)
  !  real*8 :: f(dim+1),cl(dim+1), coll_frac, integr, M0,rho_bg_z
  integer,parameter :: max_size = 240
  !  real*8 :: kw(max_size), Deltak1(max_size), Deltak(max_size)
  real*8 :: mass(max_size), sigma1(max_size), sigma_cobe(max_size)
  real*8 :: table_flat(1000,2),table_open(1000,2)
  integer(kind=4),parameter :: MSL = 512
  character(MSL) :: ifile1,ifile2,ifile3,ifile4,nfile1,nfile2,nfile3,nfile4
  character(MSL) :: ifile5,ifile6,ifile7,ifile8,mfile1,mfile2,mfile3,mfile4
  !ST params
  real*8,parameter ::  alphapar=0.707, ppar=0.3,  Apar=0.322
  real*8 :: n_ST_fine(dim_ps+1), n_ST(dim+1) 
  !
  integer :: halo_type
  real(8) :: coll_frac(3) !collapsed fractions in halos: low-mass (1),high-mass(2),minihalos (3)
  real*8,parameter :: large_halo=1.e12
  real*8,parameter :: minihalo=1.e10
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
  M_box = rho_bar*(boxsize/h)**3
  M_grid = M_box/real(n_box)**3
  !  M_grid0 =1.341e5*(boxsize/35)**3!M_solar 
  print*,'M_grid=',m_grid
  !  stop
  !
  ! Read interpolation table
  !
  open(unit=1,file='/scratch/00506/ilievit/cubepm_100620_14_5488_425Mpc/code/utils/PSvsSim/deltac_flat.dat',status='old')
  nbf=0
  do i=1,1000
     read(1,*,end=6) dummy, table_flat(i,1), table_flat(i,2)
     nbf=nbf+1
  enddo
6 close(unit=1)
  !
  open(20,file='../utils/PSvsSim/ms.z0_wmap3plus') !contains mass scale vs. sigma
  do i=1,max_size
     read(20,*) mass(i),sigma1(i)
     mass(i)=mass(i)/h!since mass is in 1/h units in sigma file from cmbfast
     !print*, i,mass(i),sigma1(i)
  end do
  !  stop
  close(20)
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
!  pause
  !  stop
  !
  !  m_min=800.*m_grid
  m_min=160.*m_grid
  m_max=1e15
  !*** Loop over mass scales
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
  !  print*,'which redshift to do?'
  !read(*,*) z
!  open(1,file='redshifts0.dat')
  open(1,file='../utils/PSvsSim/redshifts.dat')
  open(2,file='../utils/PSvsSim/z_sourcenum_fcoll.dat')
  read(1,*) num_z 
  do j=1,num_z
     !
     coll_frac=0.0d0
     !
     read(1,*) z
     !
     !  print*,'deltac',deltac(omega0,lambda0,z)
     !print*,'growth',growth(omega0,lambda0,z),'z=',z
     deltac_z = deltac(omega0,lambda0,z)*growth(omega0,lambda0,z)
     do i=2,dim_ps+1
        slope(i) = -(log(sigma(i))-log(sigma(i-1)))/(log(mms_fine(i))-log(mms_fine(i-1)))
     end do
     write(50,*)z
     do ii=2,dim_ps+1
        nu = deltac_z/sigma(ii)
        m_ave=(mms_fine(ii)+mms_fine(ii-1))/2.
        n_fine(ii)=sqrt(2./pi)*rho_bar/m_ave**2*nu*exp(-nu**2/2.)*slope(ii)!m_sun^-1 Mpc^-3
        ! n_fine(ii)=sqrt(2./pi)*rho_bar/mms_fine(ii-1)**2*nu*exp(-nu**2/2.)*slope(ii)!m_sun^-1 Mpc^-3
        n_ST_fine(ii)=Apar*(1+(sqrt(alphapar)*nu)**(-2.*ppar))*sqrt(2./pi)*rho_bar/m_ave**2* &
             nu*sqrt(alphapar)*exp(-alphapar*nu**2/2.)*slope(ii)!m_sun^-1 Mpc^-3
        write(50,*) m_ave,n_fine(ii),n_ST_fine(ii)
     end do

     print*,'done PS'
     !close(50)
     !
     !do number of halos integrals over coarse mass bins
     !
     do i=1,dim
        n(i) = 0.0d0
        n_ST(i) = 0.0d0
    end do
     do ii=1,dim
        do i=2,dim_ps       
           if(mms_fine(i)>mms(ii-1).and.mms_fine(i)<=mms(ii))then
              n(ii) = n(ii)+(n_fine(i)+n_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
              n_ST(ii) = n_ST(ii)+(n_ST_fine(i)+n_ST_fine(i-1))*(mms_fine(i)-mms_fine(i-1))/2.
          end if
        end do
        n(ii)=n(ii)/(mms(ii)-mms(ii-1))
        n_ST(ii)=n_ST(ii)/(mms(ii)-mms(ii-1))
       ! write(38,*) (mms(ii)+mms(ii-1))/2.,n(ii)
     end do
     !stop
     !------------------------------------------------------------------------------------  
     !
     if(z.ge.10.)then
        write(ifile1,'("../../results/",f6.3,"halo.dat")') z
        write(nfile1,'("../../results/PS/",f6.3,"PS_ST_sim.dat")') z
        write(mfile1,'("../../results/PS/",f6.3,"mult_ST_sim.dat")') z
     else
        write(ifile1,'("../../results/",f5.3,"halo.dat")') z
        write(nfile1,'("../../results/PS/",f5.3,"PS_ST_resc.dat")') z
        write(mfile1,'("../../results/PS/",f5.3,"mult_ST_resc.dat")') z
     end if
     !
     print*,ifile1
     print*,nfile1,mfile1
     !  stop
     open(unit=31,file=ifile1,form='formatted')
     open(unit=101,file=nfile1)
     open(unit=102,file=mfile1)
     !
!       odc=1
!     do odc=1,4
        do i=1,dim+1
           n_simul(i) = 0.0d0
           nn=0
        end do
!     end do
     !mass correction factor: set all at 130 overdensity
     correction=1. 
!!$     correction(1)=0.893 !178
!!$     correction(2)=1.!130
!!$     correction(3)=1.103 !100
!!$     correction(4)=1.792 !30
     mass_halos = 0.0d0
     m_min=1.0d20
     m_max=0.0d0
!     do lob=0,6
!        read(21) nhalo
!        print*,'number of halos=', nhalo
     do !ll=1,nhalo
        !           read(21,end=133) halo_pos(:), x_mean(:,:), v_mean(:,:), l(:,:), &
        !                v_disp(:), radius_calc(:), halo_mass(:), imass(:)
!        read(21,'(6f20.10)',end=133) halo_pos(:), halo_mass, radius_calc, v_disp !
        read(unit=31,end=133,fmt=*) &
!        read(31,'(17f20.10)',end=133) 
             halo_pos(:), x_mean(:), v_mean(:), l(:), &
             v_disp, radius_calc, halo_mass, imass, halo_mass_uncorrected
!         print*,'check',halo_pos(:), x_mean(:), v_mean(:), l(:), &
!              v_disp, radius_calc, halo_mass, imass
!        pause

        if(halo_mass > 160. .and. imass >0)then
!           nn=nn+1
           mass0=halo_mass*M_grid/correction
           !        mass0=halo_mass_uncorrected*M_grid/correction
           !        mass0=imass*M_grid/correction
           !mass_halos = mass_halos+mass0
           !     print*,'mass check',l,mass_halos,mass
           if(mass0.gt.m_max) m_max=mass0
           if(mass0.lt.m_min.and. mass0 .gt.0.) m_min=mass0        
           do ii=1,dim+1
              if(mass0>mms(ii-1).and.mass0<=mms(ii)) n_simul(ii) = n_simul(ii)+1
           end do

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

           coll_frac(halo_type)=coll_frac(halo_type)+halo_mass*M_grid           
           nn(halo_type)=nn(halo_type)+1
        endif
     end do
     !100  print*,'N=',nn
133  continue

     coll_frac=coll_frac/M_box
     print*,'coll. fracs =',coll_frac
     
101  print*,'N=',nn
     print*,'min,max',m_min,m_max
     !     m_min=1e8
     !     m_max=1e13
     do ii=2,dim+1
        m_ave=(mms(ii)+mms(ii-1))/2.
        !        write(101,*) mms(ii-1),n(ii),(n_simul(ii,odc)/(boxsize/h)**3/(mms(ii)-mms(ii-1)),odc=1,4)
        !        write(102,*)  mms(ii-1),n(ii)*mms(ii-1)**2/rho_bar,(n_simul(ii,odc)/(boxsize/h)**3/(mms(ii)-mms(ii-1))*mms(ii-1)**2/rho_bar,odc=1,4)
        write(101,10) real(mms(ii-1)),real(n(ii)),real(n_ST(ii)),n_simul(ii)/(boxsize/h)**3/(mms(ii)-mms(ii-1))
        write(102,10) mms(ii-1),real(n(ii)*mms(ii-1)**2/rho_bar),real(n_ST(ii)*mms(ii-1)**2/rho_bar),&
             real(n_simul(ii)/(boxsize/h)**3/(mms(ii)-mms(ii-1))*mms(ii-1)**2/rho_bar)
     end do
     write(2,*) z,real(nn),real(coll_frac)
     close(unit=21)
     close(unit=101)
     close(unit=102)
  end do
10   format(4(e14.5))
     stop
     
end program PS_check
