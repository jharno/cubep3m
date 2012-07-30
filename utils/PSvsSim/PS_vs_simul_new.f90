!! reads in halo catalogues and compares to the predictions of the
!! Press-Schechter approximation
!! Compile with: ifort PS_vs_simul_new.f90 deltac.f growth.f sigma_cobe_CMBfast.f spline.f90 splint.f90 -o PSvsSim
program PS_check
  implicit none

!frequently changed parameters are found in this header file:
  include '../../parameters'

!file path to halo catalogs 

  character(len=*), parameter :: halopath=output_path

!file path to redshift list of halo catalogs (from cubepm.par)

  character(len=*), parameter :: halofinds=cubepm_root//'input/halofinds_JD' 

!file path to table data
  
  character(len=*), parameter :: table_data=cubepm_root//'utils/PSvsSim/'

!size of smallest mass bin for determining mass function
! 1.0d8 - 1.0d11 in general.  Needs to be automated!

  real*8,parameter :: minimum_mass_bin=1.0d11

!these are all internal parameters and should not need modification

  integer,parameter :: max_halos=4000000
  integer,parameter :: max_size = 500
  integer,parameter :: dim = 70
  integer(kind=4),parameter :: MSL = 180

  real(8), parameter :: np = (nc/2.0)**3
  real(8), parameter :: sim_mp = (real(nc)**3)/np

  real*8,parameter :: omegabh2=0.0226,omega0=0.279,lambda0=0.721
  real*8,parameter :: h=0.701,an=0.96d0,tcmb=2.73
  real*8,parameter :: M_sun = 1.989d33, cm_to_Mpc=3.086d24
  real*8,parameter :: pi=3.1415927, grav_const=6.672d-8 !cm^3/g/s^2

  real(8), parameter :: Gkgs = grav_const/1.0d3
  real(8), parameter :: h0kgs = h*1.0d5/(cm_to_Mpc/100.0)
  real(8), parameter :: Vol = (box/h*(cm_to_Mpc/100.0))**3
  real(8), parameter :: rho_ckgs = (3.0*h0kgs**2)/(8.0*pi*Gkgs)
  real(8), parameter :: M_grid = (omega0*rho_ckgs*Vol/np)/sim_mp
  real(8), parameter :: M_grid_s = M_grid/(M_sun/1.0d3) !M_solar 

  integer :: i,j,k,nn,l,nbf, nbo, ii, halonum
  integer(4) :: num_halofinds,cur_halofind,fstat,partnum2(max_halos)
  real(4) :: dummy,halocatalogue(8,max_halos),ratio
  real*8 :: M,mass,mass_halos, m_min, m_max, phi, r, rho_max,xx,yy,zz
  real*8 :: omegab, h0, deltac, growth, sigma_cobe
  real*8 :: rho_crit_0, rho_crit_z, rho_crit_0_MpcM_sun
  real*8 :: rho_bar, deltac_z, z, nu, nu_old
  real*8 :: n(dim+1), mms(0:dim+1), length_scale(0:dim+1), slope(dim+1)
  real*8 :: n_simul(dim+1), sigma(0:dim+1), sigma_z(0:dim+1), f(dim+1)
  real*8 :: cl(dim+1), coll_frac, integr, M0,rho_bg_z
  real*8 :: kw(max_size), Deltak1(max_size), Deltak(max_size)
  real*8 :: table_flat(1000,2)
  real(4) :: z_halofind(max_size)
  character(7) :: z_s
  character(MSL) :: ifile2,nfile2
  real*4 reading_buffer(28)
#ifdef PLPLOT
  real*8 :: pkplot(3,dim)
#endif

  common /interp/ table_flat, nbo, nbf

  common/power/ kw,Deltak1,Deltak   
#ifdef PLPLOT
  common/pkplot_PS/ pkplot
#endif

!! Read in cmbfast power spectrum

  ifile2=table_data//'power_CMBfast_z0.dat'
  open(20,file=ifile2)
  do i=1,500
     read(20,*) kw(i),Deltak1(i)
  end do
  close(20)

  call spline(kw,Deltak1,max_size,1d31,1d31,Deltak)

!! Read interpolation table

  ifile2=table_data//'deltac_flat.dat'
  open(unit=1,file=ifile2,status='old')
  nbf=0
  do i=1,1000
     read(1,*,end=6) dummy, table_flat(i,1), table_flat(i,2)
     nbf=nbf+1
  enddo
6 close(unit=1)
  ifile2=''

!! Cosmological parameters.

  h0 = h*100.
  omegab = omegabh2/h**2
  rho_crit_0=1.88d-29*h**2   !g/cm^3
  rho_crit_0_MpcM_sun=rho_crit_0*cm_to_Mpc**3/M_sun

!! Read in halofinds to compare 

  open(11,file=halofinds,status='old',iostat=fstat)
  if (fstat /= 0) then
    write(*,*) 'error opening halofind list file'
    write(*,*) 'file:',halofinds
    stop
  endif
  do num_halofinds=1,max_size
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_halofind(num_halofinds)
  enddo
41  num_halofinds=num_halofinds-1
51  close(11)
  write(*,*) 'halofinds to recompose:'
  do i=1,num_halofinds
    write(*,'(f7.3)') z_halofind(i)
  enddo

do cur_halofind=1,num_halofinds
  z=z_halofind(cur_halofind)
  write(z_s,'(f7.3)') z
  z_s=adjustl(z_s)
  
  ifile2=halopath//z_s(1:len_trim(z_s))//"halo.dat"
  nfile2=halopath//z_s(1:len_trim(z_s))//"PS.dat"

  open(unit=11,file=ifile2) 
  open(unit=21,file=nfile2,status='replace')

!! Loop over mass scales

  do i=1,dim+1
     mms(i)=minimum_mass_bin*10.**(5.*(real(i-1))/dim)

!! Compute length scale in Mpc.

     length_scale(i)=1.533e-04*(3.*mms(i)/(4.*pi*omega0*h**2))**0.33333333

!! Compute sigma at z=0

     sigma(i)=sigma_cobe(omega0,omegabh2,lambda0,h0,tcmb,an,length_scale(i))
     sigma(i)=0.84/0.49*sigma(i) !normalize sigma_8 as in simulations
     write(0,*) mms(i),length_scale(i),sigma(i)
  end do

  rho_bar = rho_crit_0_MpcM_sun*omega0!*(1.+z)**3

  print *,'growth',growth(omega0,lambda0,z),'z=',z

  deltac_z = deltac(omega0,lambda0,z)*growth(omega0,lambda0,z)

  do i=2,dim+1
     slope(i) = -(log(sigma(i))-log(sigma(i-1)))/(log(mms(i))-log(mms(i-1)))
  end do

  do ii=2,dim+1
     nu = deltac_z/sigma(ii)
     n(ii)=sqrt(2./pi)*rho_bar/mms(ii)**2*nu*exp(-nu**2/2.)  !m_sun^-1 Mpc^-3
  end do

  print*,'done PS'

  do i=1,dim+1
     n_simul(i) = 0.0d0
  end do

  nn=0
  mass_halos = 0.0d0
  m_min=1.0d20
  m_max=0.0d0

  print *, 'Reading halo file'
  do l=1,max_halos
     read(11,'(28f20.4)',end=100)  reading_buffer! xx,yy,zz,M,rho_max,r
     M = reading_buffer(17)
     nn=nn+1
     mass=M*M_grid_s
     mass_halos = mass_halos+mass
     if(mass.gt.m_max) m_max=mass
     if(mass.lt.m_min) m_min=mass        
     do ii=1,dim+1
        if(mass>mms(ii-1).and.mass<=mms(ii)) n_simul(ii) = n_simul(ii)+1
     end do
  end do

100  print*,'N=',nn
     print*,'min,max',m_min,m_max

     close(11)

     j=dim
     do ii=2,dim+1
        write(21,*) mms(ii),n(ii)*slope(ii),n_simul(ii)/(box/h)**3/(mms(ii)-mms(ii-1))
#ifdef PLPLOT
        pkplot(1,ii-1)=real(mms(ii),kind=8)
        pkplot(2,ii-1)=real(n(ii)*slope(ii),kind=8)
        if (pkplot(2,ii-1)<1e-30) then
          j=ii-2
          exit
        endif
        pkplot(3,ii-1)=real(n_simul(ii)/(box/h)**3/(mms(ii)-mms(ii-1)),kind=8)
#endif
!        print*,'check', mms(ii),n(ii)*slope(ii),n_simul(ii)/(box/h)**3/(mms(ii)-mms(ii-1))
     end do
   
     close(21)

#ifdef PLPLOT
     print *,'minval',minval(pkplot(2,:))
     do ii=1,dim
       print *,'blah',pkplot(:,ii)
     enddo
     ii=3
     call plot_power(ii,j,pkplot(:,:dim),nfile2(1:len_trim(nfile2)))
#endif
 
enddo

end program PS_check
