!! reads in halo catalogues and compares to the predictions of the
!! Press-Schechter approximation
!! Compile with: ifort PS_vs_simul_new.f90 deltac.f growth.f sigma_cobe_CMBfast.f spline.f splint.f -o PSvsSim.x
program PS_check
  implicit none

  real(8), parameter :: box = 35.0 !Mpc h^-1
  integer(4), parameter :: nc = 320

  character(len=*), parameter :: file_path='/scratch/merz/cubepm_35Mpc/'
  character(len=*), parameter :: halofinds='/home/merz/codes/cubepm/input/halofinds'

  integer,parameter :: max_halos=4000000
  integer,parameter :: max_size = 500
  integer,parameter :: dim = 70
  integer(kind=4),parameter :: MSL = 180

  real(8), parameter :: np = (nc/2.0)**3
  real(8), parameter :: sim_mp = (real(nc)**3)/np

  real*8,parameter :: omegabh2=0.02156,omega0=0.27,lambda0=0.73,h=0.7,an=1.0d0,tcmb=2.73
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
  real*8 :: n(dim+1), mms(0:dim+1), length_scale(0:dim+1), slope(dim+1), n_simul(dim+1)
  real*8 :: sigma(0:dim+1), sigma_z(0:dim+1),f(dim+1),cl(dim+1), coll_frac, integr, M0,rho_bg_z
  real*8 :: kw(max_size), Deltak1(max_size), Deltak(max_size)
  real*8 :: table_flat(1000,2)
  real(4) :: z_halofind(max_size)
  character(7) :: z_s
  character(MSL) :: ifile2,nfile2

  common /interp/ table_flat, nbo, nbf

  common/power/ kw,Deltak1,Deltak   

!! Read in cmbfast power spectrum

  open(20,file='power_CMBfast_z0.dat')
  do i=1,500
     read(20,*) kw(i),Deltak1(i)
  end do
  close(20)

  call spline(kw,Deltak1,max_size,1d31,1d31,Deltak)

!! Read interpolation table

  open(unit=1,file='deltac_flat.dat',status='old')
  nbf=0
  do i=1,1000
     read(1,*,end=6) dummy, table_flat(i,1), table_flat(i,2)
     nbf=nbf+1
  enddo
6 close(unit=1)

!! Cosmological parameters.

  h0 = h*100.
  omegab = omegabh2/h**2
  rho_crit_0=1.88d-29*h**2!g/cm^3
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
    write(*,'(f5.1)') z_halofind(i)
  enddo

do cur_halofind=1,num_halofinds
  z=z_halofind(cur_halofind)
  write(z_s,'(f7.3)') z
  z_s=adjustl(z_s)
  
  ifile2=file_path//z_s(1:len_trim(z_s))//"halo.dat"
  nfile2=file_path//z_s(1:len_trim(z_s))//"PS.dat"

  open(unit=11,file=ifile2) 
  open(unit=21,file=nfile2,status='replace')

!! Loop over mass scales

  do i=1,dim+1
     mms(i)=1d10*10.**(5.*(real(i-1))/dim)

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

  do l=1,max_halos
     read(11,'(6f20.10)',end=100) xx,yy,zz,M,rho_max,r
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

     do ii=2,dim+1
        write(21,*) mms(ii),n(ii)*slope(ii),n_simul(ii)/(100./h)**3/(mms(ii)-mms(ii-1))
!        print*,'check', mms(ii),n(ii)*slope(ii),n_simul(ii)/(100./h)**3/(mms(ii)-mms(ii-1))
     end do
    
     close(21)
 
enddo

end program PS_check
