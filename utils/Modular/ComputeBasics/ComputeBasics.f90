program ComputeDipole
  use Parameters
  use Variables
  use mMPI
  use FieldIO
  use Pencil
  use Power
  use Deriv3D
  use Dipole
  implicit none

  !--------------------!
  !Variable Declaration!
  !--------------------!

  character(len=*), parameter :: odir = './' !Location to write output files to

  character(len=200) :: filename
  real, dimension(Ncells,Ncells,Ncells) :: dm,nu,ha
  real, dimension(Ncells,Ncells,Ncells,3) :: vl,vl2
  real, parameter :: On = 3*0.1/93.14/0.67**2
  real, parameter :: Oc = 0.27+0.05
  real, parameter :: fn = On/(On+Oc)
  real, parameter :: fc = Oc/(On+Oc)

  integer, parameter :: nbins = 16
  real, dimension(nbins) :: pow,kpow,kbins,kernel
  real, dimension(nbins) :: pcc,pnn,phh,pcn,phc,phn,vp,vpg,vpv

  integer, parameter :: rbins = 30
  real, dimension(rbins) :: r,xcn,xhc,xhn,scn,shc,shn

  character(len=1), dimension(3), parameter :: xyz = (/'x','y','z'/)

  integer :: i,k,stat

  !---------!
  !Main Code!
  !---------!

  call start_mpi
  call setup_pencil
  call omp_set_num_threads(Nomp)

  !Read in density fields

  !!Cold dark matter
  filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'den'//trim(adjustl(rank_s))//'.bin'
  if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
  call read_field3(dm,trim(adjustl(filename)))
  dm = dm-1.0 !Density contrast

  !!Neutrinos
  filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'den'//trim(adjustl(rank_s))//'_nu.bin'
  if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
  call read_field3(nu,trim(adjustl(filename)))
  nu = nu-1.0

  !!Halos
  filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'den'//trim(adjustl(rank_s))//'_h.bin'
  if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
  call read_field3(ha,trim(adjustl(filename)))
  ha = ha-1.0

  !Compute power spectra and cross power spectra

  !Setup k bins
  do k=1,nbins
     kbins(k) = (k-1.0)*log10(nc/2-1.0)/(nbins-1.0)
  end do
  kbins = (2.0*pi/Lbox)*(10.0**kbins)
  
  !!DM Auto
  call one_write('Computing dark matter power spectrum')
  cube = dm
  call reset_pencil
  call cp_fftw(1)
  pow = 0.0; kpow = 0.0; kernel= 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab,kernel)
  pcc = pow
  slab2 = slab

  !!NU Auto
  call one_write('Computing neutrino power spectrum')
  cube = nu
  call reset_pencil
  call cp_fftw(1)
  pow = 0.0; kpow = 0.0; kernel= 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab,kernel)  
  pnn = pow
  slab3 = slab

  !!HA Auto
  call one_write('Computing halo power spectrum')
  cube = ha
  call reset_pencil
  call cp_fftw(1)
  pow = 0.0; kpow = 0.0; kernel= 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab,kernel)
  phh = pow

  !!DMxNU
  call one_write('Computing dark matter neutrino cross power')
  pow = 0.0; kpow = 0.0; kernel= 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab2,slab3,kernel)
  pcn = pow

  !!HAxDM
  call one_write('Computing halo dark matter cross power')
  pow = 0.0; kpow = 0.0; kernel= 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab2,kernel)
  phc = pow

  !!HAxNU
  call one_write('Computing halo neutrino cross power')
  pow = 0.0; kpow = 0.0; kernel= 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab3,kernel)
  phn = pow

  !Write out power spectra
  if (rank.eq.0) then
     filename = odir//redshift//'power.txt'
     write(*,*) 'Writing power spectra to file: '//trim(adjustl(filename))
     open(unit=11,file=filename,status='replace',iostat=stat,recl=500)
     if (stat/=0) call error_routine('Could not open file: '//trim(adjustl(filename)))
     do k=1,nbins
        write(11,*) kpow(k),pcc(k),pnn(k),phh(k),pcn(k),phc(k),phn(k)
     end do
     close(11)
  end if

  !Compute gradient dipoles
  call one_write('Computing dipoles')
  call gradient(fc*dm+fn*nu,vl)
  vl2 = vl
  call one_write('dm-nu')
  call compute_dipole(dm,nu,vl,r,xcn,scn)
  call one_write('ha-dm')
  call compute_dipole(ha,dm,vl,r,xhc,shc)
  call one_write('ha-nu')
  call compute_dipole(ha,nu,vl,r,xhn,shn)

  !Write out dipoles to file
  if (rank.eq.0) then
     filename = odir//redshift//'xi_g.txt'
     write(*,*) 'Writing dipoles to file: '//trim(adjustl(filename))
     open(unit=11,file=filename,status='replace',iostat=stat,recl=500)
     if (stat/=0) call error_routine('Could not open file: '//trim(adjustl(filename)))
     do k=1,rbins
        write(11,*) r(k),xcn(k),xhc(k),xhn(k),scn(k),shc(k),shn(k)
     end do
     close(11)
  end if

  !Read in dm velocity
  vl = 0.0
  do i=1,3
     filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'vel'//xyz(i)//trim(adjustl(rank_s))//'sim.bin'
     if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
     call read_field3(vl(:,:,:,i),trim(adjustl(filename)))
  end do

  call compute_velpower(vl,vl,vpv)
  call compute_velpower(vl2,vl2,vpg)
  call compute_velpower(vl,vl2,vp)
  if (rank.eq.0) then
     write(*,*) 'DM correlation with Gr'
     do k=1,nbins
        write(*,*) kpow(k),vp(k),vpv(k),vpg(k),vp(k)/(vpv(k)*vpg(k))**0.5
     end do
  end if

  !Read in nu velocity
  vl = 0.0
  do i=1,3
     filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'vel'//xyz(i)//trim(adjustl(rank_s))//'_nusim.bin'
     if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
     call read_field3(vl(:,:,:,i),trim(adjustl(filename)))
  end do

  call compute_velpower(vl,vl,vpv)
  call compute_velpower(vl2,vl2,vpg)
  call compute_velpower(vl,vl2,vp)
  if (rank.eq.0) then
     write(*,*) 'NU correlation with Gr'
     do k=1,nbins
        write(*,*) kpow(k),vp(k),vpv(k),vpg(k),vp(k)/(vpv(k)*vpg(k))**0.5
     end do
  end if

  xcn = 0.0; xhn = 0.0; xhc = 0.0;
  scn = 0.0; shn = 0.0; shc = 0.0;
  !Read in dm-nu velocity
  vl = 0.0
  do i=1,3
     filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'vel'//xyz(i)//trim(adjustl(rank_s))//'_nu-dmsim.bin'
     if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
     call read_field3(vl(:,:,:,i),trim(adjustl(filename)))
  end do

  call compute_velpower(vl,vl,vpv)
  call compute_velpower(vl2,vl2,vpg)
  call compute_velpower(vl,vl2,vp)
  if (rank.eq.0) then
     write(*,*) 'DM-NU correlation with Gr'
     do k=1,nbins
        write(*,*) kpow(k),vp(k),vpv(k),vpg(k),vp(k)/(vpv(k)*vpg(k))**0.5
     end do
  end if

  call compute_dipole(dm,nu,vl,r,xcn,scn)  

!!$  !Read in dm-ha velocity
!!$  vl = 0.0
!!$  do i=1,3
!!$     filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'vel'//xyz(i)//trim(adjustl(rank_s))//'_dm-hasim.bin'
!!$     if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
!!$     call read_field3(vl(:,:,:,i),trim(adjustl(filename)))
!!$  end do
!!$  vl = -1.0*vl !ha-dm
  call compute_dipole(ha,dm,vl,r,xhc,shc)

!!$  !Read in nu-ha velocity
!!$  vl = 0.0
!!$  do i=1,3
!!$     filename = dir//'node'//trim(adjustl(rank_s))//'/'//redshift//'vel'//xyz(i)//trim(adjustl(rank_s))//'_nu-hasim.bin'
!!$     if (rank.eq.0) write(*,*) 'Reading file: '//trim(adjustl(filename))
!!$     call read_field3(vl(:,:,:,i),trim(adjustl(filename)))     
!!$  end do
!!$  vl = -1.0*vl !ha-nu
  call compute_dipole(ha,nu,vl,r,xhn,shn)

  !Write out dipoles to file
  if (rank.eq.0) then
     filename = odir//redshift//'xi_v.txt'
     write(*,*) 'Writing dipoles to file: '//trim(adjustl(filename))
     open(unit=11,file=filename,status='replace',iostat=stat,recl=500)
     if (stat/=0) call error_routine('Could not open file: '//trim(adjustl(filename)))
     do k=1,rbins
        write(11,*) r(k),xcn(k),xhc(k),xhn(k),scn(k),shc(k),shn(k)
     end do
     close(11)
  end if

  call end_mpi

contains

  subroutine compute_velpower(g1,g2,p)
    implicit none
    real, dimension(Ncells,Ncells,Ncells,3), intent(in) :: g1,g2
    real, dimension(nbins), intent(out) :: p
    integer :: i

    p = 0.0
    do i=1,3
       call reset_pencil
       cube = g1(:,:,:,i)
       call cp_fftw(1)
       slab2 = slab
       call reset_pencil
       cube = g2(:,:,:,i)
       call cp_fftw(1)
       pow = 0.0; kpow = 0.0; kernel= 0.0
       call crosspowerspectrum(pow,kpow,kbins,slab,slab2,kernel)    
       p = p+pow
    end do

  end subroutine compute_velpower

  subroutine error_routine(str)
    implicit none
    character(len=*), intent(in) :: str
    write(*,*) 'Received error from rank: '//trim(adjustl(rank_s))
    write(*,*) 'Message: '//str
    call mpi_abort(mpi_comm_world)
  end subroutine error_routine

end program ComputeDipole
