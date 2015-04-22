#include "./preprocessor"
program density_ps
  use Parameters
  use Variables
  use mMPI
  use FieldIO
  use Pencil
  use Power
  implicit none

  !Filename
  character(len=100) :: fnd

  !Binning info
  integer, parameter :: nbins = 64
  real, dimension(nbins) :: pow,kpow,kbins,kernel

  integer :: i

  call start_mpi

  if (rank==0) write(*,*) 'Starting program density_ps'

  call setup_pencil
  call omp_set_num_threads(Nomp)

  !Compute bins
  do i=1,nbins
    kbins(i) = (i-1.0)*log10(nc/2-1.0)/(nbins-1.0)
  end do
  kbins = (2.0*pi/Lbox)*(10.0**kbins)

  !$omp workshare
    cube = 0.0
  !$omp end workshare

  !DMxDM
  if (rank==0) write(*,*) 'Computing dark matter autospectrum'
  !!Obtain density contrast
  fnd = '/scratch2/p/pen/emberson/cubep3m/timing/cubep3m_movie_th2/node'//trim(adjustl(rank_s))//'/0.000den'//trim(adjustl(rank_s))//'.bin'
  call read_field3(cube,fnd)

  !!Compute power spectra
  call cp_fftw(1)
  kernel=0.0
  pow = 0.0
  kpow = 0.0
  call autopowerspectrum(pow,kpow,kbins,slab,kernel)

  !!Write spectrum
  fnd = dir//'ps/'//redshift//'delta_den_jd_dmdm.dat'
  if (rank==0) call write_spectrum(fnd,kbins,kpow,pow,kernel)

  !!Store dark matter for cross spectrum
  slab2 = slab

  !NuxNu
  if (rank==0) write(*,*) 'Computing neutrino autospectrum'
  call reset_pencil
  !$omp workshare
    cube=0.0
  !$omp end workshare
  !!Obtain density contrast
  fnd = '/scratch2/p/pen/emberson/cubep3m/timing/cubep3m_movie_th2/node'//trim(adjustl(rank_s))//'/0.000den'//trim(adjustl(rank_s))//'_nu.bin'
  call read_field3(cube,fnd)

  !!Compute power spectra
  call cp_fftw(1)
  kernel=0.0
  pow = 0.0
  kpow = 0.0
  call autopowerspectrum(pow,kpow,kbins,slab,kernel)

  !!Write spectrum
  fnd = dir//'ps/'//redshift//'delta_den_jd_nunu.dat'
  if (rank==0) call write_spectrum(fnd,kbins,kpow,pow,kernel)

  !DMxNU
  if (rank==0) write(*,*) 'Computing dark matter cross neutrino spectrum'
  kernel=0.0
  pow = 0.0
  kpow = 0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab2,kernel)
  !!Write spectrum
  fnd = dir//'ps/'//redshift//'delta_den_jd_nudm.dat'
  if (rank==0) call write_spectrum(fnd,kbins,kpow,pow,kernel)

  call cp_fftw(0)
  if (rank.eq.0) write(*,*) 'Finished program delta_ps'
  call end_mpi

contains

  subroutine write_spectrum(fileO,kbins,kpow,pow,kernel)
    implicit none
    character(len=*), intent(in) :: fileO
    real, dimension(:), intent(in) :: kbins, kpow, pow,kernel
    integer :: i,stat
    if (rank==0) write(*,*) 'Entering subroutine write_spectrum'
    open(unit=12,file=trim(adjustl(fileO)),status='replace',iostat=stat)
    if (stat/=0) then
      write(*,*) 'ERROR in program density_ps in subroutine write_spectrum opening file: '//trim(adjustl(fileO))
      call mpi_abort(mpi_comm_world,ierr,ierr)
    end if
    do i=1,size(pow)
      write(12,*) kbins(i), kpow(i), pow(i),kernel(i)
    end do
    close(12)
    if (rank==0) write(*,*) 'Finished subroutine write_spectrum'
  end subroutine write_spectrum

end program density_ps
