#include "./preprocessor"
program velocity_ps
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
  real, dimension(nbins) :: pow,kpow,kbins,powsum

  character(len=*), parameter :: jd = 'nu-dm'

  integer :: i

  call start_mpi

  if (rank==0) write(*,*) 'Starting program velocity_ps'

  call setup_pencil
  call omp_set_num_threads(Nomp)

  !Compute bins
  do i=1,nbins
    kbins(i) = (i-1.0)*log10(nc/2-1.0)/(nbins-1.0)
  end do
  kbins = (2.0*pi/Lbox)*(10.0**kbins)
  powsum = 0.0

  !$omp workshare
    cube = 0.0
  !$omp end workshare

  !NUxNU
  if (rank==0) write(*,*) 'Computing neutrino autospectrum'
  !!Obtain vel_x
  fnd = '/scratch2/p/pen/emberson/cubep3m/timing/cubep3m_movie_th2/node'//trim(adjustl(rank_s))//'/'//redshift//'velx'//trim(adjustl(rank_s))//'_'//jd//'rdm.bin'
  call read_field3(cube,fnd)
!  cube=cube-1.0
  !!Compute power spectrum
  call cp_fftw(1)
  pow = 0.0
  kpow = 0.0

  call crosspowerspectrum(pow,kpow,kbins,slab,slab)
  powsum = powsum+pow

  !!Obtain vel_y
  fnd = '/scratch2/p/pen/emberson/cubep3m/timing/cubep3m_movie_th2/node'//trim(adjustl(rank_s))//'/'//redshift//'vely'//trim(adjustl(rank_s))//'_'//jd//'rdm.bin'
  call read_field3(cube,fnd)

  !!Compute power spectrum
  call reset_pencil
  !cube=cube-1.0
  call cp_fftw(1)
  pow=0.0
  kpow=0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab)
  powsum = powsum+pow

  !!Obtain vel_z
  fnd = '/scratch2/p/pen/emberson/cubep3m/timing/cubep3m_movie_th2/node'//trim(adjustl(rank_s))//'/'//redshift//'velz'//trim(adjustl(rank_s))//'_'//jd//'rdm.bin'
  call read_field3(cube,fnd)

  !!Compute power spectrum
  call reset_pencil
  !cube=cube-1.0
  call cp_fftw(1)
  pow=0.0
  kpow=0.0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab)
  powsum = powsum+pow

  !!Write spectrum
  fnd = dir//'ps/'//redshift//'delta_vel_'//jd//'_jd_nunu.dat'
  if (rank==0) call write_spectrum(fnd,kbins,kpow,powsum)

  call cp_fftw(0)
  if (rank.eq.0) write(*,*) 'Finished program velocity_ps'
  call end_mpi

contains

  subroutine write_spectrum(fileO,kbins,kpow,pow)
    implicit none
    character(len=*), intent(in) :: fileO
    real, dimension(:), intent(in) :: kbins, kpow, pow
    integer :: i,stat
    if (rank==0) write(*,*) 'Entering subroutine write_spectrum'
    open(unit=12,file=trim(adjustl(fileO)),status='replace',iostat=stat)
    if (stat/=0) then
      write(*,*) 'ERROR in program velocity_ps in subroutine write_spectrum opening file: '//trim(adjustl(fileO))
      call mpi_abort(mpi_comm_world,ierr,ierr)
    end if
    do i=1,size(pow)
      write(12,*) kbins(i), kpow(i), pow(i)
    end do
    close(12)
    if (rank==0) write(*,*) 'Finished subroutine write_spectrum'
  end subroutine write_spectrum

end program velocity_ps
