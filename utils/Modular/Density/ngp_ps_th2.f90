#define VERBOSITY 1

program ngp_ps_th2
  use Parameters
  use Variables
  use mMPI
  use Zip
  use Pencil
  use Power
  implicit none

  !File names
  character(len=100) :: fn, fn2, fn3
  
  !Binning variables
  integer, parameter :: nbins = 64
  real, dimension(nbins) :: pow,kpow,kbins,kernel

  !Redshift variables
  character(len=*), parameter :: checkpoints_file = './checkpoints'
  real, dimension(100) :: checkpoints
  integer :: num_checkpoints
  character(len=100) :: checkpoint_str
  
  integer :: i,stat

  !Setup mpi, omp, pencil
  call start_mpi
#if (VERBOSITY>0)
  if (rank==0) write(*,*) 'Starting program ngp_ps_th2'
#endif
  call setup_pencil
  call omp_set_num_threads(Nomp)

  !Compute bins
  do i=1,nbins
    kbins(i) = (i-1.0)*log10(nc/2-1.0)/(nbins-1.0)
  end do
  kbins = (2.0*pi/Lbox)*(10.0**kbins)
  
  !Determine checkpoints
  checkpoints=-1
  open(unit=21,file=checkpoints_file,iostat=stat,status='old')
  if (stat/=0) call error("Could not open file: "//checkpoints_file)
  i=1
  do
     read(unit=21,end=41,fmt='(f20.10)') checkpoints(i)
     i=i+1
  end do
  41 num_checkpoints = i-1
  close(21)

#if (VERBOSITY>0)
  if (rank==0) write(*,*) "Number of checkpoints: ",num_checkpoints
#endif

  !Loop through checkpoints
  do i=1,num_checkpoints
     write(checkpoint_str,'(f7.3)') checkpoints(i)

#if (VERBOSITY>0)
     if (rank==0) write(*,*) "Computing power for checkpoint z="//trim(adjustl(checkpoint_str))
#endif
   
     call reset_pencil
     
     !Obtain density contrast for dark matter
     fn2 = dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(checkpoint_str))//'zip2_'//trim(adjustl(rank_s))//'.dat'
     fn3 = dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(checkpoint_str))//'zip3_'//trim(adjustl(rank_s))//'.dat'
     call ngp_delta_zip23(cube, fn2, fn3)

     !Fourier transform
     call cp_fftw(1)
     
     !Store slab in slab2
     !$omp workshare
     slab2 = slab
     !$omp end workshare

     !Obtain density contrast for neutrinos
     fn2 = dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(checkpoint_str))//'zip2_'//trim(adjustl(rank_s))//'_nu.dat'
     fn3 = dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(checkpoint_str))//'zip3_'//trim(adjustl(rank_s))//'_nu.dat'
     call ngp_delta_zip23(cube,fn2,fn3)

     !Fourier transform
     call reset_pencil
     call cp_fftw(1)

     !Compute DMxDM power spectrum
#if (VERBOSITY>0)
     if (rank==0) write(*,*) "Computing DMxDM power spectrum for z="//trim(adjustl(checkpoint_str))
#endif
     pow=0.0
     kpow=0.0
     kernel=0.0
     call autopowerspectrum(pow,kpow,kbins,slab2,kernel)
     
     fn = './ngp_density_delta'//trim(adjustl(checkpoint_str))//'_dmdm.dat'
     if (rank==0) call write_spectrum(fn,kbins,kpow,pow,kernel)

     !Compute DMxNU power spectrum
#if (VERBOSITY>0)
     if(rank==0) write(*,*) "Computing DMxNU power spectrum for z="//trim(adjustl(checkpoint_str))
#endif
     pow=0.0
     kpow=0.0
     kernel=0.0
     call crosspowerspectrum(pow,kpow,kbins,slab2,slab,kernel)

     fn = './ngp_density_delta'//trim(adjustl(checkpoint_str))//'_dmnu.dat'
     if (rank==0) call write_spectrum(fn,kbins,kpow,pow,kernel)

     !Compute NUxNU power spectrum
#if (VERBOSITY>0)
     if(rank==0) write(*,*) "Computing NUxNU power spectrum for z="//trim(adjustl(checkpoint_str))
#endif
     pow=0.0
     kpow=0.0
     kernel=0.0
     call autopowerspectrum(pow,kpow,kbins,slab,kernel)

     fn = './ngp_density_delta'//trim(adjustl(checkpoint_str))//'_nunu.dat'
     if (rank==0) call write_spectrum(fn,kbins,kpow,pow,kernel)

#if (VERBOSITY>0)
     if (rank==0) write(*,*) "Finished checkpoint z="//trim(adjustl(checkpoint_str))
#endif

  end do

  call cp_fftw(0)
#if (VERBOSITY>0)
  if (rank.eq.0) write(*,*) 'Finished program ngp_ps_th2'
#endif
  call end_mpi

contains

  subroutine write_spectrum(fileO,kbins,kpow,pow,kernel)
    implicit none
    character(len=*), intent(in) :: fileO
    real, dimension(:), intent(in) :: kbins, kpow, pow, kernel
    integer :: i,stat
#if (VERBOSITY>0)
    if (rank==0) write(*,*) 'Entering subroutine write_spectrum'
#endif
    open(unit=12,file=trim(adjustl(fileO)),status='replace',iostat=stat)
    if (stat/=0) call error("In subroutine write_spectrum, could not open file "//trim(adjustl(fileO)))
    do i=1,size(pow)
      write(12,*) kbins(i), kpow(i), pow(i), kernel(i)
    end do
    close(12)
#if (VERBOSITY>0)
    if (rank==0) write(*,*) 'Finished subroutine write_spectrum'
#endif
  end subroutine write_spectrum

  subroutine error(str)
    implicit none
    character(len=*), intent(in) :: str
    write(*,*) "Error in program ngp_ps_th2"
    write(*,*) "Message: "//str
    call mpi_abort(mpi_comm_world,ierr,ierr)
  end subroutine error

end program ngp_ps_th2
