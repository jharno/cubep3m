!*******************************************************
! fft_fine.f90 - subroutines 
! for fft on the fine mesh. 
! Written by Joachim Harnois-Deraps, 2013/02
!*******************************************************

subroutine cubepm_fftw2(c, thread)
   use omp_lib
   implicit none
   !include 'fftw_f77.i'
   include 'cubepm.fh'
   include 'fftw3.f03'

   character c
   integer(4) :: thread

   if (firstfftw2(thread)) then

!#ifdef _OPENMP
!#ifdef OPENMP
!      call omp_set_lock(lockfftw2)
!      write (*,*) "lockfftw2 set!"
!#endif
      !call sfftw_init_threads()
      !call sfftw_plan_with_nthreads(cores)
      !write(*,*)'using OMP in FFT'
 
      call sfftw_plan_dft_r2c_3d(fftw2_plan(thread),nf_tile,nf_tile,nf_tile,rho_f(:,:,:,thread),rho_f(:,:,:,thread),FFTW_ESTIMATE)
      call sfftw_plan_dft_c2r_3d(fftw2_iplan(thread),nf_tile,nf_tile,nf_tile,rho_f(:,:,:,thread),rho_f(:,:,:,thread),FFTW_ESTIMATE)
 
      !call rfftw3d_f77_create_plan(fftw2_plan(thread),nf_tile, &
      !      nf_tile,nf_tile, FFTW_FORWARD, FFTW_MEASURE + FFTW_IN_PLACE)

      !call rfftw3d_f77_create_plan(fftw2_iplan(thread),nf_tile, &
      !      nf_tile,nf_tile, FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE)

     firstfftw2(thread)=.false.
!#ifdef _OPENMP
!#ifdef OPENMP
!      write (*,*) "releasing lockfftw2!"
!      call omp_unset_lock(lockfftw2)
!#endif
   endif

   if (c .eq. 'f') then
     !call rfftwnd_f77_one_real_to_complex(fftw2_plan(thread),rho_f(:,:,:,thread),0)
     call sfftw_execute(fftw2_plan(thread))
   elseif (c .eq. 'b') then
     !call rfftwnd_f77_one_complex_to_real(fftw2_iplan(thread),rho_f(:,:,:,thread),0) 
     call sfftw_execute(fftw2_iplan(thread))
     rho_f(:,:,:,thread) = rho_f(:,:,:,thread) / real(nf_tile)**3
   elseif (c .eq. 'q') then

!#ifdef _OPENMP
!#ifdef OPENMP
!     call omp_set_lock(lockfftw2)
!     WRITE (*,*) "lockfftw2 set (destroy)!"
!#endif
     call sfftw_destroy_plan(fftw2_plan(thread))
     call sfftw_destroy_plan(fftw2_iplan(thread))
 
     !call rfftwnd_f77_destroy_plan(fftw2_plan(thread)) 
     !call rfftwnd_f77_destroy_plan(fftw2_iplan(thread))
!#ifdef _OPENMP
!#ifdef OPENMP
!     write (*,*) "releasing lockfftw2 (destroy)!"
!     call omp_unset_lock(lockfftw2)
!#endif
   elseif(c .eq. 'o') then
       ! Initialization of FFT plans in a thread safe manner
   else 
     write(*,*) 'FFT direction error'
   endif

   return
 end subroutine cubepm_fftw2

