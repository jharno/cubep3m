 subroutine cubepm_fftw2(c, thread)
   implicit none
   include 'fftw_f77.i'
   include 'cubepm.fh'

   character c
   integer(4) :: thread

   if (firstfftw2(thread)) then
      call rfftw3d_f77_create_plan(fftw2_plan(thread),nf_tile, &
            nf_tile,nf_tile, FFTW_FORWARD, FFTW_MEASURE + FFTW_IN_PLACE)

      call rfftw3d_f77_create_plan(fftw2_iplan(thread),nf_tile, &
            nf_tile,nf_tile, FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE)
     firstfftw2(thread)=.false.
   endif

   if (c .eq. 'f') then
     call rfftwnd_f77_one_real_to_complex(fftw2_plan(thread),rho_f(:,:,:,thread),0)
   elseif (c .eq. 'b') then
     call rfftwnd_f77_one_complex_to_real(fftw2_iplan(thread),rho_f(:,:,:,thread),0) 
     rho_f(:,:,:,thread) = rho_f(:,:,:,thread) / real(nf_tile)**3
   elseif (c .eq. 'q') then
     call rfftwnd_f77_destroy_plan(fftw2_plan(thread)) 
     call rfftwnd_f77_destroy_plan(fftw2_iplan(thread))
   else 
     write(*,*) 'FFT direction error'
   endif

   return
 end subroutine cubepm_fftw2

