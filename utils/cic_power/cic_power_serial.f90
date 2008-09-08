!compile with ifort cic_power_serial.f90 -fpp -O3 -DBINARY -DFFTW -I$MCKENZIE_FFTW_INC_PATH $MCKENZIE_FFTW_LIB_PATH/libfftw3f.a -o cic_power_serial.x
implicit none

 !! Cosmo parameters
  real, parameter :: box=100.0

  !! dir is the directory for files
  character(*), parameter :: ipl='/scratch/merz/p3m/0.000xv0.dat' !'/uber-scratch/merz/cubepm_8_320_100Mpc_060424.out/xv.ic' !'/cita/scratch/eh14/merz/cubepm_no_fmesh_scaled/0.000xv0.dat' 
  character(*), parameter :: ops='/scratch/merz/p3m/p3m_80_ps.dat'

  !! nc is the number of cells per box length
  real, parameter :: nc_orig=80.
  integer, parameter :: nc=320

  integer, parameter :: np=160**3
  real(4), dimension(6,np) :: xv
 
  real(4), dimension(2,nc) :: psdm

  real(4), dimension(nc+2,nc,nc) :: d

  integer, parameter :: nt=1
  character(len=80) :: fn
  integer :: cp,k,np_in,i
  real :: CRUD(11) !for pp
  real :: kr

#ifdef FFTW
call fft3d(d,'f')
#endif

  d=-1.0
  psdm=0.0
#ifdef BINARY
  open(11,file=ipl,form='binary')
#else
  open(11,file=ipl,form='unformatted')
#endif
  read(11) np_in
  read(11) CRUD
!  if (np/=np) then
!    print *,'np_in .ne. np'
!    stop
!  endif
!  print *,np_in
!  print *,CRUD
  read(11) xv
  close(11)
  do i=1,np
    xv(:,i)=xv(:,i)*(real(nc)/nc_orig)
  enddo
  call cicmass
#ifdef FFTW
  call fft3d(d,'f')
#else
  call fft3d(d,nc,'f')
#endif
  call powerspectrum(d,psdm)

  open(11,file=ops,recl=500)
  do k=2,nc/2+1
    kr=2*3.141592654*(k-1)/box
    write(11,*) kr,psdm(:,k)
  enddo
  close(11)

contains
 subroutine cicmass
    implicit none
    real, parameter :: ncr=nc
    real, parameter :: mp=(ncr)**3/np 

    integer i
    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    do i=1,np
       x=mod(xv(1,i)-0.5+ncr,ncr)
       y=mod(xv(2,i)-0.5+ncr,ncr)
       z=mod(xv(3,i)-0.5+ncr,ncr)

       i1=floor(x)+1
       i2=mod(i1,nc)+1
       dx1=i1-x
       dx2=1-dx1
       j1=floor(y)+1
       j2=mod(j1,nc)+1
       dy1=j1-y
       dy2=1-dy1
       k1=floor(z)+1
       k2=mod(k1,nc)+1
       dz1=k1-z
       dz2=1-dz1

       dz1=mp*dz1
       dz2=mp*dz2
       d(i1,j1,k1)=d(i1,j1,k1)+dx1*dy1*dz1
       d(i2,j1,k1)=d(i2,j1,k1)+dx2*dy1*dz1
       d(i1,j2,k1)=d(i1,j2,k1)+dx1*dy2*dz1
       d(i2,j2,k1)=d(i2,j2,k1)+dx2*dy2*dz1
       d(i1,j1,k2)=d(i1,j1,k2)+dx1*dy1*dz2
       d(i2,j1,k2)=d(i2,j1,k2)+dx2*dy1*dz2
       d(i1,j2,k2)=d(i1,j2,k2)+dx1*dy2*dz2
       d(i2,j2,k2)=d(i2,j2,k2)+dx2*dy2*dz2
    enddo
 end subroutine cicmass

 subroutine powerspectrum(delta,ps)
    implicit none
    real, parameter :: ncr=nc
    integer, parameter :: kpt=nc/nt+min(1,mod(nc,nt))
             
    real ps(2,nc),delta(nc+2,nc,nc)
             
    integer it
    integer i,j,k
    integer k1,k2
    real kr,kx,ky,kz,w1,w2,pow
    real pst(3,nc,nt)
             
    real time1,time2
    call cpu_time(time1)


    !! Compute power spectrum
    !$omp parallel do default(shared) private(it,i,j,k) &
    !$omp& private(kr,kx,ky,kz,k1,k2,w1,w2,pow)
    do it=1,nt
       pst(:,:,it)=0
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          if (k .lt. nc/2+2) then
             kz=k-1
          else
             kz=k-1-nc
          endif
          do j=1,nc
             if (j .lt. nc/2+2) then
                ky=j-1
             else
                ky=j-1-nc
             endif
             do i=1,nc+2,2
                kx=(i-1)/2
                kr=sqrt(kx**2+ky**2+kz**2)
                if (kr .ne. 0) then
                   k1=ceiling(kr)
                   k2=k1+1
                   w1=k1-kr
                   w2=1-w1
                   pow=sum((delta(i:i+1,j,k)/ncr**3)**2)
                   pst(1,k1,it)=pst(1,k1,it)+w1*pow
                   pst(2,k1,it)=pst(2,k1,it)+w1*pow**2
                   pst(3,k1,it)=pst(3,k1,it)+w1
                   pst(1,k2,it)=pst(1,k2,it)+w2*pow
                   pst(2,k2,it)=pst(2,k2,it)+w2*pow**2
                   pst(3,k2,it)=pst(3,k2,it)+w2
                endif
             enddo
          enddo
       enddo 
    enddo  
    !$omp end parallel do

    !! Divide by weights
    !! ps(1,k) stores p(k)
    !! ps(2,k) stores standard deviation
    do k=1,nc
       if (pst(3,k,1) .eq. 0) then
          ps(:,k)=0
       else  
          ps(:,k)=pst(1:2,k,1)/pst(3,k,1)
          ps(2,k)=sqrt(abs(ps(2,k)-ps(1,k)**2))
          ps(1:2,k)=4*3.141596254*(k-1)**3*ps(1:2,k)
       endif
    enddo

    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called power spectrum'
    return
  end subroutine powerspectrum

#ifdef FFTW
  subroutine fft3d(a,c)
#else
  subroutine fft3d(a,n,c)
#endif
    implicit none
#ifdef FFTW
    include 'fftw3.f'
    real a(nc+2,nc,nc)
    integer(8), save :: plan,iplan
    logical :: first_fft
#else
    integer n
    real a(n+2,n,n)
    external sfft_3d
#endif
    character c
#ifdef FFTW
    data first_fft /.true./
    if (first_fft) then
      first_fft=.false.
      call sfftw_plan_dft_r2c_3d(plan, nc, nc, nc, a, a,FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan, nc, nc, nc, a, a,FFTW_MEASURE)
    endif
#endif

    if (c .eq. 'f') then
#ifdef FFTW
      call sfftw_execute_dft_r2c(plan,a,a)
#else
      call sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
#endif
    else
#ifdef FFTW
      call sfftw_execute_dft_c2r(iplan,a,a)
      a=a/real(nc)**3
#else
      call sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
#endif
    endif
    
    return
  end subroutine fft3d

end
