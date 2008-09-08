!! init.f90 written by Hy Trac
!! Last modified on Apr 14, 2005
!! Compile with:  f90 -fast -omp -lcxmlp -o init.x init.f90
!! merz 051004 :: Modified to use FFTW3 
!! Compile with: ifort init.f90 -fpp -O3 -DFFTW -I$MCKENZIE_FFTW_INC_PATH $MCKENZIE_FFTW_LIB_PATH/libfftw3f.a -o init.x
!! Power spectra have been checked
!! Correlation between gas and dm has been checked


program main
  implicit none

  include 'parameters.h'

  real, parameter :: redshift=zi
  real, parameter :: scalefactor=1/(1+redshift)

  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc
  real, parameter    :: npr=np
  
  !! Other parameters
  real, parameter :: pi=3.14159

  !! Gas arrays
  real, dimension(5,nc,nc,nc) :: u

  !! Dark matter arrays
  real, dimension(6,np,np,np) :: xvp

  !! Power spectrum arrays
  real, dimension(5,nk) :: tf
  real, dimension(2,nc) :: pkm,pkdm,pkg

  !! Force arrays
  real, dimension(nc+2,nc,nc)     :: noise
  real, dimension(nc+2,nc,nc)     :: deltam,deltab
  real, dimension(nc+2,nc,nc)     :: force,phi
  real, dimension(hc+1,hc+1,hc+1) :: wx,wy,wz

  !! Equivalence arrays to save memory
  equivalence(force,noise)
  

  !! Common block
  common u,xvp,deltam,deltab,force,phi,wx,wy,wz


  !$ call omp_set_num_threads(nt)
  call writeparams
#ifdef FFTW
!  initialize  fftw3
  call fft3d(noise,'f')
#endif

  call initvar

  call transferfnc
  call noisemap
  call deltafield

  call dm
  call writedm
  call gas
  !call writegas
  call writepowerspectra

contains


  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nt      ', nt
    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*)
    write(*,*) 'n_s      ',ns
    write(*,*) 'sigma_8  ',s8
    write(*,*) 'omega_m  ',omegam
    write(*,*) 'omega_l  ',omegal
    write(*,*) 'omega_b  ',omegab
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*) 'redshift ',redshift
    write(*,*)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams


  subroutine writedm
    implicit none
    character*60 :: fn
    
    real time1,time2
    call cpu_time(time1)


    !! Dark matter data
    !! xvp(1:3) = x,y,z in grid units from 0 to nc
    !! xvp(4:6) = vx,vy,vz in km/s

    fn=dir//'xvp.init'
    write(*,*) 'Writing ',fn
#ifdef BINARY
    open(11,file=fn,form='binary')
#else
    open(11,file=fn,form='unformatted')
#endif
    write(11) xvp
    close(11)


    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write dm'
    return
  end subroutine writedm


  subroutine writegas
    implicit none
    character*60 :: fn

    real time1,time2
    call cpu_time(time1)

    
    !! Gas data
    !! u(1) = delta
    !! u(2:4) = vx,vy,vz in km/s

    fn=dir//'u.init'
    write(*,*) 'Writing ',fn
#ifdef BINARY
    open(11,file=fn,form='binary')
#else
    open(11,file=fn,form='unformatted')
#endif
    write(11) u
    close(11)


    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write gas'
    return
  end subroutine writegas


  subroutine writepowerspectra
    implicit none
    integer      :: k
    real         :: kr
    character*60 :: fn

    real time1,time2
    call cpu_time(time1)


    !! Output power spectrum
    !! 1st column is k
    !! 2nd is matter p(k)
    !! 3rd is standard deviation
    !! 4th is gas p(k)
    !! 5th is standard deviation
    !! 6th is dm p(k)
    !! 7th is standard deviation
    fn=dir//'pk.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       write(11,*) kr,pkm(:,k),pkg(:,k),pkdm(:,k)
    enddo
    close(11)


    !! Output cmbfast power spectrum
    fn=dir//'pk0.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,2*nc+1
       kr=2*pi*(k-1)/(2*box)
       write(11,*) kr,power(kr,1,2),power(kr,1,3)
    enddo
    close(11)
    

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write power spectra'
    return
  end subroutine writepowerspectra


!!------------------------------------------------------------------!!


  subroutine transferfnc
    implicit none
    integer :: i,k
    real    :: kr,kmax
    real*8  :: v8

    real time1,time2
    call cpu_time(time1)


    !! Get transfer function from CMBFAST
    write(*,*) 'Reading ',fntf
    open(11,file=fntf)
    read(11,*) tf
    close(11)


    !! Compute \Delta^2
    do k=1,nk
       kr     =tf(1,k)
       tf(2,k)=kr**(3+ns)*tf(2,k)**2/(2*pi**2)
       tf(3,k)=kr**(3+ns)*tf(3,k)**2/(2*pi**2)
    enddo


    !! Compute dk
    tf(4,1)=tf(1,2)/2
    do k=2,nk-1
       tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
    enddo
    tf(4,nk)=tf(1,nk)-tf(1,nk-1)


    !! Compute variance in 8 h^-1 Mpc spheres
    v8=0
    kmax=2*pi*sqrt(3.)*hc/box
    do k=1,nk
       if (tf(1,k) .gt. kmax) exit
       v8=v8+tf(2,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
    enddo


    !! Normalize to \sigma_8
    !! Include growth factor
    tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*Dgrow(scalefactor)**2


    !! tf(1,i) stores k
    !! tf(2,i) stores \Delta^2_m
    !! tf(3,i) stores \Delta^2_b
    !! tf(4,i) stores dk


    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called transfer fnc'
    return
  end subroutine transferfnc


  subroutine noisemap
    implicit none
    integer      :: i,j,k,one
    real         :: x,x1,x2
    character*60 :: fn
!    integer, dimension(1,nc) :: iseed
    integer,allocatable,dimension(:)  :: iseed
    real, dimension(nc)      :: rseed

    real time1,time2
    call cpu_time(time1)


    one=1
    if (.true.) then
       write(*,*) 'Generating seeds'
       call random_seed
       call random_seed(size=i)
       allocate(iseed(i))
       do k=1,i !nc
          call random_number(x)
          iseed(k)=int(x*huge(0))
       enddo
    else
       fn=dir//'seed.init'
       write(*,*) 'Reading ',fn
       open(11,file=fn)
       call random_seed
       call random_seed(size=i)
       do k=1,i
          read(11,*) j,iseed(k),x
       enddo
       close(11)
    endif

    call random_seed(put=iseed(1:i))
    !! Generate random reals between 0 and 1
    !$omp parallel do default(shared) private(k)
    do k=1,nc
!       call random_seed(size=one)
!       call random_seed(put=iseed(:,k))
       call random_number(noise(:,:,k))
       rseed(k)=noise(1,1,k)
    enddo
    !$omp end parallel do
 

    fn=dir//'seed.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn)
    do k=1,i
       write(11,*) k,iseed(k),rseed(k)
    enddo
    close(11)
    

    !! Generate random Gaussian numbers
    !$omp parallel do default(shared) private(i,j,k,x1,x2)
    do k=1,nc
       do j=1,nc
          do i=1,nc,2
             x1=2*pi*noise(i,j,k)
             x2=sqrt(-2*log(noise(i+1,j,k)))
             noise(i,j,k)=x2*cos(x1)
             noise(i+1,j,k)=x2*sin(x1)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT white noise field
!    call fft3d(noise,noise,nc,'f')
#ifdef FFTW
    print *,sum(noise)
    call fft3d(noise,'f')
    print *,sum(noise)
#else
    call fft3d(noise,nc,'f')
#endif

    !! Improve sampling at large scales (small k)
    !! Use only for diagnostic purposes
    if (.false.) call normalizenoise


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called noisemap'
    return
  end subroutine noisemap


  subroutine normalizenoise
    implicit none
    real, parameter :: err=0.1

    integer :: i,j,k,n
    integer :: k1,k2
    real    :: kr,kx,ky,kz
    real    :: pk,wk,norm
    real*8  :: A(2,hc,nc)

    real time1,time2
    call cpu_time(time1)


    !$omp parallel do default(shared) private(i,j,k,n) &
    !$omp& private(kr,kx,ky,kz)
    do k=1,nc
       A(:,:,k)=0
       if (k .lt. hc+2) then
          kz=k-1
       else
          kz=k-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          do i=1,nc+2,2
             kx=(i-1)/2
             kr=sqrt(kx**2+ky**2+kz**2)
             n=int(kr)
             if (n .gt. 0 .and. n .le. hc) then
                A(1,n,k)=A(1,n,k)+sum(noise(i:i+1,j,k)**2)
                A(2,n,k)=A(2,n,k)+1
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Compute normalization
    do k=2,nc
       A(:,:,1)=A(:,:,1)+A(:,:,k)
    enddo
    do k=1,hc
       norm=A(1,k,1)/A(2,k,1)/ncr**3
       A(1,k,1)=1
       if (norm .gt. 1+err) A(1,k,1)=sqrt(norm/(1+err))
       if (norm .lt. 1-err) A(1,k,1)=sqrt(norm/(1-err))
       write(*,*) k,real(norm/A(1,k,1)**2)
    enddo


    !$omp parallel do default(shared) private(i,j,k,n) &
    !$omp& private(kr,kx,ky,kz)
    do k=1,nc
       if (k .lt. hc+2) then
          kz=k-1
       else
          kz=k-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          do i=1,nc+2,2
             kx=(i-1)/2
             kr=sqrt(kx**2+ky**2+kz**2)
             n=int(kr)
             if (n .gt. 0 .and. n .le. hc) then
                noise(i:i+1,j,k)=noise(i:i+1,j,k)/A(1,n,1)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called normalize noise'
    return
  end subroutine normalizenoise


  subroutine deltafield
    implicit none
    integer :: i,j,k
    real    :: kr,kx,ky,kz
    real    :: powb,powm

    real time1,time2
    call cpu_time(time1)


    !! Interpolate \Delta^2 to grid
    !! Determine Fourier modes \delta(k)

    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(kr,kx,ky,kz,powb,powm)
    do k=1,nc
       if (k .lt. hc+2) then
          kz=k-1
       else
          kz=k-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          do i=1,nc+2,2
             kx=(i-1)/2
             kr=sqrt(kx**2+ky**2+kz**2)
             if (kr .eq. 0) then
                deltam(i:i+1,j,k)=0
                deltab(i:i+1,j,k)=0
             else
                powm=power(2*pi*kr/box,1,2)/(4*pi*kr**3)
                powb=power(2*pi*kr/box,1,3)/(4*pi*kr**3)
                deltam(i:i+1,j,k)=sqrt(powm*ncr**3)*noise(i:i+1,j,k)
                deltab(i:i+1,j,k)=sqrt(powb*ncr**3)*noise(i:i+1,j,k)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Compute matter and gas power spectrum
    call powerspectrum(deltam,pkm)
    call powerspectrum(deltab,pkg)


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called delta field'
    return
  end subroutine deltafield


!!------------------------------------------------------------------!!


  subroutine forcefield(delta)
    implicit none
    logical, parameter :: correct=.false.
    real, dimension(nc+2,nc,nc) :: delta

    integer :: i,j,k
    integer :: i1,j1,k1,ii,ir
    integer :: im,ip,jm,jp,km,kp
    real    :: r,x,y,z
    real    :: ksq,kx,ky,kz

    real time1,time2
    call cpu_time(time1)


    !! Construct potential kernel in Fourier space
    !$omp parallel do default(shared) private(i,j,k,ip) &
    !$omp& private(ksq,kx,ky,kz)
    do k=1,nc
       if (k .lt. hc+2) then
          kz=k-1
       else
          kz=k-1-nc
       endif
       kz=2*pi*kz/ncr
       !kz=2*sin(pi*kz/ncr)
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          ky=2*pi*ky/ncr
          !ky=2*sin(pi*ky/ncr)
          do i=1,nc+2,2
             ip=i+1
             kx=(i-1)/2
             kx=2*pi*kx/ncr
             !kx=2*sin(pi*kx/ncr)
             ksq=kx**2+ky**2+kz**2
             if (ksq .eq. 0) then
                phi(i:ip,j,k)=0
             else
                phi(i,j,k)=-4*pi/ksq
                phi(ip,j,k)=0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Inverse FFT potential kernel
!    call fft3d(phi,phi,nc,'b')
#ifdef FFTW
    call fft3d(phi,'b')
#else
    call fft3d(phi,nc,'b')
#endif

    !! Construct Ewald force kernel in real space 
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(im,ip,jm,jp,km,kp) &
    !$omp& private(r,x,y,z)
    do k=1,nc
       kp=mod(k,nc)+1
       km=mod(k-2+nc,nc)+1
       if (k .lt. hc+2) then
          z=k-1
       else
          z=k-1-nc
       endif
       do j=1,nc
          jp=mod(j,nc)+1
          jm=mod(j-2+nc,nc)+1
          if (j .lt. hc+2) then
             y=j-1
          else
             y=j-1-nc
          endif
          do i=1,nc
             ip=mod(i,nc)+1
             im=mod(i-2+nc,nc)+1
             if (i .lt. hc+2) then
                x=i-1
             else
                x=i-1-nc
             endif
             r=sqrt(x**2+y**2+z**2)

             force(i,j,k)=(phi(im,j,k)-phi(ip,j,k))/2
             if (correct) then
                if (r .ne. 0 .and. r .le. 8) force(i,j,k)=-x/r**3
             endif
             force(1,j,k)=0
             force(hc+1,j,k)=0
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT force kernel
!    call fft3d(force,force,nc,'f')
#ifdef FFTW
    call fft3d(force,'f')
#else
    call fft3d(force,nc,'f')
#endif

    !! Force kernel is odd and imaginary
    !$omp parallel do default(shared) private(i,j,k)
    do k=1,hc+1
       do j=1,hc+1
          do i=1,hc+1
             wx(i,j,k)=force(2*i,j,k)
             wy(j,i,k)=force(2*i,j,k)
             wz(k,j,i)=force(2*i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Complex multiply to get Fourier space forces
    !! Temporarily store forces in array u


    !! x-component
    !$omp parallel default(shared) private(i,j,k,j1,k1,ir,ii)
    !$omp do
    do k=1,hc+1
       do j=1,hc+1
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wx(i,j,k)
             force(ii,j,k)= delta(ir,j,k)*wx(i,j,k)
          enddo
       enddo
       do j=hc+1,nc
          j1=nc+2-j
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wx(i,j1,k)
             force(ii,j,k)= delta(ir,j,k)*wx(i,j1,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=hc+2,nc
       k1=nc+2-k
       do j=1,hc+1
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wx(i,j,k1)
             force(ii,j,k)= delta(ir,j,k)*wx(i,j,k1)
          enddo
       enddo
       do j=hc+2,nc
          j1=nc+2-j
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wx(i,j1,k1)
             force(ii,j,k)= delta(ir,j,k)*wx(i,j1,k1)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

 
    !! Inverse FFT force field
!    call fft3d(force,force,nc,'b')
#ifdef FFTW
    call fft3d(force,'b')
#else
    call fft3d(force,nc,'b')
#endif

    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc
       u(2,:,:,k)=force(1:nc,:,k)/(4*pi)
    enddo
    !$omp end parallel do


    !! y-component
    !$omp parallel default(shared) private(i,j,k,j1,k1,ir,ii)
    !$omp do
    do k=1,hc+1
       do j=1,hc+1
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wy(i,j,k)
             force(ii,j,k)= delta(ir,j,k)*wy(i,j,k)
          enddo
       enddo
       do j=hc+1,nc
          j1=nc+2-j
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)= delta(ii,j,k)*wy(i,j1,k)
             force(ii,j,k)=-delta(ir,j,k)*wy(i,j1,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=hc+2,nc
       k1=nc+2-k
       do j=1,hc+1
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wy(i,j,k1)
             force(ii,j,k)= delta(ir,j,k)*wy(i,j,k1)
          enddo
       enddo
       do j=hc+2,nc
          j1=nc+2-j
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)= delta(ii,j,k)*wy(i,j1,k1)
             force(ii,j,k)=-delta(ir,j,k)*wy(i,j1,k1)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel 


    !! Inverse FFT force field
!    call fft3d(force,force,nc,'b')
#ifdef FFTW
    call fft3d(force,'b')
#else
    call fft3d(force,nc,'b')
#endif

    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc
       u(3,:,:,k)=force(1:nc,:,k)/(4*pi)
    enddo
    !$omp end parallel do


    !! z-component
    !$omp parallel default(shared) private(i,j,k,j1,k1,ir,ii)
    !$omp do
    do k=1,hc+1
       do j=1,hc+1
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wz(i,j,k)
             force(ii,j,k)= delta(ir,j,k)*wz(i,j,k)
          enddo
       enddo
       do j=hc+1,nc
          j1=nc+2-j
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)=-delta(ii,j,k)*wz(i,j1,k)
             force(ii,j,k)= delta(ir,j,k)*wz(i,j1,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=hc+2,nc
       k1=nc+2-k
       do j=1,hc+1
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)= delta(ii,j,k)*wz(i,j,k1)
             force(ii,j,k)=-delta(ir,j,k)*wz(i,j,k1)
          enddo
       enddo
       do j=hc+2,nc
          j1=nc+2-j
          do i=1,hc+1
             ii=2*i
             ir=ii-1
             force(ir,j,k)= delta(ii,j,k)*wz(i,j1,k1)
             force(ii,j,k)=-delta(ir,j,k)*wz(i,j1,k1)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel 


    !! Inverse FFT force field
!    call fft3d(force,force,nc,'b')
#ifdef FFTW
    call fft3d(force,'b')
#else
    call fft3d(force,nc,'b')
#endif

    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc
       u(4,:,:,k)=force(1:nc,:,k)/(4*pi)
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called force field'
    return
  end subroutine forcefield


  subroutine potentialfield(delta)
    implicit none
    logical, parameter :: correct=.false.
    real, dimension(nc+2,nc,nc) :: delta

    integer :: i,j,k
    integer :: im,ip,jm,jp,km,kp
    real    :: r,x,y,z
    real    :: kr,ksq,kx,ky,kz
    real    :: phi8

    real time1,time2
    call cpu_time(time1)


    !! Construct uncorrected potential kernel in Fourier space
    !$omp parallel do default(shared) private(i,j,k,ip) &
    !$omp& private(ksq,kx,ky,kz)
    do k=1,nc
       if (k .lt. hc+2) then
          kz=k-1
       else
          kz=k-1-nc
       endif
       kz=2*sin(pi*kz/ncr)
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          ky=2*sin(pi*ky/ncr)
          do i=1,nc+2,2
             ip=i+1
             kx=(i-1)/2
             kx=2*sin(pi*kx/ncr)
             ksq=kx**2+ky**2+kz**2
             if (ksq .eq. 0) then
                phi(i:ip,j,k)=0
             else
                phi(i,j,k)=-4*pi/ksq
                phi(ip,j,k)=0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    if (correct) then
       !! Forward FFT potential kernel
!       call fft3d(phi,phi,nc,'b')
#ifdef FFTW
       call fft3d(phi,'b')
#else
       call fft3d(phi,nc,'b')
#endif

       !! Compute offset at r=8
       phi8=(phi(9,1,1)+phi(nc-7,1,1) &
            +phi(1,9,1)+phi(1,nc-7,1) &
            +phi(1,1,9)+phi(1,1,nc-7))/6
       

       !! Construct Ewald potential kernel in real space
       !$omp parallel do default(shared) private(i,j,k) &
       !$omp& private(r,x,y,z)
       do k=1,nc
          if (k .lt. hc+2) then
             z=k-1
          else
             z=k-1-nc
          endif
          do j=1,nc
             if (j .lt. hc+2) then
                y=j-1
             else
                y=j-1-nc
             endif
             do i=1,nc
                if (i .lt. hc+2) then
                   x=i-1
                else
                   x=i-1-nc
                endif
                r=sqrt(x**2+y**2+z**2)
                if (r .gt. 8) then
                   phi(i,j,k)=phi(i,j,k)-(phi8+1/8.)
                else
                   if (r .eq. 0) then
                      phi(i,j,k)=-2.5
                   else
                      phi(i,j,k)=-1/r
                   endif
                endif
             enddo
          enddo
       enddo
       !$omp end parallel do


       !! Forward FFT potential kernel
!       call fft3d(phi,phi,nc,'f')
#ifdef FFTW
       call fft3d(phi,'f')
#else
       call fft3d(phi,nc,'f')
#endif
    endif


    !! Complex multiply density field with potential kernel
    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc
       do j=1,nc
          do i=1,nc+2,2
             phi(i:i+1,j,k)=delta(i:i+1,j,k)*phi(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Inverse FFT potential field
!    call fft3d(phi,phi,nc,'b')
#ifdef FFTW
    call fft3d(phi,'b')
#else
    call fft3d(phi,nc,'b')
#endif

    !! Finite-difference to get displacement field
    !! Temporarily store in array u
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(im,ip,jm,jp,km,kp)
    do k=1,nc
       kp=mod(k,nc)+1
       km=mod(k-2+nc,nc)+1
       do j=1,nc
          jp=mod(j,nc)+1
          jm=mod(j-2+nc,nc)+1       
          do i=1,nc
             ip=mod(i,nc)+1
             im=mod(i-2+nc,nc)+1
             u(2,i,j,k)=(phi(im,j,k)-phi(ip,j,k))/2/(4*pi)
             u(3,i,j,k)=(phi(i,jm,k)-phi(i,jp,k))/2/(4*pi)
             u(4,i,j,k)=(phi(i,j,km)-phi(i,j,kp))/2/(4*pi)
          enddo
       enddo
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called potential field'
    return
  end subroutine potentialfield


!!------------------------------------------------------------------!!


  subroutine dm
    implicit none
    integer :: i,j,k
    integer :: i1,j1,k1
    real    :: d,dmin,dmax
    real*8  :: dsum,dvar
    real, dimension(3) :: dis

    real time1,time2
    call cpu_time(time1)


    !! Calculate force field for dark matter density field
    if (.false.) then
       call forcefield(deltam)
    else
       call potentialfield(deltam)
    endif


    !! Displace particles
    !! Displacement field is stored in u(2:4,i,j,k)
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(i1,j1,k1)
    do k=1,np
       k1=(nc/np)*(k-1)+1
       do j=1,np
          j1=(nc/np)*(j-1)+1
          do i=1,np
             i1=(nc/np)*(i-1)+1
             xvp(1,i,j,k)=mod(u(2,i1,j1,k1)+(i1-0.5)+ncr,ncr)
             xvp(2,i,j,k)=mod(u(3,i1,j1,k1)+(j1-0.5)+ncr,ncr)
             xvp(3,i,j,k)=mod(u(4,i1,j1,k1)+(k1-0.5)+ncr,ncr)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Initialized density field to be zero
    !$omp parallel do default(shared) private(k)
    do k=1,nc
       deltam(:,:,k)=0
    enddo
    !$omp end parallel do


    !! Determine velocity for each particle
    !! Assign masses to grid to compute dm power spectrum
    !$omp parallel do default(shared) private(j,k)
    do k=1,np
       do j=1,np
          call cicvelocity(xvp(:,:,j,k))
       enddo
    enddo
    !$omp end parallel do

    !! scale for cubep3m
    do k=1,np
      do j=1,np
        do i=1,np
          xvp(:,i,j,k)=xvp(:,i,j,k)/5.0
        enddo
      enddo
    enddo

    !! Convert dm density field to delta field
    dmin=0
    dmax=0
    dsum=0
    dvar=0

    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) &
    !$omp& reduction(max:dmax) &
    !$omp& reduction(+:dsum,dvar)
    do k=1,nc
       do j=1,nc
          do i=1,nc
             deltam(i,j,k)=deltam(i,j,k)-1
             d=deltam(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo
    !$omp end parallel do


    dsum=dsum/nc**3
    dvar=sqrt(dvar/nc**3)
    write(*,*)
    write(*,*) 'DM min    ',dmin
    write(*,*) 'DM max    ',dmax
    write(*,*) 'Delta sum ',real(dsum)
    write(*,*) 'Delta var ',real(dvar)
    write(*,*)


    !! Forward FFT dm delta field
!    call fft3d(deltam,deltam,nc,'f')
#ifdef FFTW
    call fft3d(deltam,'f')
#else
    call fft3d(deltam,nc,'f')
#endif

    !! Compute dm power spectrum
    call powerspectrum(deltam,pkdm)


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine dm


  subroutine cicvelocity(xvp)
    implicit none
    real, parameter       :: mp=(ncr/np)**3
    real, dimension(6,np) :: xvp

    integer :: i
    integer :: i1,i2,j1,j2,k1,k2
    real    :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2
    real    :: vf,v(3)

    vf=vfactor(scalefactor)

    do i=1,np
       x=mod(xvp(1,i)-0.5+ncr,ncr)
       y=mod(xvp(2,i)-0.5+ncr,ncr)
       z=mod(xvp(3,i)-0.5+ncr,ncr)

       i1=int(x)+1
       i2=mod(i1,nc)+1
       dx1=i1-x
       dx2=1-dx1
       j1=int(y)+1
       j2=mod(j1,nc)+1
       dy1=j1-y
       dy2=1-dy1
       k1=int(z)+1
       k2=mod(k1,nc)+1
       dz1=k1-z
       dz2=1-dz1

       v=0
       v=v+u(2:4,i1,j1,k1)*dx1*dy1*dz1
       v=v+u(2:4,i2,j1,k1)*dx2*dy1*dz1
       v=v+u(2:4,i1,j2,k1)*dx1*dy2*dz1
       v=v+u(2:4,i2,j2,k1)*dx2*dy2*dz1
       v=v+u(2:4,i1,j1,k2)*dx1*dy1*dz2
       v=v+u(2:4,i2,j1,k2)*dx2*dy1*dz2
       v=v+u(2:4,i1,j2,k2)*dx1*dy2*dz2
       v=v+u(2:4,i2,j2,k2)*dx2*dy2*dz2
       xvp(4:6,i)=v*vf

       dz1=mp*dz1
       dz2=mp*dz2
       deltam(i1,j1,k1)=deltam(i1,j1,k1)+dx1*dy1*dz1
       deltam(i2,j1,k1)=deltam(i2,j1,k1)+dx2*dy1*dz1
       deltam(i1,j2,k1)=deltam(i1,j2,k1)+dx1*dy2*dz1
       deltam(i2,j2,k1)=deltam(i2,j2,k1)+dx2*dy2*dz1
       deltam(i1,j1,k2)=deltam(i1,j1,k2)+dx1*dy1*dz2
       deltam(i2,j1,k2)=deltam(i2,j1,k2)+dx2*dy1*dz2
       deltam(i1,j2,k2)=deltam(i1,j2,k2)+dx1*dy2*dz2
       deltam(i2,j2,k2)=deltam(i2,j2,k2)+dx2*dy2*dz2
    enddo

    return
  end subroutine cicvelocity


!!--------------------------------------------------------------!!


  subroutine gas
    implicit none
    integer :: i,j,k
    real    :: vf
    real    :: d,dmin,dmax
    real*8  :: dsum,dvar

    real time1,time2
    call cpu_time(time1)


    !! Calculate force field for baryon density field
    if (.false.) then
       call forcefield(deltab)
    else
       call potentialfield(deltab)
    endif


    !! Forward FFT density field
!    call fft3d(deltab,deltab,nc,'b')
#ifdef FFTW
    call fft3d(deltab,'b')
#else
    call fft3d(deltab,nc,'b')
#endif

    !! Compute density, momentum, and total energy
    dmin=0
    dmax=0
    dsum=0
    dvar=0
    vf=vfactor(scalefactor)

    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) &
    !$omp& reduction(max:dmax) &
    !$omp& reduction(+:dsum,dvar)
    do k=1,nc
       do j=1,nc
          do i=1,nc
             d=deltab(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
             d=omegab/omegam*(1+d)
             u(1,i,j,k)=d
             u(2:4,i,j,k)=d*u(2:4,i,j,k)*vf
             u(5,i,j,k)=sum(u(2:4,i,j,k)**2)/d/2
          enddo
       enddo
    enddo
    !$omp end parallel do


    dsum=dsum/nc**3
    dvar=sqrt(dvar/nc**3)
    write(*,*)
    write(*,*) 'Gas min   ',dmin
    write(*,*) 'Gas max   ',dmax
    write(*,*) 'Delta sum ',real(dsum)
    write(*,*) 'Delta var ',real(dvar)
    write(*,*)


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called gas'
  end subroutine gas


!!--------------------------------------------------------------!!


  subroutine powerspectrum(delta,pk)
    implicit none
    integer, parameter          :: kpt=nc/nt+min(1,mod(nc,nt))
    real, dimension(2,nc)       :: pk
    real, dimension(nc+2,nc,nc) :: delta

    integer :: i,j,k
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow
    real, dimension(3,nc,nc) :: pkt

    real time1,time2
    call cpu_time(time1)


    !! Compute power spectrum
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(kr,kx,ky,kz,k1,k2,w1,w2,pow)
    do k=1,nc
       pkt(:,:,k)=0
       if (k .lt. hc+2) then
          kz=k-1
       else
          kz=k-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
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
                pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                pkt(3,k1,k)=pkt(3,k1,k)+w1
                pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                pkt(3,k2,k)=pkt(3,k2,k)+w2
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Merge power spectrum from threads
    !! Divide by weights
    !! pk(1,k) stores gas pk(k)
    !! pk(2,k) stores standard deviation
    do k=2,nc
       pkt(:,:,1)=pkt(:,:,1)+pkt(:,:,k)
    enddo
    do k=1,nc
       if (pkt(3,k,1) .eq. 0) then
          pk(:,k)=0
       else
          pk(:,k)=pkt(1:2,k,1)/pkt(3,k,1)
          pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pkt(3,k,1)-1)))
          pk(1:2,k)=4*pi*(k-1)**3*pk(1:2,k)
       endif
    enddo


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called power spectrum'
    return
  end subroutine powerspectrum


!!------------------------------------------------------------------!!


  subroutine initvar
    implicit none
    integer :: k

    real time1,time2
    call cpu_time(time1)


    !! Initialize large arrays in parallel
    !$omp parallel default(shared) private(k)
    !$omp do
    do k=1,nc
       u(:,:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,np
       xvp(:,:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc
       noise(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc
       deltam(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc
       deltab(:,:,k)=0
    enddo
    !$omp end do
    !$omp end parallel


    call cpu_time(time2)
    time2=(time2-time1)
    write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar


!!------------------------------------------------------------------!!


  function power(kr,ix,iy)
   implicit none
    real    :: kr
    integer :: ix,iy

    integer :: i,i1,i2
    real    :: x,y,x1,x2,y1,y2
    real    :: power

    i1=1
    i2=nk
    do while (i2-i1 .gt. 1)
       i=(i1+i2)/2
       if (kr .gt. tf(ix,i)) then
          i1=i
       else
          i2=i
       endif
    enddo

    x1=log(tf(ix,i1))
    y1=log(tf(iy,i1))
    x2=log(tf(ix,i2))
    y2=log(tf(iy,i2))
    x =log(kr)
    y =y1+(y2-y1)*(x-x1)/(x2-x1)
    power=exp(y)
    
    return
  end function power


  function Dgrow(a)
    implicit none
    real, parameter :: om=omegam
    real, parameter :: ol=omegal
    real a
    real Dgrow

    real g,ga,hsq,oma,ola

    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    Dgrow=a*ga/g
    return
  end function Dgrow


  function vfactor(a)
    implicit none
    real :: a

    real :: H,km,lm
    real :: vfactor

    lm=omegal/omegam
    km=(1-omegam-omegal)/omegam
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H

    return
  end function vfactor


  function tophat(x)
    implicit none
    real :: x,tophat

    if (x .ne. 0) then
       tophat=3*(sin(x)-cos(x)*x)/x**3
    else
       tophat=1
    endif

    return
  end function tophat

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
    external sfft
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

end program main
