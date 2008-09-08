
!! init.f90 written by Hy Trac on Sept 23, 2003
!! Parallelized :: Hugh Merz :: 041216
!! efc -fpp2 init_new.f90 /opt/fftw-3.0.1/lib/libfftw3f.a -openmp -tpp1 -O3 -o init.x


program main
  implicit none

  !! Cosmo parameters
  real, parameter :: box=100
  real, parameter :: s8=0.84
  real, parameter :: omegam=0.30
  real, parameter :: omegal=0.70
  real, parameter :: redshift=100
  real, parameter :: scalefactor=1/(1+redshift)

  !! nt is the number of threads
  !! dir is the directory for files
  integer, parameter :: nt=4
  character(*), parameter :: dir='./'

  !! nc is the number of cells per box length
  !! np should be set to nc (1:1) or hc (1:2)
  integer, parameter :: nc=128
  integer, parameter :: hc=nc/2
  integer, parameter :: np=hc

  !! nk is the length of the tf.init file (from cmbfast)
  character(*), parameter :: tffn='tf.wmap'
  integer, parameter :: nk=380

  !! Other parameters
  real, parameter :: pi=3.14159

  !! Power spectrum array
  real, dimension(5,nk) :: tf
  real, dimension(4,nc) :: psinit,psg,psdm

  !! Large arrays

!HM - these need to be decomposed.  Cubic?
!HM - use cubic decomposition fftw wrapper

  real, dimension(nc+2,nc,nc) :: init,noise,phi,den
  real, dimension(3,np,np,np) :: xp

  equivalence(init,noise)
  equivalence(phi,den)
  common init,phi,xp

  call OMP_SET_NUM_THREADS(nt)

  call writeparams
  call initvar

  call noisemap
  call transferfnc
  call powerspectrum
  call kernel
  call zeldovich
  call densitystatistics

  call writeps
  call writexp
  call writedelta


contains


  subroutine initvar
    implicit none
    integer k

    real time1,time2
    call cpu_time(time1)

    !! Initialize large arrays in parallel
    !$omp parallel default(shared) private(k)
    !$omp do
    do k=1,nc
       init(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc
       phi(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,np
       xp(:,:,:,k)=0
    enddo
    !$omp end do
    !$omp end parallel


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar


  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*) 'redshift ',redshift
    write(*,*) 'box      ',box
    write(*,*) 'sigma_8  ',s8
    write(*,*) 'omega_m  ',omegam
    write(*,*) 'omega_l  ',omegal
    write(*,*)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams


  subroutine writexp
    implicit none
    character*40 fn
    
    real time1,time2
    call cpu_time(time1)

    fn=dir//'xp.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,form='unformatted')
    write(11) xp
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write xp'
    return
  end subroutine writexp


  subroutine writedelta
    implicit none
    character*40 fn

    real time1,time2
    call cpu_time(time1)

    fn=dir//'od.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,form='unformatted')
    write(11) init(1:nc,:,:)
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write delta'
    return
  end subroutine writedelta


  subroutine writeps
    implicit none
    integer k
    real kr,Da,pow
    character*40 fn

    real time1,time2
    call cpu_time(time1)


    !! Calculate growth factor at scalefactor a
    Da=scalefactor*grow(scalefactor)/grow(1.)


    !! Output power spectrum
    !! First column is k
    !! Second is \Delta^2 (unsampled)
    !! Third is gas p(k)
    !! Fourth is standard deviation
    !! Fifth is dm p(k)
    !! Sixth is standard deviation

    fn=dir//'ps1.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       pow=power(kr)*Da**2
       write(11,*) kr,pow,psinit(1:2,k),psg(1:2,k),psdm(1:2,k)
    enddo
    close(11)

    fn=dir//'ps2.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       pow=power(kr)*Da**2
       write(11,*) kr,pow,psinit(3:4,k),psg(3:4,k),psdm(3:4,k)
    enddo
    close(11)


    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write ps'
    return
  end subroutine writeps


!!------------------------------------------------------------------!!


  subroutine noisemap
    implicit none
    integer, parameter :: kpt=nc/nt
    logical, parameter :: fast=.false.

    integer it
    integer i,j,k,one
    real x,x1,x2
    integer iseed(1,nt)

    real time1,time2
    call cpu_time(time1)


    one=1
    !! Generate random reals between 0 and 1
    if (fast) then
       !$omp parallel do default(shared) private(it,k,x)
       do it=1,nt
          call cpu_time(x)
          iseed(:,it)=it*ceiling(1000*x)+it
          call random_seed(size=one)
          call random_seed(put=iseed(:,it))
          call random_seed(get=iseed(:,it))
          do k=1+(it-1)*kpt,min(nc,it*kpt),2
             call random_number(noise(1:nc,:,k))
             noise(1:nc,:,k+1)=transpose(noise(1:nc,:,k))
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do default(shared) private(it,k,x)
       do it=1,nt
          call cpu_time(x)
          iseed(:,it)=it*ceiling(1000*x)+it
          call random_seed(size=one)
          call random_seed(put=iseed(:,it))
          call random_seed(get=iseed(:,it))
          do k=1+(it-1)*kpt,min(nc,it*kpt)
             call random_number(noise(1:nc,:,k))
          enddo
       enddo
       !$omp end parallel do
    endif


    !! Generate random Gaussian numbers
    !$omp parallel do default(shared) private(i,j,k,x1,x2)
    do k=1,nc,2
       do j=1,nc
          do i=1,nc,2
             x1=2*pi*noise(i,j,k)
             x2=sqrt(-2*log(noise(i+1,j,k)))
             noise(i,j,k)=x2*cos(x1)
             noise(i+1,j,k)=x2*sin(x1)
          enddo
          noise(nc+1:,j,k)=0
       enddo
       do j=1,nc
          do i=1,nc,2
             x1=2*pi*noise(i+1,j,k+1)
             x2=sqrt(-2*log(noise(i,j,k+1)))
             noise(i,j,k+1)=x2*cos(x1)
             noise(i+1,j,k+1)=x2*sin(x1)
          enddo
          noise(nc+1:,j,k+1)=0
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT white noise field
    call fft3d(noise,nc,'f')


    !! Check noise
    !! This improves sampling at large scales (small k)
    call checknoise


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called noisemap'
    return
  end subroutine noisemap


  subroutine checknoise
    implicit none
    integer, parameter :: kc=hc
    integer, parameter :: kpt=nc/nt
    real, parameter :: norm=nc**3

    integer it,i,j,k
    integer ksq,kx,ky,kz
    integer k1,k2
    real ps,ws
    real*8 A(2,(kc+1)**2,nt)

    real time1,time2
    call cpu_time(time1)


    !$omp parallel do default(shared) private(it,i,j,k) &
    !$omp& private(ksq,kx,ky,kz)
    do it=1,nt
       A(:,:,it)=0
       do k=1+(it-1)*kpt,it*kpt
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
                ksq=kx**2+ky**2+kz**2
                if (ksq .gt. 0 .and. ksq .le. kc**2) then
                   A(1,ksq,it)=A(1,ksq,it)+sum(noise(i:i+1,j,k)**2)
                   A(2,ksq,it)=A(2,ksq,it)+1
                endif
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Merge threads
    do it=2,nt
       A(:,:,1)=A(:,:,1)+A(:,:,it)
    enddo


    !! Compute normalization
    do i=1,kc
       k1=i**2
       k2=(i+1)**2-1
       ps=sum(A(1,k1:k2,1))
       ws=sum(A(2,k1:k2,1))
       if (ws .ne. 0) A(1,k1:k2,1)=sqrt(ps/ws/norm)
       write(*,*) k1,k2,sqrt(1.*k1),sqrt(1.*k2),A(1,k1,1)
    enddo


    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(ksq,kx,ky,kz)
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
             ksq=kx**2+ky**2+kz**2
             if (ksq .gt. 0 .and. ksq .le. kc**2) then
                noise(i:i+1,j,k)=noise(i:i+1,j,k)/A(1,ksq,1)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called check noise'
    return
  end subroutine checknoise


!!--------------------------------------------------------------!!


  subroutine transferfnc
    implicit none
    integer i
    real kmax,norm
    real*8 v8

    real time1,time2
    call cpu_time(time1)


    !! Get transfer function from CMBfast
    write(*,*) 'Reading ',tffn
    open(11,file=tffn)
    read(11,*) tf
    close(11)


    !! Compute \Delta^2
    do i=1,nk
       tf(3,i)=tf(1,i)**(3+1)*tf(2,i)**2/(2*pi**2)
    enddo


    !! Compute dk
    tf(4,1)=tf(1,2)/2
    do i=2,nk-1
       tf(4,i)=(tf(1,i+1)-tf(1,i-1))/2
    enddo
    tf(4,nk)=tf(1,nk)-tf(1,nk-1)


    !! Compute variance in 8 h^-1 Mpc spheres
    v8=0
    kmax=2*pi*sqrt(3.)*hc/box
    do i=1,nk
       if (tf(1,i) .gt. kmax) exit
       v8=v8+tf(3,i)*tophat(tf(1,i)*8)**2*tf(4,i)/tf(1,i)
    enddo


    !! Normalize to \sigma_8
    tf(3,:)=tf(3,:)*(s8**2/v8)


    !! tf(1,i) stores k
    !! tf(2,i) stores the transfer function
    !! tf(3,i) stores \Delta^2
    !! tf(4,i) stores dk


    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called transfer fnc'
    return
  end subroutine transferfnc


  subroutine powerspectrum
    implicit none
    real, parameter :: ncr=nc
    integer, parameter :: kpt=nc/nt

    integer it,i,j,k
    integer k1,k2
    real kx,ky,kz,kr,w1,w2
    real Da,delta,p0,p1,p2,p3,p4
    real*8 ps0(5,nc,nt),pst(5,nc,nt)

    real time1,time2
    call cpu_time(time1)


    !! Calculate growth factor at scalefactor a
    Da=scalefactor*grow(scalefactor)/grow(1.)


    !! Interpolate \Delta^2 to grid
    !! Determine Fourier modes \delta(k)
    !! Compute both sampled and unsampled power spectra

    !$omp parallel do default(shared) private(it,i,j,k) &
    !$omp& private(kr,kx,ky,kz,k1,k2,w1,w2,delta,p0,p1,p2,p3,p4)
    do it=1,nt
       pst(:,:,it)=0
       ps0(:,:,it)=0
       do k=1+(it-1)*kpt,min(nc,it*kpt)
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
                   init(i:i+1,j,k)=0
                else
                   k1=ceiling(kr)
                   k2=k1+1
                   w1=k1-kr
                   w2=1-w1
                   p0=power(2*pi*kr/box)*Da**2/(4*pi*kr**3)
                   p1=p0
                   p2=p1**2
                   p3=p1*(4*pi*kr**3)
                   p4=p3**2
                   ps0(1,k1,it)=ps0(1,k1,it)+p1*w1
                   ps0(2,k1,it)=ps0(2,k1,it)+p2*w1
                   ps0(3,k1,it)=ps0(3,k1,it)+p3*w1
                   ps0(4,k1,it)=ps0(4,k1,it)+p4*w1
                   ps0(5,k1,it)=ps0(5,k1,it)+w1
                   ps0(1,k2,it)=ps0(1,k2,it)+p1*w2
                   ps0(2,k2,it)=ps0(2,k2,it)+p2*w2
                   ps0(3,k2,it)=ps0(3,k2,it)+p3*w2
                   ps0(4,k2,it)=ps0(4,k2,it)+p4*w2
                   ps0(5,k2,it)=ps0(5,k2,it)+w2
                   p1=p0*sum(noise(i:i+1,j,k)**2)/ncr**3
                   p2=p1**2
                   p3=p1*(4*pi*kr**3)
                   p4=p3**2
                   pst(1,k1,it)=pst(1,k1,it)+p1*w1
                   pst(2,k1,it)=pst(2,k1,it)+p2*w1
                   pst(3,k1,it)=pst(3,k1,it)+p3*w1
                   pst(4,k1,it)=pst(4,k1,it)+p4*w1
                   pst(5,k1,it)=pst(5,k1,it)+w1
                   pst(1,k2,it)=pst(1,k2,it)+p1*w2
                   pst(2,k2,it)=pst(2,k2,it)+p2*w2
                   pst(3,k2,it)=pst(3,k2,it)+p3*w2
                   pst(4,k2,it)=pst(4,k2,it)+p4*w2
                   pst(5,k2,it)=pst(5,k2,it)+w2
                   delta=sqrt(p0*ncr**3)
                   init(i:i+1,j,k)=delta*noise(i:i+1,j,k)
                endif
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Merge power spectrum from threads
    do it=2,nt
       pst(:,:,1)=pst(:,:,1)+pst(:,:,it)
    enddo
    do it=2,nt
       ps0(:,:,1)=ps0(:,:,1)+ps0(:,:,it)
    enddo


    !! Divide by weights
    !! psg(1,k) stores gas p(k)
    !! psg(2,k) stores standard deviation
    do k=1,nc
       if (pst(5,k,1) .eq. 0) then
          psg(:,k)=0
       else
          psg(:,k)=pst(1:4,k,1)/pst(5,k,1)
          psg(2,k)=sqrt(abs(psg(2,k)-psg(1,k)**2))
          psg(1:2,k)=4*pi*(k-1)**3*psg(1:2,k)
          psg(4,k)=sqrt(abs(psg(4,k)-psg(3,k)**2))
       endif
    enddo
    do k=1,nc
       if (ps0(5,k,1) .eq. 0) then
          psinit(:,k)=0
       else
          psinit(:,k)=ps0(1:4,k,1)/ps0(5,k,1)
          psinit(2,k)=sqrt(abs(psinit(2,k)-psinit(1,k)**2))
          psinit(1:2,k)=4*pi*(k-1)**3*psinit(1:2,k)
          psinit(4,k)=sqrt(abs(psinit(4,k)-psinit(3,k)**2))
       endif
    enddo

    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called powerspectrum'
    return
  end subroutine powerspectrum


!!------------------------------------------------------------------!!


  subroutine kernel
    implicit none
    real, parameter :: rmax=hc

    integer i,j,k
    real r,x,y,z

    real time1,time2
    call cpu_time(time1)


    !! Calculate displacement kernels
    !! Kernel is the same as isotropic potential kernel
    !$omp parallel do default(shared) private(i,j,k,r,x,y,z)
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
             if (r .eq. 0) then
                phi(i,j,k)=-2.5
             else
                phi(i,j,k)=-1/min(r,rmax)
             endif
          enddo
          phi(nc+1:,j,k)=0
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT kernels
    call fft3d(phi,nc,'f')


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called kernel'
    return
  end subroutine kernel


  subroutine zeldovich
    implicit none
    integer i,j,k

    real time1,time2
    call cpu_time(time1)


    !! Calculate Fourier space potential
    !! Complex multiply delta field with kernel
    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc
       do j=1,nc
          do i=1,nc+2,2
             phi(i:i+1,j,k)=phi(i,j,k)*init(i:i+1,j,k)/(4*pi)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Inverse Fourier transform potential field
    call fft3d(phi,nc,'b')


    !! Calculate displacement field
    select case (np)
    case (nc)
       call displacement1
    case (hc)
       call displacement2
    end select


    !! Compute dark matter power spectrum
    call dmpowerspectrum


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called zeldovich'
    return
  end subroutine zeldovich


  subroutine displacement1
    implicit none
    real, parameter :: ncr=nc

    integer i,j,k
    integer im,ip,jm,jp,km,kp

    real time1,time2
    call cpu_time(time1)


    !! Finite difference potential to get displacement field
    !! Move particles using displacement field

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
             xp(1,i,j,k)=(phi(im,j,k)-phi(ip,j,k))/2
             xp(2,i,j,k)=(phi(i,jm,k)-phi(i,jp,k))/2
             xp(3,i,j,k)=(phi(i,j,km)-phi(i,j,kp))/2             
             xp(1,i,j,k)=mod(xp(1,i,j,k)+(i-0.5)+ncr,ncr)
             xp(2,i,j,k)=mod(xp(2,i,j,k)+(j-0.5)+ncr,ncr)
             xp(3,i,j,k)=mod(xp(3,i,j,k)+(k-0.5)+ncr,ncr)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called displacement np=nc'
    return
  end subroutine displacement1


  subroutine displacement2
    implicit none
    real, parameter :: ncr=nc

    integer i,j,k
    integer i0,j0,k0,i1,j1,k1
    integer im,ip,jm,jp,km,kp
    real*8 dis(3)

    real time1,time2
    call cpu_time(time1)


    !! Finite difference potential to get displacement field
    !! Average over 8 cells to get coarser displacement field
    !! Move particles using displacement field

    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(i0,j0,k0,i1,j1,k1,im,ip,jm,jp,km,kp,dis)
    do k=1,nc,2
       k0=(k+1)/2
       do j=1,nc,2
          j0=(j+1)/2
          do i=1,nc,2
             i0=(i+1)/2
             dis=0
             do k1=k,k+1
                kp=mod(k1,nc)+1
                km=mod(k1-2+nc,nc)+1
                do j1=j,j+1
                   jp=mod(j1,nc)+1
                   jm=mod(j1-2+nc,nc)+1
                   do i1=i,i+1
                      ip=mod(i1,nc)+1
                      im=mod(i1-2+nc,nc)+1
                      dis(1)=dis(1)+(phi(im,j1,k1)-phi(ip,j1,k1))/2
                      dis(2)=dis(2)+(phi(i1,jm,k1)-phi(i1,jp,k1))/2
                      dis(3)=dis(3)+(phi(i1,j1,km)-phi(i1,j1,kp))/2
                   enddo
                enddo
             enddo
             dis=dis/8
             xp(1,i0,j0,k0)=mod(dis(1)+i+ncr,ncr)
             xp(2,i0,j0,k0)=mod(dis(2)+j+ncr,ncr)
             xp(3,i0,j0,k0)=mod(dis(3)+k+ncr,ncr)
          enddo
       enddo
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called displacement np=hc'
    return
  end subroutine displacement2


  subroutine dmpowerspectrum
    implicit none
    real, parameter :: ncr=nc
    integer, parameter :: kpt=nc/nt
    integer, parameter :: npt=np/nt

    integer it
    integer i,j,k
    integer k1,k2
    real kr,kx,ky,kz,w1,w2
    real p1,p2,p3,p4
    real*8 pst(5,nc,nt)


    !! Mass assignment using cic to density field
    !! Note that the average density is set to be 1
    !! Therefore initial to -1 for delta field

    !$omp parallel do default(shared) private(k)
    do k=1,nc
       den(:,:,k)=-1
    enddo
    !$omp end parallel do


    !$omp parallel do default(shared) private(it,k)
    do it=1,nt
       do k=1+(it-1)*npt,it*npt
          call cicmass(k)
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT density field
    call fft3d(den,nc,'f')


    !! Compute power spectrum
    !$omp parallel do default(shared) private(it,i,j,k) &
    !$omp& private(kr,kx,ky,kz,k1,k2,w1,w2,p1,p2,p3,p4)
    do it=1,nt
       pst(:,:,it)=0
       do k=1+(it-1)*kpt,it*kpt
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
                   p1=sum((den(i:i+1,j,k)/ncr**3)**2)
                   p2=p1**2
                   p3=p1*(4*pi*kr**3)
                   p4=p3**2
                   pst(1,k1,it)=pst(1,k1,it)+p1*w1
                   pst(2,k1,it)=pst(2,k1,it)+p2*w1
                   pst(3,k1,it)=pst(3,k1,it)+p3*w1
                   pst(4,k1,it)=pst(4,k1,it)+p4*w1
                   pst(5,k1,it)=pst(5,k1,it)+w1
                   pst(1,k2,it)=pst(1,k2,it)+p1*w2
                   pst(2,k2,it)=pst(2,k2,it)+p2*w2
                   pst(3,k2,it)=pst(3,k2,it)+p3*w2
                   pst(4,k2,it)=pst(4,k2,it)+p4*w2
                   pst(5,k2,it)=pst(5,k2,it)+w2
                endif
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Merge power spectrum from threads
    do it=2,nt
       pst(:,:,1)=pst(:,:,1)+pst(:,:,it)
    enddo


    !! Divide by weights
    !! psdm(1,k) stores dm p(k)
    !! psdm(2,k) stores standard deviation
    do k=1,nc
       if (pst(5,k,1) .eq. 0) then
          psdm(:,k)=0
       else
          psdm(:,k)=pst(1:4,k,1)/pst(5,k,1)
          psdm(2,k)=sqrt(psdm(2,k)-psdm(1,k)**2)
          psdm(1:2,k)=4*pi*(k-1)**3*psdm(1:2,k)
          psdm(4,k)=sqrt(psdm(4,k)-psdm(3,k)**2)
       endif
    enddo


    return
  end subroutine dmpowerspectrum


  subroutine cicmass(k)
    implicit none
    real, parameter :: ncr=nc
    real, parameter :: mp=(ncr/np)**3
    integer k

    integer i,j
    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    do j=1,np
       do i=1,np
          x=mod(xp(1,i,j,k)-0.5+ncr,ncr)
          y=mod(xp(2,i,j,k)-0.5+ncr,ncr)
          z=mod(xp(3,i,j,k)-0.5+ncr,ncr)

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
          den(i1,j1,k1)=den(i1,j1,k1)+dx1*dy1*dz1
          den(i2,j1,k1)=den(i2,j1,k1)+dx2*dy1*dz1
          den(i1,j2,k1)=den(i1,j2,k1)+dx1*dy2*dz1
          den(i2,j2,k1)=den(i2,j2,k1)+dx2*dy2*dz1
          den(i1,j1,k2)=den(i1,j1,k2)+dx1*dy1*dz2
          den(i2,j1,k2)=den(i2,j1,k2)+dx2*dy1*dz2
          den(i1,j2,k2)=den(i1,j2,k2)+dx1*dy2*dz2
          den(i2,j2,k2)=den(i2,j2,k2)+dx2*dy2*dz2
       enddo
    enddo

    return
  end subroutine cicmass


!!------------------------------------------------------------------!!


  subroutine densitystatistics
    implicit none
    integer, parameter :: kpt=nc/nt

    integer i,j,k
    real d,dmin,dmax
    real*8 dsum,dvar

    real time1,time2
    call cpu_time(time1)


    !! Inverse Fourier transform to get real space gas density field
    call fft3d(init,nc,'b')


    dmin=0
    dmax=0
    dsum=0
    dvar=0
    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) reduction(max:dmax) reduction(+:dsum,dvar)
    do k=1,nc
       do j=1,nc
          do i=1,nc
             d=init(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo
    !$omp end parallel do

    dvar=sqrt(dvar/nc**3)
    write(*,*)
    write(*,*) 'Gas min  ',dmin
    write(*,*) 'Gas max  ',dmax
    write(*,*) 'Deltasum ',real(dsum)
    write(*,*) 'Deltavar ',real(dvar)
    write(*,*)


    !! Inverse Fourier transform to get real space dm density field
    call fft3d(den,nc,'b')

    dmin=0
    dmax=0
    dsum=0
    dvar=0
    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) reduction(max:dmax) reduction(+:dsum,dvar)
    do k=1,nc
       do j=1,nc
          do i=1,nc
             d=den(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo
    !$omp end parallel do

    dvar=sqrt(dvar/nc**3)
    write(*,*)
    write(*,*) 'DM min   ',dmin
    write(*,*) 'DM max   ',dmax
    write(*,*) 'Deltasum ',real(dsum)
    write(*,*) 'Deltavar ',real(dvar)
    write(*,*)


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called density statistics'
    return
  end subroutine densitystatistics


!!------------------------------------------------------------------!!


  function power(kr)
    implicit none
    real power,kr

    integer i,i1,i2
    real x,y,x1,x2,y1,y2

    i1=1
    i2=nk
    do
       if (i2-i1 .eq. 1) then
          exit
       else
          i=(i1+i2)/2
          if (kr .gt. tf(1,i)) then
             i1=i
          else
             i2=i
          endif
       endif
    enddo

    x1=log(tf(1,i1))
    y1=log(tf(3,i1))
    x2=log(tf(1,i2))
    y2=log(tf(3,i2))
    x=log(kr)
    y=y1+(y2-y1)*(x-x1)/(x2-x1)
    power=exp(y)

    return
  end function power


  function grow(a)
    implicit none
    real grow,a

    real hsq,oma,ola

    hsq=omegam/a**3+(1-omegam-omegal)/a**2+omegal
    oma=omegam/(a**3*hsq)
    ola=omegal/hsq
    grow=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    return
  end function grow


  function tophat(x)
    implicit none
    real tophat,x

    if (x .ne. 0) then
       tophat=3*(sin(x)-cos(x)*x)/x**3
    else
       tophat=1
    endif

    return
  end function tophat


  subroutine fft3d(a,n,c)
    implicit none
#ifndef CXML
    include '/opt/fftw-3.0.1/include/fftw3.f'
#endif
    integer n
    character c
    real a(n+2,n,n)
#ifdef CXML
    external sfft_3d
    
    if (c .eq. 'f') then
       call sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
    else
       call sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
    endif
#else
    integer(8) :: plan,iplan 
    if (c .eq. 'f') then
      call sfftw_plan_dft_r2c_3d(plan,n,n,n,a,a,FFTW_ESTIMATE)
      call sfftw_execute(plan)
      call sfftw_destroy_plan(plan)
    else 
      call sfftw_plan_dft_c2r_3d(iplan,n,n,n,a,a,FFTW_ESTIMATE)
      call sfftw_execute(iplan)
      a = a / real(n)**3
      call sfftw_destroy_plan(iplan)
    endif
#endif

    return
  end subroutine fft3d


end program main
