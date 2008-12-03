!! cicpow.f90 written by Hy Trac on Feb 23, 2003
!! Compile with f90 -fast -omp -lcxmlp -o cicpow.x cicpow.f90
!! Modified to read in cubepm checkpoints -- Hugh Merz 041119

program main
  implicit none

  !! Cosmological parameters
  real, parameter :: box=100


  !! nt is the number of threads
  !! dir is the directory for data files
  integer, parameter :: nt=4
  integer, parameter :: nout=9
  character(*), parameter :: dir=''
  real(4), dimension(nout) :: zout

  !! nc is the number of cells per length
  integer, parameter :: nc=320
  integer, parameter :: hc=nc/2


  !! np should be set to nc (1:1) or hc (1:2)
  integer :: np
  integer,parameter :: np_max=hc**3
  integer :: rs

  !! Dark matter arrays
  real, dimension(6,np_max) :: xv
  integer, dimension(np_max) :: ll
  integer, dimension(nc,nc) :: hoc
  integer, dimension(2,nc,nc,nt) :: htoc


  !! Power spectrum arrays
  real, dimension(2,nc) :: ps
  real, dimension(nc+2,nc,nc) :: d

  common d,xv,ll,htoc,hoc,ps

  call set_rs
  do rs=1,nout
    call readdm
    call cic
    call powerspectrum
    call writeps
  enddo

contains

  subroutine set_rs
    implicit none
    character*60 fn

    fn=dir//'checkpoints'
    write(*,*) 'Reading ',fn
    open(11,file=fn)
    read(11,*) zout
    close(11) 
  end subroutine set_rs

  subroutine readdm
    implicit none
    character*60 fn
    character*7 r_s 

    real time1,time2
    call cpu_time(time1)

    write(r_s,'(f7.3)') zout(rs)
    r_s=adjustl(r_s)
    !! Read particle data file
    fn=r_s(1:len_trim(r_s))//'xv.dat'
    write(*,*) 'Reading ',fn

    open(11,file=fn,form='binary')
    read(11) np
    if (np > np_max) then
      write(*,*) 'np > np_max',np,np_max
      stop
    endif
    read(11) xv(:,:np)
    close(11)


    call cpu_time(time2)
    time2=time2-time1
    write(*,*) time2, '  Called read dm'
    return
  end subroutine readdm


  subroutine writeps
    implicit none
    integer k
    character*60 fn
    character*7 r_s

    real time1,time2
    call cpu_time(time1)


    !! Output power spectrum
    !! First column is k
    !! Second column is \Delta^2
   
    write(r_s,'(f7.3)') zout(rs)
    fn=r_s(1:len_trim(r_s))//'ps.dm.dat'
    open(11,file=fn)
    do k=2,hc+1
       write(11,*) ps(:,k)
    enddo
    close(11)


    call cpu_time(time2)
    time2=time2-time1
    write(*,*) time2,'  Called write ps'
    return
  end subroutine writeps


  subroutine cic
    implicit none
    integer, parameter :: kpt=nc/nt+min(1,mod(nc,nt))
    integer :: npt

    integer it,i,j,k,ip,toe

    real time1,time2
    call cpu_time(time1)

    npt=np/nt+min(1,mod(np,nt))

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,it*npt)
          j=floor(xv(2,ip))+1
          k=floor(xv(3,ip))+1
          ll(ip)=0
          if (htoc(1,j,k,it) .eq. 0) then
             htoc(:,j,k,it)=ip
          else
             ll(htoc(2,j,k,it))=ip
             htoc(2,j,k,it)=ip
          endif
       enddo
    enddo
    !$omp end parallel do


    !! Merge chaining lists
    !$omp parallel do default(shared) private(it,j,k,toe)
    do k=1,nc
       do j=1,nc
          hoc(j,k)=0
          do it=1,nt
             if (hoc(j,k) .eq. 0) then
                hoc(j,k)=htoc(1,j,k,it)
                toe=htoc(2,j,k,it)
             else
                if (htoc(1,j,k,it) .ne. 0) then
                   ll(toe)=htoc(1,j,k,it)
                   toe=htoc(2,j,k,it)
                endif
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Initialize density field
    !$omp parallel do default(shared) private(it,k)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          d(:,:,k)=0
       enddo
    enddo
    !$omp end parallel do


    !! Add particle density to density field
    !$omp parallel do default(shared) private(it,j,k,ip)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          do j=1,nc
             ip=hoc(j,k)
             call cicmass(ip)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Convert density to overdensity
    !$omp parallel do default(shared) private(it,i,j,k)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          do j=1,nc
             do i=1,nc
                d(i,j,k)=d(i,j,k)-1
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,*) time2,'  Called cic'
    return
  end subroutine cic


  subroutine cicmass(ip)
    implicit none
    real, parameter :: ncr=nc
    real :: mp

    integer ip

    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    mp=ncr**3/np

    do
       if (ip .eq. 0) exit

       x=mod(xv(1,ip)-0.5+ncr,ncr)
       y=mod(xv(2,ip)-0.5+ncr,ncr)
       z=mod(xv(3,ip)-0.5+ncr,ncr)

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

       ip=ll(ip)
    enddo

    return
  end subroutine cicmass


!!--------------------------------------------------------------!!


  subroutine powerspectrum
    implicit none
    integer, parameter :: kpt=nc/nt+min(1,mod(nc,nt))

    integer it,i,j,k
    real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi
    real pst(2,nc,nt)

    real time1,time2
    call cpu_time(time1)


    !! Forward FFT density field
    call fft3d(d,nc,'f')


    !! Compute power spectrum
    
    !$omp parallel do default(shared) &
    !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
    do it=1,nt
       pst(:,:,it)=0
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
                if (kr .ne. 0) then
                   k1=ceiling(kr)
                   k2=k1+1
                   w1=k1-kr
                   w2=1-w1
                   pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                   pst(1,k1,it)=pst(1,k1,it)+w1
                   pst(2,k1,it)=pst(2,k1,it)+w1*pow
                   pst(1,k2,it)=pst(1,k2,it)+w2
                   pst(2,k2,it)=pst(2,k2,it)+w2*pow
                endif
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

    !! Merge power spectrum from threads
    ps=0
    do it=1,nt
       ps=ps+pst(:,:,it)
    enddo


    !! Divide by weights
    !! Convert P(k) to \Delta^2(k)
    !! Store k in ps(1,i)
    
    pi=acos(-1.)
    do k=1,nc
       if (ps(1,k) .ne. 0) then
          ps(2,k)=4*pi*((k-1)**3)*ps(2,k)/ps(1,k)
          ps(1,k)=2*pi*(k-1)/box
       endif
    enddo


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,*) time2,'  Called power spectrum'
    return
  end subroutine powerspectrum


!!--------------------------------------------------------------!!


  subroutine fft3d(a,n,c)
    implicit none
    integer n
    character c
    real a(n+2,n,n)
    external sfft
    
    if (c .eq. 'f') then
       call sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
    else
       call sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
    endif
    return
  end subroutine fft3d


end program main
