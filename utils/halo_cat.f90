! halo_cat.f90 :: construct halo catalog from recomposed pm data
! to compile :: f90 -lcxmlp halo_cat.f90 indexedsort.f90 -o halo_cat.x
 
implicit none

   real(4), parameter :: phi_threshold = -500.0
   integer(4), parameter :: min_cell_search = 18
   real(4), parameter :: overdensity = 178.0

   integer(4), parameter :: nc=320

   character(len=*),parameter :: fpath='./'
   character(len=*),parameter :: checkpoints='./checkpoints'

   character(len=80) :: ifile,catalog
   character(len=7) :: z_s
   integer(4), parameter :: max_input=100
   integer(4), parameter :: hc=nc/2
   integer(4), parameter :: np_max=hc**3
   real(4), parameter :: ncr=nc
   real(4), parameter :: hcr=hc
   real(4), parameter :: pi=3.141592654

   real(kind=4), dimension(6,np_max) :: xv
   integer(4), dimension(np_max) :: ll
   integer(4), dimension(nc,nc) :: hoc
   real(4), dimension(nc+2,nc,nc) :: rho,phi
   real(4), dimension(hc+1,hc+1,hc+1) :: wgr

   integer(4), parameter :: nsbox=32
   integer(4), parameter :: nlist=(nsbox+1)**3*5
   integer(4), parameter :: max_num_peaks=nsbox**3

   real(4) :: phimin, amass, amtot
   real(4), dimension(max_num_peaks) :: peak
   integer(4), dimension(3,max_num_peaks) :: ipeak
   integer(4), dimension(max_num_peaks) :: isortpeak
   integer(4), dimension(nlist) :: isortdist
   integer(4), dimension(3,nlist) :: idist
   real(4), dimension(nlist) :: rdist

   integer(kind=4) :: ii,kx,ky,kz,i,j,k,cp,nploc,fstat,num_checkpoints
   integer(kind=4) :: i1,j1,k1,ix0,iy0,iz0,ix,iy,iz
   integer(kind=4) :: im,ip,jm,jp,km,kp,num_halos,irtot,iloc
   real(kind=4) :: r,x,y,z,ka,kb,kc,kr
   real(kind=4), dimension(max_input) :: z_checkpoint

!! initialize halo finder 

    ii=0
    do i=-nsbox,nsbox
      do j=-nsbox,nsbox
        do k=-nsbox,nsbox
          r=sqrt(i**2+j**2+k**2*1.)
          if (r>nsbox) cycle
          ii=ii+1
          if (ii>nlist) then
            write(*,*) 'ii exceeded ',nlist
            pause
          endif
          idist(:,ii)=(/i,j,k/)
          rdist(ii)=r
        enddo
      enddo
    enddo
    irtot=ii
 
    isortdist(:ii)=(/ (i,i=1,ii) /)
    call indexedsort(ii,rdist,isortdist)
    idist(:,:ii)=idist(:,isortdist(:ii))

!! construct gravitational kernel

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
             r=sqrt(x**2+y**2+z**2*1.)
             if (r .eq. 0) then
                phi(i,j,k)=-2.5
             else
                phi(i,j,k)=-1/r
             endif
          enddo
          phi(nc+1:,j,k)=0
       enddo
    enddo

    call fft3d(phi,nc,'f')

    do k=1,hc+1
       kz=k-1
       do j=1,hc+1
          ky=j-1
          do i=1,hc+1
             kx=i-1
             kr=sqrt(1.*kx**2+ky**2+kz**2)
             if (kr .eq. 0 .or. kr .gt. 8) then
                wgr(i,j,k)=phi(2*i-1,j,k)
             else
                ka=2*sin(pi*kx/ncr)
                kb=2*sin(pi*ky/ncr)
                kc=2*sin(pi*kz/ncr)
                wgr(i,j,k)=-4*pi/(ka**2+kb**2+kc**2)
             endif
          enddo
       enddo
    enddo


!! Read in checkpoints to generate halos for 

  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) write(*,*) &
     'error opening checkpoint list file:',checkpoints
  do num_checkpoints=1,max_input
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
  enddo
41  num_checkpoints=num_checkpoints-1
51  close(11)
  write(*,*) 'checkpoints to construct halos for:'
  do i=1,num_checkpoints
    write(*,'(f5.1)') z_checkpoint(i)
  enddo

  do cp=1,num_checkpoints

   write(z_s,'(f7.3)') z_checkpoint(cp)
   z_s=adjustl(z_s)

   ifile=fpath//z_s(1:len_trim(z_s))//'xv.dat'
   catalog=fpath//z_s(1:len_trim(z_s))//'halo.dat'

   open (unit=13,file=ifile,form='binary')
   read(13) nploc
   if (nploc > np_max) then
     write(*,*) 'too many particles to store',nploc,np_max
     stop
   endif
   read(13) xv(:,:nploc)
   close(13)

!! build linked list
   hoc=0
   ll=0
   do ip=1,nploc
     i=int(xv(1,ip))+1
     j=int(xv(2,ip))+1
     ll(ip)=hoc(i,j)
     hoc(i,j)=ip
   enddo 

!! build rho
  
   rho=0.0
   
   do k=1,nc
     do j=1,nc
       ip=hoc(j,k)
       call cicmass(ip)
     enddo
   enddo  

!! build potential

   call fft3d(rho,nc,'f')

   do k=1,nc
     if (k .lt. hc+2) then
       k1=k
     else
       k1=nc+2-k
     endif
     do j=1,nc
        if (j .lt. hc+2) then
          j1=j
        else
          j1=nc+2-j
        endif
        do i=1,nc+2,2
          i1=(i-1)/2+1
          phi(i:i+1,j,k)=rho(i:i+1,j,k)*wgr(i1,j1,k1)
        enddo
     enddo
   enddo

   call fft3d(phi,nc,'b')

! find potential minima 

   num_halos = 0
   do k=1,nc
     kp=mod(k,nc)+1
     km=mod(k-2+nc,nc)+1
     do j=1,nc
       jp=mod(j,nc)+1
       jm=mod(j-2+nc,nc)+1
       do i=1,nc
         ip=mod(i,nc)+1
         im=mod(i-2+nc,nc)+1
         if (i > 1 .and. i < nc .and. j > 1 .and. j < nc &
             .and. k > 1 .and. k < nc) then
           phimin=minval(phi(i-1:i+1,j-1:j+1,k-1:k+1))
         else
           phimin=min(phi(im,jm,km),phi(im,jm,k),phi(im,jm,kp), &
                      phi(im,j,km), phi(im,j,k), phi(im,j,kp), &
                      phi(im,jp,km),phi(im,jp,k),phi(im,jp,kp), &
                      phi(i,jm,km), phi(i,jm,k), phi(i,jm,kp), &
                      phi(i,j,km),  phi(i,j,k),  phi(i,j,kp), &
                      phi(i,jp,km), phi(i,jp,k), phi(i,jp,kp),&
                      phi(ip,jm,km),phi(ip,jm,k),phi(ip,jm,kp), &
                      phi(ip,j,km), phi(ip,j,k), phi(ip,j,kp), &
                      phi(ip,jp,km),phi(ip,jp,k),phi(ip,jp,kp))
         endif
         if (phimin .eq. phi(i,j,k) .and. phimin < phi_threshold) then
            if (num_halos > max_num_peaks - 2) then    
              write(*,*) 'too many halos',num_halos,max_num_peaks-2
              stop 
            endif
            num_halos=num_halos+1
            ipeak(:,num_halos)=(/i,j,k/)
            peak(num_halos)=phimin
         endif
       enddo
     enddo
   enddo

   print *,num_halos,' local gravitational minima'

   isortpeak(:num_halos)=(/ (i,i=1,num_halos) /)
   call indexedsort(num_halos,peak,isortpeak)

   ipeak(:,:num_halos)=ipeak(:,isortpeak(:num_halos))

   call fft3d(rho,nc,'b')

   open(66,file=catalog,form='formatted')

   do iloc=1,num_halos
     ix0=ipeak(1,iloc)
     iy0=ipeak(2,iloc)
     iz0=ipeak(3,iloc)
     amtot=0
     do i=1,irtot
       ix=modulo(ix0+idist(1,i)-1,nc)+1
       iy=modulo(iy0+idist(2,i)-1,nc)+1
       iz=modulo(iz0+idist(3,i)-1,nc)+1
       amass=rho(ix,iy,iz)
       rho(ix,iy,iz)=0
       amtot=amtot+amass
       if (i > min_cell_search .and. amtot/(real(i)) < overdensity) exit
     enddo
     if (i>irtot) then
       i=irtot
       write(*,*) 'ran out of irtot, halo:',iloc
     endif
     write(66,'(3i6,3f12.1)') ix0,iy0,iz0,amtot,peak(iloc),rdist(i) 
   enddo

   close(66)

   write(*,*) 'done z=',z_checkpoint(cp)
 
  enddo

contains

  subroutine cicmass(ip)
    implicit none
    real ::  mp
    integer ip

    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    mp=ncr**3/nploc

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
       rho(i1,j1,k1)=rho(i1,j1,k1)+dx1*dy1*dz1
       rho(i2,j1,k1)=rho(i2,j1,k1)+dx2*dy1*dz1
       rho(i1,j2,k1)=rho(i1,j2,k1)+dx1*dy2*dz1
       rho(i2,j2,k1)=rho(i2,j2,k1)+dx2*dy2*dz1
       rho(i1,j1,k2)=rho(i1,j1,k2)+dx1*dy1*dz2
       rho(i2,j1,k2)=rho(i2,j1,k2)+dx2*dy1*dz2
       rho(i1,j2,k2)=rho(i1,j2,k2)+dx1*dy2*dz2
       rho(i2,j2,k2)=rho(i2,j2,k2)+dx2*dy2*dz2

       ip=ll(ip)
    enddo

    return
  end subroutine cicmass

end

  subroutine fft3d(a,n,c)
    implicit none
    integer n
    character c
    real a(n+2,n,n)
    external sfft_3d

    if (c .eq. 'f') then
       call sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
    else
       call sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
    endif

    return
  end subroutine fft3d
