
  subroutine powerspectrum
    implicit none

    include 'mpif.h'
    include 'dist_init.fh'

    real, parameter :: ncr=nc
    integer, parameter :: kpt=nc_slab/nt

    integer it,i,j,k
    integer k1,k2
    real kx,ky,kz,kr,w1,w2,kzg
    real Da,delta,p0,p1,p2,p3,p4
    real*8 ps0(5,nc,nt),pst(5,nc,nt)

    real time1,time2
    call cpu_time(time1)

    !! Calculate growth factor at scalefactor a
    Da=scalefactor*grow(scalefactor)/grow(1.)

    !! Interpolate \Delta^2 to grid
    !! Determine Fourier modes \delta(k)
    !! Compute both sampled and unsampled power spectra

!! work on this in slab coords

    !$omp parallel do default(shared) private(it,i,j,k) &
    !$omp& private(kr,kzg,kx,ky,kz,k1,k2,w1,w2,delta,p0,p1,p2,p3,p4)
    do it=1,nt
       pst(:,:,it)=0
       ps0(:,:,it)=0
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          kgz = k+nc_slab*rank
          if (kgz .lt. hc+2) then
             kz=kgz-1
          else
             kz=kgz-1-nc
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
   		   p1=p0*sum(slab(i:i+1,j,k)**2)/ncr**3
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
   		   init(i:i+1,j,k)=delta*slab(i:i+1,j,k)
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

    !! Merge to master node, scatter
    call mpi_reduce(pst(:,:,1),ps_sum,5*nc,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
  
    if (rank == 0) then

    !! Divide by weights
    !! psg(1,k) stores gas p(k)
    !! psg(2,k) stores standard deviation

      do k=1,nc
         if (ps_sum(5,k) .eq. 0) then
            psg(:,k)=0
         else
            psg(:,k)=ps_sum(1:4,k)/ps_sum(5,k)
            psg(2,k)=sqrt(abs(psg(2,k)-psg(1,k)**2))
            psg(1:2,k)=4*pi*(k-1)**3*psg(1:2,k)
            psg(4,k)=sqrt(abs(psg(4,k)-psg(3,k)**2))
         endif
      enddo

    endif

    call mpi_bcast(psg,4*nc,mpi_double_precision,0,mpi_comm_world,ierr)

    call mpi_reduce(ps0(:,:,1),ps_sum,5*nc,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

    if (rank == 0) then

      do k=1,nc
         if (ps_sum(5,k) .eq. 0) then
            psinit(:,k)=0
         else
            psinit(:,k)=ps_sum(1:4,k)/ps_sum(5,k)
            psinit(2,k)=sqrt(abs(psinit(2,k)-psinit(1,k)**2))
            psinit(1:2,k)=4*pi*(k-1)**3*psinit(1:2,k)
            psinit(4,k)=sqrt(abs(psinit(4,k)-psinit(3,k)**2))
         endif
      enddo

    endif

    call mpi_bcast(psinit,4*nc,mpi_double_precision,0,mpi_comm_world,ierr)

    call cpu_time(time2)
    time2=(time2-time1)/nt
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called powerspectrum'
    return
  end subroutine powerspectrum

