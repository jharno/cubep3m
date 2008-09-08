  subroutine dmpowerspectrum
    implicit none

    include 'mpif.h'
    include 'dist_init.fh'

    real, parameter :: ncr=nc
    integer, parameter :: kpt=nc_node_dim/nt
    integer, parameter :: npt=np_node_dim/nt

    integer it
    integer i,j,k
    integer k1,k2
    real kr,kx,ky,kz,w1,w2,kzg
    real p1,p2,p3,p4
    real*8 pst(5,nc,nt)

    !! Mass assignment using cic to density field
    !! Note that the average density is set to be 1
    !! Therefore initial to -1 for delta field

    phi(0,:,:)=0.0
    phi(nc_node_dim+1,:,:)=0.0 
    phi(:,0,:)=0.0 
    phi(:,nc_node_dim+1,:)=0.0 
    phi(:,:,0)=0.0 
    phi(:,:,nc_node_dim+1)=0.0 

    !$omp parallel do default(shared) private(k)
    do k=1,nc_node_dim
       phi(1:nc_node_dim,1:nc_node_dim,k)=-1
    enddo
    !$omp end parallel do

    !$omp parallel do default(shared) private(it,k)
    do it=1,nt
       do k=1+(it-1)*npt,it*npt
          call cicmass(k)
       enddo
    enddo
    !$omp end parallel do

    !! pass buffers to adjacent nodes and accumulate.
    !! copy to cube for processing

    call phi_buffer(1)

    cube=phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)

    !! Forward FFT density field

    call fftw(1)

    !! Compute power spectrum
    !$omp parallel do default(shared) private(it,i,j,k) &
    !$omp& private(kr,kzg,kx,ky,kz,k1,k2,w1,w2,p1,p2,p3,p4)
    do it=1,nt
       pst(:,:,it)=0
       do k=1+(it-1)*kpt,it*kpt
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
                if (kr .ne. 0) then
                   k1=ceiling(kr)
                   k2=k1+1
                   w1=k1-kr
                   w2=1-w1
                   p1=sum((slab(i:i+1,j,k)/ncr**3)**2)
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

    ps_sum=0.0
    !! Merge to master node, scatter
    call mpi_reduce(pst(:,:,1),ps_sum,5*nc,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

    if (rank == 0) then
      do k=1,nc
        if (ps_sum(5,k) .eq. 0) then
          psdm(:,k)=0
        else
          psdm(:,k)=ps_sum(1:4,k)/ps_sum(5,k)
          psdm(2,k)=sqrt(psdm(2,k)-psdm(1,k)**2)
          psdm(1:2,k)=4*pi*(k-1)**3*psdm(1:2,k)
          psdm(4,k)=sqrt(psdm(4,k)-psdm(3,k)**2)
        endif
      enddo
    endif

    call mpi_bcast(psdm,4*nc,mpi_double_precision,0,mpi_comm_world,ierr)

    return
  end subroutine dmpowerspectrum


  subroutine cicmass(k)
    implicit none

    include 'dist_init.fh'

    real, parameter :: ncr=nc
    real, parameter :: mp=(ncr/np)**3
    integer k

    integer i,j
    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    do j=1,np_node_dim
       do i=1,np_node_dim
          x=xv_grid(1,i,j,k)-0.5
          y=xv_grid(2,i,j,k)-0.5
          z=xv_grid(3,i,j,k)-0.5

          i1=floor(x)+1
          i2=i1+1
          dx1=i1-x
          dx2=1-dx1
          j1=floor(y)+1
          j2=j1+1
          dy1=j1-y
          dy2=1-dy1
          k1=floor(z)+1
          k2=k1+1
          dz1=k1-z
          dz2=1-dz1

          dz1=mp*dz1
          dz2=mp*dz2
          phi(i1,j1,k1)=phi(i1,j1,k1)+dx1*dy1*dz1
          phi(i2,j1,k1)=phi(i2,j1,k1)+dx2*dy1*dz1
          phi(i1,j2,k1)=phi(i1,j2,k1)+dx1*dy2*dz1
          phi(i2,j2,k1)=phi(i2,j2,k1)+dx2*dy2*dz1
          phi(i1,j1,k2)=phi(i1,j1,k2)+dx1*dy1*dz2
          phi(i2,j1,k2)=phi(i2,j1,k2)+dx2*dy1*dz2
          phi(i1,j2,k2)=phi(i1,j2,k2)+dx1*dy2*dz2
          phi(i2,j2,k2)=phi(i2,j2,k2)+dx2*dy2*dz2
       enddo
    enddo

    return
  end subroutine cicmass


