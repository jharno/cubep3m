 subroutine checknoise
    implicit none
 
    include 'mpif.h'
    include 'dist_init.fh'

    integer, parameter :: kc=hc
    integer, parameter :: kpt=nc_slab/nt
    real, parameter :: norm=nc**3

    integer it,i,j,k
    integer ksq,kx,ky,kz,kg
    integer k1,k2
    real ps,ws
    real*8 A(2,(kc+1)**2,nt)

    real time1,time2
    call cpu_time(time1)

    !$omp parallel do default(shared) private(it,i,j,k,kg) &
    !$omp& private(ksq,kx,ky,kz)
    do it=1,nt
       A(:,:,it)=0
       do k=1+(it-1)*kpt,it*kpt
          kg=k+nc_slab*rank
          if (kg .lt. hc+2) then
             kz=kg-1
          else
             kz=kg-1-nc
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
                   A(1,ksq,it)=A(1,ksq,it)+sum(slab(i:i+1,j,k)**2)
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

    !! Merge to master node, scatter
    call mpi_reduce(A(:,:,1),Asum,2*(hc+1)**2,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

    !! Compute normalization
    if (rank == 0) then
      do i=1,kc
         k1=i**2
         k2=(i+1)**2-1
         ps=sum(Asum(1,k1:k2))
         ws=sum(Asum(2,k1:k2))
         if (ws .ne. 0) Asum(1,k1:k2)=sqrt(ps/ws/norm)
         write(*,*) k1,k2,sqrt(1.*k1),sqrt(1.*k2),Asum(1,k1)
      enddo
    endif

    call mpi_bcast(Asum,2*(hc+1)**2,mpi_double_precision,0,mpi_comm_world,ierr)

    !$omp parallel do default(shared) private(i,j,k,kg) &
    !$omp& private(ksq,kx,ky,kz)
    do k=1,nc_slab
       kg=k+nc_slab*rank
       if (kg .lt. hc+2) then
          kz=kg-1
       else
          kz=kg-1-nc
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
                slab(i:i+1,j,k)=slab(i:i+1,j,k)/Asum(1,ksq)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do

    call cpu_time(time2)
    time2=(time2-time1)/nt
    if (rank==0) write(*,"(f8.2,a)") time2,'  Called check noise'
    return
  end subroutine checknoise

