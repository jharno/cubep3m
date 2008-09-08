!! potential is in cube

!! This version is for np=nc -- currently not used
!! This needs to be finished
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

!! This version is for np=hc -- currently only one used
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

    !! Have to buffer this array by passing with adjacent cells, 6 passes.

    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube

    !! Populate buffer (1 cell thick):

    call phi_buffer(0)
    
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(i0,j0,k0,i1,j1,k1,im,ip,jm,jp,km,kp,dis)
    do k=1,nc_node_dim,2
       k0=(k+1)/2
       do j=1,nc_node_dim,2
          j0=(j+1)/2
          do i=1,nc_node_dim,2
             i0=(i+1)/2
             dis=0
             do k1=k,k+1
                kp=k1+1
                km=k1-1
                do j1=j,j+1
                   jp=j1+1
                   jm=j1-1
                   do i1=i,i+1
                      ip=i1+1
                      im=i1-1
                      dis(1)=dis(1)+(phi(im,j1,k1)-phi(ip,j1,k1))/2
                      dis(2)=dis(2)+(phi(i1,jm,k1)-phi(i1,jp,k1))/2
                      dis(3)=dis(3)+(phi(i1,j1,km)-phi(i1,j1,kp))/2
                   enddo
                enddo
             enddo
             dis=dis/8
             xv_grid(1,i0,j0,k0)=dis(1)+i
             xv_grid(2,i0,j0,k0)=dis(2)+j
             xv_grid(3,i0,j0,k0)=dis(3)+k
          enddo
       enddo
    enddo
    !$omp end parallel do

    call cpu_time(time2)
    time2=(time2-time1)/nt
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called displacement np=hc'
    return
  end subroutine displacement2

