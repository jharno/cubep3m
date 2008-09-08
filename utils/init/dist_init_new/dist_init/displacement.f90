!! potential is in cube

#ifdef FIXED_DIS
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
    do k=1,nc_node_dim
       kp=k+1
       km=k-1
       do j=1,nc_node_dim
          jp=j+1
          jm=j-1
          do i=1,nc_node_dim
             ip=i+1
             im=i-1
             xp(1,i,j,k)=(phi(im,j,k)-phi(ip,j,k))/2
             xp(2,i,j,k)=(phi(i,jm,k)-phi(i,jp,k))/2
             xp(3,i,j,k)=(phi(i,j,km)-phi(i,j,kp))/2
             xp(1,i,j,k)=mod(xp(1,i,j,k)+(i-0.5)) !+ncr,ncr)
             xp(2,i,j,k)=mod(xp(2,i,j,k)+(j-0.5)) !+ncr,ncr)
             xp(3,i,j,k)=mod(xp(3,i,j,k)+(k-0.5)) !+ncr,ncr)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called displacement np=nc'
    return
  end subroutine displacement1
#endif

!! This version is for np=hc -- currently only one used
  subroutine displacement2
    implicit none

    include 'mpif.h'
    include 'dist_init.fh'

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
    phi=0.0
    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube

    if (rank == 7) write(*,*) '7s',cube(nc_node_dim,nc_node_dim,nc_node_dim), &
                                    phi(nc_node_dim,nc_node_dim,nc_node_dim)
    !! Populate buffer (1 cell thick):

    call phi_buffer(0)

#ifdef DEBUG  
    if (rank == 0) write(*,*) 'phi after buffer' 
    if (rank == 0) write(*,*) phi(0:3,0:3,0:3)
    call mpi_barrier(mpi_comm_world,ierr)
    if (rank == 7) write(*,*) 'sender'
    if (rank == 7) write(*,*) phi(nc_node_dim,nc_node_dim,nc_node_dim)
#endif
 
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
             if (dis(1) > 3.0) then
               write(*,*) 'dis out of bounds',rank,i0,j0,k0,dis(1),xv_grid(1,i0,j0,k0)
               call mpi_abort(mpi_comm_world,ierr,ierr)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do

    call cpu_time(time2)
    time2=(time2-time1)/nt
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called displacement np=hc'
    return
  end subroutine displacement2

