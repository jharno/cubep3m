  subroutine kernel
    implicit none
 
    include 'dist_init.fh' 

    real, parameter :: rmax=hc

    integer i,j,k,ig,jg,kg
    real r,x,y,z

    real time1,time2
    call cpu_time(time1)

    !! Calculate displacement kernels
    !! Kernel is the same as isotropic potential kernel
    !$omp parallel do default(shared) private(i,j,k,ig,jg,kg,r,x,y,z)
    do k=1,nc_node_dim
       kg=k+nc_node_dim*cart_coords(1)
       if (kg .lt. hc+2) then
          z=kg-1
       else
          z=kg-1-nc
       endif
       do j=1,nc_node_dim
          jg=j+nc_node_dim*cart_coords(2)
          if (jg .lt. hc+2) then
             y=jg-1
          else
             y=jg-1-nc
          endif
          do i=1,nc_node_dim
             ig=i+nc_node_dim*cart_coords(3)
             if (ig .lt. hc+2) then
                x=ig-1
             else
                x=ig-1-nc
             endif
             r=sqrt(x**2+y**2+z**2)
             if (r .eq. 0) then
                cube(i,j,k)=-2.5
             else
                cube(i,j,k)=-1/min(r,rmax)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT kernels
    call fftw(1)

    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called kernel'
    return
  end subroutine kernel

