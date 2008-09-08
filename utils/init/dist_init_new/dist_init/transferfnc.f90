!! calculates transfer function
  subroutine transferfnc
    implicit none
    include 'dist_init.fh'

    integer i
    real kmax,norm,tophat
    real*8 v8

    real time1,time2
    external tophat

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
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called transfer fnc'
    return
  end subroutine transferfnc
