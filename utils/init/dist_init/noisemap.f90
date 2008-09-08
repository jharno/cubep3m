
  subroutine noisemap
    implicit none

    include 'dist_init.fh'

    integer(4), parameter :: kpt=nc_node_dim/nt

    integer(4) it
    integer(4) i,j,k,one
    real(4) x,x1,x2
    integer(4) iseed(1,nt)

    real time1,time2
    call cpu_time(time1)


    one=1
    !! Generate random reals between 0 and 1
       !$omp parallel do default(shared) private(it,k,x)
       do it=1,nt
          call cpu_time(x)
          iseed(:,it)=it*ceiling(1000*x)+it
          call random_seed(size=one)
          call random_seed(put=iseed(:,it))
          call random_seed(get=iseed(:,it))
          do k=1+(it-1)*kpt,min(nc,it*kpt)
             call random_number(cube(1:nc_node_dim,:,k))
          enddo
       enddo
       !$omp end parallel do

    !! Generate random Gaussian numbers
    !$omp parallel do default(shared) private(i,j,k,x1,x2)
    do k=1,nc_dim_node,2
       do j=1,nc_dim_node
          do i=1,nc_dim_node,2
             x1=2*pi*cube(i,j,k)
             x2=sqrt(-2*log(cube(i+1,j,k)))
             cube(i,j,k)=x2*cos(x1)
             cube(i+1,j,k)=x2*sin(x1)
          enddo
       enddo
       do j=1,nc_dim_node
          do i=1,nc_dim_node,2
             x1=2*pi*cube(i+1,j,k+1)
             x2=sqrt(-2*log(cube(i,j,k+1)))
             cube(i,j,k+1)=x2*cos(x1)
             cube(i+1,j,k+1)=x2*sin(x1)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Forward FFT white noise field
    call fftw(1)

    !! Check noise
    !! This improves sampling at large scales (small k)
    call checknoise


    call cpu_time(time2)
    time2=(time2-time1)/nt
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called noisemap'
    return
  end subroutine noisemap

