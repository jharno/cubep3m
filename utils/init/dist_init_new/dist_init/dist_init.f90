
!! init.f90 written by Hy Trac on Sept 23, 2003
!! efc -fpp2 init_new.f90 /opt/fftw-3.0.1/lib/libfftw3f.a -openmp -tpp1 -O3 -o init.x
!! Parallelized :: Hugh Merz :: 041216 

program main
  implicit none

  include 'mpif.h'
  include 'dist_init.fh'

  call mpi_initialize

  !$ call OMP_SET_NUM_THREADS(nt)

  if (rank == 0)  call writeparams
  call initvar

  call noisemap

!! noise is in slab

  call transferfnc
  call powerspectrum

!! init contains delta field in fourier space!

  call kernel
!! BUG AFTER HERE

!! kernel is in slab

  call zeldovich

  call densitystatistics

!! gas density is in cube

  if (rank == 0)  call writeps
  call writexv
  call writedelta

  call fftw(0)
  call mpi_finalize(ierr)

contains

  subroutine initvar
    implicit none
    integer k

    real time1,time2
    call cpu_time(time1)

    !! Initialize large arrays in parallel
    !$omp parallel default(shared) private(k)
    !$omp do
    do k=1,nc_node_dim
       init(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc_node_dim
       phi(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,max_np
       xv(:,k)=0
    enddo
    !$omp end do
    !$omp end parallel

    !! Initialize fftw so that it generates the plans!
    firstfftw=.true.

    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar


  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nc      ', nc
    write(*,*) 'nodes   ', nodes
    write(*,*) 'np      ', np
    write(*,*) 'redshift ',redshift
    write(*,*) 'box      ',box
    write(*,*) 'sigma_8  ',s8
    write(*,*) 'omega_m  ',omegam
    write(*,*) 'omega_l  ',omegal
    write(*,*)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams


  subroutine writexv
    implicit none
    character*40 fn
    character(len=4) :: rank_s   
 
    real time1,time2
    call cpu_time(time1)

    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)

    fn=ic_path//'xv'//rank_s(1:len_trim(rank_s))//'.ic'
    !!write(*,*) 'rank',rank,'Writing ',fn
    open(11,file=fn,form='binary')
    write(11) np_local
    write(11) xv(:,1:np_local)
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called write xv'
    return
  end subroutine writexv


  subroutine writedelta
    implicit none
    character*40 fn
    character(len=4) :: rank_s

    real time1,time2
    call cpu_time(time1)

    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)

    fn=ic_path//'od'//rank_s(1:len_trim(rank_s))//'.ic'
    !!write(*,*) 'rank',rank,'Writing ',fn
    open(11,file=fn,form='binary')
    write(11) cube(:,:,:)
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called write delta'
    return
  end subroutine writedelta


  subroutine writeps
    implicit none
    integer k
    real kr,Da,pow,grow,power
    character*40 fn

    real time1,time2
    external grow,power

    call cpu_time(time1)


    !! Calculate growth factor at scalefactor a
    Da=scalefactor*grow(scalefactor)/grow(1.)


    !! Output power spectrum
    !! First column is k
    !! Second is \Delta^2 (unsampled)
    !! Third is gas p(k)
    !! Fourth is standard deviation
    !! Fifth is dm p(k)
    !! Sixth is standard deviation

    fn=ic_path//'ps1.init'
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    fn=ic_path//'ps2.init'
    write(*,*) 'Writing ',fn
    open(12,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       pow=power(kr)*Da**2
       write(11,*) kr,pow,psinit(1:2,k),psg(1:2,k),psdm(1:2,k)
       write(12,*) kr,pow,psinit(3:4,k),psg(3:4,k),psdm(3:4,k)
    enddo
    close(11)
    close(12)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write ps'
    return
  end subroutine writeps

!!------------------------------------------------------------------!!

end program main
