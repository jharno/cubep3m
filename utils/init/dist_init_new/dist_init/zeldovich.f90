  subroutine zeldovich
    implicit none

    include 'dist_init.fh'

    integer i,j,k

    real time1,time2
    call cpu_time(time1)

    !! Calculate Fourier space potential
    !! Complex multiply delta field with kernel
    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc_slab
       do j=1,nc
          do i=1,nc+2,2
             slab(i:i+1,j,k)=slab(i,j,k)*init(i:i+1,j,k)/(4*pi)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Inverse Fourier transform potential field
    call fftw(-1)

    !! potential is in cube

    if (rank == 0) write(*,*) cube(1:3,1:3,1:3)

    !! Calculate displacement field
    select case (np)
    case (nc)
       call displacement1
    case (hc)
       if (rank == 0) write(*,*) 'called displacement hc'
       call displacement2
    end select

    !! Compute dark matter power spectrum
    call dmpowerspectrum


    call cpu_time(time2)
    time2=(time2-time1)/nt
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called zeldovich'
    return
  end subroutine zeldovich


