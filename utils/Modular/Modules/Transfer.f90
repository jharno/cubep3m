module Transfer
  use Parameters
  use Variables
  use mMPI
  use Pencil
  implicit none

contains

  !Reads a transfer function file
  subroutine read_transfer(tfa,fn)
    implicit none
    character(len=*), intent(in) :: fn
    real, dimension(:,:), intent(out) :: tfa
    integer :: nk, nt, i, stat

    if (rank==0) write(*,*) 'Entering subroutine read_transfer'

    nt = size(tfa,dim=1)
    nk = size(tfa,dim=2)

    open(unit=21,file=fn,status='old',iostat=stat)
    if (stat/=0) then
       write(*,*) 'ERROR in module Transfer in subroutine read_transfer with file: '//fn
       call mpi_abort(mpi_comm_world,ierr,ierr)
    end if
    do i=1,nk
       read(21,*) tfa(:,i)
    end do
    close(21)

    if (rank==0) write(*,*) 'Finished subroutine read_transfer'

  end subroutine read_transfer

  !Applies a transfer function to field slabI.
  subroutine apply_tf(slabI,ktf,tf)
    implicit none
    real, dimension(:,:,:), intent(inout) :: slabI
    real, dimension(:), intent(in) :: ktf, tf

    integer :: k,j,i,kg,jg,ig,mg,l,ind,dx,dxy,nk
    real :: kz,ky,kx,kr,kphys,interp,w1,w2

    real, parameter :: ncr = 1.0*nc
    integer, parameter :: hc = nc/2

    if (rank==0) write(*,*) 'Entering subroutine apply_tf'

    nk = size(tf)

    ind = 0
    dx = fsize(1)
    dxy = dx * fsize(2)

    !$omp parallel default(none) &
    !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kphys,kr,interp,w1,w2,l) &
    !$omp& shared(slabI,ktf,tf,nk,dx,dxy,mypadd,fstart)
    !$omp do schedule(dynamic)
    do k = 1, nc_pen+mypadd
       ind = (k-1)*nc_node_dim*nc/2
       do j = 1, nc_node_dim
          do i = 1, nc, 2
             kg = ind / dxy
             mg = ind - kg * dxy
             jg = mg / dx
             ig = mg - jg * dx
             kg = fstart(3) + kg
             jg = fstart(2) + jg
             ig = 2 * (fstart(1) + ig) - 1
             ind = ind + 1
             if (kg < hc+2) then
                kz=kg-1
             else
                kz=kg-1-nc
             endif
             if (jg < hc+2) then
                ky=jg-1
             else
                ky=jg-1-nc
             endif
             kx = (ig-1)/2
             
             kphys = 2*pi*sqrt(kx**2+ky**2+kz**2)/Lbox
             kx = 2*sin(pi*kx/ncr)
             ky = 2*sin(pi*ky/ncr)
             kz = 2*sin(pi*kz/ncr)
             kr = sqrt(kx**2+ky**2+kz**2)
             
             if (kphys <= ktf(1)) then
                interp = tf(1)
             else if (kphys >= ktf(nk)) then
                interp = tf(nk)
             else
                do l=1,nk-1
                   if (kphys <= ktf(l+1)) then
                      w1 = ktf(l+1) - kphys
                      w2 = kphys - ktf(l)
                      interp = (w1*tf(l)+w2*tf(l+1))/(ktf(l+1)-ktf(l))
                      exit
                   endif
                enddo
             endif

!             if (kr==0) then
!                kr = 0.000001
!             end if

             slabI(i:i+1,j,k) = slabI(i:i+1,j,k)*interp

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    if (rank==0) write(*,*) 'Finished subroutine apply_tf'

  end subroutine apply_tf

end module Transfer
