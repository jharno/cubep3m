#include "./preprocessor"
module Deriv3d
  use Parameters
  use Variables
  use mMPI
  use Pencil

  real, dimension(nc,nc_node_dim,nc_pen+2) :: slabx,slaby,slabz

contains

  subroutine gradient(grid,griddim)
    implicit none
    real, dimension(Ncells,Ncells,Ncells), intent(in) :: grid
    real, dimension(Ncells,Ncells,Ncells,3), intent(out) :: griddim

    integer :: k,j,i,kg,jg,ig,mg,l,ind,dx,dxy
    real :: kz,ky,kx,kr,kphys

    real(8) :: timer_beg,timer_end

    real, parameter :: ncr = 1.0*nc
    integer, parameter :: hc = nc/2
    
    call get_time(timer_beg)

    cube = grid
    call reset_pencil
    call cp_fftw(1)

    ind = 0
    dx = fsize(1)
    dxy = dx * fsize(2)

    !$omp parallel default(none) &
    !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kphys,kr,l) &
    !$omp& shared(slabx,slaby,slabz,slab,dx,dxy,mypadd,fstart)
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
             
             kx = 2*sin(pi*kx/ncr)
             ky = 2*sin(pi*ky/ncr)
             kz = 2*sin(pi*kz/ncr)
             kr = sqrt(kx**2+ky**2+kz**2)

             if (kr.eq.0) cycle

             slabx(i+1,j,k) = -1.0*slab(i,j,k)*kx/kr**2.0
             slaby(i+1,j,k) = -1.0*slab(i,j,k)*ky/kr**2.0
             slabz(i+1,j,k) = -1.0*slab(i,j,k)*kz/kr**2.0

             slabx(i,j,k) = slab(i+1,j,k)*kx/kr**2.0
             slaby(i,j,k) = slab(i+1,j,k)*ky/kr**2.0
             slabz(i,j,k) = slab(i+1,j,k)*kz/kr**2.0

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    call force_pencil
    slab = slabx
    call cp_fftw(-1)
    griddim(:,:,:,1) = -1.0*cube
    
    call force_pencil
    slab = slaby
    call cp_fftw(-1)
    griddim(:,:,:,2) = -1.0*cube
    
    call force_pencil
    slab = slabz
    call cp_fftw(-1)
    griddim(:,:,:,3) = -1.0*cube

#if VERBOSITY>0
    if (rank.eq.0) then
       call get_time(timer_end)
       write(*,*) '[mod gradient - subroutine divergence] Time taken: ',timer_end-timer_beg,'sec'
    end if
#endif

    return
  end subroutine gradient

  subroutine divergence(griddim,grid)
    implicit none
    real, dimension(Ncells,Ncells,Ncells), intent(out) :: grid
    real, dimension(Ncells,Ncells,Ncells,3), intent(in) :: griddim
    integer :: dim

    real(8) :: timer_beg,timer_end
    
    integer :: k,j,i,kg,jg,ig,mg,l,ind,dx,dxy
    real :: kz,ky,kx,kr,kphys

    real, parameter :: ncr = 1.0*nc
    integer, parameter :: hc = nc/2

    call get_time(timer_beg)

    grid = 0.0

    call reset_pencil
    cube = griddim(:,:,:,1)
    call cp_fftw(1)
    slabx = slab

    call reset_pencil
    cube = griddim(:,:,:,2)
    call cp_fftw(1)
    slaby = slab

    call reset_pencil
    cube = griddim(:,:,:,3)
    call cp_fftw(1)
    slabz = slab

    ind = 0
    dx = fsize(1)
    dxy = dx * fsize(2)

    !$omp parallel default(none) &
    !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kphys,kr,l) &
    !$omp& shared(slabx,slaby,slabz,slab,dx,dxy,mypadd,fstart)
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

             if (kr.eq.0) cycle

             slab(i+1,j,k) = -1.0*(slabx(i,j,k)*kx + slaby(i,j,k)*ky &
                  &+ slabz(i,j,k)*kz)

             slab(i,j,k) = (slabx(i+1,j,k)*kx + slaby(i+1,j,k)*ky &
                  &+ slabz(i+1,j,k)*kz)

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    call force_pencil
    call cp_fftw(-1)
    grid = cube

#if VERBOSITY>0
    if (rank.eq.0) then
       call get_time(timer_end)
       write(*,*) '[mod gradient - subroutine divergence] Time taken: ',timer_end-timer_beg,'sec'
    end if
#endif
    return    

  end subroutine divergence

end module Deriv3d
