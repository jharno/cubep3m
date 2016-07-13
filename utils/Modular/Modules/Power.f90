#include "./preprocessor"
module Power
  use Parameters
  use Variables
  use mMPI
  use Pencil
  implicit none
  
  !Size of 4 byte float
  integer, parameter :: r4 = selected_real_kind(4)

contains
  
  subroutine autopowerspectrum(power, kaxis, kbins,slabA, kernel)
    implicit none
    real :: power(:), kaxis(:)

    !Optional argument: use if you do *NOT* want to deconvolve the particle interpolation scheme 
    !per k mode but instead at the end.  interpExponent should also be set to 0
    real, optional :: kernel(:)

    real, intent(in) :: kbins(:)
    real(8), allocatable :: powertask(:), kaxistask(:), powersum(:), kaxissum(:)
    real(8), allocatable :: kcount(:), kcounttask(:)
    real(8), allocatable :: kernelsum(:), kerneltask(:)
    integer :: Nbins
    integer :: i,j,k,ind,kg,mg,jg,ig,dx,dxy
    
    real :: kr, kx, ky, kz, sincx, sincy, sincz, sinc
    integer :: nk
    
    integer, parameter :: hc = nc/2
    real, parameter :: ncr = real(nc)
    real, parameter:: hcr = real(hc)
    
    real, dimension(nc, nc_node_dim, nc_pen+2) :: slabA
#ifndef NGP
    integer, parameter :: interpExponent = 4
#if (VERBOSITY>0)
    if (rank==0) write(*,*) "Assuming CIC particle interpolation"
#endif
#else
    integer, parameter :: interpExponent = 0
#if (VERBOSITY>0)
    if (rank==0) write(*,*) "Using NGP particle interpolation"
#endif
#endif
    
#if (VERBOSITY>0)
    if (rank==0) write(*,*) '[Mod - Power] Starting subroutine autopowerspectrum'
#endif

    Nbins = size(power)
    
    allocate(powertask(Nbins))
    allocate(powersum(Nbins))
    allocate(kaxistask(Nbins))
    allocate(kaxissum(Nbins))
    allocate(kcount(Nbins))
    allocate(kcounttask(Nbins))
    allocate(kerneltask(Nbins))
    allocate(kernelsum(Nbins))
    

    power = 0.0
    kaxis = 0.0
    if (present(kernel)) kernel=0.0
    powertask = 0.0
    kaxistask = 0.0
    powersum = 0.0
    kaxissum = 0.0
    kcount = 0
    kcounttask = 0
    kerneltask = 0
    kernelsum = 0

    dx = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
    
    !$omp parallel do &
    !$omp& default(none),shared(slabA,kbins,mypadd,fstart,dx,dxy,rank), &
    !$omp& private(k,j,i,ind,kg,mg,jg,ig,kx,ky,kz,kr,nk,sincx,sincy,sincz,sinc), &
    !$omp& reduction(+:powertask,kaxistask,kcounttask,kerneltask)
    do k=1,nc_pen+mypadd
       ind = (k-1)*nc_node_dim*hc
       do j=1, nc_node_dim
          do i=1,nc,2
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
             kr = 1.0*sqrt(kx**2+ky**2+kz**2)
             
             if(kx.eq.0 .and. ky <=0 .and. kz <=0)cycle;
             if(kx.eq.0 .and. ky >0 .and. kz <0)cycle;
             
             if (kr .gt. 0 .AND. kr < hc ) then
                sincx = merge(sin(pi*kx/nc)/(pi*kx/nc),1.0,kx/=0.0)
                sincy = merge(sin(pi*ky/nc)/(pi*ky/nc),1.0,ky/=0.0)
                sincz = merge(sin(pi*kz/nc)/(pi*kz/nc),1.0,kz/=0.0)
                sinc = (sincx*sincy*sincz)
                
                nk = nearest_k(2.0*pi*kr/Lbox, kbins)
                
                powertask(nk) = powertask(nk) + sum(slabA(i:i+1,j,k)*slabA(i:i+1,j,k))/sinc**interpExponent
                kaxistask(nk) = kaxistask(nk) + kr
                kcounttask(nk) = kcounttask(nk) + 1.0
                kerneltask(nk) = kerneltask(nk) + sinc**2.0
             end if
          end do
       end do
    end do
    !$omp end parallel do
    
    call mpi_reduce(powertask,powersum,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for powertask')
    
    call mpi_reduce(kaxistask,kaxissum,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for kaxistask')
    
    call mpi_reduce(kcounttask,kcount,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for kcounttask')
    
    call mpi_reduce(kerneltask,kernelsum,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for kerneltask')

    deallocate(powertask)
    deallocate(kaxistask)
    deallocate(kcounttask)
    deallocate(kerneltask)
    
    if (rank == 0) then
       do k=1,Nbins
          if (kcount(k)>0) then
             kaxissum(k) = 2.0*pi*kaxissum(k) / kcount(k) / Lbox
             kaxis(k) = real(kaxissum(k),kind=r4)
             powersum(k) = Lbox**3*(kaxissum(k)**3/2.0/pi**2) * powersum(k) / kcount(k) / ncr**6
             power(k) = real(powersum(k),kind=r4)
             if (present(kernel)) then
                kernel(k) = real(kernelsum(k)/kcount(k),kind=r4)
             end if
          else
             kaxis(k) = -1
             power(k) = -1
             if (present(kernel)) then
                kernel(k) = -1
             end if
          end if
       end do
    end if
    
    deallocate(kcount)
    deallocate(powersum)
    deallocate(kaxissum)
    deallocate(kernelsum)
    
    call mpi_bcast(power,Nbins,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(kaxis,Nbins,mpi_real,0,mpi_comm_world,ierr)
    if (present(kernel)) call mpi_bcast(kernel,Nbins,mpi_real,0,mpi_comm_world,ierr)

#if (VERBOSITY>0)
    if (rank==0) write(*,*) '[Mod - Power] Finished subroutine autopowerspectrum'
#endif
    
  end subroutine autopowerspectrum
  
  subroutine crosspowerspectrum(power, kaxis, kbins, slabA,slabB, kernel)
    implicit none
    real :: power(:), kaxis(:)

    !Optional argument: use if you do *NOT* want to deconvolve the particle interpolation scheme
    !per k mode but instead at the end.  interpExponent should also be set to 0
    real, optional :: kernel(:)

    real, intent(in) :: kbins(:)
    real(8), allocatable :: powertask(:), kaxistask(:), powersum(:), kaxissum(:)
    real(8), allocatable :: kcount(:), kcounttask(:)
    real(8), allocatable :: kerneltask(:), kernelsum(:)
    integer :: Nbins
    integer :: i,j,k,ind,kg,mg,jg,ig,dx,dxy
    
    real :: kr, kx, ky, kz, sincx, sincy, sincz, sinc
    integer :: nk

    integer, parameter :: hc = nc/2
    real, parameter :: ncr = real(nc)
    real, parameter:: hcr = real(hc)
    
    real, dimension(nc, nc_node_dim, nc_pen+2) :: slabA,slabB
    
#ifndef NGP
    integer, parameter :: interpExponent = 4
#if (VERBOSITY>0)
    if (rank==0) write(*,*) "Assuming CIC particle interpolation"
#endif
#else
    integer, parameter :: interpExponent = 2
#if (VERBOSITY>0)
    if (rank==0) write(*,*) "Using NGP particle interpolation"
#endif
#endif
    
#if (VERBOSITY>0)
    if (rank==0) write(*,*) '[Mod - Power] Starting subroutine crosspowerspectrum'
#endif
    
    Nbins = size(power)
    
    allocate(powertask(Nbins))
    allocate(kaxistask(Nbins))
    allocate(powersum(Nbins))
    allocate(kaxissum(Nbins))
    allocate(kcount(Nbins))
    allocate(kcounttask(Nbins))
    allocate(kerneltask(Nbins))
    allocate(kernelsum(Nbins))

    power = 0.0
    kaxis = 0.0
    powertask = 0.0
    kaxistask = 0.0
    powersum = 0.0
    kaxissum = 0.0
    kcount = 0
    kcounttask = 0
    kerneltask = 0.0
    kernelsum = 0.0

    dx = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
    
    !$omp parallel do &
    !$omp& default(none),shared(slabA, slabB, kbins,mypadd,fstart,dx,dxy,rank), &
    !$omp& private(k,j,i,ind,kg,mg,jg,ig,kx,ky,kz,kr,nk,sincx,sincy,sincz,sinc), &
    !$omp& reduction(+:powertask,kaxistask,kcounttask,kerneltask)
    do k=1,nc_pen+mypadd
       ind = (k-1)*nc_node_dim*hc
       do j=1, nc_node_dim
          do i=1,nc,2
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
             kr = 1.0*sqrt(kx**2+ky**2+kz**2)
             
             if(kx.eq.0 .and. ky <=0 .and. kz <=0)cycle;
             if(kx.eq.0 .and. ky >0 .and. kz <0)cycle;
             
             if (kr .gt. 0 .AND. kr < hc ) then
                sincx = merge(sin(pi*kx/nc)/(pi*kx/nc),1.0,kx/=0.0)
                sincy = merge(sin(pi*ky/nc)/(pi*ky/nc),1.0,ky/=0.0)
                sincz = merge(sin(pi*kz/nc)/(pi*kz/nc),1.0,kz/=0.0)
                sinc = (sincx*sincy*sincz)
                
                nk = nearest_k(2.0*pi*kr/Lbox, kbins)
                
                powertask(nk) = powertask(nk) + sum(slabA(i:i+1,j,k)*slabB(i:i+1,j,k))/sinc**interpExponent
                kaxistask(nk) = kaxistask(nk) + kr
                kcounttask(nk) = kcounttask(nk) + 1.0
                kerneltask(nk) = kerneltask(nk) + sinc**2.0
             end if
          end do
       end do
    end do
    !$omp end parallel do
    
    call mpi_reduce(powertask,powersum,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for powertask')
    
    call mpi_reduce(kaxistask,kaxissum,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for kaxistask')
    
    call mpi_reduce(kcounttask,kcount,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectrum at call to mpi_reduce &
         & for kcounttask')

    call mpi_reduce(kerneltask,kernelsum,Nbins,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    if (ierr/=0) call power_error_stop('Error in subroutine power spectruma t call to mpi_reduce &
         & for kernelsum')
    
    deallocate(powertask)
    deallocate(kaxistask)
    deallocate(kcounttask)
    deallocate(kerneltask)
    
    if (rank == 0) then
       do k=1,Nbins
          if (kcount(k)>0) then
             kaxissum(k) = 2.0*pi*kaxissum(k) / kcount(k) / Lbox
             kaxis(k) = real(kaxissum(k),kind=r4)
             powersum(k) = Lbox**3*(kaxissum(k)**3/2.0/pi**2) * powersum(k) / kcount(k) / ncr**6
             power(k) = real(powersum(k),kind=r4)
             if (present(kernel)) kernel(k) = real(kernelsum(k)/kcount(k),kind=r4)
          else
             kaxis(k) = -1
             power(k) = -1
             if (present(kernel)) kernel(k) = -1
          end if
       end do
    end if
    
    deallocate(kcount)
    deallocate(powersum)
    deallocate(kaxissum)
    deallocate(kernelsum)
    
    call mpi_bcast(power,Nbins,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(kaxis,Nbins,mpi_real,0,mpi_comm_world,ierr)
    if (present(kernel)) call mpi_bcast(kernel,Nbins,mpi_real,0,mpi_comm_world,ierr)

#if (VERBOSITY>0)
    if (rank==0) write(*,*) '[Mod - Power] Finished subroutine crosspowerspectrum'
#endif

  end subroutine crosspowerspectrum

  function nearest_k(kr, kaxis) result(i)
    implicit none
    real, intent(in) :: kr, kaxis(:)
    integer :: i
    integer :: n
    n = size(kaxis)
    i=n
    do while(kr < kaxis(i) .AND. i>1)
       i=i-1
    end do
    i = merge(i+1,i,kr - kaxis(i) > 0.5)
  end function nearest_k
  
  subroutine power_error_stop(expl)
    implicit none
    character(len=*) :: expl
    write(*,*) '[Mod - Power]'
    write(*,*) '-->'//expl
    call mpi_abort(mpi_comm_world, ierr, ierr)
  end subroutine power_error_stop
  
end Module Power
