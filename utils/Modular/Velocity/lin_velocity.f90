#include './preprocessor'
program lin_velocity
  use Parameters
  use Variables
  use mMPI
  use FieldIO
  use Pencil
  use Transfer
  implicit none

  !Density file
  character(len=100) :: fnd

  !Transfer functions
  integer, parameter :: ntf = 7
  integer, parameter :: ndtfk = 1000
  integer, parameter :: nvtfk = 1000
  real, dimension(ntf,ndtfk) :: dtf
  real, dimension(ntf,nvtfk) :: vtf
  character(len=*), parameter :: fdtf = dir//'transfer/ith2_mnu0p05_z0p05_tk.dat'
  character(len=*), parameter :: fvtf = dir//'transfer/ith2_mnu0p05_z0p05_v_tk.dat'

  character(len=*), parameter :: folder = 'nu'

  integer, parameter :: kcol = 1
  integer, parameter :: dcol = 2
  integer, parameter :: ncol = 6

  integer :: i,j,k

  call start_mpi

  if (rank==0) write(*,*) 'Starting program lin_velocity'

  call setup_pencil
  call omp_set_num_threads(Nomp)

  !$omp workshare
  cube=0.0
  !$omp end workshare

  !Read in transfer functions
  if (rank==0) write(*,*) 'Reading in transfer functions'
  call read_transfer(dtf,fdtf)
  call read_transfer(vtf,fvtf)

  vtf(2:nvtfk,:) = -1.0*vtf(2:nvtfk,:)

  !Read in dm density field and fourier transform
  if (rank==0) write(*,*) 'Reading in density fields'
  fnd = dir//'fields/dm/den/'//redshift//'den'//trim(adjustl(rank_s))//'.dat'
  call read_field3(cube, trim(adjustl(fnd)))
  call cp_fftw(1)

  !call apply_tf(slab,dtf(kcol,:),(vtf(dcol,:)-vtf(ncol,:))/dtf(dcol,:))
  call apply_tf(slab,dtf(kcol,:),(vtf(ncol,:))/dtf(dcol,:))   
#ifdef DEBUG
  if (rank==0) write(*,*) slab(1:2,1:2,1:2)
  do k=1,size(slab,dim=3)
     do j=1,size(slab,dim=2)
        do i=1,size(slab,dim=1)
           if ( isnan(slab(i,j,k)) ) write(*,*) "ISNAN @",i,j,k
        end do
     end do
  end do
#endif

  !$omp workshare
  slab2 = slab
  !$omp end workshare

  !Compute vel_x
  call compute_slab_dim(slab2,slab,1)
  call cp_fftw(-1)

  !Write vel_x
  fnd = dir//'fields/'//folder//'/vel/'//redshift//'velx'//trim(adjustl(rank_s))//'.dat'
  call write_field3(cube,fnd)

  !Compute vel_y
  call compute_slab_dim(slab2,slab,2)
  call force_pencil
  call cp_fftw(-1)

  !Write vel_y
  fnd = dir//'fields/'//folder//'/vel/'//redshift//'vely'//trim(adjustl(rank_s))//'.dat'
  call write_field3(cube,fnd)

  !Compute vel_z
  call compute_slab_dim(slab2,slab,3)
  call force_pencil
  call cp_fftw(-1)

  !Write vel_z
  fnd =dir//'fields/'//folder//'/vel/'//redshift//'velz'//trim(adjustl(rank_s))//'.dat'
  call write_field3(cube,fnd)

  call cp_fftw(0)
  if (rank==0) write(*,*) 'Finished program lin_velocity'
  call end_mpi

contains

  subroutine compute_slab_dim(slabI,slabO,dim)
    implicit none
    real, dimension(:,:,:), intent(in) :: slabI
    real, dimension(:,:,:), intent(out) :: slabO
    integer, intent(in) :: dim

    integer :: k,j,i,kg,jg,ig,mg,ind,dx,dxy
    real :: kz,ky,kx,kr

    real, parameter :: ncr = 1.0*nc
    integer, parameter :: hc = nc/2

    if (rank==0) write(*,*) 'Entering subroutine compute_slab_dim with dim',dim

    ind = 0
    dx = fsize(1)
    dxy = dx * fsize(2)

    if (dim==1) then 
       !$omp parallel default(none) &                                                                                                                                                        
       !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kr) &                                                                                                             
       !$omp& shared(slabO,slabI,dx,dxy,mypadd,fstart)                                                                                                                                   
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

                kr = (kx**2+ky**2+kz**2)**0.5
                if (kr.eq.0) kr = 0.00001

!                slabO(i:i+1,j,k) = slabI(i:i+1,j,k)*kx/kr
                slabO(i,j,k) = -1.0*slabI(i+1,j,k)*kx/kr
                slabO(i+1,j,k) = slabI(i,j,k)*kx/kr

             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
    else if (dim==2) then
       !$omp parallel default(none) &                                                                                                                                                        
       !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kr) &                                                                                                             
       !$omp& shared(slabO,slabI,dx,dxy,mypadd,fstart)                                                                                                                                   
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

                kr = (kx**2+ky**2+kz**2)**0.5
                if (kr.eq.0) kr = 0.00001

!                slabO(i:i+1,j,k) = slabI(i:i+1,j,k)*ky/kr
                slabO(i,j,k) = -1.0*slabI(i+1,j,k)*ky/kr
                slabO(i+1,j,k) = slabI(i,j,k)*ky/kr

             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
    else if (dim==3) then
       !$omp parallel default(none) &                                                                                                                                                        
       !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kr) &                                                                                                             
       !$omp& shared(slabO,slabI,dx,dxy,mypadd,fstart)                                                                                                                                   
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

                kr = (kx**2+ky**2+kz**2)**0.5
                if (kr.eq.0) kr = 0.00001

!                slabO(i:i+1,j,k) = slabI(i:i+1,j,k)*kz/kr
                slabO(i,j,k) = -1.0*slabI(i+1,j,k)*kz/kr
                slabO(i+1,j,k) = slabI(i,j,k)*kz/kr
             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
    else 
       if (rank==0) write(*,*) 'Severe error in program lin_velocity in subroutine compute_slab_dim: unknown dim'
       call mpi_abort(mpi_comm_world,ierr,ierr)
    end if

    if (rank==0) write(*,*) 'Finished subroutine compute_slab_dim'

  end subroutine compute_slab_dim

end program lin_velocity
