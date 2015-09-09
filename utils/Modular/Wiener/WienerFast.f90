program WienerFast
  use Parameters
  use Variables
  use mMPI
  use Pencil
  use Power
  use FieldIO
  use Transfer
  use Deriv3d
  implicit none

  !!Field Names
  character(len=*), parameter :: odir = '/scratch2/p/pen/dinman/HaloHalo/cubep3m_movie_th2_dm/fieldsfast/'
  character(len=*), parameter :: fdir = '/scratch2/p/pen/dinman/HaloHalo/cubep3m_movie_th2_dm/'

  character(len=*), parameter :: vsuff = '_dm-hasim'
  character(len=*), parameter :: dsuff = '_h'
!!$  character(len=*), parameter :: dsuff2 = '_h'

  character(len=100) :: fdelta,fdelta2,fvelx,fvely,fvelz,fn

  !!Power spectrum variables
  integer, parameter :: nbins = 16
  real, dimension(nbins) :: pow,kpow,kbins,kernel,pii,pri,prr,pww,pvv,pvw
  real, dimension(nbins) :: b,n2,w
  integer :: nempty

  !Fields
  real, dimension(Ncells,Ncells,Ncells,3) :: cubedim,cubedim2
  real, dimension(Ncells,Ncells,Ncells) :: cube2,cube3

  !!Useful variables
  integer :: i,j,k,stat
  real :: p,q,r
  character(len=*),dimension(3),parameter :: xyz = (/'x','y','z'/)

  !Start MPI and set field names
  call start_mpi
  call setup_pencil
  call omp_set_num_threads(Nomp)
  if (rank==0) write(*,*) 'Beginning program WienerFast'

  !Setup kbins
  do k=1,nbins
     kbins(k) = (k-1.0)*log10(nc/2-1.0)/(nbins-1.0)
  end do
  kbins = (2.0*pi/Lbox)*(10.0**kbins)

  pii = 0
  pri = 0
  prr = 0
  b = 0
  n2 = 0
  w = 0

  !Field variables
  fdelta = fdir//'node'//trim(adjustl(rank_s))//'/0.000den'//trim(adjustl(rank_s))//dsuff//'.bin'  
!!$  fdelta2 = fdir//'node'//trim(adjustl(rank_s))//'/0.000den'//trim(adjustl(rank_s))//dsuff2//'.bin'

  fvelx = fdir//'node'//trim(adjustl(rank_s))//'/0.000velx'//trim(adjustl(rank_s))//vsuff//'.bin'
  fvely = fdir//'node'//trim(adjustl(rank_s))//'/0.000vely'//trim(adjustl(rank_s))//vsuff//'.bin'
  fvelz = fdir//'node'//trim(adjustl(rank_s))//'/0.000velz'//trim(adjustl(rank_s))//vsuff//'.bin'

  !Read in delta (cdm)
  if (rank==0) write(*,*) 'Reading in density contrast: '//trim(adjustl(fdelta))
  call read_field3(cube,trim(adjustl(fdelta)))
  cube = cube - 1.0

!!$  cube2 = cube
!!$  if (rank==0) write(*,*) 'Reading in density contrast: '//trim(adjustl(fdelta2))
!!$  call read_field3(cube,trim(adjustl(fdelta2)))
!!$  cube = cube - 1.0
!!$  cube = cube2 - cube

  !Fourier transform
  call cp_fftw(1)
  !Compute power spectrum
  pow = 0; kpow = 0; kernel = 0
  call autopowerspectrum(pow,kpow,kbins,slab,kernel)
  prr = pow

  !Store in slab2
  slab2 = slab  

  !Read in velocity fields
  call read_field3(cubedim(:,:,:,1),trim(adjustl(fvelx)))
  call read_field3(cubedim(:,:,:,2),trim(adjustl(fvely)))
  call read_field3(cubedim(:,:,:,3),trim(adjustl(fvelz)))

  call divergence(cubedim,cube2)
  
  cube = cube2
  call reset_pencil
  call cp_fftw(1) !Slab stores velocity divergence

  pow = 0; kpow = 0; kernel = 0
  call autopowerspectrum(pow,kpow,kbins,slab,kernel)
  pii = pow

  pow = 0; kpow = 0; kernel = 0
  call crosspowerspectrum(pow,kpow,kbins,slab,slab2,kernel)
  pri = pow

  !Compute bias, noise(^2) and wiener
  !pii = divergence; prr = delta; pri = cross
  b = pri/prr
  n2 = pii - b**2 * prr
  w = b*(pii)/(pii+n2)

  !Compute Wiener field
  slab = slab2
  call wiener(slab,kpow,w)
  call force_pencil
  call cp_fftw(-1) 
  cube3 = cube
  !Take gradient
  call gradient(cube3,cubedim2) !cubedim2 stores grad T delta
  !Compute cross power
  pww = 0
  pvv = 0
  pvw = 0
  do i=1,3
     call reset_pencil
     cube = cubedim2(:,:,:,i)
     call cp_fftw(1)
     slab2 = slab !Stores wiener field
     pow = 0; kpow = 0; kernel = 0
     call autopowerspectrum(pow,kpow,kbins,slab2,kernel)
     pww = pww+pow

     call reset_pencil
     cube = cubedim(:,:,:,i)
     call cp_fftw(1)
     pow = 0; kpow = 0; kernel = 0
     call autopowerspectrum(pow,kpow,kbins,slab,kernel)
     pvv = pvv+pow

     pow = 0; kpow = 0; kernel = 0
     call crosspowerspectrum(pow,kpow,kbins,slab,slab2,kernel)
     pvw = pvw + pow
  end do

  !Write out Wiener Fields
  do i=1,3  
     fn = odir//'fields/vel'//xyz(i)//trim(adjustl(rank_s))//dsuff//vsuff//'.bin'
     open(unit=11,file=fn,status='replace',access='stream',iostat=stat)
     if (stat.ne.0) ERROR STOP
     write(11) cubedim2(:,:,:,i)
     close(11)
  end do

  !Write filters
  if (rank.eq.0) then
     fn = odir//'filters'//dsuff//vsuff//'.txt'
     open(unit=11,file=fn,status='replace',iostat=stat,recl=2000)
     if (stat/=0) ERROR STOP
     do k=1,nbins
        write(11,*) kpow(k),prr(k),pii(k),pri(k),b(k),n2(k),w(k),pvv(k),pww(k),pvw(k)
     end do
     close(11)
  end if
  if (rank==0) write(*,*) 'Fin'
  call end_mpi


contains

  subroutine wiener(slabio,wk,wf)
    implicit none
    real, dimension(nc, nc_node_dim, nc_pen+2), intent(inout) :: slabio
    real, dimension(nbins), intent(in) :: wk,wf

    integer :: k,j,i,kg,jg,ig,mg,l,ind,dx,dxy
    real :: kz,ky,kx,kr,kphys

    integer :: nk,nw
    real :: w,m,y1,x1

    real, parameter :: ncr = 1.0*nc
    integer, parameter :: hc = nc/2
    
    nw = size(wk)
    if (nw.ne.size(wf)) ERROR STOP

    ind = 0
    dx = fsize(1)
    dxy = dx * fsize(2)

    !$omp parallel default(none) &
    !$omp& private(k,j,i,kg,mg,jg,ig,ind,kz,ky,kx,kphys,kr,nk,m,x1,y1,w) &
    !$omp& shared(slabio,dx,dxy,mypadd,fstart,wk,wf,rank,ierr,nw)
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

             kr = sqrt(kx**2+ky**2+kz**2)
             kphys = kr*(2.0*pi/Lbox)
             
!!$             kx = 2*sin(pi*kx/ncr)
!!$             ky = 2*sin(pi*ky/ncr)
!!$             kz = 2*sin(pi*kz/ncr)
!!$             kr = sqrt(kx**2+ky**2+kz**2)

             if (kr.eq.0) cycle
             !Compute wiener filter
             if ( kphys .gt. wk(nw) ) then
                nk = nw
             else if ( kphys.lt.wk(1) ) then
                nk = 2
             else
                nk = (kphys/wk(nw)) *nc/2.0
                if (nk.ge.nw) then
                   nk = nw-1
                else if (nk.le.1) then
                   nk = 2
                end if
                do
                   if (nk.gt.nw .or. nk.le.1) then
                      write(*,*) 'ERROR in wiener subroutine!', &
                           &nk,nw,kphys,wk(nk-1),wk(nk),wf(nk-1),wf(nk),w
                      call mpi_abort(mpi_comm_world, ierr, ierr)
                   end if
                   if ( kphys .lt. wk(nk) ) then
                      if (kphys.ge.wk(nk-1)) then
                         exit
                      else
                         nk = nk-1
                      end if
                   else
                      nk = nk + 1
                   end if
                end do
             end if
             if (wf(nk).gt.0 .and. wf(nk-1).gt.0) then
                m = log10( wf(nk)/wf(nk-1) ) / log10( wk(nk)/wk(nk-1) )
                y1 = wf(nk-1)
                x1 = wk(nk-1)
            
                w = y1 * (kphys/x1)**m
             else
                w = 0.0
             end if
             if (isnan(w)) write(*,*) m,y1,x1
             
!!$             if (rank.eq.0) write(*,*) nk,nw,kphys,wk(nk-1),wk(nk),wf(nk-1),wf(nk),w

             slabio(i:i+1,j,k) = slabio(i:i+1,j,k)*w
             
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel    

  end subroutine wiener
  

end program WienerFast
