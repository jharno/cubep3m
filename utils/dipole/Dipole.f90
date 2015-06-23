! Dipole code in Coarray Fortran

!! module load intel/15.0.2 intelmpi/5.0.3.048
!! ifort -O3 -fpp -xHost -coarray=shared Dipole.f90
!! ifort -O3 -fpp -xHost -coarray=distributed Dipole.f90
!! export I_MPI_PROCESS_MANAGER=mpd
!! export FOR_COARRAY_NUM_IMAGES=nn**2
!! [use mpd.hosts]

!#define DARK
implicit none

! parameters
integer(4),parameter ::        nn = 2               ! number of node per dimension
integer(4),parameter ::        nc = 288             ! number of coarse grid per node per dimension
real(4),parameter ::           lnode = 50           ! length per node per dimension
integer(4),parameter ::        ngroup = 2           ! number of halo groups
integer(4),parameter ::        ndim = 3             ! number of dimensions
real(4),parameter ::           lbox = lnode * nn

real(4),parameter ::           r1=1, r2=sqrt(2.)/2, r3=sqrt(3.)/3
character(1),parameter ::      xyz(3)=(/'x','y','z'/)
character(*),parameter ::      redshift = '0.000'   ! file prefix
character(*),parameter ::      DirJD='/scratch2/p/pen/emberson/cubep3m/TH2/GPC/'
character(*),parameter ::      DirHR='/bgq/scratch/p/pen/haoran/tiannu/fields/'
character(*),parameter ::      DirQT='/scratch2/p/pen/quentinb/cubep3m/cubep3m_multinu/JDfields/node'
character(*),parameter ::      fn_out='./Dipole.txt' ! output file name
#ifdef DARK
  character(*),parameter ::    Sim_name='TD'
  character(*),parameter ::    dir_apdx='0'
#else
  character(*),parameter ::    Sim_name='TN'
  character(*),parameter ::    dir_apdx=''
#endif

! variables
character(len=8) ::            srank ! rank string
character(len=1000) ::         fn
integer ::                     nrank,irank,crank(3),c1,c2,c3, i,j,k,l,m,n, zip_factor, igroup
real ::                        r     ! physical length
logical ::                     head

! coarrays
real(4) :: g(nc+1,nc+1,nc+1,ngroup)[nn,nn,*]
real(4) :: v(nc+1,nc+1,nc+1,ndim)[nn,nn,*], vnu(nc+1,nc+1,nc+1,ndim)[nn,nn,*]
real(4),codimension[nn,nn,*] :: xi1,xi2,xi3
real,dimension(nn,nn,nn),codimension[nn,nn,*] ::  xi1g,xi2g,xi3g

!==================================================================================================================!

!! coarray setup
nrank=num_images()
irank=this_image()
crank=this_image(g)
c1=crank(1); c2=crank(2); c3=crank(3);
head=(irank==1)
write(srank,'(i8)') irank-1 ! irank starts from 1
!! print some information
if (head) then
  print*, 'Coarray Fortran Dipole'
  print*, 'nn =', nn
  print*, 'nc =', nc
  print*, 'lnode =', lnode
  print*, 'lbox =', lbox
  print*, 'ngroup =', ngroup
  print*, 'Output: ', fn_out
  print*, '========================='
endif
if (nrank /= nn**3) then
  if (head) then
    print*, 'nodes per dim =',nn
    print*, 'num_images() =',nrank
    print*, 'incorrect number of images'
  endif
  stop
endif
sync all

!! read fields
fn=DirJD//'node'//trim(adjustl(srank))//'/'//redshift//'den'//trim(adjustl(srank))//'_nu.bin'
!fn=DirJD//'node'//trim(adjustl(srank))//fnd1//trim(adjustl(srank))//'_h_low.bin'
!fn=DirJD//'node'//trim(adjustl(srank))//fnd1//trim(adjustl(srank))//'_h.bin'
!fn=DirHR//'particle/den'//fnd1//'_nc1728_'//trim(adjustl(srank))//'_nu.dat'
!fn=DirHR//'ha/den/'//redshift//'hden_nc1728_'//Sim_name//trim(adjustl(srank))//'_G1.dat'
call read_field3(g(1:nc,1:nc,1:nc,1),fn)

fn=DirJD//'node'//trim(adjustl(srank))//'/'//redshift//'den'//trim(adjustl(srank))//'.bin'
!fn=DirJD//'node'//trim(adjustl(srank))//fnd1//trim(adjustl(srank))//'_h_high.bin'
!fn=DirHR//'particle/den'//fnd1//'_nc1728_'//trim(adjustl(srank))//'_dm.dat'
!fn=DirHR//'ha/den/'//redshift//'hden_nc1728_'//Sim_name//trim(adjustl(srank))//'_G2.dat'
call read_field3(g(1:nc,1:nc,1:nc,2),fn)
if (nc==288) g=g-1. !! for density, convert to density contrast

do i=1,ndim !! read velocity from each x,y,z component
  !fn=DirJD//'node'//trim(adjustl(srank))//fnv1//xyz(i)//trim(adjustl(srank))//'_nurha.bin'
  !fn=DirHR//'nu/vel/0.010vel'//xyz(i)//trim(adjustl(srank))//'.dat'
  !call read_field3(vnu(1:nc,1:nc,1:nc,i),fn)
  fn=DirJD//'node'//trim(adjustl(srank))//'/'//redshift//'vel'//xyz(i)//trim(adjustl(srank))//'_nu-dmsim.bin'
  !fn=DirHR//'p'//dir_apdx//'/vel/'//redshift//'vel'//xyz(i)//trim(adjustl(srank))//'.dat'
  !fn=DirQT//trim(adjustl(srank))//'/0.000vel'//xyz(i)//trim(adjustl(srank))//'_nu-dmrhasimrha.bin'
  call read_field3(v(1:nc,1:nc,1:nc,i),fn)
enddo
!v=v-vnu

if (head) then
  print*, 'g1',g(1,1,1,1)
  print*, 'g2',g(1,1,1,2)
  print*, 'vx',v(1,1,1,1)
  print*, 'vy',v(1,1,1,2)
  print*, 'vz',v(1,1,1,3)
endif

if (head) open(unit=76,file=fn_out,status='replace',form='formatted')
sync all

m=nc ! m is the computing range; m<=nc
do ! loop over scales
  r=lnode/m; xi1=0; xi2=0; xi3=0; n=m+1
  if (head) print*, 'compute ncc =', m
  g(n,:m,:m,:)=g(1,:m,:m,:)[mod(c1,nn)+1,c2,c3]; sync all
  g(:n,n,:m,:)=g(:n,1,:m,:)[c1,mod(c2,nn)+1,c3]; sync all
  g(:n,:n,n,:)=g(:n,:n,1,:)[c1,c2,mod(c3,nn)+1]; sync all
  v(n,:m,:m,:)=v(1,:m,:m,:)[mod(c1,nn)+1,c2,c3]; sync all
  v(:n,n,:m,:)=v(:n,1,:m,:)[c1,mod(c2,nn)+1,c3]; sync all
  v(:n,:n,n,:)=v(:n,:n,1,:)[c1,c2,mod(c3,nn)+1]; sync all
  !print*, 'r =', r, r*sqrt(2.), r*sqrt(3.)
  xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,1:m,1:m,:)),(/r1,0.,0./))*1.d0)
  xi1=xi1+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,1:m,1:m,:)),(/-r1,0.,0./))*1.d0)
  xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,2:n,1:m,:)),(/0.,r1,0./))*1.d0)
  xi1=xi1+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(1:m,1:m,1:m,:)),(/0.,-r1,0./))*1.d0)
  xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,1:m,2:n,:)),(/0.,0.,r1/))*1.d0)
  xi1=xi1+sum(g(1:m,1:m,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(1:m,1:m,1:m,:)),(/0.,0.,-r1/))*1.d0)
  !! sqrt(2)
  xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,2:n,2:n,:)),(/0.,r2,r2/))*1.d0)
  xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(1:m,1:m,2:n,:)),(/0.,-r2,r2/))*1.d0)
  xi2=xi2+sum(g(1:m,2:n,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,2:n,2:n,:)+v(1:m,1:m,1:m,:)),(/0.,-r2,-r2/))*1.d0)
  xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(1:m,2:n,1:m,:)),(/0.,r2,-r2/))*1.d0)
  xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,1:m,2:n,:)),(/r2,0.,r2/))*1.d0)
  xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,1:m,2:n,:)),(/-r2,0.,r2/))*1.d0)
  xi2=xi2+sum(g(2:n,1:m,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,1:m,2:n,:)+v(1:m,1:m,1:m,:)),(/-r2,0.,-r2/))*1.d0)
  xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(2:n,1:m,1:m,:)),(/r2,0.,-r2/))*1.d0)
  xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,2:n,1:m,:)),(/r2,r2,0./))*1.d0)
  xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,2:n,1:m,:)),(/-r2,r2,0./))*1.d0)
  xi2=xi2+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,2:n,1:m,:)+v(1:m,1:m,1:m,:)),(/-r2,-r2,0./))*1.d0)
  xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(2:n,1:m,1:m,:)),(/r2,-r2,0./))*1.d0)
  !! sqrt(3)  
  xi3=xi3+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,2:n,2:n,:)),(/r3,r3,r3/))*1.d0)
  xi3=xi3+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,2:n,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,2:n,2:n,:)),(/-r3,r3,r3/))*1.d0)
  xi3=xi3+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,2:n,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(2:n,1:m,2:n,:)),(/r3,-r3,r3/))*1.d0)
  xi3=xi3+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(2:n,2:n,1:m,:)+v(1:m,1:m,2:n,:)),(/-r3,-r3,r3/))*1.d0)
  xi3=xi3+sum(g(1:m,1:m,2:n,1)*g(2:n,2:n,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(2:n,2:n,1:m,:)),(/r3,r3,-r3/))*1.d0)
  xi3=xi3+sum(g(2:n,1:m,2:n,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(2:n,1:m,2:n,:)+v(1:m,2:n,1:m,:)),(/-r3,r3,-r3/))*1.d0)
  xi3=xi3+sum(g(1:m,2:n,2:n,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,2:n,2:n,:)+v(2:n,1:m,1:m,:)),(/r3,-r3,-r3/))*1.d0)
  xi3=xi3+sum(g(2:n,2:n,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,2:n,2:n,:)+v(1:m,1:m,1:m,:)),(/-r3,-r3,-r3/))*1.d0)
  sync all
  xi1g(c1,c2,c3)[1,1,1]=xi1 ! send to the head image
  xi2g(c1,c2,c3)[1,1,1]=xi2
  xi3g(c1,c2,c3)[1,1,1]=xi3; sync all


  if (head) then
    xi1=sum(xi1g); xi2=sum(xi2g); xi3=sum(xi3g) ! sum over images
    print*,'xi1,2,3 =', xi1/m**3/nn**3/6, xi2/m**3/nn**3/12, xi3/m**3/nn**3/8
    write(unit=76,fmt='(F12.6,E18.8)') r,          xi1/m**3/nn**3/6
    write(unit=76,fmt='(F12.6,E18.8)') r*sqrt(2.), xi2/m**3/nn**3/12
    write(unit=76,fmt='(F12.6,E18.8)') r*sqrt(3.), xi3/m**3/nn**3/8
  endif
  zip_factor=f(m)
!  if (m/zip_factor==1 .or. zip_factor==-1) exit
  if (m<zip_factor .or. zip_factor==-1) exit !! or forward data to head image (  )
  do igroup=1,ngroup
    call zip(g(:m,:m,:m,igroup),zip_factor)
  enddo
  call zip(v(:m,:m,:m,1),zip_factor)
  call zip(v(:m,:m,:m,2),zip_factor)
  call zip(v(:m,:m,:m,3),zip_factor)
  m=m/zip_factor
enddo

if (head) close(76)

!==========================================================================================================!
contains

function f(x)
  integer :: x,f
  f=-1
  if (mod(x,2)==0) then
    f=2
  elseif (mod(x,3)==0) then
    f=3
  endif
endfunction

function hat(vec)
  real :: vec(:,:,:,:), hat(size(vec,1),size(vec,2),size(vec,3),size(vec,4))
  hat = vec / spread(sqrt(sum(vec**2,dim=4)),dim=4,ncopies=3)
endfunction

function vdot(vec,vn)
  real :: vec(:,:,:,:), vn(3), vdot(size(vec,1),size(vec,2),size(vec,3))
  vdot=vec(:,:,:,1)*vn(1)+vec(:,:,:,2)*vn(2)+vec(:,:,:,3)*vn(3)
endfunction

subroutine zip(vec,x)
  integer :: i,x
  real :: vec(:,:,:), vv(size(vec,1)/x,size(vec,2),size(vec,3))
  vv=0
  do i=1,x
    vv=vv+vec(i::x,:,:)
  enddo
  vec(:size(vec,1)/x,:,:)=vv
  vv=0
  do i=1,x
    vv(:,:size(vec,2)/x,:)=vv(:,:size(vec,2)/x,:)+vec(:size(vec,1)/x,i::x,:)
  enddo
  vec(:size(vec,1)/x,:size(vec,2)/x,:)=vv(:,:size(vec,2)/x,:)
  vv=0
  do i=1,x
    vv(:,:,:size(vec,3)/x)=vv(:,:,:size(vec,3)/x)+vec(:size(vec,1)/x,:size(vec,2)/x,i::x)
  enddo
  vec(:size(vec,1)/x,:size(vec,2)/x,:size(vec,3)/x)=vv(:,:size(vec,2)/x,:size(vec,3)/x)/x**3
endsubroutine

subroutine read_field3(grid,fileI)
  implicit none
  real, dimension(:,:,:), intent(out) :: grid
  character(*), intent(in) :: fileI
  open(unit=21,file=trim(adjustl(fileI)),status='old',access='stream')
  !open(unit=21,file=trim(adjustl(fileI)),status='old',access='stream',convert='big_endian')
  read(21) grid
  close(21)
  !print*, 'read complete: ', trim(adjustl(fileI))
endsubroutine

subroutine pbc3(vec)
  implicit none
  integer :: n1,n2,n3
  real :: vec(:,:,:)
  n1=size(vec,1); n2=size(vec,2); n3=size(vec,3)
  vec(n1,:,:)=vec(1,:,:)
  vec(:,n2,:)=vec(:,1,:)
  vec(:,:,n3)=vec(:,:,1)
endsubroutine

end
