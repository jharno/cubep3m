!compile with ifort den_proj_serial.f90 -fpp -O3 -xN -o den_proj_serial.x 
implicit none

  !! dir is the directory for files
  character(*), parameter :: ipl='/scratch/merz/p3m/0.000xv0.dat' !'/uber-scratch/merz/p3m_scaled/p3m/0.000xv0.dat' !'/uber-scratch/merz/cubepm_ngp_scaled_shake.1/0.000xv0.dat' !'/uber-scratch/merz/cubepm_8_320_100Mpc_060424.out/0.000xv.dat'
  character(*), parameter :: ops='/scratch/merz/p3m_fixed_' !'/scratch/merz/cic_'

  !! nc is the number of cells per box length
  real, parameter :: nc_orig=80. !320
  integer, parameter :: nc=512

  real, dimension(3),parameter :: offset=(/19.46673,-0.9722952,14.7756/) !offset=(/0.0,0.0,0.0/) 
  real :: CRUD(11) !10) 

  integer, parameter :: np=160**3
  real(4), dimension(6,np) :: xv
  real(4), dimension(nc,nc,nc) :: d
  real(4), dimension(nc,nc) :: proj
  integer(1), dimension(nc,nc) :: map
  character(len=80) :: fn
  integer :: cp,k,np_in,i,j
  real :: kr

  common /rarr/ xv,d,proj

  print *,'starting serial density projection program'
  d=0.0
  open(11,file=ipl,form='binary')
  read(11) np_in
  read(11) CRUD
  if (np/=np) then
    print *,'np_in .ne. np'
    stop
  endif
  print *,np_in
  print *,CRUD
  read(11) xv
  close(11)
  do i=1,np
    xv(:3,i)=(xv(:3,i)-offset(:))*(real(nc)/nc_orig)
  enddo
  call cicmass

  proj=0.0
  do k=1,nc
    proj=proj+d(:,:,k)
  enddo
  proj=proj-minval(proj)
  proj=proj/maxval(proj)
  map=127*sqrt(sqrt(proj))
  open(10,file=ops//'xy_proj.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc,nc
  write(10,*) 127
  close(10)
  open(10,file=ops//'xy_proj.pgm',access='append',form='binary')
  write(10) map
  close(10)

  proj=0.0
  do j=1,nc
    proj=proj+d(:,j,:)
  enddo
  proj=proj-minval(proj)
  proj=proj/maxval(proj)
  map=127*sqrt(sqrt(proj))
  open(10,file=ops//'xz_proj.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc,nc
  write(10,*) 127
  close(10)
  open(10,file=ops//'xz_proj.pgm',access='append',form='binary')
  write(10) map
  close(10)

  proj=0.0
  do i=1,nc
    proj=proj+d(i,:,:)
  enddo
  proj=proj-minval(proj)
  proj=proj/maxval(proj)
  map=127*sqrt(sqrt(proj))
  open(10,file=ops//'yz_proj.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc,nc
  write(10,*) 127
  close(10)
  open(10,file=ops//'yz_proj.pgm',access='append',form='binary')
  write(10) map
  close(10)
 
contains
 subroutine cicmass
    implicit none
    real, parameter :: ncr=nc
    real, parameter :: mp=(ncr)**3/np 

    integer i
    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    do i=1,np
       x=mod(xv(1,i)-0.5+ncr,ncr)
       y=mod(xv(2,i)-0.5+ncr,ncr)
       z=mod(xv(3,i)-0.5+ncr,ncr)

       i1=floor(x)+1
       i2=mod(i1,nc)+1
       dx1=i1-x
       dx2=1-dx1
       j1=floor(y)+1
       j2=mod(j1,nc)+1
       dy1=j1-y
       dy2=1-dy1
       k1=floor(z)+1
       k2=mod(k1,nc)+1
       dz1=k1-z
       dz2=1-dz1

       dz1=mp*dz1
       dz2=mp*dz2
       d(i1,j1,k1)=d(i1,j1,k1)+dx1*dy1*dz1
       d(i2,j1,k1)=d(i2,j1,k1)+dx2*dy1*dz1
       d(i1,j2,k1)=d(i1,j2,k1)+dx1*dy2*dz1
       d(i2,j2,k1)=d(i2,j2,k1)+dx2*dy2*dz1
       d(i1,j1,k2)=d(i1,j1,k2)+dx1*dy1*dz2
       d(i2,j1,k2)=d(i2,j1,k2)+dx2*dy1*dz2
       d(i1,j2,k2)=d(i1,j2,k2)+dx1*dy2*dz2
       d(i2,j2,k2)=d(i2,j2,k2)+dx2*dy2*dz2
    enddo
 end subroutine cicmass
end
