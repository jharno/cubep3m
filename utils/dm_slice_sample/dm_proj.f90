implicit none
  real*4, dimension(3) :: x
  integer,parameter :: nc=1728
  integer,parameter :: nc_proj=nc
  real,parameter :: zh=72.
  real,parameter :: zl=0.
  real,parameter :: eps=0.0001
  real,parameter :: zdepth=zh-zl
  integer,parameter :: bdepth=127
  real,parameter :: mag=real(nc_proj)/real(nc)
  real(4), allocatable, dimension(:,:) :: proj
  integer(1), allocatable, dimension(:,:) :: map
  integer i,j

  print *,'mag=',mag

  allocate(proj(nc_proj,nc_proj))
  allocate(map(nc_proj,nc_proj))

  proj=1.0

  open(12,file='0.000dm_slice.dat',form='binary')
  do
    read(12,end=432) x
    j=j+1
    do i=1,3
      if (abs(x(i))<eps) x(i)=sign(eps,x(i))
    enddo
    x=modulo(x,real(nc))
!    proj(floor(mag*x(1))+1,floor(mag*x(2))+1) = proj(floor(mag*x(1))+1, floor(mag*x(2))+1) + 1.0+(zdepth-int(x(3)))**-2.0
    proj(floor(mag*x(1))+1,floor(mag*x(2))+1) = proj(floor(mag*x(1))+1, floor(mag*x(2))+1) + 1.0/ (zdepth-x(3)+1.)**0.5
  enddo
432 continue
  close(12)
  print *,'read in',j,'particles' 

  print *,'min proj=',minval(proj)
  print *,'max proj=',maxval(proj)

  proj=log(proj)
  
  print *,'min log proj=',minval(proj)
  print *,'max log proj=',maxval(proj)

  proj=proj-minval(proj)
  proj=proj/maxval(proj)

!  map=bdepth*sqrt(sqrt(proj))
  map=bdepth*proj
  print *,'max map=',maxval(map)
  print *,'min map=',minval(map)
  deallocate(proj)

  open(10,file='xy_proj.pgm',status='replace')
  write(10,'(2hP5)')
  write(10,*) nc_proj,nc_proj
  write(10,*) bdepth 
  close(10)
  open(10,file='xy_proj.pgm',access='append',form='binary')
  write(10) map
  close(10)
  deallocate(map)
end
