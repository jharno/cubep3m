! plot_proj.f90 - convert cubepm density projections to portable greymap format

implicit none
integer, parameter :: n=320 	      !length of box 
character(len=*),parameter :: projections='/home/merz/codes/cubepm/input/projections'
character(len=*),parameter :: ipath='/scratch/merz/cubepm/'   !input density projection path 
character(len=*),parameter :: opath='/scratch/merz/cubepm/pgms/'      !output pgm path 
integer(4),parameter :: max_input=100
character(len=max_input) :: ofile,ifile
real(4),dimension(max_input) :: z_projection
integer(4) :: num_projections,cur_projection
integer(4) :: i,dim
integer(4) :: fstat
real(4) :: a
real(4), dimension(n,n) :: den 
integer(1), dimension(n,n) :: map

!! Read in projections to recompose

  open(11,file=projections,status='old',iostat=fstat)
  if (fstat /= 0) then
    write(*,*) 'error opening projections list file'
    write(*,*) 'file:',projections
    stop 
  endif
  do num_projections=1,max_input
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_projection(num_projections)
  enddo
41  num_projections=num_projections-1
51  close(11)
  write(*,*) 'projections to recompose:'
  do i=1,num_projections
    write(*,'(f5.1)') z_projection(i)
  enddo

do cur_projection=1,num_projections

!! Now start doing each projection
do dim=1,3

  if (dim == 1) then
    write(ifile,'(f7.3,"proj_xy.dat")') z_projection(cur_projection)
    write(ofile,'(f7.3,"proj_xy.pgm")') z_projection(cur_projection)
    print *,'calculating xy'
  elseif (dim == 2) then
    write(ifile,'(f7.3,"proj_xz.dat")') z_projection(cur_projection)
    write(ofile,'(f7.3,"proj_xz.pgm")') z_projection(cur_projection)
    print *,'calculating xz'
  else
    write(ifile,'(f7.3,"proj_yz.dat")') z_projection(cur_projection)
    write(ofile,'(f7.3,"proj_yz.pgm")') z_projection(cur_projection)
    print *,'calculating yz'
  endif
  ifile=adjustl(ifile)
  ofile=adjustl(ofile)

open(10,file=ipath//ifile,form='binary',iostat=fstat)
if (fstat /= 0) then
  print *,'error opening density file:',ipath//ifile
  stop
endif
read(10) a
read(10) den 
close(10)
print *,'scalefactor=',a,'redshift=',1.0/a-1.0
print *,'total density =',sum(den)

den=den-minval(den)
den=den/maxval(den)

! one may wish to change the mapping function 

map=127*sqrt(sqrt(den))

open(10,file=opath//ofile,status='replace',iostat=fstat)
if (fstat /= 0) then
  print *,'error opening output pgm:',opath//ofile
  stop
endif
write(10,'(2hP5)')
write(10,*) n,n
write(10,*) 127
close(10)
open(10,file=opath//ofile,access='append',form='binary',iostat=fstat)
if (fstat /= 0) then
  print *,'error opening pgm file to append:',ofile
  stop
endif
write(10) map
close(10)
enddo
enddo
end
