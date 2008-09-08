implicit none
real :: xv(6)
real*8 :: xva(6)
integer :: i,np
open(10,file='0.000xv.dat',form='binary')
read(10) np
print *,np
xva=0.0
do i=1,np
  read(10) xv(:)
  xva=xva+xv
enddo
close(10)
xva=xva/real(np)
print *,xva
end
