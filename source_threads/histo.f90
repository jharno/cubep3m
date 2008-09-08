implicit none
integer, parameter :: dp=2000
integer, parameter :: nbin=30
real pd(6,dp),sep,histo(4,nbin)
integer i,isep
open(10,file='pair_F.dat')
do i=1,dp
read(10,'(6f16.8)') pd(:,i)
enddo
histo=0.
do i=1,dp
  sep=pd(1,i)
  isep=int(sep*10.0)+1
  histo(1,isep)=histo(1,isep)+1.
  histo(2,isep)=histo(2,isep)+sep
  histo(3,isep)=histo(3,isep)+pd(2,i) 
  histo(4,isep)=histo(4,isep)+pd(3,i)
enddo
do i=1,nbin
  histo(2:4,i)=histo(2:4,i)/histo(1,i)
  print *,histo(:,i)
enddo
do i=1,nbin
!  print *,histo(2,i),(histo(4,i)-histo(3,i))/histo(4,i)
enddo
end
