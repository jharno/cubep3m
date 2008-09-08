!! This is the cubepm ic decomposer.
!! Breaks monolithic file into cubic sections

implicit none
integer, parameter 	:: nc=352    	!! number of cells total 
integer, parameter 	:: nn_dim=2  	!! number of nodes / dimension
integer, parameter 	:: nn=nn_dim**3 !! number of nodes total
real, parameter 	:: ncc=nc/nn_dim !! number of cells / cubic 
integer, parameter 	:: hc=nc/2   	!! half-grid length
integer, parameter 	:: np=hc**3  	!! number of particles total
real, dimension(6,np) 	:: xv     	!! particle list
integer 		:: nploc(nn_dim,nn_dim,nn_dim) 	!! each nodes particle count
real, dimension(6,np/nn*2,nn_dim,nn_dim,nn_dim) :: xv_out !! nodal particle list
character(len=40) fn
character(len=4) rank_s
integer i,j,k,m
real x,y,z

!! Read in particle positions

fn='xvp.init'
write(*,*) 'Reading ',fn
open(12,file=fn,form='binary')
read(12) xv
close(12)

!! Allocate particles to correct node

nploc=0
do i=1,np
  j=int(xv(1,i)/ncc)+1
  k=int(xv(2,i)/ncc)+1
  m=int(xv(3,i)/ncc)+1
  nploc(j,k,m)=nploc(j,k,m)+1
  xv_out(:,nploc(j,k,m),j,k,m)=xv(:,i)
  xv_out(1,nploc(j,k,m),j,k,m)=xv_out(1,nploc(j,k,m),j,k,m)-(j-1)*ncc
  xv_out(2,nploc(j,k,m),j,k,m)=xv_out(2,nploc(j,k,m),j,k,m)-(k-1)*ncc
  xv_out(3,nploc(j,k,m),j,k,m)=xv_out(3,nploc(j,k,m),j,k,m)-(m-1)*ncc
enddo

write(*,*) 'number of particles=',sum(nploc)

!! Write out each nodes initial conditions 

do k=1,nn_dim
  do j=1,nn_dim
    do i=1,nn_dim
      write(*,*) 'rank:',(i-1)+(j-1)*nn_dim+(k-1)*nn_dim**2,'nploc:',nploc(i,j,k)
      fn=''
      write(rank_s,'(i4)') (i-1)+(j-1)*nn_dim+(k-1)*nn_dim**2
      rank_s=adjustl(rank_s)
      fn='xv'//rank_s(1:len_trim(rank_s))//'.ic'
      open(12,file=fn,form='binary',status='new')
      write(12) nploc(i,j,k)
      write(12) xv_out(1:6,1:nploc(i,j,k),i,j,k)
      close(12)
    enddo
  enddo
enddo

end 

