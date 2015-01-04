! generate offsets, for -DREADOFFSET in cubep3m:update_position
implicit none
integer,parameter :: nts=100000
integer(4) :: seedsize, getseed(2), i,j,k,l
real(4) :: offset(3,nts)
call random_seed(size=seedsize)
call random_seed(put=(/647575,2048/))
call random_seed(get=getseed)
print*,'seedsize =',seedsize
print*,'getseed =',getseed
call random_number(harvest=offset)
open(11,file='offsets.dat',status='replace',access='stream')
write(11) offset
close(11)
print*,'wrote offsets.dat'
print*,'first 5 nts',offset(:,:5)
print*,'last  5 nts',offset(:,nts-4:)
end
