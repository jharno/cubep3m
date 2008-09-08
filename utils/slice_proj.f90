! slice_proj.f90 - samples a checkpoint file and produces a 
! particle distribution from a thin slice of the original volume

implicit none

   character(len=*),parameter :: fpath='/mnt/scratch/eh1/merz/cubepm/recomposed/'
   character(len=*),parameter :: checkpoints='./checkpoints'
   character(len=80) :: ifile,ofile
   character(len=7) :: z_s
   integer(4), parameter :: max_input=100
   real(kind=4), parameter :: z_lb = 5.0  !lower bound of slice
   real(kind=4), parameter :: z_ub = 6.0  !upper bound of slice

   real(kind=4), dimension(6) :: xv
   integer(kind=4) :: i,j,cp,nploc,fstat,num_checkpoints
   real(kind=4) :: z
   real(kind=4), dimension(max_input) :: z_checkpoint

!! Read in checkpoints to recompose

  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) write(*,*) &
     'error opening checkpoint list file:',checkpoints
  do num_checkpoints=1,max_input
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
  enddo
41  num_checkpoints=num_checkpoints-1
51  close(11)
  write(*,*) 'checkpoints to slice:'
  do i=1,num_checkpoints
    write(*,'(f5.1)') z_checkpoint(i)
  enddo

  do cp=1,num_checkpoints

   write(z_s,'(f7.3)') z_checkpoint(cp)
   z_s=adjustl(z_s)

   ifile=fpath//z_s(1:len_trim(z_s))//'xv.dat'

   open (unit=13,file=ifile,form='binary')
   read(13) nploc

   ofile=fpath//z_s(1:len_trim(z_s))//'slice.dat'
   open (unit=11,file=ofile,form='formatted')

   do j=1,nploc
     read(13) xv
     if (xv(3).ge.z_lb.and.xv(3).lt.z_ub)  &
         write(11,'(6f20.10)') xv
   enddo

   close(13)
   close(11)

   write(*,*) 'done z=',z_checkpoint(cp)
 
  enddo

end 
