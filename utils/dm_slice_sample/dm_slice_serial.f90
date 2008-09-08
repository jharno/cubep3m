! dm_slice_serial.f90 - serially samples a checkpoint file and produces a 
! particle distribution from a thin slice of the original volume
! ifort -fpp -DBINARY dm_slice_serial.f90 -o dmss.x
implicit none

! frequently used parameters are found in this header file:
  include '../../parameters'

!list of redshift to generate slices for 

   character(len=*),parameter :: checkpoints=cubepm_root//'input/checkpoints' 

!lower and upper bound of particle positions to include in slice
!currently only slices in z-dimension

   real(kind=4), parameter :: z_lb = 9.0  !lower bound of slice
   real(kind=4), parameter :: z_ub = 11.0  !upper bound of slice
   real(kind=4), parameter :: ncc = nc/nodes_dim

   integer(kind=4), parameter :: nn = nodes_dim**3

! the rest are all internal variables

   character(len=80) :: ifile,ofile
   character(len=7) :: z_s
   character(len=4) :: rank_s
   integer(4), parameter :: max_input=100

  real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p,dt_pp_acc
               
  integer(4) :: np_local,np_current, nts, sim_checkpoint, sim_projection, sim_halofind

   real(kind=4), dimension(6) :: xv
   integer(kind=4) :: rank,i,j,k,cp,nploc,fstat,num_checkpoints
   real(kind=4) :: z
   real(kind=4), dimension(max_input) :: z_checkpoint
   integer :: node_coords(nn,3)           !! coordinates of node in cube

!! Calculate node_coords
do rank=0,nn-1
do k=1,nodes_dim
  do j=1,nodes_dim
    do i=1,nodes_dim
      if (rank == (i-1)+(j-1)*nodes_dim+(k-1)*nodes_dim**2)  &
           node_coords(rank,:)=(/(i-1),(j-1),(k-1)/)
    enddo
  enddo
enddo
enddo

!! Read in list of checkpoints to recompose

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

!! create a slab for each checkpoint

  do cp=1,num_checkpoints

   write(z_s,'(f7.3)') z_checkpoint(cp)
   z_s=adjustl(z_s)

   ifile=z_s(1:len_trim(z_s))//'dm_slice.dat'
!! open output_slice
   open(unit=11,file=ifile,form='formatted',status='replace')

!! loop over each rank

   do rank=0,nn-1
   if (node_coords(rank,3)*ncc<z_lb.and.(node_coords(rank,3)+1)*ncc>z_ub) then
     write(rank_s,'(i4)') rank
   rank_s=adjustl(rank_s)
     ifile=z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
#ifdef BINARY
   open (unit=13,file=ifile,form='binary')
#else
   open (unit=13,file=ifile,form='unformatted')
   !!! THIS WILL NOT WORK PROPERLY DUE TO RECORDS :(
#endif

#ifdef PPINT
     read(13) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
                 sim_projection,sim_halofind,mass_p
#else
     read(13) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,sim_checkpoint, &
                 sim_projection,sim_halofind,mass_p
#endif
 
!     ofile=output_path//z_s(1:len_trim(z_s))//'slice.dat'
!     open (unit=11,file=ofile,form='formatted')

     do j=1,np_local
       read(13) xv
       xv(1:3)=xv(1:3)+node_coords(rank,:)*ncc
       if (xv(3).ge.z_lb.and.xv(3).lt.z_ub)  &
           write(11,'(6f20.10)') xv
     enddo

     close(13)
   endif

   enddo

   close(11)
   write(*,*) 'done z=',z_checkpoint(cp)

  enddo

end 
