!! This is the cubepm particle list recomposer.
!! Accumulates distributed particle list into monolithic file 
!! compile with: mpif77 recompose.f90 -o recompose.x

implicit none

include 'mpif.h'

integer, parameter 	:: nc=320   	!! number of cells total 
integer, parameter 	:: nn_dim=2  	!! number of nodes / dimension
character(len=*),parameter :: output_path = '/scratch/merz/cubepm_35Mpc/'
character(len=*),parameter :: checkpoints = '/home/merz/codes/cubepm/input/checkpoints'
real, parameter 	:: rnc = nc
integer, parameter      :: max_input=100!! maximum number of checkpoints
integer, parameter 	:: nn=nn_dim**3 !! number of nodes total
real, parameter 	:: ncc=nc/nn_dim!! number of cells / cubic 
integer, parameter 	:: hc=nc/2   	!! half-grid length
integer, parameter 	:: np=hc**3  	!! maximum number of particles total
integer, parameter      :: max_np=np    !! maximum number of particles / node

real, dimension(6,max_np)  :: xv     	!! particle list

integer :: node_coords(3) 		!! coordinates of node in cube	
character(len=80) ofile 
character(len=4) rank_s
character(len=7) z_s
integer i,j,k,m
real x,y,z,z_write
integer(4) :: nn_returned, ierr,rank, num_checkpoints,fstat,cur_checkpoint
real(4), dimension(max_input) :: z_checkpoint

integer(4) :: np_local,np_current, nts, sim_checkpoint, sim_projection, sim_halofind
integer(4) :: nodes_returned,tag
real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p,np_total

integer(4), dimension(0:nn-1) :: np_found
integer(4), dimension(mpi_status_size) :: status

!! Initialize MPI
 
call mpi_init(ierr)
if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
if (nodes_returned /= nn ) then
  write(*,*) 'recompose compiled for a different number of nodes'
  write(*,*) 'mpirun nodes=',nodes_returned,'recompose nodes=',nn
  call mpi_abort(mpi_comm_world,ierr,ierr)
endif
call mpi_comm_rank(mpi_comm_world,rank,ierr)
if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

!! Calculate node_coords

do k=1,nn_dim
  do j=1,nn_dim
    do i=1,nn_dim
      if (rank == (i-1)+(j-1)*nn_dim+(k-1)*nn_dim**2)  &
           node_coords(:)=(/(i-1),(j-1),(k-1)/)
    enddo
  enddo
enddo

if (rank == 0) write(*,*) 'rank, cartesian coords in cube'

do i=0,nn-1
  if (rank == i) write(*,*) rank,node_coords
  call mpi_barrier(mpi_comm_world,ierr)
enddo

!! Read in checkpoints to recompose

if (rank == 0) then
  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) then
    write(*,*) 'error opening checkpoint list file'
    write(*,*) 'rank',rank,'file:',checkpoints
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
  do num_checkpoints=1,max_input
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
  enddo
41  num_checkpoints=num_checkpoints-1
51  close(11)
  write(*,*) 'checkpoints to recompose:'
  do i=1,num_checkpoints
    write(*,'(f5.1)') z_checkpoint(i)
  enddo
endif

call mpi_bcast(num_checkpoints,1,mpi_integer,0,mpi_comm_world,ierr)

do cur_checkpoint=1,num_checkpoints

!! Read in particle positions

  np_found=0
  np_current=0  

  if (rank == 0) then
    z_write = z_checkpoint(cur_checkpoint)
    write(*,*) 'processing z=',z_write
  endif

  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)

  write(rank_s,'(i4)') rank
  rank_s=adjustl(rank_s)

  ofile=output_path//z_s(1:len_trim(z_s))//'xv'// &
  	rank_s(1:len_trim(rank_s))//'.dat'

  open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')

  if (fstat /= 0) then
    write(*,*) 'error opening checkpoint'
    write(*,*) 'rank',rank,'file:',ofile
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,sim_checkpoint, &
               sim_projection,sim_halofind,mass_p

  if (np_local > max_np) then
    write(*,*) 'too many particles to store'
    write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

!! tally up total number of particles

  call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
                       mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'number of particles =', int(np_total,4)
  if (rank == 0 .and. int(np_total,4) > max_np) then
    write(*,*) 'np total > max_np',np_total,max_np
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

  call mpi_gather(np_local,1,mpi_integer,np_found,1,mpi_integer,0, &
                mpi_comm_world,ierr)

  if (rank == 0) then
    write(*,*) 'particles / node'
    do i=0,nn-1
      write(*,*) i, np_found(i)
    enddo
  endif

  read(21) xv(:,:np_local)
  close(21)


!! prepare particles for accumulation 

  do j=1,np_local
    xv(1:3,j)=modulo(xv(1:3,j)+node_coords(:)*ncc,rnc)
  enddo

!! send particles back to master node
  
  np_current=np_local
  do i=1,nn-1
    tag=i
    if (rank == 0) then
      call mpi_recv(xv(1,np_current+1),6*np_found(i),mpi_real, &
                    i,tag,mpi_comm_world,status,ierr)
    elseif (rank == i) then
      call mpi_send(xv,6*np_local,mpi_real,0,tag,mpi_comm_world,ierr)
    endif
    call mpi_barrier(mpi_comm_world,ierr)
    if (rank == 0) np_current=np_current+np_found(i)
  enddo
 
!! write out full particle list to disk

  if (rank == 0) then
    if (np_current /= int(np_total,4)) then
      write(*,*) 'total particles in /= np_total'
      write(*,*) np_current,int(np_total,4)
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    ofile=output_path//z_s(1:len_trim(z_s))//'xv.dat'
    open(23,file=ofile,form='binary',status='replace',iostat=fstat)
    if (fstat /= 0) then
      write(*,*) 'error opening composed checkpoint'
      write(*,*) 'rank',rank,'file:',ofile
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    write(23) np_current
    write(23) xv(:,:np_current)
    close(23) 
  endif

enddo

end 

