!! This is the cubepm particle list recomposer.
!! Sep 27th 2005 :: added gas -- not finished... have to think about this one
!! Accumulates distributed particle list and gas checkpoints into a single checkpoint file 
!! compile with: mpif77 recompose.f90 -o recompose

implicit none
include 'mpif.h'

! frequently changed parameters are found in this header file:
include '../../parameters'

! list of redshift to recompose (same as in cubepm.par)

character(len=*),parameter :: checkpoints =cubepm_root//'input/checkpoints' 

real, parameter    :: rnc = nc
integer, parameter :: max_input=100 !! maximum number of checkpoints
integer, parameter :: nn=nodes_dim**3  !! number of nodes total
real, parameter    :: ncc=nc/nodes_dim !! number of cells / cubic 
integer, parameter :: hc=nc/2       !! half-grid length
integer, parameter :: max_np=6*(ncc/2)**3  !! maximum number of particles / node

real, dimension(6,max_np)  :: xv    !! particle list

#ifdef MHD_IN_PROGRESS
!.....type for dimensions, work load distribution, comm handles.........
      type comm_wld
        integer :: g !global number of zones in one direction
        integer :: r !global index of index 0 in local array
        integer :: m,n !start and end of local section without buffers
        integer :: l !dimension of array inclusive buffers
        integer, dimension(4) :: requests !communication handles
      end type comm_wld

integer :: cur_iter
real :: cur_t
type(comm_wld) :: nx,ny,nz
real :: u(5,ncc,ncc,ncc),b(3,ncc,ncc,ncc)
#endif

integer :: node_coords(3)           !! coordinates of node in cube	
character(len=80) ofile 
character(len=4) rank_s
character(len=7) z_s
integer i,j,k,m
real x,y,z,z_write
integer(4) :: nn_returned, ierr,rank, num_checkpoints,fstat,cur_checkpoint
real(4), dimension(max_input) :: z_checkpoint

integer(4) :: np_local,np_current, nts, sim_checkpoint, sim_projection, sim_halofind
integer(4) :: nodes_returned,tag
real(4) :: a,t,tau,dt_f_acc,dt_c_acc,dt_pp_acc,mass_p
integer*8 :: npl8,np_total

integer(4), dimension(0:nn-1) :: np_found
integer(4), dimension(mpi_status_size) :: status

#ifdef MHD_IN_PROGRESS
common /rarr/ xv,u,b
#else
common /rarr/ xv
#endif
!! Initialize MPI

print *,'size of major mem=',4*(6*max_np+8*ncc*ncc*ncc)/1024./1024.,'MB'
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

do k=1,nodes_dim
  do j=1,nodes_dim
    do i=1,nodes_dim
      if (rank == (i-1)+(j-1)*nodes_dim+(k-1)*nodes_dim**2)  &
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

!! Build up coherent set of checkpoint name strings

  if (rank == 0) then
    z_write = z_checkpoint(cur_checkpoint)
    write(*,*) 'processing z=',z_write
  endif

  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)
  write(rank_s,'(i4)') rank
  rank_s=adjustl(rank_s)

!! Try to open the dark matter checkpoint

  ofile=output_path//z_s(1:len_trim(z_s))//'xv'// &
  	rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
  open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
  open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif


  if (fstat /= 0) then
    write(*,*) 'error opening dark matter checkpoint'
    write(*,*) 'rank',rank,'file:',ofile
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

#ifdef MHD_IN_PROGRESS

  !! Try to open the mhd baryonic gas checkpoint

  ofile=output_path//z_s(1:len_trim(z_s))//'mhd'// &
  	rank_s(1:len_trim(rank_s))//'.dat'
#ifdef BINARY 
  open(185,file=ofile,status='old',iostat=fstat,form='binary')
#else 
  open(185,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif
  if (fstat /= 0) then
    write(*,*) 'error opening baryonic gas checkpoint'
    write(*,*) 'rank',rank,'file:',ofile
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
#endif

  !! Process the dark matter checkpoints

#ifdef PPINT
  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,sim_checkpoint, &
              sim_projection,sim_halofind,mass_p
#else
  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,sim_checkpoint, &
               sim_projection,sim_halofind,mass_p
#endif

  if (np_local > max_np) then
    write(*,*) 'too many particles to store'
    write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

!! tally up total number of particles
  npl8=int(np_local,kind=8)
  call mpi_reduce(npl8,np_total,1,mpi_double_precision, &
                       mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'number of particles =', np_total
!  if (rank == 0 .and. np_total > max_np) then
!    write(*,*) 'np total > max_np',np_total,max_np
!    call mpi_abort(mpi_comm_world,ierr,ierr)
!  endif

  call mpi_gather(np_local,1,mpi_integer,np_found,1,mpi_integer,0, &
                mpi_comm_world,ierr)

  if (rank == 0) then
    write(*,*) 'particles / node'
    do i=0,nn-1
      write(*,*) i, np_found(i)
    enddo
  endif
  read(21) xv(:,:np_local)
!  do j=1,np_local
!    read(21) xv(:,j)
!  enddo
  close(21)
!! prepare particles for accumulation 
  do j=1,np_local
    xv(1:3,j)=modulo(xv(1:3,j)+node_coords(:)*ncc,rnc)
  enddo
!! start writing
  call mpi_barrier(mpi_comm_world,ierr) !! DEBUG
  if (rank==0) then
    print *,'finished reading particles'
!! local root process writes first
    ofile=output_path//z_s(1:len_trim(z_s))//'xv.dat'
#ifdef BINARY
    open(23,file=ofile,form='binary',status='replace',iostat=fstat)
#else
    open(23,file=ofile,form='unformatted',status='replace',iostat=fstat)
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening composed checkpoint'
      write(*,*) 'rank',rank,'file:',ofile
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    np_current=sum(np_found)
    write(23) np_current
!    do i=1,np_local
!      write(23) xv(:,i)    
!    enddo
    write(23) xv(:,:np_local)
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  if (rank==0) print *,'finished initial write/modulus'
!! send particles back to master node
  
  np_current=np_local
  do i=1,nn-1
    tag=i
    if (rank == 0) then
      call mpi_recv(xv,6*np_found(i),mpi_real, &
                    i,tag,mpi_comm_world,status,ierr)
!!! ANOTHER POTENTIAL ERROR FOR FORTRAN UNFORMATTED FILES !!!
      do j=1,np_found(i)
        write(23) xv(:,j)
      enddo
    elseif (rank == i) then
      call mpi_send(xv,6*np_local,mpi_real,0,tag,mpi_comm_world,ierr)
    endif
    call mpi_barrier(mpi_comm_world,ierr)
    if (rank==0)  print *,'finished rank',i
  enddo
  if (rank==0) close(23) 

#ifdef MHD_IN_PROGRESS

  !! It's easy, send the data back and write it to disk

    read(185) cur_iter,cur_t,nx,ny,nz
    read(185) u
    read(185) b
    close(185)

#endif

enddo

end 
