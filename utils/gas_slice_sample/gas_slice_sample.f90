!! This is the cubepm gas slice sampler.
!! Reads in distributed checkpoints and generates data slices
!! compile with: mpif77 -fpp -DBINARY gas_slice_sample.f90 -o gas_slice_sample

implicit none
include 'mpif.h'

integer, parameter :: z_lb=10
integer, parameter :: slice_width=2
integer, parameter :: z_ub=z_lb+slice_width-1

! frequently changed parameters are found in this header file:
include '../../parameters'

! list of checkpoint redshift (same as in cubepm.par)
character(len=*),parameter :: checkpoints =cubepm_root//'input/checkpoints' 

real, parameter    :: rnc = nc
integer, parameter :: max_input=100 !! maximum number of checkpoints
integer, parameter :: nn=nodes_dim**3  !! number of nodes total
real, parameter    :: ncc=nc/nodes_dim !! number of cells / cubic 
integer, parameter :: hc=nc/2       !! half-grid length
real, dimension(5,nc,nc,slice_width) :: gas_slice,gas_slice_local
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
real :: u(5,ncc,ncc,ncc),b(3,ncc,ncc,ncc) !! gas and magnetic fields

integer :: node_coords(3)           !! coordinates of node in cube	
character(len=80) ofile 
character(len=4) rank_s
character(len=7) z_s
integer i,j,k,m,i0,j0
real x,y,z,z_write
integer(4) :: nn_returned, ierr,rank, num_checkpoints,fstat,cur_checkpoint
real(4), dimension(max_input) :: z_checkpoint

integer(4) :: nts, sim_checkpoint, sim_projection, sim_halofind
integer(4) :: nodes_returned,tag
real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p
integer*8 :: npl8,np_total

integer(4), dimension(mpi_status_size) :: status

common /rarr/ u,b,gas_slice,gas_slice_local

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

!! Read in gas checkpoint

  read(185) cur_iter,cur_t,nx,ny,nz
  read(185) u
  read(185) b
  close(185)

!! Build slices

  gas_slice=0.0
  gas_slice_local=0.0
  if (node_coords(3)==0) then
    i0=node_coords(1)*ncc
    j0=node_coords(2)*ncc
    gas_slice_local(:,1+i0:ncc+i0,1+j0:ncc+j0,:)=u(:,:,:,z_lb:z_ub)
  endif   

  call mpi_reduce(gas_slice_local,gas_slice,5*nc*nc*2,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

  if (rank==0) then
    ofile=output_path//z_s(1:len_trim(z_s))//'gas_slice.dat'

    open(152,file=ofile,form='formatted')
    do j=1,nc
	write(152,'(160f16.8)') gas_slice(1,:,j,1)
    enddo
    close(152)
  endif

enddo

end 
