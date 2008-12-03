!! cic_power.f90 Parallelized: Hugh Merz Jun 15, 2005
!! updated to use less memory - lots of equivalencing: Jun 12, 2007
!! mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT cic_power_mem_check.f90 -o mem_check -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
program cic_power 
  implicit none
  include 'mpif.h'

! frequently changed parameters are found in this header file:
  include '../../parameters'

  logical, parameter :: correct_kernel=.false.

  character(len=*), parameter :: checkpoints=cubepm_root//'/input/checkpoints'

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc
  real, parameter    :: npr=np

  !! internals
  integer, parameter :: max_checkpoints=100
  real, dimension(max_checkpoints) :: z_checkpoint
  integer num_checkpoints, cur_checkpoint

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
!  integer(4), parameter :: np_node_dim = np/nodes_dim
!  integer(4), parameter :: np_buffer = 5*np_node_dim**3
!  integer(4), parameter :: np_buffer = 0.35*np_node_dim**3

!  integer(4), parameter :: np_buffer = 4*np_node_dim**3
!  integer(4), parameter :: max_np = np_node_dim**3 + np_buffer

  integer(4), parameter :: max_np = density_buffer * ( ((nf_tile-2*nf_buf)*tiles_node_dim/2)**3 + &
                                  (8*nf_buf**3 + 6*nf_buf*(((nf_tile-2*nf_buf)*tiles_node_dim)**2) + &
                                  12*(nf_buf**2)*((nf_tile-2*nf_buf)*tiles_node_dim))/8.0 )
  integer(4), parameter :: np_buffer=int(2./3.*max_np)

  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! parallelization variables
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local

  integer(8) :: plan, iplan

  logical :: firstfftw

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.14159

  !! Dark matter arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(3,np_buffer) :: xp_buf
  real, dimension(3*np_buffer) :: send_buf, recv_buf
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: den 
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: den_buf 

  !! Power spectrum arrays
  real, dimension(2,nc) :: pkdm
#ifdef PLPLOT
  real*8, dimension(3,nc) :: pkplot
#endif

  !! Fourier transform arrays
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
  real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work

  !! Equivalence arrays to save memory
  equivalence (den,slab_work,recv_cube,xp_buf) 
  equivalence (xvp,slab,cube)  !! merz --  not sure if xvp is larger than slab?????

  !! Common block
#ifdef PLPLOT
!  common xvp,send_buf,slab_work,den_buf,den,cube,slab,xp_buf,recv_buf,pkdm,pkplot
  common xvp,send_buf,den_buf,den,recv_buf,pkdm,pkplot
#else
!  common xvp,send_buf,slab_work,den_buf,den,cube,slab,xp_buf,recv_buf,pkdm
  common xvp,send_buf,den_buf,den,recv_buf,pkdm
#endif

  real(4),parameter::mb=1./1024./1024.

  print *,'size xvp',4.*6.*max_np*mb
  print *,'size send_buf,recv_buf',2.*4.*3.*np_buffer*mb
  print *,'size den_buf',mb*4.*(nc_node_dim+2)**2
  print *,'size den',mb*4.*real(nc_node_dim+2,4)**3
  print *,'size pk',mb*5.0*4.*nc
  print *,'size slab (opt)',mb*4.*(nc+2)*nc*nc_slab
  
  end

