module Parameters
  implicit none
  public
 
  integer, parameter :: Nmpi = 1 !#Number of mpi tasks = nodes_dim**3
  integer, parameter :: Nomp = 8 !#Number of openmp threads

  integer, parameter :: Ncells = 576 !# cells in each node = nc / nodes_dim
  real, parameter :: Lbox = 1200.0 !Physical size of box in h/Mpc

  character(len=*), parameter :: dir = '/scratch/p/pen/dinman/Dipole_Simulations/m100mev/nodes64_h500/'
  character(len=*), parameter :: endian = 'little_endian'

  character(len=*), parameter :: redshift = '0.000'
  real, parameter :: z = 0.010
  real, parameter :: Om = 0.32+3*0.1/93.14/0.67**2
  
end module Parameters
