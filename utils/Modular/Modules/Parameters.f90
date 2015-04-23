module Parameters
  implicit none
  public
  
  integer, parameter :: Nglobal_dim = 24 ! nodes_dim in runing cubep3m. set to 24 for TianNu
  integer, parameter :: Ncells = 288 !# cells in each node = nc / nodes_dim
  real, parameter :: Lbox = 100 !Physical size of box in h/Mpc
  
  integer, parameter :: Nmpi = 8 !#Number of mpi tasks = nodes_dim**3
  integer, parameter :: Nomp = 1 !#Number of openmp threads

  character(len=*), parameter :: dir = '/scratch2/p/pen/dinman/test_th2/'
  character(len=*), parameter :: endian = 'little_endian'
  
  character(len=*), parameter :: redshift = '0.000'
  real, parameter :: z = 0 ! Redshift

  real, parameter :: Om = 0.32+0.05/93.14/0.67**2

end module Parameters
