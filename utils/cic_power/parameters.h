  !! Cosmo parameters
  real, parameter :: ns=1.0
  real, parameter :: s8=0.9
  real, parameter :: omegam=0.27
  real, parameter :: omegal=0.73
  real, parameter :: omegab=0.044

  real, parameter :: box=100.0
  real, parameter :: zi=140

  !! dir is the directory for files
  character(*), parameter :: dir='/scratch/merz/cubepm_p3m_scaled/'

  !! nk is the length of the initial power spectrum file
  integer, parameter      :: nk=408
  character(*), parameter :: fntf='~/data/cmbfast/Ilian_run/cmbfast.lcdm'
 
  !! nt is the number of threads
  integer, parameter      :: nt=1

  !! nout is the number of data outputs
  integer, parameter :: nout=300

  !! nc is the number of cells per box length
  integer, parameter :: nc=80
