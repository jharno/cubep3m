  !! Cosmo parameters
  real, parameter :: ns=1.0
  real, parameter :: s8=0.9
  real, parameter :: omegam=0.25
  real, parameter :: omegal=0.75
  real, parameter :: omegab=0.04

  real, parameter :: box=100 !182.5 !140.0
  real, parameter :: zi=200

  !! dir is the directory for files
  character(*), parameter :: dir='/scratch/merz/cubep3m/'

  !! nk is the length of the initial power spectrum file
  integer, parameter      :: nk=408
  character(*), parameter :: fntf='~/data/cmbfast/Ilian_run/cmbfast.lcdm'
 
  !! nt is the number of threads
  integer, parameter      :: nt=1

  !! nout is the number of data outputs
  integer, parameter :: nout=300

  !! nc is the number of cells per box length
  integer, parameter :: nc=200 !192 !256

  !! proj_nc is the pixel size of the projection
  integer, parameter :: proj_nc=512

  !! slab width 
  real, parameter :: slab_width=1.0/8.0 !(1/6 at box=100.0)

  !!maximum expansion rate 
  real, parameter :: ramax=0.01645
