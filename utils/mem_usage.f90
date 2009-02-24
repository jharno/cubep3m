implicit none
integer, parameter :: nf_buf=48
integer :: nf,tiles,nodes_dim,nc_slab,nc,nodes,max_np,cores
real, parameter :: mb=1024.**2
real :: density_buffer
print *,'enter nodes_dim'
read(5,*) nodes_dim
nodes=nodes_dim**3
print *,'enter cores / node'
read(5,*) cores
print *,'enter cells / coarse slab'
read(5,*) nc_slab
print *,'enter number of fine mesh tiles / node'
read(5,*) tiles
print *,'enter fraction of density to allow for buffer (1.0 == no buffer, 2.0 == 2x avg density, etc)'
read(5,*) density_buffer

nf=4*nc_slab*nodes_dim**2/tiles+nf_buf
nc=((nf-nf_buf)*tiles*nodes_dim)/4
if (mod(nc,nodes)/=0) print *,'WARNING !!!!!!!!!!!!, cannot make ICs',real(nc)/nodes
if (mod(4*nc_slab*nodes_dim**2,tiles)/=0) print *,'WARNING !!!!!!!!!!!!, cannot fit tiles into mesh evenly',real(4*nc_slab*nodes_dim**2)/tiles
if (nf<2*nf_buf) then
   print *,'WARNING !!!!!!!!!, number of fine cells / tile is too small'
   print *,' should be at least 4 times as large as the fine mesh buffer ( =',2*nf_buf,' )'
   print *,'increase cells / coarse slab or decrease fine mesh tiles / node'
endif
if ( (nf+2)*nf*nf*cores < 3*(nc+2)*nc*nc_slab ) then
  print *,'reverse equivalence for cmplx_rho_f and tmp_kern_c!'
endif


print *,'number fine cells / tile =',nf 
print *,'total simulation volume =',(nf-nf_buf)*tiles*nodes_dim
print *,'coarse mesh volume =',nc
print *,'total number of particles =',(real(nf-nf_buf)*tiles*nodes_dim/2.)**3
max_np = density_buffer * ( ((nf-nf_buf)*tiles/2)**3 + &
         (8*(nf_buf/2)**3 + 6*(nf_buf/2)*(((nf-nf_buf)*tiles)**2) + 12*((nf_buf/2)**2)*((nf-nf_buf)*tiles))/8.0 )
print *,'maximum particles in particle list =',max_np 
print *,'maximum buffer particles =', 2*max_np/6

print *,'size of fine mesh density=',8.0*cores*(nf+2)*nf*nf/mb
print *,'size of fine mesh kernel=',2.0*(nf+2)*nf*nf/mb 
print *,'size of fine mesh force=',12.0*cores*(nf-nf_buf+2)**3/mb
!print *,'size of max_np=', 3.* ((nf * tiles) / 2. )**3 

print *,'size of buffer particle list=',(4.*2.*2.*max_np)/mb
print *,'size of particle list=',(28.*max_np)/mb
print *,'size of PID list=',(8.*max_np)/mb

print *,'size of density projections=',((nf-nf_buf)*tiles*nodes_dim)**2*3.*4./mb
print *,'size of gas and magnetic field=', 32.0*real((nf-nf_buf)*tiles+6)**3/mb

if ( (nf+2)*nf*nf*cores < 3*(nc+2)*nc*nc_slab ) then
  print *,'reverse equivalence for cmplx_rho_f and tmp_kern_c!'
endif


end
