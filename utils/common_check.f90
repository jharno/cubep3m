!check common block usage and estimate memory usage
implicit none
include '../source/cubepm.par'
integer :: s_rho_f,s_rho_ft,s_kern_f,s_force_f,s_ck,s_kern_c,s_rho_c,s_force_c
integer :: s_complx_rho_c,s_force_c_buffer,s_slab,s_slab_work,s_recv_cube,s_send_buf
integer :: s_recv_buf,s_fast_buf,s_fast_pos,s_xv,s_ll,s_hoc,s_toc

s_rho_f=(nf_tile+2)*nf_tile*nf_tile
s_rho_ft=2*nf_tile*(nf_tile/2+1)*nf_tile
s_kern_f=nf_tile*(nf_tile/2+1)*nf_tile
s_force_f=3*(nf_tile-nf_buf+1-nf_buf+1)**3
s_ck=3*(nc_dim+2)*nc_dim*nc_dim
s_kern_c=3*(nc_dim/2+1)*nc_dim*nc_slab
s_rho_c=nc_node_dim**3
s_force_c=3*(nc_node_dim+2)**3
s_complx_rho_c=(nc_dim+2)*nc_dim*nc_slab
s_force_c_buffer=3*(nc_node_dim+2)**2
s_slab=nc_dim*(nc_dim+2)*nc_slab
s_slab_work=nc_dim*(nc_dim+2)*nc_slab
s_recv_cube=nc_node_dim*nc_node_dim*nc_slab*nodes_slab
s_send_buf=max_buf
s_recv_buf=max_buf
s_fast_buf=max_buf
s_fast_pos=max_buf/2
s_xv=6*max_np
s_ll=max_np
s_hoc=(hoc_nc_h-hoc_nc_l+1)**3

print *,'rho_f',s_rho_f
print *,'rho_ft',s_rho_ft
print *,'kern_f',s_kern_f
print *,'force_f',s_force_f
print *,'ck',s_ck
print *,'kern_c',s_kern_c
print *,'rho_c',s_rho_c
print *,'force_c',s_force_c
print *,'complx_rho_c',s_complx_rho_c
print *,'force_c_buffer',s_force_c_buffer
print *,'slab',s_slab
print *,'slab_work',s_slab_work
print *,'recv_cube',s_recv_cube
print *,'send_buf',s_send_buf
print *,'recv_buf',s_recv_buf
print *,'fast_buf',s_fast_buf
print *,'fast_pos',s_fast_pos
print *,'xv',s_xv
print *,'ll',s_ll
print *,'hoc',s_hoc
end
