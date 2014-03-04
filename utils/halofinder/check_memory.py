# ----------------------- parameters -----------------------

#
# These ones are usually changed
#

nodes_dim      = 12
tiles_node_dim = 16
nf_tile        = 120

ngrid_max   = 220 
nc_halo_max = 64

hvir = True

#
# These ones are usually not changed
#

cores     = 8 
nf_cutoff = 16
nf_buf    = nf_cutoff + 8

mesh_scale    = 4
np_buffer_fac = 4

# ----------------------------------------------------------

#
# Number of mesh cells per dimension
#

nc = (nf_tile - 2 * nf_buf) * tiles_node_dim * nodes_dim

#
# Parameters used in halofind.f90
#

nodes       = nodes_dim**3
tiles_node  = tiles_node_dim**3
np          = nc / 2
nc_node_dim = nc / nodes_dim
np_node_dim = np / nodes_dim
np_buffer   = np_node_dim**3 / np_buffer_fac
max_np      = np_node_dim**3 + np_buffer
nc_buf      = nf_buf / mesh_scale
nc_tile_dim = (nf_tile - 2 * nf_buf) / mesh_scale
hoc_nc_l    = 1 - nc_buf
hoc_nc_h    = nc_tile_dim*tiles_node_dim + nc_buf
max_halo_np = nc**2

#
# Array sizes
#

xvp      = 6 * max_np
xvp_buf  = 6 * np_buffer
send_buf = 6 * np_buffer
recv_buf = 6 * np_buffer

print "Particle arrays: "
print "xvp     = ", 4. * xvp / 1024.**3 
print "xvp_buf = ", 4. * xvp_buf / 1024.**3
print "send_buf + recv_buf = ", 4. * (send_buf + recv_buf) / 1024.**3 
print "TOTAL = ", (xvp + xvp_buf + send_buf + recv_buf) * 4. / 1024.**3
print

ll  = max_np
hoc = (hoc_nc_h - hoc_nc_l + 1)**3

print "Linked list arrays: "
print "ll = ", 4. * ll / 1024.**3
print "hoc = ", 4. * hoc / 1024.**3
print "TOTAL = ", (ll + hoc) * 4. / 1024.**3
print

rho_f = (nf_tile + 2) * nf_tile**2

print "rho_f: ", rho_f * 4. / 1024.**3
print

max_maxima = 5 * nc_halo_max**3
nlist      = 5 * (nc_halo_max+1)**3

isortpeak = max_maxima
isortdist = nlist
ipeak     = 3 * max_maxima
idist     = 3 * nlist
den_peak  = max_maxima
rdist     = nlist
hpart_odc = max_np
ilist_odc = max_halo_np
halo_mesh_mass = max_maxima
finegrid = ngrid_max**3

if hvir:

    hpart_vir = max_np
    ilist_vir = max_halo_np

else:

    hpart_vir = 0
    ilist_vir = 0

print "Halo arrays: "
print "isortpeak = ", isortpeak * 4. / 1024.**3
print "isortdist = ", isortdist * 4. / 1024.**3
print "ipeak     = ", ipeak * 4. / 1024.**3
print "idist     = ", idist * 4. / 1024.**3
print "den_peak  = ", den_peak * 4. / 1024.**3
print "rdist     = ", rdist * 4. / 1024.**3
print "hpart_odc = ", hpart_odc / 1024.**3
print "ilist_odc = ", ilist_odc * 4. / 1024.**3
print "hpart_vir = ", hpart_vir / 1024.**3
print "ilist_vir = ", ilist_vir * 4. / 1024.**3
print "halo_mesh_mass = ", halo_mesh_mass * 4. / 1024.**3
print "finegrid = ", finegrid * 4. / 1024.**3
print "TOTAL = ", (isortpeak + isortdist + ipeak + idist + den_peak + rdist + hpart_odc/4. + ilist_odc + hpart_vir/4. + ilist_vir + halo_mesh_mass + finegrid) * 4. / 1024.**3
print

