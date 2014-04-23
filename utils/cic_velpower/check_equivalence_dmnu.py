import numpy

# ----------------------- parameters -----------------------

# Parameters that we need
nodes_dim      = 2
tiles_node_dim = 6
nf_tile        = 176 
density_buffer = 1.5
ratio_nudm_dim = 2
nfine_buf      = 4
mesh_scale     = 2
nt             = 8

# Set this true if using the momentum_density/mass_density method
MOMENTUM = False 

# Set this true if using P3DFFT for pencil decomposition
pencil = True

# Usually won't need to change these
nf_cutoff   = 16
nf_buf      = nf_cutoff + 8

# ----------------------------------------------------------

#
# Compute values that we need
#

nc = (nf_tile - 2 * nf_buf) * tiles_node_dim * nodes_dim
nodes = nodes_dim**3
nc_node_dim = nc / nodes_dim
nc_pen = nc_node_dim / nodes_dim
nc_slab = nc / nodes
nodes_pen = nodes_dim
nodes_slab = nodes_dim**2

max_np  = int(density_buffer * (((nf_tile - 2 * nf_buf) * tiles_node_dim / 2)**3 + \
    (8 * nf_buf**3 + 6 * nf_buf * (((nf_tile - 2 * nf_buf) * tiles_node_dim)**2) + \
    12 * (nf_buf**2) * ((nf_tile - 2 * nf_buf) * tiles_node_dim)) / 8.))
np_buffer = int(2.*max_np/3.)

nc_buf = nfine_buf / mesh_scale
num_ngbhs = (2*nc_buf+1)**3
max_npart_search = int(num_ngbhs*float(max_np)/float(nc_node_dim)**2)
nm_node_dim = nc_node_dim / mesh_scale

print
print "Sizes of various variables:"
print "nc          = ", nc
print "nc_node_dim = ", nc_node_dim
print "nc_slab     = ", nc_slab
print "nc_pen      = ", nc_pen
print "max_np      = ", max_np
print "np_buffer   = ", np_buffer
print "num_ngbhs   = ", num_ngbhs
print "max_npart_search = ", max_npart_search

#
# Size of equivalence arrays declared in cic_mompower_p3dfft.f90 
#

print
print "Checking equivalence statements in cic_mompower_p3dfft.f90 ... "
print

# First statement (recv_cube, xp_buf) 
if pencil:
    recv_cube = 4 * (nc_node_dim**2 * nc_pen * nodes_pen) 
else:
    recv_cube = 4 * (nc_node_dim**2 * nc_slab * nodes_slab)
xp_buf = 4 * 6 * np_buffer

maxeq1 = max(recv_cube, xp_buf)

if maxeq1 != recv_cube:

    print "Change equivalence and common block accordingly:"
    print "recv_cube = ", recv_cube
    print "xp_buf    = ", xp_buf
    print

# Second statement (momden, send_buf)

momden = 4 * (1 * (nc_node_dim+2)**3)
send_buf = 4 * (6*np_buffer)

maxeq2 = max(momden, send_buf)

if maxeq2 != momden:

    print "Change equivalence and common block accordingly:"
    print "momden   = ", momden
    print "send_buf = ", send_buf
    print

# Third statement (only used if SLAB is true)

if pencil:
    slab_work = 0
else:
    slab_work = 4 * ((nc + 2) * nc * nc_slab)
recv_buf = 4 * 6 * np_buffer

maxeq3 = max(slab_work, recv_buf)

if not pencil and maxeq3 != slab_work:

    print "Change equivalence and common block accordingly:"
    print "slab_work = ", slab_work
    print "recv_buf  = ", recv_buf
    print

#
# Compute size of the other equivalence statements that won't change 
#

if pencil:
    slab = 4 * (nc * nc_node_dim * (nc_pen + 2))
else:
    slab = 4 * ((nc + 2) * nc * nc_slab)
massden_send_buff = 4 *(nc_node_dim+2)**2
massden_recv_buff = 4 *(nc_node_dim+2)**2

maxeq4 = slab 
maxeq5 = massden_send_buff
maxeq6 = massden_recv_buff

#
# Compute the sizes of other big arrays
#

xvp = 4 * (6 * max_np)
xvp_dm = xvp / ratio_nudm_dim**3
if MOMENTUM: GID = max_np
else: GID = 0
GID_dm = GID / ratio_nudm_dim**3 
if MOMENTUM: massden = 4 * (nc_node_dim+2)**3
else: massden = 0
if not MOMENTUM:
    hoc = 4 * ((nm_node_dim+2*nc_buf)**3) 
    hoc_dm = 4 * ((nm_node_dim+2*nc_buf)**3) 
    ll  = 4 * max_np
    ll_dm = ll / ratio_nudm_dim**3
    ipos = 4 * 2 * nt * max_npart_search
    rpos = 4 * 2 * nt * max_npart_search
else:
    hoc = 0
    hoc_dm = 0
    ll  = 0
    ll_dm = 0
    ipos = 0
    rpos = 0 

print "Total memory used by large arrays [GB]: ", (maxeq1 + maxeq2 + maxeq3 + maxeq4 + maxeq5 + maxeq6 + xvp + xvp_dm + massden + hoc + hoc_dm + ll + ll_dm + rpos + ipos + GID + GID_dm) / 1024.**3
print

