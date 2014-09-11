import numpy

# ----------------------- parameters -----------------------

# Parameters that we need
nodes_dim      = 4
tiles_node_dim = 4
nf_tile        = 176 
density_buffer = 2.0
ratio_nudm_dim = 1
max_np_h       = 1000000 

# Set this ture if using GROUPS algorithm
groups = True

# Set this true if using HALOMASS algorithm
halomass = False 

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

print
print "Sizes of various variables:"
print "nc          = ", nc
print "nc_node_dim = ", nc_node_dim
print "nc_slab     = ", nc_slab
print "nc_pen      = ", nc_pen
print "max_np      = ", max_np
print "np_buffer   = ", np_buffer

#
# Size of equivalence arrays declared in cic_power_dmnu.f90 
#

print
print "Checking equivalence statements in cic_power_dmnu.f90 ... "
print

# First statement (den, slab_work, recv_cube, xp_buf)
den = 4 * (nc_node_dim+2)**3
if pencil:
    recv_cube = 4 * (nc_node_dim**2 * nc_pen * nodes_pen) 
    slab_work = 0
else:
    recv_cube = 4 * (nc_node_dim**2 * nc_slab * nodes_slab)
    slab_work = 4 * ((nc + 2) * nc * nc_slab)
if halomass:
    xp_buf = 4 * 4 * np_buffer
else:
    xp_buf = 3 * 4 * np_buffer

maxeq1 = max(den, slab_work, recv_cube, xp_buf)

if maxeq1 != den: 

    print "Change equivalence and common block accordingly:"
    print "den       = ", den
    if not pencil: print "slab_work = ", slab_work
    print "recv_cube = ", recv_cube
    print "xp_buf    = ", xp_buf
    print

# Second statement (slab, cube) 
xvp = 4 * (3 * max_np) 
cube = 4 * nc_node_dim**3
if pencil:
    slab = 4 * nc * nc_node_dim * (nc_pen+2)
else:
    slab = 4 * ((nc + 2) * nc * nc_slab)

maxeq2 = max(slab, cube)

if maxeq2 != slab:

    print "Change equivalence and common block accordingly:"
    print "slab = ", slab
    print "cube = ", cube
    print

# Third statement (slab2, send_buf)
slab2 = slab
if halomass:
    send_buf = 4 * 4 * np_buffer
else:
    send_buf = 4 * 3 * np_buffer

maxeq3 = max(slab2, send_buf)

if maxeq3 != slab2:

    print "Change equivalence and common block accordingly:"
    print "slab2    = ", slab2
    print "send_buf = ", send_buf
    print

#
# Compute the sizes of other big arrays
#

xvp_dm = xvp / ratio_nudm_dim**3
if halomass: xvp_h  = 4 * 4 * max_np_h
else: xvp_h  = 4 * 3 * max_np_h
if groups : GID = max_np
else: GID = 0
recv_buf = send_buf
den_buf = 4 * (nc_node_dim+2)**2

print "Total memory used by large arrays [GB]: ", (maxeq1 + maxeq2 + maxeq3 + xvp + xvp_dm + xvp_h + recv_buf + GID + den_buf) / 1024.**3 
print

