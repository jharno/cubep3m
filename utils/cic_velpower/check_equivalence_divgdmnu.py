import numpy

# ----------------------- parameters -----------------------

# Parameters that we need
nodes_dim      = 2
tiles_node_dim = 6
nf_tile        = 240 
density_buffer = 1.8
ratio_nudm_dim = 2
nfine_buf      = 8 
nfine_buf_h    = 64
mesh_scale     = 1
mesh_scale_h   = 8
nt             = 16
max_np_h       = 1000000

# Set this true if using P3DFFT for pencil decomposition
pencil = True

# Set this true if you are using curl
curl = False 

# Set this true if you are using momentum density field method for neutrinos
momentum = True

# Set this true if using -DCOARSE_HACK and choose appropriate value of coarsen_factor
coarse_hack = True 
coarsen_factor = 4

# Usually won't need to change these
nf_cutoff   = 16
nf_buf      = nf_cutoff + 8

# ----------------------------------------------------------

#
# Compute values that we need
#

nc = (nf_tile - 2 * nf_buf) * tiles_node_dim * nodes_dim
if coarse_hack: nc /= coarsen_factor
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

nc_buf   = nfine_buf / mesh_scale
nc_buf_h = nfine_buf_h / mesh_scale_h
num_ngbhs   = (2*nc_buf+1)**3
num_ngbhs_h = (2*nc_buf_h+1)**3
nm_node_dim   = nc_node_dim / mesh_scale
nm_node_dim_h = nc_node_dim / mesh_scale_h

#
# Make sure domain decomposition is valid
#

if pencil:
    if nc_node_dim%nodes_dim != 0:
        print "\nERROR: nc_node_dim, nodes_dim, nc_pen = ", nc_node_dim, nodes_dim, nc_pen
        exit()
else:
    if nc_dim%nodes != 0:
        print "\nERROR: nc_dim, nodes, nc_slab = ", nc_dim, nodes, nc_slab
        exit()

#
# Print some info to screen
#

print
print "Sizes of various variables:"
print "nc          = ", nc
print "nc_node_dim = ", nc_node_dim
print "nc_slab     = ", nc_slab
print "nc_pen      = ", nc_pen
print "max_np      = ", max_np
print "np_buffer   = ", np_buffer
print "num_ngbhs   = ", num_ngbhs

#
# Size of equivalence arrays declared in cic_mompower_p3dfft.f90 
#

print
print "Checking equivalence statements in cic_crossvel_dmnu.f90 ... "
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

# Second statement (velden, send_buf)
velden   = 4 * (nc_node_dim+1)**3
send_buf = 4 * 6 * np_buffer
maxeq2 = max(velden, send_buf)

if maxeq2 != velden:

    print "Change equivalence and common block accordingly:"
    print "velden   = ", velden
    print "send_buf = ", send_buf
    print

# Third statement equivalence(velden2, recv_buf)
velden2  = velden
recv_buf = send_buf

maxeq3 = max(velden2, recv_buf)
if maxeq3 != velden2:

    print "Change equivalence and common block accordingly:"
    print "velden2  = ", velden2
    print "recv_buf = ", recv_buf
    print

# Fourth statement equivalence(hoc, hoc_dm, hoc_h)
hoc    = 4 * ((nm_node_dim+2*nc_buf)**3)
hoc_dm = hoc
hoc_h  = 4 * ((nm_node_dim_h+2*nc_buf_h)**3) 

maxeq4 = max(hoc, hoc_dm, hoc_h)
if maxeq4 != hoc:

    print "Change equivalence and common block accordingly:"
    print "hoc    = ", hoc
    print "hoc_dm = ", hoc_dm
    print "hoc_h  = ", hoc_h
    print

# Fifth statement equivalence(slab, cube) (with velden3 if curl)
if pencil:
    slab = 4 * (nc * nc_node_dim * (nc_pen + 2))
else:
    slab = 4 * ((nc + 2) * nc * nc_slab)
cube = 4 * nc_node_dim**3
velden3 = 0
if curl:
    velden3 = velden

maxeq5 = max(slab, cube, velden3)

eqtest = slab
if curl: eqtest = velden3

if maxeq5 != eqtest:

    print "Change equivalence and common block accordingly:"
    print "slab = ", slab
    print "cube = ", cube
    if curl: print "velden3 = ", velden3
    print

#
# Compute the sizes of other big arrays
#

slab_work = slab
if pencil: slab_work = 0
slab2 = slab
velden4 = 0
if momentum:
    velden4 = velden

#
# Determine memory usage of large P3DFFT arrays if applicable
#

bytes_p3dfft = 0

if pencil:

    if nodes_dim == 1:
        nm = int(nc_node_dim * nc_node_dim * (nc_node_dim+2) / 2)
    else:
        nm = int(nc_node_dim * (nc_node_dim + nodes_dim) * (nc_node_dim + 2./nodes_dim) / 2.)

    buf1 = 4 * 2 * nm
    buf2 = 4 * 2 * nm
    R    = 4 * 2 * nm
    bytes_p3dfft = buf1 + buf2 + R

xvp = 4 * (6 * max_np)
xvp_dm = xvp / ratio_nudm_dim**3
xvp_h  = 4 * 6 * max_np_h 
veldivg = 4 * nc_node_dim**3
velden_send_buff = 4 * (nc_node_dim+2)**2
velden_recv_buff = velden_send_buff

t1 = maxeq1 + maxeq2 + maxeq3 + maxeq4 + maxeq5
t2 = slab_work + slab2 + xvp + xvp_dm + xvp_h + veldivg + velden_send_buff + velden_recv_buff + velden4
print "Total memory used by large arrays [GB]: ", (t1 + t2 + bytes_p3dfft) / 1024.**3
if pencil: print "Total memory of large P3DFFT arrays [GB]: ", bytes_p3dfft / 1024.**3
print

