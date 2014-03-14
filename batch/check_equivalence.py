import numpy

# ----------------------- parameters -----------------------

#
# These ones are usually changed
#

nodes_dim      = 8
tiles_node_dim = 16
nf_tile        = 112 

density_buffer = 1.5

# Factor to reduce max_buf by in cubepm.par 
srfac = 1

# Set this true if using P3DFFT for pencil decomposition
pencil = True

# These set the total number of threads
cores          = 8
nested_threads = 2

#
# These ones are usually not changed
#

mesh_scale  = 4
nc_halo_max = 64
max_llf     = 100000
ngrid_max   = 240
nf_cutoff   = 16
nf_buf      = nf_cutoff + 8

# ----------------------------------------------------------

#
# Number of mesh cells per dimension
#

nc = (nf_tile - 2 * nf_buf) * tiles_node_dim * nodes_dim

#
# Parameters calculated in cubepm.par
#

nodes      = nodes_dim**3
tiles_node = tiles_node_dim**3

nc_tile_dim = (nf_tile - 2 * nf_buf) / mesh_scale
nc_node_dim = nc_tile_dim * tiles_node_dim
nc_dim      = nc_node_dim * nodes_dim
nodes_slab  = nodes_dim * nodes_dim
nodes_pen   = nodes_dim
nc_slab     = nc_dim / nodes
nc_pen      = nc_node_dim / nodes_dim

max_np  = int(density_buffer * (((nf_tile - 2 * nf_buf) * tiles_node_dim / 2)**3 + \
    (8 * nf_buf**3 + 6 * nf_buf * (((nf_tile - 2 * nf_buf) * tiles_node_dim)**2) + \
    12 * (nf_buf**2) * ((nf_tile - 2 * nf_buf) * tiles_node_dim)) / 8.))
max_buf = 2 * max_np / srfac

nlist       = 5 * (nc_halo_max + 1)**3
max_maxima  = 5 * nc_halo_max**3
max_halo_np = nc**2

print
print "Sizes of various variables:"
print "nc          = ", nc
print "nodes       = ", nodes
print "tiles_node  = ", tiles_node
print "nc_tile_dim = ", nc_tile_dim
print "nc_node_dim = ", nc_node_dim
print "nc_dim      = ", nc_dim
print "nodes_slab  = ", nodes_slab
print "nodes_pen   = ", nodes_pen
print "nc_slab     = ", nc_slab
print "nc_pen      = ", nc_pen
print "max_np      = ", max_np
print "max_buf     = ", max_buf
print "nlist       = ", nlist
print "max_maxima  = ", max_maxima
print "max_halo_np = ", max_halo_np
print
print "NOTE: You have chosen srfac = ", srfac

#
# Size of equivalence arrays declared in cubepm.fh
#

print
print "Checking equivalence statements in cubepm.fh ... "
print

# First statement

isortpos  = max_halo_np
isortpeak = max_maxima
isortdist = nlist

if not (isortpos >= isortpeak and isortpos >= isortdist):

    print "Reverse equivalence between isortpos, isortdist, and isortpeak !!"
    print "isortpos  = ", isortpos
    print "isortdist = ", isortdist
    print "isortpeak = ", isortpeak
    print 

# Second statement

force_f = 3 * ((nf_tile - nf_buf + 1) - (nf_buf - 1) + 1)**3 * cores 
force_c = 3 * (nc_node_dim + 2)**3

bytes_eq2 = max(force_f, force_c)

if force_c > force_f:

    print "Reverse equivalence between force_f and force_c !!"
    print "force_f = ", force_f
    print "force_c = ", force_c
    print

# Third statement

rho_f = (nf_tile + 2) * nf_tile**2 * cores
rho_c = nc_node_dim**3

bytes_eq3 = max(rho_f, rho_c)

if rho_c > rho_f:

    print "Reverse equivalence between rho_f and rho_c !!"
    print "rho_f = ", rho_f
    print "rho_c = ", rho_c
    print

# Fourth statement

if (pencil): # Slab work is different here

    send_buf       = max_buf
    force_c_buffer = 3 * (nc_node_dim + 2)**2
    fast_buf       = max_buf
    slab_work      = nc_dim * nc_node_dim * (nc_pen + 2)

else:

    send_buf       = max_buf
    force_c_buffer = 3 * (nc_node_dim + 2)**2
    fast_buf       = max_buf
    slab_work      = (nc_dim + 2) * nc_dim * nc_slab

bytes_eq4 = max(send_buf, force_c_buffer, fast_buf, slab_work)

if not (send_buf >= force_c_buffer and send_buf >= fast_buf and send_buf >= slab_work):

    print "Reverse equivalence between send_buf, force_c_buffer, fast_buf, and slab_work !!"
    print "send_buf       = ", send_buf
    print "force_c_buffer = ", force_c_buffer
    print "fast_buf       = ", fast_buf
    print "slab_work      = ", slab_work
    print

# Fifth statement

if (pencil): # tmp_kern_c is different here

    cmplx_rho_f = (nf_tile + 2) * nf_tile**2 * cores
    tmp_kern_c  = 3 * nc_dim * nc_node_dim * (nc_pen + 2)

else:

    cmplx_rho_f = (nf_tile + 2) * nf_tile**2 * cores
    tmp_kern_c  = 3 * (nc_dim + 2) * nc_dim * nc_slab

bytes_eq5 = max(cmplx_rho_f, tmp_kern_c)

if tmp_kern_c > cmplx_rho_f:

    print "Reverse equivalence between cmplx_rho_f and tmp_kern_c !!"
    print "cmplx_rho_f = ", cmplx_rho_f
    print "tmp_kern_c  = ", tmp_kern_c
    print

# Sixth statement

if (pencil): # recv_cube is different here

    recv_buf  = max_buf
    recv_cube = nc_node_dim**2 * nc_pen * nodes_pen
    fast_pos  = max_buf / 2

else:

    recv_buf  = max_buf
    recv_cube = nc_node_dim**2 * nc_slab * nodes_slab
    fast_pos  = max_buf / 2

bytes_eq6 = max(recv_buf, recv_cube, fast_pos)

if not (recv_buf >= recv_cube and recv_buf >= fast_pos):

    print "Reverse equivalence between recv_buf, recv_cube, and fast_pos !!"
    print "recv_buf  = ", recv_buf
    print "recv_cube = ", recv_cube
    print "fast_pos  = ", fast_pos
    print

# Seventh statement

xv = 6 * max_np
ck = 3 * nc_node_dim**3

bytes_eq7 = max(xv, ck)

if ck > xv:

    print "Reverse equivalence between xv and ck !!"
    print "xv = ", xv
    print "ck = ", ck
    print

# Total memory usage from these arrays

bytes = bytes_eq2 + bytes_eq3 + bytes_eq4 + bytes_eq5 + bytes_eq6 + bytes_eq7
memGB = 4. * bytes / 1024. / 1024. / 1024.

print "Total memory usage from equivalenced .fh arrays: " + str(memGB) + " GB"
print

#
# Print size information for some big arrays
#

send_buf_PID = 8*max_buf
recv_buf_PID = 8*max_buf
PID          = 8*max_np
ll           = 4*max_np
llf          = 4*max_llf*mesh_scale*mesh_scale*mesh_scale*cores*nested_threads
pp_force_accum = 4*3*max_llf*cores*nested_threads
pp_ext_force_accum = 4*max_np*cores
pos = 4*4*max_halo_np
finegrid = 4*ngrid_max**3
ilist_odc = 4*max_halo_np
ilist_vir = 4*max_halo_np
hpart_odc = 1*max_np
hpart_vir = 1*max_np
halo_mesh_mass = 4*max_maxima

print "Sizes of some of the largest arrays [GB]:"
print "PID = ", PID / 1024. / 1024. / 1024.
print "send/recv_buf_PID = ", (send_buf_PID+recv_buf_PID) / 1024. / 1024. / 1024.
print "ll = ", ll / 1024. / 1024. / 1024.
print "llf = ", llf / 1024. / 1024. / 1024.
print "pp_force_accum = ", pp_force_accum / 1024. / 1024. / 1024.
print "pp_ext_force_accum = ", pp_ext_force_accum / 1024. / 1024. / 1024.
print "pos = ", pos / 1024. / 1024. / 1024.
print "finegrid = ", finegrid / 1024. / 1024. / 1024.
print "ilist_odc = ", ilist_odc / 1024. / 1024. / 1024.
print "ilist_vir = ", ilist_vir / 1024. / 1024. / 1024.
print "hpart_odc = ", hpart_odc / 1024. / 1024. / 1024.
print "hpart_vir = ", hpart_odc / 1024. / 1024. / 1024.
print "halo_mesh_mass = ", halo_mesh_mass / 1024. / 1024. / 1024.

#
# Size of equivalence arrays declared in dist_init_dm.f90 
#

print
print "Checking equivalence statements in dist_init_dm.f90 ... "
print

nc_node_dim = nc / nodes_dim
nc_slab     = nc / nodes
nodes_slab  = nodes_dim**2

# First statement

phi       = (nc_node_dim + 2)**3
slab_work = (nc + 2) * nc * nc_slab 
recv_cube = nc_node_dim**2 * nc_slab * nodes_slab

if not (phi >= slab_work and phi >= recv_cube):

    print "Reverse equivalence between phi, slab_work, and recv_cube !!"
    print "phi       = ", phi
    print "slab_work = ", slab_work
    print "recv_cube = ", recv_cube

# Second statement

slab = (nc + 2) * nc * nc_slab 
cube = nc_node_dim**3

if cube > slab:

    print "Reverse equivalence between slab and cube !!"
    print "slab = ", slab
    print "cube = ", cube
    print

