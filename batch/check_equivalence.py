import numpy

# ---------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------

#
# These ones are usually changed
#

nodes_dim      = 2
tiles_node_dim = 6
nf_tile        = 240
density_buffer = 1.5

# Factor to reduce max_buf by in cubepm.par 
srfac = 1

# Set this true if using P3DFFT for pencil decomposition
pencil = True
pencil_ICs = True
slab_fftw3_ICs = False 

# Set this true if using ZIP checkpointing method
zip = True

# Set this true if using neutrinos
neutrinos = True
ratio_nudm_dim = 2

# Set this true if not using extended pp
no_extpp = True 

# Set this true if not using projections
no_proj = True

# Set this true if not using PIDs (ignored if neutrinos is True)
no_pid = True 

# Set this true if you are turning off openmp in ICs
no_openmp_ICs = False

# These set the total number of threads
cores          = 8
nested_threads = 2

#
# These ones are usually not changed
#

mesh_scale  = 4
part_ratio  = 2
nc_halo_max = 128
max_llf     = 100000
ngrid_max   = 650 
nf_cutoff   = 16
nf_buf      = nf_cutoff + 8
pp_range    = 2
nlist       = 5 * (nc_halo_max + 1)**3
max_maxima  = 5 * nc_halo_max**3
max_halo_np = 5 * (nc_halo_max + 1)**3

dbuf_ICs       = 1.4
mesh_scale_ICs = 4

# ---------------------------------------------------------------------
# MAIN CUBEP3M CODE
# ---------------------------------------------------------------------

#
# Simulation parameters derived from those above 
#

nc          = (nf_tile - 2 * nf_buf) * tiles_node_dim * nodes_dim
nodes       = nodes_dim**3
tiles_node  = tiles_node_dim**3
nc_tile_dim = (nf_tile - 2 * nf_buf) / mesh_scale
nc_node_dim = nc_tile_dim * tiles_node_dim
nc_dim      = nc_node_dim * nodes_dim
nodes_slab  = nodes_dim * nodes_dim
nodes_pen   = nodes_dim
nc_slab     = nc_dim / nodes
nc_pen      = nc_node_dim / nodes_dim
nc_buf      = nf_buf / mesh_scale
hoc_nc_l    = 1 - nc_buf
hoc_nc_h    = nc_node_dim + nc_buf
nf_physical_tile_dim = nf_tile - 2 * nf_buf
nf_physical_dim      = nf_physical_tile_dim * tiles_node_dim * nodes_dim 
max_np  = int(density_buffer * (((nf_tile - 2 * nf_buf) * tiles_node_dim / 2)**3 + \
    (8 * nf_buf**3 + 6 * nf_buf * (((nf_tile - 2 * nf_buf) * tiles_node_dim)**2) + \
    12 * (nf_buf**2) * ((nf_tile - 2 * nf_buf) * tiles_node_dim)) / 8.))
max_buf = 2 * max_np / srfac

#
# Consistency check
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

# ---------------------------------------------------------------------
# MAIN CUBEP3M CODE
# ---------------------------------------------------------------------


print "\n---------------------------------------------------------------------"
print "MAIN CUBEP3M CODE"
print "---------------------------------------------------------------------\n"
print "VARIOUS SIMULATION PARAMETERS" 
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

#
# Determine size of large arrays declared in cubepm.fh 
#

rho_f               = 4 * ((nf_tile+2) * nf_tile**2 * cores)
cmplx_rho_f         = 4 * ((nf_tile+2) * nf_tile**2 * cores)
kern_f              = 4 * (3 * (nf_tile/2+1) * nf_tile**2)
force_f             = 4 * (3 * ((nf_tile-nf_buf+1)-(nf_buf-1)+1)**3 * cores)
llf                 = 4 * (max_llf * mesh_scale**3 * cores * nested_threads) 
pp_force_accum      = 4 * (3 * max_llf * cores * nested_threads)
pp_ext_force_accum  = 4 * 3 * (max_np * cores)
ll_fine             = 4 * (max_np * cores)
hoc_fine            = 4 * ((nf_physical_tile_dim+2*pp_range)**3 * cores) 
ck                  = 4 * (3 * nc_node_dim**3)
kern_c              = 4 * (3 * (nc_dim/2+1) * nc_dim * nc_slab)
tmp_kern_c          = 4 * (3 * (nc_dim+2) * nc_dim * nc_slab)
rho_c               = 4 * nc_node_dim**3
force_c             = 4 * (3 * (nc_node_dim+2)**3)
cmplx_rho_c         = 4 * ((nc_dim+2) * nc_dim * nc_slab)
force_c_buffer      = 4 * (3 * (nc_node_dim+2)**2)
slab                = 4 * ((nc_dim+2) * nc_dim * nc_slab)
slab_work           = 4 * ((nc_dim+2) * nc_dim * nc_slab)
recv_cube           = 4 * (nc_node_dim**2 * nc_slab * nodes_slab)
send_buf            = 4 * max_buf
recv_buf            = 4 * max_buf
fast_buf            = 4 * max_buf
fast_pos            = 4 * (max_buf / 2)
PID                 = 8 * max_np
send_buf_PID        = 8 * max_buf
recv_buf_PID        = 8 * max_buf
xv                  = 4 * (6 * max_np)
ll                  = 4 * max_np
hoc                 = 4 * (hoc_nc_h - hoc_nc_l + 1)**3
rp_buf              = 4 * nf_physical_dim**2
rho_pxy             = 4 * nf_physical_dim**2
rho_pxz             = 4 * nf_physical_dim**2
rho_pyz             = 4 * nf_physical_dim**2
isortpos            = 4 * max_halo_np
isortpeak           = 4 * max_maxima
isortdist           = 4 * nlist
idist               = 4 * (3 * nlist)
ipeak               = 4 * (3 * max_maxima) 
den_peak            = 4 * max_maxima
pos                 = 4 * (4 * max_halo_np)
rdist               = 4 * nlist
halo_mesh_mass      = 4 * max_maxima
finegrid            = 4 * ngrid_max**3
ilist_odc           = 4 * max_halo_np
ilist_vir           = 4 * max_halo_np
hpart_odc           = 1 * max_np
hpart_vir           = 1 * max_np

# Changes if P3DFFT is used instead of slab decomposition on the coarse mesh.
if pencil:
    kern_c              = 4 * (3 * nc_dim/2 * nc_node_dim * (nc_pen+2))
    tmp_kern_c          = 4 * (3 * nc_dim * nc_node_dim * (nc_pen+2))
    cmplx_rho_c         = 4 * (nc_dim * nc_node_dim * (nc_pen+2))
    slab                = 4 * (nc_dim * nc_node_dim * (nc_pen+2)) 
    slab_work           = 4 * (nc_dim * nc_node_dim * (nc_pen+2))
    recv_cube           = 4 * (nc_node_dim**2 * nc_pen * nodes_pen)

# Changes if extended pp is not being used.
if no_extpp: 
    pp_ext_force_accum  = 0
    ll_fine             = 0
    hoc_fine            = 0

# Changes if projections are turned off
if no_proj:
    rp_buf              = 0
    rho_pxy             = 0
    rho_pxz             = 0
    rho_pyz             = 0

# Changes if neutrinos are being used.
if neutrinos:
    PID                 = 1 * max_np
    send_buf_PID        = 1 * max_buf
    recv_buf_PID        = 1 * max_buf

# Changes if PIDs are not used
if no_pid and not neutrinos:
    PID                 = 0
    send_buf_PID        = 0 
    recv_buf_PID        = 0 

#
# Determine memory usage of equivalenced arrays
#

bytes_eq1 = max(isortpos, isortpeak, isortdist)
bytes_eq2 = max(force_f, force_c)
bytes_eq3 = max(rho_f, rho_c)
bytes_eq4 = max(send_buf, force_c_buffer, fast_buf, slab_work)
bytes_eq5 = max(cmplx_rho_f, tmp_kern_c)
bytes_eq6 = max(recv_buf, recv_cube, fast_pos)
bytes_eq7 = max(xv, ck)

#
# Print memory stats to screen
#

slen = 60
print "MEMORY USAGE OF LARGE ARRAYS IN CUBEPM.FH".ljust(slen), "[GB]"
print "equiv(isortpos, isortpeak, isortdist) ".ljust(slen), bytes_eq1/1024.**3
print "equiv(force_f, force_c)".ljust(slen), bytes_eq2/1024.**3
print "equiv(rho_f, rho_c)".ljust(slen), bytes_eq3/1024.**3
print "equiv(send_buf, force_c_buffer, fast_buf, slab_work)".ljust(slen), bytes_eq4/1024.**3
print "equiv(cmplx_rho_f, tmp_kern_c)".ljust(slen), bytes_eq5/1024.**3
print "equiv(recv_buf, recv_cube, fast_pos)".ljust(slen), bytes_eq6/1024.**3
print "equiv(xv, ck)".ljust(slen), bytes_eq7/1024.**3
print "kern_f".ljust(slen), kern_f/1024.**3 
print "kern_c".ljust(slen), kern_c/1024.**3
print "cmplx_rho_c".ljust(slen), cmplx_rho_c/1024.**3
print "slab".ljust(slen), slab/1024.**3
print "PID".ljust(slen), PID/1024.**3
print "send_buf_PID+recv_buf_PID".ljust(slen), (send_buf_PID+recv_buf_PID)/1024.**3 
print "llf".ljust(slen), llf/1024.**3
print "ll".ljust(slen), ll/1024.**3
print "hoc".ljust(slen), hoc/1024.**3
print "pos".ljust(slen), pos/1024.**3
print "finegrid".ljust(slen), finegrid/1024.**3
print "idist+rdist".ljust(slen), (idist+rdist)/1024.**3
print "halo_mesh_mass+ipeak+den_peak".ljust(slen), (halo_mesh_mass+ipeak+den_peak)/1024.**3
print "ilist_odc+ilist_vir".ljust(slen), (ilist_odc+ilist_vir)/1024.**3
print "hpart_odc+hpart_vir".ljust(slen), (hpart_odc+hpart_vir)/1024.**3
print "rp_buf+rho_pxy+rho_pxz+rho_pyz".ljust(slen), (rp_buf+rho_pxy+rho_pxz+rho_pyz)/1024.**3
print "pp_force_accum".ljust(slen), pp_force_accum/1024.**3
print "pp_ext_force_accum".ljust(slen), pp_ext_force_accum/1024.**3
print "ll_fine".ljust(slen), ll_fine/1024.**3
print "hoc_fine".ljust(slen), hoc_fine/1024.**3
print 

totalbytes = bytes_eq1 + bytes_eq2 + bytes_eq3 + bytes_eq4 + bytes_eq5 + bytes_eq6 + bytes_eq7 + \
    kern_f + kern_c + cmplx_rho_c + slab + PID + send_buf_PID + recv_buf_PID + llf + ll + hoc + \
    pos + finegrid + idist + rdist + halo_mesh_mass + ipeak + den_peak + ilist_odc + ilist_vir + \
    hpart_odc + hpart_vir + rp_buf + rho_pxy + rho_pxz + rho_pyz + pp_force_accum + pp_ext_force_accum + \
    ll_fine + hoc_fine

totalparts = (nc_node_dim * mesh_scale / part_ratio)**3 
if neutrinos: totalparts += totalparts / ratio_nudm_dim**3

bytesperpart = float(totalbytes) / totalparts

print "TOTAL MEMORY USAGE: " + str(totalbytes/1024.**3) + " GB"
print "PARTICLES/NODE:     " + str(totalparts)
print "TOTAL PARTICLES:    " + str(totalparts*nodes)
print "BYTES/PARTICLE:     " + str(bytesperpart)
print

#
# Print notes to screen
#

if srfac != 1:
    print " *** Make sure to change max_buf = 2 * max_np / srfac in cubepm.par *** " 
    print "     srfac = ", srfac
    print

if no_extpp:
    print " *** Make sure to comment out pp_ext_force_accum, ll_fine, and hoc_fine in cubepm.fh ***"
    print 

if no_proj:
    print " *** Make sure to comment out rp_buf, rho_pxy, rho_pxz, and rho_pyz in cubepm.fh ***"
    print

if bytes_eq1 != isortpos:
    print " *** Adjust common block for isortpos, isortdist, and isortpeak *** "
    print "     isortpos  = ", isortpos
    print "     isortdist = ", isortdist
    print "     isortpeak = ", isortpeak
    print

if bytes_eq2 != force_f:
    print " *** Adjust common block for force_f and force_c *** "
    print "     force_f = ", force_f
    print "     force_c = ", force_c
    print

if bytes_eq3 != rho_f:
    print " *** Adjust common block for rho_f and rho_c *** "
    print "     rho_f = ", rho_f
    print "     rho_c = ", rho_c
    print

if bytes_eq4 != send_buf:
    print " *** Adjust common block for send_buf, force_c_buffer, fast_buf, and slab_work *** "
    print "     send_buf       = ", send_buf
    print "     force_c_buffer = ", force_c_buffer
    print "     fast_buf       = ", fast_buf
    print "     slab_work      = ", slab_work
    print

if bytes_eq5 != cmplx_rho_f:
    print " *** Adjust common block for cmplx_rho_f and tmp_kern_c *** "
    print "     cmplx_rho_f = ", cmplx_rho_f
    print "     tmp_kern_c  = ", tmp_kern_c
    print

if bytes_eq6 != recv_buf:
    print " *** Adjust common block for recv_buf, recv_cube, and fast_pos *** "
    print "     recv_buf  = ", recv_buf
    print "     recv_cube = ", recv_cube
    print "     fast_pos  = ", fast_pos
    print

if bytes_eq7 != xv:
    print " *** Adjust common block for xv and ck *** "
    print "     xv = ", xv
    print "     ck = ", ck
    print


# ---------------------------------------------------------------------
# INITIAL CONDITIONS
# ---------------------------------------------------------------------

print "\n---------------------------------------------------------------------"
print "INITIAL CONDITIONS" 
print "---------------------------------------------------------------------\n"

#
# Readjust parameters since ICs use the fine mesh instead of the coarse mesh.
#

nc_node_dim = nc / nodes_dim
nc_pen = nc / nodes_dim**2
nc_slab     = nc / nodes
nodes_slab  = nodes_dim**2
nodes_pen   = nodes_dim
np_node_dim = nc_node_dim / 2
num_threads_ic = cores * nested_threads 
if no_openmp_ICs: num_threads_ic = 1

#
# Determine size of largest arrays
# 

cube        = 4 * nc_node_dim**3
recv_cube   = 4 * nc_node_dim**2 * nc_slab * nodes_slab
slab        = 4 * (nc + 2) * nc * nc_slab
slab_work   = 4 * (nc + 2) * nc * nc_slab
phi         = 4 * (nc_node_dim + 2)**3
phi_buf     = 4 * (nc_node_dim + 2)**2
slab_cmplx  = 0

# Changes if P3DFFT is used instead of slab decomposition on the fine mesh.
if pencil_ICs:
    recv_cube   = 4 * nc_node_dim**2 * nc_pen * nodes_pen
    slab        = 4 * nc * nc_node_dim * (nc_pen+2)
    slab_work   = 4 * nc * nc_node_dim * (nc_pen+2)

# Consistency check
if pencil_ICs and slab_fftw3_ICs:
    print "You can only pick one of pencil_ICs/slab_fftw3_ICs !!"
    exit()

if slab_fftw3_ICs: 
    slab_cmplx = 2 * slab

if neutrinos and zip:

    #
    # Sizes of some arrays are different depending on ratio_nudm_dim 
    #
    
    np_node_dim_dm = np_node_dim / ratio_nudm_dim
    nm_node_dim = nc_node_dim / mesh_scale_ICs
    max_np    = int(dbuf_ICs*np_node_dim)*np_node_dim**2
    max_np_dm = int(dbuf_ICs*np_node_dim_dm)*np_node_dim_dm**2    
    np_buffer    = max_np - np_node_dim**3
    np_buffer_dm = max_np_dm - np_node_dim_dm**3

    xvp    = 4*6*max_np
    xvp_dm = 4*6*max_np_dm
    xp_buf    = 4*6*np_buffer
    xp_buf_dm = 4*6*np_buffer_dm
    send_buf = xp_buf ; recv_buf = xp_buf
    send_buf_dm = xp_buf_dm ; recv_buf_dm = xp_buf_dm
    hoc = 4 * nm_node_dim * nm_node_dim * nm_node_dim
    ll  = 4 * max_np
    ll_dm = 4 * max_np_dm

    slab_cache = slab

    #
    # Determine memory usage of equivalenced and other arrays
    #

    bytes_eq1    = max(phi, slab_work, recv_cube, send_buf)
    bytes_eq1_dm = max(phi, slab_work, recv_cube, send_buf_dm)
    if not slab_fftw3_ICs:
        bytes_eq2    = max(slab, cube, xp_buf, hoc)
        bytes_eq2_dm = max(slab, cube, xp_buf_dm, hoc)
    else:
        bytes_eq2    = max(cube, xp_buf, hoc)
        bytes_eq2_dm = max(cube, xp_buf_dm, hoc)
    bytes_eq3    = max(slab_cache, ll, recv_buf)
    bytes_eq3_dm = max(slab_cache, ll_dm, recv_buf_dm)
    bytes_eq4    = xvp
    bytes_eq4_dm = xvp_dm
    bytes_oth    = phi_buf
    if slab_fftw3_ICs:
        bytes_oth += slab_cmplx
        
    #
    # Determine memory usage of large P3DFFT arrays if applicable
    #

    bytes_p3dfft = 0

    if pencil_ICs:

        if nodes_dim == 1:
            nm = int(nc_node_dim * nc_node_dim * (nc_node_dim+2) / 2)
        else:
            nm = int(nc_node_dim * (nc_node_dim + nodes_dim) * (nc_node_dim + 2./nodes_dim) / 2.)

        buf1 = 4 * 2 * nm
        buf2 = 4 * 2 * nm
        R    = 4 * 2 * nm
        bytes_p3dfft = buf1 + buf2 + R

    #
    # Print memory usage
    #

    totalbytes = bytes_eq1 + bytes_eq2 + bytes_eq3 + bytes_eq4 + bytes_oth + bytes_p3dfft
    totalbytes_dm = bytes_eq1_dm + bytes_eq2_dm + bytes_eq3_dm + bytes_eq4_dm + bytes_oth + bytes_p3dfft

    print "TOTAL NU MEMORY USAGE:  " + str(totalbytes/1024.**3) + " GB"
    print "TOTAL DM MEMORY USAGE:  " + str(totalbytes_dm/1024.**3) + " GB"
    if pencil:
        print "P3DFFT MEMORY USAGE: " + str(bytes_p3dfft/1024.**3) + " GB"
    print

    #
    # Print notes to screen
    #

    for i in xrange(2):

        if i == 0:
            bb = "NEUTRINO EXECUTABLE"
            b1 = bytes_eq1
            b2 = bytes_eq2
            b3 = bytes_eq3

        if i == 1:
            bb = "DARK MATTER EXECUTABLE"
            b1 = bytes_eq1_dm
            b2 = bytes_eq2_dm
            b3 = bytes_eq3_dm

        if b1 != phi:
            print " *** Adjust common block for phi, slab_work, recv_cube, and send_buf in " + bb + " *** "
            print "     phi       = ", phi
            print "     slab_work = ", slab_work
            print "     recv_cube = ", recv_cube
            if i == 0:
                print "   send_buf_nu = ", send_buf
            else:
                print "   send_buf_dm = ", send_buf_dm
            print

        if not slab_fftw3_ICs :
            if b2 != slab:
                print " *** Adjust common block for slab, cube, hoc, and xp_buf in " + bb + " *** "
                print "     slab = ", slab
                print "     cube = ", cube
                if i == 0:
                    print "xp_buf_nu = ", xp_buf
                else:
                    print "xp_buf_dm = ", xp_buf_dm
                print "      hoc = ", hoc
                print
        else:
            if b2 != cube:
                print " *** Adjust common block for cube, hoc, and xp_buf in " + bb + " *** "
                print "     cube = ", cube
                if i == 0:
                    print "xp_buf_nu = ", xp_buf
                else:
                    print "xp_buf_dm = ", xp_buf_dm
                print "      hoc = ", hoc
                print

        if b3 != slab_cache:
            print " *** Adjust common block for slab_cache, ll, and recv_buf in " + bb + " *** " 
            print "     slab_cache = ", slab_cache
            if i == 0:
                print "          ll_nu = ", ll
                print "    recv_buf_nu = ", recv_buf 
            else:
                print "          ll_dm = ", ll_dm
                print "    recv_buf_dm = ", recv_buf_dm


else:

    xvp     = 4 * 6 * np_node_dim**2 * num_threads_ic    

    #
    # Determine memory usage of equivalenced and other arrays
    #

    bytes_eq1 = max(phi, slab_work, recv_cube)
    if not slab_fftw3_ICs:
        bytes_eq2 = max(slab, cube)
    else:
        bytes_eq2 = max(slab_cmplx, slab)
    bytes_oth = phi_buf + xvp
    if slab_fftw3_ICs:
        bytes_oth += cube

    #
    # Determine memory usage of large P3DFFT arrays if applicable
    #

    bytes_p3dfft = 0

    if pencil_ICs:

        if nodes_dim == 1:
            nm = int(nc_node_dim * nc_node_dim * (nc_node_dim+2) / 2)
        else:
            nm = int(nc_node_dim * (nc_node_dim + nodes_dim) * (nc_node_dim + 2./nodes_dim) / 2.)

        buf1 = 4 * 2 * nm
        buf2 = 4 * 2 * nm
        R    = 4 * 2 * nm
        bytes_p3dfft = buf1 + buf2 + R

    #
    # Print memory usage
    #

    totalbytes = bytes_eq1 + bytes_eq2 + bytes_oth + bytes_p3dfft

    print "TOTAL MEMORY USAGE:  " + str(totalbytes/1024.**3) + " GB"
    if pencil:
        print "P3DFFT MEMORY USAGE: " + str(bytes_p3dfft/1024.**3) + " GB"
    print

    #
    # Print notes to screen
    #

    if bytes_eq1 != phi:
        print " *** Adjust common block for phi, slab_work, and recv_cube *** " 
        print "     phi       = ", phi
        print "     slab_work = ", slab_work
        print "     recv_cube = ", recv_cube
        print

    if not slab_fftw3_ICs:
        if bytes_eq2 != slab:
            print " *** Adjust common block for slab and cube *** "
            print "     slab = ", slab
            print "     cube = ", cube
            print

