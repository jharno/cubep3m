! -------------------------------------------------------------------------------------------------------
! CubeP3M inline halofinder written by JD Emberson September 2013.
! Determines mass of halos with overdensity halo_vir (function of redshift) as
! well as halo_odc (specified in cubepm.par; independent of redshift) by  
! searching over the local particle distribution around each fine mesh density peak.
! -------------------------------------------------------------------------------------------------------
 
subroutine halofind 
    use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif
    implicit none

    include "mpif.h"
#    include "cubepm.fh"

    real(4) :: z_write
    integer(4) :: i, j, k, fstat
    integer(4), dimension(3) :: tile
    integer(4) :: num_candidates

    character (len=max_path) :: ofile
    character (len=7) :: z_s
    character (len=3) :: t_s
    character (len=5) :: r_s

    integer(4) :: pp, iloc, ii, jj 
    integer(8) :: imass_odc, i_odc
    real(4) :: mass_proxy, mass_odc, r_odc
    real(4) :: v_disp
    real(4), dimension(3) :: x_mean, x2_mean, var_x, v_mean, v2_mean, l, offset, dx, l_CM
    real(4), dimension(3) :: r_wrt_halo, v_wrt_halo, v2_wrt_halo, hpos
    integer, parameter :: N_p = 50
    real(4), dimension(N_p) :: E
#ifndef NEUTRINOS
#ifdef PID_FLAG
    integer(8), dimension(N_p) :: pid_halo
    real(4), dimension(6,N_p) ::xv_halo
#endif
#else
    real(4), dimension(3) :: x_mean_nu, v_mean_nu
    integer(4) :: n_nu
#endif
    real(4) :: dist, speed, E_tmp
    real(4), dimension(6) :: I_ij
    real(8) :: st1, st2
    integer(8) :: np_halo_local_odc, np_halo_odc, nhalo_tot, num_candidates_tot
    real(4) :: xflat, halo_vir, mass_vir, r_vir
    integer(8) :: imass_vir, i_vir
    integer(8) :: np_halo_local_vir, np_halo_vir

#ifdef NEUTRINOS
    real(4), parameter :: nu_search_radius = 2. !! Distance in Mpc/h to search for neutrino properties around each halo  
    real(4), parameter :: nu_search_radius_cells = nu_search_radius * nf_physical_dim / box 
    real(4) :: g4 = box * nf_buf / real(nf_physical_dim) 
    if (rank == 0) write(*,*) "Neutrino search radius: ", nu_search_radius_cells
    if (rank == 0 .and. g4 < nu_search_radius) write(*,*) "WARNING: Buffer size is only ", g4, " < ", nu_search_radius
#endif

    !
    ! Find halo candidates based on local overdensities for each tile
    !

    num_candidates = 0

    do i = 1, tiles_node
        tile(3) = (i-1) / (tiles_node_dim * tiles_node_dim)
        j = i - tile(3) * tiles_node_dim * tiles_node_dim
        tile(2) = (j-1) /  tiles_node_dim
        j = j - tile(2) * tiles_node_dim
        tile(1) = j - 1
        call find_halo_candidates(tile, num_candidates)
    enddo
    
    !
    ! Sort density maxima 
    !

    isortpeak(:num_candidates) = (/ (i, i=1, num_candidates) /)
    call indexedsort(num_candidates, den_peak(:), isortpeak(:))
    ipeak(:, :num_candidates)  = ipeak(:, isortpeak(:num_candidates))
    halo_mesh_mass(:num_candidates) = halo_mesh_mass(isortpeak(:num_candidates))

    !
    ! Determine total number of halo candidates
    !

    call mpi_reduce(int(num_candidates, kind=8), num_candidates_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    !
    ! Determine Delta_vir using equation (6) of Bryan et al. (1997) for a flat universe
    !

    xflat    = (omega_m / a**3) / (omega_m / a**3 + omega_l) - 1.
    halo_vir = 18.*pi**2 + 82.*xflat - 39.*xflat**2

    !
    ! Open halo file
    !

    z_write = z_halofind(cur_halofind)
    call mpi_bcast(z_write, 1, mpi_real, 0, mpi_comm_world, ierr)
    write(z_s ,"(f7.3)") z_write 
    z_s = adjustl(z_s)
    
    write(r_s, "(i5)") rank
    r_s = adjustl(r_s)

    ofile = output_path//'/node'//r_s(1:len_trim(r_s))//'/'//z_s(1:len_trim(z_s))//"halo"//r_s(1:len_trim(r_s))//".dat"

    open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")

    if (fstat /= 0) then
        write(*,*) "Error opening halo catalog for write"
        write(*,*) "rank", rank, "file:", ofile
        stop
    endif

    !
    ! Determine which candidates are to be considered halos and write their properties to file.
    !

    !! Halo ticker
    nhalo = 0

    !! Ticker recording how many particle searches end up empty handed
    search_fail = 0

    !! Node offsets
    offset(1) = cart_coords(3)*nf_physical_node_dim
    offset(2) = cart_coords(2)*nf_physical_node_dim
    offset(3) = cart_coords(1)*nf_physical_node_dim

    !! Initialize so that no particles are yet part of a halo
    hpart_odc = 0
    hpart_vir = 0

#ifdef NEUTRINOS
    !! Exclude neutrinos from being included in halos by
    !! making the halofinder think they have already been inclued in a halo
    do i = 1, np_local
#ifdef NUPID
        if (PID(i) > 0) then
#else
        if (PID(i) > 1) then
#endif
            hpart_odc(i) = 1
            hpart_vir(i) = 1
        endif
    enddo 
#endif

    !! Initialize refined mesh
    finegrid(:, :, :) = 0.

    !! Write file header (will rewrite nhalo later)
    write(12) nhalo
    write(12) halo_vir
    write(12) halo_odc

    if (rank == 0) then

        write(*,*) "Searching over ", num_candidates_tot, " halo candidates..."
        write(*,*) "halo_vir, halo_odc = ", halo_vir, halo_odc

    endif

    do iloc = num_candidates, 1, -1

        !! Fine mesh mass and position of this halo
        mass_proxy = halo_mesh_mass(iloc)
        hpos(:)    = ipeak(:, iloc)

        !! Search for particles by looking at local particle distribution
        call find_halo_particles(halo_vir, mass_proxy, hpos(:), r_vir, i_vir, 1)
        call find_halo_particles(halo_odc, mass_proxy, hpos(:), r_odc, i_odc, 2)

        !! The following conditions must pass to be considered a halo
        if (i_vir >= min_halo_particles .and. i_odc >= min_halo_particles) then

            !! Halo ticker
            nhalo = nhalo + 1

            !! Initialize halo properties
            imass_vir = 0
            x_mean = 0.
            x2_mean = 0.
            v_mean = 0.
            v2_mean = 0.
            v_wrt_halo = 0.
            v2_wrt_halo = 0.
            l = 0.
            l_CM = 0.

            !! Loop over particles in this halo to get halo properties 
            do ii = 1, i_vir
                pp = ilist_vir(ii)
                imass_vir = imass_vir + 1
                x_mean = x_mean + xv(:3, pp)
                x2_mean = x2_mean + xv(:3, pp)**2
                v_mean = v_mean + xv(4:, pp)
                v2_mean = v2_mean + xv(4:, pp)**2
                dx = hpos(:) - xv(:3, pp)
                l(1) = l(1) + (dx(3)*xv(5,pp) - dx(2)*xv(6,pp))
                l(2) = l(2) + (dx(1)*xv(6,pp) - dx(3)*xv(4,pp))
                l(3) = l(3) + (dx(2)*xv(4,pp) - dx(1)*xv(5,pp))

                !! Remove this particle from future halo candidates
                hpart_vir(pp) = 1
            enddo

            imass_odc = 0
            do ii = 1, i_odc
                pp = ilist_odc(ii)
                imass_odc = imass_odc + 1

                !! Remove this particle from future halo candidates
                hpart_odc(pp) = 1
            enddo
            mass_odc = mass_p * imass_odc

#ifdef NEUTRINOS
            call neutrino_properties(hpos(:), nu_search_radius_cells, x_mean_nu, v_mean_nu, n_nu)
            x_mean_nu = x_mean_nu + offset 
#endif

            mass_vir = mass_p * imass_vir
            hpos(:) = hpos(:) + offset
            x_mean = x_mean/real(imass_vir) + offset
            x2_mean = x2_mean/real(imass_vir)
            v_mean = v_mean/real(imass_vir)
            v2_mean = v2_mean/real(imass_vir)
            l = l/real(imass_vir)
            l_CM(1) = l(1) - (x_mean(3)*v_mean(2) - x_mean(2)*v_mean(3))
            l_CM(2) = l(2) - (x_mean(1)*v_mean(3) - x_mean(3)*v_mean(1))
            l_CM(3) = l(3) - (x_mean(2)*v_mean(1) - x_mean(1)*v_mean(2))
            v_disp = sqrt(v2_mean(1) + v2_mean(2) + v2_mean(3))
            var_x = real(imass_vir)/(real(imass_vir-1)) * (x2_mean - (x_mean-offset)**2)

#ifndef NEUTRINOS
#ifdef PID_FLAG
            pid_halo = 0
#endif
#endif
            E = 0.
            E_tmp = 0.
            I_ij = 0.

            ii = 1
            do ii = 1, i_vir
                pp = ilist_vir(ii)

                r_wrt_halo = xv(:3,pp) - (x_mean - offset)
                v_wrt_halo = xv(4:,pp) - v_mean
                v2_wrt_halo = v2_wrt_halo + v_wrt_halo(:)**2
                dist = sqrt(r_wrt_halo(1)**2 +r_wrt_halo(2)**2 + r_wrt_halo(3)**2)
                I_ij(1) = I_ij(1)+r_wrt_halo(1)*r_wrt_halo(1)
                I_ij(2) = I_ij(2)+r_wrt_halo(1)*r_wrt_halo(2)
                I_ij(3) = I_ij(3)+r_wrt_halo(1)*r_wrt_halo(3)
                I_ij(4) = I_ij(4)+r_wrt_halo(2)*r_wrt_halo(2)
                I_ij(5) = I_ij(5)+r_wrt_halo(2)*r_wrt_halo(3)
                I_ij(6) = I_ij(6)+r_wrt_halo(3)*r_wrt_halo(3)

                speed = sqrt(v_wrt_halo(1)**2 +v_wrt_halo(2)**2 + v_wrt_halo(3)**2)
                E_tmp = 0.5*(speed)**2 - imass_vir*mass_p*G/dist

                !! Find the most bound particles
                do jj = 1, N_p
                    if (E_tmp < E(jj)) then
                        E(jj+1:N_p) = E(jj:N_p-1)
                        E(jj) = E_tmp
#ifndef NEUTRINOS
#ifdef PID_FLAG
                        pid_halo(jj+1:N_p) = pid_halo(jj:N_p-1)
                        pid_halo(jj) = PID(pp)
                        xv_halo(:,jj) = xv(:,pp)
#endif
#endif
                        exit
                    endif
                enddo
            enddo

            v2_wrt_halo(:) = v2_wrt_halo(:)/real(imass_vir)

#ifdef DISP_MESH
            !! Subtract shake offset
            hpos = hpos - shake_offset
            x_mean = x_mean - shake_offset
#ifdef NEUTRINOS 
            x_mean_nu = x_mean_nu - shake_offset
#endif
#endif

#ifndef NEUTRINOS
#ifdef PID_FLAG
            if (halo_write) write(12) hpos(:), mass_vir, mass_odc, r_vir, r_odc, x_mean, v_mean, l_CM, v2_wrt_halo, var_x, pid_halo, xv_halo
#else
            if (halo_write) write(12) hpos(:), mass_vir, mass_odc, r_vir, r_odc, x_mean, v_mean, l_CM, v2_wrt_halo, var_x
#endif
#else
            if (halo_write) write(12) hpos(:), mass_vir, mass_odc, r_vir, r_odc, x_mean, v_mean, l_CM, v2_wrt_halo, var_x, I_ij, x_mean_nu, v_mean_nu, n_nu 
#endif


        endif !! i_vir/i_odc test 

    enddo !! iloc loop 

    !
    ! Rewrite nhalo in the header and close the file
    !

    rewind(12)
    write(12) nhalo
    close(12)

    !
    ! Count particles that were found within halos
    !

    np_halo_local_odc = 0
    do ii = 1, np_local
#ifdef NEUTRINOS
#ifdef NUPID
        if (PID(ii) == 0 .and. hpart_odc(ii) == 1) np_halo_local_odc = np_halo_local_odc + 1
#else
        if (PID(ii) == 1 .and. hpart_odc(ii) == 1) np_halo_local_odc = np_halo_local_odc + 1
#endif
#else
        if (hpart_odc(ii) == 1) np_halo_local_odc = np_halo_local_odc + 1
#endif
    enddo
    call mpi_reduce(np_halo_local_odc, np_halo_odc, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    np_halo_local_vir = 0
    do ii = 1, np_local
#ifdef NEUTRINOS
#ifdef NUPID
        if (PID(ii) == 0 .and. hpart_vir(ii) == 1) np_halo_local_vir = np_halo_local_vir + 1
#else
        if (PID(ii) == 1 .and. hpart_vir(ii) == 1) np_halo_local_vir = np_halo_local_vir + 1
#endif
#else
        if (hpart_vir(ii) == 1) np_halo_local_vir = np_halo_local_vir + 1
#endif
    enddo
    call mpi_reduce(np_halo_local_vir, np_halo_vir, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    call mpi_reduce(int(nhalo, kind=8), nhalo_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    !
    ! Print some info to screen
    !

    if (rank == 0) then

        write(*,*) "Done ... found ", nhalo_tot, " halos"
        write(*,*) "Total Particles found in halos = ", np_halo_vir, np_halo_odc

        if (search_fail > 0) then
            write(*,*) "CONSIDER INCREASING search_ratio ... search_fail = ", search_fail
        endif

    endif

    !
    ! Increment halofind counter 
    !

    cur_halofind  = cur_halofind + 1
    halofind_step = .false.

end subroutine halofind 

! -------------------------------------------------------------------------------------------------------

subroutine find_halo_candidates(tile, ic)
    !
    ! Finds density maxima and returns a sorted array containing potential halo candidates.
    !
    use omp_lib
    implicit none

#    include "cubepm.fh"
!    include "mpif.h"

    integer(4), dimension(3) :: offset
    integer(4), dimension(3) :: cic_l,cic_h,tile
    integer(4), intent(inout) :: ic

    integer(4) :: i, j, k, pp, thread
    integer(4) :: ii, ix, iy, iz, ix0, iy0, iz0
    real(4) :: amtot, denmax

    !! Tile offset in local coordinates

    offset = tile * nf_physical_tile_dim - nf_buf

    !! Initialize density
    thread = 1
    rho_f(:, :, :, thread) = 0.

    !! Limits for mass assignment. Ignore outermost buffer (4 fine cells)
    cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1

    !! Calculate fine mesh density for tile

    do k = cic_l(3), cic_h(3)
        do j = cic_l(2), cic_h(2)
            do i = cic_l(1), cic_h(1)
                pp = hoc(i, j, k)
#ifdef NGPH
#ifdef NEUTRINOS
                call fine_ngp_mass(pp, tile, thread, 1)
#else
                call fine_ngp_mass(pp, tile, thread)
#endif
#else
#ifdef NEUTRINOS
                call fine_cic_mass(pp, tile, thread, 1)
#else
                call fine_cic_mass(pp, tile, thread)
#endif
#endif
            enddo
        enddo
    enddo

    !
    ! Find density maxima 
    !

    do k = 1 + nf_buf, nf_buf + nf_physical_tile_dim
        do j = 1 + nf_buf, nf_buf + nf_physical_tile_dim
            do i = 1 + nf_buf, nf_buf + nf_physical_tile_dim

                !! Maximum density around nearest cells
                denmax   = maxval(rho_f(i-1:i+1, j-1:j+1, k-1:k+1, thread))

                !! Continue if this cell is a large enough local maxima 
                if (denmax == rho_f(i, j, k, thread) .and. denmax > den_peak_cutoff) then

                    if (ic > max_maxima -2) then
                        write(*,*) "Too many halos: ic, max_maxima", ic, max_maxima
                        exit
                    endif

                    !! Find the fine mesh mass of this peak
                    amtot = 0.
                    ix0 = i
                    iy0 = j
                    iz0 = k
                    do ii = 1, irtot
                        ix = ix0 + idist(1, ii)
                        if (ix < 5 .or. ix > nf_tile-4) cycle
                        iy = iy0 + idist(2, ii)
                        if (iy < 5 .or. iy > nf_tile-4) cycle
                        iz = iz0 + idist(3, ii)
                        if (iz < 5 .or. iz > nf_tile-4) cycle
                        amtot = amtot + rho_f(ix, iy, iz, thread)
                        if (complete_shell .and. rdist(ii) == rdist(ii+1)) cycle
                        if (ii > 18 .and. amtot/(real(ii)) < halo_odc) exit
                    enddo

                    !! Consider this a halo candidate if the fine mesh mass is large enough
                    if (amtot > mass_p*min_halo_particles/2.) then

                        ic = ic + 1

                        !! Store integer coordinates of the peak as well as its mesh mass
                        ipeak(:, ic) = (/real(i), real(j), real(k)/) - 0.5 + offset
                        den_peak(ic) = denmax
                        halo_mesh_mass(ic) = amtot

                    endif

                endif !! denmax

            enddo !! i loop
        enddo !! j loop
    enddo !! k loop

end subroutine find_halo_candidates

! -------------------------------------------------------------------------------------------------------

subroutine find_halo_particles(HODC, HMASS, HPOS, RODC, ITOT, DOVIR)
    !
    ! Finds refined density maxima and then sorts particles based on their
    ! distance from here, moving one particle at a time until the target radius
    ! and overdensity are reached.
    !
    use omp_lib
    implicit none

#    include "cubepm.fh"
    include "mpif.h"

    real(4), intent(in) :: HODC, HMASS
    real(4), dimension(3), intent(inout) :: HPOS
    real(4), intent(out) :: RODC
    integer(8), intent(out) :: ITOT
    integer(4), intent(in) :: DOVIR
    logical :: HALOVIR

    real :: r, dgrid, maxfinegrid, rrefine, rsearch
    integer :: pp, np_search, ii, jj, kk, i, j, k
    integer :: crbox(3, 2), csbox(3, 2), frbox(3, 2)
    integer :: ngrid(3)
    real :: dr(3), p(3)

    real :: odci, odcj, r1, r2, d1, d2, w1, w2

    integer, parameter :: search_ratio = 4
    integer, parameter :: refine_ratio = 5

    !
    ! Determine which overdensity we are trying to reach 
    !

    if (DOVIR == 1) then

        HALOVIR = .true.

    else

        HALOVIR = .false.

    endif

    !
    ! Initialize parameters
    !

    !! Use the mass proxy to guess at the size of the halo
    !! This will be the radius within which the refined density peak will be found 
    rrefine = (3. * HMASS / 4. / pi / HODC)**(1./3.)

    !! This (larger) radius will be used to store particle positions determine which are part of the halo
    rsearch = search_ratio * rrefine

    !! Coarse mesh cells within the refined region 
    crbox(:, 1) = int((HPOS(:)-rrefine)/mesh_scale) + 1
    crbox(:, 2) = int((HPOS(:)+rrefine)/mesh_scale) + 1

    !! Coarse mesh cells within the search region
    csbox(:, 1) = int((HPOS(:)-rsearch)/mesh_scale) + 1
    csbox(:, 2) = int((HPOS(:)+rsearch)/mesh_scale) + 1

    !! Boundary of the refined region in fine mesh cell units
    frbox(:, 1) = mesh_scale*(crbox(:, 1) - 1)
    frbox(:, 2) = mesh_scale*crbox(:, 2)

    !! Number of extra refined cells in this region and their spacing
    ngrid(:) = refine_ratio*(frbox(:, 2) - frbox(:, 1))
    dgrid    = 1./refine_ratio

    if (ngrid(1) > ngrid_max .or. ngrid(2) > ngrid_max .or. ngrid(3) > ngrid_max) write(*,*) "ERROR: ngrid = ", ngrid

    !! Exclude coarse mesh cells in the search region that lie within the tile buffer
    csbox(1, 1) = max(csbox(1, 1), hoc_nc_l)
    csbox(1, 2) = min(csbox(1, 2), hoc_nc_h)
    csbox(2, 1) = max(csbox(2, 1), hoc_nc_l)
    csbox(2, 2) = min(csbox(2, 2), hoc_nc_h)
    csbox(3, 1) = max(csbox(3, 1), hoc_nc_l)
    csbox(3, 2) = min(csbox(3, 2), hoc_nc_h)

    !
    ! Store particle positions within the search radius at the same time as NGP
    ! interpolating particles within the refined radius in order to find the
    ! density peak.
    !

    np_search = 0

    if (HALOVIR) then

        do k = csbox(3,1), csbox(3,2)
            do j = csbox(2,1), csbox(2,2)
                do i = csbox(1,1), csbox(1,2)
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (hpart_vir(pp) == 0) then !! particle is not yet part of a halo
                            p(:) = xv(:3, pp)
                            dr   = HPOS(:) - p(:)
                            r    = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                            if (r < rsearch) then
                                np_search = np_search + 1
                                pos(np_search, 1:3) = p(:)
                                ilist_vir(np_search) = pp
                                if (r < rrefine) then
                                    ii = int((p(1)-frbox(1,1))/dgrid) + 1
                                    jj = int((p(2)-frbox(2,1))/dgrid) + 1
                                    kk = int((p(3)-frbox(3,1))/dgrid) + 1
                                    finegrid(ii, jj, kk) = finegrid(ii, jj, kk) + 1
                                endif
                            endif
                        endif
                        pp = ll(pp)
                    enddo !! pp loop
                enddo !! i loop
            enddo !! j loop
        enddo !! k loo

    else

        !! Same as above but do not find a new centre. Use the centre determined
        !! from halo_vir which was called first.

        do k = csbox(3,1), csbox(3,2)
            do j = csbox(2,1), csbox(2,2)
                do i = csbox(1,1), csbox(1,2)
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (hpart_odc(pp) == 0) then !! particle is not yet part of a halo
                            p(:) = xv(:3, pp)
                            dr   = HPOS(:) - p(:)
                            r    = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                            if (r < rsearch) then
                                np_search = np_search + 1
                                pos(np_search, 1:3) = p(:)
                                ilist_odc(np_search) = pp
                            endif
                        endif
                        pp = ll(pp)
                    enddo !! pp loop
                enddo !! i loop
            enddo !! j loop
        enddo !! k loop

    endif

    if (HALOVIR) then

        !
        ! Find refined mesh density maximum
        !

        maxfinegrid = 0.

        do k = 1, ngrid(3)
           do j = 1, ngrid(2)
              do i = 1, ngrid(1)
                 if (finegrid(i, j, k) > maxfinegrid) then
                    maxfinegrid = finegrid(i, j, k)
                    HPOS(1) = frbox(1,1) + (i-0.5)*dgrid
                    HPOS(2) = frbox(2,1) + (j-0.5)*dgrid
                    HPOS(3) = frbox(3,1) + (k-0.5)*dgrid
                 endif
                 finegrid(i, j, k) = 0. !! Set to zero for next candidate
              enddo
           enddo
        enddo

    endif

    !
    ! Sort the particles within the search region based on their distance from the centre.
    !

    if (np_search > max_halo_np) write(*,*) "ERROR: np_search, max_halo_np = ", np_search, max_halo_np

    do i = 1, np_search
        dr = HPOS(:) - pos(i, 1:3)
        pos(i, 4) = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
    enddo

    isortpos(:np_search) = (/ (i, i=1, np_search) /)
    call indexedsort(np_search, pos(:, 4), isortpos(:))
    if (HALOVIR) then
        ilist_vir(:np_search) = ilist_vir(isortpos(:np_search))
    else
        ilist_odc(:np_search) = ilist_odc(isortpos(:np_search))
    endif

    !
    ! Move one particle at a time until we reach the radius for which the desired overdensity is obtained.
    !

    RODC = 0.
    ITOT = 0

    do j = 2, np_search !! Assume the first particle is part of the halo

        odcj = threeover4pi * j * mass_p / pos(j, 4)**3

        if (odcj <= HODC) then

            odci = threeover4pi * (j-1) * mass_p / pos(j-1, 4)**3

            !! Interpolate halo radius
            r2 = log10(pos(j, 4))
            r1 = log10(pos(j-1, 4))
            d2 = log10(odcj)
            d1 = log10(odci)
            w1 = log10(HODC) - d2
            w2 = d1 - log10(HODC)

            RODC = 10**((w1 * r1 + w2 * r2) / (d1 - d2))

            if (pos(j, 4) <= RODC) then
                ITOT = j
            else
                ITOT = j-1
            endif

            exit

        endif

    enddo

    if (ITOT == 0) search_fail = search_fail + 1

end subroutine find_halo_particles 

! -------------------------------------------------------------------------------------------------------

#ifdef NEUTRINOS
subroutine neutrino_properties(HPOS, RSEARCH, XMEAN, VMEAN, NNU) 
    !
    ! Finds the mean velocity of all neutrinos within a radius of RSEARCH from
    ! the halo centre HPOS
    !

    implicit none

#    include "cubepm.fh"
    include "mpif.h"

    real(4), dimension(3), intent(in) :: HPOS
    real(4), intent(in) :: RSEARCH
    real(4), dimension(3), intent(out) :: XMEAN, VMEAN
    integer(4), intent(out) :: NNU

    integer :: csbox(3, 2)
    integer :: pp, k, j, i
    real :: r  
    real :: dr(3), p(3)

    !! Coarse mesh cells within the search region
    csbox(:, 1) = int((HPOS(:)-RSEARCH)/mesh_scale) + 1
    csbox(:, 2) = int((HPOS(:)+RSEARCH)/mesh_scale) + 1

    !! Exclude coarse mesh cells in the search region that lie within the tile buffer
    csbox(1, 1) = max(csbox(1, 1), hoc_nc_l)
    csbox(1, 2) = min(csbox(1, 2), hoc_nc_h)
    csbox(2, 1) = max(csbox(2, 1), hoc_nc_l)
    csbox(2, 2) = min(csbox(2, 2), hoc_nc_h)
    csbox(3, 1) = max(csbox(3, 1), hoc_nc_l)
    csbox(3, 2) = min(csbox(3, 2), hoc_nc_h)
   
    !! Initialize mean velocity and particle counter
    XMEAN(:) = 0.
    VMEAN(:) = 0.
    NNU = 0

    do k = csbox(3,1), csbox(3,2)
        do j = csbox(2,1), csbox(2,2)
            do i = csbox(1,1), csbox(1,2)
                pp = hoc(i, j, k)
                do
                    if (pp == 0) exit
#ifdef NUPID
                    if (PID(pp) > 0) then !! this is a neutrino
#else
                    if (PID(pp) > 1) then !! this is a neutrino
#endif
                        p(:) = xv(:3, pp)
                        dr   = HPOS(:) - p(:)
                        r    = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                        if (r <= RSEARCH) then
                            XMEAN(:) = XMEAN(:) + p(:)
                            VMEAN(:) = VMEAN(:) + xv(4:6, pp)
                            NNU = NNU + 1
                        endif
                    endif
                    pp = ll(pp)
                enddo !! pp loop
            enddo !! i loop
        enddo !! j loop
    enddo !! k loop

    if (NNU > 0) then
        XMEAN(:) = XMEAN(:) / NNU
        VMEAN(:) = VMEAN(:) / NNU
    endif
 
end subroutine neutrino_properties 
#endif 

! -------------------------------------------------------------------------------------------------------

subroutine initialize_halofind

    use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif

    implicit none
#    include "cubepm.fh"
    include 'mpif.h'

    integer(4) :: ii, i, j, k, fstat
    real(4) :: r
    integer(4), allocatable , dimension(:,:) :: idist_tmp
    allocate(idist_tmp(3,nlist))

! Loop through a box of length 2*nc_halo_max
! if cell is within sphere of radius = box length / 2
! include distince in rdist at entry ii
! ordered bottom left to top right

    ii=0
    do i=-nc_halo_max,nc_halo_max
      do j=-nc_halo_max,nc_halo_max
        do k=-nc_halo_max,nc_halo_max
          r=sqrt(real(i)**2+real(j)**2+real(k)**2)
          if (r>nc_halo_max) cycle
          ii=ii+1
          if (ii>nlist) then
            write(*,*) 'ii exceeded ',nlist
            pause
          endif
          idist(:,ii)=(/i,j,k/)
          rdist(ii)=r
        enddo
      enddo
    enddo
    irtot=ii

! sorts the rdist array from lowest to highest radial position
! from center of sphere saves rdist array position values in idist

    isortdist(:ii)=(/ (i,i=1,ii) /)
    call indexedsort(ii,rdist,isortdist)
    if(rank==0)write(*,*)'debug flag 1'    
    !idist(:,:ii)=idist(:,isortdist(:ii))
    idist_tmp(:,:ii)=idist(:,isortdist(:ii))
    idist(:,:ii)=idist_tmp(:,:ii)
    deallocate(idist_tmp)

    if(rank==0)write(*,*)'debug flag 2'
    
end subroutine initialize_halofind

