!! add mass to coarse mesh density along nodal boundry
  subroutine coarse_cic_mass_boundry(pp)
    use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif
    implicit none

#    include "cubepm.fh"

    integer(4) :: pp
    integer(4), dimension(3) :: i1,i2
    real(4), dimension(3) :: x,dx1,dx2
#ifdef DEBUG_CCIC
    integer(4) :: i,j,k
#endif

#if defined(COARSEPROJ) && defined(NEUTRINOS)
#ifdef NUPID
    integer, parameter :: pdm = 0
#else
    integer, parameter :: pdm = 1
#endif
#endif

    do
      if (pp == 0) exit
      x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
#ifdef COARSE_NGP
      dx1(:) = 0.0
      dx2(:) = 1.0
#else
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1.0 - dx1(:)
#endif
#ifdef MHD
      dx1(1) = mass_p * dx1(1) * (1.0 - omega_b/omega_m)
      dx2(1) = mass_p * dx2(1) * (1.0 - omega_b/omega_m)
#else

#ifdef NEUTRINOS
#ifdef NUPID
      dx1(1) = mass_p * dx1(1) * mass_p_nudm_fac(nuPIDmap(PID(pp)))
      dx2(1) = mass_p * dx2(1) * mass_p_nudm_fac(nuPIDmap(PID(pp)))
#else
      dx1(1) = mass_p * dx1(1) * mass_p_nudm_fac(PID(pp)) 
      dx2(1) = mass_p * dx2(1) * mass_p_nudm_fac(PID(pp))
#endif
#else
      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)
#endif

#endif

#ifdef DEBUG_CCIC
!      if (i1(1) == 1 .and. i1(2) == 1 .and. &
!          i1(3) == 1 .or. i2(1) == 1 .and. &
!          i2(2) == 1 .and. i2(3) == 1) then
        write(*,*) 'boundary',pp,x(:),i1(:),i2(:),dx1(:),dx2(:)
!      endif

       
       do k=i1(3),i2(3)
         do j=i1(2),i2(2)
           do i=i1(1),i2(1)
             write(*,*) rank,i,j,k,rho_c(i,j,k)
           enddo
         enddo
       enddo

#endif

      if (i1(3) >= 1 .and. i1(3) <= nc_node_dim) then
        if (i1(2) >= 1 .and. i1(2) <= nc_node_dim) then

          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then

            rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) + &
                                       dx1(1) * dx1(2) * dx1(3)
          endif
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then

            rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) + &
                                       dx2(1) * dx1(2) * dx1(3)
          endif
        endif

        if (i2(2) >= 1 .and. i2(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then

            rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) + &
                                       dx1(1) * dx2(2) * dx1(3)
          endif
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then

            rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) + &
                                       dx2(1) * dx2(2) * dx1(3)
          endif
        endif
      endif

      if (i2(3) >= 1 .and. i2(3) <= nc_node_dim) then
        if (i1(2) >= 1 .and. i1(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then

            rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) + &
                                       dx1(1) * dx1(2) * dx2(3)
          endif
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then

            rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) + &
                                       dx2(1) * dx1(2) * dx2(3)
          endif
        endif
        if (i2(2) >= 1 .and. i2(2) <= nc_node_dim) then
          if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then

            rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) + &
                                       dx1(1) * dx2(2) * dx2(3)
          endif
          if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then

            rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) + &
                                       dx2(1) * dx2(2) * dx2(3)
          endif
        endif
      endif

#ifdef DEBUG_CCIC
       do k=i1(3),i2(3)
         do j=i1(2),i2(2)
           do i=i1(1),i2(1)
             write(*,*) 'after',rank,i,j,k,rho_c(i,j,k)
           enddo
         enddo
       enddo
#endif

#ifdef COARSEPROJ
      if (doCoarseProj) then
          x(:)  = x(:) - soffcproj
          i1(:) = floor(x(:)) + 1
          i2(:) = i1(:) + 1
#ifdef COARSE_NGP
          dx1(:) = 0.0
          dx2(:) = 1.0
#else
          dx1(:) = i1(:) - x(:)
          dx2(:) = 1.0 - dx1(:)
#endif
          if (all(i1 >= 0) .and. all(i2 <= nc_node_dim+1)) then
#ifdef NEUTRINOS
            if (PID(pp) == pdm) then
#endif
              crhoproj(i1(1),i1(2),i1(3)) = crhoproj(i1(1),i1(2),i1(3)) &
                                       + dx1(1) * dx1(2) * dx1(3)

              crhoproj(i2(1),i1(2),i1(3)) = crhoproj(i2(1),i1(2),i1(3)) &
                                       + dx2(1) * dx1(2) * dx1(3)

              crhoproj(i1(1),i2(2),i1(3)) = crhoproj(i1(1),i2(2),i1(3)) &
                                       + dx1(1) * dx2(2) * dx1(3)

              crhoproj(i2(1),i2(2),i1(3)) = crhoproj(i2(1),i2(2),i1(3)) &
                                       + dx2(1) * dx2(2) * dx1(3)

              crhoproj(i1(1),i1(2),i2(3)) = crhoproj(i1(1),i1(2),i2(3)) &
                                       + dx1(1) * dx1(2) * dx2(3)

              crhoproj(i2(1),i1(2),i2(3)) = crhoproj(i2(1),i1(2),i2(3)) &
                                       + dx2(1) * dx1(2) * dx2(3)

              crhoproj(i1(1),i2(2),i2(3)) = crhoproj(i1(1),i2(2),i2(3)) &
                                       + dx1(1) * dx2(2) * dx2(3)

              crhoproj(i2(1),i2(2),i2(3)) = crhoproj(i2(1),i2(2),i2(3)) &
                                       + dx2(1) * dx2(2) * dx2(3)
#ifdef NEUTRINOS
            else
              crhoproj_nu(i1(1),i1(2),i1(3)) = crhoproj_nu(i1(1),i1(2),i1(3)) &
                                       + dx1(1) * dx1(2) * dx1(3)

              crhoproj_nu(i2(1),i1(2),i1(3)) = crhoproj_nu(i2(1),i1(2),i1(3)) &
                                       + dx2(1) * dx1(2) * dx1(3)

              crhoproj_nu(i1(1),i2(2),i1(3)) = crhoproj_nu(i1(1),i2(2),i1(3)) &
                                       + dx1(1) * dx2(2) * dx1(3)

              crhoproj_nu(i2(1),i2(2),i1(3)) = crhoproj_nu(i2(1),i2(2),i1(3)) &
                                       + dx2(1) * dx2(2) * dx1(3)

              crhoproj_nu(i1(1),i1(2),i2(3)) = crhoproj_nu(i1(1),i1(2),i2(3)) &
                                       + dx1(1) * dx1(2) * dx2(3)

              crhoproj_nu(i2(1),i1(2),i2(3)) = crhoproj_nu(i2(1),i1(2),i2(3)) &
                                       + dx2(1) * dx1(2) * dx2(3)

              crhoproj_nu(i1(1),i2(2),i2(3)) = crhoproj_nu(i1(1),i2(2),i2(3)) &
                                       + dx1(1) * dx2(2) * dx2(3)

              crhoproj_nu(i2(1),i2(2),i2(3)) = crhoproj_nu(i2(1),i2(2),i2(3)) &
                                       + dx2(1) * dx2(2) * dx2(3)
            endif
#endif
          endif
      endif
#endif 

      pp = ll(pp)
    enddo

  end subroutine coarse_cic_mass_boundry
