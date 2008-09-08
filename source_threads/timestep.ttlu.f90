!! increment timestep 
  subroutine timestep

#ifdef MHD
    use mpi_tvd_mhd
#endif

    implicit none
 
    include 'mpif.h'
    include 'cubepm.fh'

    real(4) :: ra,da_1,da_2,dt_e
    integer(4) :: n

#ifdef MHD
    real(4) :: cfl,amax,cmax,D1,D2,dta,dtc
#endif

    nts = nts + 1
    if (nts /= 1) dt_old = dt

#ifdef MHD
!! calculate cfl for gas and calculate minimum gas timestep
    call calcfl(u,b,nx,ny,nz,cfl)
    call fluid_minmax(u,amax,cmax,D1,D2)
    if (rank==0) then
      print*,'cfl=',cfl,'amax=',amax,'cmax=',cmax
      print*,'min gas density =',D1,'max gas density =',D2
    endif 
    cmax=cmax*cfactor
    amax=amax*cfactor
    freeze=cmax
#ifdef DEBUG_CFL
!    print *,'cfl hardwired to 0.7'
    dta=0.7/amax
    dtc=0.7/cmax
    dt_f_acc=0.7*dt_f_acc
    dt_c_acc=0.7*dt_c_acc
#ifdef PPINT
    dt_pp_acc=0.7*dt_pp_acc
#endif
#else
    dta=cfl/amax
    dtc=cfl/cmax
    dt_f_acc=cfl*dt_f_acc
    dt_c_acc=cfl*dt_c_acc
#ifdef PPINT
    dt_pp_acc=cfl*dt_pp_acc
#endif
#endif
#endif

    if (rank == 0) then

      if (cosmo) then

        dt_e=dt_max

! restrict expansion

        n=0
        do
          n=n+1
          call expansion(a,dt_e,da_1,da_2)
          da=da_1+da_2
          ra=da/(a+da)
          if (ra.gt.ra_max) then
            dt_e=dt_e*(ra_max/ra)
          else
            exit
          endif
          if (n .gt. 10) exit
        enddo

#ifdef RESTRICT_DA
        n=0
        do
          call expansion(a,dt_e,da_1,da_2)
          da=da_1+da_2
          if (da.gt.da_max) then
            dt_e=dt_e*(da_max/da)
          else
            exit
          endif
          n=n+1
          if (n .gt. 10) exit 
        enddo
#endif

! take the minimum of all the limits 

#ifdef MHD
#ifdef PPINT
        dt = min(dt_e,dt_f_acc,dt_pp_acc,dt_c_acc,dta,dtc)
#else
        dt = min(dt_e,dt_f_acc,dt_c_acc,dta,dtc)
#endif
#else
#ifdef PPINT
        dt = min(dt_e,dt_f_acc,dt_pp_acc,dt_c_acc)
#else
        dt = min(dt_e,dt_f_acc,dt_c_acc)
#endif
#endif

        dt = dt * dt_scale

        call expansion(a,dt,da_1,da_2)

        da=da_1+da_2 

! Check to see if we are checkpointing / projecting / halofinding this step 

        checkpoint_step=.false.
        projection_step=.false.
        halofind_step=.false.

        if (a_checkpoint(cur_checkpoint) <= a_projection(cur_projection) &
            .and. a_checkpoint(cur_checkpoint) <= a_halofind(cur_halofind) &
            .and. cur_checkpoint <= num_checkpoints) then 

          if (a+da > a_checkpoint(cur_checkpoint)) then
            checkpoint_step=.true.
            dt=dt*(a_checkpoint(cur_checkpoint)-a)/da
            call expansion(a,dt,da_1,da_2)
            if (cur_checkpoint == num_checkpoints) final_step=.true.
          endif

          if (checkpoint_step.and.a_projection(cur_projection) == &
              a_checkpoint(cur_checkpoint).and.cur_projection <= &
              num_projections) projection_step=.true.

          if (checkpoint_step.and.a_halofind(cur_halofind) == &
              a_checkpoint(cur_checkpoint).and.cur_halofind <= &
              num_halofinds)  halofind_step=.true.

        elseif (a_checkpoint(cur_checkpoint) <= a_projection(cur_projection) &
            .and. cur_checkpoint <= num_checkpoints) then

          if (a+da > a_halofind(cur_halofind)) then
            halofind_step=.true.
            dt=dt*(a_halofind(cur_halofind)-a)/da
            call expansion(a,dt,da_1,da_2)
          endif

        elseif (a_checkpoint(cur_checkpoint) <= a_halofind(cur_halofind) &
            .and. cur_checkpoint <= num_checkpoints) then

          if (a+da > a_projection(cur_projection)) then
            projection_step=.true.
            dt=dt*(a_projection(cur_projection)-a)/da
            call expansion(a,dt,da_1,da_2)
          endif

        else

          if (a_projection(cur_projection) <= a_halofind(cur_halofind) &
              .and. cur_projection <= num_projections) then

            if (a+da > a_projection(cur_projection)) then
              projection_step=.true.
              dt=dt*(a_projection(cur_projection)-a)/da
              call expansion(a,dt,da_1,da_2)
            endif
        
            if (projection_step.and.a_projection(cur_projection) == &
                a_halofind(cur_halofind).and.cur_halofind <= &
                num_halofinds) halofind_step=.true.
 
          else

            if (a+da > a_halofind(cur_halofind) .and. cur_halofind <= &
                num_halofinds) then
              halofind_step=.true.
              dt=dt*(a_halofind(cur_halofind)-a)/da
              call expansion(a,dt,da_1,da_2)
            endif

          endif

        endif

!! Calculate timestep parameters to be used

        dt_gas=dt/4      
        da=da_1+da_2 
        ra=da/(a+da)
        a_mid=a+(da/2) !da_1

        write(*,*)
        write(*,*) 'Sweep number: ',nts
        write(*,*) 'Tau         : ',tau,tau+dt
        write(*,*) 'Redshift    : ',1.0/a-1.0,1.0/(a+da)-1.0
        write(*,*) 'Scale factor: ',a,a_mid,a+da
        write(*,*) 'Expansion   : ',ra
#ifdef MHD
#ifdef PPINT
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_pp_acc,dt_c_acc,dta,dtc
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_c_acc,dta,dtc
#endif
#else
#ifdef PPINT
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_pp_acc,dt_c_acc
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_c_acc
#endif
#endif

        tau=tau+dt
        t=t+dt
        a=a+da

      else
   
        a = 1.0
        a_mid = a
        da = 0.0
#ifdef PPINT
        if (pair_infall) then
!          dt = min(0.1/sqrt(G*mass_p/cur_sep**2),dt_f_acc,dt_pp_acc,dt_c_acc,dt_max_v)
          dt = min(0.05/sqrt(G*mass_p/cur_sep**2),dt_f_acc,dt_pp_acc,dt_c_acc)
        else
          dt = min(1.0,dt_f_acc,dt_pp_acc,dt_c_acc)
        endif
#else
        dt = min(1.0,dt_f_acc,dt_c_acc)
#endif
        if (pairwise_ic) dt=1.0
        if (shake_test_ic) dt=1.0
        t = t + dt
        if (rank == 0) write(*,*) 'nts=',nts,'t=',t,'dt=',dt,dt_f_acc,dt_pp_acc,dt_c_acc,dt_max_v,0.1/sqrt(G*mass_p/cur_sep**2)
  
      endif

    endif

    ! broadcast timestep variables 

    call mpi_bcast(a,1,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(a_mid,1,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(dt,1,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(dt_gas,1,mpi_real,0,mpi_comm_world,ierr)
    call mpi_bcast(checkpoint_step,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(projection_step,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(halofind_step,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(final_step,1,mpi_logical,0,mpi_comm_world,ierr)

  end subroutine timestep

!! Expansion subroutine :: Hy Trac -- trac@cita.utoronto.ca
!! Added Equation of State for Dark Energy :: Pat McDonald -- pmcdonal@cita.utoronto.ca
  subroutine expansion(a0,dt0,da1,da2)
    implicit none

    include 'cubepm.par'

    real(4) :: a0,dt0,dt_x,da1,da2
    real(8) :: a_x,adot,addot,atdot,arkm,a3rlm,omHsq

    !! Expand Friedman equation to third order and integrate
    dt_x=dt0/2
    a_x=a0
    omHsq=4.0/9.0
!    a3rlm=a_x**3*omega_l/omega_m
    a3rlm=a_x**(-3*wde)*omega_l/omega_m
    arkm=a_x*(1.0-omega_m-omega_l)/omega_m
    adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
    da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

    a_x=a0+da1
    omHsq=4.0/9.0
!    a3rlm=a_x**3*omega_l/omega_m
    a3rlm=a_x**(-3*wde)*omega_l/omega_m
    arkm=a_x*(1.0-omega_m-omega_l)/omega_m
    adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
    da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  end subroutine expansion

