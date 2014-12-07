!! increment timestep 
  subroutine timestep

#ifdef MHD
    use mpi_tvd_mhd
#endif

    implicit none
 
    include 'mpif.h'
#    include "cubepm.fh"

    real(4) :: ra,da_1,da_2,dt_e,am
    integer(4) :: n

#ifdef MHD
    real(4) :: cfl,amax,cmax,D1,D2,dta,dtc
#endif
    real(4) :: vmax, vmax_local
    real(4) :: Dx

    nts = nts + 1
    if (nts /= 1) dt_old = dt

    !! Compute maximum timestep allowed by the maximum velocity
    vmax_local = 0.
    do n = 1, np_local
        vmax_local = max(vmax_local, maxval(abs(xv(4:6,n))))
    enddo 
    call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
#ifdef DISP_MESH

#ifdef NEUTRINOS
      Dx=real(nf_buf)-0.5*mesh_scale ! 5.5 coarse grids
#else
      Dx = real(nf_buf)-4.*mesh_scale ! 2.0 coarse grids
#endif
#else
    Dx = real(nf_buf)
#endif
    !! fbuf is the fraction of the buffer that the fastest particle is allowed to move. 
    !! As long as the maximum velocity increases no more than zeta = 2/fbuf - 1 compared 
    !! to the previous time step then this method will work.
    dt_vmax = fbuf * Dx / vmax

    if (rank == 0) write(*,*) 'vmax and maximum timestep from vmax=',vmax,dt_vmax

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
    dt_vmax=0.7*dt_vmax
#ifdef PPINT
    dt_pp_acc=0.7*dt_pp_acc
#endif
#else
    dta=cfl/amax
    dtc=cfl/cmax
    dt_f_acc=cfl*dt_f_acc
    dt_c_acc=cfl*dt_c_acc
    dt_vmax=cfl*dt_vmax
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
#ifdef PP_EXT
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc,dta,dtc)
#else
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc,dta,dtc)
#endif
#else
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_c_acc,dta,dtc)
#endif
#else
#ifdef PPINT
#ifdef PP_EXT
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc)
#else
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc)
#endif
#else
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_c_acc)
#endif
#endif

        dt = dt * dt_scale

        call expansion(a,dt,da_1,da_2)

        da=da_1+da_2 

! Check to see if we are checkpointing / projecting / halofinding this step 

        checkpoint_step=.false.
        projection_step=.false.
        halofind_step=.false.
        

        am=min(a_checkpoint(cur_checkpoint),a_projection(cur_projection),a_halofind(cur_halofind))
 
#ifdef debug_timestep
        write(*,*) cur_checkpoint,a_checkpoint(cur_checkpoint),cur_projection, a_projection(cur_projection),cur_halofind, a_halofind(cur_halofind),am, a,da
#endif
        if (a_checkpoint(cur_checkpoint)==am) then
          
          if (a+da > a_checkpoint(cur_checkpoint)) then
            checkpoint_step=.true.
            dt=dt*(a_checkpoint(cur_checkpoint)-a)/da
            call expansion(a,dt,da_1,da_2)
            if (cur_checkpoint == num_checkpoints) final_step=.true.
            if (a_projection(cur_projection) == am .and. cur_projection <= num_projections) projection_step=.true. 
            if (a_halofind(cur_halofind) == am .and. cur_halofind <= num_halofinds) halofind_step=.true.
          endif

        elseif (a_projection(cur_projection)==am .and. cur_projection <= num_projections) then 

          if (a+da > a_projection(cur_projection)) then
            projection_step=.true.
            dt=dt*(a_projection(cur_projection)-a)/da
            call expansion(a,dt,da_1,da_2)
            if (a_halofind(cur_halofind) == am .and. cur_halofind <= num_halofinds) halofind_step=.true.
          endif

        elseif (a_halofind(cur_halofind) == am  .and. cur_halofind <= num_halofinds) then
  
          if (a+da > a_halofind(cur_halofind)) then
            halofind_step=.true.
            dt=dt*(a_halofind(cur_halofind)-a)/da
            call expansion(a,dt,da_1,da_2)
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
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc,dta,dtc
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_c_acc,dta,dtc
#endif
#else
#ifdef PPINT
#ifdef PP_EXT
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc
#endif
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_c_acc
#endif
#endif

        tau=tau+dt
        t=t+dt
        a=a+da

      else ! not cosmo
   
        a = 1.0
        a_mid = a
        da = 0.0
#ifdef PPINT
        if (pair_infall) then
!          dt = min(0.1/sqrt(G*mass_p/cur_sep**2),dt_f_acc,dt_pp_acc,dt_c_acc,dt_max_v)
          dt = min(0.05/sqrt(G*mass_p/cur_sep**2),dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc)
        else
#ifdef PP_EXT
           dt = min(1.0,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc)
#else
           dt = min(1.0,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc)
#endif
        endif
#else
        dt = min(1.0,dt_f_acc,dt_vmax,dt_c_acc)
#endif
        if (pairwise_ic) dt=1.0
        if (shake_test_ic) dt=1.0
        t = t + dt
        if (rank == 0) write(*,*) 'nts=',nts,'t=',t,'dt=',dt,dt_f_acc,dt_vmax,dt_pp_acc, dt_pp_ext_acc,dt_c_acc,dt_max_v,0.1/sqrt(G*mass_p/cur_sep**2)
  
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

#    include "cubepm.par"

    real(4) :: a0,dt0,dt_x,da1,da2
    real(8) :: a_x,adot,addot,atdot,arkm,a3rlm,omHsq
    real(8), parameter :: e = 2.718281828459046


#ifdef Chaplygin
    call Chaplygin(a0,dt0,da1,da2)
    return
    write(*,*) '***** IF I SEE THIS, SOMETHING IS WRONG!! *****'
#endif

    !! Expand Friedman equation to third order and integrate
    dt_x=dt0/2
    a_x=a0
    omHsq=4.0/9.0
!    a3rlm=a_x**3*omega_l/omega_m
    a3rlm=a_x**(-3*wde)*omega_l/omega_m
!    a3rlm=a_x**(-3*wde - 3*w_a)*(omega_l/omega_m)*e**(3*w_a*(a_x - 1))
    arkm=a_x*(1.0-omega_m-omega_l)/omega_m

    adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde + w_a*(a_x - 1))*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(3*w_a**2*a_x**2 + 6*w_a*(1-wde-w_a)+ (2.0-3.0*(wde+w_a))*(1.0-(wde+w_a)))*a3rlm)

    da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

    a_x=a0+da1
    omHsq=4.0/9.0
!    a3rlm=a_x**3*omega_l/omega_m
    a3rlm=a_x**(-3*wde)*omega_l/omega_m
!    a3rlm=a_x**(-3*wde - 3*w_a)*(omega_l/omega_m)*e**(3*w_a*(a_x - 1))
    arkm=a_x*(1.0-omega_m-omega_l)/omega_m

    adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde + w_a*(a_x - 1))*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(3*w_a**2*a_x**2 + 6*w_a*(1-wde-w_a)+ (2.0-3.0*(wde+w_a))*(1.0-(wde+w_a)))*a3rlm)

    da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  end subroutine expansion

!! Added Equation of State of Chaplygin gas:: Joachim Harnois-Deraps -- jharno@cita.utoronto.ca
  subroutine Chaplygin(a0,dt0,da1,da2)
    implicit none

#    include "cubepm.par"

    real(4) :: a0,dt0,dt_x,da1,da2
    real(8) :: a_x,adot,addot,atdot,arkm,omHsq, a3rchm,G_ch


    !! Expand Friedman equation to third order and integrate
    dt_x=dt0/2
    a_x=a0
    omHsq=4.0/9.0

    a3rchm=a_x**(-3)*omega_ch/omega_m
    arkm=a_x*(1.0-omega_m-omega_ch)/omega_m
    G_ch = A_ch + (1.0-A_ch)*a_x**(-3.0-3.0*alpha_ch)

    adot=sqrt(omHsq*a_x**3*(1.0 + arkm + a3rchm*(G_ch)**(1.0/(1.0+alpha_ch))))    
    addot=a_x**2*omHsq*(1.5 + 2.0*arkm + 3.0*a3rchm*A_ch*(G_ch)**(-alpha_ch/(1.0+alpha_ch)))
    atdot=a_x*adot*omHsq*(3.0 + 6.0*arkm +3.0*a3rchm*( G_ch)**(1.0/(1.0+alpha_ch) - 2.0)*(5.0*A_ch**2 + 3.0*A_ch*(1.0-A_ch)*a_x**(-3.0-3.0*alpha_ch)*(2.0+alpha_ch/2.0) + (1.0 - A_ch)**2*a_x**(-6.0-6.0*alpha_ch) ))


    da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

    a_x=a0+da1
    omHsq=4.0/9.0


    adot=sqrt(omHsq*a_x**3*(1.0 + arkm + a3rchm*(G_ch)**(1.0/(1.0+alpha_ch))))    
    addot=a_x**2*omHsq*(1.5 + 2.0*arkm + 3.0*a3rchm*A_ch*(G_ch)**(-alpha_ch/(1.0+alpha_ch)))
    atdot=a_x*adot*omHsq*(3.0 + 6.0*arkm +3.0*a3rchm*( G_ch)**(1.0/(1.0+alpha_ch) - 2.0)*(5.0*A_ch**2 + 3.0*A_ch*(1.0-A_ch)*a_x**(-3.0-3.0*alpha_ch)*(2.0+alpha_ch/2.0) + (1.0 - A_ch)**2*a_x**(-6.0-6.0*alpha_ch) ))

    da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

#ifdef debug_Chalpygin
    write(*,*) 'Called Chaplygin expansion:'
    write(*,*) 'with A_ch =',A_ch, 'alpha_ch =',alpha_ch
    write(*,*)'a_init =',a0, 'dt =',dt0
    write(*,*)'da1 =',da1, 'da2 =',da2
#endif


  end subroutine Chaplygin

