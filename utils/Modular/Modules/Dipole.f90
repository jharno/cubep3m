!Module to compute dipole
!Same as Haoran Yu's code except using 
!MPI for parallelization instead of CoArray
!Changes: 
!!now includes 1 node term when running in parallel
!!now includes estimate of variance

!Last updated: July 15, 2016
module Dipole
  use Parameters
  use Variables
  use mMPI
  
contains

  subroutine compute_monopole(denA,denB,rarr,xarr,sarr)
    implicit none
    real, dimension(Ncells,Ncells,Ncells), intent(in) :: denA, denB
    real, dimension(:), intent(out) :: rarr, xarr, sarr
    real, dimension(Ncells+1,Ncells+1,Ncells+1,2) :: g
  
    real, parameter :: r1=1.0, r2=sqrt(2.)/2.0, r3=sqrt(3.)/3.0

    integer, dimension(3), parameter :: cci = (/2,1,0/)
    integer, dimension(3) :: ncoord
    real(8) :: xi1,xi2,xi3 !dipole on node
    real(8) :: xi1g,xi2g,xi3g !mean dipole
    real(8) :: xi1s,xi2s,xi3s,xi1sm,xi2sm,xi3sm,xi1ss,xi2ss,xi3ss !error estimation

    real :: r
    integer :: i,j,k,l,m,n,nn,zip_factor,c,nume,err,lr,rr

    !MPI Stuff
    integer :: mpi_comm_sub, mpi_comm_sub_sub
    integer :: sub_rank

    g=0;

    g(:Ncells,:Ncells,:Ncells,1) = denA
    g(:Ncells,:Ncells,:Ncells,2) = denB

    rarr = 0; xarr = 0; sarr = 0;

    call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.true./),mpi_comm_sub,err )
    if (err/=mpi_success) call dipole_error_stop('Error in creating subvolume communicator')

    call mpi_comm_rank(mpi_comm_sub,sub_rank,err )
    if (err/=mpi_success) call dipole_error_stop('Error in computing subvolume slab rank')

    call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.false./),mpi_comm_sub_sub,err )
    if (err/=mpi_success) call dipole_error_stop('Error in computing sub-sub volume communicator')

    !First compute distributed piece
    nn = nodes_dim
    m = Ncells
    c = 1
    do
       r = Lbox/m/nn
       n = m+1
       xi1=0; xi2=0; xi3=0
       if (rank.eq.0) write(*,*) 'Scale: ',r,m

       !Fill boundary cells
       !!Pass x 
       !if (rank.eq.0) write(*,*) 'Pass x'
       call mpi_cart_shift( mpi_comm_cart,cci(1),1,lr,rr,err )
       if (err/=mpi_success) call dipole_error_stop('Could not determine cart x')
        
       nume = size( g(1,:m,:m,:) )
       if (nume/=size( g(n,:m,:m,:) )) call dipole_error_stop('Error: nume for g 1,n in x do not match')
       call mpi_sendrecv( g(1,:m,:m,:),nume,mpi_real,lr,rank**2,g(n,:m,:m,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv g in x')

       !if (rank.eq.0) write(*,*) 'Pass y'
       call mpi_cart_shift( mpi_comm_cart,cci(2),1,lr,rr,err )
       if (err/=mpi_success) call dipole_error_stop('Could not determine cart y')

       nume = size( g(:n,1,:m,:) )
       if (nume/=size( g(:n,n,:m,:) )) call dipole_error_stop('Error: nume for g 1,n in y do not match')
       call mpi_sendrecv( g(:n,1,:m,:),nume,mpi_real,lr,rank**2,g(:n,n,:m,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv g in y')

       !!Pass z
       !if (rank.eq.0) write(*,*) 'Pass z'
       call mpi_cart_shift( mpi_comm_cart,cci(3),1,lr,rr,err )
       if (err/=mpi_success) call dipole_error_stop('Could not determine cart z')

       nume = size( g(:n,:n,1,:) )
       if (nume/=size( g(:n,:n,n,:) )) call dipole_error_stop('Error: nume for g 1,n in z do not match')
       call mpi_sendrecv( g(:n,:n,1,:),nume,mpi_real,lr,rank**2,g(:n,:n,n,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv g in z')

       !Compute local dipole
       !! 1 displacement
       xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,1:m,2)*1.d0)
       xi1=xi1+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,1:m,2)*1.d0)

       xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,1:m,2)*1.d0)
       xi1=xi1+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,1:m,2)*1.d0)

       xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,1:m,2:n,2)*1.d0)
       xi1=xi1+sum(g(1:m,1:m,2:n,1)*g(1:m,1:m,1:m,2)*1.d0)

       !! 2 displacements
       xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,2:n,2)*1.d0)
       xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(1:m,2:n,1:m,2)*1.d0)
       xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,2:n,2)*1.d0)
       xi2=xi2+sum(g(1:m,2:n,2:n,1)*g(1:m,1:m,1:m,2)*1.d0)

       xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,1:m,2)*1.d0)
       xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,1:m,2)*1.d0)
       xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,1:m,2)*1.d0)
       xi2=xi2+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,1:m,2)*1.d0)

       xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,2:n,2)*1.d0)
       xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(2:n,1:m,1:m,2)*1.d0)
       xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,2:n,2)*1.d0)
       xi2=xi2+sum(g(2:n,1:m,2:n,1)*g(1:m,1:m,1:m,2)*1.d0)
       

       !! 3 displacements
       xi3=xi3+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,2:n,2)*1.d0)

       xi3=xi3+sum(g(1:m,1:m,2:n,1)*g(2:n,2:n,1:m,2)*1.d0)
       xi3=xi3+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,2:n,2)*1.d0)
       xi3=xi3+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,2:n,2)*1.d0)

       xi3=xi3+sum(g(1:m,2:n,2:n,1)*g(2:n,1:m,1:m,2)*1.d0)
       xi3=xi3+sum(g(2:n,1:m,2:n,1)*g(1:m,2:n,1:m,2)*1.d0)
       xi3=xi3+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,2:n,2)*1.d0)

       xi3=xi3+sum(g(2:n,2:n,2:n,1)*g(1:m,1:m,1:m,2)*1.d0)

       !Local average
       xi1=xi1/m**3/6.0
       xi2=xi2/m**3/12.0
       xi3=xi3/m**3/8.0
       
       !Reduce to one node
       xi1g = 0; xi2g = 0; xi3g = 0;

       call mpi_reduce( xi1,xi1g,1,mpi_real8,mpi_sum,0,mpi_comm_world,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1g reduction')
       call mpi_reduce( xi2,xi2g,1,mpi_real8,mpi_sum,0,mpi_comm_world,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2g reduction')
       call mpi_reduce( xi3,xi3g,1,mpi_real8,mpi_sum,0,mpi_comm_world,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3g reduction')       

       !Global average
       xi1g=xi1g/nn**3
       xi2g=xi2g/nn**3
       xi3g=xi3g/nn**3

       !Compute fluctuations
       xi1s = 0; xi2s = 0; xi3s = 0;
       call mpi_reduce( xi1,xi1s,1,mpi_real8,mpi_sum,0,mpi_comm_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1s reduction')
       call mpi_reduce( xi2,xi2s,1,mpi_real8,mpi_sum,0,mpi_comm_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2s reduction')
       call mpi_reduce( xi3,xi3s,1,mpi_real8,mpi_sum,0,mpi_comm_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3s reduction')

       !Slab average
       xi1s=xi1s/nn**2
       xi2s=xi2s/nn**2
       xi3s=xi3s/nn**2

       xi1sm = 0; xi2sm = 0; xi3sm = 0;
       call mpi_reduce( xi1s,xi1sm,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1sm reduction')
       call mpi_reduce( xi2s,xi2sm,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2sm reduction')
       call mpi_reduce( xi3s,xi3sm,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3sm reduction')
       
       !Pencil average
       xi1sm=xi1sm/nn
       xi2sm=xi2sm/nn
       xi3sm=xi3sm/nn
       if ( abs(xi1sm-xi1g) .gt. 1e-5 ) call dipole_error_stop('Error in sm g comparison')

       xi1ss = 0; xi2ss = 0; xi3ss = 0;
       call mpi_reduce( xi1s**2.0,xi1ss,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1ss reduction')
       call mpi_reduce( xi2s**2.0,xi2ss,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2ss reduction')
       call mpi_reduce( xi3s**2.0,xi3ss,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3ss reduction')
       
       !Pencil average
       xi1ss=xi1ss/nn
       xi2ss=xi2ss/nn
       xi3ss=xi3ss/nn

       !Save
       if (rank.eq.0) then
          rarr(c) = r
          xarr(c) = xi1g
          sarr(c) = (xi1ss-xi1sm**2)**0.5
          c=c+1
          rarr(c) = r*sqrt(2.0)
          xarr(c) = xi2g
          sarr(c) = (xi2ss-xi2sm**2)**0.5
          c=c+1
          rarr(c) = r*sqrt(3.0)
          xarr(c) = xi3g
          sarr(c) = (xi3ss-xi3sm**2)**0.5
          c=c+1
       end if

       zip_factor = f(m)
       if (m.lt.zip_factor .or. zip_factor.eq.-1) exit
       
       call zip(g(:m,:m,:m,1),zip_factor)
       call zip(g(:m,:m,:m,2),zip_factor)

       m=m/zip_factor

    end do

    !One node piece
    if (Nmpi.gt.1) write(*,*) 'WARNING: 1 node piece not implemented in cf_2pt yet'
    
    nume=size(rarr)
    call mpi_bcast(rarr,nume,mpi_real,0,mpi_comm_world,err)
    if (err/=mpi_success) call dipole_error_stop('Error in mpi_bcast rarr')
    call mpi_bcast(xarr,nume,mpi_real,0,mpi_comm_world,err)
    if (err/=mpi_success) call dipole_error_stop('Error in mpi_bcast xarr')
    call mpi_bcast(sarr,nume,mpi_real,0,mpi_comm_world,err)
    if (err/=mpi_success) call dipole_error_stop('Error in mpi_bcast sarr')

  end subroutine compute_monopole

  subroutine compute_dipole(denA,denB,velR,rarr,xarr,sarr)
    implicit none
    real, dimension(Ncells,Ncells,Ncells), intent(in) :: denA, denB
    real, dimension(Ncells,Ncells,Ncells,3), intent(in) :: velR
    real, dimension(:), intent(out) :: rarr,xarr,sarr
    
    real, dimension(Ncells+1,Ncells+1,Ncells+1,2) :: g
    real, dimension(Ncells+1,Ncells+1,Ncells+1,3) :: v

    real, parameter :: r1=1.0, r2=sqrt(2.)/2.0, r3=sqrt(3.)/3.0
    character(1), dimension(3), parameter :: xyz = (/'x','y','z'/)
    
    integer, dimension(3), parameter :: cci = (/2,1,0/)
    integer, dimension(3) :: ncoord
    real(8) :: xi1,xi2,xi3 !dipole on node
    real(8) :: xi1g,xi2g,xi3g !mean dipole
    real(8) :: xi1s,xi2s,xi3s,xi1sm,xi2sm,xi3sm,xi1ss,xi2ss,xi3ss !error estimation

    real :: r
    integer :: i,j,k,l,m,n,nn,zip_factor,c,nume,err,lr,rr

    !MPI Stuff
    integer :: mpi_comm_sub, mpi_comm_sub_sub
    integer :: sub_rank

    g=0; v = 0;

    g(:Ncells,:Ncells,:Ncells,1) = denA
    g(:Ncells,:Ncells,:Ncells,2) = denB

    v(:Ncells,:Ncells,:Ncells,:) = velR

    rarr = 0; xarr = 0; sarr = 0

    call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.true./),mpi_comm_sub,err )
    if (err/=mpi_success) call dipole_error_stop('Error in creating subvolume communicator')

    call mpi_comm_rank(mpi_comm_sub,sub_rank,err )
    if (err/=mpi_success) call dipole_error_stop('Error in computing subvolume slab rank')

    call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.false./),mpi_comm_sub_sub,err )
    if (err/=mpi_success) call dipole_error_stop('Error in computing sub-sub volume communicator')

    !First compute distributed piece
    nn = nodes_dim
    m = Ncells
    c = 1
    do
       r = Lbox/m/nn
       n = m+1
       xi1=0; xi2=0; xi3=0
       if (rank.eq.0) write(*,*) 'Scale: ',r,m

       !Fill boundary cells
       !!Pass x 
       !if (rank.eq.0) write(*,*) 'Pass x'
       call mpi_cart_shift( mpi_comm_cart,cci(1),1,lr,rr,err )
       if (err/=mpi_success) call dipole_error_stop('Could not determine cart x')
        
       nume = size( g(1,:m,:m,:) )
       if (nume/=size( g(n,:m,:m,:) )) call dipole_error_stop('Error: nume for g 1,n in x do not match')
       call mpi_sendrecv( g(1,:m,:m,:),nume,mpi_real,lr,rank**2,g(n,:m,:m,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv g in x')

       nume = size( v(1,:m,:m,:) )
       if (nume/=size( v(n,:m,:m,:) )) call dipole_error_stop('Error: nume for v 1,n in x do not match')
       call mpi_sendrecv( v(1,:m,:m,:),nume,mpi_real,lr,rank**2,v(n,:m,:m,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Could not send v in x')

       !!Pass y
       !if (rank.eq.0) write(*,*) 'Pass y'
       call mpi_cart_shift( mpi_comm_cart,cci(2),1,lr,rr,err )
       if (err/=mpi_success) call dipole_error_stop('Could not determine cart y')

       nume = size( g(:n,1,:m,:) )
       if (nume/=size( g(:n,n,:m,:) )) call dipole_error_stop('Error: nume for g 1,n in y do not match')
       call mpi_sendrecv( g(:n,1,:m,:),nume,mpi_real,lr,rank**2,g(:n,n,:m,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv g in y')

       nume = size( v(:n,1,:m,:) )
       if (nume/=size( v(:n,n,:m,:) )) call dipole_error_stop('Error: nume for v 1,n in y do not match')
       call mpi_sendrecv( v(:n,1,:m,:),nume,mpi_real,lr,rank**2,v(:n,n,:m,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv v in y')

       !!Pass z
       !if (rank.eq.0) write(*,*) 'Pass z'
       call mpi_cart_shift( mpi_comm_cart,cci(3),1,lr,rr,err )
       if (err/=mpi_success) call dipole_error_stop('Could not determine cart z')

       nume = size( g(:n,:n,1,:) )
       if (nume/=size( g(:n,:n,n,:) )) call dipole_error_stop('Error: nume for g 1,n in z do not match')
       call mpi_sendrecv( g(:n,:n,1,:),nume,mpi_real,lr,rank**2,g(:n,:n,n,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv g in z')

       nume = size( v(:n,:n,1,:) )
       if (nume/=size( v(:n,:n,1,:) )) call dipole_error_stop('Error: nume for v 1,n in z do not match')
       call mpi_sendrecv( v(:n,:n,1,:),nume,mpi_real,lr,rank**2,v(:n,:n,n,:),nume,mpi_real,rr,rr**2,mpi_comm_world,mpi_status_ignore,err )
       if (err/=mpi_success) call dipole_error_stop('Error in send_recv v in z')

       !Compute local dipole
       !! sqrt(1)
       !if (rank.eq.0) write(*,*) 'root 1'
       xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,1:m,1:m,:)),(/r1,0.,0./))*1.d0)
       xi1=xi1+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,1:m,1:m,:)),(/-r1,0.,0./))*1.d0)
       xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,2:n,1:m,:)),(/0.,r1,0./))*1.d0)
       xi1=xi1+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(1:m,1:m,1:m,:)),(/0.,-r1,0./))*1.d0)
       xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,1:m,2:n,:)),(/0.,0.,r1/))*1.d0)
       xi1=xi1+sum(g(1:m,1:m,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(1:m,1:m,1:m,:)),(/0.,0.,-r1/))*1.d0)
       !! sqrt(2)
       !if (rank.eq.0) write(*,*) 'root 2'
       xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,2:n,2:n,:)),(/0.,r2,r2/))*1.d0)
       xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(1:m,1:m,2:n,:)),(/0.,-r2,r2/))*1.d0)
       xi2=xi2+sum(g(1:m,2:n,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,2:n,2:n,:)+v(1:m,1:m,1:m,:)),(/0.,-r2,-r2/))*1.d0)
       xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(1:m,2:n,1:m,:)),(/0.,r2,-r2/))*1.d0)
       xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,1:m,2:n,:)),(/r2,0.,r2/))*1.d0)
       xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,1:m,2:n,:)),(/-r2,0.,r2/))*1.d0)
       xi2=xi2+sum(g(2:n,1:m,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,1:m,2:n,:)+v(1:m,1:m,1:m,:)),(/-r2,0.,-r2/))*1.d0)
       xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(2:n,1:m,1:m,:)),(/r2,0.,-r2/))*1.d0)
       xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,2:n,1:m,:)),(/r2,r2,0./))*1.d0)
       xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,2:n,1:m,:)),(/-r2,r2,0./))*1.d0)
       xi2=xi2+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,2:n,1:m,:)+v(1:m,1:m,1:m,:)),(/-r2,-r2,0./))*1.d0)
       xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(2:n,1:m,1:m,:)),(/r2,-r2,0./))*1.d0)
       !! sqrt(3)  
       !if (rank.eq.0) write(*,*) 'root 3'
       xi3=xi3+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,2:n,2:n,:)),(/r3,r3,r3/))*1.d0)
       xi3=xi3+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,2:n,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,2:n,2:n,:)),(/-r3,r3,r3/))*1.d0)
       xi3=xi3+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,2:n,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(2:n,1:m,2:n,:)),(/r3,-r3,r3/))*1.d0)
       xi3=xi3+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(2:n,2:n,1:m,:)+v(1:m,1:m,2:n,:)),(/-r3,-r3,r3/))*1.d0)
       xi3=xi3+sum(g(1:m,1:m,2:n,1)*g(2:n,2:n,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(2:n,2:n,1:m,:)),(/r3,r3,-r3/))*1.d0)
       xi3=xi3+sum(g(2:n,1:m,2:n,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(2:n,1:m,2:n,:)+v(1:m,2:n,1:m,:)),(/-r3,r3,-r3/))*1.d0)
       xi3=xi3+sum(g(1:m,2:n,2:n,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,2:n,2:n,:)+v(2:n,1:m,1:m,:)),(/r3,-r3,-r3/))*1.d0)
       xi3=xi3+sum(g(2:n,2:n,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,2:n,2:n,:)+v(1:m,1:m,1:m,:)),(/-r3,-r3,-r3/))*1.d0)

       !Local average
       xi1=xi1/m**3/6
       xi2=xi2/m**3/12
       xi3=xi3/m**3/8

       !Reduce to one node 
       xi1g = 0; xi2g = 0; xi3g = 0;
       call mpi_reduce( xi1,xi1g,1,mpi_real8,mpi_sum,0,mpi_comm_world,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1g reduction')
       call mpi_reduce( xi2,xi2g,1,mpi_real8,mpi_sum,0,mpi_comm_world,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2g reduction')
       call mpi_reduce( xi3,xi3g,1,mpi_real8,mpi_sum,0,mpi_comm_world,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3g reduction')

       !Global average
       xi1g=xi1g/nn**3
       xi2g=xi2g/nn**3
       xi3g=xi3g/nn**3

       !Compute fluctuations
       xi1s = 0; xi2s = 0; xi3s = 0;
       call mpi_reduce( xi1,xi1s,1,mpi_real8,mpi_sum,0,mpi_comm_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1s reduction')
       call mpi_reduce( xi2,xi2s,1,mpi_real8,mpi_sum,0,mpi_comm_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2s reduction')
       call mpi_reduce( xi3,xi3s,1,mpi_real8,mpi_sum,0,mpi_comm_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3s reduction')

       !Slab average
       xi1s=xi1s/nn**2
       xi2s=xi2s/nn**2
       xi3s=xi3s/nn**2

       xi1sm = 0; xi2sm = 0; xi3sm = 0;
       call mpi_reduce( xi1s,xi1sm,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1sm reduction')
       call mpi_reduce( xi2s,xi2sm,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2sm reduction')
       call mpi_reduce( xi3s,xi3sm,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3sm reduction')
       
       !Pencil average
       xi1sm=xi1sm/nn
       xi2sm=xi2sm/nn
       xi3sm=xi3sm/nn
       if ( abs(xi1sm-xi1g) .gt. 1e-5 ) call dipole_error_stop('Error in sm g comparison')

       xi1ss = 0; xi2ss = 0; xi3ss = 0;
       call mpi_reduce( xi1s**2.0,xi1ss,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi1ss reduction')
       call mpi_reduce( xi2s**2.0,xi2ss,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi2ss reduction')
       call mpi_reduce( xi3s**2.0,xi3ss,1,mpi_real8,mpi_sum,0,mpi_comm_sub_sub,err )
       if (err/=mpi_success) call dipole_error_stop('Error in xi3ss reduction')
       
       !Pencil average
       xi1ss=xi1ss/nn
       xi2ss=xi2ss/nn
       xi3ss=xi3ss/nn

       !Save
       if (rank.eq.0) then
          rarr(c) = r
          xarr(c) = xi1g
          sarr(c) = (xi1ss-xi1sm**2)**0.5
          c=c+1
          rarr(c) = r*sqrt(2.0)
          xarr(c) = xi2g
          sarr(c) = (xi2ss-xi2sm**2)**0.5
          c=c+1
          rarr(c) = r*sqrt(3.0)
          xarr(c) = xi3g
          sarr(c) = (xi3ss-xi3sm**2)**0.5
          c=c+1
       end if
       
       zip_factor = f(m)
       if (m.lt.zip_factor .or. zip_factor.eq.-1) exit
       
       call zip(g(:m,:m,:m,1),zip_factor)
       call zip(g(:m,:m,:m,2),zip_factor)

       call zip(v(:m,:m,:m,1),zip_factor)
       call zip(v(:m,:m,:m,2),zip_factor)
       call zip(v(:m,:m,:m,3),zip_factor)
       
       m=m/zip_factor

    end do

    !Now compute one node piece if parallel
    if (Nmpi.gt.1) then
       !First pass remaining pieces to head node
       if (rank .ne. 0) then
          nume = size(g(:m,:m,:m,:))
          call mpi_send(g(:m,:m,:m,:),nume,mpi_real,0,rank**2,mpi_comm_world,err)
          if (err/=mpi_success) call dipole_error_stop('Error in global send g')
          nume = size(v(:m,:m,:m,:))
          call mpi_send(v(:m,:m,:m,:),nume,mpi_real,0,rank**2,mpi_comm_world,err)
          if (err/=mpi_success) call dipole_error_stop('Error in global send v')
       else
          !write(*,*) 'Computing large scale dipole on root node'
          do n=1,nn**3-1
             call mpi_cart_coords(mpi_comm_cart,n,3,ncoord,err)
             ncoord = ncoord+1 !So it goes from 1->nodes_dim rather than 0->nodes_dim-1
             if (err/=mpi_success) call dipole_error_stop('Error in mpi_cart_coords') 
             nume = size(g(:m,:m,:m,:))
             call mpi_recv(g(1+(ncoord(1+cci(1))-1)*m:ncoord(1+cci(1))*m, &
                  & 1+(ncoord(1+cci(2))-1)*m:ncoord(1+cci(2))*m, &
                  1+(ncoord(1+cci(3))-1)*m:ncoord(1+cci(3))*m,:),&
                  & nume,mpi_real,n,n**2,mpi_comm_world,mpi_status_ignore,err)
             if (err/=mpi_success) call dipole_error_stop('Error in global recv g')
             nume = size(v(:m,:m,:m,:))
             call mpi_recv(v(1+(ncoord(1+cci(1))-1)*m:ncoord(1+cci(1))*m, &
                  & 1+(ncoord(1+cci(2))-1)*m:ncoord(1+cci(2))*m, &
                  1+(ncoord(1+cci(3))-1)*m:ncoord(1+cci(3))*m,:),&
                  & nume,mpi_real,n,n**2,mpi_comm_world,mpi_status_ignore,err)
             if (err/=mpi_success) call dipole_error_stop('Error in global recv v')
             if (n.eq.nn**3-1) m = m*ncoord(1+cci(3)) !New m            
          end do
       end if

       !Now compute dipole just on head node
       if (rank.eq.0) then
          do

             !Reduce first so we don't repeat
             zip_factor = f(m)
             if (m.lt.zip_factor .or. zip_factor.eq.-1) exit
             
             call zip(g(:m,:m,:m,1),zip_factor)
             call zip(g(:m,:m,:m,2),zip_factor)
             
             call zip(v(:m,:m,:m,1),zip_factor)
             call zip(v(:m,:m,:m,2),zip_factor)
             call zip(v(:m,:m,:m,3),zip_factor)
             
             m=m/zip_factor

             r = Lbox/m
             n = m+1
             xi1=0; xi2=0; xi3=0
             !if (rank.eq.0) write(*,*) 'Scale: ',r,m
             
             !Buffer
             g(n,:m,:m,:)=g(1,:m,:m,:)
             g(:n,n,:m,:)=g(:n,1,:m,:)
             g(:n,:n,n,:)=g(:n,:n,1,:)
             v(n,:m,:m,:)=v(1,:m,:m,:)
             v(:n,n,:m,:)=v(:n,1,:m,:)
             v(:n,:n,n,:)=v(:n,:n,1,:)

             !Compute local dipole
             !! sqrt(1)
             !if (rank.eq.0) write(*,*) 'root 1'
             xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,1:m,1:m,:)),(/r1,0.,0./))*1.d0)
             xi1=xi1+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,1:m,1:m,:)),(/-r1,0.,0./))*1.d0)
             xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,2:n,1:m,:)),(/0.,r1,0./))*1.d0)
             xi1=xi1+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(1:m,1:m,1:m,:)),(/0.,-r1,0./))*1.d0)
             xi1=xi1+sum(g(1:m,1:m,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,1:m,2:n,:)),(/0.,0.,r1/))*1.d0)
             xi1=xi1+sum(g(1:m,1:m,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(1:m,1:m,1:m,:)),(/0.,0.,-r1/))*1.d0)
             !! sqrt(2)
             !if (rank.eq.0) write(*,*) 'root 2'
             xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(1:m,2:n,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(1:m,2:n,2:n,:)),(/0.,r2,r2/))*1.d0)
             xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(1:m,1:m,2:n,:)),(/0.,-r2,r2/))*1.d0)
             xi2=xi2+sum(g(1:m,2:n,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(1:m,2:n,2:n,:)+v(1:m,1:m,1:m,:)),(/0.,-r2,-r2/))*1.d0)
             xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(1:m,2:n,1:m,:)),(/0.,r2,-r2/))*1.d0)
             xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,1:m,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,1:m,2:n,:)),(/r2,0.,r2/))*1.d0)
             xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,1:m,2:n,:)),(/-r2,0.,r2/))*1.d0)
             xi2=xi2+sum(g(2:n,1:m,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,1:m,2:n,:)+v(1:m,1:m,1:m,:)),(/-r2,0.,-r2/))*1.d0)
             xi2=xi2+sum(g(1:m,1:m,2:n,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(2:n,1:m,1:m,:)),(/r2,0.,-r2/))*1.d0)
             xi2=xi2+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,1:m,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,2:n,1:m,:)),(/r2,r2,0./))*1.d0)
             xi2=xi2+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,2:n,1:m,:)),(/-r2,r2,0./))*1.d0)
             xi2=xi2+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,2:n,1:m,:)+v(1:m,1:m,1:m,:)),(/-r2,-r2,0./))*1.d0)
             xi2=xi2+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(2:n,1:m,1:m,:)),(/r2,-r2,0./))*1.d0)
             !! sqrt(3)  
             !if (rank.eq.0) write(*,*) 'root 3'
             xi3=xi3+sum(g(1:m,1:m,1:m,1)*g(2:n,2:n,2:n,2)*vdot(hat(v(1:m,1:m,1:m,:)+v(2:n,2:n,2:n,:)),(/r3,r3,r3/))*1.d0)
             xi3=xi3+sum(g(2:n,1:m,1:m,1)*g(1:m,2:n,2:n,2)*vdot(hat(v(2:n,1:m,1:m,:)+v(1:m,2:n,2:n,:)),(/-r3,r3,r3/))*1.d0)
             xi3=xi3+sum(g(1:m,2:n,1:m,1)*g(2:n,1:m,2:n,2)*vdot(hat(v(1:m,2:n,1:m,:)+v(2:n,1:m,2:n,:)),(/r3,-r3,r3/))*1.d0)
             xi3=xi3+sum(g(2:n,2:n,1:m,1)*g(1:m,1:m,2:n,2)*vdot(hat(v(2:n,2:n,1:m,:)+v(1:m,1:m,2:n,:)),(/-r3,-r3,r3/))*1.d0)
             xi3=xi3+sum(g(1:m,1:m,2:n,1)*g(2:n,2:n,1:m,2)*vdot(hat(v(1:m,1:m,2:n,:)+v(2:n,2:n,1:m,:)),(/r3,r3,-r3/))*1.d0)
             xi3=xi3+sum(g(2:n,1:m,2:n,1)*g(1:m,2:n,1:m,2)*vdot(hat(v(2:n,1:m,2:n,:)+v(1:m,2:n,1:m,:)),(/-r3,r3,-r3/))*1.d0)
             xi3=xi3+sum(g(1:m,2:n,2:n,1)*g(2:n,1:m,1:m,2)*vdot(hat(v(1:m,2:n,2:n,:)+v(2:n,1:m,1:m,:)),(/r3,-r3,-r3/))*1.d0)
             xi3=xi3+sum(g(2:n,2:n,2:n,1)*g(1:m,1:m,1:m,2)*vdot(hat(v(2:n,2:n,2:n,:)+v(1:m,1:m,1:m,:)),(/-r3,-r3,-r3/))*1.d0)
             
             rarr(c) = r
             xarr(c) = xi1/m**3/6
             c=c+1
             rarr(c) = r*sqrt(2.0)
             xarr(c) = xi2/m**3/12
             c=c+1
             rarr(c) = r*sqrt(3.0)
             xarr(c) = xi3/m**3/8
             c=c+1
             
          end do
       end if
    
    end if

    nume=size(rarr)
    call mpi_bcast(rarr,nume,mpi_real,0,mpi_comm_world,err)
    if (err/=mpi_success) call dipole_error_stop('Error in mpi_bcast rarr')
    call mpi_bcast(xarr,nume,mpi_real,0,mpi_comm_world,err)
    if (err/=mpi_success) call dipole_error_stop('Error in mpi_bcast xarr')
    call mpi_bcast(sarr,nume,mpi_real,0,mpi_comm_world,err)
    if (err/=mpi_success) call dipole_error_stop('Error in mpi_bcast sarr')

  end subroutine compute_dipole

  function f(x)
    integer :: x,f
    f=-1
    if (mod(x,2)==0) then 
       f=2
    else if (mod(x,3)==0) then
       f=3
    end if
  end function f

  function hat(vec)
    real :: vec(:,:,:,:), hat(size(vec,1),size(vec,2),size(vec,3),size(vec,4))
    hat = vec / spread(sqrt(sum(vec**2,dim=4)),dim=4,ncopies=3)
  end function hat

  function vdot(vec,vn)
    real :: vec(:,:,:,:), vn(3), vdot(size(vec,1),size(vec,2),size(vec,3))
    vdot=vec(:,:,:,1)*vn(1)+vec(:,:,:,2)*vn(2)+vec(:,:,:,3)*vn(3)
  endfunction vdot

  subroutine zip(vec,x)
    integer :: i,x
    real :: vec(:,:,:), vv(size(vec,1)/x,size(vec,2),size(vec,3))
    vv=0
    do i=1,x
       vv=vv+vec(i::x,:,:)
    enddo
    vec(:size(vec,1)/x,:,:)=vv
    vv=0
    do i=1,x
       vv(:,:size(vec,2)/x,:)=vv(:,:size(vec,2)/x,:)+vec(:size(vec,1)/x,i::x,:)
    enddo
    vec(:size(vec,1)/x,:size(vec,2)/x,:)=vv(:,:size(vec,2)/x,:)
    vv=0
    do i=1,x
       vv(:,:,:size(vec,3)/x)=vv(:,:,:size(vec,3)/x)+vec(:size(vec,1)/x,:size(vec,2)/x,i::x)
    enddo
    vec(:size(vec,1)/x,:size(vec,2)/x,:size(vec,3)/x)=vv(:,:size(vec,2)/x,:size(vec,3)/x)/x**3
  end subroutine zip

  subroutine pbc3(vec)
    implicit none
    integer :: n1,n2,n3
    real :: vec(:,:,:)
    n1=size(vec,1); n2=size(vec,2); n3=size(vec,3)
    vec(n1,:,:)=vec(1,:,:)
    vec(:,n2,:)=vec(:,1,:)
    vec(:,:,n3)=vec(:,:,1)
  end subroutine pbc3

  subroutine dipole_error_stop(expl)
    implicit none
    character(len=*) :: expl
    write(*,*) '[Mod - Dipole]'
    write(*,*) '-->'//expl
    call mpi_abort(mpi_comm_world, ierr, ierr)
  end subroutine dipole_error_stop

end module Dipole
