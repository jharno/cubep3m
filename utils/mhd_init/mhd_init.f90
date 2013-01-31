!! cic_power.f90 Parallelized: Hugh Merz Jun 15, 2005
!! This version is used to calculate the power spectrum of the initial conditions
!! Compile with:mpif77 -shared-intel -fpp -g -O0 -CB -fpe0 -xT -DBINARY -mt_mpi mhd_init.f90 -o mhd_init -lm -ldl


program mhd_init
  implicit none
  include 'mpif.h'
  
! frequently changed parameters are found in this header file:
  include '../../parameters'
  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np= hc
  real, parameter    :: npr=np

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: np_node_dim = np/nodes_dim
  integer(4), parameter :: np_buffer = np_node_dim**3
  integer(4), parameter :: max_np = np_node_dim**3 + np_buffer
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim
  integer(4), parameter :: nc_slab = nc / nodes

  !! parallelization variables
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr
  real, dimension(2,nc) :: pkdm
  real, dimension(3,nc) :: pktsum,pksum

  integer(4) :: np_local

  integer(8) :: plan, iplan

  logical :: firstfftw

! :: simulation variables
 
  !! Other parameters
  real, parameter :: pi=3.14159

#ifdef KAISER
  real, parameter :: a = 1/(1+z_i)
#endif

  !! Dark matter arrays
  real, dimension(6,max_np) :: xvp
  real, dimension(3,np_buffer) :: xp_buf
  real, dimension(3*np_buffer) :: send_buf, recv_buf
  

!! MHD arrays
  real, dimension(5,nc_node_dim,nc_node_dim,nc_node_dim) :: u
  real, dimension(3,nc_node_dim,nc_node_dim,nc_node_dim) :: b
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: den
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: den_buf
  !real, dimension(nc+2,nc,nc_slab) :: slab, slab_work
  !real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
  !real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
  
  !! Common block

  common xvp,u, xp_buf, send_buf, recv_buf


!!---start main--------------------------------------------------------------!!

  call mpi_initialize
!  if (rank == 0) call writeparams
  firstfftw=.true.  ! initialize fftw so that it generates the plans
  call initvar
  call read_particles
  call pass_particles
  call darkmatter
  !if (rank == 0) call writepowerspectra
  !call cp_fftw(0)
  call mpi_finalize(ierr)

contains

!!---------------------------------------------------------------------------!!

  subroutine mpi_initialize
    implicit none
    
    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder
  
!! set up global mpi communicator

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (nodes_returned /= nodes ) then
      write(*,*) 'cic_init_power compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'cic_init_power nodes=',nodes 
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    if (mod(nc,nodes) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into slabs'
      write(*,*) 'nc=',nc,'nodes=',nodes,'mod(nc,nodes) != 0'
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    if (rank==0) then
      write(*,*) 'mhd_init running on',nodes,'nodes'
      write(*,*) 'using cubic distribution:',nodes_dim,'nodes per dimension'
      write(*,*) nc,'cells in mesh'
    endif

!! calculate coordinates within slab for cube processes

    slab_coord(3) = rank / nodes_slab
    slab_rank = rank - slab_coord(3) * nodes_slab
    slab_coord(2) = slab_rank / nodes_dim
    slab_coord(1) = slab_rank - slab_coord(2) * nodes_dim
    do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_neighbor(i,j) = i + j * nodes_dim + slab_coord(3) &
                           * nodes_slab
      enddo
    enddo

!! create cartesian communicator based on cubic decomposition

    dims(:) = nodes_dim
    periodic(:) = .true.
    reorder = .false.
    ndim = 3

    call mpi_cart_create(mpi_comm_world, ndim,dims, periodic, &
                       reorder, mpi_comm_cart, ierr)
    call mpi_comm_rank(mpi_comm_cart, cart_rank, ierr)
    call mpi_cart_coords(mpi_comm_cart, cart_rank, ndim,  &
                         cart_coords, ierr)

! cart_neighbor(1) -> down (negative z)
! cart_neighbor(2) -> up (positive z)
! cart_neighbor(3) -> back (negative y)
! cart_neighbor(4) -> front (positive y)
! cart_neighbor(5) -> left (negative x)
! cart_neighbor(6) -> right (positive x)

    do i = 0, ndim-1
      call mpi_cart_shift(mpi_comm_cart, i, 1, cart_neighbor(2*(i+1)-1), &
                          cart_neighbor(2*(i+1)), ierr)
    enddo

#ifdef DEBUG_LOW
  do i=0,nodes-1
    if (i==rank) write(*,'(8i4)') rank,cart_rank,cart_neighbor
    call mpi_barrier(mpi_comm_world,ierr)
  enddo
#endif

  end subroutine mpi_initialize

!!---------------------------------------------------------------------------!!

  subroutine read_particles
    implicit none
    
    real z_write,np_total
    integer j,fstat
    character*512 :: fn
    character(len=6) :: rank_s
#ifdef DEBUG
    integer :: i,pe
    real*8 :: xva(6)
    real*4 :: dmin,dmax
#endif

    if (rank==0) then
      print *,'reading initial conditions'
    endif

    write(rank_s,'(i6)') rank 

    rank_s=adjustl(rank_s)
    fn=scratch_path//'xv'//rank_s(1:len_trim(rank_s))//'.ic'
#ifdef BINARY
    open(21,file=fn,status='old',form='binary',iostat=fstat)
#else
    open(21,file=fn,status='old',form='unformatted',iostat=fstat)
#endif
     if (fstat /= 0) then
      print *,'error opening:',fn
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    read(21) np_local

    if (np_local > max_np) then
      write(*,*) 'too many particles to store'
      write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!! tally up total number of particles
    call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
                         mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'number of particles =', int(np_total,8)

#ifdef DEBUG
    xva=0.
    dmin=1000.0
    dmax=-1000.0
#endif

    do j=1,np_local
      read(21) xvp(:,j)
#ifdef DEBUG
      do i=1,3
        if (xvp(i,j)<dmin) then
          dmin=xvp(i,j)
          pe=j
        endif
        if (xvp(i,j)>dmax) dmax=xvp(i,j)
      enddo
      xva=real(xvp(:,j),kind=8)+xva
#endif
    enddo
    close(21)


#ifdef DEBUG
    do j=0,nodes-1
      if (rank==j) then
        print *,rank,'averages:',xva/real(np_local)
        print *,rank,'min',dmin,'max',dmax
        print *,rank,'bad particle',xvp(:,pe)
      endif
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif
 
  end subroutine read_particles

!!------------------------------------------------------------------!!

  subroutine darkmatter
    implicit none
    integer :: i,j,k
    integer :: i1,j1,k1,fstat
    real    :: d,dmin,dmax,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart,sum_dm_local,sum_dm
    real, dimension(3) :: dis
    character(len=4) :: rank_string
    character(len=100) :: check_name
    real time1,time2
    call cpu_time(time1)

     
    call GetU
    !! call GetB

 
    if (rank==0) then
       print *,'Wrinting u and b to file'
    endif
 
    write(rank_string,'(i4)') rank
    rank_string=adjustl(rank_string)
  
    check_name=output_path//'mhd_ic'// &
               rank_string(1:len_trim(rank_string))//'.dat'

!! open and write density file   
#ifdef BINARY
    open(unit=21,file=check_name,status='replace',iostat=fstat,form='binary')
#else
    open(unit=21,file=check_name,status='replace',iostat=fstat,form='unformatted')
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening density file'
      write(*,*) 'rank',rank,'file:',check_name
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    write(21) u
    write(21) b

!#endif


    sum_dm_local=0.0
    do k=1,nc_node_dim
      do j=1,nc_node_dim
        do i=1,nc_node_dim
          sum_dm_local=real(u(1,i,j,k),kind=8)+sum_dm_local
        enddo
      enddo
    enddo
    call mpi_reduce(sum_dm_local,sum_dm,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) print *,'DM total mass=',sum_dm

    !! Convert dm density field to delta field
    dmin=0
    dmax=0
    dsum=0
    dvar=0

    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim
             !u(i,j,k)=u(i,j,k)-1.0
             d=u(1,i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo

    call mpi_reduce(dsum,dsumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dvar,dvart,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dmin,dmint,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
    call mpi_reduce(dmax,dmaxt,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

    if (rank==0) then
      dsum=dsumt/real(nc)**3
      dvar=sqrt(dvart/real(nc)**3)
      write(*,*)
      write(*,*) 'DM min    ',dmint
      write(*,*) 'DM max    ',dmaxt
      write(*,*) 'Delta sum ',real(dsum)
      write(*,*) 'Delta var ',real(dvar)
      write(*,*)
    endif


    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine darkmatter


  subroutine pass_particles
    implicit none

    integer i,pp,np_buf,np_exit,np_final,npo,npi
    real x(3),lb,ub
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03

    lb=0.0
    ub=real(nc_node_dim)

    np_buf=0
    pp=1
    do
      if (pp > np_local) exit
      x=xvp(:3,pp)
      if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
          x(3) < lb .or. x(3) >= ub ) then
!        write (*,*) 'PARTICLE OUT',xv(:,pp)
        np_buf=np_buf+1
        if (np_buf > np_buffer) then
          print *,rank,'np_buffer =',np_buffer,'exceeded - np_buf =',np_buf
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif 
        xp_buf(:,np_buf)=xvp(:3,pp)
        xvp(:,pp)=xvp(:,np_local)
        np_local=np_local-1
        cycle 
      endif
      pp=pp+1
    enddo
 
    call mpi_reduce(np_buf,np_exit,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr) 

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'np_exit=',np_buf
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    if (rank == 0) print *,'total exiting particles =',np_exit

! pass +x

    tag=11 
    npo=0
    pp=1

    do 
      if (pp > np_buf) exit
      if (xp_buf(1,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'np_out=',npo
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(6), &
                              tag,cart_neighbor(5),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(1,np_buf+pp)=max(xp_buf(1,np_buf+pp)-ub,lb)
    enddo

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'x+ np_local=',np_local
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then

        np_local=np_local+1
        xvp(:3,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'x+ np_exit=',np_buf,np_local
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

! pass -x

    tag=12
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(1,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle 
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(5), &
                              tag,cart_neighbor(6),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(1,np_buf+pp)=min(xp_buf(1,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
  
    np_buf=np_buf+npi

! pass +y

    tag=13 
    npo=0
    pp=1
    do 
      if (pp > np_buf) exit
      if (xp_buf(2,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle 
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(4), &
                              tag,cart_neighbor(3),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)


    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(2,np_buf+pp)=max(xp_buf(2,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

! pass -y

    tag=14
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(2,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(3), &
                              tag,cart_neighbor(4),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,rrequest,rierr)

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(2,np_buf+pp)=min(xp_buf(2,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
  
    np_buf=np_buf+npi

! pass +z

    tag=15 
    npo=0
    pp=1
    do 
      if (pp > np_buf) exit
      if (xp_buf(3,pp) >= ub) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle 
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(2), &
                              tag,cart_neighbor(1),tag,mpi_comm_world, &
                              status,ierr) 

    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)

    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(3,np_buf+pp)=max(xp_buf(3,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do 
      if (pp > npi) exit 
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
   
    np_buf=np_buf+npi

! pass -z

    tag=16
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(3,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*3+1:npo*3)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(1), &
                              tag,cart_neighbor(2),tag,mpi_comm_world, &
                              status,ierr)


    call mpi_isend(send_buf,npo*3,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*3,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*3+1:pp*3)
      xp_buf(3,np_buf+pp)=min(xp_buf(3,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:3,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle 
      endif
      pp=pp+1
    enddo
  
    np_buf=np_buf+npi

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'particles left in buffer=',np_buf
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    call mpi_reduce(np_buf,np_exit,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) print *,'total buffered particles =',np_exit

    call mpi_reduce(np_local,np_final,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) then
      print *,'total particles =',np_final
      if (np_final /= np**3) then
        print *,'ERROR: total number of particles incorrect after passing'
      endif

    endif
 
!!  Check for particles out of bounds

    do i=1,np_local
      if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
          xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
          xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
        print *,'particle out of bounds',rank,i,xvp(:3,i),nc_node_dim
      endif
    enddo

  end subroutine pass_particles



  subroutine GetU
    implicit none
    real, parameter :: mp=(ncr/np)**3!*omega_b/omega_m
    real, parameter ::T_CMB = 2.725 ! in K
    real, parameter ::gamma = 5./3.
    real, parameter :: k_B = 1.38065E-23 ! in J/K
    real,parameter  :: h = 0.701
    real(8), parameter :: UnitConversion = (1+z_i)**(-5.0)*ncr**2/(box*h)**2/omega_m/4.2302E-16
    real, parameter :: a_i = 1/(1+z_i)
    real, parameter :: mu = 1.22 ! reduced mass
    real, parameter :: mproton = 1.6726E-27 ! in kg
    real :: Nprime, Ephys2sim, Econst
 
    integer :: i,i1,i2,j1,j2,k1,k2
    real    :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2,vf,v(3)
    real    :: E_thermal

    E_thermal=0.

    if (z_i > 150.) then

        !! Nprime is the number of physical particles represented by each sim
        !! particle
        !! Ephys2sim converts physical energy units (Joules) to simulation units
        !! Econst stores the remaing numerical factors from Nprime and Ephys2sim

        Econst = (4. / 9.) * 1.e-10
        Nprime = omega_b * box**3 / mu / mproton / ncr**3
        Ephys2sim = a_i**2 * ncr**5 / omega_m**2 / box**5

#ifdef CMB_coupling
        E_thermal = Econst * Nprime * k_B * T_CMB * (1. + z_i) * Ephys2sim
#endif
    endif

    do i=1,np_local
       x=xvp(1,i)!-0.5
       y=xvp(2,i)!-0.5
       z=xvp(3,i)!-0.5
       v=xvp(4:6,i)

       i1=floor(x)+1
       i2=i1+1
       dx1=i1-x
       dx2=1-dx1
       j1=floor(y)+1
       j2=j1+1
       dy1=j1-y
       dy2=1-dy1
       k1=floor(z)+1
       k2=k1+1
       dz1=k1-z
       dz2=1-dz1

       if (i1 < 1 .or. i2 > nc_node_dim .or. j1 < 1 .or. &
           j2 > nc_node_dim .or. k1 < 1 .or. k2 > nc_node_dim) then 
         print *,'particle out of bounds',i1,i2,j1,j2,k1,k2,nc_node_dim
         cycle
       endif 

       dz1=mp*dz1
       dz2=mp*dz2
      
       u(1,i1,j1,k1)=u(1,i1,j1,k1)+dx1*dy1*dz1
       u(1,i2,j1,k1)=u(1,i2,j1,k1)+dx2*dy1*dz1
       u(1,i1,j2,k1)=u(1,i1,j2,k1)+dx1*dy2*dz1
       u(1,i2,j2,k1)=u(1,i2,j2,k1)+dx2*dy2*dz1
       u(1,i1,j1,k2)=u(1,i1,j1,k2)+dx1*dy1*dz2
       u(1,i2,j1,k2)=u(1,i2,j1,k2)+dx2*dy1*dz2
       u(1,i1,j2,k2)=u(1,i1,j2,k2)+dx1*dy2*dz2
       u(1,i2,j2,k2)=u(1,i2,j2,k2)+dx2*dy2*dz2
 
       u(2,i1,j1,k1)=u(2,i1,j1,k1)+dx1*dy1*dz1*v(1)
       u(2,i2,j1,k1)=u(2,i2,j1,k1)+dx2*dy1*dz1*v(1)
       u(2,i1,j2,k1)=u(2,i1,j2,k1)+dx1*dy2*dz1*v(1)
       u(2,i2,j2,k1)=u(2,i2,j2,k1)+dx2*dy2*dz1*v(1)
       u(2,i1,j1,k2)=u(2,i1,j1,k2)+dx1*dy1*dz2*v(1)
       u(2,i2,j1,k2)=u(2,i2,j1,k2)+dx2*dy1*dz2*v(1)
       u(2,i1,j2,k2)=u(2,i1,j2,k2)+dx1*dy2*dz2*v(1)
       u(2,i2,j2,k2)=u(2,i2,j2,k2)+dx2*dy2*dz2*v(1)
 
       u(3,i1,j1,k1)=u(3,i1,j1,k1)+dx1*dy1*dz1*v(2)
       u(3,i2,j1,k1)=u(3,i2,j1,k1)+dx2*dy1*dz1*v(2)
       u(3,i1,j2,k1)=u(3,i1,j2,k1)+dx1*dy2*dz1*v(2)
       u(3,i2,j2,k1)=u(3,i2,j2,k1)+dx2*dy2*dz1*v(2)
       u(3,i1,j1,k2)=u(3,i1,j1,k2)+dx1*dy1*dz2*v(2)
       u(3,i2,j1,k2)=u(3,i2,j1,k2)+dx2*dy1*dz2*v(2)
       u(3,i1,j2,k2)=u(3,i1,j2,k2)+dx1*dy2*dz2*v(2)
       u(3,i2,j2,k2)=u(3,i2,j2,k2)+dx2*dy2*dz2*v(2)

       u(4,i1,j1,k1)=u(4,i1,j1,k1)+dx1*dy1*dz1*v(3)
       u(4,i2,j1,k1)=u(4,i2,j1,k1)+dx2*dy1*dz1*v(3)
       u(4,i1,j2,k1)=u(4,i1,j2,k1)+dx1*dy2*dz1*v(3)
       u(4,i2,j2,k1)=u(4,i2,j2,k1)+dx2*dy2*dz1*v(3)
       u(4,i1,j1,k2)=u(4,i1,j1,k2)+dx1*dy1*dz2*v(3)
       u(4,i2,j1,k2)=u(4,i2,j1,k2)+dx2*dy1*dz2*v(3)
       u(4,i1,j2,k2)=u(4,i1,j2,k2)+dx1*dy2*dz2*v(3)
       u(4,i2,j2,k2)=u(4,i2,j2,k2)+dx2*dy2*dz2*v(3)
  
       u(5,i1,j1,k1)=u(5,i1,j1,k1)+dx1*dy1*dz1*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i2,j1,k1)=u(5,i2,j1,k1)+dx2*dy1*dz1*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i1,j2,k1)=u(5,i1,j2,k1)+dx1*dy2*dz1*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i2,j2,k1)=u(5,i2,j2,k1)+dx2*dy2*dz1*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i1,j1,k2)=u(5,i1,j1,k2)+dx1*dy1*dz2*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i2,j1,k2)=u(5,i2,j1,k2)+dx2*dy1*dz2*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i1,j2,k2)=u(5,i1,j2,k2)+dx1*dy2*dz2*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
       u(5,i2,j2,k2)=u(5,i2,j2,k2)+dx2*dy2*dz2*(v(1)**2+v(2)**2+v(3)**2)/2.0  + E_thermal
    enddo
       !u=u*0.000000001

    return
  end subroutine GetU

  subroutine GetB
    implicit none
    real, parameter :: mp=(ncr/np)**3!*omega_b/omega_m
    real, parameter ::T_CMB = 2.725 ! in K
    real, parameter ::gamma = 5./3.
    real, parameter :: k_B = 1.38065E-23 ! in J/K
    real,parameter  :: h = 0.701
    real(8), parameter :: UnitConversion = (1+z_i)**(-5.0)*ncr**2/(box*h)**2/omega_m/4.2302E-16

    integer :: i,i1,i2,j1,j2,k1,k2,  ix,iy,iz
    real    :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2,vf,v(3)
    real    :: E_thermal

    E_thermal = 0! k_B*T_CMB*(1+z_i)*UnitConversion
    print*, 'Generate initial B field'
    print*, 'nc_node_dim=',nc_node_dim
    do ix=1,nc_node_dim
       do iy=1,nc_node_dim
          do iz=1,nc_node_dim
             if (sqrt((ix-nc_node_dim/2)**2.+(iy-nc_node_dim/2)**2. &
                +(iz-nc_node_dim/2)**2.)<nc_node_dim/8.) then
             b(3,ix,iy,iz)=0.1;
             endif
          enddo
       enddo
    enddo
    !b(1,:,:,:)=0.04
    !b(2,:,:,:)=0.05
    !b(3,:,:,:)=0.06
    !print*, 'B(x)=',b(1,30:31,30:31,30:31)
    !print*, 'B(y)=',b(2,30:31,30:31,30:31)
    !print*, 'B(z)=',b(3,30:31,30:31,30:31)
    print*, 'B(z)=',b(3,:,nc_node_dim/2, nc_node_dim/2)

  end subroutine GetB


!!------------------------------------------------------------------!!

  subroutine initvar
    implicit none
    integer :: k

    real time1,time2
    call cpu_time(time1)

    do k=1,max_np
       xvp(:,k)=0
    enddo
    do k=1,nc_node_dim
       u(:,:,:,k)=0
    enddo
    do k=1,nc_node_dim
       b(:,:,:,k)=0
    enddo
    
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

end program mhd_init
