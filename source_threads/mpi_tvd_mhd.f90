!=======================================================================
! 3-D Decomposition MPI TVD MHD Module
!
! written November 2001 by Ue-Li Pen, pen@cita.utoronto.ca
! debugged to pass alven test March 29, 2002
! Optimized routines 05.06.2002 Hugh Merz merz@cita.utoronto.ca
! MPI routines 05.06.2002 Hugh Merz merz@cita.utoronto.ca
! Cleanup 03.2003 Matthias Liebendoerfer liebend@cita.utoronto.ca
! Comm_module & sweep written for cubic decomposition 04.2003 Matthias
! Optimized routines customized (no performance gain) 04.2003 Matthias
! Modularized 03.2005 Hugh
! !$omp description added back 01.2012 Yu
!
! B field is stored on the left side of each cell
!
!=======================================================================
      module mpi_tvd_mhd

      implicit none

      real, parameter :: gamma = 5.0/3.0

!.....type for dimensions, work load distribution, comm handles.........
      type comm_wld
        integer :: g !global number of zones in one direction
        integer :: r !global index of index 0 in local array
        integer :: m,n !start and end of local section without buffers
        integer :: l !dimension of array inclusive buffers
        integer, dimension(4) :: requests !communication handles
      end type comm_wld

!.....mpi information...................................................
      integer :: comm_cart  !communicator used for gas processes
      integer :: comm_cart_rank  !rank of process in comm_cart
      integer, dimension(3) :: comm_cart_coord 
!coordinates of process in comm_cart -- spans (0:nodes_dim-1)
! (1)-- z dir, (2)-- y dir, (3)-- x dir
      integer, dimension(6) :: comm_cart_neighbor    
! rank of neighboring processes in comm_cart
! comm_cart_neighbor(1) -> down (negative z)
! comm_cart_neighbor(2) -> up (positive z)
! comm_cart_neighbor(3) -> back (negative y)
! comm_cart_neighbor(4) -> front (positive y)
! comm_cart_neighbor(5) -> left (negative x)
! comm_cart_neighbor(6) -> right (positive x)

!.....internal parameters...............................................
      integer, parameter :: nr=4 !number of requests per direction
      integer, parameter :: yzbuf=3 !permanent xyz overlap between cubes
      integer, parameter :: xbuf=6 !additional volatile x overlap
      integer, parameter :: nu=5 !dimension of hydrodynamical variable
      integer, parameter :: nb=3 !dimension of magnetic field

!.....message tags......................................................
      integer, parameter :: xtagdw=1
      integer, parameter :: xtagup=2
      integer, parameter :: ytagdw=3
      integer, parameter :: ytagup=4
      integer, parameter :: ztagdw=5
      integer, parameter :: ztagup=6

!.....memory structures.................................................
      type(comm_wld) :: nx,ny,nz
      real, dimension (:,:,:,:), allocatable :: u,b
      real, dimension (:), allocatable :: send0,send1,recv0,recv1 

      contains

!=======================================================================
      subroutine mpi_tvd_mhd_state_output(out_dir,cur_iter,cur_t,z_write)
      implicit none
      
      include 'mpif.h'
 
      integer, intent(in) :: cur_iter
      real, intent(in) :: cur_t
      character(*) :: out_dir
      character(len=80) :: fn1,fn2,fn3
      character(len=7) :: z_write
      integer :: fstat,ierr

      if (comm_cart_rank==0) print *,'writing out state',cur_iter,cur_t,z_write
!      write(fn1,'(i10,"state")'),cur_iter
      fn1=z_write(1:len_trim(z_write))//'mhd'
      write(fn2,'(i10,".dat")') comm_cart_rank
      fn1=adjustl(fn1)
      fn2=adjustl(fn2) 
      fn3=out_dir(1:len_trim(out_dir))//fn1(1:len_trim(fn1))//fn2(1:len_trim(fn2))
#ifdef BINARY
      open(185,file=fn3,form='binary',status='replace',iostat=fstat)
#else
      open(185,file=fn3,form='unformatted',status='replace',iostat=fstat)
#endif
      if (fstat /= 0) then
        write(*,*) 'error opening state output file'
        write(*,*) 'rank',comm_cart_rank,'file:',fn3
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      write(185) cur_iter,cur_t,nx,ny,nz
      write(185) u(:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      write(185) b(:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      close(185)
      
      end subroutine mpi_tvd_mhd_state_output
!=======================================================================
      subroutine mpi_tvd_mhd_ic(ic_dir)
      implicit none

      include 'mpif.h'

      character(len=80) :: fn1,fn2,fn3
      character(*) :: ic_dir
      integer :: fstat,ierr
    
      u=0.0
      b=0.0
      if (comm_cart_rank==0) print *,'reading in mhd initial conditions'
      fn1=ic_dir(1:len_trim(ic_dir))//'mhd_ic'
      write(fn2,'(i10,".dat")') comm_cart_rank
      fn2=adjustl(fn2)
      fn3=fn1(1:len_trim(fn1))//fn2(1:len_trim(fn2))
#ifdef BINARY
      open(186,file=fn3,form='binary',status='old',iostat=fstat)
#else
      open(186,file=fn3,form='unformatted',status='old',iostat=fstat)
#endif
      if (fstat /= 0) then
        print *,'error opening mhd initial conditions'
        print *,'rank',comm_cart_rank,'file:',fn3
        call mpi_abort(mpi_comm_world,ierr,ierr) 
      endif
      read(186) u(:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      read(186) b(:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      close(186)
     
      end subroutine mpi_tvd_mhd_ic 

!=======================================================================
      subroutine mpi_tvd_mhd_restart(ic_dir,z_s)
      implicit none

      include 'mpif.h'

      character(len=80) :: fn1,fn2,fn3
      character(*) :: ic_dir
      integer :: fstat,ierr
!      real(4) :: z_write
      character(len=7) :: z_s
      character(len=4) :: rank_s
      integer :: cur_iter
      real :: cur_t
  
      u=0.0
      b=0.0
      if (comm_cart_rank==0) print *,'reading mhd restart checkpoint'
!      if (comm_cart_rank==0) z_write = z_checkpoint(restart_checkpoint)
      !call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

!      write(z_s,'(f7.3)') z_write
!      z_s=adjustl(z_s)

      write(rank_s,'(i4)') comm_cart_rank
      rank_s=adjustl(rank_s)

      fn1=ic_dir(1:len_trim(ic_dir))//z_s(1:len_trim(z_s))//'mhd'
      write(fn2,'(i10,".dat")') comm_cart_rank
      fn2=adjustl(fn2)
      fn3=fn1(1:len_trim(fn1))//fn2(1:len_trim(fn2))
#ifdef BINARY
      open(186,file=fn3,form='binary',status='old',iostat=fstat)
#else
      open(186,file=fn3,form='unformatted',status='old',iostat=fstat)
#endif
      if (fstat /= 0) then
        print *,'error opening mhd restart checkpoint'
        print *,'rank',comm_cart_rank,'file:',fn3
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif
      read(186) cur_iter,cur_t,nx,ny,nz
      read(186) u(:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      read(186) b(:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      close(186)

      end subroutine mpi_tvd_mhd_restart 

!=======================================================================
      subroutine mpi_tvd_mhd_init(nc_in,comm_cart_in, &
                    comm_cart_rank_in,comm_cart_coord_in, &
                    comm_cart_neighbor_in,comm_cart_nodes_dim)
      implicit none

      include 'mpif.h'

      integer, intent(in), dimension(3) :: nc_in 
      integer, intent(in) :: comm_cart_in 
      integer, intent(in) :: comm_cart_rank_in 
      integer, intent(in), dimension(3) :: comm_cart_coord_in 
      integer, intent(in), dimension(6) :: comm_cart_neighbor_in
      integer, intent(in) :: comm_cart_nodes_dim

      integer :: nmax,ierr

      comm_cart=comm_cart_in
      comm_cart_rank=comm_cart_rank_in
      comm_cart_coord=comm_cart_coord_in
      comm_cart_neighbor=comm_cart_neighbor_in

      nmax=maxval(nc_in) 

!.....populate the comm_wld structures..................................
      nx%g=comm_cart_nodes_dim*nc_in(1)
      nx%r=comm_cart_coord(3)*nc_in(1)-yzbuf
      nx%m=yzbuf+1
      nx%n=yzbuf+nc_in(1)
      nx%l=2*yzbuf+nc_in(1)

      ny%g=comm_cart_nodes_dim*nc_in(2)
      ny%r=comm_cart_coord(2)*nc_in(2)-yzbuf
      ny%m=yzbuf+1
      ny%n=yzbuf+nc_in(2)
      ny%l=2*yzbuf+nc_in(2)

      nz%g=comm_cart_nodes_dim*nc_in(3)
      nz%r=comm_cart_coord(1)*nc_in(3)-yzbuf
      nz%m=yzbuf+1
      nz%n=yzbuf+nc_in(3)
      nz%l=2*yzbuf+nc_in(3)
 
!.....allocate dynamic memory structures................................
      allocate(u(nu,nx%l,ny%l,nz%l))
      allocate(b(nb,nx%l,ny%l,nz%l))
      allocate(send0((nu+nb)*(xbuf+yzbuf)*(2*yzbuf+nmax)*(2*yzbuf+nmax)))
      allocate(send1((nu+nb)*(xbuf+yzbuf)*(2*yzbuf+nmax)*(2*yzbuf+nmax)))
      allocate(recv0((nu+nb)*(xbuf+yzbuf)*(2*yzbuf+nmax)*(2*yzbuf+nmax)))
      allocate(recv1((nu+nb)*(xbuf+yzbuf)*(2*yzbuf+nmax)*(2*yzbuf+nmax)))

!.....install persistent communications in x-direction..................
      call mpi_send_init(send0,(nu+nb)*(xbuf+yzbuf)*ny%l*nz%l, &
       MPI_REAL,comm_cart_neighbor(5),xtagdw,comm_cart, &
       nx%requests(1),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_send_init(send1,(nu+nb)*(xbuf+yzbuf)*ny%l*nz%l, &
       MPI_REAL,comm_cart_neighbor(6),xtagup,comm_cart, &
       nx%requests(2),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_recv_init(recv1,(nu+nb)*(xbuf+yzbuf)*ny%l*nz%l, &
       MPI_REAL,comm_cart_neighbor(6),xtagdw,comm_cart, &
       nx%requests(3),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_recv_init(recv0,(nu+nb)*(xbuf+yzbuf)*ny%l*nz%l, &
       MPI_REAL,comm_cart_neighbor(5),xtagup,comm_cart, &
       nx%requests(4),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)

!.....install persistent communications in y-direction..................
      call mpi_send_init(send0,(nu+nb)*(xbuf+yzbuf)*nz%l*nx%l, &
       MPI_REAL,comm_cart_neighbor(3),ytagdw,comm_cart, &
       ny%requests(1),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_send_init(send1,(nu+nb)*(xbuf+yzbuf)*nz%l*nx%l, &
       MPI_REAL,comm_cart_neighbor(4),ytagup,comm_cart, &
       ny%requests(2),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_recv_init(recv1,(nu+nb)*(xbuf+yzbuf)*nz%l*nx%l, &
       MPI_REAL,comm_cart_neighbor(4),ytagdw,comm_cart, &
       ny%requests(3),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_recv_init(recv0,(nu+nb)*(xbuf+yzbuf)*nz%l*nx%l, &
       MPI_REAL,comm_cart_neighbor(3),ytagup,comm_cart,  &
       ny%requests(4),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)

!.....install persistent communications in z-direction..................
      call mpi_send_init(send0,(nu+nb)*(xbuf+yzbuf)*nx%l*ny%l, &
       MPI_REAL,comm_cart_neighbor(1),ztagdw,comm_cart, &
       nz%requests(1),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_send_init(send1,(nu+nb)*(xbuf+yzbuf)*nx%l*ny%l, &
       MPI_REAL,comm_cart_neighbor(2),ztagup,comm_cart, &
       nz%requests(2),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_recv_init(recv1,(nu+nb)*(xbuf+yzbuf)*nx%l*ny%l, &
       MPI_REAL,comm_cart_neighbor(2),ztagdw,comm_cart, &
       nz%requests(3),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)
      call mpi_recv_init(recv0,(nu+nb)*(xbuf+yzbuf)*nx%l*ny%l, &
       MPI_REAL,comm_cart_neighbor(1),ztagup,comm_cart, &
       nz%requests(4),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)

      end subroutine mpi_tvd_mhd_init
!=======================================================================
      subroutine mpi_tvd_mhd_finalize
      implicit none

      deallocate(u)
      deallocate(b)
      deallocate(send0)
      deallocate(send1)
      deallocate(recv0)
      deallocate(recv1)

      end subroutine mpi_tvd_mhd_finalize
!=======================================================================
      subroutine comm_bufferupdate(u,b,nx,ny,nz)
!     this routine updates all buffers. It is intended for
!     initialization and debugging. Buffer updates during regular
!     execution are handled in the subroutine sweep.

      type(comm_wld), intent(in) :: nx,ny,nz
      real, dimension(nu,nx%l,ny%l,nz%l), intent(out) :: u
      real, dimension(nb,nx%l,ny%l,nz%l), intent(out) :: b

      integer, parameter :: b1=yzbuf
      integer :: sru,srb
      integer, dimension(4) :: ushape,bshape

!.....copy to send buffer...............................................
      sru = nu*nx%l*ny%l*b1
      srb = nb*nx%l*ny%l*b1
      send0(1:sru)         = reshape(u(:,:,:,nz%m:nz%m+b1-1),(/ sru /))
      send0(sru+1:sru+srb) = reshape(b(:,:,:,nz%m:nz%m+b1-1),(/ srb /))
      send1(1:sru)         = reshape(u(:,:,:,nz%n-b1+1:nz%n),(/ sru /))
      send1(sru+1:sru+srb) = reshape(b(:,:,:,nz%n-b1+1:nz%n),(/ srb /))
!.....do communication..................................................
      call comm_mpistartall(nz%requests)
      call comm_mpiwaitall(nz%requests)
!.....fetch from receive buffer.........................................
      ushape = (/ nu,ny%l,nz%l,b1 /)
      bshape = (/ nb,ny%l,nz%l,b1 /)
      u(:,:,:,nz%m-b1:nz%m-1) = reshape(recv0(1:sru),ushape)
      b(:,:,:,nz%m-b1:nz%m-1) = reshape(recv0(sru+1:sru+srb),bshape)
      u(:,:,:,nz%n+1:nz%n+b1) = reshape(recv1(1:sru),ushape)
      b(:,:,:,nz%n+1:nz%n+b1) = reshape(recv1(sru+1:sru+srb),bshape)

      end subroutine comm_bufferupdate
!=======================================================================
      subroutine comm_mpistartall(requests)

      implicit none

      include 'mpif.h'

      integer, dimension(nr), intent(in) :: requests
      integer :: ierr

      call mpi_startall(nr,requests,ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)

      end subroutine comm_mpistartall
!=======================================================================
      subroutine comm_mpiwaitall(requests)

      implicit none

      include 'mpif.h'

      integer, dimension(nr), intent(in) :: requests
      integer, dimension(MPI_STATUS_SIZE,nr) :: status
      integer :: ierr

      call mpi_waitall(nr,requests,status,ierr)
      if (ierr.ne.mpi_success) call comm_mpistop(ierr)

      end subroutine comm_mpiwaitall
!=======================================================================
      subroutine comm_mpistop(istop)

      implicit none

      include 'mpif.h'

      integer :: istop,ierr

      write(*,*) comm_cart_rank,': ERROR, MPI ERROR= ',istop
      write(*,*) 'please see MPI documentation for details'

      call mpi_finalize(ierr)
      stop

      end subroutine comm_mpistop
!=======================================================================
      subroutine sweep(forward,u,b,nx,ny,nz,dt)

      logical, intent(in) :: forward
      type(comm_wld), intent(in) :: nx,ny,nz
      real, intent(in) :: dt
      real, dimension(nu,nx%l,ny%l,nz%l) :: u
      real, dimension(nb,nx%l,ny%l,nz%l) :: b

!-----------------------------------------------------------------------
!
!     actions at left boundary (right boundary is symmetric):
!
!                                 ___ = array section with original data
!                        nx%m     --- = array section with updated data
!                         |       ... = array section with corrupt data
!                         V
!     |   xbuf   | yzbuf |   xbuf   |   xbuf   |
!                                   | yzbuf |
!
!                |.......|____________________________ u _______________
!                        |____ send buf ____|
!     initiate communications
!     |..................|__ ut _______________|
!                        |..........|----------------- u ---------------
!     continue when communications complete
!     |__ receive buf ___|
!     |_____________________ ut _______________|
!     |..........|---------- ut ----|..........|
!                |------------------------------------ u ---------------
!
!-----------------------------------------------------------------------
!.....nx1,nxm,nxn,nxl match ut,bt to locations 1,nx%m,nx%n,nx%l in u,b..

      integer, parameter :: b1=xbuf+yzbuf
      integer, parameter :: b2=2*xbuf
      integer, parameter :: nx1=xbuf+1,nxm=xbuf+yzbuf+1
      integer, parameter :: nxn=5*xbuf+yzbuf,nxl=5*xbuf+2*yzbuf
      integer, parameter :: nt=6*xbuf+2*yzbuf
      
      integer :: sru,srb
      integer, dimension(4) :: ushape,bshape
      real, dimension(nu,nt,ny%l,nz%l) :: ut
      real, dimension(nb,nt,ny%l,nz%l) :: bt
      type(comm_wld) :: nbd

!.....write boundaries to send buffer...................................
      sru = nu*b1*ny%l*nz%l
      srb = nb*b1*ny%l*nz%l
      send0(1:sru)         = reshape(u(:,nx%m:nx%m+b1-1,:,:),(/ sru /))
      send0(sru+1:sru+srb) = reshape(b(:,nx%m:nx%m+b1-1,:,:),(/ srb /))
      send1(1:sru)         = reshape(u(:,nx%n-b1+1:nx%n,:,:),(/ sru /))
      send1(sru+1:sru+srb) = reshape(b(:,nx%n-b1+1:nx%n,:,:),(/ srb /))
!.....initiate communication............................................
      call comm_mpistartall(nx%requests)
!.....copy boundaries of u,b into work arrays ut,bt.....................
      ut(:,nxm:nxm+b2-1,:,:) = u(:,nx%m:nx%m+b2-1,:,:)
      bt(:,nxm:nxm+b2-1,:,:) = b(:,nx%m:nx%m+b2-1,:,:)
      ut(:,nxn-b2+1:nxn,:,:) = u(:,nx%n-b2+1:nx%n,:,:)
      bt(:,nxn-b2+1:nxn,:,:) = b(:,nx%n-b2+1:nx%n,:,:)
!.....work on central parts in x-direction..............................
      if (forward) then
        call fluidx(u,b,nx%l,ny%l,nz%l,dt)
        call advectbyzx(u,b,nx%l,ny%l,nz%l,dt)
      else
        call advectbyzx(u,b,nx%l,ny%l,nz%l,dt)
        call fluidx(u,b,nx%l,ny%l,nz%l,dt)
      endif
!.....wait for communication to complete................................
      call comm_mpiwaitall(nx%requests)
!.....fetch boundaries of ut,bt from receive buffer.....................
      ushape = (/ nu,b1,ny%l,nz%l /)
      bshape = (/ nb,b1,ny%l,nz%l /)
      ut(:,nxm-b1:nxm-1,:,:) = reshape(recv0(1:sru),ushape)
      bt(:,nxm-b1:nxm-1,:,:) = reshape(recv0(sru+1:sru+srb),bshape)
      ut(:,nxn+1:nxn+b1,:,:) = reshape(recv1(1:sru),ushape)
      bt(:,nxn+1:nxn+b1,:,:) = reshape(recv1(sru+1:sru+srb),bshape)
!.....work on boundaries in x-direction.................................
      if (forward) then
        call fluidx(ut,bt,nt,ny%l,nz%l,dt)
        call advectbyzx(ut,bt,nt,ny%l,nz%l,dt)
      else
        call advectbyzx(ut,bt,nt,ny%l,nz%l,dt)
        call fluidx(ut,bt,nt,ny%l,nz%l,dt)
      endif
!.....copy result on boundaries of u and b..............................
      u(:,1:nx%m+xbuf-1,:,:) = ut(:,nx1:nxm+xbuf-1,:,:)
      b(:,1:nx%m+xbuf-1,:,:) = bt(:,nx1:nxm+xbuf-1,:,:)
      u(:,nx%n-xbuf+1:nx%l,:,:) = ut(:,nxn-xbuf+1:nxl,:,:)
      b(:,nx%n-xbuf+1:nx%l,:,:) = bt(:,nxn-xbuf+1:nxl,:,:)

      end subroutine sweep
!=======================================================================
      subroutine advectbyzx(u,b,nx,ny,nz,dt)

      implicit none

      integer, intent(in) :: nx,ny,nz
      real, intent(in) :: dt
      real, dimension(nu,nx,ny,nz) :: u
      real, dimension(nb,nx,ny,nz) :: b

      integer j,k,jm,km,im,ip
      real, dimension(nx) :: fluxbx,b1x,vx

      !$omp parallel do private(j,jm,vx,b1x,fluxbx)
      do k=1,nz
        do j=2,ny
          jm=j-1
          vx=(u(2,:,jm,k)+u(2,:,j,k))/(u(1,:,jm,k)+u(1,:,j,k))
          b1x=b(2,:,j,k)
          call tvdb(fluxbx,b1x,vx,nx,dt)
          b(2,:,j,k)=b1x
          b1x = cshift(fluxbx,-1)
          b(1,:,j,k)=b(1,:,j,k)-b1x
          b(1,:,jm,k)=b(1,:,jm,k)+b1x
        enddo
      enddo
      do k=2,nz
        !$omp parallel do private(km,vx,b1x,fluxbx)
        do j=1,ny
          km=k-1
          vx=(u(2,:,j,km)+u(2,:,j,k))/(u(1,:,j,km)+u(1,:,j,k))
          b1x=b(3,:,j,k)
          call tvdb(fluxbx,b1x,vx,nx,dt)
          b(3,:,j,k)=b1x
          b1x = cshift(fluxbx,-1)
          b(1,:,j,k)=b(1,:,j,k)-b1x
          b(1,:,j,km)=b(1,:,j,km)+b1x
        enddo
      enddo

      end subroutine advectbyzx
!=======================================================================
      subroutine calcfl(u,b,nx,ny,nz,cfl)

      implicit none

      include 'mpif.h'
      type(comm_wld), intent(in) :: nx,ny,nz
      real, intent(in) :: u(nu,nx%l,ny%l,nz%l)
      real, intent(in) :: b(nb,nx%l,ny%l,nz%l)
      real, intent(out) :: cfl

#ifdef DEBUG_CFL
      integer :: cfl_max_loc(3)
      real :: cfl_max_u(5)
      real :: cloc
#endif

      real cmax,ierr
      real v,ps,p,c,b11,b21,b31,temp5
      integer i,j,k,ip,jp,kp

      if ((nx%l.le.nx%n).or.(ny%l.le.ny%n).or.(nz%l.le.nz%n)) then
        write(6,*) 'Error: nxyz%l<=nxyz%n in calcfl.f90!'
        call comm_mpistop(0)
      endif
      c=0
      !$omp parallel do private(j,i,kp,jp,ip,b11,b21,b31,v,ps,p,temp5) reduction(max:c)
      do k=nz%m,nz%n
        kp=k+1
        do j=ny%m,ny%n
          jp=j+1
          do i=nx%m,nx%n
            ip=i+1
            b11=(b(1,i,j,k)+b(1,ip,j,k))/2
            b21=(b(2,i,j,k)+b(2,i,jp,k))/2
            b31=(b(3,i,j,k)+b(3,i,j,kp))/2
            temp5=b11*b11+b21*b21+b31*b31
            v=maxval(abs(u(2:4,i,j,k)/u(1,i,j,k)))
            ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2)   &
                 *(gamma-1)+(2-gamma)*temp5/2
            p=ps-temp5/2
#ifdef DEBUG_CFL
            if (v+sqrt(abs((temp5*2+gamma*p)/u(1,i,j,k))) > c) then
	      cfl_max_loc=(/i,j,k/)
	      cfl_max_u=u(:,i,j,k)
	      c=v+sqrt(abs((temp5*2+gamma*p)/u(1,i,j,k)))
	    endif
#else
         !  c=max(c,v+sqrt(abs((temp5*(gamma-1)+gamma*p)/u(1,i,j,k))))
            c=max(c,v+sqrt(abs((temp5*2+gamma*p)/u(1,i,j,k))))
#endif
          end do
        end do
      end do
#ifdef DEBUG_CFL
      cloc=c
#endif
      !call mpi_reduce(c,cmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,   &
      !                comm_cart, ierr)
      !call mpi_bcast(cmax,1,MPI_DOUBLE_PRECISION,0,comm_cart, ierr)
      call mpi_reduce(c,cmax,1,MPI_REAL,MPI_MAX,0,   &
                      comm_cart, ierr)
      call mpi_bcast(cmax,1,MPI_REAL,0,comm_cart, ierr)
#ifdef DEBUG_CFL
      if (cloc == cmax) print *,'cfl MAX loc',comm_cart_rank,cloc,cfl_max_loc
      if (cloc == cmax) print *,'cfl MAX u',cfl_max_u
#endif
      cfl=1/cmax

      end subroutine calcfl
!=======================================================================
      subroutine fluid_minmax(u,amax,cmax,D1,D2)
        implicit none
   
        include 'mpif.h'
 
        real, parameter :: Emin=1e-10
        real, parameter :: gg=gamma*(gamma-1)

        real :: u(nu,nx%l,ny%l,nz%l)
        real, intent(out) :: amax,cmax,D1,D2

        integer i,j,k,ip,ierr
        real E,KE,c,amaxl,cmaxl,D1l,D2l
        real dv(3,2)

        D1=1e10
        D2=-1e10
        amax=-1e10
        cmax=-1e10

        do k=nz%m,nz%n
          do j=ny%m,ny%n 
!this needs to be corrected (?)
            dv(:,2)=u(2:4,nx%m-1,j,k)/u(1,nx%m-1,j,k)
            do i=nx%m,nx%n,2
              D1=min(D1,u(1,i,j,k))
              D2=max(D2,u(1,i,j,k))
              dv(:,1)=u(2:4,i,j,k)/u(1,i,j,k)
              amax=max(amax,maxval(abs(dv(:,2)-dv(:,1))))
              KE=u(1,i,j,k)*sum(dv(:,1)**2)/2
              E=max(Emin,u(5,i,j,k)-KE)
              u(5,i,j,k)=KE+E
              c=maxval(abs(dv(:,1)))+sqrt(gg*E/u(1,i,j,k))
              cmax=max(cmax,c)

              ip=i+1
              D1=min(D1,u(1,ip,j,k))
              D2=max(D2,u(1,ip,j,k))
              dv(:,2)=u(2:4,ip,j,k)/u(1,ip,j,k)
              amax=max(amax,maxval(abs(dv(:,2)-dv(:,1))))
              KE=u(1,ip,j,k)*sum(dv(:,2)**2)/2
              E=max(Emin,u(5,ip,j,k)-KE)
              u(5,ip,j,k)=KE+E
              c=maxval(abs(dv(:,2)))+sqrt(gg*E/u(1,ip,j,k))
              cmax=max(cmax,c)
            enddo
          enddo
        enddo

        amaxl=amax
        cmaxl=cmax
        D1l=D1
        D2l=D2  
        call mpi_reduce(amaxl,amax,1,MPI_REAL,MPI_MAX,0,comm_cart, ierr)
        call mpi_reduce(cmaxl,cmax,1,MPI_REAL,MPI_MAX,0,comm_cart, ierr)
        call mpi_reduce(D1l,D1,1,MPI_REAL,MPI_MIN,0,comm_cart, ierr)
        call mpi_reduce(D2l,D2,1,MPI_REAL,MPI_MAX,0,comm_cart, ierr)
        call mpi_bcast(amax,1,MPI_REAL,0,comm_cart, ierr)
        call mpi_bcast(cmax,1,MPI_REAL,0,comm_cart, ierr)
        call mpi_bcast(D1,1,MPI_REAL,0,comm_cart, ierr)
        call mpi_bcast(D2,1,MPI_REAL,0,comm_cart, ierr)
 
      end subroutine fluid_minmax
!=======================================================================
      subroutine fluidx(u,b,nx,ny,nz,dt)

      implicit none

      integer, intent(in) :: nx,ny,nz
      real, intent(in) :: dt
      real, dimension(nu,nx,ny,nz) :: u
      real, dimension(nb,nx,ny,nz) :: b

      real, dimension(nb,nx) :: b3x
      real, dimension(nu,nx) :: u1x
      integer j,k,jp,kp

      !$omp parallel do private(j,u1x,b3x,jp,kp)
      do k=1,nz-1
        kp = k+1
        do j=1,ny-1
          jp = j+1
          b3x=0.5*b(:,:,j,k)
          b3x(1,:)=b3x(1,:)+cshift(b3x(1,:),1)
          b3x(2,:)=b3x(2,:)+0.5*b(2,:,jp,k)
          b3x(3,:)=b3x(3,:)+0.5*b(3,:,j,kp)
          call tvd1(u(:,:,j,k),b3x,nx,dt)
        end do
      end do

      end subroutine fluidx
!=======================================================================
      subroutine mpi_tvd_mhd_init_cond(u,b,nx,ny,nz)
    
      implicit none
  
      type(comm_wld), intent(in) :: nx,ny,nz
      real, intent(out) :: u(nu,nx%l,ny%l,nz%l)
      real, intent(out) :: b(nb,nx%l,ny%l,nz%l)

      real, parameter :: p0=3.0/5.0
      real, parameter :: epsilon=0.1

      integer i,im
  
      u=0
      b=0
      u(1,:,:,:)=1
      u(5,:,:,:)=1.5*p0
  
      u(2,:,1,1)=0.0001*sin( (/ (2*3.14159*(nx%r+i)/nx%g, i=1,nx%l) /) )
      u(1,:,1,1)=u(1,:,1,1)+u(2,:,1,1)
      u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+u(2,:,:,:)**2/u(1,:,:,:)/2
      u(5,:,1,1)=u(5,:,1,1)+p0*u(2,:,1,1)*5./3./(2./3.)
!    return

! circularly polarized alven wave:
! background : rho=B_x=1 all others zero
      u=0
      b=0
      b(1,:,:,:)=1
      u(1,:,:,:)=1
      u(5,:,:,:)=0.001  ! to keep things stable
      do i=1,nx%l
        u(3,i,:,:)=epsilon*sin( 2*3.14159*(nx%r+i)/nx%g )
        u(4,i,:,:)=epsilon*cos( 2*3.14159*(nx%r+i)/nx%g )
      enddo
      b(2,:,:,:)=-u(3,:,:,:)
      b(3,:,:,:)=-u(4,:,:,:)

      if (ny%m.lt.2) write(6,*) 'Warning: ny%m<2 in init.f90'
      do i=2,ny%l
        im=i-1
        b(2,:,i,:)=(b(2,:,i,:)+b(2,:,im,:))/2
      enddo
      if (nz%m.lt.2) write(6,*) 'Warning: nz%m<2 in init.f90'
      do i=2,nz%l
        im=i-1
        b(3,:,:,i)=(b(3,:,:,i)+b(3,:,:,im))/2
      enddo
      u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+sum(u(2:4,:,:,:)**2,1) &
                 /u(1,:,:,:)/2
      return

! alven wave:
! background : rho=B_x=1 all others zero
! \dot{v_y}+\div_x(-B_y)=0
! \dot{B_y}+\div_x (       -   v_y)     = 0
! let v_y=\epsilon sin(2 \pi (x-t)/L)
! then B_y=-v_y
      u=0
      b=0
      b(1,:,:,:)=1
      u(1,:,:,:)=1
      u(5,:,:,:)=.001  ! to keep things stable
      do i=1,nx%g
        u(3,i,:,:)=0.1*sin( 2*3.14159*i/nx%g )
      enddo
      b(2,:,:,:)=-u(3,:,:,:)
      b(2,:,:,:)=(b(2,:,:,:)+cshift(b(2,:,:,:),-1,1))/2
      u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+u(3,:,:,:)**2/u(1,:,:,:)/2
      return

! magnetosonic
      u=0
      b=0
      b(2,:,:,:)=1
      u(1,:,:,:)=1
      u(5,:,:,:)=.001  ! to keep things stable
      do i=1,nx%g
        u(2,i,:,:)=0.01*sin( 2*3.14159*i/nx%g )
      enddo
      b(2,:,:,:)=b(2,:,:,:)+u(2,:,:,:)
      u(1,:,:,:)=b(2,:,:,:)
      u(5,:,:,:)=u(5,:,:,:)+sum(b**2,1)/2+sum(u(2:4,:,:,:)**2,1) &
                 /u(1,:,:,:)/2

      end subroutine mpi_tvd_mhd_init_cond
!=======================================================================
      recursive subroutine mhdflux(fr,flm,u,b)

      implicit none
      real, dimension(nu), intent(in) :: u
      real, dimension(nb), intent(in) :: b
      real, dimension(nu), intent(out) :: fr,flm

      real, dimension(nu) :: v
      real :: vx,usqr,bsqr,bdotu,ps,p,c

      vx = u(2)/u(1)
      usqr  = u(2)*u(2) + u(3)*u(3) + u(4)*u(4)
      bsqr  = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)
      bdotu = b(1)*u(2) + b(2)*u(3) + b(3)*u(4)
      ps = (u(5) - 0.5*usqr/u(1))*(gamma-1) + (2-gamma)*0.5*bsqr
      v(1) = u(2)
      v(2:4) = u(2:4)*vx - b*b(1)
      v(2) = v(2) + ps
      v(5) = (u(5)+ps)*vx - b(1)*bdotu/u(1)
      p = ps - 0.5*bsqr
      c = abs(vx) + sqrt(abs( (bsqr + gamma*p)/u(1) ))
      if (c.gt.0.) v = v/c
      fr = c*(u+v)
      flm = c*(u-v)

      end subroutine mhdflux
!=======================================================================
      subroutine transposef(u,b,nx,ny,nz)

      implicit none

      integer nx,ny,nz
      real u(nu,nx,ny,nz),b(nb,nx,ny,nz)

      integer i,j,k
      real ut(nu,ny,nz,nx),bt(nb,ny,nz,nx)

      !$omp parallel do default(none) shared(u,b,ut,bt,ny,nx,nz) private(i,j)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            ut(:,j,k,i)=u((/1,3,4,2,5/),i,j,k)
            bt(:,j,k,i)=b((/2,3,1/),i,j,k)
          end do
        end do
      end do
      call copyf(u,ut,nx*ny*nz*5)
      call copyf(b,bt,nx*ny*nz*3)

      end subroutine transposef
!=======================================================================
      subroutine transposeb(u,b,nx,ny,nz)
     
      implicit none

      integer nx,ny,nz
      real u(nu,nx,ny,nz),b(nb,nx,ny,nz)

      integer i,j,k
      real ut(nu,nz,nx,ny),bt(nb,nz,nx,ny)

      !$omp parallel do default(none) shared(u,b,ut,bt,ny,nx,nz) private(i,j)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            ut(:,k,i,j)=u((/1,4,2,3,5/),i,j,k)
            bt(:,k,i,j)=b((/3,1,2/),i,j,k)
          end do
        end do
      end do
      call copyf(u,ut,nx*ny*nz*5)
      call copyf(b,bt,nx*ny*nz*3)

      end subroutine transposeb
!=======================================================================
      subroutine copyf(a,b,n)

      implicit none

      integer n
      real a(n),b(n)

      a=b

      end subroutine copyf
!=======================================================================
      recursive subroutine tvd1(u,b,n,dt)

      implicit none

      integer, intent(in) :: n
      real, intent(in) :: dt
      real, dimension(nu,n) :: u
      real, dimension(nb,n) :: b

      real, dimension(nu) :: u4,u3,u2,u1,uu
      real, dimension(nb) :: b4,b3
      real, dimension(nu) :: gr4,gr3,gl3
      real, dimension(nu) :: fr3,fr2,fr1,fl2,fl1
      real, dimension(nu,2) :: df1,df0,d1,tmp !1 is left and 2 is right
      real, dimension(nu) :: flux3,flux2,flux1,flux0

      integer :: i

!.....initialize........................................................
      df1 = 0.
      flux1 = 0.
      u2 = 0.
      fr2 = 0.
      fl2 = 0.
      u3 = 0.
      fr3 = 0.
      flux3 = 0.
      u4 = 0.
      b4 = 0.
      gr4 = 0.
!.....down shift variables belonging to second update...................
      do i=1,n
        df0 = df1
        flux0 = flux1
        u1 = u2
        fr1 = fr2
        fl1 = fl2
        u2 = u3
        fr2 = fr3
!.....down shift variables belonging to first update....................
        flux2 = flux3
        u3 = u4
        b3 = b4
        gr3 = gr4
        u4 = u(:,i)
        b4 = b(:,i)
!.....first update......................................................
        call mhdflux(gr4,gl3,u4,b4)
        flux3 = 0.5*(gr3-gl3)
        uu = u3 - 0.5*(flux3-flux2)*dt
!.....second update.....................................................
        call mhdflux(fr3,fl2,uu,b3)
        df1(:,1) = 0.5*(fl1-fl2)
        df1(:,2) = 0.5*(fr2-fr1)
        tmp = df1*df0
        where (tmp.gt.0.)
          d1 = 2.*tmp/(df1+df0)
        elsewhere
          d1 = 0.
        end where
        flux1 = 0.5*(fr1-fl1+d1(:,2)-d1(:,1))
        if (i.gt.3) u(:,i-3) = u1 - (flux1-flux0)*dt
      enddo

      end subroutine tvd1
!=======================================================================
      recursive subroutine tvdb(flux,b,vg,n,dt)

      implicit none

      integer, intent(in) :: n
      real, intent(in) :: dt
      real, dimension(n) :: flux,b,vg

!     unlike the B field, the flux lives on the right cell boundary

      integer :: i,i4,i1
      real :: b4,b3,b2,b1,bb3,tmp
      real :: f3,f2,f1,f0,w3,w2,w1,dw2,dw1,d2,d1
      real :: vg4,vg3,vh3,vh2,vh1

!.....initialize........................................................
      f1 = 0.
      w2 = 0.
      dw2 = 0.
      d2 = 0.
      vh2 = 0.
      b2 = 0.
      vh3 = 0.
      f3 = 0.
      w3 = 0.
      b3 = 0.
      vg4 = 0.
      b4 = 0.
!.....down shift variables belonging to second update...................
      do i=1,n
        f0 = f1
        vh1 = vh2
        w1 = w2
        dw1 = dw2
        d1 = d2
        b1 = b2
        vh2 = vh3
        w2 = w3
        b2 = b3
!.....down shift variables belonging to first update....................
        f2 = f3
        vg3 = vg4
        b3 = b4
        vg4 = vg(i)
        b4 = b(i)
!.....first update......................................................
        vh3 = 0.5*(vg3+vg4)
        if (vh3.gt.0.) then
          f3 = b3*vg3
        else
          f3 = b4*vg4
        endif
        bb3 = b3 - 0.5*(f3-f2)*dt
!.....second update.....................................................
        w3 = vg3*bb3
        dw2 = 0.5*(w3-w2)
        tmp = dw1*dw2
        if (tmp.gt.0.) then
          d2 = 2.*tmp/(dw1+dw2)
        else
          d2 = 0.
        endif
        if (vh1.gt.0) then
          f1=(w1+d1)*dt
        else
          f1=(w2-d2)*dt
        end if
        if (i.gt.3) then
          flux(i-3) = f1
          b(i-3) = b1 - (f1-f0)
        endif
      end do

      end subroutine tvdb
!=======================================================================
      end module mpi_tvd_mhd
!=======================================================================
