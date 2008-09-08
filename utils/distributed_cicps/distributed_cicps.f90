!! distributed_cicps.f90 :: cubic decomposition mass power spectrum code

!! Original cicpow.f90 written by Hy Trac on Feb 23, 2003
!! Compile with mpif77 distributed_cicps.f90 -o dcicps.x 

program main
  implicit none

  include 'mpif.h'

  !! Cosmological parameters
  real, parameter :: box=100

  !! nt is the number of threads
  !! dir is the directory for data files
  integer, parameter :: nt=2

  integer, parameter      :: nc=320       !! number of cells total 
  integer, parameter 	:: hc=nc/2
  integer, parameter      :: nn_dim=2     !! number of nodes / dimension
  character(len=*),parameter :: output_path = '/scratch/merz/cubepm/160_10/'
  character(len=*),parameter :: checkpoints = '/home/merz/cubepm/inputs/checkpoints'
  real, parameter         :: rnc = nc
  integer, parameter      :: max_input=100!! maximum number of checkpoints
  integer, parameter      :: nn=nn_dim**3 !! number of nodes total
  real, parameter         :: ncc=nc/nn_dim!! number of cells / cube / node 
  integer, parameter      :: hc=nc/2      !! half-grid length
  integer, parameter      :: np=hc**3     !! maximum number of particles total
  integer, parameter      :: max_np=np    !! maximum number of particles / node

  real, dimension(6,max_np)  :: xv        !! particle list

  integer :: node_coords(3)               !! coordinates of node in cube  
  character(len=80) ofile
  character(len=4) rank_s
  character(len=7) z_s
  integer i,j,k,m
  real x,y,z,z_write
  integer(4) :: nn_returned, ierr,rank, num_checkpoints,fstat,cur_checkpoint
  real(4), dimension(max_input) :: z_checkpoint

  integer(4) :: np_local,np_current, nts, sim_checkpoint, sim_projection
  integer(4) :: nodes_returned,tag
  real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p,np_total

  integer(4), dimension(0:nn-1) :: np_found
  integer(4), dimension(mpi_status_size) :: stat

  integer :: rs

  integer, dimension(np_max) :: ll
  integer, dimension(ncc,ncc,ncc) :: hoc

  !! Power spectrum arrays
  real, dimension(2,nc) :: ps
  real, dimension(ncc,ncc,ncc) :: rho

  common rho,xv,ll,htoc,hoc,ps

!! Initialize MPI

  call mpi_init(ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
  if (nodes_returned /= nn ) then
    write(*,*) 'recompose compiled for a different number of nodes'
    write(*,*) 'mpirun nodes=',nodes_returned,'recompose nodes=',nn
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

!! Calculate node_coords

  do k=1,nn_dim
    do j=1,nn_dim
      do i=1,nn_dim
        if (rank == (i-1)+(j-1)*nn_dim+(k-1)*nn_dim**2)  &
           node_coords(:)=(/(i-1),(j-1),(k-1)/)
      enddo
    enddo
  enddo

  if (rank == 0) write(*,*) 'rank, cartesian coords in cube'

  do i=0,nn-1
    if (rank == i) write(*,*) rank,node_coords
    call mpi_barrier(mpi_comm_world,ierr)
  enddo


!! Read in checkpoints to recompose

if (rank == 0) then
  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) then
    write(*,*) 'error opening checkpoint list file'
    write(*,*) 'rank',rank,'file:',checkpoints
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
  do num_checkpoints=1,max_input
    read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
  enddo
41  num_checkpoints=num_checkpoints-1
51  close(11)
  write(*,*) 'checkpoints to recompose:'
  do i=1,num_checkpoints
    write(*,'(f5.1)') z_checkpoint(i)
  enddo
endif

call mpi_bcast(num_checkpoints,1,mpi_integer,0,mpi_comm_world,ierr)

do cur_checkpoint=1,num_checkpoints

!! Read in particle positions

  np_found=0
  np_current=0

  if (rank == 0) then
    z_write = z_checkpoint(cur_checkpoint)
    write(*,*) 'processing z=',z_write
  endif

  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)

  write(rank_s,'(i4)') rank
  rank_s=adjustl(rank_s)

  ofile=output_path//z_s(1:len_trim(z_s))//'xv'// &
        rank_s(1:len_trim(rank_s))//'.dat'

  open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')

  if (fstat /= 0) then
    write(*,*) 'error opening checkpoint'
    write(*,*) 'rank',rank,'file:',ofile
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

  read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,sim_checkpoint, &
               sim_projection,mass_p

  if (np_local > max_np) then
    write(*,*) 'too many particles to store'
    write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

!! tally up total number of particles

  call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
                       mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'number of particles =', int(np_total,4)
  if (rank == 0 .and. int(np_total,4) > max_np) then
    write(*,*) 'np total > max_np',np_total,max_np
    call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

  call mpi_gather(np_local,1,mpi_integer,np_found,1,mpi_integer,0, &
                mpi_comm_world,ierr)

  if (rank == 0) then
    write(*,*) 'particles / node'
    do i=0,nn-1
      write(*,*) i, np_found(i)
    enddo
  endif

  read(21) xv(:,:np_local)
  close(21)

!! Construct linked list

  call linked_list

!! Pass Particles

  call particle_pass


  call cic
  call powerspectrum
  call writeps

enddo

contains

  subroutine particle_pass
    implicit none

    real(4), parameter :: rnf_buf = 8
    real(4), parameter :: 
    integer(4) :: i,j,k,pp,tag
    integer(4) :: nppx,nppy,nppz,npmx,npmy,npmz,np_buf
    integer(4), dimension(mpi_status_size) :: status,sstatus,rstatus
    integer(4) :: srequest,rrequest,sierr,rierr

    tag=11
    np_buf = 0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_h - hoc_pass_depth, hoc_nc_h
          pp = hoc(i,j,k,1)
          do while (pp /= 0)
            if (xv(1,pp) >= nf_physical_node_dim - rnf_buf) then
              np_buf(1) = np_buf(1) + 1
              send_buf((np_buf(1)-1)*6+1:np_buf(1)*6)=xv(:,pp)
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo


  end subroutine particle_pass

  subroutine writeps
    implicit none
    integer k
    character*60 fn
    character*7 r_s

    real time1,time2
    call cpu_time(time1)


    !! Output power spectrum
    !! First column is k
    !! Second column is \Delta^2
   
    write(r_s,'(f7.3)') zout(rs)
    fn=r_s(1:len_trim(r_s))//'ps.dm.dat'
    open(11,file=fn)
    do k=2,hc+1
       write(11,*) ps(:,k)
    enddo
    close(11)


    call cpu_time(time2)
    time2=time2-time1
    write(*,*) time2,'  Called write ps'
    return
  end subroutine writeps


  subroutine cic
    implicit none
    integer, parameter :: kpt=nc/nt+min(1,mod(nc,nt))
    integer :: npt

    integer it,i,j,k,ip,toe

    real time1,time2
    call cpu_time(time1)

    npt=np/nt+min(1,mod(np,nt))

    !! Construct chaining lists in parallel
    !$omp parallel do default(shared) private(it,ip,j,k)
    do it=1,nt
       htoc(:,:,:,it)=0
       do ip=1+(it-1)*npt,min(np,it*npt)
          j=floor(xv(2,ip))+1
          k=floor(xv(3,ip))+1
          ll(ip)=0
          if (htoc(1,j,k,it) .eq. 0) then
             htoc(:,j,k,it)=ip
          else
             ll(htoc(2,j,k,it))=ip
             htoc(2,j,k,it)=ip
          endif
       enddo
    enddo
    !$omp end parallel do


    !! Merge chaining lists
    !$omp parallel do default(shared) private(it,j,k,toe)
    do k=1,nc
       do j=1,nc
          hoc(j,k)=0
          do it=1,nt
             if (hoc(j,k) .eq. 0) then
                hoc(j,k)=htoc(1,j,k,it)
                toe=htoc(2,j,k,it)
             else
                if (htoc(1,j,k,it) .ne. 0) then
                   ll(toe)=htoc(1,j,k,it)
                   toe=htoc(2,j,k,it)
                endif
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Initialize density field
    !$omp parallel do default(shared) private(it,k)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          d(:,:,k)=0
       enddo
    enddo
    !$omp end parallel do


    !! Add particle density to density field
    !$omp parallel do default(shared) private(it,j,k,ip)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          do j=1,nc
             ip=hoc(j,k)
             call cicmass(ip)
          enddo
       enddo
    enddo
    !$omp end parallel do


    !! Convert density to overdensity
    !$omp parallel do default(shared) private(it,i,j,k)
    do it=1,nt
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          do j=1,nc
             do i=1,nc
                d(i,j,k)=d(i,j,k)-1
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,*) time2,'  Called cic'
    return
  end subroutine cic


  subroutine cicmass(ip)
    implicit none
    real, parameter :: ncr=nc
    real :: mp

    integer ip

    integer i1,i2,j1,j2,k1,k2
    real x,y,z,dx1,dx2,dy1,dy2,dz1,dz2

    mp=ncr**3/np

    do
       if (ip .eq. 0) exit

       x=mod(xv(1,ip)-0.5+ncr,ncr)
       y=mod(xv(2,ip)-0.5+ncr,ncr)
       z=mod(xv(3,ip)-0.5+ncr,ncr)

       i1=floor(x)+1
       i2=mod(i1,nc)+1
       dx1=i1-x
       dx2=1-dx1
       j1=floor(y)+1
       j2=mod(j1,nc)+1
       dy1=j1-y
       dy2=1-dy1
       k1=floor(z)+1
       k2=mod(k1,nc)+1
       dz1=k1-z
       dz2=1-dz1

       dz1=mp*dz1
       dz2=mp*dz2
       d(i1,j1,k1)=d(i1,j1,k1)+dx1*dy1*dz1
       d(i2,j1,k1)=d(i2,j1,k1)+dx2*dy1*dz1
       d(i1,j2,k1)=d(i1,j2,k1)+dx1*dy2*dz1
       d(i2,j2,k1)=d(i2,j2,k1)+dx2*dy2*dz1
       d(i1,j1,k2)=d(i1,j1,k2)+dx1*dy1*dz2
       d(i2,j1,k2)=d(i2,j1,k2)+dx2*dy1*dz2
       d(i1,j2,k2)=d(i1,j2,k2)+dx1*dy2*dz2
       d(i2,j2,k2)=d(i2,j2,k2)+dx2*dy2*dz2

       ip=ll(ip)
    enddo

    return
  end subroutine cicmass


!!--------------------------------------------------------------!!


  subroutine powerspectrum
    implicit none
    integer, parameter :: kpt=nc/nt+min(1,mod(nc,nt))

    integer it,i,j,k
    real kr,kx,ky,kz,k1,k2,w1,w2,pow,pi
    real pst(2,nc,nt)

    real time1,time2
    call cpu_time(time1)


    !! Forward FFT density field
    call fft3d(d,nc,'f')


    !! Compute power spectrum
    
    !$omp parallel do default(shared) &
    !$omp& private(it,i,j,k,kr,kx,ky,kz,k1,k2,w1,w2,pow)
    do it=1,nt
       pst(:,:,it)=0
       do k=1+(it-1)*kpt,min(nc,it*kpt)
          if (k .lt. hc+2) then
             kz=k-1
          else
             kz=k-1-nc
          endif
          do j=1,nc
             if (j .lt. hc+2) then
                ky=j-1
             else
                ky=j-1-nc
             endif
             do i=1,nc+2,2
                kx=(i-1)/2
                kr=sqrt(kx**2+ky**2+kz**2)
                if (kr .ne. 0) then
                   k1=ceiling(kr)
                   k2=k1+1
                   w1=k1-kr
                   w2=1-w1
                   pow=(d(i,j,k)/nc**3)**2+(d(i+1,j,k)/nc**3)**2
                   pst(1,k1,it)=pst(1,k1,it)+w1
                   pst(2,k1,it)=pst(2,k1,it)+w1*pow
                   pst(1,k2,it)=pst(1,k2,it)+w2
                   pst(2,k2,it)=pst(2,k2,it)+w2*pow
                endif
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

    !! Merge power spectrum from threads
    ps=0
    do it=1,nt
       ps=ps+pst(:,:,it)
    enddo


    !! Divide by weights
    !! Convert P(k) to \Delta^2(k)
    !! Store k in ps(1,i)
    
    pi=acos(-1.)
    do k=1,nc
       if (ps(1,k) .ne. 0) then
          ps(2,k)=4*pi*((k-1)**3)*ps(2,k)/ps(1,k)
          ps(1,k)=2*pi*(k-1)/box
       endif
    enddo


    call cpu_time(time2)
    time2=(time2-time1)/nt
    write(*,*) time2,'  Called power spectrum'
    return
  end subroutine powerspectrum


!!--------------------------------------------------------------!!


  subroutine fft3d(a,n,c)
    implicit none
    integer n
    character c
    real a(n+2,n,n)
    external sfft
    
    if (c .eq. 'f') then
       call sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
    else
       call sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
    endif
    return
  end subroutine fft3d


end program main
