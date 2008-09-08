!! This is the serial cubepm particle list recomposer.
!! Accumulates distributed particle list into monolithic file 
!! compile with: ifort serial_recompose.f90 -o s_recompose

implicit none

! frequently changed parameters are found in this header file:
include '/home/merz/codes/cubepm/parameters'

! list of redshift to recompose (same as in cubepm.par)

character(len=*),parameter :: checkpoints =cubepm_root//'input/checkpoints' 
character(len=*),parameter :: serial_output_path='./'
real, parameter    :: rnc = nc
integer, parameter :: max_input=100 !! maximum number of checkpoints
integer, parameter :: nn=nodes_dim**3  !! number of nodes total
real, parameter    :: ncc=nc/nodes_dim !! number of cells / cubic 
integer, parameter :: hc=nc/2       !! half-grid length
integer, parameter :: max_np=6*(ncc/2)**3  !! maximum number of particles / node

real, dimension(6,max_np)  :: xv    !! particle list

integer :: node_coords(3,0:nn-1)           !! coordinates of node in cube	
character(len=80) ofile 
character(len=4) rank_s
character(len=7) z_s
integer i,j,k,m
real x,y,z,z_write
integer(4) :: nn_returned, ierr,rank, num_checkpoints,fstat,cur_checkpoint
real(4), dimension(max_input) :: z_checkpoint

integer(4) :: np_local,np_current, nts, sim_checkpoint, sim_projection, sim_halofind
integer(4) :: nodes_returned,tag
real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p,np_total

integer(4), dimension(0:nn-1) :: np_found
!integer(4), dimension(mpi_status_size) :: status

!! Calculate node_coords

do rank=0,nn-1
  do k=1,nodes_dim
    do j=1,nodes_dim
      do i=1,nodes_dim
        if (rank == (i-1)+(j-1)*nodes_dim+(k-1)*nodes_dim**2)  &
           node_coords(:,rank)=(/(i-1),(j-1),(k-1)/)
      enddo
    enddo
  enddo
enddo

write(*,*) 'rank, cartesian coords in cube'

do rank=0,nn-1
  write(*,*) rank,node_coords(:,rank)
enddo

!! Read in checkpoints to recompose

  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) then
    write(*,*) 'error opening checkpoint list file'
    write(*,*) 'file:',checkpoints
    stop !call mpi_abort(mpi_comm_world,ierr,ierr)
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

do cur_checkpoint=1,num_checkpoints

!! Read in particle positions

  np_found=0
  np_current=0  

  z_write = z_checkpoint(cur_checkpoint)
  write(*,*) 'processing z=',z_write

  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)

  do rank=0,nn-1
    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)
    ofile=serial_output_path//z_s(1:len_trim(z_s))//'xv'// &
    rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
    open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
    open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

    if (fstat /= 0) then
      write(*,*) 'error opening checkpoint'
      write(*,*) 'rank',rank,'file:',ofile
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

!    read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,sim_checkpoint, &
    read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,dt_pp_acc,sim_checkpoint, &
                 sim_projection,sim_halofind,mass_p

    if (np_local > max_np) then
      write(*,*) 'too many particles to store'
      write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    np_found(rank)=np_local
    close(21)
  enddo

!! tally up total number of particles

  !call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
  !                     mpi_sum,0,mpi_comm_world,ierr)
  !if (rank == 0) write(*,*) 'number of particles =', int(np_total,4)
!  if (rank == 0 .and. int(np_total,4) > max_np) then
!    write(*,*) 'np total > max_np',np_total,max_np
!    call mpi_abort(mpi_comm_world,ierr,ierr)
!  endif

  !call mpi_gather(np_local,1,mpi_integer,np_found,1,mpi_integer,0, &
  !              mpi_comm_world,ierr)

  !if (rank == 0) then
    write(*,*) 'particles / node'
    do i=0,nn-1
      write(*,*) i, np_found(i)
    enddo
  !endif

!  if (rank==0) then
!! local root process writes first
    ofile=serial_output_path//z_s(1:len_trim(z_s))//'xv.dat'
#ifdef BINARY
    open(23,file=ofile,form='binary',status='replace',iostat=fstat)
#else
    open(23,file=ofile,form='unformatted',status='replace',iostat=fstat)
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening composed checkpoint'
      write(*,*) 'rank',rank,'file:',ofile
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    np_current=sum(np_found)
    write(23) np_current

    do rank=0,nn-1
      write(rank_s,'(i4)') rank
      rank_s=adjustl(rank_s)
      ofile=serial_output_path//z_s(1:len_trim(z_s))//'xv'// &
      rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
      open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
      open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        stop !call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

!      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,sim_checkpoint, &
      read(21) np_local,a,t,tau,nts,dt_f_acc,dt_c_acc,dt_pp_acc,sim_checkpoint, &
                 sim_projection,sim_halofind,mass_p

      read(21) xv(:,:np_local)
!      do j=1,np_local
!        read(21) xv(:,j)
!      enddo
      close(21)
!! prepare particles for accumulation 
      do j=1,np_local
        xv(1:3,j)=modulo(xv(1:3,j)+node_coords(:,rank)*ncc,rnc)
      enddo
      write(23) xv(:,:np_local)
    enddo
    close(23)

enddo

end 
