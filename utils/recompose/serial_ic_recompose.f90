!! This is the serial cubepm initial condition recomposer.
!! Accumulates distributed particle list into monolithic file 
!! compile with: ifort serial_ic_recompose.f90 -o s_recompose

implicit none

! frequently changed parameters are found in this header file:
include '/home/merz/codes/cubepm/parameters'

character(len=*),parameter :: serial_output_path='./'
real, parameter    :: rnc = nc
!! Set this to change the mesh size of the recomposed output
real, parameter    :: conv_scale = 0.25
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
integer(4) :: nn_returned, ierr,rank,fstat

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

!! Read in particle positions

  np_found=0
  np_current=0  

  do rank=0,nn-1
    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)
    ofile=serial_output_path//'xv'//rank_s(1:len_trim(rank_s))//'.ic'

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

    read(21) np_local

    if (np_local > max_np) then
      write(*,*) 'too many particles to store'
      write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    np_found(rank)=np_local
    close(21)
  enddo

    write(*,*) 'particles / node'
    do i=0,nn-1
      write(*,*) i, np_found(i)
    enddo
    write(*,*) 'total number of particles=',sum(np_found)
    ofile=serial_output_path//'xv.ic'
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
      ofile=serial_output_path//'xv'//rank_s(1:len_trim(rank_s))//'.ic'

#ifdef BINARY
      open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
#else
      open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

      if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        stop 
      endif

      read(21) np_local

      read(21) xv(:,:np_local)
!      do j=1,np_local
!        read(21) xv(:,j)
!      enddo
      close(21)
!! prepare particles for accumulation 
      do j=1,np_local
        xv(1:3,j)=modulo(xv(1:3,j)+node_coords(:,rank)*ncc,rnc)
        xv(:,j)=conv_scale*xv(:,j)
      enddo
      write(23) xv(:,:np_local)
    enddo
    close(23)
end 
