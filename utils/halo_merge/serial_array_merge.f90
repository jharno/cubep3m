!! This is the serial cubep3m halo array merging program 
!! compile with: ifort -fpp -DBINARY serial_array_merge.f90 -o serial_array_merge

implicit none

! frequently changed parameters are found in this header file:
include '../../parameters'

! list of redshift to recompose (same as in cubepm.par)
!character(len=*),parameter :: checkpoints =cubepm_root//'input/halopoints' 
character(len=*),parameter :: checkpoints =cubepm_root//'input/halofinds' 
character(len=*),parameter :: serial_output_path=output_path
real, parameter    :: rnc = nc
integer, parameter :: max_input=100 !! maximum number of checkpoints
integer, parameter :: nn=nodes_dim**3  !! number of nodes total
integer(4), parameter    :: ncc=nc/nodes_dim !! number of cells / cubic 
!integer(4), parameter :: nc_node_dim = ncc
!integer(4), parameter   :: nc_dim = nc
integer(4), parameter :: nc_dim = nc/4
integer(4), parameter :: nc_node_dim = nc_dim / nodes_dim


integer, parameter :: hc=nc/2       !! half-grid length

integer :: node_coords(3,0:nn-1)           !! coordinates of node in cube	
character(len=80) ofile,fn3,fn1,fn2
character(len=4) rank_s
character(len=7) z_s
integer i,j,k,m
real x,y,z,z_write
integer(4) :: nn_returned, ierr,rank, num_checkpoints,fstat,cur_checkpoint
real(4), dimension(max_input) :: z_checkpoint

integer(4) :: nodes_returned,tag

integer(4),parameter :: ctf_scale = 4
integer(4),parameter :: fc_scale = 8 
integer(4),parameter :: rmesh = (nc_node_dim*ctf_scale)/fc_scale


real(4), dimension(:), allocatable :: rarr
real(4), dimension(:,:), allocatable :: varr

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
!    write(*,'(f5.1)') z_checkpoint(i)
    write(*,'(f5.2)') z_checkpoint(i)
  enddo

do cur_checkpoint=1,num_checkpoints

  z_write = z_checkpoint(cur_checkpoint)
  write(*,*) 'processing z=',z_write

  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)

!! first do rho_c

  do rank=0,nn-1
    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)
    ofile=serial_output_path//z_s(1:len_trim(z_s))//'rho_c'// &
    rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
    open(unit=21+rank,file=ofile,status='old',iostat=fstat,form='binary')
#else
    open(unit=21+rank,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

    if (fstat /= 0) then
      write(*,*) 'error opening coarse density'
      write(*,*) 'rank',rank,'file:',ofile
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

  enddo

! now time to do gas, lets first read in dimensions

!  print*,'nc_node_dim=',nc_node_dim
  allocate(rarr(nc_node_dim))


! now we need to read in data from each node, and write it appropriately.
! can't read in from nodes sequentially, have to jump around.

    fn1=serial_output_path(1:len_trim(serial_output_path))//z_s(1:len_trim(z_s))//'rho_c.dat'
#ifdef BINARY
    open(199,file=fn1,form='binary',status='replace',iostat=fstat)
#else 
    open(199,file=fn1,form='unformatted',status='replace',iostat=fstat)
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening global coarse mesh density'
      write(*,*) 'file:',fn1
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

! now we need to read in data from each node, and write it appropriately.
! can't read in from nodes sequentially, have to jump around.

  do k=1,nc_dim
    do j=1,nc_dim
      do i=1,nc_dim,nc_node_dim
        do rank=0,nn-1
          if ( node_coords(1,rank)==(i-1)/nc_node_dim.and. node_coords(2,rank)==(j-1)/nc_node_dim.and. &
               node_coords(3,rank)==(k-1)/nc_node_dim ) then
!             print '(7i8)',i,j,k,rank,node_coords(:,rank)
            read(21+rank) rarr
            write(199) rarr
          endif
        enddo
      enddo
    enddo
  enddo

  close(199)
  do rank=0,nn-1
    close(21+rank)
  enddo

  deallocate(rarr)

#ifdef HALO_VEL_FIELD

!! next do velocity field
!  real(4), dimension(3,nc_node_dim*mesh_scale/fine_clumping_scale,nc_node_dim*mesh_scale/fine_clumping_scale, &
!                     nc_node_dim*mesh_scale/fine_clumping_scale) :: velocity_field

  do rank=0,nn-1
    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)
    ofile=serial_output_path//z_s(1:len_trim(z_s))//'vel'// &
    rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
    open(unit=21+rank,file=ofile,status='old',iostat=fstat,form='binary')
#else
    open(unit=21+rank,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

    if (fstat /= 0) then
      write(*,*) 'error opening velocity field'
      write(*,*) 'rank',rank,'file:',ofile
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

  enddo

  allocate(varr(3,rmesh))

! now we need to read in data from each node, and write it appropriately.
! can't read in from nodes sequentially, have to jump around.

    fn1=serial_output_path(1:len_trim(serial_output_path))//z_s(1:len_trim(z_s))//'vel.dat'
#ifdef BINARY
    open(199,file=fn1,form='binary',status='replace',iostat=fstat)
#else 
    open(199,file=fn1,form='unformatted',status='replace',iostat=fstat)
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening global velocity field'
      write(*,*) 'file:',fn1
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

! now we need to read in data from each node, and write it appropriately.
! can't read in from nodes sequentially, have to jump around.

  do k=1,rmesh*nodes_dim
    do j=1,rmesh*nodes_dim
      do i=1,rmesh*nodes_dim,rmesh
        do rank=0,nn-1
          if ( node_coords(1,rank)==(i-1)/rmesh.and. node_coords(2,rank)==(j-1)/rmesh.and. &
               node_coords(3,rank)==(k-1)/rmesh) then
!             print '(7i8)',i,j,k,rank,node_coords(:,rank)
            read(21+rank) varr
            write(199) varr
          endif
        enddo
      enddo
    enddo
  enddo

  close(199)
  do rank=0,nn-1
    close(21+rank)
  enddo

  deallocate(varr)

#endif

!! next do the fine clumping array 
!  real(4), dimension(3,nc_node_dim*mesh_scale/fine_clumping_scale,nc_node_dim*mesh_scale/fine_clumping_scale, &
!                     nc_node_dim*mesh_scale/fine_clumping_scale) :: velocity_field

  do rank=0,nn-1
    write(rank_s,'(i4)') rank
    rank_s=adjustl(rank_s)
    ofile=serial_output_path//z_s(1:len_trim(z_s))//'fc'// &
    rank_s(1:len_trim(rank_s))//'.dat'

#ifdef BINARY
    open(unit=21+rank,file=ofile,status='old',iostat=fstat,form='binary')
#else
    open(unit=21+rank,file=ofile,status='old',iostat=fstat,form='unformatted')
#endif

    if (fstat /= 0) then
      write(*,*) 'error opening fine clumping checkpoint'
      write(*,*) 'rank',rank,'file:',ofile
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

  enddo

  allocate(rarr(rmesh))

! now we need to read in data from each node, and write it appropriately.
! can't read in from nodes sequentially, have to jump around.

    fn1=serial_output_path(1:len_trim(serial_output_path))//z_s(1:len_trim(z_s))//'fc.dat'
#ifdef BINARY
    open(199,file=fn1,form='binary',status='replace',iostat=fstat)
#else 
    open(199,file=fn1,form='unformatted',status='replace',iostat=fstat)
#endif
    if (fstat /= 0) then
      write(*,*) 'error opening global fine clumping file'
      write(*,*) 'file:',fn1
      stop !call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

! now we need to read in data from each node, and write it appropriately.
! can't read in from nodes sequentially, have to jump around.

  do k=1,rmesh*nodes_dim
    do j=1,rmesh*nodes_dim
      do i=1,rmesh*nodes_dim,rmesh
        do rank=0,nn-1
          if ( node_coords(1,rank)==(i-1)/rmesh.and. node_coords(2,rank)==(j-1)/rmesh.and. &
               node_coords(3,rank)==(k-1)/rmesh) then
!             print '(7i8)',i,j,k,rank,node_coords(:,rank)
            read(21+rank) rarr
            write(199) rarr
          endif
        enddo
      enddo
    enddo
  enddo

  close(199)
  do rank=0,nn-1
    close(21+rank)
  enddo

  deallocate(rarr)

enddo

end 
