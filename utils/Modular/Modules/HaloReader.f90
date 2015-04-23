!Fortran module to read HaloFiles produced by CUBEP3M
! May have to change number of parameters per halo at some point (e.g. from 35)
! May have to change header at some point

module HaloReader
  use Parameters
  use Variables
  use mMPI

  integer, parameter :: parameters_per_halo = 35
!  real(4), dimension(5) :: header
  real(4), dimension(2) :: header
contains

  !read_halo_file:
    !hfile - name of halo file to be read
    !halos - array to store the halo information (not the header)
    !max_nh - size of halos array
    !nh - actual number of halos - will be determined by read_halo_file

  subroutine read_halo_file(hfile,halos,nh)
    implicit none
    integer :: i, stat
    integer, intent(out) ::  nh    
    character(len=*),intent(in) :: hfile
    real(4), dimension(:,:), intent(out) :: halos
    integer(4) :: dummy

#if VERBOSITY > 0
    if (rank==0) write(*,*) '[Module -  HaloReader] Reading halos in file ', trim(hfile)
#endif 

    if (size(halos,dim=1) /= parameters_per_halo) call halo_error_stop('Error in subroutine read_halo_file: &
      & halo array dim 1 is incorrect')

    !Open File
    open(unit=11,file=trim(hfile),status='old',iostat=stat,access='stream')
    if (stat .NE. 0) call halo_error_stop('Error: could not open file: '//trim(hfile))

    !Check that number of halos is acceptable
    read(11) nh
    if (nh > size(halos,dim=2)) call halo_error_stop('Error: particle number exceeds max_nh')

    !Read header information (unused)
    read(11) header(:)

    !Read halos
    do i=1,nh
       read(11) halos(1:34,i)
       read(11) dummy
       halos(35,i) = real(dummy)
    end do
    close(11)

#if VERBOSITY > 0
    if (rank==0) write(*,*) '[Module - HaloReader] Finished reading halos'
#endif 
  end subroutine read_halo_file

  pure function local_coordinates(ghpos) result(lhpos)
    implicit none
    !integer :: k,j,i,index
    !integer, dimension(3) :: indices
    real(4), dimension(3), intent(in) :: ghpos
    real(4), dimension(3) :: lhpos
    lhpos(:) = ghpos(:) - slab_coord(:)*nc_node_dim
    !do k=1,nodes_dim
    !  do j=1,nodes_dim
    !    do i=1,nodes_dim
    !      index = (k-1)*nodes_dim**2 + (j-1)*nodes_dim + i - 1
    !      if (index .EQ. rank) then
    !        indices(1)=i
    !        indices(2)=j 
    !        indices(3)=k
    !      endif
    !    enddo
    !  enddo
    !enddo
    !lhpos = ghpos - (indices-1)*nc/nodes_dim
  end function local_coordinates

  subroutine halo_error_stop(expl)
    implicit none
    character(len=*), intent(in) :: expl
    write(*,*) '[Module - HaloReader]'
    write(*,*) '-->'//expl
    call mpi_abort(mpi_comm_world, ierr, ierr)
  end subroutine halo_error_stop

end module HaloReader

! -------------
! 1:3   = hpos=3
! 4     = mass_vir=1
! 5     = mass_odc=1
! 6     = r_vir=1
! 7     = r_odc=1
! 8:10  = x_mean=3
! 11:13 = v_mean=3
! 14:16 = l_CM=3
! 17:19 = v2_wrt_halo=3
! 20:22 = var_x=3
! 23:28 = I_ij=6
! 29:31 = x_mean_nu=3
! 32:34 = v_mean_nu=3
! 35    = n_nu=1
! -------------

