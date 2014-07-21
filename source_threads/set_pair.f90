!! this routine sets particle pairs for testing purposes
  subroutine set_pair
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    integer(4) :: i,j
    real(4) :: rnum,sep(3)

    np_local=0
    xv(:,:2)=0.0

    if (rank==0) then
      do i=1,3
        call random_number(rnum)
        if (rnum>=0.9999) rnum=0.995
        xv(i,1)=rnum*real(nf_physical_dim,4)
      enddo
      call random_number(rnum)
      if (.not.pair_infall) then
        cur_sep=modulo(cur_sep+sep_inc*(rnum+0.5),max_sep)
        if (cur_sep < min_sep ) cur_sep = min_sep+sep_inc*(rnum+0.5)
      endif
      print *,'current seperation=',cur_sep
      call random_number(rnum)
      sep(1)=cur_sep*2.0*(rnum-0.5)
      call random_number(rnum)
      sep(2)=sqrt(cur_sep**2 - sep(1)**2)*2.0*(rnum-0.5)
      call random_number(rnum)
      if (rnum<0.5) then
        sep(3)=-sqrt(cur_sep**2-sep(1)**2-sep(2)**2)
      else
        sep(3)=sqrt(cur_sep**2-sep(1)**2-sep(2)**2)
      endif
      do i=1,3
        xv(i,2)=modulo(xv(i,1)+sep(i),real(nf_physical_dim,4))
      enddo
if (pair_infall) then
  xv(:3,1)=(/150.0,150.0,149.75/)
  xv(:3,2)=(/150.0,150.0,149.75+cur_sep/)
endif
if (pairwise_ic) then
if (nts==1) then
xv(:3,1)=(/  34.65000153 ,   60.22747803  ,  46.03750229 /) 
xv(:3,2)=(/  34.91682053 ,   59.85746002  ,  45.87303162 /)
!xv(:3,1)=(/ 78.72641754 , 147.02244568 ,94.36387634 /)
!!xv(:3,2)=(/57.30688477 ,   151.37088013 ,87.10507202 /) 
!xv(:3,1)=(/   227.84477234  , 139.25076294  ,  35.63391876/)
!xv(:3,2)=(/ 244.55119324 ,  137.25819397 ,  32.75370789/)
!  xv(:3,1)=(/216.41,5.40,85.01/)
!  xv(:3,2)=(/216.46,5.37,85.06/)
!!  xv(:3,1)=(/30.700,30.700,30.700/)
!!  xv(:3,2)=(/30.900,30.900,30.900/)
endif
endif
  !    if (nts.eq.1) then
  !      xv(:3,1)=(/311.24270630,288.43865967,263.88589478/)
  !      xv(:3,2)=(/311.30111694,288.31353760,264.10906982/)
  !    endif

!      do j=1,2
!        if (nts.eq.1) then
!           xv(:3,1)=(/22.0,22.0,22.0/)
!           xv(:3,2)=(/62.0,22.0,22.0/)
!          xv(:3,1)=(/304.35351562,127.90506744,168.49984741/) 
!          xv(:3,2)=(/319.88720703,122.08447266,167.59216309/)
!          xv(:3,1)=(/290.64648438,145.65747070,237.38870239/)
!          xv(:3,2)=(/300.61291504,143.01974487,242.72528076/)
!        else
!          call random_seed()
!          do i=1,3
!            call random_number(rnum)
!            if (rnum.ge.0.9999) rnum=0.995
!            if (j==1) then
!              xv(i,j)=rnum*nf_physical_dim
!            else
!              xv(i,j)=mod(xv(i,1)+(rnum-0.5)*3.0*(mod(nts,14)+1.2)+ &
!                      real(nf_physical_dim,kind=4),real(nf_physical_dim,kind=4)) 
!            endif
!          enddo
!        endif
!        xv(4:,j)=0.0
!      enddo
 endif

    call mpi_bcast(xv(1:6,1:2),12,mpi_real,0,mpi_comm_world,ierr)

    if (xv(1,1) >= cart_coords(3)*nf_physical_node_dim.and.xv(1,1) &
               < (cart_coords(3)+1)*nf_physical_node_dim) then
      if (xv(2,1) >= cart_coords(2)*nf_physical_node_dim.and.xv(2,1) &
                 < (cart_coords(2)+1)*nf_physical_node_dim) then
        if (xv(3,1) >= cart_coords(1)*nf_physical_node_dim.and.xv(3,1) &
                   < (cart_coords(1)+1)*nf_physical_node_dim) then
          np_local=1
        endif
      endif
    endif        

    if (xv(1,2) >= cart_coords(3)*nf_physical_node_dim.and.xv(1,2) &
               < (cart_coords(3)+1)*nf_physical_node_dim) then
      if (xv(2,2) >= cart_coords(2)*nf_physical_node_dim.and.xv(2,2) &
                 < (cart_coords(2)+1)*nf_physical_node_dim) then
        if (xv(3,2) >= cart_coords(1)*nf_physical_node_dim.and.xv(3,2) &
                   < (cart_coords(1)+1)*nf_physical_node_dim) then
          if (np_local == 0) then
            np_local=1
            xv(:,1)=xv(:,2)
          else 
            np_local=2
          endif
        endif
      endif
    endif        
   
    do i=1,np_local
#ifdef DIAG
      write(*,*) rank, xv(:,i)
#endif
      xv(1,i)=xv(1,i)-cart_coords(3)*nf_physical_node_dim
      xv(2,i)=xv(2,i)-cart_coords(2)*nf_physical_node_dim
      xv(3,i)=xv(3,i)-cart_coords(1)*nf_physical_node_dim
    enddo

  end subroutine set_pair 
