!! generate linked list for particles based on their positions
!! within the coarse mesh
  subroutine link_list
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    integer(4) :: i,j,k,pp,pc
    integer(4) :: omp_get_thread_num,omp_get_num_threads
    external omp_get_thread_num,omp_get_num_threads

    call system_clock(count=count_i)

#ifdef DEBUG_LOW
    write(*,*) 'rank',rank,'np_local',np_local
#endif

    hoc(:,:,:)=0
    np_buf=0

     pp=1
     do
!94 continue
       if (pp > np_local) exit
       i=floor(xv(1,pp)/mesh_scale)+1
       j=floor(xv(2,pp)/mesh_scale)+1
       k=floor(xv(3,pp)/mesh_scale)+1
       if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
           j < hoc_nc_l .or. j > hoc_nc_h .or. &
           k < hoc_nc_l .or. k > hoc_nc_h) then
         print*, 'link_list: particle moved out of buffer!',xv(:,pp)
         print*, 'check timestep & update_position!'
         call mpi_abort(mpi_comm_world,ierr,ierr)
         !xv(:,pp)=xv(:,np_local)

#ifdef PID_FLAG
         PID(pp)=PID(np_local)

#ifdef DEBUG_PID         
         write(*,*) 'Moved particle', pp, 'with PID', PID(pp),'to the end'  
         !pause
#endif

#endif
         np_local=np_local-1
         np_buf=np_buf+1
!         goto 94
         cycle
       else
         ll(pp)=hoc(i,j,k)
         hoc(i,j,k)=pp
       endif
       pp=pp+1
     enddo

!    !$omp do
!    do pp=1,np_local
!      i=floor(xv(1,pp)/mesh_scale)+1
!      j=floor(xv(2,pp)/mesh_scale)+1
!      k=floor(xv(3,pp)/mesh_scale)+1
!      if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
!          j < hoc_nc_l .or. j > hoc_nc_h .or. &
!          k < hoc_nc_l .or. k > hoc_nc_h) then
!        np_buf(thrdnum) = np_buf(thrdnum) + 1
!        fast_buf((np_buf(thrdnum)-1)*6+1:np_buf(thrdnum)*6,thrdnum)=xv(:,pp)
!        fast_pos((np_buf(thrdnum)-1)*4+1:np_buf(thrdnum)*4,thrdnum)=(/pp,i,j,k/)
!      else
!        if (hoc(i,j,k,thrdnum) == 0) toc(i,j,k,thrdnum)=pp
!        ll(pp)=hoc(i,j,k,thrdnum)
!        hoc(i,j,k,thrdnum)=pp
!      endif
!    enddo
!    !$omp end do

#ifdef DEBUG_LOW
    pc=0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_l,hoc_nc_h
          pp=hoc(i,j,k)
          do while (pp .ne. 0)
            pc=pc+1
            pp=ll(pp)
          enddo 
        enddo
      enddo
    enddo
    write(*,*) 'rank',rank,'particles in linked_list=',pc
#endif

#ifdef DIAG
    if (np_buf.gt.0) write(*,*) rank,'deleted',np_buf,'particles in ll'
#endif

!! since our linked list only points forward, we cannot implement something like
!! this here.  Need to restructure linked list creation.

!    do k=1,nt
!      do j=np_buf(k),1,-1
!
!        !! fast_buf((j-1)*6+1:j*6,k) is each particle
!        !! fast_pos((j-1)*4+1:j*4,k) is the position in the particle list followed
!        !!                            by the position in the coarse mesh
!
!        !! for now we will delete these particles.
!
!        xv(:,fast_pos((j-1)*4+1,k))=xv(:,np_local)
!        np_local=np_local-1          
!
!!! NEED TO FIX LINKED LIST NOW
!
!        !! need to calculate node to go to based on fast_pos 
!        !! send and recieve all incoming particles (look at continues scatter algorithm) 
!        !! add to particle and linked list
!
!      enddo
!    enddo

#ifdef DEBUG
    write(*,*) 'rank',rank,'particles deleted=',np_buf

    pc=0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_l,hoc_nc_h
          pp=hoc(i,j,k)
          do while (pp .ne. 0)
            pc=pc+1
            pp=ll(pp)
          enddo 
        enddo
      enddo
    enddo
    write(*,*) 'rank',rank,'particles in linked_list=',pc
    call mpi_reduce(pc,k,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'total number of particles in linked list=',k 
#endif

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('linklist',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'link list finished',real(count_f-count_i)/real(count_r)
#endif


  end subroutine link_list
