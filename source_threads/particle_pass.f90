!! pass particles to adjacent nodes
  subroutine particle_pass
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    real(4), parameter :: rnf_buf = nf_buf
    integer(4) :: np_max
    integer(4) :: np_buf_max,np_buf_max_dir
    integer(4) :: np_local0

#ifdef DEBUG_PID
    real(4) :: np_total
    integer(4) :: np_local_i
#endif
    integer(4) :: i,j,k,pp,tag
    integer(4) :: nppx,nppy,nppz,npmx,npmy,npmz
    integer(4), dimension(mpi_status_size) :: status,sstatus,rstatus
    integer(4) :: srequest,rrequest,sierr,rierr
    integer(4) :: ikill, ikill_loc

    np_buf_max_dir = 0

#ifdef MPI_TIME
    call mpi_barrier(mpi_comm_world,ierr)
#endif

    call system_clock(count=count_i)

    !! Keep track of np_local at the start in case we need to checkpoint_kill
    np_local0 = np_local

    !! If this variable becomes 1 then we must checkpoint kill
    ikill_loc = 0

#ifdef DEBUG_PID
    np_local_i=np_local
    do i=0,nodes-1
      if (i == rank) print *, rank,'number of particles before pass=',np_local
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
    call mpi_reduce(real(np_local,kind=4),np_total,1,mpi_real, &
                           mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'total number of particles before pass=', int(np_total,8)
    j=0

    do i=1,np_local
       if (xv(1,i) >= (nf_physical_node_dim - rnf_buf) .or. xv(1,i) < rnf_buf .or. &
            xv(2,i) >= (nf_physical_node_dim - rnf_buf) .or. xv(2,i) < rnf_buf .or. &
            xv(3,i) >= (nf_physical_node_dim - rnf_buf) .or. xv(3,i) < rnf_buf) then
          !        print '(6f10.4)', xv(:,i)
          j=j+1
          !else
          !print '(6f10.4)', xv(:,i)
       endif
       if (xv(2,i) < 0.00001 .and. xv(2,i) > -0.00001) print *, i,rank,xv(:,i)            
       !    0.0000   12.2546   -0.0540   -0.1137   -0.0829
    enddo
    call mpi_reduce(j,k,1,mpi_integer, &
                           mpi_sum,0,mpi_comm_world,ierr)
    do i=0,nodes-1
      if (i == rank) print *, rank,'number of particles leaving=',j
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
    if (rank == 0) print *,'total out=',k
#endif

    ! Uncomment this if you want to print max number of particles in any node, before particle pass (i.e. with empty buffers):
    !call mpi_reduce(np_local,np_max,1,mpi_integer, &
    !                       mpi_max,0,mpi_comm_world,ierr)
    !if(rank==0) write(*,*) '*** max np_local (no ghosts)   = ' , np_max, ' ***'

! pass +x

    tag=11
    np_buf = 0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_h - hoc_pass_depth, hoc_nc_h 
          pp = hoc(i,j,k)
#ifdef DEBUG
          if (pp == 0 .and. i<nc_node_dim .and. j<nc_node_dim .and. k<nc_node_dim &
                      .and. i>0 .and. j>0 .and. k>0 .and. grid_ic)  &
                      write(*,*) 'no particles in coarse cell',i,j,k
#endif
          do while (pp /= 0)
            if (xv(1,pp) >= nf_physical_node_dim - rnf_buf) then
              np_buf = np_buf + 1
              send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
#ifdef PID_FLAG
              send_buf_PID(np_buf)=PID(pp)
#endif
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo

    !! Check to see if we need to checkpoint kill
    if (np_buf*6 > max_buf) then
        write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
        ikill_loc = 1 
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (+x pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    nppx = np_buf

    call mpi_reduce(nppx,np_buf_max_dir,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    if(rank==0) then
       if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
    endif

#ifdef DEBUG_PID
     write(*,*) 'Before nppx exchange'
    do i=0,nodes-1
       if (i==rank) write(*,*) 'rank',rank,'np_out=',nppx,'+x'
       call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

#ifdef DEBUG_PID_INTENSE
    do i=1,nppx
       write(*,*) i,'send_buf_PID=',send_buf_PID(i), 'send_buf=', send_buf((i-1)*6+1:i*6 - 3)
    enddo
#endif

    call mpi_sendrecv_replace(nppx,1,mpi_integer,cart_neighbor(6), &
                               tag,cart_neighbor(5),tag,mpi_comm_world, &
                               status,ierr)
#ifdef DEBUG_PID
    write(*,*) 'After nppx exchange'
    do i=0,nodes-1
      if (rank==i) write(*,*) 'rank',rank,'nppx=',nppx
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Check to see if we need to checkpoint kill
    if (np_local+nppx > max_np) then
        write(*,*) 'rank:',rank,'exceeded max_np in pass',nppx+np_local,max_np
        ikill_loc = 1
    endif 
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then 
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (+x pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,nppx*6,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,rrequest,rierr)
#ifdef PID_FLAG

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

#ifdef NEUTRINOS
#ifdef NUPID
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(6), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppx,MPI_integer8,cart_neighbor(5), &
         tag,mpi_comm_world,rrequest,rierr)
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(6), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppx,MPI_integer1,cart_neighbor(5), &
         tag,mpi_comm_world,rrequest,rierr)
#endif
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(6), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppx,MPI_integer8,cart_neighbor(5), &
         tag,mpi_comm_world,rrequest,rierr)
#endif

#endif

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

 
    do i=1,nppx
      xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
      xv(1,np_local+i)=max(xv(1,np_local+i)-nf_physical_node_dim,-rnf_buf)
#ifdef PID_FLAG
      PID(np_local+i)=recv_buf_PID(i)
#endif
    enddo

    np_local=np_local+nppx

#ifdef DEBUG_PID_INTENSE
    write(*,*) 'After x+ pass :' 
    do i = np_local-nppx+1,np_local
       write(*,*) i,'PID', PID(i),'xv',xv(1:3,i)
    enddo
#endif

! pass -x

    np_buf = 0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_l,hoc_nc_l + hoc_pass_depth
          pp = hoc(i,j,k)
          do while (pp /= 0)
            if (xv(1,pp) < rnf_buf) then
              np_buf = np_buf + 1
              send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
#ifdef PID_FLAG
              send_buf_PID(np_buf)=PID(pp)
#endif
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo

    !! Check to see if we need to checkpoint kill
    if (np_buf*6 > max_buf) then
        write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (-x pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    npmx = np_buf

    call mpi_reduce(npmx,np_buf_max_dir,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    if(rank==0) then
       if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'np_out=',npmx,'-x'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_sendrecv_replace(npmx,1,mpi_integer,cart_neighbor(5), &
                               tag,cart_neighbor(6),tag,mpi_comm_world, &
                               status,ierr)

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'npmx=',npmx
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Check to see if we need to checkpoint kill
    if (np_local+npmx > max_np) then
        write(*,*) 'rank:',rank,'exceeded max_np in pass',npmx+np_local,max_np
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (-x pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(5), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npmx*6,mpi_real,cart_neighbor(6), &
                   tag,mpi_comm_world,rrequest,rierr)
#ifdef PID_FLAG

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

#ifdef NEUTRINOS
#ifdef NUPID
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(5), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmx,MPI_integer8,cart_neighbor(6), &
         tag,mpi_comm_world,rrequest,rierr)
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(5), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmx,MPI_integer1,cart_neighbor(6), &
         tag,mpi_comm_world,rrequest,rierr)
#endif
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(5), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmx,MPI_integer8,cart_neighbor(6), &
         tag,mpi_comm_world,rrequest,rierr)
#endif

#endif
    
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)


    do i=1,npmx
      xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
#ifdef PID_FLAG
      PID(np_local+i)=recv_buf_PID(i)
#endif
      if (abs(xv(1,np_local+i)).lt.eps) then
        if (xv(1,np_local+i) < 0.0) then
          xv(1,np_local+i)=-eps
        else
          xv(1,np_local+i)=eps
        endif
      endif
      xv(1,np_local+i)=min(xv(1,np_local+i)+real(nf_physical_node_dim,4), &
                       nf_physical_node_dim+rnf_buf-eps)
    enddo

    np_local=np_local+npmx


! add additional particles to linked list
!! should add/subtract offsets here!

     pp=np_local-npmx-nppx+1
     do
91 continue
      if (pp > np_local) exit
      i=floor(xv(1,pp)/mesh_scale)+1
      j=floor(xv(2,pp)/mesh_scale)+1
      k=floor(xv(3,pp)/mesh_scale)+1
#ifdef DIAG
      if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        write(*,*) 'x-pass particle out of link list range - deleted'
        write(*,*) rank,xv(:,pp),pp,i,j,k,hoc_nc_l,hoc_nc_h
        xv(:,pp)=xv(:,np_local)
#ifdef PID_FLAG
        PID(pp)=PID(np_local)
#endif
        np_local=np_local-1
        goto 91
      endif
#endif
      ll(pp)=hoc(i,j,k)
      hoc(i,j,k)=pp
      pp=pp+1
    enddo

! pass -y

    np_buf = 0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_l,hoc_nc_l + hoc_pass_depth
        do i=hoc_nc_l,hoc_nc_h
          pp = hoc(i,j,k)
          do while (pp /= 0)
            if (xv(2,pp) < rnf_buf) then
              np_buf = np_buf + 1
              send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
#ifdef PID_FLAG
              send_buf_PID(np_buf)=PID(pp)
#endif
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo

    !! Check to see if we need to checkpoint kill
    if (np_buf*6 > max_buf) then
        write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (-y pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    npmy = np_buf

    call mpi_reduce(npmy,np_buf_max_dir,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    if(rank==0) then
       if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'np_out=',npmy,'-y'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_sendrecv_replace(npmy,1,mpi_integer,cart_neighbor(3), &
                               tag,cart_neighbor(4),tag,mpi_comm_world, &
                               status,ierr)

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'npmy=',npmy
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Check to see if we need to checkpoint kill
    if (np_local+npmy > max_np) then
        write(*,*) 'rank:',rank,'exceeded max_np in pass',npmy+np_local,max_np
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (-y pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) rank,'finished modifying buffer'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npmy*6,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,rrequest,rierr)
#ifdef PID_FLAG

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

#ifdef NEUTRINOS
#ifdef NUPID
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(3), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmy,MPI_integer8,cart_neighbor(4), &
         tag,mpi_comm_world,rrequest,rierr)
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(3), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmy,MPI_integer1,cart_neighbor(4), &
         tag,mpi_comm_world,rrequest,rierr)
#endif
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(3), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmy,MPI_integer8,cart_neighbor(4), &
         tag,mpi_comm_world,rrequest,rierr)
#endif

#endif

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) rank,'finished sending buffer'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    do i=1,npmy
      xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
#ifdef PID_FLAG
      PID(np_local+i)=recv_buf_PID(i)
#endif
      if (abs(xv(2,np_local+i)).lt.eps) then
        if (xv(2,np_local+i) < 0.0) then
          xv(2,np_local+i)=-eps
        else
          xv(2,np_local+i)=eps
        endif
      endif
      xv(2,np_local+i)=min(xv(2,np_local+i)+real(nf_physical_node_dim,4), &
                       nf_physical_node_dim+rnf_buf-eps)
    enddo
    np_local=np_local+npmy


! pass +y

   tag=11
    np_buf = 0
    do k=hoc_nc_l,hoc_nc_h
      do j=hoc_nc_h - hoc_pass_depth, hoc_nc_h 
        do i=hoc_nc_l,hoc_nc_h
          pp = hoc(i,j,k)
          do while (pp /= 0)
            if (xv(2,pp) >= nf_physical_node_dim - rnf_buf) then
              np_buf = np_buf + 1
              send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
#ifdef PID_FLAG
              send_buf_PID(np_buf)=PID(pp)
#endif  
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo

    !! Check to see if we need to checkpoint kill
    if (np_buf*6 > max_buf) then
        write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (+y pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    nppy = np_buf

    call mpi_reduce(nppy,np_buf_max_dir,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    if(rank==0) then
       if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'np_out=',nppy,'+y'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_sendrecv_replace(nppy,1,mpi_integer,cart_neighbor(4), &
                               tag,cart_neighbor(3),tag,mpi_comm_world, &
                               status,ierr)

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'nppy=',nppy
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Check to see if we need to checkpoint kill
    if (np_local+nppy > max_np) then
        write(*,*) 'rank:',rank,'exceeded max_np in pass',nppy+np_local,max_np
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (+y pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(4), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,nppy*6,mpi_real,cart_neighbor(3), &
                   tag,mpi_comm_world,rrequest,rierr)
#ifdef PID_FLAG

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

#ifdef NEUTRINOS
#ifdef NUPID
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(4), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppy,MPI_integer8,cart_neighbor(3), &
         tag,mpi_comm_world,rrequest,rierr)
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(4), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppy,MPI_integer1,cart_neighbor(3), &
         tag,mpi_comm_world,rrequest,rierr)
#endif
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(4), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppy,MPI_integer8,cart_neighbor(3), &
         tag,mpi_comm_world,rrequest,rierr)
#endif

#endif

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)


    do i=1,nppy
      xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
      xv(2,np_local+i)=max(xv(2,np_local+i)-nf_physical_node_dim,-rnf_buf)
#ifdef PID_FLAG 
      PID(np_local+i)=recv_buf_PID(i)
#endif
    enddo
    np_local=np_local+nppy



! add additional particles to linked list

     pp=np_local-npmy-nppy+1
     do
92 continue
      if (pp > np_local) exit
      i=floor(xv(1,pp)/mesh_scale)+1
      j=floor(xv(2,pp)/mesh_scale)+1
      k=floor(xv(3,pp)/mesh_scale)+1
#ifdef DIAG
      if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        write(*,*) 'y-pass particle out of link list range - deleted'
        write(*,*) rank,xv(:,pp),pp,i,j,k,hoc_nc_l,hoc_nc_h
        xv(:,pp)=xv(:,np_local)
#ifdef PID_FLAG 
        PID(pp)=PID(np_local)
#endif
        np_local=np_local-1
        goto 92
      endif
#endif
      ll(pp)=hoc(i,j,k)
      hoc(i,j,k)=pp
      pp=pp+1
    enddo

! pass +z

    tag=11
    np_buf = 0
    do k=hoc_nc_h - hoc_pass_depth, hoc_nc_h 
      do j=hoc_nc_l, hoc_nc_h
        do i=hoc_nc_l, hoc_nc_h
          pp = hoc(i,j,k)
          do while (pp /= 0)
            if (xv(3,pp) >= nf_physical_node_dim - rnf_buf) then
              np_buf = np_buf + 1
              send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
#ifdef PID_FLAG
              send_buf_PID(np_buf)=PID(pp)
#endif
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo

    !! Check to see if we need to checkpoint kill
    if (np_buf*6 > max_buf) then
        write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (+z pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    nppz = np_buf

    call mpi_reduce(nppz,np_buf_max_dir,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    if(rank==0) then
       if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'np_out=',nppz,'+z'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_sendrecv_replace(nppz,1,mpi_integer,cart_neighbor(2), &
                               tag,cart_neighbor(1),tag,mpi_comm_world, &
                               status,ierr)

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'nppz=',nppz
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Check to see if we need to checkpoint kill
    if (np_local+nppz > max_np) then
        write(*,*) 'rank:',rank,'exceeded max_np in pass',nppz+np_local,max_np
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (+z pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,nppz*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
#ifdef PID_FLAG

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

#ifdef NEUTRINOS
#ifdef NUPID
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(2), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppz,MPI_integer8,cart_neighbor(1), &
         tag,mpi_comm_world,rrequest,rierr)
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(2), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppz,MPI_integer1,cart_neighbor(1), &
         tag,mpi_comm_world,rrequest,rierr)
#endif
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(2), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,nppz,MPI_integer8,cart_neighbor(1), &
         tag,mpi_comm_world,rrequest,rierr)
#endif

#endif

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)


    do i=1,nppz
      xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
      xv(3,np_local+i)=max(xv(3,np_local+i)-nf_physical_node_dim,-rnf_buf)
#ifdef PID_FLAG
      PID(np_local+i)=recv_buf_PID(i)
#endif
    enddo

    np_local=np_local+nppz

! pass -z

    np_buf = 0
    do k=hoc_nc_l,hoc_nc_l + hoc_pass_depth
      do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_l,hoc_nc_h
          pp = hoc(i,j,k)
          do while (pp /= 0)
            if (xv(3,pp) < rnf_buf) then
              np_buf = np_buf + 1
              send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
#ifdef PID_FLAG 
              send_buf_PID(np_buf)=PID(pp)
#endif
            endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo

    !! Check to see if we need to checkpoint kill
    if (np_buf*6 > max_buf) then
        write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (-z pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    npmz = np_buf

    call mpi_reduce(npmz,np_buf_max_dir,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    if(rank==0) then
       if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
    endif

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'np_out=',npmz,'-z'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_sendrecv_replace(npmz,1,mpi_integer,cart_neighbor(1), &
                               tag,cart_neighbor(2),tag,mpi_comm_world, &
                               status,ierr)

#ifdef DEBUG
    do i=0,nodes-1
      if (i==rank) write(*,*) 'rank',rank,'npmz=',npmz
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Check to see if we need to checkpoint kill
    if (np_local+npmz > max_np) then
        write(*,*) 'rank:',rank,'exceeded max_np in pass',npmz+np_local,max_np
        ikill_loc = 1
    endif
    call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (ikill > 0) then
        if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting (-z pass) ... "
        !! Reset np_local to its starting point so that we don't write duplicates. 
        np_local = np_local0
#if defined(ZIP) || defined(ZIPDM)
        call checkpoint_kill(.false.)
#else   
        call checkpoint_kill
#endif
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npmz*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
#ifdef PID_FLAG

#ifdef NEUTRINOS
#ifdef NUPID
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(1), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmz,MPI_integer8,cart_neighbor(2), &
         tag,mpi_comm_world,rrequest,rierr)
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(1), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmz,MPI_integer1,cart_neighbor(2), &
         tag,mpi_comm_world,rrequest,rierr)
#endif
#else
    call mpi_isend(send_buf_PID,np_buf,MPI_integer8,cart_neighbor(1), &
         tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf_PID,npmz,MPI_integer8,cart_neighbor(2), &
         tag,mpi_comm_world,rrequest,rierr)
#endif

#endif

    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)


    do i=1,npmz
      xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
#ifdef PID_FLAG
      PID(np_local+i)=recv_buf_PID(i)
#endif
      if (abs(xv(3,np_local+i)).lt.eps) then
        if (xv(3,np_local+i) < 0.0) then
          xv(3,np_local+i)=-eps
        else
          xv(3,np_local+i)=eps
        endif
      endif 
      xv(3,np_local+i)=min(xv(3,np_local+i)+nf_physical_node_dim, &
                       nf_physical_node_dim+rnf_buf-eps)
    enddo

    np_local=np_local+npmz


! add additional particles to linked list

     pp=np_local-npmz-nppz+1
     do
93 continue
      if (pp > np_local) exit
      i=floor(xv(1,pp)/mesh_scale)+1
      j=floor(xv(2,pp)/mesh_scale)+1
      k=floor(xv(3,pp)/mesh_scale)+1
#ifdef DIAG
      if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        write(*,*) 'z-pass particle out of link list range - deleted'
        write(*,*) rank,xv(:,pp),pp,i,j,k,hoc_nc_l,hoc_nc_h
        xv(:,pp)=xv(:,np_local)
#ifdef PID_FLAG
        PID(pp)=PID(np_local)
#endif
        np_local=np_local-1
        goto 93
      endif
#endif
      ll(pp)=hoc(i,j,k)
      hoc(i,j,k)=pp
      pp=pp+1
    enddo

#ifdef DEBUG
    j=0
    do i=np_local_i+1,np_local
      if (xv(1,i) < nf_physical_node_dim .and. xv(1,i) >= 0.0 .and. &
          xv(2,i) < nf_physical_node_dim .and. xv(2,i) >= 0.0 .and. &
          xv(3,i) < nf_physical_node_dim .and. xv(3,i) >= 0.0) then
        j=j+1
!        print '(6f10.4)', xv(:,i)
      endif
    enddo 
    call mpi_reduce(j,k,1,mpi_integer, &
                           mpi_sum,0,mpi_comm_world,ierr)
    do i=0,nodes-1
      if (i == rank) print *, rank,'number of particles incoming=',j
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
    if (rank == 0) print *,'total in=',k
    do i=0,nodes-1
      if (rank == i) then
        write(*,*) 'rank:',rank,'np_in =',nppx+nppy+nppz+npmx+npmy+npmz, 'np_loc =', np_local
      endif
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_reduce(np_local,np_max,1,mpi_integer, &
                           mpi_max,0,mpi_comm_world,ierr)
    min_den_buf = real(np_max)*density_buffer/real(max_np)
    if(rank==0) write(*,*) '*************** Density_buffer Analysis *************'
    if(rank==0) write(*,*) '*** max np allowed             = ' , max_np, '    ***'
    if(rank==0) write(*,*) '*** max np_local (with ghosts) = ' , np_max, '    ***'
    !if(rank==0) write(*,*) '*** no density_buffer          = ', real(max_np)/density_buffer, ' ***'
    if(rank==0) write(*,*) '*** min density_buffer allowed = ', min_den_buf, ' ***'

    if(rank==0) write(*,*) '*************** SendRecv Analysis *******************'
    if(rank==0) write(*,*) '*** max np_buf                 = ' , np_buf_max, '    ***'
    if(rank==0) write(*,*) '*** max allowed                = ' , max_buf/6 , '    ***'
    if(rank==0) write(*,*) '*****************************************************'

!#ifdef WRITELOG
!    if(rank==0 .and. record_den_buf) then
!       write(unit=76,fmt='(f10.6)',advance='yes') real(np_max)*density_buffer/real(max_np)
!    endif
!#endif

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('par pass',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'particle pass finished',real(count_f-count_i)/real(count_r)
#endif

  end subroutine particle_pass

