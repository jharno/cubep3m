!! delete particles which lie outside of the physical region
  subroutine delete_particles
    implicit none
    include 'mpif.h'
#    include "cubepm.fh"

    integer(4) :: pp,i
#ifdef DIAG
    integer*8 :: npl8,np_total
#endif

    call system_clock(count=count_i)
    
    pp=1
    do
42     if (pp > np_local) exit
       if (xv(1,pp) >= nf_physical_node_dim .or. xv(1,pp) < 0.0 .or. &
            xv(2,pp) >= nf_physical_node_dim .or. xv(2,pp) < 0.0 .or. &
            xv(3,pp) >= nf_physical_node_dim .or. xv(3,pp) < 0.0) then
#ifdef DEBUG
          if (xv(1,pp) == nf_physical_node_dim .or. xv(1,pp) == 0.0 .or. &
               xv(2,pp) == nf_physical_node_dim .or. xv(2,pp) == 0.0 .or. &
               xv(3,pp) == nf_physical_node_dim .or. xv(3,pp) == 0.0) & 
               write(*,*) 'particle on boundry:',xv(:,pp)
#endif
          
#ifdef DEBUG_PID_INTENSE                 
          write(*,*) 'Moving particle', pp, 'with x =', xv(1:3,pp), ',PID=', PID(pp)
          write(*,*) 'to pp=',np_local, 'with x =', xv(1:3,np_local),',PID=',PID(np_local)   
          !pause
#endif
          
          xv(:,pp)=xv(:,np_local)
#ifdef PID_FLAG
          PID(pp)=PID(np_local)
#endif
          
#ifdef DEBUG_PID_INTENSE                 
          write(*,*) 'Moved particle', pp, 'with x =', xv(1:3,pp), ',PID=', PID(pp)
          write(*,*) 'to pp=',np_local, 'with x =', xv(1:3,np_local),',PID=',PID(np_local)   
          !pause
#endif
          
          
          
          np_local=np_local-1
          goto 42
       endif
       pp=pp+1
    enddo
    
#ifdef DEBUG
    do pp=0,nodes-1
       if (rank==pp) then
          print *, 'rank:',rank,'np_local=',np_local
          print '(6f16.8)', xv(:,np_local+1)
       endif
       call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif
#ifdef DIAG
    npl8=int(np_local,kind=8)
    call mpi_reduce(npl8,np_total,1,MPI_INTEGER8,mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'total number of particles =', np_total   
#endif
    
#ifdef DEBUG_PID
    write(*,*) 'Particles after deleting;'
    do i= 1,1000
       write(*,*) i, 'PID=', PID(i), 'xv=',xv(1:3,i)
    enddo
#endif

    
    if (shake_test_ic) print *,'P0',xv(:,1)
    
    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('del part',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'particle deletion finished',real(count_f-count_i)/real(count_r) 
#endif
    
  end subroutine delete_particles
