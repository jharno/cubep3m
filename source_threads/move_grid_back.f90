!! move grid back
  subroutine move_grid_back
    implicit none
#ifdef DISP_MESH 
    include 'mpif.h'
#endif
#    include "cubepm.fh"

    integer(4) :: i
    real(4), dimension(3) :: remove_offset
    logical :: doremove
    real(4), parameter :: maxshake = 16. 

#ifdef DISP_MESH 

    call mpi_bcast(shake_offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call system_clock(count=count_i)

    remove_offset(:) = 0.
    doremove = .false.

    if (force_grid_back) then

        doremove = .true.
        remove_offset = shake_offset

    else

        do i = 1, 3
            if (abs(shake_offset(i)) > maxshake) then
                if (shake_offset(i) > 0.) then
                    remove_offset(i) = maxshake 
                else
                    remove_offset(i) = -maxshake 
                endif
                doremove = .true.
            endif
        enddo

    endif

    if (doremove) then

        if (rank == 0) write(*,*) "remove_offset = ", remove_offset
    
        !$omp parallel do default(shared) private(i)
        do i=1,np_local
            xv(1:3,i)=xv(1:3,i)-remove_offset(:)
        enddo
        !$omp end parallel do
        shake_offset = shake_offset - remove_offset 

        call link_list
        call particle_pass
        call delete_particles

        if (rank == 0) write(*,*) "moved shake_offset = ", shake_offset

    endif

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('move_grid_back',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'move grid finished',real(count_f-count_i)/real(count_r)
#endif

#else
    if(rank==0) write(*,*)  '*** Could not move back, no off set to start with! ***'
#endif

  end subroutine move_grid_back
