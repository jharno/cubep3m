! -------------------------------------------------------------------------------------------------------
! Writes all particles within some specified subvolume of the entire simulation 
! -------------------------------------------------------------------------------------------------------
subroutine boxcheckpoint

    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    character (len=max_path) :: ofile
    character (len=5) :: step_string 
    character (len=6) :: rank_s
#ifdef NEUTRINOS
    character (len=max_path) :: ofile_nu
#endif
    real, dimension(3) :: x
    integer :: i
    integer :: np_dm, np_nu
    integer :: fstat

    if (any(rank .eq. nsubcuberanks) .and. mod(nts,writeBoxEverySteps) == 0 .and. a >= writeBoxAboveA) then

        if (rank == 0) write(*,*) "Writing boxcheckpoint for nts = ", nts

        !! Attach rank number to the output file(s) 
        write(rank_s,'(i6)') rank
        rank_s=adjustl(rank_s)

        !! Attach sweep number to the output file(s)
        write(step_string, "(i5.5)") nts

        !! Open output file(s)         
        ofile = output_path//"boxcheckpoint/node"//rank_s(1:len_trim(rank_s))//"/xv"//rank_s(1:len_trim(rank_s))//"_"//step_string//".bin"
        open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error opening boxcheckpoint file for write'
            write(*,*) 'rank',rank,'file:',ofile
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
#ifdef NEUTRINOS
        ofile_nu = output_path//"boxcheckpoint/node"//rank_s(1:len_trim(rank_s))//"/xv"//rank_s(1:len_trim(rank_s))//"_"//step_string//"_nu.bin" 
        open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error opening boxcheckpoint file for write'
            write(*,*) 'rank',rank,'file:',ofile_nu
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
#endif

#ifdef NEUTRINOS
        !! Determine how many dark matter and neutrino particles this rank has
        np_dm = 0
        np_nu = 0

        do i = 1, np_local
#ifdef NUPID
            if (PID(i) == 0) then
#else
            if (PID(i) == 1) then
#endif
                np_dm = np_dm + 1
            else
                np_nu = np_nu + 1
            endif
        enddo

        write(12) np_dm
        write(22) np_nu
#else
        write(12) np_local
#endif

        !! Now write the data
        do i = 1, np_local
#ifdef DISP_MESH
            x = xv(1:3,i) - shake_offset            
#else
            x = xv(1:3,i)
#endif            

#ifdef NEUTRINOS
            if (PID(i) == 1) then
#endif
                write(12) x(:)
                write(12) xv(4:6,i)
#ifdef NEUTRINOS
            else
                write(22) x(:)
                write(22) xv(4:6,i)
            endif
#endif
        enddo
        
        !! Close output file(s)
        close(12)
#ifdef NEUTRINOS
        close(22)
#endif

    endif 

    call mpi_barrier(mpi_comm_world, ierr)

end subroutine boxcheckpoint
