! -------------------------------------------------------------------------------------------------------
! Writes all particles within some specified cubic region on a specified node every x number of timesteps.
! -------------------------------------------------------------------------------------------------------
subroutine zoomcheckpoint

    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    integer, parameter :: rank2write = 0
    integer, parameter :: writeEverySteps = 2
    real, parameter :: width = 90.
    real, dimension(3), parameter :: centre = (/46.9507, 240.598, 190.069/)

    character (len=max_path) :: ofile
    character (len=5) :: step_string 
#ifdef NEUTRINOS
    character (len=max_path) :: ofile_nu
#endif
    real, dimension(3) :: x
    real, dimension(3) :: xmin
    real, dimension(3) :: xmax
    integer :: i
    integer :: fstat

    if (rank == rank2write .and. mod(nts,writeEverySteps) == 0) then

        write(*,*) "Writing zoomcheckpoint for nts = ", nts

        !! Store bounding values
        do i = 1, 3
            xmin(i) = centre(i) - width/2.
            xmax(i) = centre(i) + width/2.
            if (xmin(i) < 0. .or. xmax(i) > real(nf_physical_node_dim)) then
                write(*,*) "WARNING: Bad choice of centre and width in zoomcheckpoint: ", centre, width
            endif
        enddo

        !! Attach sweep number to the output file(s)
        write(step_string, "(i5.5)") nts

        !! Open output file(s)         
        ofile = output_path//"zoomcheckpoint/partdata"//step_string//".bin"
        open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error opening zoomcheckpoint file for write'
            write(*,*) 'rank',rank,'file:',ofile
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
#ifdef NEUTRINOS
        ofile_nu = output_path//"zoomcheckpoint/partdata"//step_string//"_nu.bin"
        open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
            write(*,*) 'error opening zoomcheckpoint file for write'
            write(*,*) 'rank',rank,'file:',ofile_nu
            call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
#endif

        !! Write only those particles within the region of interest
        do i = 1, np_local
#ifdef DISP_MESH
            x = xv(1:3,i) - shake_offset            
#else
            x = xv(1:3,i)
#endif            
            if (x(1) >= xmin(1) .and. x(1) <= xmax(1) .and. x(2) >= xmin(2) &
                .and. x(2) <= xmax(2) .and. x(3) >= xmin(3) .and. x(3) <= xmax(3)) then
#ifdef NEUTRINOS
                if (PID(i) == 1) then
#endif
                    write(12) x(:) - xmin(:)
                    write(12) xv(4:6,i)
#ifdef NEUTRINOS
                else
                    write(22) x(:) - xmin(:)
                    write(22) xv(4:6,i)
                endif
#endif
            endif
        enddo
        
        !! Close output file(s)
        close(12)
#ifdef NEUTRINOS
        close(22)
#endif

    endif 

    call mpi_barrier(mpi_comm_world, ierr)

end subroutine zoomcheckpoint
