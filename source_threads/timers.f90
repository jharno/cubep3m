!! start system clock timer
subroutine t_start(system_clock_count_init)

    implicit none

    integer :: system_clock_count_init

    call system_clock(count=system_clock_count_init)

  end subroutine t_start

!! calculate time elapsed since system_clock_count_init
  real(kind=4) function t_elapsed(system_clock_count_init)

    implicit none

    integer :: system_clock_count,system_clock_rate,system_clock_count_max
    integer, intent(in) :: system_clock_count_init

    call system_clock(count=system_clock_count, &
                    count_rate=system_clock_rate, count_max=system_clock_count_max)

    if (system_clock_count < system_clock_count_init) then
      t_elapsed=real(system_clock_count + system_clock_count_max - system_clock_count_init) &
               / real(system_clock_rate)
    else
      t_elapsed=real(system_clock_count - system_clock_count_init) &
               / real(system_clock_rate)
    endif

  end function t_elapsed

!! generic f90 system clock stopwatch timer
!! pass it 'i' to start, 't' to stop
real function timer(op)
    implicit none
    character(len=1), intent(in) :: op
    integer :: system_clock_count_init,system_clock_count,system_clock_rate,system_clock_count_max
    save system_clock_count_init, system_clock_rate, system_clock_count_max

    if (op=='i') then
      call system_clock(count=system_clock_count_init,count_rate=system_clock_rate,count_max=system_clock_count_max)
      timer=0.0
    elseif (op=='t') then
      call system_clock(count=system_clock_count)
      if (system_clock_count < system_clock_count_init) then
        timer=real(system_clock_count + system_clock_count_max - system_clock_count_init)/real(system_clock_rate)
      else
        timer=real(system_clock_count - system_clock_count_init)/real(system_clock_rate)
      endif
    else
      print *,'timer operation:',op,'not supported'
      timer=-42.0
    endif
end function timer

!! print out the date/time 
subroutine datestamp
  implicit none
  character(len=8) :: t
  character(len=8) :: td
  call time(t)
  call date_and_time(td)
  print *,td,' ',t
end subroutine datestamp

!! a small routine that prints the max/min/avg of 'local_time' accross multiple nodes
subroutine mpi_time_analyze(string_tag,local_time,myrank,processes)
  real :: local_time,max_time,min_time,avg_time
  integer :: processes,ierr,myrank
  character(len=8) :: string_tag
  include 'mpif.h'
  call mpi_reduce(local_time,max_time,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(local_time,min_time,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(local_time,avg_time,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
  if (myrank==0) print *,string_tag,' : ',max_time,avg_time/real(processes),min_time
end subroutine mpi_time_analyze

