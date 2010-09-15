
!! NEED TOO MUCH MEMORY ! DON'T RUN IT ON RICKY/BUBBLES !
!! modified a bit on Ilian's code

program serial_halo_merge
  !! This is the serial cubepm halo catalog merging routine.
  !! Accumulates distributed distributed halo catalogs into monolithic file 
  !! compile with: ifort indexedsort.f90 serial_halo_merge.f90 -o s_halo_merge [-fpp -DPID_FLAG optional]
  !! Can't Merge with the PID for now...

  implicit none
!  include '/smaug/data/theorie/iliev/cubep3m_threads_largerun/parameters'
  !include '/cita/d/scratch-2week/jharno/workdir/parameters'
  include '../../parameters'


  character(len=*),parameter :: halofinds ='/cita/d/scratch-2week/jharno/workdir/input/halofinds' 
  character(len=*),parameter :: serial_output_path=output_path !'/cita/d/scratch-2week/jharno/cubep3m_new_V1/out/RUN-1/' !halo_data/'
  integer, parameter :: nn=nodes_dim**3 !! number of nodes total
  integer, parameter :: max_input=100 !! maximum number of catalogs to merge 
  integer, parameter :: max_halos=10*128**3 !!maximum number of halos / catalog
  
  integer, parameter :: nhv=20  !! number of halo quantities in each catalog entry
!  integer, parameter :: nhv=17  !! number of halo quantities in each catalog entry
  
  real(4), dimension(nhv) :: halo_input_buffer
  
  
  real(4), dimension(nhv,max_halos) :: halo_list,halo_list0
  integer(4), dimension(max_halos) :: isorthalo
  real(4), dimension(max_halos) :: mass
  
  character(len=180) ofile 
  character(len=5) rank_s
  character(len=7) z_s
  integer i,j,k,m
  real x,y,z,z_write
  
  integer(4) :: thread,nn_returned, ierr,rank, num_halofinds,fstat,cur_halofind
  real(4), dimension(max_input) :: z_halofind
  
  integer(4) :: nh_local,nh_sum,nh_in
  integer(4) :: nodes_returned,tag
  real(4) :: nh_total
  
  integer(4), dimension(0:nn-1) :: nh_found
  !integer(4), dimension(mpi_status_size) :: status

#ifdef PID_FLAG    
    integer, parameter :: N_p = 10
    integer(8), dimension(N_p) :: pid_halo
    integer(8), dimension(N_p,max_halos) :: pid_list, pid_list0
#endif
 

 
  common / rarr / halo_list,mass
  common / iarr / isorthalo
  
  open(11,file=halofinds,status='old',iostat=fstat)
  if (fstat /= 0) then
     write(*,*) 'error opening halofind list file'
     write(*,*) 'rank',rank,'file:',halofinds
     stop !call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
  do num_halofinds=1,max_input
     read(unit=11,err=51,end=41,fmt='(f20.10)') z_halofind(num_halofinds)
  enddo
41 num_halofinds=num_halofinds-1
51 close(11)
  write(*,*) 'halofinds to recompose:'
  do i=1,num_halofinds
     !    write(*,'(f5.1)') z_halofind(i)
     write(*,'(f20.10)') z_halofind(i)
  enddo
  
  !call mpi_bcast(num_halofinds,1,mpi_integer,0,mpi_comm_world,ierr)
  
  do cur_halofind=1,num_halofinds
     
     !! Read in halos 
     nh_local=0
     nh_found=0
     nh_sum=0  
     
     !  if (rank == 0) then
     z_write = z_halofind(cur_halofind)
     write(*,*) 'processing z=',z_write
     !  endif
     
     !  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
     
     write(z_s,'(f7.3)') z_write
     z_s=adjustl(z_s)

     do rank=0,nn-1
        write(rank_s,'(i5)') rank
        rank_s=adjustl(rank_s)
        
        ofile=serial_output_path//z_s(1:len_trim(z_s))//'halo'// &
             rank_s(1:len_trim(rank_s))//'.dat'
        !#ifdef BINARY
        open(unit=21,file=ofile,status='old',iostat=fstat,form='binary')
        !#else
        !  open(unit=21,file=ofile,status='old',iostat=fstat,form='unformatted')
        !#endif
        !    open(unit=21,file=ofile,status='old',iostat=fstat,form='formatted')
        if (fstat /= 0) then
           write(*,*) 'error opening catalog'
           write(*,*) 'rank',rank,'file:',ofile
           stop !call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        
        read(21) nh_in
        print*,'num. halos=',nh_in
        nh_local=0
        do
           !      read(21,end=112,fmt='(6f20.10)') halo_list(:,nh_local+1)
           !      read(21,end=112,fmt='(16f20.10)') halo_list(:,nh_local+1)
#ifdef PID_FLAG
           read(21,end=112,err=113) halo_input_buffer, pid_halo
#else
           read(21,end=112,err=113) halo_input_buffer
#endif
!           if(halo_input_buffer(16)<=0.0)then
!              print*,halo_input_buffer
!              pause
!           end if
           nh_local=nh_local+1
           halo_list0(:,nh_local)=halo_input_buffer

#ifdef PID_FLAG
           !pid_list0(:,nh_local)=pid_halo(:)
#endif

113        continue
        enddo
112     close(21)
        if (nh_in /= nh_local) then
           write(*,*) 'expected',nh_in,' halos, found ',nh_local
           !    call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        if (nh_local >= max_halos) then
           write(*,*) 'too many halos to read in'
           stop !call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        nh_found(rank)=nh_local
        halo_list(:,sum(nh_found)-nh_local+1:sum(nh_found))=halo_list0(:,1:nh_local)
#ifdef PID_FLAG        
        !pid_list(:,sum(nh_found)-nh_local+1:sum(nh_found))=pid_list0(:,1:nh_local)
#endif

     enddo
     !  call mpi_gather(nh_local,1,mpi_integer,nh_found,1,mpi_integer,0, &
     !                  mpi_comm_world,ierr)
     
     !if (rank == 0) then
     do j=0,nn-1
        write(*,*) 'rank:',j,'number of halos=',nh_found(j)
        !      write(*,*) 'number of halos=',nh_local !nh_found(j)
     enddo
     nh_sum=sum(nh_found)
     write(*,*) 'total halos=',nh_sum
     if (nh_sum > max_halos) then
        print *,'too many halos to store'
        !      call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     !  endif
     
     !  do j=1,nn-1
     !    tag=j
     !    if (rank == 0) then
     !      call mpi_recv(halo_list(1,nh_local+1),6*nh_found(j),mpi_real, &
     !                    j,tag,mpi_comm_world,status,ierr)
     !    elseif (rank == j) then
     !      call mpi_send(halo_list,6*nh_local,mpi_real,0,tag,mpi_comm_world,ierr)
     !    endif
     !    call mpi_barrier(mpi_comm_world,ierr)
     !    if (rank ==0) nh_local=nh_local+nh_found(j)
     !  enddo
     
     !  if (rank == 0) then
     !    print *,'finished pass'
     !    nh_sum=nh_local
!     mass(:nh_sum)=halo_list(4,:nh_sum)
     mass(:nh_sum)=halo_list(15,:nh_sum)
     isorthalo(:nh_sum)=(/ (k,k=1,nh_sum) /)
     print*,'calling sort'
     call indexedsort(nh_sum,mass,isorthalo)
     print*,'finished sort'
     print*,'starting halo list',nh_sum
!     halo_list(:,:nh_sum)=halo_list(:,isorthalo(:nh_sum))
     print*,'done halo list',nh_sum
     ofile=serial_output_path//z_s(1:len_trim(z_s))//'halo.dat'
     open(unit=31,file=ofile,status='replace',iostat=fstat,form='formatted')
     if (fstat /= 0) then
        write(*,*) 'error opening output catalog'
        write(*,*) 'rank',rank,'file:',ofile
        stop !call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     print*,'nh_sum=',nh_sum
     do j=nh_sum,1,-1
        !      write(31,'(6f20.10)') halo_list(:,j)
        !        write(31,'(17f20.10)') halo_list(:,j)
!#ifdef PID_FLAG
!        if(halo_list(16,j)>0) write(31,*) halo_list(:,j),pid_list(:,j)
!#else        
        if(halo_list(16,j)>0) write(31,'(20f20.10)') halo_list(:,j)
!#endif
!        print*,'check',halo_list(:,j)
     enddo
!     pause
     close(31)
     !  endif
     
  enddo
  
end program serial_halo_merge
