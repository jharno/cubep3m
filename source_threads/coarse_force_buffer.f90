!! pass coarse mesh force along boundries to adjacent nodes
  subroutine coarse_force_buffer
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    integer(4) :: buffer_size
    integer(4) :: tag
    integer(4) :: status(MPI_STATUS_SIZE)


    call system_clock(count=count_i)

    buffer_size = 3 * (nc_node_dim + 2)**2

    tag=64

    force_c_buffer(:,:,:)=force_c(:,1,:,:)
  !send to node in -x 
    call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
                              cart_neighbor(5),tag,cart_neighbor(6), &
                              tag,mpi_comm_cart,status,ierr)
    force_c(:,nc_node_dim+1,:,:)=force_c_buffer(:,:,:)

    force_c_buffer(:,:,:)=force_c(:,nc_node_dim,:,:)
  !send to node in +x
    call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
                              cart_neighbor(6),tag,cart_neighbor(5), &
                              tag,mpi_comm_cart,status,ierr)
    force_c(:,0,:,:)=force_c_buffer(:,:,:)

    force_c_buffer(:,:,:)=force_c(:,:,1,:)
  !send to node in -y  
    call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)
    force_c(:,:,nc_node_dim+1,:)=force_c_buffer(:,:,:)

    force_c_buffer(:,:,:)=force_c(:,:,nc_node_dim,:)
  !send to node in +y 
    call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
                              cart_neighbor(4),tag,cart_neighbor(3), &
                              tag,mpi_comm_cart,status,ierr)
    force_c(:,:,0,:)=force_c_buffer(:,:,:)

    force_c_buffer(:,:,:)=force_c(:,:,:,1)
  !send to node in -z 
    call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)
    force_c(:,:,:,nc_node_dim+1)=force_c_buffer(:,:,:)

    force_c_buffer(:,:,:)=force_c(:,:,:,nc_node_dim)
  !send to node in +z 
    call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &
                              tag,mpi_comm_cart,status,ierr)
    force_c(:,:,:,0)=force_c_buffer(:,:,:)

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('cf  buff',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'coarse force buffer finished',real(count_f-count_i)/real(count_r)
#endif

  end subroutine coarse_force_buffer
