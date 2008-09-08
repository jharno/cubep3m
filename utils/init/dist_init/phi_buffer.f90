!! phi_buffer -- buffer cubic decomposition mesh

subroutine phi_buffer(bt)
  implicit none

  include 'mpif.h'
  include 'dist_init.fh'

  integer(4) :: bt
  integer(4) :: buffer_size
  integer(4) :: tag
  integer(4) :: status(MPI_STATUS_SIZE)

  buffer_size = (nc_node_dim + 2)**2

  tag=64

!! send to node in -x
    if (bt==1) then
      phi_buf(:,:)=phi(0,:,:)
    else
      phi_buf(:,:)=phi(1,:,:)
    endif
    call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(5),tag,cart_neighbor(6), &
                              tag,mpi_comm_cart,status,ierr)
    if (bt==1) then
      phi(nc_node_dim,:,:)=phi(nc_node_dim,:,:)+phi_buf(:,:)
    else
      phi(nc_node_dim+1,:,:)=phi_buf(:,:)
    endif

!! send to node in +x
    if (bt==1) then
      phi_buf(:,:)=phi(nc_node_dim+1,:,:)
    else
      phi_buf(:,:)=phi(nc_node_dim,:,:)
    endif
    call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(6),tag,cart_neighbor(5), &
                              tag,mpi_comm_cart,status,ierr)
    if (bt==1) then
      phi(1,:,:)=phi(1,:,:)+phi_buf(:,:)
    else
      phi(0,:,:)=phi_buf(:,:)
    endif


!! send to node in -y
    if (bt==1) then
      phi_buf(:,:)=phi(:,0,:)
    else
      phi_buf(:,:)=phi(:,1,:)
    endif
    call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)
    if (bt==1) then
      phi(:,nc_node_dim,:)=phi(:,nc_node_dim,:)+phi_buf(:,:)
    else
      phi(:,nc_node_dim+1,:)=phi_buf(:,:)
    endif

!! send to node in +y
    if (bt==1) then
      phi_buf(:,:)=phi(:,nc_node_dim+1,:)
    else
      phi_buf(:,:)=phi(:,nc_node_dim,:)
    endif
    call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(4),tag,cart_neighbor(3), &
                              tag,mpi_comm_cart,status,ierr)
    if (bt==1) then
      phi(:,1,:)=phi(:,1,:)+phi_buf(:,:)
    else
      phi(:,0,:)=phi_buf(:,:)
    endif

!! send to node in -z
    if (bt==1) then
      phi_buf(:,:)=phi(:,:,0)
    else
      phi_buf(:,:)=phi(:,:,1)
    endif
    call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)
    if (bt==1) then
      phi(:,:,nc_node_dim)=phi(:,:,nc_node_dim)+phi_buf(:,:)
    else
      phi(:,:,nc_node_dim+1)=phi_buf(:,:)
    endif

!! send to node in +z
    if (bt==1) then
      phi_buf(:,:)=phi(:,:,nc_node_dim+1)
    else
      phi_buf(:,:)=phi(:,:,nc_node_dim)
    endif
    call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &
                              tag,mpi_comm_cart,status,ierr)
    if (bt==1) then
      phi(:,:,1)=phi(:,:,1)+phi_buf(:,:)
    else
      phi(:,:,0)=phi_buf(:,:)
    endif

  end subroutine phi_buffer 
