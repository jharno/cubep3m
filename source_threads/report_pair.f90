!! this routine reports pair-wise force accuracy
  subroutine report_pair
    implicit none

    include 'mpif.h'
#    include "cubepm.fh"

    real(4), dimension(6,2) :: pair 
    integer(4), dimension(mpi_status_size) :: status
    integer(4) :: i,j,tag
    real(4) :: r(3),acc(3),v(3),F(3),magr,magF,del_F(2),frac_err(2)
    real(4) :: F_sim(3,2),magF_sim(2),magF_sim_r(2),magF_sim_t(2)
    integer(4) np_local_store
    real(4) :: xv_store(6,2)

    tag=42

    np_local_store=np_local
    xv_store(:,:)=xv(:,:2)

    if (rank==0) then
      if (np_local==2) then
        pair=xv(:,:2)
      elseif (np_local==1) then
        pair(:,1)=xv(:,1)
        call mpi_recv(pair(:,2),6,mpi_real,mpi_any_source, &
                      tag,mpi_comm_world,status,ierr)
      else
        call mpi_recv(pair(:,1),6,mpi_real,mpi_any_source, &
                      tag,mpi_comm_world,status,ierr)
        call mpi_recv(pair(:,2),6,mpi_real,mpi_any_source, &
                      tag,mpi_comm_world,status,ierr)
      endif
#ifdef DIAG
      write(*,'(6f10.5)') pair(:,1)
      write(*,'(6f10.5)') pair(:,2)
#endif
      r=pair(:3,1)-pair(:3,2)
      do j=1,3
        if (r(j) < -real(nf_physical_dim/2,kind=4)) &
            r(j) = r(j)+real(nf_physical_dim,kind=4) 
        if (r(j) > real(nf_physical_dim/2,kind=4)) &
            r(j) = r(j)-real(nf_physical_dim,kind=4)
      enddo
      dt_max_v=0.01/max(sqrt(pair(4,1)**2+pair(5,1)**2+pair(6,1)**2), &
                sqrt(pair(4,2)**2+pair(5,2)**2+pair(6,2)**2))
      magr=sqrt(r(1)**2+r(2)**2+r(3)**2)
      print *,'current seperation:',magr
      if (pair_infall) cur_sep=magr
      F=-G*r/magr**3.0
      magF=sqrt(F(1)**2+F(2)**2+F(3)**2)
      acc=F*mass_p
      v=acc*dt
      F_sim=pair(4:,:)/dt/mass_p
!! NEED velocity after fine mesh update for F_sim_fine
      do j=1,2
        magF_sim(j)=sqrt(F_sim(1,j)**2+F_sim(2,j)**2+F_sim(3,j)**2)
        magF_sim_r(j)=(F_sim(1,j)*F(1)*(-1.0)**(j) + &
                     F_sim(2,j)*F(2)*(-1.0)**(j) + &
                     F_sim(3,j)*F(3)*(-1.0)**(j))/magF
        magF_sim_t(j)=sqrt(magF_sim(j)**2-magF_sim_r(j)**2)
        del_F(j)=magF-magF_sim(j)
        frac_err(j)=del_F(j)/magF
      enddo
      open(41,file='pair_xv.dat',position='append')
      write(41,'(12f16.8)') pair 
      close(41)
      open(42,file='pair_F.dat',position='append')
      do j=1,2
        write(42,'(6f16.8)') magr, magF_sim(j), magF, frac_err(j), &
                             magF_sim_r(j),magF_sim_t(j)
      enddo
      close(42)
    else
      if (np_local==1) then
        xv(1,1)=xv(1,1)+cart_coords(3)*nf_physical_node_dim
        xv(2,1)=xv(2,1)+cart_coords(2)*nf_physical_node_dim
        xv(3,1)=xv(3,1)+cart_coords(1)*nf_physical_node_dim
        call mpi_send(xv(:,1),6,mpi_real,0,tag,mpi_comm_world,ierr)
      elseif (np_local==2) then
        do i=1,np_local
          xv(1,i)=xv(1,i)+cart_coords(3)*nf_physical_node_dim
          xv(2,i)=xv(2,i)+cart_coords(2)*nf_physical_node_dim
          xv(3,i)=xv(3,i)+cart_coords(1)*nf_physical_node_dim
        enddo
        call mpi_send(xv(:,1),6,mpi_real,0,tag,mpi_comm_world,ierr)
        call mpi_send(xv(:,2),6,mpi_real,0,tag,mpi_comm_world,ierr)
      endif
    endif  

    if (pair_infall) then
      np_local=np_local_store
      xv(:,:2)=xv_store(:,:)
    else
      np_local=0
    endif  

  end subroutine report_pair 
