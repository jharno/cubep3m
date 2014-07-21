!! calculate mass power spectrum using coarse mesh density
  subroutine coarse_power
    use omp_lib
#ifdef FFTMKL 
    use MKL_DFTI
#endif
    implicit none
  
    include 'mpif.h'
#    include "cubepm.fh"

    integer(4), parameter :: hc=nc_dim/2
    integer(4) :: i,j,k,kg,jg,ig
    integer(4) :: k1, k2
    real(4) :: kz,ky,kx,kr,w1,w2,pow,x,y,z,sync_x, sync_y, sync_z,kernel
    character(len=7) :: z_s
    integer(4) :: fstat
    character(len=max_path) :: ofile
    real(4) :: rho_c_mean, sum_od 
    integer(4) :: ind
    integer(4) :: mg, dx, dxy

    rho_c_mean=(real(nf_physical_dim/2))**3*mass_p/(real(nc_dim))**3

    !print *,'rank',rank,'rho_c_mean',rho_c_mean

!! First we have to change cmplx_rho_c to overdensity and re-transform

    slab=cmplx_rho_c
    call cubepm_fftw(-1)
    rho_c=rho_c/rho_c_mean-1.0
    !write(*,*) 'rank',rank,'sum(overdensity)=',sum(rho_c)
    call mpi_reduce(sum(rho_c),sum_od,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr) 
    if (rank==0) write(*,*) 'total(overdensity)=',sum_od
    call cubepm_fftw(1)

    ps_c=0.0

    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0

    do k = 1, nc_pen+mypadd
        do j = 1, nc_node_dim
            do i = 1, nc_dim, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc_dim
                endif
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc_dim
                endif
                kx = (ig-1)/2                
                kr=sqrt(kx**2+ky**2+kz**2)
                if(kx.eq.0 .and. ky <=0 .and. kz <=0)cycle;
                if(kx.eq.0 .and. ky >0 .and. kz <0)cycle;
                if (kr /= 0.0) then
                    k1=ceiling(kr)
                    k2=k1+1
                    w1=k1-kr
                    w2=1-w1
                    x = pi*real(kx)/nc_dim
                    y = pi*real(ky)/nc_dim
                    z = pi*real(kz)/nc_dim
                        
                    if(x==0) then 
                       sync_x = 1
                    else
                       sync_x = sin(x)/x
                    endif
                    if(y==0) then 
                       sync_y = 1
                    else
                       sync_y = sin(y)/y
                    endif
                    if(z==0) then 
                       sync_z = 1
                    else
                       sync_z = sin(z)/z
                    endif

                    kernel = sync_x*sync_y*sync_z
!#ifdef NGPPS
                    w1=1
                    w2=0
!#endif                
                    pow=(slab(i,j,k)/(real(nc_dim))**3)**2+(slab(i+1,j,k)/(real(nc_dim))**3)**2/kernel**4
                    ps_c(1,k1)=ps_c(1,k1)+w1
                    ps_c(2,k1)=ps_c(2,k1)+w1*pow
                    ps_c(1,k2)=ps_c(1,k2)+w2
                    ps_c(2,k2)=ps_c(2,k2)+w2*pow
                endif
            enddo
        enddo
    enddo
   
!! Reduce ps on master node

    ps_c_sum=0.0
    call mpi_reduce(ps_c,ps_c_sum,nc_dim*2,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
 
!! Divide by weights / convert P(k) to \delta^2(k)

    if (rank == 0) then
      do k=1,nc_dim
        if (ps_c_sum(1,k) /= 0) then
          ps_c_sum(2,k)=4.0*pi*((real(k)-1)**3)*ps_c_sum(2,k)/ps_c_sum(1,k)
          ps_c_sum(1,k)=2.0*pi*(real(k)-1)/box
        endif
      enddo

!! write power spectrum to disk 

      write(z_s,'(f7.3)') 1/a_mid - 1.0
      z_s=adjustl(z_s)
      ofile=output_path//z_s(1:len_trim(z_s))//'ps.dat'
      open(50,file=ofile,status='replace',iostat=fstat,form='formatted')

      if (fstat /= 0) then
        write(*,*) 'error:',fstat,'opening power spectrum file for write'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      do k=1,nc_dim
        write(50,'(2f20.10)') ps_c_sum(:,k)
      enddo

      close(50)

    endif

  end subroutine coarse_power     
