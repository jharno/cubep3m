!! construct fine mesh force kernel
  subroutine fine_kernel

    implicit none 

    include 'mpif.h'
#    include "cubepm.fh"
  
    integer(kind=4) :: i,j,k,temp(3)
    integer(kind=4) :: errcode,fstat
    real(kind=4) :: rtemp(2)
   
!    open(unit=18,file=kernel_path//'wfxyzf.ascii',status='old',iostat=fstat)
!    open(unit=18,file=kernel_path//'wfxyzf.1.ascii',status='old',iostat=fstat)
!    open(unit=18,file=kernel_path//'wfxyzf.2.ascii',status='old',iostat=fstat)
    open(unit=18,file=kernel_path//'wfxyzf.3.ascii',status='old',iostat=fstat)
    if (fstat /= 0) then
      write(*,*) 'rank:',rank,'error opening fine mesh kernel'
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

! first do x-direction

    rho_f=0.0

    do k=1,nf_cutoff
      do j=1,nf_cutoff
        do i=1,nf_cutoff
          read(18,'(3i4,3e16.8)') temp(1),temp(2),temp(3),rho_f(i,j,k,1),rtemp(1),rtemp(2)
          if (temp(1).ne.i.or.temp(2).ne.j.or.temp(3).ne.k) then
            write(*,*) 'error reading in fine mesh kernel'
            errcode=-1
            call mpi_abort(mpi_comm_world,errcode,ierr)
          endif
        end do
      end do
    end do
           
#ifdef PP_EXT

    if (pp_ext_force_flag) then
       if (rank == 0) write(*,*) 'Using fine mesh kernel modified for pp extended to', pp_range,'neighbours'
       !! Set the kernel to 0 for immediate neighbors
       do k=1,pp_range+1
          do j=1,pp_range+1
             do i=1,pp_range+1
                rho_f(i,j,k,1) = 0
             enddo
          enddo
       enddo
    else
       !write(*,*) 'NOT using modified fine mesh kernel (flag set to .false.)'
    endif

#endif
#ifdef FPPKCORR
    write(*,*) 'Using fine mesh kernel modified to match pp/fine mesh force'
    !! boost the fine grid neighbours for improved force plots, when averaged over multiple time steps
    do k=1,4
       do j=1,4
          do i=1,4
             if(sqrt(real(i-1)**2 + real(j-1)**2 + real(k-1)**2) .le. 1.0 )then
               rho_f(i,j,k,1) = rho_f(i,j,k,1)*1.6
             endif
          enddo
       enddo
    enddo
#endif

!! Reflect accross y/2 plane
             
      do j=2,nf_cutoff
        rho_f(1:nf_cutoff,nf_tile-j+2,1:nf_cutoff,1)=rho_f(1:nf_cutoff,j,1:nf_cutoff,1)
      enddo  
               
!! Reflect accross x/2 plane
               
      do i=2,nf_cutoff
        rho_f(nf_tile-i+2,1:nf_tile,1:nf_cutoff,1)=-rho_f(i,1:nf_tile,1:nf_cutoff,1)
      enddo  

!! Reflect accross z/2 plane

      do k=2,nf_cutoff
        rho_f(1:nf_tile,1:nf_tile,nf_tile-k+2,1)=rho_f(1:nf_tile,1:nf_tile,k,1)
      enddo

!! Transform to Fourier space

      call cubepm_fftw2('f',1)

!! Select imaginary component for force kernel

      do k=1,nf_tile
        do j=1,nf_tile
          do i=1,nf_tile/2+1
            kern_f(1,i,j,k)=rho_f(2*i,j,k,1)
          enddo
        enddo
      enddo

! now do the y-direction

      rho_f=0.0

    rewind(18)
    do k=1,nf_cutoff
      do j=1,nf_cutoff
        do i=1,nf_cutoff
          read(18,*) temp(1),temp(2),temp(3),rtemp(1),rho_f(i,j,k,1),rtemp(2)
          if (temp(1).ne.i.or.temp(2).ne.j.or.temp(3).ne.k) then
            write(*,*) 'error reading in fine mesh kernel'
            errcode=-1
            call mpi_abort(mpi_comm_world,errcode,ierr)
          endif
        end do
      end do
    end do

#ifdef PP_EXT

    if (pp_ext_force_flag) then
       !write(*,*) 'Using modified fine mesh kernel'
       !! Set the kernel to 0 for immediate neighbors
       do k=1,pp_range+1
          do j=1,pp_range+1
             do i=1,pp_range+1
                rho_f(i,j,k,1) = 0
             enddo
          enddo
       enddo
    else
       !write(*,*) 'NOT using modified fine mesh kernel (flag set to .false.)'
    endif

#endif
#ifdef FPPKCORR
    !write(*,*) 'Using fine mesh kernel modified to match pp/fine mesh force'
    !! boost the fine grid neighbours for improved force plots, when averaged over multiple time steps
    do k=1,4
       do j=1,4
          do i=1,4
             if(sqrt(real(i-1)**2 + real(j-1)**2 + real(k-1)**2) .le. 1.0 )then
                rho_f(i,j,k,1) = rho_f(i,j,k,1)*1.6
             endif
          enddo
       enddo
    enddo
#endif


!! Reflect accross y/2 plane

      do j=2,nf_cutoff
        rho_f(1:nf_cutoff,nf_tile-j+2,1:nf_cutoff,1)=-rho_f(1:nf_cutoff,j,1:nf_cutoff,1)
      enddo

!! Reflect accross x/2 plane

      do i=2,nf_cutoff
        rho_f(nf_tile-i+2,1:nf_tile,1:nf_cutoff,1)=rho_f(i,1:nf_tile,1:nf_cutoff,1)
      enddo

!! Reflect accross z/2 plane

      do k=2,nf_cutoff
        rho_f(1:nf_tile,1:nf_tile,nf_tile-k+2,1)=rho_f(1:nf_tile,1:nf_tile,k,1)
      enddo

!! Transform to Fourier space
        
      call cubepm_fftw2('f',1)
            
!! Select imaginary component for force kernel
        
      do k=1,nf_tile
        do j=1,nf_tile
          do i=1,nf_tile/2+1
            kern_f(2,i,j,k)=rho_f(2*i,j,k,1)
          enddo
        enddo
      enddo

! now do the z-direction

      rho_f=0.0

    rewind(18)
    do k=1,nf_cutoff
      do j=1,nf_cutoff
        do i=1,nf_cutoff
          read(18,*) temp(1),temp(2),temp(3),rtemp(1),rtemp(2),rho_f(i,j,k,1)
          if (temp(1).ne.i.or.temp(2).ne.j.or.temp(3).ne.k) then
            write(*,*) 'error reading in fine mesh kernel'
            errcode=-1
            call mpi_abort(mpi_comm_world,errcode,ierr)
          endif
        end do
      end do
    end do
      
#ifdef PP_EXT

    if (pp_ext_force_flag) then
       !write(*,*) 'Using modified fine mesh kernel'
       !! Set the kernel to 0 for immediate neighbors
       do k=1,pp_range+1
          do j=1,pp_range+1
             do i=1,pp_range+1
                rho_f(i,j,k,1) = 0
             enddo
          enddo
       enddo
    else
       !write(*,*) 'NOT using modified fine mesh kernel (flag set to .false.)'
    endif

#endif
#ifdef FPPKCORR
    !write(*,*) 'Using fine mesh kernel modified to match pp/fine mesh force'
    !! boost the fine grid neighbours for improved force plots, when averaged over multiple time steps
    do k=1,4
       do j=1,4
          do i=1,4
             if(sqrt(real(i-1)**2 + real(j-1)**2 + real(k-1)**2) .le. 1.0 )then
                rho_f(i,j,k,1) = rho_f(i,j,k,1)*1.6
             endif
          enddo
       enddo
    enddo
#endif


!! Reflect accross y/2 plane
  
      do j=2,nf_cutoff
        rho_f(1:nf_cutoff,nf_tile-j+2,1:nf_cutoff,1)=rho_f(1:nf_cutoff,j,1:nf_cutoff,1)
      enddo 

!! Reflect accross x/2 plane

      do i=2,nf_cutoff
        rho_f(nf_tile-i+2,1:nf_tile,1:nf_cutoff,1)=rho_f(i,1:nf_tile,1:nf_cutoff,1)
      enddo

!! Reflect accross z/2 plane

      do k=2,nf_cutoff
        rho_f(1:nf_tile,1:nf_tile,nf_tile-k+2,1)=-rho_f(1:nf_tile,1:nf_tile,k,1)
      enddo

!! Transform to Fourier space

      call cubepm_fftw2('f',1)

!! Select imaginary component for force kernel

      do k=1,nf_tile
        do j=1,nf_tile
          do i=1,nf_tile/2+1
            kern_f(3,i,j,k)=rho_f(2*i,j,k,1)
          enddo
        enddo
      enddo

    close(18)

    end subroutine fine_kernel

!-----------------------------------------------------------------------------!

!! coarse grid force kernel construction
    subroutine coarse_kernel

      implicit none

      include 'mpif.h'
#      include "cubepm.fh"

      integer(kind=4) :: i,j,k,i0,j0,k0,hc1,temp(3)
      integer(4) :: il,ih,jl,jh,kl,kh,kx,ky,kz,fstat
      real(kind=4) :: x,y,z,r,kr,wa,wb,wc,ka,kb,kc
      character(len=80) :: kfile,kfilepath
!      real(kind=4) :: ck_table(3,6,6,6)
      real(kind=4) :: ck_table(3,4,4,4)
      integer(4) :: k1,j1,i1

      integer :: mg, dxy, dx
      integer :: ind, ig, jg, kg

      if (rank==0) print *,'inside coarse kernel init'

      hc1=nc_dim/2+1
      ck=0.0

!! offsets for local volume

      il=1+nc_node_dim*cart_coords(3)
      jl=1+nc_node_dim*cart_coords(2)
      kl=1+nc_node_dim*cart_coords(1)
      ih=nc_node_dim*(cart_coords(3)+1)
      jh=nc_node_dim*(cart_coords(2)+1)
      kh=nc_node_dim*(cart_coords(1)+1)

!! construct uncorrected force kernel

      do k=kl,kh
         k1=k-kl+1
         if (k .lt. nc_dim/2+2) then
             z=k-1
         else
             z=k-1-nc_dim
         endif
         z=mesh_scale*z
         do j=jl,jh
           j1=j-jl+1
           if (j .lt. nc_dim/2+2) then
             y=j-1
           else
             y=j-1-nc_dim
           endif
           y=mesh_scale*y
            do i=il,ih
             i1=i-il+1
             if (i .lt. nc_dim/2+2) then
               x=i-1
             else
               x=i-1-nc_dim
             endif
             x=mesh_scale*x
             r=sqrt(x*x+y*y+z*z)
             if (r.eq.0.0) then
               ck(:,i1,j1,k1)=0.0
             else
               ck(1,i1,j1,k1)=-x/r**3
               ck(2,i1,j1,k1)=-y/r**3
               ck(3,i1,j1,k1)=-z/r**3
             endif
           enddo
         enddo
       enddo

      if (rank==0) print *,'finished coarse kernel'

!! Read in kernel for 2-level matching

!        open(unit=11,file=kernel_path//'wfxyzc.ascii',status='old',iostat=fstat)
!        open(unit=11,file=kernel_path//'wfxyzc.1.ascii',status='old',iostat=fstat)
        open(unit=11,file=kernel_path//'wfxyzc.2.ascii',status='old',iostat=fstat)
        if (fstat /= 0) then
          write(*,*) 'rank:',rank,'error opening coarse mesh kernel'
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        do k=1,4
          do j=1,4
            do i=1,4
              read(11,'(3i4,3e16.8)') temp(:),ck_table(:,i,j,k)
#ifdef DEBUG_VEL
              if (rank==0) print *,temp(:),ck(:,i,j,k)
#endif
            enddo
          enddo
        enddo
        close(11)

      if (rank==0) print *,'finished kernel 2-level matching'

!! copy corrections to other octants of the kernel
!! assuming nc_node_dim > 8 

      if (nodes_dim==1) then
        ck(:,:4,:4,:4)=ck_table(:,:,:,:)
        do k=2,4
          ck(1:2,:4,:4,nc_node_dim-k+2)=ck_table(1:2,:,:,k)
          ck(3,:4,:4,nc_node_dim-k+2)=-ck_table(3,:,:,k)
        enddo
        do j=2,4
          ck(1,:4,nc_node_dim-j+2,:4)=ck_table(1,:,j,:)
          ck(2,:4,nc_node_dim-j+2,:4)=-ck_table(2,:,j,:)
          ck(3,:4,nc_node_dim-j+2,:4)=ck_table(3,:,j,:)
        enddo
        do k=2,4
          do j=2,4
            ck(1,:4,nc_node_dim-j+2,nc_node_dim-k+2)=ck_table(1,:,j,k)
            ck(2:3,:4,nc_node_dim-j+2,nc_node_dim-k+2)=-ck_table(2:3,:,j,k)
          enddo
        enddo
        do i=2,4
          ck(1,nc_node_dim-i+2,:4,:4)=-ck_table(1,i,:,:)
          ck(2:3,nc_node_dim-i+2,:4,:4)=ck_table(2:3,i,:,:)
        enddo
        do k=2,4
          do i=2,4
            ck(1,nc_node_dim-i+2,:4,nc_node_dim-k+2)=-ck_table(1,i,:,k)
            ck(2,nc_node_dim-i+2,:4,nc_node_dim-k+2)=ck_table(2,i,:,k)
            ck(3,nc_node_dim-i+2,:4,nc_node_dim-k+2)=-ck_table(3,i,:,k)
          enddo
        enddo
        do j=2,4
          do i=2,4
            ck(1:2,nc_node_dim-i+2,nc_node_dim-j+2,:4)=-ck_table(1:2,i,j,:)
            ck(3,nc_node_dim-i+2,nc_node_dim-j+2,:4)=ck_table(3,i,j,:)
          enddo
        enddo
        do k=2,4
          do j=2,4
            do i=2,4
              ck(1:3,nc_node_dim-i+2,nc_node_dim-j+2,nc_node_dim-k+2)=-ck_table(:,i,j,k)
            enddo
          enddo
        enddo
      else
        if (cart_coords(3)==0.and.cart_coords(2)==0.and.cart_coords(1)==0) then
          ck(:,:4,:4,:4)=ck_table(:,:,:,:)
        elseif( cart_coords(3)==0.and.cart_coords(2)==0.and.cart_coords(1)==nodes_dim-1) then
          do k=2,4
            ck(1:2,:4,:4,nc_node_dim-k+2)=ck_table(1:2,:,:,k)
            ck(3,:4,:4,nc_node_dim-k+2)=-ck_table(3,:,:,k)
          enddo
        elseif( cart_coords(3)==0.and.cart_coords(2)==nodes_dim-1.and.cart_coords(1)==0) then
          do j=2,4
            ck(1,:4,nc_node_dim-j+2,:4)=ck_table(1,:,j,:)
            ck(2,:4,nc_node_dim-j+2,:4)=-ck_table(2,:,j,:)
            ck(3,:4,nc_node_dim-j+2,:4)=ck_table(3,:,j,:)
          enddo
        elseif( cart_coords(3)==0.and.cart_coords(2)==nodes_dim-1.and.cart_coords(1)==nodes_dim-1) then
          do k=2,4
            do j=2,4
              ck(1,:4,nc_node_dim-j+2,nc_node_dim-k+2)=ck_table(1,:,j,k)
              ck(2:3,:4,nc_node_dim-j+2,nc_node_dim-k+2)=-ck_table(2:3,:,j,k)
            enddo
          enddo
        elseif( cart_coords(3)==nodes_dim-1.and.cart_coords(2)==0.and.cart_coords(1)==0) then
          do i=2,4
            ck(1,nc_node_dim-i+2,:4,:4)=-ck_table(1,i,:,:)
            ck(2:3,nc_node_dim-i+2,:4,:4)=ck_table(2:3,i,:,:)
          enddo
        elseif( cart_coords(3)==nodes_dim-1.and.cart_coords(2)==0.and.cart_coords(1)==nodes_dim-1) then
          do k=2,4
            do i=2,4
              ck(1,nc_node_dim-i+2,:4,nc_node_dim-k+2)=-ck_table(1,i,:,k)
              ck(2,nc_node_dim-i+2,:4,nc_node_dim-k+2)=ck_table(2,i,:,k)
              ck(3,nc_node_dim-i+2,:4,nc_node_dim-k+2)=-ck_table(3,i,:,k)
            enddo
          enddo
        elseif( cart_coords(3)==nodes_dim-1.and.cart_coords(2)==nodes_dim-1.and.cart_coords(1)==0) then
          do j=2,4
            do i=2,4
              ck(1:2,nc_node_dim-i+2,nc_node_dim-j+2,:4)=-ck_table(1:2,i,j,:)
              ck(3,nc_node_dim-i+2,nc_node_dim-j+2,:4)=ck_table(3,i,j,:)
            enddo
          enddo
        elseif( cart_coords(3)==nodes_dim-1.and.cart_coords(2)==nodes_dim-1.and.cart_coords(1)==nodes_dim-1) then
          do k=2,4
            do j=2,4
              do i=2,4
                ck(1:3,nc_node_dim-i+2,nc_node_dim-j+2,nc_node_dim-k+2)=-ck_table(:,i,j,k)
              enddo
            enddo
          enddo
        endif
      endif

      if (rank==0) print *,'finished octants'

!! transform and save imaginary result in kern_c

!! For the analytically matched kernel we have to apply an additional correction

#ifdef LRCKCORR

!! Copy corrected kernel to force_c for later use

      do k=1,nc_node_dim
        do j=1,nc_node_dim
          do i=1,nc_node_dim
            force_c(:,i,j,k)=ck(:,i,j,k)
          enddo
        enddo
      enddo

!! construct uncorrected force kernel in ck

      do k=kl,kh
         k1=k-kl+1
         if (k .lt. nc_dim/2+2) then
             z=k-1
         else
             z=k-1-nc_dim
         endif
         z=mesh_scale*z
         do j=jl,jh
           j1=j-jl+1
           if (j .lt. nc_dim/2+2) then
             y=j-1
           else
             y=j-1-nc_dim
           endif
           y=mesh_scale*y
            do i=il,ih
             i1=i-il+1
             if (i .lt. nc_dim/2+2) then
               x=i-1
             else
               x=i-1-nc_dim
             endif
             x=mesh_scale*x
             r=sqrt(x*x+y*y+z*z)
             if (r.eq.0.0) then
               ck(:,i1,j1,k1)=0.0
             else
               ck(1,i1,j1,k1)=-x/r**3
               ck(2,i1,j1,k1)=-y/r**3
               ck(3,i1,j1,k1)=-z/r**3
             endif
           enddo
         enddo
       enddo

      if (rank==0) print *,'finished uncorrected force kernel'

!! transform un-corrected kernel, store in tmp_kern_c

      rho_c=ck(1,:,:,:)
      call cubepm_fftw(1)
      do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
          do i=1,nc_dim !+2
            tmp_kern_c(1,i,j,k)=slab(i,j,k)
          enddo
        enddo
      enddo
      
      if (rank==0) print *,'finished transforming uncorrected kernel'

      rho_c=ck(2,:,:,:)
      call cubepm_fftw(1)
      do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
          do i=1,nc_dim !+2
            tmp_kern_c(2,i,j,k)=slab(i,j,k)
          enddo
        enddo
      enddo
      
      if (rank==0) print *,'finished second transform'

      rho_c=ck(3,:,:,:)
      call cubepm_fftw(1)
      do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
          do i=1,nc_dim !+2
            tmp_kern_c(3,i,j,k)=slab(i,j,k)
          enddo
        enddo
      enddo

      if (rank==0) print *,'finished third transform'

!! transform corrected kernel and apply long range correction

! x

      rho_c=force_c(1,1:nc_node_dim,1:nc_node_dim,1:nc_node_dim) 
      call cubepm_fftw(1)

! -------------------------------------------------------------------
! CHECK DETAILS OF THESE LOOPS LATER TO MAKE SURE CORRECT FOR P3DFFT
! -------------------------------------------------------------------

#ifdef DEBUG_P3D_MPI
    if (rank == 0) write(*,*) "starting 1st kc loop"
#endif

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
                if (kg < nc_dim/2+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc_dim
                endif
                if (jg < nc_dim/2+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc_dim
                endif
                kx = (ig-1)/2
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr <= 8.) then
                    ka=2*sin(pi*kx/real(nc_dim))
                    kb=2*sin(pi*ky/real(nc_dim))
                    kc=2*sin(pi*kz/real(nc_dim))
                    if (kx /= 0) then
                        wa=slab(i+1,j,k)
                        wb=tmp_kern_c(1,i+1,j,k)
                        wc=4.*pi*ka/(ka**2+kb**2+kc**2)/16.
                        slab(i+1,j,k)=wa*(wc/wb)
                    endif
                endif
            enddo
        enddo
    enddo

    do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
            do i=1,nc_dim/2 !+1
                kern_c(1,i,j,k)=slab(2*i,j,k)
            enddo
        enddo
    enddo

! y

    rho_c=force_c(2,1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)
    call cubepm_fftw(1)

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
                if (kg < nc_dim/2+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc_dim
                endif
                if (jg < nc_dim/2+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc_dim
                endif
                kx = (ig-1)/2
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr <= 8.) then
                    ka=2*sin(pi*kx/real(nc_dim))
                    kb=2*sin(pi*ky/real(nc_dim))
                    kc=2*sin(pi*kz/real(nc_dim))
                    if (ky /= 0) then
                        wa=slab(i+1,j,k)
                        wb=tmp_kern_c(2,i+1,j,k)
                        wc=4.*pi*kb/(ka**2+kb**2+kc**2)/16.
                        slab(i+1,j,k)=wa*(wc/wb)
                    endif
                endif
            enddo
        enddo
    enddo

    do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
            do i=1,nc_dim/2
                kern_c(2,i,j,k)=slab(2*i,j,k)
            enddo
        enddo
    enddo

! z 

      rho_c=force_c(3,1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)
      call cubepm_fftw(1)

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
                if (kg < nc_dim/2+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc_dim
                endif
                if (jg < nc_dim/2+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc_dim
                endif
                kx = (ig-1)/2
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr <= 8.) then
                    ka=2*sin(pi*kx/real(nc_dim))
                    kb=2*sin(pi*ky/real(nc_dim))
                    kc=2*sin(pi*kz/real(nc_dim))
                    if (kz /= 0) then
                        wa=slab(i+1,j,k)
                        wb=tmp_kern_c(3,i+1,j,k)
                        wc=4.*pi*kc/(ka**2+kb**2+kc**2)/16.
                        slab(i+1,j,k)=wa*(wc/wb)
                    endif
                endif
            enddo
        enddo
    enddo

    do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
            do i=1,nc_dim/2
                kern_c(3,i,j,k)=slab(2*i,j,k)
            enddo
        enddo
    enddo

    if (rank==0) print *,'finished x y z'

! This is the non-LRCKCORR case:

#else

      rho_c=ck(1,:,:,:)
      call cubepm_fftw(1)
      do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
          do i=1,nc_dim/2
            kern_c(1,i,j,k)=slab(2*i,j,k)
          enddo
        enddo
      enddo

      rho_c=ck(2,:,:,:)
      call cubepm_fftw(1)
      do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
          do i=1,nc_dim/2
            kern_c(2,i,j,k)=slab(2*i,j,k)
          enddo
        enddo
      enddo

      rho_c=ck(3,:,:,:)
      call cubepm_fftw(1)
      do k=1,nc_pen+mypadd !nc_slab
        do j=1,nc_node_dim !nc_dim
          do i=1,nc_dim/2
            kern_c(3,i,j,k)=slab(2*i,j,k)
          enddo
        enddo
      enddo
#endif

#ifdef DEBUG_VEL
      print *,'coarse kernel',rank,'min',minval(kern_c),'max',maxval(kern_c)
#endif

      if (rank==0) print *,'finished coarse kernel' 

  end subroutine coarse_kernel
