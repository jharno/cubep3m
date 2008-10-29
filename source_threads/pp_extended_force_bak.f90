!! Computes pp force between particles of neighboring fine cells
 subroutine pp_extended_force(tile,thread)
   implicit none
   include 'cubep5m.fh'

   integer(4)               :: im, jm, km, i, j, k, ib, jb, kb
   integer(4), dimension(3) :: tile, i1, i2, cic_l, cic_h
   real(4),    dimension(3) :: x, offset
   
   !particle loop
   integer(4)               :: pp, ip1, ip2, pp1, pp2, thread
   
   
   !Extended pp particle count 
   integer ipl_ext(mesh_scale+buff_fm,mesh_scale+buff_fm,mesh_scale+buff_fm)

   !Force calculation
   real sep(3), force_pp(3), rmag, pp_force_mag
   real(4), dimension(3,max_llf,cores) :: pp1_force_ext_accum, pp2_force_ext_accum


   
   
   ! Loop only over the physical volume plus a layer of buffer in the '-' direction
   ! The '+' direction is taken care of in the coming loop
  
   cic_l(:) = nc_tile_dim * tile(:) + 1
   cic_h(:) = nc_tile_dim * (tile(:) + 1)
   
   do i = 1,3
      if(tile(i) == 0) cic_l(i) = 0
   enddo

#ifdef DEBUG
   write(*,*) 'Called pp_extended_force'
   write(*,*) 'Tile :', tile, 'thread :', thread
   write(*,*) 'Calculatied BC:', cic_l, cic_h    
#endif

   !write(*,*) 'llf_ext(1,1,1,1,thread) = ',llf_ext(1,1,1,1,thread)


   do k = cic_l(3), cic_h(3)
      do j = cic_l(2), cic_h(2)
         do i = cic_l(1), cic_h(1)
            
            
            

            !Get the core of the cube from current coarse cell       
            
            pp=hoc(i,j,k)
            
            ipl_ext = 0
            !write(*,*) 'llf_ext(1,1,1,1,thread) = ',llf_ext(1,1,1,1,thread)
           
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)                             
               i1(:) = floor(x(:)) + 1
             
               ! Get fine mesh coordinate 
               ! Could have negative values in the global  buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
               
               !Get particle count
               ipl_ext(i2(1),i2(2),i2(3))=ipl_ext(i2(1),i2(2),i2(3))+1
               if (ipl_ext(i2(1),i2(2),i2(3))>max_llf) then
                  print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                  stop
               endif
               

               ! Build a linked list for particles in that extended coarse cell
               llf_ext(ipl_ext(i2(1),i2(2),i2(3)),i2(1),i2(2),i2(3),thread)=pp
               !write(*,*) 'llf_ext', llf_ext(ipl_ext(i2(1),i2(2),i2(3)),i2(1),i2(2),i2(3),thread), pp 
             
               pp = ll(pp)
               
            enddo
#ifdef P5M_DEBUG            
            write(*,*) 'Core done on thread',thread, 'with particles', sum(ipl_ext) 
#endif
            !write(*,*) 'xy = 11, along z-axis', ipl_ext(1,1,1), ipl_ext(1,1,2), ipl_ext(1,1,3), ipl_ext(1,1,4), ipl_ext(1,1,5)
            !write(*,*) 'xy = 31, along z-axis', ipl_ext(3,1,1), ipl_ext(3,1,2), ipl_ext(3,1,3), ipl_ext(3,1,4), ipl_ext(3,1,5)
            !write(*,*) 'xy = 13, along z-axis', ipl_ext(1,3,1), ipl_ext(1,3,2), ipl_ext(1,3,3), ipl_ext(1,3,4), ipl_ext(1,3,5)
            !write(*,*) 'xy = 33, along z-axis', ipl_ext(3,3,1), ipl_ext(3,3,2), ipl_ext(3,3,3), ipl_ext(3,3,4), ipl_ext(3,3,5)



            !Get the buffer to complete the cube
            

            !**************
            !x
            !**************

            pp=hoc(i+1,j,k)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               !Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
            
               ! Verify that particle in x-buffer 
               if(i2(1) <= buff_fm)  then

                  !write(*,*) 'particle in x-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1)+mesh_scale,i2(2),i2(3))=ipl_ext(i2(1)+mesh_scale,i2(2),i2(3))+1
                  if (ipl_ext(i2(1)+mesh_scale,i2(2),i2(3))>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1)+mesh_scale,i2(2),i2(3)),i2(1)+mesh_scale,i2(2),i2(3),thread)=pp
                  !write(*,*) 'llf_ext', llf_ext(ipl_ext(i2(1)+mesh_scale,i2(2),i2(3)),i2(1)+mesh_scale,i2(2),i2(3),thread), pp 
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif

               pp = ll(pp)
               
            enddo
            
#ifdef P5M_DEBUG            
            write(*,*) 'x-buffer done on thread ', thread,'with particles', sum(ipl_ext) 
#endif
            !!***should be 'core' + 'x-buff'            

            !write(*,*) 'xy = 51, along z-axis', ipl_ext(5,1,1), ipl_ext(5,1,2), ipl_ext(5,1,3), ipl_ext(5,1,4), ipl_ext(5,1,5)
            !write(*,*) 'xy = 53, along z-axis', ipl_ext(5,3,1), ipl_ext(5,3,2), ipl_ext(5,3,3), ipl_ext(5,3,4), ipl_ext(5,3,5)


            !**************           
            !y
            !**************

            pp=hoc(i,j+1,k)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               ! Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
            
               ! Verify that particle in y-buffer 
               if(i2(2) <= buff_fm)  then

                  !write(*,*) 'particle in y-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1),i2(2)+mesh_scale,i2(3))=ipl_ext(i2(1),i2(2)+mesh_scale,i2(3))+1
                  if (ipl_ext(i2(1),i2(2)+mesh_scale,i2(3))>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1),i2(2)+mesh_scale,i2(3)),i2(1),i2(2)+mesh_scale,i2(3),thread)=pp
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif
                 
               pp = ll(pp)
               
            enddo
 
#ifdef P5M_DEBUG            
            write(*,*) 'y-buffer done on thread ', thread, 'with particles', sum(ipl_ext)
#endif            
            !!***should be 'core' + 'x-buff + y-buff'            

            !write(*,*) 'xy = 15, along z-axis', ipl_ext(1,5,1), ipl_ext(1,5,2), ipl_ext(1,5,3), ipl_ext(1,5,4), ipl_ext(1,5,5)
            !write(*,*) 'xy = 35, along z-axis', ipl_ext(3,5,1), ipl_ext(3,5,2), ipl_ext(3,5,3), ipl_ext(3,5,4), ipl_ext(3,5,5)
           

            !**************
            !z
            !**************

            pp=hoc(i,j,k+1)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               !Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                 endif
               enddo
            
               ! Verify that particle in z-buffer 
               if(i2(3) <= buff_fm)  then

                  !write(*,*) 'particle in z-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1),i2(2),i2(3)+mesh_scale)=ipl_ext(i2(1),i2(2),i2(3)+mesh_scale)+1
                  if (ipl_ext(i2(1),i2(2),i2(3)+mesh_scale)>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1),i2(2),i2(3)+mesh_scale),i2(1),i2(2),i2(3)+mesh_scale,thread)=pp
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif

               pp = ll(pp)
               
            enddo
            
#ifdef P5M_DEBUG            
            write(*,*) 'z-buffer done on thread ', thread, 'with particles', sum(ipl_ext)
#endif
            !!***should be 'core' + 'x-buff + y-buff + z-buff'            

            !write(*,*) 'xy = 11, along z-axis', ipl_ext(1,1,1), ipl_ext(1,1,2), ipl_ext(1,1,3), ipl_ext(1,1,4), ipl_ext(1,1,5), ipl_ext(1,1,6)
            !write(*,*) 'xy = 33, along z-axis', ipl_ext(3,3,1), ipl_ext(3,3,2), ipl_ext(3,3,3), ipl_ext(3,3,4), ipl_ext(3,3,5), ipl_ext(3,3,6)


            !**************
            !xy
            !**************

            pp=hoc(i+1,j+1,k)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               !Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
            
               ! Verify that particle in xy-buffer 
               if(i2(1) <= buff_fm .and. i2(2) <= buff_fm)  then

                  !write(*,*) 'particle in xy-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3))=ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3))+1
                  if (ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3))>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3)),i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3),thread)=pp
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif

               pp = ll(pp)
               
            enddo
 
            !write(*,*) 'xy = 55, along z-axis', ipl_ext(5,5,1), ipl_ext(5,5,2), ipl_ext(5,5,3), ipl_ext(5,5,4), ipl_ext(5,5,5), ipl_ext(5,5,6)

           
#ifdef P5M_DEBUG            
            write(*,*) 'xy-buffer done on thread ', thread, i,j,k
#endif
            
            !**************
            !xz
            !**************

            pp=hoc(i+1,j,k+1)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               !Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
            
               ! Verify that particle in xz-buffer 
               if(i2(1) <= buff_fm .and. i2(3) <= buff_fm)  then

                  !write(*,*) 'particle in xz-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1)+mesh_scale,i2(2),i2(3)+mesh_scale)=ipl_ext(i2(1)+mesh_scale,i2(2),i2(3)+mesh_scale)+1
                  if (ipl_ext(i2(1)+mesh_scale,i2(2),i2(3)+mesh_scale)>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1)+mesh_scale,i2(2),i2(3)+mesh_scale),i2(1)+mesh_scale,i2(2),i2(3)+mesh_scale,thread)=pp
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif
                 
               pp = ll(pp)
               
            enddo
 
            !write(*,*) 'xy = 51, along z-axis', ipl_ext(5,1,1), ipl_ext(5,1,2), ipl_ext(5,1,3), ipl_ext(5,1,4), ipl_ext(5,1,5), ipl_ext(5,1,6)

#ifdef P5M_DEBUG                       
            write(*,*) 'xz-buffer done on thread ', thread, i,j,k
#endif
            
            !**************
            !yz
            !**************

            pp=hoc(i,j+1,k+1)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               !Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
            
               ! Verify that particle in yz-buffer 
               if(i2(2) <= buff_fm .and. i2(3) <= buff_fm)  then

                  !write(*,*) 'particle in yz-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1),i2(2)+mesh_scale,i2(3)+mesh_scale)=ipl_ext(i2(1),i2(2)+mesh_scale,i2(3)+mesh_scale)+1
                  if (ipl_ext(i2(1),i2(2)+mesh_scale,i2(3)+mesh_scale)>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1),i2(2)+mesh_scale,i2(3)+mesh_scale),i2(1),i2(2)+mesh_scale,i2(3)+mesh_scale,thread)=pp
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif
                 
               pp = ll(pp)
               
            enddo

            !write(*,*) 'xy = 35, along z-axis', ipl_ext(3,5,1), ipl_ext(3,5,2), ipl_ext(3,5,3), ipl_ext(3,5,4), ipl_ext(3,5,5), ipl_ext(3,5,6)
            
#ifdef P5M_DEBUG            
            write(*,*) 'yz-buffer done on thread ', thread, i,j,k
#endif
            
            !**************
            !xyz
            !**************

            pp=hoc(i+1,j+1,k+1)
            
            do
               if (pp == 0) exit
               x(:) = xv(1:3,pp)! + offset(:)
               i1(:) = floor(x(:)) + 1
            
               !Get fine mesh coordinate 
               ! Could have negative values in the tile buffer: need to convert!
               do im=1,3
                  i2(im)=mod(i1(im)-1,mesh_scale)+1
                  if(i2(im)<1) then
                     !write(*,*) 'in the global buffer !' , i2(im), x(im)
                     i2(im) = i2(im) + mesh_scale
                  endif
               enddo
            
               ! Verify that particle in xyz-buffer 
               if(i2(1) <= buff_fm .and. i2(2) <= buff_fm .and. i2(3) <= buff_fm)  then

                  !write(*,*) 'particle in xyz-buffer' , i2(1), i2(2), i2(3)                
                  ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3)+mesh_scale)=ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3)+mesh_scale)+1
                  if (ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3)+mesh_scale)>max_llf) then
                     print *,'exceeded max_llf',max_llf,i1,i2,ipl_ext
                     stop
                  endif

                  !Add to the extended coarse cell linked list
                  llf_ext(ipl_ext(i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3)+mesh_scale),i2(1)+mesh_scale,i2(2)+mesh_scale,i2(3)+mesh_scale,thread)=pp
                  
               else
                  !write(*,*) 'particle not in buffer' , i2(1), i2(2), i2(3)
               endif
                 
               pp = ll(pp)
               
            enddo
 
            !write(*,*) 'xy = 55, along z-axis', ipl_ext(5,5,1), ipl_ext(5,5,2), ipl_ext(5,5,3), ipl_ext(5,5,4), ipl_ext(5,5,5), ipl_ext(5,5,6)
           
#ifdef P5M_DEBUG            
            write(*,*) 'xyx-buffer done on thread ', thread, i,j,k

            write(*,*) '***************************'
            write(*,*) 'Extended cube done on cell ',i,j,k, sum(ipl_ext)
            write(*,*) '***************************'
#endif
  
            !*************************************
            !*************************************
            !*************************************


            !write(*,*)'Starting fine cell loop for extended pp interaction  in coarse cell' ,i,j,k

            
            ! 1- loop over all particles in each fine cells in the non-extended coarse cell
            do km=1,mesh_scale
               do jm=1,mesh_scale
                  do im=1,mesh_scale
                                                             

                     !write(*,*) 'Stats before cycle:' ,im,jm,km, ipl_ext(im,jm,km)                     
                     
                     if(ipl_ext(im,jm,km) == 0) then
                        !write(*,*) 'No particles here!'
                        cycle
                     endif
                     
                     !write(*,*) 'Stats after cycle :' ,im,jm,km, ipl_ext(im,jm,km), llf_ext(ipl_ext(im,jm,km),im,jm,km,thread)                    

                     pp1_force_ext_accum(:,:ipl_ext(im,jm,km),thread)=0.0
                     pp2_force_ext_accum(:,:ipl_ext(im,jm,km),thread)=0.0
                      
                     do ip1=1,ipl_ext(im,jm,km)
                        pp1=llf_ext(ip1,im,jm,km,thread)
                        
                        !write(*,*) 'Inside the ip1 loop', ip1!, pp1


                        ! 2- loop over each neighboring cells 
                        do kb = 0,buff_fm
                           do jb = 0,buff_fm
                              do ib = 0,buff_fm

                                 if((ib+jb+kb) == 0) cycle ! could do PP interaction here...
                                               
                                 ! 3- loop over each particles there
                                 if(ipl_ext(im+ib,jm+jb,km+kb) >= 0) then
                                    do ip2=1,ipl_ext(im+ib,jm+jb,km+kb)
                                       pp2=llf_ext(ip2,im+ib,jm+jb,km+kb,thread)
                                       
#ifdef P5M_DEBUG
                                       write(*,*) 'close neighbors' , ip1, ip2, pp1, pp2 ,i,j,k,im,jm,km,xv(:3,pp1),xv(:3,pp2)
#endif
                                       sep=xv(:3,pp1)-xv(:3,pp2)
                                       rmag=sqrt(sep(1)*sep(1)+sep(2)*sep(2)+sep(3)*sep(3))
                                       if (rmag>rsoft) then
                                          force_pp=mass_p*(sep/(rmag*pp_bias)**3)  !mass_p divides out below
                                          pp1_force_ext_accum(:,ip1,thread)=pp_force_accum(:,ip1,thread)-force_pp
                                          pp2_force_ext_accum(:,ip2,thread)=pp_force_accum(:,ip2,thread)+force_pp
                                          if (pp_force_ext_flag) then
                                             xv(4:,pp1)=xv(4:,pp1)-force_pp*a_mid*G*dt
                                             xv(4:,pp2)=xv(4:,pp2)+force_pp*a_mid*G*dt
                                             
                                             write(*,*) 'Added extended force between  pp1 and pp2 = ',pp1 , pp2
                                             write(*,*) 'standing at' , xv(:3,pp1), xv(:3,pp2)
                                             write(*,*) 'in fine cells ',im,jm,km
                                             write(*,*) 'and',im+ib,jm+jb,km+kb
                                             write(*,*) 'in coarse cell ',i,j,k
                                             write(*,*) 'and tile', tile
                                          endif
                                       endif
                                       
                                    enddo
                                 endif

                              enddo
                           enddo
                        enddo


 
                     enddo ! ipl1

                  enddo !im 
               enddo !jm
            enddo !km
       

         enddo !cic_l(1)
      enddo !cic_l(2)
   enddo !cic_l(3)
   
   !write(*,*) 'Done tile' , tile ,'on thread' , thread







 end subroutine  pp_extended_force
 
