      function deltac(omega0,lambda0,z)
*
* This subroutine computes the PS density threshold delta_crit 
* at redshift z, for EdS models, open models, and flat lambda models.
*
*    OMEGA0   : Density parameter
*    LAMBDA0  : Cosmological constant
*
      implicit double precision (a-h,m,o-z)
      double precision lambda0

      common /interp/ table_flat(1000,2), nbo, nbf

      if(omega0.eq.1.) then
           deltac=1.6865
c$$$      else if(lambda0.eq.0.) then
c$$$           omega=omega0*(1.+z)/(1.+omega0*z)
c$$$           do i=2,nbo
c$$$                if(table_open(i,1).gt.omega) go to 1
c$$$           enddo
c$$$           stop 'error'
c$$$    1      x1=table_open(i-1,1)
c$$$           x2=table_open(i,1)
c$$$           y1=table_open(i-1,2)
c$$$           y2=table_open(i,2)
c$$$           dx=x2-x1
c$$$           deltac=((x2-omega0)*y1+(omega0-x1)*y2)/dx
      else
           omega=omega0*(1.+z)**3/(omega0*(1.+z)**3+1.-omega0)
           do i=2,nbf
                if(table_flat(i,1).gt.omega) go to 2
           enddo
           stop 'error'
    2      x1=table_flat(i-1,1)
           x2=table_flat(i,1)
           y1=table_flat(i-1,2)
           y2=table_flat(i,2)
           dx=x2-x1
           deltac=((x2-omega0)*y1+(omega0-x1)*y2)/dx
      endif

      return
      end
