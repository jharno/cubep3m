      function growth(omega0,lambda0,z)
*
* This subroutine computes the linear growth factor between redshift z
* and the present, for EdS models, open models, and flat lambda models.
*
*    OMEGA0   : Density parameter
*    LAMBDA0  : Cosmological constant
*
      implicit double precision (a-h,m,o-z)
      external integrant0
      double precision lambda0, integrant0

      if(omega0.eq.1.) then
           growth=z+1.
      else if(lambda0.eq.0.) then
           x0=1./omega0-1.
           d0=1.+3./x0+3.*dsqrt((1.+x0)/x0**3)
     +                   *dlog(dsqrt(1.+x0)-dsqrt(x0))
           x1=x0/(z+1.)
           d1=1.+3./x1+3.*dsqrt((1.+x1)/x1**3)
     +                   *dlog(dsqrt(1.+x1)-dsqrt(x1))
           growth=d0/d1
      else
           down=0.
           a0=(lambda0/omega0)**0.3333333333333
           call INTEG0(down,a0,integrant0,area)
           d0=2.5*dsqrt(1.+1./a0**3)*area
           a1=a0/(z+1.)
           call INTEG0(down,a1,integrant0,area)
           d1=2.5*dsqrt(1.+1./a1**3)*area
           growth=d0/d1
      endif

      return
      end

*====================================================================
      subroutine integ0(a,b,F,area)
*
* This subroutine computes integrals using Simpson's rule.
*
      implicit double precision (a-h,o-z)

      data tol /1.e-08/

      if(a.eq.b) then
           area=0.
           return
      endif

      areaold=0.
      sumend=F(a)+F(b)
      sumeven=0.
      sumodd=0.

      do n=1,25
           i=2.**n
           h=(b-a)/dfloat(i)
           sumodd=0.
           do j=1,i-1,2
                c=a+j*h
                sumodd=sumodd+F(c)
           enddo
           area=(sumend+4.*sumodd+2.*sumeven)*h/3.
           if(dabs((area-areaold)/area).lt.tol) return
           sumeven=sumeven+sumodd
           areaold=area
      enddo

      write(6,1000)
 1000 format(/5x,'Error, no convergence.')
      stop

      end
*=========================================================================
      function integrant0(x)
      implicit double precision (a-z)

      integrant0=(x/(1.+x**3))**1.5
      
      return
      end
