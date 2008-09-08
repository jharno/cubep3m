      function sigma_cobe(omega0,omegabh2,lambda0,h0,tcmb,an,scale,fact)
*
* This subroutine computes sigma at a given length scale,
* using a COBE-normalized CDM spectrum.
*
*    OMEGA0   : Density parameter
*    OMEGABH2 : Baryon density parameter in units of h**(-2)
*    LAMBDA0  : Cosmological constant
*    H0       : Hibble constants in km/s/Mpc
*    TCMB     : Temperature of CMB in Kelvin
*    AN       : Slope of primordial spectrum
*    SCALE    : Scale of interest
*
      implicit double precision (a-h,m,o-z)
      parameter (c=2.9979d+10)
      external integrant
      double precision lambda0, integrant
      
      integer max_size
      parameter (max_size = 500)
      real*8 k(max_size), Deltak1(max_size), Deltak(max_size)

      common /transfer/ expn, coef1, coef2, anusq, oh2, betac, t27

      common/power/ k,Deltak1,Deltak,h   

      data tol /1.e-8/

      h=h0/100.
      oh2=omega0*h**2
      omegab=omegabh2/h**2
      obo0=omegab/omega0
      t27=tcmb/2.7
      expn=an

      b1=0.313*oh2**(-0.419)*(1.+0.607*oh2**0.674)
      b2=0.238*oh2**0.223
      zd=1291.*oh2**0.251/(1.+0.659*oh2**0.828)*(1.+b1*omegabh2**b2)
      zeq=2.50e+04*oh2/t27**4
      yd=(1.+zeq)/(1.+zd)
      s=44.5*dlog(9.83/oh2)/dsqrt(1.+10.*omegabh2**0.75)

      fc=(omega0-omegab)/omega0
      fcb=1.
      fb=omegab/omega0
      fnb=fb
      pc=0.25*(5.-dsqrt(1.+24.*fc))

      betac=1./(1.-0.949*fnb)

      anu=(fc/fcb)*((5.-2.*pc)/5.)*(1.-0.553*fnb+0.126*fnb**3)
     +    *(1.+yd)**(-pc)*(1.+0.5*pc*(1.+1./(3.-4.*pc)/7.)/(1.+yd))
      anusq=dsqrt(anu)

* Compute conversion factor between x, q, and ks

      coef1=t27**2/scale
      coef2=s/scale

* Compute sigma

      down=0.
      up=1.
      call INTEG(down,up,integrant,area)
      do i=1,20
           areaold=area
           down=up
           up=2.*up
           call INTEG(down,up,integrant,area1)
           area=area+area1
           if(dabs((area-areaold)/areaold).lt.tol) go to 1
c           print*,area,areaold,(area-areaold)/areaold
      enddo
      stop 'no convergence'

c    1 ampl=DELTAH(omega0,lambda0,an)**2*(c/(1.d+07*scale*h))**(3.+an)
    1 ampl = 1.0d0
      sigma_cobe=dsqrt(ampl*area)

      return
      end
*=========================================================================
      function deltah(omega0,lambda0,an)
      implicit double precision (a-h,o-z)
      double precision lambda0

      an1=an-1.
      if(lambda0.eq.0.) then
           deltah=1.95e-05*omega0**(-0.35-0.19*dlog(omega0)-0.17*an1)
     +            *dexp(-an1-0.14*an1**2)
      else
           deltah=1.94e-05*omega0**(-0.785-0.05*dlog(omega0))
     +            *dexp(-0.95*an1-0.169*an1**2)
      endif 

      return
      end
*====================================================================
      subroutine integ(a,b,F,area)
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
      function integrant(x)
      implicit double precision (a-z)
      parameter (e=2.718281828452)

      integer max_size
      parameter (max_size = 500)
      real*8 k(max_size), Deltak1(max_size), Deltak(max_size)
      real*8 k1, t27, scale1

      common/power/ k,Deltak1,Deltak,h    

      common /transfer/ expn, coef1, coef2, anusq, oh2, betac, t27

      ks=coef2*x
      scale1 = t27**2/coef1
c
      gammaeff=oh2*(anusq+(1.-anusq)/(1.+(0.43*ks)**4))
      qeff=coef1*x/gammaeff

      if(qeff.eq.0.) then
c     t=1.
         delta=1.0d0
         w=1.
         integrant = 0.0d0
      else
         k1=h*x/scale1
         call splint(k,Deltak1,Deltak,max_size,k1,Delta)
         w=9.*(dsin(x)-x*dcos(x))**2/x**6
         integrant=1./x*Delta*w
      endif

c      integrant=x**(2.+expn)*t**2*w
c         integrant=Delta*w!/2./3.1415**2

      return
      end
