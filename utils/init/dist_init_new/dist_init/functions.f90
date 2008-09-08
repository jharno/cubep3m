! functions for dist_init.f90

  function power(kr)
    implicit none

    include 'dist_init.fh'

    real power,kr

    integer i,i1,i2
    real x,y,x1,x2,y1,y2

    i1=1
    i2=nk
    do
       if (i2-i1 .eq. 1) then
          exit
       else
          i=(i1+i2)/2
          if (kr .gt. tf(1,i)) then
             i1=i
          else
             i2=i
          endif
       endif
    enddo

    x1=log(tf(1,i1))
    y1=log(tf(3,i1))
    x2=log(tf(1,i2))
    y2=log(tf(3,i2))
    x=log(kr)
    y=y1+(y2-y1)*(x-x1)/(x2-x1)
    power=exp(y)

    return
  end function power


  function grow(a)
    implicit none
 
    include 'dist_init.par'

    real grow,a

    real hsq,oma,ola

    hsq=omegam/a**3+(1-omegam-omegal)/a**2+omegal
    oma=omegam/(a**3*hsq)
    ola=omegal/hsq
    grow=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    return
  end function grow


  function tophat(x)
    implicit none
    real tophat,x

    if (x .ne. 0) then
       tophat=3*(sin(x)-cos(x)*x)/x**3
    else
       tophat=1
    endif

    return
  end function tophat

