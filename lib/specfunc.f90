module specfunc

! Special functions:
!
! BESJ0A: Bessel function  J_0
! BESJ1A: Bessel function  J_1
! BESMS:  modified spherical Bessel functions multiplied by exponential
! LEGA:   associated Legendre functions
! SPHGEN: generalized spherical functions 
! GAMLN:  logarithmi! Gamma function 
! FACTI:  factorial function


contains


       real(8) function BESJ0A(x)

! Bessel function J0 using a polynomial approximation (accuracy
! better than 10(-7)). Version 2003-11-10.
!
! Copyright (C) 2003 Karri Muinonen

       implicit none
       real(8) :: x,y,f0,the0

       y=x/3.0d0
       if  (x.le.3.0d0) then
        y=x/3.0d0
        BESJ0A=1.0d0-2.2499997d0*y**2+1.2656208d0*y**4- &
        0.3163866d0*y**6+0.0444479d0*y**8-0.0039444d0*y**10+ &
        0.0002100d0*y**12
       else
        y=3.0d0/x
        f0=0.79788456d0-0.00000077d0*y-0.00552740d0*y**2- &
        0.00009512d0*y**3+0.00137237d0*y**4-0.00072805d0*y**5+ &
        0.00014476d0*y**6
        the0=x-0.78539816d0-0.04166397d0*y-0.00003954d0*y**2+ &
        0.00262573d0*y**3-0.00054125d0*y**4-0.00029333d0*y**5+ &
        0.00013558d0*y**6
        BESJ0A=f0*cos(the0)/sqrt(x)
       endif
       end function BESJ0A



       real(8) function BESJ1A(x)

! Bessel function J1 using a polynomial approximation (accuracy
! better than 10(-7)). Version 2003-11-10.
!
! Copyright (C) 2003 Karri Muinonen

       implicit none
       real(8) :: x,y,f1,the1

       y=x/3.0d0
       if  (x.le.3.0d0) then
        y=x/3.0d0
        BESJ1A=x*(0.5d0-0.56249985d0*y**2+0.21093573d0*y**4- &
        0.03954289d0*y**6+0.00443319d0*y**8-0.00031761d0*y**10+ &
        0.00001109d0*y**12)
       else
        y=3.0d0/x
        f1=0.79788456d0+0.00000156d0*y+0.01659667d0*y**2+ &
        0.00017105d0*y**3-0.00249511d0*y**4+0.00113653d0*y**5- &
        0.00020033d0*y**6
        the1=x-2.35619449d0+0.12499612d0*y+0.00005650d0*y**2- &
        0.00637879d0*y**3+0.00074348d0*y**4+0.00079824d0*y**5- &
        0.00029166d0*y**6
        BESJ1A=f1*cos(the1)/sqrt(x)
       endif
       end function BESJ1A



       subroutine BESMS(BESISE,x,n)

! Generates modified spherical Bessel functions multiplied by an
! exponential: i_0(x)*exp(-x),...,i_n(x)*exp(-x). Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer :: j1,n,nini
       real(8) :: BESISE(0:256),x,b,b0,b1,renorm
      
! Orders n=0 and n=1:

       if (n.le.1) then
        BESISE(0)=exp(-x)*sinh(x)/x
        BESISE(1)=exp(-x)*(-sinh(x)/x**2)+cosh(x)/x
        return
       endif

! Downward recurrence:

       nini=max(n+4,int(1.5d0*x))
       b1=0.0d0
       b0=exp(-x)*2.0d0*x
       do 10 j1=nini,n,-1
        b=(2*j1+1)*b0/x+b1
        b1=b0
        b0=b
10     end do
       BESISE(n)=b1
       BESISE(n-1)=b0
       do 20 j1=n,2,-1
        BESISE(j1-2)=(2*j1-1)*BESISE(j1-1)/x+BESISE(j1)
20     end do

! Renormalization:

       renorm=exp(-x)*(sinh(x)/x)/BESISE(0)
       do 30 j1=0,n
        BESISE(j1)=renorm*BESISE(j1)
30     end do
       end subroutine BESMS



       subroutine LEGA(LEGP,x,lmax,m)

! Computes associated Legendre functions from degree l=m
! up to l=lmax. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer :: lmax,l,m
       real(8) :: LEGP(0:256,0:256),x
       double complex GSP(0:256,0:256,-2:2),i

       i=dcmplx(0.0d0,1.0d0)

! Check degree, orders, and argument:

       if (lmax.lt.0) &
       stop 'Trouble in LEGA: degree negative.'

       if (m.gt.lmax .or. m.lt.0) &
       stop 'Trouble in LEGA: order out of range.'

       if (abs(x).gt.1.0d0) &
       stop 'Trouble in LEGA: argument out of range.'

! Compute associated Legendre functions with the help of
! the generalized spherical functions:
   
       call SPHGEN(GSP,x,lmax,m,0)

       do 10 l=m,lmax
        LEGP(l,m)=dreal(i**m* &
        sqrt(FACTI(l+m)/FACTI(l-m))*GSP(l,m,0))
10     end do
       end subroutine LEGA



       subroutine SPHGEN(GSP,x,lmax,m1,m2)

! Computes generalized spherical functions from degree max(abs(m1),abs(m2))
! up to lmax. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer :: lmax,l,m0,m1,m2,m12,p12
       real(8) :: x
       double complex GSP(0:256,0:256,-2:2),i

       i=dcmplx(0.0d0,1.0d0)

! Check degree, orders, and argument:

       if (lmax.lt.0) &
       stop 'Trouble in SPHGEN: degree negative.'

       if (abs(m1).gt.lmax .or. abs(m2).gt.min(2,lmax) .or. m1.lt.0) &
       stop 'Trouble in SPHGEN: order out of range.'

       if (abs(x).gt.1.0d0) &
       stop 'Trouble in SPHGEN: argument out of range.'

! Compute generalized spherical functions:
   
       m0=max(abs(m1),abs(m2))
       m12=abs(m1-m2)
       p12=abs(m1+m2)

       if (m0.gt.0) then

        if (m12.ne.0 .and. p12.ne.0) then
         GSP(m0,m1,m2)=(-i)**m12/2.0d0**m0* &
         sqrt(FACTI(2*m0)/(FACTI(m12)*FACTI(p12))* &
         (1.0d0-x)**m12*(1.0d0+x)**p12)
        elseif (m12.eq.0) then
         GSP(m0,m1,m2)=1.0d0/2.0d0**m0* &
        sqrt(FACTI(2*m0)/FACTI(p12)*(1.0d0+x)**p12)
        else
         GSP(m0,m1,m2)=(-i)**m12/2.0d0**m0* &
         sqrt(FACTI(2*m0)/FACTI(m12)*(1.0d0-x)**m12)
        endif

        if (m0.eq.lmax) return

        GSP(m0+1,m1,m2)=(2*m0+1)*(m0*(m0+1)*x-m1*m2)*GSP(m0,m1,m2)/ &
        (m0*sqrt(dble((m0+1)**2-m2**2))*sqrt(dble((m0+1)**2-m1**2)))

        if (m0+1.eq.lmax) return

        do 10 l=m0+1,lmax-1
         GSP(l+1,m1,m2)=((2*l+1)*(l*(l+1)*x-m1*m2)*GSP(l,m1,m2) &
         -(l+1)*sqrt(dble((l**2-m1**2)*(l**2-m2**2)))*GSP(l-1,m1,m2))/ &
         (l*sqrt(dble(((l+1)**2-m1**2)*((l+1)**2-m2**2))))
10      end do

       else

        GSP(0,0,0)=1.0d0
        if (lmax.eq.0) return
        GSP(1,0,0)=x
        if (lmax.eq.1) return

        do 20 l=m0+1,lmax-1
         GSP(l+1,0,0)=((2*l+1)*x*GSP(l,0,0)-l*GSP(l-1,0,0))/(l+1)
20      end do

       endif
       end subroutine SPHGEN



       function GAMLN(x)

! Function GAMLN computes the natural logarithm of the
! Gamma function with positive real argument. The algorithm
! is based on the approximation by Lanczos (1964: SIAM
! Journal on Numerical Analysis, ser. B, vol 1., pp. 86-96)
! Version 2002-12-16.
!
! Copyright (C) 2002 Timo Nousiainen

       implicit none
       real(8) :: GAMLN,stp,pi,x,xx,temp,arg1,arg2
       parameter (pi=3.1415926535898d0)

       if (x.le.0.0d0) then
        write (*,*) 'Trouble in GAMLN: x must be positive! Exiting.'
        stop
       endif

       xx=x-1.0d0

       stp=sqrt(2.0d0*pi)
       arg1=xx+0.5d0
       arg2=xx+5.5d0
 
       temp=1.000000000190015d0+76.18009172947146d0/(xx+1.0d0)- &
           86.50532032941677d0/(xx+2.0d0)+ &
           24.01409824083091d0/(xx+3.0d0)- &
           1.231739572450155d0/(xx+4.0d0)+ &
           0.1208650973866179d-2/(xx+5.0d0)- &
           0.5395239384953d-5/(xx+6.0d0)

       temp=dlog(temp*stp)+(arg1)*dlog(arg2)-(arg2)
       GAMLN=temp
       end function GAMLN



       function FACTI(x)

! Function FACTI returns a factorial of an integer argument x
! Values are precomputed for small values for speed. Although
! a factorial of an integer argument is also an integer, it is
! handled as real number with real(8) to handle
! large values. Version 2002-12-16.
!
! Copyright (C) 2002 Timo Nousiainen

       implicit none
       integer :: x 
       real(8) :: FACTI,xx

       if (x.lt.0) then
        write (*,*) 'Trouble in FACTI: x must be non-negative! Exiting.'
        stop
       endif

       if (x.le.14) then
        if (x.eq.0)  FACTI=1.0d0
        if (x.eq.1)  FACTI=1.0d0
        if (x.eq.2)  FACTI=2.0d0
        if (x.eq.3)  FACTI=6.0d0
        if (x.eq.4)  FACTI=24.0d0
        if (x.eq.5)  FACTI=120.0d0
        if (x.eq.6)  FACTI=720.0d0
        if (x.eq.7)  FACTI=5040.0d0
        if (x.eq.8)  FACTI=40320.0d0
        if (x.eq.9)  FACTI=362880.0d0
        if (x.eq.10) FACTI=3628800.0d0
        if (x.eq.11) FACTI=39916800.0d0
        if (x.eq.12) FACTI=479001600.0d0
        if (x.eq.13) FACTI=6227020800.0d0
        if (x.eq.14) FACTI=87178291200.0d0
       else
        xx=x*1.0d0
        FACTI=exp(GAMLN(xx+1.0d0))
       endif
       end function FACTI
end module specfunc
