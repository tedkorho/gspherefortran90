module randev

! Random deviates:
!
! RNDU: uniform distribution within [0, 1)
! RNDG: Gaussian distribution with zero mean and unit standard deviation

contains

       real(8) function RNDU(irnd)

! Quick and dirty uniform distribution [0,1) based on the TOMS
! library routine SUNIF. Initialize by calling RNDU(irnd) with 
! x = 4*n+1. Calls with negative argument return a random number.
! Version 2002-12-16.
!
! Copyright (C) 2002 Timo Nousiainen

       implicit none
       integer :: irnd
       real(8) :: r,factor,two28
       save r
       data factor /41475557.0d0/, two28 /268435456.d0/

       if (irnd.ge.0) goto 10

       r=dmod(r*factor,1.0d0)
       RNDU=dble(r)
       return

! Initialization

10     r=dble(float(irnd))/two28
       r=dmod(r*factor,1.0d0)
       RNDU=dble(r)
       irnd=-1
       return
       end function RNDU



       subroutine RNDG(r1)

! Returns a normally distributed random deviate with zero mean and 
! unit variance. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer :: flg,irnd,xrandom
       real(8) :: q1,q2,r1,r2
       save flg,r2
       data flg/0/
       common irnd

       if (flg.eq.1) then
        r1=r2
        flg=0
        return
       endif

       flg=0
10     r1=2.0d0*RNDU(irnd)-1.0d0
       r2=2.0d0*RNDU(irnd)-1.0d0
       q1=r1**2+r2**2
       if (q1.ge.1.0d0 .or. q1.le.0.0d0) goto 10

       q2=sqrt(-2.0d0*log(q1)/q1)
       r1=r1*q2
       r2=r2*q2
       flg=1
       end subroutine RNDG

end module randev

