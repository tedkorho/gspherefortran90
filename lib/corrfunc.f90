module corrfunc
! Correlation functions and related coefficients:
!
! SGSCFSTD: standard deviations for the spherical harmonics coefficients
! CSELL   : correlation length and angle, derivatives
! CSLEGP  : Legendre expansion of the correlation function
! CS1CF   : power law correlation function
! CS2CF   : modified Gaussian correlation function
! CS3CF   : input correlation function coefficients from file
        ! Modules used:
        use specfunc    ! Special functions.

contains

       subroutine SGSCFSTD(SCFSTD,CSCF,beta,lmin,lmax)

! Generates the standard deviations for the spherical harmonics 
! coefficients of the logarithmi! radial distance. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,m
       real(8) SCFSTD(0:256,0:256),CSCF(0:256),beta

       if (lmin.eq.0) then
        SCFSTD(0,0)=beta*sqrt(CSCF(0))
        do 20 l=1,lmax
         SCFSTD(l,0)=beta*sqrt(CSCF(l))
         do 10 m=1,l
          SCFSTD(l,m)=SCFSTD(l,0)* &
          sqrt(2.0d0*FACTI(l-m)/FACTI(l+m))
10       end do
20      end do
       else
        do 40 l=lmin,lmax
         SCFSTD(l,0)=beta*sqrt(CSCF(l))
         do 30 m=1,l
          SCFSTD(l,m)=SCFSTD(l,0)* &
          sqrt(2.0d0*FACTI(l-m)/FACTI(l+m))
30       end do
40      end do
       endif
       end



       subroutine CSELL(CSCF,gam,ell,cs2d,cs4d,lmin,lmax)

! Computes the correlation angle and length, and second and fourth 
! derivatives for a correlation function expressed in Legendre series with
! normalized coefficients. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax
       real(8) CSCF(0:256),gam,ell,cs2d,cs4d

       cs2d=0.0d0
       cs4d=0.0d0
       do 10 l=lmin,lmax
        cs2d=cs2d-CSCF(l)*l*(l+1.0d0)/2.0d0
        cs4d=cs4d+CSCF(l)*l*(l+1.0d0)*(3.0d0*l**2+3.0d0*l-2.0d0)/8.0d0
10     end do
       if (cs2d.eq.0.0d0) then
        ell=2.0d0
       else
        ell=1.0d0/sqrt(-cs2d)
       endif
       gam=2.0d0*asin(0.5d0*ell)
       end



       real(8) function CSLEGP(CSCF,xi,lmin,lmax)

! Computes the Legendre series expansion for the 
! correlation function. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax
       real(8) LEGP(0:256,0:256),CSCF(0:256),xi

       call LEGA(LEGP,xi,lmax,0)
       LEGP(0,0)=1.0d0

       CSLEGP=0.0d0
       do 10 l=lmin,lmax
        CSLEGP=CSLEGP+LEGP(l,0)*CSCF(l)
10     end do
       end



       subroutine CS1CF(CSCF,nu,lmin,lmax)

! Returns the Legendre coefficients for the correlation 
! function with power-law Legendre coefficients. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax
       real(8) CSCF(0:256),nu,norm

       do 10 l=0,lmin-1
        CSCF(l)=0.0d0
10     end do

       norm=0.0d0
       do 20 l=lmin,lmax
        CSCF(l)=1.0d0/l**nu
        norm=norm+CSCF(l)
20     end do

       do 30 l=lmin,lmax
        CSCF(l)=CSCF(l)/norm
30     end do
       end



       subroutine CS2CF(CSCF,ell,lmin,lmax)

! Returns the Legendre coefficients for the modified Gaussian correlation 
! function. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax
       real(8) BESISE(0:256),CSCF(0:256),ell,z,norm

       z=1.0d0/ell**2
       call BESMS(BESISE,z,lmax)

       do 10 l=0,lmin-1
        CSCF(l)=0.0d0
10     end do

       norm=0.0d0
       do 20 l=lmin,lmax
        CSCF(l)=(2*l+1)*BESISE(l)
        norm=norm+CSCF(l)
20     end do

       do 30 l=lmin,lmax
        CSCF(l)=CSCF(l)/norm
30     end do
       end



       subroutine CS3CF(CSCF,lmin,lmax)

! Inputs Legendre coefficients for the correlation 
! function from a file (and normalizes). Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,ll,lmin,lmax
       real(8) CSCF(0:256),norm

       open(unit=1, file='cscf.dat', status='old')

       do 10 l=0,lmin-1
        CSCF(l)=0.0d0
10     end do

       norm=0.0d0
       do 20 l=lmin,lmax
        read (1,*) ll,CSCF(l) 
        if (l.ne.ll) stop &
        'Trouble in CS3CF: degree inconsistency.'
        if (CSCF(ll).lt.0.0d0) stop &
        'Trouble in CS3CF: negative Legendre coefficient.'
        norm=norm+CSCF(l)
20     end do
       close(unit=1)

       if (norm.eq.0.0d0) stop &
       'Trouble in CS3CF: no nonzero Legendre coefficients.'

       do 30 l=lmin,lmax
        CSCF(l)=CSCF(l)/norm
30     end do
       end
end module corrfunc


