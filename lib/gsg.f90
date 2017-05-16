module gsg
! Gaussian sphere generator
!
! RGSSD:   discrete spherical-coordinate representation
! RGSTD:   discrete triangle representation
! RGS:     radial distance
! SGS:     logarithm of radial distance
! SGSCF:   spherical harmonics coefficient generation
        ! Modules used:
        use specfunc    ! Special functions.
        use randev      ! Random deviates.
        use voper       ! Vector operations

contains


        subroutine RGSSD(X,MU,PHI,ACF,BCF,rmax,beta, &
                        nthe,nphi,lmin,lmax)

! Discrete spherical-coordinate representation for a sample G-sphere.
! Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: nthe,nphi,lmin,lmax
       real(8),intent(in) :: MU(0:180),PHI(0:360),ACF(0:256,0:256), &
         BCF(0:256,0:256),beta
       real(8),intent(inout) :: X(0:180,0:360,3),rmax
       integer :: j1,j2
       real(8) :: r,nu

       rmax=0.0d0
       do 20 j1=0,nthe
        nu=sqrt(1.0d0-MU(j1)**2)
        do 10 j2=0,nphi
         r=RGS(ACF,BCF,MU(j1),PHI(j2),beta,lmin,lmax)
         X(j1,j2,1)=r*nu*cos(PHI(j2))
         X(j1,j2,2)=r*nu*sin(PHI(j2))
         X(j1,j2,3)=r*MU(j1)
         if (r.gt.rmax) rmax=r
10      end do
20     end do
       end subroutine RGSSD



       subroutine RGSTD(X,N,MU,PHI,ACF,BCF,rmax,beta, &
                       IT,nnod,ntri,lmin,lmax)

! Discrete triangle representation for a sample G-sphere.
! Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: nnod,ntri,lmin,lmax,IT(260000,3)
       real(8),intent(in) :: MU(130000),ACF(0:256,0:256), &
         BCF(0:256,0:256),PHI(130000),beta
       real(8),intent(inout) :: X(130000,3),N(260000,3),rmax
       integer :: j1,j2
       real(8) :: X1(3),X2(3),X3(3),r,nu

! Node coordinates:

       rmax=0.0d0
       do 10 j1=1,nnod
        r=RGS(ACF,BCF,MU(j1),PHI(j1),beta,lmin,lmax)
        nu=sqrt(1.0d0-MU(j1)**2)
        X(j1,1)=r*nu*cos(PHI(j1))
        X(j1,2)=r*nu*sin(PHI(j1))
        X(j1,3)=r*MU(j1)
        if (r.ge.rmax) rmax=r
10     end do

! Outer unit triangle normals:

       do 50 j1=1,ntri
        do 20 j2=1,3
         r=X(IT(j1,1),j2)
         X1(j2)=X(IT(j1,2),j2)-r
         X2(j2)=X(IT(j1,3),j2)-r
20      end do
        call VPRO(X3,X1,X2)
        r=sqrt(X3(1)**2+X3(2)**2+X3(3)**2)
        do 40 j2=1,3
         N(j1,j2)=X3(j2)/r
40      end do
50     end do
       end subroutine RGSTD



       real(8) function RGS(ACF,BCF,mu,phi,beta,lmin,lmax)

! Radial distance in a given direction for a sample G-sphere.
! Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in) :: ACF(0:256,0:256),BCF(0:256,0:256), &
       mu,phi,beta 
       ! real(8) :: SGS

       RGS=exp(SGS(ACF,BCF,mu,phi,lmin,lmax)-0.5d0*beta**2)
       end function RGS



       real(8) function SGS(ACF,BCF,mu,phi,lmin,lmax)

! Logarithmi! radial distance in a given direction for a sample G-sphere.
! Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in) :: ACF(0:256,0:256),BCF(0:256,0:256), &
       mu,phi 
       integer :: l,m
       real(8) :: LEGP(0:256,0:256),CPHI(256),SPHI(256)

       if (lmax.eq.0) then
        SGS=ACF(0,0)
        return
       endif

! Precomputation of sines, cosines, and associated Legendre functions:

       call LEGA(LEGP,mu,lmax,0)
       do 10 m=1,lmax
        call LEGA(LEGP,mu,lmax,m)
        CPHI(m)=cos(m*phi)
        SPHI(m)=sin(m*phi)
10     end do
       LEGP(0,0)=1.0d0

! Sum up:

       SGS=0.0d0
       do 20 l=lmin,lmax
        SGS=SGS+LEGP(l,0)*ACF(l,0)
20     end do
       do 40 m=1,lmax
        do 30 l=max(m,lmin),lmax
         SGS=SGS+ &
         LEGP(l,m)*(ACF(l,m)*CPHI(m)+BCF(l,m)*SPHI(m))
30      end do
40     end do
       end function SGS


       subroutine SGSCF(ACF,BCF,SCFSTD,lmin,lmax)

! Generates the sample spherical harmonics coefficients for the 
! logarithmi! radial distance of the G-sphere. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer :: l,lmin,lmax,m
       real(8) :: ACF(0:256,0:256),BCF(0:256,0:256), &
       SCFSTD(0:256,0:256),rn

       do 10 l=lmin,lmax
        call RNDG(rn)
        ACF(l,0)=rn*SCFSTD(l,0)
        BCF(l,0)=0.0d0
10     end do
       do 30 m=1,lmax
        do 20 l=max(m,lmin),lmax
         call RNDG(rn)
         ACF(l,m)=rn*SCFSTD(l,m)
         call RNDG(rn)
         BCF(l,m)=rn*SCFSTD(l,m)
20      end do
30     end do
       end subroutine SGSCF

end module gsg

