module gsaxg
! Gaussian axisymmetric sphere generator
!
! RGSAXSD: discrete spherical-coordinate representation
! RGSAXTD: discrete triangle representation
! RGSAX:   radial distance
! SGSAX:   logarithm of radial distance
! SGSAXCF: spherical harmonics coefficient and random orientation generator
       !
       ! Required modules:
       !
        use randev      ! Random deviates
        use voper       ! Vector operations
        use specfunc    ! Special functions

contains

        subroutine RGSAXSD(X,MU,PHI,ACF,CEU,SEU,rmax,beta, &
                          nthe,nphi,lmin,lmax)

! Discrete spherical-coordinate representation for a sample axisymmetric
! G-sphere. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: nthe,nphi,lmin,lmax
       real(8),intent(in) :: CEU(3),SEU(3),MU(0:180), &
         PHI(0:360),beta
       real(8),intent(in),allocatable :: ACF(:,:)
       real(8),intent(inout) :: X(0:180,0:360,3),rmax
       integer :: j1,j2
       real(8) :: r,nu

       rmax=0.0d0
       do 20 j1=0,nthe
        nu=sqrt(1.0d0-MU(j1)**2)
        do 10 j2=0,nphi
         r=RGSAX(ACF,CEU,SEU,MU(j1),PHI(j2),beta,lmin,lmax)
         X(j1,j2,1)=r*nu*cos(PHI(j2))
         X(j1,j2,2)=r*nu*sin(PHI(j2))
         X(j1,j2,3)=r*MU(j1)
         if (r.gt.rmax) rmax=r
10      end do
20     end do
       end subroutine RGSAXSD



       subroutine RGSAXTD(X,N,MU,PHI,ACF,CEU,SEU,rmax,beta, &
                         IT,nnod,ntri,lmin,lmax)

! Discrete triangle representation for a sample axisymmetri! G-sphere.
! Version 2017-05-17. Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax,ntri,nnod
       integer,intent(in),allocatable :: IT(:,:)
       real(8),intent(in) :: CEU(3),SEU(3),beta
       real(8),intent(in),allocatable :: MU(:),PHI(:),ACF(:,:)
       real(8),intent(inout),allocatable :: X(:,:),N(:,:)
       real(8),intent(inout) :: rmax
       integer :: j1,j2
       real(8) :: X1(3),X2(3),X3(3),r,nu
       
       if(.not. (allocated(IT) .or. allocated(MU) .or. allocated(PHI))) &
          stop 'trouble in RGSAXTD: arrays not allocated'
         
       if(size(mu).ne.nnod .or. size(phi).ne.nnod .or. size(IT).ne.3*ntri) stop &
          'trouble in RGSAXTD: array length mismatch'
       
       allocate(X(nnod,3),N(ntri,3))

! Node coordinates:

       rmax=0.0d0
       do 10 j1=1,nnod
        r=RGSAX(ACF,CEU,SEU,MU(j1),PHI(j1),beta,lmin,lmax)
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
       end subroutine RGSAXTD



       real(8) function RGSAX(ACF,CEU,SEU,mu,phi,beta, &
                                      lmin,lmax)

! Radial distance in a given direction for a sample axisymmetric G-sphere.
! Version 2002-12-16. Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in) :: CEU(3),SEU(3),mu,phi,beta
       real(8),intent(in),allocatable :: ACF(:,:)

       RGSAX=exp(SGSAX(ACF,CEU,SEU,mu,phi,lmin,lmax)-0.5d0*beta**2)
       end function RGSAX



       real(8) function SGSAX(ACF,CEU,SEU,mu,phi,lmin,lmax)

! Logarithmic radial distance in a given direction for a sample 
! axisymmetric G-sphere. Version 2003-11-07, revised from Version 2002-12-16.
! Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002, 2003 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in) :: mu,phi
       real(8),intent(in),allocatable :: ACF(:,:)
       integer :: l
       real(8) :: CEU(3),SEU(3),nu,cphi,sphi,r0,mu0,X(3)
       real(8),allocatable :: LEGP(:,:)

       if (lmax.eq.0) then
        SGSAX=ACF(0,0)
        return
       endif

! Euler rotation:

       nu=sqrt(1.0d0-mu**2)
       cphi=cos(phi)
       sphi=sin(phi)
       X(1)=nu*cphi
       X(2)=nu*sphi
       X(3)=mu
       call VROTEUT(X,CEU,SEU)
       
       r0=sqrt(X(1)**2+X(2)**2+X(3)**2)
       if (r0.lt.1.0d-12) &
       stop 'Trouble in SGSAX: radial distance too small.'
       mu0=X(3)/r0

! Precomputation of Legendre polynomials:
       allocate(LEGP(0:lmax,0:lmax))
       call LEGAA(LEGP,mu0,lmax,0)
       LEGP(0,0)=1.0d0

! Sum up:

       SGSAX=0.0d0
       do 10 l=lmin,lmax
        SGSAX=SGSAX+LEGP(l,0)*ACF(l,0)
10     end do

       deallocate(LEGP)
       
       end function SGSAX



       subroutine SGSAXCF(ACF,CEU,SEU,SCFSTD,lmin,lmax)

! Generates the sample spherical harmonics coefficients and Euler rotation
! angle sines and cosines for the logarithmi! radial distance of 
! the axisymmetri! G-sphere. Note that fixed orientation is obtained by
! using the Euler angle sines and cosines that are currently commented out.  
! Version 2002-11-07, revised from Version 2002-12-16.
! Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002, 2003 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in),allocatable :: SCFSTD(:,:)
       real(8),intent(inout) :: CEU(3),SEU(3)
       real(8),intent(inout),allocatable :: ACF(:,:)
       integer :: irnd,l
       real(8) :: alpha,gamma,rn,pi
       parameter (pi=3.1415926535898d0)
       common irnd

       do 10 l=lmin,lmax
        call RNDG(rn)
        ACF(l,0)=rn*SCFSTD(l,0)*sqrt(dble(2*l+1))
10     end do
       gamma=2.0d0*pi*RNDU(irnd)
       alpha=2.0d0*pi*RNDU(irnd)
       CEU(1)=cos(gamma)
       CEU(2)=1.0d0-2.0d0*RNDU(irnd)
       CEU(3)=cos(alpha)
       SEU(1)=sin(gamma)
       SEU(2)=sqrt(1.0d0-CEU(2)**2)
       SEU(3)=sin(alpha)
!       gamma=0.0d0
!       alpha=0.0d0
!       CEU(1)=1.0d0
!       CEU(2)=1.0d0
!       CEU(3)=1.0d0
!       SEU(1)=0.0d0
!       SEU(2)=0.0d0
!       SEU(3)=0.0d0
       end subroutine SGSAXCF

end module gsaxg

