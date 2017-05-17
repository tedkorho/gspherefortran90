module gsg
! Gaussian sphere generator
!
! RGSSD:   discrete spherical-coordinate representation
! RGSTD:   discrete triangle representation
! RGS:     radial distance
! SGS:     logarithm of radial distance
! SGSCF:   spherical harmonics coefficient generation
       !
       ! Required modules:
       !
        use specfunc    ! Special functions.
        use randev      ! Random deviates.
        use voper       ! Vector operations.

contains


        subroutine RGSSD(X,MU,PHI,ACF,BCF,rmax,beta, &
                        nthe,nphi,lmin,lmax)

! Discrete spherical-coordinate representation for a sample G-sphere.
! Version 2002-12-16. Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: nthe,nphi,lmin,lmax
       real(8),intent(in) :: MU(0:180),PHI(0:360),beta
       real(8),intent(in),allocatable :: ACF(:,:),BCF(:,:)
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
! Version 2017-05-17. Updated with dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: nnod,ntri,lmin,lmax
       integer,intent(in),allocatable :: IT(:,:)
       real(8),intent(in) :: beta
       real(8),intent(in),allocatable :: ACF(:,:),BCF(:,:),MU(:),PHI(:)
       real(8),intent(inout) :: rmax
       real(8),intent(inout),allocatable :: X(:,:),N(:,:)
       integer :: j1,j2
       real(8) :: X1(3),X2(3),X3(3),r,nu
       
       if(.not. (allocated(IT) .or. allocated(MU) .or. allocated(PHI))) &
          stop 'trouble in RGSTD: arrays not allocated'
         
       if(size(mu).ne.nnod .or. size(phi).ne.nnod .or. size(IT).ne.3*ntri) stop &
          'RGSTD: array length mismatch'
       
       allocate(X(nnod,3),N(ntri,3))

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
! Version 2002-12-16. Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in) :: mu,phi,beta
       real(8),intent(in),allocatable :: ACF(:,:),BCF(:,:)

       RGS=exp(SGS(ACF,BCF,mu,phi,lmin,lmax)-0.5d0*beta**2)
       end function RGS



       real(8) function SGS(ACF,BCF,mu,phi,lmin,lmax)

! Logarithmi! radial distance in a given direction for a sample G-sphere.
! Version 2002-12-16.
! Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in) :: mu,phi
       real(8),intent(in),allocatable :: ACF(:,:),BCF(:,:)

       integer :: l,m
       real(8) :: CPHI(256),SPHI(256)
       real(8),allocatable :: LEGP(:,:)

       if (lmax.eq.0) then
        SGS=ACF(0,0)
        return
       endif

! Precomputation of sines, cosines, and associated Legendre functions:
       allocate(LEGP(0:lmax,0:lmax))
       call LEGAA(LEGP,mu,lmax,0)
       do 10 m=1,lmax
        call LEGAA(LEGP,mu,lmax,m)
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
       
       deallocate(LEGP)
       
       end function SGS


       subroutine SGSCF(ACF,BCF,SCFSTD,lmin,lmax)

! Generates the sample spherical harmonics coefficients for the 
! logarithmi! radial distance of the G-sphere. Version 2002-12-16.
! Updated for dynamic arrays by Teo Korhonen.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       integer,intent(in) :: lmin,lmax
       real(8),intent(in),allocatable :: SCFSTD(:,:)
       real(8),intent(inout),allocatable :: ACF(:,:),BCF(:,:)
       integer :: l,m
       real(8) :: rn

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
