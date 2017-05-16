module voper

! Vector calculus:
! VROTEU:  vector rotation using Euler angles
! VROTEUT: transpose of vector rotation using Euler angles
! VROTX:   vector rotation about the x-axis
! VROTY:   vector rotation about the y-axis
! VROTZ:   vector rotation about the z-axis
! VPRO:    vector product
! SPRO:    scalar product
! VDIF:    vector difference
! VDIFN:   normalized vector difference

contains

       subroutine VROTEU(X,CA,SA)

! Vector rotation using Euler angles. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       real(8),intent(inout) :: X(3),CA(3),SA(3)

       call VROTZ(X,CA(1),SA(1))
       call VROTY(X,CA(2),SA(2))
       call VROTZ(X,CA(3),SA(3))
       end subroutine VROTEU



       subroutine VROTEUT(X,CA,SA)

! Transpose of vector rotation using Euler angles. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

       implicit none
       real(8),intent(inout) :: X(3),CA(3),SA(3)

       call VROTZ(X,CA(3),-SA(3))
       call VROTY(X,CA(2),-SA(2))
       call VROTZ(X,CA(1),-SA(1))
       end subroutine VROTEUT



       subroutine VROTX(X,c,s)

! Vector rotation about the x-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       real(8),intent(inout) :: X(3),c,s
       real(8) :: q

       q   = c*X(2)+s*X(3)
       X(3)=-s*X(2)+c*X(3)
       X(2)=q
       end subroutine VROTX



       subroutine VROTY(X,c,s)

! Vector rotation about the y-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       real(8) :: X(3),c,s
       real(8) :: q

       q   = c*X(3)+s*X(1)
       X(1)=-s*X(3)+c*X(1)
       X(3)=q
       end subroutine VROTY



       subroutine VROTZ(X,c,s)

! Vector rotation about the z-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       real(8),intent(in) :: c,s
       real(8),intent(inout) :: X(3)
       real(8) :: q

       q   = c*X(1)+s*X(2)
       X(2)=-s*X(1)+c*X(2)
       X(1)=q
       end subroutine VROTZ



       subroutine VPRO(XY,X,Y)

! Vector product. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       real(8),intent(inout) :: XY(3),X(3),Y(3)

       XY(1)=X(2)*Y(3)-X(3)*Y(2)
       XY(2)=X(3)*Y(1)-X(1)*Y(3)
       XY(3)=X(1)*Y(2)-X(2)*Y(1)    
       end subroutine VPRO



       subroutine SPRO(XY,X,Y)

! Scalar product. Version 2003-11-07.
!
! Copyright (C) 2002 Karri Muinonen

       implicit none
       real(8),intent(inout) :: XY,X(3),Y(3)

       XY=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)    
       end subroutine SPRO



       subroutine VDIF(XY,X,Y)

! Difference of two vectors. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

       implicit none
       real(8),intent(in) :: X(3),Y(3)
       real(8),intent(inout) :: XY(3)
       integer :: j1

       do 10 j1=1,3
        XY(j1)=X(j1)-Y(j1)
10     end do
       end subroutine VDIF



       subroutine VDIFN(XY,dxy,X,Y)

! Normalized difference of two vectors. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

       implicit none
       real(8),intent(in) :: X(3),Y(3)
       real(8),intent(inout) :: dxy,XY(3)
       integer :: j1
       
       dxy=0.0d0
       do 10 j1=1,3
        XY(j1)=X(j1)-Y(j1)
        dxy=dxy+XY(j1)**2
10     end do
       if (dxy.eq.0.0d0) &
       stop 'Trouble in VDIFN: zero vector.'
       dxy=sqrt(dxy)
       do 20 j1=1,3
        XY(j1)=XY(j1)/dxy
20     end do
       end subroutine VDIFN

end module voper
