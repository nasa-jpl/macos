C***********************************************************************
C	Begin file mathsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

      MODULE math_mod
      IMPLICIT NONE

      PUBLIC
      PRIVATE :: CSWAP, DROTG, DDOTLP, CNRM2, zrhqr, balanc, hqr
      CONTAINS

      SUBROUTINE RZERO(A,N)
      SREAL A(1)
      Integer :: N
      A(1:N)=0E0
      RETURN
      END SUBROUTINE RZERO
C
      SUBROUTINE DADD(A,B,C,N)
      Integer :: N
      REAL*8 A(N),B(N),C(N)
      A(1:N)=B(1:N)+C(1:N)
      RETURN
      END SUBROUTINE DADD
C
      ! Added by jzlou
      SUBROUTINE CMATSUM(A,S,M,N)
      Complex*16, dimension(1:M,1:N) :: A
      Complex*16 :: S
      Integer :: M,N,i,j
      S=DCMPLX(0d0,0d0)
      Do i=1,M
	Do j=1,N
	  S=S+A(i,j)
	End Do
      End Do
      END SUBROUTINE CMATSUM
C
      ! Sum all elements of a matrix
      ! Added by jzlou
      SUBROUTINE DMATSUM(A,S,M,N)
      Real*8, dimension(1:M,1:N) :: A
      Real*8 :: S
      Integer :: M,N,i,j
      S=0d0
      Do i=1,M
        Do j=1,N
          S=S+A(i,j)
        End Do
      End Do
      END SUBROUTINE DMATSUM
C
      FUNCTION RVECSUM(V,N)
      Real*8, dimension(1:N) :: V
      Real*8 :: RVECSUM
      Integer :: N,i
      RVECSUM=0d0
      DO i=1,N
        RVECSUM=RVECSUM+V(i)
      END DO
      END FUNCTION RVECSUM

C
      SUBROUTINE DSUB(A,B,C,N)
      REAL*8 A(N),B(N),C(N)
      Integer :: N
      A(1:N)=B(1:N)-C(1:N)
      RETURN
      END SUBROUTINE DSUB
C
      SUBROUTINE DSMPROD(A,B,S,N)
      REAL*8 A(1),B(1),S
      Integer :: N
      A(1:N)=S*B(1:N)
      RETURN
      END SUBROUTINE DSMPROD
C
      SUBROUTINE DNEGATE(A,B,N)
      REAL*8 :: A(1),B(1)
      Integer :: N
      A(1:N)=-B(1:N)
      RETURN
      END SUBROUTINE DNEGATE
C
      SUBROUTINE DEQUATE(A,B,N)
      REAL*8 A(1),B(1)
      Integer :: N
      A(1:N)=B(1:N)
      RETURN
      END SUBROUTINE DEQUATE
C
      SUBROUTINE REQUATE(A,B,N)
      SREAL A(1),B(1)
      Integer :: N
      A(1:N)=B(1:N)
      RETURN
      END SUBROUTINE REQUATE

C added by sws ...
      SUBROUTINE rdEQUATE(A,B,N)
      SREAL A(1)
      REAL*8 B(1)
      Integer :: N
      A(1:N)=B(1:N)
      RETURN
      END SUBROUTINE rdEQUATE

      SUBROUTINE drdSUB(A,B,C,N)
      REAL*8 A(N), C(N)
      SREAL B(N)
      Integer :: N
      A(1:N)=B(1:N)-C(1:N)
      RETURN
      END SUBROUTINE drdSUB
C
      SUBROUTINE DUNITVEC(YHAT,Y,DMAG,N)
      REAL*8 YHAT(1),Y(1),DMAG
      Integer :: N,I
      DMAG=0D0
      DO 1 I=1,N
	DMAG=DMAG+Y(I)**2
    1 CONTINUE
      DMAG=DSQRT(DMAG)
      IF (DMAG.GT.0d0) THEN
        DO 2 I=1,N
	  YHAT(I)=Y(I)/DMAG
    2   CONTINUE
      END IF
      RETURN
      END SUBROUTINE DUNITVEC
C
      SUBROUTINE DPERPUVEC(U,X,Y)
C  U is a unit vector in the direction of the component of X
C  perpendicular to Y: U=unit(perp(Y)*X)
C  -- assuming Y is unit vector -jzlou
      REAL*8 U(3),X(3),Y(3),V(3),S
      S=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
      V(1)=X(1)-S*Y(1)
      V(2)=X(2)-S*Y(2)
      V(3)=X(3)-S*Y(3)
      S=SQRT(V(1)*V(1)+V(2)*V(2)+V(3)*V(3))
      U(1)=V(1)/S
      U(2)=V(2)/S
      U(3)=V(3)/S
      RETURN
      END SUBROUTINE DPERPUVEC
C
      FUNCTION DMAG(X,N)
      REAL*8 X(N),DMAG
      Integer :: N,I
      DMAG=0D0
        DO 1 I=1,N
          IF (ABS(X(I)).GT.1D30) THEN
            DMAG=1D30
          ELSE
            DMAG=DMAG+X(I)*X(I)
          END IF
    1   CONTINUE
        DMAG=DSQRT(DMAG)
      RETURN
      END FUNCTION DMAG
C
      FUNCTION DDOTC(X,Y)
      REAL*8 X(3),Y(3),DDOTC
      DDOTC=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
      RETURN
      END FUNCTION DDOTC

      FUNCTION DDOTC2(X,Y)
      REAL*8 X(2),Y(2),DDOTC2
      DDOTC2=X(1)*Y(1)+X(2)*Y(2)
      RETURN
      END FUNCTION DDOTC2

      FUNCTION DDOTN(X,Y,N)
      INTEGER N,i
      REAL*8 X(N),Y(N),DDOTN
      DDOTN=0d0
      DO i=1,N
        DDOTN=DDOTN+X(i)*Y(i)
      End Do
      RETURN
      END FUNCTION DDOTN

C
      SUBROUTINE DUNITIZE(X)
      REAL*8 X(3),Y
      Y=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      IF (Y.GT.1D-10) THEN
        X(1)=X(1)/Y; X(2)=X(2)/Y; X(3)=X(3)/Y
      ELSE
        X(1:3)=0D0
      END IF
      RETURN
      END SUBROUTINE DUNITIZE
C

C ************** OLD and WRONG !!! ***************************************
C jzlou - compute matrix product A = C*B,
C         where C is NCxNB, B is NBxNA, and A is NCxNA
C         A, B, C are stored as 1-D arrays, and
C         matrices are formed row by row.
C         The above was verified on 01/04/2005 in
C         /proj/jwst/jzlou/macos-vY/macos-v3.2b/frm-lnx-taper/tst_DMPROD.F
C ************************************************************************

C
C This is the correct one - 03/16/05
C jzlou - compute matrix product A = B*C,
C         where B is NAxNB, C is NBxNC
C         and A is NAxNC

      SUBROUTINE DMPROD(A,B,C,NA,NB,NC)
      REAL*8 A(1),B(1),C(1)
      Integer :: NA,NB,NC,NAMAX,NBMAX,I,J,K,
     &           NAPTR,NBPTR,NCPTR
      NAMAX=NA*NC
      NBMAX=NA*NB
      DO I=1,NA
        J=-NB
        DO NAPTR=I,NAMAX,NA
          J=J+NB
          K=0
          A(NAPTR)=0D0
          DO NBPTR=I,NBMAX,NA
            K=K+1
            NCPTR=J+K
            A(NAPTR)=A(NAPTR)+B(NBPTR)*C(NCPTR)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE DMPROD

C
      SUBROUTINE DMPRDB(A,B,C,NA,NB,NC)
      REAL*8 A(1),B(1),C(1),BIK
      Integer :: NA,NB,NC,NAMAX,NBMAX,
     &           NAPTR,NBPTR,NCPTR,I,J,K
      NAMAX=NA*NC
      NBMAX=NA*NB
      DO 1 NAPTR=1,NAMAX
    1 A(NAPTR)=0D0
      DO 3 I=1,NA
      K=0
      DO 3 NBPTR=I,NBMAX,NA
      K=K+1
      BIK=B(NBPTR)
      IF (BIK.EQ.0D0) GO TO 3
      NCPTR=K
      DO 2 NAPTR=I,NAMAX,NA
      A(NAPTR)=A(NAPTR)+BIK*C(NCPTR)
    2 NCPTR=NCPTR+NB
    3 CONTINUE
      RETURN
      END SUBROUTINE DMPRDB
C
      SUBROUTINE RMPRDB(A,B,C,NA,NB,NC)
      SREAL A(1),B(1),C(1),BIK
      Integer :: NA,NB,NC,NAMAX,NBMAX,
     &           NAPTR,NBPTR,NCPTR,I,J,K
      NAMAX=NA*NC
      NBMAX=NA*NB
      DO 1 NAPTR=1,NAMAX
    1 A(NAPTR)=0D0
      DO 3 I=1,NA
      K=0
      DO 3 NBPTR=I,NBMAX,NA
      K=K+1
      BIK=B(NBPTR)
      IF (BIK.EQ.0D0) GO TO 3
      NCPTR=K
      DO 2 NAPTR=I,NAMAX,NA
      A(NAPTR)=A(NAPTR)+BIK*C(NCPTR)
    2 NCPTR=NCPTR+NB
    3 CONTINUE
      RETURN
      END SUBROUTINE RMPRDB
C
      SUBROUTINE DMPRDC(A,B,C,NA,NB,NC)
      REAL*8 A(1),B(1),C(1),CJK
      Integer :: NA,NB,NC,NAM1,NAMAX,NCMAX,
     &           NAPTR,NBPTR,NCPTR,I,J,K,II
      NAM1=NA-1
      NAMAX=NA*NC
      NCMAX=NB*NC
      DO 1 NAPTR=1,NAMAX
    1 A(NAPTR)=0D0
      K=-NA
      DO 3 J=1,NB
      K=K+NA
      I=-NAM1
      DO 3 NCPTR=J,NCMAX,NB
      I=I+NA
      CJK=C(NCPTR)
      IF (CJK.EQ.0D0) GO TO 3
      II=I+NAM1
      NBPTR=K
      DO 2 NAPTR=I,II
      NBPTR=NBPTR+1
    2 A(NAPTR)=A(NAPTR)+B(NBPTR)*CJK
    3 CONTINUE
      RETURN
      END SUBROUTINE DMPRDC
C
      SUBROUTINE RMPRDC(A,B,C,NA,NB,NC)
      SREAL A(1),B(1),C(1),CJK
      Integer :: NA,NB,NC,NAM1,NAMAX,NCMAX,
     &           NAPTR,NBPTR,NCPTR,I,J,K,II
      NAM1=NA-1
      NAMAX=NA*NC
      NCMAX=NB*NC
      DO 1 NAPTR=1,NAMAX
    1 A(NAPTR)=0D0
      K=-NA
      DO 3 J=1,NB
      K=K+NA
      I=-NAM1
      DO 3 NCPTR=J,NCMAX,NB
      I=I+NA
      CJK=C(NCPTR)
      IF (CJK.EQ.0D0) GO TO 3
      II=I+NAM1
      NBPTR=K
      DO 2 NAPTR=I,II
      NBPTR=NBPTR+1
    2 A(NAPTR)=A(NAPTR)+B(NBPTR)*CJK
    3 CONTINUE
      RETURN
      END SUBROUTINE RMPRDC
C
      SUBROUTINE RMPAC(A,B,C,NA,NB,NC)
      SREAL A(1),B(1),C(1),CJK
      Integer :: NA,NB,NC,NAM1,NAMAX,NCMAX,
     &           NAPTR,NBPTR,NCPTR,I,J,K,II
      NAM1=NA-1
      NAMAX=NA*NC
      NCMAX=NB*NC
      K=-NA
      DO 3 J=1,NB
      K=K+NA
      I=-NAM1
      DO 3 NCPTR=J,NCMAX,NB
      I=I+NA
      CJK=C(NCPTR)
      IF (CJK.EQ.0D0) GO TO 3
      II=I+NAM1
      NBPTR=K
      DO 2 NAPTR=I,II
      NBPTR=NBPTR+1
    2 A(NAPTR)=A(NAPTR)+B(NBPTR)*CJK
    3 CONTINUE
      RETURN
      END SUBROUTINE RMPAC
C
      SUBROUTINE DTRPOS(A,B,NA,NB)
      REAL*8 A(NA,NB),B(NB,NA)
      Integer :: NA,NB,I,J
      DO I=1,NA
        DO J=1,NB
          A(I,J)=B(J,I)
	END DO
      END DO
      RETURN
      END SUBROUTINE DTRPOS
C
      SUBROUTINE DZERO(A,N)
      REAL*8 A(1)
      Integer :: N
      A(1:N)=0D0
      RETURN
      END SUBROUTINE DZERO
C
      SUBROUTINE IZERO(A,N)
      INTEGER A(1)
      Integer :: N
      A(1:N)=0
      RETURN
      END SUBROUTINE IZERO
C
      SUBROUTINE LZERO(A,N)
      LOGICAL A(1)
      Integer :: N
      A(1:N)=.FALSE.
      RETURN
      END SUBROUTINE LZERO
C
      SUBROUTINE DOUTER(X,Y,Z)
      REAL*8 X(3,3),Y(3),Z(3)
      Integer :: I,J
      DO I=1,3
        DO J=1,3
          X(I,J)=Y(I)*Z(J)
        END DO
      END DO
      RETURN
      END SUBROUTINE DOUTER
C
      SUBROUTINE DREFLECT(R,N)
      REAL*8 R(3,3),N(3),DMAG
      DMAG=DSQRT(N(1)*N(1)+N(2)*N(2)+N(3)*N(3))
      N(1)=N(1)/DMAG
      N(2)=N(2)/DMAG
      N(3)=N(3)/DMAG
      R(1,1)=-2D0*N(1)*N(1)+1D0
      R(1,2)=-2D0*N(1)*N(2)
      R(1,3)=-2D0*N(1)*N(3)
      R(2,1)=R(1,2)
      R(2,2)=-2D0*N(2)*N(2)+1D0
      R(2,3)=-2D0*N(2)*N(3)
      R(3,1)=R(1,3)
      R(3,2)=R(2,3)
      R(3,3)=-2D0*N(3)*N(3)+1D0
      RETURN
      END SUBROUTINE DREFLECT
C
      SUBROUTINE DTSP2(A,B,NAROW,NACOL,NACAT,NECRD)
      SREAL A(NAROW,NACOL,NACAT)
      REAL*8 B(NAROW,NACOL,NACAT)
      Integer :: NAROW,NACOL,NACAT,NECRD,I,J
      DO 1 J=1,NECRD
        DO 1 I=1,NAROW
    1     A(I,J,NACAT)=B(I,J,NACAT)
      RETURN
      END SUBROUTINE DTSP2
C

      PURE SUBROUTINE DORTHOGANALIZE(X,Y,Z)
      REAL*8, intent(inout) :: X(3),Y(3),Z(3)
      REAL*8 :: D

      ! Normalize X
      D=SQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      X(1:3)=X(1:3)/D

      ! Z = X cross-product Y
      Z(1)=-X(3)*Y(2)+X(2)*Y(3)
      Z(2)= X(3)*Y(1)-X(1)*Y(3)
      Z(3)=-X(2)*Y(1)+X(1)*Y(2)

      ! normalize Z
      D=SQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))
      Z(1:3)=Z(1:3)/D

      ! Y = Z cross-product X
      Y(1)=-Z(3)*X(2)+Z(2)*X(3)
      Y(2)= Z(3)*X(1)-Z(1)*X(3)
      Y(3)=-Z(2)*X(1)+Z(1)*X(2)
      RETURN
      END SUBROUTINE DORTHOGANALIZE
C

C
C This routine is based on Eular's rotation theorem,
C see Spacecraft Attitude Dynamics (Peter C. Hughes), pg 10
C The subroutine argument 'th' is the Eular rotation input vector,
C with its direction defining the rotation axis, and its magnitude
C defining the amount of rotation. -jzlou
C
	SUBROUTINE Qform(Q,dQ,th)
	INTEGER i,j
	REAL*8 Q(3,3),dQ(3,3),th(3),qmag,qhat(3),cq,sq,omcq

	qmag=DSQRT(th(1)**2+th(2)**2+th(3)**2) ! amount (angle) of rotation

	IF (qmag.GT.1d-18) THEN
	  ! Here defines the 3x3 Eular rotation matrix
	  qhat(1:3)=th(1:3)/qmag  ! qhat = normalized theta angle
	  cq=COS(qmag)  ! cos(theta)
	  sq=SIN(qmag)  ! sin(theta)
	  omcq=1d0-cq  ! 1 - cos(theta)
	  Q(1,1)=omcq*qhat(1)**2+cq
	  Q(1,2)=omcq*qhat(1)*qhat(2)-sq*qhat(3)
	  Q(1,3)=omcq*qhat(1)*qhat(3)+sq*qhat(2)
	  Q(2,1)=omcq*qhat(2)*qhat(1)+sq*qhat(3)
	  Q(2,2)=omcq*qhat(2)**2+cq
	  Q(2,3)=omcq*qhat(2)*qhat(3)-sq*qhat(1)
	  Q(3,1)=omcq*qhat(3)*qhat(1)-sq*qhat(2)
	  Q(3,2)=omcq*qhat(3)*qhat(2)+sq*qhat(1)
	  Q(3,3)=omcq*qhat(3)**2+cq
	ELSE
	  Q(1,1)=1d0; Q(1,2)=0d0; Q(1,3)=0d0
	  Q(2,1)=0d0; Q(2,2)=1d0; Q(2,3)=0d0
	  Q(3,1)=0d0; Q(3,2)=0d0; Q(3,3)=1d0
	END IF
	dQ = Q
	dQ(1,1)=dQ(1,1)-1d0
	dQ(2,2)=dQ(2,2)-1d0
	dQ(3,3)=dQ(3,3)-1d0
	RETURN
	END SUBROUTINE Qform

C
      SUBROUTINE DCOMPACT(A,B,NAROW,NACOL,NBROW,NBCOL)
      REAL*8 A(NAROW,NACOL)
      REAL*8 B(NBROW,NBCOL)
      Integer :: NAROW,NACOL,NBROW,NBCOL,I,J
      DO 1 J=1,NACOL
      DO 1 I=1,NAROW
    1 A(I,J)=B(I,J)
      RETURN
      END SUBROUTINE DCOMPACT
C
      SUBROUTINE COMPACT(A,B,NAROW,NACOL,NBROW,NBCOL)
      SREAL A(NAROW,NACOL)
      REAL*8 B(NBROW,NBCOL)
      Integer :: NAROW,NACOL,NBROW,NBCOL,I,J
      DO 1 J=1,NACOL
      DO 1 I=1,NAROW
    1 A(I,J)=B(I,J)
      RETURN
      END SUBROUTINE COMPACT
C
	SUBROUTINE DPROJECT(P,N)
	REAL*8 P(3,3),N(3),DMAG
        DMAG=DSQRT(N(1)*N(1)+N(2)*N(2)+N(3)*N(3))
        N(1)=N(1)/DMAG
        N(2)=N(2)/DMAG
        N(3)=N(3)/DMAG
	P(1,1)=-N(1)*N(1)+1D0
	P(1,2)=-N(1)*N(2)
	P(1,3)=-N(1)*N(3)
	P(2,1)=P(1,2)
	P(2,2)=-N(2)*N(2)+1D0
	P(2,3)=-N(2)*N(3)
	P(3,1)=P(1,3)
	P(3,2)=P(2,3)
	P(3,3)=-N(3)*N(3)+1D0
	RETURN
	END SUBROUTINE DPROJECT

C***********************************************************************

	SUBROUTINE DXMAT(X,Y)
	REAL*8 X(3,3),Y(3)
	X(1,1)=0D0
	X(1,2)=-Y(3)
	X(1,3)= Y(2)
	X(2,1)= Y(3)
	X(2,2)=0D0
	X(2,3)=-Y(1)
	X(3,1)=-Y(2)
	X(3,2)= Y(1)
	X(3,3)=0D0
	RETURN
	END SUBROUTINE DXMAT
C
	SUBROUTINE DXPROD(X,Y,Z)
	REAL*8 X(3),Y(3),Z(3)
	X(1)=-Y(3)*Z(2)+Y(2)*Z(3)
	X(2)= Y(3)*Z(1)-Y(1)*Z(3)
	X(3)=-Y(2)*Z(1)+Y(1)*Z(2)
	RETURN
	END SUBROUTINE DXPROD

C***********************************************************************

      FUNCTION IDFLOOR(X)
      IMPLICIT NONE
      INTEGER IDFLOOR
      REAL*8 X
      IF (X.GE.0D0) THEN
        IDFLOOR=INT(X)
      ELSE
        IDFLOOR=-INT(1E0-X)
      END IF
      RETURN
      END FUNCTION IDFLOOR

C***********************************************************************

      SUBROUTINE INTNORM(X,Y,U00,U01,U10,U11,U,DUDX,DUDY)
      IMPLICIT NONE
      INTEGER IDFLOOR
      REAL*8 X,Y,U00,U01,U10,U11,U,DUDX,DUDY,C1,C2,C3
      C1=U10-U00
      C2=U11-U10-U01+U00
      C3=U01-U00
      DUDX=C2*Y+C1
      DUDY=C2*X+C3
      U=DUDX*X+C3*Y+U00
      RETURN
      END SUBROUTINE INTNORM

C***********************************************************************

      SUBROUTINE INTCUBIC(srfIdx,eltIdx,
     &                    x,y,GridMat,mGridMat,nGridMat,mGridSrf,
     &                    i0,j0,fh,dfdx,dfdy)
      use elt_mod, only : GridSrfdx

      IMPLICIT NONE
      Integer :: srfIdx,eltIdx,mGridMat,nGridMat,mGridSrf,i0,j0,i
      REAL*8 :: GridMat(nGridMat,nGridMat)
      REAL*8 :: x,y,dx,xstrt,fh,dfdx,dfdy
      Real*8 :: GridDat(nGridMat*nGridMat)
      REAL*8, dimension(:,:), pointer, save :: xa,ya
      Logical, dimension(:), pointer, save :: splineInit
      Logical, save :: first_entry=.true.

      If (first_entry) Then
        allocate(splineInit(mGridSrf)); splineInit=.false.
        allocate(xa(nGridMat,mGridSrf),ya(nGridMat,mGridSrf))
        first_entry=.false.
      End If

      ! GridMat(:,:): (in Fortran) first dimension is x, second dimension is y.
      ! When GridMat is passed to C, first dimension will become second dimension,
      ! and in C the second dimension is x, so everything works out fine.
#if 0
      Open(unit=12,file='cubic_GridMat.txt',status='replace')
      Do i=1,nGridMat
        Write(12,*) GridMat(1:nGridMat,i)
      End Do
      Close(12)
      Stop
#endif

      If (.not.splineInit(srfIdx)) Then
	! Initialize spline for surface 'srfIdx'
        dx=GridSrfdx(eltIdx)
        !xstrt=-dx*(nGridMat-1d0)/2d0
        xstrt=0d0
        Do i=1,nGridMat
          xa(i,srfIdx)=xstrt+(i-1)*dx; ya(i,srfIdx)=xa(i,srfIdx)
        End Do
        Call SPLIE2F2C(ya(:,srfIdx),xa(:,srfIdx),GridMat,
     &                 nGridMat,nGridMat,srfIdx)
        splineInit(srfIdx)=.true.
      End If

      ! Compute cubic spline of GridMat at x,y
      Call SPLIN2F2C(ya(:,srfIdx),xa(:,srfIdx),nGridMat,nGridMat,y,x,
     &               fh,dfdy,dfdx,j0,i0,srfIdx)
      END SUBROUTINE INTCUBIC

C***********************************************************************

      DOUBLE PRECISION FUNCTION DDOTLP(N,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      DDOTLP = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOTLP = DTEMP
      RETURN
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOTLP = DTEMP
      RETURN
      END FUNCTION DDOTLP

C***********************************************************************

      SUBROUTINE CAXPY(N,DA,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END SUBROUTINE CAXPY

C***********************************************************************

      SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      DOUBLE PRECISION X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     &        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      DOUBLE PRECISION T
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,SCALE,SHIFT,SL,SM,SN,
     &                 SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
      MAXIT = 30
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
            S(L) = CNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0D0) GO TO 10
               IF (X(L,L) .NE. 0.0D0) S(L) = DSIGN(S(L),X(L,L))
               CALL CSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0D0) GO TO 30
               T = -DDOTLP(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL CAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
            E(L) = CNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0D0) GO TO 80
               IF (E(LP1) .NE. 0.0D0) E(L) = DSIGN(E(L),E(LP1))
               CALL CSCAL(P-L,1.0D0/E(L),E(LP1),1)
               E(LP1) = 1.0D0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
               DO 90 I = LP1, N
                  WORK(I) = 0.0D0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL CAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL CAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0D0
  180       CONTINUE
            U(J,J) = 1.0D0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -DDOTLP(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL CAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL CSCAL(N-L+1,-1.0D0,U(L,L),1)
               U(L,L) = 1.0D0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0D0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0D0
  260          CONTINUE
               U(L,L) = 1.0D0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO 310 J = LP1, P
                  T = -DDOTLP(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL CAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0D0
  330       CONTINUE
            V(L,L) = 1.0D0
  340    CONTINUE
  350 CONTINUE
      MM = M
      ITER = 0
  360 CONTINUE
         IF (M .EQ. 0) GO TO 620
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
            GO TO 620
  370    CONTINUE
         DO 390 LL = 1, M
            L = M - LL
            IF (L .EQ. 0) GO TO 400
            TEST = DABS(S(L)) + DABS(S(L+1))
            ZTEST = TEST + DABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0D0
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + DABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + DABS(E(LS-1))
               ZTEST = TEST + DABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0D0
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
         GO TO (490,520,540,570), KASE
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO 530 K = L, M
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
  540    CONTINUE
            SCALE = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)),
     *                    DABS(S(L)),DABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550
               SHIFT = DSQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     *            CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
  570    CONTINUE
            IF (S(L) .GE. 0.0D0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL CSCAL(P,-1.0D0,V(1,L),1)
  580       CONTINUE
  590       IF (L .EQ. MM) GO TO 600
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     *            CALL CSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     *            CALL CSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END SUBROUTINE DSVDC

C***********************************************************************

      SUBROUTINE CSCAL(N,DA,DX,INCX)
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END SUBROUTINE CSCAL

C***********************************************************************

      SUBROUTINE DROTG(DA,DB,C,S)
      DOUBLE PRECISION DA,DB,C,S,ROE,SCALE,R,Z
      ROE = DB
      IF( DABS(DA) .GT. DABS(DB) ) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF( SCALE .NE. 0.0D0 ) GO TO 10
         C = 1.0D0
         S = 0.0D0
         R = 0.0D0
         GO TO 20
   10 R = SCALE*DSQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = DSIGN(1.0D0,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = 1.0D0
      IF( DABS(DA) .GT. DABS(DB) ) Z = S
      IF( DABS(DB) .GE. DABS(DA) .AND. C .NE. 0.0D0 ) Z = 1.0D0/C
      DA = R
      DB = Z
      RETURN
      END SUBROUTINE DROTG

C***********************************************************************

      SUBROUTINE QSORT(Arr,N)
      Implicit None
      Real*8, dimension(:) :: Arr
      Integer :: N,L,JStack,IR,IQ,I,J
      Real*8 :: FX,A
      Real*8, parameter :: FM=7875,FA=211,FC=1663,FMI=1d0/FM
      Integer, parameter :: M=7,NStack=10000
      Integer :: IStack(NStack)

      JStack=0; L=1; IR=N; FX=0d0
 10   If (IR-L<M) Then
        ! Sort by straight insertion
        Do J=L+1,IR
          A=Arr(J)
          Do I=J-1,1,-1
            If (Arr(I).LE.A) GOTO 12
            Arr(I+1)=Arr(I)
          End Do
          I=0
 12       Arr(I+1)=A
        End Do
        If (JStack.EQ.0) Return
        IR=IStack(JStack)
        L=IStack(JStack-1); JStack=JStack-2

      Else
        I=L; J=IR
        FX=MOD(FX*FA+FC,FM)
        IQ=L+(IR-L+1)*(FX*FMI)
        A=Arr(IQ)
        Arr(IQ)=Arr(L)
 20     Continue
 21     If (J.GT.0) Then
          If (A.LT.Arr(J)) Then
            J=J-1; GOTO 21
          End If
        End If
        If (J.LE.I) Then
          Arr(I)=A
          GOTo 30
        End If
        Arr(I)=Arr(J); I=I+1
 22     If (I.LE.N) Then
          If (A.GT.Arr(I)) Then
            I=I+1
            GOTO 22
          End If
        End If
        If (J.LE.I) Then
          Arr(J)=A; I=J
          GOTO 30
        End If
        Arr(J)=Arr(I); J=J-1
        GOTO 20
 30     JStack=JStack+2
        If (JStack.GT.NStack) Pause 'NStack must be made larger'
        If (IR-I.GE.I-L) Then
          IStack(JStack)=IR
          IStack(JStack-1)=I+1
          IR=I-1
        Else
          IStack(JStack)=I-1; IStack(JStack-1)=L
          L=I+1
        End If
      End If
      GOTO 10
      END SUBROUTINE QSORT


C***********************************************************************

      SUBROUTINE DROT (N,DX,INCX,DY,INCY,C,S)
      DOUBLE PRECISION DX(1),DY(1),DTEMP,C,S
      INTEGER I,INCX,INCY,IX,IY,N
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END SUBROUTINE DROT

C***********************************************************************

      DOUBLE PRECISION FUNCTION CNRM2( N, DX, INCX)
      INTEGER  NEXT,N,INCX,I,J,NN
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
      IF(N .GT. 0) GO TO 10
         CNRM2  = ZERO
         GO TO 300
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 70 TO NEXT
      GO TO 105
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
   75 SUM = (SUM * XMAX) * XMAX
   85 HITEST = CUTHI/FLOAT( N )
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      CNRM2 = DSQRT( SUM )
      GO TO 300
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
      CNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END FUNCTION CNRM2

C***********************************************************************

      SUBROUTINE DORVEC(zhat,zmag,xvec,yvec)
      IMPLICIT NONE
        REAL*8 zhat(3),zmag,xvec(3),yvec(3),zvec(3),S
	zvec(1)=-xvec(3)*yvec(2)+xvec(2)*yvec(3)
	zvec(2)= xvec(3)*yvec(1)-xvec(1)*yvec(3)
	zvec(3)=-xvec(2)*yvec(1)+xvec(1)*yvec(2)
	zmag=DSQRT(zvec(1)*zvec(1)+zvec(2)*zvec(2)+zvec(3)*zvec(3))
	IF (DABS(zmag).LT.1d-10) THEN
	  IF (yvec(3).NE.0d0) THEN
	    zvec(1)=yvec(2)
	    zvec(2)=-yvec(3)
	    zvec(3)=yvec(1)
	  ELSE IF (yvec(1).NE.0d0) THEN
	    zvec(1)=yvec(2)
	    zvec(2)=yvec(3)
	    zvec(3)=-yvec(1)
	  ELSE IF (yvec(2).NE.0d0) THEN
	    zvec(1)=-yvec(2)
	    zvec(2)=yvec(3)
	    zvec(3)=yvec(1)
	  END IF
	  zhat(1)=-yvec(3)*zvec(2)+yvec(2)*zvec(3)
	  zhat(2)= yvec(3)*zvec(1)-yvec(1)*zvec(3)
	  zhat(3)=-yvec(2)*zvec(1)+yvec(1)*zvec(2)
	  S=DSQRT(zhat(1)*zhat(1)+zhat(2)*zhat(2)+zhat(3)*zhat(3))
	  zhat(1)=zhat(1)/S
	  zhat(2)=zhat(2)/S
	  zhat(3)=zhat(3)/S
	ELSE
	  zhat(1)=zvec(1)/zmag
	  zhat(2)=zvec(2)/zmag
	  zhat(3)=zvec(3)/zmag
      END IF
      RETURN
      END SUBROUTINE DORVEC

C***********************************************************************

      SUBROUTINE ROOT4(n,x1,x2,x3,x4,y1,y2,y3,y4,y5)

C  Compute roots of polynomial of order up to 4 using algebraic formulas
C  Matlab version: FUNCTION [x]=root4(y)
C  For 4th order roots, use Laguer (method=1) or Eigensolver (method=2)

        IMPLICIT NONE
        INTEGER n,method
        REAL*8 y1,y2,y3,y4,y5,a,b,c,d,q,r,p2,yr(5),rtr(4),rti(4),
     &  po3,ao4,a2o4,p1,p0,p22,q3,r2,disc,dw(2),dw2(2)
        COMPLEX*16 x1,x2,x3,x4,rad,partA,partB,Y,RR,cR2,c1,c2,
     &  w,w2,cdisc,yc(5),x(4)
        EQUIVALENCE (dw(1),w),(dw2(1),w2)
        DATA dw/-0.5d0,0.8660254037844386d0/
        DATA dw2/-0.5d0,-0.8660254037844386d0/
        DATA method/2/

C  Check order of polynomial
        IF (y1.NE.0d0) THEN
          n=4
        ELSE IF (y2.NE.0d0) THEN
          n=3
        ELSE IF (y3.NE.0d0) THEN
          n=2
        ELSE IF (y4.NE.0d0) THEN
          n=1
        ELSE
          n=0
        END IF

C  Linear polynomial
        IF (n.EQ.1) THEN
           x1=-y5/y4

C  Quadratic polynomial
        ELSE IF (n.EQ.2) THEN
           IF (y4.EQ.0d0) THEN
              cdisc=-y5/y3
              disc=CDSQRT(cdisc)
              x1=disc
              x2=-disc
           ELSE
              cdisc=y4*y4-4d0*y3*y5
              rad=-(y4+SIGN(1d0,y4)*CDSQRT(cdisc))/2d0
              x1=rad/y3
              x2=y5/rad
           END IF

C  Cubic polynomial
        ELSE IF (n.EQ.3) THEN
           p2=y3/y2
           p1=y4/y2
           p0=y5/y2
           p22=p2*p2
           po3=p2/3d0
           q=(3d0*p1-p22)/9d0
           r=((2d0*p22-9d0*p1)*p2+27d0*p0)/54d0
           q3=q*q*q
           r2=r*r
           disc=r2+q3
           IF (disc.NE.0d0) THEN
              cdisc=disc
              rad=CDSQRT(cdisc)
              partA=(r-rad)**0.33333333333333333d0
              partB=(r+rad)**0.33333333333333333d0
              x1=w2*(-partA-partB)
              x2=-partA-w*partB
              x3=-w*partA-partB
              IF (.NOT.(isReal(x1).OR.isReal(x2).OR.isReal(x3)))
     &        THEN
                 x1=x1*w
                 x2=x2*w
                 x3=x3*w
                 IF (.NOT.(isReal(x1).OR.isReal(x2).OR.isReal(x3)))
     &           THEN
                    x1=x1*w
                    x2=x2*w
                    x3=x3*w
                 END IF
              END IF
              x1=x1-po3
              x2=x2-po3
              x3=x3-po3
           ELSE
              partA=r**0.33333333333333333d0
              x1=-2*partA-po3
              x2=x1
              x3=x1
           END IF

C  Quartic polynomial
        ELSE IF (n.EQ.4) THEN

          IF (method.EQ.1) THEN
C Call NR dzroots routine
            yc(1)=dcmplx(y5,0.d0)
            yc(2)=dcmplx(y4,0.d0)
            yc(3)=dcmplx(y3,0.d0)
            yc(4)=dcmplx(y2,0.d0)
            yc(5)=dcmplx(y1,0.d0)
            CALL dzroots(yc,4,x,.TRUE.)
            x1=x(1)
            x2=x(2)
            x3=x(3)
            x4=x(4)
          ELSE IF (method.EQ.2) THEN
            yr(1)=y5
            yr(2)=y4
            yr(3)=y3
            yr(4)=y2
            yr(5)=y1
            CALL zrhqr(yr,4,rtr,rti)
            x1=dcmplx(rtr(1),rti(1))
            x2=dcmplx(rtr(2),rti(2))
            x3=dcmplx(rtr(3),rti(3))
            x4=dcmplx(rtr(4),rti(4))
          ELSE
            WRITE(*,*) ' mathsub.f error:  ROOT4 not called correctly'
          END IF

        END IF

        RETURN
      END SUBROUTINE ROOT4

C***********************************************************************

      SUBROUTINE dzroots(a,m,roots,polish)
C     Needed for ROOT4, method=1
      INTEGER m,MAXM
      REAL*8 EPS
      COMPLEX*16 a(m+1),roots(m)
      LOGICAL polish
      PARAMETER (EPS=1.e-12,MAXM=101)
CU    USES dlaguer
      INTEGER i,j,jj,its
      COMPLEX*16 ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
        x=cmplx(0.d0,0.d0)
        call dlaguer(ad,j,x,its)
        if(dabs(dimag(x)).le.2.*EPS**2*dabs(dreal(x)))
     &    x=dcmplx(dreal(x),0.)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call dlaguer(a,m,roots(j),its)
14      continue
      endif
      do 16 j=2,m
        x=roots(j)
        do 15 i=j-1,1,-1
          if(dreal(roots(i)).le.dreal(x))goto 10
          roots(i+1)=roots(i)
15      continue
        i=0
10      roots(i+1)=x
16    continue
      return
      END SUBROUTINE dzroots

C***********************************************************************

      SUBROUTINE dlaguer(a,m,x,its)
C     Needed for ROOT4, method=1
      INTEGER m,its,MAXIT,MR,MT
      REAL*8 EPSS
      COMPLEX*16 a(m+1),x
      PARAMETER (EPSS=2.e-14,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL*8 abx,abp,abm,err,frac(MR)
      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=cdabs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=cdabs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=cdabs(b)+abx*err
11      continue
        err=EPSS*err
        if(cdabs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.d0*f/b
          sq=cdsqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=cdabs(gp)
          abm=cdabs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.d0) then
            dx=m/gp
          else
            dx=cdexp(dcmplx(dlog(1.d0+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END SUBROUTINE dlaguer

C***********************************************************************

      SUBROUTINE zrhqr(a,m,rtr,rti)
C     Eigenvalue root solver
C     Needed for ROOT4, method=2
      INTEGER m,MAXM
      REAL*8 a(m+1),rtr(m),rti(m)
      PARAMETER (MAXM=50)
CU    USES balanc,hqr
      INTEGER j,k
      REAL*8 hess(MAXM,MAXM),xr,xi
      if (m.gt.MAXM.or.a(m+1).eq.0.) pause 'bad args in zrhqr'
      do 12 k=1,m
        hess(1,k)=-a(m+1-k)/a(m+1)
        do 11 j=2,m
          hess(j,k)=0.d0
11      continue
        if (k.ne.m) hess(k+1,k)=1.d0
12    continue
      call balanc(hess,m,MAXM)
      call hqr(hess,m,MAXM,rtr,rti)
      do 14 j=2,m
        xr=rtr(j)
        xi=rti(j)
        do 13 k=j-1,1,-1
          if(rtr(k).le.xr)goto 1
          rtr(k+1)=rtr(k)
          rti(k+1)=rti(k)
13      continue
        k=0
1       rtr(k+1)=xr
        rti(k+1)=xi
14    continue
      return
      END SUBROUTINE zrhqr

C***********************************************************************

      SUBROUTINE balanc(a,n,np)
C     Needed for ROOT4, method=2
      INTEGER n,np
      REAL*8 a(np,np),RADIX,SQRDX
      PARAMETER (RADIX=2.d0,SQRDX=RADIX**2)
      INTEGER i,j,last
      REAL*8 c,f,g,r,s
1     continue
        last=1
        do 14 i=1,n
          c=0.d0
          r=0.d0
          do 11 j=1,n
            if(j.ne.i)then
              c=c+dabs(a(j,i))
              r=r+dabs(a(i,j))
            endif
11        continue
          if(c.ne.0.d0.and.r.ne.0.d0)then
            g=r/RADIX
            f=1.d0
            s=c+r
2           if(c.lt.g)then
              f=f*RADIX
              c=c*SQRDX
            goto 2
            endif
            g=r*RADIX
3           if(c.gt.g)then
              f=f/RADIX
              c=c/SQRDX
            goto 3
            endif
            if((c+r)/f.lt.0.95d0*s)then
              last=0
              g=1.d0/f
              do 12 j=1,n
                a(i,j)=a(i,j)*g
12            continue
              do 13 j=1,n
                a(j,i)=a(j,i)*f
13            continue
            endif
          endif
14      continue
      if(last.eq.0)goto 1
      return
      END SUBROUTINE balanc

C***********************************************************************

      SUBROUTINE hqr(a,n,np,wr,wi)
C     Needed for ROOT4, method=2
      INTEGER n,np
      REAL*8 a(np,np),wi(np),wr(np)
      INTEGER i,its,j,k,l,m,nn
      REAL*8 anorm,p,q,r,s,t,u,v,w,x,y,z
      anorm=0.
      do 12 i=1,n
        do 11 j=max(i-1,1),n
          anorm=anorm+dabs(a(i,j))
11      continue
12    continue
      nn=n
      t=0.d0
1     if(nn.ge.1)then
        its=0
2       do 13 l=nn,2,-1
          s=dabs(a(l-1,l-1))+dabs(a(l,l))
          if(s.eq.0.d0)s=anorm
          if(dabs(a(l,l-1))+s.eq.s)goto 3
13      continue
        l=1
3       x=a(nn,nn)
        if(l.eq.nn)then
          wr(nn)=x+t
          wi(nn)=0.d0
          nn=nn-1
        else
          y=a(nn-1,nn-1)
          w=a(nn,nn-1)*a(nn-1,nn)
          if(l.eq.nn-1)then
            p=0.5d0*(y-x)
            q=p**2+w
            z=dsqrt(dabs(q))
            x=x+t
            if(q.ge.0.d0)then
              z=p+sign(z,p)
              wr(nn)=x+z
              wr(nn-1)=wr(nn)
              if(z.ne.0.d0)wr(nn)=x-w/z
              wi(nn)=0.d0
              wi(nn-1)=0.d0
            else
              wr(nn)=x+p
              wr(nn-1)=wr(nn)
              wi(nn)=z
              wi(nn-1)=-z
            endif
            nn=nn-2
          else
	    if(its.eq.30)pause 'too many iterations in hqr'
            if(its.eq.10.or.its.eq.20)then
              t=t+x
              do 14 i=1,nn
                a(i,i)=a(i,i)-x
14            continue
              s=dabs(a(nn,nn-1))+dabs(a(nn-1,nn-2))
              x=0.75d0*s
              y=x
              w=-0.4375d0*s**2
            endif
            its=its+1
            do 15 m=nn-2,l,-1
              z=a(m,m)
              r=x-z
              s=y-z
              p=(r*s-w)/a(m+1,m)+a(m,m+1)
              q=a(m+1,m+1)-z-r-s
              r=a(m+2,m+1)
              s=dabs(p)+dabs(q)+dabs(r)
              p=p/s
              q=q/s
              r=r/s
              if(m.eq.l)goto 4
              u=dabs(a(m,m-1))*(dabs(q)+dabs(r))
              v=dabs(p)*(dabs(a(m-1,m-1))+dabs(z)+dabs(a(m+1,m+1)))
              if(u+v.eq.v)goto 4
15          continue
4           do 16 i=m+2,nn
              a(i,i-2)=0.d0
              if (i.ne.m+2) a(i,i-3)=0.d0
16          continue
            do 19 k=m,nn-1
              if(k.ne.m)then
                p=a(k,k-1)
                q=a(k+1,k-1)
                r=0.d0
                if(k.ne.nn-1)r=a(k+2,k-1)
                x=dabs(p)+dabs(q)+dabs(r)
                if(x.ne.0.d0)then
                  p=p/x
                  q=q/x
                  r=r/x
                endif
              endif
	      s=sign(dsqrt(p**2+q**2+r**2),p)
              if(s.ne.0.d0)then
                if(k.eq.m)then
                  if(l.ne.m)a(k,k-1)=-a(k,k-1)
                else
                  a(k,k-1)=-s*x
                endif
                p=p+s
                x=p/s
                y=q/s
                z=r/s
                q=q/p
                r=r/p
                do 17 j=k,nn
                  p=a(k,j)+q*a(k+1,j)
                  if(k.ne.nn-1)then
                    p=p+r*a(k+2,j)
                    a(k+2,j)=a(k+2,j)-p*z
                  endif
                  a(k+1,j)=a(k+1,j)-p*y
                  a(k,j)=a(k,j)-p*x
17              continue
                do 18 i=l,min(nn,k+3)
                  p=x*a(i,k)+y*a(i,k+1)
                  if(k.ne.nn-1)then
                    p=p+z*a(i,k+2)
                    a(i,k+2)=a(i,k+2)-p*r
                  endif
                  a(i,k+1)=a(i,k+1)-p*q
                  a(i,k)=a(i,k)-p
18              continue
              endif
19          continue
            goto 2
          endif
        endif
      goto 1
      endif
      return
      END SUBROUTINE hqr

C***********************************************************************

	FUNCTION isReal(x)
	IMPLICIT NONE
	LOGICAL isReal
	REAL*8 ix,rx,ratio
	COMPLEX*16 x
	rx=DBLE(x)
	ix=DIMAG(x)
	IF (ix.EQ.0d0) THEN
	   isReal=.TRUE.
	ELSE
	   ratio=ABS(rx/ix)
	   IF (ratio.GT.1d9) THEN
	      isReal=.TRUE.
	   ELSE
	      isReal=.FALSE.
	   END IF
	END IF
	RETURN
	END FUNCTION isReal

C***********************************************************************

C  poisson - poisson random number generator
C
	REAL*8 FUNCTION POISSON(theta, nts)
	IMPLICIT NONE

C  value returned is a pseudo-random number with Poisson distribution.
C     theta    Poisson PARAMETER, both mean and variance
C     nts      INTEGER*4 temporary storage
C              (user must dimension 2)

	INTEGER k, nts(2)
	REAL*8 h, t, theta, u

	u = EXP(theta)*RNG(nts)
	k= 0; t= 1; h= 1

 100	IF (h .GT. u) GO TO 1000
	k= k+1
	t= t*(theta/k)
	h= h+t
	GO TO 100

 1000	CONTINUE

	POISSON=k
	RETURN
	END FUNCTION POISSON

C***********************************************************************
C   rayleigh - rayleigh noise random number generator

	REAL*8 FUNCTION rayleigh(a,b, nts)

	IMPLICIT NONE

C     value returned is a pseudo-random number with rayleigh
C     distribution.

C     b: bound; for large b, we have
C           mean     = a*sqrt(pi/2)   (about 1.25 a)
C           variance = a*sqrt(2-pi/2) (about 0.655 a)
C
C     nts      INTEGER*4 temporary storage
C              (user must dimension 2)
C
	INTEGER nts(2)
	REAL*8 a,b, k, u
	REAL*8 v

	v= 0.5*(b/a)**2
	IF ( v .GT. 36) THEN
	   k= 1.0
	ELSE
	   k= 1.0 - exp(-v)
	END IF

	u = RNG (nts)
	v= k*u
C      rayleigh = a*sqrt( -2.0 * dlog(1.0 - v) )
	rayleigh = a*sqrt( dlog( 1.0/(1.0 - v)**2 ) )

	RETURN
	END FUNCTION rayleigh

C***********************************************************************

	REAL*8 FUNCTION rayleigh2(a, nts)

	IMPLICIT NONE

C     value returned is a pseudo-random number with rayleigh
C     distribution.

C           mean     = a*sqrt(pi/2)   (about 1.25 a)
C           variance = a*sqrt(2-pi/2) (about 0.655 a)
C
C     nts      INTEGER*4 temporary storage
C              (user must dimension 2)
C
	INTEGER nts(2)

	REAL*8 a,b,k,u,v

	u= a*gauss(nts)
	v= a*gauss(nts)

	rayleigh2 = sqrt( u**2 + v**2 )

	RETURN
	END FUNCTION rayleigh2

C***********************************************************************

	SUBROUTINE gauss_vec(x,sigma,n,iseed)

C generate n random states x,
C  zero mean, variance sigma

	IMPLICIT NONE
	INTEGER i,n
	INTEGER iseed(2)
	REAL*8 sigma(n), x(n)

	DO i = 1,n
	   x(i) = sigma(i)*gauss(iseed)
	END DO

	RETURN
	END SUBROUTINE gauss_vec

C***********************************************************************
C   gauss -  gaussian random number generator

	REAL*8 FUNCTION GAUSS (nts)
	IMPLICIT NONE
C
C     value returned is a pseudo-random number with gaussian
C     distribution. its mean is zero and its standard
C     deviation is unity
C
C     formal parameters
C
C     nts      INTEGER*4 temporary storage
C              (user must dimension 2)
C
	INTEGER nts(2)
C
C     DATA typing
C
	REAL*8 u, v, twopi
	DATA twopi / 6.28318530717959 /
C
	u = RNG (nts)
	v = RNG (nts)
	GAUSS = SQRT (-2d0 * LOG(u)) * COS(twopi * v)
C
	RETURN
	END FUNCTION GAUSS


C*****************************************************************
C     RNG   -  uniform random number generator

	REAL*8 FUNCTION RNG (nts)

C   uniform distribution on [0,1]
C   two linear congruence generators are used in a vernam
C   encription connection

C   MUST CALL RNGIN TO INITIALIZE THE SEQUENCE!!

	IMPLICIT NONE
	INTEGER  nts(2)

C     lcg bases
	REAL*8     base1, base2, magic
	PARAMETER  (base1 = 32717.0, base2 = 32713.0, magic = 3125.0)
	INTEGER    ibase1, ibase2
	PARAMETER  (ibase1 = 32717, ibase2 = 32713)

C     local variables
	INTEGER    k
	REAL*8     x

C     check sequence initialization
	IF ((nts(1) .LE. 0) .OR. (nts(1) .GT. ibase1)) THEN
	  nts(1) = 123
	END IF
	IF ((nts(2) .LE. 0) .OR. (nts(2) .GT. ibase2)) THEN
	  nts(2) = 123
	END IF

C     advance the lcgs
	x = FLOAT(nts(1))*magic
	k = INT(x/base1)
	nts(1) = INT(x - DBLE(k)*base1 + 0.5)
	x = FLOAT(nts(2))*magic
	k = INT(x/base2)
	nts(2) = INT(x - DBLE(k)*base2 + 0.5)
	k = ieor (nts(1), nts(2))
	IF (k .EQ. 0) THEN
	  k = 1
	  END IF
	RNG = DBLE(k)/32768.0

	RETURN
	END FUNCTION RNG

C*****************************************************************

C     RNGIN -  uniform random number generator initialization

	SUBROUTINE RNGIN (seed, nts)

C     seed value is transferred to temporary storage
C     seed       INTEGER to start sequence
C     nts        temporary sequence storage (dim 2, INTEGER)

	IMPLICIT NONE
	INTEGER    nts(2), seed

	IF ((seed .GT. 0) .AND. (seed .LT. 32700)) THEN
	  nts(1:2) = seed
	ELSE
	  nts(1) = 345; nts(2) = 567
	END IF
	RETURN
	END SUBROUTINE RNGIN

C***********************************************************************

      FUNCTION POIDEV(xm,idum)
      SREAL poidev,xm,PI
      INTEGER idum
      PARAMETER (PI=3.141592654)
CU    USES gammln,ran1
      SREAL alxm,em,g,oldm,sq,t,y
      SAVE alxm,g,oldm,sq
      DATA oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
        t=t*ran1(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(PI*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
      endif
      poidev=em
      return
      END FUNCTION POIDEV

#if 0
C Below is replaced by sab's updated version above

	FUNCTION POIDEV_OLD(XM,IDUM)
	PARAMETER (PI=3.141592654)
	DATA OLDM /-1./
	IF (XM.LT.12.)THEN
	  IF (XM.NE.OLDM) THEN
	    OLDM=XM
	    G=EXP(-XM)
	  ENDIF
	  EM=-1
	  T=1.
2	 EM=EM+1.
	  T=T*RAN1(IDUM)
	  IF (T.GT.G) GO TO 2
	ELSE
	  IF (XM.NE.OLDM) THEN
	    OLDM=XM
	    SQ=SQRT(2.*XM)
	    ALXM=ALOG(XM)
	    G=XM*ALXM-GAMMLN(XM+1.)
	  ENDIF
1	 Y=TAN(PI*RAN1(IDUM))
	  EM=SQ*Y+XM
	  IF (EM.LT.0.) GO TO 1
	  EM=INT(EM)
	  T=0.9*(1.+Y**2)*EXP(EM*ALXM-GAMMLN(EM+1.)-G)
	  IF (RAN1(IDUM).GT.T) GO TO 1
	ENDIF
	POIDEV=EM
	RETURN
	END FUNCTION POIDEV_OLD
#endif

C***********************************************************************

      FUNCTION gammln(xx)
      SREAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END FUNCTION gammln

#if 0
C Below is replaced by sab's updated version above

	FUNCTION GAMMLN_OLD(XX)
	REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
	DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
	DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
	X=XX-ONE
	TMP=X+FPF
	TMP=(X+HALF)*LOG(TMP)-TMP
	SER=ONE
	DO 11 J=1,6
	  X=X+ONE
	  SER=SER+COF(J)/X
11	CONTINUE
	GAMMLN=TMP+LOG(STP*SER)
	RETURN
	END FUNCTION GAMMLN_OLD
#endif

C***********************************************************************

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      SREAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END FUNCTION ran1

#if 0
C Below is replaced by sab's updated version above

	FUNCTION RAN1_OLD(IDUM)
	DIMENSION R(97)
	PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
	PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
	PARAMETER (M3=243000,IA3=4561,IC3=51349)
	DATA IFF /0/
	IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
	  IFF=1
	  IX1=MOD(IC1-IDUM,M1)
	  IX1=MOD(IA1*IX1+IC1,M1)
	  IX2=MOD(IX1,M2)
	  IX1=MOD(IA1*IX1+IC1,M1)
	  IX3=MOD(IX1,M3)
	  DO 11 J=1,97
	    IX1=MOD(IA1*IX1+IC1,M1)
	    IX2=MOD(IA2*IX2+IC2,M2)
	    R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11	CONTINUE
	  IDUM=1
	ENDIF
	IX1=MOD(IA1*IX1+IC1,M1)
	IX2=MOD(IA2*IX2+IC2,M2)
	IX3=MOD(IA3*IX3+IC3,M3)
	J=1+(97*IX3)/M3
	IF(J.GT.97.OR.J.LT.1)PAUSE
	RAN1=R(J)
	R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
	RETURN
	END FUNCTION RAN1_OLD
#endif

C***********************************************************************

C   Compute image centroid - added by jzlou
        Subroutine Centroid(arr2d,nrow,ncol,xc,yc)
        Real*8, dimension(nrow,ncol) :: arr2d
        Real*8 :: xc,yc,s, mx, my
        Integer :: nrow,ncol, ic,ir

        Call DMATSUM(arr2d,s,nrow,ncol) ! sum of eles of arr2d

        ! Default centroid values
        xc=floor(Real(nrow/2))+1; yc=floor(Real(ncol/2))+1

        If (s/=0) then
          mx=0d0; my=0d0
          Do ic=1,ncol
            ! Sum up values in each column, multiplied by ic,
            ! added up in mx
            my = my + RVECSUM(arr2d(:,ic),nrow)*ic
          End Do
          Do ir=1,nrow
            ! Sum up values in each row, multiplied by ir,
            ! added up in my
            mx = mx + RVECSUM(arr2d(ir,:),ncol)*ir
          End Do
          xc=mx/s; yc=my/s
        End If
        End Subroutine Centroid

C***********************************************************************

C    Calculate x^n - jzlou
        Function PowerF(x,n)
        Real*8 :: x, PowerF
        Integer :: n,i
        PowerF=1d0
        if (n==0 .or. x==0d0) then
          return
        else
          Do i=1,abs(n)
            PowerF=PowerF*x
          End Do
        end if
        if (n<0) PowerF=1d0/PowerF
        End Function PowerF

C***********************************************************************

C    Calculate N!  - jzlou
        Function Factorial(N)
        Integer N,i,Factorial
        Factorial=1
        if (N==0) then
          return
        else
          Do i=2,N
            Factorial=Factorial*i
          End Do
        end if
        End Function Factorial

C***********************************************************************

C Calculate RMS of 'map' using only nonzero elements - jzlou
        Function RMS2(map,m,n)
        Implicit none
        Real*8 , dimension(m,n) :: map, tmp_map
        Real*8 :: s,mean,RMS2
        Integer :: m,n,im,in, tot_nonzeros

        s=0d0; tot_nonzeros=0
        Do im=1,m
          Do in=1,n
            if (map(im,in)/=0) then
              s=s+map(im,in)
              tot_nonzeros=tot_nonzeros+1
            endif
          End Do
        End Do
        if (tot_nonzeros==0) then
          RMS2=-1d0
          return  ! error return
        end if
        mean=s/tot_nonzeros

        tmp_map = map-mean

        s=0d0
        Do im=1,m
          Do in=1,n
            if (map(im,in)/=0) then
              s=s+tmp_map(im,in)*tmp_map(im,in)
            end if
          End Do
        End Do
        s=s/tot_nonzeros; RMS2=sqrt(s)
        End Function RMS2

C***********************************************************************

	Subroutine FindPV(map,m,n,p,v)
        Implicit none
        Real*8 , dimension(m,n) :: map
        Real*8 :: p,v
        Integer :: m,n,im,in

        p=-1d12; v=1d12
        Do im=1,m
          Do in=1,n
            if (p<map(im,in)) p=map(im,in)
            if (v>map(im,in)) v=map(im,in)
          End Do
        End Do
        End Subroutine FindPV

C***********************************************************************

        Function Med(A,sz)
        Implicit none
        Integer :: sz,i,j
        Real*8 :: A(sz),B(sz),tmp,Med

        B=A
#if 0
        Do j=1,sz-1
          Do i=j+1,sz
            If (B(i)>B(j)) Then
              tmp=B(i); B(i)=B(j); B(j)=tmp
            End If
          End Do
        End Do
#endif

        Call QSort(B,sz)
        i=MOD(sz,2)
        If (i==0) Then
          Med=B(sz/2)
        Else
          Med=B(sz/2+1)
        End If
        End Function Med


C***********************************************************************
#if 1
C Bicubic spline routines

C  Bicubic spline driver routine - find interpolated values on a set
C  of untabulated locations defined by xloc and yloc, given tabulated
C  grid defined by x1a and x2a and function values by ya.

	Subroutine bcubic_int(x1a,x2a,ya,xloc,yloc,val)
	Implicit none
	LREAL, DIMENSION(:), INTENT(IN) :: x1a,x2a,xloc,yloc
        LREAL, DIMENSION(:,:), INTENT(IN) :: ya
	LREAL, DIMENSION(:), INTENT(OUT) :: val
	LREAL, DIMENSION(size(x1a),size(x2a)) :: y2a
	Integer :: iv

	! First call splie2 once
	Call splie2(x1a,x2a,ya,y2a)

	Do iv=1,size(val)
	  val(iv)=splin2(x1a,x2a,ya,y2a,yloc(iv),xloc(iv))
	End Do
	End Subroutine bcubic_int


C x1a,x2a provide tabulated locations of grid function values
C stored in ya.

	FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
        USE nn_util, ONLY : assert_eq
        IMPLICIT NONE
        ! Given x1a, x2a, ya as described in splie2 and y2a as produced by
        ! that routine; and given a desired interpolating point x1,x2; this
        ! routine returns an interpolated function value by bicubic spline
        ! interpolation.

        LREAL, DIMENSION(:), INTENT(IN) :: x1a,x2a
        LREAL, DIMENSION(:,:), INTENT(IN) :: ya,y2a
        LREAL, INTENT(IN) :: x1,x2
        LREAL :: splin2
        Integer :: j,m,ndum
        LREAL, DIMENSION(size(x1a)) :: yytmp,y2tmp2

        m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
        ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
        Do j=1,m
          yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
          ! Perform m evaluations of the row splines constructed by splie2,
          ! using the one-dimensional
          ! spline evaluator splint.
        End Do
        Call spline(x1a,yytmp,1.0d30,1.0d30,y2tmp2)
          ! Construct the one-dimensional column spline and evaluate it.
        splin2=splint(x1a,yytmp,y2tmp2,x1)
        END FUNCTION splin2


	SUBROUTINE splie2(x1a,x2a,ya,y2a)
	USE nn_util, ONLY : assert_eq
	IMPLICIT NONE

	! Given an M x N tabulated function ya, and N tabulated independent
	! variables x2a, this routine constructs one-dimensional natural cubic
        ! splines of the rows of ya and returns the second derivatives in the
	! M x N array y2a. (The array x1a is included in the argument
	! list merely for consistency with routine splin2.)
	LREAL, DIMENSION(:), INTENT(IN) :: x1a,x2a
	LREAL, DIMENSION(:,:), INTENT(IN) :: ya
	LREAL, DIMENSION(:,:), INTENT(OUT) :: y2a
	Integer :: j,m,ndum

	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
	Do j=1,m
	  Call spline(x2a,ya(j,:),1.0d30,1.0d30,y2a(j,:))
   	   ! Values 1*10^30 signal a natural spline.
   	End Do
	END SUBROUTINE splie2

	FUNCTION splint(xa,ya,y2a,x)
        USE nn_util, ONLY : assert_eq,nerror
	IMPLICIT NONE
	! Given the arrays xa and ya, which tabulate a function (with
	! the xaxis in increasing or decreasing order), and given the
	! array y2a, which is the output from spline above, and
	! given a value of x, this routine returns a cubic-spline interpolated
        ! value. The arrays xa, ya and y2a are all of the same size.
	LREAL, DIMENSION(:), INTENT(IN) :: xa,ya,y2a

	LREAL, INTENT(IN) :: x
	LREAL :: splint
	Integer :: khi,klo,n
	LREAL :: a,b,h
	n=assert_eq(size(xa),size(ya),size(y2a),'splint')
	klo=max(min(locate(xa,x),n-1),1)
	! We will find the right place in the table by means of locate's
	! bisection algorithm. This is optimal if sequential calls to this
	! routine are at random values of x. If sequential calls are in
	! order, and closely spaced, one would do better to store previous
	! values of klo and khi and test if they remain appropriate on the
	! next call.
	khi=klo+1 ! klo and khi now bracket the input value of x.
	h=xa(khi)-xa(klo)
	If (h == 0.0) call nerror('bad xa input in splint')
	   ! The xaxis must be distinct.
	a=(xa(khi)-x)/h  ! Cubic spline polynomial is now evaluated.
	b=(x-xa(klo))/h
	splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)
     &         +(b**3-b)*y2a(khi))*(h**2)/6.0d0
	END FUNCTION splint


	SUBROUTINE spline(x,y,yp1,ypn,y2)
        USE nn_util, ONLY : assert_eq
	IMPLICIT NONE
	! Given arrays x and y of length N containing a tabulated function,
	! i.e., yi = f(xi), with x1 < x2 < .. . < xN, and given values yp1
	! and ypn for the first derivative of the interpolating
	! function at points 1 and N, respectively, this routine returns an
	! array y2 of length N that contains the second derivatives of the
	! interpolating function at the tabulated points xi. If yp1 and/or ypn
	! are equal to 1 \327 1030 or larger, the routine is signaled to set the
	! corresponding boundary condition for a natural spline, with zero second
	! derivative on that boundary.

	LREAL, DIMENSION(:), INTENT(IN) :: x,y
	LREAL, INTENT(IN) :: yp1,ypn
	LREAL, DIMENSION(:), INTENT(OUT) :: y2
	Integer :: n
	LREAL, DIMENSION(size(x)) :: a,b,c,r
	n=assert_eq(size(x),size(y),size(y2),'spline')
	c(1:n-1)=x(2:n)-x(1:n-1)  ! Set up the tridiagonal equations.
	r(1:n-1)=6.0d0*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0d0*(c(2:n-1)+a(2:n-1))
	b(1)=1.0; b(n)=1.0
	If (yp1 > 0.99d30) Then
	  ! The lower boundary condition is set either to be 'natural'
	  r(1)=0.0; c(1)=0.0
	Else ! or else to have a specified first derivative.
	  r(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	  c(1)=0.5
	End If
        If (ypn > 0.99d30) Then
          ! The upper boundary condition is set either to be 'natural'
	  r(n)=0.0; a(n)=0.0
        Else
          ! Or else to have a specified first derivative.
	  r(n)=(-3.0d0/(x(n)-x(n-1)))
     &         *((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
	  a(n)=0.5
	End If
	Call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	END SUBROUTINE spline


	SUBROUTINE tridag(a,b,c,r,u)
	USE nn_util, ONLY : assert_eq,nerror
	IMPLICIT NONE
	! Solves for a vector u of size N the tridiagonal linear set
	! given by equation (2.4.1) using a serial algorithm. Input vectors b
	! (diagonal elements) and r (right-hand sides) have size N,
	! while a and c (off-diagonal elements) are size N-1.

	LREAL, DIMENSION(:), INTENT(IN) :: a,b,c,r
	LREAL, DIMENSION(:), INTENT(OUT) :: u
	LREAL, DIMENSION(size(b)) :: gam ! One vector of workspace, gam is needed.
	Integer :: n,j
	LREAL :: bet
	n=assert_eq((/size(a)+1,size(b),size(c)+1,
     &              size(r),size(u)/),'tridag')
	bet=b(1)
	If (bet == 0.0) call nerror('tridag: Error at code stage 1')
	! If this happens then you should rewrite your equations as
	! a set of order N-1, with u2 trivially eliminated.
	u(1)=r(1)/bet
	Do j=2,n ! Decomposition and forward substitution.
	  gam(j)=c(j-1)/bet
	  bet=b(j)-a(j-1)*gam(j)
	  if (bet == 0.0) then ! Algorithm fails; see below routine in Vol. 1.
	    call nerror('tridag: Error at code stage 2')
	  end if
	  u(j)=(r(j)-a(j-1)*u(j-1))/bet
	End Do
	Do j=n-1,1,-1  ! Backsubstitution.
	  u(j)=u(j)-gam(j+1)*u(j+1)
	End Do
	END SUBROUTINE tridag


	FUNCTION locate(xx,x)
	! Given an array xx(1:N), and given a value x, returns a value
	! j such that x is between xx(j) and xx(j+1). xx must be monotonic,
	! either increasing or decreasing. j=0 or j=N is returned to indicate
	! that x is out of range.
	Implicit none

	SREAL, dimension(:), intent(in) :: xx
	SREAL, intent(in) :: x
	Integer :: locate
	Integer :: n,jl,jm,ju
	Logical :: ascnd

	n=size(xx)
	ascnd=(xx(n) >= xx(1))
	jl=0; ju=n+1

	Do
	  if (ju-jl<=1) exit
	  jm=(ju+jl)/2
	  if (ascnd .eqv. (x >= xx(jm))) then
	    jl=jm
	  else
	    ju=jm
	  end if
	End Do

	If (x==xx(1)) Then
	  locate=1
	Else If (x==xx(n)) Then
	  locate=n-1
	Else
	  locate=jl
	End If
        END FUNCTION locate
#endif

	Logical Function inTrig(A,B,C,xin,yin,xp,yp)
        Real*8 :: A(3),B(3),C(3),xin,yin,xp,yp,
     &            vin(3),vp(3)
        Integer :: i

        Do i=1,3
          vin(i)=A(i)*xin+B(i)*yin+C(i)
          vp(i)=A(i)*xp+B(i)*yp+C(i)
        End Do
        inTrig=.false.
        If (((vin(1)*vp(1)).GE.0) .AND. ((vin(2)*vp(2)).GE.0)
     &      .AND. ((vin(3)*vp(3)).GE.0)) Then
          inTrig=.true.
        End If
        End Function inTrig

C ********************************************************
        Real*8 Function Calc_RMS2(datArr,ndim)
        Real*8 :: datArr(ndim),mv,t1
        Integer :: ndim,i

        If (ndim==0) Then
          Calc_RMS2 = 0d0
        Else
          mv = sum(datArr)/ndim
          t1 = 0d0
          Do i=1,ndim
            t1 = t1 + (datArr(i)-mv)**2
          End Do
          Calc_RMS2 = dsqrt(t1/ndim)
        End If
        End Function Calc_RMS2

C ********************************************************

        SUBROUTINE  CSWAP (N,DX,INCX,DY,INCY)
          IMPLICIT NONE
          DOUBLE PRECISION DX(1),DY(1),DTEMP
          INTEGER I,INCX,INCY,IX,IY,M,MP1,N
          IF(N.LE.0)RETURN
          IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
          IX = 1
          IY = 1
          IF(INCX.LT.0)IX = (-N+1)*INCX + 1
          IF(INCY.LT.0)IY = (-N+1)*INCY + 1
          DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
          END DO
          RETURN

   20     M = MOD(N,3)
          IF( M .EQ. 0 ) GO TO 40
          DO I = 1,M
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
          END DO
          IF( N .LT. 3 ) RETURN

   40     MP1 = M + 1
          DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I + 1)
            DX(I + 1) = DY(I + 1)
            DY(I + 1) = DTEMP
            DTEMP = DX(I + 2)
            DX(I + 2) = DY(I + 2)
            DY(I + 2) = DTEMP
          END DO
          RETURN
        END SUBROUTINE  CSWAP


        END MODULE math_mod


C***********************************************************************
C****** End file mathsub.F *********************************************
C***********************************************************************


