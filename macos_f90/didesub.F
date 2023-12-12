C***********************************************************************
C	Begin file didesub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

	MODULE didesub

	PRIVATE
	PUBLIC :: IntSrf, dides

	CONTAINS


	SUBROUTINE IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,
     &	  pin,ihat,XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	  wVec,pr,L,pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	use math_mod

	IMPLICIT NONE

	INTEGER NDP,IERROR,IWORK(*),IMODE,IDNDP(31),mDP
	REAL*8 XYZSRF(3,NDP),DZSRF(2,NDP),xIDES,yIDES,wVec(3),
     &	DWORK(*)
	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &  pr(3),pmon(3),xhat(3),yhat(3),zhat(3),L,Nhat(3)
	REAL*8 prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &  Nvec(3),Nmag,dNdp(3,3),Kc,Kr,
     &  k2,itpsi,psitprel,L1,L2,L3,tol,
     &	zc(3),dfdrho(3),d2fdrho2(3,3),
     &	rho(3),psitzc,fh,dzcdrho(3,3),dNvecdp(3,3),
     &	dxdrho(3),dydrho(3),rhoc(3),rhom(3),zmtzc,psitzm

	REAL*8 D1(9),D2(9),S1,S2

	INTEGER i,j

	LOGICAL LROK,ifBuild,ifLNsrf

C	Compute incidence point of ray on the base conic surface

	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
	  rhoc(i)=pmon(i)-pv(i)
   1	CONTINUE
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d10) THEN
	  mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
	  a=1d0+Kc*itpsi*itpsi
	  b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &	  ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
	  c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
	ELSE
	  a=0d0
	  b=itpsi
	  c=psitprel
	END IF

C	Compute intersection point and normal vector

	IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d10)) THEN
	  L=-c/b
	  IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98
	ELSE
	  k2=b*b-4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  IF (b.GE.0d0) THEN
	    ksqrt=-0.5d0*(b+DSQRT(k2))
	  ELSE
	    ksqrt=-0.5d0*(b-DSQRT(k2))
	  END IF
	  kp=ksqrt/a
	  km=c/ksqrt
	  IF(ifLNsrf) THEN
	    IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	      L=kp
	    ELSE
	      L=km
	    END IF
	  ELSE
	    IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	      IF (kp.GT.0D0) THEN
	        L=kp
	      ELSE
	        L=km
	      END IF
	    ELSE
	      IF (km.GT.0D0) THEN
	        L=km
	      ELSE
	        L=kp
	      END IF
	    END IF
	    IF (L.LT.0d0) GO TO 98
	  END IF
	END IF

C	Use L as a starting point, find incidence point on full surface

	L1=0.999d0*L
	L2=1.001d0*L
	CALL DZPB(L1,L2,Kc,Kr,psi,pv,prel,ihat,XYZSRF,DZSRF,
     &	IERROR,IWORK,DWORK,NDP,mDP,wVec,pmon,rhoc,xhat,yhat,zhat,L)
	tol=1d-14
	L3=DZPSolve(L1,L2,tol,Kc,Kr,psi,pv,prel,ihat,XYZSRF,DZSRF,
     &	IERROR,IWORK,DWORK,NDP,mDP,wVec,pmon,rhoc,xhat,yhat,zhat,L)
	L=L3

C	Compute normal vector

	DO 2 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  2	CONTINUE
	xIDES=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	yIDES=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

	IMODE=6
	CALL DIDES(IMODE,NDP,mDP,xIDES,yIDES,wVec,XYZSRF,DZSRF,
     &	IWORK,DWORK,IERROR)
	IF (IERROR.NE.0) CALL DZERO(wVec,3)
	fh=wVec(1)

	DO 5 i=1,3
	  dxdrho(i)=xhat(i)
	  dydrho(i)=yhat(i)
	  dFdrho(i)=wVec(2)*dxdrho(i)+wVec(3)*dydrho(i)
	  zc(i)=rho(i)-fh*zhat(i)
  5	CONTINUE

	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	psitzm=psi(1)*zhat(1)+psi(2)*zhat(2)+psi(3)*zhat(3)
	zmtzc=zhat(1)*zc(1)+zhat(2)*zc(2)+zhat(3)*zc(3)
	S1=Kr+Kc*psitzc
	S2=psitzm*Kr+psitzm*psitzc*Kc+zmtzc

	DO 6 i=1,3
	  pr(i)=pin(i)+L*ihat(i)
	  Nvec(i)=zc(i)+S1*psi(i)-S2*dFdrho(i)
   6	CONTINUE
	Nmag=DSQRT(Nvec(1)*Nvec(1)
     &       +Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  DO 7 i=1,3
	    Nhat(i)=Nvec(i)/Nmag
   7	  CONTINUE
	ELSE
	  DO 8 i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
   8	  CONTINUE
	END IF

C	Compute partials of the normal

	IF (ifBuild) THEN
	  IF (ABS(Kr).LT.1d10) THEN
	    CALL DSMPROD(D1,psi,Kc,3)
	    CALL DOUTER(M,D1,psi)
	    M(1,1)=M(1,1)+1d0
	    M(2,2)=M(2,2)+1d0
	    M(3,3)=M(3,3)+1d0
	  ELSE
	    CALL DZERO(M,9)
	  END IF

CDIDES	  S1=wVec(4,1,1)
CDIDES	  S2=wVec(5,1,1)
CDIDES	  S3=wVec(6,1,1)
	  DO 12 j=1,3
	    DO 11 i=1,3
CDIDES	      d2Fdrho2(i,j)=S1*dxdrho(i)*dxdrho(j)+
CDIDES     &	      S2*(dxdrho(i)*dydrho(j)+dydrho(i)*dxdrho(j))+
CDIDES     &	      S3*dydrho(i)*dydrho(j)
	      d2Fdrho2(i,j)=0d0
	      dzcdrho(i,j)=-zhat(i)*dFdrho(j)
  11	    CONTINUE
	    dzcdrho(j,j)=dzcdrho(j,j)+1d0
  12	  CONTINUE
	  CALL DMPROD(D1,dzcdrho,M,3,3,3)
	  CALL DMPROD(D2,D1,dzcdrho,3,3,3)
	  S1=Kr+Kc*psitzc
	  CALL DSMPROD(D1,d2Fdrho2,S1,9)
	  CALL DSUB(dNvecdp,D2,D1,9)

	  S1=1d0/Nmag
	  CALL DPROJECT(D1,Nhat)
	  CALL DMPROD(D2,D1,dNvecdp,3,3,3)
	  CALL DSMPROD(dNdp,D2,S1,9)
	END IF

	RETURN

C	Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.
	RETURN
	END SUBROUTINE IntSrf

C***********************************************************************

	FUNCTION SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	  XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	  wVec,pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

	IMPLICIT NONE

	INTEGER NDP,IERROR,IWORK(*),IMODE,IDNDP(31),mDP
	REAL*8 XYZSRF(3,NDP),DZSRF(2,NDP),xIDES,yIDES,wVec(3),DWORK(*)
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),xhat(3),
     &	yhat(3),zhat(3),pmon(3),L,SIntSrf
	REAL*8 rho(3),rhom(3),rhoc(3),fh,zc(3),zc2,psitzc
	INTEGER i
cdebug
c	REAL*8 Sval
c	INTEGER icnt
c	DATA icnt/0/
c	icnt=icnt+1
cdebug

C  Compute point-of-incidence coordinates

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	xIDES=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	yIDES=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Compute deformed surface height

	IMODE=6
	CALL DIDES(IMODE,NDP,mDP,xIDES,yIDES,wVec,XYZSRF,DZSRF,
     &	IWORK,DWORK,IERROR)
	IF (IERROR.NE.0) CALL DZERO(wVec,3)
	fh=wVec(1)

C  Compute surface function S

	DO 4 i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
  4	CONTINUE
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SIntSrf=zc2+(Kc*psitzc+2d0*Kr)*psitzc

cdebug
c	SVal=SIntSrf
c 100	FORMAT(1P,' SIntSrf: S=',d12.5,' L=',d12.5,' i=',i4)
c	WRITE(*,100)SIntSrf,L,icnt
cdebug

	RETURN
	END FUNCTION SIntSrf

C***********************************************************************

      FUNCTION DZPSolve(X1,X2,TOL,Kc,Kr,psi,pv,prel,ihat,XYZSRF,DZSRF,
     &	IERROR,IWORK,DWORK,NDP,mDP,wVec,pmon,rhoc,xhat,yhat,zhat,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NDP,IERROR,IWORK(*),mDP
      REAL*8 XYZSRF(3,NDP),DZSRF(2,NDP),wVec(3),DWORK(*)
      REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),L,!SIntSrf,
     &	DZPSolve,pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3)
      PARAMETER (ITMAX=400,EPS=1d-14)
      A=X1
      B=X2

      FA=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	wVec,pmon,rhoc,xhat,yhat,zhat,A)
      FB=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	wVec,pmon,rhoc,xhat,yhat,zhat,B)
      IF (FB*FA.GT.0.) PAUSE 'DZPB failed.'
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A; FC=FA; D=B-A; E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B; B=C; C=A
          FA=FB; FB=FC; FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          DZPSolve=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D; D=P/Q
          ELSE
            D=XM; E=D
          ENDIF
        ELSE
          D=XM; E=D;
        ENDIF
        A=B; FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &		   XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &		   wVec,pmon,rhoc,xhat,yhat,zhat,B)
11    CONTINUE
      PAUSE 'DZPSolve exceeding maximum iterations.'
      DZPSolve=B
      RETURN
      END FUNCTION DZPSolve

C***********************************************************************

	SUBROUTINE DZPB(X1,X2,Kc,Kr,psi,pv,prel,ihat,XYZSRF,DZSRF,
     &	IERROR,IWORK,DWORK,NDP,mDP,wVec,pmon,rhoc,xhat,yhat,zhat,L)
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (gold=1.6,NTRY=50)
	INTEGER NDP,IERROR,IWORK(*),mDP
	REAL*8 XYZSRF(3,NDP),DZSRF(2,NDP),wVec(3),DWORK(*)
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),L,!SIntSrf,
     &	pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3)
	F1=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	wVec,pmon,rhoc,xhat,yhat,zhat,X1)
	F2=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	wVec,pmon,rhoc,xhat,yhat,zhat,X2)
	DO 1 J=1,NTRY
	  IF(F1*F2.LT.0.)RETURN
	  IF(ABS(F1).LT.ABS(F2))THEN
	    X1=X1+gold*(X1-X2)
	    F1=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	    XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	    wVec,pmon,rhoc,xhat,yhat,zhat,X1)
	  ELSE
	    X2=X2+gold*(X2-X1)
	    F2=SIntSrf(Kc,Kr,psi,pv,prel,ihat,
     &	    XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	    wVec,pmon,rhoc,xhat,yhat,zhat,X2)
	  ENDIF
1	CONTINUE
	RETURN
	END SUBROUTINE DZPB

C***********************************************************************
C*******Start DIDES routine library*************************************
C***********************************************************************

      SUBROUTINE DIDES(MD,NDP,NDPMAX,XI,YI,ZI,XYZDAT,DZ,IWK,WRK,IERROR)
C
C     This SUBROUTINE performs smooth surface fitting when the
C     projections of the DATA points in the X-Y plane are
C     irregularly distributed in the plane.
C
C     The basis of the DIDES software is from
C
C     "C0 and C1 Surface Interpolation to Scattered DATA", in Chapter
C     12.4 of the JPL MATH77 package.  These codes are due to Charles
C     L. Lawson.
C
C     This code written by Fred T. Krogh, JPL, July 25, 1996
C
C *********************** Calling Sequence *****************************
C
C DZ     DZ(1, i) gives the value of the partial of the surface with
C        with respect to x at the i-th DATA point, and DZ(2, i) is
C        defined similarly for the partial with respect to y.  This may
C        be input, or may be computed on request.
C IERROR Returns result information.  When MD < 3, this is set
C        as follows.  IF IERROR is not zero, an error message will have
C        been printed.
C   IERROR  = 0 Normal termination, the triangular grid is complete and
C               optimized.
C           = 1 A complete triangular grid has been constructed, but the
C               process of optimizing the shapes of the triangular cells
C               did not complete.  This should not happen.
C           = 2 Either all given points are colinear or ELSE some
C               pair of the points are duplicates.
C           = 3 Some pair of the points are duplicates.
C           = 4 The dimension MB is not large enough.  Setting MB =
C               NDP + 1 will suffice.
C           = 5 The dimension MT is not large enough.  Setting MT =
C               12 *(NDP-1) will suffice.
C       IERROR is set to 0 when a CALL with MD = 3 or 4 is made.
C       When MD = 5 or 6 IERROR is set as follows.  IF IERROR is not 0,
C       an error message is printed and no interpolaton is done.
C   IERROR = 0  The point (XI, YI) is within the convex hull of the
C                triangular grid, or possibly outside by a small
C                tolerance.
C           = 1  The point (XI, YI) is outside the convex hull of the
C                triangular grid by more than the built-in small
C                tolerance.
C           = 2  Error condition.  STGFND has tested the point (XI,YI)
C                against NT triangles without resolving its status.
C                This would indicate an error in the definition of the
C                triangular grid or in the design or coding of the
C                SUBROUTINE.
C        Errors which can occur for any value of MD, all of which PRINT
C        a message are flagged with the following values of IERROR.
C   IERROR = 6  NDP < 3.
C   IERROR = 7  NDP changed from the last CALL.
C   IERROR = 8  MD < 1 or MD > 6.
C IWK    INTEGER work space of size at least 17*NDP - 8.  The first
C        12*NDP - 12 are used for the array called TRIANG in Chapter
C        12.4 of the MATH77 documentation.  The next 4*NDP + 4 locations
C        are used for the array called B in Chapter 12.4 of the MATH77
C        documentation, and the rest of IWR (NDP locations) are use as
C        working space.
C MD     Type of computation desired, indicated as follows.
C    = 1   Specify new or original values for XYZDAT(*,1) and
C          XYZDAT(*,2).  Computes a new triangular grid.
C    = 2   As for 1, but also computes DZ, the partials.
C    = 3   The values of Z have changed, discard informations based on
C          the old values.
C    = 4   As for 3, but also computes DZ, the partials.
C    = 5   Compute an interpolated value from XI and YI, and output the
C          result in ZI(1).
C    = 6   As for MD = 6, but in addition store the partial derivative
C          of the surface at (XI, YI) with respect to X in ZI(2), and
C          that with respect to Y in ZI(3).
C MT     Dimension of TRIANG, must be > 12*NDP - 6*NB - 13.  IF NB is
C        not known, IF MT > 12*NDP - 13 you are safe.
C NDP    Number of DATA points.  Must have NDP > 3.
C TRIANG Defines the triangular grid.  Each triangle is defined by six
C        INTEGER indices identifying up to three neighboring triangles
C        and three vertex points.  The number of triangles will be
C        NT = 2 * (NDP-1) - NB, where NB is the number of points on the
C        boundary of the convex hull.
C WRK    Floating point work array with dimesion .GE. max(NDP, 30).  The
C        first 28 locations are used for SAVWRK as described in Chapter
C        12.4 of the MATH77 documentation, the next two locations are
C        used to remember NT and NDP respectively.  The first NDP
C        locations are used as work space when definging the grid.
C XI,YI  Defines point where interpolation is to be done.  Only used
C        when MD > 2.
C XYZDAT Array defining the surface.
C    XYZDAT(i, 1) The x coordinate for the i-th DATA point.
C    XYZDAT(i, 2) The y coordinate for the i-th DATA point.
C    XYZDAT(i, 3) The value of the surface at the i-th DATA point.
C ZI     When MD > 4, gives results as indicated above for MD = 5 & 6.
C
C
C *********************** Variable Definitions *************************
C
C DZ     Formal argument, see above.
C IERROR Formal argument, see above.
C INFO   Passed to DTGRD to get error and status info.  IERROR is set
C        to INFO(1). See Chapter 12.4 of the MATH77 documentation for
C        details on INFO(2) and INFO(3).
C IWK    Formal argument, see above.
C MB     Space required for the array called B in Chapter 12.4 of the
C        MATH77 documentation.
C MD     Formal argument, see above.
C MD0    Internal variable = MD.
C MT     Space required for the array called TRIANG in Chapter 12.4 of
C        the MATH77 documentation.
C NDP    Formal argument, see above.
C NDP0   Saved internal variable = NDP.
C NCONT  PARAMETER set to 1 to get a continuous interpolant.  Set to 0
C        IF you don't need or want a continuous interpolant.
C NB     Only referenced in comments.  The number of boundary points in
C        the convex hull of the triangularation.  This number is in
C        INFO(2) on a RETURN after a CALL with MD = 1.
C NT     Number of triangles in the grid.
C SAVWRK Array used to save current state of interpolations.
C WRK    Formal argument, see above.
C XI     Formal argument, see above.
C XYZDAT Formal argument, see above.
C YI     Formal argument, see above.
C ZI     Formal argument, see above.
C
C
C *********************** Variable Declarations ************************
C
C Formals
      INTEGER MD, NDP, NDPMAX
      INTEGER IERROR, IWK(*)
      DOUBLE PRECISION DZ(2,*), WRK(*), XYZDAT(NDPMAX,3), XI, YI, ZI(3)
C Others
C
C              Set to NCONT 1 for continuous interpolation, 0 otherwise.
      INTEGER NCONT
      PARAMETER (NCONT = 1)

      INTEGER INFO(3), MB, MD0, MT, NDP0, NT
      DOUBLE PRECISION Q(2)
C
C *********************** Start of Executable Code *********************
C
      MD0 = MD
      NDP0 = NDP
      MT = 12 * (NDP0 - 1)
      MB = NDP0 + 1
      IF (NDP .LT. 3) THEN
         IERROR = 6
         CALL ERROUT(IERROR)
         RETURN
      ELSE IF ((MD0 .LE. 0) .OR. (MD0 .GT. 6)) THEN
         IERROR = 8
         CALL ERROUT(IERROR)
      ELSE IF (MD0 .LE. 2) THEN
C            Initialize the triangular grid for new XYZDAT(*,1:2)
         CALL DTGGRD(XYZDAT(1,1), XYZDAT(1,2), NDP0, IWK(MT+4*MB+1),
     &      WRK, IWK, MT, IWK(MT), MB, NT, INFO)
         WRK(29) = NT
         WRK(30) = NDP0
         IERROR = INFO(1)
         IF (IERROR .NE. 0) THEN
            CALL ERROUT(IERROR)
            RETURN
         END IF
      ELSE IF (NDP0 .NE. int(WRK(30))) THEN
         IERROR = 7
         CALL ERROUT(IERROR)
         RETURN
      END IF
      IF (MD0 .LE. 4) THEN
         WRK(1) = 0.D0
         NT = WRK(29)
         IF ((MD0 .EQ. 2) .OR. (MD0 .EQ. 4)) THEN
C                                            Get new partials.
            CALL DTGPD(XYZDAT(1,1),XYZDAT(1,2),XYZDAT(1,3),
     &         	       DZ,NDP0,IWK,NT,IWK(MT+4*MB+1))
         END IF
      ELSE

C	 ====================
C        DO the interpolation
C        ====================
         Q(1) = XI
         Q(2) = YI
         NT = WRK(29)
         CALL DTGFI(XYZDAT(1,1),XYZDAT(1,2),XYZDAT(1,3),DZ,NDP0,
     &          IWK,NT,NCONT,Q,ZI(1),MD0.EQ.6,ZI(2:3),IERROR,WRK)
         IF (IERROR .NE. 0) CALL ERROUT(IERROR-2)
       END IF
       RETURN
       END SUBROUTINE DIDES

C***********************************************************************

      SUBROUTINE ERROUT (INDEX)
C        Called by DIDES to PRINT error message text selected by INDEX.
      INTEGER INDEX
      CHARACTER ERRTXT(-1:8)*60
      DATA ERRTXT /
     1 'The point (XI, YI) is outside the convex hull.',
     2 'Some kind of error in the triangular grid.',
     3 'Couldn''t complete optimization of triangular grid.',
     4 'Duplicate points or some points colinear.',
     5 'Some pair of the points are duplicates.',
     6 'The dimension MB is too small, MB = NDP+1 will suffice.',
     7 'The dimension MT is too small, MT = 12*(NDP-1) will suffice.',
     8 ' NDP < 3.',
     9 ' NDP changed from last CALL.',
     A ' MD < 1 or MD > 6.' /
C
C D. Redding -- modified to prevent massive error outputs
      IF (INDEX.NE.-1) WRITE (*, '(1X,A60)') ERRTXT(INDEX)
      RETURN
      END SUBROUTINE ERROUT

C***********************************************************************
C***********************************************************************

      SUBROUTINE DTGC0(S, ZOUT, WANTDZ, DZOUT)
C>> 1996-02-02 CLL
C>> 1996-01-11 CLL
C>> 1995-11-01 CLL
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C     C.L.LAWSON, JPL, 1976 DEC 7
C  This subr interpolates over a triangle using linear interpolation.
C  This method gives C0 continuity with neighboring triangles, i.e.,
C  continuity of the value but generally not continuity of first and
C  higher order partial derivatives.
C  Optionally the subr also computes first partial derivatives.
C
C   (S(1:3, 1:4) [inout]  Columns 1-4 contain DATA set by the user to
C          specify the interpolation problem.  Col 1 contains
C          unnormalized barycentric coordinates of the interpolation
C          point.  The other cols contain DATA depending only on the
C          triangle and its vertex DATA and not on the interpolation
C          point.
C
C         S( ,1) = UNNORMALIZED BARYCENTRIC COORDS OF INTERP POINT.
C         S( ,2) = U = X COORD OF EDGE VECTOR
C         S( ,3) = V = Y COORD OF EDGE VECTOR
C         S( ,4) = Z = FCN VALUE AT VERTEX
C
C   ZOUT [out]  INTERPOLATED VALUE COMPUTED BY SUBR.
C
C   WANTDZ [in]  =.TRUE. MEANS COMPUTE DZOUT() AS WELL AS ZOUT.
C              =.FALSE. MEANS COMPUTE ONLY ZOUT AND NOT DZOUT().
C
C   DZOUT(1:2) [out]  First partial derivs w.r.t. X and Y of the
C      interpolated surface at the interpolation point.
C      Note that since this SUBROUTINE computes only a linear
C      interpolant over a triangle, the partial derivs w.r.t. X and Y
C      will be the same at every point of the triangle.  In general
C      these derivatives will jump to different values when moving to
C      an adjacent triangle.
C     ------------------------------------------------------------------
C         Details of the contents of TRI(1:7) and S(1:3, 1:4).
C
C  We assume the point Q is in the triangle indexed by INDTRI.
C  The indices of the vertices of this triangle, in counter-
C  clockwise, order are P(1), P(2), and P(3) which may be obtained as
C  P(i) = TRI(3+i) for i = 1, 2, and 3.  TRI(7) contains the same
C  value as TRI(4).  The triangle adjacent to this triangle across
C  the edge from P(i) to P(i+1) is indexed by TRI(i) for i = 1, 2,
C  and 3.  IF there is no adjacent triangle across this edge THEN
C  TRI(i) = 0.
C
C  For descriptive convenience we regard the subscript of P() and the
C  1st subsubscript of S(,) as always being reduced modulo 3 to 1, 2,
C  or 3.  Also for convenience we shall WRITE P(i) to mean the vertex
C  indexed by P(i).
C
C  The unnormalized barycentric coordinate that is zero along the edge
C  from P(i) to P(i+1) and has a positive value at P(i+2) is stored in
C  S(i,1).  The (x,y) coordinates of the vector from P(i) to P(i+1) are
C  stored in (S(i,2), S(i,3)).
C  The FUNCTION value Z at P(i+2) is stored in S(i,4) for i = 1, 2, & 3.
C     ------------------------------------------------------------------
C--D replaces "?": ?TGC0
C     ------------------------------------------------------------------
      DOUBLE PRECISION delz(2), denom, dzout(2), S(3,4), zout
      LOGICAL  WANTDZ
C     ------------------------------------------------------------------
C        To compute zout:
C        Combine the vertex FUNCTION values using the unnormalized
C        barycentric coordinates as weights, and normalize by the
C        sum of the unnormalized barycentric coordinates.
C
      zout = (s(1,1)*s(1,4) +s(2,1)*s(2,4) +s(3,1)*s(3,4)) /
     &            (S(1,1)+S(2,1)+S(3,1))
      IF(wantdz) THEN
         delz(1) = s(3,4) - s(2,4)
         delz(2) = s(1,4) - s(3,4)
         denom = s(1,2) * s(2,3) - s(2,2) * s(1,3)
         IF(denom .NE. 0.0d0) THEN
            dzout(1) = (delz(1) * s(2,3) - delz(2) * s(1,3)) / denom
            dzout(2) = (s(1,2) * delz(2) - s(2,2) * delz(1)) / denom
         ELSE
            dzout(1) = 0.0d0
            dzout(2) = 0.0d0
         END IF
      END IF
      RETURN
      END SUBROUTINE DTGC0

C***********************************************************************
C***********************************************************************

      SUBROUTINE DTGC1(NEWTRI, S, ZOUT, WANTDZ, DZOUT)
C>> 1996-02-02 CLL
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C     C.L.Lawson, JPL, 1976 Dec 7
C     This subr interpolates over a triangle using the piecewise
C     bicubic formulas of  Clough and Tocher (1965) as formulated for
C     computation by C.L.Lawson in JPL Tech. Memorandum 33-770, May,1976
C     This method gives C1 continuity with neighboring triangles, i.e.,
C     continuity of value and first partial derivatives.
C     Generally there will be jumps in second and higher order
C     derivatives across edges between triangles, and along implicit
C     edges that partition each triangle into three smaller triangles.
C     Optionally the subr also computes first partial derivatives.
C
C   NEWTRI [in]  This argument allows certain computations to be
C          skipped, saving 15 multiplies, 3 divides, and 9 adds,
C          when the user has made no changes to the contents of
C          columns 2, 3, 5, 6, 7, 8, & 9 of the S() array since a
C          previous CALL.  (Thus it is allowable to have changed the
C          contents of columns 1 and 4.)
C          Setting NEWTRI = .TRUE. causes these computations to be done,
C          while .FALSE. cause them to be skipped.
C
C   (S(1:3, 1:9) [inout]  Columns 1-6 contain DATA set by the user to
C          specify the interpolation problem.  Col 1 contains
C          unnormalized barycentric coordinates of the interpolation
C          point.  The other cols contain DATA depending only on the
C          triangle and its vertex DATA and not on the interpolation
C          point.  Cols 7-9 are computed by this SUBROUTINE, and in
C          appropriate cases can be saved from one CALL to the next.
C          See NEWTRI above.
C
C         S( ,1) = UNNORMALIZED BARYCENTRIC COORDS OF INTERP POINT.
C         S( ,2) = U = X COORD OF EDGE VECTOR
C         S( ,3) = V = Y COORD OF EDGE VECTOR
C         S( ,4) = Z = FCN VALUE AT VERTEX
C         S( ,5) = ZX = PARTIAL DERIV OF Z W.R.T. X
C         S( ,6) = ZY = PARTIAL DERIV OF Z W.R.T. Y
C         S( ,7) = HTILDA = SCALED P.D. TANGENTIAL TO EDGE AT LEFT END.
C         S( ,8) = KTILDA = SCALED P.D. TANGENTIAL TO EDGE AT RIGHT END.
C         S( ,9) = LFAC   = FACTOR INVOLVING LENGTHS OF EDGES.
C
C   ZOUT [out]  INTERPOLATED VALUE COMPUTED BY SUBR.
C
C   WANTDZ [in]  =.TRUE. MEANS COMPUTE DZOUT() AS WELL AS ZOUT.
C              =.FALSE. MEANS COMPUTE ONLY ZOUT AND NOT DZOUT().
C
C   DZOUT(1:2) [out]  FIRST PARTIAL DERIVS W.R.T. X AND Y OF THE
C              INTERPOLATED SURFACE AT THE INTERPOLATION POINT.
C     ------------------------------------------------------------------
C         Details of the contents of TRI(1:7) and S(1:3, 1:6).
C
C  We assume the point Q is in the triangle indexed by INDTRI.
C  The indices of the vertices of this triangle, in counter-
C  clockwise, order are P(1), P(2), and P(3) which may be obtained as
C  P(i) = TRI(3+i) for i = 1, 2, and 3.  TRI(7) contains the same
C  value as TRI(4).  The triangle adjacent to this triangle across
C  the edge from P(i) to P(i+1) is indexed by TRI(i) for i = 1, 2,
C  and 3.  IF there is no adjacent triangle across this edge THEN
C  TRI(i) = 0.
C
C  For descriptive convenience we regard the subscript of P() and the
C  1st subsubscript of S(,) as always being reduced modulo 3 to 1, 2,
C  or 3.  Also for convenience we shall WRITE P(i) to mean the vertex
C  indexed by P(i).
C
C  The unnormalized barycentric coordinate that is zero along the edge
C  from P(i) to P(i+1) and has a positive value at P(i+2) is stored in
C  S(i,1).  The (x,y) coordinates of the vector from P(i) to P(i+1) are
C  stored in (S(i,2), S(i,3)).
C  The FUNCTION value Z at P(i+2) is stored in S(i,4) for i = 1, 2, & 3.
C  ZX, the partial deriv of Z w.r.t. X at P(i+2) is stored in S(i,5)
C     for i = 1, 2, & 3.
C  ZY, the partial deriv of Z w.r.t. Y at P(i+2) is stored in S(i,6)
C     for i = 1, 2, & 3.
C     ------------------------------------------------------------------
C--D replaces "?": ?TGC1
C     ------------------------------------------------------------------
      INTEGER          ADD1(3), i, im1, ip1, j, m, mm1, mp1, SUB1(3)
      DOUBLE PRECISION A, B, C
      DOUBLE PRECISION DA, DB, DC, DGTILD(3), DPHI(3)
      DOUBLE PRECISION DR(3), DRHO(3), DZOUT(2)
      DOUBLE PRECISION fac, GTILDA(3), LSQ(3)
      DOUBLE PRECISION r(3),  RHO(3), PHI(3), S(3,9), sum
      DOUBLE PRECISION zout
      LOGICAL  WANTDZ,NEWTRI
      DATA ADD1(1),ADD1(2),ADD1(3)/ 2, 3, 1 /
      DATA SUB1(1),SUB1(2),SUB1(3)/ 3, 1, 2 /
C     ------------------------------------------------------------------
C              COMPUTE R() = NORMALIZED BARYCENTRIC COORDINATES.
C              AND LSQ()   = SQUARED EDGE LENGTHS
C
      FAC = 1.D00/(S(1,1)+S(2,1)+S(3,1))
      DO 10 I=1,3
         R(I) = FAC*S(I,1)
         LSQ(I)= S(I,2)**2 + S(I,3)**2
   10 CONTINUE
C                             COMPUTE PHI(I) = R(I+1) * R(I-1)
      PHI(1) = R(2)*R(3)
      PHI(2) = R(3)*R(1)
      PHI(3) = R(1)*R(2)
C
C              COMPUTE RHO() = CLOUGH-TOCHER PIECEWISE BICUBIC
C                              CORRECTION FUNCTION
C
      IF ( R(1) .LE. R(2)) THEN
      IF ( R(1) .LE. R(3)) THEN
            M = 1
      ELSE
            M = 3
      END IF
      ELSE
      IF ( R(2) .LE. R(3)) THEN
            M = 2
      ELSE
            M = 3
      END IF
      END IF
C
      A = 0.5D0 * R(M)**2
      B = (1.0D0/3.0D0)*R(M)
      C = PHI(M) + (5.0D0/3.0D0)*A
      RHO(M) = R(M) * C - A
      MP1 = ADD1(M)
      MM1 = SUB1(M)
      RHO(MP1) = A*(R(MM1) - B)
      RHO(MM1) = A*(R(MP1) - B)
C
      SUM = 0.0d0
        DO 20 I=1,3
         IP1 = ADD1(I)
         IM1 = SUB1(I)
C
         IF (NEWTRI ) THEN
C                                                               H TILDA
            S(I,7) = S(I,2)*S(IP1,5) + S(I,3)*S(IP1,6)
C                                                               K TILDA
            S(I,8) = S(I,2)*S(IM1,5) + S(I,3)*S(IM1,6)
C                                                               LFAC
            S(I,9) = 3.0D0*(LSQ(IP1)-LSQ(IM1)) / LSQ(I)
         END IF
C                                                               G TILDA
         GTILDA(I)= (R(IP1)-R(IM1))*PHI(I) + S(I,9)*RHO(I)
     &              -RHO(IP1)+RHO(IM1)
         SUM = SUM + S(I,7)*(GTILDA(I)+PHI(I))
     &             + S(I,8)*(GTILDA(I)-PHI(I))
   20 CONTINUE
      ZOUT = 0.5D0 * SUM
C
      DO 30 I=1,3
         IP1 = ADD1(I)
         IM1 = SUB1(I)
         ZOUT=ZOUT + S(I,4)*(R(I)+GTILDA(IM1)-GTILDA(IP1))
   30 CONTINUE
C                                  FINISHED COMPUTING ZOUT.
C
C                                  NOW COMPUTE DZOUT() IF REQUESTED.
      IF ( WANTDZ) THEN
         DO 100 J=1,2
            IF(j .EQ. 1) THEN
               DO 40 I=1,3
                  DR(I) = -FAC*S(I,3)
   40          CONTINUE
            ELSE
               DO 50 I=1,3
                  DR(I) = +FAC*S(I,2)
   50          CONTINUE
            END IF
C
            DPHI(1) = R(2)*DR(3) + DR(2)*R(3)
            DPHI(2) = R(3)*DR(1) + DR(3)*R(1)
            DPHI(3) = R(1)*DR(2) + DR(1)*R(2)
C
            DA = R(M)*DR(M)
            DB = (1.D0/3.D0)*DR(M)
            DC = DPHI(M) + (5.D0/3.D0)*DA
            DRHO(M) = R(M)*DC + DR(M)*C - DA
            DRHO(MP1) = A*(DR(MM1) - DB) + DA*(R(MM1) - B)
            DRHO(MM1) = A*(DR(MP1) - DB) + DA*(R(MP1) - B)
C
            SUM = 0.0d0
            DO 60 I=1,3
               IP1=ADD1(I)
               IM1=SUB1(I)
               DGTILD(I)  =(R(IP1)-R(IM1))*DPHI(I)
     &                  +(DR(IP1) - DR(IM1))*PHI(I)
     &                  +S(I,9)*DRHO(I) - DRHO(IP1) + DRHO(IM1)
               SUM = SUM +S(I,7)*(DGTILD(I)+DPHI(I))
     &                +S(I,8)*(DGTILD(I)-DPHI(I))
   60       CONTINUE
            SUM = 0.5D0 * SUM
            DO 70 I=1,3
               IP1=ADD1(I)
               IM1=SUB1(I)
               SUM = SUM + S(I,4)*(DR(I) + DGTILD(IM1) - DGTILD(IP1))
   70       CONTINUE
            DZOUT(J) = SUM
  100    CONTINUE
      END IF
      RETURN
      END SUBROUTINE DTGC1

C***********************************************************************
C***********************************************************************

      SUBROUTINE DTGFI(X,Y,Z,DZ,NP,TRIANG,NT,NCONT,Q,ZOUT,
     &                 WANTDZ,DZOUT,MODE,SAVWRK)
C>> 1996-02-02 CLL
C>> 1996-01-11 CLL
C>> 1995-10-30 CLL
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C  This routine uses _TRFND to DO lookup, and either _TRC1 or _TRC2
C  to DO interpolation in a triangular grid.
C     ------------------------------------------------------------------
C                         SUBROUTINE Arguments
C
C   X(1:NP), Y(1:NP) [in]  (x,y) coordinates of vertices of the
C                          triangular grid.
C
C   Z(1:NP) [in]  Z(i) is the value at (X(i),Y(i)) of the DATA to
C                 be interpolated.
C
C   DZ(1:2, 1:NP) [in]  DZ(1:2, i) are the values of the partial
C         derivatives the interpolation FUNCTION at (X(i),Y(i)) with
C         respect to x and y, respectively.
C
C   NP [in]  No. of vertices in the triangular grid.
C
C   TRIANG(1:6*NT) [in]  Array of INTEGER pointers defining the
C            connectivity of the triangular grid.
C
C   NT [in]  No. of triangles in the triangular grid.
C
C   NCONT [in]  = 0 or 1 to request either C0 or C1 continuity.
C
C   Q(1:2) [in]  The (x,y) coordinates of the point for which this
C             subr will attempt to find an enclosing triangle and THEN
C             DO interpolation.
C
C   ZOUT [out]  Interpolated value computed by this SUBROUTINE.
C
C   WANTDZ [in]  =.TRUE. means compute DZOUT() as well as ZOUT.
C              =.FALSE. means compute only ZOUT and not DZOUT().
C
C   DZOUT(1:2) [out]  First partial derivs w.r.t. x and y of the
C              interpolated surface at the interpolation point.
C
C   MODE [out]
C          =0  Means OK.  Q is interior or almost so.
C          =1  Q is exterior by more than the built-in tolerance.
C          =2  Bad.  SUBROUTINE is cycling in the search.  This should
C              never happen.
C
C   SAVWRK(1:28) [inout]  This array is used as work space and to save
C      quantities that may be reusable on a subsequent CALL.
C      Let us say the settings of X(), Y(), Z(), DZ(), NP, TRIANG(),
C      and NT define a "grid and DATA".  On the first CALL to this
C      SUBROUTINE for use of a particular "grid and DATA" the user must
C      set SAVWRK(1) = 0.  This SUBROUTINE will generally reset
C      SAVWRK(1) to a nonzero value on RETURN, indicating there are
C      saved quantities in SAVWRK() relating to this "grid and DATA".
C      IF the user will be switching back and forth between different
C      "grid and DATA" specifications, some saving of computation time
C      can be achieved by using distinct SAVWRK() arrays associated with
C      distinct "grid and DATA" specifications.  It is the user's
C      responsibility to set SAVWRK(1) = 0. on any CALL for which the
C      current SAVWRK() array is not the one that was used on the most
C      recent previous CALL involving the current "grid and DATA".
C
C     ------------------------------------------------------------------
C              How the array SAVWRK() is used internally.
C
C     The 1st location of SAVWRK() is regarded as containing a value
C     of INDSGN, and the remaining 27 locations contain S(1:3, 1:9).
C
C     We start by setting       INDSGN = SAVWRK(1)
C     and                       INDTRI = abs(INDSGN)
C     We start the current lookup at INDTRI.  IF INDTRI is positive this
C     causes the current lookup to start in the triangle at which the
C     previous lookup for this "grid and DATA" ended, which in some
C     cases significantly reduces the search time.
C
C     IF INDSGN > 0, and NCONT = 2, and the value of INDTRI does not
C     change during the lookup, we assume there is useful saved info in
C     S(1:3,7:9).  Otherwise we compute this info as needed.
C
C                    Contents of S(,)
C
C   (s(1:3, 1:9) [inout]  Columns 1-6 contain DATA specifying the
C              interpolation problem.  Col 1 contains unnormalized
C              barycentric coordinates of the interpolation point.
C              The other cols contain DATA depending only on the
C              triangle and its vertex DATA and not on the
C              interpolation point.  Cols 7-9 are used to save computed
C              quantities depending only on the triangle and the
C              partial deriv values from one CALL to the next.
C         s( ,1) = Unnormalized barycentric coords of interp point.
C         s( ,2) = u = x coord of edge vector
C         s( ,3) = v = y coord of edge vector
C         s( ,4) = z = Fcn value at vertex
C         s( ,5) = zx = Partial deriv of z w.r.t. x
C         s( ,6) = zy = Partial deriv of z w.r.t. y
C         s( ,7) = htilda = Scaled p.d. tangential to edge at left END.
C         s( ,8) = ktilda = Scaled p.d. tangential to edge at right END.
C         s( ,9) = lfac   = Factor involving lengths of edges.
C
C  Let the indices of the vertices of this triangle, in counter-
C  clockwise order be denoted by P(1), P(2), and P(3).
C  For descriptive convenience we regard the subscript of P() and the
C  1st subsubscript of S(,) as always being reduced modulo 3 to 1, 2,
C  or 3.  Also for convenience we shall WRITE P(i) to mean the vertex
C  indexed by P(i).
C
C  The unnormalized barycentric coordinate that is zero along the edge
C  from P(i) to P(i+1) and has a positive value at P(i+2) is stored in
C  S(i,1).  The (x,y) coordinates of the vector from P(i) to P(i+1) are
C  stored in (S(i,2), S(i,3)).
C
C  The FUNCTION value Z at P(i+2) is stored in S(i,4) for i = 1, 2, & 3.
C  ZX, the partial deriv of Z w.r.t. X at P(i+2) is stored in S(i,5)
C     for i = 1, 2, & 3.
C  ZY, the partial deriv of Z w.r.t. Y at P(i+2) is stored in S(i,6)
C     for i = 1, 2, & 3.
C     ------------------------------------------------------------------
C                      Internal Variables
C
C  NEWTRI [in]  = .FALSE. means that saved quantities are
C              present in  s(1:3, 7:9).9).
C              IF newtri=.TRUE. the subr will compute s(1:3, 7:9).
C
C  INDTRI  On entry to _TRFND, indtri designates the triangle at which
C      the search will begin.  On RETURN from _TRFND, indtri is the
C      index of the last triangle tested.
C      This is the triangle containing q IF mode=0.
C      IF mode = 1 this is a boundary triangle and q is outside
C      a boundary edge of this triangle.
C
C  TRI(1:8)  INTEGER array of pointers defining one triangle.  It is
C     never assumed to contain any saved values on entry to
C     _TRFND.  On RETURN from _TRFND, TRI(1:7) defines the triangle
C     whose index is the returned value of INDTRI.
C
C     The indices of the vertices of this triangle, in counter-
C     clockwise order, are P(1), P(2), and P(3), which may be obtained
C     as P(i) = TRI(3+i) for i = 1, 2, and 3.  TRI(7) contains the same
C     value as TRI(4).  TRI(8) will be set in this subr to equal TRI(5).
C     The triangle adjacent to this triangle across the edge from P(i)
C     to P(i+1) is indexed by TRI(i) for i = 1, 2, and 3.  IF there is
C     no adjacent triangle across this edge THEN TRI(i) = 0.
C
C
C     ------------------------------------------------------------------
C--D replaces "?": ?TGFI, ?TGFND, ?TGC0, ?TGC1
C     ------------------------------------------------------------------
      INTEGER indsgn, indtri, mode, ncont, np, nt
      INTEGER tri(8), triang(1:*)
      DOUBLE PRECISION dz(2,np), dzout(2), q(2), savwrk(28)
      DOUBLE PRECISION x(np), y(np), z(np), zero, zout
      LOGICAL newtri, wantdz
      PARAMETER(zero = 0.0d0)
C     ------------------------------------------------------------------
      indsgn = savwrk(1)
      indtri = abs(indsgn)
      CALL DTGFND(x,y,np, triang,nt,q,indtri,tri,savwrk(2),mode)
      IF (mode .EQ. 0) THEN
C
C           The point q(1:2) was found ok.
C           tri(1:7) and s(1:3, 1:3) contain info set by DTRFND.
C           Next move info into s(1:3, 4:6) and CALL DTRC0 or DTRC1 for
C           interpolation.
C           Note that savwrk(-2+i + 3*j) contains s(i,j).
C
C           S(1,4) = Z(TRI(6))
C           S(2,4) = Z(TRI(4))
C           S(3,4) = Z(TRI(5))
C
         savwrk(1+10) = Z(TRI(6))
         savwrk(1+11) = Z(TRI(4))
         savwrk(1+12) = Z(TRI(5))
         IF(ncont .EQ. 0) THEN
            CALL DTGC0(savwrk(2),zout,wantdz,dzout)
            indsgn = -indtri
         ELSE
C
C           S(1,5) = DZ(1,TRI(6))
C           S(2,5) = DZ(1,TRI(4))
C           S(3,5) = DZ(1,TRI(5))
C
C           S(1,6) = DZ(2,TRI(6))
C           S(2,6) = DZ(2,TRI(4))
C           S(3,6) = DZ(2,TRI(5))

            savwrk(1+13) = DZ(1,TRI(6))
            savwrk(1+14) = DZ(1,TRI(4))
            savwrk(1+15) = DZ(1,TRI(5))

            savwrk(1+16) = DZ(2,TRI(6))
            savwrk(1+17) = DZ(2,TRI(4))
            savwrk(1+18) = DZ(2,TRI(5))

            newtri = indtri .NE. indsgn
            CALL DTGC1(newtri, savwrk(2), zout, wantdz ,dzout)
            indsgn = indtri
         END IF
      ELSE
C           Here mode = 1 or 2.
C           1 means xy is outside the triangular grid.
C           2 means cycling has happened in DTGFND. this
C           should not happen.  Error message will have been printed
C           from DTGFND.
C
         zout = zero
         IF (wantdz) THEN
            dzout(1) = zero
            dzout(2) = zero
         END IF
         indsgn = -indtri
      END IF
      savwrk(1) = indsgn
      RETURN
      END SUBROUTINE DTGFI

C***********************************************************************
C***********************************************************************

      SUBROUTINE DTGFND(X,Y,NP, TRIANG,NT,Q,INDTRI,TRI,S,MODE)
C>> 1996-03-30 DTGFND  Krogh  Removed Fortran 90 comments.
C>> 1996-02-02 CLL
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C
C     C.L.LAWSON, JPL, 1976 DEC 3
C         LOOK UP POINT Q IN TRIANGULAR GRID.
C
C   X(1:NP), Y(1:NP) [in]  (x,y) coordinates of vertices of the
C                          triangular grid.
C
C   NP [in]  No. of vertices.
C
C   TRIANG(1:6*NT) [in]  Array of INTEGER pointers defining the
C            connectivity of the triangular grid.
C
C   NT [in]  No. of triangles in the triangular grid.
C
C   Q(1:2) [in]  The (x,y) coordinates of the point for which this
C             subr will attempt to find an enclosing triangle.
C
C   INDTRI [in]  DESIGNATES TRIANGLE AT WHICH SEARCH WILL BEGIN.
C              ON RETURN INDTRI IS THE INDEX OF THE LAST TRIANGLE
C              TESTED.  THIS IS THE TRIANGLE CONTAINING Q IF MODE=0.
C              IF MODE = 1 THIS IS A BOUNDARY TRIANGLE AND Q IS OUTSIDE
C              A BOUNDARY EDGE OF THIS TRIANGLE.
C
C   TRI(1:7) [out]  INTEGER ARRAY OF POINTERS ASSOCIATED
C              WITH TRIANGLE INDTRI.
C
C   S(1:3, 1:3) [out]  IF MODE=0 on RETURN, S(,) CONTAINS
C              S(:,1) = UNNORMALIZED BARYCENTRIC COORDINATES OF Q.
C              S(:,2) = X COORDINATES OF EDGE VECTORS.
C              S(:,3) = Y COORDINATES OF EDGE VECTORS.
C
C   MODE [out]
C          =0  MEANS OK.  Q IS INTERIOR OR ALMOST SO.
C          =1  Q IS EXTERIOR BY MORE THAN THE BUILT-IN TOLERANCE.
C          =2  BAD.  SUBR IS CYCLING IN THE SEARCH.  THIS SHOULD NEVER
C              HAPPEN.
C
C     ------------------------------------------------------------------
C         Details of the contents of TRI(1:7) and S(1:3, 1:3).
C
C  TRI() contains INTEGER pointers defining the triangle indexed by
C  INDTRI.  The indices of the vertices of this triangle, in counter-
C  clockwise order are P(1), P(2), and P(3) which may be obtained as
C  P(i) = TRI(3+i) for i = 1, 2, and 3.  TRI(7) contains the same
C  value as TRI(4).  The triangle adjacent to this triangle across
C  the edge from P(i) to P(i+1) is indexed by TRI(i) for i = 1, 2,
C  and 3.  IF there is no adjacent triangle across this edge THEN
C  TRI(i) = 0.
C
C  For descriptive convenience we regard the subscript of P() and the
C  1st subsubscript of S(,) as always being reduced modulo 3 to 1, 2,
C  or 3.  Also for convenience we shall WRITE P(i) to mean the vertex
C  indexed by P(i).
C
C  The unnormalized barycentric coordinate that is zero along the edge
C  from P(i) to P(i+1) and has a positive value at P(i+2) is stored in
C  S(i,1).  The (x,y) coordinates of the vector from P(i) to P(i+1) are
C  stored in (S(i,2), S(i,3)).
C     ------------------------------------------------------------------
C--D replaces "?": ?TGFND, ?TGGET
C     ------------------------------------------------------------------
      INTEGER icount, indtri, iskip, j, j1, j2, j2save, MODE, NP, NT
      INTEGER TRI(7), TRIANG(1:*)
      DOUBLE PRECISION p1, p2, Q(2), S(3,3), X(NP), Y(NP)
C     ------------------------------------------------------------------
      IF( INDTRI .LE. 0 .OR. INDTRI .GT. NT) INDTRI=1
      MODE = 2
      ISKIP  = 0
      J2SAVE = 0
      DO 50 ICOUNT = 1, NT
         CALL DTGGET(INDTRI,TRI,TRIANG)
         TRI(7)=TRI(4)

            DO 10 J=1,3
               J1= TRI(J+3)
               IF ( J1 .EQ. J2SAVE) THEN
                  ISKIP=J
               ELSE
                  J2 = TRI(J+4)
                  S(J,2) = X(J2)-X(J1)
                  P2 = S(J,2) * (Q(2)-Y(J1))
                  S(J,3) = Y(J2)-Y(J1)
                  P1 = S(J,3) * (Q(1)-X(J1))
                  S(J,1) = -P1+P2
C
                  IF (S(J,1) .LT. -1.0D-4*(ABS(P1)+ABS(P2))) THEN
C
C                        Q IS OUTSIDE TRIANGLE INDTRI BY MORE THAN
C                        TOLERANCE.  MOVE TO NEIGHBORING TRIANGLE IF
C                        THERE IS ONE.
C
                     IF (TRI(J) .NE. 0) THEN
                        J2SAVE = J2
                        INDTRI = TRI(J)
                        GO TO 40
                     ELSE
C D. Redding -- modified to allow extrapolation
C                        MODE = 0
                        MODE = 1
                        GO TO 80
                     END IF
                  END IF
               END IF
   10       CONTINUE
C                                  Q IS INSIDE OR ALMOST INSIDE
C                                  THIS TRIANGLE
            MODE = 0
            GO TO 90
   40    CONTINUE
   50 CONTINUE
   80 CONTINUE
   90 CONTINUE
C                             COMPUTE SKIPPED INFO FOR USE IN
C                             INTERPOLATION SUBR.
      IF (ISKIP .NE. 0) THEN
         J1=TRI(ISKIP+3)
         J2=TRI(ISKIP+4)
         S(ISKIP,2)= X(J2)-X(J1)
         S(ISKIP,3)= Y(J2)-Y(J1)
         S(ISKIP,1)= -S(ISKIP,3)*(Q(1)-X(J1))
     &               +S(ISKIP,2)*(Q(2)-Y(J1))
      END IF
      RETURN
      END SUBROUTINE DTGFND

C***********************************************************************

C  File: dtggrd.[f|for] contains DTGGRD, DTGANG, DTGADJ.
C
      SUBROUTINE DTGGRD(X,Y,NP,IP,W,TRIANG,MT,B,MB,NT, INFO)
C>> 1996-03-30 DTGGRD Krogh Added EXTERNAL statement, => generics.
C>> 1996-03-04 CLL Changed specifications of INFO().
C>> 1996-02-02 CLL
C>> 1996-01-12 CLL
C>> 1996-01-11 CLL
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C>> 1991-11-20 CLL Using Math77 SUBROUTINE DSORTP for sorting.
C
C  THIS SUBR CONSTRUCTS A LIST OF POINTERS IN THE ARRAY
C  TRIANG() DEFINING A TRIANGULAR GRID HAVING THE GIVEN (X,Y)
C  DATA AS VERTICES.
C  This will be a Delaunay triangulation.  In such a triangulation
C  the smallest interior angle in any triangle of the grid will be
C  maximal, i.e., the max-min angle property.  Also the circumcircle of
C  each triangle will contain no nodes of the grid in its interior.
C  A Delaunay triangulation is dual to a tessellation made up of
C  proximity regions associated with the names Dirichlet, Voronoi, and
C  Thiessen.
C
C     C. L. LAWSON, JPL, 1976 NOV 28
C     C.L.L., 1979 MAR 3. CONVERTED TO SFTRAN3.
C     C.L.L., 1979 JUL 22. ADDED CALL TO DTGSIZ.
C
C  X(1:NP),Y(1:NP) [in]    (X,Y) DATA POINTS.
C  NP [in]                    NO. OF DATA POINTS.  REQUIRE NP .GE. 3.
C
C  IP() [scratch]   WORK SPACE OF LENGTH at least NP.
C  W()  [scratch]   WORK SPACE OF LENGTH at least NP.
C  TRIANG() [out]  ARRAY IN WHICH THIS SUBR WILL BUILD A LIST DEFINING
C             A triangular grid.  access to this list is via the
C             four subrs DTGGET, DTGPUT, DTGSIZ, and DTGSET.
C  MT [in]   DIMENSION OF ARRAY TRIANG().  NO. OF TRIANGLES IT
C            CAN HOLD, NTLIM, IS DETERMINED BY USE OF SUBROUTINE
C            DTGSIZ.    NTLIM WILL BE MT/6 IF POINTERS ARE NOT
C            PACKED.  NTLIM = (3*MT)/6 IF POINTERS ARE PACKED
C            THREE PER WORD.
C  B(,) [out]  WORK SPACE OF DIMENSION (MB,4).
C  MB [in]   FIRST DIMENSION OF B(,).  SETTING MB = NP+1 WILL
C            ALWAYS BE ADEQUATE, HOWEVER A MUCH SMALLER VALUE
C            WILL OFTEN SUFFICE.  SUGGEST TRYING MB = 6 * (cube root
C            of NP) WHEN NP .GE. 64.
C  NT [out]  OUTPUT..  NO. OF TRIANGLES
C
C  INFO(1:3) [out, INTEGER]  Termination status information.
C
C     INFO(1) INDICATES STATUS ON TERMINATION.
C              0  NORMAL TERMINATION.  The triangular grid is complete
C                 and successfully optimized.
C              1  ERROR. The triangular grid is complete but not
C                 successfully optimized.  Apparent looping in DTGADJ.
C              2  ERROR.  ALL GIVEN POINTS ARE COLINEAR.  NO TRIANGLES
C                 CONSTRUCTED.  (OR ELSE THERE ARE SOME DUPLICATE PTS)
C              2  ERROR.  DUPLICATE POINTS.  TRIANGULAR GRID NOT
C                 COMPLETED.
C              4  ERROR.  NOT ENOUGH SPACE IN B(,).  USER MUST
C                 INCREASE MB.
C              5  NOT ENOUGH SPACE IN TRIANG().  USER MUST MAKE
C                 MT BIGGER.
C     INFO(2)  NO. OF BOUNDARY POINTS.
C     INFO(3)  SMALLEST SETTING OF MB THAT WOULD SUFFICE FOR THIS DATA.
C
C        [Relations to internal variables or previously used variable
C         names:  INFO(1) = MODE, INFO(2) = NB, INFO(3) = MBUSED,
C         The value of INFO(1) = MODE will be  set to 0 when
C         OK1 = .TRUE. and 1 when OK1 = .FALSE.]
C     ------------------------------------------------------------------
C--D replaces "?": ?TGGRD, ?TGANG, ?TGSET, ?TGGET, ?TGPUT, ?TGSIZ,
C--&               ?TGADJ, ?SORTP
C     ------------------------------------------------------------------
      INTEGER NP, mb, mt
      INTEGER AVAIL, B(MB,4), ext, TRI(6)
      INTEGER i, i1, i2, i3, i3prev, index, INFO(3), ip(np), isave
      INTEGER jfix, jhi, jlo, jnew
      INTEGER k, kb1, kbfix, kbhi, kblo, kbnew
      INTEGER loc004, loc005
      INTEGER mode
      INTEGER n20025
      INTEGER nb, nt, ntlim, mbused,  TRIANG(1:MT)
      DOUBLE PRECISION aref, cross, dx12, dy12, pangle
      DOUBLE PRECISION W(np), X(NP), xmean, Y(NP), ymean
      LOGICAL  OK, OK1
C     ------------------------------------------------------------------
      CALL DTGSIZ( MT, NTLIM )
      MODE = 0
      OK1 = .TRUE.
C
C         FIND AN EXTREME POINT (X(EXT),Y(EXT))
C
      EXT=1
      DO 20002 I=1,NP
            IF(     X(I) .LT. X(EXT)
     &        .OR. (X(I) .EQ. X(EXT)  .AND.  Y(I) .LT. Y(EXT) )) EXT=I
20002 CONTINUE
C
C              COMPUTE SQUARED DISTANCES FROM (X(EXT),Y(EXT)).
C              ALSO INITIALIZE IP().
      DO 20005 I=1,NP
         W(I) = (X(I)-X(EXT))**2 + (Y(I)-Y(EXT))**2
20005 CONTINUE
C
C        Use _SORTP from the JPL MATH77 library to
C        determine a sorted order for the contents of W() and set
C        indices indicating the sorted order in IP().
C        IP() will be set so that W(IP(1)) is the smallest number in
C        W(), W(IP(2)) is the next larger number, etc.  The contents
C        of IP() will be used as indices into the X() and Y() arrays.
C        The contents of W() will not be reordered by _SORTP, however
C        this does not matter here since the next time we reference W()
C        will be to store other things into it.
C
      CALL DSORTP(W, 1, NP, IP)
C
C              BUILD THE FIRST TRIANGLE USING POINTS IP(1), IP(2), AND
C              THE NEXT POINT THAT IS NOT COLINEAR WITH THESE TWO.
C
      I1=IP(1)
      I2=IP(2)
      DX12= X(I2)-X(I1)
      DY12= Y(I2)-Y(I1)
      I3PREV = 0
      DO 20009 I=3,NP
      I3 = IP(I)
      IP(I) = I3PREV
         CROSS =-(X(I3)-X(I2))*DY12 + (Y(I3)-Y(I2))*DX12
      IF(CROSS .NE. 0.)THEN
      IF(CROSS .LT. 0.)THEN
C
C                                  SWAP I2 AND I3
               ISAVE=I2
               I2   =I3
               I3   =ISAVE
      END IF
      GO TO 20008
      END IF
         I3PREV=I3
20009 CONTINUE
C
C              ERROR IF DROP THRU TO HERE.  ALL NP POINTS ARE COLINEAR
C              OR ELSE POINTS I1 AND I2 ARE IDENTICAL.
C
      MODE = 2
      GO TO 999
C
C
C                   HERE THE VERTICES OF THE FIRST TRIANGLE ARE I1,I2,I3
C
C                   COMPUTE MEAN OF VERTICES OF FIRST TRIANGLE.
C
20008 XMEAN = (X(I1)+X(I2)+X(I3))/3.
      YMEAN = (Y(I1)+Y(I2)+Y(I3))/3.
C
C              COMPUTE THE PSEUDOANGLE OF VECTOR FROM (XMEAN,YMEAN)
C              TO I1.
C
      AREF = DTGANG(X(I1)-XMEAN, Y(I1)-YMEAN, 0.D0)
C
C         INITIALIZE BOUNDARY POINTER STRUCTURE B(,).
C         B(K,1) = FWD POINTER
C         B(K,2) = BACKWARD POINTER
C         B(K,3) = A BOUNDARY POINT
C         B(K,4) = A BOUNDARY TRIANGLE
C         W(K)   = PSEUDOANGLE OF POINT B(K,3)
C
      DO 20016 K=1,MB
         B(K,1) = K+1
20016 CONTINUE
20017 B(MB,1) = 0
      B(4,1)  = 0
      AVAIL   = 5
      DO 20019 K=1,4
         B(K,2) = K-1
         B(K,4) = 1
20019 CONTINUE
20020 B(4,4) = 0
C
      B(1,3) = I1
      B(2,3) = I2
      B(3,3) = I3
      B(4,3) = I1
C
      W(1) = 0.
      W(2) = DTGANG(X(I2)-XMEAN,Y(I2)-YMEAN,AREF)
      W(3) = DTGANG(X(I3)-XMEAN,Y(I3)-YMEAN,AREF)
      W(4) = 8.
      NB = 3
      MBUSED = 3
C
C                   BUILD FIRST TRIANGLE IN TRIANG().
C
      IF(NTLIM .LE. 0)THEN
         MODE = 5
         GO TO 999
      END IF
      CALL DTGSET(1,0,0,0,I1,I2,I3,TRIANG,MT)
      NT=1
C     ******************************************************************
C
C                   BEGIN MAIN LOOP.  ADD ONE NEW POINT AT A TIME INTO
C                   THE TRIANGULAR GRID STRUCTURE.
C
C     ******************************************************************
      KBNEW=2
      INDEX=4
      N20025=NP
      GO TO 20026
20024 INDEX=INDEX+1
20026 IF ((N20025-INDEX).LT.0) GO TO 20025
         JNEW = IP(INDEX)
         PANGLE = DTGANG(X(JNEW)-XMEAN, Y(JNEW)-YMEAN, AREF)
C
C                   USE  PANGLE  AS SEARCH KEY IN THE TABLE OF
C                   PSEUDO ANGLES OF BOUNDARY POINTS.
C                   START SEARCH AT PREVIOUS NEW BOUNDARY POINT.
C
         KB1= KBNEW
      IF( PANGLE .GT. W(KB1))THEN
            KBHI = KB1
      GO TO 20031
20029 IF ( PANGLE .LE. W(KBHI) ) GO TO 20030
20031          KBLO = KBHI
               KBHI = B(KBLO,1)
      GO TO 20029
20030 CONTINUE
      ELSE
C                             HERE WE HAVE  PANGLE .LE. W(KB1)
C
            KBLO = KB1
      GO TO 20034
20032 IF ( PANGLE .GE. W(KBLO) ) GO TO 20033
20034       KBHI = KBLO
               KBLO = B(KBHI,2)
      GO TO 20032
20033 CONTINUE
      END IF
C
         JLO = B(KBLO,3)
         JHI = B(KBHI,3)
C
C                             TEST FOR DUPLICATE POINTS
C
      IF( PANGLE .EQ. W(KBLO) )THEN
         IF( X(JNEW) .EQ. X(JLO) .AND. Y(JNEW).EQ.Y(JLO))THEN
               MODE = 3
               GO TO 999
         END IF
      ELSE
      IF(PANGLE .EQ. W(KBHI) )THEN
      IF(X(JNEW).EQ.X(JHI) .AND. Y(JNEW).EQ.Y(JHI))THEN
                  MODE = 3
                  GO TO 999
      END IF
      END IF
      END IF
C
C                             ATTACH POINT JNEW TO JHI AND JLOW.
C
      GO TO 30001
C
C                             THE ABOVE PROCEDURE ASSIGNS JNEW THE
C                             INDEX KBNEW IN THE BOUNDARY LIST.
C
C              LOOP BACKWARD ALONG BOUNDARY FROM KBNEW TO FIND OTHER
C              POINTS TO WHICH POINT JNEW CAN BE CONNECTED.
C
20043    KBHI= KBLO
         KBLO = B(KBHI,2)
20044 IF (.NOT.(KBLO.NE.0)) GO TO 20045
            JHI = B(KBHI,3)
            JLO = B(KBLO,3)
            CROSS = -(Y(JHI)-Y(JNEW)) * (X(JLO)-X(JHI))
     &              +(X(JHI)-X(JNEW)) * (Y(JLO)-Y(JHI))
      IF (CROSS .LE. 0.) GO TO 20045
C
C                   CONNECT POINT JNEW TO JLO
C
      GO TO 30002
20046       KBHI = KBLO
            KBLO = B(KBHI,2)
      GO TO 20044
C
C              LOOP FORWARD ALONG BOUNDARY FROM KBNEW TO FIND OTHER
C              POINTS TO WHICH POINT JNEW CAN BE CONNECTED.
C
20045    KBLO = B(KBNEW,1)
         KBHI = B(KBLO,1)
20047 IF (.NOT.(KBHI.NE.0)) GO TO 20048
            JHI= B(KBHI,3)
            JLO= B(KBLO,3)
            CROSS = -(Y(JLO)-Y(JHI)) * (X(JNEW)-X(JLO))
     &              +(X(JLO)-X(JHI)) * (Y(JNEW)-Y(JLO))
      IF (CROSS .LE. 0) GO TO 20048
C
C                                  CONNECT POINT JNEW TO JHI
C
      GO TO 30003
20049       KBLO = KBHI
            KBHI = B(KBLO,1)
      GO TO 20047
20048 GO TO 20024
C
20025 MBUSED = MBUSED + 1
C
  999 CONTINUE
      IF(MODE .EQ. 0 .AND. .NOT. OK1) MODE = 1
      INFO(1) = MODE
      INFO(2) = NB
      INFO(3) = MBUSED
      RETURN
C
C     ******************************************************************
C
C                             END MAIN LOOP
C
C     ******************************************************************
C
C     ***********************************
C                                       *
C     PROCEDURE ( ATTACH-JNEW-JHI-JLO )
C                                       *
C     ***********************************
C
30001 IF(AVAIL .EQ. 0)THEN
C
C                                  ERROR..  INSUFFICIENT STORAGE IN B()
C                                  TO REPRESENT THE BOUNDARY.
C
            MODE = 4
            GO TO 999
      END IF
      IF(NT   .GE. NTLIM)THEN
C
C                                  ERROR..  INSUFFICIENT STORAGE IN
C                                  TRIANG()
C
            MODE = 5
            GO TO 999
      END IF
C
C                                  UPDATE TRIANG().
C
         NT=NT+1
         KBFIX = KBLO
      ASSIGN 20054 TO LOC004
      GO TO 30004
20054    CALL DTGSET(NT,B(KBLO,4),0,0,JHI,JLO,JNEW,TRIANG,MT)
C
C                                  UPDATE B() AND W().
         NB = NB+1
         IF (NB .GT. MBUSED) MBUSED = NB
         KBNEW = AVAIL
         AVAIL = B(AVAIL,1)
         B(KBNEW,1) = KBHI
         B(KBNEW,2) = KBLO
         B(KBNEW,3) = JNEW
         B(KBNEW,4) = NT
         W(KBNEW)   = PANGLE
         B(KBLO,1)  = KBNEW
         B(KBHI,2)  = KBNEW
         B(KBLO,4)  = NT
C
      ASSIGN 20055 TO LOC005
      GO TO 30005
20055 GO TO 20043
C
C     *********************************
C                                     *
C     PROCEDURE ( CONNECT-JNEW-JLO )
C                                     *
C     *********************************
C
30002 IF(NT .GE. NTLIM)THEN
C
C                        ERROR..  INSUFFICIENT SPACE IN TRIANG().
C
            MODE = 5
            GO TO 999
      END IF
C
C                                  UPDATE TRIANG()
C
         NT=NT+1
         KBFIX=KBLO
      ASSIGN 20058 TO LOC004
      GO TO 30004
20058    KBFIX=KBHI
      ASSIGN 20059 TO LOC004
      GO TO 30004
20059    CALL DTGSET(NT,B(KBLO,4),0,B(KBHI,4),JHI,JLO,JNEW,TRIANG,MT)
C
C                                  UPDATE B()
         B(KBLO,1) = KBNEW
         B(KBNEW,2)= KBLO
         B(KBLO,4) = NT
         NB=NB-1
         B(KBHI,1) = AVAIL
         AVAIL     = KBHI
C
      ASSIGN 20060 TO LOC005
      GO TO 30005
20060 GO TO 20046
C
C     *********************************
C                                     *
C     PROCEDURE ( CONNECT-JNEW-JHI )
C                                     *
C     *********************************
C
30003 IF(NT .GE. NTLIM)THEN
C
C                        ERROR..  INSUFFICIENT SPACE IN TRIANG().
            MODE = 5
            GO TO 999
      END IF
C                                  UPDATE TRIANG()
C
         NT=NT+1
         KBFIX= KBLO
      ASSIGN 20063 TO LOC004
      GO TO 30004
20063    KBFIX= KBNEW
      ASSIGN 20064 TO LOC004
      GO TO 30004
20064    CALL DTGSET(NT,B(KBLO,4),B(KBNEW,4),0,JHI,JLO,JNEW,TRIANG,MT)
C
C                                  UPDATE B()
         B(KBNEW,1) = KBHI
         B(KBHI,2)  = KBNEW
         B(KBNEW,4) = NT
         NB=NB-1
         B(KBLO,1)  = AVAIL
         AVAIL      = KBLO
C
      ASSIGN 20065 TO LOC005
      GO TO 30005
20065 GO TO 20049
C
C     *
C     **********************************
C                                      *
C     PROCEDURE ( UPDATE-USING-KBFIX )
C                                      *
C     **********************************
C
30004    JFIX= B(KBFIX,3)
         CALL DTGGET(B(KBFIX,4),TRI,TRIANG)
      DO 20066 I=1,3
      IF(TRI(I+3) .EQ. JFIX)THEN
               TRI(I) = NT
      GO TO 20067
      END IF
20066 CONTINUE
20067    CALL DTGPUT(B(KBFIX,4),TRI,TRIANG,MT)
      GO TO LOC004,(20054,20058,20059,20063,20064)
C
C     ****************************
C                                *
C     PROCEDURE ( CALL-DTGADJ )
C                                *
C     ****************************
C
30005    CALL DTGADJ(NT,1,X,Y,TRIANG,MT,B,MB, KBNEW, OK)
         OK1 = OK1 .AND. OK
      GO TO LOC005,(20055,20060,20065)
C
      END SUBROUTINE DTGGRD

C***********************************************************************

      DOUBLE PRECISION FUNCTION DTGANG(XX,YY,AREF)
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C
C     Compute the pseudoangle  DTGANG  between the reference direction
C     whose pseudoangle is  ARDF  and the vector  (XX,YY).
C     Pseudoangle is measured counterclockwise with period 8.
C     Require  0. .LE. AREF .LT. 8. on input.
C     On RETURN DTGANG will be in the interval 0. .LE. DTGANG .LT. 8.
C
C     C.L.LAWSON, JPL, 1976 NOV 8
C     ------------------------------------------------------------------
      DOUBLE PRECISION xx,yy,aref,a
      DOUBLE PRECISION eight, four, two, zero
      PARAMETER(eight = 8.0d0, four = 4.0d0, two = 2.0d0, zero = 0.0d0)
C     ------------------------------------------------------------------
      IF(abs(xx) .GT. abs(yy))THEN
         a = yy/xx
      IF(xx .LT. zero)THEN
            a=a+four
      ELSE
      IF(yy .LT. zero)THEN
               a=a+ eight
      END IF
      END IF
      ELSE
      IF(yy .EQ. zero)THEN
            a = zero
      ELSE
            a = two - xx/yy
            IF(yy .LT. zero) a = a+four
      END IF
      END IF
C     .                            A IS NOW THE PSEUDOANGLE OF (XX,YY)
C     .                            RELATIVE TO (1.,0.) AND LYING IN
C     .                            THE INTERVAL   0. .LE. A .LT. 8.
      a = a - aref
      IF(a .LT. zero) a = a + eight
      DTGANG = a
      RETURN
      END FUNCTION DTGANG

C***********************************************************************

      SUBROUTINE DTGADJ(TIN,NIN,X,Y,TRIANG,MT,B,MB,KBNEW, OK)
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C
C     C.L.LAWSON, JPL, 1976 NOV 30, Changed 1977 APR 7
C     Changed 1977 June 16 to get rid of internal arrays
C     used to stack triangles to be tested.
C
C     Given triangle TIN and the neighboring triangle that is in
C     position NIN relative to TIN,  Test for possible exchange of
C     their COMMON edge.
C     Use the criterion of maximizing the smallest angle, which is the
C     same as the criterion of the empty circumcircle.
C     IF the swap is made THEN certain other edges must be tested.
C     This subr does the tests, swaps edges as needed, and does
C     tests and swaps on other edges as needed.
C     The counter COUNT is used to guard against infinite looping
C     which should never happen anyway.  IF count exceeds 100 the
C     subr sets OK = .FALSE. indicating that the triangulation
C     produced may not be optimal.
C     ------------------------------------------------------------------
C                       SUBROUTINE Arguments
C
C  TIN [INTEGER, in]
C  NIN [INTEGER, in]
C  X(),Y() [floating, in]
C  TRIANG() [INTEGER, inout]
C  MT [INTEGER, in]
C  B(,1:4) [INTEGER, inout]
C  MB [INTEGER, in]
C  KBNEW [INTEGER, in]  Index into the boundary array, B(,).
C  OK [LOGICAL, out]
C     ------------------------------------------------------------------
      INTEGER mb
      INTEGER  ADD1(3)
      INTEGER  B(MB,4), COUNT
      INTEGER i, jpt, k, kbkwd,KBNEW, kfwd
      INTEGER more, mt
      INTEGER nabor1, nabor2, nin
      INTEGER  QV1,QV2,QV3,QV4,QV1S(3),QV2S(3),QV3S(3),SUB1(3)
      INTEGER  t1, T11,T12,t2, T21,T22, tfix, tin
      INTEGER  TNEW, TOLD, TRI1(6),TRI2(6), TRIANG(1:MT)
      DOUBLE PRECISION a123, a234, a341, a412
      DOUBLE PRECISION diag13, diag24
      DOUBLE PRECISION dx12, dx23, dx34, dx41, dy12, dy23, dy34, dy41
      DOUBLE PRECISION s12, s23, s34, s41
      DOUBLE PRECISION tau13, tau24
      DOUBLE PRECISION X(1:*),Y(1:*)
      LOGICAL  OK, SWAP
      DATA QV1S(1),QV1S(2),QV1S(3)/ 5, 6, 4/
      DATA QV2S(1),QV2S(2),QV2S(3)/ 6, 4, 5/
      DATA QV3S(1),QV3S(2),QV3S(3)/ 4, 5, 6/
      DATA ADD1(1),ADD1(2),ADD1(3)/ 2, 3, 1/
      DATA SUB1(1),SUB1(2),SUB1(3)/ 3, 1, 2/
C     ------------------------------------------------------------------
      T1=TIN
      NABOR1 = NIN
      COUNT = 1
      MORE = 1
      OK = .TRUE.
20002 IF (.NOT.(MORE.GT.0.AND.OK)) GO TO 20003
      GO TO 30001
20004 IF (.NOT.( SWAP)) GO TO 20006
      GO TO 30002
20007    MORE =MORE+1
         COUNT = COUNT+1
         OK = COUNT .LE. 100
      GO TO 20005
20006    MORE = MORE-1
      IF (.NOT.(MORE .GT. 0)) GO TO 20009
      GO TO 30003
20010 CONTINUE
20009 CONTINUE
20005 GO TO 20002
20003 RETURN
C
C *****************************( TEST FOR SWAP )************************
C
C     PROCEDURE ( TEST FOR SWAP )
30001    CALL DTGGET(T1,TRI1,TRIANG)
         T2 = TRI1(NABOR1)
      IF (.NOT.( T2 .NE. 0)) GO TO 20012
C                                  T2=0 MEANS SIDE NABOR OF TRIANGLE T1
C                                  IS A BOUNDARY EDGE SO THERE IS NO
C                                  TRIANGLE T2.
C                                  HERE WE HAVE  T2 .NE. 0
            CALL DTGGET(T2,TRI2,TRIANG)
C
C                             SET QV1,...,QV4 AS INDICES OF THE 4
C                             VERTICES OF THE QUADRILATERAL FORMED BY
C                             THE 2 TRIANGLES T1 AND T2
C
            I=QV1S(NABOR1)
            QV1=TRI1(I)
            I=QV2S(NABOR1)
            QV2=TRI1(I)
            I=QV3S(NABOR1)
            QV3=TRI1(I)
      DO 20013 NABOR2=1,3
      IF(TRI2(NABOR2) .EQ. T1)THEN
                  I=QV2S(NABOR2)
                  QV4=TRI2(I)
      GO TO 20014
      END IF
20013 CONTINUE
20014 GO TO 30004
20018 GO TO 20011
20012     SWAP = .FALSE.
20011 GO TO 20004
C
C*****************************( SWAP THE EDGE )*************************
C
C     PROCEDURE ( SWAP THE EDGE )
C
C                             REMOVE THE EDGE CONNECTING
C                             QV1 AND QV3.  INSERT NEW EDGE CONNECTING
C                             QV2 AND QV4.
C                             REDEFINE TRIANGLE T1 TO HAVE VERTICES QV1,
C                             QV2, AND QV4.  REDEFINE TRIANGLE T2 TO
C                             HAVE VERTICES QV2, QV3, AND QV4.
C
30002          I = ADD1(NABOR1)
               T11 = TRI1(I)
               I = ADD1(I)
               T12 = TRI1(I)
C
               I   = ADD1(NABOR2)
               T21 = TRI2(I)
               I   = ADD1(I)
               T22 = TRI2(I)
C
               CALL DTGSET(T1,T11,T2,T22,QV1,QV2,QV4,TRIANG,MT)
               CALL DTGSET(T2,T12,T21,T1,QV2,QV3,QV4,TRIANG,MT)
C
               TNEW=T2
      IF (.NOT.( T12 .NE. 0)) GO TO 20020
                  TFIX= T12
                  TOLD= T1
      ASSIGN 20021 TO LOC005
      GO TO 30005
20021 GO TO 20019
C
C           The following test of (QV2 .EQ. B(KBNEW,3)) will
C           always be TRUE IF the outer algorithm calling this
C           SUBROUTINE is the one originally planned,
C           always adding new points from
C           the exterior, not from the interior.  The FALSE
C           branch is provided to handle the more general case.
C
20020 IF (.NOT.(QV2 .EQ. B(KBNEW,3) )) GO TO 20023
                     B(KBNEW,4) = TNEW
      GO TO 20022
C1001             FORMAT(/' ADJUST..  USING FIX-B WITH QV2 =',I5,5X,
C    *            'B(KBNEW,3) =',I5)
C                 WRITE(*,1001) QV2,B(KBNEW,3)
20023                JPT = QV2
      ASSIGN 20024 TO LOC006
      GO TO 30006
20024 CONTINUE
20022 CONTINUE
C
20019          TNEW=T1
      IF (.NOT.( T22 .NE. 0)) GO TO 20026
            TFIX=T22
            TOLD = T2
            ASSIGN 20027 TO LOC005
            GO TO 30005
20027       GO TO 20025
C1002       FORMAT(/' ADJUST..  USING FIX-B WITH QV4 =',I5)
C           WRITE(*,1002) QV4
20026    JPT = QV4
         ASSIGN 20028 TO LOC006
         GO TO 30006
20028 CONTINUE
20025         NABOR1 = 3
      GO TO 20007
C
C********************************( FIX-TFIX )***************************
C
C     PROCEDURE ( FIX-TFIX )
C
30005    CALL DTGGET( TFIX,TRI2, TRIANG)
      DO 20029 I=1,3
      IF(TRI2(I) .EQ. TOLD)THEN
               TRI2(I) = TNEW
      GO TO 20030
      END IF
20029 CONTINUE
20030    CALL DTGPUT( TFIX,TRI2, TRIANG, MT)
      GO TO LOC005,(20021,20027)
C
C********************************( FIX-B )******************************
C
C     PROCEDURE ( FIX-B )
C
C         SEARCH ALTERNATELY FWD AND BKWD FROM KBNEW IN THE BOUNDARY
C         LIST B() FOR THE POINT JPT.   WHEN FOUND SET THE ASSOCIATED
C         TRIANGLE TO TNEW.
C
30006    KFWD = KBNEW
         KBKWD= B(KFWD ,2)
20034 IF(KFWD .NE. 0)THEN
      IF(B(KFWD,3) .EQ. JPT)THEN
                  B(KFWD,4) = TNEW
      GO TO 20035
      END IF
               KFWD = B(KFWD,1)
      END IF
C
      IF(KBKWD .NE. 0)THEN
      IF(B(KBKWD,3) .EQ. JPT)THEN
                  B(KBKWD,4) = TNEW
      GO TO 20035
      END IF
               KBKWD = B(KBKWD,2)
      END IF
      GO TO 20034
20035 GO TO LOC006,(20024,20028)
C
C********************************( MAXMIN ANGLE TEST )******************
C
C     PROCEDURE ( MAXMIN ANGLE TEST )
C
C         THE TEST IS MADE USING THE MAX-MIN ANGLE CRITERION AS
C     DESCRIBED IN JPL INTERNAL TECH. MEMO. 914-299, FEB., 1972
C     BY C. L. LAWSON.
C
30004 DX23 = X(QV3)-X(QV2)
      DY23 = Y(QV3)-Y(QV2)
      DX34 = X(QV4)-X(QV3)
      DY34 = Y(QV4)-Y(QV3)
      A234 = DX23*DY34 - DY23*DX34
      IF(A234 .LE. 0.)THEN
         SWAP = .FALSE.
      ELSE
         DX12 = X(QV2)-X(QV1)
         DY12 = Y(QV2)-Y(QV1)
         DX41 = X(QV1)-X(QV4)
         DY41 = Y(QV1)-Y(QV4)
         A412 = DX41*DY12 - DY41*DX12
      IF(A412 .LE. 0.)THEN
            SWAP = .FALSE.
      ELSE
            A341 = DX34*DY41 - DY34*DX41
            A123 = A234 + A412 - A341
            S12 = DX12**2 + DY12**2
            S23 = DX23**2 + DY23**2
            S34 = DX34**2 + DY34**2
            S41 = DX41**2 + DY41**2
            DIAG13 = (X(QV3)-X(QV1))**2 + (Y(QV3)-Y(QV1))**2
            DIAG24 = (X(QV4)-X(QV2))**2 + (Y(QV4)-Y(QV2))**2
            TAU13 = MIN( A123**2 / MAX(S12,S23),
     &                     A341**2 / MAX(S34,S41)) / DIAG13
            TAU24 = MIN( A234**2 / MAX(S23,S34),
     &                     A412**2 / MAX(S41,S12)) / DIAG24
            SWAP = TAU24 .GT. TAU13
      END IF
      END IF
      GO TO 20018
C
C***************************( SET NEXT T1 AND NABOR1 )******************
C
C     PROCEDURE ( SET NEXT T1 AND NABOR1 )
30003 TOLD = T1
      I = SUB1(NABOR1)
      T1=TRI1(I)
      CALL DTGGET(T1,TRI1,TRIANG)
      DO 20048 K=1,3
      IF(TRI1(K) .EQ. TOLD)THEN
            NABOR1 = SUB1(K)
      GO TO 20049
      END IF
20048 CONTINUE
20049 GO TO 20010
C
      END SUBROUTINE DTGADJ

C***********************************************************************

C  File: DTGPD.[f|for] contains DTGPD, DTGMOR, DTGLS.
      SUBROUTINE DTGPD(X, Y, Z, DZ, NP, TRIANG, NT, IWORK)
C>> 1996-05-10 DTGPD Krogh  MW removed from dim. declarator in C version
C>> 1996-03-30 DTGPD Krogh  MIN0 => MIN
C>> 1996-02-02 CLL
C>> 1996-01-15 CLL
C>> 1996-01-11 CLL
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C
C  THIS SUBR ESTIMATES FIRST PARTIAL DERIVS AT THE GIVEN DATA POINTS.
C  C.L.LAWSON, JPL, 1976 DEC 21.  EDITED COMMENTS 1979 MAR 5.
C
C  ---------------------------------------------------------------------
C                     Method
C
C      THE PARTIAL DERIVS AT POINT  JP  WILL BE ESTIMATED BY FITTING
C  A QUADRATIC IN X AND Y TO A SET OF  MFIT+1  POINTS CONSISTING OF
C  POINT  JP  AND   MFIT NEARBY POINTS IN THE TRIANGULAR GRID.
C         THE VALUE OF MFIT WILL IN GENERAL BE DIFFERENT AT DIFFERENT
C  POINTS.  MFIT WILL BE IN THE RANGE FROM MFMIN TO MFMAX.  LET NNEB
C  DENOTE THE NO. OF IMMEDIATE NEIGHBORING POINTS TO POINT JB.  IF
C  MFMIN .LE. NNEB .LE. MFMAX  THEN MFIT=NNEB.
C      IF NNEB .LT. MFMIN THEN MFIT = MFMIN.
C      IF NNEB .GT. MFMAX THEN MFIT = MFMAX.
C  ---------------------------------------------------------------------
C                 SUBROUTINE Arguments
C
C  X(), Y(), Z()     [in]
C  DZ(2,)            [out]
C  NP                [in]
C  TRIANG(), NT  [in]
C  (IWORK(I),I=1,NP) [scratch]  INITIALLY SET TO ZERO BY DTGPD.
C       DTGPD WILL SET IWORK(JP)= 1 WHEN IT ESTIMATES THE PARTIAL
C       DERIVATIVES AT POINT JP.
C
C  ---------------------------------------------------------------------
C                     COMMON BLOCK
C
C  /DTGCM1/ provides communication between DTGPD and DTGMOR.
C  ---------------------------------------------------------------------
C                 Internal variables
C
C  AVAIL  POINTS TO AVAILABLE SPACE IN KE().
C  MFMAX  Dimensioning PARAMETER used in the following declarations:
C           INTEGER   KE(MFMAX  ,4), JUSED(MFMAX)
C           DOUBLE PRECISION W(MFMAX+5,6), DSQ(MFMAX)
C  MFMIN  LOWER BOUND FOR MFIT.  SET TO A CONSTANT VALUE BY DTGPD.
C  MFIT   NO. OF NEIGHBORING POINTS USED WITH EACH POINT JP
C         TO ESTIMATE PARTIALS AT POINT JP.  GENERALLY DIFFERS
C         FOR EACH POINT JP.
C  FIRSTL,LASTL     POINTERS TO FIRST AND LAST ELEMENTS OF THE
C                   LIST CONTAINED IN KE().
C  KE(,) IS A DOUBLY  LINKED LIST DESCRIBING A CYCLE OF EDGES
C       SURROUNDING POINT JP.  AN EDGE IS DEFINED BY ITS TWO ADJACENT
C       TRIANGLES.
C           KE(L,1)   FWD POINTER  (COUNTERCLOCKWISE).
C           KE(L,2)   BKWD POINTER (CLOCKWISE).
C           KE(L,3)   TRIANGLE INSIDE EDGE L.
C           KE(L,4)   TRIANGLE OUTSIDE EDGE L,  OR ELSE ZERO IF EDGE
C                     ON THE BOUNDARY.
C  JUSED()   AN UNORDERED SET OF INDICES OF ROINTS USED IN THE LOCAL
C            FIT THAT ARE NOT IMMEDIATE NEIGHBORS OF POINT JP.
C  X0,Y0,Z0  COORDINATES OF POINT JP.
C  DSQ(L)    SQUARED DISTANCE FROM POINT JP TO THE MIDPOINT OF SIDE L
C  NEWPT     LOGICAL FLAG.  FOR A NEW POINT JP DTGPD SETS
C            NEWPT TRUE.  THIS CAUSES DTGLS TO TRIANGULARIZE THE
C            LEAST SQUARES SYSTEM FOR JP. DTGLS THEN TESTS THE
C            CONDITION OF THE SYSTEM. IF IT IS BAD DTGLS SETS NEWPT
C            FALSE AND RETURNS. ESTPT THEN USES DTGMOR TO GET ONE
C            MORE POINT AND CALLS DTGLS AGAIN, LEAVING NEWPT
C            FALSE.  DTGLS THEN ACCUMULATES THE ADDITIONAL POINT
C            AND AGAIN TESTS THE CONDITION.  EVENTUALLY EITHER
C            ENOUGH POINTS ARE ADDED TO OBTAIN ADEQUATE CONDITION
C            OR ELSE DTGLS APPLIES AN ARBITRARY STABILIZATION AND
C            SO IN EITHER CASE DTGLS FINALLY RETURNS A SOLUTION.
C
C     ------------------------------------------------------------------
C--D replaces "?": ?TGPD, ?TGGET, ?TGMOR, ?TGLS, ?TGCM1
C--&               ?TGMOR, ?TGCM1, ?TGLS,  ?ROT,  ?ROTG
C     ------------------------------------------------------------------
      COMMON  /DTGCM1/ X0, Y0, Z0, DSQ,
     &         AVAIL, MFMIN, MFIT, FIRSTL, LASTL, KE, JUSED, NEWPT
C Note -- IF nfmax is changed, mw changes and C code must change
C everywhere mw is used in a dimension declarator.
      INTEGER mfmax
      PARAMETER(mfmax = 16)
      DOUBLE PRECISION X0, Y0, Z0, DSQ(mfmax)
      INTEGER AVAIL,MFMIN,MFIT, FIRSTL, LASTL, KE(mfmax,4), JUSED(mfmax)
      LOGICAL NEWPT
      INTEGER mw
      PARAMETER(mw = mfmax+5)
      INTEGER NP
      INTEGER ADD1(3), i, IWORK(np), ip1, ipcent, it, iv1, iv2, ivert
      INTEGER j1, j2, jj1, jj2, jp
      INTEGER kp, l, LIMIT
      INTEGER LOC003, LOC004
      INTEGER lolim
      INTEGER n20006
      INTEGER NT
      INTEGER SUB1(3)
      INTEGER TNOW, TRI(6), TRIANG(1:*), TRSAVE(6)
C++ Code for .C. is INACTIVE
C      DOUBLE PRECISION dsfcn, DZ(2,NP), W(21,6)
C++ Code for ~.C. is ACTIVE
      DOUBLE PRECISION dsfcn, DZ(2,NP), W(mw,6)
C++ END
      DOUBLE PRECISION X(NP), x1, x2, Y(NP), y1, y2, Z(NP)
      LOGICAL HITBDY
C
      DATA ADD1(1),ADD1(2),ADD1(3),SUB1(1),SUB1(2),SUB1(3)/2,3,1, 3,1,2/
C     ------------------------------------------------------------------
      DSFCN(X1,Y1,X2,Y2) = ((X1+X2)*.5 - X0)**2 + ((Y1+Y2)*.5 - Y0)**2
C     ------------------------------------------------------------------
      MFMIN = 6
      LIMIT = min( MFMAX, NP-1 )
      LOLIM = min( MFMIN, NP-1 )
      NEWPT = .TRUE.
C                        SET IWORK() = 0
      DO 20002 I=1,NP
         IWORK(I)=0
20002 CONTINUE
C                                  MAIN LOOP THRU TRIANGLES
C
      IT=1
      N20006=NT
      GO TO 20007
20005 IT=IT+1
20007 IF ((N20006-IT).LT.0) GO TO 20006
         CALL DTGGET(IT,TRSAVE, TRIANG)
C                                       LOOP THRU VERTICES OF A TRIANGLE
      IVERT=1
      GO TO 20010
20008 IVERT=IVERT+1
20010 IF ((3-IVERT).LT.0) GO TO 20009
            JP= TRSAVE(IVERT+3)
      IF (.NOT.(IWORK(JP) .EQ. 0)) GO TO 20012
               IWORK(JP)=1
C1001 FORMAT(39H0ESTPD..  ESTIMATING PARTIALS AT POINT ,I5/1X)
C     THE FORTRAN END FOLLOWS..
C              WRITE(*,1001) JP
      GO TO 30001
20013 CONTINUE
20012 GO TO 20008
20009 GO TO 20005
20006 RETURN
C                        END OF TOP LEVEL CODE.   PROCEDURES FOLLOW..
C
C     *************************************
C                                         *
C     PROCEDURE ( PARTIALS FOR POINT JP )
C                                         *
C      ************************************
C
30001    X0=X(JP)
         Y0=Y(JP)
         Z0=Z(JP)
C                                  SET POINTERS FOR AVAILABLE SPACE.
      DO 20014 I=2,MFMAX
            KE(I-1,1)=I
20014 CONTINUE
         AVAIL=1
         KE(MFMAX,1)=0
C                                            BUILD FITTING EQUATIONS.
      GO TO 30002
20017 IF(MFIT .LT. LOLIM)THEN
            CALL DTGMOR ( X,Y,Z,NP, TRIANG,W)
      END IF
C                                            SOLVE FITTING EQUATIONS.
C
         CALL DTGLS(W,MW,MFIT,NEWPT, LIMIT, DZ(1,JP),DZ(2,JP))
C                     GO FOR MORE POINTS IF NEEDED TO IMPROVE CONDITION.
20020 IF(.NOT.NEWPT)THEN
           CALL DTGMOR (X,Y,Z,NP, TRIANG,W)
           CALL DTGLS(W,MW,MFIT,NEWPT, LIMIT, DZ(1,JP),DZ(2,JP))
      GO TO 20020
      END IF
C
      GO TO 20013
C
C     *****************************************
C                                             *
C     PROCEDURE ( PROCESS RING OF NEIGHBORS )
C                                             *
C      ****************************************
30002 DO 20022 I=1,6
            TRI(I)=TRSAVE(I)
20022 CONTINUE
      ASSIGN 20025 TO LOC003
      GO TO 30003
C                             BUILD FIRST EDGE IN RING OF NEIGHBORS.
20025    FIRSTL=L
         LASTL =L
         KE(L,1)=0
         KE(L,2)=0
         KE(L,3)=IT
         IV1=   ADD1(IVERT)
         KE(L,4)= TRI(IV1)
         J1 = TRI(IV1+3)
         IV2=   SUB1(IVERT)
         J2 = TRI(IV2+3)
         DSQ(L) = DSFCN(X(J1),Y(J1),X(J2),Y(J2))
C
C                            BUILD EQUAS FOR FIRST TWO NEIGHBORING PTS.
         MFIT=1
         KP=J1
      ASSIGN 20026 TO LOC004
      GO TO 30004
20026    MFIT=2
         KP=J2
      ASSIGN 20027 TO LOC004
      GO TO 30004
C
C                             MOVE LEFT AND THEN RIGHT AROUND POINT JP
C                             TO BUILD EQUAS FOR IMMEDIATELY
C                             NEIGHBORING POINTS.
C
20027 GO TO 30005
20028 IF (.NOT.( HITBDY)) GO TO 20030
      DO 20031 I=1,6
               TRI(I)=TRSAVE(I)
20031 CONTINUE
      GO TO 30006
20034 CONTINUE
20030 GO TO 20017
C     ***********************
C                           *
C     PROCEDURE (BUILD ROW)
C                           *
C      **********************
C1000 FORMAT(25H BUILD ROW..  CENTER PT =,I5,12H,  EQUA NO =,I3,
C    *13H  NEARBY PT =,I5)
C        WRITE(*,1000) JP,MFIT,KP
30004    JUSED(MFIT) = KP
         W(MFIT,4) = X(KP) - X0
         W(MFIT,5) = Y(KP) - Y0
         W(MFIT,6) = Z(KP) - Z0
         W(MFIT,1) = W(MFIT,4)**2
         W(MFIT,2) = W(MFIT,4)*W(MFIT,5)
         W(MFIT,3) = W(MFIT,5)**2
      GO TO LOC004,(20026,20027,20045,20052)
C
C     **************************
C                              *
C     PROCEDURE ( MOVE LEFT )
C                              *
C      *************************
C                                     LEFT MEANS COUNTERCLOCKWISE.
30005    JJ2 =  J2
         IPCENT=IVERT
         HITBDY = .FALSE.
20035 IF (.NOT.(MFIT.LT.MFMAX)) GO TO 20036
            I     = SUB1(IPCENT)
            TNOW = TRI(I)
      IF (.NOT.(TNOW .EQ. 0)) GO TO 20038
C                                            HAVE HIT BOUNDARY
               HITBDY  = .TRUE.
      GO TO 20036
20038          CALL DTGGET(TNOW,TRI, TRIANG)
      DO 20039 IPCENT=1,3
      IF ( TRI(IPCENT+3).EQ.JP) GO TO 20040
20039 CONTINUE
C
C                                       IDENTIFY NEXT POINT
20040          I = SUB1(IPCENT)
               KP= TRI(I+3)
C
C                   RECORD INFORMATION ABOUT  THE NEW EDGE IN
C                   THE LIST STRUCTURE.
C
      ASSIGN 20042 TO LOC003
      GO TO 30003
20042          KE(L,2)= LASTL
               KE(L,3)= TNOW
               I      = ADD1(IPCENT)
               KE(L,4)= TRI(I)
               KE(LASTL,1) = L
               LASTL  = L
               JJ1 = JJ2
               JJ2 = KP
               DSQ(L) = DSFCN(X(JJ1),Y(JJ1),X(JJ2),Y(JJ2) )
      IF (.NOT.( KP .NE. J1)) GO TO 20044
C                                       KP IS A NEW NEIGHBORING POINT.
                  MFIT = MFIT +1
      ASSIGN 20045 TO LOC004
      GO TO 30004
20045             KE(L,1) = 0
      GO TO 20043
C
C                   KP IS NOT A NEW POINT.  WE HAVE CYCLED ALL THE WAY
C                   AROUND TO THE INITIAL POINT J1.  THIS MEANS THE
C                   ENTIRE RING OF NEIGHBORS HAS BEEN FOUND.
C
20044             KE(L,1) = FIRSTL
                  KE(FIRSTL,2) = LASTL
      GO TO 20036
20043 CONTINUE
20037 GO TO 20035
C
C         ABOVE LOOP IS EXITED FOR ONE OF THREE REASONS..
C              (1) MFIT .EQ. MFMAX,
C              (2) HITBDY = .TRUE.,
C           OR (3) COMPLETE RING OF NEIGHBORS HAS BEEN BUILT
C
20036 GO TO 20028
C
C     **************************
C                              *
C     PROCEDURE ( MOVE RIGHT )
C                              *
C      *************************
C                                      RIGHT MEANS CLOCKWISE.
30006    JJ1 = J1
         IPCENT = IVERT
20046 IF (.NOT.(MFIT.LT.MFMAX)) GO TO 20047
            TNOW = TRI(IPCENT)
C
C              THE FOLLOWING TEST CAUSES AN EXIT FROM THE LOOP IF THE
C              BOUNDARY ON THE RIGHT HAS BEEN ENCOUNTERED.
C
      IF ( TNOW .EQ. 0 ) GO TO 20047
            CALL DTGGET(TNOW,TRI, TRIANG)
      DO 20048 IPCENT=1,3
      IF ( TRI(IPCENT+3) .EQ. JP) GO TO 20049
20048 CONTINUE
C                                       IDENTIFY NEXT POINT
20049       IP1=ADD1(IPCENT)
            KP = TRI(IP1+3)
C
C                   RECORD INFO ABOUT THE NEW EDGE IN THE LIST STRUCTURE
C
      ASSIGN 20051 TO LOC003
      GO TO 30003
20051       KE(L,1) = FIRSTL
            KE(L,2) = 0
            KE(L,3) = TNOW
            KE(L,4) = TRI(IP1)
            KE(FIRSTL,2) = L
            FIRSTL  = L
            JJ2 = JJ1
            JJ1 = KP
            DSQ(L) = DSFCN(X(JJ1),Y(JJ1),X(JJ2),Y(JJ2) )
C
C                                       KP IS A NEW NEIGHBORING POINT.
            MFIT = MFIT+1
      ASSIGN 20052 TO LOC004
      GO TO 30004
20052 GO TO 20046
C
C                   ABOVE LOOP IS EXITED WHEN
C                      (1) MFIT .EQ. MFMAX
C                   OR (2) HIT BOUNDARY ON RIGHT
C
20047 GO TO 20034
C
C     *****************************
C                                 *
C     PROCEDURE ( POP AVAIL TO L )
C                                 *
C      ****************************
C
30003    L = AVAIL
         AVAIL = KE(AVAIL,1)
      GO TO LOC003,(20025,20042,20051)
C
      END SUBROUTINE DTGPD

C***********************************************************************

      SUBROUTINE DTGMOR(X,Y,Z,NP,TRIANG,W)
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C
C     C.L.LAWSON, JPL, 1977 JAN 4
C         CHANGED NOV 30, 1977    TT
C         THIS SUBR IS CALLED BY DTGPD AS FOLLOWS:
C      (1) IF NEWPT = .TRUE.,
C            THE SUBR EXTENDS THE RING OF NEARBY POINTS AROUND
C          POINT (X0,Y0) UP TO A TOTAL OF MIN(MFMIN,NP-1) NEARBY POINTS
C          TO BE USED IN COMPUTING A LOCAL FIT TO THE SURFACE NEAR POINT
C          (X0,Y0).
C      (2) IF NEWPT = .FALSE.,
C            THE SUBR GETS ONLY ONE MORE NEARBY POINT AND
C          RETURNS.
C     ------------------------------------------------------------------
      COMMON  /DTGCM1/ X0, Y0, Z0, DSQ,
     &         AVAIL, MFMIN, MFIT, FIRSTL, LASTL, KE, JUSED, NEWPT
      INTEGER mfmax
      PARAMETER(mfmax = 16)
      DOUBLE PRECISION X0, Y0, Z0, DSQ(mfmax)
      INTEGER AVAIL,MFMIN,MFIT, FIRSTL, LASTL, KE(mfmax,4), JUSED(mfmax)
      LOGICAL NEWPT
C
      LOGICAL TSTLP1,TSTLP2,TSTLM1,TSTLM2,AVOK, TEST, USED
      INTEGER ADD1(3), i, im1, ip1, k, kp
      INTEGER l, limit2, lm1, lm2, lmin, lnew
      INTEGER LOC001, LOC002
      INTEGER lp1, lp2, lp3
      INTEGER NP
      INTEGER SUB1(3),TNOW, TRI(6),TRIANG(1:*),V1,V2
C++ Code for .C. is INACTIVE
C      DOUBLE PRECISION dmin, dsfcn, W(21,6)
C++ Code for ~.C. is ACTIVE
      INTEGER mw
      PARAMETER(mw = mfmax+5)
      DOUBLE PRECISION dmin, dsfcn, W(mw,6)
C++ END
      DOUBLE PRECISION X(NP), x1, x2, Y(NP), y1, y2, Z(NP)
C
      DATA ADD1(1),ADD1(2),ADD1(3),SUB1(1),SUB1(2),SUB1(3)/2,3,1, 3,1,2/
C
C     NAMELIST/DEBUG/AVAIL,MFMIN,MFIT,FIRSTL,LASTL, JUSED,
C    *X0,Y0,Z0,DSQ
C     ------------------------------------------------------------------
      DSFCN(X1,Y1,X2,Y2) = ((X1+X2)*.5 - X0)**2 + ((Y1+Y2)*.5 - Y0)**2
C     ------------------------------------------------------------------
C
      IF(NEWPT)THEN
         LIMIT2 = MIN( MFMIN, NP-1 )
      ELSE
         LIMIT2 = MFIT + 1
C                        WHEN DTGMOR IS CALLED WITH
C                        NEWPT = .FALSE. MFIT MUST BE LESS THAN MFMIN
C                        AND NP-1 SO THAT LIMIT2 CAN SAFELY BE SET
C                        TO MFIT+1 AS ABOVE.
      END IF
      AVOK = .TRUE.
20004 IF (.NOT.(MFIT.LT.LIMIT2.AND.AVOK)) GO TO 20005
C        WRITE(*,DEBUG)
C1004    FORMAT(/' KE(,)='/(7X,5I5))
C        WRITE(*,1004) (I,(KE(I,J), J=1,4), I=1,16)
C
C                   FIND EDGE WHOSE MIDPOINT IS CLOSEST TO POINT JP.
C                   SKIP BOUNDARY EDGES.
         L = FIRSTL
         LMIN=0
20006 IF( KE(L,4) .NE. 0)THEN
      IF(LMIN .EQ. 0)THEN
                  TEST = .TRUE.
      ELSE
                  TEST=DSQ(L) .LT. DMIN
      END IF
      IF(TEST)THEN
                  DMIN = DSQ(L)
                  LMIN = L
      END IF
      END IF
      IF (L .EQ. LASTL) GO TO 20007
            L=KE(L,1)
      GO TO 20006
20007    L=LMIN
C
C         EDGE L IS THE NEAREST TO POINT JP.
C         INVESTIGATE TRIANGLE ON OPPOSIDE SIDE OF EDGE L.
C
         TNOW = KE(L,4)
         LP1 = KE(L,1)
         TSTLP1 = LP1 .NE. 0
         IF(TSTLP1) TSTLP1 = KE(LP1,4) .EQ. TNOW
         LM1 = KE(L,2)
         TSTLM1 = LM1 .NE. 0
         IF(TSTLM1) TSTLM1 = KE(LM1,4) .EQ. TNOW
      IF (.NOT.(TSTLP1)) GO TO 20015
            LP2 = KE(LP1,1)
            TSTLP2 = LP2 .NE. 0
            IF(TSTLP2) TSTLP2 = KE(LP2,4) .EQ. TNOW
      IF (.NOT.(TSTLP2)) GO TO 20017
      ASSIGN 20018 TO LOC001
      GO TO 30001
20018 GO TO 20016
20017 IF (.NOT.(TSTLM1)) GO TO 20019
               LP2=LP1
               LP1=L
               L  =LM1
      ASSIGN 20020 TO LOC001
      GO TO 30001
20020 GO TO 20016
20019 ASSIGN 20021 TO LOC002
      GO TO 30002
20021 CONTINUE
20016 GO TO 20014
20015 IF (.NOT.(TSTLM1)) GO TO 20022
            LM2 = KE(LM1,2)
            TSTLM2 = LM2 .NE. 0
            IF(TSTLM2) TSTLM2 = KE(LM2,4) .EQ. TNOW
      IF (.NOT.(TSTLM2)) GO TO 20024
               LP2 = L
               LP1 = LM1
               L   = LM2
      ASSIGN 20025 TO LOC001
      GO TO 30001
20025 GO TO 20023
20024          LP1=L
               L  =LM1
      ASSIGN 20026 TO LOC002
      GO TO 30002
20026 CONTINUE
20023 GO TO 20014
20022 GO TO 30003
20027 GO TO 30004
20028 CONTINUE
20014 GO TO 20004
20005 RETURN
C                        END OF TOP LEVEL CODE.  PROCEDURES FOLLOW..
C
C     *******************************************
C                                               *
C     PROCEDURE ( DELETE EDGES L,LP1, AND LP2 )
C                                               *
C      ******************************************
C
30001    LM1 = KE(L,2)
         LP3 = KE(LP2,1)
         IF(LM1 .NE. 0) KE(LM1,1) = LP3
         IF(LP3 .NE. 0) KE(LP3,2) = LM1
      IF( L .EQ. FIRSTL)THEN
            FIRSTL = LP3
      ELSEIF( LP1 .EQ. FIRSTL  .OR.  LP2 .EQ. FIRSTL)THEN
            FIRSTL = LP3
            LASTL  = LM1
      ELSEIF( LP2 .EQ. LASTL)THEN
            LASTL = LM1
      END IF
C                   TRANSFER L, LP1, AND LP2 TO AVAILABLE SPACE LIST
C
         KE(LP2,1) = AVAIL
         AVAIL = L
 1002    FORMAT(60X,'MORPTS..  DELETING THREE EDGES.',
     &   '  TNOW=', I5)
         WRITE(*,1002) TNOW
      GO TO LOC001,(20018,20020,20025)
C
C     ****************************************************
C                                                        *
C     PROCEDURE ( DELETE EDGES L AND LP1.  ADD NEW EDGE)
C                                                        *
C      ***************************************************
C
30002    LP2 = KE(LP1,1)
         KE(L,1) = LP2
         IF(LP2 .NE. 0) KE(LP2,2)=L
C
      IF(LP1 .EQ.  FIRSTL)THEN
            FIRSTL = LP2
      ELSEIF(LP1 .EQ. LASTL)THEN
            LASTL = L
      END IF
C
         KE(LP1,1) = AVAIL
         AVAIL = LP1
         CALL DTGGET( TNOW,TRI, TRIANG)
C
      DO 20036 I=1,3
      IF ( TRI(I) .EQ. KE(L,3)) GO TO 20037
20036 CONTINUE
20037    KE(L,3) = TNOW
         IP1 = ADD1(I)
         KE(L,4) = TRI(IP1)
      IF(KE(L,4) .NE. 0)THEN
            IM1 = SUB1(I)
            V1  = TRI(IP1+3)
            V2  = TRI(IM1+3)
            DSQ(L) = DSFCN(X(V1),Y(V1),  X(V2),Y(V2) )
      END IF
C1003    FORMAT(60X,'MORPTS..  DELETING TWO EDGES, ADDING ONE.',
C    *   ' TNOW=', I5)
C        WRITE(*,1003) TNOW
      GO TO LOC002,(20021,20026)
C
C     *************************************************************
C                                                                 *
C     PROCEDURE (DELETE EDGE L.  ADD TWO NEW EDGES. /KP,AVOK/ )
C                                                                 *
C      ************************************************************
C
30003 GO TO 30005
C
C
20041 IF(AVOK)THEN
            KE(L,1) = LNEW
            KE(LNEW,1) = LP1
            KE(LNEW,2) = L
            KE(LNEW,3) = TNOW
            IF(LP1 .NE. 0) KE(LP1,2) = LNEW
            IF( L .EQ. LASTL) LASTL = LNEW
            CALL DTGGET( TNOW,TRI, TRIANG)
      DO 20044 I=1,3
      IF ( TRI(I) .EQ. KE(L,3) ) GO TO 20045
20044 CONTINUE
20045       IP1 = ADD1(I)
            KE(L,4) = TRI(IP1)
            KE(L,3) = TNOW
            IM1 = SUB1(I)
            KE(LNEW,4) = TRI(IM1)
            KP         = TRI(IM1+3)
            V1=TRI(IP1+3)
            DSQ(L)    = DSFCN( X(V1),Y(V1), X(KP),Y(KP) )
            V1=TRI(I+3)
            DSQ(LNEW) = DSFCN( X(V1),Y(V1), X(KP),Y(KP) )
      END IF
      GO TO 20027
C
C     **********************************
C                                      *
C     PROCEDURE ( PROCESS NEW VERTEX )
C                                      *
C      *********************************
C
C         HERE THE TRIANGLE TNOW HAS ONLY EDGE L IN COMMON WITH THE
C         CURRENT RING  AROUND POINT JP.  THE VERTEX OPPOSITE SIDE L
C         IN TRIANGLE TNOW IS KP.
C            TEST TO SEE IF POINT KP HAS ALREADY BEEN USED TO
C         BUILD AN EQUATION.
C
30004    USED = .FALSE.
      DO 20047 K=1,MFIT
      IF( JUSED(K) .EQ. KP)THEN
               USED = .TRUE.
      GO TO 20048
      END IF
20047 CONTINUE
C
20048 IF (.NOT.(.NOT. USED)) GO TO 20053
            MFIT = MFIT+1
      GO TO 30006
20054 CONTINUE
20053 GO TO 20028
C
C     ******************************************
C                                              *
C     PROCEDURE ( POP AVAIL TO LNEW. SET AVOK )
C                                              *
C      *****************************************
30005 IF( AVAIL .EQ. 0)THEN
            AVOK = .FALSE.
C                                            ERROR CONDITION
C
 1001    FORMAT(/' SUBR DTGMOR..  WARNING..  DIMENSION OF KE(,) NOT',
     &   ' LARGE ENOUGH TO FIND MFMIN NEARBY POINTS FOR ESTIMATING',
     &   ' PARTIAL DERIVS.'/
     &   ' PROGRAM WILL CONTINUE USING FEWER NEARBY POINTS.')
         WRITE(*,1001)
      ELSE
         LNEW = AVAIL
         AVAIL = KE(AVAIL,1)
      END IF
      GO TO 20041
C
C     ***********************
C                           *
C     PROCEDURE (BUILD ROW)
C                           *
C      **********************
C1000    FORMAT(' BUILD ROW IN DTGMOR..   ',5X,''   EQUA NO =',I3,
C    *   '  NEARBY PT =',I5)
C        WRITE(*,1000)    MFIT,KP
30006    JUSED(MFIT) = KP
         W(MFIT,4) = X(KP) - X0
         W(MFIT,5) = Y(KP) - Y0
         W(MFIT,6) = Z(KP) - Z0
         W(MFIT,1) = W(MFIT,4)**2
         W(MFIT,2) = W(MFIT,4)*W(MFIT,5)
         W(MFIT,3) = W(MFIT,5)**2
      GO TO 20054
C
      END SUBROUTINE DTGMOR

C***********************************************************************

      SUBROUTINE DTGLS (W,MW,MFIT,NEWPT, LIMIT, DZ1,DZ2)
C>> 1995-09-26 CLL Editing for inclusion into MATH77.
C
C     WEIGHT, SCALE, STABILIZE, AND SOLVE A LEAST SQRS PROBLEM.
C
C     C.L.LAWSON, JPL, 1976 FEB 27
C         CHANGED NOV 28, 1977    TT
C     EDITED CODE. NO CHANGE IN ALGORITHM. CLL. 1979 MAR 5.
C
C     ON INPUT, WITH NEWPT TRUE, COLUMNS OF W(,) CONTAIN TRANSLATED
C     VALUES OF X**2, X*Y, Y**2, X, Y, Z.
C
C          THE NUMBER, STAB, IS USED AS A MARQUARDT STABILIZATION
C     VALUE TO DAMP DOWN THE VALUES OF THE THREE SECOND PARTIAL
C     DERIVATIVES WHEN THE AVAILABLE DATA DOES NOT PRODUCE A
C     REASONABLY WELL-CONDITIONED LEAST SQUARES PROBLEM.
C     A LARGER VALUE OF STAB DAMPS THE SECOND PARTIALS MORE.
C
C     THE OUTPUT IS DZ1 AND DZ2.
C     ------------------------------------------------------------------
      DOUBLE PRECISION AVE,WT,WT1,WT2
      save             AVE,WT,WT1,WT2
      INTEGER limit
      INTEGER LOC001, LOC003
C
      INTEGER i, IROW, J, JJ, KP, LEN, MFIT, MW
      INTEGER N20008
      DOUBLE PRECISION AVESQ, C1, C2, COSINE, DMIN, DZ1, DZ2
C++ Code for .C. is INACTIVE
C      DOUBLE PRECISION SINE, STAB, SUM, W(21,6)
C++ Code for ~.C. is ACTIVE
      DOUBLE PRECISION SINE, STAB, SUM, W(MW,6)
C++ END
      LOGICAL   NEWPT
C     ------------------------------------------------------------------
      STAB = 1.D0
      IF (.NOT.(NEWPT)) GO TO 20003
       SUM = 0.0d0
       KP = MFIT
      DO 20004 I=1,KP
          SUM=SUM + W(I,1) + W(I,3)
20004 CONTINUE
C                          KP  will be converted to floating point.
       AVESQ = SUM/KP
       AVE   = SQRT(AVESQ)
C                          WE HAVE EXPERIMENTED WITH A VARIETY OF
C                          WEIGHTING FORMULAS AND FOUND NONE SIG-
C                          NIFICANTLY BETTER THAN THE UNIFORM
C                          WEIGHTING..  WT = 1.
         WT = 1.0d0
         WT1 = WT/AVE
         WT2 = WT/AVESQ
      IROW=1
      N20008=KP
      GO TO 20009
20007 IROW=IROW+1
20009 IF ((N20008-IROW).LT.0) GO TO 20008
      ASSIGN 20010 TO LOC001
      GO TO 30001
20010 GO TO 20007
C
C
20008 GO TO 30002
C
20011 GO TO 20002
C                                   SCALE NEW ROW
20003    IROW = MFIT
      ASSIGN 20012 TO LOC001
      GO TO 30001
C
20012 ASSIGN 20013 TO LOC003
      GO TO 30003
C
C
20013 CONTINUE
20002 IF(KP .GE. 5)THEN
      DMIN = min( min( min(abs(W(1,1)), abs(W(2,2))),
     &                 min(abs(W(3,3)), abs(W(4,4)))),
     &            abs(W(5,5)) )
      ELSE
          DMIN = 0.0d0
      END IF
C               TEST CONDITION OF SYSTEM.
      IF (.NOT.(DMIN .LT. 0.01d0)) GO TO 20017
C
C               SYSTEM IS ILL-CONDITIONED.
C
      IF(MFIT .LT. LIMIT)THEN
C
C               RETURN TO ADD ONE MORE POINT TO THE FITTING PROBLEM
C               IN HOPES THAT THIS WILL IMPROVE THE CONDITION.
C
          NEWPT = .FALSE.
          RETURN
      END IF
C
C               CANNOT ADD ANY MORE POINTS.  STABILIZE BY DAMPING
C               THE SECOND PARTIAL DERIVATIVES.
C
      JJ=1
      GO TO 20022
20020 JJ=JJ+1
20022 IF ((3-JJ).LT.0) GO TO 20021
            IROW = MFIT + JJ
      DO 20023 J = 1,6
               W(IROW,J) = 0.0d0
20023 CONTINUE
            W(IROW,JJ) = STAB
      ASSIGN 20026 TO LOC003
      GO TO 30003
20026 GO TO 20020
20021 CONTINUE
C
C                 SOLVE FOR DZ1 AND DZ2
C
20017 C2 = W(5,6)/W(5,5)
      C1 = (W(4,6) - C2*W(4,5)) / W(4,4)
      DZ1 = C1 / AVE
      DZ2 = C2 / AVE
      NEWPT = .TRUE.
C
      RETURN
C
C
C     PROCEDURE( SCALE ROW IROW )
C
30001  W(IROW,1) = W(IROW,1) * WT2
       W(IROW,2) = W(IROW,2) * WT2
       W(IROW,3) = W(IROW,3) * WT2
       W(IROW,4) = W(IROW,4) * WT1
       W(IROW,5) = W(IROW,5) * WT1
       W(IROW,6) = W(IROW,6) * WT
      GO TO LOC001,(20010,20012)
C
C      *******************************************
C
C     PROCEDURE (  GIVENS TRIANGULARIZATION  )
C
C     *******************************************
30002 DO 20027 J=1,6
      IF(J+1 .LE. KP)THEN
            LEN = 6-J
      DO 20032 I=J+1,KP
      IF(W(I,J) .NE. 0.0d0)THEN
                  CALL KDROTG (W(J,J),W(I,J),COSINE,SINE)
                  W(I,J) = 0.0d0
                  IF (LEN .GT. 0)
     &             CALL KDROT(LEN,W(J,J+1),MW,W(I,J+1),MW,COSINE,SINE)
      END IF
20032 CONTINUE
      END IF
20027 CONTINUE
      GO TO 20011
C
C     *************************************************
C
C     PROCEDURE (ADD ROW IROW TO TRIANGULARIZATION)
C
C      ***************************************************
C             THE NEW ROW WILL BE ACCUMULATED INTO THE TRIANGLE
C             USING GIVENS ROTATIONS.
30003 DO 20037 J=1,6
        LEN = 6 - J
      IF(W(IROW,J) .NE. 0.0d0)THEN
          CALL KDROTG (W(J,J), W(IROW,J), COSINE,SINE)
          W(IROW,J) = 0.0d0
          IF (LEN .GT. 0)
     &    CALL KDROT (LEN,W(J,J+1),MW,W(IROW,J+1),MW,COSINE,SINE)
      END IF
20037 CONTINUE
C
      GO TO LOC003,(20013,20026)
      END SUBROUTINE DTGLS

C***********************************************************************

C  File: DTGSET.[F|FOR] CONTAINS DTGSET, DTGGET, DTGPUT, DTGSIZ.
C
C     THE FOUR SUBROUTINES DTGGET, DTGPUT, DTGSET, AND DTGSIZ MUST
C     BE MUTUALLY COMPATIBLE WITH REGARD TO THE METHOD OF PACKING
C     POINTERS IN TRIANG().
C
C     IN THIS VERSION ONLY ONE POINTER IS STORED PER WORD FOR MAXIMUM
C     PORTABILITY.
C
      SUBROUTINE DTGSET( J,T1,T2,T3,P1,P2,P3, TRIANG,MT)
C>> 1996-02-02 CLL
C>> 1995-01-11 CLL Editing for inclusion into MATH77.
C     INPUT..  J AND THE INTEGERS T1,T2,T3,P1,P2,P3 DEFINING THE
C              JTH TRIANGLE.
C     THIS SUBR PUTS THIS INFO INTO THE ARRAY TRIANG().
C     C.L.LAWSON, JPL, 1979 JULY 22
C     ------------------------------------------------------------------
C--D replaces "?": ?TGSET, ?TGGET, ?TGPUT, ?TGSIZ
C     ------------------------------------------------------------------
      INTEGER  J,T1,T2,T3,P1,P2,P3, MT,TRIANG(MT)
C     ------------------------------------------------------------------
      TRIANG(6*J-5) = T1
      TRIANG(6*J-4) = T2
      TRIANG(6*J-3) = T3
      TRIANG(6*J-2) = P1
      TRIANG(6*J-1) = P2
      TRIANG(6*J) = P3
      RETURN
      END SUBROUTINE DTGSET

C***********************************************************************

      SUBROUTINE DTGGET( J, TJ, TRIANG )
C     INPUT: J, TRIANG().
C     THIS SUBR FETCHES THE JTH TRIANGLE FROM TRIANG() AND STORES IT
C     IN TJ(1:6).
C     C.L.LAWSON, JPL, 1979 JULY 22
C     ------------------------------------------------------------------
      INTEGER  J, TJ(6), TRIANG(1:*)
C     ------------------------------------------------------------------
      TJ(1) = TRIANG(6*J-5)
      TJ(2) = TRIANG(6*J-4)
      TJ(3) = TRIANG(6*J-3)
      TJ(4) = TRIANG(6*J-2)
      TJ(5) = TRIANG(6*J-1)
      TJ(6) = TRIANG(6*J)
      RETURN
      END SUBROUTINE DTGGET

C***********************************************************************

      SUBROUTINE DTGPUT( J, TJ, TRIANG, MT )
C     INPUT MT=DIMENSION OF TRIANG()
C     INPUT J AND TJ(1:6).    THIS SUBR PUTS THE CONTENTS OF
C     TJ(1:6) INTO THE ARRAY TRIANG() TO DEFINE THE JTH TRIANGLE.
C     C.L.LAWSON, JPL, 1979 JULY 22
C     ------------------------------------------------------------------
      INTEGER  J, TJ(1:6), MT, TRIANG(1:MT)
C     ------------------------------------------------------------------
      TRIANG(6*J-5) = TJ(1)
      TRIANG(6*J-4) = TJ(2)
      TRIANG(6*J-3) = TJ(3)
      TRIANG(6*J-2) = TJ(4)
      TRIANG(6*J-1) = TJ(5)
      TRIANG(6*J) = TJ(6)
      RETURN
      END SUBROUTINE DTGPUT

C***********************************************************************

      SUBROUTINE DTGSIZ( MT, NTLIM )
C     INPUT..     MT = DIMENSION OF ARRAY TRIANG()
C     OUTPUT..    NTLIM = MAX NUMBER OF TRIANGLES WHOSE DESCRIPTION
C                 POINTERS CAN BE STORED IN TRIANG().  THERE ARE SIX
C                 POINTERS PER TRIANGLE.
C     C.L.LAWSON, JPL, 1979 JULY 22
C     ------------------------------------------------------------------
      INTEGER  MT, NTLIM
C     ------------------------------------------------------------------
      NTLIM = MT/6
      RETURN
      END SUBROUTINE DTGSIZ

C***********************************************************************

      SUBROUTINE KDROT(N,X,INCX,Y,INCY,C,S)
C>> 1994-11-11 DROT  Krogh   Declared all vars.
C>> 1994-10-20 DROT   Krogh  Changes to use M77CON
C>> 1985-08-02 DROT   Lawson  Initial code.
C--D replaces "?": ?ROT
C
C     MULTIPLY THE 2 X 2 MATRIX  ( C S) TIMES THE 2 X N MATRIX (X**T)
C                                (-S C)                        (Y**T)
C     WHERE **T INDICATES TRANSPOSE.    THE ELEMENTS OF X ARE IN
C     X(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
C     LX = (-INCX)*N, AND SIMILARLY FOR Y USING LY AND INCY.
      INTEGER N, INCX, INCY, NSTEPS, I, KX, KY
      DOUBLE PRECISION X,Y,C,S,ZERO,ONE,W,Z
      DIMENSION X(*),Y(*)
C
      DATA ZERO,ONE/0.D0,1.D0/
      IF(N .LE. 0 .OR. (S .EQ. ZERO .AND. C .EQ. ONE)) GO TO 40
      IF(.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
           NSTEPS=INCX*N
           DO 10 I=1,NSTEPS,INCX
                W=X(I)
                Z=Y(I)
                X(I)=C*W+S*Z
                Y(I)=-S*W+C*Z
   10           CONTINUE
           GO TO 40
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF(INCX .LT. 0) KX=1-(N-1)*INCX
           IF(INCY .LT. 0) KY=1-(N-1)*INCY
C
           DO 30 I=1,N
                W=X(KX)
                Z=Y(KY)
                X(KX)=C*W+S*Z
                Y(KY)=-S*W+C*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END SUBROUTINE KDROT

C***********************************************************************

      SUBROUTINE KDROTG(A,B,C,S)
C>> 1994-10-20 DROTG  Krogh  Changes to use M77CON
C>> 1994-04-19 DROTG  Krogh   Changed to use generic intrinsics.
C>> 1985-08-02 DROTG  Lawson  Initial code.
C--D replaces "?": ?ROTG, ?ROT
C
C     DESIGNED BY C.L.LAWSON, JPL, 1977 SEPT 08
C
C
C     CONSTRUCT THE GIVENS TRANSFORMATION
C
C         ( C  S )
C     G = (      ) ,    C**2 + S**2 = 1 ,
C         (-S  C )
C
C     WHICH ZEROS THE SECOND ENTRY OF THE 2-VECTOR  (A,B)**T .
C
C     THE QUANTITY R = (+/-)SQRT(A**2 + B**2) OVERWRITES A IN
C     STORAGE.  THE VALUE OF B IS OVERWRITTEN BY A VALUE Z WHICH
C     ALLOWS C AND S TO BE RECOVERED BY THE FOLLOWING ALGORITHM:
C           IF Z=1  SET  C=0.D0  AND  S=1.D0
C           IF ABS(Z) .LT. 1  SET  C=SQRT(1-Z**2)  AND  S=Z
C           IF ABS(Z) .GT. 1  SET  C=1/Z  AND  S=SQRT(1-C**2)
C
C     NORMALLY, THE SUBPROGRAM KDROT(N,X,INCX,Y,INCY,C,S) WILL
C     NEXT BE CALLED TO APPLY THE TRANSFORMATION TO A 2 BY N MATRIX.
C
C ------------------------------------------------------------------
C
      DOUBLE PRECISION  A, B, C, S, U, V, R
      IF (ABS(A) .LE. ABS(B)) GO TO 10
C
C *** HERE ABS(A) .GT. ABS(B) ***
C
      U = A + A
      V = B / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF A
C
      R = SQRT(.25D0 + V**2) * U
C
C     NOTE THAT C IS POSITIVE
C
      C = A / R
      S = V * (C + C)
      B = S
      A = R
      RETURN
C
C *** HERE ABS(A) .LE. ABS(B) ***
C
   10 IF (B .EQ. 0.D0) GO TO 20
      U = B + B
      V = A / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF B
C     (R IS IMMEDIATELY STORED IN A)
C
      A = SQRT(.25D0 + V**2) * U
C
C     NOTE THAT S IS POSITIVE
C
      S = B / A
      C = V * (S + S)
      IF (C .EQ. 0.D0) GO TO 15
      B = 1.D0 / C
      RETURN
   15 B = 1.D0
      RETURN
C
C *** HERE A = B = 0.D0 ***
C
   20 C = 1.D0
      S = 0.D0
      RETURN
      END SUBROUTINE KDROTG

C***********************************************************************

      SUBROUTINE DSORTP (A, M, N, P)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  All rights reserved.  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1995-11-15 DSORTP  Krogh  SFTRAN => Fortran, removed mult. entry.
C>> 1994-10-19 DSORTP  Krogh  Changes to use M77CON
C>> 1992-11-23 DSORTP  Snyder  Add entry DSORTQ.
C>> 1991-04-02 DSORTP  Snyder  Repair no permutation vector IF m-n < 10
C>> 1988-11-22 DSORTP  Snyder  Initial code.
C--D replaces "?": ?SORTP, ?SORTQ
C
C     Sort the M:N-vector A.
C     A is not disturbed.  P is set so that A(P(J)) is the J'th element
C     of the sorted sequence.
C     Enter at DSORTQ to use pre-specified permutation vector.
C
C     To sort an array A' into descending order, let A = -A'
C     To sort an array A' into ascending order according to the
C     absolute value of the elements let A = ABS(A').
C     To sort an array A' into decending order according to the
C     absolute value of the elements let A = -ABS(A').
C
      INTEGER M, N, P(*)
C--D Next line special: I
      DOUBLE PRECISION A(*)
      INTEGER CL
C                      Get permutation vector for sorting
      DO 20 cl = m, n
         p(cl)=cl
   20 CONTINUE
      CALL DSORTQ (A, M, N, P)
      RETURN
      END SUBROUTINE DSORTP


C***********************************************************************

      SUBROUTINE DSORTQ (A, M, N, P)
      INTEGER M, N
C--D Next line special: I
      DOUBLE PRECISION A(*)
      INTEGER P(*)
C
C     *****     Local Variables     ************************************
C
C BL      is the left bound of the sub-array to be sorted at the next
C         step.
C BR      is the right bound of the sub-array to be sorted at the next
C         step.
C CL      is the current left bound of the unsorted sub-array.
C CR      is the current right bound of the unsorted sub-array.
C PARTN   is the partition element.
C PTEMP   holds elements of P during exchanges.
C STACKL  keeps track of the left bounds of sub-arrays waiting to be
C         sorted.
C STACKR  keeps track of the right bounds of sub-arrays waiting to be
C         sorted.
C STKTOP  keeps track of the top of the stacks.
C
      INTEGER BL,BR,CL,CR
C--D Next line special: I
      DOUBLE PRECISION PARTN
      INTEGER PTEMP,STACKL(32),STACKR(32),STKTOP
C
C     *****     Executable Statements     ******************************
C
      IF (n-m.GE.10) THEN
         stktop=1
         stackl(1)=m
         stackr(1)=n
   40    CONTINUE
            bl=stackl(stktop)
            br=stackr(stktop)
            stktop=stktop-1
C           Choose a partitioning element.  Use the median of the first,
C           middle and last elements.  Sort them so the extreme elements
C           can serve as sentinels during partitioning.
            cl=(bl+br)/2
            ptemp=p(cl)
            IF (a(p(bl)).GT.a(ptemp)) THEN
               p(cl)=p(bl)
               p(bl)=ptemp
               ptemp=p(cl)
            END IF
            IF (a(p(bl)).GT.a(p(br))) THEN
               cr=p(bl)
               p(bl)=p(br)
               p(br)=cr
            END IF
            IF (a(ptemp).GT.a(p(br))) THEN
               p(cl)=p(br)
               p(br)=ptemp
               ptemp=p(cl)
            END IF
            p(cl)=p(br-1)
            p(br-1)=ptemp
            partn=a(ptemp)
C           Partition the sub-array around PARTN.  Exclude the above
C           considered elements from partitioning because they're al-
C           ready in the correct subfiles.  STOP scanning on equality to
C           prevent files containing equal values from causing a loop.
            cl=bl
            cr=br-1
   80       CONTINUE
  100          cl=cl+1
               IF (a(p(cl)) .LT. partn) GO TO 100
  120          cr=cr-1
               IF (a(p(cr)) .GT. partn) GO TO 120
               IF (cl.GT.cr) GO TO 150
               ptemp=p(cl)
               p(cl)=p(cr)
               p(cr)=ptemp
               GO TO 80
  150       CONTINUE
C           Put sub-arrays on the stack IF they're big enough.  Put the
C           larger under the smaller, so the smaller will be done next.
C           This makes the upper bound of the stack depth log2 (n-m+1).
C           (The "Hibbard" modification of quicksort).
            IF (cl-bl .GT. br-cr) THEN
               IF (cl-bl.GT.10) THEN
                  stktop=stktop+1
                  stackl(stktop)=bl
                  stackr(stktop)=cr
               END IF
               IF (br-cr.GT.10) THEN
                  stktop=stktop+1
                  stackl(stktop)=cl
                  stackr(stktop)=br
               END IF
            ELSE
               IF (br-cr.GT.10) THEN
                  stktop=stktop+1
                  stackl(stktop)=cl
                  stackr(stktop)=br
               END IF
               IF (cl-bl.GT.10) THEN
                  stktop=stktop+1
                  stackl(stktop)=bl
                  stackr(stktop)=cr
               END IF
            END IF
         IF (stktop .NE.0) GO TO 40
      END IF
C     Clean up small subfiles using insertion sort on everything.
      DO 200 cr = m+1, n
         ptemp=p(cr)
         partn=a(ptemp)
         cl=cr
  180    IF (a(p(cl-1)) .GT. partn) THEN
            p(cl)=p(cl-1)
            cl=cl-1
            IF (cl .GT. m) GO TO 180
         END IF
         p(cl)=ptemp
  200   CONTINUE
        RETURN
        END SUBROUTINE DSORTQ

        END MODULE didesub