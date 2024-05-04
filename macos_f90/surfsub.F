C***********************************************************************
C	Begin file surfsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2010, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

        MODULE surfsub

        PUBLIC
        PRIVATE :: Calc_ToricNormal, SAsphere, AZPSolve, AZPB, SMonomial,
     &             MZPSolve, MZPB, UDSZPSolve, UDSZPB, GSZPSolve, SGSrf
     &             NGSrf
        CONTAINS

	SUBROUTINE FlatSrf(ifBuild,ifLNsrf,LROK,
     &                     psi,pv,prot,pin,ihat,
     &		           pr,L,Nvec,Nmag,Nhat,dNdp)
	use math_mod

	IMPLICIT NONE

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &         pr(3),L,Nhat(3)
	REAL*8 prel(3),b,c,Nvec(3),Nmag,dNdp(3,3),S1
	INTEGER i
	LOGICAL LROK,ifBuild,ifLNsrf

C	Compute incidence point of ray on the surface

	Do i=1,3
	  prel(i)=pin(i)-pv(i)
	  Nvec(i)=psi(i)
        End Do

	b=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	IF (b.EQ.0d0) GO TO 98
	c=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)
	L=-c/b  ! distance from previous ray point to intersection point on flat surf, -jzlou

	IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98

	Do i=1,3
	  pr(i)=pin(i)+L*ihat(i)
    	End Do

	Nmag=1d0
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  Do i=1,3
	    Nhat(i)=Nvec(i)
    	  End Do
	ELSE
	  DO i=1,3
	    Nhat(i)=-Nvec(i)
    	  End Do
	END IF

C	Compute partials of the normal

	IF (ifBuild) THEN
	  CALL DZERO(dNdp,9)
	END IF

	RETURN

C	Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.
	RETURN
	END SUBROUTINE FlatSrf

C***********************************************************************

! Input arguments: pin, ihat, ifLNsrf
! Output argments: pr, Nvec, Nmag, Nhat

	SUBROUTINE ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,
     &		          pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	use macos_debug
	use math_mod

	IMPLICIT NONE

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),pr(3),L,Nhat(3),
     &	       prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &         Nvec(3),Nmag,dNdp(3,3),Kc,Kr,k2,itpsi,psitprel,
     &	       D1(9),D2(9),S1,S2
	INTEGER i
	LOGICAL LROK,ifBuild,ifLNsrf

! ifLNsrf is a flag indicating whether to allow ray length 'L' to be negative -jzlou
! pr is ray incident point on current surface (to be computed) -jzlou

! ihat = RayDir(1:3,iRay) = direction of incident ray
! Nhat = surface normal at ray incident point (jzlou guessed);
!        Nhat will be used to compute direction of reflected ray (rhat) in calling routine -jzlou


C	Compute incidence point of ray on the surface

	Do i=1,3
	  prel(i)=pin(i)-pv(i) ! current Ray position (on previous surface) minus current VptElt
    	End Do
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3) ! = ihat projected onto psi
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3) ! = psi projected onto prel

	IF (ABS(Kr).LT.1d22) THEN
	  mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
	  a=1d0+Kc*itpsi*itpsi
	  b=2d0*(Kr*itpsi+Kc*itpsi*psitprel
     &	    +ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
	  c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
	ELSE
	  a = 0d0
	  b = itpsi
	  c = psitprel
	END IF

C	Compute intersection point and normal vector

	IF ((ABS(a).LE.1d-20).OR.(ABS(Kr).GE.1d20)) THEN
	  L=-c/b
	  IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98
	ELSE
	  k2 = b*b - 4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  IF (b.GE.0d0) THEN
	    ksqrt = -0.5d0*(b+SQRT(k2))
	  ELSE
	    ksqrt = -0.5d0*(b-SQRT(k2))
	  END IF
	  kp = ksqrt/a
	  km = c/ksqrt
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



	IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98

	S1 = Kr + Kc*(psitprel+L*itpsi)
        DO i=1,3
	  S2=L*ihat(i)
	  pr(i)=pin(i)+S2  ! pr is ray incident point on current surface -jzlou
	  Nvec(i)=S1*psi(i)+prel(i)+S2 ! direction of surface normal at ray incident point -jzlou
    	END DO

	Nmag = SQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	S1 = Nvec(1)*ihat(1) + Nvec(2)*ihat(2) + Nvec(3)*ihat(3)
	IF (S1 .LE. 0d0) THEN
	  DO i=1,3
	    Nhat(i)=Nvec(i)/Nmag
          END DO
	ELSE
	  DO i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
          END DO
	END IF

C	Compute partials of the normal
C  Compute partials of the local normal (Nhat) wrt psi -jzlou

	IF (ifBuild) THEN
	  IF (ABS(Kr).LT.1d10) THEN
	    CALL DSMPROD(D1,psi,Kc,3) ! D1 = psi*Kc
	    CALL DOUTER(M,D1,psi)
	    M(1,1)=M(1,1)+1d0
	    M(2,2)=M(2,2)+1d0
	    M(3,3)=M(3,3)+1d0
	  ELSE
	    CALL DZERO(M,9)
	  END IF
	  S1=1d0/Nmag
	  CALL DPROJECT(D1,Nhat)
	  CALL DMPROD(D2,D1,M,3,3,3)
	  CALL DSMPROD(dNdp,D2,S1,9) ! dNdp is derivative of Nhat wrt psi, a 3x3 matrix -jzlou
                                     ! Nhat is surface normal at ray incident point
          ! Since Nhat is computed from an analytical function/equation of psi, its derivative
          ! wrt psi can be analytically computed as well -jzlou
	END IF

	RETURN

        ! Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.

	RETURN
	END SUBROUTINE ConSrf

C***********************************************************************

	SUBROUTINE NSConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,
     &	  pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	use math_mod

	IMPLICIT NONE

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &         pr(3),L,Nhat(3),prel(3),M(3,3),a,b,c,
     &         Nvec(3),Nmag,dNdp(3,3),Kc,Kr,k2,itpsi,psitprel,
     &         ksqrt,kp,km,mpr,D1(9),D2(9),S1,S2
	INTEGER i
	LOGICAL LROK,ifBuild,ifLNsrf


C	Compute incidence point of ray on the surface

	DO i=1,3
	  prel(i)=pin(i)-pv(i)
        END DO
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d22) THEN
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

        ! Compute intersection point and normal vector

	IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d22)) THEN
	  if (.false.) then !(abs(b) .le. 1d-10) then
	    print*,'NSConSrf(): b = ', b
	    stop
	  end if
	  L=-c/b
	ELSE
	  k2=b*b-4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  ksqrt=DSQRT(k2)
	  S1=ABS(ABS(b)-ksqrt)
	  IF (S1.LT.1d-4) THEN
	    kp=(ksqrt-b)/(2d0*a)
	    km=(-ksqrt-b)/(2d0*a)
	  ELSE
	    kp=(-2d0*c)/(b+ksqrt)
	    km=(-2d0*c)/(b-ksqrt)
	  END IF
c	  IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	  IF (kp.LT.km) THEN
	    IF (kp.GT.1d-12) THEN
	      L=kp
	    ELSE
	      L=km
	    END IF
	  ELSE
	    IF (km.GT.1d-12) THEN
	      L=km
	    ELSE
	      L=kp
	    END IF
	  END IF
	END IF

cLnegOK -- Nonsequential surfaces must have positive L!!
	IF (L.LT.1d-12) GO TO 98

	S1=Kr+Kc*(psitprel+L*itpsi)
	DO 2 i=1,3
	  S2=L*ihat(i)
	  pr(i)=pin(i)+S2
	  Nvec(i)=S1*psi(i)+prel(i)+S2
   2	CONTINUE
	Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	IF (Nmag.LE.0d0) GO TO 98
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  DO 3 i=1,3
	    Nhat(i)=Nvec(i)/Nmag
   3	  CONTINUE
	ELSE
	  DO 4 i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
   4	  CONTINUE
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
	  S1=1d0/Nmag
	  CALL DPROJECT(D1,Nhat)
	  CALL DMPROD(D2,D1,M,3,3,3)
	  CALL DSMPROD(dNdp,D2,S1,9)
	END IF

	RETURN

C	Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.
	RETURN
	END
C	End of NSConSrf

C***********************************************************************

	SUBROUTINE AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,prot,
     &			  pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	use math_mod

	IMPLICIT NONE

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &  pr(3),L,xhat(3),yhat(3),Nhat(3),AnaCoef(4)

	REAL*8 prel(3),a,c,ac,
     &  Nvec(3),Nmag,dNdp(3,3),Kx,Ky,Cx,Cy,itpsi,psitprel,rho(3)

	REAL*8 D1(9),D2(9),S1,S2

	REAL*8 ix,iy,psix,psiy,px,py,ix2,iy2,px2,py2,pxxi,pyyi,L2,
     &	  opkrx,opkry,f0,f1,p0,p1,p2,c0,c1,c2,a0,a1,a2,q0,q1,q2,q3,q4,
     &	  z,z2,rhox,rhoy,rhotpsi,dSxo2,dSyo2,dSpsio2,rhomag,rhomag2

	COMPLEX*16 CL(4)

	INTEGER i,j,Order

	LOGICAL LROK,ifBuild,ifLNsrf

C  Compute incidence point of ray on the surface
C    Set surface coefficients:
	Cx=-1d0/AnaCoef(1)
	Cy=-1d0/AnaCoef(2)
	Kx=AnaCoef(3)
	Ky=AnaCoef(4)
	CALL DORTHOGANALIZE(psi,xhat,yhat)
	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
   1	CONTINUE
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	ix=ihat(1)*xhat(1)+ihat(2)*xhat(2)+ihat(3)*xhat(3)
	iy=ihat(1)*yhat(1)+ihat(2)*yhat(2)+ihat(3)*yhat(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)
c	psix=psi(1)*xhat(1)+psi(2)*xhat(2)+psi(3)*xhat(3)
c	psiy=psi(1)*yhat(1)+psi(2)*yhat(2)+psi(3)*yhat(3)
	px=prel(1)*xhat(1)+prel(2)*xhat(2)+prel(3)*xhat(3)
	py=prel(1)*yhat(1)+prel(2)*yhat(2)+prel(3)*yhat(3)
	ix2=ix*ix
	iy2=iy*iy
	px2=px*px
	py2=py*py
	pxxi=2d0*px*ix
	pyyi=2d0*py*iy
	opkrx=(1d0+Kx)*Cx*Cx
	opkry=(1d0+Ky)*Cy*Cy

	f0=psitprel
	f1=itpsi
	p0=f0*f0
	p1=2d0*f0*f1
	p2=f1*f1
	c0=Cx*px2+Cy*py2
	c1=Cx*pxxi+Cy*pyyi
	c2=Cx*ix2+Cy*iy2
	a0=opkrx*px2+opkry*py2
	a1=opkrx*pxxi+opkry*pyyi
	a2=opkrx*ix2+opkry*iy2
	q0=a0*p0-2d0*c0*f0+c0*c0
	q1=a0*p1+a1*p0-2d0*(c1*f0+c0*f1-c0*c1)
	q2=a0*p2+a1*p1+a2*p0-2d0*(c2*f0+c1*f1-c0*c2)+c1*c1
	q3=a1*p2+a2*p1-2d0*(c2*f1-c1*c2)
	q4=a2*p2+c2*c2

C    Find surface roots and select L:
        CALL ROOT4(Order,CL(1),CL(2),CL(3),CL(4),q4,q3,q2,q1,q0)
        IF (Order.EQ.0) THEN
          L = -f0/f1
        ELSE
          L      = 1d22
          rhomag = 1d22
          DO i=1,Order
            IF (isReal(CL(i))) THEN
              L2  = DBLE(CL(i))
              pr  = pin + L2*ihat
              rho = pr - pv

              rhomag2=rho(1)*rho(1)+rho(2)*rho(2)+rho(3)*rho(3)
              IF ((L2.GE.0d0).AND.(rhomag2.LT.rhomag)) THEN
                rhomag = rhomag2
                L      = L2
              END IF
            END IF
          END DO
        END IF

cLnegOK -- Anamorphic surfaces temporarily require positive L
	IF (L.EQ.1d22) GO TO 98

C  Compute normal
C    Set coefficients
	DO i=1,3
	  pr(i)=pin(i)+L*ihat(i)
	  rho(i)=pr(i)-pv(i)
 	END DO
	IF (Order.NE.0d0) THEN
	  L2=L*L
c	  z=f0+f1*L
	  z2=p0+p1*L+p2*L2
	  c=c0+c1*L+c2*L2
	  ac=a0+a1*L+a2*L2
	  a=ac/c

cdebug -- check this formula:
	  S1=1d0-2d0*a*z2/c

	  rhox=rho(1)*xhat(1)+rho(2)*xhat(2)+rho(3)*xhat(3)
	  rhoy=rho(1)*yhat(1)+rho(2)*yhat(2)+rho(3)*yhat(3)
	  rhotpsi=rho(1)*psi(1)+rho(2)*psi(2)+rho(3)*psi(3)
	  dSxo2=(z2*opkrx/c+S1*Cx)*rhox
	  dSyo2=(z2*opkry/c+S1*Cy)*rhoy
	  dSpsio2=a*rhotpsi-1d0
	  DO i=1,3
	    Nvec(i)=dSxo2*xhat(i)+dSyo2*yhat(i)+dSpsio2*psi(i)
	  END DO
	ELSE
	  CALL DEQUATE(Nvec,psi,3)
	END IF

	Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  DO 3 i=1,3
	    Nhat(i)=Nvec(i)/Nmag
   3	  CONTINUE
	ELSE
	  DO 4 i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
   4	  CONTINUE
	END IF

C	Compute partials of the normal

	IF (ifBuild) THEN
	END IF

	RETURN

C	Here if ray misses surface

  98	CONTINUE

!	print*,'** AnaSrf: ray lost !!'
	LROK=.FALSE.
	RETURN
	END SUBROUTINE AnaSrf

C***********************************************************************

        SUBROUTINE ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,prot,
     &  		    pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	! argument 'prot' not used
	use math_mod

        IMPLICIT NONE
c
        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &         pr(3),L,xhat(3),yhat(3),Nhat(3), AnaCoef(4)

c
        REAL*8 prel(3),Acoef,Bcoef,Ccoef,Nvec(3),Nvl(3),Nmag,
     &         dNdp(3,3),Krx,Kry,Cx,psitprel,rho(3),
     &         alpha, beta, gamma, nvroot
c
        REAL*8 ix,iy,itpsi,psix,psiy,px,py,ix2,iy2,px2,py2,
     &         pxxi,pyyi,L2,q0,q1,q2,q3,q4,xpt,ypt,zpt
c
        REAL*8  D1(9),D2(9),S1,S2
c
        COMPLEX*16 CL(4)
        INTEGER i,j,Order
        LOGICAL LROK,ifBuild,ifLNsrf
	Logical :: isCYL, useReddingToricNormal
	Logical :: prtL=.false.  ! debugging only, -jzlou

	if (prtL) print*,'** input psi =',psi
	if (prtL) print*,'** input ihat =',ihat

        ! flag for whether use Redding's or jzlou's toric normal
        ! calculation
	useReddingToricNormal=.true.

C
C  Compute incidence point of ray on the surface
C  Set surface coefficients:
C
	! Krx and Kry are radius in X and Y.
	! AnaCoef(4) is assumed zero for toric surface -jzlou
        ! Cx is conic constant in X direction, Cy is assumed zero (circle in
        ! Y direction)
        Krx=-AnaCoef(1); Kry=-AnaCoef(2); Cx=AnaCoef(3)

	! yhat is computed from psi cross xhat
	! xhat, yhat and psi orthogonal to each other
        CALL DORTHOGANALIZE(psi,xhat,yhat)

        prel(1:3)=pin(1:3)-pv(1:3)

	! Convert ray direction cosine from global
        ! to element-local coord system -jzlou
        itpsi = ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
        ix = ihat(1)*xhat(1)+ihat(2)*xhat(2)+ihat(3)*xhat(3)
        iy = ihat(1)*yhat(1)+ihat(2)*yhat(2)+ihat(3)*yhat(3)

	! Convert ray starting point to element-local coord system -jzlou
        psitprel = psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)
        px   = prel(1)*xhat(1)+prel(2)*xhat(2)+prel(3)*xhat(3)
        py   = prel(1)*yhat(1)+prel(2)*yhat(2)+prel(3)*yhat(3)

	! psix and psiy never used later in this routine.
	! Since psi, xhat, yhat are orthogonal to one another now,
	! psix and psiy must be zero. -jzlou
        psix = psi(1)*xhat(1)+psi(2)*xhat(2)+psi(3)*xhat(3)
        psiy = psi(1)*yhat(1)+psi(2)*yhat(2)+psi(3)*yhat(3)


	!
	! This 'if' section was added by J.Lou. 11/11/04
	!
	if (abs(Kry) .gt. 1d07) then
	  isCYL=.true.
	else
	  isCYL=.false.
	end if
	!
	if (isCYL) then

	   ! Compute coeffs for cylindrical equation
	   q4=0.0; q3=0.0
	   q2 = ix*ix+(1.0+Cx)*itpsi*itpsi
	   q1 = 2.0*((1.0+Cx)*psitprel*itpsi-Krx*itpsi+px*ix)
           q0 = (1.0+Cx)*psitprel*psitprel-2.0*Krx*psitprel+px*px

	   if (prtL) print*,'B4 ROOT4 (cyl): q2,q1,q0=',
     &                       q2,q1,q0

	   ! Solve for path length
	   ! for cylinder surface, there are at most two real roots
	   CALL ROOT4(Order,CL(1),CL(2),CL(3),CL(4),q4,q3,q2,q1,q0)

	   if (prtL) print*,'In ToricSrf(cyl): CL(1:4)=', CL(1:4)

	   if (isReal(CL(1)) .and. dble(CL(1)).ge.0 .and.
     &         isReal(CL(2)) .and. dble(CL(2)).ge.0) then
	     if (AnaCoef(1).ge.0) then
	       ! Take 'first' intesection point, which is the one
	       ! closer to the ray starting point, i.e. the smaller
	       ! positive 'L' solution.
	       if (dble(CL(1))<dble(CL(2))) then
	         L=dble(CL(1))
	       else
	         L=dble(CL(2))
	       end if
	     else ! AnaCoef(1).lt.0
	       ! Take 'second' intesection point, which is the one
	       ! farther to the ray starting point, i.e. the larger
	       ! positive 'L' solution.
	       if (dble(CL(1))>dble(CL(2))) then
                 L=dble(CL(1))
               else
                 L=(CL(2))
               end if
	     end if
	   else if ((isReal(CL(1)) .and. dble(CL(1)).ge.0) .or.
     &              (isReal(CL(2)) .and. dble(CL(2)).ge.0)) then
	     if (isReal(CL(1)) .and. dble(CL(1)).ge.0) then
	       L=dble(CL(1))
	     else
	       L=dble(CL(2))
	     end if
	   else
	     ! No valid ray intersection with cylindrical surface
	   end if

	   if (prtL) print*,'***In cyl: B4 goto 95, L=', L
	   goto 95  ! skip true Toric case
	end if ! isCyL

        ix2  = ix*ix; iy2  = iy*iy
        px2  = px*px; py2  = py*py
        pxxi = px*ix; pyyi = py*iy

        If((Cx.eq.-1.0) .and. (.not.isCYL)) Cx=Cx-2.0e-6

        alpha = Kry - Krx/(1.0+Cx)
	!print*,'alpha =', alpha

        beta  = alpha*alpha + Krx*Krx/(1.0+Cx)/(1.0+Cx)
        gamma = 2.0*alpha/(1.0+Cx)
	if (prtL) print*,'gamma =', gamma

        Acoef = itpsi*itpsi+ix2/(1.0+Cx)+iy2
	if (prtL) print*,'Acoef =', Acoef
        Bcoef = 2.0 * (itpsi*(psitprel-Kry) + pxxi/(1.0+Cx) + pyyi)
	if (prtL) print*,'Bcoef =', Bcoef
        Ccoef = (Kry-psitprel)*(Kry-psitprel)
     &         + px2/(1.0+Cx) + py2 - beta
	if (prtL) print*,'Ccoef =', Ccoef

        IF(dabs(alpha).gt.1.d-9) Then
           Acoef = Acoef/gamma
           Bcoef = Bcoef/gamma
           Ccoef = Ccoef/gamma

	   if (prtL) print*,'2: Acoef,Bcoef,Ccoef=',Acoef,Bcoef,Ccoef

           q4 = Acoef*Acoef
           q3 = 2.0*Acoef*Bcoef
           q2 = Bcoef*Bcoef + 2.0*Acoef*Ccoef + (1.0+Cx)*ix2
           q1 = 2.0 * (Bcoef*Ccoef+(1.0+Cx)*pxxi)
           q0 = Ccoef*Ccoef-Krx*Krx + (1.0+Cx)*px2

           if(dabs(q4).lt.1d-16) q4=0d0
           if(dabs(q3).lt.1d-16) q3=0d0
           if(dabs(q2).lt.1d-16) q2=0d0
           if(dabs(q1).lt.1d-16) q1=0d0
           if(dabs(q0).lt.1d-16) q0=0d0

           !
           ! Find surface roots and select L:
           !

	   if (prtL) print*,'B4 ROOT4: q4,q3,q2,q1,q0=',q4,q3,q2,q1,q0

           CALL ROOT4(Order,CL(1),CL(2),CL(3),CL(4),q4,q3,q2,q1,q0)

           if (prtL) print*,'In ToricSrf(): Order=', Order
           if (prtL) print*,'In ToricSrf(): CL(1:4)=', CL(1:4)

           IF( Order.EQ.0 ) THEN
              L= -psitprel/itpsi
	      if (prtL) print*,'1: L= -psitprel/itpsi=',L
!          ELSE IF (Order.EQ.2) THEN
! need to fill in something here ... -jzlou

	   ELSE

	     if (.false.) then

	     ! dcr's original version
	     ! L is taken as minimum of roots when Cx < 1.0. -jzlou
             L= 1d22
             DO i=1,Order
               IF (isReal(CL(i))) THEN
                  L2=DBLE(CL(i))
                  IF(L2.GE.0.0) THEN
                    IF(L.EQ.1D22.or.Cx.LE.1.0)  THEN
                       L=MIN(L,L2)
                    ELSE
                       L=MAX(L,L2)
                    ENDIF
                  END IF
               END IF
             END DO

	     end if  ! dcr version

	     if (prtL) then
	        print*,'===>ToricSrf (ReddingCode): L=', L
	        print*,'Quit ...'
	        stop
	     end if

              if (.true.) then
C       Fix by S.A.B.  maybe need abs(Krx) and abs(Kyr) ???  doesn't always work
              IF (Krx.LT.Kry) THEN
                L= CL(4)
	        if (prtL) print*,'2: L= CL(4)=',L
              ELSE
                L= CL(3)
	        if (prtL) print*,'3: L= CL(3)=',L
              END IF
	      end if ! sab version

           END IF
        ELSE
           IF( dabs(Acoef).gt.1.d-14 ) Then
              CL(1) = (-Bcoef + dsqrt(Bcoef*Bcoef-4.0*Acoef*Ccoef))
     &                / (2.0*Acoef)
              CL(2) = (-Bcoef - dsqrt(Bcoef*Bcoef-4.0*Acoef*Ccoef))
     &                / (2.0*Acoef)
              DO i=1,2
                 L= 1d22
                 IF (isReal(CL(i))) THEN
                    L2=DBLE(CL(i))
                    IF (L2.GE.0d0) THEN
                       IF(L.EQ.1D22.or.Cx.LE.0.0)  THEN
                          L=MIN(L,L2)
	                  if (prtL) print*,'4: L= MIN(L,L2)=',L
                       ELSE
                          L=MAX(L,L2)
	                  if (prtL) print*,'5: L= MAX(L,L2)=',L
                       ENDIF
                    END IF
                 END IF
              END DO
           ELSE IF (dabs(Bcoef).gt.1.d-12) Then
              L= -Ccoef/Bcoef
	      if (prtL) print*,'6: L= -Ccoef/Bcoef=',L
           ELSE
              L= -psitprel/itpsi
              if (prtL) print*,'7: L= -psitprel/itpsi=',L
           END IF
        END IF

	if (prtL) print*,'**** In ToricSurf(): OPL = ', L

        IF (L.EQ.1d22) GOTO 98

	! for cylindrical case to jump
 95     continue

        !
        !  Compute normal. First set coefficients.
        !
        DO i=1,3
          pr(i)  = pin(i)+L*ihat(i)
          rho(i) = pr(i)-pv(i)  ! rho is local, never used after
        END DO

        IF( Order.NE.0d0 ) THEN
           xpt = px + L*ix
           ypt = py + L*iy
           zpt = psitprel + L*itpsi

	   if (isCYL) then
	     ! cylinder, added by jzlou
	     Nvl(1) =2.0*xpt
	     Nvl(2) =0.0
	     Nvl(3) =2.0*((1.0+Cx)*zpt-Krx)
	   else
	     ! true torus, not cylinder

	     if (useReddingToricNormal) then
               nvroot = dsqrt(Krx*Krx-(1.0+Cx)*xpt*xpt)
               Nvl(1) =-0.5*xpt/(Kry-zpt)*(1.0/(1.0+Cx)+alpha/nvroot)
               Nvl(2) =-0.5*ypt/(Kry-zpt)
               Nvl(3) = 0.5d0
	     else
	       ! use torus normal computed by Lou
	       Call Calc_ToricNormal(Kry,Krx,Cx,xpt,ypt,zpt,Nvl)
	     end if
	   end if

           Nvec(1)=xhat(1)*Nvl(1)+yhat(1)*Nvl(2)+psi(1)*Nvl(3)
           Nvec(2)=xhat(2)*Nvl(1)+yhat(2)*Nvl(2)+psi(2)*Nvl(3)
           Nvec(3)=xhat(3)*Nvl(1)+yhat(3)*Nvl(2)+psi(3)*Nvl(3)
        ELSE
           CALL DEQUATE(Nvec,psi,3)
        END IF

        if (prtL) print*,'Toric: Nvec(1:3) = ', Nvec(1:3)

        Nmag=dsqrt(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
        S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)

        IF( S1.LE.0d0 ) THEN
          Nhat(1:3) = Nvec(1:3)/Nmag
        ELSE
          Nhat(1:3) =-Nvec(1:3)/Nmag
        END IF
        !
        !       Compute partials of the normal
        !
        IF (ifBuild) THEN
        END IF

        RETURN

        ! Here if ray misses surface
  98    CONTINUE

        LROK=.FALSE.
        RETURN
        END SUBROUTINE ToricSrf


	SUBROUTINE Calc_ToricNormal(Ry,Rx,Kx,x,y,z,nv)
        ! Kx = conic constant in X
        Implicit none
        Real*8 :: Ry,Rx,Kx,x,y,z,nv(3)
        Real*8 :: x2,tt1,tt2,Zx

        tt1=(1.0+Kx)/(Rx*Rx)
        x2=x*x; tt2=sqrt(1.0-tt1*x2)
	Zx=x2/(1.0+tt2)/Rx

        nv(1)=2*x*(1.0+tt2)+tt1*x2*x/tt2
        nv(1)=nv(1)/Rx/((1.0+tt2)*(1.0+tt2))
        nv(1)=2.0*nv(1)*(Ry-Zx)
        nv(2)=2.0*y; nv(3)=2.0*(z-Ry)
        END SUBROUTINE Calc_ToricNormal

C***********************************************************************

	SUBROUTINE AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,coef,psi,pv,prot,
     &		       pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,nAsphCoef)
	use math_mod

	IMPLICIT NONE

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),Kc,Kr,
     &         pr(3),L,Nvec(3),Nmag,Nhat(3),dNdp(3,3)
	Real*8 :: coef(nAsphCoef)
	Integer :: nAsphCoef

	REAL*8 prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &    k2,itpsi,psitprel,L1,L2,L3,tol,
     &	  z,rho2,h2,h4,h6,h8,h10,zc(3),dfdrho(3),d2fdrho2(3,3),
     &	  rho(3),psitzc,psitrho,fh,Ppsi(3,3),ImpsidF(3,3),
     &    dNvecdp(3,3),hh1,hh2
	REAL*8 D1(9),D2(9),S1,S2
	Real*8, save :: sc(9,2)

	INTEGER i,j
	LOGICAL LROK,ifBuild,ifLNsrf
	Logical, save :: first_entry=.true.

	If (first_entry) then
	  sc(1,1)=4
	  Do i=1,9
	    if (i.ge.2) sc(i,1)=sc(i-1,1)+2
	    sc(i,2)=sc(i,1)*(sc(i,1)-2)
	  End Do
	  first_entry=.false.
	End If

C	Compute incidence point of ray on the surface

	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
   1	CONTINUE
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d22) THEN
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
	  if (.false.) then ! (abs(b) .le. 1d-10) then
	    print*,'AsphSrf(): b = ', b
	    stop
	  end if
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

C	Use L as a starting point, find incidence pt on full surface

        IF (ABS(L).GT.1e-7) THEN
          L1=0.999d0*L
          L2=1.001d0*L
        ELSE
          L1=-1e-7
          L2=1e-7
        END IF

	CALL AZPB(L1,L2,Kc,Kr,psi,pv,prel,ihat,coef,L,nAsphCoef)
	tol=1d-14
	L3=AZPSolve(L1,L2,tol,Kc,Kr,psi,pv,prel,ihat,coef,L,nAsphCoef)
	L=L3

C	Compute normal vector

	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
   	End Do
	psitrho=psi(1)*rho(1)+psi(2)*rho(2)+psi(3)*rho(3)

	z=psi(1)*rho(1)+psi(2)*rho(2)+psi(3)*rho(3)
	rho2=rho(1)*rho(1)+rho(2)*rho(2)+rho(3)*rho(3)
	h2=rho2-z*z
	h4=h2*h2
	h6=h4*h2
	h8=h6*h2
	h10=h8*h2

#if 0
! Original
        fh=coef(1)*h4+coef(2)*h6+coef(3)*h8+coef(4)*h10
        S1=4d0*coef(1)*h2+6d0*coef(2)*h4+8d0*coef(3)*h6+10d0*coef(4)*h8
#else
! jzlou
        hh1=h4; hh2=h2
        fh=0d0; S1=0d0
        Do i=1,nAsphCoef
          fh=fh+coef(i)*hh1
          S1=S1+sc(i,1)*coef(i)*hh2
          hh1=hh1*h2; hh2=hh2*h2
        End Do
#endif

	Do i=1,3
	  zc(i)=rho(i)-fh*psi(i)
	  dFdrho(i)=S1*(rho(i)-psitrho*psi(i))
   	End Do
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)

	S1=Kr+Kc*psitzc
	S2=Kr+(1d0+Kc)*psitzc

	Do i=1,3
	  pr(i)=pin(i)+L*ihat(i)
	  Nvec(i)=zc(i)+S1*psi(i)-S2*dFdrho(i)
    	End Do

	Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  Do i=1,3
	    Nhat(i)=Nvec(i)/Nmag
    	  End Do
	ELSE
	  Do i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
    	  End Do
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

#if 0
! Original
!	  S1=8d0*coef(1)+24d0*coef(2)*h4+48d0*coef(3)*h6+80d0*coef(4)*h8
!	  S2=4d0*coef(1)*h2+6d0*coef(2)*h2+8d0*coef(3)*h4+1d1*coef(4)*h6

          ! should be -jzlou
          S1=8d0*coef(1)+24d0*coef(2)*h2+48d0*coef(3)*h4+80d0*coef(4)*h6
          S2=4d0*coef(1)*h2+6d0*coef(2)*h4+8d0*coef(3)*h6+1d1*coef(4)*h8
          ! S1=(dS2/dh)/h  -jzlou
#else
! jzlou
	  hh1=1d0; hh2=h2
	  S1=0d0; S2=0d0
	  Do i=2,nAsphCoef
            S1=S1+sc(i,2)*coef(i)*hh1
	    S2=S2+sc(i,1)*coef(i)*hh2
            hh1=hh1*h2; hh2=hh2*h2
          End Do
#endif

	  CALL DPROJECT(Ppsi,psi)
	  CALL DMPROD(D1,Ppsi,rho,3,3,1)

	  Do j=1,3
	    Do i=1,3
	      d2Fdrho2(i,j)=S1*D1(i)*D1(j)+S2*Ppsi(i,j)
	      ImpsidF(i,j)=-psi(i)*dFdrho(j)
   	    End Do
	    ImpsidF(j,j)=ImpsidF(j,j)+1d0
   	  End Do
	  CALL DMPROD(D1,ImpsidF,M,3,3,3)
	  CALL DMPROD(D2,D1,ImpsidF,3,3,3)
	  S1=Kr+(1d0+Kc)*psitzc
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
	END SUBROUTINE AsphSrf

C***********************************************************************

        SUBROUTINE AsphGridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,coef,psi,pv,
     &                 prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,nAsphCoef,
     &                 nGridMat,mGridMat,GridMat,GridSrfdx,pmon,
     &                 xhat,yhat,zhat)
	use elt_mod, only : iEltToGridSrf,ifGridDataDefined
        use traceutil_mod, only : crnt_Elt
        use math_mod

        IMPLICIT NONE

        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),Kc,Kr,
     &         pr(3),L,Nvec(3),Nmag,Nhat(3),dNdp(3,3),
     &         Nvec_Asph(3),coef(nAsphCoef),dNdp_asph(3,3),
     &         dNdp2(3,3),L_conic,L_asph
        Integer :: nAsphCoef

        REAL*8 prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &    k2,itpsi,psitprel,L1,L2,L3,tol,
     &    z,rho2,h2,h4,h6,h8,h10,zc(3),dFdrho(3),
     &    d2Fdrho2_asph(3,3),d2Fdrho2_grid(3,3),
     &    rho(3),psitzc,psitrho,fh,Ppsi(3,3),ImpsidF(3,3),
     &    dNvecdp(3,3),hh1,hh2,zc_asph(3),rhoc(3),zmtzc,psitzm,
     &    psitzc_asph,dzcdrho(3,3),grdm_rms2
        REAL*8 D1(9),D2(9),S1,S2,S1_asph,S2_asph
        Real*8, save :: sc(9,2)

        INTEGER i,j
        LOGICAL LROK,ifBuild,ifLNsrf
        Logical, save :: first_entry=.true.

        INTEGER :: nGridMat,mGridMat
        REAL*8 :: GridMat(mGridMat,mGridMat),GridSrfdx,pmon(3),
     &            xhat(3),yhat(3),zhat(3)


!	 CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,coef,psi,pv,
!    &            prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,nAsphCoef)
!        return


#if 0
	 If (.not. ifGridDataDefined(iEltToGridSrf(crnt_Elt))) Then
           ! if grid data are not initialized/provided for this element, just
 	   ! call subroutine AsphSrf() and return
	   CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,coef,psi,pv,
     &            prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,nAsphCoef)
	   !print*,'iElt =',crnt_Elt,', GridMat NOT initialized!'
           return
	 Else
           !print*,'iElt =',crnt_Elt,', GridMat initialized!'
         End If
#endif


        If (first_entry) Then
          sc(1,1)=4
          Do i=1,9
            if (i.ge.2) sc(i,1)=sc(i-1,1)+2
            sc(i,2)=sc(i,1)*(sc(i,1)-2)
          End Do
          first_entry=.false.
        End If

C       Compute incidence point of ray on the surface

        Do i=1,3
          prel(i)=pin(i)-pv(i)
          rhoc(i)=pmon(i)-pv(i)
        End Do
        itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
        psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

        IF (ABS(Kr).LT.1d22) THEN
          mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
          a=1d0+Kc*itpsi*itpsi
          b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &    ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
          c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
        ELSE
          a=0d0
          b=itpsi
          c=psitprel
        END IF

C       Compute intersection point and normal vector

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

        L_conic=L

C       Use L as a starting point, find incidence pt on full (aspherical) surface

        IF (ABS(L).GT.1e-7) THEN
          L1=0.999d0*L
          L2=1.001d0*L
        ELSE
          L1=-1e-7
          L2=1e-7
        END IF

        CALL AZPB(L1,L2,Kc,Kr,psi,pv,prel,ihat,coef,L,nAsphCoef)
        tol=1d-14
        L3=AZPSolve(L1,L2,tol,Kc,Kr,psi,pv,prel,ihat,coef,L,nAsphCoef)
        L=L3
        L_asph=L  ! L_asph is L on top of aspherical surface (conic+pure_aspherical)

C      Compute normal vector

        Do i=1,3
          rho(i)=prel(i)+L*ihat(i)
        End Do
        psitrho=psi(1)*rho(1)+psi(2)*rho(2)+psi(3)*rho(3)

        z=psi(1)*rho(1)+psi(2)*rho(2)+psi(3)*rho(3)
        rho2=rho(1)*rho(1)+rho(2)*rho(2)+rho(3)*rho(3)
        h2=rho2-z*z
        h4=h2*h2
        h6=h4*h2
        h8=h6*h2
        h10=h8*h2
        !
        hh1=h4; hh2=h2
        fh=0d0; S1=0d0
        Do i=1,nAsphCoef
          fh=fh+coef(i)*hh1
          S1=S1+sc(i,1)*coef(i)*hh2
          hh1=hh1*h2; hh2=hh2*h2
        End Do

        Do i=1,3
          zc(i)=rho(i)-fh*psi(i)
          zc_asph(i) = -fh*psi(i)
          dFdrho(i)=S1*(rho(i)-psitrho*psi(i)) ! dFdrho is from pure aspherical
        End Do
        psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
        psitzc_asph
     &    =psi(1)*zc_asph(1)+psi(2)*zc_asph(2)+psi(3)*zc_asph(3)

        S1=Kr+Kc*psitzc
        S2=Kr+(1d0+Kc)*psitzc
        !S1_asph=Kr+Kc*psitzc_asph
        S1_asph=Kc*psitzc_asph
        S2_asph=Kr+(1d0+Kc)*psitzc

#if 0
 ! Test only - use apherical part of Nvec/Nhat  -jzlou
        Do i=1,3
          pr(i)=pin(i)+L*ihat(i)
          Nvec(i)=zc(i)+S1*psi(i)-S2*dFdrho(i)
        End Do

        Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
        S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
        IF (S1.LE.0d0) THEN
          Do i=1,3
            Nhat(i)=Nvec(i)/Nmag
          End Do
        ELSE
          Do i=1,3
            Nhat(i)=-Nvec(i)/Nmag
          End Do
        END IF
#endif

        ! Separate conic and pure_aspherical parts of Nvec

        ! Compute Nvec of pure-aspherical part - Nvec_Asph
        Do i=1,3
          rho(i)=prel(i)+L*ihat(i)    ! rho updated with new 'L' on aspherical -jzlou
          Nvec_Asph(i)=zc_asph(i)+S1_asph*psi(i)-S2_asph*dFdrho(i)
            ! Nvec_Asph is "local surface normal from "pure aspherical" -jzlou
        End Do

#if 0
        ! Take only conic part
        zc(1:3)=rho(1:3) ! only conic part
        psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
        S1=Kr+Kc*psitzc
        DO i=1,3
          pr(i)=pin(i)+L_asph*ihat(i)
          Nvec(i)=Nvec_Asph(i) + zc(i)+S1*psi(i)
        END DO

        Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
        S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
        IF (S1.LE.0d0) THEN
          Do i=1,3
            Nhat(i)=Nvec(i)/Nmag
          End Do
        ELSE
          Do i=1,3
            Nhat(i)=-Nvec(i)/Nmag
          End Do
        END IF
#endif


#if 1
!
! Now add grid surface contribution on top of aspherical surface -jzlou
!

C       grdm_rms2 =
C    &    Calc_RMS2(GridMat(1:nGridMat,1:nGridMat),nGridMat*nGridMat)
C       !print*,'** gridm_rms2 =',grdm_rms2

        ! 'TRUE' = complete code, with grid part included
        ! 'FALSE' = skip the grid part, just use aspherical part
        If (.TRUE.) Then
! ----
          IF (ABS(L).GT.1e-7) THEN
            L1=0.999d0*L_conic
            L2=1.001d0*L_conic
          ELSE
            L1=-1e-7
            L2=1e-7
          END IF

          ! This function updates L1 and L2 -jzlou
          CALL GSZPB(L1,L2,Kc,Kr,psi,pv,prel,ihat,
     &           nGridMat,mGridMat,GridMat,GridSrfdx,
     &           pmon,rhoc,xhat,yhat,zhat,L_conic)

          ! This function updates 'L' for grid surface -jzlou
          tol=1d-10
          L3=GSZPSolve(L1,L2,tol,Kc,Kr,psi,pv,prel,ihat,
     &         nGridMat,mGridMat,GridMat,GridSrfdx,
     &         pmon,rhoc,xhat,yhat,zhat,L_conic)

          L = L_asph+L3-L_conic

C  Compute normal
          CALL NGSrf(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &               nGridMat,mGridMat,GridMat,GridSrfdx,pmon,rhoc,
     &               xhat,yhat,zhat,L,zc,fh,dFdrho,d2Fdrho2_grid)

          psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
          psitzm=psi(1)*zhat(1)+psi(2)*zhat(2)+psi(3)*zhat(3)
          zmtzc=zhat(1)*zc(1)+zhat(2)*zc(2)+zhat(3)*zc(3)
          S1=Kr+Kc*psitzc
          S2=psitzm*Kr+psitzm*psitzc*Kc+zmtzc
! ----
          ! Compute total Nvec = Nvec of pure-aspherical + Nvec of conic+grid
          DO i=1,3
            pr(i)=pin(i)+L*ihat(i)  ! final ray footprint on surface -jzlou
            Nvec(i)=Nvec_Asph(i) + zc(i)+S1*psi(i)-S2*dFdrho(i)
             ! Nvec_Asph(i) is "pure aspherical" contribution to local surface
             ! normal,
             ! "zc(i)+S1*psi(i)-S2*dFdrho(i)" is conic and grid surface
             ! contributions to local surface normal. -jzlou
          End Do
          !
        Else
          L=L_asph  ! skip grid part

          ! Compute Nvec = Nvec of pure-aspherical + Nvec of conic part
          zc(1:3)=rho(1:3) ! only conic part
          psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
          S1=Kr+Kc*psitzc
          DO i=1,3
            pr(i)=pin(i)+L_asph*ihat(i)
            Nvec(i)=Nvec_Asph(i) + zc(i)+S1*psi(i)
          END DO
        End If  ! If .TRUE. or .FALSE.

        Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
        S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
        IF (S1.LE.0d0) THEN
          Do i=1,3
            Nhat(i)=Nvec(i)/Nmag
          End Do
        ELSE
          Do i=1,3
            Nhat(i)=-Nvec(i)/Nmag
          End Do
        END IF


C      Compute partials of the normal

        IF (ifBuild) THEN
          ! Since local surface normal Nvec = Nvec_pure_asph + Nvec_conic + Nvec_grid,
          ! we have
          ! dNdp = dNdp_pure_asph + dNdp_conic + dNdp+grid
          ! -jzlou

          IF (ABS(Kr).LT.1d10) THEN
            CALL DSMPROD(D1,psi,Kc,3)
            CALL DOUTER(M,D1,psi)
            M(1,1)=M(1,1)+1d0
            M(2,2)=M(2,2)+1d0
            M(3,3)=M(3,3)+1d0
          ELSE
            CALL DZERO(M,9)
          END IF

          ! Compute dNdp_pure_asph
          hh1=1d0; hh2=h2
          S1=0d0; S2=0d0
          Do i=2,nAsphCoef
            S1=S1+sc(i,2)*coef(i)*hh1
            S2=S2+sc(i,1)*coef(i)*hh2
            hh1=hh1*h2; hh2=hh2*h2
          End Do

          CALL DPROJECT(Ppsi,psi)
          CALL DMPROD(D1,Ppsi,rho,3,3,1)

          Do j=1,3
            Do i=1,3
              d2Fdrho2_asph(i,j)=S1*D1(i)*D1(j)+S2*Ppsi(i,j)
              ImpsidF(i,j)=-psi(i)*dFdrho(j)
            End Do
            !ImpsidF(j,j)=ImpsidF(j,j)+1d0
            ImpsidF(j,j)=ImpsidF(j,j)
          End Do
          CALL DMPROD(D1,ImpsidF,M,3,3,3)
          CALL DMPROD(D2,D1,ImpsidF,3,3,3)
          S1=Kr+(1d0+Kc)*psitzc
          CALL DSMPROD(D1,d2Fdrho2_asph,S1,9) ! D1 = S1*d2Fdrho2_asph
          CALL DSUB(dNvecdp,D2,D1,9) ! dNvecdp = D2-D1

          S1=1d0/Nmag
          CALL DPROJECT(D1,Nhat)
          CALL DMPROD(D2,D1,dNvecdp,3,3,3)
          CALL DSMPROD(dNdp_asph,D2,S1,9)

          !
          ! Now add conic+grid dNdp to the total dNdp
          !
          Do j=1,3
            Do i=1,3
              dzcdrho(i,j)=-zhat(i)*d2Fdrho2_grid(i,j)
            End Do
            dzcdrho(j,j)=dzcdrho(j,j)+1d0
          End Do
          CALL DMPROD(D1,dzcdrho,M,3,3,3)
          CALL DMPROD(D2,D1,dzcdrho,3,3,3)
          S1=Kr+Kc*psitzc
          CALL DSMPROD(D1,d2Fdrho2_grid,S1,9)
          CALL DSUB(dNvecdp,D2,D1,9)

          S1=1d0/Nmag
          CALL DPROJECT(D1,Nhat)
          CALL DMPROD(D2,D1,dNvecdp,3,3,3)
          CALL DSMPROD(dNdp2,D2,S1,9)

          dNdp = dNdp_asph+dNdp2
        END IF ! ifBuild
#endif

        RETURN

C       Here if ray misses surface

  98    CONTINUE
        LROK=.FALSE.
        RETURN
        END SUBROUTINE AsphGridSrf

C***********************************************************************

	FUNCTION SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,L,nAsphCoef)
	IMPLICIT NONE
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),L,SAsphere
	REAL*8 rho(3),h2,h4,h6,h8,h10,fh,z,rho2,zc,zc2,psitzc,
     &         zcv(3),coef(nAsphCoef),hh
	Integer :: nAsphCoef,i

#if 0
	if (nAsphCoef /= 4) then
          print*,'*** SAsphere() Error: nAsphCoef = ',nAsphCoef
          stop
	else
	  print*,'** SAsphere(): AsphCoef=',coef(1:nAsphCoef)
          stop
        end if
#endif

	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
        END DO

	z=psi(1)*rho(1)+psi(2)*rho(2)+psi(3)*rho(3)
	rho2=rho(1)*rho(1)+rho(2)*rho(2)+rho(3)*rho(3)
	h2=rho2-z*z
	h4=h2*h2
	h6=h4*h2
	h8=h6*h2
	h10=h8*h2
#if 0
! Original
        fh=coef(1)*h4+coef(2)*h6+coef(3)*h8+coef(4)*h10
#else
        hh=h4
        fh=coef(1)*hh
        Do i=2,nAsphCoef
          hh=hh*h2
          fh=fh+coef(i)*hh
        End Do
#endif

	DO i=1,3
	  zcv(i)=rho(i)-fh*psi(i)
        END DO
	zc2=zcv(1)*zcv(1)+zcv(2)*zcv(2)+zcv(3)*zcv(3)
	psitzc=psi(1)*zcv(1)+psi(2)*zcv(2)+psi(3)*zcv(3)
	SAsphere=zc2+Kc*psitzc*psitzc+2d0*Kr*psitzc

	RETURN
	END FUNCTION SAsphere

C***********************************************************************

      REAL*8 FUNCTION AZPSolve(X1,X2,TOL,Kc,Kr,psi,pv,prel,ihat,
     &                           coef,L,nAsphCoef)
        USE Constants, ONLY: EPS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(nAsphCoef),L
      Integer :: nAsphCoef
      PARAMETER (ITMAX=400)
      A=X1
      B=X2
      FA=SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,A,nAsphCoef)
      FB=SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,B,nAsphCoef)
      IF(FB*FA.GT.0.) PAUSE 'AZPB failed.'
      FC=FB
      DO ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2d0*EPS*ABS(B)+0.5d0*TOL
        XM=.5*(C-B)
C        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
        IF((ABS(XM).LE.TOL1).OR.(FB.EQ.0.) .OR.(ITER.GE.ITMAX))THEN
          AZPSolve=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2d0*XM*S
            Q=1d0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2d0*XM*Q*(Q-R)-(B-A)*(R-1d0))
            Q=(Q-1d0)*(R-1d0)*(S-1d0)
          ENDIF
          IF(P.GT.0d0) Q=-Q
          P=ABS(P)
          IF(2d0*P .LT. MIN(3d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,B,nAsphCoef)
      END DO
      PAUSE 'AZPSolve exceeding maximum iterations.'
      AZPSolve=B
      RETURN
      END FUNCTION AZPSolve

C***********************************************************************

	SUBROUTINE AZPB(X1,X2,Kc,Kr,psi,pv,prel,ihat,coef,L,nAsphCoef)
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (gold=1.6d0,NTRY=50)
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(nAsphCoef),L
	Integer :: nAsphCoef
	F1=SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,X1,nAsphCoef)
	F2=SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,X2,nAsphCoef)
	DO J=1,NTRY
	  IF(F1*F2.LT.0.)RETURN
	  IF(ABS(F1).LT.ABS(F2))THEN
	    X1 = X1 + gold*(X1-X2)
	    F1 = SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,X1,nAsphCoef)
	  ELSE
	    X2 = X2 + gold*(X2-X1)
	    F2 = SAsphere(Kc,Kr,psi,pv,prel,ihat,coef,X2,nAsphCoef)
	  ENDIF
        END DO
	RETURN
	END SUBROUTINE AZPB

C***********************************************************************
        ! ---------------------------------------------------------------------
        !           ANSI    n   m   Z(R,theta)      Z(x,y)
        ! zc(1)  |    0  |  0   0 | 1                 1
        ! zc(2)  |    1  |  1  -1 | R*sin(Theta)      y
        ! zc(3)  |    2  |  1   1 | R*cos(theta)      x
        ! zc(4)  |    3  |  2  -2 | R^2*sin(2*theta)  2*x*y
        ! zc(5)  |    4  |  2   0 | 2*R^2 - 1         2(x^2+y^2) - 1
        ! zc(6)  |    5  |  2   2 | R^2*cos(2*theta)  x^2 - y^2

        !---------------------------------|---------------------------------|---------------------------------|
        !   idx  | ANSI  |        |  Norm |   idx  | ANSI  |        |  Norm |   idx  | ANSI  |        |  Norm |
        !        |  idx  |  n   m | Square|        |  idx  |  n   m | Square|        |  idx  |  n   m | Square|
        ! -------|-------|--------|-------| -------|-------|--------|-------| -------|-------|--------|-------|
        !     1  |    0  |  0   0 |    1  |    22  |   21  |  6  -6 |   14  |    46  |   45  |  9  -9 |   20  |
        !---------------------------------|    23  |   22  |  6  -4 |   14  |    47  |   46  |  9  -7 |   20  |
        !     2  |    1  |  1  -1 |    4  |    24  |   23  |  6  -2 |   14  |    48  |   47  |  9  -5 |   20  |
        !     3  |    2  |  1   1 |    4  |    25  |   24  |  6   0 |    7  |    49  |   48  |  9  -3 |   20  |
        !---------------------------------|    26  |   25  |  6   2 |   14  |    50  |   49  |  9  -1 |   20  |
        !     4  |    3  |  2  -2 |    6  |    27  |   26  |  6   4 |   14  |    51  |   50  |  9   1 |   20  |
        !     5  |    4  |  2   0 |    3  |    28  |   27  |  6   6 |   14  |    52  |   51  |  9   3 |   20  |
        !     6  |    5  |  2   2 |    6  |---------------------------------|    53  |   52  |  9   5 |   20  |
        !---------------------------------|    29  |   28  |  7  -7 |   16  |    54  |   53  |  9   7 |   20  |
        !     7  |    6  |  3  -3 |    8  |    30  |   29  |  7  -5 |   16  |    55  |   54  |  9   9 |   20  |
        !     8  |    7  |  3  -1 |    8  |    31  |   30  |  7  -3 |   16  |---------------------------------|
        !     9  |    8  |  3   1 |    8  |    32  |   31  |  7  -1 |   16  |    56  |   55  | 10 -10 |   22  |
        !    10  |    9  |  3   3 |    8  |    33  |   32  |  7   1 |   16  |    57  |   56  | 10  -8 |   22  |
        !---------------------------------|    34  |   33  |  7   3 |   16  |    58  |   57  | 10  -6 |   22  |
        !    11  |   10  |  4  -4 |   10  |    35  |   34  |  7   5 |   16  |    59  |   58  | 10  -4 |   22  |
        !    12  |   11  |  4  -2 |   10  |    36  |   35  |  7   7 |   16  |    60  |   59  | 10  -2 |   22  |
        !    13  |   12  |  4   0 |    5  |---------------------------------|    61  |   60  | 10   0 |   11  |
        !    14  |   13  |  4   2 |   10  |    37  |   36  |  8  -8 |   18  |    62  |   61  | 10   2 |   22  |
        !    15  |   14  |  4   4 |   10  |    38  |   37  |  8  -6 |   18  |    63  |   62  | 10   4 |   22  |
        !---------------------------------|    39  |   38  |  8  -4 |   18  |    64  |   63  | 10   6 |   22  |
        !    16  |   15  |  5  -5 |   12  |    40  |   39  |  8  -2 |   18  |    65  |   64  | 10   8 |   22  |
        !    17  |   16  |  5  -3 |   12  |    41  |   40  |  8   0 |    9  |    66  |   65  | 10  10 |   22  |
        !    18  |   17  |  5  -1 |   12  |    42  |   41  |  8   2 |   18  |---------------------------------|
        !    19  |   18  |  5   1 |   12  |    43  |   42  |  8   4 |   18  |
        !    20  |   19  |  5   3 |   12  |    44  |   43  |  8   6 |   18  |
        !    21  |   20  |  5   5 |   12  |    45  |   44  |  8   8 |   18  |
        !---------------------------------|---------------------------------|

        !      x^i y^(i-j)
        !     --------------------------
        !     ( 0,  0 ) = ( 1   *  1  )
        !     --------------------------
        !     ( 0,  1 ) = ( 1   *  y  )
        !     ( 1,  0 ) = ( x   *  1  )
        !     --------------------------
        !     ( 0,  2 ) = ( 1   *  y^2)
        !     ( 1,  1 ) = ( x   *  y  )
        !     ( 2,  0 ) = ( x^2 *  1  )
        ! ---------------------------------------------------------------------

        ! Mapping ANSI Zernike Ordering to Monomials
        SUBROUTINE ZerntoMon1(zc_in, mc, ie)
          use elt_mod ! , only : ZernTypeL,ZnNmFac_m,mZernModes,mMonCoef,NORM_RMS_PARAM_ANSI

          IMPLICIT NONE
          REAL*8, intent(in)  :: zc_in(mZernModes)   ! Zernike Coefficients (ANSI ordering)
          REAL*8, intent(out) :: mc(mMonCoef)        ! Monomial Coefficients
          INTEGER             :: ie                  ! Element ID

          REAL*8  :: zc(mZernModes)
          ! ----------------------------------------
          zc = zc_in

          ! Apply rms Normalisation Factor
          IF (ANY(ZernType_rms_normalised == ZernTypeL(ie))) zc = zc * NORM_RMS_PARAM_ANSI


          ! ANSI Zernike (index shifted by 1) to Monomials
          mc( 1) = zc(1)-zc(5)+zc(13)-zc(25)+zc(41)-zc(61)

          mc( 2) = zc(2)-2d0*zc(8)+3d0*zc(18)-4d0*zc(32)+5d0*zc(50)
          mc( 3) = zc(3)-2d0*zc(9)+3d0*zc(19)-4d0*zc(33)+5d0*zc(51)

          mc( 4) = 2d0*zc(5)-zc(6)-6d0*zc(13)+3d0*zc(14)+12d0*zc(25)
     &             -6d0*zc(26)-20d0*zc(41)+10d0*zc(42)+30d0*zc(61)-15d0*zc(62)
          mc( 5) = 2d0*zc(4)-6d0*zc(12)+12d0*zc(24)-20d0*zc(40)+30d0*zc(60)
          mc( 6) = 2d0*zc(5)+zc(6)-6d0*zc(13)-3d0*zc(14)+12d0*zc(25)+6d0*zc(26)
     &            -20d0*zc(41)-10d0*zc(42)+30d0*zc(61)+15d0*zc(62)
          mc( 7) = -zc(7)+3d0*zc(8)+4d0*zc(17)-12d0*zc(18)-10d0*zc(31)+30d0*zc(32)+20d0*zc(49)-60d0*zc(50)
          mc( 8) = 3d0*zc(9)-3d0*zc(10)-12d0*zc(19)+12d0*zc(20)+30d0*zc(33)-30d0*zc(34)-60d0*zc(51)+60d0*zc(52)
          mc( 9) = 3d0*zc(7)+3d0*zc(8)-12d0*zc(17)-12d0*zc(18)+30d0*zc(31)+30d0*zc(32)-60d0*zc(49)-60d0*zc(50)
          mc(10) = 3d0*zc(9)+zc(10)-12d0*zc(19)-4d0*zc(20)+30d0*zc(33)+10d0*zc(34)-60d0*zc(51)-20d0*zc(52)

          mc(11) = 6d0*zc(13)-4d0*zc(14)+zc(15)-30d0*zc(25)+20d0*zc(26)-5d0*zc(27)
     &             +90d0*zc(41)-60d0*zc(42)+15d0*zc(43)-210d0*zc(61)+140d0*zc(62)-35d0*zc(63)
          mc(12) = -4d0*zc(11)+8d0*zc(12)+20d0*zc(23)-40d0*zc(24)-60d0*zc(39)+120d0*zc(40)+140d0*zc(59)-280d0*zc(60)
          mc(13) = 12d0*zc(13)-6d0*zc(15)-60d0*zc(25)+30d0*zc(27)+180d0*zc(41)-90d0*zc(43)-420d0*zc(61)+210d0*zc(63)
          mc(14) = 4d0*zc(11)+8d0*zc(12)-20d0*zc(23)-40d0*zc(24)+60d0*zc(39)+120d0*zc(40)-140d0*zc(59)-280d0*zc(60)
          mc(15) = 6d0*zc(13)+4d0*zc(14)+zc(15)-30d0*zc(25)-20d0*zc(26)-5d0*zc(27)+90d0*zc(41)+60d0*zc(42)
     &             +15d0*zc(43)-210d0*zc(61)-140d0*zc(62)-35d0*zc(63)
          mc(16) = zc(16)-5d0*zc(17)+10d0*zc(18)-6d0*zc(30)+30d0*zc(31)
     &             -60d0*zc(32)+21d0*zc(48)-105d0*zc(49)+210d0*zc(50)
          mc(17) = 10d0*zc(19)-15d0*zc(20)+5d0*zc(21)-60d0*zc(33)+90d0*zc(34)
     &            -30d0*zc(35)+210d0*zc(51)-315d0*zc(52)+105d0*zc(53)
          mc(18) = -10d0*zc(16)+10d0*zc(17)+20d0*zc(18)+60d0*zc(30)-60d0*zc(31)
     &             -120d0*zc(32)-210d0*zc(48)+210d0*zc(49)+420d0*zc(50)
          mc(19) = 20d0*zc(19)-10d0*zc(20)-10d0*zc(21)-120d0*zc(33)+60d0*zc(34)
     &             +60d0*zc(35)+420d0*zc(51)-210d0*zc(52)-210d0*zc(53)
          mc(20) = 5d0*zc(16)+15d0*zc(17)+10d0*zc(18)-30d0*zc(30)-90d0*zc(31)
     &             -60d0*zc(32)+105d0*zc(48)+315d0*zc(49)+210d0*zc(50)
          mc(21) = 10d0*zc(19)+5d0*zc(20)+zc(21)-60d0*zc(33)-30d0*zc(34)
     &             -6d0*zc(35)+210d0*zc(51)+105d0*zc(52)+21d0*zc(53)
          mc(22) = 20d0*zc(25)-15d0*zc(26)+6d0*zc(27)-zc(28)-140d0*zc(41)+105d0*zc(42)
     &             -42d0*zc(43)+7d0*zc(44)+560d0*zc(61)-420d0*zc(62)+168d0*zc(63)-28d0*zc(64)
          mc(23) = 6d0*zc(22)-24d0*zc(23)+30d0*zc(24)-42d0*zc(38)+168d0*zc(39)
     &            -210d0*zc(40)+168d0*zc(58)-672d0*zc(59)+840d0*zc(60)
          mc(24) = 60d0*zc(25)-15d0*zc(26)-30d0*zc(27)+15d0*zc(28)-420d0*zc(41)+105d0*zc(42)
     &             +210d0*zc(43)-105d0*zc(44)+1680d0*zc(61)-420d0*zc(62)-840d0*zc(63)+420d0*zc(64)
          mc(25) = -20d0*zc(22)+60d0*zc(24)+140d0*zc(38)-420d0*zc(40)-560d0*zc(58)+1680d0*zc(60)
          mc(26) = 60d0*zc(25)+15d0*zc(26)-30d0*zc(27)-15d0*zc(28)-420d0*zc(41)-105d0*zc(42)
     &             +210d0*zc(43)+105d0*zc(44)+1680d0*zc(61)+420d0*zc(62)-840d0*zc(63)-420d0*zc(64)
          mc(27) = 6d0*zc(22)+24d0*zc(23)+30d0*zc(24)-42d0*zc(38)-168d0*zc(39)-210d0*zc(40)
     &             +168d0*zc(58)+672d0*zc(59)+840d0*zc(60)
          mc(28) = 20d0*zc(25)+15d0*zc(26)+6d0*zc(27)+zc(28)-140d0*zc(41)-105d0*zc(42)-42d0*zc(43)
     &             -7d0*zc(44)+560d0*zc(61)+420d0*zc(62)+168d0*zc(63)+28d0*zc(64)
          mc(29) = -zc(29)+7d0*zc(30)-21d0*zc(31)+35d0*zc(32)+8d0*zc(47)-56d0*zc(48)+168d0*zc(49)-280d0*zc(50)
          mc(30) = 35d0*zc(33)-63d0*zc(34)+35d0*zc(35)-7d0*zc(36)-280d0*zc(51)+504d0*zc(52)-280d0*zc(53)+56d0*zc(54)
          mc(31) = 21d0*zc(29)-63d0*zc(30)+21d0*zc(31)+105d0*zc(32)-168d0*zc(47)+504d0*zc(48)-168d0*zc(49)-840d0*zc(50)
          mc(32) = 105d0*zc(33)-105d0*zc(34)-35d0*zc(35)+35d0*zc(36)-840d0*zc(51)+840d0*zc(52)+280d0*zc(53)-280d0*zc(54)
          mc(33) = -35d0*zc(29)-35d0*zc(30)+105d0*zc(31)+105d0*zc(32)+280d0*zc(47)+280d0*zc(48)-840d0*zc(49)-840d0*zc(50)
          mc(34) = 105d0*zc(33)-21d0*zc(34)-63d0*zc(35)-21d0*zc(36)-840d0*zc(51)+168d0*zc(52)+504d0*zc(53)+168d0*zc(54)
          mc(35) = 7d0*zc(29)+35d0*zc(30)+63d0*zc(31)+35d0*zc(32)-56d0*zc(47)-280d0*zc(48)-504d0*zc(49)-280d0*zc(50)
          mc(36) = 35d0*zc(33)+21d0*zc(34)+7d0*zc(35)+zc(36)-280d0*zc(51)-168d0*zc(52)-56d0*zc(53)-8d0*zc(54)

          mc(37) = 70d0*zc(41)-56d0*zc(42)+28d0*zc(43)-8d0*zc(44)+zc(45)-630d0*zc(61)
     &             +504d0*zc(62)-252d0*zc(63)+72d0*zc(64)-9d0*zc(65)
          mc(38) = -8d0*zc(37)+48d0*zc(38)-112d0*zc(39)+112d0*zc(40)+72d0*zc(57)
     &             -432d0*zc(58)+1008d0*zc(59)-1008d0*zc(60)
          mc(39) = 280d0*zc(41)-112d0*zc(42)-112d0*zc(43)+112d0*zc(44)-28d0*zc(45)
     &             -2520d0*zc(61)+1008d0*zc(62)+1008d0*zc(63)-1008d0*zc(64)+252d0*zc(65)
          mc(40) = 56d0*zc(37)-112d0*zc(38)-112d0*zc(39)+336d0*zc(40)
     &             -504d0*zc(57)+1008d0*zc(58)+1008d0*zc(59)-3024d0*zc(60)
          mc(41) = 420d0*zc(41)-280d0*zc(43)+70d0*zc(45)-3780d0*zc(61)+2520d0*zc(63)-630d0*zc(65)
          mc(42) = -56d0*zc(37)-112d0*zc(38)+112d0*zc(39)+336d0*zc(40)
     &             +504d0*zc(57)+1008d0*zc(58)-1008d0*zc(59)-3024d0*zc(60)
          mc(43) = 280d0*zc(41)+112d0*zc(42)-112d0*zc(43)-112d0*zc(44)-28d0*zc(45)
     &             -2520d0*zc(61)-1008d0*zc(62)+1008d0*zc(63)+1008d0*zc(64)+252d0*zc(65)
          mc(44) = 8d0*zc(37)+48d0*zc(38)+112d0*zc(39)+112d0*zc(40)-72d0*zc(57)
     &             -432d0*zc(58)-1008d0*zc(59)-1008d0*zc(60)
          mc(45) = 70d0*zc(41)+56d0*zc(42)+28d0*zc(43)+8d0*zc(44)+zc(45)-630d0*zc(61)
     &             -504d0*zc(62)-252d0*zc(63)-72d0*zc(64)-9d0*zc(65)

          mc(46) = zc(46)-9d0*zc(47)+36d0*zc(48)-84d0*zc(49)+126d0*zc(50)
          mc(47) = 126d0*zc(51)-252d0*zc(52)+180d0*zc(53)-63d0*zc(54)+9d0*zc(55)
          mc(48) = -36d0*zc(46)+180d0*zc(47)-288d0*zc(48)+504d0*zc(50)
          mc(49) = 504d0*zc(51)-672d0*zc(52)+252d0*zc(54)-84d0*zc(55)
          mc(50) = 126d0*zc(46)-126d0*zc(47)-504d0*zc(48)+504d0*zc(49)+756d0*zc(50)
          mc(51) = 756d0*zc(51)-504d0*zc(52)-504d0*zc(53)+126d0*zc(54)+126d0*zc(55)
          mc(52) = -84d0*zc(46)-252d0*zc(47)+672d0*zc(49)+504d0*zc(50)
          mc(53) = 504d0*zc(51)-288d0*zc(53)-180d0*zc(54)-36d0*zc(55)
          mc(54) = 9d0*zc(46)+63d0*zc(47)+180d0*zc(48)+252d0*zc(49)+126d0*zc(50)
          mc(55) = 126d0*zc(51)+84d0*zc(52)+36d0*zc(53)+9d0*zc(54)+zc(55)

          mc(56) = 252d0*zc(61)-210d0*zc(62)+120d0*zc(63)-45d0*zc(64)+10d0*zc(65)-zc(66)
          mc(57) = 10d0*zc(56)-80d0*zc(57)+270d0*zc(58)-480d0*zc(59)+420d0*zc(60)
          mc(58) = 1260d0*zc(61)-630d0*zc(62)-360d0*zc(63)+585d0*zc(64)-270d0*zc(65)+45d0*zc(66)
          mc(59) = -120d0*zc(56)+480d0*zc(57)-360d0*zc(58)-960d0*zc(59)+1680d0*zc(60)
          mc(60) = 2520d0*zc(61)-420d0*zc(62)-1680d0*zc(63)+630d0*zc(64)+420d0*zc(65)-210d0*zc(66)
          mc(61) = 252d0*zc(56)-1260d0*zc(58)+2520d0*zc(60)
          mc(62) = 2520d0*zc(61)+420d0*zc(62)-1680d0*zc(63)-630d0*zc(64)+420d0*zc(65)+210d0*zc(66)
          mc(63) = -120d0*zc(56)-480d0*zc(57)-360d0*zc(58)+960d0*zc(59)+1680d0*zc(60)
          mc(64) = 1260d0*zc(61)+630d0*zc(62)-360d0*zc(63)-585d0*zc(64)-270d0*zc(65)-45d0*zc(66)
          mc(65) = 10d0*zc(56)+80d0*zc(57)+270d0*zc(58)+480d0*zc(59)+420d0*zc(60)
          mc(66) = 252d0*zc(61)+210d0*zc(62)+120d0*zc(63)+45d0*zc(64)+10d0*zc(65)+zc(66)

          ! set undefined values to zero
          mc((mZernModes+1):) = 0d0

        END SUBROUTINE ZerntoMon1

C***********************************************************************

        ! ---------------------------------------------------------------------
        ! Mapping from Born & Wolf => ANSI Ordering => Monomial Coefficients
        ! ---------------------------------------------------------------------
        SUBROUTINE ZerntoMon2(zc,mc,ie)
          use elt_mod, only : mZernModes, mMonCoef

          IMPLICIT NONE
          REAL*8 zc(mZernModes),mc(mMonCoef),z2c(mZernModes)
          Integer :: ie
          ! ----------------------------------------
                              ! | "ANSI"|  B&W  | (n,m)  |  Zrn
                              ! | Order | Order | n ,  m |
                              ! |  Idx  |  Idx  |  Idx   |
          z2c( 1) = zc( 1)    ! |    0  |    1  |  0   0 | 1

          z2c( 2) = zc( 3)    ! |    1  |    3  |  1  -1 | R*sin(Theta)
          z2c( 3) = zc( 2)    ! |    2  |    2  |  1   1 | R*cos(theta)

          z2c( 4) = zc( 6)    ! |    3  |    6  |  2  -2 | R^2*sin(2*theta)
          z2c( 5) = zc( 5)    ! |    4  |    5  |  2   0 | 2*R^2 - 1
          z2c( 6) = zc( 4)    ! |    5  |    4  |  2   2 | R^2*cos(2*theta)

          z2c( 7) = zc(10)    ! |    6  |   10  |  3  -3 |
          z2c( 8) = zc( 9)    ! |    7  |    9  |  3  -1 |
          z2c( 9) = zc( 8)    ! |    8  |    8  |  3   1 |
          z2c(10) = zc( 7)    ! |    9  |    7  |  3   3 |

          z2c(11) = zc(15)    ! |   10  |   15  |  4  -4 |
          z2c(12) = zc(14)    ! |   11  |   14  |  4  -2 |
          z2c(13) = zc(13)    ! |   12  |   13  |  4   0 |
          z2c(14) = zc(12)    ! |   13  |   12  |  4   2 |
          z2c(15) = zc(11)    ! |   14  |   11  |  4   4 |

          z2c(16) = zc(21)    ! |   15  |   21  |  5  -5 |
          z2c(17) = zc(20)    ! |   16  |   20  |  5  -3 |
          z2c(18) = zc(19)    ! |   17  |   19  |  5  -1 |
          z2c(19) = zc(18)    ! |   18  |   18  |  5   1 |
          z2c(20) = zc(17)    ! |   19  |   17  |  5   3 |
          z2c(21) = zc(16)    ! |   20  |   16  |  5   5 |

          z2c(22) = zc(28)    ! |   21  |   28  |  6  -6 |
          z2c(23) = zc(27)    ! |   22  |   27  |  6  -4 |
          z2c(24) = zc(26)    ! |   23  |   26  |  6  -2 |
          z2c(25) = zc(25)    ! |   24  |   25  |  6   0 |
          z2c(26) = zc(24)    ! |   25  |   24  |  6   2 |
          z2c(27) = zc(23)    ! |   26  |   23  |  6   4 |
          z2c(28) = zc(22)    ! |   27  |   22  |  6   6 |

          z2c(29) = zc(36)    ! |   28  |   36  |  7  -7 |
          z2c(30) = zc(35)    ! |   29  |   35  |  7  -5 |
          z2c(31) = zc(34)    ! |   30  |   34  |  7  -3 |
          z2c(32) = zc(33)    ! |   31  |   33  |  7  -1 |
          z2c(33) = zc(32)    ! |   32  |   32  |  7   1 |
          z2c(34) = zc(31)    ! |   33  |   31  |  7   3 |
          z2c(35) = zc(30)    ! |   34  |   30  |  7   5 |
          z2c(36) = zc(29)    ! |   35  |   29  |  7   7 |

          z2c(37) = zc(45)    ! |   36  |   45  |  8  -8 |
          z2c(38) = zc(44)    ! |   37  |   44  |  8  -6 |
          z2c(39) = zc(43)    ! |   38  |   43  |  8  -4 |
          z2c(40) = zc(42)    ! |   39  |   42  |  8  -2 |
          z2c(41) = zc(41)    ! |   40  |   41  |  8   0 |
          z2c(42) = zc(40)    ! |   41  |   40  |  8   2 |
          z2c(43) = zc(39)    ! |   42  |   39  |  8   4 |
          z2c(44) = zc(38)    ! |   43  |   38  |  8   6 |
          z2c(45) = zc(37)    ! |   44  |   37  |  8   8 |

          z2c(46) = zc(55)    ! |   45  |   55  |  9  -9 |
          z2c(47) = zc(54)    ! |   46  |   54  |  9  -7 |
          z2c(48) = zc(53)    ! |   47  |   53  |  9  -5 |
          z2c(49) = zc(52)    ! |   48  |   52  |  9  -3 |
          z2c(50) = zc(51)    ! |   49  |   51  |  9  -1 |
          z2c(51) = zc(50)    ! |   50  |   50  |  9   1 |
          z2c(52) = zc(49)    ! |   51  |   49  |  9   3 |
          z2c(53) = zc(48)    ! |   52  |   48  |  9   5 |
          z2c(54) = zc(47)    ! |   53  |   47  |  9   7 |
          z2c(55) = zc(46)    ! |   54  |   46  |  9   9 |

          z2c(56) = zc(66)    ! |   55  |   66  | 10 -10 |
          z2c(57) = zc(65)    ! |   56  |   65  | 10  -8 |
          z2c(58) = zc(64)    ! |   57  |   64  | 10  -6 |
          z2c(59) = zc(63)    ! |   58  |   63  | 10  -4 |
          z2c(60) = zc(62)    ! |   59  |   62  | 10  -2 |
          z2c(61) = zc(61)    ! |   60  |   61  | 10   0 |
          z2c(62) = zc(60)    ! |   61  |   60  | 10   2 |
          z2c(63) = zc(59)    ! |   62  |   59  | 10   4 |
          z2c(64) = zc(58)    ! |   63  |   58  | 10   6 |
          z2c(65) = zc(57)    ! |   64  |   57  | 10   8 |
          z2c(66) = zc(56)    ! |   65  |   56  | 10  10 |

          CALL ZerntoMon1(z2c,mc,ie)

        END SUBROUTINE ZerntoMon2

C***********************************************************************

        ! ---------------------------------------------------------------------
        ! Mapping from Extended Fringe => ANSI Ordering => Monomial Coefficients
        ! ---------------------------------------------------------------------
        SUBROUTINE ZerntoMon_ZFE(zc,mc,ie)
          use elt_mod, only : mZernModes, mMonCoef

          IMPLICIT NONE
          REAL*8 zc(mZernModes),mc(mMonCoef),z2c(mZernModes)
          Integer :: ie
          ! ----------------------------------------
          z2c(:) = 0d0

                               !           Extended
                               ! | "ANSI"| Fringe | (n,m)  |
                               ! | Order | Order  | n ,  m |
                               ! |  Idx  |  Idx   |  Idx   |
          z2c( 1) = zc(  1)    ! |    0  |     1  |  0   0 |

          z2c( 2) = zc(  3)    ! |    1  |     3  |  1  -1 |
          z2c( 3) = zc(  2)    ! |    2  |     2  |  1   1 |

          z2c( 4) = zc(  6)    ! |    3  |     6  |  2  -2 |
          z2c( 5) = zc(  4)    ! |    4  |     4  |  2   0 |
          z2c( 6) = zc(  5)    ! |    5  |     5  |  2   2 |

          z2c( 7) = zc( 11)    ! |    6  |    11  |  3  -3 |
          z2c( 8) = zc(  8)    ! |    7  |     8  |  3  -1 |
          z2c( 9) = zc(  7)    ! |    8  |     7  |  3   1 |
          z2c(10) = zc( 10)    ! |    9  |    10  |  3   3 |

          z2c(11) = zc( 18)    ! |   10  |    18  |  4  -4 |
          z2c(12) = zc( 13)    ! |   11  |    13  |  4  -2 |
          z2c(13) = zc(  9)    ! |   12  |     9  |  4   0 |
          z2c(14) = zc( 12)    ! |   13  |    12  |  4   2 |
          z2c(15) = zc( 17)    ! |   14  |    17  |  4   4 |

          z2c(16) = zc( 27)    ! |   15  |    27  |  5  -5 |
          z2c(17) = zc( 20)    ! |   16  |    20  |  5  -3 |
          z2c(18) = zc( 15)    ! |   17  |    15  |  5  -1 |
          z2c(19) = zc( 14)    ! |   18  |    14  |  5   1 |
          z2c(20) = zc( 19)    ! |   19  |    19  |  5   3 |
          z2c(21) = zc( 26)    ! |   20  |    26  |  5   5 |

          z2c(22) = zc( 38)    ! |   21  |    38  |  6  -6 |
          z2c(23) = zc( 29)    ! |   22  |    29  |  6  -4 |
          z2c(24) = zc( 22)    ! |   23  |    22  |  6  -2 |
          z2c(25) = zc( 16)    ! |   24  |    16  |  6   0 |
          z2c(26) = zc( 21)    ! |   25  |    21  |  6   2 |
          z2c(27) = zc( 28)    ! |   26  |    28  |  6   4 |
          z2c(28) = zc( 37)    ! |   27  |    37  |  6   6 |

          z2c(29) = zc( 51)    ! |   28  |    51  |  7  -7 |
          z2c(30) = zc( 40)    ! |   29  |    40  |  7  -5 |
          z2c(31) = zc( 31)    ! |   30  |    31  |  7  -3 |
          z2c(32) = zc( 24)    ! |   31  |    24  |  7  -1 |
          z2c(33) = zc( 23)    ! |   32  |    23  |  7   1 |
          z2c(34) = zc( 30)    ! |   33  |    30  |  7   3 |
          z2c(35) = zc( 39)    ! |   34  |    39  |  7   5 |
          z2c(36) = zc( 50)    ! |   35  |    50  |  7   7 |

          z2c(37) = zc( 66)    ! |   36  |    66  |  8  -8 |
          z2c(38) = zc( 53)    ! |   37  |    53  |  8  -6 |
          z2c(39) = zc( 42)    ! |   38  |    42  |  8  -4 |
          z2c(40) = zc( 33)    ! |   39  |    33  |  8  -2 |
          z2c(41) = zc( 25)    ! |   40  |    25  |  8   0 |
          z2c(42) = zc( 32)    ! |   41  |    32  |  8   2 |
          z2c(43) = zc( 41)    ! |   42  |    41  |  8   4 |
          z2c(44) = zc( 52)    ! |   43  |    52  |  8   6 |
          z2c(45) = zc( 65)    ! |   44  |    65  |  8   8 |

                               ! |   45  |    83  |  9  -9 |
                               ! |   46  |    68  |  9  -7 |
          z2c(48) = zc( 55)    ! |   47  |    55  |  9  -5 |
          z2c(49) = zc( 44)    ! |   48  |    44  |  9  -3 |
          z2c(50) = zc( 35)    ! |   49  |    35  |  9  -1 |
          z2c(51) = zc( 34)    ! |   50  |    34  |  9   1 |
          z2c(52) = zc( 43)    ! |   51  |    43  |  9   3 |
          z2c(53) = zc( 54)    ! |   52  |    54  |  9   5 |
                               ! |   53  |    67  |  9   7 |
                               ! |   54  |    82  |  9   9 |

                               ! |   55  |   102  | 10 -10 |
                               ! |   56  |    85  | 10  -8 |
                               ! |   57  |    70  | 10  -6 |
          z2c(59) = zc( 57)    ! |   58  |    57  | 10  -4 |
          z2c(60) = zc( 46)    ! |   59  |    46  | 10  -2 |
          z2c(61) = zc( 36)    ! |   60  |    36  | 10   0 |
          z2c(62) = zc( 45)    ! |   61  |    45  | 10   2 |
          z2c(63) = zc( 56)    ! |   62  |    56  | 10   4 |
                               ! |   63  |    69  | 10   6 |
                               ! |   64  |    84  | 10   8 |
                               ! |   65  |   101  | 10  10 |

          CALL ZerntoMon1(z2c,mc,ie)

        END SUBROUTINE ZerntoMon_ZFE

C***********************************************************************

        ! ---------------------------------------------------------------------
        ! Mapping from Fringe => ANSI Ordering => Monomial Coefficients
        ! ---------------------------------------------------------------------
        SUBROUTINE ZerntoMon3(zc,mc,ie)
          use elt_mod,   only: mZernModes, mMonCoef, ZernType_rms_normalised, ZernTypeL
          use Constants, only: EPS

	  IMPLICIT NONE
          REAL*8,  INTENT(IN)   :: zc(mZernModes)  ! Zernike Coefficients: Fringe
          REAL*8,  INTENT(INOUT):: mc(mMonCoef)    ! Monomial Coefficients: relay to Monomial Mapping
          INTEGER, INTENT(IN)   :: ie              ! Element ID

          REAL*8 :: z2c(mZernModes)
          ! ----------------------------------------
          z2c(:) = 0d0   ! reset

                              ! | "ANSI"| Fringe| (n,m)  |  Zrn
                              ! | Order | Order | n ,  m |
                              ! |  Idx  |  Idx  |  Idx   |
          z2c( 1) = zc( 1)    ! |    0  |    1  |  0   0 |

          z2c( 2) = zc( 3)    ! |    1  |    3  |  1  -1 |
          z2c( 3) = zc( 2)    ! |    2  |    2  |  1   1 |

          z2c( 4) = zc( 6)    ! |    3  |    6  |  2  -2 |
          z2c( 5) = zc( 4)    ! |    4  |    4  |  2   0 |
          z2c( 6) = zc( 5)    ! |    5  |    5  |  2   2 |

          z2c( 7) = zc(11)    ! |    6  |   11  |  3  -3 |
          z2c( 8) = zc( 8)    ! |    7  |    8  |  3  -1 |
          z2c( 9) = zc( 7)    ! |    8  |    7  |  3   1 |
          z2c(10) = zc(10)    ! |    9  |   10  |  3   3 |

          z2c(11) = zc(18)    ! |   10  |   18  |  4  -4 |
          z2c(12) = zc(13)    ! |   11  |   13  |  4  -2 |
          z2c(13) = zc( 9)    ! |   12  |    9  |  4   0 |
          z2c(14) = zc(12)    ! |   13  |   12  |  4   2 |
          z2c(15) = zc(17)    ! |   14  |   17  |  4   4 |

          z2c(16) = zc(27)    ! |   15  |   27  |  5  -5 |
          z2c(17) = zc(20)    ! |   16  |   20  |  5  -3 |
          z2c(18) = zc(15)    ! |   17  |   15  |  5  -1 |
          z2c(19) = zc(14)    ! |   18  |   14  |  5   1 |
          z2c(20) = zc(19)    ! |   19  |   19  |  5   3 |
          z2c(21) = zc(26)    ! |   20  |   26  |  5   5 |

                              ! |   21  |       |  6  -6 |
          z2c(23) = zc(29)    ! |   22  |   29  |  6  -4 |
          z2c(24) = zc(22)    ! |   23  |   22  |  6  -2 |
          z2c(25) = zc(16)    ! |   24  |   16  |  6   0 |
          z2c(26) = zc(21)    ! |   25  |   21  |  6   2 |
          z2c(27) = zc(28)    ! |   26  |   28  |  6   4 |
                              ! |   27  |       |  6   6 |

                              ! |   28  |       |  7  -7 |
                              ! |   29  |       |  7  -5 |
          z2c(31) = zc(31)    ! |   30  |   31  |  7  -3 |
          z2c(32) = zc(24)    ! |   31  |   24  |  7  -1 |
          z2c(33) = zc(23)    ! |   32  |   23  |  7   1 |
          z2c(34) = zc(30)    ! |   33  |   30  |  7   3 |
                              ! |   34  |       |  7   5 |
                              ! |   35  |       |  7   7 |

                              ! |   36  |       |  8  -8 |
                              ! |   37  |       |  8  -6 |
                              ! |   38  |       |  8  -4 |
          z2c(40) = zc(33)    ! |   39  |   33  |  8  -2 |
          z2c(41) = zc(25)    ! |   40  |   25  |  8   0 |
          z2c(42) = zc(32)    ! |   41  |   32  |  8   2 |
                              ! |   42  |       |  8   4 |
                              ! |   43  |       |  8   6 |
                              ! |   44  |       |  8   8 |

                              ! |   45  |       |  9  -9 |
                              ! |   46  |       |  9  -7 |
                              ! |   47  |       |  9  -5 |
                              ! |   48  |       |  9  -3 |
          z2c(50) = zc(35)    ! |   49  |   35  |  9  -1 |
          z2c(51) = zc(34)    ! |   50  |   34  |  9   1 |
                              ! |   51  |       |  9   3 |
                              ! |   52  |       |  9   5 |
                              ! |   53  |       |  9   7 |
                              ! |   54  |       |  9   9 |

                              ! |   55  |       | 10 -10 |
                              ! |   56  |       | 10  -8 |
                              ! |   57  |       | 10  -6 |
                              ! |   58  |       | 10  -4 |
                              ! |   59  |       | 10  -2 |
          z2c(61) = zc(36)    ! |   60  |   36  | 10   0 |
                              ! |   61  |       | 10   2 |
                              ! |   62  |       | 10   4 |
                              ! |   63  |       | 10   6 |
                              ! |   64  |       | 10   8 |
                              ! |   65  |       | 10  10 |
                              ! ...
                              ! |   84  |   37  | 12   0 | <<< Special

	  CALL ZerntoMon1(z2c,mc,ie)

	  ! Special: update Monomials for Z37  (radial n=12, azimuthal m=0)
          !   ==> CODE V Manual, p751: LensSystemSetupRM.pdf
          !
          ! ANSI: j=84  n=12  m=0
          !  --> Z84 = C84 * (1-42 R^2 + 420*R^4-1680*R^6 + 3150*R^8 - 2772*R^10 + 924*R^12)

          IF (ABS(zc(37)) > EPS) THEN

            z2c(1) = zc(37)

            ! check if rms normalisation is to be applied
            IF (ANY(ZernType_rms_normalised == ZernTypeL(ie))) z2c(1) = z2c(1)*sqrt(13d0)

                                                !  x^j * y^(i-j)
            mc( 1) = mc( 1) +         z2c(1)    !  { 0 ,  0 }

            mc( 4) = mc( 4) -    42d0*z2c(1)    !  { 0 ,  2 }
            mc( 6) = mc( 6) -    42d0*z2c(1)    !  { 2 ,  0 }

            mc(11) = mc(11) +   420d0*z2c(1)    !  { 0 ,  4 }
            mc(13) = mc(13) +   840d0*z2c(1)    !  { 2 ,  2 }
            mc(15) = mc(15) +   420d0*z2c(1)    !  { 4 ,  0 }

            mc(22) = mc(22) -  1680d0*z2c(1)    !  { 0 ,  6 }
            mc(24) = mc(24) -  5040d0*z2c(1)    !  { 2 ,  4 }
            mc(26) = mc(26) -  5040d0*z2c(1)    !  { 4 ,  2 }
            mc(28) = mc(28) -  1680d0*z2c(1)    !  { 6 ,  0 }

            mc(37) = mc(37) +  3150d0*z2c(1)    !  { 0 ,  8 }
            mc(39) = mc(39) + 12600d0*z2c(1)    !  { 2 ,  6 }
            mc(41) = mc(41) + 18900d0*z2c(1)    !  { 4 ,  4 }
            mc(43) = mc(43) + 12600d0*z2c(1)    !  { 6 ,  2 }
            mc(45) = mc(45) +  3150d0*z2c(1)    !  { 8 ,  0 }

            mc(56) = mc(56) -  2772d0*z2c(1)    !  { 0 , 10 }
            mc(58) = mc(58) - 13860d0*z2c(1)    !  { 2 ,  8 }
            mc(60) = mc(60) - 27720d0*z2c(1)    !  { 4 ,  6 }
            mc(62) = mc(62) - 27720d0*z2c(1)    !  { 6 ,  4 }
            mc(64) = mc(64) - 13860d0*z2c(1)    !  { 8 ,  2 }
            mc(66) = mc(66) -  2772d0*z2c(1)    !  {10 ,  0 }

            mc(79) = mc(79) +   924d0*z2c(1)    !  { 0,  12 }
            mc(81) = mc(81) +  5544d0*z2c(1)    !  { 2,  10 }
            mc(83) = mc(83) + 13860d0*z2c(1)    !  { 4,   8 }
            mc(85) = mc(85) + 18480d0*z2c(1)    !  { 6,   6 }
            mc(87) = mc(87) + 13860d0*z2c(1)    !  { 8,   4 }
            mc(89) = mc(89) +  5544d0*z2c(1)    !  {10,   2 }
            mc(91) = mc(91) +   924d0*z2c(1)    !  {12,   0 }

          END IF

	END SUBROUTINE ZerntoMon3

C***********************************************************************

       Subroutine ZerntoMon4(Zc,mc)
       use elt_mod  !, only : ZernTypeL
       Implicit none

       Real*8 Zc(1:mZernModes),mc(mMonCoef)
       Integer, parameter :: mzern=45,mpolyterm=15
       Integer, save :: nterm(mzern),hzi(mzern,mpolyterm,2),
     &                  mij2k(0:14,0:14),iz,i,j,k,jm,ii,jj
       Real*8, save :: hz(mpolyterm,mzern)
       Logical, save :: first_entry=.true.

       nterm=(/1,1,1,3,1,2,3,3,2,2,6,4,3,5,3,
     &         6,6,5,5,6,6,10,6,9,6,9,3,10,10,10,
     &         9,9,10,10,10,10,15,13,10,14,10,15,7,14,10/)

       If (first_entry) Then
#include "hzdef.inc"
          k=0; jm=0
          Do i=0,14
            ii=0; jj=jm
            Do j=0,i
              k=k+1
              mij2k(ii,jj)=k
              ii=ii+1; jj=jj-1
            End Do
            jm=jm+1
          End Do
          first_entry=.false.
        End If

        mc=0d0
        Do i=1,mzern
          Do j=1,nterm(i)
            k=mij2k(hzi(i,j,1),hzi(i,j,2))
            mc(k)=mc(k)+Zc(i)*hz(j,i)
          End Do
        End Do
        Return
	End Subroutine ZerntoMon4

C***********************************************************************

	Subroutine ZerntoMon5(Zc,mc)
        Implicit None

        REAL*8 :: Zc(4),mc(120)

        Integer, parameter :: mzern=4,maxpolyterms=28
        Integer, save :: nterm(mzern),fzi(mzern,maxpolyterms,2),
     &                   mij2k(0:14,0:14),i,j,k,jm,ii,jj,iz,ioff
        Real*8, save :: fz(maxpolyterms,mzern)
        Logical, save :: first_entry=.true.

        nterm=(/15,15,21,28/)

        If (first_entry) Then
C ---
#include "fzdef.inc"
C ---
	  k=0; jm=0
          Do i=0,14
            ii=0; jj=jm
            Do j=0,i
              k=k+1
              mij2k(ii,jj)=k
              ii=ii+1; jj=jj-1
            End Do
            jm=jm+1
          End Do
          first_entry=.false.
        End If

	Do i=1,mzern
          Do j=1,nterm(i)
            k=mij2k(fzi(i,j,1),fzi(i,j,2))
            mc(k)=mc(k)+Zc(i)*fz(j,i)
          End Do
        End Do
        End Subroutine ZerntoMon5

C***********************************************************************

        ! ---------------------------------------------------------------------
        ! Zernike-Order Mapping: Noll ==> ANSI => Monomials
        ! ---------------------------------------------------------------------
        SUBROUTINE ZerntoMon6(zc,mc,ie)
          use elt_mod, only: mZernModes, mMonCoef

          IMPLICIT NONE
          REAL*8,  INTENT(IN)   :: zc(mZernModes)  ! Zernike Coefficients: Noll
          REAL*8,  INTENT(INOUT):: mc(mMonCoef)    ! Monomial Coefficients: relay to Monomial Mapping
          INTEGER, INTENT(IN)   :: ie              ! Element ID

          REAL*8 :: z2c(mZernModes)
          ! ----------------------------------------

                              ! | "ANSI"|  Noll | (n,m)  |
                              ! | Order | Order | n ,  m |
                              ! |  Idx  |  Idx  |  Idx   |
          z2c( 1) = zc( 1)    ! |    0  |    1  |  0   0 |

          z2c( 2) = zc( 3)    ! |    1  |    3  |  1  -1 |
          z2c( 3) = zc( 2)    ! |    2  |    2  |  1   1 |

          z2c( 4) = zc( 5)    ! |    3  |    5  |  2  -2 |
          z2c( 5) = zc( 4)    ! |    4  |    4  |  2   0 |
          z2c( 6) = zc( 6)    ! |    5  |    6  |  2   2 |

          z2c( 7) = zc( 9)    ! |    6  |    9  |  3  -3 |
          z2c( 8) = zc( 7)    ! |    7  |    7  |  3  -1 |
          z2c( 9) = zc( 8)    ! |    8  |    8  |  3   1 |
          z2c(10) = zc(10)    ! |    9  |   10  |  3   3 |

          z2c(11) = zc(15)    ! |   10  |   15  |  4  -4 |
          z2c(12) = zc(13)    ! |   11  |   13  |  4  -2 |
          z2c(13) = zc(11)    ! |   12  |   11  |  4   0 |
          z2c(14) = zc(12)    ! |   13  |   12  |  4   2 |
          z2c(15) = zc(14)    ! |   14  |   14  |  4   4 |

          z2c(16) = zc(21)    ! |   15  |   21  |  5  -5 |
          z2c(17) = zc(19)    ! |   16  |   19  |  5  -3 |
          z2c(18) = zc(17)    ! |   17  |   17  |  5  -1 |
          z2c(19) = zc(16)    ! |   18  |   16  |  5   1 |
          z2c(20) = zc(18)    ! |   19  |   18  |  5   3 |
          z2c(21) = zc(20)    ! |   20  |   20  |  5   5 |

          z2c(22) = zc(27)    ! |   21  |   27  |  6  -6 |
          z2c(23) = zc(25)    ! |   22  |   25  |  6  -4 |
          z2c(24) = zc(23)    ! |   23  |   23  |  6  -2 |
          z2c(25) = zc(22)    ! |   24  |   22  |  6   0 |
          z2c(26) = zc(24)    ! |   25  |   24  |  6   2 |
          z2c(27) = zc(26)    ! |   26  |   26  |  6   4 |
          z2c(28) = zc(28)    ! |   27  |   28  |  6   6 |

          z2c(29) = zc(35)    ! |   28  |   35  |  7  -7 |
          z2c(30) = zc(33)    ! |   29  |   33  |  7  -5 |
          z2c(31) = zc(31)    ! |   30  |   31  |  7  -3 |
          z2c(32) = zc(29)    ! |   31  |   29  |  7  -1 |
          z2c(33) = zc(30)    ! |   32  |   30  |  7   1 |
          z2c(34) = zc(32)    ! |   33  |   32  |  7   3 |
          z2c(35) = zc(34)    ! |   34  |   34  |  7   5 |
          z2c(36) = zc(36)    ! |   35  |   36  |  7   7 |

          z2c(37) = zc(45)    ! |   36  |   45  |  8  -8 |
          z2c(38) = zc(43)    ! |   37  |   43  |  8  -6 |
          z2c(39) = zc(41)    ! |   38  |   41  |  8  -4 |
          z2c(40) = zc(39)    ! |   39  |   39  |  8  -2 |
          z2c(41) = zc(37)    ! |   40  |   37  |  8   0 |
          z2c(42) = zc(38)    ! |   41  |   38  |  8   2 |
          z2c(43) = zc(40)    ! |   42  |   40  |  8   4 |
          z2c(44) = zc(42)    ! |   43  |   42  |  8   6 |
          z2c(45) = zc(44)    ! |   44  |   44  |  8   8 |

          z2c(46) = zc(55)    ! |   45  |   55  |  9  -9 |
          z2c(47) = zc(53)    ! |   46  |   53  |  9  -7 |
          z2c(48) = zc(51)    ! |   47  |   51  |  9  -5 |
          z2c(49) = zc(49)    ! |   48  |   49  |  9  -3 |
          z2c(50) = zc(47)    ! |   49  |   47  |  9  -1 |
          z2c(51) = zc(46)    ! |   50  |   46  |  9   1 |
          z2c(52) = zc(48)    ! |   51  |   48  |  9   3 |
          z2c(53) = zc(50)    ! |   52  |   50  |  9   5 |
          z2c(54) = zc(52)    ! |   53  |   52  |  9   7 |
          z2c(55) = zc(54)    ! |   54  |   54  |  9   9 |

          z2c(56) = zc(65)    ! |   55  |   65  | 10 -10 |
          z2c(57) = zc(63)    ! |   56  |   63  | 10  -8 |
          z2c(58) = zc(61)    ! |   57  |   61  | 10  -6 |
          z2c(59) = zc(59)    ! |   58  |   59  | 10  -4 |
          z2c(60) = zc(57)    ! |   59  |   57  | 10  -2 |
          z2c(61) = zc(56)    ! |   60  |   56  | 10   0 |
          z2c(62) = zc(58)    ! |   61  |   58  | 10   2 |
          z2c(63) = zc(60)    ! |   62  |   60  | 10   4 |
          z2c(64) = zc(62)    ! |   63  |   62  | 10   6 |
          z2c(65) = zc(64)    ! |   64  |   64  | 10   8 |
          z2c(66) = zc(66)    ! |   65  |   66  | 10  10 |

          CALL ZerntoMon1(z2c,mc,ie)

        END SUBROUTINE ZerntoMon6

C***********************************************************************

       Subroutine ZerntoMon7(Zc,mc,r)
        use elt_mod, only: mZernModes, mMonCoef
        Implicit none

       Real*8 :: Zc(mZernModes),mc(mMonCoef),r
       Integer, parameter :: mzern=45,mpolyterm=15
       Integer, save :: nterm(mzern),azi(mzern,mpolyterm,2),
     &                  mij2k(0:14,0:14),iz,i,j,k,jm,ii,jj
       Real*8, save :: az(mpolyterm,mzern)
       Logical, save :: first_entry=.true.

       nterm=(/1,1,1,3,1,2,3,3,2,2,6,4,3,3,2,6,6,5,5,3,
     &         3,10,6,8,4,7,3,4,10,10,9,9,7,7,4,4,15,12,
     &         10,12,8,8,7,5,4/)

       If (first_entry) Then
#include "azdef.inc"
          k=0; jm=0
          Do i=0,14
            ii=0; jj=jm
            Do j=0,i
              k=k+1
              mij2k(ii,jj)=k
              ii=ii+1; jj=jj-1
            End Do
            jm=jm+1
          End Do
          first_entry=.false.
        End If

        mc=0d0
        Do i=1,mzern
          Do j=1,nterm(i)
            k=mij2k(azi(i,j,1),azi(i,j,2))
            mc(k)=mc(k)+Zc(i)*az(j,i)
          End Do
        End Do
        Return
        End Subroutine ZerntoMon7

C***********************************************************************

        SUBROUTINE MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,coef,psi,pv,prot,
     &    pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
        use math_mod
        use Constants, only: EPS

        IMPLICIT NONE

        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),coef(120),
     &  pr(3),rad,pmon(3),xhat(3),yhat(3),zhat(3),L,Nhat(3),yi,yimj,xj

        REAL*8 prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &  Nvec(3),Nmag,dNdp(3,3),Kc,Kr,
     &  k2,itpsi,psitprel,L1,L2,L3,tol,
     &  zc(3),dfdrho(3),d2fdrho2(3,3),
     &  rho(3),psitzc,fh,dzcdrho(3,3),dNvecdp(3,3),x,
     &  y,dxdrho(3),dydrho(3),rhoc(3),rhom(3),zmtzc,psitzm, fh_dz

        REAL*8 D1(9),D2(9),S1,S2,S3

        INTEGER i,j,k

        LOGICAL LROK,ifBuild,ifLNsrf

        ! Compute incidence point of ray on the surface

        DO i=1,3
          prel(i)=pin(i)-pv(i)
          rhoc(i)=pmon(i)-pv(i)
        END DO
        itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
        psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

        IF (ABS(Kr).LT.1d10) THEN
          mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
          a=1d0+Kc*itpsi*itpsi
          b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &      ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
          c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
        ELSE
          a=0d0
          b=itpsi
          c=psitprel
        END IF

        ! Compute intersection point and normal vector

        IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d10)) THEN
          L=-c/b
          IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98
        ELSE
          k2=b*b-4d0*a*c
          IF (k2.LT.0d0) GO TO 98
          IF (b.GE.0d0) THEN
            ksqrt=-0.5d0*(b+SQRT(k2))
          ELSE
            ksqrt=-0.5d0*(b-SQRT(k2))
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

         ! Use L as a starting point, find incidence pt on full surface

        IF (ABS(L).GT.1d-7) THEN
          L1=0.999d0*L
          L2=1.001d0*L
        ELSE
          L1=-1d-7
          L2=1d-7
        END IF

        CALL MZPB(L1,L2,Kc,Kr,psi,pv,prel,ihat,coef,rad,
     &            pmon,rhoc,xhat,yhat,zhat,L)
        tol=2d-15
        L3=MZPSolve(L1,L2,tol,Kc,Kr,psi,pv,prel,ihat,coef,
     &              rad,pmon,rhoc,xhat,yhat,zhat,L)
        L=L3

        ! Compute normal vector
        rho  = prel + L*ihat
        rhom = rho - rhoc

        x = (xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))/rad
        y = (yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))/rad

        k  = 1
        S1 = 0d0
        S2 = 0d0
        fh = COEF(1)

        IF ((ABS(x) <= EPS).AND.(ABS(y) <= EPS)) THEN

          S1 = S1 + coef(3)
          S2 = S2 + coef(2)

        ELSE IF (ABS(x) <= EPS) THEN

            yi=y
            DO CONCURRENT (i=1:14)
              DO CONCURRENT (j=0:i)

                k=k+1
                IF (j>1) CYCLE

              IF (ABS(COEF(k)) > EPS) THEN

                  IF (j.EQ.0) THEN
                    fh_dz = coef(k)*yi
                    fh    = fh + fh_dz
                    S2    = S2 + DBLE(i)*fh_dz/y
                  END IF

                IF (j.EQ.1) S1 = S1 + coef(k)*yi/y

                END IF
                END DO
              yi=yi*y
            END DO

          ELSE IF (ABS(y) <= EPS) THEN

            DO CONCURRENT (i=1:14)
              xj=1d0
              DO CONCURRENT (j=0:i)
                k=k+1

              IF (ABS(COEF(k)) > EPS) THEN

                IF (i.EQ.j) THEN
                  fh_dz = coef(k)*xj
                  fh    = fh + fh_dz
                  S1    = S1 + DBLE(j)*fh_dz/x
                END IF

                IF (i-j.EQ.1) S2 = S2 + coef(k)*xj

                END IF
                xj=xj*x
              END DO
            END DO

          ELSE

            yi=y
            DO CONCURRENT (i=1:14)
              xj=1d0
              yimj=yi
              DO CONCURRENT (j=0:i)
                k=k+1
                IF (ABS(COEF(k)) > EPS) THEN
                  IF (j.EQ.0) THEN
                    fh_dz = coef(k)*yi
                  fh    = fh + fh_dz
                    S2    = S2 + DBLE(i)*fh_dz/y
                  ELSE IF (i.EQ.j) THEN
                    fh_dz = coef(k)*xj
                  fh    = fh + fh_dz
                    S1    = S1 + DBLE(j)*fh_dz/x
                  ELSE
                    fh_dz = coef(k)*xj*yimj
                  fh    = fh + fh_dz
                    S1    = S1 + DBLE(j)*fh_dz/x
                    S2    = S2 + DBLE(i-j)*fh_dz/y
                  END IF
                END IF
                yimj=yimj/y
                xj=xj*x
              END DO
              yi=yi*y
          END DO
        END IF

        dxdrho  = xhat / rad
        dydrho  = yhat / rad
        dFdrho  = S1*dxdrho + S2*dydrho
        zc      = rho - fh*zhat

        psitzc = psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
        psitzm = psi(1)*zhat(1)+psi(2)*zhat(2)+psi(3)*zhat(3)
        zmtzc  = zhat(1)*zc(1)+zhat(2)*zc(2)+zhat(3)*zc(3)

        S1 = Kr + Kc*psitzc
        S2 = psitzm*Kr + psitzm*psitzc*Kc + zmtzc

        pr   = pin + L*ihat
        Nvec = zc + S1*psi - S2*dFdrho

        Nmag = SQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
        S1   = Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
        IF (S1.LE.0d0) THEN
          Nhat = Nvec/Nmag
        ELSE
          Nhat = -Nvec/Nmag
        END IF

        ! Compute partials of the normal

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

          k=1
          S1=0d0
          S2=0d0
          S3=0d0
          DO i=2,14
            DO j=0,i
              k=k+1
              IF (coef(k).NE.0d0) THEN
                IF (j.GT.1) THEN
                  S1=S1+coef(k)*DBLE(j*(j-1))*(x**(j-2))*(y**(i-j))
                END IF
                IF (j.LT.(i-1)) THEN
                  S2=S2+coef(k)*DBLE((i-j)*(i-j-1))*(x**j)*(y**(i-j-2))
                END IF
                IF ((j.GT.0).AND.(j.LT.i)) THEN
                  S3=S3+coef(k)*DBLE(j*(i-j))*(x**(j-1))*(y**(i-j-1))
                END IF
              END IF
            END DO
          END DO

          DO j=1,3
            DO i=1,3
              d2Fdrho2(i,j)=S1*dxdrho(i)*dxdrho(j)+
     &        S2*dydrho(i)*dydrho(j)+S3*(dxdrho(i)*dydrho(j)+
     &        dydrho(i)*dxdrho(j))
              dzcdrho(i,j)=-zhat(i)*dFdrho(j)
            END DO
            dzcdrho(j,j)=dzcdrho(j,j)+1d0
          END DO
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
        ! Here if ray misses surface

  98    CONTINUE
	LROK=.FALSE.
	RETURN
	END SUBROUTINE MonSrf

C***********************************************************************

        REAL*8 FUNCTION SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,
     &                            rhoc,xhat,yhat,zhat,L)
          use Constants, only: EPS
          IMPLICIT NONE
          REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(120),xhat(3),
     &           yhat(3),zhat(3),rad,pmon(3),L,x,y,yi,yimj,xj
          REAL*8 rho(3),rhom(3),rhoc(3),fh,zc(3),zc2,psitzc

          INTEGER i,j,k

          ! Compute normalized point-of-incidence coordinates

          !                   Monomials
          !  k   i   j   i-j   x     y
          !  ===========================
          !  1   0   0    0    1  *  1     i==j
          !  ---------------------------
          !  2   1   0    1    1  *  y     j==0
          !  3   1   1    0    x  *  1     i==j, j==1
          !  ---------------------------
          !  4   2   0    2    1   * y^2   j==0
          !  5   2   1    1    x   * y           j==1
          !  6   2   2    0    x^2 * 1     i==j
          !  ---------------------------
          !  7   3   0    3    1   * y^3   j==0
          !  8   3   1    2    x   * y^2         j==1
          !  9   3   2    1    x^2 * y
          ! 10   3   3    0    x^3 * 1     i==j
          !  ---------------------------

          rho  = prel + L*ihat
          rhom = rho - rhoc

          x = (xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))/rad
          y = (yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))/rad

          k  = 0
          fh = 0d0

          IF ((ABS(x)<=EPS).AND.(ABS(y)<=EPS)) THEN

            fh = coef(1)

          ELSE IF (ABS(x) <= EPS) THEN

            yi=1d0
            DO CONCURRENT (i=0:14)
              DO CONCURRENT (j=0:i)
                k=k+1
                IF ((j.EQ.0).AND.(ABS(COEF(k)) > EPS)) fh=fh+coef(k)*yi
              END DO
              yi=yi*y
            END DO

          ELSE IF (ABS(y) <= EPS) THEN

            DO CONCURRENT (i=0:14)
              xj=1d0
              DO CONCURRENT (j=0:i)
                k=k+1
                IF ((ABS(coef(k)) > EPS).AND.(i.EQ.j)) fh=fh+coef(k)*xj
                xj=xj*x
              END DO
            END DO

          ELSE

            yi=1d0
            DO CONCURRENT (i=0:14)
              xj=1d0
              yimj=yi
              DO CONCURRENT (j=0:i)
                k=k+1
                fh=fh+coef(k)*xj*yimj
                xj=xj*x
                yimj=yimj/y
              END DO
              yi=yi*y
            END DO
          END IF

          zc        = rho - fh*zhat
          zc2       = zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
          psitzc    = psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
          SMonomial = zc2 + (Kc*psitzc + 2d0*Kr)*psitzc

          RETURN
        END FUNCTION SMonomial

C***********************************************************************
! Check: Numerical Recipies: pg 404

        REAL*8 FUNCTION MZPSolve(X1,X2,TOL,Kc,Kr,psi,pv,prel,ihat,coef,
     &                           rad,pmon,rhoc,xhat,yhat,zhat,L)
          use Constants, only: EPS

          IMPLICIT REAL*8 (A-H,O-Z)
          REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(120),L,
     &           rad,pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3)
          PARAMETER (ITMAX=400) !,EPS=1d-14)
          A=X1
          B=X2
          FA=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &                 xhat,yhat,zhat,A)
          FB=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &                 xhat,yhat,zhat,B)
          IF(FB*FA.GT.0d0) PAUSE 'MZPB failed.'
          FC=FB
          DO ITER=1,ITMAX
            IF(FB*FC.GT.0d0) THEN
              C=A
              FC=FA
              D=B-A
              E=D
            ENDIF
            IF(ABS(FC).LT.ABS(FB)) THEN
              A=B
              B=C
              C=A
              FA=FB
              FB=FC
              FC=FA
            ENDIF
            TOL1=2d0*EPS*ABS(B)+0.5d0*TOL
            !   TOL1=2d0*EPS*ABS(B)+0.5d0*TOL
            !   TOL1 = EPS_def
              XM=0.5d0*(C-B)
            IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0d0) THEN
              MZPSolve=B
              RETURN
            ENDIF
            IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S=FB/FA
              IF(A.EQ.C) THEN
                P=2d0*XM*S
                Q=1d0-S
              ELSE
                Q=FA/FC
                R=FB/FC
                P=S*(2d0*XM*Q*(Q-R)-(B-A)*(R-1d0))
                Q=(Q-1d0)*(R-1d0)*(S-1d0)
              ENDIF
              IF(P.GT.0d0) Q=-Q
              P=ABS(P)
              IF(2d0*P .LT. MIN(3d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                E=D
                D=P/Q
              ELSE
                D=XM
                E=D
              ENDIF
            ELSE
              D=XM
              E=D
            ENDIF
            A=B
            FA=FB
            IF(ABS(D) .GT. TOL1) THEN
              B=B+D
            ELSE
              B=B+SIGN(TOL1,XM)
            ENDIF
            FB=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &                   xhat,yhat,zhat,B)
          END DO
          PAUSE 'MZPSolve exceeding maximum iterations.'
          MZPSolve=B
          RETURN
        END FUNCTION MZPSolve

C***********************************************************************

        SUBROUTINE MZPB(X1,X2,Kc,Kr,psi,pv,prel,ihat,coef,rad,
     &                  pmon,rhoc,xhat,yhat,zhat,L)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (gold=1.6d0,NTRY=50)
        REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(120),L,
     &         rad,pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3)
        IF(X1.EQ.X2)PAUSE 'You have to guess an initial range'
        F1=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &               xhat,yhat,zhat,X1)
        F2=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &               xhat,yhat,zhat,X2)
        DO J=1,NTRY
          IF(F1*F2.LT.0d0)RETURN
          IF(ABS(F1).LT.ABS(F2))THEN
            X1=X1+gold*(X1-X2)
            F1=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &                   xhat,yhat,zhat,X1)
          ELSE
            X2=X2+gold*(X2-X1)
            F2=SMonomial(Kc,Kr,psi,pv,prel,ihat,coef,rad,pmon,rhoc,
     &                   xhat,yhat,zhat,X2)
          ENDIF
        END DO
        RETURN
        END SUBROUTINE MZPB

C***********************************************************************

	SUBROUTINE UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,coef,
     &	  	         param,psi,pv,prot,pin,ihat,pr,L,rad,
     &	  	         pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	use math_mod

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifLNsrf
	INTEGER i,j,UDSrfType
	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),coef(*),
     &   pr(3),rad,pmon(3),xhat(3),yhat(3),zhat(3),L,Nhat(3),
     &   prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &   Nvec(3),Nmag,dNdp(3,3),Kc,Kr,
     &   k2,itpsi,psitprel,L1,L2,L3,tol,
     &	 zc(3),dFdrho(3),d2Fdrho2(3,3),
     &	 psitzc,fh,dzcdrho(3,3),dNvecdp(3,3),rhoc(3),zmtzc,psitzm,
     &   D1(9),D2(9),S1,S2,param(6)

C	Compute incidence point of ray on the surface

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

        ! Compute intersection point and normal vector

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

        ! Use L as a starting point, find incidence pt on full surface

        IF (ABS(L).GT.1e-7) THEN
          L1=0.999d0*L
          L2=1.001d0*L
        ELSE
          L1=-1d-7
          L2=1d-7
        END IF

	CALL UDSZPB(L1,L2,UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,
     &		    param,rad,pmon,rhoc,xhat,yhat,zhat,L)
	tol=1d-10
	L3=UDSZPSolve(L1,L2,tol,UDSrfType,Kc,Kr,psi,pv,prel,ihat,
     &	              coef,param,rad,pmon,rhoc,xhat,yhat,zhat,L)
	L=L3

        ! Compute normal

	CALL UDSrfNorm(ifBuild,UDSrfType,Kc,Kr,psi,pv,prel,ihat,
     &	       coef,param,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,
     &	       fh,dFdrho,d2Fdrho2)

	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	psitzm=psi(1)*zhat(1)+psi(2)*zhat(2)+psi(3)*zhat(3)
	zmtzc=zhat(1)*zc(1)+zhat(2)*zc(2)+zhat(3)*zc(3)
	S1=Kr+Kc*psitzc
	S2=psitzm*Kr+psitzm*psitzc*Kc+zmtzc

	DO i=1,3
	  pr(i)=pin(i)+L*ihat(i)
	  Nvec(i)=zc(i)+S1*psi(i)-S2*dFdrho(i)
        END DO
	Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  DO i=1,3
	    Nhat(i)=Nvec(i)/Nmag
          END DO
	ELSE
	  DO i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
          END DO
	END IF

        ! Compute partials of the normal

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

	  DO j=1,3
	    DO i=1,3
	      dzcdrho(i,j)=-zhat(i)*d2Fdrho2(i,j)
            END DO
	    dzcdrho(j,j)=dzcdrho(j,j)+1d0
          END DO
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
	END SUBROUTINE UDSrf

C***********************************************************************


      FUNCTION UDSZPSolve(X1,X2,TOL,UDSrfType,Kc,Kr,psi,pv,prel,ihat,
     &	coef,param,rad,pmon,rhoc,xhat,yhat,zhat,L)

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER UDSrfType
      REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),L,SUDSrf,
     &	UDSZPSolve,rad,pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3),param(6)
      PARAMETER (ITMAX=400,EPS=1d-10)
      A=X1
      B=X2
      FA=SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &          rad,pmon,rhoc,xhat,yhat,zhat,A)
      FB=SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &          rad,pmon,rhoc,xhat,yhat,zhat,B)
      IF(FB*FA.GT.0d0) PAUSE 'UDSZPB failed.'
      FC=FB
      DO ITER=1,ITMAX
        IF(FB*FC.GT.0d0) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2d0*EPS*ABS(B)+0.5d0*TOL
        XM=0.5d0*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0d0)THEN
          UDSZPSolve=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2d0*XM*S
            Q=1d0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2d0*XM*Q*(Q-R)-(B-A)*(R-1d0))
            Q=(Q-1d0)*(R-1d0)*(S-1d0)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2d0*P .LT. MIN(3d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &		  rad,pmon,rhoc,xhat,yhat,zhat,B)
      END DO
      PAUSE 'UDSZPSolve exceeding maximum iterations.'
      UDSZPSolve=B
      RETURN
      END FUNCTION UDSZPSolve

C***********************************************************************

	SUBROUTINE UDSZPB(X1,X2,UDSrfType,Kc,Kr,psi,pv,prel,ihat,
     &                    coef,param,rad,pmon,rhoc,xhat,yhat,zhat,L)

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (gold=1.6d0,NTRY=50)
	INTEGER UDSrfType
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),L,
     &         rad,pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3),param(6)
	IF(X1.EQ.X2)PAUSE 'You have to guess an initial range'
	F1=SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &            rad,pmon,rhoc,xhat,yhat,zhat,X1)
	F2=SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &            rad,pmon,rhoc,xhat,yhat,zhat,X2)
	DO J=1,NTRY
	  IF(F1*F2.LT.0d0)RETURN
	  IF(ABS(F1).LT.ABS(F2))THEN
	    X1 = X1 + gold*(X1-X2)
	    F1 = SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &                  rad,pmon,rhoc,xhat,yhat,zhat,X1)
	  ELSE
	    X2 = X2 + gold*(X2-X1)
	    F2 = SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,coef,param,
     &                  rad,pmon,rhoc,xhat,yhat,zhat,X2)
	  ENDIF
        END DO
	RETURN
	END SUBROUTINE UDSZPB

C***********************************************************************

 	SUBROUTINE GridInit(GridFile,jGridSrf,nGridMatS,ifGridInit,
     &                      isBinFil)
	use param_mod
        use elt_mod
	use math_mod

	IMPLICIT NONE

	CHARACTER(len=*) :: GridFile
	LOGICAL exist,ifGridInit
	Logical :: isBinFil
	INTEGER jGridSrf,nGridMatS,i,j,iDat,ICLEN

  103	FORMAT(' Grid surface array is of dimension',i5,' by ',i5)

C   If no file, zero matrix and return
        IF ((GridFile(1:1).EQ.' ') .OR.
     &      (GridFile(1:ICLEN(GridFile)).EQ.'none') .OR.
     &      (GridFile(1:ICLEN(GridFile)).EQ.'NONE') .OR.
     &      (GridFile(1:ICLEN(GridFile)).EQ.'None')) THEN
          GridFile='none'
          WRITE(*,*)' No grid matrix provided'
          CALL DZERO(GridMat(1,1,jGridSrf),mGridMat*mGridMat)
          ifGridInit=.FALSE.
          RETURN
        END IF

        ! Check file existence

	!print*,'B4 INQ: GridFile=',GridFile(1:ICLEN(GridFile))

	INQUIRE (FILE=GridFile,EXIST=exist)
	IF (exist) THEN
	  if (.not.isBinFil) then
            OPEN (UNIT=94,FILE=GridFile,STATUS='OLD',
     &            FORM='FORMATTED')
	  else
	    ! binary input
	    OPEN (UNIT=94,FILE=GridFile,STATUS='OLD',
     &            FORM='UNFORMATTED')
	  end if
	ELSE
	  WRITE(*,*)' File ',GridFile(1:ICLEN(GridFile)),
     &	  ' does not exist -sub GridInit in surfsub.f'
	  CALL DZERO(GridMat(1,1,jGridSrf),(mGridMat*mGridMat))
	  ifGridInit=.FALSE.
	  RETURN
	END IF

C  Read input data
	WRITE(*,*)' Reading surface grid data file ',
     &            GridFile(1:ICLEN(GridFile))
	WRITE(*,*)' mGridMat,nGridMatS =', mGridMat,nGridMatS
        iDat=0
	!DO j=1,mGridMat
	DO j=1,nGridMatS
	  iDat=iDat+1
	  if (.not.isBinFil) then
	    READ(94,*,ERR=2,END=1)
     &        (GridMat(i,j,jGridSrf),i=1,nGridMatS)
	  else
            ! read binary data
	    print*,'* Bin read: reading a line'
	    !READ(94,ERR=2,END=1)
	    READ(94,END=1)
     &         GridMat(1:nGridMatS,j,jGridSrf)
	  end if
	END DO

  1	CONTINUE
	ifGridInit=.TRUE.
	WRITE(*,103) nGridMatS,iDat
	CLOSE(94)
	!print*,'**** GridMat=',GridMat(1:96,1:96,jGridSrf)
	RETURN

  2	CONTINUE
	ifGridInit=.FALSE.
	WRITE(*,*) ' Error during read GridData file!'
	CALL DZERO(GridMat(1,1,jGridSrf),(mGridMat*mGridMat))
	CLOSE(94)
	RETURN
	END SUBROUTINE GridInit

C***********************************************************************

	SUBROUTINE AmplInit(AmplFile,jAmplSrf,nAmplMatS,ifAmplInit)
	use param_mod
	use elt_mod
	use math_mod

	IMPLICIT NONE

	CHARACTER(len=*) :: AmplFile
        LOGICAL exist,ifAmplInit,LCMP
        INTEGER jAmplSrf,nAmplMatS,i,j,iDat,ICLEN

	If (LCMP(AmplFile(1:1),' ',1) .OR.
     &      LCMP(AmplFile(1:ICLEN(AmplFile)),'none',4)) Then
	  WRITE(*,*)' No amplitude matrix provided'
	  CALL DZERO(AmplMat(1,1,jAmplSrf),mAmplMat*mAmplMat)
          ifAmplInit=.FALSE.
          RETURN
	End If

	INQUIRE (FILE=AmplFile,EXIST=exist)
	IF (exist) THEN
          OPEN (UNIT=94,FILE=AmplFile,STATUS='OLD',
     &          FORM='FORMATTED')
        ELSE
          WRITE(*,*)' File ',AmplFile(1:ICLEN(AmplFile)),
     &      ' does not exist - subroutine AmplInit in surfsub.f'
          CALL DZERO(AmplMat(1,1,jAmplSrf),(mAmplMat*mAmplMat))
          ifAmplInit=.FALSE.
          RETURN
        END IF

C Read amplitude data
	WRITE(*,*)' Reading surface amplitude file ',
     &  	  AmplFile(1:ICLEN(AmplFile))
        WRITE(*,*)' mAmplMat,nAmplMatS =',mAmplMat,nAmplMatS

	iDat=0
        DO j=1,nAmplMatS
          iDat=iDat+1
	  READ(94,*,ERR=2,END=1)
     &        (AmplMat(i,j,jAmplSrf),i=1,nAmplMatS)
        END DO

  1     CONTINUE
        ifAmplInit=.TRUE.
        WRITE(*,104) nAmplMatS,iDat-1

 104    FORMAT(' Amplitude array is of dimension',i5,' by ',i5)
        CLOSE(94)
        RETURN

  2     CONTINUE
        ifAmplInit=.FALSE.
        WRITE(*,*) ' Error during read surface amplitude file!'
	CALL DZERO(AmplMat(1,1,jAmplSrf),(mAmplMat*mAmplMat))
        CLOSE(94)
        RETURN
	END SUBROUTINE AmplInit

C***********************************************************************

	SUBROUTINE GridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,
     &	  psi,pv,prot,pin,ihat,pr,L,
     &	  nGridMat,mGridMat,GridMat,GridSrfdx,
     &	  pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	use macos_debug
	use math_mod

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifLNsrf
	INTEGER i,j,nGridMat,mGridMat
	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &         pr(3),pmon(3),xhat(3),yhat(3),zhat(3),L,Nhat(3),
     &         prel(3),M(3,3),a,b,c,ksqrt,kp,km,mpr,
     &         Nvec(3),Nmag,dNdp(3,3),Kc,Kr,
     &         k2,itpsi,psitprel,L1,L2,L3,tol,
     &	       GridMat(mGridMat,mGridMat),GridSrfdx,
     &	       zc(3),dFdrho(3),d2Fdrho2(3,3),psitzc,fh,
     &	       dzcdrho(3,3),dNvecdp(3,3),rhoc(3),zmtzc,psitzm,
     &         D1(9),D2(9),S1,S2

	if (.false.) then
	  print*,'**Begin GridSrf: GridSrfdx = ',GridSrfdx
	  print*,'**Begin GridSrf: pmon = ',pmon(1:3)
	  print*,'**Begin GridSrf: nGridMat = ',nGridMat
	end if

        ! Compute incidence point of ray on the surface

	DO i=1,3
	  prel(i)=pin(i)-pv(i)
	  rhoc(i)=pmon(i)-pv(i)
        END DO
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

        ! Compute intersection point and normal vector

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

	if (.false.) then !(iRay_debug == 2526 .and.
C    &      iElt_debug == 90) then
          print*,' ****==> iRay = 2526, L = ', L,
     &      'a,b,c,kp,km,Kc = ', a, b, c, kp, km,Kc
          print*,'****==>2526, ihat, psi, itpsi =',
     &      ihat, psi, itpsi
        end if

        !  Use L as a starting point, find incidence pt on full surface

        IF (ABS(L).GT.1e-7) THEN
          L1=0.999d0*L
          L2=1.001d0*L
        ELSE
          L1=-1d-7
          L2=1d-7
        END IF

	CALL GSZPB(L1,L2,Kc,Kr,psi,pv,prel,ihat,
     &	       nGridMat,mGridMat,GridMat,GridSrfdx,
     &	       pmon,rhoc,xhat,yhat,zhat,L)
	tol=1d-10
	L3=GSZPSolve(L1,L2,tol,Kc,Kr,psi,pv,prel,ihat,
     &	     nGridMat,mGridMat,GridMat,GridSrfdx,
     &	     pmon,rhoc,xhat,yhat,zhat,L)
	L=L3

C  Compute normal
	CALL NGSrf(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  nGridMat,mGridMat,GridMat,GridSrfdx,
     &	  pmon,rhoc,xhat,yhat,zhat,L,zc,
     &	  fh,dFdrho,d2Fdrho2)

	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	psitzm=psi(1)*zhat(1)+psi(2)*zhat(2)+psi(3)*zhat(3)
	zmtzc=zhat(1)*zc(1)+zhat(2)*zc(2)+zhat(3)*zc(3)
	S1=Kr+Kc*psitzc
	S2=psitzm*Kr+psitzm*psitzc*Kc+zmtzc

	DO i=1,3
	  pr(i)=pin(i)+L*ihat(i)
	  Nvec(i)=zc(i)+S1*psi(i)-S2*dFdrho(i)
        END DO
	Nmag=DSQRT(Nvec(1)*Nvec(1)+Nvec(2)*Nvec(2)+Nvec(3)*Nvec(3))
	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.LE.0d0) THEN
	  DO i=1,3
	    Nhat(i)=Nvec(i)/Nmag
   	  END DO
	ELSE
	  DO i=1,3
	    Nhat(i)=-Nvec(i)/Nmag
   	  END DO
	END IF


        ! Compute partials of the normal

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

	  DO j=1,3
	    DO i=1,3
	      dzcdrho(i,j)=-zhat(i)*d2Fdrho2(i,j)
  	    END DO
	    dzcdrho(j,j)=dzcdrho(j,j)+1d0
  	  END DO
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

        ! Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.
	RETURN
	END SUBROUTINE GridSrf

C***********************************************************************

	FUNCTION GSZPSolve(X1,X2,TOL,Kc,Kr,psi,pv,prel,ihat,
     &	  nGridMat,mGridMat,GridMat,GridSrfdx,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

	IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER nGridMat,mGridMat
	REAL*8 GridMat(mGridMat,mGridMat),GridSrfdx
	INTEGER UDSrfType
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),L,
     &	  GSZPSolve,rad,pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3)
	  PARAMETER (ITMAX=400,EPS=1d-10)
	A=X1
	B=X2
	FA=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &		 nGridMat,mGridMat,GridMat,GridSrfdx,
     &		 pmon,rhoc,xhat,yhat,zhat,A)
	FB=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &		 nGridMat,mGridMat,GridMat,GridSrfdx,
     &		 pmon,rhoc,xhat,yhat,zhat,B)
	IF(FB*FA.GT.0.) PAUSE 'GSZPB failed.'
	FC=FB
	DO ITER=1,ITMAX
	  IF(FB*FC.GT.0d0) THEN
	    C=A; FC=FA; D=B-A; E=D
	  ENDIF
	  IF(ABS(FC).LT.ABS(FB)) THEN
	    A=B; B=C; C=A
	    FA=FB; FB=FC; FC=FA
	  ENDIF
	  TOL1=2d0*EPS*ABS(B)+0.5d0*TOL
	  XM=0.5d0*(C-B)
	  IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0d0) THEN
	    GSZPSolve=B
	    RETURN
	  ENDIF
	  IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
	    S=FB/FA
	    IF(A.EQ.C) THEN
	      P=2d0*XM*S
	      Q=1d0-S
	    ELSE
	      Q=FA/FC
	      R=FB/FC
	      P=S*(2d0*XM*Q*(Q-R)-(B-A)*(R-1d0))
	      Q=(Q-1d0)*(R-1d0)*(S-1d0)
	    ENDIF
	    IF(P.GT.0.) Q=-Q
	    P=ABS(P)
	    IF(2d0*P .LT. MIN(3d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
	      E=D
	      D=P/Q
	    ELSE
	      D=XM
	      E=D
	    ENDIF
	  ELSE
	    D=XM
	    E=D
	  ENDIF
	  A=B
	  FA=FB
	  IF(ABS(D) .GT. TOL1) THEN
	    B=B+D
	  ELSE
	    B=B+SIGN(TOL1,XM)
	  ENDIF
	  FB=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &	  	   nGridMat,mGridMat,GridMat,GridSrfdx,
     &	  	   pmon,rhoc,xhat,yhat,zhat,B)
        END DO
	PAUSE 'GSZPSolve exceeding maximum iterations.'
	GSZPSolve=B
	RETURN
	END FUNCTION GSZPSolve

C***********************************************************************

	SUBROUTINE GSZPB(X1,X2,Kc,Kr,psi,pv,prel,ihat,
     &	  nGridMat,mGridMat,GridMat,GridSrfdx,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

	  IMPLICIT REAL*8 (A-H,O-Z)
	  INTEGER nGridMat,mGridMat
	  REAL*8 GridMat(mGridMat,mGridMat),GridSrfdx
	  PARAMETER (gold=1.6,NTRY=50)
	  INTEGER UDSrfType
	  REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),L,
     &	         pmon(3),xhat(3),yhat(3),zhat(3),rhoc(3)

	IF(X1.EQ.X2)PAUSE 'GSZPB initial range error'
	F1=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &	         nGridMat,mGridMat,GridMat,GridSrfdx,
     &	         pmon,rhoc,xhat,yhat,zhat,X1)
	F2=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &	         nGridMat,mGridMat,GridMat,GridSrfdx,
     &	         pmon,rhoc,xhat,yhat,zhat,X2)
	DO J=1,NTRY
	  IF(F1*F2.LT.0.)RETURN
	  IF(ABS(F1).LT.ABS(F2))THEN
	    X1=X1+gold*(X1-X2)
	    F1=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &	             nGridMat,mGridMat,GridMat,GridSrfdx,
     &	             pmon,rhoc,xhat,yhat,zhat,X1)
	  ELSE
	    X2=X2+gold*(X2-X1)
	    F2=SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &	             nGridMat,mGridMat,GridMat,GridSrfdx,
     &	             pmon,rhoc,xhat,yhat,zhat,X2)
	  ENDIF
        END DO
	RETURN
	END SUBROUTINE GSZPB

C***********************************************************************

	FUNCTION SGSrf(Kc,Kr,psi,pv,prel,ihat,
     &	               nGridMat,mGridMat,GridMat,dAct,
     &		       pmon,rhoc,xhat,yhat,zhat,L)
        use param_mod, only : mGridSrf
        use elt_mod, only : iEltToGridSrf,GridSrfOrder
        use traceutil_mod, only : crnt_Elt
	use math_mod

C  Computes surface function for gridded data surface

	IMPLICIT NONE

	INTEGER nGridMat,mGridMat,gridSrfIdx
	REAL*8 GridMat(mGridMat,mGridMat)

	INTEGER i,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),xhat(3),
     &	  yhat(3),zhat(3),pmon(3),L,SGSrf,rho(3),rhom(3),rhoc(3),
     &	  fh,dAct,x,y,zc(3),zc2,psitzc,tiltx,tilty,xi,yj,diCtr

        ! Compute point-of-incidence in surface coords

	diCtr=INT((nGridMat+1)/2)
	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
	END DO

        ! xi is non-integer coordinate in 1-99 frame, similar for yi
	xi=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))/dAct+diCtr
	yj=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))/dAct+diCtr

        ! Find grid coordinates of ray incidence point

	i0=IDFLOOR(xi) ! floor integer part of xi
	i1=i0+1
	j0=IDFLOOR(yj)
	j1=j0+1
	x=xi-DBLE(i0)  ! fractional part, ray position from i0
	y=yj-DBLE(j0)


        ! Compute excess surface height

	IF ((i0.LT.1).OR.(i1.GT.nGridMat).OR.(j0.LT.1).OR.
     &	    (j1.GT.nGridMat)) THEN
	  fh=0d0
	  tiltx=0d0
          tilty=0d0
	ELSE
          If (GridSrfOrder(crnt_Elt)==1) Then
            ! bi-linear interpolation of GridMat
	    CALL INTNORM(x,y,GridMat(i0,j0),GridMat(i0,j1),GridMat(i1,j0),
     &	    	         GridMat(i1,j1),fh,tiltx,tilty)
          Else IF (GridSrfOrder(crnt_Elt)==3) Then
	    ! cubic-spline interpolation of GridMat
            gridSrfIdx=iEltToGridSrf(crnt_Elt)
!           print*,'******* mGridSrf,crnt_Elt,gridSrfIdx =',
!    &	           mGridSrf,crnt_Elt,gridSrfIdx
            x=xi*dAct; y=yj*dAct
            CALL INTCUBIC(gridSrfIdx,crnt_Elt,
     &                    x,y,GridMat(1:nGridMat,1:nGridMat),
     &                    mGridMat,nGridMat,mGridSrf,
     &                    i0,j0,fh,tiltx,tilty)
	  Else
	    WRITE(*,*)'*** Function SGSrf: Invalid SrfOrder: ',
     &                GridSrfOrder(crnt_Elt)
            Stop
          End If
	END IF

        ! Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SGSrf=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SGSrf

C***********************************************************************

	SUBROUTINE NGSrf(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  nGridMat,mGridMat,GridMat,dAct,
     &	  pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use macos_debug
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER nGridMat,mGridMat
	REAL*8 GridMat(mGridMat,mGridMat)

	LOGICAL ifBuild
	INTEGER i,iAxis,jAxis,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),xhat(3),
     &	       yhat(3),zhat(3),pmon(3),L,rho(3),rhom(3),rhoc(3),
     &	       fh,x,y,zc(3),dAct,xi,yj,
     &	       dFdrho(3),d2Fdrho2(3,3),tiltx,tilty,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=INT((nGridMat+1)/2)
	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
	END DO
	xi=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))/dAct+diCtr
	yj=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))/dAct+diCtr

C  Find actuator coordinates of ray incidence point

	i0=IDFLOOR(xi)
	i1=i0+1
	j0=IDFLOOR(yj)
	j1=j0+1
	x=xi-DBLE(i0)
	y=yj-DBLE(j0)

C  Compute excess surface height

	IF ((i0.LT.1).OR.(i1.GT.nGridMat).OR.(j0.LT.1).OR.
     &	(j1.GT.nGridMat)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,GridMat(i0,j0),GridMat(i0,j1),GridMat(i1,j0),
     &	    GridMat(i1,j1),fh,tiltx,tilty)
	END IF

	tiltx=tiltx/dAct
	tilty=tilty/dAct

	DO jAxis=1,3
 	  dFdrho(jAxis)=tiltx*xhat(jAxis)+tilty*yhat(jAxis)
 	  DO iAxis=1,3
   	    d2Fdrho2(iAxis,jAxis)=0d0
   	  END DO
   	END DO

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NGSrf

C***********************************************************************

#if 0
	SUBROUTINE FindConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,
     &                        pin,ihat,pr,L)
	use math_mod

	IMPLICIT NONE

	INTEGER i
	LOGICAL LROK,ifBuild,ifLNsrf
	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &  pr(3),L,S1
	REAL*8 prel(3),a,b,c,ksqrt,kp,km,mpr,Kc,Kr,
     &  k2,itpsi,psitprel

        ! Compute incidence point of ray on the surface

	DO i=1,3
	  prel(i)=pin(i)-pv(i)
        END DO
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d22) THEN
	  mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
	  a=1d0+Kc*itpsi*itpsi
	  b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &    ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
	  c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
	ELSE
	  a=0d0
	  b=itpsi
	  c=psitprel
	END IF

        ! Compute intersection point

	IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d10)) THEN
	  L=-c/b
	ELSE
	  k2=b*b-4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  ksqrt=DSQRT(k2)
	  S1=ABS(ABS(b)-ksqrt)
	  IF (S1.LT.1d-4) THEN
	    kp=(ksqrt-b)/(2d0*a)
	    km=(-ksqrt-b)/(2d0*a)
	  ELSE
	    kp=(-2d0*c)/(b+ksqrt)
	    km=(-2d0*c)/(b-ksqrt)
	  END IF
	  IF (kp.LT.km) THEN
	    IF (kp.GT.1d-12) THEN
	      L=kp
	    ELSE
	      L=km
	    END IF
	  ELSE
	    IF (km.GT.1d-12) THEN
	      L=km
	    ELSE
	      L=kp
	    END IF
	  END IF
	END IF

        ! LnegOK -- Nonsequential surfaces must have positive L!!
	IF (L.LT.1d-12) GO TO 98

	RETURN

        ! Here if ray misses surface

  98    CONTINUE
	LROK=.FALSE.
	RETURN
	END SUBROUTINE FindConSrf

#endif
C***********************************************************************

        END MODULE surfsub
