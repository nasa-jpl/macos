C***********************************************************************
C	Begin file usersub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************


#include "realtype.h"

	SUBROUTINE UDSinit(ReadOK,iElt)
C  UDSinit provides a shell for user-defined surface function 
C  intialization routines. Which routine is invoked depends on variable 
C  UDSrfType. Current offerings are:
C
C	1: 1313-actuator deformable mirror.
C	2: 349-actuator deformable mirror.
C	3: 397-actuator triangular array deformable mirror.
C	4: 349-actuator deformable mirror with squared-out
C	   influence function.
C	5: 349-actuator deformable mirror with ideal actuator
C	   influence functions (bilinear approximation).
C	6: 349-actuator deformable mirror with squared-out
C	   influence function.
C	7: 149-by-149 x-y grid for Hubble Primary Mirror.
C	8: 115-by-115 x-y grid for Hubble Secondary Mirror.
C	9: Radial-actuator DM.
C      10: SQAna-by-SQAna DM with damped-sine influence fcns.
C      12: 185-by-185 x-y grid

	use param_mod
	use cfiles_mod
	use elt_mod
	use math_mod

	IMPLICIT NONE

	LOGICAL ReadOK,IHSTPM,IHSTSM,IDMRAD,IGridSrf,ISQARRAYDM,EXIST
	INTEGER i,j,iElt,SQAna
	REAL*8 S1,S2

	IF (UDSrfType(iElt).EQ.7) THEN
	  ReadOK=IHSTPM(UDSrfFile(iElt))

	ELSE IF (UDSrfType(iElt).EQ.8) THEN
	  ReadOK=IHSTSM(UDSrfFile(iElt))

	ELSE IF (UDSrfType(iElt).EQ.9) THEN
	  ReadOK=IDMRAD(UDSrfFile(iElt))

	ELSE IF (UDSrfType(iElt).EQ.10) THEN
	  SQAna=UDSrfParam(3,iElt)
	  ReadOK=ISQARRAYDM(UDSrfFile(iElt),SQAna,UDSrfCoef(1,iElt))

	  write(*,*)'SQAna =',SQAna

	ELSE IF (UDSrfType(iElt).EQ.11) THEN
	  ReadOK=IGridSrf(UDSrfFile(iElt))

	ELSE IF ((UDSrfType(iElt).GE.1).AND.(UDSrfType(iElt).LT.9)) 
     &  THEN
	  INQUIRE (FILE=UDSrfFile(iElt),EXIST=EXIST)
	  CALL DZERO(UDSrfCoef(1,iElt),mUDScoef)
	  IF (EXIST) THEN
	    OPEN (UNIT=27,FILE=UDSrfFile(iElt),STATUS='OLD',
     &	    FORM='FORMATTED')
	    DO i=1,mUDScoef   ! mUDScoef defined in macos_param.txt
	      READ(27,*,END=2,ERR=1)S1,S2
	      j=S1
	      UDSrfCoef(j,iElt)=S2
	    END DO
	    GO TO 2

 1	    CONTINUE
	    WRITE(*,*)' Error reading UDSrfFile file'
	    ReadOK=.FALSE.
	    RETURN

 2	    CONTINUE
	    CLOSE(27)
	    ReadOK=.TRUE.
	  END IF

	ELSE
	  WRITE(*,*)' Unknown user-defined surface type.'
	  ReadOK=.FALSE.
	END IF

	RETURN
	END SUBROUTINE UDSinit

C***********************************************************************

	FUNCTION SUDSrf(UDSrfType,Kc,Kr,psi,pv,prel,ihat,
     &	  coef,param,rad,pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

C  SUDSrf provides a shell for user-defined surface function routines.
C  Which routine is invoked depends on variable UDSrfType. Current
C  offerings are:
C
C	1: 1313-actuator deformable mirror.
C	2: 349-actuator deformable mirror.
C	3: 397-actuator triangular array deformable mirror.
C	4: 349-actuator deformable mirror with squared-out
C	   influence function.
C	5: 349-actuator deformable mirror with ideal actuator
C	   influence functions (bilinear approximation).
C	6: 349-actuator deformable mirror with squared-out
C	   influence function.
C	7: 149-by-149 x-y grid for Hubble Primary Mirror.
C	8: 115-by-115 x-y grid for Hubble Secondary Mirror.
C	9: Radial-actuator DM.
C      10: SQAna-by-SQAna DM with damped-sine influence fcns.
C      12: 185-by-185 x-y grid

C  Note that beta*dAct=1.79519221165970

	IMPLICIT NONE
	INTEGER UDSrfType,SQAna
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),rad,
     &	  pmon(3),rhoc(3),xhat(3),yhat(3),zhat(3),L,param(6),
     &	  SUDSrf,SDM349,SDM355,SDM349I,SDM349B,SHSTPM,SHSTSM,SDMRAD,
     &	  SGridSrf,SDM1313,SSQARRAY

	IF (UDSrfType.EQ.1) THEN
c	  dAct=3.5d0
c	  beta=0.5129120604742d0
	  SUDSrf=SDM1313(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                 coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.2) THEN
c	  dAct=0.0036d0
c	  beta=4.986645032d2
C	  dAct=0.00375d0
C	  beta=4.7872e+02
	  SUDSrf=SDM349(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.3) THEN
C	  beta=2.5645654315019d0 
C	  dAct=0.7d0
	  SUDSrf=SDM355(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.4) THEN
C	  beta=2.5645654315019d0 
C	  dAct=0.7d0
	  SUDSrf=SDM349B(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                 coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.5) THEN
c	  dAct=0.0036d0
c	  beta=4.986645032d2
C	  dAct=0.00375d0
C	  beta=4.7872e+02
	  SUDSrf=SDM349I(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                 coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.6) THEN
C	  beta=6.51399619601475d0
C	  dAct=2.7559d-01
	  SUDSrf=SDM349B(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                 coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.7) THEN
C	  beta=0d0
C	  dAct=1.683d0
	  SUDSrf=SHSTPM(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	  	        coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.8) THEN
C	  beta=0d0
C	  dAct=0.302d0
	  SUDSrf=SHSTSM(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	  	        coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.9) THEN
C	  drad=0.85666666666667d0
C	  dq=6.28318530717959d0/16
C	  minrad=15.87166666666667d0
	  SUDSrf=SDMRAD(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                param(3),coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.10) THEN
C	  beta=2.5645654315019d0 
C	  dAct=0.7d0
C	  UDSna=63d0
	  SUDSrf=SSQARRAY(Kc,Kr,psi,pv,prel,ihat,param(1),
     &	  param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	ELSE IF (UDSrfType.EQ.11) THEN
C	  beta=0d0
C	  dAct=0.57432432432432d0
	  SUDSrf=SGridSrf(Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	                  coef,rad,pmon,rhoc,xhat,yhat,zhat,L)
	END IF

	RETURN
	END  FUNCTION SUDSrf

C***********************************************************************

	SUBROUTINE UDSrfNorm(ifBuild,UDSrfType,Kc,Kr,psi,pv,prel,ihat,
     &	 coef,param,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,dFdrho,d2Fdrho2)
	use math_mod

C  UDSrfNorm provides a shell for user-defined surface normal routines.
C  Which routine is invoked depends on variable UDSrfType. Current
C  offerings are:
C
C	1: 1313-actuator deformable mirror.
C	2: 349-actuator deformable mirror.
C	3: 397-actuator triangular array deformable mirror.
C	4: 349-actuator deformable mirror with squared-out
C	   influence function.
C	5: 349-actuator deformable mirror with ideal actuator
C	   influence functions (bilinear approximation).
C	6: 349-actuator deformable mirror with squared-out
C	   influence function.
C	7: 149-by-149 x-y grid for Hubble Primary Mirror.
C	8: 115-by-115 x-y grid for Hubble Secondary Mirror.
C	9: Radial-actuator DM.
C      12: 185-by-185 x-y grid

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER UDSrfType
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),rad,
     &	pmon(3),rhoc(3),xhat(3),yhat(3),zhat(3),L,zc(3),fh,
     &	dFdrho(3),d2Fdrho2(3,3),param(6)

C	Param order = beta,dact or drad,dq,minrad

	IF (UDSrfType.EQ.1) THEN
	  CALL NDM1313(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.2) THEN
	  CALL NDM349(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.3) THEN
	  CALL NDM355(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.4) THEN
	  CALL NDM349B(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.5) THEN
	  CALL NDM349I(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.6) THEN
	  CALL NDM349B(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.7) THEN
	  CALL NHSTPM(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.8) THEN
	  CALL NHSTSM(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.9) THEN
	  CALL NDMRAD(ifBuild,Kc,Kr,psi,pv,prel,ihat,param(1),param(2),
     &	    param(3),coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	    dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.10) THEN
	  CALL NSQARRAY(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	ELSE IF (UDSrfType.EQ.11) THEN
	  CALL NGridSrf(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	    param(1),param(2),coef,rad,pmon,rhoc,xhat,yhat,zhat,
     &	    L,zc,fh,dFdrho,d2Fdrho2)
	END IF

	RETURN
	END SUBROUTINE UDSrfNorm

C***********************************************************************

	FUNCTION IGridSrf(UDSFile)
	use math_mod

	IMPLICIT NONE

	INTEGER ArchFlag,Sun,Mac,DEC,RecordLength
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	INTEGER UDSna,UDSctr
	CHARACTER(len=*) :: UDSFile
	PARAMETER (UDSna=185, UDSctr=93)
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm,inVec(UDSna)
	COMMON /GridSrf1/ UDSmat,UDSnorm
	SAVE /GridSrf1/

	LOGICAL exist,IGridSrf
	INTEGER i,j,iDat,ICLEN

  103	FORMAT(' Gridded surface array is of dimension',i4,' by',i4)

	UDSnorm=0.6328d-4

C   Set record length depending on target machine
	IF (ArchFlag.EQ.Sun) THEN
	  RecordLength=8*UDSna
	ELSE IF (ArchFlag.EQ.Mac) THEN
	  RecordLength=UDSna
	ELSE IF (ArchFlag.EQ.DEC) THEN
	  RecordLength=UDSna
	END IF

C   Check file existence
	IF (UDSFile.EQ.' ') UDSFile='GridSrf.txt'
	INQUIRE (FILE=UDSfile,EXIST=exist)
	IF (exist) THEN
          OPEN (UNIT=94,FILE=UDSfile,STATUS='OLD',FORM='FORMATTED')
	ELSE
	  WRITE(*,*)' File ',UDSFile(1:ICLEN(UDSFile)),
     &	  ' does not exist (func IGridSrf in usersub).'
	  IGridSrf=.FALSE.
	  RETURN
	END IF

C  Read input data
        iDat=0
	DO j=1,UDSna
	  iDat=iDat+1
	  READ(94,*,ERR=2,END=1) (UDSmat(i,j),i=1,UDSna)
	END DO

  1	CONTINUE
	IGridSrf=.TRUE.
	WRITE(*,103) UDSna,iDat
	CLOSE(94)
	RETURN

  2	CONTINUE
	IGridSrf=.FALSE.
	WRITE(*,*) ' Error during read: IGridSrf failed!!!'
	CLOSE(94)
	RETURN

	END FUNCTION IGridSrf  

C***********************************************************************

	FUNCTION SGridSrf(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

C  Computes surface function for gridded data surface

	IMPLICIT NONE

	INTEGER UDSna,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=185, UDSctr=93)
	PARAMETER (udsfile='GridSrf.bin')
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm
	COMMON /GridSrf1/ UDSmat,UDSnorm
	SAVE /GridSrf1/

	INTEGER i,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	 yhat(3),zhat(3),rad,pmon(3),L,SGridSrf,rho(3),rhom(3),rhoc(3),
     &	 fh,beta,dAct,x,y,zc(3),zc2,psitzc,tiltx,tilty,xi,yj,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
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

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSna)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
	END IF

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SGridSrf=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SGridSrf

C***********************************************************************

	SUBROUTINE NGridSrf(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSctr
	PARAMETER (UDSna=185, UDSctr=93)
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm
	COMMON /GridSrf1/ UDSmat,UDSnorm
	SAVE /GridSrf1/

	LOGICAL ifBuild
	INTEGER i,iAxis,jAxis,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),
     &	fh,x,y,zc(3),beta,dAct,xi,yj,
     &	dFdrho(3),d2Fdrho2(3,3),tiltx,tilty,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
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

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSna)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
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
	END SUBROUTINE NGridSrf

C***********************************************************************
C***********************************************************************

	FUNCTION IDMRAD(UDSFile)

	IMPLICIT NONE

	INTEGER ArchFlag,Sun,Mac,DEC,RecordLength
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	INTEGER UDSna,UDSnb,UDSctr
	CHARACTER(len=*) :: UDSFile
	PARAMETER (UDSna=31, UDSnb=16, UDSctr=0)
	REAL*8 UDSmat(UDSna,UDSnb),UDSnorm,inVec(UDSna)
	COMMON /UDSdata3/ UDSmat,UDSnorm
	SAVE /UDSdata3/

	LOGICAL exist,IDMRAD
	INTEGER i,j,iDat,ICLEN

  103	FORMAT(' DMRAD array is of dimension',i4,' by',i4)

	UDSnorm=0.6328d-4

C   Check file existence
	IF (UDSFile.EQ.' ') UDSFile='dmrad.waffle.txt'
	INQUIRE (FILE=UDSfile,EXIST=exist)
	IF (exist) THEN
	  OPEN (UNIT=94,FILE=UDSfile,FORM='FORMATTED',STATUS='OLD')
	ELSE
	  WRITE(*,*)' File ',UDSFile(1:ICLEN(UDSFile)),
     &	  ' does not exist (func IDMRAD in usersub.f).'
	  IDMRAD=.FALSE.
	  RETURN
	END IF

C  Read input data
        iDat=0
	DO i=1,UDSna
	  iDat=iDat+1
	  READ(94,*,ERR=2,END=1) (UDSmat(i,j),j=1,UDSnb)
	END DO

  1	CONTINUE
	IDMRAD=.TRUE.
	WRITE(*,103) iDat,UDSnb
	CLOSE(94)
	RETURN

  2	IDMRAD=.FALSE.
	WRITE(*,*) ' Error during read: IDMRAD failed!!!'
	RETURN

	END FUNCTION IDMRAD

C***********************************************************************
C***********************************************************************

	FUNCTION IHSTPM(UDSFile)

	IMPLICIT NONE

	INTEGER ArchFlag,Sun,Mac,DEC,RecordLength
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	INTEGER UDSna,UDSctr
	CHARACTER(len=*) :: UDSFile
	PARAMETER (UDSna=149, UDSctr=75)
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm,inVec(UDSna)
	COMMON /UDSdata1/ UDSmat,UDSnorm
	SAVE /UDSdata1/

	LOGICAL exist,IHSTPM
	INTEGER i,j,iDat,ICLEN

  103	FORMAT(' HST PM array is of dimension',i4,' by',i4)

	UDSnorm=0.6328d-4

C   Set record length depending on target machine
	IF (ArchFlag.EQ.Sun) THEN
	  RecordLength=8*UDSna
	ELSE IF (ArchFlag.EQ.Mac) THEN
	  RecordLength=UDSna
	ELSE IF (ArchFlag.EQ.DEC) THEN
	  RecordLength=UDSna
	END IF

C   Check file existence
	IF (UDSFile.EQ.' ') UDSFile='hstpm.bin'
	INQUIRE (FILE=UDSfile,EXIST=exist)
	IF (exist) THEN
	  OPEN (UNIT=94,FILE=UDSfile,ACCESS='DIRECT',
     &	  RECL=RecordLength,STATUS='OLD')
	ELSE
	  WRITE(*,*)' File ',UDSFile(1:ICLEN(UDSFile)),
     &	  ' does not exist (func IHSTPM in usersub.f).'
	  IHSTPM=.FALSE.
	  RETURN
	END IF

C  Read input data
        iDat=0
	DO j=1,UDSna
	  iDat=iDat+1
c	  READ(94,REC=j,ERR=2,END=1) (UDSmat(i,j),i=1,UDSna)
	  READ(94,REC=j,ERR=2) (UDSmat(i,j),i=1,UDSna)
	END DO

  1	CONTINUE
	IHSTPM=.TRUE.
	WRITE(*,103) UDSna,iDat
	CLOSE(94)
	RETURN

  2	CONTINUE
	IHSTPM=.FALSE.
	WRITE(*,*) ' Error during read: IHSTPM failed!!!'
	RETURN

	END FUNCTION IHSTPM

C***********************************************************************
C***********************************************************************

	FUNCTION IHSTSM(UDSFile)

	IMPLICIT NONE

	INTEGER ArchFlag,Sun,Mac,DEC,RecordLength
	PARAMETER(Sun=0,Mac=1,DEC=2,ArchFlag=Sun)

	INTEGER UDSna,UDSctr
	CHARACTER(len=*) :: UDSFile
	PARAMETER (UDSna=115, UDSctr=58)
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm,inVec(UDSna)
	COMMON /UDSdata2/ UDSmat,UDSnorm
	SAVE /UDSdata2/

	LOGICAL exist,IHSTSM
	INTEGER i,j,iDat,ICLEN

  103	FORMAT(' HST SM array is of dimension',i4,' by',i4)

	UDSnorm=0.6328d-4

C   Set record length depending on target machine
	IF (ArchFlag.EQ.Sun) THEN
	  RecordLength=8*UDSna
	ELSE IF (ArchFlag.EQ.Mac) THEN
	  RecordLength=UDSna
	ELSE IF (ArchFlag.EQ.DEC) THEN
	  RecordLength=UDSna
	END IF

C   Chec file existence
	IF (UDSFile.EQ.' ') UDSFile='hstsm.bin'
	INQUIRE (FILE=UDSfile,EXIST=exist)
	IF (exist) THEN
	  OPEN (UNIT=94,FILE=UDSfile,ACCESS='DIRECT',
     &	  RECL=RecordLength,STATUS='OLD')
	ELSE
	  WRITE(*,*)' File ',UDSFile(1:ICLEN(UDSFile)),
     &	  ' does not exist (func IHSTSM in usersub.f).'
	  IHSTSM=.FALSE.
	  RETURN
	END IF

C  Read input data
        iDat=0
	DO j=1,UDSna
	  iDat=iDat+1
c	  READ(94,REC=j,ERR=2,END=1) (UDSmat(i,j),i=1,UDSna)
	  READ(94,REC=j,ERR=2) (UDSmat(i,j),i=1,UDSna)
	END DO

  1	CONTINUE
	IHSTSM=.TRUE.
	WRITE(*,103) UDSna,iDat
	CLOSE(94)
	RETURN

  2	CONTINUE
	IHSTSM=.FALSE.
	WRITE(*,*) ' Error during read: IHSTSM failed!!!'
	RETURN

	END FUNCTION IHSTSM

C***********************************************************************
C***********************************************************************

	FUNCTION SDM349(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	pmon,rhoc,xhat,yhat,zhat,L)

C  Computes surface function for 349 actuator, square array DM,
C  units are inches

	IMPLICIT NONE
	INTEGER i,j,actMap(-10:10,-10:10),actIndex(2,349),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SDM349,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),zc2,psitzc,czc

	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/10/
	DATA pio4/0.78539816339745d0/
	DATA (actMap(-10,i),i=-10,10)/0,0,0,0,0,0,0,1,2,3,4,5,6,
     &	7,0,0,0,0,0,0,0/
	DATA (actMap(-9,i),i=-10,10)/0,0,0,0,0,8,9,10,11,12,13,14,15,16,
     &	17,18,0,0,0,0,0/
	DATA (actMap(-8,i),i=-10,10)/0,0,0,0,19,20,21,22,23,24,25,26,27,
     &	28,29,30,31,0,0,0,0/
	DATA (actMap(-7,i),i=-10,10)/0,0,0,32,33,34,35,36,37,38,39,40,
     &	41,42,43,44,45,46,0,0,0/
	DATA (actMap(-6,i),i=-10,10)/0,0,47,48,49,50,51,52,53,54,55,
     &	56,57,58,59,60,61,62,63,0,0/
	DATA (actMap(-5,i),i=-10,10)/0,64,65,66,67,68,69,70,71,72,73,
     &	74,75,76,77,78,79,80,81,82,0/
	DATA (actMap(-4,i),i=-10,10)/0,83,84,85,86,87,88,89,90,91,92,
     &	93,94,95,96,97,98,99,100,101,0/
	DATA (actMap(-3,i),i=-10,10)/102,103,104,105,106,107,108,109,
     &	110,111,112,113,114,115,116,117,118,119,120,121,122/
	DATA (actMap(-2,i),i=-10,10)/123,124,125,126,127,128,129,130,
     &	131,132,133,134,135,136,137,138,139,140,141,142,143/
	DATA (actMap(-1,i),i=-10,10)/144,145,146,147,148,149,150,151,
     &	152,153,154,155,156,157,158,159,160,161,162,163,164/
	DATA (actMap(0,i),i=-10,10)/165,166,167,168,169,170,171,172,
     &	173,174,175,176,177,178,179,180,181,182,183,184,185/
	DATA (actMap(1,i),i=-10,10)/186,187,188,189,190,191,192,193,
     &	194,195,196,197,198,199,200,201,202,203,204,205,206/
	DATA (actMap(2,i),i=-10,10)/207,208,209,210,211,212,213,214,
     &	215,216,217,218,219,220,221,222,223,224,225,226,227/
	DATA (actMap(3,i),i=-10,10)/228,229,230,231,232,233,234,235,
     &	236,237,238,239,240,241,242,243,244,245,246,247,248/
	DATA (actMap(4,i),i=-10,10)/0,249,250,251,252,253,254,255,
     &	256,257,258,259,260,261,262,263,264,265,266,267,0/
	DATA (actMap(5,i),i=-10,10)/0,268,269,270,271,272,273,274,
     &	275,276,277,278,279,280,281,282,283,284,285,286,0/
	DATA (actMap(6,i),i=-10,10)/0,0,287,288,289,290,291,292,
     &	293,294,295,296,297,298,299,300,301,302,303,0,0/
	DATA (actMap(7,i),i=-10,10)/0,0,0,304,305,306,307,308,309,
     &	310,311,312,313,314,315,316,317,318,0,0,0/
	DATA (actMap(8,i),i=-10,10)/0,0,0,0,319,320,321,322,323,
     &	324,325,326,327,328,329,330,331,0,0,0,0/
	DATA (actMap(9,i),i=-10,10)/0,0,0,0,0,332,333,334,335,336,
     &	337,338,339,340,341,342,0,0,0,0,0/
	DATA (actMap(10,i),i=-10,10)/0,0,0,0,0,0,0,343,344,345,346,
     &	347,348,349,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,349)/-3,-2,-1, 0, 1, 2, 3,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-6,-5,-4,-3,-2,-1, 0,
     &	1, 2, 3, 4, 5, 6,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-10,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,
     &	-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-7,-6,-5,-4,-3,-2,-1, 0,
     &	 1, 2, 3, 4, 5, 6, 7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-3,-2,-1, 0, 1, 2, 3/
	DATA (actIndex(1,i),i=1,349)/-10,-10,
     &	-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-8,-8,-8,
     &	-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,
     &	-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,
     &	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-3,-3,-3,-3,
     &	-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
     &	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 
     &	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
     &	 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
     &	 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
     &	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     &	 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
     &	 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     &	 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     &	 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10/ 

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add influence functions of neighboring actuators

	imin=MAX((iRho-idLim),-nzc)
	imax=MIN((iRho+idLim), nzc)
	jmin=MAX((jRho-idLim),-nzc)
	jmax=MIN((jRho+idLim), nzc)
	fh=0d0
	DO j=jmin,jmax
   	  DO i=imin,imax
   	    nearAct=actMap(i,j)
   	    IF (nearAct.GT.0) THEN
   	      IF (coef(nearAct).NE.0d0) THEN
   	        xAct=dAct*actIndex(1,nearAct)
   	        yAct=dAct*actIndex(2,nearAct)
   	        bd=beta*SQRT((x-xAct)**2+(y-yAct)**2)
   	        IF (bd.LT.czc) 
     &	          fh=fh+coef(nearAct)*EXP(-bd)*SIN(bd+pio4)
   	      END IF
   	    END IF
   	  END DO
   	END DO

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SDM349=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SDM349

C***********************************************************************
C***********************************************************************

	SUBROUTINE NDM349(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	dFdrho,d2Fdrho2)

C  Computes surface function for 349 actuator, square array DM

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER i,j,actMap(-10:10,-10:10),actIndex(2,349),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,iAxis,jAxis,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),rhoAx,
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),
     &	dFdrho(3),d2Fdrho2(3,3),esbd,ecbd,ebd,beta2,d,bod,dVec(3),
     &	CdF,CdF2,dx,dy,czc

	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/10/
	DATA pio4/0.78539816339745d0/
	DATA (actMap(-10,i),i=-10,10)/0,0,0,0,0,0,0,1,2,3,4,5,6,
     &	7,0,0,0,0,0,0,0/
	DATA (actMap(-9,i),i=-10,10)/0,0,0,0,0,8,9,10,11,12,13,14,15,16,
     &	17,18,0,0,0,0,0/
	DATA (actMap(-8,i),i=-10,10)/0,0,0,0,19,20,21,22,23,24,25,26,27,
     &	28,29,30,31,0,0,0,0/
	DATA (actMap(-7,i),i=-10,10)/0,0,0,32,33,34,35,36,37,38,39,40,
     &	41,42,43,44,45,46,0,0,0/
	DATA (actMap(-6,i),i=-10,10)/0,0,47,48,49,50,51,52,53,54,55,
     &	56,57,58,59,60,61,62,63,0,0/
	DATA (actMap(-5,i),i=-10,10)/0,64,65,66,67,68,69,70,71,72,73,
     &	74,75,76,77,78,79,80,81,82,0/
	DATA (actMap(-4,i),i=-10,10)/0,83,84,85,86,87,88,89,90,91,92,
     &	93,94,95,96,97,98,99,100,101,0/
	DATA (actMap(-3,i),i=-10,10)/102,103,104,105,106,107,108,109,
     &	110,111,112,113,114,115,116,117,118,119,120,121,122/
	DATA (actMap(-2,i),i=-10,10)/123,124,125,126,127,128,129,130,
     &	131,132,133,134,135,136,137,138,139,140,141,142,143/
	DATA (actMap(-1,i),i=-10,10)/144,145,146,147,148,149,150,151,
     &	152,153,154,155,156,157,158,159,160,161,162,163,164/
	DATA (actMap(0,i),i=-10,10)/165,166,167,168,169,170,171,172,
     &	173,174,175,176,177,178,179,180,181,182,183,184,185/
	DATA (actMap(1,i),i=-10,10)/186,187,188,189,190,191,192,193,
     &	194,195,196,197,198,199,200,201,202,203,204,205,206/
	DATA (actMap(2,i),i=-10,10)/207,208,209,210,211,212,213,214,
     &	215,216,217,218,219,220,221,222,223,224,225,226,227/
	DATA (actMap(3,i),i=-10,10)/228,229,230,231,232,233,234,235,
     &	236,237,238,239,240,241,242,243,244,245,246,247,248/
	DATA (actMap(4,i),i=-10,10)/0,249,250,251,252,253,254,255,
     &	256,257,258,259,260,261,262,263,264,265,266,267,0/
	DATA (actMap(5,i),i=-10,10)/0,268,269,270,271,272,273,274,
     &	275,276,277,278,279,280,281,282,283,284,285,286,0/
	DATA (actMap(6,i),i=-10,10)/0,0,287,288,289,290,291,292,
     &	293,294,295,296,297,298,299,300,301,302,303,0,0/
	DATA (actMap(7,i),i=-10,10)/0,0,0,304,305,306,307,308,309,
     &	310,311,312,313,314,315,316,317,318,0,0,0/
	DATA (actMap(8,i),i=-10,10)/0,0,0,0,319,320,321,322,323,
     &	324,325,326,327,328,329,330,331,0,0,0,0/
	DATA (actMap(9,i),i=-10,10)/0,0,0,0,0,332,333,334,335,336,
     &	337,338,339,340,341,342,0,0,0,0,0/
	DATA (actMap(10,i),i=-10,10)/0,0,0,0,0,0,0,343,344,345,346,
     &	347,348,349,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,349)/-3,-2,-1, 0, 1, 2, 3,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-6,-5,-4,-3,-2,-1, 0,
     &	1, 2, 3, 4, 5, 6,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-10,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,
     &	-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-7,-6,-5,-4,-3,-2,-1, 0,
     &	 1, 2, 3, 4, 5, 6, 7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-3,-2,-1, 0, 1, 2, 3/
	DATA (actIndex(1,i),i=1,349)/-10,-10,
     &	-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-8,-8,-8,
     &	-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,
     &	-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,
     &	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-3,-3,-3,-3,
     &	-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
     &	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 
     &	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
     &	 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
     &	 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
     &	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     &	 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
     &	 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     &	 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     &	 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10/ 

	beta2=beta*beta

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add neighboring actuator influence functions and derivatives

	imin=MAX((iRho-idLim),-nzc)
	imax=MIN((iRho+idLim), nzc)
	jmin=MAX((jRho-idLim),-nzc)
	jmax=MIN((jRho+idLim), nzc)
	fh=0d0
	IF (ifBuild) THEN
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
  	    DO iAxis=1,3
	      d2Fdrho2(iAxis,jAxis)=0d0
   	    END DO
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          d=SQRT((x-xAct)**2+(y-yAct)**2)
   	          bd=beta*d
   	          IF (bd.LT.czc) THEN
   	            DO iAxis=1,3
   	              rhoAx=xAct*xhat(iAxis)+
     &                yAct*yhat(iAxis)-rhoc(iAxis)
c 	              dVec(iAxis)=rhom(iAxis)-rhoAx
  	              dVec(iAxis)=(1d0-zhat(iAxis))*rhom(iAxis)-rhoAx
   	            END DO
   	            IF (d.lt.1d-12) THEN
   	              ebd=coef(nearAct)
    	              esbd=ebd*0.70710678118655d0
    	              ecbd=ebd*0.70710678118655d0
   	              fh=fh+esbd
     	              CdF2=-2d0*beta*beta*ecbd
 	              DO iAxis=1,3
   	                d2Fdrho2(iAxis,iAxis)=d2Fdrho2(iAxis,iAxis)+CdF2
   	              END DO
    	            ELSE
   	              bod=beta/d
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
   	              fh=fh+esbd
   	              CdF=bod*(ecbd-esbd)
   	              CdF2=-2d0*bod*bod*ecbd
   	              DO jAxis=1,3
   	                dFdrho(jAxis)=dFdrho(jAxis)+CdF*dVec(jAxis)
   	                d2Fdrho2(jAxis,jAxis)=d2Fdrho2(jAxis,jAxis)+CdF
   	                DO iAxis=1,3
   	                  d2Fdrho2(iAxis,jAxis)=d2Fdrho2(iAxis,jAxis)+
     &   	          (CdF2-CdF/d/d)*dVec(iAxis)*dVec(jAxis)
   	                END DO
   	              END DO
   	            END IF
   	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	ELSE
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          dx=x-xAct
   	          dy=y-yAct
   	          d=SQRT(dx*dx+dy*dy)
   	          IF (d.lt.1d-12) THEN
   	            ebd=coef(nearAct)
    	            esbd=ebd*0.70710678118655d0
  	            fh=fh+esbd
   	          ELSE
   	            bd=beta*d
    	            IF (bd.LT.czc) THEN
  	              bod=beta/d
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
  	              fh=fh+esbd
	              CdF=bod*(ecbd-esbd)
 	              DO jAxis=1,3
 	                dFdrho(jAxis)=dFdrho(jAxis)
     & 	                +CdF*(dx*xhat(jAxis)+dy*yhat(jAxis))
   	              END DO
  	            END IF
  	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	END IF

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NDM349

C***********************************************************************
C***********************************************************************

	FUNCTION SDM355(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	pmon,rhoc,xhat,yhat,zhat,L)

C  Computes surface function for 355 actuator, hex array DM

	IMPLICIT NONE
	INTEGER i,j,actMap(-21:21,-11:11),actIndex(2,355),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SDM355,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),zc2,psitzc,
     &	h2,w2,sin60,cos60

	DATA idLim/9/
	DATA pio4/0.78539816339745d0/
	DATA sin60/0.86602540378444d0/
	DATA cos60/0.5d0/
	DATA (actMap(-21,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,1,0,2,0,
     &	0,0,0,0,0,0,0,0,0/
	DATA (actMap(-20,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,0,3,0,0,
     &	0,0,0,0,0,0,0,0,0/
	DATA (actMap(-19,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,4,0,5,0,
     &	0,0,0,0,0,0,0,0,0/
	DATA (actMap(-18,i),i=-11,11)/0,0,0,0,0,0,0,6,0,7,0,8,0,9,
     &	0,10,0,0,0,0,0,0,0/
	DATA (actMap(-17,i),i=-11,11)/0,0,0,0,0,0,0,0,11,0,12,0,13,
     &	0,14,0,0,0,0,0,0,0,0/
	DATA (actMap(-16,i),i=-11,11)/0,0,0,0,0,0,0,15,0,16,0,17,0,
     &	18,0,19,0,0,0,0,0,0,0/
	DATA (actMap(-15,i),i=-11,11)/0,0,0,0,20,0,21,0,22,0,23,0,
     &	24,0,25,0,26,0,27,0,0,0,0/
	DATA (actMap(-14,i),i=-11,11)/0,0,0,0,0,28,0,29,0,30,0,31,
     &	0,32,0,33,0,34,0,0,0,0,0/
	DATA (actMap(-13,i),i=-11,11)/0,0,0,0,35,0,36,0,37,0,38,0,
     &	39,0,40,0,41,0,42,0,0,0,0/
	DATA (actMap(-12,i),i=-11,11)/0,43,0,44,0,45,0,46,0,47,0,
     &	48,0,49,0,50,0,51,0,52,0,53,0/
	DATA (actMap(-11,i),i=-11,11)/0,0,54,0,55,0,56,0,57,0,58,
     &	0,59,0,60,0,61,0,62,0,63,0,0/
	DATA (actMap(-10,i),i=-11,11)/0,64,0,65,0,66,0,67,0,68,0,
     &	69,0,70,0,71,0,72,0,73,0,74,0/
	DATA (actMap(-9,i),i=-11,11)/75,0,76,0,77,0,78,0,79,0,80,
     &	0,81,0,82,0,83,0,84,0,85,0,86/
	DATA (actMap(-8,i),i=-11,11)/0,87,0,88,0,89,0,90,0,91,0,
     &	92,0,93,0,94,0,95,0,96,0,97,0/
	DATA (actMap(-7,i),i=-11,11)/0,0,98,0,99,0,100,0,101,0,
     &	102,0,103,0,104,0,105,0,106,0,107,0,0/
	DATA (actMap(-6,i),i=-11,11)/0,108,0,109,0,110,0,111,0,
     &	112,0,113,0,114,0,115,0,116,0,117,0,118,0/
	DATA (actMap(-5,i),i=-11,11)/0,0,119,0,120,0,121,0,122,0,
     &	123,0,124,0,125,0,126,0,127,0,128,0,0/
	DATA (actMap(-4,i),i=-11,11)/0,129,0,130,0,131,0,132,0,
     &	133,0,134,0,135,0,136,0,137,0,138,0,139,0/
	DATA (actMap(-3,i),i=-11,11)/140,0,141,0,142,0,143,0,144,
     &	0,145,0,146,0,147,0,148,0,149,0,150,0,151/
	DATA (actMap(-2,i),i=-11,11)/0,152,0,153,0,154,0,155,0,
     &	156,0,157,0,158,0,159,0,160,0,161,0,162,0/
	DATA (actMap(-1,i),i=-11,11)/0,0,163,0,164,0,165,0,166,0,
     &	167,0,168,0,169,0,170,0,171,0,172,0,0/
	DATA (actMap(0,i),i=-11,11)/0,173,0,174,0,175,0,176,0,
     &	177,0,178,0,179,0,180,0,181,0,182,0,183,0/
	DATA (actMap(1,i),i=-11,11)/0,0,184,0,185,0,186,0,187,0,
     &	188,0,189,0,190,0,191,0,192,0,193,0,0/
	DATA (actMap(2,i),i=-11,11)/0,194,0,195,0,196,0,197,0,
     &	198,0,199,0,200,0,201,0,202,0,203,0,204,0/
	DATA (actMap(3,i),i=-11,11)/205,0,206,0,207,0,208,0,209,
     &	0,210,0,211,0,212,0,213,0,214,0,215,0,216/
	DATA (actMap(4,i),i=-11,11)/0,217,0,218,0,219,0,220,0,
     &	221,0,222,0,223,0,224,0,225,0,226,0,227,0/
	DATA (actMap(5,i),i=-11,11)/0,0,228,0,229,0,230,0,231,0,
     &	232,0,233,0,234,0,235,0,236,0,237,0,0/
	DATA (actMap(6,i),i=-11,11)/0,238,0,239,0,240,0,241,0,
     &	242,0,243,0,244,0,245,0,246,0,247,0,248,0/
	DATA (actMap(7,i),i=-11,11)/0,0,249,0,250,0,251,0,252,0,
     &	253,0,254,0,255,0,256,0,257,0,258,0,0/
	DATA (actMap(8,i),i=-11,11)/0,259,0,260,0,261,0,262,0,
     &	263,0,264,0,265,0,266,0,267,0,268,0,269,0/
	DATA (actMap(9,i),i=-11,11)/270,0,271,0,272,0,273,0,274,
     &	0,275,0,276,0,277,0,278,0,279,0,280,0,281/
	DATA (actMap(10,i),i=-11,11)/0,282,0,283,0,284,0,285,0,
     &	286,0,287,0,288,0,289,0,290,0,291,0,292,0/
	DATA (actMap(11,i),i=-11,11)/0,0,293,0,294,0,295,0,296,0,
     &	297,0,298,0,299,0,300,0,301,0,302,0,0/
	DATA (actMap(12,i),i=-11,11)/0,303,0,304,0,305,0,306,0,
     &	307,0,308,0,309,0,310,0,311,0,312,0,313,0/
	DATA (actMap(13,i),i=-11,11)/0,0,0,0,314,0,315,0,316,0,
     &	317,0,318,0,319,0,320,0,321,0,0,0,0/
	DATA (actMap(14,i),i=-11,11)/0,0,0,0,0,322,0,323,0,324,
     &	0,325,0,326,0,327,0,328,0,0,0,0,0/
	DATA (actMap(15,i),i=-11,11)/0,0,0,0,329,0,330,0,331,0,
     &	332,0,333,0,334,0,335,0,336,0,0,0,0/
	DATA (actMap(16,i),i=-11,11)/0,0,0,0,0,0,0,337,0,338,0,
     &	339,0,340,0,341,0,0,0,0,0,0,0/
	DATA (actMap(17,i),i=-11,11)/0,0,0,0,0,0,0,0,342,0,343,
     &	0,344,0,345,0,0,0,0,0,0,0,0/
	DATA (actMap(18,i),i=-11,11)/0,0,0,0,0,0,0,346,0,347,0,
     &	348,0,349,0,350,0,0,0,0,0,0,0/
	DATA (actMap(19,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,351,0,
     &	352,0,0,0,0,0,0,0,0,0,0/
	DATA (actMap(20,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,0,353,
     &	0,0,0,0,0,0,0,0,0,0,0/
	DATA (actMap(21,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,354,0,
     &	355,0,0,0,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,355)/-1,1,0,-1,1,-4,-2,0,2,4,-3,-1,1,3,
     &	-4,-2,0,2,4,-7,-5,-3,-1,1,3,5,7,-6,-4,-2,0,2,4,6,-7,-5,-3,
     &	-1,1,3,5,7,-10,-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,
     &	5,7,9,-10,-8,-6,-4,-2,0,2,4,6,8,10,-11,-9,-7,-5,-3,-1,1,3,5,
     &	7,9,11,-10,-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,
     &	-10,-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,
     &	-8,-6,-4,-2,0,2,4,6,8,10,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,-10,
     &	-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,
     &	-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,
     &	0,2,4,6,8,10,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,-10,-8,-6,-4,-2,
     &	0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,0,2,4,
     &	6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,0,2,4,6,8,10,
     &	-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,-10,-8,-6,-4,-2,0,2,4,6,8,10,
     &	-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,0,2,4,6,8,10,-7,-5,
     &	-3,-1,1,3,5,7,-6,-4,-2,0,2,4,6,-7,-5,-3,-1,1,3,5,7,-4,-2,0,2,
     &	4,-3,-1,1,3,-4,-2,0,2,4,-1,1,0,-1,1/
	DATA (actIndex(1,i),i=1,355)/-21,-21,-20,-19,-19,-18,-18,-18,
     &	-18,-18,-17,-17,-17,-17,-16,-16,-16,-16,-16,-15,-15,-15,
     &	-15,-15,-15,-15,-15,-14,-14,-14,-14,-14,-14,-14,-13,-13,
     &	-13,-13,-13,-13,-13,-13,-12,-12,-12,-12,-12,-12,-12,-12,
     &	-12,-12,-12,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-10,
     &	-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,
     &	-9,-9,-9,-9,-9,-9,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,
     &	-7,-7,-7,-7,-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,-4,-4,-4,-4,-4,-4,-4,
     &	-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,
     &	0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,
     &	3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,
     &	5,5,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,
     &	8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,
     &	10,10,10,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,
     &	12,12,12,12,12,13,13,13,13,13,13,13,13,14,14,14,14,14,14,
     &	14,15,15,15,15,15,15,15,15,16,16,16,16,16,17,17,17,17,18,
     &	18,18,18,18,19,19,20,21,21/

	w2=dAct*sin60
	h2=dAct*cos60

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/h2)
	jRho=IDNINT(y/w2)

C  Add influence functions of neighboring actuators

	imin=MAX((iRho-idLim),-21)
	imax=MIN((iRho+idLim), 21)
	jmin=MAX((jRho-idLim),-11)
	jmax=MIN((jRho+idLim), 11)
	fh=0d0
	DO j=jmin,jmax
   	  DO i=imin,imax
   	    nearAct=actMap(i,j)
   	    IF (nearAct.GT.0) THEN
   	      IF (coef(nearAct).NE.0d0) THEN
   	        xAct=h2*actIndex(1,nearAct)
   	        yAct=w2*actIndex(2,nearAct)
   	        bd=beta*SQRT((x-xAct)**2+(y-yAct)**2)
   	        fh=fh+coef(nearAct)*EXP(-bd)*SIN(bd+pio4)
   	      END IF
   	    END IF
   	  END DO
   	END DO

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SDM355=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SDM355

C***********************************************************************
C***********************************************************************

	SUBROUTINE NDM355(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	dFdrho,d2Fdrho2)

C  Computes surface function for 355 actuator, hex array DM

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER i,j,actMap(-21:21,-11:11),actIndex(2,355),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,iAxis,jAxis
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),rhoAx,
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),
     &	dFdrho(3),d2Fdrho2(3,3),esbd,ecbd,ebd,beta2,d,bod,dVec(3),
     &	CdF,CdF2,h2,w2,sin60,cos60,dx,dy

	DATA idLim/9/
	DATA pio4/0.78539816339745d0/
	DATA sin60/0.86602540378444d0/
	DATA cos60/0.5d0/
	DATA (actMap(-21,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,1,0,2,0,
     &	0,0,0,0,0,0,0,0,0/
	DATA (actMap(-20,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,0,3,0,0,
     &	0,0,0,0,0,0,0,0,0/
	DATA (actMap(-19,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,4,0,5,0,
     &	0,0,0,0,0,0,0,0,0/
	DATA (actMap(-18,i),i=-11,11)/0,0,0,0,0,0,0,6,0,7,0,8,0,9,
     &	0,10,0,0,0,0,0,0,0/
	DATA (actMap(-17,i),i=-11,11)/0,0,0,0,0,0,0,0,11,0,12,0,13,
     &	0,14,0,0,0,0,0,0,0,0/
	DATA (actMap(-16,i),i=-11,11)/0,0,0,0,0,0,0,15,0,16,0,17,0,
     &	18,0,19,0,0,0,0,0,0,0/
	DATA (actMap(-15,i),i=-11,11)/0,0,0,0,20,0,21,0,22,0,23,0,
     &	24,0,25,0,26,0,27,0,0,0,0/
	DATA (actMap(-14,i),i=-11,11)/0,0,0,0,0,28,0,29,0,30,0,31,
     &	0,32,0,33,0,34,0,0,0,0,0/
	DATA (actMap(-13,i),i=-11,11)/0,0,0,0,35,0,36,0,37,0,38,0,
     &	39,0,40,0,41,0,42,0,0,0,0/
	DATA (actMap(-12,i),i=-11,11)/0,43,0,44,0,45,0,46,0,47,0,
     &	48,0,49,0,50,0,51,0,52,0,53,0/
	DATA (actMap(-11,i),i=-11,11)/0,0,54,0,55,0,56,0,57,0,58,
     &	0,59,0,60,0,61,0,62,0,63,0,0/
	DATA (actMap(-10,i),i=-11,11)/0,64,0,65,0,66,0,67,0,68,0,
     &	69,0,70,0,71,0,72,0,73,0,74,0/
	DATA (actMap(-9,i),i=-11,11)/75,0,76,0,77,0,78,0,79,0,80,
     &	0,81,0,82,0,83,0,84,0,85,0,86/
	DATA (actMap(-8,i),i=-11,11)/0,87,0,88,0,89,0,90,0,91,0,
     &	92,0,93,0,94,0,95,0,96,0,97,0/
	DATA (actMap(-7,i),i=-11,11)/0,0,98,0,99,0,100,0,101,0,
     &	102,0,103,0,104,0,105,0,106,0,107,0,0/
	DATA (actMap(-6,i),i=-11,11)/0,108,0,109,0,110,0,111,0,
     &	112,0,113,0,114,0,115,0,116,0,117,0,118,0/
	DATA (actMap(-5,i),i=-11,11)/0,0,119,0,120,0,121,0,122,0,
     &	123,0,124,0,125,0,126,0,127,0,128,0,0/
	DATA (actMap(-4,i),i=-11,11)/0,129,0,130,0,131,0,132,0,
     &	133,0,134,0,135,0,136,0,137,0,138,0,139,0/
	DATA (actMap(-3,i),i=-11,11)/140,0,141,0,142,0,143,0,144,
     &	0,145,0,146,0,147,0,148,0,149,0,150,0,151/
	DATA (actMap(-2,i),i=-11,11)/0,152,0,153,0,154,0,155,0,
     &	156,0,157,0,158,0,159,0,160,0,161,0,162,0/
	DATA (actMap(-1,i),i=-11,11)/0,0,163,0,164,0,165,0,166,0,
     &	167,0,168,0,169,0,170,0,171,0,172,0,0/
	DATA (actMap(0,i),i=-11,11)/0,173,0,174,0,175,0,176,0,
     &	177,0,178,0,179,0,180,0,181,0,182,0,183,0/
	DATA (actMap(1,i),i=-11,11)/0,0,184,0,185,0,186,0,187,0,
     &	188,0,189,0,190,0,191,0,192,0,193,0,0/
	DATA (actMap(2,i),i=-11,11)/0,194,0,195,0,196,0,197,0,
     &	198,0,199,0,200,0,201,0,202,0,203,0,204,0/
	DATA (actMap(3,i),i=-11,11)/205,0,206,0,207,0,208,0,209,
     &	0,210,0,211,0,212,0,213,0,214,0,215,0,216/
	DATA (actMap(4,i),i=-11,11)/0,217,0,218,0,219,0,220,0,
     &	221,0,222,0,223,0,224,0,225,0,226,0,227,0/
	DATA (actMap(5,i),i=-11,11)/0,0,228,0,229,0,230,0,231,0,
     &	232,0,233,0,234,0,235,0,236,0,237,0,0/
	DATA (actMap(6,i),i=-11,11)/0,238,0,239,0,240,0,241,0,
     &	242,0,243,0,244,0,245,0,246,0,247,0,248,0/
	DATA (actMap(7,i),i=-11,11)/0,0,249,0,250,0,251,0,252,0,
     &	253,0,254,0,255,0,256,0,257,0,258,0,0/
	DATA (actMap(8,i),i=-11,11)/0,259,0,260,0,261,0,262,0,
     &	263,0,264,0,265,0,266,0,267,0,268,0,269,0/
	DATA (actMap(9,i),i=-11,11)/270,0,271,0,272,0,273,0,274,
     &	0,275,0,276,0,277,0,278,0,279,0,280,0,281/
	DATA (actMap(10,i),i=-11,11)/0,282,0,283,0,284,0,285,0,
     &	286,0,287,0,288,0,289,0,290,0,291,0,292,0/
	DATA (actMap(11,i),i=-11,11)/0,0,293,0,294,0,295,0,296,0,
     &	297,0,298,0,299,0,300,0,301,0,302,0,0/
	DATA (actMap(12,i),i=-11,11)/0,303,0,304,0,305,0,306,0,
     &	307,0,308,0,309,0,310,0,311,0,312,0,313,0/
	DATA (actMap(13,i),i=-11,11)/0,0,0,0,314,0,315,0,316,0,
     &	317,0,318,0,319,0,320,0,321,0,0,0,0/
	DATA (actMap(14,i),i=-11,11)/0,0,0,0,0,322,0,323,0,324,
     &	0,325,0,326,0,327,0,328,0,0,0,0,0/
	DATA (actMap(15,i),i=-11,11)/0,0,0,0,329,0,330,0,331,0,
     &	332,0,333,0,334,0,335,0,336,0,0,0,0/
	DATA (actMap(16,i),i=-11,11)/0,0,0,0,0,0,0,337,0,338,0,
     &	339,0,340,0,341,0,0,0,0,0,0,0/
	DATA (actMap(17,i),i=-11,11)/0,0,0,0,0,0,0,0,342,0,343,
     &	0,344,0,345,0,0,0,0,0,0,0,0/
	DATA (actMap(18,i),i=-11,11)/0,0,0,0,0,0,0,346,0,347,0,
     &	348,0,349,0,350,0,0,0,0,0,0,0/
	DATA (actMap(19,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,351,0,
     &	352,0,0,0,0,0,0,0,0,0,0/
	DATA (actMap(20,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,0,353,
     &	0,0,0,0,0,0,0,0,0,0,0/
	DATA (actMap(21,i),i=-11,11)/0,0,0,0,0,0,0,0,0,0,354,0,
     &	355,0,0,0,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,355)/-1,1,0,-1,1,-4,-2,0,2,4,-3,-1,1,3,
     &	-4,-2,0,2,4,-7,-5,-3,-1,1,3,5,7,-6,-4,-2,0,2,4,6,-7,-5,-3,
     &	-1,1,3,5,7,-10,-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,
     &	5,7,9,-10,-8,-6,-4,-2,0,2,4,6,8,10,-11,-9,-7,-5,-3,-1,1,3,5,
     &	7,9,11,-10,-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,
     &	-10,-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,
     &	-8,-6,-4,-2,0,2,4,6,8,10,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,-10,
     &	-8,-6,-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,
     &	-4,-2,0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,
     &	0,2,4,6,8,10,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,-10,-8,-6,-4,-2,
     &	0,2,4,6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,0,2,4,
     &	6,8,10,-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,0,2,4,6,8,10,
     &	-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,-10,-8,-6,-4,-2,0,2,4,6,8,10,
     &	-9,-7,-5,-3,-1,1,3,5,7,9,-10,-8,-6,-4,-2,0,2,4,6,8,10,-7,-5,
     &	-3,-1,1,3,5,7,-6,-4,-2,0,2,4,6,-7,-5,-3,-1,1,3,5,7,-4,-2,0,2,
     &	4,-3,-1,1,3,-4,-2,0,2,4,-1,1,0,-1,1/
	DATA (actIndex(1,i),i=1,355)/-21,-21,-20,-19,-19,-18,-18,-18,
     &	-18,-18,-17,-17,-17,-17,-16,-16,-16,-16,-16,-15,-15,-15,
     &	-15,-15,-15,-15,-15,-14,-14,-14,-14,-14,-14,-14,-13,-13,
     &	-13,-13,-13,-13,-13,-13,-12,-12,-12,-12,-12,-12,-12,-12,
     &	-12,-12,-12,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-10,
     &	-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,
     &	-9,-9,-9,-9,-9,-9,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,
     &	-7,-7,-7,-7,-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,-4,-4,-4,-4,-4,-4,-4,
     &	-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,
     &	0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,
     &	3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,
     &	5,5,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,
     &	8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,
     &	10,10,10,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,
     &	12,12,12,12,12,13,13,13,13,13,13,13,13,14,14,14,14,14,14,
     &	14,15,15,15,15,15,15,15,15,16,16,16,16,16,17,17,17,17,18,
     &	18,18,18,18,19,19,20,21,21/

	beta2=beta*beta
	w2=dAct*sin60
	h2=dAct*cos60

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/h2)
	jRho=IDNINT(y/w2)

C  Add neighboring actuator influence functions and derivatives

	imin=MAX((iRho-idLim),-22)
	imax=MIN((iRho+idLim), 22)
	jmin=MAX((jRho-idLim),-11)
	jmax=MIN((jRho+idLim), 11)
	fh=0d0
	IF (ifBuild) THEN
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
  	    DO iAxis=1,3
	      d2Fdrho2(iAxis,jAxis)=0d0
   	    END DO
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=h2*actIndex(1,nearAct)
   	          yAct=w2*actIndex(2,nearAct)
   	          DO iAxis=1,3
   	            rhoAx=xAct*xhat(iAxis)+yAct*yhat(iAxis)-rhoc(iAxis)
c 	            dVec(iAxis)=rhom(iAxis)-rhoAx
  	            dVec(iAxis)=(1d0-zhat(iAxis))*rhom(iAxis)-rhoAx
   	          END DO
   	          d=SQRT((x-xAct)**2+(y-yAct)**2)
   	          bd=beta*d
   	          IF (d.lt.1d-12) THEN
   	            ebd=coef(nearAct)
    	            esbd=ebd*0.70710678118655d0
    	            ecbd=ebd*0.70710678118655d0
   	            fh=fh+esbd
     	            CdF2=-2d0*beta*beta*ecbd
 	            DO iAxis=1,3
   	              d2Fdrho2(iAxis,iAxis)=d2Fdrho2(iAxis,iAxis)+CdF2
   	            END DO
    	          ELSE
   	            bod=beta/d
   	            ebd=coef(nearAct)*EXP(-bd)
    	            esbd=ebd*SIN(bd+pio4)
    	            ecbd=ebd*COS(bd+pio4)
   	            fh=fh+esbd
   	            CdF=bod*(ecbd-esbd)
   	            CdF2=-2d0*bod*bod*ecbd
   	            DO jAxis=1,3
   	              dFdrho(jAxis)=dFdrho(jAxis)+CdF*dVec(jAxis)
   	              d2Fdrho2(jAxis,jAxis)=d2Fdrho2(jAxis,jAxis)+CdF
   	              DO iAxis=1,3
   	                d2Fdrho2(iAxis,jAxis)=d2Fdrho2(iAxis,jAxis)+
     &   	        (CdF2-CdF/d/d)*dVec(iAxis)*dVec(jAxis)
   	              END DO
   	            END DO
   	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	ELSE
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=h2*actIndex(1,nearAct)
   	          yAct=w2*actIndex(2,nearAct)
   	          dx=x-xAct
   	          dy=y-yAct
   	          d=SQRT(dx*dx+dy*dy)
   	          IF (d.lt.1d-12) THEN
   	            ebd=coef(nearAct)
    	            esbd=ebd*0.70710678118655d0
  	            fh=fh+esbd
   	          ELSE
   	            bd=beta*d
   	            bod=beta/d
   	            ebd=coef(nearAct)*EXP(-bd)
    	            esbd=ebd*SIN(bd+pio4)
    	            ecbd=ebd*COS(bd+pio4)
  	            fh=fh+esbd
	            CdF=bod*(ecbd-esbd)
 	            DO jAxis=1,3
 	              dFdrho(jAxis)=dFdrho(jAxis)
     & 	              +CdF*(dx*xhat(jAxis)+dy*yhat(jAxis))
   	            END DO
  	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	END IF

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NDM355

C***********************************************************************
C***********************************************************************

	FUNCTION SDM349B(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	pmon,rhoc,xhat,yhat,zhat,L)

C  Computes surface function for 349 actuator, square array DM,
C  units are inches

	IMPLICIT NONE
	INTEGER i,j,actMap(-10:10,-10:10),actIndex(2,349),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,shapeOrder,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SDM349B,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),zc2,psitzc,
     &	shapeOrdInv,dn,dp,dx,dx2,dx4,dx8,dy,dy2,dy4,dy8,czc

	DATA shapeOrder/8/
	DATA shapeOrdInv/0.125d0/
	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/10/
	DATA pio4/0.78539816339745d0/
	DATA (actMap(-10,i),i=-10,10)/0,0,0,0,0,0,0,1,2,3,4,5,6,
     &	7,0,0,0,0,0,0,0/
	DATA (actMap(-9,i),i=-10,10)/0,0,0,0,0,8,9,10,11,12,13,14,15,16,
     &	17,18,0,0,0,0,0/
	DATA (actMap(-8,i),i=-10,10)/0,0,0,0,19,20,21,22,23,24,25,26,27,
     &	28,29,30,31,0,0,0,0/
	DATA (actMap(-7,i),i=-10,10)/0,0,0,32,33,34,35,36,37,38,39,40,
     &	41,42,43,44,45,46,0,0,0/
	DATA (actMap(-6,i),i=-10,10)/0,0,47,48,49,50,51,52,53,54,55,
     &	56,57,58,59,60,61,62,63,0,0/
	DATA (actMap(-5,i),i=-10,10)/0,64,65,66,67,68,69,70,71,72,73,
     &	74,75,76,77,78,79,80,81,82,0/
	DATA (actMap(-4,i),i=-10,10)/0,83,84,85,86,87,88,89,90,91,92,
     &	93,94,95,96,97,98,99,100,101,0/
	DATA (actMap(-3,i),i=-10,10)/102,103,104,105,106,107,108,109,
     &	110,111,112,113,114,115,116,117,118,119,120,121,122/
	DATA (actMap(-2,i),i=-10,10)/123,124,125,126,127,128,129,130,
     &	131,132,133,134,135,136,137,138,139,140,141,142,143/
	DATA (actMap(-1,i),i=-10,10)/144,145,146,147,148,149,150,151,
     &	152,153,154,155,156,157,158,159,160,161,162,163,164/
	DATA (actMap(0,i),i=-10,10)/165,166,167,168,169,170,171,172,
     &	173,174,175,176,177,178,179,180,181,182,183,184,185/
	DATA (actMap(1,i),i=-10,10)/186,187,188,189,190,191,192,193,
     &	194,195,196,197,198,199,200,201,202,203,204,205,206/
	DATA (actMap(2,i),i=-10,10)/207,208,209,210,211,212,213,214,
     &	215,216,217,218,219,220,221,222,223,224,225,226,227/
	DATA (actMap(3,i),i=-10,10)/228,229,230,231,232,233,234,235,
     &	236,237,238,239,240,241,242,243,244,245,246,247,248/
	DATA (actMap(4,i),i=-10,10)/0,249,250,251,252,253,254,255,
     &	256,257,258,259,260,261,262,263,264,265,266,267,0/
	DATA (actMap(5,i),i=-10,10)/0,268,269,270,271,272,273,274,
     &	275,276,277,278,279,280,281,282,283,284,285,286,0/
	DATA (actMap(6,i),i=-10,10)/0,0,287,288,289,290,291,292,
     &	293,294,295,296,297,298,299,300,301,302,303,0,0/
	DATA (actMap(7,i),i=-10,10)/0,0,0,304,305,306,307,308,309,
     &	310,311,312,313,314,315,316,317,318,0,0,0/
	DATA (actMap(8,i),i=-10,10)/0,0,0,0,319,320,321,322,323,
     &	324,325,326,327,328,329,330,331,0,0,0,0/
	DATA (actMap(9,i),i=-10,10)/0,0,0,0,0,332,333,334,335,336,
     &	337,338,339,340,341,342,0,0,0,0,0/
	DATA (actMap(10,i),i=-10,10)/0,0,0,0,0,0,0,343,344,345,346,
     &	347,348,349,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,349)/-3,-2,-1, 0, 1, 2, 3,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-6,-5,-4,-3,-2,-1, 0,
     &	1, 2, 3, 4, 5, 6,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-10,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,
     &	-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-7,-6,-5,-4,-3,-2,-1, 0,
     &	 1, 2, 3, 4, 5, 6, 7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-3,-2,-1, 0, 1, 2, 3/
	DATA (actIndex(1,i),i=1,349)/-10,-10,
     &	-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-8,-8,-8,
     &	-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,
     &	-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,
     &	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-3,-3,-3,-3,
     &	-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
     &	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 
     &	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
     &	 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
     &	 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
     &	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     &	 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
     &	 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     &	 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     &	 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10/ 

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add influence functions of neighboring actuators

	imin=MAX((iRho-idLim),-nzc)
	imax=MIN((iRho+idLim), nzc)
	jmin=MAX((jRho-idLim),-nzc)
	jmax=MIN((jRho+idLim), nzc)
	fh=0d0
	DO j=jmin,jmax
   	  DO i=imin,imax
   	    nearAct=actMap(i,j)
   	    IF (nearAct.GT.0) THEN
   	      IF (coef(nearAct).NE.0d0) THEN
   	        xAct=dAct*actIndex(1,nearAct)
   	        yAct=dAct*actIndex(2,nearAct)
   	        dx=x-xAct
   	        dx2=dx*dx
   	        dx4=dx2*dx2
   	        dx8=dx4*dx4
   	        dy=y-yAct
   	        dy2=dy*dy
   	        dy4=dy2*dy2
   	        dy8=dy4*dy4
   	        dn=(dx8+dy8)
   	        dp=dn**shapeOrdInv
c   	        dn=((x-xAct)**8+(y-yAct)**8)
c   	        dp=SQRT(SQRT(SQRT(dn)))
c   	        dn=((x-xAct)**2+(y-yAct)**2)
c   	        dp=SQRT(dn)
   	        bd=beta*dp
   	        IF (bd.LT.czc) 
     &	          fh=fh+coef(nearAct)*EXP(-bd)*SIN(bd+pio4)
   	      END IF
   	    END IF
   	  END DO
   	END DO

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SDM349B=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SDM349B

C***********************************************************************
C***********************************************************************

	SUBROUTINE NDM349B(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	dFdrho,d2Fdrho2)

C  Computes surface function for 349 actuator, square array DM

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER i,j,actMap(-10:10,-10:10),actIndex(2,349),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,iAxis,jAxis,shapeOrder,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),rhoAx,
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),dx2,dx4,dx7,dx8,dy2,
     &	dy4,dy7,dy8,dFdrho(3),d2Fdrho2(3,3),esbd,ecbd,ebd,beta2,d,bod,
     &	dVec(3),CdF,CdF2,dx,dy,shapeOrdInv,dn,czc

	DATA shapeOrder/8/
	DATA shapeOrdInv/0.125d0/
	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/10/
	DATA pio4/0.78539816339745d0/
	DATA (actMap(-10,i),i=-10,10)/0,0,0,0,0,0,0,1,2,3,4,5,6,
     &	7,0,0,0,0,0,0,0/
	DATA (actMap(-9,i),i=-10,10)/0,0,0,0,0,8,9,10,11,12,13,14,15,16,
     &	17,18,0,0,0,0,0/
	DATA (actMap(-8,i),i=-10,10)/0,0,0,0,19,20,21,22,23,24,25,26,27,
     &	28,29,30,31,0,0,0,0/
	DATA (actMap(-7,i),i=-10,10)/0,0,0,32,33,34,35,36,37,38,39,40,
     &	41,42,43,44,45,46,0,0,0/
	DATA (actMap(-6,i),i=-10,10)/0,0,47,48,49,50,51,52,53,54,55,
     &	56,57,58,59,60,61,62,63,0,0/
	DATA (actMap(-5,i),i=-10,10)/0,64,65,66,67,68,69,70,71,72,73,
     &	74,75,76,77,78,79,80,81,82,0/
	DATA (actMap(-4,i),i=-10,10)/0,83,84,85,86,87,88,89,90,91,92,
     &	93,94,95,96,97,98,99,100,101,0/
	DATA (actMap(-3,i),i=-10,10)/102,103,104,105,106,107,108,109,
     &	110,111,112,113,114,115,116,117,118,119,120,121,122/
	DATA (actMap(-2,i),i=-10,10)/123,124,125,126,127,128,129,130,
     &	131,132,133,134,135,136,137,138,139,140,141,142,143/
	DATA (actMap(-1,i),i=-10,10)/144,145,146,147,148,149,150,151,
     &	152,153,154,155,156,157,158,159,160,161,162,163,164/
	DATA (actMap(0,i),i=-10,10)/165,166,167,168,169,170,171,172,
     &	173,174,175,176,177,178,179,180,181,182,183,184,185/
	DATA (actMap(1,i),i=-10,10)/186,187,188,189,190,191,192,193,
     &	194,195,196,197,198,199,200,201,202,203,204,205,206/
	DATA (actMap(2,i),i=-10,10)/207,208,209,210,211,212,213,214,
     &	215,216,217,218,219,220,221,222,223,224,225,226,227/
	DATA (actMap(3,i),i=-10,10)/228,229,230,231,232,233,234,235,
     &	236,237,238,239,240,241,242,243,244,245,246,247,248/
	DATA (actMap(4,i),i=-10,10)/0,249,250,251,252,253,254,255,
     &	256,257,258,259,260,261,262,263,264,265,266,267,0/
	DATA (actMap(5,i),i=-10,10)/0,268,269,270,271,272,273,274,
     &	275,276,277,278,279,280,281,282,283,284,285,286,0/
	DATA (actMap(6,i),i=-10,10)/0,0,287,288,289,290,291,292,
     &	293,294,295,296,297,298,299,300,301,302,303,0,0/
	DATA (actMap(7,i),i=-10,10)/0,0,0,304,305,306,307,308,309,
     &	310,311,312,313,314,315,316,317,318,0,0,0/
	DATA (actMap(8,i),i=-10,10)/0,0,0,0,319,320,321,322,323,
     &	324,325,326,327,328,329,330,331,0,0,0,0/
	DATA (actMap(9,i),i=-10,10)/0,0,0,0,0,332,333,334,335,336,
     &	337,338,339,340,341,342,0,0,0,0,0/
	DATA (actMap(10,i),i=-10,10)/0,0,0,0,0,0,0,343,344,345,346,
     &	347,348,349,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,349)/-3,-2,-1, 0, 1, 2, 3,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-6,-5,-4,-3,-2,-1, 0,
     &	1, 2, 3, 4, 5, 6,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-10,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,
     &	-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-7,-6,-5,-4,-3,-2,-1, 0,
     &	 1, 2, 3, 4, 5, 6, 7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-3,-2,-1, 0, 1, 2, 3/
	DATA (actIndex(1,i),i=1,349)/-10,-10,
     &	-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-8,-8,-8,
     &	-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,
     &	-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,
     &	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-3,-3,-3,-3,
     &	-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
     &	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 
     &	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
     &	 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
     &	 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
     &	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     &	 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
     &	 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     &	 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     &	 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10/ 

	beta2=beta*beta

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add neighboring actuator influence functions and derivatives

	imin=MAX((iRho-idLim),-10)
	imax=MIN((iRho+idLim), 10)
	jmin=MAX((jRho-idLim),-10)
	jmax=MIN((jRho+idLim), 10)
	fh=0d0
	IF (ifBuild) THEN

	  PAUSE 'Build not supported for this DM'

   	ELSE
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          dx=x-xAct
   	          dx2=dx*dx
   	          dx4=dx2*dx2
   	          dx8=dx4*dx4
   	          dy=y-yAct
   	          dy2=dy*dy
   	          dy4=dy2*dy2
   	          dy8=dy4*dy4
   	          dn=(dx8+dy8)
   	          d=dn**shapeOrdInv
c   	          dn=(dx*dx+dy*dy)
c   	          d=SQRT(dn)
   	          IF (d.lt.1d-12) THEN
   	            ebd=coef(nearAct)
    	            esbd=ebd*0.70710678118655d0
  	            fh=fh+esbd
   	          ELSE
   	            bd=beta*d
    	            IF (bd.LT.czc) THEN
   	              bod=beta*(d**(1-shapeOrder))
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
  	              fh=fh+esbd
	              CdF=bod*(ecbd-esbd)
   	              dx7=dx4*dx2*dx
   	              dy7=dy4*dy2*dy
 	              DO jAxis=1,3
 	                dFdrho(jAxis)=dFdrho(jAxis)
     & 	                +CdF*(dx7*xhat(jAxis)+dy7*yhat(jAxis))
   	              END DO
  	            END IF
  	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	END IF

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NDM349B

C***********************************************************************
C***********************************************************************

	FUNCTION SDM349I(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

C  Computes surface function for 349 actuator, square array DM,
C  units are inches

	IMPLICIT NONE
	INTEGER i,actMap(-10:10,-10:10),actIndex(2,349),iAct,
     &	i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SDM349I,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,x,y,zc(3),zc2,psitzc,tiltx,tilty,u00,u01,u10,u11,
     &	xi,yj

	DATA (actMap(-10,i),i=-10,10)/0,0,0,0,0,0,0,1,2,3,4,5,6,
     &	7,0,0,0,0,0,0,0/
	DATA (actMap(-9,i),i=-10,10)/0,0,0,0,0,8,9,10,11,12,13,14,15,16,
     &	17,18,0,0,0,0,0/
	DATA (actMap(-8,i),i=-10,10)/0,0,0,0,19,20,21,22,23,24,25,26,27,
     &	28,29,30,31,0,0,0,0/
	DATA (actMap(-7,i),i=-10,10)/0,0,0,32,33,34,35,36,37,38,39,40,
     &	41,42,43,44,45,46,0,0,0/
	DATA (actMap(-6,i),i=-10,10)/0,0,47,48,49,50,51,52,53,54,55,
     &	56,57,58,59,60,61,62,63,0,0/
	DATA (actMap(-5,i),i=-10,10)/0,64,65,66,67,68,69,70,71,72,73,
     &	74,75,76,77,78,79,80,81,82,0/
	DATA (actMap(-4,i),i=-10,10)/0,83,84,85,86,87,88,89,90,91,92,
     &	93,94,95,96,97,98,99,100,101,0/
	DATA (actMap(-3,i),i=-10,10)/102,103,104,105,106,107,108,109,
     &	110,111,112,113,114,115,116,117,118,119,120,121,122/
	DATA (actMap(-2,i),i=-10,10)/123,124,125,126,127,128,129,130,
     &	131,132,133,134,135,136,137,138,139,140,141,142,143/
	DATA (actMap(-1,i),i=-10,10)/144,145,146,147,148,149,150,151,
     &	152,153,154,155,156,157,158,159,160,161,162,163,164/
	DATA (actMap(0,i),i=-10,10)/165,166,167,168,169,170,171,172,
     &	173,174,175,176,177,178,179,180,181,182,183,184,185/
	DATA (actMap(1,i),i=-10,10)/186,187,188,189,190,191,192,193,
     &	194,195,196,197,198,199,200,201,202,203,204,205,206/
	DATA (actMap(2,i),i=-10,10)/207,208,209,210,211,212,213,214,
     &	215,216,217,218,219,220,221,222,223,224,225,226,227/
	DATA (actMap(3,i),i=-10,10)/228,229,230,231,232,233,234,235,
     &	236,237,238,239,240,241,242,243,244,245,246,247,248/
	DATA (actMap(4,i),i=-10,10)/0,249,250,251,252,253,254,255,
     &	256,257,258,259,260,261,262,263,264,265,266,267,0/
	DATA (actMap(5,i),i=-10,10)/0,268,269,270,271,272,273,274,
     &	275,276,277,278,279,280,281,282,283,284,285,286,0/
	DATA (actMap(6,i),i=-10,10)/0,0,287,288,289,290,291,292,
     &	293,294,295,296,297,298,299,300,301,302,303,0,0/
	DATA (actMap(7,i),i=-10,10)/0,0,0,304,305,306,307,308,309,
     &	310,311,312,313,314,315,316,317,318,0,0,0/
	DATA (actMap(8,i),i=-10,10)/0,0,0,0,319,320,321,322,323,
     &	324,325,326,327,328,329,330,331,0,0,0,0/
	DATA (actMap(9,i),i=-10,10)/0,0,0,0,0,332,333,334,335,336,
     &	337,338,339,340,341,342,0,0,0,0,0/
	DATA (actMap(10,i),i=-10,10)/0,0,0,0,0,0,0,343,344,345,346,
     &	347,348,349,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,349)/-3,-2,-1, 0, 1, 2, 3,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-6,-5,-4,-3,-2,-1, 0,
     &	1, 2, 3, 4, 5, 6,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-10,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,
     &	-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-7,-6,-5,-4,-3,-2,-1, 0,
     &	 1, 2, 3, 4, 5, 6, 7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-3,-2,-1, 0, 1, 2, 3/
	DATA (actIndex(1,i),i=1,349)/-10,-10,
     &	-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-8,-8,-8,
     &	-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,
     &	-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,
     &	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-3,-3,-3,-3,
     &	-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
     &	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 
     &	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
     &	 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
     &	 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
     &	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     &	 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
     &	 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     &	 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     &	 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10/ 

C  Compute point-of-incidence in surface coords

	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
	END DO
	xi=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))/dAct
	yj=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))/dAct

C  Find actuator coordinates of ray incidence point

	i0=IDFLOOR(xi)
	i1=i0+1
	j0=IDFLOOR(yj)
	j1=j0+1
	x=xi-DBLE(i0)
	y=yj-DBLE(j0)

C  Compute excess surface height

	IF ((i0.GE.-10).AND.(i0.LE.10).AND.(j0.GE.-10).AND.(j0.LE.10))
     &	THEN
   	  iAct=actMap(i0,j0)
   	  IF (iAct.GT.0) THEN
	    u00=coef(iAct)
	  ELSE
	    u00=0d0
	  END IF
	ELSE
	  u00=0d0
	END IF

	IF ((i1.GE.-10).AND.(i1.LE.10).AND.(j0.GE.-10).AND.(j0.LE.10))
     &	THEN
   	  iAct=actMap(i1,j0)
   	  IF (iAct.GT.0) THEN
	    u10=coef(iAct)
	  ELSE
	    u10=0d0
	  END IF
	ELSE
	  u10=0d0
	END IF

	IF ((i0.GE.-10).AND.(i0.LE.10).AND.(j1.GE.-10).AND.(j1.LE.10))
     &	THEN
   	  iAct=actMap(i0,j1)
   	  IF (iAct.GT.0) THEN
	    u01=coef(iAct)
	  ELSE
	    u01=0d0
	  END IF
	ELSE
	  u01=0d0
	END IF

	IF ((i1.GE.-10).AND.(i1.LE.10).AND.(j1.GE.-10).AND.(j1.LE.10))
     &	THEN
   	  iAct=actMap(i1,j1)
   	  IF (iAct.GT.0) THEN
	    u11=coef(iAct)
	  ELSE
	    u11=0d0
	  END IF
	ELSE
	  u11=0d0
	END IF

	CALL INTNORM(x,y,u00,u01,u10,u11,fh,tiltx,tilty)

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SDM349I=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SDM349I

C***********************************************************************
C***********************************************************************

	SUBROUTINE NDM349I(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for 349 actuator, square array DM

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER i,actMap(-10:10,-10:10),actIndex(2,349),
     &	iAxis,jAxis,iAct,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),
     &	fh,x,y,zc(3),beta,dAct,xi,yj,
     &	dFdrho(3),d2Fdrho2(3,3),tiltx,tilty,u00,u01,u10,u11

	DATA (actMap(-10,i),i=-10,10)/0,0,0,0,0,0,0,1,2,3,4,5,6,
     &	7,0,0,0,0,0,0,0/
	DATA (actMap(-9,i),i=-10,10)/0,0,0,0,0,8,9,10,11,12,13,14,15,16,
     &	17,18,0,0,0,0,0/
	DATA (actMap(-8,i),i=-10,10)/0,0,0,0,19,20,21,22,23,24,25,26,27,
     &	28,29,30,31,0,0,0,0/
	DATA (actMap(-7,i),i=-10,10)/0,0,0,32,33,34,35,36,37,38,39,40,
     &	41,42,43,44,45,46,0,0,0/
	DATA (actMap(-6,i),i=-10,10)/0,0,47,48,49,50,51,52,53,54,55,
     &	56,57,58,59,60,61,62,63,0,0/
	DATA (actMap(-5,i),i=-10,10)/0,64,65,66,67,68,69,70,71,72,73,
     &	74,75,76,77,78,79,80,81,82,0/
	DATA (actMap(-4,i),i=-10,10)/0,83,84,85,86,87,88,89,90,91,92,
     &	93,94,95,96,97,98,99,100,101,0/
	DATA (actMap(-3,i),i=-10,10)/102,103,104,105,106,107,108,109,
     &	110,111,112,113,114,115,116,117,118,119,120,121,122/
	DATA (actMap(-2,i),i=-10,10)/123,124,125,126,127,128,129,130,
     &	131,132,133,134,135,136,137,138,139,140,141,142,143/
	DATA (actMap(-1,i),i=-10,10)/144,145,146,147,148,149,150,151,
     &	152,153,154,155,156,157,158,159,160,161,162,163,164/
	DATA (actMap(0,i),i=-10,10)/165,166,167,168,169,170,171,172,
     &	173,174,175,176,177,178,179,180,181,182,183,184,185/
	DATA (actMap(1,i),i=-10,10)/186,187,188,189,190,191,192,193,
     &	194,195,196,197,198,199,200,201,202,203,204,205,206/
	DATA (actMap(2,i),i=-10,10)/207,208,209,210,211,212,213,214,
     &	215,216,217,218,219,220,221,222,223,224,225,226,227/
	DATA (actMap(3,i),i=-10,10)/228,229,230,231,232,233,234,235,
     &	236,237,238,239,240,241,242,243,244,245,246,247,248/
	DATA (actMap(4,i),i=-10,10)/0,249,250,251,252,253,254,255,
     &	256,257,258,259,260,261,262,263,264,265,266,267,0/
	DATA (actMap(5,i),i=-10,10)/0,268,269,270,271,272,273,274,
     &	275,276,277,278,279,280,281,282,283,284,285,286,0/
	DATA (actMap(6,i),i=-10,10)/0,0,287,288,289,290,291,292,
     &	293,294,295,296,297,298,299,300,301,302,303,0,0/
	DATA (actMap(7,i),i=-10,10)/0,0,0,304,305,306,307,308,309,
     &	310,311,312,313,314,315,316,317,318,0,0,0/
	DATA (actMap(8,i),i=-10,10)/0,0,0,0,319,320,321,322,323,
     &	324,325,326,327,328,329,330,331,0,0,0,0/
	DATA (actMap(9,i),i=-10,10)/0,0,0,0,0,332,333,334,335,336,
     &	337,338,339,340,341,342,0,0,0,0,0/
	DATA (actMap(10,i),i=-10,10)/0,0,0,0,0,0,0,343,344,345,346,
     &	347,348,349,0,0,0,0,0,0,0/
	DATA (actIndex(2,i),i=1,349)/-3,-2,-1, 0, 1, 2, 3,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-6,-5,-4,-3,-2,-1, 0,
     &	1, 2, 3, 4, 5, 6,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-10,-9,-8,-7,-6,-5,-4,
     &	-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,-6,
     &	-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,-7,
     &	-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,-9,
     &	-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &	-9,-8,-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-9,-8,
     &	-7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-8,-7,-6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-7,-6,-5,-4,-3,-2,-1, 0,
     &	 1, 2, 3, 4, 5, 6, 7,-6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6,-5,
     &	-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,-3,-2,-1, 0, 1, 2, 3/
	DATA (actIndex(1,i),i=1,349)/-10,-10,
     &	-10,-10,-10,-10,-10,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-8,-8,-8,
     &	-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,
     &	-7,-7,-7,-7,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,
     &	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-4,-4,
     &	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-3,-3,-3,-3,
     &	-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2,
     &	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
     &	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 
     &	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     &	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
     &	 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
     &	 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
     &	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     &	 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
     &	 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     &	 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
     &	 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10/ 

C  Compute point-of-incidence in surface coords

	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
	END DO
	xi=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))/dAct
	yj=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))/dAct

C  Find actuator coordinates of ray incidence point

	i0=IDFLOOR(xi)
	i1=i0+1
	j0=IDFLOOR(yj)
	j1=j0+1
	x=xi-DBLE(i0)
	y=yj-DBLE(j0)

C  Compute excess surface height

	IF ((i0.GE.-10).AND.(i0.LE.10).AND.(j0.GE.-10).AND.(j0.LE.10))
     &	THEN
   	  iAct=actMap(i0,j0)
   	  IF (iAct.GT.0) THEN
	    u00=coef(iAct)
	  ELSE
	    u00=0d0
	  END IF
	ELSE
	  u00=0d0
	END IF

	IF ((i1.GE.-10).AND.(i1.LE.10).AND.(j0.GE.-10).AND.(j0.LE.10))
     &	THEN
   	  iAct=actMap(i1,j0)
   	  IF (iAct.GT.0) THEN
	    u10=coef(iAct)
	  ELSE
	    u10=0d0
	  END IF
	ELSE
	  u10=0d0
	END IF

	IF ((i0.GE.-10).AND.(i0.LE.10).AND.(j1.GE.-10).AND.(j1.LE.10))
     &	THEN
   	  iAct=actMap(i0,j1)
   	  IF (iAct.GT.0) THEN
	    u01=coef(iAct)
	  ELSE
	    u01=0d0
	  END IF
	ELSE
	  u01=0d0
	END IF

	IF ((i1.GE.-10).AND.(i1.LE.10).AND.(j1.GE.-10).AND.(j1.LE.10))
     &	THEN
   	  iAct=actMap(i1,j1)
   	  IF (iAct.GT.0) THEN
	    u11=coef(iAct)
	  ELSE
	    u11=0d0
	  END IF
	ELSE
	  u11=0d0
	END IF

	CALL INTNORM(x,y,u00,u01,u10,u11,fh,tiltx,tilty)
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
	END SUBROUTINE NDM349I

C***********************************************************************

	FUNCTION SHSTPM(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=149, UDSctr=75, udsfile='hstpm.bin')
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm
	COMMON /UDSdata1/ UDSmat,UDSnorm
	SAVE /UDSdata1/

	INTEGER i,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SHSTPM,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,x,y,zc(3),zc2,psitzc,tiltx,tilty,xi,yj,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
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

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSna)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
	END IF

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SHSTPM=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SHSTPM

C***********************************************************************

	SUBROUTINE NHSTPM(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=149, UDSctr=75, udsfile='hstpm.bin')
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm
	COMMON /UDSdata1/ UDSmat,UDSnorm
	SAVE /UDSdata1/

	LOGICAL ifBuild
	INTEGER i,iAxis,jAxis,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),
     &	fh,x,y,zc(3),beta,dAct,xi,yj,
     &	dFdrho(3),d2Fdrho2(3,3),tiltx,tilty,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
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

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSna)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
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
	END SUBROUTINE NHSTPM

C***********************************************************************
C***********************************************************************

	FUNCTION SHSTSM(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=115, UDSctr=58, udsfile='hstsm.bin')
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm
	COMMON /UDSdata2/ UDSmat,UDSnorm
	SAVE /UDSdata2/

	INTEGER i,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SHSTSM,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,x,y,zc(3),zc2,psitzc,tiltx,tilty,xi,yj,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
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

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSna)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
	END IF

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SHSTSM=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SHSTSM

C***********************************************************************
C***********************************************************************

	SUBROUTINE NHSTSM(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=115, UDSctr=58, udsfile='hstsm.bin')
	REAL*8 UDSmat(UDSna,UDSna),UDSnorm
	COMMON /UDSdata2/ UDSmat,UDSnorm
	SAVE /UDSdata2/

	LOGICAL ifBuild
	INTEGER i,iAxis,jAxis,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),
     &	fh,x,y,zc(3),beta,dAct,xi,yj,
     &	dFdrho(3),d2Fdrho2(3,3),tiltx,tilty,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
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

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSna)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTNORM(x,y,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
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
	END SUBROUTINE NHSTSM

C***********************************************************************

	FUNCTION SDMRAD(Kc,Kr,psi,pv,prel,ihat,drad,dq,minrad,coef,rad,
     &	                pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSnb,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=31, UDSnb=16, UDSctr=0)
	REAL*8 UDSmat(UDSna,UDSnb),UDSnorm
	COMMON /UDSdata3/ UDSmat,UDSnorm
	SAVE /UDSdata3/

	INTEGER i,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SDMRAD,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,x,y,zc(3),zc2,psitzc,tiltx,tilty,xi,yj,diCtr,
     &	r,q,ri,qj,drad,qq,dq,qhat(3),rhat(3),drq,rdq,minrad

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
	END DO
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))
	rad=SQRT(x*x+y*y)
	qq=ATAN2(y,-x)
	IF (qq.LT.0d0) qq=qq+6.28318530717959d0

	IF (rad.LT.minrad) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	  DO i=1,3
	    zc(i)=rho(i)
   	  END DO
	  zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	  psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	  SDMRAD=zc2+(Kc*psitzc+2d0*Kr)*psitzc
	  RETURN
	END IF

	IF (rad.GE.1d-22) THEN
	  DO i=1,3
	    rhat(i)=rhom(i)/rad
	  END DO
	ELSE
	  DO i=1,3
	    rhat(i)=xhat(i)
	  END DO
	END IF
	CALL DXPROD(qhat,zhat,rhat)

C  Find actuator coordinates of ray incidence point

	ri=(rad-minrad)/drad
	qj=qq/dq    ! ri and qj are polar coords of ray at surface - jzlou
	i0=IDFLOOR(ri)
	i1=i0+1
	j0=IDFLOOR(qj)
	j1=j0+1
	q=qj-DBLE(j0)
	IF (j0.EQ.0) j0=UDSnb
	r=ri-DBLE(i0)

C  Compute excess surface height

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSnb)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTRAD(r,q,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	              UDSmat(i1,j1),fh,tiltx,tilty)
	END IF

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SDMRAD=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SDMRAD

C***********************************************************************

	SUBROUTINE NDMRAD(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	 drad,dq,minrad,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	 dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for HST primary mirror

	IMPLICIT NONE

	INTEGER UDSna,UDSnb,UDSctr
	CHARACTER*12 udsfile
	PARAMETER (UDSna=31, UDSnb=16, UDSctr=0)
	REAL*8 UDSmat(UDSna,UDSnb),UDSnorm,inVec(UDSna)
	COMMON /UDSdata3/ UDSmat,UDSnorm
	SAVE /UDSdata3/

	LOGICAL ifBuild
	INTEGER i,iAxis,jAxis,i1,i0,j1,j0
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),
     &	fh,x,y,zc(3),beta,dAct,xi,yj,
     &	r,q,ri,qj,drad,qq,dq,qhat(3),rhat(3),drq,rdq,minrad,
     &	dFdrho(3),d2Fdrho2(3,3),tiltx,tilty,diCtr

C  Compute point-of-incidence in surface coords

	diCtr=UDSctr
	DO i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
	END DO
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))
	rad=SQRT(x*x+y*y)
	qq=ATAN2(y,-x)
	IF (qq.LT.0d0) qq=qq+6.28318530717959d0

	IF (rad.LT.minrad) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	  DO jAxis=1,3
	    zc(jAxis)=rho(jAxis)
 	    dFdrho(jAxis)=0d0
 	    DO iAxis=1,3
   	      d2Fdrho2(iAxis,jAxis)=0d0
   	    END DO
   	  END DO
	  RETURN
	END IF

	IF (rad.GE.1d-22) THEN
	  DO i=1,3
	    rhat(i)=rhom(i)/rad
	  END DO
	ELSE
	  DO i=1,3
	    rhat(i)=xhat(i)
	  END DO
	END IF
	CALL DXPROD(qhat,zhat,rhat)
	rdq=rad*dq

C  Find actuator coordinates of ray incidence point

	ri=(rad-minrad)/drad
	qj=qq/dq
	i0=IDFLOOR(ri)
	i1=i0+1
	j0=IDFLOOR(qj)
	j1=j0+1
	q=qj-DBLE(j0)
	IF (j0.EQ.0) j0=UDSnb
	r=ri-DBLE(i0)

C  Compute excess surface height

	IF ((i0.LT.1).OR.(i1.GT.UDSna).OR.(j0.LT.1).OR.
     &	(j1.GT.UDSnb)) THEN
	  fh=0d0
	  tiltx=0d0
	  tilty=0d0
	ELSE
	  CALL INTRAD(r,q,UDSmat(i0,j0),UDSmat(i0,j1),UDSmat(i1,j0),
     &	  UDSmat(i1,j1),fh,tiltx,tilty)
	END IF

	tiltx=tiltx/drad
	IF (rdq.GT.1d-10) THEN
	  tilty=tilty/rdq
	ELSE
	  tilty=0
	END IF

	DO jAxis=1,3
 	  dFdrho(jAxis)=tiltx*rhat(jAxis)+tilty*qhat(jAxis)
 	  DO iAxis=1,3
   	    d2Fdrho2(iAxis,jAxis)=0d0
   	  END DO
   	END DO

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NDMRAD

C***********************************************************************

      SUBROUTINE INTRAD(X,Y,U00,U01,U10,U11,U,DUDX,DUDY)
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
      END SUBROUTINE INTRAD

C***********************************************************************

	FUNCTION SDM1313(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,rad,
     &	  pmon,rhoc,xhat,yhat,zhat,L)
	use math_mod

c	dAct=3.5d0
c	beta=0.5129120604742d0

C  Computes surface function for 1313 actuator, square array DM

	IMPLICIT NONE
	INTEGER i,j,actMap(-20:20,-20:20),actIndex(2,1313),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SDM1313,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),zc2,psitzc,czc
	COMMON /DM1313/ actIndex,actMap
	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/20/
	DATA pio4/0.78539816339745d0/

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add influence functions of neighboring actuators

	imin=MAX((iRho-idLim),-nzc)
	imax=MIN((iRho+idLim), nzc)
	jmin=MAX((jRho-idLim),-nzc)
	jmax=MIN((jRho+idLim), nzc)
	fh=0d0
	DO j=jmin,jmax
   	  DO i=imin,imax
   	    nearAct=actMap(i,j)
   	    IF (nearAct.GT.0) THEN
   	      IF (coef(nearAct).NE.0d0) THEN
   	        xAct=dAct*actIndex(1,nearAct)
   	        yAct=dAct*actIndex(2,nearAct)
   	        bd=beta*SQRT((x-xAct)**2+(y-yAct)**2)
   	        IF (bd.LT.czc) 
     &	          fh=fh+coef(nearAct)*EXP(-bd)*SIN(bd+pio4)
   	      END IF
   	    END IF
   	  END DO
   	END DO

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SDM1313=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SDM1313

C***********************************************************************
C***********************************************************************

	SUBROUTINE NDM1313(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for 1313 actuator, square array DM

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER i,j,actMap(-20:20,-20:20),actIndex(2,1313),nearAct,
     &	iRho,jRho,imin,imax,jmin,jmax,idLim,iAxis,jAxis,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),rhoAx,
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),
     &	dFdrho(3),d2Fdrho2(3,3),esbd,ecbd,ebd,beta2,d,bod,dVec(3),
     &	CdF,CdF2,dx,dy,czc
	COMMON /DM1313/ actIndex,actMap

	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/20/
	DATA pio4/0.78539816339745d0/

	beta2=beta*beta

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add neighboring actuator influence functions and derivatives

	imin=MAX((iRho-idLim),-nzc)
	imax=MIN((iRho+idLim), nzc)
	jmin=MAX((jRho-idLim),-nzc)
	jmax=MIN((jRho+idLim), nzc)
	fh=0d0
	IF (ifBuild) THEN
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
  	    DO iAxis=1,3
	      d2Fdrho2(iAxis,jAxis)=0d0
   	    END DO
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          d=SQRT((x-xAct)**2+(y-yAct)**2)
   	          bd=beta*d
   	          IF (bd.LT.czc) THEN
   	            DO iAxis=1,3
   	              rhoAx=xAct*xhat(iAxis)+
     &                yAct*yhat(iAxis)-rhoc(iAxis)
c 	              dVec(iAxis)=rhom(iAxis)-rhoAx
  	              dVec(iAxis)=(1d0-zhat(iAxis))*rhom(iAxis)-rhoAx
   	            END DO
   	            IF (d.lt.1d-12) THEN
   	              ebd=coef(nearAct)
    	              esbd=ebd*0.70710678118655d0
    	              ecbd=ebd*0.70710678118655d0
   	              fh=fh+esbd
     	              CdF2=-2d0*beta*beta*ecbd
 	              DO iAxis=1,3
   	                d2Fdrho2(iAxis,iAxis)=d2Fdrho2(iAxis,iAxis)+CdF2
   	              END DO
    	            ELSE
   	              bod=beta/d
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
   	              fh=fh+esbd
   	              CdF=bod*(ecbd-esbd)
   	              CdF2=-2d0*bod*bod*ecbd
   	              DO jAxis=1,3
   	                dFdrho(jAxis)=dFdrho(jAxis)+CdF*dVec(jAxis)
   	                d2Fdrho2(jAxis,jAxis)=d2Fdrho2(jAxis,jAxis)+CdF
   	                DO iAxis=1,3
   	                  d2Fdrho2(iAxis,jAxis)=d2Fdrho2(iAxis,jAxis)+
     &   	          (CdF2-CdF/d/d)*dVec(iAxis)*dVec(jAxis)
   	                END DO
   	              END DO
   	            END IF
   	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	ELSE
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          dx=x-xAct
   	          dy=y-yAct
   	          d=SQRT(dx*dx+dy*dy)
   	          IF (d.lt.1d-12) THEN
   	            ebd=coef(nearAct)
    	            esbd=ebd*0.70710678118655d0
  	            fh=fh+esbd
   	          ELSE
   	            bd=beta*d
    	            IF (bd.LT.czc) THEN
  	              bod=beta/d
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
  	              fh=fh+esbd
	              CdF=bod*(ecbd-esbd)
 	              DO jAxis=1,3
 	                dFdrho(jAxis)=dFdrho(jAxis)
     & 	                +CdF*(dx*xhat(jAxis)+dy*yhat(jAxis))
   	              END DO
  	            END IF
  	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	END IF

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NDM1313

C***********************************************************************
C***********************************************************************

	BLOCK DATA DM1313BD
	IMPLICIT NONE
	INTEGER aiv(2626),amv(1681),i
	COMMON /DM1313/aiv,amv
	DATA (aiv(i),i=1,160)/
     &	-20,-4,-20,-3,-20,-2,-20,-1,-20,0,-20,1,-20,2,-20,3,
     &	-20,4,-19,-7,-19,-6,-19,-5,-19,-4,-19,-3,-19,-2,-19,-1,
     &	-19,0,-19,1,-19,2,-19,3,-19,4,-19,5,-19,6,-19,7,
     &	-18,-9,-18,-8,-18,-7,-18,-6,-18,-5,-18,-4,-18,-3,-18,-2,
     &	-18,-1,-18,0,-18,1,-18,2,-18,3,-18,4,-18,5,-18,6,
     &	-18,7,-18,8,-18,9,-17,-11,-17,-10,-17,-9,-17,-8,-17,-7,
     &	-17,-6,-17,-5,-17,-4,-17,-3,-17,-2,-17,-1,-17,0,-17,1,
     &	-17,2,-17,3,-17,4,-17,5,-17,6,-17,7,-17,8,-17,9,
     &	-17,10,-17,11,-16,-12,-16,-11,-16,-10,-16,-9,-16,-8,-16,-7,
     &	-16,-6,-16,-5,-16,-4,-16,-3,-16,-2,-16,-1,-16,0,-16,1/
	DATA (aiv(i),i=161,320)/
     &	-16,2,-16,3,-16,4,-16,5,-16,6,-16,7,-16,8,-16,9,
     &	-16,10,-16,11,-16,12,-15,-13,-15,-12,-15,-11,-15,-10,-15,-9,
     &	-15,-8,-15,-7,-15,-6,-15,-5,-15,-4,-15,-3,-15,-2,-15,-1,
     &	-15,0,-15,1,-15,2,-15,3,-15,4,-15,5,-15,6,-15,7,
     &	-15,8,-15,9,-15,10,-15,11,-15,12,-15,13,-14,-14,-14,-13,
     &	-14,-12,-14,-11,-14,-10,-14,-9,-14,-8,-14,-7,-14,-6,-14,-5,
     &	-14,-4,-14,-3,-14,-2,-14,-1,-14,0,-14,1,-14,2,-14,3,
     &	-14,4,-14,5,-14,6,-14,7,-14,8,-14,9,-14,10,-14,11,
     &	-14,12,-14,13,-14,14,-13,-15,-13,-14,-13,-13,-13,-12,-13,-11,
     &	-13,-10,-13,-9,-13,-8,-13,-7,-13,-6,-13,-5,-13,-4,-13,-3/
	DATA (aiv(i),i=321, 480)/
     &	-13,-2,-13,-1,-13,0,-13,1,-13,2,-13,3,-13,4,-13,5,
     &	-13,6,-13,7,-13,8,-13,9,-13,10,-13,11,-13,12,-13,13,
     &	-13,14,-13,15,-12,-16,-12,-15,-12,-14,-12,-13,-12,-12,-12,-11,
     &	-12,-10,-12,-9,-12,-8,-12,-7,-12,-6,-12,-5,-12,-4,-12,-3,
     &	-12,-2,-12,-1,-12,0,-12,1,-12,2,-12,3,-12,4,-12,5,
     &	-12,6,-12,7,-12,8,-12,9,-12,10,-12,11,-12,12,-12,13,
     &	-12,14,-12,15,-12,16,-11,-17,-11,-16,-11,-15,-11,-14,-11,-13,
     &	-11,-12,-11,-11,-11,-10,-11,-9,-11,-8,-11,-7,-11,-6,-11,-5,
     &	-11,-4,-11,-3,-11,-2,-11,-1,-11,0,-11,1,-11,2,-11,3,
     &	-11,4,-11,5,-11,6,-11,7,-11,8,-11,9,-11,10,-11,11/
	DATA (aiv(i),i=481, 640)/
     &	-11,12,-11,13,-11,14,-11,15,-11,16,-11,17,-10,-17,-10,-16,
     &	-10,-15,-10,-14,-10,-13,-10,-12,-10,-11,-10,-10,-10,-9,-10,-8,
     &	-10,-7,-10,-6,-10,-5,-10,-4,-10,-3,-10,-2,-10,-1,-10,0,
     &	-10,1,-10,2,-10,3,-10,4,-10,5,-10,6,-10,7,-10,8,
     &	-10,9,-10,10,-10,11,-10,12,-10,13,-10,14,-10,15,-10,16,
     &	-10,17,-9,-18,-9,-17,-9,-16,-9,-15,-9,-14,-9,-13,-9,-12,
     &	-9,-11,-9,-10,-9,-9,-9,-8,-9,-7,-9,-6,-9,-5,-9,-4,
     &	-9,-3,-9,-2,-9,-1,-9,0,-9,1,-9,2,-9,3,-9,4,
     &	-9,5,-9,6,-9,7,-9,8,-9,9,-9,10,-9,11,-9,12,
     &	-9,13,-9,14,-9,15,-9,16,-9,17,-9,18,-8,-18,-8,-17/
	DATA (aiv(i),i=641, 800)/
     &	-8,-16,-8,-15,-8,-14,-8,-13,-8,-12,-8,-11,-8,-10,-8,-9,
     &	-8,-8,-8,-7,-8,-6,-8,-5,-8,-4,-8,-3,-8,-2,-8,-1,
     &	-8,0,-8,1,-8,2,-8,3,-8,4,-8,5,-8,6,-8,7,
     &	-8,8,-8,9,-8,10,-8,11,-8,12,-8,13,-8,14,-8,15,
     &	-8,16,-8,17,-8,18,-7,-19,-7,-18,-7,-17,-7,-16,-7,-15,
     &	-7,-14,-7,-13,-7,-12,-7,-11,-7,-10,-7,-9,-7,-8,-7,-7,
     &	-7,-6,-7,-5,-7,-4,-7,-3,-7,-2,-7,-1,-7,0,-7,1,
     &	-7,2,-7,3,-7,4,-7,5,-7,6,-7,7,-7,8,-7,9,
     &	-7,10,-7,11,-7,12,-7,13,-7,14,-7,15,-7,16,-7,17,
     &	-7,18,-7,19,-6,-19,-6,-18,-6,-17,-6,-16,-6,-15,-6,-14/
	DATA (aiv(i),i=801, 960)/
     &	-6,-13,-6,-12,-6,-11,-6,-10,-6,-9,-6,-8,-6,-7,-6,-6,
     &	-6,-5,-6,-4,-6,-3,-6,-2,-6,-1,-6,0,-6,1,-6,2,
     &	-6,3,-6,4,-6,5,-6,6,-6,7,-6,8,-6,9,-6,10,
     &	-6,11,-6,12,-6,13,-6,14,-6,15,-6,16,-6,17,-6,18,
     &	-6,19,-5,-19,-5,-18,-5,-17,-5,-16,-5,-15,-5,-14,-5,-13,
     &	-5,-12,-5,-11,-5,-10,-5,-9,-5,-8,-5,-7,-5,-6,-5,-5,
     &	-5,-4,-5,-3,-5,-2,-5,-1,-5,0,-5,1,-5,2,-5,3,
     &	-5,4,-5,5,-5,6,-5,7,-5,8,-5,9,-5,10,-5,11,
     &	-5,12,-5,13,-5,14,-5,15,-5,16,-5,17,-5,18,-5,19,
     &	-4,-20,-4,-19,-4,-18,-4,-17,-4,-16,-4,-15,-4,-14,-4,-13/
	DATA (aiv(i),i=961,1120)/
     &	-4,-12,-4,-11,-4,-10,-4,-9,-4,-8,-4,-7,-4,-6,-4,-5,
     &	-4,-4,-4,-3,-4,-2,-4,-1,-4,0,-4,1,-4,2,-4,3,
     &	-4,4,-4,5,-4,6,-4,7,-4,8,-4,9,-4,10,-4,11,
     &	-4,12,-4,13,-4,14,-4,15,-4,16,-4,17,-4,18,-4,19,
     &	-4,20,-3,-20,-3,-19,-3,-18,-3,-17,-3,-16,-3,-15,-3,-14,
     &	-3,-13,-3,-12,-3,-11,-3,-10,-3,-9,-3,-8,-3,-7,-3,-6,
     &	-3,-5,-3,-4,-3,-3,-3,-2,-3,-1,-3,0,-3,1,-3,2,
     &	-3,3,-3,4,-3,5,-3,6,-3,7,-3,8,-3,9,-3,10,
     &	-3,11,-3,12,-3,13,-3,14,-3,15,-3,16,-3,17,-3,18,
     &	-3,19,-3,20,-2,-20,-2,-19,-2,-18,-2,-17,-2,-16,-2,-15/
	DATA (aiv(i),i=1121,1280)/
     &	-2,-14,-2,-13,-2,-12,-2,-11,-2,-10,-2,-9,-2,-8,-2,-7,
     &	-2,-6,-2,-5,-2,-4,-2,-3,-2,-2,-2,-1,-2,0,-2,1,
     &	-2,2,-2,3,-2,4,-2,5,-2,6,-2,7,-2,8,-2,9,
     &	-2,10,-2,11,-2,12,-2,13,-2,14,-2,15,-2,16,-2,17,
     &	-2,18,-2,19,-2,20,-1,-20,-1,-19,-1,-18,-1,-17,-1,-16,
     &	-1,-15,-1,-14,-1,-13,-1,-12,-1,-11,-1,-10,-1,-9,-1,-8,
     &	-1,-7,-1,-6,-1,-5,-1,-4,-1,-3,-1,-2,-1,-1,-1,0,
     &	-1,1,-1,2,-1,3,-1,4,-1,5,-1,6,-1,7,-1,8,
     &	-1,9,-1,10,-1,11,-1,12,-1,13,-1,14,-1,15,-1,16,
     &	-1,17,-1,18,-1,19,-1,20,0,-20,0,-19,0,-18,0,-17/
	DATA (aiv(i),i=1281,1440)/
     &	0,-16,0,-15,0,-14,0,-13,0,-12,0,-11,0,-10,0,-9,
     &	0,-8,0,-7,0,-6,0,-5,0,-4,0,-3,0,-2,0,-1,
     &	0,0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,
     &	0,8,0,9,0,10,0,11,0,12,0,13,0,14,0,15,
     &	0,16,0,17,0,18,0,19,0,20,1,-20,1,-19,1,-18,
     &	1,-17,1,-16,1,-15,1,-14,1,-13,1,-12,1,-11,1,-10,
     &	1,-9,1,-8,1,-7,1,-6,1,-5,1,-4,1,-3,1,-2,
     &	1,-1,1,0,1,1,1,2,1,3,1,4,1,5,1,6,
     &	1,7,1,8,1,9,1,10,1,11,1,12,1,13,1,14,
     &	1,15,1,16,1,17,1,18,1,19,1,20,2,-20,2,-19/
	DATA (aiv(i),i=1441,1600)/
     &	2,-18,2,-17,2,-16,2,-15,2,-14,2,-13,2,-12,2,-11,
     &	2,-10,2,-9,2,-8,2,-7,2,-6,2,-5,2,-4,2,-3,
     &	2,-2,2,-1,2,0,2,1,2,2,2,3,2,4,2,5,
     &	2,6,2,7,2,8,2,9,2,10,2,11,2,12,2,13,
     &	2,14,2,15,2,16,2,17,2,18,2,19,2,20,3,-20,
     &	3,-19,3,-18,3,-17,3,-16,3,-15,3,-14,3,-13,3,-12,
     &	3,-11,3,-10,3,-9,3,-8,3,-7,3,-6,3,-5,3,-4,
     &	3,-3,3,-2,3,-1,3,0,3,1,3,2,3,3,3,4,
     &	3,5,3,6,3,7,3,8,3,9,3,10,3,11,3,12,
     &	3,13,3,14,3,15,3,16,3,17,3,18,3,19,3,20/
	DATA (aiv(i),i=1601,1760)/
     &	4,-20,4,-19,4,-18,4,-17,4,-16,4,-15,4,-14,4,-13,
     &	4,-12,4,-11,4,-10,4,-9,4,-8,4,-7,4,-6,4,-5,
     &	4,-4,4,-3,4,-2,4,-1,4,0,4,1,4,2,4,3,
     &	4,4,4,5,4,6,4,7,4,8,4,9,4,10,4,11,
     &	4,12,4,13,4,14,4,15,4,16,4,17,4,18,4,19,
     &	4,20,5,-19,5,-18,5,-17,5,-16,5,-15,5,-14,5,-13,
     &	5,-12,5,-11,5,-10,5,-9,5,-8,5,-7,5,-6,5,-5,
     &	5,-4,5,-3,5,-2,5,-1,5,0,5,1,5,2,5,3,
     &	5,4,5,5,5,6,5,7,5,8,5,9,5,10,5,11,
     &	5,12,5,13,5,14,5,15,5,16,5,17,5,18,5,19/
	DATA (aiv(i),i=1761,1920)/
     &	6,-19,6,-18,6,-17,6,-16,6,-15,6,-14,6,-13,6,-12,
     &	6,-11,6,-10,6,-9,6,-8,6,-7,6,-6,6,-5,6,-4,
     &	6,-3,6,-2,6,-1,6,0,6,1,6,2,6,3,6,4,
     &	6,5,6,6,6,7,6,8,6,9,6,10,6,11,6,12,
     &	6,13,6,14,6,15,6,16,6,17,6,18,6,19,7,-19,
     &	7,-18,7,-17,7,-16,7,-15,7,-14,7,-13,7,-12,7,-11,
     &	7,-10,7,-9,7,-8,7,-7,7,-6,7,-5,7,-4,7,-3,
     &	7,-2,7,-1,7,0,7,1,7,2,7,3,7,4,7,5,
     &	7,6,7,7,7,8,7,9,7,10,7,11,7,12,7,13,
     &	7,14,7,15,7,16,7,17,7,18,7,19,8,-18,8,-17/
	DATA (aiv(i),i=1921,2080)/
     &	8,-16,8,-15,8,-14,8,-13,8,-12,8,-11,8,-10,8,-9,
     &	8,-8,8,-7,8,-6,8,-5,8,-4,8,-3,8,-2,8,-1,
     &	8,0,8,1,8,2,8,3,8,4,8,5,8,6,8,7,
     &	8,8,8,9,8,10,8,11,8,12,8,13,8,14,8,15,
     &	8,16,8,17,8,18,9,-18,9,-17,9,-16,9,-15,9,-14,
     &	9,-13,9,-12,9,-11,9,-10,9,-9,9,-8,9,-7,9,-6,
     &	9,-5,9,-4,9,-3,9,-2,9,-1,9,0,9,1,9,2,
     &	9,3,9,4,9,5,9,6,9,7,9,8,9,9,9,10,
     &	9,11,9,12,9,13,9,14,9,15,9,16,9,17,9,18,
     &	10,-17,10,-16,10,-15,10,-14,10,-13,10,-12,10,-11,10,-10/
	DATA (aiv(i),i=2081,2240)/
     &	10,-9,10,-8,10,-7,10,-6,10,-5,10,-4,10,-3,10,-2,
     &	10,-1,10,0,10,1,10,2,10,3,10,4,10,5,10,6,
     &	10,7,10,8,10,9,10,10,10,11,10,12,10,13,10,14,
     &	10,15,10,16,10,17,11,-17,11,-16,11,-15,11,-14,11,-13,
     &	11,-12,11,-11,11,-10,11,-9,11,-8,11,-7,11,-6,11,-5,
     &	11,-4,11,-3,11,-2,11,-1,11,0,11,1,11,2,11,3,
     &	11,4,11,5,11,6,11,7,11,8,11,9,11,10,11,11,
     &	11,12,11,13,11,14,11,15,11,16,11,17,12,-16,12,-15,
     &	12,-14,12,-13,12,-12,12,-11,12,-10,12,-9,12,-8,12,-7,
     &	12,-6,12,-5,12,-4,12,-3,12,-2,12,-1,12,0,12,1/
	DATA (aiv(i),i=2241,2400)/
     &	12,2,12,3,12,4,12,5,12,6,12,7,12,8,12,9,
     &	12,10,12,11,12,12,12,13,12,14,12,15,12,16,13,-15,
     &	13,-14,13,-13,13,-12,13,-11,13,-10,13,-9,13,-8,13,-7,
     &	13,-6,13,-5,13,-4,13,-3,13,-2,13,-1,13,0,13,1,
     &	13,2,13,3,13,4,13,5,13,6,13,7,13,8,13,9,
     &	13,10,13,11,13,12,13,13,13,14,13,15,14,-14,14,-13,
     &	14,-12,14,-11,14,-10,14,-9,14,-8,14,-7,14,-6,14,-5,
     &	14,-4,14,-3,14,-2,14,-1,14,0,14,1,14,2,14,3,
     &	14,4,14,5,14,6,14,7,14,8,14,9,14,10,14,11,
     &	14,12,14,13,14,14,15,-13,15,-12,15,-11,15,-10,15,-9/
	DATA (aiv(i),i=2401,2560)/
     &	15,-8,15,-7,15,-6,15,-5,15,-4,15,-3,15,-2,15,-1,
     &	15,0,15,1,15,2,15,3,15,4,15,5,15,6,15,7,
     &	15,8,15,9,15,10,15,11,15,12,15,13,16,-12,16,-11,
     &	16,-10,16,-9,16,-8,16,-7,16,-6,16,-5,16,-4,16,-3,
     &	16,-2,16,-1,16,0,16,1,16,2,16,3,16,4,16,5,
     &	16,6,16,7,16,8,16,9,16,10,16,11,16,12,17,-11,
     &	17,-10,17,-9,17,-8,17,-7,17,-6,17,-5,17,-4,17,-3,
     &	17,-2,17,-1,17,0,17,1,17,2,17,3,17,4,17,5,
     &	17,6,17,7,17,8,17,9,17,10,17,11,18,-9,18,-8,
     &	18,-7,18,-6,18,-5,18,-4,18,-3,18,-2,18,-1,18,0/
	DATA (aiv(i),i=2561,2626)/
     &	18,1,18,2,18,3,18,4,18,5,18,6,18,7,18,8,
     &	18,9,19,-7,19,-6,19,-5,19,-4,19,-3,19,-2,19,-1,
     &	19,0,19,1,19,2,19,3,19,4,19,5,19,6,19,7,
     &	20,-4,20,-3,20,-2,20,-1,20,0,20,1,20,2,20,3,
     &	20,4/
	DATA (amv(i),i=1,180)/
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,473,514,555,596,637,678,719,760,
     &	801,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,356,395,434,474,515,556,
     &	597,638,679,720,761,802,842,881,920,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,282,319,357,
     &	396,435,475,516,557,598,639,680,721,762,803,843,
     &	882,921,959,996,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	212,247,283,320,358,397,436,476,517,558,599,640,
     &	681,722,763,804,844,883,922,960,997,1033,1068,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,179,213,248,284,321,359,398,437/
	DATA (amv(i),i=181, 360)/
     &	477,518,559,600,641,682,723,764,805,845,884,923,
     &	961,998,1034,1069,1103,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,148,180,214,249,
     &	285,322,360,399,438,478,519,560,601,642,683,724,
     &	765,806,846,885,924,962,999,1035,1070,1104,1136,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	119,149,181,215,250,286,323,361,400,439,479,520,
     &	561,602,643,684,725,766,807,847,886,925,963,1000,
     &	1036,1071,1105,1137,1167,0,0,0,0,0,0,0,
     &	0,0,0,0,92,120,150,182,216,251,287,324,
     &	362,401,440,480,521,562,603,644,685,726,767,808,
     &	848,887,926,964,1001,1037,1072,1106,1138,1168,1196,0,
     &	0,0,0,0,0,0,0,0,67,93,121,151,
     &	183,217,252,288,325,363,402,441,481,522,563,604,
     &	645,686,727,768,809,849,888,927,965,1002,1038,1073/
	DATA (amv(i),i=361, 540)/
     &	1107,1139,1169,1197,1223,0,0,0,0,0,0,0,
     &	44,68,94,122,152,184,218,253,289,326,364,403,
     &	442,482,523,564,605,646,687,728,769,810,850,889,
     &	928,966,1003,1039,1074,1108,1140,1170,1198,1224,1248,0,
     &	0,0,0,0,0,45,69,95,123,153,185,219,
     &	254,290,327,365,404,443,483,524,565,606,647,688,
     &	729,770,811,851,890,929,967,1004,1040,1075,1109,1141,
     &	1171,1199,1225,1249,0,0,0,0,0,25,46,70,
     &	96,124,154,186,220,255,291,328,366,405,444,484,
     &	525,566,607,648,689,730,771,812,852,891,930,968,
     &	1005,1041,1076,1110,1142,1172,1200,1226,1250,1271,0,0,
     &	0,0,26,47,71,97,125,155,187,221,256,292,
     &	329,367,406,445,485,526,567,608,649,690,731,772,
     &	813,853,892,931,969,1006,1042,1077,1111,1143,1173,1201,
     &	1227,1251,1272,0,0,0,10,27,48,72,98,126/
	DATA (amv(i),i=541, 720)/
     &	156,188,222,257,293,330,368,407,446,486,527,568,
     &	609,650,691,732,773,814,854,893,932,970,1007,1043,
     &	1078,1112,1144,1174,1202,1228,1252,1273,1290,0,0,11,
     &	28,49,73,99,127,157,189,223,258,294,331,369,
     &	408,447,487,528,569,610,651,692,733,774,815,855,
     &	894,933,971,1008,1044,1079,1113,1145,1175,1203,1229,1253,
     &	1274,1291,0,0,12,29,50,74,100,128,158,190,
     &	224,259,295,332,370,409,448,488,529,570,611,652,
     &	693,734,775,816,856,895,934,972,1009,1045,1080,1114,
     &	1146,1176,1204,1230,1254,1275,1292,0,1,13,30,51,
     &	75,101,129,159,191,225,260,296,333,371,410,449,
     &	489,530,571,612,653,694,735,776,817,857,896,935,
     &	973,1010,1046,1081,1115,1147,1177,1205,1231,1255,1276,1293,
     &	1305,2,14,31,52,76,102,130,160,192,226,261,
     &	297,334,372,411,450,490,531,572,613,654,695,736/
	DATA (amv(i),i=721, 900)/
     &	777,818,858,897,936,974,1011,1047,1082,1116,1148,1178,
     &	1206,1232,1256,1277,1294,1306,3,15,32,53,77,103,
     &	131,161,193,227,262,298,335,373,412,451,491,532,
     &	573,614,655,696,737,778,819,859,898,937,975,1012,
     &	1048,1083,1117,1149,1179,1207,1233,1257,1278,1295,1307,4,
     &	16,33,54,78,104,132,162,194,228,263,299,336,
     &	374,413,452,492,533,574,615,656,697,738,779,820,
     &	860,899,938,976,1013,1049,1084,1118,1150,1180,1208,1234,
     &	1258,1279,1296,1308,5,17,34,55,79,105,133,163,
     &	195,229,264,300,337,375,414,453,493,534,575,616,
     &	657,698,739,780,821,861,900,939,977,1014,1050,1085,
     &	1119,1151,1181,1209,1235,1259,1280,1297,1309,6,18,35,
     &	56,80,106,134,164,196,230,265,301,338,376,415,
     &	454,494,535,576,617,658,699,740,781,822,862,901,
     &	940,978,1015,1051,1086,1120,1152,1182,1210,1236,1260,1281/
	DATA (amv(i),i=901,1080)/
     &	1298,1310,7,19,36,57,81,107,135,165,197,231,
     &	266,302,339,377,416,455,495,536,577,618,659,700,
     &	741,782,823,863,902,941,979,1016,1052,1087,1121,1153,
     &	1183,1211,1237,1261,1282,1299,1311,8,20,37,58,82,
     &	108,136,166,198,232,267,303,340,378,417,456,496,
     &	537,578,619,660,701,742,783,824,864,903,942,980,
     &	1017,1053,1088,1122,1154,1184,1212,1238,1262,1283,1300,1312,
     &	9,21,38,59,83,109,137,167,199,233,268,304,
     &	341,379,418,457,497,538,579,620,661,702,743,784,
     &	825,865,904,943,981,1018,1054,1089,1123,1155,1185,1213,
     &	1239,1263,1284,1301,1313,0,22,39,60,84,110,138,
     &	168,200,234,269,305,342,380,419,458,498,539,580,
     &	621,662,703,744,785,826,866,905,944,982,1019,1055,
     &	1090,1124,1156,1186,1214,1240,1264,1285,1302,0,0,23,
     &	40,61,85,111,139,169,201,235,270,306,343,381/
	DATA (amv(i),i=1081,1260)/
     &	420,459,499,540,581,622,663,704,745,786,827,867,
     &	906,945,983,1020,1056,1091,1125,1157,1187,1215,1241,1265,
     &	1286,1303,0,0,24,41,62,86,112,140,170,202,
     &	236,271,307,344,382,421,460,500,541,582,623,664,
     &	705,746,787,828,868,907,946,984,1021,1057,1092,1126,
     &	1158,1188,1216,1242,1266,1287,1304,0,0,0,42,63,
     &	87,113,141,171,203,237,272,308,345,383,422,461,
     &	501,542,583,624,665,706,747,788,829,869,908,947,
     &	985,1022,1058,1093,1127,1159,1189,1217,1243,1267,1288,0,
     &	0,0,0,43,64,88,114,142,172,204,238,273,
     &	309,346,384,423,462,502,543,584,625,666,707,748,
     &	789,830,870,909,948,986,1023,1059,1094,1128,1160,1190,
     &	1218,1244,1268,1289,0,0,0,0,0,65,89,115,
     &	143,173,205,239,274,310,347,385,424,463,503,544,
     &	585,626,667,708,749,790,831,871,910,949,987,1024/
	DATA (amv(i),i=1261,1440)/
     &	1060,1095,1129,1161,1191,1219,1245,1269,0,0,0,0,
     &	0,0,66,90,116,144,174,206,240,275,311,348,
     &	386,425,464,504,545,586,627,668,709,750,791,832,
     &	872,911,950,988,1025,1061,1096,1130,1162,1192,1220,1246,
     &	1270,0,0,0,0,0,0,0,91,117,145,175,
     &	207,241,276,312,349,387,426,465,505,546,587,628,
     &	669,710,751,792,833,873,912,951,989,1026,1062,1097,
     &	1131,1163,1193,1221,1247,0,0,0,0,0,0,0,
     &	0,0,118,146,176,208,242,277,313,350,388,427,
     &	466,506,547,588,629,670,711,752,793,834,874,913,
     &	952,990,1027,1063,1098,1132,1164,1194,1222,0,0,0,
     &	0,0,0,0,0,0,0,0,147,177,209,243,
     &	278,314,351,389,428,467,507,548,589,630,671,712,
     &	753,794,835,875,914,953,991,1028,1064,1099,1133,1165,
     &	1195,0,0,0,0,0,0,0,0,0,0,0/
	DATA (amv(i),i=1441,1620)/
     &	0,0,178,210,244,279,315,352,390,429,468,508,
     &	549,590,631,672,713,754,795,836,876,915,954,992,
     &	1029,1065,1100,1134,1166,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,211,245,280,316,
     &	353,391,430,469,509,550,591,632,673,714,755,796,
     &	837,877,916,955,993,1030,1066,1101,1135,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,246,281,317,354,392,431,470,510,551,592,
     &	633,674,715,756,797,838,878,917,956,994,1031,1067,
     &	1102,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,318,355,393,
     &	432,471,511,552,593,634,675,716,757,798,839,879,
     &	918,957,995,1032,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,394,433,472,512,553,594,635,676/
	DATA (amv(i),i=1621,1681)/
     &	717,758,799,840,880,919,958,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	513,554,595,636,677,718,759,800,841,0,0,0,
     &	0,0,0,0,0,0,0,0,0,0,0,0,
     &	0/

	END BLOCK DATA

C***********************************************************************
C***********************************************************************

	FUNCTION ISQARRAYDM(SQAFile,SQAna,coef)

	IMPLICIT NONE

	INTEGER i,j,k
	CHARACTER(len=*) :: SQAFile

	INTEGER mSQAna,mSQActr,SQAna,SQActr,mSQA2,mlo,mhi
	PARAMETER (mSQAna=63,mSQActr=32,mlo=-31,mhi=31,
     &	mSQA2=mSQAna*mSQAna)
	INTEGER actMap(mlo:mhi,mlo:mhi),actIndex(2,mSQA2),Nact,ilo,ihi
	COMMON /SQARRAYDM1/ SQActr,actMap,actIndex,Nact,ilo,ihi
	SAVE /SQARRAYDM1/

	LOGICAL exist,ISQARRAYDM
	INTEGER iDat,ICLEN
	REAL*8 coef(*),S1,S2,x,y,r2,di

  103	FORMAT(' Actuator array is ',i4,' by',i4,
     &	', with ',i5,' actuators')

C   Check that DM doesn't have too many actuators
	IF (SQAna.GT.mSQAna) THEN
	  WRITE(*,*)' WARNING: Too many actuators for SQARRAY!'
	  WRITE(*,*)' SQAna reset to mSQAna'
	  SQAna=mSQAna
	END IF

C   Set up maps
	ilo=-(SQAna-1)/2
	ihi=ilo+SQAna-1
	di=0.5d0*DBLE(ilo+ihi)
	r2=(SQAna/2)**2
	k=0
	DO j=mlo,mhi
	  DO i=mlo,mhi
	    k=k+1
	    actMap(i,j)=0
	    actIndex(1,k)=0
	    actIndex(2,k)=0
	  END DO
	END DO
	k=0
	DO j=ilo,ihi
	  y=DBLE(j)-di
	  DO i=ilo,ihi
	    x=DBLE(i)-di
	    IF ((x*x+y*y).LE.r2) THEN
	      k=k+1
	      actMap(i,j)=k
	      actIndex(1,k)=i
	      actIndex(2,k)=j
	    END IF
	  END DO
	END DO
	Nact=k

C   Check file existence
	IF (SQAFile.EQ.' ') SQAFile='SQARRAYDM.txt'
	INQUIRE (FILE=SQAfile,EXIST=exist)
	IF (exist) THEN
	  WRITE(*,*)' Opening file ',SQAFile(1:ICLEN(SQAFile))
          OPEN (UNIT=94,FILE=SQAfile,STATUS='OLD',FORM='FORMATTED')
	ELSE
	  WRITE(*,*)' File ',SQAFile(1:ICLEN(SQAFile)),
     &	  ' does not exist (func ISQARRAYDM in usersub.f).'
	  ISQARRAYDM=.FALSE.
	  RETURN
	END IF

C  Read input data
        iDat=0
	DO i=1,Nact
	  iDat=iDat+1
	  READ(94,*,ERR=2,END=1) S1, S2
	  j=S1
	  coef(j)=S2
	END DO

  1	CONTINUE
	ISQARRAYDM=.TRUE.
	WRITE(*,103) SQAna,SQAna,Nact
	CLOSE(94)	
	RETURN

  2	CONTINUE
	ISQARRAYDM=.FALSE.
	WRITE(*,*) ' Error during read: ISQARRAYDM failed!!!'
	CLOSE(94)
	RETURN

	END FUNCTION ISQARRAYDM

C***********************************************************************

	FUNCTION SSQARRAY(Kc,Kr,psi,pv,prel,ihat,beta,dAct,coef,
     &	rad,pmon,rhoc,xhat,yhat,zhat,L)

C  Computes surface function for square array DM

	IMPLICIT NONE
	INTEGER i,j,nearAct,iRho,jRho,imin,imax,jmin,jmax,idLim,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,SSQARRAY,rho(3),rhom(3),rhoc(3),
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),zc2,psitzc,czc

	INTEGER mSQAna,mSQActr,SQAna,SQActr,mSQA2,mlo,mhi
	PARAMETER (mSQAna=63,mSQActr=32,mlo=-31,mhi=31,
     &	mSQA2=mSQAna*mSQAna)
	INTEGER actMap(mlo:mhi,mlo:mhi),actIndex(2,mSQA2),Nact,ilo,ihi
	COMMON /SQARRAYDM1/ SQActr,actMap,actIndex,Nact,ilo,ihi
	SAVE /SQARRAYDM1/

	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/5/
	DATA pio4/0.78539816339745d0/

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add influence functions of neighboring actuators

	imin=MAX((iRho-idLim), ilo)
	imax=MIN((iRho+idLim), ihi)
	jmin=MAX((jRho-idLim), ilo)
	jmax=MIN((jRho+idLim), ihi)
	fh=0d0
	DO j=jmin,jmax
   	  DO i=imin,imax
   	    nearAct=actMap(i,j)
   	    IF (nearAct.GT.0) THEN
   	      IF (coef(nearAct).NE.0d0) THEN
   	        xAct=dAct*actIndex(1,nearAct)
   	        yAct=dAct*actIndex(2,nearAct)
   	        bd=beta*SQRT((x-xAct)**2+(y-yAct)**2)
   	        IF (bd.LT.czc) 
     &	          fh=fh+coef(nearAct)*EXP(-bd)*SIN(bd+pio4)
   	      END IF
   	    END IF
   	  END DO
   	END DO

C  Compute surface function

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO
	zc2=zc(1)*zc(1)+zc(2)*zc(2)+zc(3)*zc(3)
	psitzc=psi(1)*zc(1)+psi(2)*zc(2)+psi(3)*zc(3)
	SSQARRAY=zc2+(Kc*psitzc+2d0*Kr)*psitzc

	RETURN
	END FUNCTION SSQARRAY

C***********************************************************************

	SUBROUTINE NSQARRAY(ifBuild,Kc,Kr,psi,pv,prel,ihat,
     &	  beta,dAct,coef,rad,pmon,rhoc,xhat,yhat,zhat,L,zc,fh,
     &	  dFdrho,d2Fdrho2)
	use math_mod

C  Computes surface function for 349 actuator, square array DM

	IMPLICIT NONE
	LOGICAL ifBuild
	INTEGER i,j,nearAct,iRho,jRho,imin,imax,jmin,jmax,idLim,
     &	iAxis,jAxis,nzc
	REAL*8 Kc,Kr,psi(3),pv(3),prel(3),ihat(3),coef(*),xhat(3),
     &	yhat(3),zhat(3),rad,pmon(3),L,rho(3),rhom(3),rhoc(3),rhoAx,
     &	fh,beta,dAct,xAct,yAct,x,y,bd,pio4,zc(3),
     &	dFdrho(3),d2Fdrho2(3,3),esbd,ecbd,ebd,beta2,d,bod,dVec(3),
     &	CdF,CdF2,dx,dy,czc

	INTEGER mSQAna,mSQActr,SQAna,SQActr,mSQA2,mlo,mhi
	PARAMETER (mSQAna=63,mSQActr=32,mlo=-31,mhi=31,
     &	mSQA2=mSQAna*mSQAna)
	INTEGER actMap(mlo:mhi,mlo:mhi),actIndex(2,mSQA2),Nact,ilo,ihi
	COMMON /SQARRAYDM1/ SQActr,actMap,actIndex,Nact,ilo,ihi
	SAVE /SQARRAYDM1/

	DATA idLim/3/
	DATA czc/5.4977871437821d0/
	DATA nzc/10/
	DATA pio4/0.78539816339745d0/

	beta2=beta*beta

C  Compute point-of-incidence in surface coords

	DO 1 i=1,3
	  rho(i)=prel(i)+L*ihat(i)
	  rhom(i)=rho(i)-rhoc(i)
  1	CONTINUE
	x=(xhat(1)*rhom(1)+xhat(2)*rhom(2)+xhat(3)*rhom(3))
	y=(yhat(1)*rhom(1)+yhat(2)*rhom(2)+yhat(3)*rhom(3))

C  Find actuator coordinates of ray incidence point

	iRho=IDNINT(x/dAct)
	jRho=IDNINT(y/dAct)

C  Add neighboring actuator influence functions and derivatives

	imin=MAX((iRho-idLim), ilo)
	imax=MIN((iRho+idLim), ihi)
	jmin=MAX((jRho-idLim), ilo)
	jmax=MIN((jRho+idLim), ihi)
	fh=0d0
	IF (ifBuild) THEN
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
  	    DO iAxis=1,3
	      d2Fdrho2(iAxis,jAxis)=0d0
   	    END DO
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          d=SQRT((x-xAct)**2+(y-yAct)**2)
   	          bd=beta*d
   	          IF (bd.LT.czc) THEN
   	            DO iAxis=1,3
   	              rhoAx=xAct*xhat(iAxis)+
     &                yAct*yhat(iAxis)-rhoc(iAxis)
c 	              dVec(iAxis)=rhom(iAxis)-rhoAx
  	              dVec(iAxis)=(1d0-zhat(iAxis))*rhom(iAxis)-rhoAx
   	            END DO
   	            IF (d.lt.1d-12) THEN
   	              ebd=coef(nearAct)
    	              esbd=ebd*0.70710678118655d0
    	              ecbd=ebd*0.70710678118655d0
   	              fh=fh+esbd
     	              CdF2=-2d0*beta*beta*ecbd
 	              DO iAxis=1,3
   	                d2Fdrho2(iAxis,iAxis)=d2Fdrho2(iAxis,iAxis)+CdF2
   	              END DO
    	            ELSE
   	              bod=beta/d
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
   	              fh=fh+esbd
   	              CdF=bod*(ecbd-esbd)
   	              CdF2=-2d0*bod*bod*ecbd
   	              DO jAxis=1,3
   	                dFdrho(jAxis)=dFdrho(jAxis)+CdF*dVec(jAxis)
   	                d2Fdrho2(jAxis,jAxis)=d2Fdrho2(jAxis,jAxis)+CdF
   	                DO iAxis=1,3
   	                  d2Fdrho2(iAxis,jAxis)=d2Fdrho2(iAxis,jAxis)+
     &   	          (CdF2-CdF/d/d)*dVec(iAxis)*dVec(jAxis)
   	                END DO
   	              END DO
   	            END IF
   	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	ELSE
	  DO jAxis=1,3
 	    dFdrho(jAxis)=0d0
   	  END DO
	  DO j=jmin,jmax
   	    DO i=imin,imax
   	      nearAct=actMap(i,j)
   	      IF (nearAct.GT.0) THEN
   	        IF (coef(nearAct).NE.0d0) THEN
   	          xAct=dAct*actIndex(1,nearAct)
   	          yAct=dAct*actIndex(2,nearAct)
   	          dx=x-xAct
   	          dy=y-yAct
   	          d=SQRT(dx*dx+dy*dy)
   	          IF (d.lt.1d-12) THEN
   	            ebd=coef(nearAct)
    	            esbd=ebd*0.70710678118655d0
  	            fh=fh+esbd
   	          ELSE
   	            bd=beta*d
    	            IF (bd.LT.czc) THEN
  	              bod=beta/d
   	              ebd=coef(nearAct)*EXP(-bd)
    	              esbd=ebd*SIN(bd+pio4)
    	              ecbd=ebd*COS(bd+pio4)
  	              fh=fh+esbd
	              CdF=bod*(ecbd-esbd)
 	              DO jAxis=1,3
 	                dFdrho(jAxis)=dFdrho(jAxis)
     & 	                +CdF*(dx*xhat(jAxis)+dy*yhat(jAxis))
   	              END DO
  	            END IF
  	          END IF
  	        END IF
  	      END IF
   	    END DO
   	  END DO
   	END IF

	DO i=1,3
	  zc(i)=rho(i)-fh*zhat(i)
   	END DO

	RETURN
	END SUBROUTINE NSQARRAY

