C***********************************************************************
C       Begin file srtrace.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C 
C  Subroutine SRTRACE traces a single ray through a set of optical surfaces, 
C  either in forward or backward order.
C

C ***** For references only *******************************************
C	SUBROUTINE CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
C    &    ifLNegOK,Cmatrix,iStartElt,iEndElt,iStartRay,iEndRay,
C    &    RayMap,SegMap,ifRayPrint,ifIntSrfSet,L1,nBadRays,
C    &    StartIndRef,StartExtinc,DrawRay,DrawElt,nDrawElt,xDraw,yDraw)
C *********************************************************************

#include "realtype.h"

	Subroutine SRTRACE_Test(iStartElt,iEndElt,RayPos_in,RayDir_in)
	use param_mod
        use elt_mod
	use elemsub_mod
        use traceutil_mod
	use tracesub_mod  ! for accessing getIndex() function

        use macos_debug

        IMPLICIT NONE

	! Subroutine arguments
	Integer :: iStartElt,iEndElt
	Real*8 :: RayPos_in(3), RayDir_in(3)

	! Local variables
	Logical :: ifRayElt(0:mElt),ifIntSrfOK(mElt),if,ifDraw,
     &             ifTransmit,ifLNegOK,ifLNsrf,LRTdummy,ifRayToSeg,
     &             ifUseBaseSrf,ifReturn,L1(1),LRayTrans

	Integer :: i,j,iElt,iRay,iNode,nEm1,nEm2,RayIndex(2,mRay),
     &    PrevNonSeg,NextNonSeg,itst,imin,iFirstNonSeq,iLastNonSeq,
     &    iStartRay,iEndRay,IERROR,nDrawElt,
     &    RayMap(mpts,mpts),SegMap(mpts,mpts),DrawElt(mDrawElt)

	Integer :: NextElt(0:mElt),PrevElt(0:mElt),iminPrev,
     &    jGridSrf(mElt),
     &    iSecondElt,AdjRay(4),LastElt,iEltLoop,NextLast,
     &    iIntSrf,NDP3,jElt,iNextSeq,nBadRays,iNextElt,CrntNSCnt,
     &    iStep,ii

	SREAL Cmatrix(7,mCm,bRay),DrawRay(2,mDrawElt),
     &         dxndxi(7,7,0:mElt),dxndui(7,6,mElt),spot(3)

	REAL*8 :: Lmin,na,nb,RayPosTst(3),RayLTst,PrevDir(3),PrevPos(3),
     &    xa,xb,CurIndRef,PrevIndRef,NextIndRef,StartIndRef,
     &    Intxyz,CurExtinc,PrevExtinc,NextExtinc,StartExtinc,xDraw(3),
     &    yDraw(3),CWU

	! These were arguments to CTRACE.
	! For now, these are "test" parameters
	Logical,parameter :: ifPol=.false.,ifRayPrint=.true.,
     &      ifBuild=.false.,ifDMBuild=.false.
      	Logical :: ifIntSrfSet(1:mElt)

	DATA ifUseBaseSrf/.FALSE./

	tAsphCoef=maxAsphCoef  ! added by jzlou

        print*,'iStartElt,iEndElt=',iStartElt,iEndElt

        ifIntSrfSet(1:mElt)=.false.
C
	! Define optical element indexes in the order to be traced 
	if (iStartElt<iEndElt) then
	  iStep=1
	else
	  iStep=-1
	end if
	iSecondElt=iStartElt+iStep

	! Set RayPos and RayDir iStartElt
	RayPos(1:3,iStartElt)=RayPos_in
        RayDir(1:3,iStartElt)=RayDir_in

C       Find start and end points for segmented elements
        PrevNonSeg=0
        iElt=1
        DO WHILE (iElt.LE.nElt)
          IF (EltID(iElt).NE.11) THEN
            NextElt(PrevNonSeg)=iElt
            PrevElt(iElt)=PrevNonSeg
            PrevNonSeg=iElt
          ELSE
            PrevElt(iElt)=PrevNonSeg
            NextNonSeg=iElt+1
            DO WHILE (EltID(NextNonSeg).EQ.11)
              NextNonSeg=NextNonSeg+1
            END DO
            NextElt(iElt)=NextNonSeg
          END IF
          iElt=iElt+1
        END DO
! Now, for iElt, PrevElt(iElt) gives index of previous non-segment element,
! NextElt(iElt) gives index of next non-segment element.

	! Trace single ray only, setting ray id = 1
	iRay=1  

	! Initialize Index of Reflection etc.
        StartIndRef=IndRef(iStartElt)
        StartExtinc=Extinc(iStartElt)
        CumRayL(iRay)=0d0
        LRayOK(iRay)=.TRUE.
        PrevNonSeg=iStartElt


C       Set wavelength in base units; set index of refraction for Glass;
C       convert Zernike coefficients to monomials; initialize LEltOK;
C       check interpolated surfaces; compute obscuration coords

        CWU=CWM*1e6
        CWB=CWM/CBM
        WaveBU=CWB*Wavelen

	iFirstNonSeq=iStartElt
        iLastNonSeq=iEndElt

        DO 17 iElt=iSecondElt,iEndElt,iStep
C Set index
          IF (LGlass(iElt)) THEN
            IndRef(iElt)=getIndex(CWU*Wavelen,GlassCoef(1,iElt),
     &                            GlassType(iElt))
          END IF
C Set Zerns
          IF (SrfType(iElt).EQ.8) THEN
            IF (ZernTypeL(iElt).EQ.1 .OR.
     &          ZernTypeL(iElt).EQ.4) THEN
C Malacara's
              CALL ZerntoMon1(ZernCoef(1,iElt),MonCoef(1,iElt))
            ELSEIF (ZernTypeL(iElt).EQ.2 .OR.
     &              ZernTypeL(iElt).EQ.5) THEN
C Noll's
              CALL ZerntoMon2(ZernCoef(1,iElt),MonCoef(1,iElt))
            ELSEIF (ZernTypeL(iElt).EQ.3 .OR.
     &              ZernTypeL(iElt).EQ.6) THEN
C Fringe's
              CALL ZerntoMon3(ZernCoef(1,iElt),MonCoef(1,iElt))
            END IF
          END IF

	  iEltToiWF(iElt)=0
          LEltOK(iElt)=.TRUE.
          IF (SrfType(iElt).EQ.5) THEN
            IF (ifIntSrfSet(iElt)) THEN
              iIntSrf=iEltToIntSrf(iElt)
              ifIntSrfOK(iElt)=.TRUE.
            ELSE
              SrfType(iElt)=2
              iIntSrf=1
              ifIntSrfOK(iElt)=.FALSE.
            END IF
          ELSE
            iIntSrf=1
            ifIntSrfOK(iElt)=.TRUE.
          END IF
          IF ((EltID(iElt).EQ.7).OR.(nObs(iElt).GT.0).OR.
     &        (ApType(iElt).GT.0)) THEN
            LEltObs(iElt)=.TRUE.

            CALL DEQUATE(zObs(1,iElt),psiElt(1,iElt),3)
            CALL DXPROD(yObs(1,iElt),zObs(1,iElt),xObs(1,iElt))
            CALL DUNITIZE(yObs(1,iElt))
            CALL DXPROD(xObs(1,iElt),yObs(1,iElt),zObs(1,iElt))
          ELSE
            LEltObs(iElt)=.FALSE.
          END IF
          IF (SrfType(iElt).EQ.9) THEN
            jGridSrf(iElt)=iEltToGridSrf(iElt)
          ELSE
            jGridSrf(iElt)=1
          END IF
 17     CONTINUE


C  Begin iRay loop to trace a single ray and compute its partials
	
C       Loop through each element

        iElt=iSecondElt
        CurIndRef=StartIndRef
        PrevIndRef=CurIndRef
        CurExtinc=StartExtinc
        PrevExtinc=CurExtinc
        ifReturn=.FALSE.

	DO 10 iEltLoop=iSecondElt,iEndElt,iStep
          IF (iElt.GT.iEndElt) GO TO 10
          IF (iElt.LT.nElt) THEN
            iNextElt=iElt+1
          ELSE
            iNextElt=nElt
          END IF

	  IF ((SrfType(iElt).EQ.5).AND.ifIntSrfOK(iElt)) THEN
            iIntSrf=iEltToIntSrf(iElt)
          END IF

          IF (ifRayPrint) THEN
            CALL DEQUATE(PrevDir,RayDir(1,iRay),3)
            CALL DEQUATE(PrevPos,RayPos(1,iRay),3)
          END IF

	  IF (ifLNegOK) THEN
            IF ((EltID(iElt).EQ.3).OR.(EltID(iElt).EQ.7).OR.
     &        (EltID(iElt).EQ.8).OR.(EltID(PrevNonSeg).EQ.3).OR.
     &        (EltID(PrevNonSeg).EQ.7).OR.(EltID(PrevNonSeg).EQ.8))
     &      THEN
              ifLNsrf=.TRUE.
            ELSE
              ifLNsrf=.FALSE.
            END IF
          END IF


          crnt_Elt=iElt  ! crnt_Elt is shared with
                         ! ChkRayTrans function in elemsub.F

C       Trace rays using appropriate element routines
C                                               ! Reflector
          IF (EltID(iElt).EQ.1) THEN
	      print*,' --> Trace through a Reflector'
              PrevNonSeg=iElt

              CALL Reflector(
     &          ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &          LEltObs(iElt),LRayOK(iRay),LRayTrans,KcElt(iElt),
     &          KrElt(iElt),AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &          AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &          UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &          GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &          psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &          RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &          lMon(iElt),pMon(1,iElt),
     &          xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &          CurIndRef,IndRef(iElt),
     &          CurExtinc,Extinc(iElt),WaveBU,
     &          ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &          ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &          RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &          XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &          IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &	        EltCoat(iElt),EltCoatThk(:,iElt),
     &          IndRefArr(:,iElt),ExtincArr(:,iElt))

	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
              CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)

              PrevIndRef=CurIndRef
              PrevExtinc=CurExtinc
              ifRayElt(iElt)=.TRUE.

C                                               ! Focal Plane
          ELSE IF (EltID(iElt).EQ.2) THEN
	      print*,' --> Trace through a Focal Plane'
              PrevNonSeg=iElt

              ifReturn=.FALSE.

              CALL FocalPln(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &          LRayOK(iRay),LRayTrans,psiElt(1,iElt),
     &          VptElt(1,iElt),RptElt(1,iElt),RayPos(1,iRay),
     &          RayDir(1,iRay),RayE(1,iRay),CurIndRef,
     &          IndRef(iElt),CurExtinc,Extinc(iElt),WaveBU,
     &          ApType(iElt),ApVec(1,iElt),nObs(iElt),
     &          ObsType(1,iElt),ObsVec(1,1,iElt),
     &          xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &          RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),
     &          RayE(1,iRay),spot,
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt))

              IF (ifReturn) RayL(iRay)=-RayL(iRay)
              IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
              CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
              PrevIndRef=CurIndRef
              PrevExtinc=CurExtinc
              ifRayElt(iElt)=.TRUE.
	  !
C                                               ! Reference Surface 
            ELSE IF (EltID(iElt).EQ.3) THEN
	      print*,' --> Trace through a Reference Surface'
              PrevNonSeg=iElt

              CALL RefSrf(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &          LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &          psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &          RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &          CurIndRef,IndRef(iElt),
     &          CurExtinc,Extinc(iElt),WaveBU,
     &          ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &          ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),
     &          zObs(1,iElt),
     &          RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt))

	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
              IF (ifReturn) RayL(iRay)=-RayL(iRay)
              CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
              PrevIndRef=CurIndRef
              PrevExtinc=CurExtinc
              ifRayElt(iElt)=.TRUE.
	      !
C                                               ! Reflector Segment
            ELSE IF (EltID(iElt).EQ.11) THEN
	      print*,' --> Trace through a Reflector Segment'
              ifRayToSeg=.FALSE.
              DO i=1,nRayToSeg(iRay)
                IF (RayToSegMap(iRay,i).EQ.EltToSegMap(iElt))
     &            ifRayToSeg=.TRUE.
              END DO
              IF (ifRayToSeg) THEN
c               NextNonSeg=PrevNonSeg+nSeg
c               PrevElt(iElt)=PrevNonSeg
c               NextElt(PrevNonSeg)=iElt
c               NextElt(iElt)=NextNonSeg
                PrevNonSeg=iElt
                NextNonSeg=NextElt(iElt)
                PrevElt(NextNonSeg)=iElt

	        CALL Reflector(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &            LEltObs(iElt),
     &            LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &            AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &            AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &            UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &            GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &            psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &            RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &            lMon(iElt),pMon(1,iElt),
     &            xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &            CurIndRef,IndRef(iElt),
     &            CurExtinc,Extinc(iElt),WaveBU,
     &            ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &            ObsVec(1,1,iElt),
     &            xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &            RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &            dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &            XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &            IWORK(1,iIntSrf),DWORK(1,iIntSrf),
     &            NDP(iIntSrf),mDP,
     &	   	  EltCoat(iElt),EltCoatThk(:,iElt),
     &            IndRefArr(:,iElt),ExtincArr(:,iElt))

	        IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
                CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
                PrevIndRef=CurIndRef
                PrevExtinc=CurExtinc
                ifRayElt(iElt)=.TRUE.
              ELSE
                PrevElt(iElt)=0
                NextElt(iElt)=0
                ifRayElt(iElt)=.FALSE.
              END IF
	        
	  END IF ! Element selection 
	  iElt=iElt+1
 10	CONTINUE 
  
	End Subroutine SRTRACE_Test
