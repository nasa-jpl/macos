C***********************************************************************
C	Begin file propsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2008, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C***********************************************************************
C	Changes pending for propagation:
C	1. Check dx calc and conditionals for sfprop etc.
C	2. Replace (Z2-Z1) by (Z1-Z2) everywhere in NFprop.
C	3. Reinstate S=Z1/Z2 everywhere in NFprop.
C	4. Shift of x and y is different? [-w/2+dx,w/2]
C	5. Change applyfac routines (done)
C	6. Make alpha user definable (dftsub)
C***********************************************************************

#include "realtype.h"

        Module propsub_mod
	use math_mod
        use dft_mod

	PUBLIC :: RunProp, CPROPAGATE, NFPROP, FFPROP, FRPROP

        PRIVATE  ! everythinb below is module-private

	Logical :: ifFnCalc=.FALSE.  ! turned off now

        Contains

C
C       This routine is run when executing 'INtensity' or 'PIx' commands.
C       It calls CPROPAGATE to perform actual propagation. All rays are
C       traced to one element at a time.  -jzlou
C

	SUBROUTINE RunProp(iElt,npts,ifPol,ifGrid,ifPropagate,
     &	  ifBuild,ifTrace,ifRegrid,ifVecDif3,ifIntSrfSet,ifLNegOK,
     &	  iEndElt,iCurWFElt,iCurRayElt,iWF,Dicr,Djcr,dxSource,
     &	  iObsOpt,nBadRays,StartIndRef,StartExtinc,DWF,OPDMat)

        use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
	use traceutil_mod

        use sourcsub_mod
        use tracesub_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifLNegOK,ifBuild,ifTrace,
     &	        ifRegrid(mElt),ifVecDif3,ifIntSrfSet(mElt)
	INTEGER iElt,npts,ICLEN,iCurWFElt,iCurRayElt,iWF,
     &	        iObsOpt,nBadRays,iStartElt,iEndElt
	SREAL StartTime,DTime,EndTime
	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &         OPDMat(mpts,mpts)
	COMPLEX*16 DWF(mdttl,mdttl)

 116    FORMAT(' Compute time was',f10.4,' sec')
 120	FORMAT(' Tracing ',i9,' rays and propagating ',i9,
     &	' grid points...')


	tAsphCoef=maxAsphCoef  ! added by jzlou

	IF ((.NOT.ifPropagate).OR.(iElt.NE.iCurWFElt)) THEN
	  iEndElt=iElt
	  IF (ifPropagate.AND.(iCurWFElt.LE.iEndElt)) THEN
	    iStartElt=iCurWFElt
	  ELSE
	    iStartElt=0
	  END IF
	  ifBuild=.FALSE.
	  ifPropagate=.TRUE.
	  ifTrace=.FALSE.
	  IF ((iStartElt.EQ.0).OR.(iStartElt.LT.iEndElt)) THEN
	    CALL RTIME(StartTime)
	    IF (iStartElt.EQ.0) THEN
	      taperFactor=1.0  ! reset tapering factor -jzlou
	      IF (ABS(zSource).GE.1d10) THEN
	        CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	                       ifPropagate)
	      ELSE
	        CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	                      ifPropagate)
	      END IF
	      ifGrid=.TRUE.
	    END IF
	    WRITE(*,120) nRay, md2
	    ifPropagate=.TRUE.
	    CALL CPROPAGATE(RayIndex,iWF,ifRegrid,ifPol,ifVecDif3,
     &	      iObsOpt,RayMap,SegMap,iCurRayElt,iStartElt,iEndElt,
     &	      ifIntSrfSet,ifLNegOK,BeamType,rxBeam,ryBeam,Dicr,Djcr,
     &	      dxSource,DWF,OPDMat,npts,nBadRays,StartIndRef,StartExtinc)
	      iCurWFElt=iEndElt
	      CALL RTIME(EndTime)
	      DTime=EndTime-StartTime
	      WRITE(*,116)DTime
	  END IF
	  IF (nBadRays.GT.0) CALL WARN(LEltOK,iEndElt,nBadRays)
	END IF
	iWF=iEltToiWF(iElt)

	RETURN
	END SUBROUTINE RunProp

C***********************************************************************

	SUBROUTINE CPROPAGATE(RayIndex,iWF,ifRegrid,ifPol,ifVecDif3,
     &	  iObsOpt,RayMap,SegMap,iCurRayElt,iStartElt,iEndElt,ifIntSrfSet,
     &    ifLNegOK,BeamType,rxBeam,ryBeam,Dicr,Djcr,dxSource,DWF,OPDMat,
     &    npts,nBadRays,StartIndRef,StartExtinc)

	use param_mod
        use dft_mod, only: NFPropDFT, FFPropDFT
        use elt_mod
        use sourcsub_mod
        use surfsub, only: ZerntoMon1, ZerntoMon2, ZerntoMon3,
     &                     ZerntoMon4
        use tracesub_mod
	use elemsub_mod
	use traceutil_mod

	IMPLICIT NONE

	LOGICAL ifBuild,LROKdummy,ifRegrid(mElt),ifAccumOPD,LRayTrans,
     &	 ifOPDCalc,ifIntSrfSet(mElt),ifIntSrfOK(mElt),ifDMBuild,ifPol,
     &	 ifVecDif3,ifReturn,ifTransmit,ifLNegOK,ifLNsrf,LRTdummy,
     &	 ifRayToSeg,ifUseBaseSrf,ifReturnElt(mElt)
cdebug     &	,ifRayPrint

	INTEGER i,j,k,iElt,iSeg,iRay,iEm1,RayIndex(2,mRay),
     &    PrevNonSeg(mRay),NextNonSeg,itst,imin,iFirstNonSeq,
     &	  iLastNonSeq,iWF,nEstart,iLstart,iIntSrf,iminPrev,iObsOpt,
     &	  iLend,RayMap(mpts,mpts),SegMap(mpts,mpts),PrevElt(0:mElt),
     &	  iCurRayElt,iStartElt,iEndElt,iSecondElt,iEltLoop,IERROR,
     &	  BeamType,iNextSeq,npts,nBadRays,iNextElt,NextElt(0:mElt),
     &	  iPrevNonSeg,jGridSrf(mElt),CrntNSCnt,nSides,iObs,nPassRays
cdebug  &, pRay

	SREAL dxdx(7,7),dxdu(7,6),spot(3)

	REAL*8 dxMin(0:mElt),dxMax(0:mElt),dxAvg(0:mElt),
     &         dxDev(0:mElt),OPDMat(mpts,mpts)

	! CumLStart needs be saved so when iStartElt /= 0 upon re-entering
	! this routine, its previous values can be used. -jzlou
	REAL*8, dimension(:), allocatable, save :: CumLStart

	REAL*8 Lmin,TPL,DPHI,Dicr,Djcr,dxSource,
     &	       Cphi, Lin1,Lin2,RayPosTst(3),RayLTst,dxsf,xa,xb,
     &	       rxBeam,ryBeam,na,nb,CWU,SumL,RefCumRayL,OPDScal,
     &	       CurIndRef,NextIndRef,PrevIndRef,StartIndRef,
     &	       zStart,zEnd,alpha,z1,dx1p,z2,dx2p,dzp,DAvgL,
     &         PNSIndRef,PNSExtinc,CurExtinc,PrevExtinc,
     &	       NextExtinc,StartExtinc,ddumy(2),OPD_in(mdttl,mdttl)
cdebug     &	,z1z2,PrevDir(3),PrevPos(3),Intxyz

	COMPLEX*16 CD,DWF(mdttl,mdttl),Evec(3),C1

	DATA ifUseBaseSrf/.FALSE./

	If (.not.allocated(CumLStart)) allocate(CumLStart(mRay))

! format stmts
#include "propsub.inc"


	ifBuild=.FALSE.; ifDMBuild=.FALSE.
	ifAccumOPD=.FALSE.; ifOPDCalc=.FALSE.
	ifLNsrf=.FALSE.

cdebug
c	ifRayPrint=.TRUE.
c	pRay=1
cdebug


C  1. Initialize raytrace parameters

	EltID(0)=0
	iSecondElt=iStartElt+1

C	Find start and end points for segmented elements

	iSeg=0
	iPrevNonSeg=0
	iElt=1
	DO WHILE (iElt.LE.nElt)
	  IF (EltID(iElt).NE.11) THEN
	    NextElt(iPrevNonSeg)=iElt
	    PrevElt(iElt)=iPrevNonSeg
	    iPrevNonSeg=iElt
	  ELSE
	    PrevElt(iElt)=iPrevNonSeg
	    NextNonSeg=iElt+1
	    DO WHILE (EltID(NextNonSeg).EQ.11)
	      NextNonSeg=NextNonSeg+1
	    END DO
	    NextElt(iElt)=NextNonSeg
	  END IF
	  iElt=iElt+1
	END DO

C	Initialize input beam if first element is source

	IF (iStartElt.EQ.0) THEN
	  StartIndRef=IndRef(0)
	  StartExtinc=Extinc(0)
	  DO 3 iRay=1,nRay
	    LRayOK(iRay)=.TRUE.
	    CumRayL(iRay)=0d0
	    CumLStart(iRay)=0d0
	    LRayPass(iRay)=.TRUE.
 3        CONTINUE
	END IF

C	Set wavelength in base units; set index of refraction for Glass;
C	convert Zernike coefficients to monomials; initialize LEltOK;
C	check interpolated surfaces; compute obscuration coords

	CWU=CWM*1e6
	CWB=CWM/CBM
     	WaveBU=CWB*Wavelen

	DO 16 iElt=iSecondElt,iEndElt
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
	      CALL ZerntoMon1(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)
	    ELSEIF (ZernTypeL(iElt).EQ.2 .OR.
     &              ZernTypeL(iElt).EQ.5) THEN
C Noll's
	      CALL ZerntoMon2(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)
	    ELSEIF (ZernTypeL(iElt).EQ.3 .OR.
     &              ZernTypeL(iElt).EQ.6) THEN
C Fringe's
	      CALL ZerntoMon3(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)
	    ELSEIF (ZernTypeL(iElt).EQ.7) THEN
              CALL ZerntoMon4(ZernCoef(1,iElt),MonCoef(1,iElt))
	    ELSEIF (ZernTypeL(iElt).EQ.8) THEN
              CALL ZerntoMon6(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)
	    ELSEIF (ZernTypeL(iElt).EQ.9) THEN
              CALL ZerntoMon7(ZernCoef(1,iElt),MonCoef(1,iElt),
     &                        ZernAnnuRatio(iElt))
	    END IF
	  END IF
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
	    ifIntSrfOK(iElt)=.TRUE.
	    iIntSrf=1
	  END IF
	  IF (SrfType(iElt).EQ.9) THEN
	    jGridSrf(iElt)=iEltToGridSrf(iElt)
	  ELSE
	    jGridSrf(iElt)=1
	  END IF
	  IF ((EltID(iElt).EQ.7).OR.(nObs(iElt).GT.0).OR.
     &	      (ApType(iElt).GT.0)) THEN
	    ! EltID=7 is of element type "Obscuring"  -jzlou
	    LEltObs(iElt)=.TRUE.
	    CALL DEQUATE(zObs(1,iElt),psiElt(1,iElt),3)
	    CALL DXPROD(yObs(1,iElt),zObs(1,iElt),xObs(1,iElt))
	    CALL DUNITIZE(yObs(1,iElt))
	    CALL DXPROD(xObs(1,iElt),yObs(1,iElt),zObs(1,iElt))
	  ELSE
	    LEltObs(iElt)=.FALSE.
	  END IF
 16	CONTINUE

	DO 43 iRay=1,nRay
	  PrevNonSeg(iRay)=iStartElt
  43	CONTINUE

C  2. Initialize diffraction parameters

	TPL=6.28318530717959D0/WaveBU
	iElt=0
	nEstart=iStartElt
	iLstart=iStartElt
	IF (iStartElt.EQ.0) THEN
	  iWF=1
	  CALL DXCALC(iElt,nEstart,dxMin(iElt),dxMax(iElt),
     &	         dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
	END IF

C      Set ifReturnElt variable for each element

	ifReturn=.FALSE.
	DO iElt=1,nElt
	  ifReturnElt(iElt)=ifReturn
	  IF (EltID(iElt).EQ.8) ifReturn=.NOT.ifReturn
	END DO


C  3. Begin iElt loop to trace beams between each pair of elements

	IF (iElt.LT.nElt) THEN
	  iNextElt=iElt+1
	ELSE
	  iNextElt=nElt
	END IF
	nBadRays=0
	iElt=iSecondElt  ! iElt = next Elt to trace to -jzlou
	iEm1=iElt-1
	CurIndRef=StartIndRef
	PrevIndRef=CurIndRef
	CurExtinc=StartExtinc
	PrevExtinc=CurExtinc

	!
	! iElt is defined as next Elt to trace to in the following
        ! loop. -jzlou
	DO 14 iEltLoop=iSecondElt,iEndElt
	  crnt_Elt=iElt  ! crnt_Elt, defined in traceutil_mod, is shared
                         ! with ChkRayTrans function in elemsub.F

	  iEm1=iElt-1
	  iEltToiWF(iElt)=0
	  IF ((SrfType(iElt).EQ.5).AND.ifIntSrfOK(iElt)) THEN
	    iIntSrf=iEltToIntSrf(iElt)
	  END IF

C	Initialize non-sequential element search parameters

	  PNSIndRef=CurIndRef
	  PNSExtinc=CurExtinc

	  IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.12)) THEN
	    iFirstNonSeq=iElt
	    iNextSeq=iElt
 706	    CONTINUE
	      IF ((EltID(iNextSeq).EQ.9).OR.
     &	      (EltID(iNextSeq).EQ.12)) THEN
	        iNextSeq=iNextSeq+1
	        GO TO 706
	      END IF
	      iLastNonSeq=iNextSeq-1
	  END IF

C	Set conditions for accumulating OPD

	  IF ((PropType(iElt).NE.1).OR.((PropType(iEm1).NE.3)
     &	      .AND.(PropType(iEm1).NE.5).AND.(PropType(iEm1).NE.6)
     &	      .AND.(PropType(iEm1).NE.7).AND.(PropType(iEm1).LT.10)))
     &	  THEN
	    IF ((PropType(iElt).EQ.11).AND.(PropType(iEm1).EQ.10)) THEN
	      ifAccumOPD=.FALSE.
	    ELSE
	      ifAccumOPD=.TRUE.
	    END IF
	  ELSE
	    ifAccumOPD=.FALSE.
	  END IF

	  IF (ifLNegOK) THEN
	    IF ((EltID(iElt).EQ.3).OR.(EltID(iElt).EQ.7).OR.
     &	        (EltID(iElt).EQ.8).OR.(EltID(iEm1).EQ.3).OR.
     &	        (EltID(iElt).EQ.12).OR.
     &	        (EltID(iEm1).EQ.7).OR.(EltID(iEm1).EQ.8)) THEN
	      ifLNsrf=.TRUE.
	    ELSE
	      ifLNsrf=.FALSE.
	    END IF
	  END IF

C  4. Begin iRay loop to trace each ray to the current element

	  DO 9 iRay=1,nRay
	    crnt_traced_ray=iRay

	    IF (LRayOK(iRay)) THEN

C	Trace rays using appropriate element routines
C						! Reflector
	      IF (EltID(iElt).EQ.1) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt

	        CALL Reflector(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	         LEltObs(iElt),
     &	         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	         UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &	         GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &	         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	         lMon(iElt),pMon(1,iElt),
     &	         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	         CurIndRef,IndRef(iElt),
     &	         CurExtinc,Extinc(iElt),WaveBU,
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),XYZSRF(1,1,iIntSrf),
     &	         DZSRF(1,1,iIntSrf),IERROR,IWORK(1,iIntSrf),
     &	         DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &           EltCoat(iElt),EltCoatThk(:,iElt),
     &           IndRefArr(:,iElt),ExtincArr(:,iElt))

	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3).AND.
     &	             (PropType(PrevElt(iElt)).NE.10).AND.
     &	             (PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD) THEN
    	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
    	          IF (ifPol) THEN
    	            RayE(1,iRay)=Evec(1)
    	            RayE(2,iRay)=Evec(2)
    	            RayE(3,iRay)=Evec(3)
    	          END IF
    	        END IF
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Focal Plane
	      ELSE IF (EltID(iElt).EQ.2) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
c	        IF ((EltID(PrevElt(iElt)).EQ.8).AND.
c     &	        (EltID(PrevElt(iElt)-1).NE.8)) THEN
c	          ifReturn=.TRUE.
c	        ELSE
c	          ifReturn=.FALSE.
c	        END IF

	        CALL FocalPln(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &	          LRayOK(iRay),LRayTrans,
     &	          psiElt(1,iElt),VptElt(1,iElt),
     &	          RptElt(1,iElt),RayPos(1,iRay),RayDir(1,iRay),
     &	          RayE(1,iRay),CurIndRef,IndRef(iElt),
     &	          CurExtinc,Extinc(iElt),WaveBU,ApType(iElt),
     &	          ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	          ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),
     &	          zObs(1,iElt),RayPos(1,iRay),RayDir(1,iRay),
     &            RayL(iRay),Evec,spot,dxdx,dxdu)

	        IF (ifReturnElt(iElt)) RayL(iRay)=-RayL(iRay)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3).AND.
     &	          (PropType(PrevElt(iElt)).NE.10).AND.
     &	          (PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD) THEN
	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
    	          IF (ifPol) THEN
    	            RayE(1,iRay)=Evec(1)
    	            RayE(2,iRay)=Evec(2)
    	            RayE(3,iRay)=Evec(3)
    	          END IF
    	        END IF
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Reference Surface
	      ELSE IF (EltID(iElt).EQ.3) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        CALL RefSrf(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &	        LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	        CurIndRef,IndRef(iElt),
     &	        CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &	        dxdx,dxdu)
	        IF (ifReturnElt(iElt)) RayL(iRay)=-RayL(iRay)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3).AND.
     &	          (PropType(PrevElt(iElt)).NE.10).AND.
     &	          (PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD) THEN
	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
    	          IF (ifPol) THEN
    	            RayE(1,iRay)=Evec(1)
    	            RayE(2,iRay)=Evec(2)
    	            RayE(3,iRay)=Evec(3)
    	          END IF
    	        END IF
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Conic Mirror Segment
	      ELSE IF (EltID(iElt).EQ.11) THEN
	        ifRayToSeg=.FALSE.
	        DO i=1,nRayToSeg(iRay)
	          IF (RayToSegMap(iRay,i).EQ.EltToSegMap(iElt))
     &	            ifRayToSeg=.TRUE.
	        END DO
	        IF (ifRayToSeg) THEN
		  NextNonSeg=NextElt(iElt)
c		  PrevElt=PrevNonSeg(iRay)
		  PrevNonSeg(iRay)=iElt
	          CALL Reflector(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	              LEltObs(iElt),
     &	              LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	              AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	              AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	              UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &	              GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &	              psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	              RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	              lMon(iElt),pMon(1,iElt),
     &	              xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	              CurIndRef,IndRef(iElt),
     &	              CurExtinc,Extinc(iElt),WaveBU,
     &	           ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	           ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),
     &	           zObs(1,iElt),
     &	           RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &             dxdx,dxdu,SrfType(iElt),XYZSRF(1,1,iIntSrf),
     &	           DZSRF(1,1,iIntSrf),IERROR,IWORK(1,iIntSrf),
     &	           DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &             EltCoat(iElt),EltCoatThk(:,iElt),
     &             IndRefArr(:,iElt),ExtincArr(:,iElt))
	          IF (.NOT.LRayTrans) THEN
	            LRayPass(iRay)=.FALSE.
	            IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	            .AND.(PropType(PrevElt(iElt)).NE.10)
     &	            .AND.(PropType(PrevElt(iElt)).NE.12))
     &	              WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	          END IF
	          IF (ifAccumOPD) THEN
	            CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
    	            IF (ifPol) THEN
    	              RayE(1,iRay)=Evec(1)
    	              RayE(2,iRay)=Evec(2)
    	              RayE(3,iRay)=Evec(3)
    	            END IF
    	          END IF
	          PrevIndRef=CurIndRef
	          PrevExtinc=CurExtinc
	        ELSE
	          PrevElt(iElt)=0
	        END IF
C				! Holographic Optical Element (HOE)
	      ELSE IF (EltID(iElt).EQ.4) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        CALL HOE(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	         LEltObs(iElt),
     &	         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	         UDSrfParam(1,iElt),
     &	         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	         lMon(iElt),pMon(1,iElt),
     &	         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	         CurIndRef,IndRef(iElt),
     &	         CurExtinc,Extinc(iElt),WaveBU,
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),
     &	         ObsType(1,iElt),ObsVec(1,1,iElt),
     &	         xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),
     &	         h1HOE(1,iElt),h2HOE(1,iElt),OrderHOE(iElt),
     &	         CWB*WaveHOE(iElt),XYZSRF(1,1,iIntSrf),
     &	         DZSRF(1,1,iIntSrf),IERROR,IWORK(1,iIntSrf),
     &	         DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD)
     &	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Ruled Grating
	      ELSE IF (EltID(iElt).EQ.5) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        CALL Grating(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	         LEltObs(iElt),
     &	         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	         UDSrfParam(1,iElt),
     &	         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	         lMon(iElt),pMon(1,iElt),
     &	         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	         CurIndRef,IndRef(iElt),
     &	         CurExtinc,Extinc(iElt),WaveBU,
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),
     &	         OrderHOE(iElt),RuleWidth(iElt),h1HOE(1,iElt),
     &	         XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	         IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD)
     &	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Ruled TrGrating
	      ELSE IF (EltID(iElt).EQ.13) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        CALL TrGrating(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	         LEltObs(iElt),
     &	         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	         UDSrfParam(1,iElt),
     &	         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	         lMon(iElt),pMon(1,iElt),
     &	         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	         CurIndRef,IndRef(iElt),
     &	         CurExtinc,Extinc(iElt),WaveBU,
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),
     &	         OrderHOE(iElt),RuleWidth(iElt),h1HOE(1,iElt),
     &	         XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	         IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD)
     &	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc

C                                    ! Diffractive Optical Element (DOE)
              ELSE IF (EltID(iElt).EQ.17) THEN
c               PrevElt=PrevNonSeg(iRay)
                PrevNonSeg(iRay)=iElt
                CALL DoeTrGrating(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &           LEltObs(iElt),
     &           LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &           AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &           AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &           UDSrfParam(1,iElt),
     &           psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &           RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &           lMon(iElt),pMon(1,iElt),
     &           xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &           CurIndRef,IndRef(iElt),
     &           CurExtinc,Extinc(iElt),WaveBU,
     &      ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &           ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),
     &           zObs(1,iElt),
     &           RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),OrderHOE(iElt),
     &           XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &           IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
                IF (.NOT.LRayTrans) THEN
                  LRayPass(iRay)=.FALSE.
                  IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &            .AND.(PropType(PrevElt(iElt)).NE.10)
     &            .AND.(PropType(PrevElt(iElt)).NE.12))
     &              WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
                END IF
                IF (ifAccumOPD)
     &            CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
                PrevIndRef=CurIndRef
                PrevExtinc=CurExtinc

C						! Refractor
	      ELSE IF (EltID(iElt).EQ.6) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt

	        CALL Refractor(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	         LEltObs(iElt),
     &	         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	         UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &           GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &	         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	         lMon(iElt),pMon(1,iElt),
     &	         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	         CurIndRef,IndRef(iElt),LNewIndRef,NewIndRef,
     &	         CurExtinc,Extinc(iElt),WaveBU,
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),
     &           zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),
     &	         XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	         IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &		 EltCoat(iElt),EltCoatThk(:,iElt),
     &           IndRefArr(:,iElt),ExtincArr(:,iElt),
     &           IsVarIndRefElt(iElt),IndRefGridArr(iElt),
     &           IsGrinRefElt(iElt),GrinRefArr(iElt))

	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD) THEN
	          If (.not. IsGrinRefElt(iElt)) Then
  	            CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
                  Else
	            ! Graded index lens
                    CumRayL(iRay)=CumRayL(iRay)+RayL(iRay)
	          End If
    	          IF (ifPol) THEN
    	            RayE(1,iRay)=Evec(1)
    	            RayE(2,iRay)=Evec(2)
    	            RayE(3,iRay)=Evec(3)
    	          END IF
    	        END IF
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Conic Obscuration
	      ELSE IF (EltID(iElt).EQ.7) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        CALL ObsSrf(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &	        LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),
     &	        RptElt(1,iElt),RayPos(1,iRay),
     &	        RayDir(1,iRay),RayE(1,iRay),
     &	        CurIndRef,IndRef(iElt),CurExtinc,
     &	        Extinc(iElt),WaveBU,ApType(iElt),
     &	        ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),
     &	        RayDir(1,iRay),RayL(iRay),Evec,dxdx,dxdu)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD) THEN
	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
    	          IF (ifPol) RayE(1:3,iRay)=Evec(1:3)
    	        END IF
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C						! Return element
	      ELSE IF (EltID(iElt).EQ.8) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        ifTransmit=.FALSE.
c	        IF (EltID(PrevElt(iElt)).EQ.8) THEN
c	          ifReturn=.TRUE.
c	        ELSE
c	          ifReturn=.FALSE.
c	        END IF
	        CALL ReturnSrf(ifBuild,ifPol,ifReturn,ifTransmit,
     &	        ifLNsrf,LEltObs(iElt),LRayOK(iRay),LRayTrans,
     &	        KcElt(iElt),KrElt(iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	        CurIndRef,IndRef(iElt),
     &	        CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &	        dxdx,dxdu)
	        IF (ifReturnElt(iElt)) RayL(iRay)=-RayL(iRay)
	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD) THEN
	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
    	          IF (ifPol) RayE(1:3,iRay)=Evec(1:3)
    	        END IF
	        PrevIndRef=CurIndRef
	        PrevExtinc=CurExtinc
C					! Nonsequential Elements
	      ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.12))
     &	      THEN
	        CurIndRef=PNSIndRef
	        CurExtinc=PNSExtinc
	        iminPrev=iFirstNonSeq-1

C	Loop until ray hits next sequential element
	        CrntNSCnt=0
 708	        CONTINUE
	          CrntNSCnt=CrntNSCnt+1
cdebug
c	          IF (ifRayPrint) THEN
c	            CALL DEQUATE(PrevDir,RayDir(1,iRay),3)
c	            CALL DEQUATE(PrevPos,RayPos(1,iRay),3)
c	          END IF
cdebug
C	Find next element as least ray length
	          Lmin=1D22
	          imin=iNextSeq
	          DO 7 itst=iFirstNonSeq,iNextSeq

	            ! test -jzlou
                    crnt_Elt=itst  ! crnt_Elt is shared with
                      ! ChkRayTrans function in elemsub.F

	            LROKdummy=.TRUE.
	            LRTdummy=.TRUE.
c	            IF (ifUseBaseSrf) THEN
c	              CALL FindConSrf(ifBuild,ifLNsrf,LROKdummy,
c    &	              KcElt(itst),KrElt(itst),
c    &	              psiElt(1,itst),VptElt(1,itst),RptElt(1,itst),
c    &	              RayPos(1,iRay),RayDir(1,iRay),
c    &	              RayPosTst,RayLTst)
c	            ELSE
	              CALL FindSrf(ifLNsrf,LEltObs(itst),
     &	              LROKdummy,LRTdummy,KcElt(itst),KrElt(itst),
     &	              AsphCoef(1:tAsphCoef,itst),MonCoef(1,itst),
     &	              AnaCoef(1,itst),UDSrfCoef(1,itst),UDSrfType(itst),
     &	              UDSrfParam(1,itst),nGridMat(itst),mGridMat,
     &	              GridMat(1,1,jGridSrf(itst)),GridSrfdx(itst),
     &	              psiElt(1,itst),VptElt(1,itst),RptElt(1,itst),
     &	              RayPos(1,iRay),RayDir(1,iRay),
     &	              lMon(itst),pMon(1,itst),
     &	              xMon(1,itst),yMon(1,itst),zMon(1,itst),
     &	              ApType(itst),ApVec(1,itst),nObs(itst),
     &	              ObsType(1,itst),ObsVec(1,1,itst),
     &	              xObs(1,itst),yObs(1,itst),zObs(1,itst),
     &	              RayPosTst,RayLTst,SrfType(itst),
     &	              XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	              IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),
     &	              mDP,maxAsphCoef)
c	            END IF
C note: RayLTst changed from 1d-12 to 1d-22
	            IF ((SrfType(iminPrev).LE.2).OR.(iminPrev.NE.itst))
     &	            THEN
	              IF (LROKdummy.AND.(RayLTst.LT.Lmin).AND.LRTdummy.AND.
     &	              (RayLTst.GT.1d-8).AND.((SrfType(iminPrev).LE.2)
     &	              .OR.(iminPrev.NE.itst))) THEN
	                Lmin=RayLTst
	                imin=itst
	              END IF
	            END IF
  7	          CONTINUE
C	Here when element number is determined...
	          IF (Lmin.LT.1d22) THEN
	            IF (imin.NE.iNextSeq) THEN
	              LROKdummy=.TRUE.
	              LRTdummy=.TRUE.
	              IF (EltID(imin).EQ.12) THEN
	                IF ((SrfType(imin).EQ.5).AND.ifIntSrfOK(imin))
     &	                  iIntSrf=iEltToIntSrf(iElt)

	                ! test -jzlou
                        crnt_Elt=imin  ! crnt_Elt is shared with
                           ! ChkRayTrans function in elemsub.F

	                CALL Reflector(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	                 LEltObs(imin),
     &	                 LROKdummy,LRTdummy,KcElt(imin),KrElt(imin),
     &	                 AsphCoef(1:tAsphCoef,imin),MonCoef(1,imin),
     &	                 AnaCoef(1,imin),UDSrfCoef(1,imin),
     &	                 UDSrfType(imin),UDSrfParam(1,imin),
     &	                 nGridMat(imin),mGridMat,
     &	                 GridMat(1,1,jGridSrf(imin)),GridSrfdx(imin),
     &	                 psiElt(1,imin),VptElt(1,imin),RptElt(1,imin),
     &	                 RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	                 lMon(imin),pMon(1,imin),
     &	                 xMon(1,imin),yMon(1,imin),zMon(1,imin),
     &	                 CurIndRef,IndRef(imin),
     &	                 CurExtinc,Extinc(imin),WaveBU,
     &	                 ApType(imin),ApVec(1,imin),nObs(imin),
     &	                 ObsType(1,imin),ObsVec(1,1,imin),
     &	                 xObs(1,imin),yObs(1,imin),zObs(1,imin),
     &	                 RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),
     &	                 Evec,dxidxim1(1,1,imin),dxidui(1,1,imin),
     &	                 SrfType(imin),XYZSRF(1,1,iIntSrf),
     &	                 DZSRF(1,1,iIntSrf),IERROR,IWORK(1,iIntSrf),
     &	                 DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &                   EltCoat(imin),EltCoatThk(:,imin),
     &          	 IndRefArr(:,imin),ExtincArr(:,imin))
	                NextIndRef=CurIndRef
	                NextExtinc=CurExtinc
	              ELSE
	                IF (EltID(PrevNonSeg(iRay)).NE.9) THEN
	                  na=CurIndRef
	                  nb=IndRef(imin)
	                  xa=CurExtinc
	                  xb=Extinc(imin)
	                ELSE IF ((ABS(KrElt(imin)).LT.
     &	                ABS(KrElt(PrevNonSeg(iRay))))) THEN
	                  na=CurIndRef
	                  nb=IndRef(imin)
	                  xa=CurExtinc
	                  xb=Extinc(imin)
	                ELSE
	                  na=IndRef(imin)
	                  nb=IndRef(imin-1)
	                  xa=Extinc(imin)
	                  xb=Extinc(imin-1)
	                END IF
	                CurIndRef=na
	                NextIndRef=nb
	                CurExtinc=xa
	                NextExtinc=xb

                        ! test -jzlou
                        crnt_Elt=imin  ! crnt_Elt is shared with
                           ! ChkRayTrans function in elemsub.F

	                CALL NSRefractor(ifBuild,ifDMBuild,ifPol,
     &	                 ifLNsrf,LEltObs(imin),LROKdummy,LRTdummy,
     &	                 KcElt(imin),KrElt(imin),
     &	                 AsphCoef(1:tAsphCoef,imin),MonCoef(1,imin),
     &	                 AnaCoef(1,imin),UDSrfCoef(1,imin),
     &	                 UDSrfType(imin),UDSrfParam(1,imin),
     &                   nGridMat(imin),mGridMat,
     &                   GridMat(1,1,jGridSrf(imin)),GridSrfdx(imin),
     &	                 psiElt(1,imin),VptElt(1,imin),RptElt(1,imin),
     &	                 RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	                 lMon(imin),pMon(1,imin),
     &	                 xMon(1,imin),yMon(1,imin),zMon(1,imin),
     &	                 na,nb,xa,xb,WaveBU,nCoatElt(imin),
     &	                 CoatIndxElt(1,imin),CoatThkElt(1,imin),
     &	                 ApType(imin),ApVec(1,imin),nObs(imin),
     &	                 ObsType(1,imin),ObsVec(1,1,imin),
     &	                 xObs(1,imin),yObs(1,imin),zObs(1,imin),
     &	                 RayPos(1,iRay),
     &	                 RayDir(1,iRay),RayL(iRay),Evec,
     &	                 dxdx,dxdu,SrfType(imin),
     &	                 XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	                 IWORK(1,iIntSrf),DWORK(1,iIntSrf),
     &	                 NDP(iIntSrf),mDP)
	              END IF

	              IF (.NOT.LRTdummy) THEN
	                LRayPass(iRay)=.FALSE.
	                IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt))
     &	                  .NE.3).AND.(PropType(PrevElt(iElt)).NE.10)
     &	                  .AND.(PropType(PrevElt(iElt)).NE.12))
     &	                WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	              END IF
	              IF (.NOT.LROKdummy) THEN
	                LRayOK(iRay)=.FALSE.
	                PrevElt(iElt)=PrevNonSeg(iRay)
	              ELSE
	                PrevElt(iElt)=PrevNonSeg(iRay)
	                PrevNonSeg(iRay)=imin
	                IF (ifAccumOPD) THEN
	                  CumRayL(iRay)=CumRayL(iRay)+
     &	                  CurIndRef*RayL(iRay)
    	                  IF (ifPol) THEN
    	                    RayE(1,iRay)=Evec(1)
    	                    RayE(2,iRay)=Evec(2)
    	                    RayE(3,iRay)=Evec(3)
    	                  END IF
    	                END IF
	                PrevIndRef=CurIndRef
	                CurIndRef=NextIndRef
	                PrevExtinc=CurExtinc
	                CurExtinc=NextExtinc
cdebug
c	                IF ((ifRayPrint).AND.(iRay.EQ.pRay)) THEN
c	                  WRITE(*,115)pRay,PrevElt,
c     &	                  EltName(PrevElt),imin,EltName(imin),
c     &	                  (PrevPos(i),i=1,3),
c     &	                  (RayPos(i,pRay),i=1,3),
c     &	                  (PrevDir(i),i=1,3),RayL(iRay),
c     &	                  CumRayL(pRay)
c	                  WRITE(*,118)PrevIndRef,PrevExtinc,CurIndRef,
c     &	                  CurExtinc
c	                  IF (ifPol) THEN
c	                    Intxyz=RayE(1,iRay)*DCONJG(RayE(1,iRay))+
c     &	                    RayE(2,iRay)*DCONJG(RayE(2,iRay))+
c     &	                    RayE(3,iRay)*DCONJG(RayE(3,iRay))
c	                    WRITE(*,119)RayE(1,iRay),
c     &	                    RayE(2,iRay),RayE(3,iRay),Intxyz
c	                  END IF
c	                END IF
cdebug
	              END IF
	            END IF
	          ELSE
	            LRayOK(iRay)=.FALSE.
	            PrevElt(iElt)=PrevNonSeg(iRay)
	          END IF
	          IF (.NOT.LRayOK(iRay)) THEN
	            nBadRays=nBadRays+1
	            LEltOK(imin)=.FALSE.
	            IF (IRAY.GT.1) THEN
	              RayMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	              SegMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	            END IF
	            GO TO 9
	          END IF

#if 1
                  ! Tapering complex amplitude for this ray on non-sequential
                  ! element if needed
                  If (LTaperedElt(imin).and.LRayOK(iRay)) Then
                    !i=RayIndex(1,iRay); j=RayIndex(2,iRay)
                    !i=i+(mdttl-npts)/2+1; j=j+(mdttl-npts)/2+1
                    i=RayWFMap(1,iRay); j=RayWFMap(2,iRay)
	            !print*,'*** Ray ',iRay,' is tapered by ',taperFactor(iRay)
                    WFElt(i,j,iWF)=taperFactor(iRay)*WFElt(i,j,iWF)
                  End if
#endif

	          iminPrev=imin

	          if (NSCnt(iElt).gt.0 .and.
     &                CrntNSCnt.eq.NSCnt(iElt)) then
                     imin=iNextSeq
                  end if

	          IF (imin.NE.iNextSeq) THEN
	            GO TO 708
	          END IF
	        iElt=iLastNonSeq
C					! Lens Array
	      ELSE IF (EltID(iElt).EQ.10) THEN
c	        PrevElt=PrevNonSeg(iRay)
	        PrevNonSeg(iRay)=iElt
	        CALL LensArray(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	         LEltObs(iElt),
     &	         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	         UDSrfParam(1,iElt),
     &	         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	         lMon(iElt),pMon(1,iElt),
     &	         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	         CurIndRef,IndRef(iElt),LNewIndRef,NewIndRef,
     &	         CurExtinc,Extinc(iElt),WaveBU,
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec,
     &           dxdx,dxdu,SrfType(iElt),
     &	         LensArrayType(iElt),LensArrayWidth(iElt),
     &           IsVarIndRefElt(iElt),TotalLensletsArray(iElt),iElt,
     &	         XYZSRF(1,1,iIntSrf),
     &	         DZSRF(1,1,iIntSrf),IERROR,IWORK(1,iIntSrf),
     &	         DWORK(1,iIntSrf),NDP(iIntSrf),mDP)

	        IF (.NOT.LRayTrans) THEN
	          LRayPass(iRay)=.FALSE.
	          IF ((iRay.GT.1).AND.(PropType(PrevElt(iElt)).NE.3)
     &	          .AND.(PropType(PrevElt(iElt)).NE.10)
     &	          .AND.(PropType(PrevElt(iElt)).NE.12))
     &	            WFElt(RayWFMap(1,iRay),RayWFMap(2,iRay),iWF)=0d0
	        END IF
	        IF (ifAccumOPD)
     &	          CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	        PrevIndRef=CurIndRef; PrevExtinc=CurExtinc

C  End of element section
	      END IF

	      IF (.NOT.LRayOK(iRay)) THEN
	        nBadRays=nBadRays+1
	        LEltOK(iElt)=.FALSE.
	        IF (IRAY.GT.1) THEN
	          RayMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	          SegMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	        END IF
	      END IF
	    END IF

C	End of iRay loop

	! End of DO 9 loop, all rays traced to current element, -jzlou
  9	CONTINUE

!	print*,'***** Right after DO 9 loop in propsub, iElt,iEndElt =',
!    &         iElt,iEndElt


C    Reset values for next elt loop pass

	IF ((EltID(iElt).EQ.6).OR.(EltID(iElt).EQ.10)) THEN
	  CurIndRef=IndRef(iElt)
	  CurExtinc=Extinc(iElt)
	END IF

	StartIndRef=CurIndRef
	StartExtinc=CurExtinc
	IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.12))
     &	    iElt=iLastNonSeq

cdebug
c	IF (ifRayPrint.AND. (EltType(iElt).NE.13)
c     &	    .AND.(EltType(iElt).NE.24)) THEN
c	  WRITE(*,115) pRay,PrevElt,
c     &	                 EltName(PrevElt),iElt,EltName(iElt),
c     &	                 (PrevPos(i),i=1,3),
c     &	                 (RayPos(i,pRay),i=1,3),
c     &	                 (PrevDir(i),i=1,3),RayL(pRay),
c     &	  CumRayL(pRay)
c	  WRITE(*,118)PrevIndRef,PrevExtinc,CurIndRef,CurExtinc
c	  IF (ifPol) THEN
c	    Intxyz=RayE(1,iRay)*DCONJG(RayE(1,iRay))+
c     &	    RayE(2,iRay)*DCONJG(RayE(2,iRay))+
c     &	    RayE(3,iRay)*DCONJG(RayE(3,iRay))
c	    WRITE(*,119)RayE(1,iRay),
c     &	    RayE(2,iRay),RayE(3,iRay),Intxyz
c	  END IF
c	END IF
cdebug

	If (.false.) Then
	  ! test only
          print*,'**** At element: ',iElt,', LTaperedElt(iElt) = ',
     &      LTaperedElt(iElt)
	  !stop
	End If


	! Apodization effects added here because the user eventually would want
        ! to see the result of diffraction propagation.
	! For tapered (apodized) aperture or tapered obscuration, add tapering
	! effect, which scales the complex E field, but no change to the phase.
	! Here is for a sequential element.
	If ((EltID(iElt).NE.9).AND.(EltID(iElt).NE.12).AND.
     &      LTaperedElt(iElt)) Then
!	  print*,'=> LTaperedElt(iElt), iElt=',iElt,', taperFactor(iRay) =',
!    &             taperFactor(iRay)
          ! Multiply tapering factor to each ray. -jzlou
          Do j=1,mdttl
            Do i=1,mdttl
              iRay=RayID(i,j)
              ! Scale WFElt(i,j,iWF) by tapering factor
              WFElt(i,j,iWF)=taperFactor(iRay)*WFElt(i,j,iWF)
            End Do
          End Do
        End If

#if 0
	If (EltID==17 .and. (GratingMaterial(iElt)==Dielectric .or.
     &        GratingMaterial(iElt)==Metal)) Then
	    ! is FzpTrGrating, and need to include transimission coefficient effect
	    ! to compute diffraction efficiency.

	End If
#endif


C    Begin diffraction propagation section

C      Find starting surface for propagation
	IF ((PropType(iElt).NE.1).AND.
     &	    (PropType(iElt).NE.PropType(iEm1))) THEN
	  IF ((PropType(iEm1).EQ.3).OR.(PropType(iEm1).EQ.10)
     &	      .OR.(PropType(iEm1).GE.12)) THEN
	    j=nEstart
	ELSE
	  j=iElt
	END IF
	  !print*,'**B4 call DXCALC before diffraction propatation 1'
          !print*,'iElt =',iElt
	  !print*,'*nRay =',nRay
	  CALL DXCALC(iElt,j,dxMin(iElt),dxMax(iElt),
     &	         dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
	  !print*,'** dxAvg(iElt) =',dxAvg(iElt)
	  iLend=iElt
	END IF

C      Update parameters during propagation if type 1, 2, 4, 8 or 9
	  IF ((PropType(iElt).EQ.2).OR.(PropType(iElt).EQ.4)
     &	    .OR.(PropType(iElt).EQ.8).OR.(PropType(iElt).EQ.9)) THEN
	    iLend=iElt
	  END IF
	  IF ((PropType(iElt).EQ.1).AND.(PropType(iEm1).EQ.1))
     &	   THEN
	    iLend=iElt
	  END IF

C	Find last surface for propagation
	  IF (((PropType(iEm1).NE.1).AND.
     &	      (PropType(iElt).NE.PropType(iEm1))).OR.
     &	      (iElt.EQ.iEndElt)) THEN
	    iEltToiWF(iElt)=iWF

	    !print*,'**B4 call DXCALC before diffraction propatation 2'
            !print*,'iElt =',iElt
	    !print*,'*nRay =',nRay

C	    Compute scaling
	    CALL DXCALC(iElt,nEstart,dxMin(iElt),dxMax(iElt),
     &	           dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)


C  For reference only -jzlou
C       PropTypeName(1)='Geometric'
C       PropTypeName(2)='NFS1surf'
C       PropTypeName(3)='FarField'
C       PropTypeName(4)='NFP1surf'
C       PropTypeName(5)='NFSpherical'
C       PropTypeName(6)='NFPlane'
C       PropTypeName(7)='SF1surf'
C       PropTypeName(8)='SpatialFilter'
C       PropTypeName(9)='GeomUpdate'
C       PropTypeName(10)='NF1'
C       PropTypeName(11)='NF2'
C       PropTypeName(12)='Fresnel'
C       PropTypeName(13)='NFdftPlane'
C       PropTypeName(14)='NFdftSphere'
C       PropTypeName(15)='FFdft'

C		Geometric propagation
cscale	    z1z2=dxElt(nEstart)/dxElt(iElt)
	    IF (ifPol) THEN
	      IF (PropType(nEstart).NE.11) THEN
	        DO 211 j=1,mdttl
	          DO 210 i=1,mdttl
	            iRay=RayID(i,j)
	            IF ((iRay.GT.nRay).OR.(iRay.LE.0)) THEN
	              IF (PropType(nEstart).NE.11) THEN
	                IF (ifVecDif3) THEN
	                  WFElt(i,j,1)=0D0
	                  WFElt(i,j,2)=0D0
	                  WFElt(i,j,3)=0D0
	                ELSE
	                  WFElt(i,j,iWF)=0D0
	                ENDIF
	              ENDIF
	            ELSE IF (LRayOK(iRay)) THEN
	              IF (ifVecDif3) THEN
	                DO k=1,3
	                  ! For vector diffraction, E fields in X,Y and Z directions
			  ! diffracted separately. -jzlou
	                  WFElt(i,j,k)=RayE(k,iRay)
cscale	                  WFElt(i,j,k)=-z1z2*RayE(k,iRay)
     	                END DO
	              ELSE
	                ! Polarization case:
	                ! Here a complex optical field is formed in
			! WFElt(*,*,iWF) array, with its phase equal
		        ! to OPD and its amplitude equal to the 2-norm
			! of RayE. WFElt array is used in a diffraction
                        ! propagation below. For a more complete picture, read
	                ! comments in the non-polarization section, -jzlou
	                Cphi=CumRayL(iRay)-CumLStart(iRay)
	                ! in the following,
	                ! TPL*Cphi = (2*pi/lambda)*OPL -jzlou
	                DPHI=MOD((TPL*Cphi),6.28318530717959D0)
cscale	                C1=z1z2*SQRT(DBLE(...))
	                C1=SQRT(DBLE(
     &	                	RayE(1,iRay)*DCONJG(RayE(1,iRay))+
     &	                	RayE(2,iRay)*DCONJG(RayE(2,iRay))+
     &	                	RayE(3,iRay)*DCONJG(RayE(3,iRay))))
	                WFElt(i,j,iWF)=C1*EXP(DCMPLX(0D0,DPHI))
	                CumLStart(iRay)=CumRayL(iRay)
	              END IF
	            ELSE
	              IF (ifVecDif3) THEN
	                WFElt(i,j,1:3)=0D0
	              ELSE
	                WFElt(i,j,iWF)=0D0
	              ENDIF
	              CumLStart(iRay)=0D0
	              CumRayL(iRay)=0D0
	            END IF
  210	          CONTINUE
  211	        CONTINUE
	      END IF
	    ELSE
!	      print*,'***** 2: After DO 9 loop in propsub, iElt,iEndElt =',
!    &         iElt,iEndElt
!	      print*,'***** 2: nEstart =',nEstart
!	      print*,'***** 2: PropType(nEstart) =',PropType(nEstart)

	      ! Non-polarization propagation. -jzlou
	      IF (PropType(nEstart).NE.11) THEN   ! 11 = NF2 prop-type
	        DO 11 j=1,mdttl
	          DO 10 i=1,mdttl
	            iRay=RayID(i,j)
	            IF ((iRay.GT.nRay).OR.(iRay.LE.0)) THEN
	              CD=0D0
	            ELSE IF (LRayOK(iRay)) THEN
	              ! This section is reached when
	              !  1) (iElt-1) to iElt is NOT geometric propagation, or
	              !  2) iElt is the last surface of propagation.
		      ! Here the complex optical field in WFElt(*,*,iWF) is updated.
	              ! Geometric propagation has so far reached iElt.
	              ! If propagation from (iElt-1) to iElt is geometric, CumRayL(iRay)
	              ! has been accumulated upto iElt. If propagation from (iElt-1) to
		      ! iElt is physical, CumRayL(iRay) has been accumulated only to (iElt-1).
	              ! So if the propagation from (iElt-1) to iElt is physical, the complex
		      ! field WElt formed below is the input E field for the physical
		      ! propagation below. If the propagation from (iElt-1) to iElt is
                      ! geometric, (which has been completed above), the complex WElt field
                      ! formed below will be the final E field at iElt.  -jzlou
	              !
	              Cphi=CumRayL(iRay)-CumLStart(iRay)
	              DPHI=MOD((TPL*Cphi),6.28318530717959D0)

	              if (iRay==-2) then
	                print*,'***** iRay = 2, i,j,DPHI =',i,j,DPHI
	                print*,'*****  WFElt(i,j,iWF) =',WFElt(i,j,iWF)
	              end if

		      ! OPD is combined with existing complex E field below, -jzlou
cscale	              CD=z1z2*EXP(DCMPLX(0D0,DPHI))*WFElt(i,j,iWF)
	              CD=EXP(DCMPLX(0D0,DPHI))*WFElt(i,j,iWF)
	              CumLStart(iRay)=CumRayL(iRay)
	            ELSE
	              CD=0D0
	              CumLStart(iRay)=0D0
	              CumRayL(iRay)=0D0
	            END IF
	            ! Complex E field is stored in WFElt here
	            WFElt(i,j,iWF)=CD
   10	          END DO
   11	        END DO


#if 0
! To utilize MACOS FFPROP to do far-field diffraction for an externally input
! exit-pupil wavefront (OPD) array. -jzlou, 01/2011
                 Open(Unit=12,file='OPD_in.txt',Status='OLD')
                 Do j=1,mdttl/2
                   Read(12,*) OPD_in(1:mdttl/2,j)
	         End Do
	         Do j=1,mdttl/2
	           Do i=1,mdttl/2
	             DPHI=MOD((TPL*OPD_in(i,j)),6.28318530717959D0)
	             CD=EXP(DCMPLX(0D0,DPHI))
	             WFElt(i+mdttl/4,j+mdttl/4,iWF)=CD
	           End Do
	         End Do
 	         print*,' ***After forming WFElt(i,j,iWF), CumRayL(1)=',
     &                  CumRayL(1)
! When external wavefront is inserted here, it should work for near-field diffraction
! propagation as well - the type of diffraction is selected below. -jzlou, 01/2011
#endif


#if 1
	        ! If at last propagating surface, compute OPDMat - added by jzlou
	        If (iElt==iEndElt) Then
		  IF (LRayOK(1)) THEN
                    RefCumRayL=CumRayL(1)
                  ELSE
                    RefCumRayL=0d0
                  END IF
	          !
	          nPassRays=0
	          Do iRay=2,nRay
		    IF (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &                ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &                ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) THEN
                      nPassRays=nPassRays+1
            	      OPDScal=CumRayL(iRay)-RefCumRayL
                      SumL=SumL+OPDScal
                    END IF
	          End Do
	          !
	          IF (nPassRays.GT.0) THEN
	            DAvgL=SumL/dble(nPassRays)
	            IF (LUseChfRayIfOK .and. LRayOK(1)) Then
                      ! use chief ray OPL
                      Do iRay=2,nRay
                        If (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &          	  ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &          	  ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) Then
                	  OPDScal=CumRayL(iRay)-RefCumRayL
                	  OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))=OPDScal
		        End If
		      End Do
	 	    ELSE
	              ! use average OPL
            	      DO iRay=2,nRay
              	        IF (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &         	  ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &         	  ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) THEN
               	  OPDScal=(CumRayL(iRay)-RefCumRayL)-DAvgL
                  OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))=OPDScal
                        END IF
            	      END DO
	            END IF
		  END IF  ! If nPassRays

	          ! Dump OPD for debugging
                  !Open(Unit=15,File='Opd_propsub.txt')
                  !Do j=1,npts+1
                  !  Write(15,*) OPDMat(1:npts+1,j)
                  !End Do
                  !Close(15)
	        End If  ! If iElt=iEndElt
#endif
	      END IF
	    END IF

#if 1
          ! For computing MTF at exit pupil only -jzlou
          if (ifInMTFCalc .and. iElt==nElt) then
            ! WFElt (complex wavefront) has been computed at exit
            ! pupil (nElt-1), which is used in computing MTF;
            ! no need to diffraction propagation to focal plane
            return
          end if
#endif

C       Geometric propagation has been performed upto iElt.
C
C       Now physically propagate optical field from nEstart to iElt
C       (nEstart and iEm1 are (iElt-1)). The input to physical propagation
C       is the complex E field WFElt formed above (see comments above), which
C       is updated in the physical propagation below. So after physical propagation,
C       both OPD map and WFElt field exist at iElt. Upon next physical propagation,
C       WFElt will be combined with (CumRayL(iRay)-CumLStart(iRay)) to form the
C       needed input field for diffiraction. -jzlou

C       Since geometric propagation has reached iElt, after the physical propagation,
C       all the rays' positions and directions are still determined by the geometric
C       propagation from (iElt-1) to iElt, unless ReGrid is called. -jzlou

	    iLstart=iElt

C					Near-field sphere-to-sphere
	    IF (PropType(iEm1).EQ.2) THEN
	      dxElt(nEstart)=ABS(dxElt(iElt)*zElt(nEstart)/zElt(iElt))
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)

              print*,
     &            '=====> Near-field sphere-to-plane diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL NFPROP(WFElt(1,1,iWF),dxElt(iElt),WFElt(1,1,iWF),
     &	        dxElt(nEstart),zElt(iElt),zElt(nEstart),
     &	        WaveBU,mdttl,DWF)

	      WRITE(*,502)nEstart,iElt,iWF,zElt(nEstart),dxElt(nEstart),
     &	       zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	       dxDev(iElt),Lin2
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					Far-field sphere-to-plane
	    ELSE IF (PropType(iEm1).EQ.3) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)
	      IF (ifVecDif3) THEN
                print*,
     &  '==> 1: Far-field (vector) sphere-to-plane diffraction ...'
		If (ifFnCalc)
     &            Fn = FnCalc_DZ(zElt(nEstart),WaveBU,1)

	        ! Far-field diffraction for polarized light,
	        CALL PFFPROP(WFElt(1,1,1),dxElt(iElt),
     &	               dxElt(nEstart),zElt(nEstart),WaveBU,mdttl)

	        WRITE(*,509)nEstart,iElt,iWF,zElt(nEstart),
     &	          dxElt(nEstart),dxMin(nEstart),dxMax(nEstart),
     &	          dxDev(nEstart),Lin1,zElt(iElt),dxElt(iElt)
	      ELSE
                print*,
     &     '===> 2: Far-field sphere-to-plane diffraction ...'
	        If (ifFnCalc)
     &            Fn = FnCalc_DZ(zElt(nEstart),WaveBU,1)

	        CALL FFPROP(WFElt(1,1,iWF),dxElt(iElt),
     &	          WFElt(1,1,iWF),dxElt(nEstart),zElt(nEstart),
     &	          WaveBU,mdttl,DWF)

	        WRITE(*,503)nEstart,iElt,iWF,zElt(nEstart),
     &	          dxElt(nEstart),dxMin(nEstart),dxMax(nEstart),
     &	          dxDev(nEstart),Lin1,zElt(iElt),dxElt(iElt)
	      END IF

	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C    Now compute obscurations in diffracted light
	      IF (EltID(iElt).EQ.7)  ! 'Obscuring' element
     & 	        CALL FFObscure(psiElt(1,iElt),VptElt(1,iElt),
     &	   RptElt(1,iElt),RayPos(1,1),dxElt(iElt),(npts+1),
     &	   ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	   ObsVec(1,1,iElt),xObs(1,iElt),WFElt(1,1,iWF),mdttl)

C					Near-field plane-to-plane
	    ELSE IF (PropType(iEm1).EQ.4) THEN
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)
	      dxElt(nEstart)=dxElt(iElt)

           print*,'===> Near-field plane-to-plane diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL PPPROP(WFElt(1,1,iWF),dxElt(iElt),WFElt(1,1,iWF),
     &	             dxElt(nEstart),zElt(iElt),zElt(nEstart),
     &	             WaveBU,mdttl,DWF)

	      WRITE(*,504)
     &          nEstart,iElt,iWF,zElt(nEstart),dxElt(nEstart),
     &	        zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	        dxDev(iElt),Lin2
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					Near-field using 2 RS
	    ELSE IF (PropType(iEm1).EQ.5) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)

              print*,'===> Near-field sphere-to-sphere diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL NFPROP(WFElt(1,1,iWF),dxElt(iElt),WFElt(1,1,iWF),
     &	        	  dxElt(nEstart),zElt(iElt),zElt(nEstart),
     &	        	  WaveBU,mdttl,DWF)

	      WRITE(*,505)nEstart,iElt,iWF,zElt(nEstart),dxElt(nEstart),
     &	       dxMin(nEstart),dxMax(nEstart),dxDev(nEstart),Lin1,
     &	       zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	       dxDev(iElt),Lin2

	      ! Optionally re-create OPD grid at iElt if system aberration
	      ! is large.
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C				Near-field plane-to-plane using 2 RS
	    ELSE IF (PropType(iEm1).EQ.6) THEN
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)
	      dxElt(nEstart)=dxElt(iElt)

              print*,'===> Near-field plane-to-plane diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL PPPROP(WFElt(1,1,iWF),dxElt(iElt),WFElt(1,1,iWF),
     &	                  dxElt(nEstart),zElt(iElt),zElt(nEstart),
     &	                  WaveBU,mdttl,DWF)

	      WRITE(*,504)nEstart,iElt,iWF,zElt(nEstart),dxElt(nEstart),
     &	       zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	       dxDev(iElt),Lin2
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

#if 0
	      print*,' ****** nEstart,iElt =',nEstart,iElt
	      print*,' ****** RayPos(1:3,1) =',RayPos(1:3,1)
	      IF (ifRegrid(iElt)) THEN
                 print*,'Calling ReGrid ...'
                 CALL ReGrid(iElt)
              ELSE
                 print*,'NOT Calling ReGrid ...'
              END IF
              stop
#endif

C		Near-field sphere-to-sphere with pin-hole spatial filter
	    ELSE IF (PropType(iEm1).EQ.7) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)

              print*,
     &  '===> Near-field sphere-to-sphere (pin hole) diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL SFPROP(WFElt(1,1,iWF),dxElt(iElt),WFElt(1,1,iWF),
     &	        dxElt(nEstart),zElt(iElt),zElt(nEstart),
     &	        WaveBU,mdttl,DWF,PinHole(nEstart),dxsf)

	      WRITE(*,507)nEstart,iElt,iWF,PinHole(nEstart),dxsf,
     &	       zElt(nEstart),dxElt(nEstart),
     &	       dxMin(nEstart),dxMax(nEstart),dxDev(nEstart),Lin1,
     &	       zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	       dxDev(iElt),Lin2
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C		Near-field sphere-to-sphere with pin-hole spatial filter
	    ELSE IF (PropType(iEm1).EQ.8) THEN
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)

              print*,
     &          '=====> Near-field sphere-to-sphere diffraction ...'
    	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL SFPROP(WFElt(1,1,iWF),dxElt(iElt),WFElt(1,1,iWF),
     &	        dxElt(nEstart),zElt(iElt),zElt(nEstart),
     &	        WaveBU,mdttl,DWF,PinHole(nEstart),dxsf)


	      WRITE(*,508)nEstart,iElt,iWF,PinHole(nEstart),dxsf,
     &	        zElt(nEstart),dxElt(nEstart),
     &	        zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	        dxDev(iElt),Lin2
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)


C  Far-field sphere-to-sphere: first part of S2S calc
	    ELSE IF (PropType(iEm1).EQ.10) THEN
C    NOTE: set dx, z according to nf prop surfaces!!!
	      zStart=zElt(nEstart)
	      zEnd=zElt(iElt+1)
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)

              print*,'===> Near-field NF1/NF2 diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc_DZ(zStart-zEnd,WaveBU,1)

	      CALL SPH2PL(WFElt(1,1,iWF),WFElt(1,1,iWF),
     &	                  dxElt(nEstart),zEnd,zStart,
     &	                  WaveBU,mdttl,DWF)


	      WRITE(*,510)nEstart,iElt,iWF,
     &	       zElt(nEstart),dxElt(nEstart),dxMin(nEstart),
     &	       dxMax(nEstart),dxDev(nEstart),Lin1,
     &	       zElt(iElt),dxElt(iElt)
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C    Now compute obscurations in diffracted light
	      IF (EltID(iElt).EQ.7)
     & 	        CALL FFObscure(psiElt(1,iElt),VptElt(1,iElt),
     &	         RptElt(1,iElt),RayPos(1,1),dxElt(iElt),(npts+1),
     &	         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),WFElt(1,1,iWF),mdttl)

C  Far-field plane-to-sphere: second part of S2S calc
	    ELSE IF (PropType(iEm1).EQ.11) THEN
C    NOTE: set dx, z according to nf prop surfaces!!!
	      Lin2=1d2*(dxElt(iElt)-dxDev(iElt))/dxElt(iElt)

	      CALL PL2SPH(WFElt(1,1,iWF),WFElt(1,1,iWF),mdttl,DWF)

	      WRITE(*,511)nEstart,iElt,iWF,
     &	      zElt(nEstart),dxElt(nEstart),
     &	      zElt(iElt),dxElt(iElt),dxMin(iElt),dxMax(iElt),
     &	      dxDev(iElt),Lin2
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					Fresnel sphere-to-sphere
	    ELSE IF (PropType(iEm1).EQ.12) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)

              print*,'===> Fresnel sphere-to-sphere diffraction ...'
	      If (ifFnCalc)
     &          Fn = FnCalc(zElt(iElt),zElt(nEstart),WaveBU)

	      CALL FRPROP(WFElt(1,1,iWF),dxElt(iElt),
     &	        WFElt(1,1,iWF),dxElt(nEstart),zElt(iElt),
     &	        zElt(nEstart),WaveBU,mdttl,DWF)

	      WRITE(*,512)nEstart,iElt,iWF,zElt(nEstart),
     &	       dxElt(nEstart),dxMin(nEstart),dxMax(nEstart),
     &	       dxDev(nEstart),Lin1,zElt(iElt),dxElt(iElt)
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					NF DFT plane-to-plane
	    ELSE IF (PropType(iEm1).EQ.13) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)

c	      CALL nfpropdft(uin,npt,dx1,utst,npt,dx2,lambda,z)
Cmore efficient	-- CALL WFrepack(DWF,WFElt(1,1,iWF),nGridPts)

	      CALL NFPropDFT(WFElt(1,1,iWF),mdttl,dxElt(nEstart),
     &	             WFElt(1,1,iWF),mdttl,dxElt(iElt),WaveBU,
     &	             zElt(nEstart)-zElt(iElt))

	      WRITE(*,513)nEstart,iElt,iWF,zElt(nEstart),
     &	       dxElt(nEstart),dxMin(nEstart),dxMax(nEstart),
     &	       dxDev(nEstart),Lin1,zElt(iElt),dxElt(iElt)
	       IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					NF DFT sphere-to-sphere
	    ELSE IF (PropType(iEm1).EQ.14) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)

C  Siegman-Sziklas coord xform:
	      alpha=1d0
	      z1=zElt(nEstart)
	      dx1p=dxElt(nEstart)
	      z2=zElt(iElt)
	      dx2p=alpha*dxElt(nEstart)
	      dzp=alpha**2*(z2-z1)*z1/z2

	      CALL NFPropDFT(WFElt(1,1,iWF),mdttl,dx1p,
     &               WFElt(1,1,iWF),mdttl,dx2p,WaveBU,dzp)

	      WRITE(*,514)nEstart,iElt,iWF,zElt(nEstart),
     &	       dxElt(nEstart),dxMin(nEstart),dxMax(nEstart),
     &	       dxDev(nEstart),Lin1,zElt(iElt),dxElt(iElt)
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					FF DFT
	    ELSE IF (PropType(iEm1).EQ.15) THEN
	      Lin1=1d2*(dxElt(nEstart)-dxDev(nEstart))/dxElt(nEstart)
	      CALL FFPropDFT(WFElt(1,1,iWF),mdttl,dxElt(nEstart),
     &	        WFElt(1,1,iWF),mdttl,dxElt(iElt),WaveBU,
     &	        zElt(nEstart)-zElt(iElt))
	      WRITE(*,515)nEstart,iElt,iWF,zElt(nEstart),
     &	       dxElt(nEstart),dxMin(nEstart),dxMax(nEstart),
     &	       dxDev(nEstart),Lin1,zElt(iElt),dxElt(iElt)
	      IF (ifRegrid(iElt)) CALL ReGrid(iElt)

C					Geometric prop
	    ELSE IF ((PropType(iEm1).EQ.9).OR.(iElt.EQ.iEndElt)) THEN
	       WRITE(*,506) nEstart,iElt,iWF

	    END IF

C       Reset RayE following a NF propagation
C           IF (ifPol) THEN
C           END IF
          END IF


C     If starting surface for propagation, compute dx, set next start
	  IF ((PropType(iElt).NE.1).AND.
     &	    ((PropType(iElt).NE.PropType(iEm1)).OR.(iEm1.EQ.0))) THEN
	    nEstart=iElt
c	    CALL DXCALC(iElt,iElt,dxMin(iElt),dxMax(iElt),
c     &	      dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
c	      iLend=iElt

            ! Calc beam diameter at starting surface
	    ddumy=0d0
	    If (ifFnCalc)
     &        Fn=FnCalc_DZ(ddumy(1),ddumy(2),0) ! first two args are place holder only
	  END IF

C	End of propagation section

	  IF ((iElt.GE.nElt).OR.(iElt.GE.iEndElt)) GO TO 27
	  iElt=iElt+1

C	End of iElt loop

 14     CONTINUE
 27     CONTINUE
	IF (iLend.LT.iEndElt) THEN
	  iCurRayElt=iLend
	ELSE
	  iCurRayElt=iEndElt
	END IF
c	DO 15 iElt=iSecondElt,iEndElt
c	  IF (.NOT.ifIntSrfOK(iElt)) THEN
c	    SrfType(iElt)=5
c	    WRITE(*,124)iElt
c	  END IF
 15	CONTINUE
	RETURN
	END SUBROUTINE CPROPAGATE



C***************************************************************************
!	This is the Near Field (Fresnel) sphere-to sphere propagation
!       routine. Z1 and Z2 and are distances from the two (source and image)
!       reference planes to the center of curvature of the spherical
!	source surface (see MACOS manual). When specifiying Z1 and Z2 in Rx,
!       the origin is set at the center of curvature of the spherical source
!       surface, so Z1 and Z2 are both negative, |Z2| < |Z1|, and (Z2-Z1) is
!       positive.
!	Algorithm: Propagation of angular spectrum between planes (Goodman),
!       plus idea of spherical-to-plane transformation (Sziklas & Siegman)
!	which is reflected in the factor Z1/Z2 in computing S.
!       -jzlou
C***************************************************************************

	SUBROUTINE NFPROP(WF2,DX2,WF1,DX1,Z2,Z1,LAMBDA,NFFT,DWF)

	IMPLICIT NONE
	INTEGER I,J,II,JJ,NN(2),C3,NFFT
	REAL*8 DX1,DX2,Z1,Z2,S,LAMBDA,NSQ
	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),N2,M2,
     &		   C4,C5

	NN(1)=NFFT
	NN(2)=NFFT
	NSQ=NFFT*NFFT
	C3=-NFFT/2-1
	S=1D0
	C4=DCMPLX(S)
	S=-3.14159265358979*LAMBDA*(Z2-Z1)*(Z1/Z2)*
     &	  (1E0/DX1/DBLE(NFFT))**2
	C5=DCMPLX(0D0,S)  ! make a pure complex (imaginary) value

	! Shift DC to the center ? -jzlou
	CALL DSWAP2(DWF,WF1,NFFT)
	CALL DFOURN(DWF,NN,2,1)
	CALL DSWAP2(WF2,DWF,NFFT)

	DO 2 J=1,NFFT
	  JJ=C3+J
	  N2=DCMPLX(JJ*JJ)
	  DO 1 I=1,NFFT
	    II=C3+I
	    M2=DCMPLX(II*II)
	    WF2(I,J)=WF2(I,J)*EXP(C5*(M2+N2))/NSQ
    1	  CONTINUE
    2	CONTINUE

	CALL DSWAP2(DWF,WF2,NFFT)
	CALL DFOURN(DWF,NN,2,-1)
	CALL DSWAP2(WF2,DWF,NFFT)

	RETURN
	END SUBROUTINE NFPROP

C ********************************************************************************

C*********************************************************************************
!        Near field plane-to-plane propagation. The code is almost the
!        same as NFPROP except for the lack of Sziklas sphere-to-plane
!        transformation.
!        More about the angular spectrum approach of propagating plane waves
!        between two planes (Goodman):
!        1) This is a general approach applicable to both near and far-field
!           propagations, but with far-field, FFT of input field is simpler.
!        2) Angular spectrum <==> Fourier transform of a field.
!        3) Angular spectrum approach for plane wave propagation:
!           a) Find the angular spectrum of source plane wave - FFT(source-wave).
!           b) "Propagating" source spectrum by distance (Z2-Z1) - multiply
!              the source spectrum by a complex phase factor. The result is
!              is the angular spectrum of the destination wave field.
!           c) Inverse Fourier transform the destination angular spectrum.
!        4) The implemention below assumes the wave propagates in Z direction,
!           a special case of the general formula in Goodman's book.
!        -jzlou
C********************************************************************************

	SUBROUTINE PPPROP(WF2,DX2,WF1,DX1,Z2,Z1,LAMBDA,NFFT,DWF)

	IMPLICIT NONE
	INTEGER I,J,II,JJ,NN(2),C3,NFFT
	REAL*8 DX1,DX2,Z1,Z2,S,LAMBDA,NSQ
	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),N2,M2,C5

	NN(1)=NFFT
	NN(2)=NFFT
	NSQ=NFFT*NFFT
	C3=-NFFT/2-1
	S=-3.14159265358979*LAMBDA*(Z2-Z1)*(1E0/DX1/DBLE(NFFT))**2
	C5=DCMPLX(0D0,S)

	CALL DSWAP2(DWF,WF1,NFFT)
	CALL DFOURN(DWF,NN,2,1)
	CALL DSWAP2(WF2,DWF,NFFT)

	DO 2 J=1,NFFT
	  JJ=C3+J
	  N2=DCMPLX(JJ*JJ)
	  DO 1 I=1,NFFT
	    II=C3+I
	    M2=DCMPLX(II*II)
	    WF2(I,J)=WF2(I,J)*EXP(C5*(M2+N2))/NSQ
    1	  CONTINUE
    2	CONTINUE

	CALL DSWAP2(DWF,WF2,NFFT)
	CALL DFOURN(DWF,NN,2,-1)
	CALL DSWAP2(WF2,DWF,NFFT)

	RETURN
	END SUBROUTINE PPPROP



C***********************************************************************
!	  Far field diffraction propagation -jzlou
C***********************************************************************

	SUBROUTINE FFPROP(WF2,DX2,WF1,DX1,DZ,LAMBDA,NFFT,DWF)
        use Constants, only: PI
        use dft_mod, only: applyfac2
	IMPLICIT NONE
	INTEGER I,J,NN(2),NFFT
	REAL*8 DX1,DX2,DZ,S,LAMBDA,TWOPI,C2
	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),C1

	NN(1)=NFFT
	NN(2)=NFFT
	S=1D0/DBLE(NFFT)
	C1=DCMPLX(S)

	CALL DSWAP2(DWF,WF1,NFFT)
	CALL DFOURN(DWF,NN,2,1)
	CALL DSWAP2(WF2,DWF,NFFT)

	DO J=1,NFFT
	  DO I=1,NFFT
	    !WF2(I,J)=C1*WF2(I,J)
	    WF2(I,J)=WF2(I,J)  ! updated by jzlou, when DFOURN is NR based
          END DO
        END DO

	C1 = 1d0/DCMPLX(0d0,lambda*dz)
	C2 = pi/lambda/dz
	CALL applyfac2(wf2,NFFT,NFFT,C1*dx1*dx1,C2,dx2)

	!wf2 = CDEXP(DCMPLX(0d0,(2d0*pi/lambda)*dz))*wf2   ! added by jzlou

	RETURN
	END SUBROUTINE FFPROP


C***********************************************************************

	SUBROUTINE SFPROP(WF2,DX2,WF1,DX1,Z2,Z1,LAMBDA,NFFT,DWF,PHRAD,
     &	                  DXFF)

	IMPLICIT NONE
	LOGICAL IFPOS
	INTEGER I,J,II,JJ,NN(2),C3,NFFT,II2,JJ2,IRAD2,IPHRAD,IPHRAD2
	REAL*8 DX1,DX2,Z1,Z2,S,LAMBDA,NSQ,PHRAD,DXFF,PHRAD2
 	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),N2,
     &	           M2,C4,C5

	PHRAD2 = 0  ! jzlou

	NN(1)=NFFT
	NN(2)=NFFT
	NSQ=NFFT*NFFT
	C3=-NFFT/2-1
cscale	S=Z1/Z2
	S=1D0
	C4=DCMPLX(S)
	S=-3.14159265358979*LAMBDA*(Z2-Z1)*(Z1/Z2)*
     &	  (1E0/DX1/DBLE(NFFT))**2
	C5=DCMPLX(0D0,S)
	DXFF=LAMBDA*(Z1/DX1)/DBLE(NFFT)

	IF (PHRAD.GE.0D0) THEN
	  IFPOS=.TRUE.
	ELSE
	  IFPOS=.FALSE.
	END IF

	IPHRAD2=(PHRAD*NFFT/2)**2
	IF (PHRAD2.LT.0) IPHRAD2=0

	CALL DSWAP2(DWF,WF1,NFFT)
	CALL DFOURN(DWF,NN,2,1)
	CALL DSWAP2(WF2,DWF,NFFT)

	DO 2 J=1,NFFT
	  JJ=C3+J
	  JJ2=JJ*JJ
	  N2=DCMPLX(JJ2)
	  DO 1 I=1,NFFT
	    II=C3+I
	    II2=II*II
	    IRAD2=II2+JJ2
	    IF (IFPOS) THEN
	      IF (IRAD2.LE.IPHRAD2) THEN
	        M2=DCMPLX(II2)
	        WF2(I,J)=WF2(I,J)*EXP(C5*(M2+N2))/NSQ
	      ELSE
	        WF2(I,J)=0D0
	      END IF
	    ELSE
	      IF (IRAD2.GT.IPHRAD2) THEN
	        M2=DCMPLX(II2)
	        WF2(I,J)=WF2(I,J)*EXP(C5*(M2+N2))/NSQ
	      ELSE
	        WF2(I,J)=0D0
	      END IF
	    END IF
    1	  CONTINUE
    2	CONTINUE

	CALL DSWAP2(DWF,WF2,NFFT)
	CALL DFOURN(DWF,NN,2,-1)
	CALL DSWAP2(WF2,DWF,NFFT)

	DO 4 J=1,NFFT
	  DO 3 I=1,NFFT
	    WF2(I,J)=C4*WF2(I,J)
    3	  CONTINUE
    4	CONTINUE

	RETURN
	END SUBROUTINE SFPROP



C***********************************************************************
!	Near field (Fresnel) diffraction, direct FFT based (not
!       angular spectrum).
!	-jzlou
C***********************************************************************

	SUBROUTINE FRPROP(WF2,DX2,WF1,DX1,Z2,Z1,LAMBDA,NFFT,DWF)

	IMPLICIT NONE
	INTEGER NN(2),NFFT,II2,JJ2, NSQ, I,J
	REAL*8 DX1,DX2,S,SP,LAMBDA,Z1,Z2,N2,M2,PI,C3,II,JJ
	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT)
        COMPLEX*16 C4,C5
C

	NN(1) = NFFT
	NN(2) = NFFT
	NSQ = NFFT*NFFT
	C3  = -NFFT/2-1
c
	S  = 3.14159265358979/LAMBDA/(Z2-Z1)*DX1*DX1*Z1/Z2
	C5 = DCMPLX(0D0,S)
	SP = 3.14159265358979/LAMBDA/(Z2-Z1)*DX2*DX2*Z1/Z2
        C4 = DCMPLX(0D0,SP)
C
	DO 2 J=1,NFFT
	   JJ=C3+J
	   JJ2=JJ*JJ
	   N2=DCMPLX(JJ2,0.0)
	   DO 1 I=1,NFFT
	      II=C3+I
	      II2=II*II
	      M2=DCMPLX(II2,0.0)
	      WF1(I,J)=WF1(I,J)*EXP(C5*(M2+N2))/NFFT
    1	   CONTINUE
    2	CONTINUE
C
	CALL DSWAP2(DWF,WF1,NFFT)
	CALL DFOURN(DWF,NN,2,1)
	CALL DSWAP2(WF2,DWF,NFFT)
C
	DO J=1,NFFT
           JJ=C3+J
	   JJ2=JJ*JJ
	   N2=DCMPLX(JJ2,0.0)
	   DO I=1,NFFT
	      II=C3+I
	      II2=II*II
	      M2=DCMPLX(II2,0.0)
	      WF2(I,J)=WF2(I,J)*EXP(C4*(M2+N2))*DX1/DX2
	   END DO
    	END DO
C
	RETURN
	END SUBROUTINE FRPROP


C***********************************************************************
!	Far-field sphere-to-plane diffraction
!	-jzlou
C***********************************************************************

	SUBROUTINE SPH2PL(WF2,WF1,DX1,Z2,Z1,LAMBDA,NFFT,DWF)

	IMPLICIT NONE
	LOGICAL IFPOS
	INTEGER I,J,II,JJ,NN(2),C3,NFFT,II2,JJ2
	REAL*8 DX1,Z1,Z2,S,LAMBDA,NSQ
	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),N2,M2,
     &		   C4,C5

	NN(1)=NFFT
	NN(2)=NFFT
	C4=1D0/DBLE(NFFT)
	C3=-NFFT/2-1
cscale	S=Z1/Z2
	S=-3.14159265358979*LAMBDA*(Z2-Z1)*(Z1/Z2)*
     &	 (1E0/DX1/DBLE(NFFT))**2
	 C5=DCMPLX(0D0,S)

	CALL DSWAP2(DWF,WF1,NFFT)
	CALL DFOURN(DWF,NN,2,1)
	CALL DSWAP2(WF2,DWF,NFFT)

	DO 2 J=1,NFFT
	  JJ=C3+J
	  JJ2=JJ*JJ
	  N2=DCMPLX(JJ2)
	  DO 1 I=1,NFFT
	    II=C3+I
	    II2=II*II
	    M2=DCMPLX(II2)
	    WF2(I,J)=WF2(I,J)*C4*EXP(C5*(M2+N2))
    1	  CONTINUE
    2	CONTINUE

	RETURN
	END SUBROUTINE SPH2PL

C***********************************************************************

	SUBROUTINE PL2SPH(WF2,WF1,NFFT,DWF)

	IMPLICIT NONE
	INTEGER I,J,NN(2),NFFT
	COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),C1

        !print*,'******* IN PL2SPH '
        !stop

	NN(1)=NFFT
	NN(2)=NFFT
	C1=1D0/DBLE(NFFT)

	DO 2 J=1,NFFT
	  DO 1 I=1,NFFT
	    WF2(I,J)=C1*WF2(I,J)
    1	  CONTINUE
    2	CONTINUE

	CALL DSWAP2(DWF,WF2,NFFT)
	CALL DFOURN(DWF,NN,2,-1)
	CALL DSWAP2(WF2,DWF,NFFT)

	RETURN
	END SUBROUTINE PL2SPH

C***********************************************************************

	SUBROUTINE PFFPROP(WF,DX1,DX2,DZ,LAMBDA,NFFT)
	IMPLICIT NONE
	INTEGER I,J,K,NN(2),NFFT
	REAL*8 DX1,DX2,DZ,S,LAMBDA
	COMPLEX*16 WF(NFFT,NFFT,3),C1

	NN(1)=NFFT
	NN(2)=NFFT
	S=1D0/DBLE(NFFT)
	C1=DCMPLX(S)

	DO 3 K=1,3
	  CALL DSWAP3(WF(1,1,K),NFFT)
	  CALL DFOURN(WF(1,1,K),NN,2,1)
	  CALL DSWAP3(WF(1,1,K),NFFT)

	  DO 2 J=1,NFFT
	    DO 1 I=1,NFFT
	      WF(I,J,K)=C1*WF(I,J,K)
    1	    CONTINUE
    2	  CONTINUE
    3	CONTINUE

	RETURN
	END SUBROUTINE PFFPROP

C***********************************************************************

	SUBROUTINE FFObscure(psi,pv,prot,pCR,dx,nGridPts,
     &		    ApType,ApVec,nObs,ObsType,ObsVec,xGrid,WF,NFFT)

C  Turn off wavefront points that are obscured.

	IMPLICIT NONE
	INTEGER i,ApType,nObs,ObsType(nObs),iObs,NFFT,iWF,jWF,
     &	nGridPts,NFCTR
	LOGICAL ifHitsPosObs,ifHitsNegObs,ifHitsAp
	REAL*8 pv(3),prot(3),pCR(3),psi(3),dx,xCell,yCell
	REAL*8 prel(3),pxCR,pyCR,dxCtr,
     &	xGrid(3),yGrid(3),zGrid(3),px,py,ApVec(6),ObsVec(6,nObs)
	REAL*8 TwoPi,Kd,K123,q1xq2y,q2xq3y,q3xq1y,
     &	q1yq2x,q2yq3x,q3yq1x,q1x,q1y,q2x,q2y,q3x,q3y,rc,
     &	xcmin,xcmax,ycmin,ycmax,yt,R2,A2,A3

	COMPLEX*16 WF(NFFT,NFFT)

	DATA TwoPi/6.28318530717959D0/

	NFCTR=NFFT/2+1

C  Compute element coordinates
	CALL DEQUATE(zGrid,psi,3)
	CALL DUNITIZE(zGrid)
	CALL DXPROD(yGrid,zGrid,xGrid)
	CALL DUNITIZE(yGrid)
	CALL DXPROD(xGrid,yGrid,zGrid)

C  Find center of wf on element
	IF (MOD(nGridPts,2).EQ.0) THEN
	  dxCtr=0.5d0
	ELSE
	  dxCtr=0d0
	END IF
	DO 1 i=1,3
	  prel(i)=pCR(i)-pv(i)
   1	CONTINUE
	pxCR=DDOTC(xGrid,prel)
	pyCR=DDOTC(yGrid,prel)

C Loop through each grid pt and check if obscured
	DO jWF=1,NFFT
	  yCell=(DBLE(jWF-NFCTR)+dxCtr)*dx
	  py=pyCR+yCell
	  DO iWF=1,NFFT

C  Compute grid pt location on element
	    xCell=(DBLE(iWF-NFCTR)+dxCtr)*dx
	    px=pxCR+xCell

C  Compute aperture and obscurations
	    ifHitsAp=.TRUE.
	    IF (ApType.EQ.1) THEN
	      r2=(ApVec(2)-px)**2+(ApVec(3)-py)**2
	      A2=ApVec(1)**2
	      IF (r2.GT.A2) ifHitsAp=.FALSE.
	    ELSE IF (ApType.EQ.2) THEN
	      IF ((px.LT.ApVec(1)).OR.(px.GT.ApVec(2)).OR.
     &	      (py.LT.ApVec(3)).OR.(py.GT.ApVec(4))) ifHitsAp=.FALSE.
	    ELSE IF (ApType.EQ.3) THEN
	      IF ((px.LT.ApVec(1)).OR.(px.GT.ApVec(2)).OR.
     &	      (py.LT.ApVec(3)).OR.(py.GT.ApVec(4))) ifHitsAp=.FALSE.
	      rc=ApVec(5)
	      xcmin=ApVec(1)+rc
	      xcmax=ApVec(2)-rc
	      ycmin=ApVec(3)+rc
	      ycmax=ApVec(4)-rc
	      IF ((px.LT.xcmin).AND.(py.LT.ycmin)) THEN
	        yt=ycmin-SQRT(rc*rc-(px-xcmin)**2)
	        IF (py.LT.yt) ifHitsAp=.FALSE.
	      ELSE IF ((px.GT.xcmax).AND.(py.LT.ycmin)) THEN
	        yt=ycmin-SQRT(rc*rc-(px-xcmax)**2)
	        IF (py.LT.yt) ifHitsAp=.FALSE.
	      ELSE IF ((px.LT.xcmin).AND.(py.GT.ycmax)) THEN
	        yt=ycmax+SQRT(rc*rc-(px-xcmin)**2)
	        IF (py.GT.yt) ifHitsAp=.FALSE.
	      ELSE IF ((px.GT.xcmax).AND.(py.GT.ycmax)) THEN
	        yt=ycmax+SQRT(rc*rc-(px-xcmax)**2)
	        IF (py.GT.yt) ifHitsAp=.FALSE.
	      END IF
	    END IF

C  Check if grid pt hits an obscuration
	    ifHitsPosObs=.FALSE.
	    ifHitsNegObs=.FALSE.
	    IF (ifHitsAp) THEN
	      DO 10 iObs=1,nObs
	        IF (ABS(ObsType(iObs)).EQ.1) THEN
	          r2=(ObsVec(2,iObs)-px)**2+(ObsVec(3,iObs)-py)**2
	          A2=ObsVec(1,iObs)**2
	          IF (r2.LT.A2) THEN
	            IF (ObsType(iObs).GT.0) THEN
	              ifHitsPosObs=.TRUE.
	            ELSE
	              ifHitsNegObs=.TRUE.
	            END IF
	          END IF
	        ELSE IF (ObsType(iObs).EQ.2) THEN
	          IF ((px.GT.ObsVec(1,iObs)).AND.(px.LT.ObsVec(2,iObs))
     &	          .AND.(py.GT.ObsVec(3,iObs)).AND.
     &	          (py.LT.ObsVec(4,iObs))) THEN
	            IF (ObsType(iObs).GT.0) THEN
	              ifHitsPosObs=.TRUE.
	            ELSE
	              ifHitsNegObs=.TRUE.
	            END IF
 	          END IF
	        ELSE IF (ObsType(iObs).EQ.3) THEN
	          r2=(ObsVec(2,iObs)-px)**2+(ObsVec(3,iObs)-py)**2
	          A2=ObsVec(1,iObs)**2
	          A3=ObsVec(4,iObs)**2
	          IF ((r2.LT.A2).AND.(r2.GT.A3)) THEN
	            IF (ObsType(iObs).GT.0) THEN
	              ifHitsPosObs=.TRUE.
	            ELSE
	              ifHitsNegObs=.TRUE.
	            END IF
 	          END IF
	        ELSE IF (ObsType(iObs).EQ.4) THEN
	          q1x=ObsVec(1,iObs)
	          q1y=ObsVec(2,iObs)
	          q2x=ObsVec(3,iObs)
	          q2y=ObsVec(4,iObs)
	          q3x=ObsVec(5,iObs)
	          q3y=ObsVec(6,iObs)
	          q1xq2y=q1x*q2y
	          q2xq3y=q2x*q3y
	          q3xq1y=q3x*q1y
	          q1yq2x=q1y*q2x
	          q2yq3x=q2y*q3x
	          q3yq1x=q3y*q1x
	          K123=DABS(q1xq2y+q2xq3y+q3xq1y-q1yq2x-q2yq3x-q3yq1x)
	          Kd=DABS(px*q2y+q2xq3y+q3x*py-py*q2x-q2yq3x-q3y*px)
	          IF (Kd.GT.K123) THEN
	            GO TO 10
	          ELSE
	           Kd=Kd+DABS(q1x*py+px*q3y+q3xq1y-q1y*px-py*q3x-q3yq1x)
	            IF (Kd.GT.K123) THEN
	              GO TO 10
	            ELSE
	           Kd=Kd+DABS(q1xq2y+q2x*py+px*q1y-q1yq2x-q2y*px-py*q1x)
	              IF (Kd.GT.K123) THEN
	                GO TO 10
	              ELSE
	                IF (ObsType(iObs).GT.0) THEN
	                  ifHitsPosObs=.TRUE.
	                ELSE
	                  ifHitsNegObs=.TRUE.
	                END IF
	              END IF
	            END IF
	          END IF
	        ELSE IF (ABS(ObsType(iObs)).EQ.5) THEN
	          r2=((ObsVec(3,iObs)-px)/ObsVec(1,iObs))**2
     &            +((ObsVec(4,iObs)-py)/ObsVec(2,iObs))**2
	          IF (r2.LT.1.d0) THEN
	            IF (ObsType(iObs).GT.0) THEN
	              ifHitsPosObs=.TRUE.
	            ELSE
	              ifHitsNegObs=.TRUE.
	            END IF
	          END IF
	        END IF
  10	      CONTINUE
	    END IF

C  Reset wf if grid pt is obscured
	    IF (.NOT.ifHitsAp.OR.(ifHitsPosObs.AND..NOT.ifHitsNegObs))
     &	      WF(iWF,jWF)=0d0

	  END DO
	END DO

	RETURN
	END SUBROUTINE FFObscure

C *****************************************************************************

        Real*8 Function FnCalc(Z2,Z1,Lambda)
        Implicit none

        Real*8 :: Z2,Z1,Lambda,DZ
        Real*8, parameter :: eps_fncalc = 1.0D-16

        !FnCalc=0d0 ! debug
        !return

        if (abs(Z2-Z1) < eps_fncalc) then
          print*,
     &     'Fncalc() Warning: Diffraction Distance (Z2-Z1) Almost Zero!'
          return
        end if

	DZ=abs(Z2-Z1)
        FnCalc = FnCalc_DZ(DZ,Lambda,1)
        End Function FnCalc

C ***************************************************************************

        Real*8 Function FnCalc_DZ(DZ,Lambda,op_flag)
        use param_mod
        use elt_mod

        Implicit none

        Integer :: op_flag, i
        Real*8 :: DZ, Lambda, wv_len
        Real*8 :: dd, RayOne(3), dx(3)
        Real*8, save :: beam_dm = 0.0D0

	FnCalc_DZ=0d0

        If (op_flag == 0) Then
          RayOne(1:3) = RayPos(1:3,1)
          beam_dm = 0.0D0
          if (nRay > 1) then
            Do i=2, nRay
              dx(1) = RayPos(1,i)-RayOne(1)
              dx(2) = RayPos(2,i)-RayOne(2)
              dx(3) = RayPos(3,i)-RayOne(3)
              dd = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)
              if (dd > beam_dm) beam_dm = dd
            End Do
          else if (nRay == 1) then
              print*,
     &          'Fncalc_DZ(): only found a single ray in the beam!'
              return
          end if
	  print*,'=====> Beam diameter estimate = ', dd
          return  ! compute beam diameter only
        End If

        DZ = abs(DZ)
        wv_len = Lambda*1.0  ! no conversion

        dd = 2.0*dsqrt(beam_dm)  ! estimate of beam diameter
        FnCalc_DZ = 4.*beam_dm/(DZ*wv_len)  ! Frensel # = D/(Z*lambda)

        print*,'=====> Beam diameter estimate = ', dd
        print*,'=====> DZ = ', DZ
        print*,'=====> Wave length = ', wv_len
        print*,'=====> Fresnel number of propagation =', FnCalc_DZ
        End Function FnCalc_DZ

        End Module propsub_mod
