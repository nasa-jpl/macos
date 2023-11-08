C***********************************************************************
C	Begin file tracesub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module tracesub_mod
	use math_mod
	use dopt_mod

	Logical, save :: isDOP=.false.

        Contains

! npts=nGridpts-1, and nGridpts is read in from Rx in command-line macos. 

	SUBROUTINE OPD(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,OPDMat,WFE,WFEPV,DAvgL,
     &	  L1,Cmatrix,DrawRayVec,DrawEltVec)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
        use sourcsub_mod
        use macos_mod, only : cUnit
	use loh_I1, only : BaseUnits_FLG
	use traceutil_mod, only : 
     &        LUseChfRayIfOK,opdRayMask,LCalcWFElt,
     &        LShowDTime,nValidRays,OPDRefRayLen_FLG,
     &        OPDRefRayLen,ifCtraceCalled,LRayOKHist,
     &        RayPosHist,SaveRayPosHist_Flg,SaveVis3d_dat_Flg
#ifdef CMACOS
        use macos_vars_mod, only : DAvgCumL
#else
        use smacos_vars_mod, only : DAvgCumL
#endif


	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifOPDCalc,ifRayPrint,
     &    ifLNegOK,ifBuild,ifTrace,ifRegrid(mElt),ifVecDif3,
     &    ifIntSrfSet(mElt),ifDMBuild,ifDraw,L1(mRay),ifCalWFElt

	INTEGER npts,i,j,ICLEN,iCurWFElt,iCurRayElt,ie,
     &	  nBadRays,iStartElt,iEndElt,nPassRays,iObsOpt,
     &    iStartRay,iEndRay,nDrawElt,iElt,iEm1,iRay,mpt2,
     &    DrawElt(mDrawElt),DrawEltVec(mDrawElt,mDrawRay)

	SREAL StartTime,DTime,EndTime,
     &	  Cmatrix(7,mCm,bRay),DrawRay(2,mDrawElt),
     &	  DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,Dn,
     &    WFE,WFEPV,MinL,MaxL,DAvgL,SumL,SumL2,SumCumL,
     &    OPDScal,OPDMat(mpts,mpts),dxMin,dxAvg,dxMax,dxDev,
     &    xDraw(3),yDraw(3),RefCumRayL,WFmax,WFmin
! jzlou: RefCumRayL = RayL for chief ray

 116    FORMAT(' Compute time was',f10.4,' sec')
 117	FORMAT(' Tracing ',i9,' rays...')

	iElt=iEndElt
	IF (ifTrace .AND. (iCurRayElt.LE.iEndElt)) THEN
	  iStartElt=iCurRayElt
	ELSE IF (ifPropagate .AND. (iCurRayElt.EQ.iEndElt)) THEN
	  iStartElt=iCurRayElt
	ELSE
	  iStartElt=0
	END IF

	! If an OPD command has been previously run, the next OPD
        ! to the same element will not get into the following 'if'
	! section because iStartElt=iEndElt, so redundant calls to
        ! CTRACE routine is avoided. -jzlou

	IF ((iStartElt.EQ.0).OR.(iStartElt.LT.iEndElt)) THEN
	  CALL RTIME(StartTime)
	  IF (iStartElt.EQ.0) THEN
	    IF (ABS(zSource).GE.1d10) THEN
              ! source at infinity - jzlou
	      CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	                     ifPropagate)
	    ELSE
	      CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	                    ifPropagate)
	    END IF
	    ifGrid=.TRUE.
	  END IF
	  if (.not. isDOP)
     &      WRITE(*,*) 'OPD: Rays being traced = ', nRay
	  iStartRay=1
	  iEndRay=nRay

          If (SaveRayPosHist_Flg .OR. SaveVis3d_dat_Flg) Then
            ! save ray positions at source, all source rays are assumed valid.
            Do i=1,nRay
              LRayOKHist(i,0) = .true. 
              RayPosHist(1:3,i,0)=RayPos(1:3,i) 
            End Do
          End If

          LRayOKHist(:,1:nElt)=.FALSE. ! default, 09/2014 

	  CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	    ifLNegOK,Cmatrix,
     &	    iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	    StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

	  iCurRayElt=iEndElt
	  ifTrace=.TRUE.
          ifCtraceCalled=.TRUE.
	  ifOPDCalc=.TRUE.
	  ifBuild=.FALSE.
	  ifDMBuild=.FALSE.
	  ifPropagate=.FALSE.
	  ifRayPrint=.FALSE.
	  CALL RTIME(EndTime)
	  DTime=EndTime-StartTime
	  if ((.not.isDOP) .and. (LShowDTime)) WRITE(*,116) DTime
        ELSE
          ifCtraceCalled=.FALSE.
	END IF

	nGridPts=npts+1
	SumL=0d0
	SumL2=0d0
        SumCumL=0d0
	MaxL=-1d38
	MinL=1d38
	mpt2=mpts*mpts
	CALL DZERO(OPDMat,mpt2)

	! Following code calculates RMS WFE, in which chief ray
	! is not really needed (would be cancelled in RMS formula) -jzlou
cstartscottfix
	IF (LRayOK(1)) THEN
	   RefCumRayL=CumRayL(1)
	ELSE IF (OPDRefRayLen_FLG) THEN
           RefCumRayL=OPDRefRayLen
	ELSE
	   RefCumRayL=0d0
	END IF
cendscottfix

#if 1
	! Define ray mask: 1 for ray not-lost AND passed, 0 for otherwise,
	! for the purpose of sensitivity calculations
	opdRayMask=0
	Do iRay=2,nRay
	  If (LRayOK(iRay).AND.LRayPass(iRay)) 
     &	    opdRayMask(RayIndex(1,iRay),RayIndex(2,iRay))=1
	End Do
#endif

! iObsOpt 
!  = 0 to plot all rays (unobscured and obscured) on spot diagram
!  = 1 to plot unobscured rays only
!  = 2 to plot obscured rays only

	nPassRays=0
	DO 7 iRay=2,nRay
	  IF (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &	    ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &	    ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) THEN
	    nPassRays=nPassRays+1
            passRayIds(nPassRays)=iRay 
	    OPDScal=CumRayL(iRay)-RefCumRayL
	    SumL=SumL+OPDScal
	    SumCumL=SumCumL+CumRayL(iRay)
            !print*,'iRay =',iRay,',  CumRayL(iRay) =',CumRayL(iRay)
	  END IF
  7	CONTINUE
        nValidRays=nPassRays

	IF (nPassRays.GT.0) THEN
	  Dn=nPassRays
	  DAvgL=SumL/Dn  ! SumOfOPD/numPassRays=averageOPD, -jzlou
          DAvgCumL=SumCumL/Dn
          If (BaseUnits_FLG) Then
	    WRITE(*,11) nPassRays,DAvgCumL,cUnit
 11         FORMAT(2P,' OPD: nPassRays =',I8,',','    DAvgCumL =',
     &             D17.9,A3)
	  Else
            WRITE(*,12) nPassRays,DAvgCumL
 12         FORMAT(2P,' OPD: nPassRays =',I8,',','    DAvgCumL =',D17.9) 
	  End If
	  WFE=0d0; WFmax=-1d22; WFmin=1d22
	  IF (LUseChfRayIfOK .and. LRayOK(1)) Then
	    ! use chief ray OPL
            !print*,'***** OPD: using chfRay OPL, RefCumRayL =',RefCumRayL 
	    DO iRay=2,nRay
              IF (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &          ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &          ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) THEN
	        !print*,'**** ChfRay Len =',RefCumRayL
	        OPDScal=CumRayL(iRay)-RefCumRayL 
	        WFE=WFE+(CumRayL(iRay)-RefCumRayL-DAvgl)**2
	        OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))=OPDScal
	        If (OPDScal>WFmax) WFmax=OPDScal
	        If (OPDScal<WFmin) WFmin=OPDScal
	      END IF
	    END DO
	  ELSE IF (OPDRefRayLen_FLG) THEN
            ! use Rx provided chief ray OPL - added by jzlou
            DO iRay=2,nRay
              IF (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &          ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &          ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) THEN
                OPDScal=CumRayL(iRay)-RefCumRayL
                WFE=WFE+(CumRayL(iRay)-RefCumRayL-DAvgl)**2
                OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))=OPDScal
                If (OPDScal>WFmax) WFmax=OPDScal
                If (OPDScal<WFmin) WFmin=OPDScal
              END IF
            END DO
            WFEPV=WFmax-WFmin
	  ELSE
	    ! use average OPL
            !print*,'***** OPD: using averaged OPL ...'
            DO iRay=2,nRay
	      IF (LRayOK(iRay).AND.((iObsOpt.EQ.0).OR.
     &          ((iObsOpt.EQ.1).AND.(LRayPass(iRay))).OR.
     &          ((iObsOpt.EQ.2).AND.(.NOT.LRayPass(iRay))))) THEN
      	        OPDScal=(CumRayL(iRay)-RefCumRayL)-DAvgl
	        WFE=WFE+OPDScal*OPDScal
	        OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))=OPDScal
	        If (OPDScal>WFmax) WFmax=OPDScal
	        If (OPDScal<WFmin) WFmin=OPDScal
	      END IF
	    END DO
	  END IF
	  WFEPV=WFmax-WFmin
	  IF (nPassRays.GT.1) THEN
            WFE=SQRT(WFE/(Dn-REAL(1d0,KIND(WFE))))
          ELSE
            WFE=SQRT(WFE/Dn)
          END IF
	  !WFE=SQRT(WFE/Dn)

	  If (LCalcWFElt) Then

	  End If
	ELSE
	  WRITE(*,*)' All rays are obscured or lost.'
	  DAvgl=99999d32
	  WFE=99999d32; WFEPV=0d0
	END IF
	END SUBROUTINE OPD


C***********************************************************************

	SUBROUTINE FEX(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,
     &	  Cmatrix,DrawRayVec,DrawEltVec,
     &	  iElt,zp,psip,CrossPt,ifCentroid,
     &	  CRSpot,CentroidSpot,dCRCentroid)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod

        use sourcsub_mod
	use traceutil_mod, only : Rx_FEXCentrFlg

	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifRayPrint,ifLNegOK,
     &	  ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),ifDMBuild,
     &	  ifDraw,ifCentroid

	INTEGER i,iSpot,iRay,mpt2,npts,iCurWFElt,iCurRayElt,
     &	  nBadRays,iStartElt,iEndElt,iObsOpt,iStartRay,
     &	  iEndRay,nDrawElt,DrawEltVec(mDrawElt,mDrawRay),
     &	  iElt,iEm1

	SREAL Cmatrix(7,mCm,bRay),DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &	  xDraw(3),yDraw(3),CentroidSpot(3),dCRCentroid(3),CRSpot(3),
     &	  DV1(9),S1,th(3),
     &	  oldChfRayDir(3),oldChfRayPos(3),
     &	  zp,psip(3),CrossPt(3),cr1dir(3),cr1pos(3),cr2dir(3),cr2pos(3),
     &	  cr2indir(3),cr2inpos(3),CrossDist

C         ***Find exit pupil (with option to minimize tilt)***

	IF (ABS(zSource).GE.1d10) THEN
	  CALL ColSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	    ifPropagate)
	ELSE
	  CALL PtSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	    ifPropagate)
	END IF

	iEm1=iElt-1
	iStartRay=1
	iEndRay=1
	iStartElt=0
	iEndElt=iEm1 
	CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &		ifLNegOK,Cmatrix,
     &		iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &		ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &		StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

	CALL DEQUATE(cr1dir,RayDir(1,1),3)
	CALL DEQUATE(cr1pos,RayPos(1,1),3)
        !print*,' ** FEX at focal-return: cr1pos =', cr1pos
        !print*,' ** FEX at focal-return: cr1dir =', cr1dir
c	th(1)=3.5355d-6
c	th(2)=3.5355d-6
c	th(3)=3.5355d-6
	!CALL DSMPROD(th,xGrid,1d-5,3)  ! jzlou trying a smaller dangle
	CALL DSMPROD(th,xGrid,5d-6,3)
	CALL DEQUATE(cr2indir,ChfRayDir,3)
	CALL DEQUATE(cr2inpos,ChfRayPos,3)
	CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
	CALL DEQUATE(RayDir(1,1),cr2indir,3)
	CALL DEQUATE(RayPos(1,1),cr2inpos,3)

	CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	  	 ifLNegOK,Cmatrix,
     &	  	 iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	  	 ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	         StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
	CALL DEQUATE(cr2dir,RayDir(1,1),3)
	CALL DEQUATE(cr2pos,RayPos(1,1),3)
        !print*,' ** FEX at focal-return: cr2pos =', cr2pos
        !print*,' ** FEX at focal-return: cr2dir =', cr2dir

	CALL FindCrossPt(cr1pos,cr1dir,cr2pos,cr2dir,zp,
     &		         CrossDist,CrossPt)

	print*,'***** FEX: CrossPt =',CrossPt
        print*,'***** FEX: zp =',zp

	! Rx_FEXCentrFlg below is controlled in optical Rx, which is defaulted 
        ! to TRUE if FEXCentroid flag is not referenced in Rx.  -jzlou

C  Centroid or chief ray option tested
	IF (ifCentroid .AND. Rx_FEXCentrFlg) THEN
          print*,'***** FEX: ifCentroid used!\n'

	  IF (ifTrace) THEN
	    i=iCurRayElt
	  ELSE IF (ifPropagate) THEN
	    i=iCurWFElt
	  ELSE 
	    i=0
	  END IF
	  iEndElt=iElt+1
	  IF (ifTrace.AND.(iCurRayElt.LE.iEndElt)) THEN
	    iStartElt=iCurRayElt
	  ELSE IF (ifPropagate.AND.(iCurWFElt.EQ.iEndElt)) THEN
	    iStartElt=iCurWFElt
	  ELSE
	    iStartElt=0
	  END IF

	  IF ((iStartElt.EQ.0).OR.(iStartElt.LT.iEndElt)) THEN
	    ifRayPrint=.FALSE.
	    IF (iStartElt.EQ.0) THEN
	      IF (ABS(zSource).GE.1d10) THEN
	        CALL ColSource(nTestPts,Dicr,Djcr,dxSource,ifPol,
     &	        ifPropagate)
	      ELSE
	        CALL PtSource(nTestPts,Dicr,Djcr,dxSource,ifPol,
     &	        ifPropagate)
	      END IF
	      ifGrid=.TRUE.
	    END IF
	    !WRITE(*,*)' Tracing ',nRay,' rays...'
	    iStartRay=1
	    iEndRay=nRay

	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	      ifLNegOK,Cmatrix,
     &	      iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	      ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	      StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

	    iCurRayElt=iEndElt
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifPropagate=.FALSE.
	  END IF

	  iSpot=0
	  nGridPts=npts+1
	  mpt2=mpts*mpts
	  DO 23 iRay=1,nRay
	    IF (LRayOK(iRay).AND.((iRay.EQ.1).OR.(iObsOpt.EQ.0).OR.
     &	        (LRayPass(iRay).AND.(iObsOpt.EQ.1)))) THEN
	      IF (iRay.EQ.1) THEN
	        CALL DEQUATE(CRSpot,RayPos(1,1),3)
	        CALL DZERO(DV1,3)
	      END IF
	      IF ((iRay.NE.1).OR.
     &	        ((iRay.EQ.1).AND.(nGridpts.EQ.1))) THEN
                iSpot=iSpot+1
                CALL DADD(DV1,DV1,RayPos(1,iRay),3)
	      END IF
	    END IF
 23	  CONTINUE
	  IF (iSpot.GT.0) THEN
            S1=1d0/DBLE(iSpot)
	    CALL DSMPROD(CentroidSpot,DV1,S1,3)
	    CALL DSUB(dCRCentroid,CRSpot,CentroidSpot,3)
	  ELSE
	    WRITE(*,*)' All rays were obscured or lost'
	  END IF

	  IF (nBadRays.GT.0) CALL WARN(LEltOK,iEndElt,nBadRays)

	  CALL DSUB(psip,CentroidSpot,CrossPt,3)
	  CALL DUNITIZE(psip)
	ELSE
	  CALL DNEGATE(psip,cr1dir,3)
	  print*,'***** FEX: Centroid NOT used!'
	END IF
	END SUBROUTINE FEX


C***********************************************************************

	SUBROUTINE FDPN(dpElt,npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &    ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &    ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &    Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
	use sourcsub_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
	use traceutil_mod

	IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: ANS 
        LOGICAL :: ifPol,ifGrid,ifPropagate,ifRayPrint,
     &    ifLNegOK,ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &    ifDMBuild,ifDraw,ifLocalCoord

        INTEGER :: dpElt,npts,ICLEN,iCurWFElt,iCurRayElt,
     &             nBadRays,iStartElt,iEndElt,iObsOpt,
     &             DrawEltVec(mDrawElt,mDrawRay)

        SREAL :: Cmatrix(7,mCm,bRay),
     &           DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 :: Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &          ChfRayDirSv(3),ChfRayPosSv(3),cr1Pos(3),cr1Dir(3),
     &          cr2Pos(3),cr2Dir(3),zp,th(3),CrsDist,DplnPos(3),
     &          DplnPosSag(3),DplnPosTan(3)

	IF (ABS(zSource).GE.1d10) THEN
	  CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                   ifPropagate)
        ELSE
	  CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                  ifPropagate) 
        END IF

	ChfRayDirSv(1:3)=ChfRayDir(1:3)
        ChfRayPosSv(1:3)=ChfRayPos(1:3)	

	! Trace nominal chief ray 
	iStartElt=0; iEndElt=dpElt-1  ! from source to the elt before detector
	!ifRayPrint=.TRUE.
	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)
        cr1Pos=CRIncidPos(1:3,iEndElt); cr1Dir=CROutgoDir(1:3,iEndElt)

	!print*,'cr1Dir =', cr1Dir
	!print*,'cr1Pos =', cr1Pos

	IF (ABS(zSource).GE.1d10) THEN
	  ! Collimated source
          ! Trace a differential chief ray displaced in -yGrid 
          CALL DSMPROD(th,-yGrid,1d-5,3)
          CALL DADD(ChfRayPos,ChfRayPosSv,th,3)
	ELSE
          ! Point source
	  ! Trace a differential chief ray rotated wrt xGrid
          CALL DSMPROD(th,xGrid,5d-6,3)
          CALL RotChfRay(th,ChfRayDir,ChfRayPos,ChfRayPosSv)	  	    
	END IF
	RayDir(1:3,1)=ChfRayDir(1:3)
	RayPos(1:3,1)=ChfRayPos(1:3)

	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)
	cr2Pos=CRIncidPos(1:3,iEndElt); cr2Dir=CROutgoDir(1:3,iEndElt)
        !cr2Dir=CRIncidDir(1:3,iEndElt)
	ifTrace=.FALSE.; ifBuild=.FALSE.; ifPropagate=.FALSE.

	! Find crossing point of 2 chief rays
        CALL FindCrossPt(cr1Pos,cr1Dir,cr2Pos,cr2Dir,
     &                   zp,CrsDist,DplnPosSag)

	IF (ABS(zSource).GE.1d10) THEN
          ! Collimated source
          ! Trace a differential chief ray displaced in -yGrid
          CALL DSMPROD(th,-xGrid,1d-5,3)
          CALL DADD(ChfRayPos,ChfRayPosSv,th,3)
        ELSE
          ! Point source
          ! Trace a differential chief ray rotated wrt xGrid
          CALL DSMPROD(th,yGrid,5d-6,3)
          CALL RotChfRay(th,ChfRayDir,ChfRayPos,ChfRayPosSv)
        END IF
        RayDir(1:3,1)=ChfRayDir(1:3)
        RayPos(1:3,1)=ChfRayPos(1:3)

	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)
        cr2Pos=CRIncidPos(1:3,iEndElt); cr2Dir=CROutgoDir(1:3,iEndElt)
	ifTrace=.FALSE.; ifBuild=.FALSE.; ifPropagate=.FALSE.

	CALL FindCrossPt(cr1Pos,cr1Dir,cr2Pos,cr2Dir,
     &                   zp,CrsDist,DplnPosTan)

	DplnPos(1:3)=0.5d0*(DplnPosSag+DplnPosTan)	

        WRITE(*,*)' New Detector Position =',DplnPos(1:3)
	cr1Dir(1:3)=-cr1Dir(1:3)
        WRITE(*,*)' New Detector Normal   =',cr1Dir(1:3)
	if (dpElt>nElt) return  ! dpElt does not exist in Rx

	CALL CACCEPT(ANS,'YES', 'Accept the New Detector Position? ')
        IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	  VptElt(1:3,dpElt)=DplnPos(1:3)
	  RptElt(1:3,dpElt)=DplnPos(1:3)
!	  psiElt(1:3,dpElt)=-cr1Dir(1:3)
 	  psiElt(1:3,dpElt)=cr1Dir(1:3)
     &       *sign(1d0,DDOTC(cr1Dir(1:3),psiElt(1:3,dpElt)))
	  If (dpElt>2 .and. (LnkElt(dpElt-2)==dpElt)) Then
	    ! Elt (dpElt-2) is focal return surface
	    VptElt(1:3,dpElt-2)=VptElt(1:3,dpElt)
            RptElt(1:3,dpElt-2)=RptElt(1:3,dpElt)
            psiElt(1:3,dpElt-2)=psiElt(1:3,dpElt)
	  End If
   	  WRITE(*,*)' --Detector Element Position ( Elt',dpElt,') Updated!'
	END IF
        WRITE(*,*)' '

	ChfRayPos=ChfRayPosSv; ChfRayDir=ChfRayDirSv
	END SUBROUTINE FDPN


C***********************************************************************

#if 0
	SUBROUTINE FBMSZ(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &    ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &    ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &    Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec,BmSiz)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
        use sourcsub_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
        use traceutil_mod

        IMPLICIT NONE

	LOGICAL :: ifPol,ifGrid,ifPropagate,ifRayPrint,
     &    ifLNegOK,ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &    ifDMBuild,ifDraw,ifLocalCoord

        INTEGER :: npts,ICLEN,iCurWFElt,iCurRayElt,
     &             nBadRays,iStartElt,iEndElt,iObsOpt,
     &             DrawEltVec(mDrawElt,mDrawRay)

        SREAL :: Cmatrix(7,mCm,bRay),
     &           DrawRayVec(2,mDrawElt,mDrawRay)

        REAL*8 :: Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
		  xObs_loc(3),yObs_loc(3),zObs_loc(3)

	If (iEndElt==0) Then
	  ! Is source, return source aperture size
	  BmSiz=Aperture
	  return 
	Else
	  If (.not.LxObs(iEndElt)) Then
	    ! need to define xObs (arbitary) first
	    Call CalcOrthoVec(psiElt(:,iEndElt),xObs_loc) 
	    Call DUNITIZE(xObs_loc)
	  Else
	    xObs_loc=xObs(1:3,iEndElt)
	  End If
	  zObs_loc=psiElt(1:3,iEndElt)
	  CALL DXPROD(yObs_loc,zObs_loc,xObs_loc)

	End If	
	END SUBROUTINE FBMSZ
#endif

C***********************************************************************

! The 'STOP' subroutine centers beam at a specifed point on a specified
! surface, by repositioning the source in XY plane. Source is NOT rotated
! in this beam center operation, as evidenced in the routine where only
! sensitivity wrt ChfRayPos in X and Y are computed. -jzlou

! This routine is called for both 'CENTER' command and 'STOP' command.
! For 'STOP' command, it is called ONLY WHEN the stop is speicied on an 
! element surface. 

	SUBROUTINE STOP(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,Cmatrix,DrawRayVec,
     &	  DrawEltVec,iLocal,xTarg,k,sqErr,ifBeamCentered)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifOPDCalc,ifRayPrint,ifLNegOK,
     &	  ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),ifDMBuild,
     &	  ifDraw,ifBeamCentered,ifLocalCoord

	INTEGER k,npts,ICLEN,iCurWFElt,iCurRayElt,
     &	        nBadRays,iStartElt,iEndElt,iObsOpt,iStartRay,
     &	        iEndRay,nDrawElt,
     &	        DrawEltVec(mDrawElt,mDrawRay),iLocal

	SREAL Cmatrix(7,mCm,bRay),
     &	       DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &	 xDraw(3),yDraw(3),
     &	 DV1(9),DV2(9),S1,
     &	 xTarg(3),xInt(2),xInt0(2),thx,thy,dthPt,ddxPt(2),
     &	 xLocal(3),yLocal(3),zLocal(3),
     &	 newChfRayDir(3),newChfRayPos(3),oldChfRayDir(3),oldChfRayPos(3),
     &	 crLocalPos(3),dxPt(2),sqErr,GPt(2,2),GPtInv(2,2),
     &	 del(3),errTol

C         ***Center Beam in System Stop***

#if 0
	IF (ABS(zSource).LE.1d10) THEN
	  ! With point source (source at finite distance), 
	  ! enforce ChfRay thru STOP position at element by tip/tilt 
          ! chief ray in source plane. -jzlou
#include "stop_set.inc"
	  return 
	End If
#endif

c Note: still to be done...
c 1. Replace "enter stop position" querys with auto calc
c 2. Include stop kinematics in linear source perturbation models
c 3. Include local source coordinates in linear source models

	dthPt=1d-3
	iStartRay=1
	iEndRay=1
	iStartElt=0
	iEndElt=iLocal  ! = id of STOP element
	nMinPts=3
	errTol=1d-22
	CALL DEQUATE(oldChfRayDir,ChfRayDir,3)
	CALL DEQUATE(oldChfRayPos,ChfRayPos,3)
	CALL DEQUATE(newChfRayPos,ChfRayPos,3)

C Find nominal incidence point and local coords
C (Need to define local frame axes (xLocal and yLocal) at
C  STOP element. -jzlou) 
	CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	 ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &	 ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,
     &	 dxSource,nBadRays,StartIndRef,StartExtinc,
     &	 Cmatrix,DrawRayVec,DrawEltVec,
     &	 iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)

	IF (.NOT.ifLocalCoord) THEN
	  WRITE(*,*)'Local coordinates could not be computed'
	  ifBeamCentered=.FALSE.
	END IF

C  Nominal (unperturbed) incidence pt in local coords
	CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
	CALL DSUB(DV1,crLocalPos,VptElt(1,iLocal),3)
	xInt(1)=DDOTC(xLocal,DV1)
	xInt(2)=DDOTC(yLocal,DV1)
	xInt0(1:2)=xInt(1:2)

C  Loop until desired pt is matched (max 100)
	DO 388 k=1,100

C    Check for convergence
	  CALL DSUB(dxPt,xInt,xTarg,2)  ! xTarg is STOP offset from VptElt in
                             ! local frame, also is target position in local frame
	  sqErr=dxPt(1)*dxPt(1)+dxPt(2)*dxPt(2)
	  IF (sqErr.LE.errTol) GO TO 386

C    Compute x-gradient
	  CALL DSMPROD(del,xGrid,dthPt,3) ! dx projected in xGrid
	  CALL DADD(ChfRayPos,newChfRayPos,del,3) ! ChfRayPos perturbed along xGrid

	  ! Trace chief ray from source to STOP element (=iEndElt)
	  CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	   ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	   ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	   iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	   nBadRays,StartIndRef,StartExtinc,
     &	   Cmatrix,DrawRayVec,DrawEltVec)

	  CALL DSUB(DV1,RayPos,crLocalPos,3) ! diff of ChfRayPos footprint at STOP 
	  ddxPt(1)=DDOTC(DV1,xLocal)
	  ddxPt(2)=DDOTC(DV1,yLocal)
	  GPt(1,1)=ddxPt(1)/dthPt
	  GPt(2,1)=ddxPt(2)/dthPt

C    Compute y-gradient	  
	  CALL DSMPROD(del,yGrid,dthPt,3)
	  CALL DADD(ChfRayPos,newChfRayPos,del,3)

	  CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	   ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	   ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	   iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	   nBadRays,StartIndRef,StartExtinc,
     &	   Cmatrix,DrawRayVec,DrawEltVec)

	  CALL DSUB(DV1,RayPos,crLocalPos,3)
	  ddxPt(1)=DDOTC(DV1,xLocal)
	  ddxPt(2)=DDOTC(DV1,yLocal)
	  GPt(1,2)=ddxPt(1)/dthPt
	  GPt(2,2)=ddxPt(2)/dthPt

C    Compute next step
	  S1=1d0/(GPt(1,1)*GPt(2,2)-GPt(1,2)*GPt(2,1))
	  GPtInv(1,1)=S1*GPt(2,2)
	  GPtInv(1,2)=-S1*GPt(1,2)
	  GPtInv(2,1)=-S1*GPt(2,1)
	  GPtInv(2,2)=S1*GPt(1,1)
	  thx=-(GPtInv(1,1)*dxPt(1)+GPtInv(1,2)*dxPt(2))
	  thy=-(GPtInv(2,1)*dxPt(1)+GPtInv(2,2)*dxPt(2))
	  CALL DSMPROD(DV1,xGrid,thx,3)
	  CALL DSMPROD(DV2,yGrid,thy,3)
	  CALL DADD(del,DV1,DV2,3)
	  CALL DADD(ChfRayPos,newChfRayPos,del,3)

	  CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	   ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	   ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	   iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	   nBadRays,StartIndRef,StartExtinc,
     &	   Cmatrix,DrawRayVec,DrawEltVec)

	  IF (nBadRays.NE.0) THEN
	    WRITE(*,*)'Chief ray becomes undefined during search'
	    ifBeamCentered=.FALSE.
	  END IF

	  CALL DEQUATE(newChfRayPos,ChfRayPos,3)

	  CALL DEQUATE(crLocalPos,RayPos,3)
	  CALL DSUB(DV1,RayPos,VptElt(1,iLocal),3)
	  xInt(1)=DDOTC(xLocal,DV1)
	  xInt(2)=DDOTC(yLocal,DV1)
	  S1=SQRT(thx*thx+thy*thy)
	  IF ((dthPt.GT.S1).AND.(dthPt.GT.5d-12)) dthPt=0.2*dthPt
 388	CONTINUE

C   Here if did not converge

	WRITE(*,*)' WARNING: Stop calculation failed to converge!'
	CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
	ifBeamCentered=.FALSE.
	
 386	CONTINUE

C   Here if converged
	ifBeamCentered=.TRUE.
	k=k-1
	END SUBROUTINE STOP

C***********************************************************************

	SUBROUTINE LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	  ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &	  ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,
     &	  dxSource,nBadRays,StartIndRef,StartExtinc,
     &	  Cmatrix,DrawRayVec,DrawEltVec,
     &	  iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifRayPrint,ifLNegOK,
     &	  ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &    ifDMBuild,ifDraw,ifLocalCoord

	INTEGER npts,iCurWFElt,iCurRayElt,
     &	  nBadRays,iStartElt,iEndElt,iObsOpt,iStartRay,
     &	  iEndRay,DrawEltVec(mDrawElt,mDrawRay),
     &	  iLocal

	SREAL Cmatrix(7,mCm,bRay),
     &	  DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &	  delMag,thMag,del(3),th(3),xLocal(3),yLocal(3),zLocal(3),
     &	  DMAG1,DMAG2,oldChfRayDir(3),oldChfRayPos(3),crLocalPos(3),
     &	  v1(3),v2(3)

	delMag=1d-6
	thMag=1d-6
	iStartRay=1
	iEndRay=1
	iStartElt=0
	iEndElt=iLocal
	CALL DEQUATE(oldChfRayDir,ChfRayDir,3)
	CALL DEQUATE(oldChfRayPos,ChfRayPos,3)

C  Find nominal incidence point
        !ifRayPrint=.TRUE. ! test only
	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	  ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	  ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	  iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,
     &	  Cmatrix,DrawRayVec,DrawEltVec)
	IF (nBadRays.NE.0) THEN
	  WRITE(*,*)'Nominal chief ray becomes undefined'
	  ifLocalCoord=.FALSE.
	END IF

	CALL DEQUATE(crLocalPos,RayPos,3)

C  Find stop elt X-axis
C    First try translating differential ray:
	CALL DSMPROD(del,xGrid,delMag,3)
	CALL DADD(ChfRayPos,oldChfRayPos,del,3)

	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	  ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	  ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	  iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,
     &	  Cmatrix,DrawRayVec,DrawEltVec)

	CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
	IF (nBadRays.NE.0) THEN
	  WRITE(*,*)'Differential ray (x-t) becomes undefined'
	  ifLocalCoord=.FALSE.
	END IF
	CALL DSUB(xLocal,RayPos,crLocalPos,3)
	CALL DUNITVEC(v1,xLocal,DMAG1,3)

C    Didn't work, try rotating differential ray:
	IF (DMAG1.LT.1d-3) THEN
	  CALL DSMPROD(th,xGrid,thMag,3)
	  CALL DADD(ChfRayDir,oldChfRayDir,th,3)
	  CALL DUNITIZE(ChfRayDir)
	  CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	    nBadRays,StartIndRef,StartExtinc,
     &	    Cmatrix,DrawRayVec,DrawEltVec)
	  CALL DEQUATE(ChfRayDir,oldChfRayDir,3)
	  IF (nBadRays.NE.0) THEN
	    WRITE(*,*)'Differential ray (x-r) becomes undefined'
	    ifLocalCoord=.FALSE.
	  END IF
	  CALL DSUB(xLocal,RayPos,crLocalPos,3)
	  CALL DUNITVEC(v2,xLocal,DMAG2,3)
	  IF (DMAG2.GT.DMAG1) THEN
	    CALL DEQUATE(xLocal,v2,3)
	  ELSE
	    CALL DEQUATE(xLocal,v1,3)
	  END IF
	ELSE
	  CALL DEQUATE(xLocal,v1,3)
	END IF

C  Find stop elt Y-axis
C    First try translating differential ray:
	CALL DSMPROD(del,yGrid,delMag,3)
	CALL DADD(ChfRayPos,oldChfRayPos,del,3)
	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	  ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	  ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	  iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,
     &	  Cmatrix,DrawRayVec,DrawEltVec)
	CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
	IF (nBadRays.NE.0) THEN
	  WRITE(*,*)'Tilted differential ray becomes undefined'
	  ifLocalCoord=.FALSE.
	END IF
	CALL DSUB(yLocal,RayPos,crLocalPos,3)
	CALL DUNITVEC(v1,yLocal,DMAG1,3)

C    Didn't work, try rotating differential ray:
	IF (DMAG1.LT.1d-3) THEN
	  CALL DSMPROD(th,yGrid,thMag,3)
	  CALL DADD(ChfRayDir,oldChfRayDir,th,3)
	  CALL DUNITIZE(ChfRayDir)
	  CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	    nBadRays,StartIndRef,StartExtinc,
     &	    Cmatrix,DrawRayVec,DrawEltVec)
	  CALL DEQUATE(ChfRayDir,oldChfRayDir,3)
	  IF (nBadRays.NE.0) THEN
	    WRITE(*,*)'Rotated differential ray becomes undefined'
	    ifLocalCoord=.FALSE.
	  END IF
	  CALL DUNITVEC(v2,yLocal,DMAG2,3)
	  IF (DMAG2.GT.DMAG1) THEN
	    CALL DEQUATE(yLocal,v2,3)
	  ELSE
	    CALL DEQUATE(yLocal,v1,3)
	  END IF
	ELSE
	  CALL DEQUATE(yLocal,v1,3)
	END IF

C  Successful termination
	ifLocalCoord=.TRUE.
	CALL DXPROD(zLocal,xLocal,yLocal)
	CALL DUNITIZE(zLocal)
	CALL DXPROD(yLocal,zLocal,xLocal)
	END SUBROUTINE LocalCoord

C***********************************************************************

	SUBROUTINE CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &	  ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,
     &    ifRegrid,ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &	  iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,
     &	  Cmatrix,DrawRayVec,DrawEltVec)

        use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod

        use sourcsub_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifRayPrint,ifLNegOK,
     &	  ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &    ifDMBuild,ifDraw

	INTEGER npts,iCurWFElt,iCurRayElt,
     &	        nBadRays,iStartElt,iEndElt,iObsOpt,iStartRay,
     &	        iEndRay,nDrawElt,DrawEltVec(mDrawElt,mDrawRay)

	SREAL Cmatrix(7,mCm,bRay),
     &	      DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &	       xDraw(3),yDraw(3)

	ifDraw=.FALSE.
	iStartElt=0
	iStartRay=1
	iEndRay=1
	nMinPts=3

	IF (ABS(zSource).GE.1d10) THEN
	  CALL ColSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	                 ifPropagate)
	ELSE
	  CALL PtSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	                ifPropagate)
	END IF
	CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &		ifLNegOK,Cmatrix,
     &	 iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	 ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	 StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
	END SUBROUTINE CRTrace

C***********************************************************************

	SUBROUTINE SRTrace(xoff,yoff,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec,rayS)

	use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
	use traceutil_mod, only : RayStateStr,thisRay

        use sourcsub_mod

        IMPLICIT NONE

        LOGICAL ifPol,ifGrid,ifPropagate,ifRayPrint,ifLNegOK,
     &   ifDraw,ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &   ifDMBuild
	INTEGER iCurWFElt,iCurRayElt,nBadRays,ie,
     &          iStartElt,iEndElt,iObsOpt,iStartRay,iEndRay,
     &          nDrawElt,DrawEltVec(mDrawElt,mDrawRay),ICLEN
        SREAL Cmatrix(7,mCm,bRay),
     &        DrawRayVec(2,mDrawElt,mDrawRay)
        REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &         xDraw(3),yDraw(3),xoff,yoff
	Type(RayStateStr), pointer :: rayS

	ifDraw=.FALSE.; nMinPts=3
        iStartElt=0 
	! ray 1 is still chief ray; 
        ! ray 2 is the selected ray to trace if (xoff,yoff)/=(0d0,0d0)
	If (xoff==0d0 .and. yoff==0d0) Then
	  iStartRay=1; iEndRay=1
        Else
          iStartRay=2; iEndRay=2
	End If

	! Test only
	!iEndElt=nElt; xoff=0.3; yoff=0.2

	IF (ABS(zSource).GE.1d10) THEN
	  Call SetColSrcRay(nMinPts,Dicr,Djcr,dxSource,
     &                      xoff,yoff,ifPol,ifPropagate)
	ELSE
	  Call SetPtSrcRay(nMinPts,Dicr,Djcr,dxSource,
     &                     xoff,yoff,ifPol,ifPropagate)
	END If 

	CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &   ifLNegOK,Cmatrix,
     &   iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &   ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &   StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
	rayS=>thisRay

#ifdef CSMACOS
        ! Dumpe ray state to a file
        Open(12,file=filnam(1:ICLEN(filnam)),status="replace")
        Do ie=1,iEndElt
          Write(12,*) 
     &      rayS%InciPos(1:3,ie),rayS%InciDir(1:3,ie),
     &      rayS%InciSrfN(1:3,ie),
     &      rayS%PrevIndRef(ie),rayS%CurIndRef(ie)
        End Do
        close(12)
#endif
	END SUBROUTINE SRTrace

C***********************************************************************

	SUBROUTINE CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,
     &    ifDraw,ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &    iStartRay,iEndRay,
     &	  RayMap,SegMap,ifRayPrint,ifIntSrfSet,L1,nBadRays,
     &	  StartIndRef,StartExtinc,DrawRay,DrawElt,nDrawElt,
     &    xDraw,yDraw)

! L1 arg = LRayPass array defined in elt_mod.F -jzlou

c History
c  24Jul97 loh:  added ZernTypeL stuff

	use param_mod
        use elt_mod
	use elemsub_mod
	use traceutil_mod
	use dopt_mod, only : beamOpt
	use macos_debug

	IMPLICIT NONE

	! Arguments 
	LOGICAL ifBuild,ifDMBuild,LROKdummy,ifRayPrint,
     &      LRayTrans,ifPol,
     &	    ifLNegOK,L1(mRay),ifIntSrfSet(mElt),ifDraw
! ifIntSrfSet: related to interpolated surface -jzlou

	Integer :: iStartElt,iEndElt,iStartRay,iEndRay,
     &       	   nBadRays,nDrawElt,RayMap(mpts,mpts),
     &             SegMap(mpts,mpts),DrawElt(mDrawElt)

	SREAL :: DrawRay(2,mDrawElt),Cmatrix(7,mCm,bRay)
	Real*8 :: StartIndRef,StartExtinc,xDraw(3),yDraw(3)

        ! Local variables
	Logical ifRayElt(0:mElt),ifIntSrfOK(mElt),
     &	        ifTransmit,ifLNsrf,LRTdummy,ifRayToSeg,
     &	        ifUseBaseSrf,ifReturn

	INTEGER i,j,iElt,iRay,iNode,nEm1,nEm2,RayIndex(2,mRay),
     &     PrevNonSeg,NextNonSeg,itst,imin,iFirstNonSeq,
     &     iLastNonSeq,IERROR

        Integer :: NextElt(0:mElt),PrevElt(0:mElt),iminPrev,
     &      jGridSrf(mElt),
     &	    iSecondElt,AdjRay(4),LastElt,iEltLoop,NextLast,
     &	    iIntSrf,NDP3,jElt,iNextSeq,iNextElt,CrntNSCnt

	SREAL :: dxndxi(7,7,0:mElt),dxndui(7,6,mElt),spot(3)

	REAL*8 Lmin,na,nb,RayPosTst(3),RayLTst,PrevDir(3),
     &         PrevPos(3),
     &	       xa,xb,CurIndRef,PrevIndRef,NextIndRef,
     &	       Intxyz,CurExtinc,PrevExtinc,NextExtinc,CWU

	DATA ifUseBaseSrf/.FALSE./

 115    FORMAT(1P,' Ray ',i4,' segment from Element ',i4,' (',a12,
     &  ') to Element ',i4,' (',a12,'):'/
     &  '   Starting point=',3D22.14/'        End point=',3D22.14/
     &  '        Direction=',3D22.14/'   Segment Length=',D22.14/
     &	'     Total Length=',D22.14)
 116    FORMAT(' Ray',i6,' has aperture coordinates',i5,' and',i5,
     &	'.'/' WF coordinates are ',i4,' and ',i4,'.'/
     &	' Adjacent rays are ',3(i4,', '),i4,'.')
 126    FORMAT(' Ray',i6,' has aperture coordinates',i5,' and',i5,
     &  '.'/' WF coordinates are ',i4,' and ',i4,'.')
 117    FORMAT(' Ray',i6,' is the chief ray.')
 118	FORMAT(1P,'       Index(i-1)=',2d17.9/
     &	'         Index(i)=',2d17.9)
 119	FORMAT(1P,15X,'Ex=',2d17.9/15X,'Ey=',2d17.9/
     &	15X,'Ez=',2d17.9/13X,'Ixyz=',d17.9)
 122	FORMAT(/' Ray',i6,' becomes undefined after element',i3)
 123	FORMAT(/' Ray',i6,' is obscured by element',i3)
 124	FORMAT(' Interpolated surface',i4,' traced using base ',
     &	'conicoid only.'/
     &  ' Use SINT command to set interpolation data.')

C  0. Initialize parameters

	PrevElt=0  ! default, added by jzlou to avoid code crash in building partials

	tAsphCoef=maxAsphCoef  ! added by jzlou

	!ChfRayHist%keepHist=1
	If (ChfRayHist%keepHist) Then
	  ! Need to keep chief ray history during ray-tracing
	  If (.not. allocated(ChfRayHist%SrfId)) Then
	    allocate(ChfRayHist%SrfId(nElt),
     &               ChfRayHist%RayPosAtSrf(3,nElt),
     &               ChfRayHist%RayDirAtSrf(3,nElt))
          End If
	  ChfRayHist%tot_srfs=0
	End If

	nBadRays=0
	EltID(0)=0
	iSecondElt=iStartElt+1   ! = 0 when run 'opd' commnd, -jzlou

C	Find start and end points for segmented elements

	PrevNonSeg=0
	iElt=1
	DO WHILE (iElt.LE.nElt)
	  IF (EltID(iElt).NE.11) THEN  ! iElt is non-segment -jzlou
	    NextElt(PrevNonSeg)=iElt
	    PrevElt(iElt)=PrevNonSeg
	    PrevNonSeg=iElt
	  ELSE  ! iElt is segment -jzlou
	    PrevElt(iElt)=PrevNonSeg
	    NextNonSeg=iElt+1
	    DO WHILE (EltID(NextNonSeg).EQ.11)
	      NextNonSeg=NextNonSeg+1
	    END DO
	    NextElt(iElt)=NextNonSeg
	  END IF
	  iElt=iElt+1
	END DO

	! For each iElt, PrevElt(iElt) is previous non-segment, NextElt(iElt) 
        ! is next non-segment. -jzlou

C	Initialize input beam if first element is source

	IF (iStartElt.EQ.0) THEN
	  StartIndRef=IndRef(0)
	  StartExtinc=Extinc(0)
	  DO 3 iRay=1,nRay
	    CumRayL(iRay)=0d0
	    LRayOK(iRay)=.TRUE.
	    L1(iRay)=.TRUE.
 3        CONTINUE
	  PrevNonSeg=0
	ELSE
	  PrevNonSeg=iStartElt
	END IF

C	Set wavelength in base units; set index of refraction for Glass;
C	convert Zernike coefficients to monomials; initialize LEltOK; 
C	check interpolated surfaces; compute obscuration coords

	CWU=CWM*1e6  ! scaling wavelength to um, only for glass table -jzlou 
	CWB=CWM/CBM
     	WaveBU=CWB*Wavelen  
        ! WaveBU is the value of wavelength converted to base unit,
        ! for later use, -jzlou.

	! L1(1:nRay), once set to .FALSE. (meaining blocked), will not become
        ! .TRUE. at subsequent elements, in the loop below. -jzlou

	iFirstNonSeq=0
	iLastNonSeq=nElt
	DO 17 iElt=iSecondElt,iEndElt

C Set index
	  IF (LGlass(iElt)) THEN

	    !print*,'** CTRACE, iElt=', iElt
	    !print*,'** GlassType(iElt)=', GlassType(iElt)

	    IndRef(iElt)=getIndex(CWU*Wavelen,GlassCoef(1,iElt),
     &                            GlassType(iElt))
	  END IF
C Set Zerns
	  IF (SrfType(iElt).EQ.8) THEN
	    IF (ZernTypeL(iElt).EQ.1 .OR.
     &          ZernTypeL(iElt).EQ.4) THEN
C Malacara's
	      CALL ZerntoMon1(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)
	    ELSE IF (ZernTypeL(iElt).EQ.2 .OR.
     &              ZernTypeL(iElt).EQ.5) THEN
C Noll's
	      CALL ZerntoMon2(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)
	    ELSE IF (ZernTypeL(iElt).EQ.3 .OR.
     &               ZernTypeL(iElt).EQ.6) THEN
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
     &	      (ApType(iElt).GT.0)) THEN
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
 17	CONTINUE

C	Initialize DM partials data
C	i=7*mDMdof*bRay
C	IF (ifDMBuild) CALL RZERO(CmDM,i)

C	Begin iRay loop to trace each ray and compute its partials

	DO 14 iRay=iStartRay,iEndRay
	  crnt_traced_ray=iRay; iRay_debug=iRay

	  IF (ifRayPrint) THEN
	    IF (iRay.EQ.1) THEN
	      if (.not.LBWKnoPrint) WRITE(*,117) iRay
	    ELSE
	      i=RayIndex(1,iRay)
	      j=RayIndex(2,iRay)
	      IF ((i.GT.1).AND.(i.LT.mpts)) THEN
	        AdjRay(1)=RayMap((i-1),j)
	        AdjRay(2)=RayMap((i+1),j)
	      ELSE IF (i.EQ.1) THEN
	        AdjRay(1)=0
	        AdjRay(2)=RayMap((i+1),j)
	      ELSE IF (i.EQ.mpts) THEN
	        AdjRay(1)=RayMap((i-1),j)
	        AdjRay(2)=0
	      END IF
	      IF ((j.GT.1).AND.(j.LT.mpts)) THEN
	        AdjRay(3)=RayMap(i,(j-1))
	        AdjRay(4)=RayMap(i,(j+1))
	      ELSE IF (j.EQ.1) THEN
	        AdjRay(3)=0
	        AdjRay(4)=RayMap(i,(j+1))
	      ELSE IF (j.EQ.mpts) THEN
	        AdjRay(3)=RayMap(i,(j-1))
	        AdjRay(4)=0
	      END IF
	      If (.not.ifSRT) Then
	        WRITE(*,116) iRay,RayIndex(1,iRay),RayIndex(2,iRay),
     &	          RayWFMap(1,iRay),RayWFMap(2,iRay),(AdjRay(i),i=1,4)
	      Else
	        WRITE(*,126) iRay,RayIndex(1,iRay),RayIndex(2,iRay),
     &            RayWFMap(1,iRay),RayWFMap(2,iRay)
	      End If
	    END IF
	  END IF

	  IF (ifDraw) THEN
	    nDrawElt=1
	    DrawElt(nDrawElt)=iStartElt
	    DrawRay(1,nDrawElt)=DDOTC(xDraw,RayPos(1,iRay))
	    DrawRay(2,nDrawElt)=DDOTC(yDraw,RayPos(1,iRay))
	  ELSE
	    nDrawElt=0
	  END IF

c	  PrevElt(0)=0
c	  NextElt(0)=1
	  ifRayElt(0)=.TRUE.
	  LRayTrans=.TRUE.
	  ifLNsrf=.FALSE.


C       Loop through each element

	  CurIndRef=StartIndRef

	  iElt=iSecondElt
	  PrevIndRef=CurIndRef
	  CurExtinc=StartExtinc
	  PrevExtinc=CurExtinc
	  ifReturn=.FALSE.

#if 1
	  If (IsVarIndRefElt(iStartElt)) Then
            ! StartElt has variable IndRef
	    CurIndRef=VarStrtIndRef(iRay)
	  Else
	    CurIndRef=StartIndRef
          End If
#endif

	  DO 10 iEltLoop=iSecondElt,iEndElt
            iElt_debug=iEltLoop

	    IF (iElt.GT.iEndElt) GO TO 10
	    IF (iElt.LT.nElt) THEN
	      iNextElt=iElt+1
	    ELSE
	      iNextElt=nElt
	    END IF
	    IF ((SrfType(iElt).EQ.5).AND.ifIntSrfOK(iElt)) THEN
	      iIntSrf=iEltToIntSrf(iElt)
	    END IF

	    IF (ifRayPrint .OR. beamOpt%ifOptBeamDir 
     &          .OR. beamOpt%ifOptBeamPos 
     &          .OR. beamOpt%ifOptBeamRefRayDir
     &          .OR. iRay==1) THEN
	      CALL DEQUATE(PrevDir,RayDir(1,iRay),3)
	      CALL DEQUATE(PrevPos,RayPos(1,iRay),3)
	    END IF

	    IF (ifLNegOK) THEN
	      IF ((EltID(iElt).EQ.3).OR.(EltID(iElt).EQ.7).OR.
     &	        (EltID(iElt).EQ.8).OR.(EltID(PrevNonSeg).EQ.3).OR.
     &	        (EltID(PrevNonSeg).EQ.7).OR.(EltID(PrevNonSeg).EQ.8)
     &          .OR.(EltID(iElt).EQ.10)) THEN ! added LensArray, jzlou 
	        ifLNsrf=.TRUE.
	      ELSE
	        ifLNsrf=.FALSE.
	      END IF
	    END IF

#if 0 
            ! added by jzlou 10/2009
	    If ((EltID(PrevNonSeg).EQ.3).AND.
     &          ABS(KrElt(PrevNonSeg))<1d22) Then
	      ! Doesn't allow negative L when the previous reference 
              ! surface is 'powered'
              ifLNsrf=.FALSE.
            End If  
#endif	    
 
	    crnt_Elt=iElt  ! crnt_Elt, defined in traceutil_mod, is shared 
                           ! with ChkRayTrans function in elemsub.F

C	Trace rays using appropriate element routines
C						! Reflector
	    IF (EltID(iElt).EQ.1) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt
	      CALL Reflector(
     &          ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	        LEltObs(iElt),LRayOK(iRay),LRayTrans,KcElt(iElt),
     &          KrElt(iElt),AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	        AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	        UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &	        GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	        lMon(iElt),pMon(1,iElt),
     &	        xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	        CurIndRef,IndRef(iElt),
     &	        CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	        XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	        IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &          EltCoat(iElt),EltCoatThk(:,iElt),
     &          IndRefArr(:,iElt),ExtincArr(:,iElt))

 	      if (.false. .and. iEltLoop==21) then
 	        print*,'*** Elt 21 Aft Reflector(),LRayOK(iRay)=',
     &			LRayOK(iRay)
       	      end if

	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)

	      if (.false.) then !if (iRay == (iEndRay-45)) then
	        print*,'**After Reflector RayDir(:,iRay) = ', 
     &            RayDir(:,iRay)
	        !print*,'==>iSecondElt,iEndElt = ',iSecondElt,iEndElt
	        print*,'===> CumRayL(iEndRay-45)=', 
     &            CumRayL(iEndRay-45),
     &            ',RayL(iEndRay-45)=',RayL(iEndRay-45),
     &            'SrfType(iElt)=',SrfType(iElt),
     &            'iElt = ', iElt
	      end if

	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.
C						! Focal Plane
	    ELSE IF (EltID(iElt).EQ.2) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt

c	      IF ((EltID(iElt-1).EQ.8).AND.
c     &	      (EltID(iElt-2).NE.8)) THEN
c	        ifReturn=.TRUE.
c	      ELSE
c	        ifReturn=.FALSE.
c	      END IF

	      ifReturn=.FALSE.

	      CALL FocalPln(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &	        LRayOK(iRay),LRayTrans,psiElt(1,iElt),
     &	        VptElt(1,iElt),RptElt(1,iElt),RayPos(1,iRay),
     &	        RayDir(1,iRay),RayE(1,iRay),CurIndRef,
     &	        IndRef(iElt),CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),
     &          ObsType(1,iElt),ObsVec(1,1,iElt),
     &	        xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),
     &	        RayE(1,iRay),spot,
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt))

	      IF (ifReturn) RayL(iRay)=-RayL(iRay)
	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.
C						! Reference Surface
c	    ELSE IF ((EltID(iElt).EQ.3).AND.
c     &	    (ifRayPrint.OR.ifDraw.OR.(iElt.EQ.iEndElt))) THEN
	    ELSE IF (EltID(iElt).EQ.3) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt

	      CALL RefSrf(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &	        LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	        CurIndRef,IndRef(iElt),
     &	        CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),
     &          zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt))

c	      IndRef(iElt)=IndRef(PrevElt(iElt))
	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      IF (ifReturn) RayL(iRay)=-RayL(iRay)
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.
C						! Reflector Segment
	    ELSE IF (EltID(iElt).EQ.11) THEN
	      ifRayToSeg=.FALSE.
	      DO i=1,nRayToSeg(iRay)
	        IF (RayToSegMap(iRay,i).EQ.EltToSegMap(iElt))
     &	          ifRayToSeg=.TRUE.
	      END DO
	      IF (ifRayToSeg) THEN
c		NextNonSeg=PrevNonSeg+nSeg
c		PrevElt(iElt)=PrevNonSeg
c		NextElt(PrevNonSeg)=iElt
c		NextElt(iElt)=NextNonSeg
		PrevNonSeg=iElt
		NextNonSeg=NextElt(iElt)
		PrevElt(NextNonSeg)=iElt
	        CALL Reflector(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	          LEltObs(iElt),
     &	          LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	          AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	          AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	          UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &	          GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &	          psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	          RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	          lMon(iElt),pMon(1,iElt),
     &	          xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	          CurIndRef,IndRef(iElt),
     &	          CurExtinc,Extinc(iElt),WaveBU,
     &	          ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &            dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	          XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	          IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &		  EltCoat(iElt),EltCoatThk(:,iElt),
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
C				! Holographic Optical Element (HOE)
	    ELSE IF (EltID(iElt).EQ.4) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt
	      CALL HOE(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	       LEltObs(iElt),
     &	       LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	       AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	       AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	       UDSrfParam(1,iElt),
     &	       psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	       lMon(iElt),pMon(1,iElt),
     &	       xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	       CurIndRef,IndRef(iElt),
     &	       CurExtinc,Extinc(iElt),WaveBU,
     &	       ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	       ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &         dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	       h1HOE(1,iElt),h2HOE(1,iElt),OrderHOE(iElt),
     &	       CWB*WaveHOE(iElt),
     &	       XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	       IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.
C						! Reflective Grating
	    ELSE IF (EltID(iElt).EQ.5) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt
	      CALL Grating(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	       LEltObs(iElt),
     &	       LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	       AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	       AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	       UDSrfParam(1,iElt),
     &	       psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	       lMon(iElt),pMon(1,iElt),
     &	       xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	       CurIndRef,IndRef(iElt),
     &	       CurExtinc,Extinc(iElt),WaveBU,
     &	       ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	       ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &         dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	       OrderHOE(iElt),RuleWidth(iElt),h1HOE(1,iElt),
     &	       XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	       IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
	       IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	       CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.
C						! Transmissive Grating
	    ELSE IF (EltID(iElt).EQ.13) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt
	      CALL TrGrating(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	       LEltObs(iElt),
     &	       LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	       AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	       AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	       UDSrfParam(1,iElt),
     &	       psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	       lMon(iElt),pMon(1,iElt),
     &	       xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	       CurIndRef,IndRef(iElt),
     &	       CurExtinc,Extinc(iElt),WaveBU,
     &	       ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	       ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &         dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	       OrderHOE(iElt),RuleWidth(iElt),h1HOE(1,iElt),
     &	       XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	       IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      CurIndRef=IndRef(iElt)
	      CurExtinc=Extinc(iElt)
	      ifRayElt(iElt)=.TRUE.
 
C                                ! Diffractive Optical Element (DOE)
            ELSE IF (EltID(iElt).EQ.17) THEN
c             PrevElt(iElt)=PrevNonSeg
c             NextElt(PrevNonSeg)=iElt
              PrevNonSeg=iElt
              CALL DoeTrGrating(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &         LEltObs(iElt),
     &         LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &         AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &         AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &         UDSrfParam(1,iElt),
     &         psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &         RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &         lMon(iElt),pMon(1,iElt),
     &         xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &         CurIndRef,IndRef(iElt),
     &         CurExtinc,Extinc(iElt),WaveBU,
     &         ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &         ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &         RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &         dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &         OrderHOE(iElt),
     &         XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &         IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)
              IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
              CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
              PrevIndRef=CurIndRef
              PrevExtinc=CurExtinc
              CurIndRef=IndRef(iElt)
              CurExtinc=Extinc(iElt)
              ifRayElt(iElt)=.TRUE.


C						! Refractor
	    ELSE IF (EltID(iElt).EQ.6) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt

!	      if (iElt==44) then
!	      print*,'****In CTRACE,B4 Refractor,AsphCoef(1,iElt)=',
!    &          AsphCoef(1:4,iElt),',  iElt=',iElt
!	      end if

	      CALL Refractor(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	        LEltObs(iElt),
     &	        LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	        AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	        AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	        UDSrfParam(1,iElt),nGridMat(iElt),mGridMat,
     &          GridMat(1,1,jGridSrf(iElt)),GridSrfdx(iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	        lMon(iElt),pMon(1,iElt),
     &	        xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	        CurIndRef,IndRef(iElt),LNewIndRef,NewIndRef,
     &	        CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	        XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	        IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP,
     &          EltCoat(iElt),EltCoatThk(:,iElt),
     &          IndRefArr(:,iElt),ExtincArr(:,iElt),
     &		IsVarIndRefElt(iElt),IndRefGridArr(iElt),
     &	        IsGrinRefElt(iElt),GrinRefArr(iElt))
	        ! RayPos(1,iRay) is both input and output variable to above routine 

	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      If (.not. IsGrinRefElt(iElt)) Then
  	        CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      Else
	        ! Graded index lens
	        CumRayL(iRay)=CumRayL(iRay)+RayL(iRay)
	      End If
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc

	      If (LNewIndRef) Then
                CurIndRef=NewIndRef
                VarCurIndRef(iRay)=NewIndRef
              Else
                CurIndRef=IndRef(iElt)
                VarCurIndRef(iRay)=IndRef(iElt)
              End If
	      CurExtinc=Extinc(iElt)
	      ifRayElt(iElt)=.TRUE.
C						! Conic Obscuration
	    ELSE IF ((EltID(iElt).EQ.7).AND.(ifDraw.OR.
     &	        ifRayPrint.OR.(iElt.EQ.iEndElt).OR.(.NOT.ifBuild))) 
     &	     THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt
	      CALL ObsSrf(ifBuild,ifPol,ifLNsrf,LEltObs(iElt),
     &	       LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	       psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	       CurIndRef,IndRef(iElt),CurExtinc,
     &	       Extinc(iElt),WaveBU,
     &	       ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	       ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),
     &         RayL(iRay),RayE(1,iRay),
     &	       dxidxim1(1,1,iElt),dxidui(1,1,iElt))
	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      IndRef(iElt)=CurIndRef
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.
C						! Return element
	    ELSE IF (EltID(iElt).EQ.8) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      ifTransmit=.FALSE.
c	      IF (EltID(PrevNonSeg).EQ.8) THEN
c	        ifReturn=.TRUE.
c	      ELSE
c	        ifReturn=.FALSE.
c	      END IF

	      PrevNonSeg=iElt
	      CALL ReturnSrf(ifBuild,ifPol,ifReturn,ifTransmit,
     &	       ifLNsrf,LEltObs(iElt),LRayOK(iRay),LRayTrans,
     &	       KcElt(iElt),KrElt(iElt),
     &	       psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	       CurIndRef,IndRef(iElt),
     &	       CurExtinc,Extinc(iElt),WaveBU,
     &	       ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	       ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	       RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &         dxidxim1(1,1,iElt),dxidui(1,1,iElt))
	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      IndRef(iElt)=CurIndRef
	      IF (ifReturn) RayL(iRay)=-RayL(iRay)

	      ifReturn=.NOT.ifReturn

	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      ifRayElt(iElt)=.TRUE.

C					! Nonsequential Elements
	    ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.12))
     &	    THEN

	     if (.false.) then
	       print*,'Nonseq elem: iElt = ', iElt
	       print*,'EltID=', EltID(iElt)
	       !print*,'ApVec(1:4,iElt) = ', ApVec(1:4,iElt) 
	     end if


C	Initialize search parameters
	      iFirstNonSeq=iElt
	      iNextSeq=iElt
 706	      CONTINUE
	        IF ((EltID(iNextSeq).EQ.9).OR.
     &	        (EltID(iNextSeq).EQ.12)) THEN
	          ifRayElt(iNextSeq)=.FALSE.
	          iNextSeq=iNextSeq+1
	          GO TO 706
	        END IF
	        iLastNonSeq=iNextSeq-1
	        iminPrev=iElt-1

C	Loop until ray hits next sequential element
	      CrntNSCnt=0
 708	      CONTINUE
	        CrntNSCnt=CrntNSCnt+1
	        IF (ifRayPrint .OR. beamOpt%ifOptBeamDir
     &              .OR. beamOpt%ifOptBeamPos 
     &              .OR. beamOpt%ifOptBeamRefRayDir
     &              .OR. iRay==1) THEN
	          CALL DEQUATE(PrevDir,RayDir(1,iRay),3)
	          CALL DEQUATE(PrevPos,RayPos(1,iRay),3)
	        END IF
C	Find next element as least ray length
	        Lmin=1D22
	        imin=iNextSeq
	        DO 707 itst=iFirstNonSeq,iNextSeq

	          crnt_Elt=itst  ! crnt_Elt, defined in traceutil_mod, is shared 
                                 ! with ChkRayTrans function in elemsub.F

	          LROKdummy=.TRUE.
	          LRTdummy=.TRUE.
c	          IF (ifUseBaseSrf) THEN
c	            CALL FindConSrf(ifBuild,ifLNsrf,LROKdummy,
c    &	            KcElt(itst),KrElt(itst),
c    &	            psiElt(1,itst),VptElt(1,itst),RptElt(1,itst),
c    &	            RayPos(1,iRay),RayDir(1,iRay),
c    &	            RayPosTst,RayLTst)
c	          ELSE
	            CALL FindSrf(ifLNsrf,LEltObs(itst),
     &	             LROKdummy,LRTdummy,KcElt(itst),KrElt(itst),
     &	             AsphCoef(1:tAsphCoef,itst),MonCoef(1,itst),
     &	             AnaCoef(1,itst),UDSrfCoef(1,itst),UDSrfType(itst),
     &	             UDSrfParam(1,itst),nGridMat(itst),mGridMat,
     &	             GridMat(1,1,jGridSrf(itst)),GridSrfdx(itst),
     &	             psiElt(1,itst),VptElt(1,itst),RptElt(1,itst),
     &	             RayPos(1,iRay),RayDir(1,iRay),
     &	             lMon(itst),pMon(1,itst),
     &	             xMon(1,itst),yMon(1,itst),zMon(1,itst),
     &	             ApType(itst),ApVec(1,itst),nObs(itst),
     &	             ObsType(1,itst),ObsVec(1,1,itst),
     &	             xObs(1,itst),yObs(1,itst),zObs(1,itst),
     &	             RayPosTst,RayLTst,SrfType(itst),
     &	             XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	             IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),
     &               mDP,maxAsphCoef)
c	          END IF
C note: RayLTst changed from 1d-12 to 1d-8
	          IF (LROKdummy.AND.(RayLTst.LT.Lmin).AND.
     &	          (RayLTst.GT.1d-8).AND.LRTdummy.AND.
     &	          ((SrfType(iminPrev).LE.2).OR.(iminPrev.NE.itst))) THEN
	            Lmin=RayLTst
	            imin=itst
	          END IF
 707	        CONTINUE

C	Here when element number is determined...
	        IF (Lmin.LT.1d22) THEN
	          IF (imin.NE.iNextSeq) THEN
	            LROKdummy=.TRUE.
	            LRTdummy=.TRUE.
	            IF (EltID(imin).EQ.12) THEN
	              IF ((SrfType(imin).EQ.5).AND.ifIntSrfOK(imin)) 
     &	                iIntSrf=iEltToIntSrf(iElt)

                      crnt_Elt=imin  ! crnt_Elt is shared with
                           	     ! ChkRayTrans function in elemsub.F

	              CALL Reflector(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	               LEltObs(imin),
     &	               LROKdummy,LRTdummy,KcElt(imin),KrElt(imin),
     &	               AsphCoef(1:tAsphCoef,imin),MonCoef(1,imin),
     &	               AnaCoef(1,imin),UDSrfCoef(1,imin),UDSrfType(imin),
     &	               UDSrfParam(1,imin),nGridMat(imin),mGridMat,
     &	               GridMat(1,1,jGridSrf(imin)),GridSrfdx(imin),
     &	               psiElt(1,imin),VptElt(1,imin),RptElt(1,imin),
     &	               RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	               lMon(imin),pMon(1,imin),
     &	               xMon(1,imin),yMon(1,imin),zMon(1,imin),
     &	               CurIndRef,IndRef(imin),
     &	               CurExtinc,Extinc(imin),WaveBU,
     &	               ApType(imin),ApVec(1,imin),nObs(imin),
     &	               ObsType(1,imin),ObsVec(1,1,imin),
     &	               xObs(1,imin),yObs(1,imin),zObs(1,imin),
     &	               RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),
     &	               RayE(1,iRay),dxidxim1(1,1,imin),dxidui(1,1,imin),
     &	               SrfType(imin),
     &	               XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &                 IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),
     &                 mDP,EltCoat(imin),EltCoatThk(:,imin),
     &                 IndRefArr(:,imin),ExtincArr(:,imin))
	              NextIndRef=CurIndRef
	              NextExtinc=CurExtinc
	            ELSE
	              IF (EltID(PrevNonSeg).NE.9) THEN
	                na=CurIndRef
	                nb=IndRef(imin)
	                xa=CurExtinc
	                xb=Extinc(imin)
	              ELSE IF ((ABS(KrElt(imin)).LT.
     &	                        ABS(KrElt(PrevNonSeg)))) THEN
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

	              crnt_Elt=imin  ! crnt_Elt is shared with
                           	     ! ChkRayTrans function in elemsub.F

	              CALL NSRefractor(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	                LEltObs(imin),LROKdummy,LRTdummy,
     &	                KcElt(imin),KrElt(imin),
     &	                AsphCoef(1:tAsphCoef,imin),MonCoef(1,imin),
     &	                AnaCoef(1,imin),UDSrfCoef(1,imin),
     &	                UDSrfType(imin),UDSrfParam(1,imin),
     &                  nGridMat(imin),mGridMat,
     &                  GridMat(1,1,jGridSrf(imin)),GridSrfdx(imin),
     &	                psiElt(1,imin),VptElt(1,imin),RptElt(1,imin),
     &	                RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	                lMon(imin),pMon(1,imin),
     &	                xMon(1,imin),yMon(1,imin),zMon(1,imin),
     &	                na,nb,xa,xb,WaveBU,nCoatElt(imin),
     &	                CoatIndxElt(1,imin),CoatThkElt(1,imin),
     &	                ApType(imin),ApVec(1,imin),nObs(imin),
     &	                ObsType(1,imin),ObsVec(1,1,imin),
     &	                xObs(1,imin),yObs(1,imin),zObs(1,imin),
     &	                RayPos(1,iRay),
     &	                RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &	                dxidxim1(1,1,imin),dxidui(1,1,imin),
     &	                SrfType(imin),
     &	                XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	                IWORK(1,iIntSrf),DWORK(1,iIntSrf),
     &	                NDP(iIntSrf),mDP)
	            END IF

	            IF (.NOT.LRTdummy) THEN
	              LRayTrans=.FALSE.
	              L1(iRay)=.FALSE.
	            END IF

	            IF (.NOT.LROKdummy) THEN
	              LRayOK(iRay)=.FALSE.
	              PrevElt(iFirstNonSeq)=PrevNonSeg
	            ELSE
	              PrevElt(imin)=PrevNonSeg
	              NextElt(PrevNonSeg)=imin
	              PrevNonSeg=imin
	              CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	              PrevIndRef=CurIndRef
	              PrevExtinc=CurExtinc
	              CurIndRef=NextIndRef
	              CurExtinc=NextExtinc
	              ifRayElt(imin)=.TRUE.
	              IF (ifDraw) THEN
	                nDrawElt=nDrawElt+1
	                DrawElt(nDrawElt)=imin
	                DrawRay(1,nDrawElt)=DDOTC(xDraw,RayPos(1,iRay))
	                DrawRay(2,nDrawElt)=DDOTC(yDraw,RayPos(1,iRay))
	              END IF

                      ! Save incident chief ray state:
                      ! direction, intersection point, surf nomal, -jzlou
                      If (iRay==1) Then
                        CRIncidDir(1:3,iElt)=PrevDir(1:3)
                        CROutgoDir(1:3,iElt)=RayDir(1:3,iRay)
                        CRIncidPos(1:3,iElt)=RayPos(1:3,iRay)
                        CRNhat(1:3,iElt)=-NhatP(1:3) ! CRNhat is outward from surface
                      End If

	              IF (ifRayPrint) THEN
	                if (.not.LBWKnoPrint) then
	                  WRITE(*,115) iRay, PrevElt(imin),
     &	                    EltName(PrevElt(imin)),imin,EltName(imin),
     &	                    (PrevPos(i),i=1,3),
     &	                    (RayPos(i,iRay),i=1,3),
     &	                    (PrevDir(i),i=1,3),RayL(iRay),
     &	                    CumRayL(iRay)
	                  WRITE(*,118)PrevIndRef,PrevExtinc,CurIndRef,
     &	                              CurExtinc
	                end if

	                IF (ifPol) THEN
	                  Intxyz=RayE(1,iRay)*DCONJG(RayE(1,iRay))+
     &	                  RayE(2,iRay)*DCONJG(RayE(2,iRay))+
     &	                  RayE(3,iRay)*DCONJG(RayE(3,iRay))
	                  WRITE(*,119)RayE(1,iRay),
     &	                  RayE(2,iRay),RayE(3,iRay),Intxyz
	                END IF
	              END IF

	              If (beamOpt%ifOptBeamDir .and. 
     &                    imin==beamOpt%beamDirElt .and. iRay==1) Then 
                          ! Get chief ray direction at BeamDir Elt, added by jzlou
	                 ChfRayDirAtElt(1:3)=PrevDir(1:3) 
	              End If
	              If (beamOpt%ifOptBeamPos .and.
     &                    imin==beamOpt%beamPosElt .and. iRay==1) Then
                          ! Get chief ray position at beamPosElt, added by jzlou
                         ChfRayPosAtElt(1:3)=RayPos(1:3,iRay)
                      End If
	              If (beamOpt%ifOptBeamRefRayDir .and.
     &                    (imin==beamOpt%beamRefRayElt)) Then
	                CurRefRayDir(1:3)=PrevDir(1:3)
	              End If
	              If (ChfRayHist%keepHist .and. iRay==1) Then
                        ! Keep chief ray history
                        ChfRayHist%tot_srfs=ChfRayHist%tot_srfs+1
                        ChfRayHist%SrfId(ChfRayHist%tot_srfs)=imin
                        ChfRayHist%RayPosAtSrf(1:3,ChfRayHist%tot_srfs)
     &                    =RayPos(1:3,iRay)
                        ChfRayHist%RayDirAtSrf(1:3,ChfRayHist%tot_srfs)
     &                    =PrevDir(1:3)
                      End If
                      !
	            END IF
	          END IF
	        ELSE
	          LRayOK(iRay)=.FALSE.
	          PrevElt(iFirstNonSeq)=PrevNonSeg
	        END IF
	        IF (.NOT.LRayOK(iRay)) THEN
	          nBadRays=nBadRays+1
	          LastElt=PrevElt(iElt)
	          LEltOK(imin)=.FALSE.
	          IF (IRAY.GT.1) THEN
	            RayMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	            SegMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	          END IF
	          IF (ifRayPrint) THEN
	            WRITE(*,122)iRay,LastElt
	          END IF
	          GO TO 14
	        END IF
	        iElt=imin
	        iminPrev=imin

                ! **************************************************************
                ! Adding saving RayPosHist(1:3,jRay,iElt) for nonseq elt? -jzlou
                ! **************************************************************
                If (SaveRayPosHist_Flg .OR. SaveVis3d_dat_Flg) Then
                  LRayOKHist(iRay,iElt) = .false.
                  If (LRayOK(iRay) .and. ifRayElt(iElt)) Then
                    LRayOKHist(iRay,iElt) = .true.
                    RayPosHist(1:3,iRay,iElt)=RayPos(1:3,iRay)
                  End If
                End If

	        If (NSCnt(iElt).gt.0 .and. 
     &              CrntNSCnt.eq.NSCnt(iElt)) Then
	          ! This is for optimizing ray-tracing speed with segmented mirror;
	          ! NSCount should NOT be defined in Rx for CC -jzlou
                  imin=iNextSeq
                End If

	        IF (imin.NE.iNextSeq) THEN
                  PrevElt(iNextSeq)=imin ! In case this imin is the last NS elt hit,
                                         ! 01/06/2009
	          GO TO 708
	        END IF
	      iElt=iLastNonSeq

C					! Lens Array
	    ELSE IF (EltID(iElt).EQ.10) THEN
c	      PrevElt(iElt)=PrevNonSeg
c	      NextElt(PrevNonSeg)=iElt
	      PrevNonSeg=iElt
	      CALL LensArray(ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &	        LEltObs(iElt),
     &	        LRayOK(iRay),LRayTrans,KcElt(iElt),KrElt(iElt),
     &	        AsphCoef(1:tAsphCoef,iElt),MonCoef(1,iElt),
     &	        AnaCoef(1,iElt),UDSrfCoef(1,iElt),UDSrfType(iElt),
     &	        UDSrfParam(1,iElt),
     &	        psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay),
     &	        lMon(iElt),pMon(1,iElt),
     &	        xMon(1,iElt),yMon(1,iElt),zMon(1,iElt),
     &	        CurIndRef,IndRef(iElt),LNewIndRef,NewIndRef, 
     &	        CurExtinc,Extinc(iElt),WaveBU,
     &	        ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	        ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	        RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),RayE(1,iRay),
     &          dxidxim1(1,1,iElt),dxidui(1,1,iElt),SrfType(iElt),
     &	        LensArrayType(iElt),LensArrayWidth(iElt),
     &          IsVarIndRefElt(iElt),TotalLensletsArray(iElt),iElt,
     &	        XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),IERROR,
     &	        IWORK(1,iIntSrf),DWORK(1,iIntSrf),NDP(iIntSrf),mDP)

	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      PrevIndRef=CurIndRef; PrevExtinc=CurExtinc

	      If (LNewIndRef) Then
	        CurIndRef=NewIndRef
	        VarCurIndRef(iRay)=NewIndRef
	      Else
	        CurIndRef=IndRef(iElt)  
	        VarCurIndRef(iRay)=IndRef(iElt)
	      End If
              CurExtinc=Extinc(iElt)
	      ifRayElt(iElt)=.TRUE.

C  End of element section
	    END IF

#if 0
	    if (iStartRay==1738) then
	      print*,'** iStartRay=1738: iElt =',iElt,
     &         ', ifRayElt(iElt) =',ifRayElt(iElt),
     &         ', EltID(iElt) =',EltID(iElt),', LRayOK(iRay)=',LRayOK(iRay),
     &         ', ifDraw =',ifDraw
	    end if

	    if (iElt==25) then
              print*,'** Right after End of element section,LRayOK(iRay)=',
     &	 	LRayOK(iRay), ',  iRay=', iRay
	      print*,'** ifRayElt(iElt) =',ifRayElt(iElt)
	      print*,'EltID(iElt) =',EltID(iElt)
	      print*,'beamOpt%ifOptBeamPos =',beamOpt%ifOptBeamPos
	      print*,'beamOpt%beamPosElt =',beamOpt%beamPosElt 
	    end if
#endif

	    IF (ifRayElt(iElt).AND.
     &          (EltID(iElt).NE.9).AND.(EltID(iElt).NE.12)) THEN
	      ! not for non-sequential reflectors and refractors 
	      If ((beamOpt%ifOptBeamDir).AND.(iElt==beamOpt%beamDirElt)
     &            .AND.(iRay==1)) Then 
                ! Get chief ray direction at beamDirElt, added by jzlou
                ChfRayDirAtElt(1:3)=PrevDir(1:3)
                !print*,'Ctrace 2*: ChfRayDirAtElt(1:3)=',ChfRayDirAtElt(1:3)
              End If
	      If ((beamOpt%ifOptBeamPos).AND.(iElt==beamOpt%beamPosElt)
     &            .AND.(iRay==1)) Then
                ! Get chief ray position at beamPosElt, added by jzlou
                ChfRayPosAtElt(1:3)=RayPos(1:3,iRay)	
              End If
	      If ((beamOpt%ifOptBeamRefRayDir) .and. 
     &            (iElt==beamOpt%beamRefRayElt)) Then
                CurRefRayDir(1:3)=PrevDir(1:3)
!	        print*,' **** iRay,CurRefRayDir(1:3) =',iRay,
!    &                 CurRefRayDir(1:3)
              End If
	      If ((iRay==1) .AND. ChfRayHist%keepHist) Then
	        ! Keep chief ray history
	        ChfRayHist%tot_srfs=ChfRayHist%tot_srfs+1
                ChfRayHist%SrfId(ChfRayHist%tot_srfs)=iElt
                ChfRayHist%RayPosAtSrf(1:3,ChfRayHist%tot_srfs)
     &            =RayPos(1:3,iRay)
                If (.false. .and. iElt==7) Then
	         print*,'***** MACOS at elt 7: RayPos =',RayPos(1:3,iRay)  
                End If
	        ChfRayHist%RayDirAtSrf(1:3,ChfRayHist%tot_srfs)
     &            =PrevDir(1:3)
	      End If
	    END IF

	    ! Note that LRayOK(iRay) is false only when iRay missed the surface.
            ! Obscurations and Apertures don't affect LRayOK(iRay), but they 
            ! affect L1(iRay) (= LRayPass) -jzlou   

	    IF (.NOT.LRayOK(iRay)) THEN
	      nBadRays=nBadRays+1
	      LEltOK(iElt)=.FALSE.
	      LastElt=PrevElt(iElt)
	      IF (IRAY.GT.1) THEN
	        RayMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	        SegMap(RayIndex(1,iRay),RayIndex(2,iRay))=0
	      END IF
	      IF (ifRayPrint) THEN
!	        print*,'** After End of element section ...'
	        WRITE(*,122)iRay,LastElt
	      END IF
	      GO TO 14
	    ELSE IF (ifRayPrint.AND.ifRayElt(iElt).AND.
     &	      (EltID(iElt).NE.9).AND.(EltID(iElt).NE.12)) THEN
	      if (.not.LBWKnoPrint) then
	        WRITE(*,115) iRay, PrevElt(iElt),
     &            EltName(PrevElt(iElt)),iElt,EltName(iElt),
     &            (PrevPos(i),i=1,3),
     &            (RayPos(i,iRay),i=1,3),
     &	          (PrevDir(i),i=1,3),RayL(iRay),
     &	          CumRayL(iRay)
	        WRITE(*,118)PrevIndRef,PrevExtinc,CurIndRef,CurExtinc
	      end if

	      IF (ifPol) THEN
	        Intxyz=RayE(1,iRay)*DCONJG(RayE(1,iRay))+
     &	        RayE(2,iRay)*DCONJG(RayE(2,iRay))+
     &	        RayE(3,iRay)*DCONJG(RayE(3,iRay))
	        WRITE(*,119)RayE(1,iRay),
     &	        RayE(2,iRay),RayE(3,iRay),Intxyz
	      END IF
	      IF (.NOT.LRayTrans) THEN
	        WRITE(*,123)iRay,iElt
	        LRayTrans=.TRUE.
	      END IF
	    ELSE IF (ifDraw.AND.ifRayElt(iElt).AND.
     &	      (EltID(iElt).NE.9).AND.(EltID(iElt).NE.12)) THEN
	      nDrawElt=nDrawElt+1
	      DrawElt(nDrawElt)=iElt
	      DrawRay(1,nDrawElt)=DDOTC(xDraw,RayPos(1,iRay))
	      DrawRay(2,nDrawElt)=DDOTC(yDraw,RayPos(1,iRay))
	    END IF

	    ! Save incident chief ray state:
            ! direction, intersection point, surf nomal, -jzlou
            If (iRay==1) Then
              CRIncidDir(1:3,iElt)=PrevDir(1:3)
              CROutgoDir(1:3,iElt)=RayDir(1:3,iRay)
              CRIncidPos(1:3,iElt)=RayPos(1:3,iRay)
              CRNhat(1:3,iElt)=-NhatP(1:3) ! CRNhat is outward from surface
            End If

	    ! Save chief ray path length at this element
	    If (iRay==1) CumCRL(iElt)=CumRayL(iRay)
	    If (LSavRayL .and. LRayOK(iRay) .and. ifRayElt(iElt)) Then 
              CumRayLsav(iElt)=CumRayL(iRay)
!	     print*,'******svL iElt =', iElt,': LRayOK(iRay),LRayPass(iRay)=',
!    &         LRayOK(iRay),LRayPass(iRay)  
	    End If

	    If (ifSRT) Then
              ! save ray state at this element
              thisRay%InciPos(1:3,iElt)=RayPos(1:3,crnt_traced_ray)
              thisRay%InciDir(1:3,iElt)=RayDir(1:3,crnt_traced_ray)
              thisRay%InciSrfN(1:3,iElt)=Nsrfm(1:3)
              thisRay%PrevIndRef(iElt)=PrevIndRef
              thisRay%CurIndRef(iElt)=CurIndRef
            End If

            ! ******************************
            ! Adding saving RayPosHist(1:3,jRay,iElt) here? -jzlou
            ! ******************************
            If (SaveRayPosHist_Flg .OR. SaveVis3d_dat_Flg) Then
              LRayOKHist(crnt_traced_ray,iElt) = .false.
              If (LRayOK(iRay) .and. ifRayElt(iElt)) Then
                LRayOKHist(crnt_traced_ray,iElt) = .true.
                RayPosHist(1:3,crnt_traced_ray,iElt) =
     &                              RayPos(1:3,crnt_traced_ray)
              End If
            End If

	    iElt=iElt+1

C	  End of iElt loop
 10     CONTINUE


C  5. Compute composite beam train partials:

	  IF (ifBuild) THEN
	    nEm1=iEndElt-1
	    nEm2=nEm1-1
	    NextLast=PrevElt(iEndElt)
	    CALL REQUATE(dxndxi(1,1,NextLast),dxidxim1(1,1,iEndElt),49)
	    CALL REQUATE(dxndui(1,1,iEndElt),dxidui(1,1,iEndElt),42)
	    IF (iEndElt.GT.1) THEN
	      jElt=PrevElt(iEndElt)

	      DO 11 iElt=nEm2,0,-1
!	        print*,' Begin Loop11: jElt=',jElt, 
!    &		  ', PrevElt(jElt)=', PrevElt(jElt)
	        jElt=PrevElt(jElt)
	        IF (ifRayElt(jElt).AND.(jElt.NE.NextLast).AND.
     &	          (EltID(jElt).NE.3).AND.(EltID(jElt).NE.7)) THEN
c	        IF (ifRayElt(jElt).AND.(jElt.NE.NextLast)) THEN
		  CALL RMPRDC(dxndxi(1,1,jElt),
     &            	      dxndxi(1,1,NextElt(jElt)),
     &            	      dxidxim1(1,1,NextElt(jElt)),7,7,7)
	        END IF
 11	      CONTINUE
	      DO 12 iElt=1,nEm1
	        IF (ifRayElt(iElt).AND.(EltID(iElt).NE.3)
     &	            .AND.(EltID(iElt).NE.7)) THEN
c	        IF (ifRayElt(iElt)) THEN
		  CALL RMPRDC(dxndui(1,1,iElt),dxndxi(1,1,iElt),
     &            dxidui(1,1,iElt),7,7,6)
	        ELSE
		  CALL RZERO(dxndui(1,1,iElt),42)
	        END IF
 12	      CONTINUE
	    END IF

C  6. Load Cmatrix:

	    CALL REQUATE(Cmatrix(1,1,iRay),dxndxi(1,1,0),49)
	    DO 13 iElt=1,iEndElt
	      i=8+(iElt-1)*6
	      CALL REQUATE(Cmatrix(1,i,iRay),dxndui(1,1,iElt),42)
	      IF (ifRayElt(iElt)) THEN
	        LRayToElt(iElt,iRay)=.TRUE.
	      ELSE
	        LRayToElt(iElt,iRay)=.FALSE.
	      END IF
 13	    CONTINUE
	  END IF  ! ifBuild

C  7. Compute interpolated surface partials and load CmDM:
C	  IF (ifBuild.AND.ifDMBuild) THEN
C	  END IF

#if 0
	  If (LFastRayTrace) Then
	    if (LRayOK(iRay).AND.LRayPass(iRay)) then
               FastRayTrace(iRay)=.TRUE.
	    else
	      FastRayTrace(iRay)=.FALSE.
	    end if
	  End If
#endif

C  9. End of iRay loop

 14     CONTINUE

C 10. Close up

#if 1
	If (IsVarIndRefElt(iEndElt)) Then
          ! iEndElt has variable IndRef
          VarStrtIndRef(1:iEndRay)=VarCurIndRef(1:iEndRay)
        Else
  	  StartIndRef=CurIndRef
        End If
#endif
	
	StartIndRef=CurIndRef
	StartExtinc=CurExtinc
c	DO 18 iElt=iSecondElt,iEndElt
c	  IF (.NOT.ifIntSrfOK(iElt)) THEN
c	    SrfType(iElt)=5
c	    IF (.NOT.ifDraw) THEN
c	      WRITE(*,124)iElt
c	    END IF
c	  END IF
c 18	CONTINUE
	ifRayPrint=.FALSE.

	RETURN
	END SUBROUTINE CTRACE

C***********************************************************************

	SUBROUTINE CRSOPTIMIZE(iElt)
	use param_mod
        use elt_mod

#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
	IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: ANS
	INTEGER i,iElt,iEm1,nCalls
	REAL*8 Vpt(3),Rpt(3),psi(3),fRSmin,RMSmin,
     &	       TOL,df,fA,fB,fC,RA,RB,RC,zp,L

 108	FORMAT(1P,' Reference surface optimization results:'/
     &	' TTL OPD=',D17.9,'  nCalls=',i4/
     &	'   Old f=',D17.9,'   New f=',D17.9/
     &	'   Old z=',D17.9,'   New z=',D17.9/
     &	' Old psi=',D17.9,' New psi=',D17.9,2(/9x,D17.9,9x,D17.9)/
     &	' Old Vpt=',D17.9,' New Vpt=',D17.9,2(/9x,D17.9,9x,D17.9))

C  2. Initialize various parameters

	nCalls=0
	iEm1=iElt-1
	TOL=CWB*Wavelen/1d2
c	CALL DACCEPT(TOL,TOL,1,'Enter tolerance: ')
	df=fElt(iElt)/1d2
c	CALL DACCEPT(df,df,1,'Enter test increment: ')
	fA=fElt(iElt)
	fB=fA+df

C  3. Compute minimum RMSOPD

c********debug
c	CALL RSTRACE(fA,iElt,RMSOPD,nCalls,Vpt,Rpt,psi,L)
c	CALL RSTRACE(fB,iElt,RMSOPD,nCalls,Vpt,Rpt,psi,L)
c	PAUSE 'Hit return to optimize'
c	write(*,*)('Calling MNZPB')
c	write(*,*)('Calling ZPSolve')
c********debug

	CALL MNZPB(fA,fB,fC,RA,RB,RC,iElt,nCalls,Vpt,Rpt,psi,L)

	RMSmin= 
     &    ZPSolve(fA,fB,fC,iElt,TOL,fRSmin,nCalls,Vpt,Rpt,psi,L)

	IF ((EltID(iElt).EQ.8).AND.(EltID(iEm1).EQ.8)) THEN
	  zp=-SIGN(1d0,DDOTC(RayDir(1,1),psi(1)))*fRSmin
	ELSE
	  zp=SIGN(1d0,DDOTC(RayDir(1,1),psi(1)))*fRSmin
	END IF

	IF (ABS(fElt(iElt)).GE.1d21) THEN
	  fRSmin=fElt(iElt)
	  zp=zElt(iElt)
	END IF

C  4. Print the solution

	WRITE(*,108)RMSmin,nCalls,fElt(iElt),fRSmin,zElt(iElt),zp,
     &	 psiElt(1,iElt),psi(1),psiElt(2,iElt),psi(2),psiElt(3,iElt),
     &	 psi(3),VptElt(1,iElt),Vpt(1),VptElt(2,iElt),Vpt(2),
     &	 VptElt(3,iElt),Vpt(3)

	CALL CACCEPT(ANS,'YES',
     &  'Accept the new element data? ')
	IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	  fElt(iElt)=fRSmin
	  KrElt(iElt)=-(1d0+eElt(iElt))*fElt(iElt)
	  zElt(iElt)=zp
	  DO 1 i=1,3
	    psiElt(i,iElt)=psi(i)
	    VptElt(i,iElt)=Vpt(i)
	    RptElt(i,iElt)=Rpt(i)
  1	  CONTINUE
	END IF

	RETURN
	END SUBROUTINE CRSOPTIMIZE

C***********************************************************************

      SUBROUTINE MNZPB(AX,BX,CX,FA,FB,FC,iElt,nCalls,Vpt,Rpt,psi,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20, MaxCalls=200)
      REAL*8 Vpt(3),Rpt(3),psi(3),L
c     FA=FUNC(AX)
c     FB=FUNC(BX)
      CALL RSTRACE(AX,iElt,FA,nCalls,Vpt,Rpt,psi,L)
      CALL RSTRACE(BX,iElt,FB,nCalls,Vpt,Rpt,psi,L)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
c     FC=FUNC(CX)
      CALL RSTRACE(CX,iElt,FC,nCalls,Vpt,Rpt,psi,L)
1     IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
c         FU=FUNC(U)
          CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BX)
c         FU=FUNC(U)
          CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
c         FU=FUNC(U)
          CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
c           FU=FUNC(U)
            CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
c         FU=FUNC(U)
          CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
        ELSE
          U=CX+GOLD*(CX-BX)
c         FU=FUNC(U)
          CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        IF (nCalls.LT.MaxCalls) GO TO 1
      ENDIF
      RETURN
      END SUBROUTINE MNZPB

C***********************************************************************

      FUNCTION ZPSolve(AX,BX,CX,iElt,TOL,XMIN,nCalls,Vpt,Rpt,psi,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      REAL*8 Vpt(3),Rpt(3),psi(3),L
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
c     FX=F(X)
      CALL RSTRACE(X,iElt,FX,nCalls,Vpt,Rpt,psi,L)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. 
     *        P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
C       FU=F(U)
        CALL RSTRACE(U,iElt,FU,nCalls,Vpt,Rpt,psi,L)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      PAUSE 'ZPSolve exceeded maximum iterations.'
3     XMIN=X
      ZPSolve=FX
      RETURN
      END FUNCTION ZPSolve

C***********************************************************************

	SUBROUTINE RSTRACE(fRS,iElt,RMSOPD,nCalls,Vpt,Rpt,psi,L)
	use param_mod
        use elt_mod
	use elemsub_mod

	IMPLICIT NONE

	LOGICAL :: noBuild,noPol,noObs,LNsrf,ifTransmit,ifReturn,
     &             ifRetSrf,LRT
	INTEGER :: i,iElt,iRay,iEm1,nCalls,jRay,jElt
	SREAL :: dxdx(7,7),dxdu(7,6)
	REAL*8 :: Vpt(3),Rpt(3),psi(3),pr(3),rhat(3),L,fRS,RMSOPD,CumL,
     &	          L1,signpsi,SIGN,KrRS,S1,SumL,SumL2,OPDScal,Dn,DAvgL
	COMPLEX*16 :: C1(3)

C  1. Initialize various parameters

	iEm1=iElt-1
	noBuild=.FALSE.
	noPol=.FALSE.
	noObs=.FALSE.
	LNsrf=.TRUE.
	IF (EltID(iElt).EQ.8) THEN
	  ifRetSrf=.TRUE.
	  ifTransmit=.FALSE.
c	  IF (EltID(iEm1).EQ.8) THEN
c	    ifReturn=.TRUE.
c	  ELSE
c	    ifReturn=.FALSE.
c	  END IF
	ELSE
	  ifRetSrf=.FALSE.
	END IF
	ifReturn=.FALSE.
	DO jElt=1,iEm1
	  IF (EltID(jElt).EQ.8) ifReturn=.NOT.ifReturn
	END DO
	S1=1d0
	nCalls=nCalls+1
	KcElt(iElt)=-eElt(iElt)*eElt(iElt)
	KrRS=-(1d0+eElt(iElt))*fRS

C  2. Set new element location and attitude using chief ray

	IF (.NOT.ifRetSrf) THEN
	  CALL RefSrf(noBuild,noPol,LNsrf,noObs,
     &	  LRayOK(1),LRT,KcElt(iElt),KrRS,
     &	  psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	  RayPos(1,1),RayDir(1,1),C1,S1,S1,S1,S1,S1,
     &	  ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	  ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	  pr,rhat,L,C1,dxdx,dxdu)

	  signpsi=SIGN(1D0,DDOTC(RayDir(1,1),psiElt(1,iElt)))
	ELSE
	  CALL ReturnSrf(noBuild,noPol,ifReturn,ifTransmit,
     &	  LNsrf,noObs,LRayOK(1),LRT,
     &	  KcElt(iElt),KrRS,
     &	  psiElt(1,iElt),VptElt(1,iElt),RptElt(1,iElt),
     &	  RayPos(1,1),RayDir(1,1),C1,S1,S1,S1,S1,S1,
     &	  ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	  ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	  pr,rhat,L,C1,dxdx,dxdu)
	  IF (ifReturn) THEN
	    L=-L
	    signpsi=-SIGN(1D0,DDOTC(RayDir(1,1),psiElt(1,iElt)))
	  END IF
	END IF

	DO 1 i=1,3
	  psi(i)=signpsi*rhat(i)
	  Vpt(i)=pr(i)
	  Rpt(i)=pr(i)
  1	CONTINUE

C  3. Begin iRay loop to trace each ray from previous surface

	SumL=0d0
	SumL2=0d0
	RMSOPD=0d0
	jRay=0
	L1=CumRayL(1)+IndRef(iEm1)*L
	DO 14 iRay=1,nRay
	  IF (LRayOK(iRay).AND.LRayPass(iRay)) THEN
	    jRay=jRay+1
	    IF (.NOT.ifRetSrf) THEN
	      CALL RefSrf(noBuild,noPol,LNsrf,noObs,
     &	      LRayOK(iRay),LRT,KcElt(iElt),KrRS,
     &	      psi,Vpt,Rpt,
     &	      RayPos(1,iRay),RayDir(1,iRay),C1,S1,S1,S1,S1,S1,
     &	      ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	      ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	      pr,rhat,L,C1,dxdx,dxdu)
	    ELSE
	      CALL ReturnSrf(noBuild,noPol,ifReturn,ifTransmit,
     &	      LNsrf,noObs,LRayOK(iRay),LRT,
     &	      KcElt(iElt),KrRS,
     &	      psi,Vpt,Rpt,
     &	      RayPos(1,iRay),RayDir(1,iRay),C1,S1,S1,S1,S1,S1,
     &	      ApType(iElt),ApVec(1,iElt),nObs(iElt),ObsType(1,iElt),
     &	      ObsVec(1,1,iElt),xObs(1,iElt),yObs(1,iElt),zObs(1,iElt),
     &	      pr,rhat,L,C1,dxdx,dxdu)
	      IF (ifReturn) L=-L
	    END IF
	    CumL=CumRayL(iRay)+IndRef(iEm1)*L
	    IF (iRay.EQ.1) THEN
	      L1=CumL
	    ELSE
	      OPDScal=CumL-L1
	      SumL=SumL+OPDScal
	      SumL2=SumL2+OPDScal*OPDScal
	      RMSOPD=RMSOPD+(CumL-L1)**2
	    END IF
	  END IF
 14     CONTINUE

C  4. Compute RMS OPD

	RMSOPD=SQRT(RMSOPD/FLOAT(jRay-1))

c	Dn=jRay-1
c	DAvgL=SumL/Dn
c	RMSOPD=SQRT(DABS((SumL2-2d0*DAvgL*SumL+Dn*DAvgL*DAvgL)/Dn))

C****debug
c	WRITE(*,901)nCalls,fRS,RMSOPD
c 901	FORMAT(1P,'nCalls=',i6,' Focal Length fRS=',D17.9,
c     &	' OPD=',D17.9)
C****debug

	RETURN
	END SUBROUTINE RSTRACE

C***********************************************************************

	SUBROUTINE RotChfRay(th,crDir,crPos,StopPos)
	IMPLICIT NONE

	REAL*8 th(3),Q(3,3),dQ(3,3),D1(9),D2(9),crDir(3),crPos(3),
     &	       StopPos(3)

C Perturb chief ray

	CALL Qform(Q,dQ,th)           ! 'th' is Eular vector of rotation
	CALL DMPROD(D1,Q,crDir,3,3,1) ! ChfRay rotated around 'th' by amt |th|
	CALL DEQUATE(crDir,D1,3)      ! crDir is rotated ChfRay

	CALL DSUB(D1,crPos,StopPos,3)
	CALL DMPROD(D2,Q,D1,3,3,1)
	CALL DADD(crPos,StopPos,D2,3) ! crPos is rotated ChfPos
 
	RETURN
	END SUBROUTINE RotChfRay

C***********************************************************************

	SUBROUTINE FindCrossPt(p1,r1,p2,r2,v,CrossDist,CrossPt)

	IMPLICIT NONE

	REAL*8 p1(3),p2(3),r1(3),r2(3),r1r2i,Cp1(3),Cp2(3),
     &	       rho(3),CrossPt(3),del(3),CrossDist,CL,v

c   r1r2i=1/(r1'*r2);
	r1r2i=1d0/DDOTC(r1,r2)
c   del=p2-p1;
	CALL DSUB(del,p2,p1,3)
c   CL=(r1'*del)*r1r2i;
	CL=DDOTC(r1,del)*r1r2i
c   Cp1=CL*r2-del
	Cp1(1)=CL*r2(1)-del(1)
	Cp1(2)=CL*r2(2)-del(2)
	Cp1(3)=CL*r2(3)-del(3)
c   Cp2=-r1+r2*r1r2i
	Cp2(1)=-r1(1)+r2(1)*r1r2i
	Cp2(2)=-r1(2)+r2(2)*r1r2i
	Cp2(3)=-r1(3)+r2(3)*r1r2i
c   v=Cp1'*Cp2/(Cp2'*Cp2)
	v=DDOTC(Cp1,Cp2)/DDOTC(Cp2,Cp2)
c   rho=Cp1-Cp2*v
	rho(1)=Cp1(1)-Cp2(1)*v
	rho(2)=Cp1(2)-Cp2(2)*v
	rho(3)=Cp1(3)-Cp2(3)*v
c   rhosq=rho'*rho
	CrossDist=SQRT(DDOTC(rho,rho))
c   p1v=p1+v*r1
	CrossPt(1)=p1(1)+v*r1(1)
	CrossPt(2)=p1(2)+v*r1(2)
	CrossPt(3)=p1(3)+v*r1(3)
 
	RETURN
	END SUBROUTINE FindCrossPt

C***********************************************************************

	SUBROUTINE WARN(LEltOK,nElt,nBadRays)
	IMPLICIT NONE
	INTEGER nElt,iElt,nBadRays
	LOGICAL LEltOK(nElt)

  124	FORMAT(' One or more rays become undefined at element',i4)
  125	FORMAT(' A total of ',i8,' rays were lost')
	DO 1 iElt=1,nElt
	  IF (.NOT.LEltOK(iElt)) THEN
	    WRITE(*,124)iElt
	  END IF
   1	CONTINUE
	WRITE(*,125)nBadRays
	RETURN
	END SUBROUTINE WARN

C***********************************************************************

! glassProp argument is the six glass coefs of that glass type;
! Given wavelength 'lamda', index of refraction can then be computed -jzlou
        FUNCTION getIndex(lamda,glassProp,eltGlassType)
	use param_mod

C       Compute index of refraction at wavelengths in microns

	IMPLICIT NONE

        REAL*8 :: A1,A2,A3,glassProp(6),lamda,
     &            p2,p4,p6,p8,getIndex
	INTEGER :: eltGlassType

        ! There are only *two types* of glasses, Sellmeier type
        ! or Schott type, which have different formulas for computing
        ! indx given lamda and glass coefs (or glassProp) -jzlou 
	If (eltGlassType==SellmeierGlass) Then
C The Sellmeier formula
C  Index*Index-1=B1*lamda*lamda/(lamda*lamda-C1)
C               +B2*lamda*lamda/(lamda*lamda-C2)
C               +B3*lamda*lamda/(lamda*lamda-C3)
          p2=lamda*lamda
          A1=glassProp(1)*p2/(p2-glassProp(4))
          A2=glassProp(2)*p2/(p2-glassProp(5))
          A3=glassProp(3)*p2/(p2-glassProp(6))
          getIndex=sqrt(1+A1+A2+A3)
	Else If (eltGlassType==SchottGlass) then
	  p2=lamda*lamda;p4=p2*p2;p6=p2*p2*p2;p8=p4*p4
	  getIndex=glassProp(1)+glassProp(2)*p2+
     &             glassProp(3)/p2+glassProp(4)/p4+
     &             glassProp(5)/p6+glassProp(6)/p8
	  getIndex=sqrt(getIndex)
	Else
	  print*,'getIndex(): eltGlassType=', eltGlassType
	  print*,'getIndex(): unknown glass type, quitting ...'
	  stop
	End If

        RETURN
        END FUNCTION getIndex

	End Module tracesub_mod
