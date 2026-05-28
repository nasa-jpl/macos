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
     &    ifLNegOK,ifBuild,ifTrace,ifRegrid(mElt),
     &    ifIntSrfSet(mElt),ifDMBuild,ifDraw,L1(mRay)

	INTEGER npts,i,iCurWFElt,iCurRayElt,
     &	  nBadRays,iStartElt,iEndElt,nPassRays,iObsOpt,
     &    iStartRay,iEndRay,nDrawElt,iElt,iRay,mpt2,
     &    DrawEltVec(mDrawElt,mDrawRay)

	SREAL StartTime,DTime,EndTime,
     &	  Cmatrix(7,mCm,bRay),
     &	  DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,Dn,
     &    WFE,WFEPV,MinL,MaxL,DAvgL,SumL,SumL2,SumCumL,
     &    OPDScal,OPDMat(mpts,mpts),
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
	  IF (iStartElt.EQ.0)
     &      CALL SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
     &                            ifPropagate,ifGrid)

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
	  ifTrace        = .TRUE.
          ifCtraceCalled = .TRUE.
	  ifOPDCalc      = .TRUE.
	  ifBuild        = .FALSE.
	  ifDMBuild      = .FALSE.
	  ifPropagate    = .FALSE.
	  ifRayPrint     = .FALSE.
	  CALL RTIME(EndTime)
	  DTime=EndTime-StartTime
	  if ((.not.isDOP) .and. (LShowDTime)) WRITE(*,116) DTime
        ELSE
          ifCtraceCalled=.FALSE.
	END IF

	nGridPts = npts+1
	SumL     = 0d0
	SumL2    = 0d0
        SumCumL  = 0d0
	MaxL     =-1d38
	MinL     = 1d38
	mpt2     = mpts*mpts
	CALL DZERO(OPDMat,mpt2)

	! Following code calculates RMS WFE, in which chief ray
	! is not really needed (would be cancelled in RMS formula) -jzlou
cstartscottfix
	IF (LRayOK(1) .and. (.not. OPDRefRayLen_FLG)) THEN
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
     &	  zp,psip(3),CrossPt(3),cr1dir(3),cr1pos(3),cr2dir(3),cr2pos(3),
     &	  cr2indir(3),cr2inpos(3),CrossDist

C         ***Find exit pupil (with option to minimize tilt)***

	CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate)

	iEm1      = iElt-1
	iStartRay = 1
	iEndRay   = 1
	iStartElt = 0
	iEndElt   = iEm1
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
          !   print*,'***** FEX: ifCentroid used!\n'

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
	    IF (iStartElt.EQ.0)
     &        CALL SetSourceRayGrid(nTestPts,Dicr,Djcr,dxSource,ifPol,
     &                              ifPropagate, ifGrid)

            iStartRay = 1
            iEndRay   = nRay

	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &        ifLNegOK,Cmatrix,
     &        iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &        ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &        StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

	    iCurRayElt  = iEndElt
	    ifTrace     = .FALSE.
	    ifBuild     = .FALSE.
	    ifPropagate = .FALSE.
	  END IF

	  iSpot    = 0
	  nGridPts = npts+1
	  mpt2     = mpts*mpts
	  DO iRay=1,nRay
	    IF (LRayOK(iRay).AND.((iRay.EQ.1).OR.(iObsOpt.EQ.0).OR.
     &          (LRayPass(iRay).AND.(iObsOpt.EQ.1)))) THEN
	      IF (iRay.EQ.1) THEN
	        CALL DEQUATE(CRSpot,RayPos(1,1),3)
	        CALL DZERO(DV1,3)
	      END IF
	      IF ((iRay.NE.1).OR.
     &          ((iRay.EQ.1).AND.(nGridpts.EQ.1))) THEN
                iSpot=iSpot+1
                CALL DADD(DV1,DV1,RayPos(1,iRay),3)
	      END IF
	    END IF
	  END DO
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
	  !   print*,'***** FEX: Centroid NOT used!'
	END IF
	END SUBROUTINE FEX


C***********************************************************************
C       Set eXit Pupil (SXP)
C
C       Clone of FEX with one geometry fix per Dave's spec (2026-05-19):
C       the returned zp -- which the caller uses as the EP's focal
C       length / radius of curvature -- is the chief-ray distance from
C       the EP (CrossPt) to the FP (assumed iElt+1 in a typical imager
C       Rx: ... -> EP -> FP), instead of FEX's legacy distance from
C       iEm1 to the EP.  Makes the EP radius sensitive to FP
C       perturbations through Tz_FP (focus) -- the original FEX is
C       insensitive because iEm1 is upstream and doesn't move when FP
C       moves.  Lateral FP perturbations (Tx/Ty) and FP rotations are
C       NOT captured here; they need the EP vertex / orientation to
C       follow, not just the radius (use the dw/dx FocalPlaneChannel
C       "track" mode for full coverage).
C
C       FEX is left untouched to avoid breaking existing prescriptions,
C       GMI workflows, and journals that rely on its current behavior.
C***********************************************************************

	SUBROUTINE SXP(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
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
     &	  iElt,iEm1,iFP_loc

	SREAL Cmatrix(7,mCm,bRay),DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &	  xDraw(3),yDraw(3),CentroidSpot(3),dCRCentroid(3),CRSpot(3),
     &	  DV1(9),S1,th(3),
     &	  zp,psip(3),CrossPt(3),cr1dir(3),cr1pos(3),cr2dir(3),cr2pos(3),
     &	  cr2indir(3),cr2inpos(3),CrossDist,
     &	  dvec(3),denom,t_FP,dot_dvec,zp_iEm1

C       Standard FEX chief-ray pair trace (source -> iEm1) and find
C       the two-ray crossing point -- that's the EP location.

	CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate)

	iEm1      = iElt-1
	iStartRay = 1
	iEndRay   = 1
	iStartElt = 0
	iEndElt   = iEm1
	CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &		ifLNegOK,Cmatrix,
     &		iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &		ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &		StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

	CALL DEQUATE(cr1dir,RayDir(1,1),3)
	CALL DEQUATE(cr1pos,RayPos(1,1),3)
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

	CALL FindCrossPt(cr1pos,cr1dir,cr2pos,cr2dir,zp_iEm1,
     &		         CrossDist,CrossPt)

C       Replace the FEX zp (iEm1-to-EP distance) with the chief-ray
C       distance from the EP (CrossPt) to the FP (iElt+1).
C       Post-EP chief ray direction = -cr1dir for a normally-incident
C       chief ray reflecting off the EP Return surface at its vertex.
C       r(t) = CrossPt - t * cr1dir; intersect with FP plane.

	iFP_loc = iElt + 1
	IF (iFP_loc .LE. nElt) THEN
	  dvec(1) = VptElt(1,iFP_loc) - CrossPt(1)
	  dvec(2) = VptElt(2,iFP_loc) - CrossPt(2)
	  dvec(3) = VptElt(3,iFP_loc) - CrossPt(3)
	  denom = cr1dir(1)*psiElt(1,iFP_loc) +
     &            cr1dir(2)*psiElt(2,iFP_loc) +
     &            cr1dir(3)*psiElt(3,iFP_loc)
	  IF (DABS(denom) .GT. 1d-15) THEN
	    dot_dvec = dvec(1)*psiElt(1,iFP_loc) +
     &                 dvec(2)*psiElt(2,iFP_loc) +
     &                 dvec(3)*psiElt(3,iFP_loc)
	    t_FP = -dot_dvec / denom
	    zp = t_FP
	  ELSE
	    zp = zp_iEm1
	  END IF
	ELSE
	  zp = zp_iEm1
	END IF

	print*,'***** SXP: CrossPt =',CrossPt
        print*,'***** SXP: zp_iEm1 =',zp_iEm1,
     &         '  zp_FP =',zp

C       Centroid / chief-ray option -- unchanged from FEX.

	IF (ifCentroid .AND. Rx_FEXCentrFlg) THEN
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
	    IF (iStartElt.EQ.0)
     &        CALL SetSourceRayGrid(nTestPts,Dicr,Djcr,dxSource,ifPol,
     &                              ifPropagate, ifGrid)

            iStartRay = 1
            iEndRay   = nRay

	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &        ifLNegOK,Cmatrix,
     &        iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &        ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &        StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

	    iCurRayElt  = iEndElt
	    ifTrace     = .FALSE.
	    ifBuild     = .FALSE.
	    ifPropagate = .FALSE.
	  END IF

	  iSpot    = 0
	  nGridPts = npts+1
	  mpt2     = mpts*mpts
	  DO iRay=1,nRay
	    IF (LRayOK(iRay).AND.((iRay.EQ.1).OR.(iObsOpt.EQ.0).OR.
     &          (LRayPass(iRay).AND.(iObsOpt.EQ.1)))) THEN
	      IF (iRay.EQ.1) THEN
	        CALL DEQUATE(CRSpot,RayPos(1,1),3)
	        CALL DZERO(DV1,3)
	      END IF
	      IF ((iRay.NE.1).OR.
     &          ((iRay.EQ.1).AND.(nGridpts.EQ.1))) THEN
                iSpot=iSpot+1
                CALL DADD(DV1,DV1,RayPos(1,iRay),3)
	      END IF
	    END IF
	  END DO
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
	END IF
	END SUBROUTINE SXP


C***********************************************************************
        ! Find Detector (Focal) Plane Position and Orientation

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
     &    ifDMBuild,ifDraw

        INTEGER :: dpElt,npts,iCurWFElt,iCurRayElt,
     &             nBadRays,iStartElt,iEndElt,iObsOpt,
     &             DrawEltVec(mDrawElt,mDrawRay)

        SREAL :: Cmatrix(7,mCm,bRay),
     &           DrawRayVec(2,mDrawElt,mDrawRay)

	REAL*8 :: Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &          ChfRayDirSv(3),ChfRayPosSv(3),cr1Pos(3),cr1Dir(3),
     &          cr2Pos(3),cr2Dir(3),zp,th(3),CrsDist,DplnPos(3),
     &          DplnPosSag(3),DplnPosTan(3)

	CALL SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate)

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

	cr1Pos = CRIncidPos(1:3,iEndElt)
        cr1Dir = CROutgoDir(1:3,iEndElt)

        IF (isPointSource()) THEN
          ! Point source
          ! Trace a differential chief ray rotated wrt xGrid
          CALL DSMPROD(th,xGrid,5d-6,3)
          CALL RotChfRay(th,ChfRayDir,ChfRayPos,ChfRayPosSv)
        ELSE
          ! Collimated source
          ! Trace a differential chief ray displaced in -yGrid
          CALL DSMPROD(th,-yGrid,1d-5,3)
          CALL DADD(ChfRayPos,ChfRayPosSv,th,3)
        END IF
	RayDir(1:3,1) = ChfRayDir(1:3)
	RayPos(1:3,1) = ChfRayPos(1:3)

	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)
	cr2Pos      = CRIncidPos(1:3,iEndElt)
        cr2Dir      = CROutgoDir(1:3,iEndElt)
        !cr2Dir      = CRIncidDir(1:3,iEndElt)
	ifTrace     = .FALSE.;
        ifBuild     = .FALSE.;
        ifPropagate = .FALSE.

	! Find crossing point of 2 chief rays
        CALL FindCrossPt(cr1Pos,cr1Dir,cr2Pos,cr2Dir,
     &                   zp,CrsDist,DplnPosSag)

        IF (isPointSource()) THEN
          ! Point source
          ! Trace a differential chief ray rotated wrt xGrid
          CALL DSMPROD(th,yGrid,5d-6,3)
          CALL RotChfRay(th,ChfRayDir,ChfRayPos,ChfRayPosSv)
        ELSE
          ! Collimated source
          ! Trace a differential chief ray displaced in -yGrid
          CALL DSMPROD(th,-xGrid,1d-5,3)
          CALL DADD(ChfRayPos,ChfRayPosSv,th,3)
        END IF
        RayDir(1:3,1) = ChfRayDir(1:3)
        RayPos(1:3,1) = ChfRayPos(1:3)

	CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)

        cr2Pos      = CRIncidPos(1:3,iEndElt)
        cr2Dir      = CROutgoDir(1:3,iEndElt)
	ifTrace     = .FALSE.
        ifBuild     = .FALSE.
        ifPropagate = .FALSE.

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

#if 1


        SUBROUTINE STOP(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,Cmatrix,DrawRayVec,
     &	  DrawEltVec,iLocal,xTarg,sqErr,ifBeamCentered)

          use param_mod
          use src_mod
          use cfiles_mod
          use elt_mod
          use sourcsub_mod

          IMPLICIT NONE

          ! INPUT arguments
          INTEGER,  INTENT(IN)    :: npts              ! Number of points in ray grid
          INTEGER,  INTENT(IN)    :: iObsOpt           ! Observation option
          INTEGER,  INTENT(INOUT) :: iEndElt           ! End element for trace (defined by iLocal)
          INTEGER,  INTENT(IN)    :: iCurWFElt         ! Current wavefront element
          INTEGER,  INTENT(IN)    :: iCurRayElt        ! Current ray element
          INTEGER,  INTENT(IN)    :: iLocal            ! Local element ID (STOP)

          LOGICAL,  INTENT(IN)    :: ifPol             ! Polarization flag
          LOGICAL,  INTENT(INOUT) :: ifGrid            ! Grid flag
          LOGICAL,  INTENT(IN)    :: ifPropagate       ! Propagation flag
          LOGICAL,  INTENT(IN)    :: ifDraw            ! Draw flag
          LOGICAL,  INTENT(IN)    :: ifLNegOK          ! Allow negative L flag
          LOGICAL,  INTENT(IN)    :: ifBuild           ! Build flag
          LOGICAL,  INTENT(IN)    :: ifDMBuild         ! DM build flag
          LOGICAL,  INTENT(IN)    :: ifTrace           ! Trace flag
          LOGICAL,  INTENT(IN)    :: ifRayPrint        ! Ray print flag
          LOGICAL,  INTENT(IN)    :: ifRegrid(mElt)    ! Regrid flags per element
          LOGICAL,  INTENT(IN)    :: ifIntSrfSet(mElt) ! Interface surface flags

          REAL(8),  INTENT(INOUT) :: Dicr              ! Grid increment i
          REAL(8),  INTENT(INOUT) :: Djcr              ! Grid increment j
          REAL(8),  INTENT(INOUT) :: dxSource          ! Source offset
          REAL(8),  INTENT(IN)    :: StartIndRef       ! Starting index of refraction
          REAL(8),  INTENT(IN)    :: StartExtinc       ! Starting extinction
          REAL(8),  INTENT(IN)    :: xTarg(3)          ! Target position in local coords

          SREAL,    INTENT(INOUT) :: Cmatrix(7,mCm,bRay)  ! Coating matrix

          INTEGER,  INTENT(INOUT) :: DrawEltVec(mDrawElt,mDrawRay)   ! Draw element vector (mod by LocalCoord)
          SREAL,    INTENT(INOUT) :: DrawRayVec(2,mDrawElt,mDrawRay) ! Draw ray vector     (mod by LocalCoord)

          ! OUTPUT arguments
          INTEGER,  INTENT(OUT)   :: nBadRays          ! Number of failed rays
          REAL(8),  INTENT(OUT)   :: sqErr             ! Squared error at convergence
          LOGICAL,  INTENT(OUT)   :: ifBeamCentered    ! Success flag

          ! Local Variables
          LOGICAL :: local_coord_ok
          REAL(8) :: xLocal(3), yLocal(3), zLocal(3), xInt(2)
          REAL(8) :: crLocalPos(3)

          ! -------------------------------------------------------------------
          local_coord_ok = .false.

          ! ============================================================
          ! COMPUTE LOCAL COORDINATE SYSTEM AT STOP ELEMENT
          ! Must be done BEFORE point source setup since xLocal/yLocal
          ! are used as perturbation directions in point source mode
          ! ============================================================
          if (nECoord(iLocal)>0) then

            ! use TElt when it is defined
            xLocal = TElt(1:3, 1, iLocal)
            yLocal = TElt(1:3, 2, iLocal)
            zLocal = TElt(1:3, 3, iLocal)
            local_coord_ok = .True.

          else

            ! define CSYS  (this should be defined when reading Rx  -- ToDo)
            call define_local_csys(psiElt(1:3,iLocal), xLocal, yLocal, zLocal, local_coord_ok)

            if (.not.local_coord_ok)
     &        CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &                    ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &                    ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,
     &                    dxSource,nBadRays,StartIndRef,StartExtinc,
     &                    Cmatrix,DrawRayVec,DrawEltVec,
     &                    iLocal,xLocal,yLocal,zLocal,crLocalPos,local_coord_ok)

            if (.not.local_coord_ok) then
              WRITE(*,*)'ERROR: Local coordinates could not be computed'
              return
            end if

            ! update TElt
            TElt(1:3, 1, iLocal) = xLocal
            TElt(1:3, 2, iLocal) = yLocal
            TElt(1:3, 3, iLocal) = zLocal

            TElt(4:6, 4, iLocal) = xLocal
            TElt(4:6, 5, iLocal) = yLocal
            TElt(4:6, 6, iLocal) = zLocal

            LUpdateTElt_FLG(iLocal) = .True.
            nECoord(iLocal)         = 6

          end if


          CALL ChiefRayAiming(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,Cmatrix,DrawRayVec,
     &	  DrawEltVec,iLocal,xTarg,sqErr,ifBeamCentered,xInt)

          if (ifBeamCentered) ifBeamCentered = UpdSrcGrid()

        END SUBROUTINE STOP


#define ChiefRayAiming_Debug 0

        SUBROUTINE ChiefRayAiming(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,Cmatrix,DrawRayVec,
     &	  DrawEltVec,iLocal,xTarg,sqErr,ifBeamCentered,xInt)

          use param_mod
          use src_mod
          use cfiles_mod
          use elt_mod
          use sourcsub_mod
          use Constants, only: EPS, EPS_CBRT, EPS_SQRT, TOL_ZERO, TOL_GEOM2, TOL_GEOM

          IMPLICIT NONE

          ! INPUT arguments
          INTEGER,  INTENT(IN)    :: npts              ! Number of points in ray grid
          INTEGER,  INTENT(IN)    :: iObsOpt           ! Observation option
          INTEGER,  INTENT(INOUT) :: iEndElt           ! End element for trace (defined by iLocal)
          INTEGER,  INTENT(IN)    :: iCurWFElt         ! Current wavefront element
          INTEGER,  INTENT(IN)    :: iCurRayElt        ! Current ray element
          INTEGER,  INTENT(IN)    :: iLocal            ! Local element ID (STOP)

          LOGICAL,  INTENT(IN)    :: ifPol             ! Polarization flag
          LOGICAL,  INTENT(INOUT) :: ifGrid            ! Grid flag
          LOGICAL,  INTENT(IN)    :: ifPropagate       ! Propagation flag
          LOGICAL,  INTENT(IN)    :: ifDraw            ! Draw flag
          LOGICAL,  INTENT(IN)    :: ifLNegOK          ! Allow negative L flag
          LOGICAL,  INTENT(IN)    :: ifBuild           ! Build flag
          LOGICAL,  INTENT(IN)    :: ifDMBuild         ! DM build flag
          LOGICAL,  INTENT(IN)    :: ifTrace           ! Trace flag
          LOGICAL,  INTENT(IN)    :: ifRayPrint        ! Ray print flag
          LOGICAL,  INTENT(IN)    :: ifRegrid(mElt)    ! Regrid flags per element
          LOGICAL,  INTENT(IN)    :: ifIntSrfSet(mElt) ! Interface surface flags

          REAL(8),   INTENT(INOUT) :: Dicr              ! Grid increment i
          REAL(8),   INTENT(INOUT) :: Djcr              ! Grid increment j
          REAL(8),   INTENT(INOUT) :: dxSource          ! Source offset
          REAL(8),   INTENT(IN)    :: StartIndRef       ! Starting index of refraction
          REAL(8),   INTENT(IN)    :: StartExtinc       ! Starting extinction
          REAL(8),   INTENT(IN)    :: xTarg(3)          ! Target position in local coords

          SREAL,    INTENT(INOUT) :: Cmatrix(7,mCm,bRay)  ! Coating matrix

          INTEGER,  INTENT(INOUT) :: DrawEltVec(mDrawElt,mDrawRay)  ! Draw element vector (mod by LocalCoord)
          SREAL,    INTENT(INOUT) :: DrawRayVec(2,mDrawElt,mDrawRay) ! Draw ray vector    (mod by LocalCoord)

          ! OUTPUT arguments
          INTEGER,  INTENT(OUT)   :: nBadRays          ! Number of failed rays
          REAL(8),  INTENT(OUT)   :: sqErr             ! Squared error at convergence
          LOGICAL,  INTENT(OUT)   :: ifBeamCentered    ! Success flag
          REAL(8),  INTENT(OUT)   :: xInt(2)           ! Final target position

          ! LOCAL variables
          INTEGER :: k, iStartElt, iStartRay, iEndRay, nDrawElt

          LOGICAL :: ifDraw_, ifBuild_, ifLocalCoord

          REAL(8) :: DV1(3), DV2(3), S1
          REAL(8) :: thx, thy, dthPt, ddxPt(2)
          REAL(8) :: xLocal(3), yLocal(3), zLocal(3)

          ! State management: separate accepted state from gradient baseline
          REAL(8) :: acceptedChfRayDir(3), acceptedChfRayPos(3)  ! Last accepted state
          REAL(8) :: baselineChfRayDir(3), baselineChfRayPos(3)  ! Baseline for gradient FD
          REAL(8) :: crLocalPos(3), dxPt(2), GPt(2,2), GPtInv(2,2)
          REAL(8) :: del(3)
          REAL(8) :: xDraw(3), yDraw(3)

          ! Point source support variables
          LOGICAL :: ifPointSrc
          REAL(8) :: saveSrcPos(3), fixedChfRayPos(3)
          REAL(8) :: ChfRayDir_pert(3), norm_dir, det_GPt
          REAL(8) :: sqErr_prev, step_norm, step_limit
          REAL(8) :: dthPt_fd  ! Finite difference perturbation (separate from step control)
          REAL(8) :: dthPt_fd_max  ! Maximum safe perturbation
          REAL(8) :: convergence_rate
          REAL(8) :: grad_x_mag, grad_y_mag  ! Gradient magnitudes
          REAL(8) :: grad_threshold          ! Adaptive gradient threshold (LM)
          LOGICAL :: grad_x_ok, grad_y_ok  ! Gradient quality flags
          INTEGER :: n_error_increase  ! Count consecutive error increases
          INTEGER :: n_grad_increase  ! Count gradient increase attempts

          ! Levenberg-Marquardt variables
          REAL(8) :: lambda_lm        ! LM damping parameter (adaptive)
          REAL(8) :: JtJ(2,2)         ! J^T * J (Hessian approximation)
          REAL(8) :: Jt_f(2)          ! J^T * f (gradient)
          REAL(8) :: A_lm(2,2)        ! Augmented system matrix
          REAL(8) :: det_lm           ! Determinant of augmented system
          REAL(8) :: sqErr_trial      ! Trial error for step acceptance
          LOGICAL :: step_accepted   ! Step acceptance flag
          LOGICAL :: ray_success

          ! Condition number estimation and adaptive regularization
          REAL(8) :: cond_est, reg_factor

          ! Stagnation and failure detection
          INTEGER :: n_stagnant      ! Consecutive iterations with minimal progress
          INTEGER :: n_ray_failures  ! Consecutive ray failures in inner loop
          REAL(8) :: sqErr_best       ! Best error achieved so far

          ! Source CSYS (keep original -- stable)
          REAL(8) :: xSrcGrid(3), ySrcGrid(3)

          ! Control flow variables
          LOGICAL :: converged      !, local_coord_ok

          INTEGER, PARAMETER :: MAXITERRATIONS = 100
          INTEGER, PARAMETER :: MAX_BACKTRACK = 20

          ! Tolerances (connected to constants.f90 for consistency)
          real(8), parameter :: errTol             = TOL_GEOM2,      ! ~2.7e-22 mm^2 (geometric tolerance squared)
     &                          dthPt_min          = EPS_SQRT,       ! ~1.49e-8 (minimum for finite differences)
     &                          det_GPt_min        = TOL_ZERO*1d-2,  ! ~6e-22 (more lenient Jacobian threshold)
     &                          alpha_min          = EPS_CBRT,       ! ~6e-6 (minimum backtracking fraction)
     &                          dthPt_adaptive_min = EPS_SQRT        ! ~1.49e-8 (adaptive step threshold)

          ! LM regularization options
          LOGICAL, PARAMETER :: ifUseScaling = .FALSE.  ! Use scaled condition number estimation

          ! ------------------------------------------------
          ! Initialize control flags
          ifBeamCentered = .FALSE.  ! Assume failure until proven otherwise
          converged      = .FALSE.

          ! Ray-Trace Setup parameters
          ifDraw_   = .FALSE.
          ifBuild_  = .FALSE.
          iStartRay = 1
          iEndRay   = 1
          iStartElt = 0
          iEndElt   = iLocal  ! = id of STOP element
          nMinPts   = 3

          ! Tolerances & Tracking
          dthPt              = 1d-3   ! Initial perturbation: 0.001 rad or mm
          n_error_increase   = 0      ! Track consecutive error increases
          n_grad_increase    = 0      ! Track gradient increase attempts
          sqErr_prev         = 1d10  ! Initialize for error monitoring

          ! Adaptive finite difference perturbation strategy
          ! Start small to avoid ray failures, increase if gradients insufficient
          ifPointSrc = isPointSource()
          IF (ifPointSrc) THEN
            dthPt_fd     = 1d-6  ! Start small to avoid ray failures
            dthPt_fd_max = 1d-3  ! Maximum safe angular perturbation
          ELSE
            dthPt_fd     = 1d-6  ! Collimated: position perturbation
            dthPt_fd_max = 1d-3  ! Maximum safe position perturbation
          END IF

          ! Levenberg-Marquardt initialization
          lambda_lm = 1d-3  ! Initial damping (will adapt during iteration)
          n_stagnant = 0    ! No stagnation yet
          sqErr_best = 1d10 ! Initialize best error

#if ChiefRayAiming_Debug
            WRITE(*,*)' '
            WRITE(*,*)' =========================================='
            WRITE(*,*)' USING LEVENBERG-MARQUARDT ALGORITHM'
            WRITE(*,*)' =========================================='
            WRITE(*,*)' Initial damping: lambda =', lambda_lm
#endif ! ChiefRayAiming_Debug

          ! Initialize accepted state with initial configuration
          acceptedChfRayDir(1:3) = ChfRayDir(1:3)
          acceptedChfRayPos(1:3) = ChfRayPos(1:3)

          ! require TElt to be defined
          if (nECoord(iLocal) < 0) return

          ! use TElt
          xLocal = TElt(1:3, 1, iLocal)
          yLocal = TElt(1:3, 2, iLocal)
          zLocal = TElt(1:3, 3, iLocal)

#if ChiefRayAiming_Debug
            !
            ! DISPLAY SETUP BASED ON SOURCE TYPE
            !
            IF (ifPointSrc) THEN
              ! Point source: will modify direction, keep source position fixed
              saveSrcPos = SourcePos()

              WRITE(*,*)' '
              WRITE(*,*)' ========================================'
              WRITE(*,*)' STOP SETTING: Point Source Mode'
              WRITE(*,*)' ========================================'
              WRITE(*,*)' Will modify ChfRayDir (direction)'
              WRITE(*,*)' Fixed SourcePos =', saveSrcPos
              WRITE(*,*)' Initial ChfRayDir =', ChfRayDir
              WRITE(*,*)' STOP surface axes (from optical prescription):'
              WRITE(*,*)' - xLocal =', xLocal
              WRITE(*,*)' - yLocal =', yLocal
            ELSE
              ! Collimated source: will modify position (current behavior)
              WRITE(*,*)' '
              WRITE(*,*)' ========================================'
              WRITE(*,*)' STOP SETTING: Collimated Source Mode'
              WRITE(*,*)' ========================================'
              WRITE(*,*)' Will modify ChfRayPos (position)'
              WRITE(*,*)' Initial ChfRayPos =', ChfRayPos
              WRITE(*,*)' Initial ChfRayDir =', ChfRayDir
            END IF
            WRITE(*,*)' Target (local): xTarg =', xTarg(1:2)
            WRITE(*,*)' Tolerance: errTol =', errTol
            WRITE(*,*)' '

#else ! ChiefRayAiming_Debug

            ! Point source: will modify direction, keep source position fixed
            IF (ifPointSrc) saveSrcPos = SourcePos()

#endif ! ChiefRayAiming_Debug


          ! Trace chief ray from source to STOP element (=iEndElt)
          CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,ifPropagate, ifGrid)
          CALL CTRACE(RayIndex,ifBuild_,ifDMBuild,ifPol,ifDraw_,
     &                    ifLNegOK,Cmatrix,
     &                    iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &                    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &                    StartExtinc,DrawRayVec,DrawEltVec,nDrawElt, xDraw,yDraw)

          crLocalPos = RayPos(:3, 1)  ! chief ray

          ! Source CSYS (retained from initial configuration)
          ! => zGrid, xGrid, yGrid are updated after SetSourceRayGrid call
          xSrcGrid(:) = xGrid(:)
          ySrcGrid(:) = yGrid(:)

          ! Check if initial ray trace succeeded
          IF (.NOT. LRayOK(1)) THEN

#if ChiefRayAiming_Debug
              WRITE(*,*)' '
              WRITE(*,*)' ========================================='
              WRITE(*,*)' ERROR: Initial ray trace FAILED'
              WRITE(*,*)' ========================================='
              WRITE(*,*)' The target position may be unreachable due to:'
              WRITE(*,*)' - Aperture limits or vignetting'
              WRITE(*,*)' - Total internal reflection'
              WRITE(*,*)' - Geometric constraints of the optical system'
              WRITE(*,*)' '
              WRITE(*,*)' Cannot proceed with STOP setting.'
#endif ! ChiefRayAiming_Debug

            RETURN
          END IF

          ! Nominal (unperturbed) incidence pt in local coords
          DV1(1:3) = crLocalPos(1:3) - VptElt(1:3, iLocal)  ! CALL DSUB(DV1,crLocalPos,VptElt(1,iLocal),3)
          xInt(1)  = DDOTC(xLocal,DV1)
          xInt(2)  = DDOTC(yLocal,DV1)

          ! Initialize sqErr_best with actual initial error for proper stagnation detection
          dxPt(1) = xInt(1) - xTarg(1)
          dxPt(2) = xInt(2) - xTarg(2)
          sqErr_best = dxPt(1)*dxPt(1) + dxPt(2)*dxPt(2)

#if ChiefRayAiming_Debug
            !
            ! INITIAL STATE DIAGNOSTICS
            !
            WRITE(*,*)' '
            WRITE(*,*)' === INITIAL STATE ==='
            WRITE(*,*)' RayPos         =', RayPos(:3,1)
            WRITE(*,*)' d(RayPos)      =', DV1(:3)
            WRITE(*,*)' xInt (current) =', xInt(1:2)
            WRITE(*,*)' xTarg (target) =', xTarg(1:2)
            WRITE(*,*)' Initial error  =', xInt(1:2) - xTarg(1:2)
            WRITE(*,*)' Initial sqErr  =', (xInt(1)-xTarg(1))**2 + (xInt(2)-xTarg(2))**2
            WRITE(*,*)' sqErr_best     =', sqErr_best
            WRITE(*,*)' '

#endif ! ChiefRayAiming_Debug

          ! -------------------------------------------------------------------
          ! MAIN ITERATION LOOP
          ! -------------------------------------------------------------------
          ! State management:
          ! - Accepted state  = last successfully accepted configuration
          ! - Baseline state  = state at iteration start for gradient computation
          !
          ! At each iteration:
          !   1. Save current state -> baseline
          !   2. Compute gradients by perturbing from baseline
          !   3. Apply step from baseline
          !   4. If step accepted: update accepted state
          !   5. If step rejected: restore to baseline (= accepted state)
          ! -------------------------------------------------------------------
          DO k=1, MAXITERRATIONS

            ! Check for convergence  (xTarg is STOP offset from VptElt in local frame)
            dxPt(:2) = xInt(:2) - xTarg(:2)               ! CALL DSUB(dxPt,xInt,xTarg,2)
            sqErr    = dxPt(1)*dxPt(1) + dxPt(2)*dxPt(2)  ! Squared Error

#if ChiefRayAiming_Debug
              !
              ! ITERATION DIAGNOSTICS
              !
              WRITE(*,*)' '
              WRITE(*,*)' ========== ITERATION', k, ' =========='
              WRITE(*,*)' Current position:'
              WRITE(*,*)' - xInt (actual)  =', xInt(1:2)
              WRITE(*,*)' - xTarg (target) =', xTarg(1:2)
              WRITE(*,*)' - dxPt (error)   =', dxPt(1:2)
              WRITE(*,*)' - sqErr          =', sqErr
              WRITE(*,*)' - SQRT(sqErr)    =', SQRT(sqErr)
              WRITE(*,*)' - errTol         =', errTol
              IF (ifPointSrc) THEN
                WRITE(*,*)' ChfRayDir =', ChfRayDir(1:3)
                WRITE(*,*)' ChfRayPos =', ChfRayPos(1:3)
                WRITE(*,*)' SourcePos =', SaveSrcPos(1:3), '(fixed)'
              ELSE
                WRITE(*,*)' ChfRayPos =', ChfRayPos(1:3)
                WRITE(*,*)' ChfRayDir =', ChfRayDir(1:3)
              END IF
              WRITE(*,*)' Perturbation: dthPt_fd =', dthPt_fd

#endif ! ChiefRayAiming_Debug

            ! Save baseline state for this iteration (for gradient computation)
            baselineChfRayDir(1:3) = ChfRayDir(1:3)
            baselineChfRayPos(1:3) = ChfRayPos(1:3)

            ! Check for convergence
            IF (sqErr.LE.errTol) THEN

#if ChiefRayAiming_Debug
                WRITE(*,*)' '
                WRITE(*,*)' *** CONVERGED: sqErr <= errTol ***'
#endif ! ChiefRayAiming_Debug

              converged = .TRUE.
              EXIT  ! Exit main iteration loop
            END IF

            ! LM: Additional practical convergence check
            ! If we've reached numerical precision limits, accept convergence
            ! sqErr < TOL_GEOM means position error < sqrt(TOL_GEOM) ≈ 1.3e-6 mm
            ! TOL_GEOM = EPS^(2/3) accounts for accumulated roundoff in geometric ops
            IF (sqErr <= errTol) THEN

#if ChiefRayAiming_Debug
                WRITE(*,*)' '
                WRITE(*,*)' *** LM CONVERGED: Reached numerical precision limit ***'
                WRITE(*,*)' sqErr =', sqErr, ' (< Tol =', errTol, ')'
                WRITE(*,*)' Position error ~', SQRT(sqErr), ' mm'
                WRITE(*,*)' (Target errTol =', errTol, ' is beyond precision)'
#endif ! ChiefRayAiming_Debug

              converged = .TRUE.
              EXIT  ! Exit main iteration loop
            END IF

            ! -----------------------------------------------------------------
            ! Compute x-gradient (with ray failure recovery)
            ! -----------------------------------------------------------------
            CALL ComputeGradient(xSrcGrid, 1, "X", grad_x_ok, GPt(:,1))
            IF (.NOT. grad_x_ok) EXIT  ! Exit main iteration loop

            ! -----------------------------------------------------------------
            ! Compute y-gradient (with ray failure recovery)
            ! -----------------------------------------------------------------
            CALL ComputeGradient(ySrcGrid, 2, "Y", grad_y_ok, GPt(:,2))
            IF (.NOT. grad_y_ok) EXIT  ! Exit main iteration loop

            !
            ! JACOBIAN ANALYSIS
            !
            det_GPt = GPt(1,1)*GPt(2,2) - GPt(1,2)*GPt(2,1)

#if ChiefRayAiming_Debug
              WRITE(*,*)' '
              WRITE(*,*)' === JACOBIAN MATRIX ==='
              WRITE(*,*)' GPt = [', GPt(1,1), GPt(1,2), ']'
              WRITE(*,*)' GPt = [', GPt(2,1), GPt(2,2), ']'
              WRITE(*,*)' det(GPt) =', det_GPt
              WRITE(*,*)' |det(GPt)| =', ABS(det_GPt)
              WRITE(*,*)' det_GPt_min =', det_GPt_min
#endif ! ChiefRayAiming_Debug

            ! -----------------------------------------------------------------
            ! GRADIENT QUALITY CHECK
            ! Check if gradients are too small (insufficient perturbation)
            ! -----------------------------------------------------------------
            grad_x_mag = SQRT(GPt(1,1)**2 + GPt(2,1)**2)
            grad_y_mag = SQRT(GPt(1,2)**2 + GPt(2,2)**2)

            ! Gradient is "good" if magnitude is reasonable relative to error
            ! For LM: Only check if gradients are above machine epsilon
            ! LM's adaptive regularization handles ill-conditioning
            ! Accept any gradient that's distinguishable from roundoff error
            grad_x_ok = (grad_x_mag > EPS)
            grad_y_ok = (grad_y_mag > EPS)

#if ChiefRayAiming_Debug
              WRITE(*,*)' '
              WRITE(*,*)' === GRADIENT QUALITY ==='
              WRITE(*,*)' X-gradient magnitude =', grad_x_mag, ' OK =', grad_x_ok
              WRITE(*,*)' Y-gradient magnitude =', grad_y_mag, ' OK =', grad_y_ok
              WRITE(*,*)' dthPt_fd (current)   =', dthPt_fd
#endif ! ChiefRayAiming_Debug


            ! If gradients are insufficient and we haven't tried too many times
            IF ((.NOT. grad_x_ok .OR. .NOT. grad_y_ok) .AND.
     &          n_grad_increase < 3 .AND. dthPt_fd < dthPt_fd_max) THEN

              dthPt_fd = min(dthPt_fd_max, dthPt_fd * 10.0d0)
              n_grad_increase = n_grad_increase + 1

#if ChiefRayAiming_Debug
                WRITE(*,*)' '
                WRITE(*,*)' *** Gradients too small - increasing perturbation ***'
                WRITE(*,*)' New dthPt_fd =', dthPt_fd
                WRITE(*,*)' Will retry gradient computation...'
#endif ! ChiefRayAiming_Debug

              CYCLE  ! Retry this iteration with larger perturbation
            END IF

            ! Reset counter if gradients are good
            IF (grad_x_ok .AND. grad_y_ok) n_grad_increase = 0

#if ChiefRayAiming_Debug
              ! Check Jacobian conditioning (robustness safeguard)
              ! Only fail if EXTREMELY singular (more lenient than before)
              IF (ABS(det_GPt) < det_GPt_min) THEN
                WRITE(*,*)' '
                WRITE(*,*)' *** WARNING: Jacobian nearly singular! ***'
                WRITE(*,*)' det(GPt) =', det_GPt, ' < threshold =', det_GPt_min
                WRITE(*,*)' Attempting to continue with regularization...'
                ! Don't reduce dthPt_fd or cycle - let regularization handle it
              END IF
#endif ! ChiefRayAiming_Debug

            ! -----------------------------------------------------------------
            ! LEVENBERG-MARQUARDT METHOD
            ! Compute J^T*J, J^T*f, condition number, and regularization
            ! -----------------------------------------------------------------
            reg_factor = DetLMRegularization(GPt, dxPt, ifUseScaling, JtJ, Jt_f, cond_est)

            ! LM damping adjustment loop - recompute step when lambda changes
            step_accepted = .FALSE.
            n_ray_failures = 0  ! Reset consecutive failure counter
            DO  ! LM inner loop
              ! Augmented system with adaptive regularization:
              ! A = J^T*J + lambda*diag(J^T*J) + reg_factor*I
              A_lm(1,1) = JtJ(1,1) * (1.0d0 + lambda_lm) + reg_factor
              A_lm(1,2) = JtJ(1,2)
              A_lm(2,1) = JtJ(2,1)
              A_lm(2,2) = JtJ(2,2) * (1.0d0 + lambda_lm) + reg_factor

              det_lm = A_lm(1,1)*A_lm(2,2) - A_lm(1,2)*A_lm(2,1)

#if ChiefRayAiming_Debug
                WRITE(*,*)' '
                WRITE(*,*)' === LEVENBERG-MARQUARDT ==='
                WRITE(*,*)' lambda (damping) =', lambda_lm
                WRITE(*,*)' det(A_lm) =', det_lm
#endif ! ChiefRayAiming_Debug

              ! Solve using Cramer's rule: A_lm * step = -Jt_f
              thx = -( A_lm(2,2)*Jt_f(1) - A_lm(1,2)*Jt_f(2)) / det_lm
              thy = -(-A_lm(2,1)*Jt_f(1) + A_lm(1,1)*Jt_f(2)) / det_lm

#if ChiefRayAiming_Debug
                WRITE(*,*)' thx (trial) =', thx
                WRITE(*,*)' thy (trial) =', thy
                WRITE(*,*)' step_norm =', SQRT(thx*thx + thy*thy)
                WRITE(*,*)' '
                WRITE(*,*)' === STEP APPLICATION DEBUG ==='
                IF (ifPointSrc) THEN
                  WRITE(*,*)' baselineChfRayDir =', baselineChfRayDir(1:3)
                ELSE
                  WRITE(*,*)' baselineChfRayPos =', baselineChfRayPos(1:3)
                END IF
#endif ! ChiefRayAiming_Debug

              ! Apply trial step from baseline state
              IF (ifPointSrc) THEN
                ! Apply angular rotations: thx around xSrcGrid, thy around ySrcGrid (Source CSYS)
                ! Rotation by angle θ around axis a: d_new ≈ d + θ*(a × d)
                ! Combined rotation: d_new = d + thx*(xSrcGrid × d) + thy*(ySrcGrid × d)

                ! Cross product: xSrcGrid × baselineChfRayDir
                DV1(1) = xSrcGrid(2)*baselineChfRayDir(3) - xSrcGrid(3)*baselineChfRayDir(2)
                DV1(2) = xSrcGrid(3)*baselineChfRayDir(1) - xSrcGrid(1)*baselineChfRayDir(3)
                DV1(3) = xSrcGrid(1)*baselineChfRayDir(2) - xSrcGrid(2)*baselineChfRayDir(1)

                ! Cross product: ySrcGrid × baselineChfRayDir
                DV2(1) = ySrcGrid(2)*baselineChfRayDir(3) - ySrcGrid(3)*baselineChfRayDir(2)
                DV2(2) = ySrcGrid(3)*baselineChfRayDir(1) - ySrcGrid(1)*baselineChfRayDir(3)
                DV2(3) = ySrcGrid(1)*baselineChfRayDir(2) - ySrcGrid(2)*baselineChfRayDir(1)

                ! Apply combined rotation from baseline
                ChfRayDir_pert(1:3) = baselineChfRayDir(1:3) +
     &                                      thx*DV1(1:3) + thy*DV2(1:3)

                ! Normalize
                norm_dir = SQRT(ChfRayDir_pert(1)*ChfRayDir_pert(1) +
     &                          ChfRayDir_pert(2)*ChfRayDir_pert(2) +
     &                          ChfRayDir_pert(3)*ChfRayDir_pert(3))
                ChfRayDir(1:3) = ChfRayDir_pert(1:3) / norm_dir
                ChfRayPos(1:3) = saveSrcPos(1:3) - zSource*ChfRayDir(1:3)
              ELSE
                ChfRayPos(1:3) = baselineChfRayPos(1:3) +
     &                                 thx*xSrcGrid(1:3) + thy*ySrcGrid(1:3)
              END IF

              ! Trace ray to evaluate trial step
              CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,ifPropagate, ifGrid)
              CALL CTRACE(RayIndex,ifBuild_,ifDMBuild,ifPol,ifDraw_,
     &                    ifLNegOK,Cmatrix,
     &                    iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &                    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &                    StartExtinc,DrawRayVec,DrawEltVec,nDrawElt, xDraw,yDraw)

              IF (LRayOK(1)) THEN  ! IF (nBadRays.EQ.0) THEN
                ! ----------------------------------------------
                ! Ray succeeded - compute trial error
                ! ----------------------------------------------
                DV1(1:3) = RayPos(1:3,1) - VptElt(1:3,iLocal) ! CALL DSUB(DV1,RayPos,VptElt(1,iLocal),3)
                xInt(1) = DDOTC(DV1, xLocal)                  ! DOT_PRODUCT(xLocal, DV1)
                xInt(2) = DDOTC(DV1, yLocal)                  ! DOT_PRODUCT(yLocal, DV1)
                dxPt(1) = xInt(1) - xTarg(1)
                dxPt(2) = xInt(2) - xTarg(2)
                sqErr_trial = dxPt(1)*dxPt(1) + dxPt(2)*dxPt(2)

#if ChiefRayAiming_Debug
                  WRITE(*,*)' '
                  WRITE(*,*)' === TRIAL STEP RESULT ==='
                  IF (ifPointSrc) THEN
                    WRITE(*,*)' ChfRayDir (trial) =', ChfRayDir(1:3)
                    WRITE(*,*)' ChfRayPos (trial) =', ChfRayPos(1:3)
                  ELSE
                    WRITE(*,*)' ChfRayPos (trial) =', ChfRayPos(1:3)
                  END IF
                  WRITE(*,*)' RayPos at STOP    =', RayPos(1:3,1)
                  WRITE(*,*)' xInt (trial)      =', xInt(1:2)
                  WRITE(*,*)' dxPt (trial err)  =', dxPt(1:2)
                  WRITE(*,*)' sqErr_trial       =', sqErr_trial
                  WRITE(*,*)' sqErr (current)   =', sqErr
                  WRITE(*,*)' Delta sqErr       =', sqErr_trial - sqErr
#endif ! ChiefRayAiming_Debug

                IF (sqErr_trial < sqErr) THEN

                  ! Accept step - error decreased
                  step_accepted = .TRUE.
                  lambda_lm = lambda_lm / 10.0d0  ! Reduce damping

                  ! Update best error and check for stagnation
                  IF (sqErr_trial < sqErr_best * 0.999d0) THEN
                    ! Significant progress (>0.1% improvement)
                    sqErr_best = sqErr_trial
                    n_stagnant = 0
                  ELSE
                    ! Minimal progress - increment stagnation counter
                    n_stagnant = n_stagnant + 1
                  END IF

#if ChiefRayAiming_Debug
                    WRITE(*,*)' '
                    WRITE(*,*)' === LM STEP ACCEPTED ==='
                    WRITE(*,*)' sqErr: ', sqErr, ' -> ', sqErr_trial
                    WRITE(*,*)' Reduced lambda to', lambda_lm
#endif ! ChiefRayAiming_Debug

                  EXIT  ! Exit LM inner loop

                ELSE

                  ! Reject step - error increased
#if ChiefRayAiming_Debug
                    WRITE(*,*)' '
                    WRITE(*,*)' === LM STEP REJECTED ==='
                    WRITE(*,*)' sqErr would increase: ', sqErr, ' -> ', sqErr_trial
#endif ! ChiefRayAiming_Debug

                  ! Check if we're at numerical precision limit
                  IF (sqErr < errTol) THEN

#if ChiefRayAiming_Debug
                      WRITE(*,*)' Already at numerical precision (sqErr < errTol)'
                      WRITE(*,*)' Accepting current state as converged'
#endif ! ChiefRayAiming_Debug

                    step_accepted = .TRUE.
                    EXIT
                  END IF

                  ! Increase damping
                  lambda_lm = lambda_lm * 10.0d0

#if ChiefRayAiming_Debug
                    WRITE(*,*)' Increased lambda to', lambda_lm
#endif ChiefRayAiming_Debug

                  IF (lambda_lm > 1d10) THEN
#if ChiefRayAiming_Debug
                      WRITE(*,*)' Lambda too large - convergence stalled'
#endif ! ChiefRayAiming_Debug

                    step_accepted = .FALSE.
                    EXIT
                  END IF

                  ! Restore to baseline and CYCLE to recompute with new lambda
                  IF (ifPointSrc) ChfRayDir(1:3) = baselineChfRayDir(1:3)
                  ChfRayPos(1:3) = baselineChfRayPos(1:3)

                  CYCLE  ! Go back to top of loop, recompute step
                END IF

              ELSE
                ! ----------------------------------------------
                ! Ray failed - increase damping and retry
                ! ----------------------------------------------
                n_ray_failures = n_ray_failures + 1

#if ChiefRayAiming_Debug
                  WRITE(*,*)' '
                  WRITE(*,*)' *** Ray trace failed with LM step ***'
#endif ! ChiefRayAiming_Debug

                ! Check if we're at numerical precision limit
                IF (sqErr < TOL_GEOM) THEN
#if ChiefRayAiming_Debug
                    WRITE(*,*)' Already at numerical precision (sqErr < TOL_GEOM)'
                    WRITE(*,*)' Accepting current state as converged'
#endif ! ChiefRayAiming_Debug

                  step_accepted = .TRUE.
                  EXIT
                END IF

                lambda_lm = lambda_lm * 10.0d0
#if ChiefRayAiming_Debug
                  WRITE(*,*)' Increased lambda to', lambda_lm
#endif ! ChiefRayAiming_Debug

                ! Abort if too many consecutive ray failures
                IF (n_ray_failures >= 5) THEN

#if ChiefRayAiming_Debug
                    WRITE(*,*)' '
                    WRITE(*,*)' *** TOO MANY CONSECUTIVE RAY FAILURES ***'
                    WRITE(*,*)' Target position appears unreachable from current state'
                    WRITE(*,*)' sqErr =', sqErr, ' (best achieved: ', sqErr_best, ')'
#endif ! ChiefRayAiming_Debug

                  step_accepted = .FALSE.
                  EXIT
                END IF

                IF (lambda_lm > 1d10) THEN
#if ChiefRayAiming_Debug
                    WRITE(*,*)' Lambda too large - cannot find valid step'
#endif ! ChiefRayAiming_Debug

                  step_accepted = .FALSE.
                  EXIT
                END IF

                ! Restore to baseline
                IF (ifPointSrc) ChfRayDir(1:3) = baselineChfRayDir(1:3)
                ChfRayPos(1:3) = baselineChfRayPos(1:3)

                CYCLE  ! Recompute with larger lambda
              END IF
            END DO  ! LM inner loop

            IF (.NOT. step_accepted) THEN

              ! Check if we failed due to numerical precision limits
              IF (sqErr < TOL_GEOM) THEN

#if ChiefRayAiming_Debug
                  WRITE(*,*)' '
                  WRITE(*,*)' LM reached numerical precision limit'
                  WRITE(*,*)' sqErr =', sqErr, ' (excellent convergence)'
#endif ! ChiefRayAiming_Debug

                converged = .TRUE.
                EXIT  ! Exit main iteration loop with success

              ELSE
#if ChiefRayAiming_Debug
                  WRITE(*,*)' '
                  WRITE(*,*)' ========================================'
                  WRITE(*,*)' ERROR: LM algorithm failed to find acceptable step'
                  WRITE(*,*)' ========================================'
                  WRITE(*,*)' Current sqErr  =', sqErr
                  WRITE(*,*)' Best sqErr     =', sqErr_best
                  WRITE(*,*)' Target errTol  =', errTol
                  WRITE(*,*)' '
                  WRITE(*,*)' Likely causes:'
                  WRITE(*,*)' - Target position is beyond physical aperture limits'
                  WRITE(*,*)' - Ray encounters vignetting or TIR'
                  WRITE(*,*)' - Geometric constraints prevent reaching target'
#endif ! ChiefRayAiming_Debug

                EXIT  ! Exit main iteration loop with failure
              END IF
            END IF

            ! Check for stagnation - too many iterations without progress
            IF (n_stagnant >= 10) THEN

#if ChiefRayAiming_Debug
                WRITE(*,*)' '
                WRITE(*,*)' ========================================'
                WRITE(*,*)' WARNING: Algorithm stagnated'
                WRITE(*,*)' ========================================'
                WRITE(*,*)' No significant progress in', n_stagnant, 'iterations'
                WRITE(*,*)' Current sqErr =', sqErr
                WRITE(*,*)' Best sqErr    =', sqErr_best
                WRITE(*,*)' Target errTol =', errTol
                WRITE(*,*)' '
#endif ! ChiefRayAiming_Debug

              IF (sqErr < TOL_GEOM) THEN

#if ChiefRayAiming_Debug
                  WRITE(*,*)' Accepting current state (within numerical precision)'
#endif ! ChiefRayAiming_Debug

                converged = .TRUE.

              ELSE

#if ChiefRayAiming_Debug
                  WRITE(*,*)' Target may be unreachable - aborting'
#endif ! ChiefRayAiming_Debug

                converged = .FALSE.

              END IF
              EXIT  ! Exit main iteration loop

            END IF

            ! -----------------------------------------------------------------
            ! Update accepted state after successful step
            ! -----------------------------------------------------------------
            IF (ifPointSrc) THEN

              ! Ensure ChfRayDir is normalized before saving as accepted state
              norm_dir = ChfRayDir(1)*ChfRayDir(1) +
     &                   ChfRayDir(2)*ChfRayDir(2) +
     &                   ChfRayDir(3)*ChfRayDir(3)
              IF (ABS(norm_dir - 1.0d0) > 2d0*EPS_SQRT) THEN
                ChfRayDir(1:3) = ChfRayDir(1:3) / sqrt(norm_dir)

#if ChiefRayAiming_Debug
                  WRITE(*,*)' WARNING: ChfRayDir not normalized! |ChfRayDir| =', norm_dir
                  WRITE(*,*)' Normalized ChfRayDir =', ChfRayDir(1:3)
#endif ! ChiefRayAiming_Debug

              END IF
              acceptedChfRayDir(:) = ChfRayDir(:)
              acceptedChfRayPos(:) = saveSrcPos(1:3) - zSource*ChfRayDir(:)
              ChfRayPos(:)         = acceptedChfRayPos(:)
            ELSE
              acceptedChfRayPos(:) = ChfRayPos(:)
            END IF

            crLocalPos(:) = RayPos(:,1)
            DV1(1:3)      = RayPos(1:3,1) - VptElt(1:3, iLocal) ! CALL DSUB(DV1,RayPos,VptElt(1,iLocal),3)
            xInt(1)       = DDOTC(xLocal, DV1)                   ! DOT_PRODUCT(xLocal, DV1)
            xInt(2)       = DDOTC(yLocal, DV1)                   ! DOT_PRODUCT(yLocal, DV1)

#if ChiefRayAiming_Debug
              !
              ! CONVERGENCE METRICS AND ADAPTIVE CONTROL
              !
              WRITE(*,*)' '
              WRITE(*,*)' === CONVERGENCE METRICS ==='
              WRITE(*,*)' sqErr (new)  =', sqErr
              WRITE(*,*)' sqErr (prev) =', sqErr_prev
#endif ! ChiefRayAiming_Debug


            ! LM handles adaptation internally - minimal post-processing
            sqErr_prev = sqErr

          END DO  ! Main iteration loop (k=1, MAXITERRATIONS)

          ! -----------------------------------------------------------------
          ! Finalize: Check convergence and set output flag
          ! -----------------------------------------------------------------
          IF (converged) THEN
            ifBeamCentered = .TRUE.

#if ChiefRayAiming_Debug
              WRITE(*,*)' '
              WRITE(*,*)' SUCCESS: Stop calculation converged!'
              WRITE(*,*)' Final sqErr =', sqErr
              WRITE(*,*)' Iterations =', k-1
#endif ! ChiefRayAiming_Debug

          ELSE
#if ChiefRayAiming_Debug
              WRITE(*,*)' '
              WRITE(*,*)' ========================================'
              WRITE(*,*)' *** CONVERGENCE FAILURE DIAGNOSTICS ***'
              WRITE(*,*)' ========================================'
              WRITE(*,*)' Final sqErr     =', sqErr
              WRITE(*,*)' Target errTol   =', errTol
              WRITE(*,*)' Ratio sqErr/tol =', sqErr/errTol
              WRITE(*,*)' SQRT(sqErr)     =', SQRT(sqErr)
              WRITE(*,*)' Iterations used =', k-1
              WRITE(*,*)' Max iterations  =', MAXITERRATIONS
              WRITE(*,*)' '
              WRITE(*,*)' Final position:'
              WRITE(*,*)' - xInt (actual)  =', xInt(1:2)
              WRITE(*,*)' - xTarg (target) =', xTarg(1:2)
              WRITE(*,*)' - dxPt (error)   =', dxPt(1:2)
              WRITE(*,*)' '
              WRITE(*,*)' Final perturbation: dthPt =', dthPt
              WRITE(*,*)' Minimum allowed: dthPt_min =', dthPt_min
              WRITE(*,*)' '
              IF (ifPointSrc) THEN
                WRITE(*,*)' Final ChfRayDir =', ChfRayDir(1:3)
              ELSE
                WRITE(*,*)' Final ChfRayPos =', ChfRayPos(1:3)
              END IF
              WRITE(*,*)' '
              WRITE(*,*)' Possible causes:'
              IF (sqErr > 1d-10) THEN
                WRITE(*,*)' - Large residual: algorithm diverged or stalled'
              END IF
              IF (dthPt < dthPt_min * 10) THEN
                WRITE(*,*)' - Small perturbation: Jacobian may be ill-conditioned'
              END IF
              IF (k >= MAXITERRATIONS) THEN
                WRITE(*,*)' - Max iterations reached: slow convergence'
              END IF
              WRITE(*,*)' '
#endif ! ChiefRayAiming_Debug

            ! Restore original state
            ChfRayPos(:) = acceptedChfRayPos(:)
            ChfRayDir(:) = acceptedChfRayDir(:)
            ! ifBeamCentered already FALSE
          END IF

          CONTAINS

          ! -----------------------------------------------------------------------
          !  SUBROUTINE: ComputeGradient
          !
          !  PURPOSE: Compute gradient column by finite difference perturbation
          !           Handles both point source (direction) and collimated (position)
          !
          !  INPUTS:
          !    pertAxis(3)  - Perturbation axis (xSrcGrid or ySrcGrid)
          !    gradCol      - Gradient column index (1 for x, 2 for y)
          !    axisName     - Axis name for debug messages ("X" or "Y")
          !
          !  OUTPUTS:
          !    gradient_ok  - .TRUE. if gradient computed successfully
          !    GPt_col(2)   - Computed gradient column GPt(:,gradCol)
          !
          !  SIDE EFFECTS:
          !    - Modifies dthPt_fd (may reduce if ray fails)
          !    - Traces rays via CTRACE
          !    - Restores baseline state before returning
          ! -----------------------------------------------------------------------
          SUBROUTINE ComputeGradient(pertAxis, gradCol, axisName,
     &                               gradient_ok, GPt_col)
            IMPLICIT NONE
            REAL(8),      INTENT(IN)  :: pertAxis(3)
            INTEGER,      INTENT(IN)  :: gradCol
            CHARACTER(*), INTENT(IN)  :: axisName
            LOGICAL,      INTENT(OUT) :: gradient_ok
            REAL(8),      INTENT(OUT) :: GPt_col(2)

            ! Local variables (ChfRayDir_pert accessed from parent via host association)
            LOGICAL :: ray_success
            REAL(8) :: ddxPt(3)

            gradient_ok = .FALSE.
            ray_success = .FALSE.

            ! Retry loop with perturbation reduction
            DO
              IF (ifPointSrc) THEN
                ! Point source: perturb direction by rotating around pertAxis
                ! For small angle dthPt_fd, rotation of ChfRayDir around pertAxis:
                ! ChfRayDir_new ≈ ChfRayDir + dthPt_fd * (pertAxis × ChfRayDir)
                ! Cross product: pertAxis × ChfRayDir
                ChfRayDir_pert(1) = pertAxis(2)*ChfRayDir(3) - pertAxis(3)*ChfRayDir(2)
                ChfRayDir_pert(2) = pertAxis(3)*ChfRayDir(1) - pertAxis(1)*ChfRayDir(3)
                ChfRayDir_pert(3) = pertAxis(1)*ChfRayDir(2) - pertAxis(2)*ChfRayDir(1)

                ! Apply rotation: d_new = d + θ*(axis × d)
                ChfRayDir(1:3) = ChfRayDir(1:3) + dthPt_fd * ChfRayDir_pert(1:3)

                ! Normalize to ensure unit vector
                norm_dir = SQRT(ChfRayDir(1)*ChfRayDir(1) +
     &                          ChfRayDir(2)*ChfRayDir(2) +
     &                          ChfRayDir(3)*ChfRayDir(3))
                ChfRayDir(1:3) = ChfRayDir(1:3) / norm_dir
                ChfRayPos(1:3) = saveSrcPos(1:3) - zSource*ChfRayDir(1:3)

              ELSE
                ! Collimated source: perturb position along pertAxis
                ChfRayPos(1:3) = ChfRayPos(1:3) + dthPt_fd*pertAxis(1:3)
              END IF

              ! Trace chief ray from source to STOP element
              CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,ifPropagate, ifGrid)
              CALL CTRACE(RayIndex,ifBuild_,ifDMBuild,ifPol,ifDraw_,
     &                    ifLNegOK,Cmatrix,
     &                    iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &                    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &                    StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

              IF (LRayOK(1)) THEN
                ray_success = .TRUE.
                EXIT  ! Success - exit retry loop
              END IF

              ! Ray failed - reduce perturbation and retry
#if ChiefRayAiming_Debug
                WRITE(*,*)' WARNING: Ray failed during ',axisName,'-gradient computation'
                WRITE(*,*)' Reducing perturbation: dthPt_fd =', dthPt_fd, '->',dthPt_fd*0.5d0
#endif ! ChiefRayAiming_Debug

              dthPt_fd = dthPt_fd * 0.5d0
              IF (dthPt_fd < dthPt_min) THEN
#if ChiefRayAiming_Debug
                  WRITE(*,*)' Perturbation too small (< dthPt_min), giving up'
#endif ! ChiefRayAiming_Debug
                EXIT  ! Exit retry loop
              END IF

              ! Restore to baseline and retry
              IF (ifPointSrc) ChfRayDir(1:3) = baselineChfRayDir(1:3)
              ChfRayPos(1:3) = baselineChfRayPos(1:3)

            END DO  ! Retry loop

            IF (.NOT. ray_success) THEN
#if ChiefRayAiming_Debug
                WRITE(*,*)' ERROR: ',axisName,'-gradient computation failed'
#endif ! ChiefRayAiming_Debug
              RETURN
            END IF

            ! Compute gradient from displacement
            DV1(1:3)   = RayPos(1:3,1) - crLocalPos(1:3)
            ddxPt(1)   = DDOTC(DV1,xLocal)
            ddxPt(2)   = DDOTC(DV1,yLocal)

            GPt_col(1) = ddxPt(1)/dthPt_fd
            GPt_col(2) = ddxPt(2)/dthPt_fd

#if ChiefRayAiming_Debug
              WRITE(*,*)' ',axisName,'-gradient: ddxPt =', ddxPt(1:2),
     &                 ' GPt(:,',gradCol,') =', GPt_col(1:2)
#endif ! ChiefRayAiming_Debug

            ! Restore to baseline for next perturbation
            IF (ifPointSrc) ChfRayDir(1:3) = baselineChfRayDir(1:3)
            ChfRayPos(1:3) = baselineChfRayPos(1:3)

            gradient_ok = .TRUE.

          END SUBROUTINE ComputeGradient

          ! -----------------------------------------------------------------------
          !  FUNCTION: DetLMRegularization
          !
          !  PURPOSE: Compute Levenberg-Marquardt system matrices and regularization
          !           Computes J^T*J, J^T*f, condition number, and adaptive regularization
          !
          !  INPUTS:
          !    GPt(2,2)     - Jacobian matrix
          !    dxPt(2)      - Error vector
          !    ifUseScaling - Use scaled condition number estimation
          !
          !  OUTPUTS (via arguments):
          !    JtJ(2,2)     - J^T * J (Hessian approximation)
          !    Jt_f(2)      - J^T * f (gradient)
          !    cond_est     - Condition number estimate
          !
          !  RETURNS:
          !    reg_factor   - Adaptive regularization factor
          ! -----------------------------------------------------------------------
          FUNCTION DetLMRegularization(GPt, dxPt, ifUseScaling,
     &                                     JtJ, Jt_f, cond_est) RESULT(reg_factor)
            IMPLICIT NONE
            REAL(8), INTENT(IN)  :: GPt(2,2), dxPt(2)
            LOGICAL, INTENT(IN)  :: ifUseScaling
            REAL(8), INTENT(OUT) :: JtJ(2,2), Jt_f(2), cond_est
            REAL(8)              :: reg_factor

            ! Local variables (moved from parent scope)
            REAL(8) :: trace_JtJ, det_JtJ, discriminant, sqrt_disc
            REAL(8) :: lambda_max, lambda_min, b2, amc
            REAL(8) :: scale, as, bs, cs

            ! ============================================================
            ! LEVENBERG-MARQUARDT METHOD
            ! Solve: (J^T*J + lambda*diag(J^T*J)) * step = -J^T*f
            ! ============================================================

            ! Compute J^T * J (Hessian approximation)
            JtJ(1,1) = GPt(1,1)*GPt(1,1) + GPt(2,1)*GPt(2,1)  ! dthx row
            JtJ(1,2) = GPt(1,1)*GPt(1,2) + GPt(2,1)*GPt(2,2)  ! cross term
            JtJ(2,1) = JtJ(1,2)                               ! symmetric
            JtJ(2,2) = GPt(1,2)*GPt(1,2) + GPt(2,2)*GPt(2,2)  ! dthy row

            b2 = JtJ(1,2)*JtJ(1,2)
            amc = JtJ(1,1) - JtJ(2,2)

            ! Compute J^T * f (gradient)
            Jt_f(1) = GPt(1,1)*dxPt(1) + GPt(2,1)*dxPt(2)
            Jt_f(2) = GPt(1,2)*dxPt(1) + GPt(2,2)*dxPt(2)

#if ChiefRayAiming_Debug
              WRITE(*,*)' '
              WRITE(*,*)' === JACOBIAN & GRADIENT DEBUG ==='
              WRITE(*,*)' GPt (Jacobian) = [', GPt(1,1), GPt(1,2), ']'
              WRITE(*,*)' GPt (Jacobian) = [', GPt(2,1), GPt(2,2), ']'
              WRITE(*,*)' dxPt (error)   = [', dxPt(1), dxPt(2), ']'
              WRITE(*,*)' JtJ(1,1) =', JtJ(1,1)
              WRITE(*,*)' JtJ(1,2) =', JtJ(1,2)
              WRITE(*,*)' JtJ(2,2) =', JtJ(2,2)
              WRITE(*,*)' Jt_f (gradient) = [', Jt_f(1), Jt_f(2), ']'
              WRITE(*,*)' scale =', scale
#endif ! ChiefRayAiming_Debug

            IF (ifUseScaling) THEN
              ! ---------------------------------------------------------------
              ! SCALED CONDITION NUMBER ESTIMATION
              ! ---------------------------------------------------------------
              ! Find scale (largest magnitude entry)
              scale = max(abs(JtJ(1,1)), abs(JtJ(1,2)), abs(JtJ(2,2)))

              if (scale <= EPS) then
                cond_est = 1d20
              else
                ! Scale matrix
                as = JtJ(1,1) / scale
                bs = JtJ(1,2) / scale
                cs = JtJ(2,2) / scale

                amc = as - cs
                b2  = bs * bs

                ! Estimate condition number and apply adaptive regularization
                ! For ill-conditioned systems, add extra damping to weak directions
                trace_JtJ = as + cs
                det_JtJ   = max(0d0, as*cs - b2) ! Clamp small negative due to roundoff

                ! Condition estimate: ratio of largest to smallest eigenvalue
                ! For 2x2: λ_max ≈ (trace + sqrt(trace²-4det))/2
                !          λ_min ≈ (trace - sqrt(trace²-4det))/2
                ! Use numerically stable formulation to avoid catastrophic cancellation

                ! Stable discriminant
                discriminant = amc*amc + 4.0d0*b2
                sqrt_disc = SQRT(discriminant)

                ! Eigenvalues (scaled: no need to scale back for "cond_est")
                lambda_max = 0.5d0 * (trace_JtJ + sqrt_disc)
                lambda_min = det_JtJ / lambda_max

                IF (lambda_min > TOL_ZERO) THEN
                  cond_est = lambda_max / lambda_min
                ELSE
                  cond_est = 1d20
                END IF

#if ChiefRayAiming_Debug
                  WRITE(*,*)' cond_est =', cond_est
#endif ! ChiefRayAiming_Debug

              end if

            ELSE
              ! ---------------------------------------------------------------
              ! UNSCALED CONDITION NUMBER ESTIMATION
              ! ---------------------------------------------------------------

              ! Estimate condition number and apply adaptive regularization
              ! For ill-conditioned systems, add extra damping to weak directions
              trace_JtJ = JtJ(1,1) + JtJ(2,2)
              det_JtJ = max(0d0, JtJ(1,1)*JtJ(2,2) - b2)  ! symmetric

              ! Condition estimate: ratio of largest to smallest eigenvalue
              ! For 2x2: λ_max ≈ (trace + sqrt(trace²-4det))/2
              !          λ_min ≈ (trace - sqrt(trace²-4det))/2
              ! Use numerically stable formulation to avoid catastrophic cancellation
              IF (det_JtJ > TOL_ZERO .AND. trace_JtJ > TOL_ZERO) THEN
                ! Stable discriminant (no cancellation)
                discriminant = amc*amc + 4d0*b2
                sqrt_disc = SQRT(discriminant)

                ! Largest eigenvalue (safe)
                lambda_max = 0.5d0 * (trace_JtJ + sqrt_disc)

                ! Smallest eigenvalue via product relation
                lambda_min = det_JtJ / lambda_max

                IF (lambda_min > TOL_ZERO) THEN
                  cond_est = lambda_max / lambda_min
                ELSE
                  cond_est = 1d20  ! Essentially singular
                END IF
              ELSE
                cond_est = 1d20  ! Essentially singular
              END IF

#if ChiefRayAiming_Debug
                WRITE(*,*)' cond_est =', cond_est
#endif ! ChiefRayAiming_Debug

            END IF

            ! Adaptive regularization factor based on conditioning
            IF (cond_est > 1d9) THEN
              reg_factor = SQRT(cond_est) * TOL_ZERO

#if ChiefRayAiming_Debug
                WRITE(*,*)' WARNING: Ill-conditioned Jacobian detected'
                WRITE(*,*)' Condition estimate =', cond_est
                WRITE(*,*)' Applying enhanced regularization =', reg_factor
#endif ! ChiefRayAiming_Debug

            ELSE IF (cond_est > 1d6) THEN
              ! Moderately ill-conditioned: apply regularization silently
              reg_factor = SQRT(cond_est) * TOL_ZERO

            ELSE
              reg_factor = 0.0d0
            END IF

          END FUNCTION DetLMRegularization

      END SUBROUTINE ChiefRayAiming


#else

        SUBROUTINE STOP(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	  ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	  iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	  nBadRays,StartIndRef,StartExtinc,Cmatrix,DrawRayVec,
     &	  DrawEltVec,iLocal,xTarg,sqErr,ifBeamCentered)

        use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
        use sourcsub_mod

	IMPLICIT NONE

      	LOGICAL ifPol,ifGrid,ifPropagate,ifRayPrint,ifLNegOK,
     &    ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),ifDMBuild,
     &    ifDraw,ifBeamCentered,ifLocalCoord,
     &    ifDraw_, ifBuild_

        INTEGER k,npts,iCurWFElt,iCurRayElt,
     &          nBadRays,iStartElt,iEndElt,iObsOpt,iStartRay,
     &          iEndRay,
     &          DrawEltVec(mDrawElt,mDrawRay),iLocal,
     &          nDrawElt

        SREAL Cmatrix(7,mCm,bRay),
     &	       DrawRayVec(2,mDrawElt,mDrawRay)

        REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &   DV1(9),DV2(9),S1,
     &   xTarg(3),xInt(2),thx,thy,dthPt,ddxPt(2),
     &   xLocal(3),yLocal(3),zLocal(3),
     &   newChfRayPos(3),oldChfRayDir(3),oldChfRayPos(3),
     &   crLocalPos(3),dxPt(2),sqErr,GPt(2,2),GPtInv(2,2),
     &   del(3),errTol,
     &   xDraw(3),yDraw(3)

        INTEGER, PARAMETER :: MaxIterrations = 100

        ! ------------------------------------------------
        ifDraw_   = .FALSE.
        ifBuild_  = .FALSE.
        dthPt     = 1d-3
        iStartRay = 1
        iEndRay   = 1
        iStartElt = 0
        iEndElt   = iLocal  ! = id of STOP element
        nMinPts   = 3
        errTol    = 1d-22

        oldChfRayDir(:) = ChfRayDir(:)  ! CALL DEQUATE(oldChfRayDir,ChfRayDir,3)
        oldChfRayPos(:) = ChfRayPos(:)  ! CALL DEQUATE(oldChfRayPos,ChfRayPos,3)
        newChfRayPos(:) = ChfRayPos(:)  ! CALL DEQUATE(newChfRayPos,ChfRayPos,3)

        ! Find nominal incidence point and local coords
        ! (Need to define local frame axes (xLocal and yLocal) at
        !  STOP element. -jzlou)

        CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw_,
     &                  ifLNegOK,ifBuild_,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &                  ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,
     &                  dxSource,nBadRays,StartIndRef,StartExtinc,
     &                  Cmatrix,DrawRayVec,DrawEltVec,
     &                  iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)

        IF (.NOT.ifLocalCoord) THEN
          WRITE(*,*)'Local coordinates could not be computed'
          go to 385
        END IF

        ! Nominal (unperturbed) incidence pt in local coords
        CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
        CALL DSUB(DV1,crLocalPos,VptElt(1,iLocal),3)
        xInt(1) = DDOTC(xLocal,DV1)
        xInt(2) = DDOTC(yLocal,DV1)

        ! Loop until desired pt is matched (max 100)
        DO 388 k=1, MaxIterrations

          ! Check for convergence  (xTarg is STOP offset from VptElt in local frame)
          dxPt(:2) = xInt(:2) - xTarg(:2)               ! CALL DSUB(dxPt,xInt,xTarg,2)
          sqErr    = dxPt(1)*dxPt(1) + dxPt(2)*dxPt(2)  ! Squared Error
          IF (sqErr.LE.errTol) GO TO 386

          ! Compute x-gradient
          CALL DSMPROD(del,xGrid,dthPt,3)         ! dx projected in xGrid
          CALL DADD(ChfRayPos,newChfRayPos,del,3) ! ChfRayPos perturbed along xGrid

          ! Trace chief ray from source to STOP element (=iEndElt)
          CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,ifPropagate, ifGrid)
          CALL CTRACE(RayIndex,ifBuild_,ifDMBuild,ifPol,ifDraw_,
     &                ifLNegOK,Cmatrix,
     &                iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &                ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &                StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
          IF (nBadRays.NE.0) GO TO 385

          CALL DSUB(DV1,RayPos,crLocalPos,3) ! diff of ChfRayPos footprint at STOP
          ddxPt(1)=DDOTC(DV1,xLocal)
          ddxPt(2)=DDOTC(DV1,yLocal)
          GPt(1,1)=ddxPt(1)/dthPt
          GPt(2,1)=ddxPt(2)/dthPt

          ! Compute y-gradient
          CALL DSMPROD(del,yGrid,dthPt,3)
          CALL DADD(ChfRayPos,newChfRayPos,del,3)

          CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,ifPropagate, ifGrid)
          CALL CTRACE(RayIndex,ifBuild_,ifDMBuild,ifPol,ifDraw_,
     &                ifLNegOK,Cmatrix,
     &                iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &                ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &                StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
          IF (nBadRays.NE.0) GO TO 385

          CALL DSUB(DV1,RayPos,crLocalPos,3)
          ddxPt(1)=DDOTC(DV1,xLocal)
          ddxPt(2)=DDOTC(DV1,yLocal)
          GPt(1,2)=ddxPt(1)/dthPt
          GPt(2,2)=ddxPt(2)/dthPt

          ! Compute next step
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

          CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,ifPropagate, ifGrid)
          CALL CTRACE(RayIndex,ifBuild_,ifDMBuild,ifPol,ifDraw_,
     &                ifLNegOK,Cmatrix,
     &                iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &                ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &                StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
          IF (nBadRays.NE.0) GO TO 385


          newChfRayPos(:) = ChfRayPos(:)           ! CALL DEQUATE(newChfRayPos,ChfRayPos,3)
          crLocalPos(:)   = RayPos(:,1)            ! CALL DEQUATE(crLocalPos,RayPos,3)
          CALL DSUB(DV1,RayPos,VptElt(1,iLocal),3)
          xInt(1)=DDOTC(xLocal,DV1)
          xInt(2)=DDOTC(yLocal,DV1)
          S1=SQRT(thx*thx+thy*thy)
          IF ((dthPt.GT.S1).AND.(dthPt.GT.5d-12)) dthPt=0.2d0*dthPt
 388    CONTINUE

        ! Here if did not converge
 385    CONTINUE
        WRITE(*,*)' WARNING: Stop calculation failed to converge!'
        ChfRayPos(:) = oldChfRayPos(:)  ! CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
        ChfRayDir(:) = oldChfRayDir(:)  ! CALL DEQUATE(ChfRayPos,oldChfRayPos,3)
        ifBeamCentered=.FALSE.
        return

        ! Here if converged
 386    CONTINUE
        ifBeamCentered=.TRUE.
        ! k=k-1
        END SUBROUTINE STOP

#endif


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

        ! Input arguments (control flags)
        LOGICAL, INTENT(IN) :: ifPol             ! Polarization flag
        LOGICAL, INTENT(IN) :: ifPropagate       ! Propagation flag
        LOGICAL, INTENT(IN) :: ifDraw            ! Draw flag
        LOGICAL, INTENT(IN) :: ifLNegOK          ! Allow negative L flag
        LOGICAL, INTENT(IN) :: ifBuild           ! Build flag
        LOGICAL, INTENT(IN) :: ifDMBuild         ! DM build flag
        LOGICAL, INTENT(IN) :: ifTrace           ! Trace flag
        LOGICAL, INTENT(IN) :: ifRayPrint        ! Ray print flag
        LOGICAL, INTENT(IN) :: ifRegrid(mElt)    ! Regrid flags per element
        LOGICAL, INTENT(IN) :: ifIntSrfSet(mElt) ! Interface surface flags

        ! Input arguments (ray trace parameters)
        INTEGER, INTENT(IN)    :: npts           ! Number of points in ray grid
        INTEGER, INTENT(IN)    :: iObsOpt        ! Observation option
        INTEGER, INTENT(IN)    :: iCurWFElt      ! Current wavefront element
        INTEGER, INTENT(IN)    :: iCurRayElt     ! Current ray element
        INTEGER, INTENT(IN)    :: iLocal         ! Local element ID (STOP) - target for trace
        INTEGER, INTENT(INOUT) :: iEndElt        ! End element for trace (not used, overridden by iLocal but update used in calling code)

        REAL*8, INTENT(IN)     :: StartIndRef    ! Starting index of refraction
        REAL*8, INTENT(IN)     :: StartExtinc    ! Starting extinction

        ! Input/Output arguments (modified by ray tracing)
        LOGICAL, INTENT(INOUT) :: ifGrid         ! Grid flag (modified by SetSourceRayGrid)

        REAL*8, INTENT(INOUT) :: Dicr            ! Grid increment i (modified by SetSourceRayGrid)
        REAL*8, INTENT(INOUT) :: Djcr            ! Grid increment j (modified by SetSourceRayGrid)
        REAL*8, INTENT(INOUT) :: dxSource        ! Source offset (modified by SetSourceRayGrid)

        SREAL, INTENT(INOUT) :: Cmatrix(7,mCm,bRay)              ! Ray matrix (modified by CRTrace)
        SREAL, INTENT(INOUT) :: DrawRayVec(2,mDrawElt,mDrawRay)  ! Draw ray vector (modified by CRTrace)
        INTEGER, INTENT(INOUT) :: DrawEltVec(mDrawElt,mDrawRay)  ! Draw element vector (modified by CRTrace)

        ! Output arguments (computed local coordinate system)
        INTEGER, INTENT(OUT) :: nBadRays         ! Number of failed rays

        REAL*8, INTENT(OUT) :: xLocal(3)         ! Local x-axis at STOP element
        REAL*8, INTENT(OUT) :: yLocal(3)         ! Local y-axis at STOP element
        REAL*8, INTENT(OUT) :: zLocal(3)         ! Local z-axis at STOP element (normal)
        REAL*8, INTENT(OUT) :: crLocalPos(3)     ! Chief ray position at STOP element

        LOGICAL, INTENT(OUT) :: ifLocalCoord     ! Success flag (TRUE if local coords computed)

        ! Local variables
        INTEGER :: iStartElt, iStartRay, iEndRay
        REAL*8 :: delMag, thMag, del(3), th(3)
        REAL*8 :: DMAG1, DMAG2, oldChfRayDir(3), oldChfRayPos(3)
        REAL*8 :: v1(3), v2(3)

	delMag            = 1d-6
	thMag             = 1d-6
	iStartRay         = 1
	iEndRay           = 1
	iStartElt         = 0
	iEndElt           = iLocal
        oldChfRayDir(1:3) = ChfRayDir(1:3)  ! CALL DEQUATE(oldChfRayDir,ChfRayDir,3)
        oldChfRayPos(1:3) = ChfRayPos(1:3)  ! CALL DEQUATE(oldChfRayPos,ChfRayPos,3)

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
          ChfRayDir(:) = oldChfRayDir(:)
          ChfRayPos(:) = oldChfRayPos(:)
          return
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
          ChfRayDir(:) = oldChfRayDir(:)
          ChfRayPos(:) = oldChfRayPos(:)
          return
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
            ChfRayDir(:) = oldChfRayDir(:)
            ChfRayPos(:) = oldChfRayPos(:)
            return
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
          ChfRayDir(:) = oldChfRayDir(:)
          ChfRayPos(:) = oldChfRayPos(:)
          return
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
            ChfRayDir(:) = oldChfRayDir(:)
            ChfRayPos(:) = oldChfRayPos(:)
            return
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
        ChfRayDir(:) = oldChfRayDir(:)
        ChfRayPos(:) = oldChfRayPos(:)

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

	ifDraw    = .FALSE.
	iStartElt = 0
	iStartRay = 1
	iEndRay   = 1
	nMinPts   = 3

	! define source with npts = nMinPts (note: it will modify nGridPts)
        !    -- Huge Overhead for simple chief-ray trace
        CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate, ifGrid)

        CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &		ifLNegOK,Cmatrix,
     &	 iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	 ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	 StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

        ! Restore  (ToDo: Needed?)
!         CALL SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
!      &                        ifPropagate,ifGrid)

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
        CALL SetSourceRay(nMinPts,Dicr,Djcr,dxSource,
     &                    xoff,yoff,ifPol,ifPropagate)

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
        use surfsub
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
     &	        ifReturn

	INTEGER i,j,iElt,iRay,nEm1,nEm2,RayIndex(2,mRay),
     &     PrevNonSeg,NextNonSeg,itst,imin,iFirstNonSeq,
     &     iLastNonSeq,IERROR

        Integer :: NextElt(0:mElt),PrevElt(0:mElt),iminPrev,
     &      jGridSrf(mElt),
     &	    iSecondElt,AdjRay(4),LastElt,iEltLoop,NextLast,
     &	    iIntSrf,jElt,iNextSeq,iNextElt,CrntNSCnt

	SREAL :: dxndxi(7,7,0:mElt),dxndui(7,6,mElt),spot(3)

	REAL*8 Lmin,na,nb,RayPosTst(3),RayLTst,PrevDir(3),
     &         PrevPos(3),
     &	       xa,xb,CurIndRef,PrevIndRef,NextIndRef,
     &	       Intxyz,CurExtinc,PrevExtinc,NextExtinc,CWU

	! DATA ifUseBaseSrf/.FALSE./

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
!-->  Reset per-ray status for this trace
	RayStatus(:) = RayStat_OK
	RayFailElt(:) = 0
	RayFailMsg(:) = ' '
	nZPFailMsg = 0
!<--
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

	    IndRef(iElt)=getIndex(CWU*Wavelen,GlassCoef(:,iElt),
     &                            GlassType(iElt))
	  END IF
C Set Zerns
	  IF ((SrfType(iElt).EQ.SrfType_Zernike).OR.(SrfType(iElt).EQ.SrfType_ZrnGrData)) THEN
            ! Zernike: ANSI Ordering
	    IF (ZernTypeL(iElt).EQ.ZernType_ANSI .OR.
     &          ZernTypeL(iElt).EQ.ZernType_NormANSI) THEN
	      CALL ZerntoMon1(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Zernike: Born & Wolf Ordering
	    ELSEIF (ZernTypeL(iElt).EQ.ZernType_BornWolf .OR.
     &              ZernTypeL(iElt).EQ.ZernType_NormBornWolf) THEN
	      CALL ZerntoMon2(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Zernike: Fringe Ordering
	    ELSEIF (ZernTypeL(iElt).EQ.ZernType_Fringe .OR.
     &              ZernTypeL(iElt).EQ.ZernType_NormFringe) THEN
	      CALL ZerntoMon3(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Special Zernike: Hex
	    ELSEIF (ZernTypeL(iElt).EQ.ZernType_NormHex) THEN
              CALL ZerntoMon4(ZernCoef(1,iElt),MonCoef(1,iElt))

            ! Zernike: Noll Ordering
            ELSEIF ((ZernTypeL(iElt).EQ.ZernType_NormNoll).OR.
     &              (ZernTypeL(iElt).EQ.ZernType_Noll)) THEN
              CALL ZerntoMon6(ZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Special Zernike: Annular Zernikes using Noll Ordering
	    ELSEIF (ZernTypeL(iElt).EQ.ZernType_NormAnnularNoll) THEN
              CALL ZerntoMon7(ZernCoef(1,iElt),MonCoef(1,iElt),
     &                        ZernAnnuRatio(iElt))
	    ELSE
	      WRITE(*,*)'** tracesub CTRACE: Elt ',iElt,
     &           ' SrfType=Zernike but ZernTypeL=',ZernTypeL(iElt),
     &           ' has no ZerntoMon implementation;',
     &           ' Zernike coefficients are being silently ignored.',
     &           ' Likely causes: ZernType= ExtFringe (no converter),',
     &           ' or Rx did not specify ZernType for this element.'
	    END IF
	  END IF

C Set FF Zernike coefficients for FreeForm surface
	  IF (SrfType(iElt).EQ.SrfType_FreeForm) THEN
            ! FF Zernike: ANSI Ordering
	    IF (FFZernTypeL(iElt).EQ.ZernType_ANSI .OR.
     &          FFZernTypeL(iElt).EQ.ZernType_NormANSI) THEN
	      CALL ZerntoMon1(FFZernCoef(1,iElt),FFCoef(1,iElt),iElt)

            ! FF Zernike: Born & Wolf Ordering
	    ELSEIF (FFZernTypeL(iElt).EQ.ZernType_BornWolf .OR.
     &              FFZernTypeL(iElt).EQ.ZernType_NormBornWolf) THEN
	      CALL ZerntoMon2(FFZernCoef(1,iElt),FFCoef(1,iElt),iElt)

            ! FF Zernike: Fringe Ordering
	    ELSEIF (FFZernTypeL(iElt).EQ.ZernType_Fringe .OR.
     &              FFZernTypeL(iElt).EQ.ZernType_NormFringe) THEN
	      CALL ZerntoMon3(FFZernCoef(1,iElt),FFCoef(1,iElt),iElt)

            ! Special FF Zernike: Hex
	    ELSEIF (FFZernTypeL(iElt).EQ.ZernType_NormHex) THEN
              CALL ZerntoMon4(FFZernCoef(1,iElt),FFCoef(1,iElt))

            ! FF Zernike: Noll Ordering
            ELSEIF ((FFZernTypeL(iElt).EQ.ZernType_NormNoll).OR.
     &              (FFZernTypeL(iElt).EQ.ZernType_Noll)) THEN
              CALL ZerntoMon6(FFZernCoef(1,iElt),FFCoef(1,iElt),iElt)

            ! Special FF Zernike: Annular Zernikes using Noll Ordering
	    ELSEIF (FFZernTypeL(iElt).EQ.ZernType_NormAnnularNoll) THEN
              CALL ZerntoMon7(FFZernCoef(1,iElt),FFCoef(1,iElt),
     &                        ZernAnnuRatio(iElt))
	    ELSE
	      WRITE(*,*)'** tracesub CTRACE: Elt ',iElt,
     &           ' SrfType=FreeForm but FFZernTypeL=',FFZernTypeL(iElt),
     &           ' has no ZerntoMon implementation;',
     &           ' FF Zernike coefficients are being silently ignored.'
	    END IF
	  END IF

C Set Mon Zernike coefficients for FreeForm surface
	  IF (SrfType(iElt).EQ.SrfType_FreeForm) THEN
            ! Mon Zernike: ANSI Ordering
	    IF (MonZernTypeL(iElt).EQ.ZernType_ANSI .OR.
     &          MonZernTypeL(iElt).EQ.ZernType_NormANSI) THEN
	      CALL ZerntoMon1(MonZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Mon Zernike: Born & Wolf Ordering
	    ELSEIF (MonZernTypeL(iElt).EQ.ZernType_BornWolf .OR.
     &              MonZernTypeL(iElt).EQ.ZernType_NormBornWolf) THEN
	      CALL ZerntoMon2(MonZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Mon Zernike: Fringe Ordering
	    ELSEIF (MonZernTypeL(iElt).EQ.ZernType_Fringe .OR.
     &              MonZernTypeL(iElt).EQ.ZernType_NormFringe) THEN
	      CALL ZerntoMon3(MonZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Special Mon Zernike: Hex
	    ELSEIF (MonZernTypeL(iElt).EQ.ZernType_NormHex) THEN
              CALL ZerntoMon4(MonZernCoef(1,iElt),MonCoef(1,iElt))

            ! Mon Zernike: Noll Ordering
            ELSEIF ((MonZernTypeL(iElt).EQ.ZernType_NormNoll).OR.
     &              (MonZernTypeL(iElt).EQ.ZernType_Noll)) THEN
              CALL ZerntoMon6(MonZernCoef(1,iElt),MonCoef(1,iElt),iElt)

            ! Special Mon Zernike: Annular Zernikes using Noll Ordering
	    ELSEIF (MonZernTypeL(iElt).EQ.ZernType_NormAnnularNoll) THEN
              CALL ZerntoMon7(MonZernCoef(1,iElt),MonCoef(1,iElt),
     &                        ZernAnnuRatio(iElt))
	    ELSE
	      WRITE(*,*)'** tracesub CTRACE: Elt ',iElt,
     &        ' SrfType=FreeForm but MonZernTypeL=',MonZernTypeL(iElt),
     &        ' has no ZerntoMon implementation;',
     &        ' Mon Zernike coefficients are being silently ignored.'
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
	  IF (nGridMat(iElt).GT.0) THEN
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

        ! ----------------------------------------------
        ! Trace rays using appropriate element routines
        ! ----------------------------------------------

            ! ----------------------------------
            ! (1) Reflector
            ! ----------------------------------

            IF (EltID(iElt).EQ.ReflectorElt) THEN
              ! PrevElt(iElt)=PrevNonSeg
              ! NextElt(PrevNonSeg)=iElt
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
     &	        pData(1,iElt),xData(1,iElt),yData(1,iElt),zData(1,iElt),
     &	        lFF(iElt),pFF(1,iElt),xFF(1,iElt),yFF(1,iElt),
     &	        zFF(1,iElt),FFCoef(1,iElt),
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

!  	      if (.false. .and. iEltLoop==21) then
!  	        print*,'*** Elt 21 Aft Reflector(),LRayOK(iRay)=',
!      &			LRayOK(iRay)
!        	      end if

	      IF (.NOT.LRayTrans) L1(iRay)=.FALSE.
	      CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)

! 	      if (.false.) then !if (iRay == (iEndRay-45)) then
! 	        print*,'**After Reflector RayDir(:,iRay) = ',
!      &            RayDir(:,iRay)
! 	        !print*,'==>iSecondElt,iEndElt = ',iSecondElt,iEndElt
! 	        print*,'===> CumRayL(iEndRay-45)=',
!      &            CumRayL(iEndRay-45),
!      &            ',RayL(iEndRay-45)=',RayL(iEndRay-45),
!      &            'SrfType(iElt)=',SrfType(iElt),
!      &            'iElt = ', iElt
! 	      end if

              PrevIndRef=CurIndRef
              PrevExtinc=CurExtinc
              ifRayElt(iElt)=.TRUE.

            ! ----------------------------------
            ! (2) Focal Plane
            ! ----------------------------------

	    ELSE IF (EltID(iElt).EQ.FocalPlaneElt) THEN
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

            ! ----------------------------------
            ! (3) Reference Surface
            ! ----------------------------------

	    ELSE IF (EltID(iElt).EQ.ReferenceElt) THEN
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

            ! ----------------------------------
            ! (11) Reflector Segment
            ! ----------------------------------

	    ELSE IF (EltID(iElt).EQ.SegmentElt) THEN
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
     &	          pData(1,iElt),xData(1,iElt),yData(1,iElt),zData(1,iElt),
     &	          lFF(iElt),pFF(1,iElt),xFF(1,iElt),yFF(1,iElt),
     &	          zFF(1,iElt),FFCoef(1,iElt),
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

            ! -------------------------------------
            ! (4) Holographic Optical Element (HOE)
            ! -------------------------------------

	    ELSE IF (EltID(iElt).EQ.HOEElt) THEN
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

            ! -------------------------------------
            ! (5) Reflective Grating
            ! -------------------------------------

	    ELSE IF (EltID(iElt).EQ.GratingElt) THEN
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
             ! CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
              CumRayL(iRay)=CumRayL(iRay)+RayL(iRay)
              PrevIndRef=CurIndRef
              PrevExtinc=CurExtinc
              ifRayElt(iElt)=.TRUE.

            ! -------------------------------------
            ! (13) Transmissive Grating
            ! -------------------------------------

	    ELSE IF (EltID(iElt).EQ.TrGratingElt) THEN
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
	!       CumRayL(iRay)=CumRayL(iRay)+CurIndRef*RayL(iRay)
	      CumRayL(iRay)=CumRayL(iRay)+RayL(iRay)
	      PrevIndRef=CurIndRef
	      PrevExtinc=CurExtinc
	      CurIndRef=IndRef(iElt)
	      CurExtinc=Extinc(iElt)
	      ifRayElt(iElt)=.TRUE.

            ! --------------------------------------
            ! (17) Diffractive Optical Element (DOE)
            ! --------------------------------------

            ELSE IF (EltID(iElt).EQ.DoeTrGratingElt) THEN
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

            ! --------------------------------------
            ! Refractor
            ! --------------------------------------

	    ELSE IF (EltID(iElt).EQ.RefractorElt) THEN
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
     &	        pData(1,iElt),xData(1,iElt),
     &	        yData(1,iElt),zData(1,iElt),
     &	        lFF(iElt),pFF(1,iElt),xFF(1,iElt),
     &	        yFF(1,iElt),zFF(1,iElt),FFCoef(1,iElt),
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

            ! --------------------------------------
            ! (7) Conic Obscuration
            ! --------------------------------------

	    ELSE IF ((EltID(iElt).EQ.ObscuringElt).AND.(ifDraw.OR.
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

            ! --------------------------------------
            ! (8) Return element
            ! --------------------------------------

	    ELSE IF (EltID(iElt).EQ.ReturnElt) THEN
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

            ! --------------------------------------
            ! (9 or 12) Nonsequential Elements
            ! --------------------------------------

            ELSE IF ((EltID(iElt).EQ.NSRefractorElt).OR.
     &               (EltID(iElt).EQ.NSReflectorElt)) THEN

        !      if (.false.) then
        !        print*,'Nonseq elem: iElt = ', iElt
        !        print*,'EltID=', EltID(iElt)
        !        !print*,'ApVec(1:4,iElt) = ', ApVec(1:4,iElt)
        !      end if


              ! Initialize search parameters
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

              ! Loop until ray hits next sequential element
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

                ! Find next element as least ray length
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
     &	             pData(1,itst),xData(1,itst),yData(1,itst),zData(1,itst),
     &	             lFF(itst),pFF(1,itst),xFF(1,itst),yFF(1,itst),
     &	             zFF(1,itst),FFCoef(1,itst),
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
     &	               pData(1,iElt),xData(1,iElt),yData(1,iElt),zData(1,iElt),
     &	               lFF(imin),pFF(1,imin),xFF(1,imin),yFF(1,imin),
     &	               zFF(1,imin),FFCoef(1,imin),
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
!-->  Record failure reason (non-sequential).  Only count this ray
!     in nBadRays if SetRayFail actually recorded a status (i.e. the
!     ray was previously OK); see elt_mod%SetRayFail.
	          IF (LZPFailed) THEN
	            IF (SetRayFail(iRay,RayStat_Bracket,iElt,
     &                  'root-find bracket/iter fail'))
     &              nBadRays=nBadRays+1
	            LZPFailed = .FALSE.
	          ELSE
	            IF (SetRayFail(iRay,RayStat_Miss,iElt,
     &                  'surface miss'))
     &              nBadRays=nBadRays+1
	          END IF
!<--
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

            ! --------------------------------------
            ! (10) Lens Array
            ! --------------------------------------

	    ELSE IF (EltID(iElt).EQ.LensArrayElt) THEN
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

! #if 0
! 	    if (iStartRay==1738) then
! 	      print*,'** iStartRay=1738: iElt =',iElt,
!      &         ', ifRayElt(iElt) =',ifRayElt(iElt),
!      &         ', EltID(iElt) =',EltID(iElt),', LRayOK(iRay)=',LRayOK(iRay),
!      &         ', ifDraw =',ifDraw
! 	    end if

! 	    if (iElt==25) then
!               print*,'** Right after End of element section,LRayOK(iRay)=',
!      &	 	LRayOK(iRay), ',  iRay=', iRay
! 	      print*,'** ifRayElt(iElt) =',ifRayElt(iElt)
! 	      print*,'EltID(iElt) =',EltID(iElt)
! 	      print*,'beamOpt%ifOptBeamPos =',beamOpt%ifOptBeamPos
! 	      print*,'beamOpt%beamPosElt =',beamOpt%beamPosElt
! 	    end if
! #endif

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
                ! If (.false. .and. iElt==7) Then
	        !  print*,'***** MACOS at elt 7: RayPos =',RayPos(1:3,iRay)
                ! End If
	        ChfRayHist%RayDirAtSrf(1:3,ChfRayHist%tot_srfs)
     &            =PrevDir(1:3)
	      End If
	    END IF

	    ! Note that LRayOK(iRay) is false only when iRay missed the surface.
            ! Obscurations and Apertures don't affect LRayOK(iRay), but they
            ! affect L1(iRay) (= LRayPass) -jzlou

	    IF (.NOT.LRayOK(iRay)) THEN
!-->  Record failure reason.  Only bump nBadRays when SetRayFail
!     actually recorded the status (first-failure wins).
	      IF (LZPFailed) THEN
	        IF (SetRayFail(iRay,RayStat_Bracket,iElt,
     &              'root-find bracket/iter fail'))
     &          nBadRays=nBadRays+1
	        LZPFailed = .FALSE.
	      ELSE
	        IF (SetRayFail(iRay,RayStat_Miss,iElt,
     &              'surface miss'))
     &          nBadRays=nBadRays+1
	      END IF
!<--
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

! #if 0
! 	  If (LFastRayTrace) Then
! 	    if (LRayOK(iRay).AND.LRayPass(iRay)) then
!                FastRayTrace(iRay)=.TRUE.
! 	    else
! 	      FastRayTrace(iRay)=.FALSE.
! 	    end if
! 	  End If
! #endif

!-->  Record obscuration (if ray passed surfaces but was blocked).
!     Bump nBadRays only on first-failure (SetRayFail is also called
!     from CPROPAGATE; an already-obscured ray must not be counted
!     again at every subsequent propagation pass).
	    IF (.NOT.L1(iRay)) THEN
	      IF (SetRayFail(iRay,RayStat_Obscured,iElt,
     &            'obscured/vignetted'))
     &        nBadRays=nBadRays+1
	    END IF
!<--

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
      WRITE(*,*) 'ZPSolve exceeded maximum iterations.'
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
     &	          L1,signpsi,SIGN,KrRS,S1,SumL,SumL2,OPDScal,DAvgL
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

          ! -------------------------------------------------------------------
          ! PURPOSE:
          !   Find the closest crossing point between two skew lines in 3D space.
          !   Computes the point on line 1 that is nearest to line 2, the minimum
          !   distance between the lines, and the parameter v along line 1.
          !   Used for entrance/exit pupil calculations and ray intersection analysis.
          !
          ! INPUTS:
          !   p1(3)      - Position vector of point on line 1
          !   r1(3)      - Direction vector of line 1 (unit vector)
          !   p2(3)      - Position vector of point on line 2
          !   r2(3)      - Direction vector of line 2 (unit vector)
          !
          ! OUTPUTS:
          !   v          - Parameter along line 1 where crossing point occurs
          !   CrossDist  - Minimum distance between the two lines
          !   CrossPt(3) - Crossing point coordinates (on line 1)
          ! -------------------------------------------------------------------

	IMPLICIT NONE

	! INPUT arguments
	REAL*8, INTENT(IN)  :: p1(3)       ! Position on line 1
	REAL*8, INTENT(IN)  :: r1(3)       ! Direction of line 1
	REAL*8, INTENT(IN)  :: p2(3)       ! Position on line 2
	REAL*8, INTENT(IN)  :: r2(3)       ! Direction of line 2

	! OUTPUT arguments
	REAL*8, INTENT(OUT) :: v           ! Parameter along line 1
	REAL*8, INTENT(OUT) :: CrossDist   ! Minimum distance between lines
	REAL*8, INTENT(OUT) :: CrossPt(3)  ! Crossing point on line 1

	! LOCAL variables
	REAL*8 :: r1r2i, Cp1(3), Cp2(3)
	REAL*8 :: rho(3), del(3), CL

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
	use elt_mod, only: RayStatus, RayFailElt,
     &    RayStat_OK, RayStat_Obscured, RayStat_Miss,
     &    RayStat_Bracket, RayStat_MaxIter, RayStat_Undef, nRay
	use surfsub, only: nZPFailMsg, mZPFailMsg
	IMPLICIT NONE
	INTEGER nElt,iElt,nBadRays,iRay
	INTEGER :: nMiss,nObs,nBrkt,nOther
	LOGICAL LEltOK(nElt)

  124	FORMAT(' One or more rays become undefined at element',i4)
  125	FORMAT(' A total of ',i8,' rays were lost')
  126	FORMAT('   Surface miss:   ',i8)
  127	FORMAT('   Obscured:       ',i8)
  128	FORMAT('   Bracket/iter:   ',i8)
  129	FORMAT('   Other:          ',i8)

	DO iElt=1,nElt
	  IF (.NOT.LEltOK(iElt)) WRITE(*,124) iElt
	END DO
	WRITE(*,125) nBadRays

!-->  Per-status summary
	nMiss=0; nObs=0; nBrkt=0; nOther=0
	DO iRay=1,nRay
	  SELECT CASE (RayStatus(iRay))
	  CASE (RayStat_Miss)
	    nMiss = nMiss + 1
	  CASE (RayStat_Obscured)
	    nObs = nObs + 1
	  CASE (RayStat_Bracket, RayStat_MaxIter)
	    nBrkt = nBrkt + 1
	  CASE (RayStat_Undef)
	    nOther = nOther + 1
	  END SELECT
	END DO
	IF (nMiss.GT.0)  WRITE(*,126) nMiss
	IF (nObs.GT.0)   WRITE(*,127) nObs
	IF (nBrkt.GT.0)  WRITE(*,128) nBrkt
	IF (nOther.GT.0) WRITE(*,129) nOther
	IF (nZPFailMsg.GT.mZPFailMsg)
     &    WRITE(*,130) mZPFailMsg, nZPFailMsg
  130	FORMAT('   (',i4,' of ',i6,
     &    ' bracket/iteration messages shown)')
!<--

	RETURN
	END SUBROUTINE WARN

C***********************************************************************

        ! ---------------------------------------------------------------------
        ! Compute index of refraction at wavelength lambda [um]
        ! ---------------------------------------------------------------------
         real(8) FUNCTION getIndex(lamda,glassProp,eltGlassType)
          implicit none
          real(8), intent(in) :: lamda,          ! Wavlength [um]
     &                           glassProp(:)    ! Glass Model Coefficients
          integer, intent(in) :: eltGlassType    ! Glass Model ID

          real(8) :: p,L2
          ! -----------------------------------------------

          SELECT CASE(eltGlassType)

            CASE (SellmeierGlass)
              ! ---------------------------------------------------------------
              ! GlassType=  Sellmeier
              ! -------------------------------
              ! Sellmeier Glass Model:
              ! n(\lambda)^2 - 1 = \sum_i frac{A_i \lambda^2}{\lambda^2 - B_i^2}
              !                                               for i=0,...,3
              !
              ! Caution:
              !  Some authors use a linear term (e.g., A) in the denominator of
              !  the Sellmeier terms and others use a squared term (e.g., A^2).
              !  The term 'A' carries the dimension of wavelength and represents
              !  the wavelength of a resonance.
              !
              ! Schott's Sellmeier Dispersion Formula:
              !
              !  n(\lambda)^2 - 1 = B1*L2/(L2-C1) + B2*L2/(L2-C2) + B3*L2/(L2-C3)
              !      with L2 = \lambda^2
              !
              ! Parameter Order: B1, B2, B3, C1, C2, C3
              ! ---------------------------------------------------------------
              L2 = lamda*lamda
              getIndex = sqrt(1d0 + glassProp(1)*L2 / (L2-glassProp(4)) +
     &                              glassProp(2)*L2 / (L2-glassProp(5)) +
     &                              glassProp(3)*L2 / (L2-glassProp(6)))


            CASE (SELLMEIER_2TA_GLASS)
              ! ---------------------------------------------------------------
              ! GlassType=  Sellm2TA
              ! -------------------------------
              ! Glass Model: n(\lamba)^2 = A + B1*L2 / (L2 - C1^2) +
              !                                B2*L2 / (L2 - C2^2)
              !
              ! Parameter Order: A, B1, B2, C1, C2
              ! ---------------------------------------------------------------
              L2 = lamda*lamda
              getIndex = sqrt(glassProp(1) +
     &                        glassProp(2)*L2 / (L2 - glassProp(4)**2) +
     &                        glassProp(3)*L2 / (L2 - glassProp(5)**2))


            CASE (KINGSLAKE1_GLASS)
              ! ---------------------------------------------------------------
              ! GlassType=  Kingslake
              ! -------------------------------
              ! Glass Model: Kettler-Drude or Helmholzz- Drude used by
              ! Li for ZnS & ZnSe
              !
              !  ==> n(lamba)^2 = A + B1 / (L2 - C1^2) + B2 / (L2 - C2^2)
              !
              ! Example:
              !  https://refractiveindex.info/?shelf=main&book=ZnS&page=Debenham
              !
              ! Parameter Order: A, B1, B2, C1, C2
              ! ---------------------------------------------------------------
              L2 = lamda*lamda
              getIndex = sqrt(glassProp(1) +
     &                       glassProp(2) / (L2 - glassProp(4)**2) +
     &                       glassProp(3) / (L2 - glassProp(5)**2))


            CASE (SchottGlass)
              ! ---------------------------------------------------------------
              ! GlassType=  Schott
              ! -------------------------------
              ! Schott = Laurent model
              !            --> "Old Schott" Model based on Laurent Series
              !
              ! SCHOTT2X4:
              !  => n(L)^2 = A + B*L^2 + C/L^2 + D/L^4 + E/L^6 + F/L^8
              !
              ! Parameter Order: A, B, C, D, E, F
              ! ---------------------------------------------------------------
              L2=lamda*lamda
              getIndex = glassProp(1) + glassProp(2)*L2 +
     &                      (glassProp(3) + (glassProp(4)+
     &                      (glassProp(5) + glassProp(6)/L2)/L2)/L2)/L2
              getIndex=sqrt(getIndex)


            CASE (CauchyGlass)
              ! ---------------------------------------------------------------
              ! GlassType= Cauchy
              ! -------------------------------
              !
              ! Cauchy Glass Model
              !
              !  => n(\lamba) = A + B/L^2 + C/L^4
              !
              ! Parameter Order: A, B, C
              ! ---------------------------------------------------------------

              L2 = lamda*lamda;
              getIndex = glassProp(1) + (glassProp(2) + glassProp(3)/L2)/L2


            CASE (HerzbergerGlass)
              ! ---------------------------------------------------------------
              ! GlassType= Herzberger
              ! -------------------------------
              !
              ! Note:
              !  HERZBERGER devised a useful equation for fitting over wide
              !  wavelength ranges. The first is a 4-term Herzberger equation.
              !  The constant 0.028 represents a short wavelength corresponding
              !  to the high frequency absorption edge, which is valid for a
              !  wide class of materials.
              !
              !  Herzberger: HERZBRGR4X2
              !
              !  n(lambda) = A + B*L^2 + C*L^4 + D*L^6 + E*p + F*p^2
              !     where p = 1/(L^2 - 0.027989)
              !
              !  --------------------------------
              ! Herzberger Formula: (using Horner's method)
              !
              !  L2 = L^2
              !  p = 1/(L^2 - 0.027989)
              !
              !  n(lambda) = A + (B + (C + D*L^2)*L^2)*L^2 + (E + F*p)*p
              !
              ! Parameter Order: A, B, C, D, E, F
              ! ---------------------------------------------------------------
              L2 = lamda*lamda
              p  = 1d0/(L2 - 0.027989d0)

              getIndex = glassProp(1) + (glassProp(2) + (glassProp(3) +
     &                                   glassProp(4)*L2)*L2)*L2 +
     &                                  (glassProp(5) + glassProp(6)*p)*p


            CASE DEFAULT
              ! => should never be reached
              getIndex = 1d22

          END SELECT

        END FUNCTION getIndex

	End Module tracesub_mod
