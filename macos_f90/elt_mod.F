C***********************************************************************
C       Begin file elt_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module elt_mod
	use param_mod
	
	Implicit none
	Save

	Integer, parameter :: mZernType=9,maxZernTerms=100,mCoat=10,
     &           	      mPrvGlass=5, ! max private glasses defined in Rx
     &                        mPrvGlassEntry=20, ! max wavlen/index entries for 
					 	 ! each glass in catalog
     &                        mES=10  ! max edge-sensors per element
	Integer, parameter :: mAsphCoef=9,  ! max terms = 9, up to term h^20
     &                        mSrfType=11   ! number of surface types as defined in prop_defs.inc 
	Integer, parameter :: mMetSrf=20,  ! max surfaces with metrology
     &                        mMetBeam=48, ! max met beams on each met srf. 
     &                        mMetPos=mMetBeam,  ! max metrology points on each met srf.
     &                        mSysMetBeam=128,    ! max total system met measurements
     &                        mRefRay=5,   ! max # reference rays, chief + 4 marginal
     &                        mPupRay=9,   ! max # pupil rays, chief + 8 marginal
     &                        mIndRef=5,   ! 
     & 			      mDOE=6       ! max # DOE elements 
        Integer :: maxAsphCoef=4  ! default to 4 terms

	Logical :: ifInMTFCalc,BWK_FLG,LStopAtEltAutoSet,LSavRayL,SrcLF_FLG,
     &             CalcSrcZFrame_FLG,DupElt_FLG
        Logical, dimension(:), allocatable :: LFastRayTrace,LRayOK,
     &                  LEltOK,LRayPass,LEltObs,ifComputeDZ,LGlass,
     &                  LTaperedElt,LApTypeDefined,LApVecDefined,
     &	                LEltMRP_FLG,LUpdateTElt_FLG,LxObsDef,
     &                  ifAmplSrf,ifOPDModGrid,ifGridDataDefined
	Logical, dimension(:,:), allocatable :: LRayToElt

        Character*32 :: EltTypeName(20), SrfTypeName(20),
     &                  ApTypeName(0:15), ObsTypeName(-10:10),
     &                  PropTypeName(15), GridTypeName(6),
     &                  ZernTypeNameL(mZernType),BaseUnits,WaveUnits

        Character*32, dimension(:), allocatable :: EltName,
     &                GlassElt, GlassName,EltMRP

        Integer :: nElt,nRay,nOutCord,nSeg,spcOption,nRefRay,nDOE,
     &             SegApType,SegApVec(3),nMetMeas

        Integer, dimension(:,:,:), allocatable :: metBeamFlg

        Integer, dimension(:,:), allocatable :: SegCoord,
     &           RayToSegMap,RayID,RayWFMap,ObsType,
     &           IWORK,EltGrp,ObsTypeRx,tMetSrfPos

        Integer, dimension(:), allocatable :: EltType,LnkElt,
     &         nECoord,EltToSegMap,SegToEltMap,PropType,ApType,
     &         nObs,iEltToiWF,SrfType,OpType,iEltToIntSrf, 
     &         nDP,LensArrayType,nCoatElt,UDSrfType, EltID,
     &         ZernaType,ZernTypeL,iEltToGridSrf,nGridMat,GridSrfOrder,
     &	       GlassType,EltCoat,nObsRx,nAmplMat,passRayIds,
     &         iEltToAmplSrf,iEltToMetSrf,nMetPos,tMetSrf,ntMetPos

	Integer, dimension(:), allocatable :: nRayToSeg,NSCnt,EdgeSensor 

        SREAL, dimension(:,:,:), allocatable :: dxidxim1,dxidui
        SREAL, dimension(:,:), allocatable :: MWFFT
        SREAL, dimension(:), allocatable :: XFP,YFP

        Real*8 :: Wavelen,Flux, zSource,Tout(7,7),focalRatio,plateScale,
     &            WavelenUM,WavelenConv,WaveBU,CBM,CWM,CWB,WavelenRx,
     &            ToutRx(7,7),sysLOS(2),sysLOSAng(2),CurRefRayDir(3)

        Real*8, dimension(:,:,:), allocatable :: TElt,ObsVec,XYZSRF,
     &          DZSRF,GridMat,ApHexBound,PolyApVtx,PolyApEqn,ObsVecRx,
     &          TEltRx,AmplMat,SrfMetPos,EdgeSensVec

	Real*8, dimension(:,:,:,:), allocatable :: ObsHexBound,
     &          PolyObsVtx,PolyObsEqn   

        Real*8, dimension(:,:), allocatable :: vptElt,
     &        RptElt, psiElt, RayPos, RayDir, ApVec, xObs,
     &        AsphCoef,MonCoef, ZernCoef,pMon,xMon,yMon,zMon,
     &        ZernCent,ZernXGrid,ZernYGrid,
     &        h1HOE,h2HOE, CoatIndxElt,CoatThkElt,AnaCoef,
     &        UDSrfCoef,yObs, zObs,DWORK,UDSrfParam,vptRx,psiRx,RptRx,
     &        GlassCoef, GlassTable,MTFMat,MTFMat2,pMonRx,xMonRx,
     &        yMonRx,zMonRx,IndRefArr,ExtincArr,EltCoatThk,ZernCoefRx,
     &	      ArrIndRef,SrfMetMea
	! IndRefArr is for coating layers, ArrIndRef is for lens element 
	! having several IndRef for different wavelengths

        Real*8, dimension(:), allocatable :: fElt,eElt,ZernRad,
     &          RayL,dxElt,IndRef,CumRayL,CumCRL,zElt,KcElt,KrElt,
     &          lMon, LensArrayWidth, OrderHOE,WaveHOE,PinHole,
     &          Extinc,RuleWidth,GridSrfdx,AmplSrfdx,KrRx,KcRx,
     &	        ZernAnnuRatio,CumRayLsav,metMeasBuf

        Real*8 :: Fn,  ! Fresnel number 
     &	          OPDdx, OPDdy, ! exit pupil OPD grid spacing
     &            RayDirBuff(3,5),rayHt,costh,efl, ! used in 'EFL' command
     &            ChfRayDirAtElt(3),ChfRayPosAtElt(3),BeamSizeAtElt,
     &	          RefRayDir(3,mRefRay)

	! Zernike normallizing coefs (Malacara) 
	Real*8 :: ZnNmFac_m(maxZernTerms)

	! Diffractive Optical Element (DOE) related arrays
	Integer, dimension(:), allocatable :: DoeId 
	Real*8 :: DoeWL(mDOE)  ! DOE design wavelength
	Real*8, dimension(:), allocatable :: DoeFl
	Real*8, dimension(:,:), allocatable :: DoePhProf  ! phase profile

        Complex*16, dimension(:,:,:), allocatable :: WFElt,JmatElt, 
     &                                               ObsMat
        Complex*16, dimension(:,:), allocatable :: RayE,WFbuff,OTFMat
	
	Integer, parameter :: maxLensGrp = 10
	Integer :: m_iCenterPMSeg,m_iFirstPMSeg,m_iLastPMSeg,
     &             m_iPM,m_iSM,m_iTM,m_iFSM,m_iPickoff, 
     &             m_iFM(2),m_iExitPupil,m_iFocalPlane, 
     &             m_nLensGrp,m_LensGrp(2,maxLensGrp), 
     &             m_LensGrpElts(maxLensGrp),refRay(mPupRay),
     &             m_err_elt
	
	Integer, parameter :: mPolySide=64  ! changed to 64 from 32 02/04/2015

	Integer, parameter :: ReflectorElt=1, FocalPlaneElt=2,
     &           ReferenceElt=3, HOEElt=4, GratingElt=5,
     &           RefractorElt=6, ObscuringElt=7,
     &           ReturnElt=8, NSRefractorElt=9,
     &           LensArrayElt=10, SegmentElt=11,
     &           NSReflectorElt=12, TrGratingElt=13,
     &           RfPolarizerElt=14, TrPolarizerElt=15
	
	Type OptBeamStr
          Logical :: ifOptBeamDir,ifOptBeamPos,ifOptBeamSize,
     &               ifOptBeamRefRayDir
          Real*8 ::  nomBeamDir(3),nomBeamPos(3),nomBeamSize,
     &               nomBeamRefRayDir(3)
	  Integer :: beamDirElt,beamPosElt,beamSizeElt,beamRefRayElt
        End Type OptBeamStr

	Type PrvGlass
	  Character(len=10) :: glsName
   	  Integer :: nEntry
	  LREAL :: IndRef(mPrvGlassEntry)
	  LREAL :: WavLen(mPrvGlassEntry)
	End Type

	Type(PrvGlass) :: PrvGlassCatalog(mPrvGlass)
	Integer :: nPrvGlass  ! actual private glasses defined in Rx


	Contains

	Subroutine elt_mod_init

        ! Dealloocate all dynamic arrays if already allocated.
        If (allocated(LRayToElt)) Then
          deallocate(LRayToElt,LRayOK,LEltOK,LRayPass,LEltObs,
     &               ifComputeDZ,LGlass,LTaperedElt,LFastRayTrace,
     &               LApTypeDefined,LApVecDefined,
     &               ifGridDataDefined,stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(logical): deallocate failed!')

          deallocate(EltName,GlassElt,GlassName,stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(char): deallocate failed!')

          deallocate(SegCoord,RayToSegMap,RayID,RayWFMap,
     &               ObsType,IWORK,LnkElt,EltMRP, LEltMRP_FLG,
     &               EltType,nECoord,EltToSegMap,SegToEltMap,
     &               PropType,ApType,nObs,iEltToiWF,SrfType,
     &		     LUpdateTElt_FLG,nObsRx,ObsTypeRx,ifAmplSrf,
     &               ifOPDModGrid,passRayIds,metBeamFlg,
     &               stat=m_err_elt)
	  if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(int): deallocate failed!')

	  deallocate(OpType, iEltToIntSrf, nDP,LensArrayType,
     &               nCoatElt,UDSrfType,EltID,ZernaType,
     &               ZernTypeL,iEltToGridSrf,nGridMat,GridSrfOrder,
     &               iEltToMetSrf,
     &               nMetPos,tMetSrf,tMetSrfPos,iEltToAmplSrf,
     &               ntMetPos,nAmplMat,stat=m_err_elt)
	  if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(int): deallocate failed!')

	  deallocate(nRayToSeg,NSCnt,EdgeSensor,stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(int): deallocate failed!')

          deallocate(dxidxim1,dxidui,MWFFT,XFP,YFP,
     &               stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(SREAL): deallocate failed!')

          deallocate(TElt,ObsVec,XYZSRF,DZSRF,GridMat,AmplMat,
     &               vptElt,RptElt,psiElt,RayPos,TEltRx,
     &               vptRx, RptRx, psiRx, ObsVecRx,LxObsDef,
     &               RayDir, ApVec, xObs, ApHexBound,ZernCoefRx,
     &               AsphCoef,MonCoef, ZernCoef,EdgeSensVec,
     &		     ZernCent,ZernXGrid,ZernYGrid,ZernRad,
     &               pMon,xMon,yMon,zMon,pMonRx,xMonRx,yMonRx,
     &               zMonRx, h1HOE,h2HOE, CoatIndxElt, CoatThkElt,
     &               PolyApVtx,PolyApEqn,PolyObsVtx,PolyObsEqn,
     &               ObsHexBound,DoeId,DoePhProf,stat=m_err_elt)
	  if (m_err_elt /= 0)
     &      call macos_memory_failure
     &             ('elt_mod_init(real*8): deallocate failed!')
     
          deallocate(AnaCoef,UDSrfCoef,yObs, zObs,DWORK,UDSrfParam,
     &               GlassCoef, GlassTable, MTFMat,MTFMat2,fElt, eElt,
     &               RayL,dxElt,IndRef,CumRayL,CumCRL,zElt,KcElt,KrElt,
     &               KrRx,KcRx,ZernAnnuRatio,CumRayLsav,
     &               metMeasBuf,stat=m_err_elt)
	  if (m_err_elt /= 0)
     &      call macos_memory_failure
     &             ('elt_mod_init(real*8): deallocate failed!')

	  deallocate(lMon,LensArrayWidth,OrderHOE,WaveHOE,PinHole,
     &               Extinc,RuleWidth,GridSrfdx,GlassType,EltGrp, 
     &               EltCoat,IndRefArr,ExtincArr,EltCoatThk,
     &               AmplSrfdx,SrfMetPos,SrfMetMea,
     &               ArrIndRef,DoeFl,stat=m_err_elt)
            if (m_err_elt /= 0)
     &      call macos_memory_failure
     &             ('elt_mod_init(real*8): deallocate failed!')

          deallocate(WFElt, JmatElt, ObsMat, RayE, WFbuff, OTFMat,
     &               stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &             ('elt_mod_init(complex*16): deallocate failed!')
        End If
	

	! Allocations

	! Logical arrays
        allocate(LRayToElt(0:mElt,bRay),LRayOK(mRay),LEltOK(mElt),
     &           LRayPass(mRay),LEltObs(mElt),ifComputeDZ(mIntSrf),
     &           LGlass(mElt),LTaperedElt(mElt),LEltMRP_FLG(0:mElt), 
     &           LApTypeDefined(mElt),LApVecDefined(mElt),
     &           LUpdateTElt_FLG(0:mElt),LFastRayTrace(mRay), 
     &           ifGridDataDefined(mGridSrf),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init: allocate(logical) failed!')
	LRayToElt = .false.; LRayOK = .false.; LEltOK = .false.
	LRayPass = .false.; LEltObs = .false.; ifComputeDZ = .false.
	LGlass = .false.; LTaperedElt=.false.
        LApTypeDefined=.false.; LApVecDefined=.false. 
        LFastRayTrace=.false.; ifGridDataDefined=.false.

        ! Character arrays
        allocate(EltName(0:mElt), GlassElt(mElt), GlassName(mGlass),
     & 	         EltMRP(0:mElt),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init: allocate(char) failed!')


        ! Integer arrays
        allocate(SegCoord(3,mElt),RayToSegMap(mRay,mSegMirs),
     &           RayID(mdttl,mdttl),RayWFMap(2,mRay),
     &           ObsType(mObs,mElt),IWORK(LIWORK,mIntSrf),
     &           EltType(0:mElt),nECoord(0:mElt),LnkElt(0:mElt),
     &           EltToSegMap(mElt),SegToEltMap(mEm1),
     &           ObsTypeRx(mObs,mElt),PropType(0:mElt),
     &           ApType(mElt),nObs(mElt),ifAmplSrf(mElt),
     &           iEltToiWF(mElt),SrfType(mElt),OpType(mElt),
     &           GlassType(mElt),EltGrp(0:mElt,0:mElt),
     &           EltCoat(0:mElt),nObsRx(mElt),ifOPDModGrid(mElt),
     &           passRayIds(mRay),DoeId(mElt),
     &           metBeamFlg(mMetPos,mMetPos,mMetSrf),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(int-1): allocate failed!')
	SegCoord=0; RayToSegMap=0; RayID=0; DoeId=0
        RayWFMap=0; ObsType=0; IWORK=0
	EltType=0; nECoord=0; EltToSegMap=0
	SegToEltMap=0; PropType=0; ApType=0
	nObs=0; iEltToiWF=0; SrfType=0; OpType=0
        LnkElt=-1; EltCoat=0; ifAmplSrf=.FALSE.; ifOPDModGrid=.FALSE.

        allocate(iEltToIntSrf(mElt),nDP(mIntSrf),
     &           LensArrayType(mElt),nCoatElt(mElt), 
     &           UDSrfType(0:mElt),EltID(0:mElt),
     &           ZernaType(mElt),ZernTypeL(mElt), 
     &           iEltToGridSrf(mElt),nGridMat(mElt),
     &           GridSrfOrder(mElt), 
     &           nRayToSeg(mRay),NScnt(mElt),nAmplMat(mElt),
     &           iEltToAmplSrf(mElt),iEltToMetSrf(mElt),
     &           nMetPos(mElt),tMetSrf(mElt),ntMetPos(mElt),
     &           EdgeSensor(mElt),stat=m_err_elt)
	if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(int-2): allocate failed!')
	iEltToIntSrf=0; nDP=0; LensArrayType=0
        nCoatElt=0; UDSrfType=0; nAmplMat=0
	EltID=0; ZernaType=0; ZernTypeL=0
	iEltToGridSrf=0; nGridMat=0; nRayToSeg=0;GridSrfOrder=1
	iEltToAmplSrf=0; NScnt=0; EdgeSensor=0

        ! SREAL arrays
        allocate(dxidxim1(7,7,mElt), dxidui(7,6,mElt),
     &           MWFFT(mdttl,mdttl), XFP(mdttl), YFP(mdttl),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(SREAL): allocate failed!')
	dxidxim1 = fzero_; dxidui = fzero_
        MWFFT = fzero_; XFP = fzero_; YFP = fzero_

	! Temporary set, need to update later ...
	mAmplMat=mGridMat; mAmplSrf=mGridSrf

        ! REAL*8 arrays
        allocate(TElt(6,6,0:mElt), ObsVec(6,mObs,mElt),
     &           XYZSRF(mDP,3,mIntSrf), DZSRF(2,mDP,mIntSrf),
     &           GridMat(mGridMat,mGridMat,mGridSrf),
     &           AmplMat(mAmplMat,mAmplMat,mAmplSrf),
     &           ApHexBound(3,6,mElt),ObsVecRx(6,mObs,mElt), 
     &           ObsHexBound(3,6,mObs,mElt),EdgeSensVec(9,mES,mElt),
     &           PolyApVtx(2,mPolySide,mElt),TEltRx(6,6,mElt),
     &           PolyApEqn(3,mPolySide,mElt),
     &           PolyObsVtx(2,mPolySide,mObs,mElt),
     &           PolyObsEqn(3,mPolySide,mObs,mElt),
     &           vptElt(3,mElt),vptRx(3,mElt),RptRx(3,mElt),
     &           psiRx(3,mElt),DoePhProf(6,mDOE),stat=m_err_elt)
	if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(real*8-0.5): allocate failed!')
	ApHexBound = -1d0; EdgeSensVec=-1d0

	allocate(RptElt(3,mElt), psiElt(3,mElt), RayPos(3,mRay),
     &           RayDir(3,mRay), ApVec(6,mElt), xObs(3,mElt),
     &           AsphCoef(mAsphCoef,mElt),MonCoef(120,mElt),
     &           ZernCoef(45,mElt),ZernCent(3,mElt),ZernRad(mElt),
     &		 ZernXGrid(3,mElt),ZernYGrid(3,mElt),
     &           pMon(3,mElt),pMonRx(3,mElt),ZernCoefRx(45,mElt),
     &           LxObsDef(mElt),tMetSrfPos(mMetPos,mElt),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(real*8-1): allocate failed!')

	TElt = dzero_; ObsVec = dzero_  
        XYZSRF = dzero_; DZSRF = dzero_
	GridMat = dzero_; vptElt = dzero_; RptElt = dzero_
	psiElt = dzero_; RayPos = dzero_; RayDir = dzero_ 
 	ApVec = dzero_; xObs = dzero_; AsphCoef = dzero_
	MonCoef = dzero_; ZernCoef = dzero_
        pMon = dzero_; AmplMat = dzero_; DoePhProf=dzero_

        allocate(xMon(3,mElt),yMon(3,mElt),zMon(3,mElt),
     &           xMonRx(3,mElt),yMonRx(3,mElt),zMonRx(3,mElt),
     &           h1HOE(3,mElt),h2HOE(3,mElt),
     &           CoatIndxElt(20,mElt),CoatThkElt(20,mElt),
     &           AnaCoef(4,mElt),
     &           UDSrfCoef(mUDScoef,mElt),yObs(3,mElt),
     &           zObs(3,mElt),DWORK(mDP,mIntSrf),UDSrfParam(10,mElt),
     &           GlassCoef(6,mElt),GlassTable(6,mGlass),
     &           MTFMat(mdttl,mdttl),MTFMat2(2*mdttl,2*mdttl),
     &           fElt(mElt), eElt(mElt),
     &           RayL(mRay), dxElt(0:mElt), IndRef(0:mElt),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(real*8-2): allocate failed!')
	xMon=dzero_; yMon = dzero_; zMon=dzero_; h1HOE=dzero_ 
	h2HOE=dzero_; CoatIndxElt=dzero_; CoatThkElt=dzero_
	AnaCoef=dzero_; UDSrfCoef=dzero_; yObs=dzero_
	zObs=dzero_; DWORK=dzero_; UDSrfParam=dzero_
	GlassCoef=dzero_; GlassTable=dzero_; fElt=dzero_ 
	eElt=dzero_; RayL=dzero_; dxElt=0; IndRef=dzero_
        MTFMat=dzero_  
	
        allocate(CumRayL(mRay),CumCRL(mElt),zElt(0:mElt),KcElt(mElt),
     &        KrElt(mElt),lMon(mElt), LensArrayWidth(mElt),
     &        OrderHOE(mElt),WaveHOE(mElt),KrRx(mElt),KcRx(mElt),
     &        PinHole(mElt),Extinc(0:mElt),AmplSrfdx(mGridSrf),
     &        RuleWidth(mElt),GridSrfdx(mGridSrf),ZernAnnuRatio(mElt),
     &        IndRefArr(0:mCoat,mElt),ExtincArr(0:mCoat,mElt),
     &        EltCoatThk(mCoat,mElt),SrfMetPos(3,mMetPos,mMetSrf),
     &        SrfMetMea(mMetBeam,mMetSrf),CumRayLsav(mElt),
     &        DoeFl(mDOE),ArrIndRef(mIndRef,mElt),
     &        metMeasBuf(mSysMetBeam),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(real*8-3): allocate failed!')
	CumRayL=dzero_; CumCRL=dzero_; zElt=dzero_; KcElt=dzero_
	KrElt = dzero_; lMon = dzero_; LensArrayWidth = dzero_
	OrderHOE = dzero_; WaveHOE = dzero_; PinHole = dzero_
	Extinc = dzero_;  RuleWidth = dzero_; GridSrfdx = dzero_
	IndRefArr=dzero_; ExtincArr=dzero_; EltCoatThk=dzero_
	AmplSrfdx=dzero_; CumRayLsav=dzero_; DoeFl=-1d0
        ArrIndRef=dzero_

        ! Complex*16 arrays
        allocate(WFElt(mdttl,mdttl,mWF), JmatElt(2,2,mElt),
     &           ObsMat(mObsMat,mObsMat,mCplxObs),
     &           RayE(3,mRay), WFbuff(mdttl,mdttl),
     &           OTFMat(mdttl,mdttl), stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(complex*16): allocate failed!')
	WFElt = czero_; JmatElt = czero_; ObsMat = czero_
	RayE = czero_; WFbuff = czero_; OTFMat=czero_

	LSavRayL=.false. 

	! Initialize Zernike normalizing factors only
	Call SetZernNormCoef(0)
        End Subroutine elt_mod_init

C  *****************************************************************

	Subroutine SetZernNormCoef(arg_flag)
        Integer :: arg_flag
	Logical, save :: first_entry=.true.

	if (arg_flag==0) first_entry=.true.  ! needed by MRESET command

	if (first_entry) then
          ! Define normalizing constants for circular Zernike terms,
          ! in Malacara order
          ZnNmFac_m(1)=1; ZnNmFac_m(2)=2; ZnNmFac_m(3)=2
          ZnNmFac_m(4)=sqrt(6d0); ZnNmFac_m(5)=sqrt(3d0)
          ZnNmFac_m(6)=sqrt(6d0); ZnNmFac_m(7)=2.*sqrt(2d0)
          ZnNmFac_m(8:10)=ZnNmFac_m(7); ZnNmFac_m(11:12)=sqrt(10d0)
          ZnNmFac_m(13)=sqrt(5d0); ZnNmFac_m(14:15)=sqrt(10d0)
          ZnNmFac_m(16:21)=2.*sqrt(3d0); ZnNmFac_m(22:24)=sqrt(14d0)
          ZnNmFac_m(25)=sqrt(7d0); ZnNmFac_m(26:28)=sqrt(14d0)
          ZnNmFac_m(29:36)=4d0; ZnNmFac_m(37:40)=3.*sqrt(2d0)
          ZnNmFac_m(41)=3d0; ZnNmFac_m(42:45)=3.*sqrt(2d0)
          ZnNmFac_m(46:66)=0.0 ! undefined beyond first 45 terms; MACOS now 
	  		       ! supports zern coef upto first 45 terms
	  first_entry=.false.
	  return
	end if
        End Subroutine SetZernNormCoef

	End Module elt_mod
