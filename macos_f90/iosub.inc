C***********************************************************************
C	Begin file iosub.inc
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.      	      |
C    +----------------------------------------------------------------+
C***********************************************************************

	SUBROUTINE LOHIN1()
	use loh_I1

	IMPLICIT NONE

C       This routine initializes the Entrance Condition Flags to .FALSE.

	ChfRayDir_FLG = .FALSE.; ChfRayPos_FLG = .FALSE.
	IndRef_FLG = .FALSE.
	Extinc_FLG = .FALSE.
	WaveUnits_FLG = .FALSE.
     	BaseUnits_FLG = .FALSE.
	Wavelen_FLG = .FALSE.
	Flux_FLG = .FALSE.
	GridType_FLG = .FALSE.
	Aperture_FLG = .FALSE.
	Obscratn_FLG = .FALSE.
	nGridpts_FLG = .FALSE.
	xGrid_FLG = .FALSE.; yGrid_FLG = .FALSE.
	nElt_FLG = .FALSE.
	nSeg_FLG = .FALSE.
	width_FLG = .FALSE.
	gap_FLG = .FALSE.
	nPetals_FLG = .FALSE.
	radCtr_FLG = .FALSE.
	SegCoord_FLG = .FALSE.
	SegXgrid_FLG = .FALSE.
	zSource_FLG = .FALSE.

	RETURN
	END SUBROUTINE LOHIN1

C***********************************************************************

	SUBROUTINE LOHIN2()
	use param_mod
	use src_mod
	use loh_I2

	IMPLICIT NONE

	INTEGER I

C       This routine initializes the Element Flags to .FALSE.

	EltName_FLG = .FALSE.
	EltType_FLG = .FALSE.
	EltID_FLG = .FALSE.
	SrfType_FLG = .FALSE.
	fElt_FLG = .FALSE.
	eElt_FLG = .FALSE.
	AsphCoef_FLG = .FALSE.
	AnaCoef_FLG = .FALSE.
	ZernCoef_FLG = .FALSE.
	ZernTypeL_FLG = .FALSE.
	MonCoef_FLG = .FALSE.
	UDSrfCoef_FLG = .FALSE.; UDSrfType_FLG = .FALSE.
	UDSrfParam_FLG = .FALSE.
	pMon_FLG = .FALSE.
	xMon_FLG = .FALSE.; yMon_FLG = .FALSE.; zMon_FLG = .FALSE.
	lMon_FLG = .FALSE.
	KrElt_FLG=.FALSE.
	KcElt_FLG=.FALSE.
	psiElt_FLG = .FALSE.
	VptElt_FLG = .FALSE.
	RptElt_FLG = .FALSE.
	IndRefE_FLG = .FALSE.
	GlassElt_FLG = .FALSE.
	GlassCoef_FLG = .TRUE.
	ExtincE_FLG = .FALSE.
	nCoat_FLG = .FALSE.
	CoatIndx_FLG = .FALSE.
	CoatThk_FLG = .FALSE.
	zElt_FLG = .FALSE.
	PropType_FLG = .FALSE.
	PinHole_FLG = .FALSE.
	nECoord_FLG = .FALSE.
	ApType_FLG = .FALSE.
	ApVec_FLG = .FALSE.
	nObs_FLG = .FALSE.
	ObsType_FLG(1:mObs) = .FALSE.
	ObsVec_FLG(1:mObs) = .FALSE.

	xObs_FLG = .FALSE.
	TElt_FLG = .FALSE.
	LensArrayType_FLG = .FALSE.
	LensArrayWidth_FLG = .FALSE.
	h1HOE_FLG = .FALSE.
	h2HOE_FLG = .FALSE.
	OrderHOE_FLG = .FALSE.
	WaveHOE_FLG = .FALSE.
	RuleWidth_FLG = .FALSE.
	nGridMat_FLG = .FALSE.
	GridMat_FLG = .FALSE.
	GridSrfdx_FLG = .FALSE.

	ZernCent_FLG  = .FALSE.

	RETURN
	END SUBROUTINE LOHIN2

C***********************************************************************

	SUBROUTINE LOHIN3()
	use loh_I3

	IMPLICIT NONE

C   This routine initializes the End Condition Flags to .FALSE.
	nOutCord_FLG = .FALSE.
	Tout_FLG = .FALSE.

	RETURN
	END SUBROUTINE LOHIN3

C***********************************************************************

	SUBROUTINE ChkDf1(npts,ENTR_OK)
	use param_mod
	use src_mod
	use elt_mod
	use loh_I1

	IMPLICIT NONE

	INTEGER npts
	REAL*8 zaxis(3)

C       Setup for new input scheme

	LOGICAL ENTR_OK

  111	FORMAT(A,' must be specified by user')
  112	FORMAT(' Default used for ',A)

	ENTR_OK = .TRUE.

C       Check for things user must specify
C       Set defaults, if necessary

	IF (ChfRayDir_FLG) THEN
	  CONTINUE
	ELSE
	  ChfRayDir(:)  = (/0d0, 0d0, 1d0/)
	  ChfRayDir_FLG = .TRUE.
	  WRITE(*,112)'ChfRayDir'
	END IF

	IF (ChfRayPos_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(ChfRayPos,3)
	  ChfRayPos_FLG=.TRUE.
	  WRITE(*,112)'ChfRayPos'
	END IF

	IF (IndRef_FLG) THEN
	  CONTINUE
	ELSE
	  IndRef(0)=1.0D0
	  IndRef_FLG=.TRUE.
	  WRITE(*,112)'IndRef(0)'
	END IF

	IF (Extinc_FLG) THEN
	  CONTINUE
	ELSE
	  Extinc(0)=0.0D0
	  Extinc_FLG=.TRUE.
c	  WRITE(*,112)'Extinc(0)'
	END IF

	IF (BaseUnits_FLG) THEN
	  CONTINUE
	ELSE
	  BaseUnits='(default-unit)'
	  BaseUnits_FLG=.TRUE.
	  WRITE(*,112)'BaseUnits'
	END IF

	IF (WaveUnits_FLG) THEN
	  CONTINUE
	ELSE
	  WaveUnits='none'
	  WaveUnits_FLG=.TRUE.
	  WRITE(*,112)'WaveUnits'
	END IF

	IF (Wavelen_FLG) THEN
	  CONTINUE
	ELSE
	  Wavelen=6.0D-7
	  Wavelen_FLG=.TRUE.
	  WRITE(*,112)'Wavelen'
	END IF

	IF (Flux_FLG) THEN
	  CONTINUE
	ELSE
	  Flux=1.0D0
	  Flux_FLG=.TRUE.
	  WRITE(*,112)'Flux'
	END IF

	IF (GridType_FLG) THEN
	  CONTINUE
	ELSE
	  GridType=1
	  GridType_FLG=.TRUE.
	  WRITE(*,112)'GridType'
	END IF

	IF (Aperture_FLG) THEN
	  CONTINUE
	ELSE
	  Aperture=1.0D0
	  Aperture_FLG=.TRUE.
	  WRITE(*,112)'Aperture'
	END IF

	IF (Obscratn_FLG) THEN
	  CONTINUE
	ELSE
	  Obscratn=0.0D0
	  Obscratn_FLG=.TRUE.
	  WRITE(*,112)'Obscratn'
	END IF

	IF (nGridpts_FLG) THEN
	  CONTINUE
	ELSE
	  nGridpts=mpts/2
	  npts=nGridpts-1
	  nGridpts_FLG=.TRUE.
	  WRITE(*,112)'nGridpts'
	END IF

	IF (xGrid_FLG) THEN
	  CONTINUE
	ELSE
	  IF (ChfRayDir_FLG) THEN
	    xGrid(1)=ChfRayDir(3)
	    xGrid(2)=ChfRayDir(1)
	    xGrid(3)=ChfRayDir(2)
	    xGrid_FLG=.TRUE.
	    WRITE(*,112)' xGrid'
	  ELSE
	    WRITE(*,111)' xGrid'
	    ENTR_OK=.FALSE.
	    RETURN
	  END IF
	END IF

	IF (yGrid_FLG) THEN
	  CONTINUE
	ELSE
	  IF (xGrid_FLG) THEN
	    CALL DNEGATE(zaxis,ChfRayDir,3)
	    CALL DXPROD(yGrid,zaxis,xGrid)
	    CALL DUNITIZE(yGrid)
	    CALL DXPROD(xGrid,yGrid,zaxis)
	    yGrid_FLG=.TRUE.
	    WRITE(*,112)' yGrid'
	  ELSE
	    WRITE(*,111)' yGrid'
	    ENTR_OK=.FALSE.
            RETURN
	  END IF
	END IF

	IF (nElt_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' nElt'
	  ENTR_OK=.FALSE.
          RETURN
	END IF

	IF (nSeg_FLG) THEN
	  CONTINUE
	ELSE
	  IF ((GridType.EQ.3).OR.(GridType.EQ.4)) THEN
	    WRITE(*,111)' nSeg'
	    ENTR_OK=.FALSE.
            RETURN
	  ELSE
	    nSeg=0
	    nSeg_FLG=.TRUE.
	    WRITE(*,112)'nSeg'
	  END IF
	END IF

	IF (width_FLG) THEN
	  CONTINUE
	ELSE
	  width=Aperture
	  width_FLG=.TRUE.
	  WRITE(*,112)'width'
	END IF

	IF (gap_FLG) THEN
	  CONTINUE
	ELSE
	  gap=0.0D0
	  gap_FLG=.TRUE.
	  WRITE(*,112)'gap'
	END IF

	IF (nPetals_FLG) THEN
	  CONTINUE
	ELSE
	  nPetals=nSeg-1
	  nPetals_FLG=.TRUE.
	  WRITE(*,112)'nPetals'
	END IF

	IF (radCtr_FLG) THEN
	  CONTINUE
	ELSE
	  radCtr=Aperture/3d0
	  radCtr_FLG=.TRUE.
	  WRITE(*,112)'radCtr'
	END IF

	IF (SegCoord_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' SegCoord'
	  ENTR_OK=.FALSE.
          RETURN
	END IF

	IF (SegXgrid_FLG) THEN
	  CONTINUE
	ELSE
	  IF (xGrid_FLG) THEN
	    CALL DEQUATE(SegXgrid(1),xGrid(1),3)
	    SegXgrid_FLG=.TRUE.
	    WRITE(*,112)' SegXgrid'
	  END IF
	END IF

	IF (zSource_FLG) THEN
	  CONTINUE
	ELSE
	  zSource=1D19
	  zSource_FLG=.TRUE.
	  WRITE(*,112)'zSource'
	END IF

	RETURN
	END SUBROUTINE ChkDf1

C***********************************************************************

	SUBROUTINE ChkDf2(iElt,ELEM_OK)

cHistory
c  24Jul97 loh:  Added flag to indicate zernike type
c  24Jul97 loh:  Fixed bug so that when the xObs default is used, the
c		 warning message is printed correctly, and fixed it so
c		 you don't get the irritating "Default Used for xObs"
c		 warning messages when you have nObs=0

	use param_mod
	use src_mod
	use elt_mod
	use loh_I2

	IMPLICIT NONE

	INTEGER i,iElt

C   Check for things user must specify
C   Setup for new input scheme

	LOGICAL ELEM_OK

  111	FORMAT(A,'(',i3,') must be specified by user')
  112	FORMAT(' Default used for ',A,'(',i4,')')
  113	FORMAT(' Too many obscurations specified -- truncated to',i3)

	ELEM_OK = .TRUE.

C       Set defaults, if necessary

	IF (EltName_FLG) THEN
	  CONTINUE
	ELSE
	  EltName(iElt)='Elt'
	  EltName_FLG=.TRUE.
	  WRITE(*,112)'EltName',iElt
	END IF

	IF (EltType_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' EltType',iElt
	  ELEM_OK = .FALSE.
	END IF

	IF (EltID_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' Element',iElt
	  ELEM_OK = .FALSE.
	END IF

	IF (SrfType_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' Surface',iElt
	  ELEM_OK = .FALSE.
	END IF

	IF (fElt_FLG) THEN
	  CONTINUE
	ELSE
	  fElt(iElt)=1D22
	  fElt_FLG=.TRUE.
	  WRITE(*,112)'fElt',iElt
	END IF

	IF (eElt_FLG) THEN
	  CONTINUE
	ELSE
	  eElt(iElt)=0.0D0
	  eElt_FLG=.TRUE.
	  WRITE(*,112)'eElt',iElt
	END IF

	IF (psiElt_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DNEGATE(psiElt(1,iElt),ChfRayDir(1),3)
	  psiElt_FLG=.TRUE.
	  WRITE(*,112)'psiElt',iElt
	END IF

	IF (AsphCoef_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(AsphCoef(1,iElt),4)
	  AsphCoef_FLG=.TRUE.
	  WRITE(*,112)'AsphCoef',iElt
	END IF

	IF (AnaCoef_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(AnaCoef(1,iElt),4)
	  AnaCoef_FLG=.TRUE.
	  WRITE(*,112)'AnaCoef',iElt
	END IF

	IF (ZernCoef_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(ZernCoef(1,iElt),45)
	  ZernCoef_FLG=.TRUE.
	  WRITE(*,112)'ZernCoef',iElt
	END IF

	IF (ZernTypeL_FLG) THEN
	  CONTINUE
	ELSE
	  ZernTypeL_FLG=.TRUE.
	  ZernTypeL(iElt)=1
	  WRITE(*,112)'ZernType',iElt
	END IF

	IF (MonCoef_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(MonCoef(1,iElt),120)
	  MonCoef_FLG=.TRUE.
	  WRITE(*,112)'MonCoef',iElt
	END IF

	IF (UDSrfCoef_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(UDSrfCoef(1,iElt),mUDScoef)
	  UDSrfCoef_FLG=.TRUE.
	  WRITE(*,112)'UDSrfCoef',iElt
	END IF

	IF (UDSrfType_FLG) THEN
	  CONTINUE
	ELSE
	  IF (SrfType(iElt).EQ.7) THEN
	    UDSrfType(iElt)=1
	    UDSrfType_FLG=.TRUE.
	    WRITE(*,112)'UDSrfType',iElt
	  ELSE
	    UDSrfType(iElt)=0
	    UDSrfType_FLG=.TRUE.
	  END IF
	END IF

	IF (UDSrfParam_FLG) THEN
	  CONTINUE
	ELSE
	  DO i=1,6
	    UDSrfParam(i,iElt)=0d0
	  END DO
	  IF (SrfType(iElt).EQ.7) THEN
	    IF (UDSrfType(iElt).EQ.1) THEN
C	Param order = beta,dact or drad,dq,minrad
	      UDSrfParam(1,iElt)=6.51399619601475d0
	      UDSrfParam(2,iElt)=2.7559d-01
	    ELSE IF (UDSrfType(iElt).EQ.2) THEN
	      UDSrfParam(1,iElt)=0.00375d0
	      UDSrfParam(2,iElt)=4.7872e+02
	    ELSE IF (UDSrfType(iElt).EQ.3) THEN
	      UDSrfParam(1,iElt)=2.5645654315019d0
	      UDSrfParam(2,iElt)=0.7d0
	    ELSE IF (UDSrfType(iElt).EQ.4) THEN
	      UDSrfParam(1,iElt)=2.5645654315019d0
	      UDSrfParam(2,iElt)=0.7d0
	    ELSE IF (UDSrfType(iElt).EQ.5) THEN
	      UDSrfParam(1,iElt)=0.00375d0
	      UDSrfParam(2,iElt)=4.7872d2
	    ELSE IF (UDSrfType(iElt).EQ.6) THEN
	      UDSrfParam(1,iElt)=6.51399619601475d0
	      UDSrfParam(2,iElt)=2.7559d-01
	    ELSE IF (UDSrfType(iElt).EQ.7) THEN
	      UDSrfParam(1,iElt)=0d0
	      UDSrfParam(2,iElt)=1.683d0
	    ELSE IF (UDSrfType(iElt).EQ.8) THEN
	      UDSrfParam(1,iElt)=0d0
	      UDSrfParam(2,iElt)=0.302d0
	    ELSE IF (UDSrfType(iElt).EQ.9) THEN
	      UDSrfParam(1,iElt)=0.85666666666667d0
	      UDSrfParam(2,iElt)=6.28318530717959d0/16
	      UDSrfParam(3,iElt)=15.87166666666667d0
	    ELSE IF (UDSrfType(iElt).EQ.10) THEN
	      UDSrfParam(1,iElt)=2.5645654315019d0
	      UDSrfParam(2,iElt)=0.7d0
	      UDSrfParam(3,iElt)=21
	    ELSE IF (UDSrfType(iElt).EQ.11) THEN
	      UDSrfParam(1,iElt)=0d0
	      UDSrfParam(2,iElt)=1.380435d0
	    END IF
	    WRITE(*,112)'UDSrfParam',iElt
	  END IF
	  UDSrfParam_FLG=.TRUE.
	END IF

	IF (KrElt_FLG) THEN
	  CONTINUE
	ELSE
	  KrElt(iElt)=-1D22
	  KrElt_FLG=.TRUE.
	  WRITE(*,112)'KrElt',iElt
	END IF

	IF (KcElt_FLG) THEN
	  CONTINUE
	ELSE
	  KcElt(iElt)=0.0D0
	  KcElt_FLG=.TRUE.
	  WRITE(*,112)'KcElt',iElt
	END IF

	IF (pMon_FLG) THEN
	  CONTINUE
	ELSE
	  IF (VptElt_FLG) THEN
	    CALL DEQUATE(pMon(1,iElt),VptElt(1,iElt),3)
	    pMon_FLG=.TRUE.
	    WRITE(*,112)'pMon',iElt
	  ELSE
	    ELEM_OK=.FALSE.
	  END IF
	END IF

	IF (xMon_FLG) THEN
	  CONTINUE
	ELSE
	  xMon(1,iElt)=psiElt(3,iElt)
	  xMon(2,iElt)=psiElt(1,iElt)
	  xMon(3,iElt)=psiElt(2,iElt)
	  xMon_FLG=.TRUE.
	  WRITE(*,112)'xMon',iElt
	END IF

	IF (yMon_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DXPROD(yMon(1,iElt),psiElt(1,iElt),xMon)
	  yMon_FLG=.TRUE.
	  WRITE(*,112)'yMon',iElt
	END IF

	IF (zMon_FLG) THEN
	  CONTINUE
	ELSE
	  IF (psiElt_FLG) THEN
	    CALL DEQUATE(zMon(1,iElt),psiElt(1,iElt),3)
	    zMon_FLG=.TRUE.
	    WRITE(*,112)'zMon',iElt
	  ELSE
	    ELEM_OK=.FALSE.
	  END IF
	END IF

	IF (lMon_FLG) THEN
	  CONTINUE
	ELSE
	  lMon(iElt)=Aperture
	  lMon_FLG=.TRUE.
	  WRITE(*,112)'lMon',iElt
	END IF

	IF (LensArrayType_FLG) THEN
	  CONTINUE
	ELSE
	  LensArrayType(iElt)=1
	  LensArrayType_FLG=.TRUE.
	  WRITE(*,112)'LensArrayType',iElt
	END IF

	IF (LensArrayWidth_FLG) THEN
	  CONTINUE
	ELSE
	  LensArrayWidth(iElt)=1
	  LensArrayWidth_FLG=.TRUE.
	  WRITE(*,112)'LensArrayWidth',iElt
	END IF

	IF (h1HOE_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' h1HOE',iElt
	  ELEM_OK = .FALSE.
	END IF

	IF (h2HOE_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' h2HOE',iElt
	  ELEM_OK = .FALSE.
	END IF

	IF (OrderHOE_FLG) THEN
	  CONTINUE
	ELSE
	  OrderHOE(iElt)=1.0
	  OrderHOE_FLG=.TRUE.
	  WRITE(*,112)' OrderHOE',iElt
	END IF

	IF (RuleWidth_FLG) THEN
	  CONTINUE
	ELSE
	  RuleWidth(iElt)=1.0
	  RuleWidth_FLG=.TRUE.
	  WRITE(*,112)'RuleWidth',iElt
	END IF

	IF (WaveHOE_FLG) THEN
	  CONTINUE
	ELSE
	  WaveHOE(iElt)=1d-7
	  WaveHOE_FLG=.TRUE.
	  WRITE(*,112)'WaveHOE',iElt
	END IF

	IF (VptElt_FLG) THEN
	  CONTINUE
	ELSE
	  WRITE(*,111)' VptElt',iElt
	  ELEM_OK = .FALSE.
	END IF

	IF (RptElt_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DEQUATE(RptElt(1,iElt),VptElt(1,iElt),3)
	  RptElt_FLG=.TRUE.
	  WRITE(*,112)'RptElt',iElt
	END IF

	IF (IndRefE_FLG) THEN
	  CONTINUE
	ELSE
	  IndRef(iElt)=1.0D0
	  IndRefE_FLG=.TRUE.
	  WRITE(*,112)'IndRef',iElt
	END IF

	IF (GlassElt_FLG) THEN
	  CONTINUE
	ELSE
	  GlassElt(iElt)='Air'
	  GlassElt_FLG=.TRUE.
	  WRITE(*,112)'GlassElt',iElt
	END IF

	IF (GlassCoef_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(GlassCoef(1,iElt),6)
	  GlassCoef_FLG=.TRUE.
	  WRITE(*,112)'GlassCoef',iElt
	END IF

	IF (ExtincE_FLG) THEN
	  CONTINUE
	ELSE
	  IF ((EltID(iElt).EQ.6).OR.(EltID(iElt).EQ.9)) THEN
	    Extinc(iElt)=0.0D0
	  ELSE
	    Extinc(iElt)=1D22
	  END IF
	  ExtincE_FLG=.TRUE.
c	  WRITE(*,112)'Extinc',iElt
	END IF

	IF (nCoat_FLG) THEN
	  CONTINUE
	ELSE
	  nCoatElt(iElt)=0
	  nCoat_FLG=.TRUE.
	END IF

	IF (CoatIndx_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(CoatIndxElt(1,iElt),20)
	  CoatIndx_FLG=.TRUE.
	END IF

	IF (CoatThk_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(CoatThkElt(1,iElt),20)
	  CoatIndx_FLG=.TRUE.
	END IF

	IF (zElt_FLG) THEN
	  CONTINUE
	ELSE
	  zElt(iElt)=fElt(iElt)
	  zElt_FLG=.TRUE.
	  WRITE(*,112)'zElt',iElt
	END IF

	IF (PropType_FLG) THEN
	  CONTINUE
	ELSE
	  IF ((EltID(iElt).EQ.3).OR.(EltID(iElt).EQ.8)) THEN
	    PropType(iElt)=3
	  ELSE
	    PropType(iElt)=1
	  END IF
	  PropType_FLG=.TRUE.
	  WRITE(*,112)'PropType',iElt
	END IF

	IF (PinHole_FLG) THEN
	  CONTINUE
	ELSE
	  PinHole_FLG=.TRUE.
	  PinHole(iElt)=1d22
	  IF ((PropType(iElt).EQ.7).OR.(PropType(iElt).EQ.8))
     &	  THEN
	    WRITE(*,112)'PinHole',iElt
	  END IF
	END IF

	IF (nECoord_FLG) THEN
	  CONTINUE
	ELSE
	  nECoord(iElt)=-6
	  nECoord_FLG=.TRUE.
	  WRITE(*,112)'nECoord',iElt
	END IF

	IF (ApType_FLG) THEN
	  CONTINUE
	ELSE
	  ApType(iElt)=0
	  ApType_FLG=.TRUE.
	  ApVec_FLG=.TRUE.
C	  WRITE(*,112)'ApType',iElt
	END IF

	IF (ApType(iElt).NE.0) THEN
	  IF (ApVec_FLG) THEN
	    CONTINUE
	  ELSE
	    IF (ApType(iElt).EQ.2) THEN
	      ApVec(1,iElt)=ChfRayPos(1)-Aperture/2
	      ApVec(2,iElt)=ChfRayPos(1)+Aperture/2
	      ApVec(3,iElt)=ChfRayPos(2)-Aperture/2
	      ApVec(4,iElt)=ChfRayPos(2)+Aperture/2
	    ELSE
	      ApVec(1,iElt)=Aperture/2
	      ApVec(2,iElt)=ChfRayPos(1)
	      ApVec(3,iElt)=ChfRayPos(2)
	    END IF
	    ApVec_FLG=.TRUE.
  	  END IF
  	END IF

	IF (nObs_FLG) THEN
	  CONTINUE
	ELSE
	  nObs(iElt)=0
	  nObs_FLG=.TRUE.
C	  WRITE(*,112)'nObs',iElt
	END IF

	IF (nObs(iElt) .NE. 0) THEN
	  IF (nObs(iElt).GT.mObs) THEN
	    nObs(iElt)=mObs
	    WRITE(*,113)mObs
	  END IF
	  DO 11 I=1,nObs(iElt)
	    IF (ObsType_FLG(I)) THEN
	      CONTINUE
	    ELSE
	      ObsType(I,iElt)=1
	      ObsType_FLG(I) =.TRUE.
	      WRITE(*,112)'ObsType',iElt
	    END IF
   11     CONTINUE
	END IF

	IF (nObs(iElt) .NE. 0) THEN
	    DO 12 I=1,nObs(iElt)
	    IF (ObsVec_FLG(I)) THEN
	      CONTINUE
	    ELSE
	      CALL DZERO(ObsVec(1,I,iElt),4)
	      ObsVec_FLG(I) =.TRUE.
	      WRITE(*,112)'ObsVec',iElt
	    END IF
   12     CONTINUE
	END IF

	IF ((nObs(iElt).NE.0).OR.(ApType(iElt).NE.0)) THEN
	  IF (xObs_FLG) THEN
	    CONTINUE
	  ELSE
	    IF (psiElt_FLG) THEN
	      xObs(1,iElt)=psiElt(3,iElt)
	      xObs(2,iElt)=psiElt(1,iElt)
	      xObs(3,iElt)=psiElt(2,iElt)
	      xObs_FLG=.TRUE.
	      WRITE(*,112)' xObs',iElt
	    ELSE
	      WRITE(*,111)' xObs',iElt
	      ELEM_OK=.FALSE.
	    END IF
	  END IF
	END IF

	IF (TElt_FLG) THEN
	  CONTINUE
	ELSE
	  CALL DZERO(TElt(1,1,iElt),36)
	  DO 13 I=1,6
	      TElt(I,I,iElt)=1.0D0
   13     CONTINUE
	  TElt_FLG=.TRUE.
	  WRITE(*,112)'TElt',iElt
	END IF

	IF ((fElt(iElt) .EQ. 0.0) .AND. (eElt(iElt) .EQ. 0.0)) THEN
	  IF (KcElt(iElt) .GT. 0.0) THEN
	    eElt(iElt) = 0d0
	    fElt(iElt) = ABS(KrElt(iElt))
	  ELSE
	    eElt(iElt) = DSQRT(-KcElt(iElt))
	    fElt(iElt) = -KrElt(iElt)/(1.0+eElt(iElt))
	  END IF
	ELSE
	  KcElt(iElt)=-eElt(iElt)*eElt(iElt)
	  KrElt(iElt)=-(1.0 + eElt(iElt))*fElt(iElt)
	END IF

	RETURN
	END SUBROUTINE ChkDf2

C***********************************************************************

	SUBROUTINE ChkDf3(EXIT_OK)
	use param_mod
	use src_mod
	use elt_mod
	use loh_I3

	IMPLICIT NONE


C       Setup for new input scheme

	LOGICAL EXIT_OK
	INTEGER iCoord

  111	FORMAT(A,' must be specified by user')
  112	FORMAT(' Default used for ',A)

	EXIT_OK = .TRUE.

C       Set defaults, if necessary

	IF (nOutCord_FLG) THEN
	  CONTINUE
	ELSE
	  nOutCord=5
	  nOutCord_FLG=.TRUE.
	  WRITE(*,112)'nOutCord'
	END IF

	IF (Tout_FLG) THEN
	  CONTINUE
	ELSE
	  IF (nOutCord.EQ.5) THEN
	    CALL DZERO(Tout(1,1),49)
	    Tout(1,1)=-psiElt(3,nElt)
	    Tout(1,2)=-psiElt(1,nElt)
	    Tout(1,3)=-psiElt(2,nElt)
	    Tout(3,4)=-psiElt(3,nElt)
	    Tout(3,5)=-psiElt(1,nElt)
	    Tout(3,6)=-psiElt(2,nElt)
	    Tout(2,1)=-psiElt(2,nElt)
	    Tout(2,2)=-psiElt(3,nElt)
	    Tout(2,3)=-psiElt(1,nElt)
	    Tout(4,4)=-psiElt(2,nElt)
	    Tout(4,5)=-psiElt(3,nElt)
	    Tout(4,6)=-psiElt(1,nElt)
	    Tout(5,7)=1
	  ELSE
 	    WRITE(*,111)' Tout'
	    EXIT_OK = .FALSE.
	  END IF
	END IF

	RETURN
	END SUBROUTINE ChkDf3

C***********************************************************************

	SUBROUTINE SumData(GridType,npts)
	use param_mod
	use elt_mod

	IMPLICIT NONE

	INTEGER GridType,npts,iElt,nGridPts

 118	FORMAT(' Optical Train Summary:'/
     &	'  Tracing rays past',i4,' elements.'/
     &	'  Aperture grid type = ',A,'with ',i4,' grid points.')
 119	FORMAT('  Aperture is segmented into',i4,' segments.')
 120	FORMAT(1P,'  Elt',i4,': ',a12,': ',a8,' with Kc=',d11.4,
     &	', Kr=',d11.4,', nECoords =',i4)
 121	FORMAT('  Number of output coordinates is',i4)

	nGridPts=npts+1

	WRITE (*,118) nElt,GridTypeName(GridType)(1:10),nGridPts
	IF (GridType.GE.3) THEN
	  WRITE (*,119) nSeg
	END IF
	DO 2 iElt=1,nElt
	  WRITE (*,120) iElt,EltTypeName(EltID(iElt)),
     &	  SrfTypeName(SrfType(iElt)),
     &	  KcElt(iElt),KrElt(iElt),nECoord(iElt)
 2	CONTINUE
	WRITE (*,121) nOutCord

	RETURN
	END SUBROUTINE SumData

C***********************************************************************

	SUBROUTINE ViewData(npts,FU)
	use param_mod
	use src_mod
	use cfiles_mod
	use elt_mod
	use traceutil_mod, only : PolyVtx3DSav,Poly3DApVecL,ifSave3DApVec

	IMPLICIT NONE

	INTEGER iElt,npts,i,j,k,iSeg,FU,iObs

C  input data

#include "viewdata_formats.inc"

C	Screen output:

	IF (FU.EQ.0) THEN

C  input data
	  WRITE(*,100)(ChfRayDir(j),j=1,3),(ChfRayPos(j),j=1,3),zSource,
     &	   IndRef(0),Extinc(0),BaseUnits,WaveUnits,Wavelen,Flux,
     &	   GridTypeName(GridType),Aperture,Obscratn,
     &	   nGridPts,(xGrid(j),j=1,3),(yGrid(j),j=1,3),nElt
	  IF ((GridType.EQ.3).OR.(GridType.EQ.4)) THEN
	    WRITE(*,117)nSeg,width,gap,(SegXgrid(j),j=1,3),
     &                  (SegCoord(i,1),i=1,3)
	    IF (nSeg.GT.1) THEN
	      DO iSeg=2,nSeg
	        WRITE(*,116)(SegCoord(i,iSeg),i=1,3)
	      END DO
	    END IF
	  ELSE IF (GridType.EQ.5) THEN
	    WRITE(*,*)'    gap= ',gap
	    WRITE(*,*)' nRings= ',nRings
	    WRITE(*,*)' nPetals= ',nPetals(0:nRings)
	    WRITE(*,*)' RingRads= ',RingRads(1:nRings)
	    WRITE(*,118)nSeg,radCtr,(SegXgrid(j),j=1,3),
     &                  SegCoord(1:2,1)
	    IF (nSeg.GT.1) THEN
	      DO iSeg=2,nSeg
	        WRITE(*,116)SegCoord(1:2,iSeg)
	      END DO
	    END IF
	  END IF

C  element data

 	  DO 3 iElt=1,nElt
C Basic params & optical params & conic params

	    IF ((SrfType(iElt).EQ.6).OR.(SrfType(iElt).EQ.10)) THEN
	      WRITE(*,98)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	      SrfTypeName(SrfType(iElt)),
     &	      (psiElt(j,iElt),j=1,3),
     &        (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	    ELSE IF (KcElt(iElt).LE.0) THEN
	      WRITE(*,101)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	      SrfTypeName(SrfType(iElt)),fElt(iElt),eElt(iElt),
     &	      KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &        (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	    ELSE
	      WRITE(*,99)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	      SrfTypeName(SrfType(iElt)),
     &	      KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &        (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	    END IF

	    If (EltGrp(0,iElt)/=0) Then
              ! Element group
              WRITE(*,133)
              WRITE(*,*) (EltGrp(i,iElt),i=0,EltGrp(0,iElt))
            End If

	    IF (LGlass(iElt)) THEN
	      WRITE(*,97)GlassElt(iElt),
     &	      (GlassCoef(j,iElt),j=1,6),Extinc(iElt)
	    ELSE
	      WRITE(*,95)IndRef(iElt),Extinc(iElt)
	    ENDIF
 	    IF (SrfType(iElt).EQ.3) THEN
	      WRITE(*,*)'AsphCoef=   ',
     &                (AsphCoef(j,iElt),j=1,maxAsphCoef)
 	    ELSE IF (SrfType(iElt).EQ.6) THEN
	      WRITE(*,104)(AnaCoef(j,iElt),j=1,4)
 	    ELSE IF (SrfType(iElt).EQ.4) THEN
	      WRITE(*,105)(MonCoef(j,iElt),j=1,120)
 	    ELSE IF (SrfType(iElt).EQ.7) THEN
	      WRITE(*,106)UDSrfType(iElt),UDSrfFile(iElt),
     &	      (UDSrfParam(j,iElt),j=1,6)
 	    ELSE IF (SrfType(iElt).EQ.8
     &               .OR. SrfType(iElt).EQ.9) THEN
	      WRITE(*,131) ZernTypeNameL(ZernTypeL(iElt))
	      WRITE(*,107)(ZernCoef(j,iElt),j=1,45)
 	    ELSE IF (SrfType(iElt).EQ.9) THEN
	      WRITE(*,132) GridFile(iElt),nGridMat(iElt),
     &	      GridSrfdx(iElt)
 	    ELSE IF (SrfType(iElt).EQ.10) THEN
	      WRITE(*,104)(AnaCoef(j,iElt),j=1,4)
	    END IF

	    IF ((SrfType(iElt).GE.4).OR.(EltID(iElt).EQ.10)) THEN
	      WRITE(*,108)(pMon(j,iElt),j=1,3),(xMon(j,iElt),j=1,3),
     &	        (yMon(j,iElt),j=1,3),(zMon(j,iElt),j=1,3)
              IF ((SrfType(iElt).EQ.4).OR.(SrfType(iElt).EQ.8)) THEN
		WRITE(*,109)lMon(iElt)
	      END IF
	    END IF

	    ! Added by J.Lou, 11-17-03
	    IF (SrfType(iElt)==9) THEN
	      ! grid surface type
              WRITE(*,108)(pMon(j,iElt),j=1,3),(xMon(j,iElt),j=1,3),
     &          (yMon(j,iElt),j=1,3),(zMon(j,iElt),j=1,3)
            END IF

	    IF (EltID(iElt).EQ.4) THEN
	      WRITE(*,110)(h1HOE(j,iElt),j=1,3),(h2HOE(j,iElt),j=1,3),
     &	      OrderHOE(iElt),WaveHOE(iElt)
	    ELSE IF ((EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13)) THEN
	      WRITE(*,111)(h1HOE(j,iElt),j=1,3),
     &	      OrderHOE(iElt),RuleWidth(iElt)
	    ELSE IF (EltID(iElt).EQ.10) THEN
	      WRITE(*,112)LensArrayType(iElt),LensArrayWidth(iElt)
	    END IF

C  Basic parameter
	    WRITE(*,102)nObs(iElt)

C  obscuration params
	    IF (nObs(iElt).GT.0) THEN
	      DO 21 iObs=1,nObs(iElt)
	        WRITE(*,124)ObsTypeName(ObsType(iObs,iElt))
	        IF (ABS(ObsType(iObs,iElt)).EQ.1) THEN
	          WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,3)
	  	ELSE IF (ABS(ObsType(iObs,iElt)).EQ.4) THEN
	          WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,6)
	        ELSE IF (ABS(ObsType(iObs,iElt)).EQ.8) THEN
                  ! sws polygonal
                  i = ObsVec(3,iObs,iElt)
                  WRITE(*,139) ObsVec(1:2,iObs,iElt),i
                  Do k=1,i
                    WRITE(*,136)(PolyObsVtx(j,k,iObs,iElt),j=1,2)
                  End Do
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.9) THEN
                  WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.10) THEN
                  WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
		ELSE
	          WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,4)
	        END IF
 21	      CONTINUE
	      WRITE(*,126)(xObs(j,iElt),j=1,3)
	    END IF

	    WRITE(*,121)ApTypeName(ApType(iElt))
	    IF (ApType(iElt).NE.0) THEN
	      IF (ApType(iElt).EQ.1) THEN
	        WRITE(*,122)(ApVec(j,iElt),j=1,3)
	      ELSE IF (ApType(iElt).EQ.2) THEN
	        WRITE(*,122)(ApVec(j,iElt),j=1,4)
	      ELSE IF (ApType(iElt).EQ.6) THEN ! hex
                WRITE(*,122)(ApVec(j,iElt),j=1,4)
	      ELSE IF (ApType(iElt).EQ.7) THEN ! polygonal
                i=ApVec(3,iElt)
                If (ifSave3DApVec.and.Poly3DApVecL(iElt)) Then
                  WRITE(*,127) i
                  Do k=1,i
                    WRITE(*,128)(PolyVtx3DSav(j,k,iElt),j=1,3)
                  End Do
                Else
                  ! 2D polygon
                  WRITE(*,129) ApVec(1:2,iElt),i
                  Do k=1,i
                    WRITE(*,136)(PolyApVtx(j,k,iElt),j=1,2)
                  End Do
                End If
	      ELSE
	        WRITE(*,122)(ApVec(j,iElt),j=1,5)
	      END IF
	    END IF

C  basic params
	    WRITE(*,113)zElt(iElt),PropTypeName(PropType(iElt))
C  pinhole
	    IF ((PropType(iElt).EQ.7).OR.(PropType(iElt).EQ.8))
     &	    THEN
	      WRITE(*,151)PinHole(iElt)
	    END IF

C  coatings
	    IF (nCoatElt(iElt).GT.0) THEN
	      WRITE(*,141)nCoatElt(iElt),(CoatIndxElt(j,iElt),
     &	      j=1,nCoatElt(iElt))
	      WRITE(*,142)(CoatThkElt(j,iElt),j=1,nCoatElt(iElt))
	    END IF

C  basic params
	    WRITE(*,114)nECoord(iElt)
	    IF (nECoord(iElt).GT.0) THEN
	      WRITE(*,115)(TElt(1,j,iElt),j=1,nECoord(iElt))
	      DO 2 k=2,6
	        WRITE(*,162)(TElt(k,j,iElt),j=1,nECoord(iElt))
 2	      CONTINUE
	    END IF
 3	  CONTINUE

C  Output Variables
	  WRITE(*,161)nOutCord,((Tout(k,j),j=1,7),k=1,nOutCord)

	ELSE

C	File output:

C  input data

	  If (ifSave3DApVec) WRITE(FU,171)

	  WRITE(FU,100)(ChfRayDir(j),j=1,3),(ChfRayPos(j),j=1,3),
     &	  zSource,IndRef(0),Extinc(0),BaseUnits,WaveUnits,Wavelen,Flux,
     &	  GridTypeName(GridType),Aperture,Obscratn,
     &	  nGridPts,(xGrid(j),j=1,3),(yGrid(j),j=1,3),nElt
	  IF ((GridType.EQ.3).OR.(GridType.EQ.4)) THEN
	    WRITE(FU,117)nSeg,width,gap,(SegXgrid(j),j=1,3),
     &      (SegCoord(i,1),i=1,3)
	    IF (nSeg.GT.1) THEN
	      DO iSeg=2,nSeg
	        WRITE(FU,116)(SegCoord(i,iSeg),i=1,3)
 	      END DO
	    END IF
	  ELSE IF (GridType.EQ.5) THEN
	    WRITE(FU,*)'    gap= ',gap
	    WRITE(FU,*)' nRings= ',nRings
            WRITE(FU,*)' nPetals= ',nPetals(0:nRings)
            WRITE(FU,*)' RingRads= ',RingRads(1:nRings)
	    WRITE(FU,118)nSeg,radCtr,(SegXgrid(j),j=1,3),
     &                   SegCoord(1:2,1)
	    IF (nSeg.GT.1) THEN
	      DO iSeg=2,nSeg
	        WRITE(FU,116)SegCoord(1:2,iSeg)
	      END DO
	    END IF
	  END IF

C  element data

 	  DO 993 iElt=1,nElt
C Basic params & optical params & conic params
	    IF ((SrfType(iElt).EQ.6).OR.(SrfType(iElt).EQ.10)) THEN
	      WRITE(FU,98)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	      SrfTypeName(SrfType(iElt)),
     &	      (psiElt(j,iElt),j=1,3),
     &        (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	    ELSE
	      WRITE(FU,99)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	      SrfTypeName(SrfType(iElt)),
     &	      KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &        (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	    END IF

	    If (EltGrp(0,iElt)/=0) Then
              ! Element group
              WRITE(FU,133)
              WRITE(FU,*) (EltGrp(i,iElt),i=0,EltGrp(0,iElt))
            End If

	    IF (LGlass(iElt)) THEN
	      WRITE(FU,97)GlassElt(iElt),
     &	      (GlassCoef(j,iElt),j=1,6),Extinc(iElt)
	    ELSE
	      WRITE(FU,95)IndRef(iElt),Extinc(iElt)
	    ENDIF
 	    IF (SrfType(iElt).EQ.3) THEN
	      WRITE(FU,*)'AsphCoef=   ',
     &                (AsphCoef(j,iElt),j=1,maxAsphCoef)
 	    ELSE IF (SrfType(iElt).EQ.6) THEN
	      WRITE(FU,104)(AnaCoef(j,iElt),j=1,4)
 	    ELSE IF (SrfType(iElt).EQ.4) THEN
	      WRITE(FU,105)(MonCoef(j,iElt),j=1,120)
 	    ELSE IF (SrfType(iElt).EQ.7) THEN
	      WRITE(FU,106)UDSrfType(iElt),UDSrfFile(iElt),
     &	      (UDSrfParam(j,iElt),j=1,6)
 	    ELSE IF (SrfType(iElt).EQ.8) THEN
	      WRITE(FU,131) ZernTypeNameL(ZernTypeL(iElt))
	      If (ZernTypeL(iElt)==9) Then
                WRITE(FU,*)'ZernAnnualRatio= ',ZernAnnuRatio(iElt)
              End If
	      WRITE(FU,107)(ZernCoef(j,iElt),j=1,45)
 	    ELSE IF (SrfType(iElt).EQ.9) THEN
	      WRITE(FU,132) GridFile(iElt),nGridMat(iElt),
     &	      GridSrfdx(iElt)
 	    ELSE IF (SrfType(iElt).EQ.10) THEN
	      WRITE(FU,104)(AnaCoef(j,iElt),j=1,4)
	    END IF

	    IF ((SrfType(iElt).GE.4).OR.(EltID(iElt).EQ.10)) THEN
	      WRITE(FU,108)(pMon(j,iElt),j=1,3),(xMon(j,iElt),j=1,3),
     &	      (yMon(j,iElt),j=1,3),(zMon(j,iElt),j=1,3)
              IF ((SrfType(iElt).EQ.4).OR.(SrfType(iElt).EQ.8)) THEN
		WRITE(FU,109)lMon(iElt)
	      END IF
	    END IF

	    IF (EltID(iElt).EQ.4) THEN
	      WRITE(FU,110)(h1HOE(j,iElt),j=1,3),(h2HOE(j,iElt),j=1,3),
     &	      OrderHOE(iElt),WaveHOE(iElt)
	    ELSE IF ((EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13)) THEN
	      WRITE(FU,111)(h1HOE(j,iElt),j=1,3),
     &	      OrderHOE(iElt),RuleWidth(iElt)
	    ELSE IF (EltID(iElt).EQ.10) THEN
	      WRITE(FU,112)LensArrayType(iElt),LensArrayWidth(iElt)
	    END IF

C  Basic parameter
	    WRITE(FU,102)nObs(iElt)

C  obscuration params
	    IF (nObs(iElt).GT.0) THEN
	      DO 9921 iObs=1,nObs(iElt)
	        WRITE(FU,124)ObsTypeName(ObsType(iObs,iElt))
	        IF (ABS(ObsType(iObs,iElt)).EQ.1) THEN
	          WRITE(FU,125)(ObsVec(j,iObs,iElt),j=1,3)
	  	ELSE IF (ABS(ObsType(iObs,iElt)).EQ.4) THEN
	          WRITE(FU,125)(ObsVec(j,iObs,iElt),j=1,6)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.7) THEN
                   WRITE(FU,125)(ObsVec(j,iObs,iElt),j=1,5)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.8) THEN
	           ! sws polygonal
                   i = ObsVec(3,iObs,iElt)
                   WRITE(FU,139) ObsVec(1:2,iObs,iElt),i
                   Do k=1,i
                     WRITE(FU,136)(PolyObsVtx(j,k,iObs,iElt),j=1,2)
                   End Do
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.9) THEN
                   WRITE(FU,125)(ObsVec(j,iObs,iElt),j=1,5)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.10) THEN
                   WRITE(FU,125)(ObsVec(j,iObs,iElt),j=1,5)
		ELSE
	          WRITE(FU,125)(ObsVec(j,iObs,iElt),j=1,4)
	        END IF
9921	      CONTINUE
	    END IF

	    IF ((nObs(iElt).GT.0).OR.(ApType(iElt).NE.0)) THEN
	      WRITE(FU,126)(xObs(j,iElt),j=1,3)
	    END IF

	    WRITE(FU,121)ApTypeName(ApType(iElt))
	    IF (ApType(iElt).NE.0) THEN
	      IF (ApType(iElt).EQ.1) THEN
	        WRITE(FU,122)(ApVec(j,iElt),j=1,3)
	      ELSE IF (ApType(iElt).EQ.2) THEN
	        WRITE(FU,122)(ApVec(j,iElt),j=1,4)
	      ELSE IF (ApType(iElt).EQ.6) THEN ! hex
                WRITE(FU,122)(ApVec(j,iElt),j=1,4)
	      ELSE IF (ApType(iElt).EQ.7) THEN ! polygonal
                i=ApVec(3,iElt)
                If (ifSave3DApVec.and.Poly3DApVecL(iElt)) Then
                  WRITE(FU,127) i
                  Do k=1,i
                    WRITE(FU,128)(PolyVtx3DSav(j,k,iElt),j=1,3)
                  End Do
                Else
                  ! 2D polygon
                  WRITE(FU,129) ApVec(1:2,iElt),i
                  Do k=1,i
                    WRITE(FU,136)(PolyApVtx(j,k,iElt),j=1,2)
                  End Do
                End If
	      ELSE
	        WRITE(FU,122)(ApVec(j,iElt),j=1,5)
	      END IF
	    END IF

C  basic params
	    WRITE(FU,113)zElt(iElt),PropTypeName(PropType(iElt))
C  pinhole
	    IF ((PropType(iElt).EQ.7).OR.(PropType(iElt).EQ.8))
     &	    THEN
	      WRITE(FU,151)PinHole(iElt)
	    END IF

C  coatings
	    IF (nCoatElt(iElt).GT.0) THEN
	      WRITE(FU,141)nCoatElt(iElt),(CoatIndxElt(j,iElt),
     &	      j=1,nCoatElt(iElt))
	      WRITE(FU,142)(CoatThkElt(j,iElt),j=1,nCoatElt(iElt))
	    END IF

C  basic params
	    WRITE(FU,114)nECoord(iElt)
	    IF (nECoord(iElt).GT.0) THEN
	      WRITE(FU,115)(TElt(1,j,iElt),j=1,nECoord(iElt))
	      DO 992 k=2,6
	        WRITE(FU,162)(TElt(k,j,iElt),j=1,nECoord(iElt))
992	      CONTINUE
	    END IF
993	  CONTINUE

C  Output Variables
	  WRITE(FU,161)nOutCord,((Tout(k,j),j=1,7),k=1,nOutCord)

	END IF

	RETURN
	END SUBROUTINE ViewData

C***********************************************************************

#if 0
! This is the 'old' version, before Sirlin made his update on the ShowData function.
	SUBROUTINE ShowData(npts,iElt)
	use param_mod
	use src_mod
	use cfiles_mod
	use elt_mod

	IMPLICIT NONE

	INTEGER iElt,npts,i,j,k,iSeg,iObs

#include "viewdata_formats.inc"

c  input data
 100	FORMAT(1P,'ChfRayDir=',3D23.15/'ChfRayPos=',3D23.15/
     &	' zSource= ',D23.15/'  IndRef= ',D23.15/'  Extinc= ',D23.15/
     &	'BaseUnits= ',A16/'WaveUnits= ',A16/' Wavelen= ',D23.15/
     &	'    Flux= ',D23.15/'GridType=   ',A16/'Aperture= ',D23.15/
     &	'Obscratn= ',D23.15/'nGridpts= ',i4/'   xGrid= ',3D23.15/
     &	'   yGrid= ',3D23.15/'    nElt= ',i4)
 116	FORMAT(1P,10x,i5,5x,i5,5x,i5)
 117	FORMAT(1P,'    nSeg= ',i4/'   width= ',D23.15/
     &	'     gap= ',D23.15/'SegXgrid= ',D23.15,2(2x,D23.15)/
     &  'SegCoord= ',i5,2(5x,i5))
 118	FORMAT(1P,'    nSeg= ',i4/
     &	'  radCtr= ',D23.15/'SegXgrid= ',D23.15,2(2x,D23.15)/
     &  'SegCoord= ',i5,5x,i5)

C Basic params & optical params & conic params
  98	FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8/
     &	'  psiElt= ',3D23.15/
     &	'  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
  99	FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8/
     &	'   KrElt= ',D23.15/'   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &	'  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
 101	FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8/
     &	'    fElt= ',D23.15/'    eElt= ',D23.15/
     &	'   KrElt= ',D23.15/'   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &	'  VptElt= ',3D23.15/'  RptElt= ',3D23.15)

 130    FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8)
 133    FORMAT(1P,'    psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
 134    FORMAT(1P,'    fElt= ',D23.15/
     &  '    eElt= ',D23.15/
     &  '   KrElt= ',D23.15/'   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
 135    FORMAT(1P,'    KrElt= ',D23.15/
     &  '   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)

  95	FORMAT(1P,'  IndRef= ',D23.15/'  Extinc= ',D23.15)
  96	FORMAT(1P,'GlassElt= ',A16/'  Extinc= ',D23.15)
  97	FORMAT(1P,'GlassElt= ',A16/'GlassCoef= ',6D23.15/
     &	'  Extinc= ',D23.15)
 102	FORMAT(1P,'    nObs= ',i4)
 121	FORMAT(1P,'  ApType=   ',A16)

C  asphere params
 103	FORMAT(1P,'AsphCoef= ',4D23.15)

C  anamorph params
 104	FORMAT(1P,' AnaCoef= ',4D23.15)

C  monomial params
 105	FORMAT(1P,' MonCoef= ',6D23.15/	19(10x,6D23.15/))

C  user-defined params
 106	FORMAT(1P,'UDSrfType= ',i4/'UDSrfFile= ',A/
     &	'UDSrfParam= ',6D23.15)

C  zernike params
 107	FORMAT(1P,'ZernCoef= ',6D23.15/ 6(10x,6D23.15/) (10x,3D23.15))

C  zernike type
 131	FORMAT(1P,'ZernType=   ',A16)

C  grid surface params
 132	FORMAT(1P,'GridFile= ',A/'nGridMat= ',i4/
     &	'GridSrfdx= ',D23.15)

C  surface coord params
 108	FORMAT(1P,'    pMon= ',3D23.15/'    xMon= ',3D23.15/
     &	'    yMon= ',3D23.15/'    zMon= ',3D23.15)
 109	FORMAT(1P,'    lMon= ',D23.15)

C  HOE params & Grating params
 110	FORMAT(1P,'   h1HOE= ',3D23.15/'   h2HOE= ',3D23.15/
     &	'OrderHOE= ',D23.15/' WaveHOE= ',D23.15)
 111	FORMAT(1P,'   h1HOE= ',3D23.15/
     &	'OrderHOE= ',D23.15/'RuleWidth= ',D23.15)

C  lens array params
 112	FORMAT(1P,'LensArrayType= ',i2/'LensArrayWidth= ',D23.15)

C  obscuration params
 122	FORMAT(1P,'   ApVec= ',6D23.15)
 124	FORMAT(1P,' ObsType=   ',A16)
 125	FORMAT(1P,'  ObsVec= ',6D23.15)
 126	FORMAT(1P,'    xObs= ',3D23.15)

C  coatings
 141	FORMAT(1P,'   nCoat= ',i4/'CoatIndx= ',5D23.15/3(10X,5D23.15))
 142	FORMAT(1P,'CoatThk= ',5D23.15/3(10X,5D23.15))

C  pinhole
 151	FORMAT(1P,' PinHole= ',D23.15)

C  basic params
 113	FORMAT(1P,'    zElt= ',D23.15/'PropType=   ',A16)
 114	FORMAT(1P,' nECoord= ',i4)
 115	FORMAT(1P,'    TElt= ',6D23.15)
 162	FORMAT(1P,10x,6D23.15)

C  Output Variables
 161	FORMAT(1P/'nOutCord= ',i4/'    Tout= ',7D23.15,6(/10x,7D23.15))

C	Screen output:

	IF (iElt.EQ.0) THEN

C  input data
	  WRITE(*,100)(ChfRayDir(j),j=1,3),(ChfRayPos(j),j=1,3),zSource,
     &	  IndRef(0),Extinc(0),BaseUnits,WaveUnits,Wavelen,Flux,
     &	  GridTypeName(GridType),Aperture,Obscratn,
     &	  nGridPts,(xGrid(j),j=1,3),(yGrid(j),j=1,3),nElt
	  IF ((GridType.EQ.3).OR.(GridType.EQ.4)) THEN
	    WRITE(*,117)nSeg,width,gap,(SegXgrid(j),j=1,3),
     &      (SegCoord(i,1),i=1,3)
	    IF (nSeg.GT.1) THEN
	      DO iSeg=2,nSeg
	        WRITE(*,116)(SegCoord(i,iSeg),i=1,3)
	      END DO
	    END IF
	  ELSE IF (GridType.EQ.5) THEN
	    WRITE(*,118)nSeg,nPetals,radCtr,(SegXgrid(j),j=1,3),
     &      SegCoord(1,1)
	    IF (nSeg.GT.1) THEN
	      DO iSeg=2,nSeg
	        WRITE(*,116)SegCoord(1,iSeg)
	      END DO
	    END IF
	  END IF
	ELSE
C  Element Data
C Basic params & optical params & conic params
	    IF (EltGrp(0,iElt)==0) THEN
	      IF ((SrfType(iElt).EQ.6).OR.(SrfType(iElt).EQ.10)) THEN
	        WRITE(*,98)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	          SrfTypeName(SrfType(iElt)),
     &	          (psiElt(j,iElt),j=1,3),
     &            (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	      ELSE IF (KcElt(iElt).LE.0) THEN
	        WRITE(*,101)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	        SrfTypeName(SrfType(iElt)),fElt(iElt),eElt(iElt),
     &	        KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &          (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	      ELSE
	        WRITE(*,99)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &	        SrfTypeName(SrfType(iElt)),
     &	        KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &          (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	      END IF
	    ELSE
	      ! With EltGrp defined at iElt
	      IF ((SrfType(iElt).EQ.6).OR.(SrfType(iElt).EQ.10)) THEN
                WRITE(*,130)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &                      SrfTypeName(SrfType(iElt))
                WRITE(*,*)'   EltGrp= ', EltGrp(1:EltGrp(0,iElt),iElt)
                WRITE(*,133)(psiElt(j,iElt),j=1,3),
     &                (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              ELSE IF (KcElt(iElt).LE.0) THEN
	        WRITE(*,130)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &            	    SrfTypeName(SrfType(iElt))
                WRITE(*,*)'   EltGrp= ', EltGrp(1:EltGrp(0,iElt),iElt)
	        WRITE(*,134)fElt(iElt),eElt(iElt),
     &       	    KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &              (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              ELSE
	        WRITE(*,130)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &                      SrfTypeName(SrfType(iElt))
                WRITE(*,*)'   EltGrp= ', EltGrp(1:EltGrp(0,iElt),iElt)
	        WRITE(*,135)
     &            KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &            (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
	      END IF
	    END IF

	    IF (LGlass(iElt)) THEN
	      WRITE(*,97)GlassElt(iElt),
     &	      (GlassCoef(j,iElt),j=1,6),Extinc(iElt)
	    ELSE
	      WRITE(*,95)IndRef(iElt),Extinc(iElt)
	    ENDIF

 	  IF (SrfType(iElt).EQ.3) THEN
	    !WRITE(*,103)(AsphCoef(j,iElt),j=1,maxAsphCoef)
	    WRITE(*,*)'AsphCoef=   ',
     &                (AsphCoef(j,iElt),j=1,maxAsphCoef)
 	  ELSE IF (SrfType(iElt).EQ.6) THEN
	    WRITE(*,104)(AnaCoef(j,iElt),j=1,4)
 	  ELSE IF (SrfType(iElt).EQ.4) THEN
	    WRITE(*,105)(MonCoef(j,iElt),j=1,120)
 	  ELSE IF (SrfType(iElt).EQ.7) THEN
	    WRITE(*,106)UDSrfType(iElt),UDSrfFile(iElt),
     &	    (UDSrfParam(j,iElt),j=1,6)
 	  ELSE IF (SrfType(iElt).EQ.8) THEN
	    WRITE(*,131) ZernTypeNameL(ZernTypeL(iElt))
	    If (ZernTypeL(iElt)==9) Then
              WRITE(*,*)'ZernAnnualRatio= ',ZernAnnuRatio(iElt)
            End If
	    WRITE(*,107)(ZernCoef(j,iElt),j=1,45)
 	  ELSE IF (SrfType(iElt).EQ.9) THEN
	    WRITE(*,132) GridFile(iElt),nGridMat(iElt),
     &	    GridSrfdx(iElt)
 	  ELSE IF (SrfType(iElt).EQ.10) THEN
	    WRITE(*,104)(AnaCoef(j,iElt),j=1,4)
	  END IF

	  IF ((SrfType(iElt).GE.4).OR.(EltID(iElt).EQ.10)) THEN
	    WRITE(*,108)(pMon(j,iElt),j=1,3),(xMon(j,iElt),j=1,3),
     &	    (yMon(j,iElt),j=1,3),(zMon(j,iElt),j=1,3)
            IF ((SrfType(iElt).EQ.4).OR.(SrfType(iElt).EQ.8)) THEN
	      WRITE(*,109)lMon(iElt)
	    END IF
	  END IF

	  IF (EltID(iElt).EQ.4) THEN
	    WRITE(*,110)(h1HOE(j,iElt),j=1,3),(h2HOE(j,iElt),j=1,3),
     &	    OrderHOE(iElt),WaveHOE(iElt)
	  ELSE IF ((EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13)) THEN
	    WRITE(*,111)(h1HOE(j,iElt),j=1,3),
     &	    OrderHOE(iElt),RuleWidth(iElt)
	  ELSE IF (EltID(iElt).EQ.10) THEN
	    WRITE(*,112)LensArrayType(iElt),LensArrayWidth(iElt)
	  END IF

C  Basic parameter
	  WRITE(*,102)nObs(iElt)

C  obscuration params
	  IF (nObs(iElt).GT.0) THEN
	    DO 21 iObs=1,nObs(iElt)
	      WRITE(*,124)ObsTypeName(ObsType(iObs,iElt))
	      IF (ABS(ObsType(iObs,iElt)).EQ.1) THEN
	        WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,3)
	      ELSE IF (ABS(ObsType(iObs,iElt)).EQ.4) THEN
	        WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,6)
              ELSE IF (ABS(ObsType(iObs,iElt)).EQ.7) THEN
                WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
              ELSE IF (ABS(ObsType(iObs,iElt)).EQ.8) THEN
                ! sws polygonal
                i = ObsVec(3,iObs,iElt)
                WRITE(*,139) ObsVec(1:2,iObs,iElt),i
                Do k=1,i
                  WRITE(*,136)(PolyObsVtx(j,k,iObs,iElt),j=1,2)
                End Do
              ELSE IF (ABS(ObsType(iObs,iElt)).EQ.9) THEN
                WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
              ELSE IF (ABS(ObsType(iObs,iElt)).EQ.10) THEN
                WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
	      ELSE
	        WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,4)
	      END IF
 21	    CONTINUE
	    WRITE(*,126)(xObs(j,iElt),j=1,3)
	  END IF

	  WRITE(*,121)ApTypeName(ApType(iElt))
	  IF (ApType(iElt).NE.0) THEN
	    IF (ApType(iElt).EQ.1) THEN
	      WRITE(*,122)(ApVec(j,iElt),j=1,3)
	    ELSE IF (ApType(iElt).EQ.2) THEN
	      WRITE(*,122)(ApVec(j,iElt),j=1,4)
	    ELSE IF (ApType(iElt).EQ.6) THEN ! hex
              WRITE(*,122)(ApVec(j,iElt),j=1,4)
	    ELSE
	      WRITE(*,122)(ApVec(j,iElt),j=1,5)
	    END IF
	  END IF

C  basic params
	  WRITE(*,113)zElt(iElt),PropTypeName(PropType(iElt))
C  pinhole
	  IF ((PropType(iElt).EQ.7).OR.(PropType(iElt).EQ.8))
     &	  THEN
	    WRITE(*,151)PinHole(iElt)
	  END IF

C  coatings
	  IF (nCoatElt(iElt).GT.0) THEN
	    WRITE(*,141)nCoatElt(iElt),(CoatIndxElt(j,iElt),
     &	    j=1,nCoatElt(iElt))
	    WRITE(*,142)(CoatThkElt(j,iElt),j=1,nCoatElt(iElt))
	  END IF

C  basic params
	  WRITE(*,114)nECoord(iElt)
	  IF (nECoord(iElt).GT.0) THEN
	    WRITE(*,115)(TElt(1,j,iElt),j=1,nECoord(iElt))
	    DO 2 k=2,6
	      WRITE(*,162)(TElt(k,j,iElt),j=1,nECoord(iElt))
 2	    CONTINUE
	  END IF

	  IF (iElt.EQ.nElt) THEN
C  Output Variables
	    WRITE(*,161)nOutCord,((Tout(k,j),j=1,7),k=1,nOutCord)
	  END IF
	END IF

	RETURN
	END SUBROUTINE ShowData
#endif


C***********************************************************************
!  This is Sirlin's updated version ShowData routine
!  06/25/2009
        SUBROUTINE ShowData(npts,iElt)
        use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
        use traceutil_mod, only : PolyVtx3DSav,Poly3DApVecL,
     &                            ifSave3DApVec

        IMPLICIT NONE

        INTEGER iElt,npts,i,j,k,iSeg,iObs

#include "viewdata_formats.inc"


C       Screen output:

        IF (iElt.EQ.0) THEN

C  input data
          WRITE(*,100)(ChfRayDir(j),j=1,3),(ChfRayPos(j),j=1,3),zSource,
     &     IndRef(0),Extinc(0),BaseUnits,WaveUnits,Wavelen,Flux,
     &     GridTypeName(GridType),Aperture,Obscratn,
     &     nGridPts,(xGrid(j),j=1,3),(yGrid(j),j=1,3),nElt
           IF ((GridType.EQ.3).OR.(GridType.EQ.4)) THEN
             WRITE(*,117)nSeg,width,gap,(SegXgrid(j),j=1,3),
     &                   (SegCoord(i,1),i=1,3)
            IF (nSeg.GT.1) THEN
              DO iSeg=2,nSeg
                WRITE(*,116)(SegCoord(i,iSeg),i=1,3)
              END DO
            END IF
          ELSE IF (GridType.EQ.5) THEN
	    WRITE(*,*)'    gap= ',gap
	    WRITE(*,*)' nRings= ',nRings
            WRITE(*,*)' nPetals= ',nPetals(0:nRings)
            WRITE(*,*)' RingRads= ',RingRads(1:nRings)
            WRITE(*,118)nSeg,radCtr,(SegXgrid(j),j=1,3),
     &                  SegCoord(1:2,1)
            IF (nSeg.GT.1) THEN
              DO iSeg=2,nSeg
                WRITE(*,116)SegCoord(1:2,iSeg)
              END DO
            END IF
          END IF
        ELSE

C  Element Data
C Basic params & optical params & conic params
            IF (EltGrp(0,iElt)==0) THEN
              IF ((SrfType(iElt).EQ.6).OR.(SrfType(iElt).EQ.10)) THEN
                WRITE(*,98)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &            SrfTypeName(SrfType(iElt)),
     &            (psiElt(j,iElt),j=1,3),
     &            (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              ELSE IF (KcElt(iElt).LE.0) THEN
                WRITE(*,101)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &          SrfTypeName(SrfType(iElt)),fElt(iElt),eElt(iElt),
     &          KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &          (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              ELSE
                WRITE(*,99)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &          SrfTypeName(SrfType(iElt)),
     &          KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &          (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              END IF
            ELSE
              ! With EltGrp defined at iElt
              IF ((SrfType(iElt).EQ.6).OR.(SrfType(iElt).EQ.10)) THEN
                WRITE(*,230)iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &                      SrfTypeName(SrfType(iElt))
                WRITE(*,*)'   EltGrp= ', EltGrp(1:EltGrp(0,iElt),iElt)
                WRITE(*,233)(psiElt(j,iElt),j=1,3),
     &            (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              ELSE IF (KcElt(iElt).LE.0) THEN
                WRITE(*,230)iElt,EltName(iElt),
     &                      EltTypeName(EltID(iElt)),
     &                      SrfTypeName(SrfType(iElt))
                WRITE(*,*)'   EltGrp= ', EltGrp(1:EltGrp(0,iElt),iElt)
                WRITE(*,234)fElt(iElt),eElt(iElt),
     &            KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &            (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              ELSE
                WRITE(*,230)
     &            iElt,EltName(iElt),EltTypeName(EltID(iElt)),
     &            SrfTypeName(SrfType(iElt))
                WRITE(*,*)'   EltGrp= ',
     &                    EltGrp(1:EltGrp(0,iElt),iElt)
                WRITE(*,235)
     &            KrElt(iElt),KcElt(iElt),(psiElt(j,iElt),j=1,3),
     &            (VptElt(j,iElt),j=1,3),(RptElt(j,iElt),j=1,3)
              END IF
            END IF

            IF (LGlass(iElt)) THEN
              WRITE(*,97)GlassElt(iElt),
     &        (GlassCoef(j,iElt),j=1,6),Extinc(iElt)
            ELSE
              WRITE(*,95)IndRef(iElt),Extinc(iElt)
            ENDIF

          IF (SrfType(iElt).EQ.3) THEN
            !WRITE(*,103)(AsphCoef(j,iElt),j=1,maxAsphCoef)
            WRITE(*,*)'AsphCoef=   ',
     &                (AsphCoef(j,iElt),j=1,maxAsphCoef)
          ELSE IF (SrfType(iElt).EQ.6) THEN
            WRITE(*,104)(AnaCoef(j,iElt),j=1,4)
          ELSE IF (SrfType(iElt).EQ.4) THEN
            WRITE(*,105)(MonCoef(j,iElt),j=1,120)
          ELSE IF (SrfType(iElt).EQ.7) THEN
            WRITE(*,106)UDSrfType(iElt),UDSrfFile(iElt),
     &      (UDSrfParam(j,iElt),j=1,6)
          ELSE IF (SrfType(iElt).EQ.8) THEN
            WRITE(*,131) ZernTypeNameL(ZernTypeL(iElt))
            If (ZernTypeL(iElt)==9) Then
              WRITE(*,*)'ZernAnnualRatio= ',ZernAnnuRatio(iElt)
            End If
            WRITE(*,107)(ZernCoef(j,iElt),j=1,45)
          ELSE IF (SrfType(iElt).EQ.9) THEN
            WRITE(*,132) GridFile(iElt),nGridMat(iElt),
     &      GridSrfdx(iElt)
          ELSE IF (SrfType(iElt).EQ.10) THEN
            WRITE(*,104)(AnaCoef(j,iElt),j=1,4)
          END IF

          IF ((SrfType(iElt).GE.4).OR.(EltID(iElt).EQ.10)) THEN
            WRITE(*,108)(pMon(j,iElt),j=1,3),(xMon(j,iElt),j=1,3),
     &      (yMon(j,iElt),j=1,3),(zMon(j,iElt),j=1,3)
            IF ((SrfType(iElt).EQ.4).OR.(SrfType(iElt).EQ.8)) THEN
              WRITE(*,109)lMon(iElt)
            END IF
          END IF

          IF (EltID(iElt).EQ.4) THEN
            WRITE(*,110)(h1HOE(j,iElt),j=1,3),(h2HOE(j,iElt),j=1,3),
     &      OrderHOE(iElt),WaveHOE(iElt)
          ELSE IF ((EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13)) THEN
            WRITE(*,111)(h1HOE(j,iElt),j=1,3),
     &      OrderHOE(iElt),RuleWidth(iElt)
          ELSE IF (EltID(iElt).EQ.10) THEN
            WRITE(*,112)LensArrayType(iElt),LensArrayWidth(iElt)
          END IF

C  Basic parameter
          WRITE(*,102)nObs(iElt)

C  obscuration params
          IF (nObs(iElt).GT.0) THEN
             DO 21 iObs=1,nObs(iElt)
                WRITE(*,124)ObsTypeName(ObsType(iObs,iElt))
                IF (ABS(ObsType(iObs,iElt)).EQ.1) THEN
                   WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,3)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.4) THEN
                   WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,6)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.7) THEN
                   WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.8) THEN
! sws polygonal
                   i = ObsVec(3,iObs,iElt)
                   WRITE(*,139) ObsVec(1:2,iObs,iElt),i
                   Do k=1,i
                      WRITE(*,136)(PolyObsVtx(j,k,iObs,iElt),j=1,2)
                   End Do
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.9) THEN
                   WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
                ELSE IF (ABS(ObsType(iObs,iElt)).EQ.10) THEN
                   WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,5)
                ELSE
                   WRITE(*,125)(ObsVec(j,iObs,iElt),j=1,4)
                END IF
 21          CONTINUE
             WRITE(*,126)(xObs(j,iElt),j=1,3)
          END IF

          WRITE(*,121)ApTypeName(ApType(iElt))
          IF (ApType(iElt).NE.0) THEN
             IF (ApType(iElt).EQ.1) THEN
                WRITE(*,122)(ApVec(j,iElt),j=1,3)
             ELSE IF (ApType(iElt).EQ.2) THEN
                WRITE(*,122)(ApVec(j,iElt),j=1,4)
             ELSE IF (ApType(iElt).EQ.6) THEN ! hex
                WRITE(*,122)(ApVec(j,iElt),j=1,4)
             ELSE IF (ApType(iElt).EQ.7) THEN ! polygonal
                i=ApVec(3,iElt)
                If (ifSave3DApVec.and.Poly3DApVecL(iElt)) Then
                   WRITE(*,127) i
                   Do k=1,i
                      WRITE(*,128)(PolyVtx3DSav(j,k,iElt),j=1,3)
                   End Do
                Else
! 2D polygon
                   WRITE(*,129) ApVec(1:2,iElt),i
                   Do k=1,i
                      WRITE(*,136)(PolyApVtx(j,k,iElt),j=1,2)
                   End Do
                End If
             ELSE
                WRITE(*,122)(ApVec(j,iElt),j=1,5)
             END IF
          END IF

C  basic params
          WRITE(*,113)zElt(iElt),PropTypeName(PropType(iElt))
C  pinhole
          IF ((PropType(iElt).EQ.7).OR.(PropType(iElt).EQ.8))
     &    THEN
            WRITE(*,151)PinHole(iElt)
          END IF

C  coatings
          IF (nCoatElt(iElt).GT.0) THEN
            WRITE(*,141)nCoatElt(iElt),(CoatIndxElt(j,iElt),
     &      j=1,nCoatElt(iElt))
            WRITE(*,142)(CoatThkElt(j,iElt),j=1,nCoatElt(iElt))
          END IF

C  basic params
          WRITE(*,114)nECoord(iElt)
          IF (nECoord(iElt).GT.0) THEN
            WRITE(*,115)(TElt(1,j,iElt),j=1,nECoord(iElt))
            DO 2 k=2,6
              WRITE(*,162)(TElt(k,j,iElt),j=1,nECoord(iElt))
 2          CONTINUE
          END IF

          IF (iElt.EQ.nElt) THEN
C  Output Variables
            WRITE(*,161)nOutCord,((Tout(k,j),j=1,7),k=1,nOutCord)
          END IF
        END IF

        RETURN
        END SUBROUTINE ShowData


C***********************************************************************
C       - Save a new MACOS prescription  - jzlou
C***********************************************************************

	SUBROUTINE CSAVE(npts)
	use param_mod
	use src_mod
	use cfiles_mod
	use elt_mod

	IMPLICIT NONE

	INTEGER l,npts,ICLEN
	CHARACTER(len=MacosCharLen) :: ANS
	!CHARACTER*132 ANS
	LOGICAL exist

	nGridPts=npts+1

C	Get file name

 43	CONTINUE
	CALL CACCEPT(filnam,' ','Enter new input file name: ')
	l=ICLEN(filnam)
	infil=filnam
	outfil=filnam
	cmatfil=filnam
	plotfil=filnam
	mfil=filnam
	pertfil=filnam
	nasfil=filnam
	hdrfil=filnam
	imgfil=filnam
	infil(l+1:l+3)='.in'
	outfil(l+1:l+4)='.out'
	cmatfil(l+1:l+5)='.cmat'
	plotfil(l+1:l+5)='.plot'
	mfil(l+1:l+2)='.m'
	pertfil(l+1:l+5)='.pert'
	nasfil(l+1:l+4)='.nas'
	imgfil(l+1:l+4)='.img'
	hdrfil(l+1:l+4)='.hdr'
	INQUIRE (file=infil,exist=exist)
	IF (exist) THEN
	  L=ICLEN(infil)
	  CALL CACCEPT(ANS,'YES',
     &	  'Input file '//infil(1:L)//' exists.  Replace?')
	  IF ((ans(1:1).EQ.'y').OR.(ans(1:1).EQ.'Y')) THEN
	    OPEN (UNIT=1,FILE=infil(1:L),STATUS='OLD')
	    CLOSE (UNIT=1,STATUS='DELETE')
	  ELSE
	    GO TO 43
	  END IF
	END IF

C	Create new input file and write input data

	WRITE(*,*) ' File '//INFIL(1:ICLEN(infil))//' being created.'
	OPEN (UNIT=1,FILE=INFIL,STATUS='NEW',FORM='FORMATTED')
	CALL ViewData(npts,1)
	CLOSE (UNIT=1)

	RETURN
	END SUBROUTINE CSAVE

C***********************************************************************

	SUBROUTINE SetIntSrf(ifIntSrfSet)

C       A routine to read interpolated surface data.

C	A binary surface data file is expected to have the name
C	'filnam.srfN.bin', where filnam is the .in-file root name and N
C	is the element number. A text surface data file is expected to
C	have the name 'filnam.srfN.txt'

C	If no file with the expected name is found, the user is
C	prompted for a file name.

C	The datafile has numbers nDP and iComputeDZ on the first line.
C	nDP gives the number of data points. iComputeDZ=1 indicates that
C	slopes are to be computed by the DIDES code; iComputeDZ=0
C	indicates that slopes are provided in the data file. Each
C	subsequent line has 3 or 5 real numbers, giving the x and y
C	coordinates of the data point, the surface height at that point,
C	and (if iComputeDZ=0) the x and y slopes at that point.

	use param_mod
	use cfiles_mod
	use elt_mod

	IMPLICIT NONE

	CHARACTER*8 CInteger
	CHARACTER(len=MacosCharLen) :: newFilNam
	!CHARACTER*132 newFilNam
	LOGICAL exist,ifIntSrfSet(mElt),ifBin,LCMP
	INTEGER i,j,iElt,iSrf,iIntSrf,L,nDatPts,iDP,jDP,ICLEN,IERROR,
     &	        IDNDP(31),iComputeDZ,nRead,IMODE
	REAL*8 xIDES,yIDES,wVec(6,3,30),dDatPts,SURFACE(mDP,5)

C  Formats

 100	FORMAT(2x,i10)
 101	FORMAT(2x,5D23.15)
 102	FORMAT(' Interpolated surface data not correctly set for ',
     &	'element',i4)

C  Loop through each element to find interpolated surfaces

	iIntSrf=0
	DO 18 iElt=1,nElt
	  IF (SrfType(iElt).NE.5) THEN
	    iEltToIntSrf(iElt)=0
	    ifIntSrfSet(iElt)=.FALSE.
	  ELSE

C  Obtain file name for surface data

	    CALL IntToChar(Cinteger,iElt,i)
	    newFilNam=filnam
	    j=ICLEN(newFilNam)
	    SrfFil=newFilNam(1:j)//'.srf'//Cinteger(1:i)//'.bin'
  1	    CONTINUE
	      INQUIRE (FILE=SrfFil,EXIST=exist)
	      IF (.NOT.exist) THEN
	        L=ICLEN(SrfFil)
	        SrfFil=SrfFil(1:(L-4))//'.txt'
	        INQUIRE (FILE=SrfFil,EXIST=exist)
	        IF (.NOT.exist) THEN
	          WRITE(*,*)'File '//SrfFil(1:L)//' does not exist.'
	          CALL CACCEPT(newFilNam,'QUIT',
     &	          'Enter new .srf-file root name or QUIT to exit:')
	          IF (LCMP(newFilNam,'QUIT',4)) GO TO 2
	          j=ICLEN(newFilNam)
	          SrfFil=newFilNam(1:j)//'.srf'//Cinteger(1:i)//'.bin'
	          GO TO 1
	        ELSE
	          ifBin=.FALSE.
	        END IF
	      ELSE
	        ifBin=.TRUE.
	      END IF

C  Read surface datafile

	    IF (ifBin) THEN
	      OPEN (43,FILE=SrfFil,STATUS='OLD',FORM='UNFORMATTED')
	      READ(43,ERR=2,END=2)nDatPts,iComputeDZ
	      iSrf=iIntSrf+1
	      IF (iComputeDZ.GT.0) THEN
	        ifComputeDZ(iSrf)=.TRUE.
	        nRead=3
	      ELSE
	        ifComputeDZ(iSrf)=.FALSE.
	        nRead=5
	      END IF
	      IF (nDatPts.GT.mDP) THEN
	        WRITE(*,*)' Too many data points.'
	        GO TO 2
	      END IF
	      IF (nDatPts.LT.3) THEN
	        WRITE(*,*)' Too few data points.'
	        GO TO 2
	      END IF
	      jDP=0
	      DO iDP=1,nDatPts
	        jDP=jDP+1
	        READ(43,ERR=2,END=4)(SURFACE(iDP,i),i=1,nRead)
	      END DO
	    ELSE
	      OPEN (43,FILE=SrfFil,STATUS='OLD',FORM='FORMATTED')
	      READ(43,*,ERR=2,END=2)nDatPts,iComputeDZ
	      iSrf=iIntSrf+1
	      IF (iComputeDZ.GT.0) THEN
	        ifComputeDZ(iSrf)=.TRUE.
	        nRead=3
	      ELSE
	        ifComputeDZ(iSrf)=.FALSE.
	        nRead=5
	      END IF
	      IF (nDatPts.GT.mDP) THEN
	        WRITE(*,*)' Too many data points.'
	        GO TO 2
	      END IF
	      IF (nDatPts.LT.3) THEN
	        WRITE(*,*)' Too few data points.'
	        GO TO 2
	      END IF
	      jDP=0
	      DO iDP=1,nDatPts
	        jDP=jDP+1
	        READ(43,*,ERR=2,END=4)(SURFACE(iDP,i),i=1,nRead)
  	      END DO
	    END IF

C  Here if data successfully read

  4	    CONTINUE
	    CLOSE(43)
	    iIntSrf=iIntSrf+1
	    nDP(iIntSrf)=jDP
	    iEltToIntSrf(iElt)=iIntSrf

	    DO i=1,3
	      DO iDP=1,jDP
	        XYZSRF(iDP,i,iIntSrf)=SURFACE(iDP,i)
	      END DO
	    END DO

	    IF (ifComputeDZ(iIntSrf)) THEN
	      IMODE=2
	    ELSE
	      ! surface slopes are provided as input -jzlou
	      IMODE=1
	      DO i=1,2
	        j=i+3
	        DO iDP=1,jDP
	          DZSRF(i,iDP,iIntSrf)=SURFACE(iDP,j)
	        END DO
	      END DO
	    END IF

	    CALL DIDES(IMODE,NDP(iIntSrf),mDP,xIDES,yIDES,wVec,
     &	    	   XYZSRF(1,1,iIntSrf),DZSRF(1,1,iIntSrf),
     &	    	   IWORK(1,iIntSrf),DWORK(1,iIntSrf),IERROR)

	    IF (IERROR.EQ.0) THEN
	      ifIntSrfSet(iElt)=.TRUE.
	      GO TO 18
	    END IF

C  Here if data not correctly set

  2	    CONTINUE
	    CLOSE(43)
	    WRITE(*,102)iElt
	    iEltToIntSrf(iElt)=0
	    ifIntSrfSet(iElt)=.FALSE.
	  END IF
 18	CONTINUE

	RETURN
	END SUBROUTINE SetIntSrf

C***********************************************************************

	SUBROUTINE GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG)

C     THIS ROUTINE TAKES THE STRING, AND DOES SOME PARSING ON IT
C     AND RETURNS THE VARIABLE NAME, THE INDEX LOCATIONS IF THE
C     VARIABLE IS AN ARRAY, AND THE VALUE OF THE VARIABLE NAME.

	IMPLICIT NONE

        CHARACTER(len=MacosCharLen) :: STRING  ! to match those in macosio.F, -jzlou
!	CHARACTER*130 STRING  ! original
!	CHARACTER*180 STRING  ! -jzlou,  to match those in macosio.F
        CHARACTER(len=MacosVarNamLen) :: VAR_NAM
!	CHARACTER*120 VALUE   ! original
 	CHARACTER(len=MacosValLen) :: VALUE   ! -jzlou,  to match those in macosio.F
	Integer, parameter :: maxStringLen=MacosCharLen,
     &                        maxValueLen=MacosValLen

	CHARACTER*4 IND(3)
	INTEGER :: I,J,K,NVTOK,ntok
        LOGICAL :: IND_FND,CLS_IND,ARR_FLG,NO_LET,VAL_LET,SP_FLG
	LOGICAL :: LInCommentMode,LCMP,in_tok

C     i=counter for string, j=counter for output variable being filled
C     k=counter for # of indexes


  123	FORMAT(A)

	! Test only
	!STRING='  MACOS  =an optical % analysis program'

	SP_FLG=.FALSE.  ! set to false on entrance; will be set to true later i
			! when a blank or comment line is found.

	VAR_NAM=' '
	VALUE=' '
	IND(1)='    '; IND(2)='    '; IND(3)='    '
	LInCommentMode=.FALSE.

   20   CONTINUE

	IF (SP_FLG) THEN
	  READ(1,123,END=354) STRING
	  !print*,'GET_EQ() SECOND STRING = ',STRING
 354	  CONTINUE
	END IF

	J=1; K=1; ntok=0
	IND_FND=.FALSE.
	CLS_IND=.FALSE.
	ARR_FLG=.FALSE.
	NO_LET=.TRUE.
	VAL_LET=.TRUE.

C       BEGIN PARSING THE STRING
	in_tok=.false.
	DO 12 I=1,maxStringLen
	  if (J .EQ. maxValueLen) exit ! reached max size of VALUE, exit loop !-jzlou

	  IF (STRING(I:I) .EQ. ' ') THEN
	    IF (NO_LET) THEN
	      CONTINUE
	    ELSE IF (.NOT.IND_FND) THEN
	      VAR_NAM(J:J)=STRING(I:I)
	      !print*,'VAR_NAM(J:J)=',VAR_NAM(J:J)
	      J=J+1
	      ! end of a previous token, increment token count
	      if (in_tok) ntok=ntok+1
	    ELSE IF ( (IND_FND) .AND. (CLS_IND) ) THEN
	      IF (VAL_LET) THEN
	        CONTINUE
	      ELSE
	        VALUE(J:J)=STRING(I:I)
	        J=J+1
	      END IF
	      if (in_tok) ntok=ntok+1
	    END IF
	    in_tok=.false.
	  ELSE IF (STRING(I:I) .EQ. '[') THEN
	    in_tok=.false.
	    CONTINUE
	  ELSE IF (STRING(I:I) .EQ. ']') THEN
	    in_tok=.false.
	    CONTINUE
	  ELSE IF (STRING(I:I) .EQ. '%') THEN
	    in_tok=.false.
	    IF (NO_LET) THEN
	      SP_FLG=.TRUE.
	      GO TO 20
	    ELSE
	      GO TO 15
	    END IF
	  ELSE IF (STRING(I:I) .EQ. '(') THEN
	    in_tok=.false.
	    IND_FND=.TRUE.
	    J=1
	  ELSE IF (STRING(I:I) .EQ. ')') THEN
	    in_tok=.false.
	    CLS_IND=.TRUE.
	    J=1
	  ELSE IF (STRING(I:I) .EQ. ';') THEN
	    in_tok=.false.; ntok=0
	    CONTINUE
	  ELSE IF (STRING(I:I) .EQ. '=') THEN
	    IND_FND=.TRUE.
	    CLS_IND=.TRUE.
	    J=1
	    in_tok=.false.; ntok=0
	  ELSE
	    ! STRING(I:I) is NOT a special character or blank
	    IF (.NOT.IND_FND) THEN
	      VAR_NAM(J:J)=STRING(I:I)
	      !print*,'VAR_NAM(J:J)=',VAR_NAM(J:J)
	      J=J+1
	      NO_LET=.FALSE.
	      in_tok=.true.
	    ELSE IF (.NOT.CLS_IND) THEN
	      IF (STRING(I:I) .EQ. ':') THEN
	        ARR_FLG=.TRUE.
	      ELSE IF (STRING(I:I) .EQ. ',') THEN
 	        J=1
	        K=K+1
	        in_tok=.false.
	      ELSE
	        IND(K)(J:J)=STRING(I:I)
	        J=J+1
	        in_tok=.true.
	      END IF
	    ELSE
	      VALUE(J:J)=STRING(I:I)
	      J=J+1
	      VAL_LET=.FALSE.
	      in_tok=.true.
	    END IF
     	  END IF
   12   CONTINUE

   15   CONTINUE

	IF (VAR_NAM(1:1) .EQ. ' ') THEN
	  SP_FLG=.TRUE.
	  GO TO 20
	END IF

	IF (LCMP(VAR_NAM,'CommentBegin',12) .or.
     &      LCMP(VAR_NAM,'/*',2)) THEN
	  LInCommentMode=.TRUE.
	  VAR_NAM(1:1)=' '; SP_FLG=.TRUE.
	ELSE IF (LCMP(VAR_NAM,'CommentEnd',10) .OR.
     &           LCMP(VAR_NAM,'*/',2)) THEN
	  LInCommentMode=.FALSE.
	  VAR_NAM(1:1)=' '; SP_FLG=.TRUE.
          GO TO 20
	END IF
	IF (LInCommentMode) GO TO 20

	NVTOK=ntok
	RETURN
	END SUBROUTINE GET_EQ

C***********************************************************************

	SUBROUTINE MB_LOH(task)
	use lohpars_mod

C      THIS SUBROUTINE OPENS AND CLOSES THE BATCH FILE

	IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: filnamloh
	!CHARACTER*132 filnamloh

	LOGICAL exist
	!INTEGER ustack,l,task,ICLEN
	INTEGER l,task,ICLEN

	IF (task .EQ. 1) THEN

c new unit number
	   ustack = ustack+1

	   CALL CACCEPT(filnamloh,' ','Enter file name: ')
	   l=ICLEN(filnamloh)
	   filnamloh(l+1:l+4)='.jou'
	   l = l+4
	   INQUIRE(file=filnamloh(1:l),exist=exist)

	   IF (exist) THEN
	      WRITE(*,*) 'Reading from journal file (',
     &		   filnamloh(1:l),')'
	      OPEN(ustack,FILE=filnamloh(1:l),STATUS='OLD')
	   ELSE
	      WRITE(*,*) '[mb_loh] Couldn''t find journal file (',
     &		   filnamloh(1:l),')'
	   ENDIF
	ELSE IF (task .EQ. 2) THEN
	   CLOSE (unit=ustack,STATUS='KEEP')
	ENDIF

	RETURN
	END SUBROUTINE MB_LOH

C***********************************************************************

	 SUBROUTINE JOU_LOH(task)

C      THIS SUBROUTINE OPENS AND CLOSES THE JOU FILE

	 IMPLICIT NONE

	 CHARACTER(len=MacosCharLen) :: filnamjou
	 LOGICAL exist
	 INTEGER ICLEN,l,task

	 IF (task .EQ. 1) THEN
	 CALL CACCEPT(filnamjou,' ','Enter journal file name: ')
	   l=ICLEN(filnamjou)
	   filnamjou(l+1:l+4)='.jou'

	   INQUIRE(file=filnamjou,exist=exist)

	   IF (exist) THEN
	     OPEN(27,FILE=filnamjou,STATUS='OLD')
	     CLOSE (unit=27,STATUS='DELETE')
	   END IF
	   OPEN(27,FILE=filnamjou,STATUS='NEW')
c sws
	   WRITE (*,*) 'Writing journal file ',
     &	      filnamjou(1:ICLEN(filnamjou)), '...'

	 ELSE IF (task .EQ. 2) THEN
	   CLOSE (unit=27,STATUS='KEEP')
	   WRITE (*,*) 'Closing journal file ',
     &	      filnamjou(1:ICLEN(filnamjou)), '.'
	 END IF
c sws
	 RETURN
	 END SUBROUTINE JOU_LOH

C***********************************************************************

	SUBROUTINE MFil1(nElt,nCMCoords,nOutCord,
     &		         nFirstRay,nLastRay,nRay)

	use param_mod

	IMPLICIT NONE

	INTEGER nFirstRay,nLastRay,nRay,nCMCoords(0:mElt),nElt,lRay,i,k,
     &	        nOutCord

 130	FORMAT(/'%   This M-file was created by the MACOS program. It'/
     &	'%   contains a linear optical model.'/)
 140	FORMAT(A12,I6,';')
 141	FORMAT(A12,I3,1x,$)
 142	FORMAT(I3,1X,$)
 143	FORMAT(I3,'];')

	WRITE(8,130)
	WRITE(8,140)'      nElt= ',nElt
	WRITE(8,140)'      nRay= ',nRay
	WRITE(8,141)'nCmatCoord=[',nCMCoords(0)
	k=nElt-1
	DO 1 i=1,k
	  WRITE(8,142)nCMCoords(i)
 1	CONTINUE
	WRITE(8,143)nCMCoords(nElt)

	RETURN
	END SUBROUTINE MFil1

C***********************************************************************

	SUBROUTINE NasFil1(nElt,ifInPar,nOutCord,nEC,
     &                     tot_row,tot_col)
        Integer :: nElt,nOutCord,nEC(nElt), tot_row,tot_col,ib
	Logical :: ifInPar

        tot_row=nOutCord ! total rows

        ! Compute total columns
	if (ifInPar) then
          tot_col=7   ! always 7 columns for src when included
	else
	  tot_col=0
	end if
        Do ib=1,nElt
          tot_col=tot_col+iabs(nEC(ib))
        End Do

        !WRITE(8,*) 'iabs(nEC)=',iabs(nEC(1:nElt))

        ! Write out DMI matrix header
	WRITE(8,301)'DMI','C','0','2','2','2',' ',tot_row,tot_col
 301    FORMAT(A3,5X,A8,A8,A8,A8,A8,A8,I8,I8)
        END SUBROUTINE NasFil1

C***********************************************************************

	SUBROUTINE NasFil1_old(nElt,nOutCord,nEC,
     &                     tot_row,tot_col)
	Integer :: nElt,nOutCord,nEC(nElt),tot_row,tot_col, ib

	tot_row=nOutCord ! total rows

	! Compute total columns
	tot_col=7   ! always 7 columns for src
	Do ib=1,nElt
          tot_col=tot_col+iabs(nEC(ib))
	End Do

	!WRITE(8,*) 'iabs(nEC)=',iabs(nEC(1:nElt))

	! Write out DMI matrix header
	WRITE(8,*)'DMI  C  0  2  2  2',tot_row,tot_col
	END SUBROUTINE NasFil1_old

C***********************************************************************

	! Save or write out C-matrix in Nastran fixed long format.
        ! 1) Write out a sensitivity matrix block, for one surface.
        !    Each block is dimension nOutCord x nEcrd
        ! 2) The block is written out in column-oriented sparse matrix
        !    form, non-zero values only.

        SUBROUTINE NasFil2(Cout,nOutCord,nEcrd,irow,icol)
        IMPLICIT NONE

        SREAL Cout(nOutCord,nEcrd)
        Integer :: nOutCord,nEcrd,irow,icol,ir,ic,i
        Logical :: ColStrt,NzBlkStrt
        Integer :: icol_loc,tot_col_nz,curr_col_nz,off

        icol_loc=icol

	! The input block is stored in cout(1:nOutCord,1:nEcrd)
        Do ic=1,nEcrd
          ! Check whether it's an empty col
          tot_col_nz=0
          Do ir=1,nOutCord
            if (cout(ir,ic)/=0d0) then
              tot_col_nz=tot_col_nz+1
            end if
          End Do
          if (tot_col_nz==0) then
	    icol_loc=icol_loc+1
            go to 211  ! col is empty
	  end if
          ColStrt=.TRUE.; NzBlkStrt=.TRUE.
          off=0
          WRITE(8,201)'DMI*','C',icol_loc
          off=off+8+16+16
          icol_loc=icol_loc+1   ! global column index
          curr_col_nz=0
          Do ir=1,nOutCord
            if (cout(ir,ic)/=0d0 .and. ColStrt) then
              ColStrt=.FALSE.
            end if
            if (cout(ir,ic)/=0d0 .and. NzBlkStrt) then
              WRITE(8,202)ir  ! row id
              off=off+16
              NzBlkStrt=.FALSE.
            else if (cout(ir,ic)==0d0 .and. (.not.NzBlkStrt)) then
              NzBlkStrt=.TRUE.
              ! if cout(ir,ic)==0d0 and NzBlkStrt, then do nothing
            end if
            if (cout(ir,ic)/=0d0) then
              WRITE(8,203)cout(ir,ic)
              off=off+16;
              curr_col_nz=curr_col_nz+1  ! total non-zeros already written
            end if
            !if ((off==72 .OR. off==56) .AND.
            if (off==72 .AND.
     &          curr_col_nz < tot_col_nz) then
	      WRITE(8,*)' '  ! end a line
              ! add continuation '*' at beginning of next line
              WRITE(8,204)'*'
              WRITE(8,205)' '  ! skip 7 chars
              off=8
            end if
          End Do  ! writing column loop
          WRITE(8,*)' '  ! new line
 211      Continue
        End Do
 201    FORMAT(A4,4X,A16,I16,$)
 202    FORMAT(I16,$)  ! row index
 203    FORMAT(D16.9,$) ! data value
 204    FORMAT(A1,$)
 205    FORMAT(A7,$)
        END SUBROUTINE NasFil2

C***********************************************************************

	! Save or write out C-mtx in Nastran fixed long format.
        ! 1) Write out a sensitivity matrix block, for one surface.
        !    Each block is dimension nOutCord x nEcrd
        ! 2) The block is written out in column-oriented sparse matrix
        !    form, non-zero values only.

	SUBROUTINE NasFil2_old_2(Cout,nOutCord,nEcrd,irow,icol)
        IMPLICIT NONE

        SREAL Cout(nOutCord,nEcrd)
        Integer :: nOutCord,nEcrd,irow,icol,ir,ic,i
        Logical :: ColStrt,NzBlkStrt
        Integer :: icol_loc,tot_col_nz,curr_col_nz,off

        icol_loc=icol

	! The input block is stored in cout(1:nOutCord,1:nEcrd)
        Do ic=1,nEcrd
          ! Check whether it's an empty col
          tot_col_nz=0
          Do ir=1,nOutCord
            if (cout(ir,ic)/=0d0) then
              tot_col_nz=tot_col_nz+1
            end if
          End Do
          if (tot_col_nz==0) go to 211  ! col is empty
          ColStrt=.TRUE.; NzBlkStrt=.TRUE.
	  off=0
          WRITE(8,201)'DMI* ','C',icol_loc,' '
	  off=off+8+16+16
          icol_loc=icol_loc+1   ! global column index
          curr_col_nz=0
          Do ir=1,nOutCord
            if (cout(ir,ic)/=0d0 .and. ColStrt) then
              ColStrt=.FALSE.
            end if
            if (cout(ir,ic)/=0d0 .and. NzBlkStrt) then
              WRITE(8,202)ir,' '  ! row id
	      off=off+16
              NzBlkStrt=.FALSE.
            else if (cout(ir,ic)==0d0 .and. (.not.NzBlkStrt)) then
              NzBlkStrt=.TRUE.
              ! if cout(ir,ic)==0d0 and NzBlkStrt, then do nothing
            end if
            if (cout(ir,ic)/=0d0) then
              WRITE(8,203)cout(ir,ic),' '
	      off=off+17;
              curr_col_nz=curr_col_nz+1  ! total non-zeros already written
            end if
	    if ((off==73 .OR. off==59) .AND.
     &          curr_col_nz < tot_col_nz) then
	      ! add continuation '*'
	      do i=1,(80-off-1)
	        WRITE(8,204)' '
	      end do
	      WRITE(8,204)'*'
              WRITE(8,*)' '  ! new line
	      WRITE(8,205)' '  ! skip first 8 chars
	      off=8
	    end if
          End Do  ! writing column loop
          WRITE(8,*)' '  ! new line
 211      Continue
        End Do
 201    FORMAT(1X,A7,A16,I15,A1,$)
 202    FORMAT(I15,A1,$)  ! row index
 203    FORMAT(D16.9,A1,$) ! data value
 204    FORMAT(A1,$)
 205    FORMAT(A8,$)
        END SUBROUTINE NasFil2_old_2

C***********************************************************************

	SUBROUTINE MFil2(Cout,nOutCord,nEcrd,irow,icol,iRay,iElt,
     &		         iExMode,jRay)

	IMPLICIT NONE

	INTEGER nOutCord,nEcrd,irow,icol,i,j,nline,jrow,jcol,iRay,iElt,
     &	iExMode,jRay
	LOGICAL isZero
	SREAL Cout(nOutCord,nEcrd)

 108	FORMAT(1P,7D23.15)
 151	FORMAT('CMAT(',i6,',',i4,':',i4,')=[',$)
 152	FORMAT(1X,E15.7,$)
 153	FORMAT('];')
 154	FORMAT('% iElt=',i3,' iRay=',i6)

	IF (iExMode.EQ.1) THEN
	  WRITE(8)((Cout(i,j),j=1,nEcrd),i=1,nOutCord)
	ELSE IF (iExMode.EQ.2) THEN
	  DO 1 j=1,nEcrd
	    WRITE(8,108)(Cout(i,j),i=1,nOutCord)
  1	  CONTINUE
	ELSE IF (iExMode.EQ.3) THEN
	  isZero=.TRUE.
	  DO 3 i=1,nOutCord
	    DO 2 j=1,nEcrd
	      IF (Cout(i,j).NE.0e0) isZero=.FALSE.
  2	    CONTINUE
  3	  CONTINUE
	  IF (.NOT.isZero) THEN
	    WRITE(8,154)iElt,jRay
	    jrow=irow
	    DO 5 i=1,nOutCord
	      jcol=icol
	      DO 4 j=1,nEcrd
	        IF (j.EQ.1) THEN
	          nline=icol+nEcrd-1
	          WRITE (8,151)jrow,icol,nline
	        END IF
	        WRITE (8,152)Cout(i,j)
	        IF (j.EQ.nEcrd) THEN
	          WRITE (8,153)
	        END IF
	        jcol=jcol+1
  4	      CONTINUE
	      jrow=jrow+1
  5	    CONTINUE
	  END IF
	END IF
	RETURN
	END SUBROUTINE MFil2

C***********************************************************************

	SUBROUTINE MFil3(Cout,nOutCord,nEcrd,irow,icol,iRay,iElt,
     &	iExMode,iNode,jRay)

	IMPLICIT NONE

	INTEGER nOutCord,nEcrd,irow,icol,i,j,nline,jrow,jcol,iRay,iElt,
     &	iExMode,iNode,jRay

	SREAL Cout(nOutCord,nEcrd)

 108	FORMAT(1P,7D23.15)
 152	FORMAT(1X,E15.7,$)
 153	FORMAT('];')
 154	FORMAT('% iNode=',i6,',iRay=',i6,', iElt=',i6,':')

	IF (iExMode.EQ.1) THEN
	  WRITE(8)((Cout(i,j),j=1,nEcrd),i=1,nOutCord)
	ELSE IF (iExMode.EQ.2) THEN
	  DO 3 j=1,nEcrd
	    WRITE(8,108)(Cout(i,j),i=1,nOutCord)
  3	  CONTINUE
	END IF
	RETURN
	END SUBROUTINE MFil3

C***********************************************************************

	SUBROUTINE CFiltRead(ifFilterData)

C       A routine to read filter data.

	use param_mod
	use src_mod
	use cfiles_mod

	IMPLICIT NONE

	LOGICAL ifFilterData,exist
	INTEGER i,j,ICLEN
	CHARACTER*132 FiltName

 100	FORMAT(' Filter includes information for ',i3,
     &	' wavelengths and ',i2,' star types.')

C  Read .filt file
	CALL CACCEPT(filtfil,' ','Enter filter file name: ')
	filtfil=filtfil(1:ICLEN(filtfil))//'.filt'
	INQUIRE (file=filtfil,exist=exist)
	IF (.NOT.exist) THEN
	  WRITE(*,*)' Filter file '//filtfil(1:ICLEN(filtfil))
     &	   //' does not exist.'
	  ifFilterData=.FALSE.
	  RETURN
	ELSE
	  OPEN (44,FILE=filtfil,STATUS='OLD',FORM='FORMATTED')
	  READ (44,*,ERR=3,END=3)FiltName,nWavelen
	  IF (nWavelen.GT.mWavelen) GO TO 2
	  nStarType=0
	  READ(44,*,ERR=3,END=3)(FiltWavelen(j),j=1,nWavelen)
	  DO 1 i=1,mStarType
	    READ(44,*,ERR=3,END=2)(FiltFlux(j,i),j=1,nWavelen)
	    nStarType=nStarType+1
 1	  CONTINUE
 2	  CONTINUE
	  WRITE(*,100)nWavelen,nStarType
	  ifFilterData=.TRUE.
	  CLOSE(44)
	  RETURN
	END IF

C  Here if read failed
 3	CONTINUE
	WRITE(*,*)
     &	' Filter file not in expected format. Data not accepted.'
	ifFilterData=.FALSE.
	CLOSE(44)
	RETURN

	END SUBROUTINE CFiltRead

C***********************************************************************

	SUBROUTINE CFiltLoad(ifFilterData)

C       A routine to load filter data interactively.

	use param_mod
	use src_mod
	use cfiles_mod

	IMPLICIT NONE

	LOGICAL ifFilterData,exist
	INTEGER i,j,ICLEN
	CHARACTER(len=MacosCharLen) :: ANS,FiltName

 100	FORMAT(' Filter includes information for ',i3,
     &	' wavelengths and ',i2,' star types.')

C  Load filter data
	CALL CACCEPT(filtfil,filnam(1:ICLEN(filnam)),
     &	'Enter filter name: ')
	filtfil=filtfil(1:ICLEN(filtfil))//'.filt'

	nStarType=1
	nWavelen=0
	DO 4 i=1,mWavelen
	  CALL DACCEPT_S(FiltWavelen(i),0d0,1,
     &	  'Enter next wavelength in base units (0 when done):')
	  IF (FiltWavelen(i).LE.0d0) GO TO 5
	    nWavelen=nWavelen+1
	    CALL DACCEPT_S(FiltFlux(i,1),0d0,1,
     &	    'Enter flux value at this wavelength:')
 4 	CONTINUE
 5 	CONTINUE
	IF (nWavelen.GT.0) THEN
	  ifFilterData=.TRUE.
	ELSE
	  ifFilterData=.FALSE.
	END IF
	WRITE(*,100)nWavelen,nStarType

	RETURN
	END SUBROUTINE CFiltLoad

C***********************************************************************

	SUBROUTINE CFiltSave(ifFilterData)

C       A routine to save filter data.

	use param_mod
	use src_mod
	use cfiles_mod

	IMPLICIT NONE

	LOGICAL ifFilterData,exist
	INTEGER i,j,ICLEN
	CHARACTER(len=MacosCharLen) :: ANS,FiltName

 100	FORMAT(' Filter includes information for ',i3,
     &	' wavelengths and ',i2,' star types.')

C  Save filter data
	IF (.NOT.ifFilterData) THEN
	  WRITE(*,*)' No filter data loaded.'
	ELSE
	  CALL CACCEPT(filtfil,filnam(1:ICLEN(filnam)),
     &	  'Enter new filter file name: ')
	  filtfil=filtfil(1:ICLEN(filtfil))//'.filt'
	  INQUIRE (file=filtfil,exist=exist)
	  IF (exist) THEN
	    WRITE(*,*)' Filter file '//filtfil(1:ICLEN(filtfil))
     &	    //' being overwritten.'
	    OPEN (UNIT=44,FILE=filtfil,STATUS='OLD')
	    CLOSE (UNIT=44,STATUS='DELETE')
	  END IF

	  FiltName=filtfil
	  OPEN (44,FILE=filtfil,STATUS='NEW',FORM='FORMATTED')
	  WRITE (44,*,ERR=3)FiltName,nWavelen
	  nStarType=0
	  WRITE(44,*,ERR=3)(FiltWavelen(j),j=1,nWavelen)
	  DO 1 i=1,nStarType
	    WRITE(44,*,ERR=3)(FiltFlux(j,i),j=1,nWavelen)
 1	  CONTINUE
	  WRITE(*,*)' Filter file created.'
	  CLOSE(44)
	  RETURN
	END IF

C  Here if read failed
 3	CONTINUE
	WRITE(*,*)' Write errors prevented creation of filter file.'
	CLOSE(44)
	RETURN

	END SUBROUTINE CFiltSave

C***********************************************************************

	SUBROUTINE EltConvert(EltType,EltID,SrfType)

	IMPLICIT NONE
	INTEGER EltType,EltID,SrfType

	IF (EltType.EQ.1) THEN
	  EltID=1
	  SrfType=2
	ELSE IF (EltType.EQ.2) THEN
	  EltID=1
	  SrfType=1
	ELSE IF (EltType.EQ.3) THEN
	  EltID=2
	  SrfType=1
	ELSE IF (EltType.EQ.4) THEN
	  EltID=3
	  SrfType=2
	ELSE IF (EltType.EQ.5) THEN
	  EltID=11
	  SrfType=2
	ELSE IF (EltType.EQ.6) THEN
	  EltID=4
	  SrfType=2
	ELSE IF (EltType.EQ.7) THEN
	  EltID=5
	  SrfType=2
	ELSE IF (EltType.EQ.8) THEN
	  EltID=6
	  SrfType=2
	ELSE IF (EltType.EQ.9) THEN
	  EltID=7
	  SrfType=2
	ELSE IF (EltType.EQ.10) THEN
	  EltID=8
	  SrfType=2
	ELSE IF (EltType.EQ.11) THEN
	  EltID=1
	  SrfType=2
	ELSE IF (EltType.EQ.12) THEN
	  EltID=1
	  SrfType=3
	ELSE IF (EltType.EQ.13) THEN
	  EltID=12
	  SrfType=2
	ELSE IF (EltType.EQ.14) THEN
	  EltID=1
	  SrfType=8
	ELSE IF (EltType.EQ.15) THEN
	  EltID=11
	  SrfType=8
	ELSE IF (EltType.EQ.16) THEN
	  EltID=1
	  SrfType=4
	ELSE IF (EltType.EQ.17) THEN
	  EltID=11
	  SrfType=4
	ELSE IF (EltType.EQ.18) THEN
	  EltID=1
	  SrfType=5
	ELSE IF (EltType.EQ.19) THEN
	  EltID=11
	  SrfType=5
	ELSE IF (EltType.EQ.20) THEN
	  EltID=6
	  SrfType=3
	ELSE IF (EltType.EQ.21) THEN
	  EltID=6
	  SrfType=8
	ELSE IF (EltType.EQ.22) THEN
	  EltID=10
	  SrfType=2
	ELSE IF (EltType.EQ.23) THEN
	  EltID=10
	  SrfType=3
	ELSE IF (EltType.EQ.24) THEN
	  EltID=9
	  SrfType=2
	ELSE IF (EltType.EQ.25) THEN
	  EltID=1
	  SrfType=6
	ELSE IF (EltType.EQ.26) THEN
	  EltID=1
	  SrfType=7
	END IF

	RETURN
	END SUBROUTINE EltConvert

C***********************************************************************

	SUBROUTINE EltSurfCompat(iElt,EltID_FLG,SrfType_FLG)
	use param_mod
	use elt_mod

	IMPLICIT NONE

	INTEGER iElt
	LOGICAL EltID_FLG,SrfType_FLG

	IF (EltID(iElt).EQ.2) THEN
	  IF (SrfType(iElt).NE.1) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.3) THEN
	  ! extended to allow reference to have GridData, 04/2009, -jzlou
	  IF ((SrfType(iElt).GT.2) .AND. (SrfType(iElt)/=9)) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.4) THEN
	  IF (SrfType(iElt).EQ.5) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.5) THEN
	  IF (SrfType(iElt).EQ.5) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.6) THEN
	  IF (SrfType(iElt).EQ.5) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.7) THEN
	  IF ((SrfType(iElt).NE.1).AND.(SrfType(iElt).NE.2)) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.8) THEN
	  IF ((SrfType(iElt).NE.1).AND.(SrfType(iElt).NE.2)) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.9) THEN
	  IF (SrfType(iElt).EQ.5) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	ELSE IF (EltID(iElt).EQ.10) THEN
	  IF ((SrfType(iElt).NE.2).AND.(SrfType(iElt).NE.3)) THEN
	    EltID_FLG=.FALSE.
	    SrfType_FLG=.FALSE.
	  END IF
	END IF

	RETURN
	END SUBROUTINE EltSurfCompat

C***********************************************************************

	Subroutine RecordElementID(iElt,eName,lastCall)
        use param_mod
        use elt_mod

        implicit none

        Integer :: iElt
        Character(len=*) eName
 	Logical :: lastCall

        Logical :: LCMP, LCMP2
        Logical, save :: firstSegFound=.false., inSeg=.false.,
     &                   first_entry = .true.

	if (first_entry) then
	  m_LensGrp(1,1:maxLensGrp) = -1
          m_LensGrpElts = 0
	  m_nLensGrp = 0
	  first_entry = .false.
	end if

	if (lastCall) then
	  ! final call to this routine
	  if (m_nLensGrp > 0) then
	    m_LensGrp(2,1:m_nLensGrp)
     &	     = m_LensGrp(1,1:m_nLensGrp)
     &         +m_LensGrpElts(1:m_nLensGrp)-1
	  end if
	  return
	end if

        If (LCMP(eName,'Seg',3)) then
          if (.not. firstSegFound) then
             m_iFirstPMSeg = iElt
             firstSegFound = .true.
             inSeg = .true.
          end if
	  !
	Else If (LCMP(eName,'Lens',4)) then
	   if (inSeg) then
            m_iLastPMSeg = iElt-1
            inSeg = .false.
          end if

	  ! Currently supports a max of four lens groups
	  if (LCMP2(eName,'Lens-1',6)) then
	    if (m_LensGrp(1,1)==-1) then
	      m_LensGrp(1,1) = iElt
	      m_nLensGrp = m_nLensGrp + 1
	    end if
	    m_LensGrpElts(1) = m_LensGrpElts(1) + 1
	    !
	  else if (LCMP2(eName,'Lens-2',6)) then
	    if (m_LensGrp(1,2)==-1) then
	      m_LensGrp(1,2) = iElt
	      m_nLensGrp = m_nLensGrp + 1
	    end if
	    m_LensGrpElts(2) = m_LensGrpElts(2) + 1
	    !
	  else if (LCMP2(eName,'Lens-3',6)) then
            if (m_LensGrp(1,3)==-1) then
              m_LensGrp(1,3) = iElt
              m_nLensGrp = m_nLensGrp + 1
            end if
            m_LensGrpElts(3) = m_LensGrpElts(3) + 1
	  else if (LCMP2(eName,'Lens-4',6)) then
            if (m_LensGrp(1,4)==-1) then
              m_LensGrp(1,4) = iElt
              m_nLensGrp = m_nLensGrp + 1
            end if
            m_LensGrpElts(4) = m_LensGrpElts(4) + 1
	  end if
	  !
        Else  ! not 'Seg' and not 'Lens'
          if (inSeg) then
            m_iLastPMSeg = iElt-1
            inSeg = .false.
          end if

	  if (LCMP(eName,'PM',2) .or.
     &        LCMP(eName,'PrimaryMirror',13)) then
	      m_iPM = iElt
          else if (LCMP(eName,'CtrSeg',6) .or.
     &             LCMP(eName,'CenterSegment',13)) then
            m_iCenterPMSeg = iElt
          else if (LCMP(eName,'SM',2) .or.
     &             LCMP(eName,'SecondaryMirror',15)) then
            m_iSM = iElt
          else if (LCMP(eName,'TM',2) .or.
     &             LCMP(eName,'TertiaryMirror',15)) then
            m_iTM = iElt
          else if (LCMP(eName,'FSM',3)) then
            m_iFSM = iElt
	  else if (LCMP(eName,'Pickoff',7)) then
	    m_iPickoff = iElt
          else if (LCMP(eName,'ExitPupil',9)) then
            m_iExitPupil = iElt
          else if (LCMP(eName,'FocalPlane',10) .or.
     &             LCMP(eName,'ImagePlane',10)) then
            m_iFocalPlane = iElt
          end if
        End If
        End Subroutine RecordElementID

C***********************************************************************
!
! Parse ValBuf to extract numerical values and put them into ValArr
!
	Subroutine Get_Values(ValBuf,ValArr,nVal)
        Implicit none

        CHARACTER(len=MacosValLen) :: ValBuf
        Real*8, dimension(:) :: ValArr
        Character(len=100):: tok
        Integer :: nVal,c,i,j
        Logical :: nextTok

        nextTok=.true. ! when true, 'state' is outside of a token
        tok(1:100)=' '; c=1; j=1; nVal=0
        Do
          If (c==MacosCharLen) exit  ! done
          If (ValBuf(c:c)==' ' .and. c<MacosCharLen) Then
            If (.not.nextTok) Then
              ! a token is complete
              !print*,'** tok = ',tok
              Read(tok,*) ValArr(j)
              !print*,'** ValArr(j) =',ValArr(j)
              j=j+1; tok(1:100)=' '; nextTok=.true.
	      nVal=nVal+1
            End If
            If (c<MacosCharLen) c=c+1
          Else
            If (nextTok) Then
              i=1; nextTok=.false.
            End If
            !print*,'** i =',i
            tok(i:i)=ValBuf(c:c); i=i+1; c=c+1
          End If
        End Do
        End Subroutine Get_Values

C***********************************************************************
C*******End file iosub.inc**********************************************
C***********************************************************************


