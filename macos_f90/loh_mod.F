C***********************************************************************
C       Begin file loh_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C
C  Data modules
C

	Module loh_I1
	Implicit none

	LOGICAL ChfRayDir_FLG,ChfRayPos_FLG,IndRef_FLG,Wavelen_FLG,
     &    Flux_FLG,GridType_FLG,Aperture_FLG,Obscratn_FLG,nGridpts_FLG,
     &    xGrid_FLG,nElt_FLG,nSeg_FLG,width_FLG,gap_FLG,SegCoord_FLG,
     &    zSource_FLG,yGrid_FLG,Extinc_FLG,SegXgrid_FLG,WaveUnits_FLG,
     &    BaseUnits_FLG,nPetals_FLG,radCtr_FLG,nRings_FLG,RingRad_FLG
	End Module loh_I1


	Module loh_I2
	Implicit none

	LOGICAL EltName_FLG,EltType_FLG,fElt_FLG,eElt_FLG,
     &    AsphCoef_FLG,ZernCoef_FLG,nZernCoef_FLG,ZernModes_FLG,
     &    pMon_FLG,xMon_FLG,yMon_FLG,zMon_FLG,
     &    lMon_FLG,psiElt_FLG,VptElt_FLG,RptElt_FLG,IndRefE_FLG,
     &    zElt_FLG,PropType_FLG,nECoord_FLG,ApType_FLG,ApVec_FLG,
     &    nObs_FLG,xObs_FLG,TElt_FLG,SegApType_FLG,ZernCent_FLG,
     &    KrElt_FLG,KcElt_FLG,LensArrayType_FLG,LensArrayWidth_FLG,
     &    h1HOE_FLG,h2HOE_FLG,OrderHOE_FLG,WaveHOE_FLG,PinHole_FLG,
     &    ExtincE_FLG,nCoat_FLG,CoatIndx_FLG,CoatThk_FLG,AnaCoef_FLG,
     &    MonCoef_FLG,RuleWidth_FLG,UDSrfCoef_FLG,UDSrfType_FLG,
     &    EltID_FLG,SrfType_FLG,UDSrfParam_FLG,ZernTypeL_FLG,
     &    nGridMat_FLG,GridMat_FLG,GridSrfdx_FLG,GlassCoef_FLG,
     &    GlassElt_FLG,nAmplMat_FLG,AmplSrfdx_FLG,LApSunBlades
	End Module loh_I2


	Module loh_I3
	Implicit none

	LOGICAL nOutCord_FLG,Tout_FLG
	End Module loh_I3
