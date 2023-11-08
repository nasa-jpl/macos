C***********************************************************************
C       Begin file traceutil_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2008, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

	Module traceutil_mod
	use param_mod, only : mpts,mDrawRay,czero_
	use elt_mod, only : mPolySide,mRefRay,mPupRay

	Implicit none
        Save

	Integer :: crnt_Elt,crnt_traced_ray,m_err_ray,tAsphCoef=0,
     &    nValidRays,iDrawRay_global,pgplotDrawPlane,numGradInd,
     &    RayTgtElt
	Integer, dimension(:,:), allocatable :: opdRayMask,gradInd_kmax
        Real*8 :: RRayPosHis(3,mPupRay)
	Real*8, dimension(:), allocatable :: taperFactor
        Real*8 :: NewIndRef,PolyVtx3D(3,mPolySide),UDxDrawGrid(3),
     &            UDyDrawGrid(3),UDzDrawGrid(3),OPDRefRayLen
        Real*8, dimension(:), allocatable :: VarStrtIndRef,VarCurIndRef
	Real*8, dimension(:,:), allocatable :: CRIncidDir,CROutgoDir,
     &       CRIncidPosNom,CRIncidPos,CRIncidSurfN,CRXLoc,CRYLoc,CRNhat 
	Real*8, dimension(:,:,:), allocatable :: PolyVtx3DSav,RayPosHist 
        Real*8, dimension(:,:,:,:), allocatable :: gradIndRayPath
	  ! note that gradIndRayPath now only supports two gradInd elems
	  ! see allocate for gradIndRayPath below 
        Logical, dimension(:,:), allocatable :: LRayOKHist
	Logical, dimension(:), allocatable :: Poly3DApVecL
        Logical :: LBWKnoPrint,LNewIndRef,ifSave3DApVec,ifSRT,
     &             LCalcWFElt,ifCtraceCalled,Rx_FEXCentrFlg,
     &             LShowDTime,LudDrawGrid,OPDRefRayLen_FLG,LRayTgtElt,
     &             ShowMetData_FLG,saveOPD_Flg,SaveRayPosHist_Flg,
     &             DumpRayPosHist_Flg,SaveVis3d_dat_Flg 
	
	! Supporting variable IndRef, etc.
	Integer, parameter :: mLenslet=250,mGradIndRaySegs=256
	Integer :: nHexLenlets=91,   ! hex lens array is hardwired,
     &             nRecLenlets=2601  ! 51*51 rec lens array by default

	! Supporting FZP grating models
	Integer :: Dielectric=1,  ! dielectric lens (grating)
     &             Metal=2	  ! metal lens (grating)
        Integer, dimension(:), allocatable :: GratingMaterial 
	Complex*16, dimension(:), allocatable :: rayTC ! transmission coefficient for rays

	Type LensletStr
          Integer :: Ltri,Rtri,Xtri,RecXtri,RecYtri
          Real*8 :: IndRef,Xctr,Yctr,NgbDist
        End Type LensletStr
        !
        Logical, dimension(:), allocatable :: IsVarIndRefElt,
     &                                        IsGrinRefElt
	Logical :: LUseChfRayIfOK,IsOverDraw
        Type(LensletStr), dimension(:,:), allocatable :: LensletStrArray
        Integer, dimension(:), allocatable :: TotalLensletsArray,
     &                         HexLensletIdx,RecLensletIdx,iElt2iGradInd

	Type IndRefGridStr
	  Real*8, dimension(:,:), allocatable :: IndRefGrid  
	  Real*8 :: width
	End Type IndRefGridStr
	Type(IndRefGridStr), dimension(:), allocatable :: IndRefGridArr 

        Type GrinRefStr
          Integer :: modelType
          Real*8, dimension(10) :: modelCoef
	  Real*8 :: lensDepth,dt
          Real*8, dimension(:,:), allocatable :: raypt  ! ray points array 
	End Type GrinRefStr 
	Type(GrinRefStr), dimension(:), allocatable :: GrinRefArr

	Type RayStateStr
          Integer :: rayId
          Real*8, dimension(:,:), allocatable :: InciPos,InciDir,
     &                                           InciSrfN
          Real*8, dimension(:), allocatable :: PrevIndRef,CurIndRef
        End Type RayStateStr
        Type(RayStateStr), target :: thisRay

	Type ChfRayHistStr
          Logical :: keepHist
	  Integer :: tot_srfs
	  Integer, dimension(:), allocatable :: SrfId
	  Real*8, dimension(:,:), allocatable :: RayPosAtSrf,RayDirAtSrf
	End Type ChfRayHistStr
        Type(ChfRayHistStr) :: ChfRayHist

	! Support grating with variable rulewidth
	Integer, parameter :: nRuleWidth=1000
	Real*8 :: ruleWidthArr(nRuleWidth)

	Contains

	Subroutine ray_mod_init(max_elt,max_ray)
	Integer :: max_elt,max_ray

	If (allocated(taperFactor)) Then
          deallocate(taperFactor,CRIncidDir,CRIncidPosNom,CRIncidPos,
     &          CROutgoDir,CRIncidSurfN,CRXLoc,CRYLoc,CRNhat,
     &          VarStrtIndRef,VarCurIndRef,IndRefGridArr,GrinRefArr,
     &          IsVarIndRefElt,LensletStrArray,TotalLensletsArray,
     &          HexLensletIdx,RecLensletIdx,Poly3DApVecL,
     &          IsGrinRefElt,PolyVtx3DSav,RayPosHist,LRayOKHist,
     &          gradIndRayPath,opdRayMask,gradInd_kmax,iElt2iGradInd,
     &          rayTC,GratingMaterial,stat=m_err_ray)
          if (m_err_ray /= 0)
     &      call macos_memory_failure
     &               ('ray_mod_init: deallocate failed!')
        End If

	allocate(taperFactor(max_ray),CRIncidDir(3,max_elt),
     &           CROutgoDir(3,max_elt),CRIncidPosNom(3,max_elt),
     &           CRIncidPos(3,max_elt),CRIncidSurfN(3,max_elt),
     &           CRXLoc(3,max_elt),CRYLoc(3,max_elt),CRNhat(3,max_elt),
     &           VarStrtIndRef(max_ray),VarCurIndRef(max_ray),
     &           IndRefGridArr(max_elt),IsVarIndRefElt(0:max_elt),
     &           TotalLensletsArray(max_elt),Poly3DApVecL(max_elt),
     &           LensletStrArray(mLensLet,max_elt),GrinRefArr(max_elt),
     &           HexLensletIdx(mLensLet*3),RecLensletIdx(mLenslet*2),
     &		 PolyVtx3DSav(3,mPolySide,max_elt),
     &           RayPosHist(3,max_ray,0:max_elt),
     &           LRayOKHist(max_ray,0:max_elt),
     &           opdRayMask(mpts,mpts),
     &           IsGrinRefElt(0:max_elt),gradInd_kmax(mDrawRay,2),
     &           gradIndRayPath(3,mGradIndRaySegs,mDrawRay,2),
     &           iElt2iGradInd(max_elt),rayTC(max_ray),
     &           GratingMaterial(max_elt),stat=m_err_ray)
	if (m_err_ray /= 0)
     &      call macos_memory_failure
     &                 ('ray_mod_init: deallocate failed!')

	taperFactor=1d0  ! default value, no tapering
	IsVarIndRefElt=.FALSE.; TotalLensletsArray=0
        LensletStrArray%Ltri=0; LensletStrArray%Rtri=0
        LensletStrArray%Xtri=0; IsGrinRefElt=.FALSE.
        LensletStrArray%RecXtri=0; LensletStrArray%RecYtri=0
        LensletStrArray%IndRef=0d0
	CRIncidDir=0d0; CROutgoDir=0d0; CRIncidPosNom=0d0
	CRIncidPos=0d0; CRIncidSurfN=0d0; iElt2iGradInd=0
	CRXLoc=0d0; CRYLoc=0d0; CRNhat=0d0
	VarStrtIndRef=0d0; VarCurIndRef=0d0
        HexLensletIdx=0; RecLensletIdx=0; LShowDTime=.TRUE.
	iDrawRay_global=0; IsOverDraw=.FALSE.; rayTC=czero_
        GratingMaterial=0
        ChfRayHist%keepHist=.FALSE.; ChfRayHist%tot_srfs=0
   	Rx_FEXCentrFlg=.TRUE.  ! default 
        SaveRayPosHist_Flg=.FALSE. ! default 
        SaveVis3d_dat_Flg=.FALSE.  ! default
	End Subroutine ray_mod_init


#if 1
        ! For matlab/mex interface
        Subroutine SetIndRefGrid(iElt,gridData,width)
        Implicit none
        Save

        Integer :: iElt,nx,ny,i,j
        Real*8, dimension(:,:) :: gridData
        Real*8 :: width

        IsVarIndRefElt(iElt)=.TRUE.
        ny=size(gridData,1); nx=size(gridData,2)
        !print*,'SetIndRefGrid: nx,ny=',nx,ny

        ! allocate or re-allocate IndRefGridArr(iElt)%IndRefGrid
        !If (.NOT.allocated(IndRefGridArr(iElt)%IndRefGrid)) Then

        If (.NOT.allocated(IndRefGridArr(iElt)%IndRefGrid)
     &      .OR. (size(IndRefGridArr(iElt)%IndRefGrid,1)/=nx)
     &      .OR. (size(IndRefGridArr(iElt)%IndRefGrid,2)/=ny))
     &  Then
          if (allocated(IndRefGridArr(iElt)%IndRefGrid)) then
            deallocate(IndRefGridArr(iElt)%IndRefGrid)
          end if
          allocate(IndRefGridArr(iElt)%IndRefGrid(ny+1,nx+1))
        End If

        IndRefGridArr(iElt)%IndRefGrid=1d0
        IndRefGridArr(iElt)%IndRefGrid(1:ny,1:nx)=gridData(1:ny,1:nx)
        IndRefGridArr(iElt)%width=width
        End Subroutine SetIndRefGrid
#endif
	End Module traceutil_mod

