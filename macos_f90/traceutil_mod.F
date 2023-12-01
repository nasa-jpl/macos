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
	! Supporting variable IndRef, etc.
	Integer, parameter :: mLenslet        = 250,
     &                        mGradIndRaySegs = 256

     	! Support grating with variable rulewidth
	Integer, parameter :: nRuleWidth = 1000

        integer :: m_err_ray
        Save

        integer :: nHexLenlets,    ! hex lens array is hardwired,
     &             nRecLenlets     ! 51*51 rec lens array by default

        ! Supporting FZP grating models
        integer :: Dielectric,     ! dielectric lens (grating)
     &             Metal           ! metal lens (grating)

        Integer :: crnt_Elt,crnt_traced_ray,tAsphCoef,
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

        Integer, dimension(:), allocatable :: GratingMaterial
	Complex*16, dimension(:), allocatable :: rayTC ! transmission coefficient for rays

	Type LensletStr
          Integer :: Ltri,Rtri,Xtri,RecXtri,RecYtri
          Real*8 :: IndRef,Xctr,Yctr,NgbDist
        End Type LensletStr

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

! >>>> ------------- used in MACOS SRT cmd -- but memory is never allocated -- FIX (ToDo)

	Type RayStateStr
          Integer :: rayId
          Real*8, dimension(:,:), allocatable :: InciPos,InciDir,
     &                                           InciSrfN
          Real*8, dimension(:), allocatable :: PrevIndRef,CurIndRef
        End Type RayStateStr
        Type(RayStateStr), target :: thisRay

! <<<< -------------

	Type ChfRayHistStr
          Logical :: keepHist
	  Integer :: tot_srfs
	  Integer, dimension(:), allocatable :: SrfId
	  Real*8, dimension(:,:), allocatable :: RayPosAtSrf,RayDirAtSrf
	End Type ChfRayHistStr
        Type(ChfRayHistStr) :: ChfRayHist

	Real*8 :: ruleWidthArr(nRuleWidth)

        contains

        subroutine ray_mod_init(max_elt,max_ray)
          implicit none
          integer, intent(in) :: max_elt,max_ray

          integer :: i
          ! ----------------------------------------- IndRefGridArr
          ! de-allocate allocations inside type definition
          if (allocated(IndRefGridArr)) then
            do concurrent (i=1:size(IndRefGridArr,1))
              if (allocated(IndRefGridArr(i)%IndRefGrid))
     &          deallocate(IndRefGridArr(i)%IndRefGrid)
            end do
          end if

          ! ----------------------------------------- GrinRefArr
          if (allocated(GrinRefArr)) then
            do concurrent (i=1:size(GrinRefArr,1))
              if (allocated(GrinRefArr(i)%raypt))
     &          deallocate(GrinRefArr(i)%raypt)
            end do
          end if

          ! -----------------------------------------

	  if (allocated(taperFactor)) then
            deallocate(taperFactor,CRIncidDir,CRIncidPosNom,CRIncidPos,
     &            CROutgoDir,CRIncidSurfN,CRXLoc,CRYLoc,CRNhat,
     &            VarStrtIndRef,VarCurIndRef,IndRefGridArr,GrinRefArr,
     &            IsVarIndRefElt,LensletStrArray,TotalLensletsArray,
     &            HexLensletIdx,RecLensletIdx,Poly3DApVecL,
     &            IsGrinRefElt,PolyVtx3DSav,RayPosHist,LRayOKHist,
     &            gradIndRayPath,opdRayMask,gradInd_kmax,iElt2iGradInd,
     &            rayTC,GratingMaterial,stat=m_err_ray)
            if (m_err_ray /= 0)
     &        call macos_memory_failure
     &                 ('ray_mod_init: deallocate failed!')
          end if

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
     &        call macos_memory_failure
     &                   ('ray_mod_init: deallocate failed!')

          ! ----------------------------------------- ChfRayHist
          if (allocated(ChfRayHist%SrfId)) then
            deallocate(ChfRayHist%SrfId,
     &                 ChfRayHist%RayPosAtSrf,
     &                 ChfRayHist%RayDirAtSrf,stat=m_err_ray)
            if (m_err_ray /= 0)
     &        call macos_memory_failure
     &                   ('ray_mod_init: deallocate failed!')
          end if

          ! initialisations
          call ray_mod_init_vars()

        End Subroutine ray_mod_init


        subroutine ray_mod_init_vars()
          implicit none
          integer :: i

          ! ----------------------------------------- ChfRayHist
          if (allocated(ChfRayHist%SrfId)) then
            deallocate(ChfRayHist%SrfId,
     &                 ChfRayHist%RayPosAtSrf,
     &                 ChfRayHist%RayDirAtSrf,stat=m_err_ray)
            if (m_err_ray /= 0)
     &        call macos_memory_failure
     &                   ('ray_mod_init: deallocate failed!')
          end if

          ChfRayHist = ChfRayHistStr(keepHist    = .FALSE.,
     &                               tot_srfs    = 0)

          ! ----------------------------------------- GrinRefArr
          if (allocated(GrinRefArr)) then
            do concurrent (i=1:size(GrinRefArr,1))
              if (allocated(GrinRefArr(i)%raypt))
     &          deallocate(GrinRefArr(i)%raypt)
            end do

            GrinRefArr(:) = GrinRefStr(modelType = 0,
     &                                 modelCoef = 0d0,
     &                                 lensDepth = 0d0,
     &                                 dt        = 0d0)
          end if

          ! ----------------------------------------- IndRefGridArr
          if (allocated(IndRefGridArr)) then
            do concurrent (i=1:size(IndRefGridArr,1))
              if (allocated(IndRefGridArr(i)%IndRefGrid))
     &          deallocate(IndRefGridArr(i)%IndRefGrid)
            end do

            IndRefGridArr(:) = IndRefGridStr(width = 0d0)
          end if
          ! -----------------------------------------


        ruleWidthArr(:) = 0d0   ! grating with variable rule width

        crnt_Elt                = 0
        crnt_traced_ray         = 0
        m_err_ray               = 0
        tAsphCoef               = 0
        nValidRays              = 0
        iDrawRay_global         = 0
        pgplotDrawPlane         = 0
        numGradInd              = 0
        RayTgtElt               = 0

        nHexLenlets             =   91     ! hex lens array is hardwired,
        nRecLenlets             = 2601     ! 51*51 rec lens array by default

        ! Supporting FZP grating models
        Dielectric              = 1        ! dielectric lens (grating)
        Metal                   = 2        ! metal lens (grating)


        RRayPosHis(:,:)         = 0d0

        NewIndRef               = 0d0
        PolyVtx3D(:,:)          = 0d0
        UDxDrawGrid(:)          = 0d0
        UDyDrawGrid(:)          = 0d0
        UDzDrawGrid(:)          = 0d0
        OPDRefRayLen            = 0d0

        LBWKnoPrint             = .FALSE.
        LNewIndRef              = .FALSE.
        ifSave3DApVec           = .FALSE.
        ifSRT                   = .FALSE.
        LCalcWFElt              = .FALSE.
        ifCtraceCalled          = .FALSE.
        Rx_FEXCentrFlg          = .TRUE.
        LShowDTime              = .TRUE.
        LudDrawGrid             = .FALSE.
        OPDRefRayLen_FLG        = .FALSE.
        LRayTgtElt              = .FALSE.
        ShowMetData_FLG         = .FALSE.
        saveOPD_Flg             = .FALSE.
        SaveRayPosHist_Flg      = .FALSE.
        DumpRayPosHist_Flg      = .FALSE.
        SaveVis3d_dat_Flg       = .FALSE.
        IsOverDraw              = .FALSE.
        LUseChfRayIfOK          = .FALSE.
        ! ------------------

        opdRayMask(:,:)         = 0
        gradInd_kmax(:,:)       = 0
        GratingMaterial(:)      = 0

        PolyVtx3DSav(:,:,:)     = 0d0
        RayPosHist(:,:,:)       = 0d0
        gradIndRayPath(:,:,:,:) = 0d0

        TotalLensletsArray(:)   = 0
        HexLensletIdx(:)        = 0
        RecLensletIdx(:)        = 0
        iElt2iGradInd(:)        = 0

        LRayOKHist(:,:)         = .FALSE.
        Poly3DApVecL(:)         = .FALSE.

        IsVarIndRefElt(:)       = .FALSE.
        IsGrinRefElt(:)         = .FALSE.

        CRIncidDir(:,:)         = 0d0
        CROutgoDir(:,:)         = 0d0
        CRIncidPosNom(:,:)      = 0d0
        CRIncidPos(:,:)         = 0d0
        CRIncidSurfN(:,:)       = 0d0
        CRXLoc(:,:)             = 0d0
        CRYLoc(:,:)             = 0d0
        CRNhat(:,:)             = 0d0
        VarStrtIndRef           = 0d0
        VarCurIndRef            = 0d0
        rayTC(:)                = (0d0, 0d0)    ! transmission coefficient for rays

        taperFactor(:)          = 1d0           ! default: no tapering

        LensletStrArray(:,:) = LensletStr(Ltri    = 0,
     &                                    Rtri    = 0,
     &                                    Xtri    = 0,
     &                                    RecXtri = 0,
     &                                    RecYtri = 0,
     &                                    IndRef  = 0d0,
     &                                    Xctr    = 0d0,
     &                                    Yctr    = 0d0,
     &                                    NgbDist = 0d0)


        end Subroutine ray_mod_init_vars


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

