C***********************************************************************
C       Begin file elt_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

        Module elt_mod
        use Constants, only:  SQRT_3,  SQRT_5,  SQRT_6,  SQRT_7,  SQRT_8, SQRT_10,
     &                       SQRT_11, SQRT_12, SQRT_14, SQRT_18, SQRT_20, SQRT_22
        use param_mod

        Implicit none
        Save

        !
        ! --------- Basic Surface Definition ---------
        !

        INTEGER,  PARAMETER:: maxEltNameLength = 32

        REAL*8,   DIMENSION(:,:),  ALLOCATABLE :: VptElt,RptElt,psiElt,VptRx,psiRx,RptRx
        REAL*8,   DIMENSION(:),    ALLOCATABLE :: fElt,eElt,zElt,KcElt,KrElt,KrRx,KcRx
        REAL*8,   DIMENSION(:,:,:),ALLOCATABLE :: TElt
        INTEGER,  DIMENSION(:),    ALLOCATABLE :: nECoord,SrfType,EltType,EltID
        CHARACTER(LEN=maxEltNameLength),
     &            DIMENSION(:),    ALLOCATABLE :: EltName
        LOGICAL,  DIMENSION(:),    ALLOCATABLE :: LEltOK

      !                                         EltID
        INTEGER, PARAMETER :: ReflectorElt    =   1,
     &                        FocalPlaneElt   =   2,
     &                        ReferenceElt    =   3,
     &                        HOEElt          =   4,
     &                        GratingElt      =   5,
     &                        RefractorElt    =   6,
     &                        ObscuringElt    =   7,
     &                        ReturnElt       =   8,
     &                        NSRefractorElt  =   9,
     &                        LensArrayElt    =  10,
     &                        SegmentElt      =  11,  ! Reflector
     &                        NSReflectorElt  =  12,
     &                        TrGratingElt    =  13,
     &                        RfPolarizerElt  =  14,
     &                        TrPolarizerElt  =  15,
     &                        CGHNullPlateElt =  16,
     &                        DoeTrGratingElt =  17,
     &                        mEltTypes       =  17   ! defines available Element Types

        CHARACTER(LEN=12),PARAMETER:: EltTypeName(mEltTypes)=
     &                     [CHARACTER(LEN=12) ::
     &                       'Reflector',         ! 1
     &                       'FocalPlane',        ! 2
     &                       'Reference',         ! 3
     &                       'HOE',               ! 4
     &                       'Grating',           ! 5
     &                       'Refractor',         ! 6
     &                       'Obscuring',         ! 7
     &                       'Return',            ! 8
     &                       'NSRefractor',       ! 9
     &                       'LensArray',         !10
     &                       'Segment',           !11
     &                       'NSReflector',       !12
     &                       'TrGrating',         !13
     &                       'RfPolarizer',       !14
     &                       'TrPolarizer',       !15
     &                       'CGHNullPlate',      !16
     &                       'DoeTrGrating']      !17  Diffractive Optical Element (DOE)

        INTEGER, PARAMETER :: SrfType_Flat              =  1,
     &                        SrfType_Conic             =  2,
     &                        SrfType_Aspheric          =  3,
     &                        SrfType_Monomial          =  4,
     &                        SrfType_Interpolated      =  5,
     &                        SrfType_Anamorphic        =  6,
     &                        SrfType_UserDefined       =  7,
     &                        SrfType_Zernike           =  8,
     &                        SrfType_GridData          =  9,
     &                        SrfType_Toric             = 10,
     &                        SrfType_AsGrData          = 11,
     &                        SrfType_MonGrData         = 12,
     &                        SrfType_ZrnGrData         = 13,
     &                        SrfType_FreeForm          = 14,
     &                        mSrfType                  = 14

!-->  Per-ray status codes
        INTEGER, PARAMETER :: RayStat_OK       = 0
        INTEGER, PARAMETER :: RayStat_Obscured = 1
        INTEGER, PARAMETER :: RayStat_Miss     = 2
        INTEGER, PARAMETER :: RayStat_Bracket  = 3
        INTEGER, PARAMETER :: RayStat_MaxIter  = 4
        INTEGER, PARAMETER :: RayStat_Undef    = 5
!<--

        ! Surface types with Gridded Surface Parameters
        INTEGER, PARAMETER :: GridTypeAll(7) = (/SrfType_GridData,
     &                                           SrfType_AsGrData,
     &                                           SrfType_MonGrData,
     &                                           SrfType_ZrnGrData,
     &                                           SrfType_FreeForm,
     &                                           SrfType_Monomial,
     &                                           SrfType_Zernike/)

        ! Surfaces having a Local Coord. System defined (required)
        INTEGER, PARAMETER :: SrfTypes_LCSYS(8) = [SrfType_Monomial,
     &                                             SrfType_Zernike,
     &                                             SrfType_AsGrData,
     &                                             SrfType_MonGrData,
     &                                             SrfType_ZrnGrData,
     &                                             SrfType_GridData,
     &                                             SrfType_FreeForm,
     &                                             SrfType_Anamorphic]

        CHARACTER(LEN=12),PARAMETER:: SrfTypeName(mSrfType)=
     &                     [CHARACTER(LEN=12) ::
     &                       'Flat',               ! 1
     &                       'Conic',              ! 2
     &                       'Aspheric',           ! 3
     &                       'Monomial',           ! 4
     &                       'Interpolated',       ! 5
     &                       'Anamorphic',         ! 6
     &                       'UserDefined',        ! 7
     &                       'Zernike',            ! 8
     &                       'GridData',           ! 9
     &                       'Toric',              !10
     &                       'AsGrData',           !11
     &                       'MonGrData',          !12
     &                       'ZrnGrData',          !13
     &                       'FreeForm']           !14

        ! Old Element Definitions:
        INTEGER, PARAMETER, DIMENSION(26,2), PRIVATE:: EltType2EltIDSrfType=
                     !  EltType           SrfType                           EltID       Element                        EltType         SurfType
     &      RESHAPE([ ReflectorElt,   SrfType_Conic,          ! [[ 1, 2], !  1   <==  Conic Mirror                  1 (Reflector)   |  2 (ConSrf)
     &                ReflectorElt,   SrfType_Flat,           !  [ 1, 1], !  2   <==  Flat Mirror                   1 (Reflector)   |  1 (FlatSrf)
     &                FocalPlaneElt,  SrfType_Flat,           !  [ 2, 1], !  3   <==  Focal Plane                   2 (FocalPln)    |  1 (FlatSrf)
     &                ReferenceElt,   SrfType_Conic,          !  [ 3, 2], !  4   <==  Reference Surface             3 (RefSrf)      |  2 (ConSrf)
     &                SegmentElt,     SrfType_Conic,          !  [11, 2], !  5   <==  Mirror Segment               11 (Segment)     |  2 (ConSrf)
     &                HOEElt,         SrfType_Conic,          !  [ 4, 2], !  6   <==  Holographic Optical Elmnt     4 (HOE)         |  2 (ConSrf)
     &                GratingElt,     SrfType_Conic,          !  [ 5, 2], !  7   <==  Straight-Ruled Grating        5 (Grating)     |  2 (ConSrf)
     &                RefractorElt,   SrfType_Conic,          !  [ 6, 2], !  8   <==  Conic Refractor               6 (Refractor)   |  2 (ConSrf)
     &                ObscuringElt,   SrfType_Conic,          !  [ 7, 2], !  9   <==  Obscuring Surface             7 (ObsSrf)      |  2 (ConSrf)
     &                ReturnElt,      SrfType_Conic,          !  [ 8, 2], ! 10   <==  Return Surface                8 (ReturnSrf)   |  2 (ConSrf)
     &                ReflectorElt,   SrfType_Conic,          !  [ 1, 2], ! 11   <==  General Conic Reflector       1 (Reflector)   |  2 (ConSrf)
     &                ReflectorElt,   SrfType_Aspheric,       !  [ 1, 3], ! 12   <==  Aspheric Reflector            1 (Reflector)   |  3 (AsphSrf)
     &                NSReflectorElt, SrfType_Conic,          !  [12, 2], ! 13   <==  Non-Sequential Reflector     12 (Reflector)   |  2 (ConSrf)
     &                ReflectorElt,   SrfType_Zernike,        !  [ 1, 8], ! 14   <==  Zernike Reflector             1 (Reflector)   |  8 (MonSrf)
     &                SegmentElt,     SrfType_Zernike,        !  [11, 8], ! 15   <==  Zernike Segment              11 (Segment)     |  8 (MonSrf)
     &                ReflectorElt,   SrfType_Monomial,       !  [ 1, 4], ! 16   <==  Monomial Reflector            1 (Reflector)   |  4 (MonSrf)
     &                SegmentElt,     SrfType_Monomial,       !  [11, 4], ! 17   <==  Monomial Segment             11 (Segment)     |  4 (MonSrf)
     &                ReflectorElt,   SrfType_Interpolated,   !  [ 1, 5], ! 18   <==  Interpolated Reflector        1 (Reflector)   |  5 (IntSrf)
     &                SegmentElt,     SrfType_Interpolated,   !  [11, 5], ! 19   <==  Interpolated Segment         11 (Segment)     |  5 (IntSrf)
     &                RefractorElt,   SrfType_Aspheric,       !  [ 6, 3], ! 20   <==  Aspheric Refractor            6 (Refractor)   |  3 (AsphSrf)
     &                RefractorElt,   SrfType_Zernike,        !  [ 6, 8], ! 21   <==  Zernike Refractor             6 (Refractor)   |  4 (MonSrf)
     &                LensArrayElt,   SrfType_Conic,          !  [10, 2], ! 22   <==  Conic Lens Array             10 (LensArray)   |  2 (ConSrf)
     &                LensArrayElt,   SrfType_Aspheric,       !  [10, 3], ! 23   <==  Aspheric Lens Array          10 (LensArray)   |  3 (AsphSrf)
     &                NSRefractorElt, SrfType_Conic,          !  [ 9, 2], ! 24   <==  Non-Sequential Refractor      9 (NSRefractor) |  2 (ConSrf)
     &                ReflectorElt,   SrfType_Anamorphic,     !  [ 1, 6], ! 25   <==  Anamorphic Reflector          1 (Reflector)   |  6 (AnaSrf)
     &                ReflectorElt,   SrfType_UserDefined ],  !  [ 1, 7]] ! 26   <==  User-Defined Reflector        1 (Reflector)   |  7 (UDSrf)
     &              [26,2], ORDER=[2,1])


        ! Surface types that require a Srf. Coord. Systems
        INTEGER, PARAMETER :: SrfType_RequireSrfCSYS(7)=
     &                      [SrfType_Zernike,
     &                       SrfType_Monomial,
     &                       SrfType_GridData,
     &                       SrfType_AsGrData,
     &                       SrfType_MonGrData,
     &                       SrfType_ZrnGrData,
     &                       SrfType_FreeForm]

        !
        ! --- Surfaces: Asphere, Zernikes, Monomials, Grid, UDS, AsGrid, MonGrid, FreeForm -------
        !
        LOGICAL,  DIMENSION(:),     ALLOCATABLE:: ifComputeDZ                              ! interpolated Srf. Data (Grid)
        REAL*8,   DIMENSION(:),     ALLOCATABLE:: ZernRad, lMon, GridSrfdx, ZernAnnuRatio
        REAL*8,   DIMENSION(:,:),   ALLOCATABLE:: AsphCoef,MonCoef,ZernCoef,pMon,
     &                                            xMon,yMon,zMon,
     &                                            AnaCoef,ZernCent,ZernXGrid,ZernYGrid,
     &                                            UDSrfCoef,UDSrfParam,
     &                                            xMonRx,yMonRx,zMonRx,pMonRx,ZernCoefRx
        REAL*8,   DIMENSION(:,:,:), ALLOCATABLE:: GridMat

        ! --- Added for multisurface (grid data coord frame):
        REAL*8,   DIMENSION(:),     ALLOCATABLE:: lData,lDataRx
        REAL*8,   DIMENSION(:,:),   ALLOCATABLE:: pData,xData,yData,zData,
     &                                            pDataRx,xDataRx,yDataRx,zDataRx

        ! --- FreeForm surface (FF monomial term: separate from design Mon term) ---
        !     lFF         : FF scale length (sentinel: lFF <= 0 means FF term inactive)
        !     pFF         : FF data center (3 x mElt)
        !     xFF,yFF     : FF local x,y axes (3 x mElt each)
        !     zFF         : FF surface axis (3 x mElt)
        !     FFZernCoef  : FF Zernike coefficients as entered (mFFCoef x mElt)
        !     FFZernTypeL : FF Zernike type index per element (uses ZernType_* constants)
        !     FFCoef      : FF monomial coefficients after ZerntoMon conversion (mFFCoef x mElt)
        !     MonZernCoef : Mon Zernike coefficients as entered (mMonCoef x mElt)
        !     MonZernTypeL: Mon Zernike type index per element (uses ZernType_* constants)
        !     ifMon       : Mon term active for this element
        !     ifFF        : FF term active for this element
        !     ifGridTerm  : Grid term active for this element
        !     (all three set at prescription-read time from SrfType and sentinels)

        INTEGER,  PARAMETER :: mFFCoef = 120   ! max FF monomial coefficients

        REAL*8,   DIMENSION(:),     ALLOCATABLE:: lFF
        REAL*8,   DIMENSION(:,:),   ALLOCATABLE:: pFF,xFF,yFF,zFF,
     &                                            FFZernCoef,FFCoef,
     &                                            MonZernCoef
        INTEGER,  DIMENSION(:),     ALLOCATABLE:: FFZernTypeL,MonZernTypeL
        LOGICAL,  DIMENSION(:),     ALLOCATABLE:: ifMon,ifFF,ifGridTerm
        ! --- end FreeForm ---

        INTEGER,  DIMENSION(:),     ALLOCATABLE:: UDSrfType,ZernaType,ZernTypeL,
     &                                            nGridMat,iEltToGridSrf,
     &                                            nDP                                      ! interpolated Srf. Data (Grid)


        INTEGER, PARAMETER :: mZernType          =  11,
     &                        mZernModes         =  66,  ! max. # of Zernike Coefficients (8th radial mode)
     &                        nZernModes_Default =  45,  ! default if "nZernCoef" is not defined in Rx
     &                        mAsphCoef          =   9,  ! max terms = 9, up to term h^20
     &                        nAsphCoef_Default  =   4,  ! default if "nAsphCoef" is not defined in Rx
     &                        mAnaCoef           =   4,
     &                        mMonCoef           = 120,
     &                        mUDSrfParam        =  10

        INTEGER, PARAMETER ::
     &          ZernType_ANSI            =  1,
     &          ZernType_BornWolf        =  2,
     &          ZernType_Fringe          =  3,   ! incl. CODE V's 37th mode
     &          ZernType_NormANSI        =  4,
     &          ZernType_NormBornWolf    =  5,
     &          ZernType_NormFringe      =  6,
     &          ZernType_NormHex         =  7,   ! ??? ToDo: Implemented?
     &          ZernType_NormNoll        =  8,
     &          ZernType_NormAnnularNoll =  9,
     &          ZernType_Noll            = 10,
     &          ZernType_ExtFringe       = 11

        INTEGER, PARAMETER :: ZernType_MaxMode(mZernType) =
     &          [mZernModes,  ! ZernType_ANSI
     &           mZernModes,  ! ZernType_BornWolf
     &                   37,  ! ZernType_Fringe   (CODE V notation: #37 = 12th radial mode)
     &           mZernModes,  ! ZernType_NormANSI
     &           mZernModes,  ! ZernType_NormBornWolf
     &                   37,  ! ZernType_NormFringe
     &                    0,  ! ZernType_NormHex          implemented ???
     &           mZernModes,  ! ZernType_NormNoll
     &           mZernModes,  ! ZernType_NormAnnularNoll
     &           mZernModes,  ! ZernType_Noll
     &           mZernModes]  ! ZernType_ExtFringe



        INTEGER, DIMENSION(5), PARAMETER ::
     &          ZernType_rms_normalised = (/ZernType_NormANSI,
     &                                      ZernType_NormBornWolf,
     &                                      ZernType_NormFringe,
     &                                      ZernType_NormNoll,
     &                                      ZernType_NormAnnularNoll/)

        INTEGER, DIMENSION(5), PARAMETER ::
     &           ZernType_not_normalised = (/ZernType_ANSI,        !  1,
     &                                       ZernType_BornWolf,    !  2,
     &                                       ZernType_Fringe,      !  3,
     &                                       ZernType_Noll,        ! 10,
     &                                       ZernType_ExtFringe/)  ! 11

        CHARACTER(LEN=15), PARAMETER :: ZernTypeNameL(mZernType)=
     &         [CHARACTER(LEN=15) ::
     &           'ANSI',             !  1   was Malacara
     &           'BornWolf',         !  2
     &           'Fringe',           !  3
     &           'NormANSI',         !  4   was Malacara
     &           'NormBornWolf',     !  5
     &           'NormFringe',       !  6
     &           'NormHex',          !  7
     &           'NormNoll',         !  8
     &           'NormAnnularNoll',  !  9   Annular Zernike
     &           'Noll',             ! 10   Noll Ordering
     &           'ExtFringe']        ! 11   Extended Fringe Ordering

        REAL(8), PARAMETER :: NORM_RMS_PARAM_ANSI(mZernModes) =
     &             (/    1d0,      2d0,      2d0,  SQRT_6,   SQRT_3,   SQRT_6,
     &               SQRT_8,   SQRT_8,   SQRT_8,   SQRT_8,   SQRT_10,  SQRT_10,
     &               SQRT_5,   SQRT_10,  SQRT_10,  SQRT_12,  SQRT_12,  SQRT_12,
     &               SQRT_12,  SQRT_12,  SQRT_12,  SQRT_14,  SQRT_14,  SQRT_14,
     &               SQRT_7,   SQRT_14,  SQRT_14,  SQRT_14,      4d0,      4d0,
     &                   4d0,      4d0,      4d0,      4d0,      4d0,      4d0,
     &               SQRT_18,  SQRT_18,  SQRT_18,  SQRT_18,      3d0,  SQRT_18,
     &               SQRT_18,  SQRT_18,  SQRT_18,  SQRT_20,  SQRT_20,  SQRT_20,
     &               SQRT_20,  SQRT_20,  SQRT_20,  SQRT_20,  SQRT_20,  SQRT_20,
     &               SQRT_20,  SQRT_22,  SQRT_22,  SQRT_22,  SQRT_22,  SQRT_22,
     &               SQRT_11,  SQRT_22,  SQRT_22,  SQRT_22,  SQRT_22,  SQRT_22/)

        !
        ! --------- Basic Surface Definition: Material ---------
        !

        Integer, parameter :: HartmannGlass       = 0,  ! Not used
     &                        SellmeierGlass      = 1,  ! Schott's definition (SELLMEIER3TA)
     &                        SchottGlass         = 2,  ! Schott = Laurent model   ("Old Schott" Model)
     &                        HerzbergerGlass     = 3,  ! HERZBRGR4X2
     &                        SELLMEIER_2TA_GLASS = 4,  ! Sellmeier, 2 Terms with Constant A (5 Parameters)
     &                        KINGSLAKE1_GLASS    = 5,  ! Kingslake with 2 Terms (5 Parameters)
     &                        CauchyGlass         = 6,  ! Cauchy
     &                        MAX_GLASS_MODELS    = 6

        Integer, parameter :: nGLASS_COEF(0:6) = [0,    ! HartmannGlass
     &                                            6,    ! SellmeierGlass
     &                                            6,    ! SchottGlass
     &                                            6,    ! HerzbergerGlass
     &                                            5,    ! SELLMEIER_2TA_GLASS
     &                                            5,    ! KINGSLAKE1_GLASS
     &                                            3]    ! CauchyGlass

        CHARACTER(LEN=10), PARAMETER :: GlassTypeName(0:6) =
     &                         [CHARACTER(LEN=10) ::
     &                          'Hartmann', 'Sellmeier',
     &                          'Schott', 'Herzberger', 'Sellm2TA',
     &                          'Kingslake', 'Cauchy']

        Integer, parameter :: mCoat=10,
     &           	      mPrvGlass=5, ! max private glasses defined in Rx
     &                        mPrvGlassEntry=20, ! max wavlen/index entries for
					 	 ! each glass in catalog
     &                        mES=10  ! max edge-sensors per element
!-->  2026-07-12 Sprint 2D S3: segmented-primary met trusses need a met
!     surface on EVERY segment (19 segs + M2 + M3 = 21 > 20) and 6
!     beams/segment (2-ring 19*6+6 = 126 ~ 128; 3-ring 228 > 128).
        Integer, parameter :: mMetSrf=64,  ! max surfaces with metrology
     &                        mMetBeam=48, ! max met beams on each met srf.
     &                        mMetPos=mMetBeam,  ! max metrology points on each met srf.
     &                        mSysMetBeam=512,    ! max total system met measurements
!<--
     &                        mRefRay=5,   ! max # reference rays, chief + 4 marginal
     &                        mPupRay=9,   ! max # pupil rays, chief + 8 marginal
     &                        mIndRef=5,   !
     & 			      mDOE=6       ! max # DOE elements
        Integer :: maxAsphCoef  ! default to 4 terms

        Logical :: ifInMTFCalc,BWK_FLG,LStopAtEltAutoSet,LSavRayL,SrcLF_FLG,
     &             CalcSrcZFrame_FLG,DupElt_FLG
        Logical, dimension(:), allocatable :: LFastRayTrace,LRayOK,
     &                  LRayPass,LEltObs,LGlass,
     &                  LTaperedElt,LApTypeDefined,LApVecDefined,
     &                  LEltMRP_FLG,LUpdateTElt_FLG,LxObsDef,
     &                  ifAmplSrf,ifOPDModGrid,ifGridDataDefined,
!-->  GridFile= parsed before nGridMat= for this element: the grid read
!     is DEFERRED (GridInit sizes its read from nGridMat) and resolved
!     by the nGridMat= handler or by ChkDf2.  GridFileBin remembers the
!     ascii/binary dispatch decided at parse time.
     &                  GridFilePending,GridFileBin
!<--
        Logical, dimension(:,:), allocatable :: LRayToElt

!-->  Per-ray status tracking arrays
        Integer, dimension(:), allocatable :: RayStatus, RayFailElt
        Character(len=40), dimension(:), allocatable :: RayFailMsg
!<--

!-->  Whole-line Rx comment round-trip (PLAN section 0 item 3):
!     comment LINES are captured at load time (GET_EQ) and re-emitted
!     by SAVE (ViewData writers).  Anchor = element context at
!     capture: 0 = header, iElt = within that element's block.
!     LRxCommentCapture gates capture to Rx loads only -- GET_EQ is
!     also called by the MOD dialog and must not capture there.
        Integer, parameter :: mRxComment = 500
        Character(len=MacosCharLen) :: RxCommentLine(mRxComment)
        Integer :: RxCommentAnchor(mRxComment)
        Integer :: nRxComment = 0
        Integer :: RxCommentCtx = 0
        Logical :: LRxCommentCapture = .FALSE.
!<--

        Character*32 :: ApTypeName(0:15), ObsTypeName(-10:10),
     &                  PropTypeName(15),
     &                  BaseUnits,WaveUnits

        Character*32, dimension(:), allocatable ::
     &                GlassElt, GlassName,EltMRP

        Integer :: nElt,nRay,nOutCord,nSeg,spcOption,nRefRay,nDOE,
     &             SegApType,SegApVec(3),nMetMeas

        Integer, dimension(:,:,:), allocatable :: metBeamFlg

        Integer, dimension(:,:), allocatable :: SegCoord,
     &           RayToSegMap,RayID,RayWFMap,ObsType,
     &           IWORK,EltGrp,ObsTypeRx,tMetSrfPos,
     &           nPolyObsVtx3DPending

        Integer, dimension(:), allocatable :: LnkElt,
     &         EltToSegMap,SegToEltMap,PropType,ApType,
     &         nObs,iEltToiWF,OpType,iEltToIntSrf,
     &         LensArrayType,nCoatElt,
     &         GridSrfOrder,
     &         GlassType,EltCoat,nObsRx,nAmplMat,passRayIds,
     &         iEltToAmplSrf,iEltToMetSrf,nMetPos,tMetSrf,ntMetPos

        Integer, dimension(:), allocatable :: nRayToSeg,NSCnt,EdgeSensor

        SREAL, dimension(:,:,:), allocatable :: dxidxim1,dxidui
        SREAL, dimension(:,:), allocatable :: MWFFT
        SREAL, dimension(:), allocatable :: XFP,YFP

        Real*8 :: Wavelen,Flux, zSource,Tout(7,7),focalRatio,plateScale,
     &            WavelenUM,WavelenConv,WaveBU,CBM,CWM,CWB,WavelenRx,
     &            ToutRx(7,7),sysLOS(2),sysLOSAng(2),CurRefRayDir(3)

        Real*8, dimension(:,:,:), allocatable :: ObsVec,XYZSRF,
     &          DZSRF,ApHexBound,PolyApVtx,PolyApEqn,ObsVecRx,
     &          TEltRx,AmplMat,SrfMetPos,EdgeSensVec

        Real*8, dimension(:,:,:,:), allocatable :: ObsHexBound,
     &          PolyObsVtx,PolyObsEqn,PolyObsVtx3DPending

        Real*8, dimension(:,:), allocatable ::
     &        RayPos, RayDir, ApVec, xObs,
     &        h1HOE,h2HOE, CoatIndxElt,CoatThkElt,
     &        yObs, zObs,DWORK,
     &        GlassCoef, GlassTable,MTFMat,MTFMat2,
     &        IndRefArr,ExtincArr,EltCoatThk,
     &	      ArrIndRef,SrfMetMea
	! IndRefArr is for coating layers, ArrIndRef is for lens element
	! having several IndRef for different wavelengths

        Real*8, dimension(:), allocatable ::
     &          RayL,dxElt,IndRef,CumRayL,CumCRL,
     &          LensArrayWidth, OrderHOE,WaveHOE,PinHole,
     &          Extinc,RuleWidth,AmplSrfdx,
     &          CumRayLsav,metMeasBuf

        Real*8 :: Fn,  ! Fresnel number
     &	          OPDdx, OPDdy, ! exit pupil OPD grid spacing
     &            RayDirBuff(3,5),rayHt,costh,efl, ! used in 'EFL' command
     &            ChfRayDirAtElt(3),ChfRayPosAtElt(3),BeamSizeAtElt,
     &            RefRayDir(3,mRefRay)

	! Zernike normallizing coefs (Malacara)
	Real*8 :: ZnNmFac_m(mZernModes)

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

        Integer, parameter :: mPolySide=128

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

C***********************************************************************
C  ParseZernType -- shared parser for ZernType / FFZernType / MonZernType
C  prescription values.
C
C  Maps a user-supplied typename string (e.g. "Noll", "NormNoll",
C  "BornWolf", "NormAnnularNoll 0.3") to:
C    typeL        -- the corresponding ZernType_* index (1..mZernType)
C    found        -- .TRUE. if the typename was recognized
C    annuRatio    -- annular pupil ratio (only meaningful when typeL == 9)
C    hasRatio     -- .TRUE. if a ratio was successfully extracted
C
C  Match policy:
C    - Norm-prefixed names (NormANSI, NormBornWolf, NormFringe, NormHex,
C      NormNoll, NormAnnularNoll) match by the first 7 characters.
C    - Un-normalized names (ANSI, BornWolf, Fringe, Noll, ExtFringe)
C      match by the first 3 characters; this preserves backward
C      compatibility with shortened forms like "Bor" for BornWolf.
C    - Earlier MACOS versions silently rewrote "Noll" -> "NormNoll" as
C      a "backward compatibility" shim.  That shim is removed: "Noll"
C      now correctly resolves to ZernType_Noll (index 10).
C
C  Centralized so all parser sites (msmacosio.inc, macosio.F) share one
C  source of truth -- a previous bug where the rewrite shim existed at
C  one site only is what motivated this consolidation.
C***********************************************************************

        Subroutine ParseZernType(value, typeL, found, annuRatio, hasRatio)
        Implicit none

        Character(*), intent(in)  :: value
        Integer,      intent(out) :: typeL
        Logical,      intent(out) :: found
        Real*8,       intent(out) :: annuRatio
        Logical,      intent(out) :: hasRatio

        Integer :: i, k, val_len
        Logical :: LCMP        ! external

        typeL     = 0
        found     = .FALSE.
        annuRatio = 0d0
        hasRatio  = .FALSE.

        If (LCMP(value, 'Norm', 4)) Then
          ! Norm-prefixed types: NormANSI(4), NormBornWolf(5),
          ! NormFringe(6), NormHex(7), NormNoll(8), NormAnnularNoll(9)
          DO i = 4, 9
            IF (LCMP(value, ZernTypeNameL(i), 7)) THEN
              typeL = i
              found = .TRUE.
            END IF
          END DO
          ! For NormAnnularNoll, look for an additional whitespace-
          ! separated numeric token after the typename and parse it.
          If (typeL == 9) Then
            val_len = LEN(value)
            ! Skip past the typename
            i = 1
            DO WHILE (i <= val_len)
              IF (value(i:i) == ' ') EXIT
              i = i + 1
            END DO
            ! Skip whitespace
            DO WHILE (i <= val_len)
              IF (value(i:i) /= ' ') EXIT
              i = i + 1
            END DO
            ! Try to parse a real number from the remainder
            If (i < val_len) Then
              READ(value(i:), *, IOSTAT=k) annuRatio
              hasRatio = (k == 0)
            End If
          End If
        Else
          ! Un-normalized types.  Indices 1..3 (ANSI, BornWolf,
          ! Fringe) plus the later additions Noll(10) and ExtFringe(11).
          ! 3-character prefix is unambiguous: no other typename starts
          ! with "ANS", "Bor", "Fri", "Nol", or "Ext".
          DO i = 1, 3
            IF (LCMP(value, ZernTypeNameL(i), 3)) THEN
              typeL = i
              found = .TRUE.
            END IF
          END DO
          DO i = 10, 11
            IF (LCMP(value, ZernTypeNameL(i), 3)) THEN
              typeL = i
              found = .TRUE.
            END IF
          END DO
          ! Legacy alias: "Malacara" was the pre-modernization name
          ! for ANSI (index 1).  Old JPL Rx files still
          ! use it; honour the user's intent rather than silently
          ! dropping the field and falling back to ZernTypeL=0.
          IF (.NOT.found .AND. LCMP(value, 'Mala', 4)) THEN
            typeL = 1
            found = .TRUE.
          END IF
        End If

        End Subroutine ParseZernType

C***********************************************************************

        Subroutine elt_mod_init

        ! Deallocate all dynamic arrays if already allocated.
        If (allocated(LRayToElt)) Then
          deallocate(LRayToElt,LRayOK,LEltOK,LRayPass,LEltObs,
     &               ifComputeDZ,LGlass,LTaperedElt,LFastRayTrace,
     &               LApTypeDefined,LApVecDefined,
     &               ifGridDataDefined,
     &               GridFilePending,GridFileBin,
     &               ifMon,ifFF,ifGridTerm,
     &               stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(logical): deallocate failed!')

          deallocate(EltName,GlassElt,GlassName,
     &               RayFailMsg,stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(char): deallocate failed!')

          deallocate(RayStatus,RayFailElt,stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(RayStat): deallocate failed!')

          deallocate(SegCoord,RayToSegMap,RayID,RayWFMap,
     &               ObsType,IWORK,LnkElt,EltMRP, LEltMRP_FLG,
     &               EltType,nECoord,EltToSegMap,SegToEltMap,
     &               PropType,ApType,nObs,iEltToiWF,SrfType,
     &               LUpdateTElt_FLG,nObsRx,ObsTypeRx,ifAmplSrf,
     &               ifOPDModGrid,passRayIds,metBeamFlg,
     &               stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(int): deallocate failed!')

          deallocate(OpType, iEltToIntSrf, nDP,LensArrayType,
     &               nCoatElt,UDSrfType,EltID,ZernaType,
     &               ZernTypeL,FFZernTypeL,MonZernTypeL,
     &               iEltToGridSrf,nGridMat,GridSrfOrder,
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
     &               ZernCent,ZernXGrid,ZernYGrid,ZernRad,
     &               pMon,xMon,yMon,zMon,pMonRx,xMonRx,yMonRx,
     &               pData,xData,yData,zData,pDataRx,xDataRx,
     &               yDataRx,zDataRx,lData,
     &               zMonRx, h1HOE,h2HOE, CoatIndxElt, CoatThkElt,
     &               PolyApVtx,PolyApEqn,PolyObsVtx,PolyObsEqn,
     &               PolyObsVtx3DPending,nPolyObsVtx3DPending,
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

          ! Deallocate FreeForm arrays
          deallocate(lFF,pFF,xFF,yFF,zFF,FFZernCoef,FFCoef,
     &               MonZernCoef,stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &             ('elt_mod_init(FreeForm real*8): deallocate failed!')

          deallocate(WFElt, JmatElt, ObsMat, RayE, WFbuff, OTFMat,
     &               stat=m_err_elt)
          if (m_err_elt /= 0)
     &      call macos_memory_failure
     &               ('elt_mod_init(complex*16): deallocate failed!')
        End If


        ! Allocations

        ! Logical arrays
        allocate(LRayToElt(0:mElt,bRay),LRayOK(mRay),LEltOK(mElt),
     &           LRayPass(mRay),LEltObs(mElt),ifComputeDZ(mIntSrf),
     &           LGlass(mElt),LTaperedElt(mElt),LEltMRP_FLG(0:mElt),
     &           LApTypeDefined(mElt),LApVecDefined(mElt),
     &           LUpdateTElt_FLG(0:mElt),LFastRayTrace(mRay),
     &           ifGridDataDefined(mGridSrf),
     &           GridFilePending(mElt),GridFileBin(mElt),
     &           ifMon(mElt),ifFF(mElt),ifGridTerm(mElt),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init: allocate(logical) failed!')

        ! Character arrays
        allocate(EltName(0:mElt), GlassElt(mElt), GlassName(mGlass),
     &           EltMRP(0:mElt),RayFailMsg(mRay),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init: allocate(char) failed!')

!-->  Per-ray status arrays
        allocate(RayStatus(mRay),RayFailElt(mRay),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init: allocate(RayStat) failed!')
!<--

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

        allocate(iEltToIntSrf(mElt),nDP(mIntSrf),
     &           LensArrayType(mElt),nCoatElt(mElt),
     &           UDSrfType(0:mElt),EltID(0:mElt),
     &           ZernaType(mElt),ZernTypeL(mElt),FFZernTypeL(mElt),MonZernTypeL(mElt),
     &           iEltToGridSrf(mElt),nGridMat(mElt),
     &           GridSrfOrder(mElt),
     &           nRayToSeg(mRay),NScnt(mElt),nAmplMat(mElt),
     &           iEltToAmplSrf(mElt),iEltToMetSrf(mElt),
     &           nMetPos(mElt),tMetSrf(mElt),ntMetPos(mElt),
     &           EdgeSensor(mElt),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(int-2): allocate failed!')

        ! SREAL arrays
        allocate(dxidxim1(7,7,mElt), dxidui(7,6,mElt),
     &           MWFFT(mdttl,mdttl), XFP(mdttl), YFP(mdttl),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(SREAL): allocate failed!')

        mAmplMat = mGridMat
        mAmplSrf = mGridSrf

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
     &           PolyObsVtx3DPending(3,mPolySide,mObs,mElt),
     &           nPolyObsVtx3DPending(mObs,mElt),
     &           vptElt(3,mElt),vptRx(3,mElt),RptRx(3,mElt),
     &           psiRx(3,mElt),DoePhProf(6,mDOE),stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(real*8-0.5): allocate failed!')

        allocate(RptElt(3,mElt), psiElt(3,mElt), RayPos(3,mRay),
     &           RayDir(3,mRay), ApVec(6,mElt), xObs(3,mElt),
     &           AsphCoef(mAsphCoef,mElt),MonCoef(mMonCoef,mElt),
     &           ZernCoef(mZernModes,mElt),ZernCent(3,mElt),
     &           ZernRad(mElt),
     &           ZernXGrid(3,mElt),ZernYGrid(3,mElt),
     &           pMon(3,mElt),pMonRx(3,mElt),
     &           pData(3,mElt),pDataRx(3,mElt),
     &           ZernCoefRx(mZernModes,mElt),
     &           LxObsDef(mElt),tMetSrfPos(mMetPos,mElt),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(real*8-1): allocate failed!')

        allocate(xMon(3,mElt),yMon(3,mElt),zMon(3,mElt),
     &           xMonRx(3,mElt),yMonRx(3,mElt),zMonRx(3,mElt),
     &           xData(3,mElt),yData(3,mElt),zData(3,mElt),
     &           xDataRx(3,mElt),yDataRx(3,mElt),zDataRx(3,mElt),
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

        allocate(CumRayL(mRay),CumCRL(mElt),zElt(0:mElt),KcElt(mElt),
     &        KrElt(mElt),lMon(mElt), lData(mElt), LensArrayWidth(mElt),
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

        ! FreeForm arrays
        allocate(lFF(mElt),
     &           pFF(3,mElt),xFF(3,mElt),yFF(3,mElt),zFF(3,mElt),
     &           FFZernCoef(mFFCoef,mElt),FFCoef(mFFCoef,mElt),
     &           MonZernCoef(mMonCoef,mElt),
     &           stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(FreeForm real*8): allocate failed!')

        ! Complex*16 arrays
        allocate(WFElt(mdttl,mdttl,mWF), JmatElt(2,2,mElt),
     &           ObsMat(mObsMat,mObsMat,mCplxObs),
     &           RayE(3,mRay), WFbuff(mdttl,mdttl),
     &           OTFMat(mdttl,mdttl), stat=m_err_elt)
        if (m_err_elt /= 0)
     &    call macos_memory_failure(
     &          'elt_mod_init(complex*16): allocate failed!')

        End Subroutine elt_mod_init

C***********************************************************************

        SUBROUTINE elt_mod_init_vars()
          USE Kinds
          IMPLICIT NONE

          LRayToElt(:,:)       = .FALSE.
          LRayOK(:)            = .FALSE.
          LEltOK(:)            = .FALSE.
          LRayPass(:)          = .FALSE.

!-->  Per-ray status reset
          RayStatus(:)         = RayStat_OK
          RayFailElt(:)        = 0
          RayFailMsg(:)        = ' '
!<--
!-->  Rx comment round-trip reset
          nRxComment           = 0
          RxCommentCtx         = 0
          LRxCommentCapture    = .FALSE.
!<--
          LEltObs(:)           = .FALSE.
          ifComputeDZ(:)       = .FALSE.
          LGlass(:)            = .FALSE.
          LTaperedElt(:)       = .FALSE.
          LApTypeDefined(:)    = .FALSE.
          LApVecDefined(:)     = .FALSE.
          LEltMRP_FLG(:)       = .FALSE.
          LUpdateTElt_FLG(:)   = .FALSE.
          LxObsDef(:)          = .FALSE.
          ifAmplSrf(:)         = .FALSE.
          ifOPDModGrid(:)      = .FALSE.
          LFastRayTrace(:)     = .FALSE.
          ifGridDataDefined(:) = .FALSE.
          GridFilePending(:)   = .FALSE.
          GridFileBin(:)       = .FALSE.

          ! FreeForm activity flags -- set at prescription-read time
          ! via SetFreeFormFlags(ie) based on SrfType and sentinels
          ifMon(:)             = .FALSE.
          ifFF(:)              = .FALSE.
          ifGridTerm(:)            = .FALSE.

          ifInMTFCalc          = .FALSE.
          BWK_FLG              = .FALSE.
          LStopAtEltAutoSet    = .FALSE.
          LSavRayL             = .FALSE.
          SrcLF_FLG            = .FALSE.    ! set when "xFrame" or "yFrame" for Source Prb. is set  (Why here?) -> macos_init.inc, macos_mod.F   <- system flag
          CalcSrcZFrame_FLG    = .FALSE.    ! set when "zFrame" is defined in Rx                    (Why here?) -> macos_init.inc, macos_mod.F
          DupElt_FLG           = .FALSE.

          EltName(:)           = ''
          GlassElt(:)          = ''
          GlassName(:)         = ''
          EltMRP(:)            = ''

          nElt                 =  0          ! total number of optical elements                (nElt <= mElt)
          nRay                 =  0
          nOutCord             =  0
          nSeg                 =  0          ! total number of segments (requires seg. source) (nSeg <= mSeg)
          spcOption            =  0

          maxAsphCoef          =  4          ! default to 4 terms  (Heritage => affects MACOS 'NEW' cmd)

          nRefRay              =  0

          nDOE                 =  0          ! Diffractive Optical Element (DOE) related arrays
          DoeId(:)             =  0
          DoeWL(:)             =  0e0_pr     ! DOE design wavelength
          DoeFl(:)             = -1e0_pr
          DoePhProf(:,:)       =  0e0_pr     ! phase profile

          EdgeSensor(:)        =  0
          EdgeSensVec(:,:,:)   = -1e0_pr

          nMetMeas             =  0          ! Metrology
          metBeamFlg(:,:,:)    =  0
          tMetSrfPos(:,:)      =  0
          iEltToMetSrf(:)      =  0
          nMetPos(:)           =  0
          tMetSrf(:)           =  0
          ntMetPos(:)          =  0
          SrfMetPos(:,:,:)     =  0e0_pr
          SrfMetMea(:,:)       =  0e0_pr
          metMeasBuf(:)        =  0e0_pr

          sysLOS(:)            =  0e0_pr
          sysLOSAng(:)         =  0e0_pr
          CurRefRayDir(:)      =  0e0_pr

          nAmplMat(:)          =  0
          AmplMat(:,:,:)       =  0e0_pr
          AmplSrfdx(:)         =  0e0_pr
          iEltToAmplSrf(:)     =  0

          SegCoord(:,:)        =  0
          RayToSegMap(:,:)     =  0
          passRayIds(:)        =  0
          RayID(:,:)           =  0
          RayWFMap(:,:)        =  0
          IWORK(:,:)           =  0
          EltGrp(:,:)          =  0
          ObsType(:,:)         =  0
          ObsTypeRx(:,:)       =  0

          EltType(:)           =  0
          LnkElt(:)            = -1
          nECoord(:)           =  0
          EltToSegMap(:)       =  0
          SegToEltMap(:)       =  0
          PropType(:)          =  0
          ApType(:)            =  0
          nObs(:)              =  0
          iEltToiWF(:)         =  0
          SrfType(:)           =  0
          OpType(:)            =  0
          iEltToIntSrf(:)      =  0
          nDP(:)               =  0
          LensArrayType(:)     =  0
          nCoatElt(:)          =  0
          UDSrfType(:)         =  0
          EltID(:)             =  0
          ZernaType(:)         =  0
          ZernTypeL(:)         =  0
          iEltToGridSrf(:)     =  0
          nGridMat(:)          =  0
          GridSrfOrder(:)      =  1
          GlassType(:)         =  0
          EltCoat(:)           =  0
          nObsRx(:)            =  0

          nRayToSeg(:)         =  0
          NSCnt(:)             =  0

          SegApType            =  0
          SegApVec(:)          =  0e0_pr

          dxidxim1(:,:,:)      = 0e0_pr
          dxidui(:,:,:)        = 0e0_pr
          MWFFT(:,:)           = 0e0_pr
          XFP(:)               = 0e0_pr
          YFP(:)               = 0e0_pr

          Wavelen              = 0e0_pr
          Flux                 = 0e0_pr
          zSource              = 0e0_pr
          Tout(:,:)            = 0e0_pr
          WavelenUM            = 0e0_pr
          WavelenConv          = 0e0_pr
          WaveBU               = 0e0_pr
          CBM                  = 0e0_pr
          CWM                  = 0e0_pr
          CWB                  = 0e0_pr
          WavelenRx            = 0e0_pr
          ToutRx(:,:)          = 0e0_pr
          BeamSizeAtElt        = 0e0_pr


          TElt(:,:,:)          =  0e0_pr
          TEltRx(:,:,:)        =  0e0_pr
          ObsVec(:,:,:)        =  0e0_pr
          XYZSRF(:,:,:)        =  0e0_pr
          DZSRF(:,:,:)         =  0e0_pr
          GridMat(:,:,:)       =  0e0_pr
          ApHexBound(:,:,:)    = -1e0_pr
          PolyApVtx(:,:,:)     =  0e0_pr
          PolyApEqn(:,:,:)     =  0e0_pr
          ObsVecRx(:,:,:)      =  0e0_pr
          ObsHexBound(:,:,:,:) =  0e0_pr
          PolyObsVtx(:,:,:,:)  =  0e0_pr
          PolyObsEqn(:,:,:,:)  =  0e0_pr
          PolyObsVtx3DPending(:,:,:,:) = 0e0_pr
          nPolyObsVtx3DPending(:,:)    = 0

          VptElt(:,:)          =  0e0_pr
          RptElt(:,:)          =  0e0_pr
          psiElt(:,:)          =  0e0_pr
          RayPos(:,:)          =  0e0_pr
          RayDir(:,:)          =  0e0_pr
          ApVec(:,:)           =  0e0_pr
          xObs(:,:)            =  0e0_pr
          AsphCoef(:,:)        =  0e0_pr
          MonCoef(:,:)         =  0e0_pr
          ZernCoefRx(:,:)      =  0e0_pr
          ZernCoef(:,:)        =  0e0_pr
          pMon(:,:)            =  0e0_pr
          xMon(:,:)            =  0e0_pr
          yMon(:,:)            =  0e0_pr
          zMon(:,:)            =  0e0_pr
          pData(:,:)           =  0e0_pr
          xData(:,:)           =  0e0_pr
          yData(:,:)           =  0e0_pr
          zData(:,:)           =  0e0_pr
          ZernCent(:,:)        =  0e0_pr
          ZernXGrid(:,:)       =  0e0_pr
          ZernYGrid(:,:)       =  0e0_pr
          h1HOE(:,:)           =  0e0_pr
          h2HOE(:,:)           =  0e0_pr
          CoatIndxElt(:,:)     =  0e0_pr
          CoatThkElt(:,:)      =  0e0_pr
          AnaCoef(:,:)         =  0e0_pr
          UDSrfCoef(:,:)       =  0e0_pr
          yObs(:,:)            =  0e0_pr
          zObs(:,:)            =  0e0_pr
          DWORK(:,:)           =  0e0_pr
          UDSrfParam(:,:)      =  0e0_pr
          vptRx(:,:)           =  0e0_pr
          psiRx(:,:)           =  0e0_pr
          RptRx(:,:)           =  0e0_pr
          GlassCoef(:,:)       =  0e0_pr
          GlassTable(:,:)      =  0e0_pr
          MTFMat(:,:)          =  0e0_pr
          MTFMat2(:,:)         =  0e0_pr
          pMonRx(:,:)          =  0e0_pr
          xMonRx(:,:)          =  0e0_pr
          yMonRx(:,:)          =  0e0_pr
          zMonRx(:,:)          =  0e0_pr
          pDataRx(:,:)         =  0e0_pr
          xDataRx(:,:)         =  0e0_pr
          yDataRx(:,:)         =  0e0_pr
          zDataRx(:,:)         =  0e0_pr
          IndRefArr(:,:)       =  0e0_pr
          ExtincArr(:,:)       =  0e0_pr
          EltCoatThk(:,:)      =  0e0_pr

          fElt(:)              =  0e0_pr
          eElt(:)              =  0e0_pr
          ZernRad(:)           =  0e0_pr
          RayL(:)              =  0e0_pr
          dxElt(:)             =  0e0_pr
          IndRef(:)            =  0e0_pr
          ArrIndRef(:,:)       =  0e0_pr
          CumRayL(:)           =  0e0_pr
          CumRayLsav(:)        =  0e0_pr
          CumCRL(:)            =  0e0_pr
          zElt(:)              =  0e0_pr
          KcElt(:)             =  0e0_pr
          KrElt(:)             =  0e0_pr
          lMon(:)              =  0e0_pr
          lData(:)             =  0e0_pr
          LensArrayWidth(:)    =  0e0_pr
          OrderHOE(:)          =  0e0_pr
          WaveHOE(:)           =  0e0_pr
          PinHole(:)           =  0e0_pr
          Extinc(:)            =  0e0_pr
          RuleWidth(:)         =  0e0_pr
          GridSrfdx(:)         =  0e0_pr
          KrRx(:)              =  0e0_pr
          KcRx(:)              =  0e0_pr
          ZernAnnuRatio(:)     =  0e0_pr

          Fn                   =  0e0_pr  ! Fresnel number
          OPDdx                =  0e0_pr  ! exit pupil OPD grid spacing
          OPDdy                =  0e0_pr
          RayDirBuff(:,:)      =  0e0_pr
          RefRayDir(:,:)       =  0e0_pr
          rayHt                =  0e0_pr
          costh                =  0e0_pr
          efl                  =  0e0_pr
          focalRatio           =  0e0_pr
          plateScale           =  0e0_pr
          ChfRayDirAtElt(:)    =  0e0_pr
          ChfRayPosAtElt(:)    =  0e0_pr

          ZnNmFac_m(:)         =  0e0_pr

          ! FreeForm arrays
          lFF(:)               =  0e0_pr
          pFF(:,:)             =  0e0_pr
          xFF(:,:)             =  0e0_pr
          yFF(:,:)             =  0e0_pr
          zFF(:,:)             =  0e0_pr
          FFZernCoef(:,:)      =  0e0_pr
          FFCoef(:,:)          =  0e0_pr
          FFZernTypeL(:)       =  0
          MonZernCoef(:,:)     =  0e0_pr
          MonZernTypeL(:)      =  0

          WFElt(:,:,:)         =  0e0_prc
          JmatElt(:,:,:)       =  0e0_prc
          ObsMat(:,:,:)        =  0e0_prc
          RayE  (:,:)          =  0e0_prc
          WFbuff(:,:)          =  0e0_prc
          OTFMat(:,:)          =  0e0_prc

          m_iCenterPMSeg       =  0
          m_iFirstPMSeg        =  0
          m_iLastPMSeg         =  0
          m_iPM                =  0
          m_iSM                =  0
          m_iTM                =  0
          m_iFSM               =  0
          m_iPickoff           =  0
          m_iFM(:)             =  0
          m_iExitPupil         =  0
          m_iFocalPlane        =  0
          m_nLensGrp           =  0
          m_LensGrp(:,:)       =  0
          m_LensGrpElts(:)     =  0
          refRay(:)            =  0
          m_err_elt            =  0

          ! PrvGlass
          nPrvGlass          =   0      ! actual private glasses defined in Rx
          PrvGlassCatalog(:) = PrvGlass(glsName =  '',
     &                                  nEntry  =   0,
     &                                  IndRef  =  0e0_pr,
     &                                  WavLen  =  0e0_pr)

        END SUBROUTINE elt_mod_init_vars

C***********************************************************************

        SUBROUTINE EltConvert(EltType,EltID,SrfType)
          IMPLICIT NONE
          INTEGER, intent(in) :: EltType   ! old Encoding ID
          INTEGER, intent(out):: EltID
          INTEGER, intent(out):: SrfType
          ! ----------------------------------
          EltID   = EltType2EltIDSrfType(EltType,1)  ! => EltTypeName
          SrfType = EltType2EltIDSrfType(EltType,2)  ! => SrfTypeName

        END SUBROUTINE EltConvert

C***********************************************************************
C  SetFreeFormFlags
C
C  Called at prescription-read time (after SrfType, lMon, lFF, nGridMat
C  are all set for element ie) to derive the activity flags ifMon, ifFF,
C  ifGridTerm.  These flags are used at ray-trace time to avoid sentinel
C  checks inside the Brent solver loop.
C
C  Sentinels:
C    Mon active  : lMon(ie) > 0
C    FF  active  : lFF(ie)  > 0
C    Grid active : nGridMat(ie) > 0
C
C  All three are only meaningful when SrfType(ie) == SrfType_FreeForm.
C***********************************************************************

        SUBROUTINE SetFreeFormFlags(ie)
          IMPLICIT NONE
          INTEGER, intent(in) :: ie

          IF (SrfType(ie) .EQ. SrfType_FreeForm) THEN
            ifMon(ie)  = (lMon(ie)    .GT. 0d0)
            ifFF(ie)   = (lFF(ie)     .GT. 0d0)
            ifGridTerm(ie) = (nGridMat(ie) .GT. 0)
          ELSE
            ifMon(ie)  = .FALSE.
            ifFF(ie)   = .FALSE.
            ifGridTerm(ie) = .FALSE.
          END IF

        END SUBROUTINE SetFreeFormFlags

C***********************************************************************
C  SetRayFail -- record the first failure for a ray
C
C  Only stores the first failure (subsequent calls for the same ray
C  are ignored), so the root cause is preserved.  Returns .TRUE. iff
C  this call actually recorded a status (i.e. the ray was previously
C  OK).  Callers gate per-failure tallies (nBadRays) on the return so
C  the same ray isn't counted more than once across propagation
C  passes that re-observe the obscured/missed state.
C***********************************************************************

!-->
        LOGICAL FUNCTION SetRayFail(iRay, iStat, iElt, cMsg)
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: iRay, iStat, iElt
          CHARACTER(len=*), INTENT(IN) :: cMsg
          IF (RayStatus(iRay) .EQ. RayStat_OK) THEN
            RayStatus(iRay)  = iStat
            RayFailElt(iRay) = iElt
            RayFailMsg(iRay) = cMsg
            SetRayFail = .TRUE.
          ELSE
            SetRayFail = .FALSE.
          END IF
        END FUNCTION SetRayFail
!<--

C***********************************************************************

        End Module elt_mod