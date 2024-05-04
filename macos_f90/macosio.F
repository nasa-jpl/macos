C***********************************************************************
C	Begin file macosio.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module macosio_mod
	use math_mod

	Implicit none

	Logical :: RxStopSet,EltStopSet,PrtPrinRots,RxNoStopSet_Flg
	REAL*8 :: StopOffset(2)

	Interface
          SUBROUTINE GridInit(GridFile,jGridSrf,nGridMatS,
     &                        ifGridInit,binFil)
          IMPLICIT NONE
          CHARACTER(len=*) :: GridFile
          LOGICAL ifGridInit
          Logical :: binFil
          INTEGER jGridSrf,nGridMatS
          End SUBROUTINE GridInit
        End Interface

	Contains

	SUBROUTINE MBFile6(npts,ifOld,ifFid,LOAD_SUCCESS)
	use Constants, only: EPS
	use param_mod
        use src_mod
        use elt_mod
        use surfsub, only: GridInit, AMPLINIT
        use cfiles_mod
        use macos_mod, only: ifStopSet,ifObjImgSet,ObjIntFn,zernUsrOpt,
     &                       zcoType,ZernCircular,ZernAnnular,ZernHex,
     &                       cUnit,ifRegrid,ifPGColor
	use lohpars_mod, only : cmdArg,RX_Flg
	use traceutil_mod
	use loh_I1
	use loh_I2
	use loh_I3
	use dopt_mod
        use macos_IO, only: SetPrtToFileID
! --------------------------------------
        IMPLICIT NONE

        INTEGER, PARAMETER :: RX_FID = 77
        LOGICAL :: isRxOpen

	CHARACTER*8 CNUM
	CHARACTER(len=MacosCharLen) :: ANS,TypeString
	LOGICAL exist,ifOld,LCMP,IHSTPM,IHSTSM,ifGridInit,
     &          ifAmplInit,ifFid,LOAD_SUCCESS
	INTEGER i,j,k,l,g,iaxis,iElt,npts,iopt,iCoord,IDEF(7),ICLEN,
     &          iSeg,iObs,jGridSrf,eltGrpStrt,eltGrpEnd,grpRange(20),
     &          nValsPerRow,fid,status,jAmplSrf,jMetSrf,
     &          dupElt,nTok,GetNumToken
	REAL*8 invec6(6),invec7(7),DDEF(7)
	REAL*8 fseg,eseg,psiseg(3),Vseg(3),S1,S2

C       Setup for new input scheme

	CHARACTER*4 IND(3)
        CHARACTER(len=MacosCharLen) :: STRING
        CHARACTER(len=MacosVarNamLen) :: VAR_NAM
        CHARACTER(len=MacosValLen) :: VALUE
!	CHARACTER*160 VALUE
	LOGICAL ARR_FLG,opt_npts_flg
	LOGICAL ENTR_OK, ELEM_OK, EXIT_OK, ApInputErr
	INTEGER ians, nVertex,StrLen

	REAL*8 D1(9)

        INTEGER :: nAsphCoef              ! Set to Default when reading in values to make it backwards compatible
        INTEGER :: nZernCoef              ! Ditto
        INTEGER :: ZernModes(mZernModes)  ! need a temp space for Zernike Single-Index mapping (hash-tbl)


c***loh<
C	File setup and initial conditions read

 102	FORMAT(1P,10x,D17.9)
 103	FORMAT(1P,10x,3D17.9)
 104	FORMAT(1P,10x,i4)
 105	FORMAT(1P,10x,i4)
 106	FORMAT(1P,10x,a8)
 107	FORMAT(' Enter coordinate ',i1,' 6-vector (Qx,Qy,Qz,Dx,Dy,Dz):')
 108	FORMAT(1P,5(10x,6D17.9/),10x,6D17.9)
 109	FORMAT(1P,6(10x,7D17.9/),10x,7D17.9)
 110	FORMAT(' Enter coordinate',i4,' 7-vector (x,y,z,x,y,z,dL):')
 111	FORMAT(1P,'nOutCord= ',i4/'    Tout= ',7D17.9,6(/10x,7D17.9))
 112	FORMAT(1P,10x,6D17.9)
 113	FORMAT(1P,'    TElt= ',6D17.9)
 114	FORMAT(' Enter Element',i4,' Data:')
 115	FORMAT(' Enter 3 Segment ',i4,' Hex Coords:')
 116	FORMAT(10x,i5,5x,i5,5x,i5)
 118	FORMAT(1P,10x,e17.10)
 120	FORMAT(1P,10x,4D17.9)
 121	FORMAT('  Too many grid points. Resetting npts to',i5)
 122	FORMAT(1P,A10,D17.9)

 123	FORMAT(A)
 124	FORMAT(A,' must have non-zero magnitude')
 126	FORMAT(1P,' PinHole= ',D17.9)
 128	FORMAT(1P,'   nCoat= ',i4/'CoatIndx= ',5D17.9/3(10X,5D17.9))
 129	FORMAT(1P,'CoatThk= ',5D17.9/3(10X,5D17.9))
 130	FORMAT(1P,' **Invalid Element Type, Element=',A)
 131	FORMAT(1P,' **Invalid Element/Surface Combination '/
     &  'Element=',A,',Surface=',A)
 132	FORMAT(1P,' **Invalid Surface Type, Surface=',A)
 133	FORMAT(1P,' **Invalid Zernike Type, Surface=',A)
 140	FORMAT(' +--------------------------------------------+'/
     &	' | Element types:                             |'/
     &	' |   Reflector      NSReflector    Segment    |'/
     &	' |   Refractor      NSRefractor    LensArray  |'/
     &	' |   FocalPlane     Reference      Return     |'/
     &	' |   HOE            Grating        TrGrating  |'/
     &	' |   RfPolarizer    TrPolarizer    Obscuring  |'/
     &	' +--------------------------------------------+')
 141	FORMAT(' +--------------------------------------------+'/
     &	' | Surface types:                             |'/
     &	' |   Flat           Conic          Aspheric   |'/
     &	' |   Anamorphic     Zernike        Monomial   |'/
     &	' |   Interpolated   UserDefined    Grid data  |'/
     &	' |   Toric          InfluenceFcn              |'/
     &	' +--------------------------------------------+')
 142	FORMAT(' +--------------------------------------------+'/
     &	' | Propagation types:                         |'/
     &	' |   Geometric      GeomUpdate     FarField   |'/
     &	' |   NFSpherical    NFS1surf                  |'/
     &	' |   NFPlane        NFP1surf                  |'/
     &	' |   NF1            NF2                       |'/
     &	' |   SpatialFilter  SF1surf                   |'/
     &	' +--------------------------------------------+')
 143	FORMAT(' +--------------------------------------------+'/
     &	' | Source ray grid types:                     |'/
     &	' |   Circular       Square         Hex        |'/
     &	' |   Pie            Flower                    |'/
     &	' +--------------------------------------------+')
 144	FORMAT(' +--------------------------------------------+'/
     &	' | Obscuration types:                         |'/
     &	' |   Circle         NegCircle                 |'/
     &	' |   Rectangle      NegRectangle              |'/
     &	' |   Annulus        NegAnnulus                |'/
     &  ' |   Triangle       NegTriangle               |'/
     &  ' |   Ellipse        NegEllipse                |'/
     &	' +--------------------------------------------+')
 145	FORMAT(' +--------------------------------------------+'/
     &	' | Aperture types:                            |'/
     &	' |   Circular       Rectangular               |'/
     &	' +--------------------------------------------+')
 146	FORMAT(' +--------------------------------------------+'/
     &	' | Zernike types:                             |'/
     &	' |   Malacara       Noll           Fringe     |'/
     &	' +--------------------------------------------+')

        ! Flags
        LOAD_SUCCESS = .FALSE.        ! so far, Rx is not loaded

        CALL reinitialise_variables() ! Reset ALL: see "call macos_init_all()" in macos_init.F
        CALL LOHIN1()                 ! def in iosub.inc
        CALL LOHIN2()                 ! def in iosub.inc

        ! Initialize constants

        DDEF(1:7)=0D0
        IDEF(1:7)=0D0
        CWM=1d0; CBM=1d0
        opt_npts_flg=.false.

C	Get file name

 43	CONTINUE
	If (.not.ifFid) Then
	  CALL CACCEPT(filnam,' ','Enter file name: ')
	Else
	  if (.not.RX_Flg) then
	    Call LSFILES(0,cmdArg,0)
            RX_Flg=.TRUE.
	  end if
	  CALL IACCEPT_S(fid,1,1,'Enter file id: ')
	  filnam(:)=' '
	  CALL GETRXFN(fid,filnam,status)
	  if (status==0) then
	    WRITE(*,*)' ** Invalid file id; can not load file'
	    WRITE(*,*)' '
	    LOAD_SUCCESS=.false.; return
	  end if
	End If

	L=ICLEN(filnam)
	IF ((L.EQ.1).AND.((filnam.EQ.'q').OR.(filnam.EQ.'Q'))) RETURN
	infil=filnam
	outfil=filnam
	cmatfil=filnam
	plotfil=filnam
	mfil=filnam
	pertfil=filnam

        ! Updated so will accept '.in' file extension -jzlou
	if (.not.LCMP(infil(l-2:l),'.in',3)) then
	  infil(l+1:l+4)='.in '
        end if
	outfil(l+1:l+5)='.out '
	cmatfil(l+1:l+6)='.cmat '
	plotfil(l+1:l+6)='.plot '
	mfil(l+1:l+3)='.m '
	pertfil(l+1:l+6)='.pert '
	L=ICLEN(infil)
	INQUIRE (file=infil,exist=exist)
	IF (exist) THEN
	  IF (ifOld) THEN
	    WRITE (*,*) ' Input file ',infil(1:L),' being loaded.'
	  ELSE
	    WRITE (*,*) ' Input file ',infil(1:L),
     &	    ' already exists. Please use a new name or use OLD.'
	    GO TO 43
	  END IF
	  iopt=1
	ELSE
	  IF (ifOld) THEN
	    WRITE (*,*) ' Input file ',infil(1:L),
     &	    ' does not exist. Check name/path or use NEW.'
	    GO TO 43
	  END IF
	  iopt=0
	END IF

	! Comment:  -jzlou
	! When iopt=1 here, read in an existing Rx file;
        ! when iopt=0 here, write system data to a new Rx.

 29	CONTINUE
	EltName(0)  = 'InputRay'
	PropType(0) = 1
	zElt(0)     = 1d22

C
C	Read existing input file
C
	IF (IOPT.EQ.1) THEN
	  OPEN (UNIT=RX_FID, FILE=INFIL, STATUS='OLD', ACTION='READ')
	  REWIND(RX_FID)

#include "msmacosio.inc"

   70     CONTINUE
	  If (.not.LOAD_SUCCESS) GOTO 99

	  ! Check ApType and ApVec consistency for all elements
	  ApInputErr=.false.
	  Do iElt=1,nElt
	    if ((LApTypeDefined(iElt) .and.
     &           (.not.LApVecDefined(iElt))) .OR.
     &          ((.not.LApTypeDefined(iElt)) .and.
     &          LApVecDefined(iElt))) then
	      ApInputErr=.true.
	      WRITE(*,*)
     &       'ApType and ApVec are inconsistent for element',iElt
	    end if
	  End Do
	  if (ApInputErr) then
	    WRITE(*,*)'Input file not properly loaded'
	    GO TO 99
	  end if

	  CALL LOHIN3()

	  READ(VALUE,*) nOutCord
	  nOutCord_FLG=.TRUE.

	  READ(RX_FID, 123, END=99) STRING
	  CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,RX_FID)

	  IF (LCMP(VAR_NAM,'Tout',4)) THEN
	    READ(VALUE,*)(Tout(1,j),j=1,7)
	    IF (nOutCord .GT. 1) THEN
	      DO k=2,nOutCord
	 	READ(RX_FID, 123, END=99)STRING
	        READ(STRING,*)(Tout(k,j),j=1,7)
              END DO
	    END IF
	    Tout_FLG=.TRUE.
	  END IF

	  CALL ChkDf3(EXIT_OK)

	  IF (EXIT_OK) THEN
	    GO TO 80
	  ELSE
	    WRITE(*,*)'  Input file not properly loaded'
	    GO TO 99
	  END IF

   80 	  CONTINUE

	  CALL SumData(GridType,npts)

	  IF (nGridpts.GT.mpts) THEN
#if 0
	    WRITE(*,*)'Warning: nGridpts (Rx) =',nGridpts,
     &        ', mpts (macos_param) =',mpts
	    CALL CACCEPT(ANS,'Y','Use nGridpts value for grid size? ')
            If (LCMP(ANS,'Y',1)) Then
	      WRITE(*,*)'Grid size set to ',nGridpts
	      mpts=nGridpts
	      mRay=mpts**2+1
	    Else
#endif
	    WRITE(*,121) mpts
            !print*,'nGridpts, mpts=', nGridpts, mpts
	    nGridpts=mpts
	    npts=nGridpts-1
	    if (.not.opt_npts_flg) then
	      opt_npts=nGridpts/2-1
	      opt_npts_flg=.true.
	    end if
	  ELSE
	    WRITE(*,*)' Grid size set to ',nGridpts
	  END IF

! ToDo: this will cause indexing faults
          IF ((GridType.GE.3).AND.(nSeg.GT.mElt)) THEN
	    WRITE(*,*) ' Too many segments.'
	    WRITE(*,*) '  MODIFY or use another version. '
	  END IF

C	Create new input file

	ELSE
	  WRITE(*,*) 'File '//infil(1:ICLEN(infil))//' being created.'
	  OPEN (1,FILE=infil,STATUS='NEW')

C	Specify input ray grid parameters

	  CALL DZERO(ChfRayDir,3)
	  ChfRayDir(3)=1.0
	  CALL DACCEPT(ChfRayDir,ChfRayDir,3,
     &	  'Enter chief ray input direction vector (x,y,z) ')
	  CALL DUNITIZE(ChfRayDir)
	  CALL DACCEPT(ChfRayPos,DDEF,3,
     &	  'Enter chief ray location vector (x,y,z) ')
	  CALL DACCEPT_S(zSource,1d22,1,
     &	  'Enter distance from source ')
	  CALL DACCEPT_S(IndRef(0),1d0,1,
     &	  'Enter chief ray index of refraction ')
	  CALL DACCEPT_S(Extinc(0),0d0,1,
     &	  'Enter chief ray extinction coefficient ')
	  CALL CACCEPT(BaseUnits,' ',
     &	  'Enter base units (m, cm, mm, in, none) :')
	  IF (LCMP(BaseUnits,'cm',2)) THEN
	    CBM=1d-2
	  ELSE IF (LCMP(BaseUnits,'mm',2)) THEN
	    CBM=1d-3
	  ELSE IF (LCMP(BaseUnits,'m',1)) THEN
	    CBM=1d0
	  ELSE IF (LCMP(WaveUnits,'in',2)) THEN
	    CBM=2.54d-2
	  ELSE IF (LCMP(WaveUnits,'none',2)) THEN
	    CBM=1d0
	  ELSE
	    WRITE(*,*)' Unknown units, unitless operation assumed'
	    BaseUnits='none'
	    WaveUnits='none'
	    CBM=1d0
	  END IF
	  IF (.NOT.LCMP(BaseUnits,'none',2))
     &	  CALL CACCEPT(WaveUnits,BaseUnits,
     &	  'Enter wavelength units (m, cm, mm, um, nm, A, in, none) :')
	  IF (LCMP(WaveUnits,'cm',2)) THEN
	    DDEF(1)=6.328d-5
	    CWM=1d-2
	  ELSE IF (LCMP(WaveUnits,'mm',2)) THEN
	    DDEF(1)=6.328d-4
	    CWM=1d-3
	  ELSE IF (LCMP(WaveUnits,'um',2)) THEN
	    DDEF(1)=6.328d-1
	    CWM=1d-6
	  ELSE IF (LCMP(WaveUnits,'nm',2)) THEN
	    DDEF(1)=6.328d2
	    CWM=1d-9
	  ELSE IF (LCMP(WaveUnits,'m',1)) THEN
	    DDEF(1)=6.328d-7
	    CWM=1d0
	  ELSE IF (LCMP(WaveUnits,'A',1)) THEN
	    DDEF(1)=6.328d3
	    CWM=1d-10
	  ELSE IF (LCMP(WaveUnits,'in',2)) THEN
	    DDEF(1)=2.4913386d-5
	    CWM=2.54d-2
	  ELSE IF (LCMP(WaveUnits,'none',2)) THEN
	    DDEF(1)=6.328d-7
	    CWM=1d0
	  ELSE
	    WRITE(*,*)' Unknown units, unitless operation assumed'
	    WaveUnits='none'
	    DDEF(1)=6.328d-7
	    CWM=1d0
	  END IF
	  CWB=CWM/CBM
	  CALL DACCEPT_S(Wavelen,DDEF(1),1,
     &	  'Enter monochromatic wavelength ')
     	  WaveBU=CWB*Wavelen
	  CALL DACCEPT_S(Flux,1d0,1,
     &	  'Enter intensity (flux over aperture area) ')
	  WRITE(*,143)
	  CALL CACCEPT(TypeString,'Circular','Enter source grid type ')
	  DO i=1,5
	    IF (LCMP(TypeString,GridTypeName(i),2)) THEN
	      GridType=i
	      GridType_FLG=.TRUE.
	    END IF
	  END DO
	  IF (.NOT.GridType_FLG) THEN
	    WRITE(*,*)' Unknown grid type'
	    GridType=1
	  END IF
	  CALL DACCEPT_S(Aperture,1d0,1,
     &	  'Enter aperture diameter or angle ')
C	  CALL DACCEPT_S(Obscratn,DDEF(1),1,
C     &	  'Enter input central obscuration:')
	  Obscratn=0d0
	  CALL IACCEPT_S(nGridpts,mpts,1,
     &	  		 'Enter number of grid points ')
	  npts=nGridpts-1
	  IF (nGridpts.GT.mpts) THEN
	    WRITE(*,*) 'Too many grid points.'
	    WRITE(*,*) '  MODIFY or use another version. '
	  END IF
	  xGrid(1)=ChfRayDir(3)
	  xGrid(2)=ChfRayDir(1)
	  xGrid(3)=ChfRayDir(2)
	  CALL DACCEPT(xGrid,xGrid,3,
     &	  'Enter input plane x-axis vector (x,y,z) ')
	  CALL DXPROD(yGrid,ChfRayDir,xGrid)
	  CALL DACCEPT(yGrid,yGrid,3,
     &	  'Enter input plane y-axis vector (x,y,z) ')
	  IF (GridType.GE.3) THEN
	    CALL IACCEPT_S(nSeg,IDEF(1),-1,'Enter number of segments ')
	    IF (nSeg.GT.mEm1) THEN
	      WRITE(*,*) 'Too many segments.'
	      WRITE(*,*) '  MODIFY or use another version. '
	    END IF
	    IF ((GridType.EQ.3).OR.(GridType.EQ.4)) THEN
	      CALL DACCEPT_S(width,Aperture,1,'Enter segment width:')
	      CALL DACCEPT_S(gap,0d0,1,'Enter gap between segments:')
	      DO iSeg=1,nSeg
	        CALL IntToChar(Cnum,iSeg,i)
	        CALL IACCEPT(SegCoord(1:3,iSeg),IDEF,-3,
     &	        'Enter 3 segment '//Cnum(1:i)//' hex coords ')
	      END DO
	    ELSE
	      WRITE(*,*)' The flower pattern is defined by:'
	      WRITE(*,*)'     nRings = number of rings around center'
	      WRITE(*,*)'    nPetals = number of petals around center'
	      WRITE(*,*)'     radCtr = minimum radius of center segment'
	      WRITE(*,*)' Petal segment order is ascending in SegXgrid'
	      WRITE(*,*)' coordinates.'
	      CALL IACCEPT_S(nRings,1,1,'Enter number of rings')
	      CALL IACCEPT(nPetals(1:nRings),IDEF,-nRings,
     &                     'Enter number of petals')
	      CALL DACCEPT_S(radCtr,(Aperture/3d0),1,
     &	      'Enter center segment radius:')
	      DO iSeg=1,nSeg
	        CALL IntToChar(Cnum,iSeg,i)
	        CALL IACCEPT_S(SegCoord(1,iSeg),IDEF(1),1,
     &	        'Enter 1 petal coord for segment '//Cnum(1:i)//
     &	        ' (iPetal):')
	      END DO
	    END IF
	    CALL DACCEPT(SegXgrid,xGrid,3,
     &	    'Enter segment plane x-axis vector (x,y,z) ')
	    CALL DXPROD(D1,ChfRayDir,SegXgrid)
	    CALL DXPROD(SegXgrid,D1,ChfRayDir)
	    CALL DUNITIZE(SegXgrid)
	  END IF


C	Loop to specify each element

	  iSeg=0
	  CALL IACCEPT_S(nElt,IDEF(1),-1,
     &	                 'Enter number of elements ')
	  IF (nElt.GT.mElt) THEN
	    WRITE(*,*) 'Too many elements.  Reset and recompile. '
	  END IF

	  DO 15 iElt=1,nElt
 11	    CONTINUE
	    WRITE(*,114) iElt

C  Assign EltID

	    CALL CACCEPT(EltName(iElt),'Elt',
     &	    	   'Enter element name (no blanks) ')
	    WRITE(*,140)
	    CALL CACCEPT(TypeString,' ','Enter element type :')
	    DO i=1,15
	      IF (LCMP(TypeString,EltTypeName(i),6)) THEN
		EltID(iElt)=i
		EltID_FLG=.TRUE.
	      END IF
	    END DO
	    IF (.NOT.EltID_FLG) THEN
	      WRITE(*,130) TypeString
	    END IF

C  Assign SrfType

	    WRITE(*,141)
	    CALL CACCEPT(TypeString,' ','Enter surface type :')
	    DO i=1,10
	      IF (LCMP(TypeString,SrfTypeName(i),3)) THEN
	        SrfType(iElt)=i
	        SrfType_FLG=.TRUE.
	      END IF
	    END DO
	    IF (.NOT.SrfType_FLG) THEN
	      WRITE(*,132) TypeString
	    END IF

C  Check to make sure no Element/Surface incompatability exists

	    CALL EltSurfCompat(iElt,EltID_FLG,SrfType_FLG)
	    IF (EltID_FLG.AND.SrfType_FLG) THEN
	      EltType_FLG=.TRUE.
	    ELSE
	      WRITE(*,131) EltTypeName(EltID(iElt)),
     &                     SrfTypeName(EltID(iElt))
	      GO TO 11
	    END IF

C  Assign default values where needed

	    IF (SrfType(iElt).EQ.1) THEN
	      KcElt(iElt)=0d0
	      KrElt(iElt)=-1d22
	      fElt(iElt)=1d22
	      eElt(iElt)=0d0
	    END IF

C  Keep track of Segments...

	    IF (EltID(iElt).NE.11) THEN
	      iSeg=0
	    END IF
	    IF (EltID(iElt).EQ.11) THEN
	      iSeg=iSeg+1
	      EltToSegMap(iElt)=iSeg
	      SegToEltMap(iSeg)=iElt
	    ELSE
	      EltToSegMap(iElt)=0
	    END IF

C  Get Basic Parameters

	    IF (iSeg.LE.1) THEN
	      CALL DACCEPT(VptElt(1:3,iElt),DDEF,-3,
     &	      'Enter element vertex location (x,y,z) ')
	    ELSE
	      CALL DEQUATE(VptElt(1,iElt),VptElt(1,iElt-1),3)
	    END IF

	    CALL DACCEPT(RptElt(1:3,iElt),VptElt(1:3,iElt),3,
     &	    'Enter element rotation point (x,y,z) ')

	    IF (iSeg.LE.1) THEN
	      CALL DNEGATE(psiElt(1,iElt),ChfRayDir(1),3)
	      CALL DACCEPT(psiElt(1:3,iElt),psiElt(1:3,iElt),3,
     &	      'Enter element principal axis (x,y,z) ')
	      CALL DUNITIZE(psiElt(1,iElt))
	    ELSE
	      CALL DEQUATE(psiElt(1,iElt),psiElt(1,iElt-1),3)
	    END IF

C  Interupt Basic Parameter Set with Conic Parameter Set

	    IF (((SrfType(iElt).NE.1).AND.(SrfType(iElt).NE.6))
     &          .AND.(iSeg.LE.1)) THEN
	      CALL IACCEPT_S(ians,1,1,
     &	                     'Enter e&f or Kc&Kr (1=e&f,2=Kc&Kr)? ')
	      IF (ians .EQ. 1) THEN
		CALL DACCEPT_S(fElt(iElt),1d22,1,
     &		'Enter element focal length ')
		CALL DACCEPT_S(eElt(iElt),0d0,1,
     &	        'Enter element eccentricity ')
		KcElt(iElt)=-eElt(iElt)*eElt(iElt)
		KrElt(iElt)=-(1.0+eElt(iElt))*fElt(iElt)
	      ELSE
		CALL DACCEPT_S(KcElt(iElt),0d0,1,
     &		'Enter Conic Constant ')
		CALL DACCEPT_S(KrElt(iElt),-1d22,1,
     &	        'Enter Radius ')
		IF(KrElt(iElt).GT.0d0)KrElt(iElt)=-KrElt(iElt)
		IF (KcElt(iElt).LE.0) THEN
		  eElt(iElt) = DSQRT(-KcElt(iElt))
		  fElt(iElt) = -KrElt(iElt)/(1d0+eElt(iElt))
		ELSE
		  fElt(iElt) = 0.0
		END IF
	      END IF
	    ELSE IF (iSeg.GT.1) THEN
	      fElt(iElt)=fElt(iElt-1)
	      eElt(iElt)=eElt(iElt-1)
	      KcElt(iElt)=KcElt(iElt-1)
	      KrElt(iElt)=KrElt(iElt-1)
	    ELSE
	      eElt(iElt)=0d0
	      fElt(iElt)=1d22
	      KcElt(iElt)=0d0
	      KrElt(iElt)=-1d22
	    END IF

C  End Conic Parameter Set
C  Back to Basic Parameter Set

	    CALL CACCEPT(ANS,'NO',
     &	      'Do you want an aperture on this element?')
	    IF (LCMP(ANS,'YES',1)) THEN
	      WRITE(*,145)
	      CALL CACCEPT(TypeString,'Circular','Enter aperture type ')
	      DO i=1,2
	        IF (LCMP(TypeString,ApTypeName(i),2)) THEN
	          ApType(iElt)=i
	          ApType_FLG=.TRUE.
	        END IF
	      END DO
	      IF (.NOT.ApType_FLG) THEN
	        WRITE(*,*)' Unknown aperture type'
	        ApType(iElt)=0
	      END IF
	    ELSE
	      ApType(iElt)=0
	    END IF

C  Get Aperture Parameters

	    IF (ApType(iElt).NE.0) THEN
	      IF (ApType(iElt).EQ.1) THEN
	 	ApVec(1,iElt)=Aperture/2
	 	ApVec(2,iElt)=ChfRayPos(1)
	 	ApVec(3,iElt)=ChfRayPos(2)
	 	CALL DACCEPT(ApVec(:,iElt),ApVec(:,iElt),3,
     &	        'Enter aperture radius, x, y ')
	      ELSE IF (ApType(iElt).EQ.2) THEN
	 	ApVec(1,iElt)=ChfRayPos(1)-Aperture/2
	 	ApVec(2,iElt)=ChfRayPos(1)+Aperture/2
	 	ApVec(3,iElt)=ChfRayPos(2)-Aperture/2
	 	ApVec(4,iElt)=ChfRayPos(2)+Aperture/2
	 	CALL DACCEPT(ApVec(:,iElt),ApVec(:,iElt),4,
     &	        'Enter aperture xmin,xmax,ymin,ymax ')
	      ELSE
	 	CALL DACCEPT(ApVec(:,iElt),DDEF(:),5,
     &	        'Enter aperture xmin,xmax,ymin,ymax,radius ')
	      END IF
	    END IF

C  End Aperture Parameters
C  Back to Basic Parameter Set

	    CALL IACCEPT_S(nObs(iElt),0,1,
     &	    'Enter number of obscurations ')
	    IF (nObs(iElt).GT.0) THEN
	      IF (nObs(iElt).GT.mObs) THEN
	        WRITE(*,*)' This version limited to ',mObs,
     &	        ' obscurations'
	        nObs(iElt)=mObs
	      END IF
	    END IF

C  Get Obscuration Parameters

	    IF (nObs(iElt).NE.0) THEN
	      WRITE(*,144)
	      DO 20 iObs=1,nObs(iElt)
	        CALL CACCEPT(TypeString,'Circle',
     &	        'Enter obscuration type ')
	        DO i=-6,6
	          IF (LCMP(TypeString,ObsTypeName(i),4)) THEN
	            ObsType(iObs,iElt)=i
	            ObsType_FLG(iObs)=.TRUE.
	          END IF
	        END DO
	        IF (.NOT.ObsType_FLG(iObs)) THEN
	          WRITE(*,*)' Unknown obscuration type'
	          ObsType(iObs,iElt)=0
	        END IF

	        IF (ABS(ObsType(iObs,iElt)).EQ.1) THEN
		  CALL DACCEPT(ObsVec(:,iObs,iElt),DDEF,-3,
     &	          'Enter obscuration radius, x, y ')
		ELSE IF (ABS(ObsType(iObs,iElt)).EQ.2) THEN
		  CALL DACCEPT(ObsVec(:,iObs,iElt),DDEF,-4,
     &            'Enter obscuration xmin, xmax, ymin, ymax ')
		ELSE IF (ABS(ObsType(iObs,iElt)).EQ.3) THEN
		  CALL DACCEPT(ObsVec(:,iObs,iElt),DDEF,-4,
     &            'Enter annulus outer radius, x, y, inner radius ')
		ELSE IF (ABS(ObsType(iObs,iElt)).EQ.4) THEN
		  CALL DACCEPT(ObsVec(:,iObs,iElt),DDEF,-6,
     &	          'Enter triangle as q1x, q1y, q2x, q2y, q3x, q3y ')
	        ELSE IF (ABS(ObsType(iObs,iElt)).EQ.5) THEN
	          CALL DACCEPT(ObsVec(:,iObs,iElt),DDEF,-4,
     &            'Enter ellipse as xrad, yrad, x, y ')
	        ELSE IF (ABS(ObsType(iObs,iElt)).EQ.6) THEN
	          CALL DACCEPT(ObsVec(:,iObs,iElt),DDEF,-4,
     &            'Enter matrix px, py, dx, dy ')

	          CALL IACCEPT_S(nGridMat(iElt),mGridMat,1,
     &	          'Enter dimension of obscuration matrix ')
	          CALL IntToChar(Cnum,iElt,i)
	          CALL CACCEPT(GridFile(iElt),
     &	            filnam(1:ICLEN(filnam))//'.Obs'//Cnum(1:i),
     &              'Enter obscuration matrix file name:')
	          INQUIRE (FILE=GridFile(iElt),EXIST=EXIST)
	          CALL GridInit(GridFile(iElt),
     &                   jGridSrf,nGridMat(iElt),ifGridInit,.FALSE.)
	          ifGridDataDefined(iEltToGridSrf(iElt)) = ifGridInit
	        END IF
 20	      CONTINUE
	    END IF

C  End Obscuration Parameters
C  Get xObs if aperture or obscuration

	    IF ((ApType(iElt).NE.0).OR.(nObs(iElt).NE.0)) THEN
	      xObs(1,iElt)=psiElt(3,iElt)
	      xObs(2,iElt)=psiElt(1,iElt)
	      xObs(3,iElt)=psiElt(2,iElt)
	      CALL DACCEPT(xObs(:,iElt),xObs(:,iElt),3,
     &	        'Enter element x-axis vector (x,y,z) ')
	    END IF
C  End xObs

C  Back to Basic Parameter Set
	    CALL DACCEPT_S(zElt(iElt),fElt(iElt),1,
     &	    'Enter Fresnel propagation distance ')
	    IF ((EltID(iElt).EQ.3).OR.(EltID(iElt).EQ.8)) THEN
	      PropType(iElt)=3
	    ELSE
	      PropType(iElt)=1
	    END IF
	    WRITE(*,142)
	    CALL CACCEPT(TypeString,'Geometric',
     &	    'Enter propagation type ')
	    DO i=1,15
	      IF (LCMP(TypeString,PropTypeName(i),6)) THEN
	        PropType(iElt)=i
	        PropType_FLG=.TRUE.
	      END IF
	    END DO
	    IF (.NOT.PropType_FLG) THEN
	      WRITE(*,*)' Unknown propagation type'
	      PropType(iElt)=1
	    END IF
 	    IF ((PropType(iElt).EQ.7).OR.(PropType(iElt).EQ.8)) THEN
      	      CALL DACCEPT_S(PinHole(iElt),1d22,1,
     &	      'Enter pinhole radius ')
	    END IF

	    CALL CACCEPT(ANS,'NO',
     &	    'Do you wish to use element-based coordinates?:')
	    CALL DZERO(TElt(1,1,iElt),36)
	    IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	      CALL IACCEPT_S(nECoord(iElt),IDEF(1),-1,
     &	      'Enter number of element coordinates (up to 6) ')
	      IF (nECoord(iElt).NE.0) THEN
	        CALL CACCEPT(ANS,'YES',
     &	        'DUNITIZE coordinate vectors?:')
	        DO 12 iCoord=1,nECoord(iElt)
	          CALL IntToChar(Cnum,iCoord,i)
		  CALL DACCEPT(invec6,DDEF,-6,'Enter coordinate '
     &  	  //Cnum(1:i)//' 6-vector (Rx,Ry,Rz,Tx,Ty,Tz) ')
	          IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	            CALL DUNITIZE(invec6)
	            CALL DUNITIZE(invec6(4))
	          END IF
	          CALL DEQUATE(TElt(1,iCoord,iElt),invec6,6)
 12             CONTINUE
	      END IF
	    ELSE
	      nECoord(iElt)=-6
	      DO 13 iCoord=1,6
	        TElt(iCoord,iCoord,iElt)=1d0
 13	      CONTINUE
	    END IF
C  End Basic Parameter Set

C  Get Optical Parameters
	    IF ((EltID(iElt).NE.2).AND.(EltID(iElt).NE.3).AND.
     &         (EltID(iElt).NE.7).AND.(EltID(iElt).NE.8)) THEN
	      CALL CACCEPT(GlassElt(iElt),'FixedIndex',
     &	      'Enter glass type (FixedIndex):')
	      LGlass(iElt)=.FALSE.
	      IF (LCMP(GlassElt(iElt),'FixedIndex',3)) THEN
	        CALL DACCEPT_S(IndRef(iElt),1d0,1,
     &	        'Enter element index of refraction ')
	        IF ((EltID(iElt).EQ.1).OR.(EltID(iElt).EQ.4).OR.
     &	        (EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13).OR.
     &	        (EltID(iElt).EQ.11).OR.(EltID(iElt).EQ.12)) THEN
	          S1=1d22
	        ELSE
	          S1=0d0
	        END IF
	        CALL DACCEPT_S(Extinc(iElt),S1,1,
     &	        'Enter element extinction coefficient ')
	      ELSE
	        i=1
	        DO WHILE (((.NOT.LGlass(iElt))).AND.(i.LE.mGlass))
	          IF (LCMP(GlassElt(iElt),GlassName(i),
     &	          ICLEN(GlassElt(iElt)))) THEN
	            LGlass(iElt)=.TRUE.
	            CALL DEQUATE(GlassCoef(1,iElt),GlassTable(1,i),6)
	          END IF
     	          i=i+1
	        END DO
	        IF (.NOT.LGlass(iElt)) THEN
		  CALL DZERO(DDEF,9)
		  DDEF(1)=1.25d0
	          CALL DACCEPT(GlassCoef(:,iElt),DDEF,6,
     &	          'Unknown glass, enter Sellmeier coefficients ')
	          LGlass(iElt)=.TRUE.
	        END IF
	      END IF
	    ELSE
	      IndRef(iElt)=1d0
	      Extinc(iElt)=0d0
	    END IF

C  End Optical Parameters
C  Get Aspheric Parameters

	    IF (SrfType(iElt).EQ.3) THEN
	      If (maxAsphCoef==4) then
	        CALL DACCEPT(AsphCoef(:,iElt),DDEF,maxAsphCoef,
     &	         'Enter aspheric coefficients A,B,C,D')
	      Else
	        CALL DACCEPT(AsphCoef(:,iElt),DDEF,maxAsphCoef,
     &	         'Enter aspheric coefficients A,B,C,D,E,F,G,H,J')
	      End If
	    ELSE
	      DO I=1,maxAsphCoef
		AsphCoef(I,iElt)=0d0
	      END DO
	    END IF

C  End Aspheric Parameters
C  Get Monomial Parameters

	    CALL DZERO(DDEF,7)
	    IF (SrfType(iElt).EQ.4) THEN
	      CALL DACCEPT(MonCoef(1:6,iElt),DDEF,6,
     &	      'Enter monomial coefficients 1-6 ')
	      CALL DACCEPT(MonCoef(7:12,iElt),DDEF,6,
     &	      'Enter monomial coefficients 7-12 ')
	      CALL DACCEPT(MonCoef(13:18,iElt),DDEF,6,
     &	      'Enter monomial coefficients 13-18 ')
	      CALL DACCEPT(MonCoef(19:24,iElt),DDEF,6,
     &	      'Enter monomial coefficients 19-24 ')
	      CALL DACCEPT(MonCoef(25:30,iElt),DDEF,6,
     &	      'Enter monomial coefficients 25-30 ')
	      CALL DACCEPT(MonCoef(31:36,iElt),DDEF,6,
     &	      'Enter monomial coefficients 31-36 ')
	      DO i=37,120
	        MonCoef(i,iElt)=0d0
	      END DO
	    ELSE
	      DO I=1,120
	        MonCoef(I,iElt)=0d0
	      END DO
	    END IF

C  End Monomial Parameters
C  Get Anamorph Parameters

	    IF (SrfType(iElt).EQ.6) THEN
	      CALL DACCEPT(AnaCoef(:,iElt),DDEF,4,
     &	      'Enter anamorphic coefficients Krx,Kry,Kcx,Kcy ')
	    ELSE IF (SrfType(iElt).EQ.10) THEN
	      CALL DACCEPT(AnaCoef(:,iElt),DDEF,3,
     &	      'Enter toric coefficients Krx,Kry,Kcx,Kcy ')
	    ELSE
	      DO I=1,4
		AnaCoef(I,iElt)=0d0
	      END DO
	    END IF

C  End Anamorph Parameters

C  Get Zernike Parameters
	    IF (SrfType(iElt).EQ.8) THEN
	      WRITE(*,146)
	      CALL CACCEPT(TypeString,'Malacara',
     &	             'Enter Zernike Type :')
	      if (LCMP(TypeString,'Norm',4)) then
                ! normalized Zernike coeff
	        DO i=4,6
                  IF (LCMP(TypeString,ZernTypeNameL(i),3)) THEN
                    ZernTypeL(iElt)=i
                    ZernTypeL_FLG=.TRUE.
                  END IF
                END DO
	      else
	        ! un-normalized Zernike coeff
	        DO i=1,3
	          IF (LCMP(TypeString,ZernTypeNameL(i),3)) THEN
		    ZernTypeL(iElt)=i
		    ZernTypeL_FLG=.TRUE.
	          END IF
	        END DO
	      end if
	      IF (.NOT.ZernTypeL_FLG) THEN
	        WRITE(*,*)' Unknown Zernike type'
	        ZernTypeL(iElt)=1
	      END IF
	    ELSE
              ZernTypeL(iElt)=1
              ZernTypeL_FLG=.TRUE.
	    END IF

	    ! This part is for entering Zernike coeff in command-lind
	    ! mode. -jzlou
	    IF (SrfType(iElt).EQ.8) THEN
	      CALL DACCEPT(ZernCoef(1:6,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 1-6 ')
	      CALL DACCEPT(ZernCoef(7:12,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 7-12 ')
	      CALL DACCEPT(ZernCoef(13:18,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 13-18 ')
	      CALL DACCEPT(ZernCoef(19:24,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 19-24 ')
	      CALL DACCEPT(ZernCoef(25:30,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 25-30 ')
	      CALL DACCEPT(ZernCoef(31:36,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 31-36 ')
	      CALL DACCEPT(ZernCoef(37:42,iElt),DDEF,6,
     &	        'Enter Zernike coefficients 37-42 ')
	      CALL DACCEPT(ZernCoef(43:45,iElt),DDEF,3,
     &	        'Enter Zernike coefficients 43-45 ')

            ELSE
              ZernCoef(:,iElt)=0d0
            END IF
C  End Zernike Parameters

C  Get User-defined Parameters & Data

	    IF (SrfType(iElt).EQ.7) THEN
	      CALL IACCEPT_S(UDSrfType(iElt),1,1,
     &	      	     'Enter user defined surface type ')
	      DO i=1,6
	        UDSrfParam(i,iElt)=0d0
	      END DO
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
	      CALL DACCEPT(UDSrfParam(:,iElt),UDSrfParam(:,iElt),6,
     &	      'Enter user defined surface parameters 1-6')
	      CALL IntToChar(Cnum,iElt,i)
	      CALL CACCEPT(UDSrfFile(iElt),
     &	      filnam(1:ICLEN(filnam))//'.Srf'//Cnum(1:i),
     &       'Enter user defined surface coefficients file name:')
	      CALL UDSinit(EXIST,iElt)
	    END IF
C  End User-defined Parameters & Data

C  Get GridSurf Parameters & Data
	    IF (SrfType(iElt).EQ.9) THEN
	      CALL IACCEPT_S(nGridMat(iElt),mGridMat,1,
     &	        'Enter dimension of gridded data matrix ')
	      CALL DACCEPT_S(GridSrfdx(iElt),1d0,1,
     &	        'Enter spacing of data grids ')
	      CALL IntToChar(Cnum,iElt,i)
	      CALL CACCEPT(GridFile(iElt),
     &	        filnam(1:ICLEN(filnam))//'.Srf'//Cnum(1:i),
     &          'Enter gridded data surface file name:')
	      INQUIRE (FILE=GridFile(iElt),EXIST=EXIST)
	      CALL GridInit(GridFile(iElt),jGridSrf,
     &               nGridMat(iElt),ifGridInit,.FALSE.)
	    END IF
C  End GridSurf Parameters & Data

C  Get Surface Coordinate Parameters
	    ! EltID(*) indicates what type of element (e.g. Reflector,
            ! Refractor etc.) that element is.
	    ! EltID(iElt)=10 says the element is a 'LensArray'
	    ! (see prop_defs.inc). -jzlou
	    !
	    IF ((SrfType(iElt).GE.4).OR.(EltID(iElt).EQ.10)) THEN
	      CALL DEQUATE(pMon(1,iElt),VptElt(1,iElt),3)
	      CALL DACCEPT(pMon(:,iElt),pMon(:,iElt),3,
     &	        'Enter surface reference point location (x,y,z) ')
	      xMon(1,iElt)=psiElt(3,iElt)
	      xMon(2,iElt)=psiElt(1,iElt)
	      xMon(3,iElt)=psiElt(2,iElt)
	      CALL DACCEPT(xMon(:,iElt),xMon(:,iElt),3,
     &	        'Enter surface x-axis unit vector (x,y,z) ')
	      CALL DXPROD(yMon(1,iElt),psiElt(1,iElt),xMon(1,iElt))
	      CALL DACCEPT(yMon(:,iElt),yMon(:,iElt),3,
     &	        'Enter surface y-axis unit vector (x,y,z) ')
	      CALL DEQUATE(zMon(1,iElt),psiElt(1,iElt),3)
	      CALL DACCEPT(zMon(:,iElt),zMon(:,iElt),3,
     &	        'Enter surface z-axis unit vector (x,y,z) ')
	      CALL DORTHOGANALIZE(xMon(1,iElt),yMon(1,iElt),
     &	        ZMon(1,iElt))
	      WRITE(*,*)' Orthoganalized zMon=',(zMon(i,iElt),i=1,3)

	      if (.false.) then
	        print*,'iElt=',iElt
	        print*,'In computing xMon,yMon,zMon ...'
	        stop
	      end if
	    END IF

	    IF ((SrfType(iElt).EQ.4).OR.(SrfType(iElt).EQ.8)) THEN
	      CALL DACCEPT_S(lMon(iElt),Aperture,1,
     &	      'Enter scale length ')
	    ELSE
	      lMon(iElt)=1d0
	    END IF
C  End Surface Coordinate Parameters

C  Get HOE/Grating Parameters
	    IF (EltID(iElt).EQ.4) THEN
	      CALL DACCEPT(h1HOE(:,iElt),DDEF,-3,
     &	      'Enter input reference point (x,y,z) ')
	      CALL DACCEPT(h2HOE(:,iElt),DDEF,-3,
     &	      'Enter output reference point (x,y,z) ')
	      CALL DACCEPT_S(OrderHOE(iElt),1d0,1,
     &	      'Enter diffraction order ')
	      CALL DACCEPT_S(WaveHOE(iElt),1d-7,1,
     &	      'Enter reference wavelength ')
	    END IF
	    IF ((EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13)) THEN
	      CALL DACCEPT_S(OrderHOE(iElt),1d0,1,
     &	      'Enter diffraction order ')
	      CALL DACCEPT_S(RuleWidth(iElt),1d0,1,
     &	      'Enter spacing of grating ruling ')
	      CALL DACCEPT(h1HOE(:,iElt),DDEF,-3,
     &	      'Enter vector perp to grating ruling (x,y,z) ')
	    END IF

C  End HOE/Grating Parameters
C  Get LensArray Parameters

	    IF (EltID(iElt).EQ.10) THEN
	      CALL DACCEPT_S(LensArrayWidth(iElt),1d0,1,
     &	      'Enter lenslet width ')
	      CALL IACCEPT_S(LensArrayType(iElt),1,1,
     &	      'Enter lens array type (1=hex, 2=square) ')
	    END IF

C  End LensArray Parameters

	    CALL ShowData(iElt)

	    CALL CACCEPT(ANS,'YES','Is this correct?')
	    IF ((ANS(1:1).EQ.'n').OR.(ANS(1:1).EQ.'N')) GO TO 11
 15	  CONTINUE

C	Specify output coordinates

	  nOutCord=5

	  invec7(1)=-psiElt(3,nElt)
	  invec7(2)=-psiElt(1,nElt)
	  invec7(3)=-psiElt(2,nElt)
	  CALL DACCEPT(invec7,invec7,3,'Enter xOut (x,y,z) ')

	  CALL DXPROD(invec6,psiElt,invec7)
	  CALL DNEGATE(invec6,invec6,3)
	  CALL DACCEPT(invec6,invec6,3,'Enter yOut (x,y,z) ')

	  CALL DZERO(Tout(1,1),49)
	  DO iCoord=1,3
	    Tout(1,iCoord)=invec7(iCoord)
	    Tout(3,3+iCoord)=invec7(iCoord)
	    Tout(2,iCoord)=invec6(iCoord)
	    Tout(4,3+iCoord)=invec6(iCoord)
	  END DO
	  Tout(5,7)=1

	  CALL ViewData()   ! output to SCreen
	END IF

	! Finalise Configuration
	npts=nGridpts-1

        ! Create Rx
        IF (iopt.EQ.0) THEN
          CALL SetPrtToFileID(RX_FID)  ! re-direct output to File
          CALL ViewData()
          CALL SetPrtToFileID(0)       ! re-direct to screen
        END IF

	iopt=0

        ! Close input file
        CLOSE (UNIT=RX_FID)


	! Setup ApType and ApVec for all segments if an aperture
        ! is defined.
        if (SegApType_FLG) Call SetSegAp

	! Final call
	Call RecordElementID(-1,"",.true.)

	if (.not. opt_npts_flg) then
	  opt_npts=nGridpts/2-1
	  opt_npts_flg=.true.
	end if

        LOAD_SUCCESS = .TRUE.
        RETURN

        !
        ! Rx loading failed
        !
   99   CONTINUE
        LOAD_SUCCESS = .FALSE.

        INQUIRE(UNIT=RX_FID, OPENED=isRxOpen)
        IF (isRxOpen) CLOSE(UNIT=RX_FID)

	RETURN
	END SUBROUTINE MBFile6


C***********************************************************************

	SUBROUTINE SetSegAp
        use param_mod
        use src_mod
        use elt_mod

	Integer :: iElt,iSeg
	Real*8 :: h,a,theta
	Logical :: LCMP
	Real*8, parameter :: pi=3.1415926

	iSeg=0; a=dsqrt(3.d0)/3.d0
        h=width/2

        theta=-pi/2.d0  ! -90 degree

	Do iElt=1,nElt
	  if (EltID(iElt)==11 .and.
     &        (.not.LCMP(EltName(iElt),'CenterSeg',9))) then
	    iSeg=iSeg+1
	    ApType(iElt)=SegApType
	    ApVec(1,iElt)=SegApVec(1) ! size of aperture
	    ApVec(4,iElt)=SegApVec(2) ! orientation

	    ! Compute segment center coords
	    ApVec(3,iElt)=SegCoord(1,iSeg)*h  ! y coord
	    ApVec(2,iElt)=a*(-SegCoord(1,iSeg)+2.*SegCoord(2,iSeg))*h ! x coord

	    !ApVec(2,iElt)=SegCoord(1,iSeg)*h  ! x coord
	    !ApVec(3,iElt)=a*(-SegCoord(1,iSeg)+2.*SegCoord(2,iSeg))*h ! y coord


	    Call SetHexBound(ApVec(1,iElt),iElt)

C	    print*,'SegCoord(1:2,iSeg),ApType,ApVec(1:4)=',
C    &           SegCoord(1:2,iSeg),ApType(iElt),ApVec(1:4,iElt)

	    ! Define xGrid axis for segment, needed in the aperture and
	    ! obscuration routine ChkRayTrans.
	    ! zObs is defined as psiElt, the normal to element surface,
	    ! in the beginning of CPROPGATE of propsub.F,
	    ! and yObs is cross product of xObs and zObs, also computed
	    ! in CPROPGATE. -jzlou
 	    !xObs(1:3,iElt)=0.d0; xObs(1,iElt)=-1
 	    !xObs(1:3,iElt)=0.d0
            !xObs(1,iElt)=-0.7071; xObs(2,iElt)=0.7071

	    xObs(1,iElt)=SegXgrid(1)*cos(theta)+SegXgrid(2)*sin(theta)
	    xObs(2,iElt)=-SegXgrid(1)*sin(theta)+SegXgrid(2)*cos(theta)

 	  else if (EltID(iElt)==11 .and.
     &             LCMP(EltName(iElt),'CenterSeg',9)) then
	    ! Skip center segment
	    iSeg=iSeg+1
	  end if
	End Do
	END SUBROUTINE SetSegAp


C***********************************************************************
C  SUBROUTINE MOD_LOH is called only for 'MODIFY' command -jzlou
C***********************************************************************

	SUBROUTINE MOD_LOH(npts)
	use param_mod
	use src_mod
	use elt_mod
	use cfiles_mod
        use surfsub, only: GridInit

	implicit none

cHistory
c  24Jul97 loh:  Added flag to indicate zernike type

	CHARACTER*4 IND(3)
        CHARACTER(len=MacosCharLen) :: STRING
	!CHARACTER*128 STRING
	CHARACTER(len=MacosVarNamLen) :: VAR_NAM
        CHARACTER(len=MacosValLen) :: VALUE
	!CHARACTER*128 VALUE
	LOGICAL ARR_FLG,LCMP,EXIST,LNOTSET,ifGridInit,ifUDSinit
	INTEGER MB_I(3),MB_J,npts,i,j,ICLEN
	REAL*8 S1,S2

	CHARACTER(len=MacosCharLen) :: TypeString

  10	FORMAT(1P,A,' = ',3D17.9)
  11	FORMAT(1P,A,' = ',D17.9)
  12	FORMAT(A,' = ',i4)
  13 	FORMAT(1P,A,'(',i3,',',i3,') = ',i4)
  17	FORMAT(1P,A,'(',i3,') = ',3D17.9)
  18	FORMAT(1P,A,'(',i3,') = ',D17.9)
  19	FORMAT(A,'(',i3,') = ',i4)
  20  	FORMAT(A,'(',i3,') = ',A8)
  21 	FORMAT(1P,A,'(',i3,') = ',4D17.9)
  22	FORMAT(1P,A,'(1:6,',i3,') = ',6D17.9)
  23	FORMAT(1P,A,'(7:12,',i3,') = ',6D17.9)
  24	FORMAT(1P,A,'(13:18,',i3,') = ',6D17.9)
  25	FORMAT(1P,A,'(19:24,',i3,') = ',6D17.9)
  26	FORMAT(1P,A,'(25:30,',i3,') = ',6D17.9)
  27	FORMAT(1P,A,'(31:36,',i3,') = ',6D17.9)
  271	FORMAT(1P,A,'(37:42,',i3,') = ',6D17.9)
  272	FORMAT(1P,A,'(43:45,',i3,') = ',3D17.9)
  28 	FORMAT(1P,A,'(',i3,',',i3,') = ',3D17.9)
  29 	FORMAT(1P,A,'(',i3,',',i3,') = ',4D17.9)
  30 	FORMAT(1P,A,'(',i3,',',i3,') = ',6D17.9)
  31	FORMAT(1P,A,'(',i3,') = ',7D17.9)
  32	FORMAT(A,'(',i3,') = ',3i4)
  33 	FORMAT(1P,A,'(',i3,') = ',5D17.9)
  34  	FORMAT(A,'(',i3,') = ',A)
  35  	FORMAT(A,'(',i3,',',i3,') = ',A)
  36  	FORMAT(A,' = ',A)
 100    FORMAT(/'COMP MODify Function'/
     &	'MODify Syntax Examples:'/
     &	'  Flux=3'/
     &	'  ChfRayPos(1)=2-MOD 1st element of vector'/
     &	'  ChfRayPos(:)=2,0,0-MOD all elements of vector'/
     &	'  RptElt(3,2)=3-MOD 3rd element of vector for iElt=2'/
     &  '  RptElt(:,2)=1.5,3,0-MOD all elements of '/
     &  '  vector for iElt=2'/
     &	'Quit: quit Modify function.')
 130	FORMAT(1P,'Invalid Element Type, Element=',A)
 132	FORMAT(1P,'Invalid Surface Type, Surface=',A)
 133	FORMAT(1P,'Invalid Zernike Type, Surface=',A)

C       read in command

	STRING(:)=' '
  998   CALL CACCEPT(STRING,'q','MOD:')
	STRING(66:128)=' '

	IF (LCMP(STRING,'QUIT',1)) THEN
	  RETURN
	ELSE IF (LCMP(STRING,'HELP',2)) THEN
	  WRITE(*,100)
	END IF

C       parse the statement


	CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,0)  ! String parsing and not parsing from a Rx file: RX_FID=0

	IF (LCMP(VAR_NAM,'A',1)) THEN
	    IF (LCMP(VAR_NAM,'Aperture',8)) THEN
	      READ(VALUE,*,ERR=999,END=999)Aperture
	    WRITE(*,11)' Aperture',Aperture
	    ELSE IF (LCMP(VAR_NAM,'AnaCoef',7)) THEN
	      IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (AnaCoef(MB_J,MB_I(2)),MB_J=1,4)
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        AnaCoef(MB_I(1),MB_I(2))
	      END IF
	      WRITE(*,21)' AnaCoef',MB_I(2),
     &	      (AnaCoef(MB_J,MB_I(2)),MB_J=1,4)
	    ELSE IF (LCMP(VAR_NAM,'AsphCoef',8)) THEN
	      IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (AsphCoef(MB_J,MB_I(2)),MB_J=1,maxAsphCoef)
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        AsphCoef(MB_I(1),MB_I(2))
	      END IF
	      WRITE(*,21)' AsphCoef',MB_I(2),
     &	      (AsphCoef(MB_J,MB_I(2)),MB_J=1,maxAsphCoef)
	    ELSE IF (LCMP(VAR_NAM,'ApType',6)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
c	      READ(VALUE,*,ERR=999,END=999)ApType(MB_I(1))
	      DO i=1,2
	        IF (LCMP(VALUE,ApTypeName(i),2)) ApType(MB_I(1))=i
	      END DO
 	      WRITE(*,34)' ApType',MB_I(1),ApTypeName(ApType(MB_I(1)))
	    ELSE IF (LCMP(VAR_NAM,'ApVec',5)) THEN
	      IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (ApVec(MB_J,MB_I(2)),MB_J=1,3)
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)ApVec(MB_I(1),MB_I(2))
	      END IF
	      IF (ApType(MB_I(2)).EQ.1) THEN
	        WRITE(*,17)' ApVec',MB_I(2),(ApVec(MB_J,MB_I(2)),MB_J=1,3)
	      ELSE IF (ApType(MB_I(2)).EQ.1) THEN
	        WRITE(*,21)' ApVec',MB_I(2),(ApVec(MB_J,MB_I(2)),MB_J=1,4)
	      ELSE
	        WRITE(*,33)' ApVec',MB_I(2),(ApVec(MB_J,MB_I(2)),MB_J=1,5)
	      END IF
	    END IF
	END IF
	IF (LCMP(VAR_NAM,'C',1)) THEN
	    IF (LCMP(VAR_NAM,'ChfRayDir',9)) THEN
	    IF (ARR_FLG) THEN
	        READ(VALUE,*,ERR=999,END=999)(ChfRayDir(MB_J),MB_J=1,3)
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(VALUE,*,ERR=999,END=999)ChfRayDir(MB_I(1))
	    END IF
	    CALL DUNITIZE(ChfRayDir)
	    WRITE(*,10)' ChfRayDir',(ChfRayDir(MB_J),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'ChfRayPos',9)) THEN
	    IF (ARR_FLG) THEN
	        READ(VALUE,*,ERR=999,END=999)(ChfRayPos(MB_J),MB_J=1,3)
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(VALUE,*,ERR=999,END=999)ChfRayPos(MB_I(1))
	    END IF
	    WRITE(*,10)' ChfRayPos',(ChfRayPos(MB_J),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'CoatIndxElt',8)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (CoatIndxElt(MB_J,MB_I(2)),MB_J=1,3)
  	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        CoatIndxElt(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' CoatIndxElt',MB_I(2),
     &	    (CoatIndxElt(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'CoatThkElt',7)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (CoatThkElt(MB_J,MB_I(2)),MB_J=1,3)
  	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        CoatThkElt(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' CoatThkElt',MB_I(2),
     &	    (CoatThkElt(MB_J,MB_I(2)),MB_J=1,3)
	  END IF

	ELSE IF (LCMP(VAR_NAM,'E',1)) THEN
	    IF (LCMP(VAR_NAM,'EltName',7)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)EltName(MB_I(1))
 	      WRITE(*,20)' EltName',MB_I(1),EltName(MB_I(1))
	    ELSE IF (LCMP(VAR_NAM,'EltType',7)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)EltType(MB_I(1))
  	      WRITE(*,19)' EltType',MB_I(1),EltType(MB_I(1))
	      CALL EltConvert(EltType(MB_I(1)),
     &	      EltID(MB_I(1)),SrfType(MB_I(1)))
	      WRITE(*,*)' Element='//EltTypeName(EltID(MB_I(1)))(1:12)
	      WRITE(*,*)' Surface='//SrfTypeName(SrfType(MB_I(1)))(1:12)
	    ELSE IF (LCMP(VAR_NAM,'Element',7)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)TypeString
	      LNOTSET=.TRUE.
	      DO i=1,15
	        IF (LCMP(TypeString,EltTypeName(i),6)) THEN
	          EltID(MB_I(1))=i
	          LNOTSET=.FALSE.
	        END IF
	      END DO
	      IF (LNOTSET) THEN
		WRITE(*,130) TypeString
	      ELSE
 	        WRITE(*,34)
     & 	        ' Element',MB_I(1),EltTypeName(EltID(MB_I(1)))
	      END IF

	    ELSE IF (LCMP(VAR_NAM,'eElt',4)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)eElt(MB_I(1))
 	    WRITE(*,18)' eElt',MB_I(1),eElt(MB_I(1))
	    KcElt(MB_I(1))=-eElt(MB_I(1))*eElt(MB_I(1))
	    KrElt(MB_I(1))=-(1.0+eElt(MB_I(1)))*fElt(MB_I(1))
 	    WRITE(*,18)' KcElt',MB_I(1),KcElt(MB_I(1))
 	    WRITE(*,18)' KrElt',MB_I(1),KrElt(MB_I(1))
	  ELSE IF (LCMP(VAR_NAM,'Extinc',6)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)Extinc(MB_I(1))
 	    WRITE(*,18)' Extinc',MB_I(1),Extinc(MB_I(1))
	  END IF

	ELSE IF (LCMP(VAR_NAM,'F',1)) THEN
	    print*,'Debug: VAR_NAM = ', VAR_NAM
	    IF (LCMP(VAR_NAM,'FLUX',4)) THEN
	      print*,'Debug: VALUE = ', VALUE
	      READ(VALUE,*)Flux
 	      WRITE(*,11)' Flux',Flux
	    ELSE IF (LCMP(VAR_NAM,'fElt',4)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)fElt(MB_I(1))
 	      WRITE(*,18)' fElt',MB_I(1),fElt(MB_I(1))
	      KrElt(MB_I(1))=-(1.0+eElt(MB_I(1)))*fElt(MB_I(1))
 	      WRITE(*,18)' KrElt',MB_I(1),KrElt(MB_I(1))
	    END IF

	ELSE IF (LCMP(VAR_NAM,'G',1)) THEN
	  IF (LCMP(VAR_NAM,'GridType',8)) THEN
c	    READ(VALUE,*)GridType
	    DO i=1,5
	      IF (LCMP(VALUE,GridTypeName(i),2)) GridType=i
	    END DO
 	    WRITE(*,36)' GridType',GridTypeName(GridType)
	  ELSE IF (LCMP(VAR_NAM,'GlassElt',6)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    LNOTSET=.TRUE.
	    i=1
	    DO WHILE (LNOTSET.AND.(i.LE.mGlass))
	      IF (LCMP(VALUE,GlassName(i),
     &	      ICLEN(GlassName(i)))) THEN
	        LNOTSET=.FALSE.
	        CALL DEQUATE(GlassCoef(1,MB_I(1)),GlassTable(1,i),6)
	      END IF
     	      i=i+1
	    END DO
	    IF (LNOTSET) THEN
	      WRITE(*,*)
     &	        ' Unknown glass, MOD GlassCoef to enter properties'
	    ELSE
	      GlassElt(MB_I(1))=VALUE
 	      WRITE(*,20)' GlassElt',MB_I(1),GlassElt(MB_I(1))
	    END IF
	  ELSE IF (LCMP(VAR_NAM,'GlassCoef',6)) THEN
	    IF (ARR_FLG) THEN
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	        (GlassCoef(MB_J,MB_I(2)),MB_J=1,6)
	    ELSE
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	        GlassCoef(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,31)' GlassCoef',MB_I(2),
     &	    (GlassCoef(MB_J,MB_I(2)),MB_J=1,6)
	  ELSE IF (LCMP(VAR_NAM,'gap',3)) THEN
	    READ(VALUE,*)gap
 	    WRITE(*,11)' gap',gap
	  ELSE IF (LCMP(VAR_NAM,'GridFile',5)) THEN
	    READ(IND(1),*,ERR=999,END=999) MB_I(1)
	    READ(VALUE,*,ERR=999,END=999) GridFile(MB_I(1))
	    WRITE(*,*)' GridFile(',MB_I(1),')=',GridFile(MB_I(1))
	    CALL GridInit(GridFile(MB_I(1)),iEltToGridSrf(MB_I(1)),
     &	    	   nGridMat(MB_I(1)),ifGridInit,.FALSE.)
	  ELSE IF (LCMP(VAR_NAM,'GridSrfdx',5)) THEN
	    READ(VALUE,*) GridSrfdx
 	    WRITE(*,11)' GridSrfdx',GridSrfdx
	  END IF

	ELSE IF (LCMP(VAR_NAM,'H',1)) THEN
 	  IF (LCMP(VAR_NAM,'h1HOE',5)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (h1HOE(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        h1HOE(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' h1HOE',MB_I(2),(h1HOE(MB_J,MB_I(2)),MB_J=1,3)
 	  ELSE IF (LCMP(VAR_NAM,'h2HOE',5)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (h2HOE(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        h2HOE(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' h2HOE',MB_I(2),(h2HOE(MB_J,MB_I(2)),MB_J=1,3)
 	  END IF

	ELSE IF (LCMP(VAR_NAM,'IndRef',6)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)IndRef(MB_I(1))
 	  WRITE(*,18)' IndRef',MB_I(1),IndRef(MB_I(1))

	ELSE IF (LCMP(VAR_NAM,'K',1)) THEN
	  IF (LCMP(VAR_NAM,'KrElt',5)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)KrElt(MB_I(1))
	    IF (KrElt(MB_I(1)).GT.0d0) KrElt(MB_I(1))=-KrElt(MB_I(1))
 	    WRITE(*,18)' KrElt',MB_I(1),KrElt(MB_I(1))
	    IF (KcElt(MB_I(1)).LT.0) THEN
 	      fElt(MB_I(1)) = ABS(KrElt(MB_I(1)))/(1d0+eElt(MB_I(1)))
 	      WRITE(*,18)' fElt',MB_I(1),fElt(MB_I(1))
	    END IF
	  ELSE IF (LCMP(VAR_NAM,'KcElt',5)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)KcElt(MB_I(1))
 	    WRITE(*,18)' KcElt',MB_I(1),KcElt(MB_I(1))
	    IF (KcElt(MB_I(1)) .GT. 0.0) THEN
	      IF ((SrfType(MB_I(1)).EQ.1).OR.
     &        (SrfType(MB_I(1)).EQ.6)) THEN
	        WRITE(*,*) ' Oblate spheroid not valid for this ',
     &	        'surface type.'
	      ELSE
	        fElt(MB_I(1)) = 0.0
	      END IF
	    ELSE
	      eElt(MB_I(1)) = DSQRT(-KcElt(MB_I(1)))
	      fElt(MB_I(1)) = -KrElt(MB_I(1))/(1.0+eElt(MB_I(1)))
 	      WRITE(*,18)' fElt',MB_I(1),fElt(MB_I(1))
 	      WRITE(*,18)' eElt',MB_I(1),eElt(MB_I(1))
	    END IF
	  END IF

	ELSE IF (LCMP(VAR_NAM,'L',1)) THEN
	  IF (LCMP(VAR_NAM,'lMon',4)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)lMon(MB_I(1))
 	    WRITE(*,18)' lMon',MB_I(1),lMon(MB_I(1))
 	  ELSE IF (LCMP(VAR_NAM,'LensArrayType',12)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)LensArrayType(MB_I(1))
 	    WRITE(*,12)' LensArrayType',MB_I(1),LensArrayType(MB_I(1))
 	  ELSE IF (LCMP(VAR_NAM,'LensArrayWidth',12)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)LensArrayWidth(MB_I(1))
 	    WRITE(*,18)' LensArrayWidth',MB_I(1),LensArrayWidth(MB_I(1))
 	  END IF

	ELSE IF (LCMP(VAR_NAM,'M',1)) THEN
	    IF (LCMP(VAR_NAM,'MonCoef',8)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)MonCoef(MB_I(1),MB_I(2))
	    IF (MB_I(1) .LE. 6) THEN
	      WRITE(*,22)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=1,6)
	    ELSE IF ((MB_I(1) .GT. 6) .AND. (MB_I(1) .LE. 12)) THEN
	      WRITE(*,23)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=7,12)
	    ELSE IF ((MB_I(1) .GT. 12) .AND. (MB_I(1) .LE. 18)) THEN
	      WRITE(*,24)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=13,18)
	    ELSE IF ((MB_I(1) .GT. 18) .AND. (MB_I(1) .LE. 24)) THEN
	      WRITE(*,25)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=19,24)
	    ELSE IF ((MB_I(1) .GT. 24) .AND. (MB_I(1) .LE. 30)) THEN
	      WRITE(*,26)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=25,30)
	    ELSE IF ((MB_I(1) .GT. 31) .AND. (MB_I(1) .LE. 36)) THEN
	      WRITE(*,27)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=31,36)
	    ELSE IF ((MB_I(1) .GT. 36) .AND. (MB_I(1) .LE. 42)) THEN
	      WRITE(*,23)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=37,42)
	    ELSE IF ((MB_I(1) .GT. 42) .AND. (MB_I(1) .LE. 48)) THEN
	      WRITE(*,24)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=43,48)
	    ELSE IF ((MB_I(1) .GT. 48) .AND. (MB_I(1) .LE. 54)) THEN
	      WRITE(*,25)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=49,54)
	    ELSE IF ((MB_I(1) .GT. 54) .AND. (MB_I(1) .LE. 60)) THEN
	      WRITE(*,26)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=55,60)
	    ELSE IF ((MB_I(1) .GT. 61) .AND. (MB_I(1) .LE. 66)) THEN
	      WRITE(*,27)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=61,66)
	    ELSE IF ((MB_I(1) .GT. 66) .AND. (MB_I(1) .LE. 72)) THEN
	      WRITE(*,23)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=67,72)
	    ELSE IF ((MB_I(1) .GT. 72) .AND. (MB_I(1) .LE. 78)) THEN
	      WRITE(*,24)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=73,78)
	    ELSE IF ((MB_I(1) .GT. 78) .AND. (MB_I(1) .LE. 84)) THEN
	      WRITE(*,25)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=79,84)
	    ELSE IF ((MB_I(1) .GT. 84) .AND. (MB_I(1) .LE. 90)) THEN
	      WRITE(*,26)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=85,90)
	    ELSE IF ((MB_I(1) .GT. 91) .AND. (MB_I(1) .LE. 96)) THEN
	      WRITE(*,27)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=91,96)
	    ELSE IF ((MB_I(1) .GT. 96) .AND. (MB_I(1) .LE. 102)) THEN
	      WRITE(*,23)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=97,102)
	    ELSE IF ((MB_I(1) .GT. 102) .AND. (MB_I(1) .LE. 108)) THEN
	      WRITE(*,24)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=103,108)
	    ELSE IF ((MB_I(1) .GT. 108) .AND. (MB_I(1) .LE. 114)) THEN
	      WRITE(*,25)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=109,114)
	    ELSE IF ((MB_I(1) .GT. 115) .AND. (MB_I(1) .LE. 120)) THEN
	      WRITE(*,26)' MonCoef',MB_I(2),(MonCoef(MB_J,MB_I(2)),
     &	       MB_J=115,120)
	    END IF
	  END IF

	ELSE IF (LCMP(VAR_NAM,'N',1)) THEN
	  IF (LCMP(VAR_NAM,'npts',4)) THEN
	    READ(VALUE,*)npts
	    nGridpts=npts+1
  	    WRITE(*,12)' npts',npts
	  ELSE IF (LCMP(VAR_NAM,'nGridpts',8)) THEN
	    READ(VALUE,*)nGridpts
	    IF (nGridpts.GT.mpts) nGridpts=mpts
	    npts=nGridpts-1
  	    WRITE(*,12)' nGridpts',nGridpts
	  ELSE IF (LCMP(VAR_NAM,'nSeg',4)) THEN
	    READ(VALUE,*)nSeg
  	    WRITE(*,12)' nSeg',nSeg
	  ELSE IF (LCMP(VAR_NAM,'nPetals',5)) THEN
	    READ(VALUE,*)nPetals
  	    WRITE(*,12)' nPetals',nPetals
	  ELSE IF (LCMP(VAR_NAM,'nECoord',7)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)nECoord(MB_I(1))
  	    WRITE(*,19)' nECoord',MB_I(1),nECoord(MB_I(1))
	  ELSE IF (LCMP(VAR_NAM,'nGridMat',6)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)nGridMat(MB_I(1))
  	    WRITE(*,19)' nGridMat',MB_I(1),nGridMat(MB_I(1))
	  ELSE IF (LCMP(VAR_NAM,'nElt',4)) THEN
	    READ(VALUE,*)nElt
  	    WRITE(*,12)' nElt',nElt
	  ELSE IF (LCMP(VAR_NAM,'nCoatElt',5)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)nCoatElt(MB_I(1))
 	    WRITE(*,19)' nCoatElt',MB_I(1),nCoatElt(MB_I(1))
	  ELSE IF (LCMP(VAR_NAM,'nObs',4)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)nObs(MB_I(1))
  	    WRITE(*,19)' nObs',MB_I(1),nObs(MB_I(1))
	  ELSE IF (LCMP(VAR_NAM,'nOutCord',8)) THEN
	    READ(VALUE,*)nOutCord
  	    WRITE(*,12)' nOutCord',nOutCord
	  END IF

	ELSE IF (LCMP(VAR_NAM,'O',1)) THEN
	    IF (LCMP(VAR_NAM,'Obscratn',8)) THEN
	      READ(VALUE,*)Obscratn
  	    WRITE(*,11)' Obscratn',Obscratn
	    ELSE IF (LCMP(VAR_NAM,'ObsType',7)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
c	      READ(VALUE,*,ERR=999,END=999)ObsType(MB_I(1),MB_I(2))
	      DO i=-3,3
	        IF (LCMP(VALUE,ObsTypeName(i),4))
     &	          ObsType(MB_I(1),MB_I(2))=i
	      END DO
 	      WRITE(*,35)' ObsType',
     & 	      MB_I(1),MB_I(2),ObsTypeName(ObsType(MB_I(1),MB_I(2)))
	    ELSE IF (LCMP(VAR_NAM,'ObsVec',6)) THEN
	      IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(IND(3),*,ERR=999,END=999)MB_I(3)
	        IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 1) THEN
	          READ(VALUE,*,ERR=999,END=999)
     &	          (ObsVec(MB_J,MB_I(2),MB_I(3)),MB_J=1,3)
	        ELSE IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 4) THEN
	          READ(VALUE,*,ERR=999,END=999)
     &	          (ObsVec(MB_J,MB_I(2),MB_I(3)),MB_J=1,6)
	        ELSE
	          READ(VALUE,*,ERR=999,END=999)
     &	          (ObsVec(MB_J,MB_I(2),MB_I(3)),MB_J=1,4)
	        END IF
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(IND(3),*,ERR=999,END=999)MB_I(3)
	        READ(VALUE,*,ERR=999,END=999)
     &	        ObsVec(MB_I(1),MB_I(2),MB_I(3))
	      END IF
	      IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 1) THEN
	        WRITE(*,28)' ObsVec',MB_I(2),MB_I(3),
     &	        (ObsVec(MB_J,MB_I(2),MB_I(3)),
     &	        MB_J=1,3)
	      ELSE IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 4) THEN
	        WRITE(*,30)' ObsVec',MB_I(2),MB_I(3),
     &	        (ObsVec(MB_J,MB_I(2),MB_I(3)),
     &	        MB_J=1,6)
	      ELSE IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 2) THEN
	        WRITE(*,29)' ObsVec',MB_I(2),MB_I(3),
     &	        (ObsVec(MB_J,MB_I(2),MB_I(3)),
     &		MB_J=1,4)
	      ELSE IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 3) THEN
	        WRITE(*,29)' ObsVec',MB_I(2),MB_I(3),
     &	        (ObsVec(MB_J,MB_I(2),MB_I(3)),
     &		MB_J=1,4)
	      ELSE IF (ABS(ObsType(MB_I(2),MB_I(3))) .EQ. 5) THEN
	        WRITE(*,29)' ObsVec',MB_I(2),MB_I(3),
     &          (ObsVec(MB_J,MB_I(2),MB_I(3)),
     &          MB_J=1,4)
	      END IF
 	  ELSE IF (LCMP(VAR_NAM,'OrderHOE',8)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)OrderHOE(MB_I(1))
 	    WRITE(*,18)' OrderHOE',MB_I(1),OrderHOE(MB_I(1))
	  END IF

	ELSE IF (LCMP(VAR_NAM,'P',1)) THEN
	    IF (LCMP(VAR_NAM,'psiElt',6)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (psiElt(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        psiElt(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' psiElt',MB_I(2),(psiElt(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'pMon',4)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (pMon(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        pMon(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' pMon',MB_I(2),(pMon(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'PropType',8)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
c	      READ(VALUE,*,ERR=999,END=999)PropType(MB_I(1))
	      DO i=1,15
	        IF (LCMP(VALUE,PropTypeName(i),5)) PropType(MB_I(1))=i
	      END DO
 	      WRITE(*,34)' PropType',
     & 	      MB_I(1),PropTypeName(PropType(MB_I(1)))
	    ELSE IF (LCMP(VAR_NAM,'PinHole',7)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)PinHole(MB_I(1))
	    WRITE(*,18)' PinHole',MB_I(1),PinHole(MB_I(1))
	  END IF

	ELSE IF (LCMP(VAR_NAM,'R',1)) THEN
	  IF (LCMP(VAR_NAM,'RptElt',6)) THEN
	    IF (ARR_FLG) THEN
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      (RptElt(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      RptElt(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' RptElt',MB_I(2),(RptElt(MB_J,MB_I(2)),MB_J=1,3)
 	  ELSE IF (LCMP(VAR_NAM,'RuleWidth',9)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)RuleWidth(MB_I(1))
 	    WRITE(*,18)' RuleWidth',MB_I(1),RuleWidth(MB_I(1))
 	  ELSE IF (LCMP(VAR_NAM,'RadCtr',5)) THEN
	    READ(VALUE,*)RadCtr
 	    WRITE(*,11)' RadCtr',RadCtr
	  END IF

	ELSE IF (LCMP(VAR_NAM,'S',1)) THEN
	  IF (LCMP(VAR_NAM,'SegCoord',8)) THEN
	    IF (ARR_FLG) THEN
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      (SegCoord(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      SegCoord(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,32)' SegCoord=',MB_I(2),
     &	      (SegCoord(MB_J,MB_I(2)),MB_J=1,3)
	  ELSE IF (LCMP(VAR_NAM,'SegXgrid',8)) THEN
	    IF (ARR_FLG) THEN
	        READ(VALUE,*,ERR=999,END=999)(SegXgrid(MB_J),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(VALUE,*,ERR=999,END=999)SegXgrid(MB_I(1))
	    END IF
	    WRITE(*,10)' SegXgrid',(SegXgrid(MB_J),MB_J=1,3)

	    ELSE IF (LCMP(VAR_NAM,'SrfType',6)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)SrfType(MB_I(1))
  	      WRITE(*,19)' SrfType',MB_I(1),SrfType(MB_I(1))

	    ELSE IF (LCMP(VAR_NAM,'Surface',7)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)TypeString
	      DO i=1,8
	        IF (LCMP(TypeString,SrfTypeName(i),3)) THEN
		  SrfType(MB_I(1))=i
	          LNOTSET=.FALSE.
	        END IF
	      END DO
	      IF (LNOTSET) THEN
		WRITE(*,132) TypeString
	      ELSE
 	        WRITE(*,34)
     & 	        ' SrfType',MB_I(1),SrfTypeName(SrfType(MB_I(1)))
	      END IF
	  END IF

	ELSE IF (LCMP(VAR_NAM,'T',1)) THEN
	  IF (LCMP(VAR_NAM,'TElt',4)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(IND(3),*,ERR=999,END=999)MB_I(3)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (TElt(MB_J,MB_I(2),MB_I(3)),MB_J=1,6)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(IND(3),*,ERR=999,END=999)MB_I(3)
	        READ(VALUE,*,ERR=999,END=999)
     &	        TElt(MB_I(1),MB_I(2),MB_I(3))
	    END IF
	      WRITE(*,30)' TElt',MB_I(2),MB_I(3),
     &	      (TElt(MB_J,MB_I(2),MB_I(3)),MB_J=1,6)
	  ELSE IF (LCMP(VAR_NAM,'Tout',4)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (Tout(MB_I(1),MB_J),MB_J=1,7)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        Tout(MB_I(1),MB_I(2))
	    END IF
	      WRITE(*,31)' Tout ',MB_I(1),(Tout(MB_I(1),MB_J),MB_J=1,7)
	  END IF

	ELSE IF (LCMP(VAR_NAM,'U',1)) THEN
	  IF (LCMP(VAR_NAM,'UDSrfFile',7)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)UDSrfFile(MB_I(1))
	    WRITE(*,*)' UDSrfFile(',MB_I(1),')=',UDSrfFile(MB_I(1))
	    CALL UDSinit(ifUDSinit,MB_I(1))

 	  ELSE IF (LCMP(VAR_NAM,'UDSrfType',7)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*,ERR=999,END=999)UDSrfType(MB_I(1))
 	    WRITE(*,19)' UDSrfType',MB_I(1),UDSrfType(MB_I(1))

 	  ELSE IF (LCMP(VAR_NAM,'UDSrfParam',7)) THEN
	    IF (ARR_FLG) THEN
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      (UDSrfParam(MB_J,MB_I(2)),MB_J=1,6)
	    ELSE
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      UDSrfParam(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,10)' UDSrfParam',(xGrid(MB_J),MB_J=1,6)
	  END IF

	ELSE IF (LCMP(VAR_NAM,'VptElt',6)) THEN
	  IF (ARR_FLG) THEN
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      (VptElt(MB_J,MB_I(2)),MB_J=1,3)
	  ELSE
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)
	      READ(VALUE,*,ERR=999,END=999)
     &	      VptElt(MB_I(1),MB_I(2))
	  END IF
	  WRITE(*,17)' VptElt',MB_I(2),(VptElt(MB_J,MB_I(2)),MB_J=1,3)

	ELSE IF (LCMP(VAR_NAM,'W',1)) THEN
	  IF (LCMP(VAR_NAM,'Wavelen',7)) THEN
	    READ(VALUE,*)Wavelen
	    WRITE(*,11)' Wavelen',Wavelen
	  ELSE IF (LCMP(VAR_NAM,'width',5)) THEN
	    READ(VALUE,*)width
	    WRITE(*,11)' width',width
 	  ELSE IF (LCMP(VAR_NAM,'WaveHOE',7)) THEN
	    READ(IND(1),*,ERR=999,END=999)MB_I(1)
	    READ(VALUE,*)WaveHOE(MB_I(1))
 	    WRITE(*,18)' WaveHOE',MB_I(1),WaveHOE(MB_I(1))
	  END IF

	ELSE IF (LCMP(VAR_NAM,'X',1)) THEN
	    IF (LCMP(VAR_NAM,'xGrid',5)) THEN
	    IF (ARR_FLG) THEN
	        READ(VALUE,*,ERR=999,END=999)(xGrid(MB_J),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(VALUE,*,ERR=999,END=999)xGrid(MB_I(1))
	    END IF
	    WRITE(*,10)' xGrid',(xGrid(MB_J),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'xMon',4)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (xMon(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        xMon(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' xMon',MB_I(2),(xMon(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'xObs',4)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (xObs(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        xObs(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' xObs',MB_I(2),(xObs(MB_J,MB_I(2)),MB_J=1,3)
	  END IF

	ELSE IF (LCMP(VAR_NAM,'Y',1)) THEN
	    IF (LCMP(VAR_NAM,'yGrid',5)) THEN
	    IF (ARR_FLG) THEN
	        READ(VALUE,*,ERR=999,END=999)(yGrid(MB_J),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(VALUE,*,ERR=999,END=999)yGrid(MB_I(1))
	    END IF
	    WRITE(*,10)' yGrid',(yGrid(MB_J),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'yMon',4)) THEN
	    IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (yMon(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        yMon(MB_I(1),MB_I(2))
	    END IF
	    WRITE(*,17)' yMon',MB_I(2),(yMon(MB_J,MB_I(2)),MB_J=1,3)
	  END IF

	ELSE IF (LCMP(VAR_NAM,'Z',1)) THEN
	    IF (LCMP(VAR_NAM,'ZernCoef',8)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(IND(2),*,ERR=999,END=999)MB_I(2)

	      if (MB_I(1)>mZernModes) then
	        ! ZernCoef index out of range
                ! ToDo --> FIX
	      end if

	      ! MB_I(1)=ZernCoef index, MB_I(2)=iElt index. -jzlou
	      READ(VALUE,*,ERR=999,END=999)ZernCoef(MB_I(1),MB_I(2))

	      IF (MB_I(1) .LE. 6) THEN
	       WRITE(*,22)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=1,6)
	      ELSE IF ((MB_I(1) .GT. 6) .AND. (MB_I(1) .LE. 12)) THEN
	       WRITE(*,23)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=7,12)
	      ELSE IF ((MB_I(1) .GT. 12) .AND. (MB_I(1) .LE. 18)) THEN
	       WRITE(*,24)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=13,18)
	      ELSE IF ((MB_I(1) .GT. 18) .AND. (MB_I(1) .LE. 24)) THEN
	       WRITE(*,25)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=19,24)
	      ELSE IF ((MB_I(1) .GT. 24) .AND. (MB_I(1) .LE. 30)) THEN
	       WRITE(*,26)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=25,30)
	      ELSE IF ((MB_I(1) .GT. 31) .AND. (MB_I(1) .LE. 36)) THEN
	       WRITE(*,27)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=31,36)
	      ELSE IF ((MB_I(1) .GT. 37) .AND. (MB_I(1) .LE. 42)) THEN
	       WRITE(*,271)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=37,42)
	      ELSE IF ((MB_I(1) .GT. 43) .AND. (MB_I(1) .LE. 45)) THEN
	       WRITE(*,272)' ZernCoef',MB_I(2),(ZernCoef(MB_J,MB_I(2)),
     &	         MB_J=43,45)
	      END IF
	    ELSE IF (LCMP(VAR_NAM,'zMon',4)) THEN
	      IF (ARR_FLG) THEN
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)
     &	        (zMon(MB_J,MB_I(2)),MB_J=1,3)
	      ELSE
	        READ(IND(1),*,ERR=999,END=999)MB_I(1)
	        READ(IND(2),*,ERR=999,END=999)MB_I(2)
	        READ(VALUE,*,ERR=999,END=999)zMon(MB_I(1),MB_I(2))
	      END IF
	      WRITE(*,17)' zMon',MB_I(2),(zMon(MB_J,MB_I(2)),MB_J=1,3)
	    ELSE IF (LCMP(VAR_NAM,'zElt',4)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)zElt(MB_I(1))
	      WRITE(*,17)' zElt',MB_I(1),zElt(MB_I(1))
	    ELSE IF (LCMP(VAR_NAM,'zSource',7)) THEN
	      READ(VALUE,*)zSource
	      WRITE(*,11)' zSource',zSource
	    ELSE IF (LCMP(VAR_NAM,'ZernType',8)) THEN
	      READ(IND(1),*,ERR=999,END=999)MB_I(1)
	      READ(VALUE,*,ERR=999,END=999)TypeString
	      if (LCMP(TypeString,'Norm',4)) then
	        ! normalized Zernike coeff
	        DO i=4,6
                  IF (LCMP(TypeString,ZernTypeNameL(i),7)) THEN
                    ZernTypeL(MB_I(1))=i
                    LNOTSET=.FALSE.
                  END IF
                END DO
	      else
	        DO i=1,3
	          IF (LCMP(TypeString,ZernTypeNameL(i),3)) THEN
		    ZernTypeL(MB_I(1))=i
	            LNOTSET=.FALSE.
	          END IF
	        END DO
	      end if
	      IF (LNOTSET) THEN
		WRITE(*,133) TypeString
	      ELSE
 	        WRITE(*,34)
     & 	        ' ZernType',MB_I(1),ZernTypeNameL(ZernTypeL(MB_I(1)))
	      END IF
	    END IF
	END IF

	GO TO 998

 999	CONTINUE
	  WRITE(*,*)' Subscript error. Please reenter.'
	  GO TO 998
	END SUBROUTINE MOD_LOH

C***********************************************************************

	SUBROUTINE READ_LOH(MODE)

C       THIS ROUTINE READS THE COMMAND STRING, PARSES IT, AND
C       PASSES THE APPROPRIATE COMMAND BACK TO THE ROUTINE
C       THAT CALLED IT.  A MAXIMUM OF 10 ENTRIES MAY BE PUT ON
C       ONE LINE.  ENTRIES ARE DELIMITED BY A SEMICOLON.  COMMENTS
C       ARE DELIMITED BY A PERCENT SIGN.  A COMMAND LINE MAY BE
C       72 CHARACTERS LONG.

c sws
	use lohpars_mod
c sws
	INTEGER MODE,L
	INTEGER ICLEN

	LOGICAL LCMP

 5	FORMAT(' ',A,$)
 6	FORMAT(' ',A)

C       Are there inputs that have been read, but not used?

	IF (read_len(pstack) .EQ. 0) THEN
	   CALL PARSE_LOH(MODE)
	END IF

c sws
c   Only do this if we got something
	IF (read_len(pstack) .NE. 0) THEN
	   CVAR = var_store(read_cur(pstack),pstack)

C          Write to screen
	   IF ( (read_cur(pstack) .GT. 1) .OR. (bat_task .EQ. 2) )THEN
	      L=ICLEN(CVAR)
	      IF (L.LE.0) L=1
	      WRITE(*,6) CVAR(1:L)
	   END IF

C  Write out journal data to file, if in journal mode
c sws
	   IF (jou_task .EQ. 2) THEN
	      IF (.not. LCMP(CVAR,'jou',3)) THEN
		 WRITE (27,6) CVAR(1:ICLEN(CVAR))
	      END IF
	   END IF
	ELSE
c got nothing, return a blank
	   CVAR = ' '

	END IF

C       Is this the last command saved?

	IF (read_cur(pstack) .EQ. read_len(pstack)) THEN
	   read_len(pstack)=0
	   read_cur(pstack)=0
	ELSE
	   read_cur(pstack)=read_cur(pstack)+1
	END IF

	!print*,'pstack =',pstack
	!print*,'read_len(pstack) =',read_len(pstack)

	RETURN
	END SUBROUTINE READ_LOH

C***********************************************************************

#if 1

	SUBROUTINE PARSE_LOH(Mode)

C       THIS ROUTINE DOES THE ACTUAL PARSING.  A MAXIMUM OF 10
C       ENTRIES MAY BE PUT ON ONE LINE.  ENTRIES ARE DELIMITED
C       BY A SEMICOLON or a blank.  COMMENTS ARE DELIMITED BY A PERCENT
C       SIGN. A COMMAND LINE MAY BE len_max CHARACTERS LONG.
C
        use lohpars_mod

        CHARACTER*130 VAR

        INTEGER Mode,i,j,k,l
        INTEGER len_max
        LOGICAL last_flag,cma
	Integer, parameter :: maxBlanks=10


C    if last flag = true, then the semi colon is the last character
C    on the line, (besides spaces)

c       maximum length of a line
        data len_max /130/

    6   FORMAT(A)

C sws
C  added for jump from popping the stack
 10     continue

c Initialize command stack
        var_store(1:10,pstack)=' '

C       Read new statement to be parsed
   98   CONTINUE

        IF (bat_task == 1) THEN
	   ! command-line mode, -jzlou
	   VAR=' '
#ifdef READLINE_LIBRARY
	   If (Mode /= 0) Then
	     Call MHIST(' ',VAR)
	   Else
	     Call MHIST('MACOS',VAR)
	   End If
#else
           READ(*,6,END=200, ERR=300) VAR
#endif
        ELSE
	  READ(ustack,6,END=100, ERR=300) VAR
#ifdef READLINE_LIBRARY
          Call MADDH(VAR)
#endif
        END IF

	! At this point, 'VAR' holds the entire command line text,
        ! such as "old j18sd_s"  -jzlou
	cmdline=' '; cmdline=VAR

C       i is counter for location in array VAR, j is counter for
C       number of tokens in line, k is counter for var_store,
C       last_flag is set to true at the end of a command (multiple
C       commands can be on same line) -jzlou

        i=1; j=1; k=1
        last_flag=.TRUE.
	cma=.false.

C  Check for % sign, semicolon, blank space...

    1   CONTINUE
        IF (i .LE. len_max) THEN
          IF (VAR(i:i) .EQ. '%') THEN
C           Comment follows, no more input data
            GO TO 2
          ELSE IF (VAR(i:i) .EQ. ';') THEN
C  End of command - multiple commands can be on same line
            j=j+1; k=1; i=i+1
            last_flag=.TRUE.
            GO TO 1
          ELSE IF (VAR(i:i) .EQ. ' ') THEN
            IF (i.GT.1 .and. .not.cma) THEN
              l=i-1
              IF ((var(l:l).NE.' ').AND.(var(l:l).NE.';')) THEN
C  End of command
	        ! start a new line in var_store
                j=j+1; k=1
                i=i+1
                last_flag=.TRUE.
	        GO TO 1
	      ELSE
                i=i+1
	      END IF
              GO TO 1
            END IF
C  Go to next character
            i=i+1
            GO TO 1
          ELSE
	    if (VAR(i:i)==',') then
              cma=.true.
	    else
              cma=.false.
	    end if
            var_store(j,pstack)(k:k)=VAR(i:i)
            k=k+1; i=i+1
            last_flag=.FALSE.
            GO TO 1
          END IF
        END IF
    2   CONTINUE

        IF (last_flag) j=j-1

        !print*,'i,j, pstack=', i,j, pstack
        !print*,'var_store(j,pstack)= ',var_store(j,pstack)

        read_len(pstack)=j
        read_cur(pstack)=1

	! At this point, var_store(1:j,pstack)(:) stores 'j'
        ! tokens read from an input line. -jzlou
#if 0
	print*,'read_len(pstack) =',read_len(pstack)
	print*,'var_store(1,pstack)=',var_store(1,pstack)
#endif
        RETURN

c sws
c current journal file empty
  100   CALL MB_LOH(2)

        IF (pstack .GT. 2) THEN
c still have more .jou files in file stack
c pop the file stack
           pstack = pstack-1
           ustack = ustack-1
c now check current read_len
c IF empty, try to fill it
           IF ( read_len(pstack) .EQ. 0 ) GO TO 10
c ELSE CONTINUE
        ELSE IF (pstack .EQ. 2) THEN
c no more .jou files in the file stack
c pop the file stack to initial values
           pstack = pstack-1
           ustack = ustack-1
c go back to console
           bat_task=1
	   var_store(1:10,pstack)=' ' ! added by jzlou to clean up cmd arg stack after
                                        ! running journal file, 01/19/07
        ELSE
c shouldn't get here
           WRITE(*,*) '[parse_loh]',
     &        ' error, pstack = ',pstack, ' bat_task = ',bat_task
c go back to console
           bat_task=1
        END IF

c STACK empty
 200    CONTINUE

        GO TO 98
c       RETURN

c sws    error RETURN in reading input
 300    CONTINUE
        WRITE(*,*) '[parse_loh]',
     &     ' file error during READ, bat_task = ',bat_task
        RETURN
        END SUBROUTINE PARSE_LOH

#endif

C *********************************************************************

	SUBROUTINE PARSE_LOH_SAVE()

C       THIS ROUTINE DOES THE ACTUAL PARSING.  A MAXIMUM OF 10
C       ENTRIES MAY BE PUT ON ONE LINE.  ENTRIES ARE DELIMITED
C       BY A SEMICOLON or a blank.  COMMENTS ARE DELIMITED BY A PERCENT
C       SIGN. A COMMAND LINE MAY BE len_max CHARACTERS LONG.
C
	use lohpars_mod

	CHARACTER*130 VAR

	INTEGER i,j,k,l
	INTEGER len_max
	LOGICAL last_flag


C    if last flag = true, then the semi colon is the last character
C    on the line, (besides spaces)

c       maximum length of a line
	data len_max /130/

    6   FORMAT(A)

C sws
C  added for jump from popping the stack
 10     continue

c Initialize command stack
	var_store(1:10,pstack)=' '

C       Read new statement to be parsed
   98   CONTINUE

	IF (bat_task .EQ. 1) THEN
c sws
	   READ(*,6,END=200, ERR=300) VAR
	ELSE
	  READ(ustack,6,END=100, ERR=300) VAR
	END IF
c sws

C       i is index into array VAR, j is counter for number of
C       commands in line, k is counter for var_store

	i=1; j=1; k=1
	last_flag=.TRUE.

C  Check for % sign, semicolon, blank space...
    1   CONTINUE
	IF (i .LE. len_max) THEN
	  IF (VAR(i:i) .EQ. '%') THEN
C           Comment follows, no more input data
	    GO TO 2
	  ELSE IF (VAR(i:i) .EQ. ';') THEN
C  End of command
	    j=j+1
	    k=1
	    i=i+1
	    last_flag=.TRUE.
	    GO TO 1
	  ELSE IF (var(i:i) .EQ. ' ') THEN
	    IF (i.GT.1) THEN
	      l=i-1
	      IF ((var(l:l).NE.' ').AND.(var(l:l).NE.';')) THEN
	        ! this test, when true, signifies the end of
	        ! a 'word' in a command, when there is a space
                ! before next input value.
C  End of command
	        j=j+1
	        k=1
	        i=i+1
	        last_flag=.TRUE.
	        GO TO 1
	      END IF
	    END IF
C  Go to next character
	    i=i+1
	    GO TO 1
	  ELSE
	    var_store(j,pstack)(k:k)=var(i:i)
	    k=k+1
	    i=i+1
	    last_flag=.FALSE.
	    GO TO 1
	  END IF
	END IF
    2   CONTINUE

	IF (last_flag) THEN
	  j=j-1
	END IF

	if (.false.) then
	print*,'i,j, pstack,last_flag=', i,j, pstack,last_flag
        print*,'var_store(j-1,pstack)= ',var_store(j-1,pstack)
        print*,'var_store(j,pstack)= ',var_store(j,pstack)
	end if

	read_len(pstack)=j  ! number of commands read
	read_cur(pstack)=1

	RETURN

c sws
c current journal file empty
  100   CALL MB_LOH(2)
	  IF (pstack .GT. 2) THEN
c still have more .jou files in file stack
c pop the file stack
	     pstack = pstack-1
	     ustack = ustack-1
c now check current read_len
c IF empty, try to fill it
	     IF ( read_len(pstack) .EQ. 0 ) GO TO 10
c ELSE CONTINUE
	  ELSE IF (pstack .EQ. 2) THEN
c no more .jou files in the file stack
c pop the file stack to initial values
	     pstack = pstack-1
	     ustack = ustack-1
c go back to console
	     bat_task=1
	  ELSE
c shouldn't get here
	     WRITE(*,*) '[parse_loh]',
     &	        ' error, pstack = ',pstack, ' bat_task = ',bat_task
c go back to console
	     bat_task=1
	  END IF

c STACK empty
 200    CONTINUE

	  GO TO 98
c	  RETURN

c sws    error RETURN in reading input
 300    CONTINUE
	  WRITE(*,*) '[parse_loh]',
     &	     ' file error during READ, bat_task = ',bat_task
	  RETURN

	END SUBROUTINE PARSE_LOH_SAVE

C***********************************************************************

	SUBROUTINE PROMPT(COUT,TEXT)
c sws
	use lohpars_mod
c sws
	CHARACTER*(*) TEXT
	CHARACTER(len=MacosCharLen) :: COUT
	!CHARACTER*132 COUT
	INTEGER ICLEN

 5    FORMAT(' ',A,$)
 6    FORMAT(A)

c sws
c   jump here if we need a new line from the console
 100  CONTINUE

#ifndef READLINE_LIBRARY
	WRITE(*,5) TEXT
#endif

c sws
c   jump here if we need a new line from a batch file
 200	CONTINUE
	CVAR=' '
	CALL READ_LOH(0)

	IF (CVAR.EQ.' ') THEN
	   IF ( bat_task .EQ. 1 ) GO TO 100
	   GO TO 200
	END IF

	COUT(1:ICLEN(CVAR))=CVAR(1:ICLEN(CVAR))

	RETURN
	END SUBROUTINE PROMPT

C***********************************************************************

	SUBROUTINE CACCEPT(COUT,CDEF,TEXT)
	use lohpars_mod

C       CHARACTER INPUT ROUTINE CACCEPT
C       TYPES OUT TEXT
C       IF INPUT IS CR, CVAR=CDEF
C       OTHERWISE, CVAR=INPUT

	CHARACTER*(*) TEXT,CDEF
	CHARACTER(len=*) COUT
	INTEGER ICLEN

 1	FORMAT(' ',A,' [',A,']: ',$)
 4	FORMAT(9A)
 5	FORMAT(' ',A,$)

 100	CONTINUE
	CVAR=' '
	COUT=' '

c if cdef is blank, read till we get something
	IF (CDEF.EQ.' ')THEN
	   WRITE(*,5) TEXT
	   CALL READ_LOH(1)
	   IF(CVAR.EQ.' ') GO TO 100
	ELSE
c show default, use it if we get a blank
	   WRITE(*,1)TEXT,CDEF
	   CALL READ_LOH(1)
	   IF(CVAR.EQ.' ') THEN
	      CVAR=CDEF
	      IF (JOU_TASK .EQ. 2) THEN
c we're building a journal file, write out the default
	         WRITE(27,*) '   ',CDEF
	      END IF
	   END IF
	END IF
	COUT(1:ICLEN(CVAR))=CVAR(1:ICLEN(CVAR))

	RETURN
	END SUBROUTINE CACCEPT

C***********************************************************************

	SUBROUTINE DACCEPT(DVAR,DDEF,N,TEXT)
	use lohpars_mod

C       DOUBLE-PRECISION N-VECTOR INPUT ROUTINE DACCEPT
C       TYPES OUT TEXT AND DEFAULT VALUES IN PROMPT LINE
C       IF INPUT IS CR, VAR=DEF
C       OTHERWISE, VAR=NPUT

	LOGICAL ifDefault
	CHARACTER*(*) TEXT
	INTEGER I,N,ABSN
	REAL*8 DVAR(N),DDEF(N)

 100	FORMAT(' ',A,$)
 101	FORMAT(1P,6D18.10)

C  Write prompt line:
 1	CONTINUE
	WRITE(*,100)TEXT
	ABSN=ABS(N)
	IF (N.LT.0) THEN
	  WRITE(*,100)'[no default]: '
	  ifDefault=.FALSE.
	ELSE
	  CALL DScreenDef(DDEF,ABSN)
	  ifDefault=.TRUE.
	END IF

C  Read response:
	CALL READ_LOH(2)
	IF (CVAR.NE.' ')THEN
	   READ(CVAR,*,ERR=1)(DVAR(I),I=1,ABSN)
	ELSE IF (ifDefault) THEN
C    No data entered, use the default
	   DO I=1,ABSN
	      DVAR(I)=DDEF(I)
  	   END DO
C    Write to read buffer
	   WRITE(CVAR,101,ERR=1)(DDEF(I),I=1,ABSN)
C    If building a journal file, write out the default
	   IF (JOU_TASK .EQ. 2) THEN
	      CALL DPrintJou(DDEF,ABSN,27)
	   END IF
	ELSE
	   GO TO 1
	END IF

	RETURN
	END SUBROUTINE DACCEPT

	SUBROUTINE DACCEPT_S(DVAR,DDEF,N,TEXT)
        CHARACTER*(*) TEXT
        REAL*8 DVAR,DDEF, dvar_arr(1), ddef_arr(1)
	INTEGER :: N
        dvar_arr(1) = DVAR; ddef_arr(1) = DDEF
        CALL DACCEPT(dvar_arr,ddef_arr,N,TEXT)
        DVAR = dvar_arr(1)
        END SUBROUTINE DACCEPT_S


C***********************************************************************

	SUBROUTINE DPrintJou(DVEC,N,FC)
	INTEGER I,N,FC
	REAL*8 DVEC(9)
 100	FORMAT(1P,4X,D18.10)
 101	FORMAT(1P,4X,D18.10)
 102	FORMAT(1P,4X,D18.10,$)
 103	FORMAT(1P,4X,D18.10,$)
 104	FORMAT(1P,',',D18.10,$)
 105	FORMAT(1P,',',D18.10,$)
 106	FORMAT(' ')
	IF (N.EQ.1) THEN
	  IF (DVEC(1).GE.0d0) THEN
	    WRITE(FC,100)DVEC(1)
	  ELSE
	    WRITE(FC,101)DVEC(1)
	  END IF
	  RETURN
	ELSE
	  IF (DVEC(1).GE.0d0) THEN
	    WRITE(FC,102)DVEC(1)
	  ELSE
	    WRITE(FC,103)DVEC(1)
	  END IF
	  DO I=2,N
	    IF (DVEC(1).GE.0d0) THEN
	      WRITE(FC,104)DVEC(I)
	    ELSE
	      WRITE(FC,105)DVEC(I)
	    END IF
	  END DO
	  WRITE(FC,106)
	  RETURN
	END IF
	END SUBROUTINE DPrintJou

C***********************************************************************

	SUBROUTINE DScreenDef(DVEC,N)
	CHARACTER*23 CNUM
	INTEGER I,N,LEN
	REAL*8 DVEC(N)
 100	FORMAT(' [',$)
 101	FORMAT(A,',',$)
 102	FORMAT(A,']:',$)
	WRITE(*,100)
	IF (N.EQ.1) THEN
	  CALL D2Screen(DVEC(1),CNUM,LEN)
	  WRITE(*,102)CNUM(1:LEN)
	ELSE
	  DO I=1,(N-1)
	    CALL D2Screen(DVEC(I),CNUM,LEN)
	    WRITE(*,101)CNUM(1:LEN)
	  END DO
	  CALL D2Screen(DVEC(N),CNUM,LEN)
	  WRITE(*,102)CNUM(1:LEN)
	END IF
	RETURN
	END SUBROUTINE DScreenDef

C***********************************************************************

	SUBROUTINE D2Screen(D,C,LEN)

	CHARACTER*23 C,C1
	INTEGER nstart,nend,LEN
	REAL*8 D

 100	FORMAT(1P,G14.6)
 101	FORMAT(G14.6)
	IF ((D.LT.0.1d0).OR.(D.GT.1d6)) THEN
	  WRITE(C1,100)D
	ELSE
	  WRITE(C1,101)D
	END IF

C  Find first non-blank
	nstart=0
 1      CONTINUE
	  nstart=nstart+1
	  IF ((C1(nstart:nstart).EQ.' ').AND.(nstart.LT.23)) GO TO 1

C  Find last non-blank
	nend=nstart
 2      CONTINUE
	  nend=nend+1
	  IF ((C1(nend:nend).NE.' ').AND.(nend.LT.23)) GO TO 2
	nend=nend-1

C  Write to character output
	C=C1(nstart:nend)
	LEN=nend-nstart+1
	RETURN
	END SUBROUTINE D2Screen

C***********************************************************************

	SUBROUTINE RACCEPT(RVAR,RDEF,N,TEXT)
C       DOUBLE-PRECISION N-VECTOR INPUT ROUTINE DACCEPT
C       TYPES OUT TEXT AND DEFAULT VALUES IN PROMPT LINE
C       IF INPUT IS CR, VAR=DEF
C       OTHERWISE, VAR=NPUT

	use lohpars_mod

	LOGICAL ifDefault
	CHARACTER*(*) TEXT
	INTEGER I,N,ABSN
	SREAL RVAR(N),RDEF(N)


 100	FORMAT(' ',A,$)
 101	FORMAT(1P,6E16.8)

C  Write prompt line:
 1	CONTINUE
	WRITE(*,100)TEXT
	ABSN=ABS(N)
	IF (N.LT.0) THEN
	  WRITE(*,100)'[no default]: '
	  ifDefault=.FALSE.
	ELSE
	  CALL RScreenDef(RDEF,ABSN)
	  ifDefault=.TRUE.
	END IF

C  Read response:
	CALL READ_LOH(3)
	IF (CVAR.NE.' ')THEN
	   READ(CVAR,*,ERR=1)(RVAR(I),I=1,ABSN)
	ELSE IF (ifDefault) THEN
C    No data entered, use the default
	   DO 9 I=1,ABSN
	      RVAR(I)=RDEF(I)
 9	   CONTINUE
C    Write to read buffer
	   WRITE(CVAR,101,ERR=1)(RDEF(I),I=1,ABSN)
C    If building a journal file, write out the default
	   IF (JOU_TASK .EQ. 2) THEN
	      CALL RPrintJou(RDEF,ABSN,27)
	   END IF
	ELSE
	   GO TO 1
	END IF

	RETURN
	END SUBROUTINE RACCEPT

	SUBROUTINE RACCEPT_S(RVAR,RDEF,N,TEXT)
        CHARACTER*(*) TEXT
        SREAL RVAR,RDEF, rvar_arr(1), rdef_arr(1)
	INTEGER :: N
        rvar_arr(1) = RVAR; rdef_arr(1) = RDEF
        CALL RACCEPT(rvar_arr,rdef_arr,N,TEXT)
        RVAR = rvar_arr(1)
        END SUBROUTINE RACCEPT_S

C***********************************************************************

	SUBROUTINE RPrintJou(RVEC,N,FC)
	INTEGER I,N,FC
	SREAL RVEC(9)
 100	FORMAT(1P,4X,E16.8)
 101	FORMAT(1P,4X,E16.8)
 102	FORMAT(1P,4X,E16.8,$)
 103	FORMAT(1P,4X,E16.8,$)
 104	FORMAT(1P,',',E16.8,$)
 105	FORMAT(1P,',',E16.8,$)
 106	FORMAT(' ')
	IF (N.EQ.1) THEN
	  IF (RVEC(1).GE.0E0) THEN
	    WRITE(FC,100)RVEC(1)
	  ELSE
	    WRITE(FC,101)RVEC(1)
	  END IF
	ELSE
	  IF (RVEC(1).GE.0E0) THEN
	    WRITE(FC,102)RVEC(1)
	  ELSE
	    WRITE(FC,103)RVEC(1)
	  END IF
	  DO I=2,N
	    IF (RVEC(1).GE.0E0) THEN
	      WRITE(FC,104)RVEC(I)
	    ELSE
	      WRITE(FC,105)RVEC(I)
	    END IF
	  END DO
	  WRITE(FC,106)
	END IF
	RETURN
	END SUBROUTINE RPrintJou

C***********************************************************************

	SUBROUTINE RScreenDef(RVEC,N)
	CHARACTER*23 CNUM
	INTEGER I,N,LEN
	SREAL RVEC(9)
 100	FORMAT(' [',$)
 101	FORMAT(A,',',$)
 102	FORMAT(A,']:',$)
	WRITE(*,100)
	IF (N.EQ.1) THEN
	  CALL R2Screen(RVEC(1),CNUM,LEN)
	  WRITE(*,102)CNUM(1:LEN)
	ELSE
	  DO I=1,(N-1)
	    CALL R2Screen(RVEC(I),CNUM,LEN)
	    WRITE(*,101)CNUM(1:LEN)
	  END DO
	  CALL R2Screen(RVEC(N),CNUM,LEN)
	  WRITE(*,102)CNUM(1:LEN)
	END IF
	RETURN
	END SUBROUTINE RScreenDef

C***********************************************************************

	SUBROUTINE R2Screen(R,C,LEN)

	CHARACTER*23 C,C1
	INTEGER nstart,nend,LEN
	SREAL R

 100	FORMAT(1P,G14.6)
 101	FORMAT(G14.6)
	IF ((R.LT.0.1d0).OR.(R.GT.1e6)) THEN
	  WRITE(C1,100)R
	ELSE
	  WRITE(C1,101)R
	END IF

C  Find first non-blank
	nstart=0
 1      CONTINUE
	  nstart=nstart+1
	  IF ((C1(nstart:nstart).EQ.' ').AND.(nstart.LT.23)) GO TO 1

C  Find last non-blank
	nend=nstart
 2      CONTINUE
	  nend=nend+1
	  IF ((C1(nend:nend).NE.' ').AND.(nend.LT.23)) GO TO 2
	nend=nend-1

C  Write to character output
	C=C1(nstart:nend)
	LEN=nend-nstart+1
	RETURN
	END SUBROUTINE R2Screen

C***********************************************************************

	SUBROUTINE IACCEPT(IVAR,IDEF,N,TEXT)
C       INTEGER N-VECTOR INPUT ROUTINE IACCEPT
C       TYPES OUT TEXT AND DEFAULT VALUES IN PROMPT LINE
C       IF INPUT IS CR, VAR=DEF
C       OTHERWISE, VAR=INPUT

	use lohpars_mod

	LOGICAL ifDefault
	CHARACTER*(*) TEXT
	INTEGER I, N, ABSN, IVAR(N), IDEF(N)


 100	FORMAT(' ',A,$)
 101	FORMAT(8I6)

C  Write prompt line:
 1	CONTINUE
	WRITE(*,100) TEXT
	ABSN=ABS(N)
	IF (N.LT.0) THEN
	  WRITE(*,100)'[no default]: '
	  ifDefault=.FALSE.
	ELSE
	  CALL IScreenDef(IDEF,ABSN)
	  ifDefault=.TRUE.
	END IF

	!print*,'IACCEPT: ifDefault=',ifDefault

C  Read response:
	CALL READ_LOH(4)
	IF (CVAR.NE.' ') THEN
	   READ(CVAR,*,ERR=1)(IVAR(I),I=1,ABSN)
	ELSE IF (ifDefault) THEN
     	   ! No data entered, use the default
	   !print*,'**** IACCEPT: use default ...'
	   DO 9 I=1,ABSN
	      IVAR(I)=IDEF(I)
 9	   CONTINUE
     	   ! Write to read buffer
	   !print*,'B4 WRITE to CVAR Buffer'
	   WRITE(CVAR,101,ERR=1)(IDEF(I),I=1,ABSN)
     	   ! If building a journal file, write out the default
	   IF (JOU_TASK .EQ. 2) THEN
	      CALL IPrintJou(IDEF,ABSN,27)
	   END IF
	ELSE
	   GO TO 1
	END IF

	RETURN
	END SUBROUTINE IACCEPT


	SUBROUTINE IACCEPT_S(IVAR,IDEF,N,TEXT)
        CHARACTER*(*) TEXT
        INTEGER IVAR,IDEF, ivar_arr(1), idef_arr(1)
	INTEGER N
        ivar_arr(1) = IVAR; idef_arr(1) = IDEF
        CALL IACCEPT(ivar_arr,idef_arr,N,TEXT)
        IVAR = ivar_arr(1)
        END SUBROUTINE IACCEPT_S

C***********************************************************************

	SUBROUTINE IPrintJou(IVEC,N,FC)
	CHARACTER*16 CNUM
	INTEGER I,N,LEN,FC,IVEC(9)
 100	FORMAT('    ',$)
 101	FORMAT(A,',',$)
 102	FORMAT(A,' ')
	WRITE(FC,100)
	IF (N.EQ.1) THEN
	  CALL I2Screen(IVEC(1),CNUM,LEN)
	  WRITE(FC,102)CNUM(1:LEN)
	ELSE
	  DO I=1,(N-1)
	    CALL I2Screen(IVEC(I),CNUM,LEN)
	    WRITE(FC,101)CNUM(1:LEN)
	  END DO
	  CALL I2Screen(IVEC(N),CNUM,LEN)
	  WRITE(FC,102)CNUM(1:LEN)
	END IF
	RETURN
	END SUBROUTINE IPrintJou

C***********************************************************************

	SUBROUTINE IScreenDef(IVEC,N)
	CHARACTER*16 CNUM
	INTEGER I,N,LEN,IVEC(N)
 100	FORMAT(' [',$)
 101	FORMAT(A,',',$)
 102	FORMAT(A,']: ',$)
	WRITE(*,100)
	IF (N.EQ.1) THEN
	  CALL I2Screen(IVEC(1),CNUM,LEN)
	  WRITE(*,102)CNUM(1:LEN)
	ELSE
	  DO I=1,(N-1)
	    CALL I2Screen(IVEC(I),CNUM,LEN)
	    WRITE(*,101)CNUM(1:LEN)
	  END DO
	  CALL I2Screen(IVEC(N),CNUM,LEN)
	  WRITE(*,102)CNUM(1:LEN)
	END IF
	RETURN
	END SUBROUTINE IScreenDef

C***********************************************************************

	SUBROUTINE I2Screen(I,C,LEN)

	CHARACTER*16 C,C1
	INTEGER nstart,nend,LEN,I

 100	FORMAT(I14)

	WRITE(C1,100)I

C  Find first non-blank
	nstart=0
 1      CONTINUE
	  nstart=nstart+1
	  IF ((C1(nstart:nstart).EQ.' ').AND.(nstart.LT.16)) GO TO 1

C  Find last non-blank
	nend=nstart
 2      CONTINUE
	  nend=nend+1
	  IF ((C1(nend:nend).NE.' ').AND.(nend.LT.16)) GO TO 2
	nend=nend-1

C  Write to character output
	C=C1(nstart:nend)
	LEN=nend-nstart+1
	RETURN
	END SUBROUTINE I2Screen

#include "iosub.inc"
	End Module macosio_mod

C***********************************************************************
C****** End file macosio.F *********************************************
C***********************************************************************
