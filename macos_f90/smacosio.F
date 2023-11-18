C***********************************************************************
C	Begin file smacosio.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module smacosio_mod
	use math_mod

	Implicit none

	Logical :: RxStopSet, EltStopSet,RxNoStopSet_Flg
        REAL*8 :: StopOffset(2)

	Contains

	SUBROUTINE MBFile6(npts,ifOld,LOAD_SUCCESS)
	use param_mod
	use src_mod
	use elt_mod
        use surfsub, only: GridInit
 	use cfiles_mod
 	use macos_mod, only: ifStopSet,ifObjImgSet,ObjIntFn,zernUsrOpt,
     &                       zcoType,ZernCircular,ZernAnnular,ZernHex,
     &                       cUnit,ifRegrid
	use traceutil_mod

	use loh_I1
	use loh_I2
	use loh_I3

	use dopt_mod

	IMPLICIT NONE

	Integer :: npts
	Logical :: ifOld

	INTEGER, PARAMETER :: RX_FID = 77
        LOGICAL :: isRxOpen
c	CHARACTER*8 CDEF
	CHARACTER(len=MacosCharLen) :: TypeString
	LOGICAL exist,LCMP,IHSTPM,IHSTSM,ifGridInit,ifAmplInit,
     &          LnOutCordPass,LOAD_SUCCESS
	INTEGER i,j,k,l,iaxis,iElt,iopt,iCoord,IDEF(7),ICLEN,
     &          iSeg,iObs,jGridSrf,eltGrpStrt,eltGrpEnd,grpRange(20),
     &          nValsPerRow,g,jAmplSrf,jMetSrf,
     &          dupElt,nTok,GetNumToken
	REAL*8 invec6(6),invec7(7),DDEF(7)
	REAL*8 fseg,eseg,psiseg(3),Vseg(3),S1,S2

C       Setup for new input scheme

	CHARACTER*4 IND(3)
	CHARACTER(len=MacosCharLen) :: STRING
 	!CHARACTER*180 STRING
	CHARACTER(MacosVarNamLen) :: VAR_NAM
        CHARACTER(len=MacosValLen) :: VALUE
	!CHARACTER*36 VAR_NAM
 	!CHARACTER*160 VALUE
	LOGICAL ARR_FLG, opt_npts_flg
	LOGICAL ENTR_OK, ELEM_OK, EXIT_OK
	INTEGER ians, nVertex

	REAL*8 D1(9)

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
 130	FORMAT(1P,'Invalid Element Type, Element=',A)
 131	FORMAT(1P,'Invalid Element/Surface Combination '/
     &  'Element=',A,',Surface=',A)
 132	FORMAT(1P,'Invalid Surface Type, Surface=',A)
 133	FORMAT(1P,'Invalid Zernike Type, Surface=',A)
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


C	Initialize constants

	DDEF(1:7)=0D0
	IDEF(1:7)=0D0
 	CWM=1d0; CBM=1d0
	opt_npts_flg=.false.

C	Get file name

 43	CONTINUE
	CALL CACCEPT(filnam,' ','Enter file name: ')
	l=ICLEN(filnam)
	IF ((l.EQ.1).AND.((filnam.EQ.'q').OR.(filnam.EQ.'Q'))) RETURN
	infil=filnam(1:l)//'.in'
	outfil=filnam(1:l)//'.out'
	cmatfil=filnam(1:l)//'.cmat'
	plotfil=filnam(1:l)//'.plot'
	mfil=filnam(1:l)//'.m'
	pertfil=filnam(1:l)//'.pert'
	L=ICLEN(infil)
	INQUIRE (file=infil(1:L),exist=exist)
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
	  iopt=0   ! ifOLD = .false. - jzlou
	END IF

 29	CONTINUE
	EltName(0)='InputRay'
	PropType(0)=1
	zElt(0)=1d19

C
C	Read existing input file
C
	IF (IOPT.EQ.1) THEN
	  OPEN (UNIT=RX_FID, FILE=infil(1:L), STATUS='OLD', ACTION='READ')
	  REWIND(RX_FID)

#include "msmacosio.inc"
   70     CONTINUE
	  If (.not.LOAD_SUCCESS) Then
            WRITE(*,*)' ** Rx input file load failed'
            CLOSE(UNIT=RX_FID); RETURN
          End If

	  CALL LOHIN3()

	  READ(VALUE,*)nOutCord
	  nOutCord_FLG=.TRUE.

	  READ(RX_FID, 123)STRING
	  CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,RX_FID)

	  IF (LCMP(VAR_NAM,'Tout',4)) THEN
	    READ(VALUE,*)(Tout(1,j),j=1,7)
	    IF (nOutCord .GT. 1) THEN
	      DO 71, k=2,nOutCord
	 	READ(RX_FID, 123)STRING
	        READ(STRING,*)(Tout(k,j),j=1,7)
   71	      CONTINUE
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
	    WRITE(*,121)mpts
	    nGridpts=mpts
	    npts=nGridpts-1
	    if (.not.opt_npts_flg) then
	      opt_npts=nGridpts/2-1
	      opt_npts_flg=.true.
	    end if
	  END IF
	  IF (nElt.GT.mElt) THEN
	    WRITE(*,*) ' Too many elements.'
	    WRITE(*,*) '  MODIFY or use another version. '
	  END IF
	  IF ((GridType.GE.3).AND.(nSeg.GT.mElt)) THEN
	    WRITE(*,*) ' Too many segments.'
	    WRITE(*,*) '  MODIFY or use another version. '
	  END IF

C	Create new input file
C
	ELSE
	  WRITE(*,*)' SMACOS does not support the NEW command.'
	END IF


        ! Close input file
        INQUIRE(UNIT=RX_FID, OPENED=isRxOpen)
        IF (isRxOpen) CLOSE(UNIT=RX_FID)


        ! Finalise Configuration
	npts=nGridpts-1

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

C           print*,'SegCoord(1:2,iSeg),ApType,ApVec(1:4)=',
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
        End SUBROUTINE SetSegAp


C***********************************************************************
C***********************************************************************

	SUBROUTINE READ_LOH

C       THIS ROUTINE READS THE COMMAND STRING, PARSES IT, AND
C       PASSES THE APPROPRIATE COMMAND BACK TO THE ROUTINE
C       THAT CALLED IT.  A MAXIMUM OF 10 ENTRIES MAY BE PUT ON
C       ONE LINE.  ENTRIES ARE DELIMITED BY A SEMICOLON.  COMMENTS
C       ARE DELIMITED BY A PERCENT SIGN.  A COMMAND MAY BE
C       256 CHARACTERS LONG.

	use lohpars_mod

	IMPLICIT NONE
CSCOMP
c sws
	LOGICAL LCMP
c sws
CSCOMP
	INTEGER l
	INTEGER ICLEN
	LOGICAL ifEcho
	COMMON /SCIO/ifEcho

 5	FORMAT(' ',A,$)
 6	FORMAT(' ',A)

C       Are there inputs that have been read, but not used?

	IF (read_len(pstack) .EQ. 0) THEN
	   CALL PARSE_LOH()
	END IF

c sws
c   only do this if we got something
	IF (read_len(pstack) .NE. 0) THEN
cdbg	   CVAR=var_store(read_cur(pstack),pstack)
	   CVAR(1:cvarlen) = var_store(read_cur(pstack),pstack)

C       Write to screen
	   IF ( (read_cur(pstack) .GT. 1) .OR. (bat_task .EQ. 2) )THEN
	      l=ICLEN(CVAR)
	      IF (l.LE.0) l=1
	      IF (ifEcho) WRITE(*,6) CVAR(1:l)
	   END IF

C  Write out journal data to file, if in journal mode
c sws
	   IF (jou_task .EQ. 2) THEN
	      IF (.NOT. LCMP(CVAR,'jou',3)) THEN
	         WRITE (27,6) CVAR(1:ICLEN(CVAR))
	      END IF
	   END IF

	ELSE
c got nothing, RETURN a blank
cdbg	   CVAR = ' '
	   CVAR(1:cvarlen) = ' '

	END IF

C       Is this the last command saved?

	IF (read_cur(pstack) .EQ. read_len(pstack)) THEN
	   read_len(pstack)=0
	   read_cur(pstack)=0
	ELSE
	   read_cur(pstack)=read_cur(pstack)+1
	END IF

	RETURN
	END SUBROUTINE READ_LOH

C***********************************************************************
C***********************************************************************

	SUBROUTINE PARSE_LOH()

C       THIS ROUTINE DOES THE ACTUAL PARSING.  A MAXIMUM OF 10
C       ENTRIES MAY BE PUT ON ONE LINE.  ENTRIES ARE DELIMITED
C       BY A SEMICOLON or a blank.  COMMENTS ARE DELIMITED BY A PERCENT
C       SIGN. A COMMAND LINE MAY BE 1000 CHARACTERS LONG.
C
	use lohpars_mod

	IMPLICIT NONE

CSCOMP
	CHARACTER*1000 VAR
	INTEGER ls
	CHARACTER*1000 STACK
	COMMON /SCCOM/ls,STACK
CSCOMP
	  INTEGER i,j,k,l
	INTEGER len_max
	  LOGICAL last_flag

C    IF last flag = true, THEN the semi colon is the last CHARACTER
C    on the line, (besides spaces)


    6   FORMAT(A)

c sws
c  added for jump from poping the stack
 10     CONTINUE

CSCOMP
c	  len_max=ls
	  IF (bat_task .EQ. 1) THEN
	    len_max=ls
	  ELSE
	    len_max=130
	  END IF

CSCOMP
	  DO i=1,10
	    var_store(i,pstack)=' '
	  END DO

C       READ new statement to be parsed

   98     CONTINUE

	IF (bat_task .EQ. 1) THEN
c sws
 99	  READ(STACK,6,END=200, ERR=300) VAR
	 ELSE
	     READ(ustack,6,END=100, ERR=300) VAR
	END IF
c sws

C       i is counter for location in array VAR, j is counter for
C       number of commands in line, k is counter for var_store

	  i=1
	  j=1
	  k=1
	  last_flag=.TRUE.

C  Check for % sign, semicolon, blank space...
    1   CONTINUE
	  IF (i .LE. len_max) THEN
	    IF (var(i:i) .EQ. '%') THEN
C           Comment follows, no more input DATA
csws        add a blank
               var_store(j,pstack)(k:k) = ' '
               k = k + 1
	      GO TO 2
	    ELSE IF (var(i:i) .EQ. ';') THEN
C  END of command
               var_store(j,pstack)(k:k) = ' '
c  now prepare for next command
	      j=j+1
	      k=1
	      i=i+1
	      last_flag=.TRUE.
	      GO TO 1
	    ELSE IF (var(i:i) .EQ. ' ') THEN
	      IF (i.GT.1) THEN
	        l=i-1
	        IF ((var(l:l).NE.' ').AND.(var(l:l).NE.';')) THEN
C  END of command
                   var_store(j,pstack)(k:k) = ' '
c  now prepare for next command
	          j=j+1
	          k=1
	          i=i+1
	          last_flag=.TRUE.
	          GO TO 1
	        END IF
	      END IF
C  GO TO next CHARACTER
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
	read_len(pstack)=j
	read_cur(pstack)=1

	RETURN

c sws
c current journal file empty
 100	CALL MB_LOH(2)
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
c go back to STACK
	     bat_task=1
c scomp stack holds only one command
c stack command must have been exe
c we've done that, so just return
	  GO TO 200
	ELSE
c shouldn't get here
	  WRITE(*,*) '[parse_loh]',
     &	  ' error, pstack = ',pstack, ' bat_task = ',bat_task
c go back to console
	  bat_task=1
	END IF

c now try to read again, without initializing
	GO TO 98

c input STACK empty, quit
 200	CONTINUE

	RETURN

c sws    error RETURN in reading input
 300	CONTINUE
	WRITE(*,*) '[parse_loh]',
     &	' file error during READ, bat_task = ',bat_task
	RETURN

	END SUBROUTINE PARSE_LOH

C***********************************************************************

	SUBROUTINE PROMPT(COUT,TEXT)
	use lohpars_mod

	CHARACTER*(*) TEXT
	CHARACTER(len=MacosCharLen) :: COUT

	INTEGER ICLEN,L
	LOGICAL ifEcho
	COMMON /SCIO/ifEcho

  5     FORMAT(' ',A,$)
  6	FORMAT(A)

c sws   old version barfs...
c	COUT=' '

	L = ICLEN(COUT)
	COUT(1:L) = ' '

c sws
c   jump here IF we need a new line from the console
100     CONTINUE
	IF (ifEcho) WRITE(*,5)TEXT

c sws
c   jump here IF we need a new line from a batch file
 200	CONTINUE
csws	CVAR=' '
	CVAR(1:cvarlen) = ' '
	CALL READ_LOH

CSCOMP
	IF ((CVAR(1:6).NE.'RETURN').AND.(ifEcho))
     &	  WRITE (*,6) CVAR(1:ICLEN(CVAR))
CSCOMP

	IF(CVAR.EQ.' ') THEN
	   IF ( bat_task .EQ. 1 ) GO TO 100
	   GO TO 200
	END IF

c sws seem to need restriction to keep from clobbering memory
c	COUT = CVAR
c	COUT(1:64) = CVAR(1:64)
c	COUT(1) = CVAR
c	COUT(1:132) = CVAR(1:132)
	L = ICLEN(CVAR)
	if (L .lt. MacosCharLen) then
cloh	   COUT(1:L) = CVAR(1:L)
c sws	   COUT = CVAR(1:L)
	   COUT(1:L) = CVAR(1:L)

	ELSE
	   L=20
cloh	   COUT(1:L) = CVAR(1:L)
csws	   COUT = CVAR(1:L)
	   COUT(1:L) = CVAR(1:L)
	END IF

	COUT(L+1:L+1)=' '

	RETURN
	END SUBROUTINE PROMPT

C***********************************************************************

	SUBROUTINE CACCEPT(COUT,CDEF,TEXT)
C       CHARACTER INPUT ROUTINE CACCEPT
C       TYPES OUT TEXT
C       IF INPUT IS CR, CVAR=CDEF
C       OTHERWISE, CVAR=INPUT

	use lohpars_mod

CSCOMP
	CHARACTER*(*) TEXT,CDEF
	CHARACTER(len=*) :: COUT
	LOGICAL ifEcho
	INTEGER L,ICLEN
	COMMON /SCIO/ifEcho

 1	FORMAT(' ',A,' (',A,'): ',$)
 4	FORMAT(9A)
 5	FORMAT(' ',A,$)

	COUT(1:MacosCharLen)=' '

 100  CONTINUE
	CVAR(1:cvarlen)=' '

c if cdef is blank, read till we get something
	IF (CDEF.EQ.' ')THEN
	   IF (ifEcho) WRITE(*,5)TEXT
	   CALL READ_LOH
	   IF(CVAR.EQ.' ') GO TO 100
	ELSE
c show default, use it if we get a blank
	   IF (ifEcho) WRITE(*,1)TEXT,CDEF
	   CALL READ_LOH
	   IF(CVAR.EQ.' ') THEN
	    CVAR=CDEF
	    IF (JOU_TASK .EQ. 2) THEN
c we're building a journal file, write out the default
	         WRITE(27,*) CDEF
	    END IF
	   END IF
	END IF

c sws
c	COUT = CVAR
c	COUT(1:64) = CVAR(1:64)
c	COUT(1:132) = CVAR(1:132)
	L = ICLEN(CVAR)
	IF (L .lt. MacosCharLen) THEN
	   COUT = CVAR(1:L)
	ELSE
	   L=20
	   COUT = CVAR(1:L)
	END IF

	COUT(L+1:L+1)=' '
	RETURN
	END SUBROUTINE CACCEPT

C***********************************************************************

	SUBROUTINE DACCEPT(DVAR,DDEF,N,TEXT)
C       DOUBLE-PRECISION N-VECTOR INPUT ROUTINE DACCEPT
C       TYPES OUT TEXT
C       IF INPUT IS CR, DVAR=DDEF
C       OTHERWISE, DVAR=INPUT

	use lohpars_mod

	INTEGER N,I
	CHARACTER*(*) TEXT
	REAL*8 DVAR(N),DDEF(N)
	LOGICAL ifEcho
	COMMON /SCIO/ifEcho

 2	FORMAT(' ',A,$)
 3	FORMAT(BN,9F31.0)
 4	FORMAT(A)
 5	FORMAT(' ',F15.8,$)
 6	FORMAT(1P,6D18.10)

 100  IF (ifEcho) WRITE(*,2)TEXT

	CALL READ_LOH
	!print*,'**** CVAR =',CVAR

	IF (CVAR.NE.' ')THEN
c        got something
	   READ(CVAR,*,ERR=100)(DVAR(I),I=1,N)
	ELSE
c         no data entered, use the default
	  DO I=1,N
	    DVAR(I)=DDEF(I)
          END DO
c            write to buffer
	  WRITE(CVAR,6,ERR=100)(DDEF(I),I=1,N)
c            write to screen
	  IF (ifEcho) WRITE(*,*) ' [', (DDEF(I),I=1,N), ']'
	  IF (JOU_TASK .EQ. 2) THEN
c we're building a journal file, write out the default
	     WRITE(27,*)(DDEF(I),I=1,N)
	   END IF
	END IF
	RETURN
	END SUBROUTINE DACCEPT

C ************************************************************

	SUBROUTINE DACCEPT_S(DVAR,DDEF,N,TEXT)
        CHARACTER*(*) TEXT
	INTEGER N
        REAL*8 DVAR,DDEF, dvar_arr(1), ddef_arr(1)
        dvar_arr(1) = DVAR; ddef_arr(1) = DDEF
        CALL DACCEPT(dvar_arr,ddef_arr,1,TEXT)
        DVAR = dvar_arr(1)
        END SUBROUTINE DACCEPT_S


C***********************************************************************

	SUBROUTINE RACCEPT(RVAR,RDEF,N,TEXT)
C       SREAL N-VECTOR INPUT ROUTINE DACCEPT
C       TYPES OUT TEXT
C       IF INPUT IS CR, RVAR=RDEF
C       OTHERWISE, RVAR=INPUT

	use lohpars_mod

	INTEGER N,I
	CHARACTER*(*) TEXT
	SREAL RVAR(N),RDEF(N)
	LOGICAL ifEcho
	COMMON /SCIO/ifEcho

 2	FORMAT(' ',A,$)
 3	FORMAT(BN,F31.0)
 4	FORMAT(A)
 5	FORMAT(' ',F15.8,$)
 6	FORMAT(1P,6E16.8)

 100  IF (ifEcho) WRITE(*,2)TEXT

	CALL READ_LOH

	IF (CVAR(1:1) .NE. ' ')THEN
c        got something
	   READ(CVAR,*,ERR=100) (RVAR(I),I=1,N)
	ELSE
c        no data entered, use the default
	   DO 14 I=1,N
	      RVAR(I)=RDEF(I)
 14      CONTINUE
c write to buffer
	   WRITE(CVAR,6,ERR=100)(RDEF(I),I=1,N)
c write to screen
	   IF (ifEcho) WRITE(6,*) ' [', (RDEF(I),I=1,N), ']'
	   IF (JOU_TASK .EQ. 2) THEN
c we're building a journal file, write out the default
	      WRITE(27,*)(RDEF(I),I=1,N)
	   END IF
	END IF
	RETURN
	END SUBROUTINE RACCEPT

	SUBROUTINE RACCEPT_S(RVAR,RDEF,N,TEXT)
	INTEGER N
        CHARACTER*(*) TEXT
        SREAL RVAR,RDEF, rvar_arr(1), rdef_arr(1)
        rvar_arr(1) = RVAR; rdef_arr(1) = RDEF
        CALL RACCEPT(rvar_arr,rdef_arr,1,TEXT)
        RVAR = rvar_arr(1)
        END SUBROUTINE RACCEPT_S

C*****************************************************************************

	SUBROUTINE IACCEPT(IVAR,IDEF,N,TEXT)
C       INTEGER N-VECTOR INPUT ROUTINE IACCEPT
C       TYPES OUT TEXT
C       IF INPUT IS CR, IVAR=IDEF,  i.e. IDEF = default values for IVAR, jzlou
C       OTHERWISE, IVAR=INPUT

	use lohpars_mod

	CHARACTER*(*) TEXT
	INTEGER N,IVAR(N),IDEF(N),I
	LOGICAL ifEcho
	COMMON /SCIO/ifEcho

 1	FORMAT(1X,A,'[',G23.16,'] ')
 2	FORMAT(' ',A,$)
 3	FORMAT(BN,9I8)
 4	FORMAT(A)
 5	FORMAT(' ',I8,$)
 6	FORMAT(6I8)

 100  IF (ifEcho) WRITE(*,2)TEXT

	CALL READ_LOH

	IF (CVAR.NE.' ')THEN
c got something
	   READ(CVAR,*,ERR=100)(IVAR(I),I=1,N)
	ELSE
c        no data entered, use the default
	   DO 14 I=1,N
	      IVAR(I)=IDEF(I)
 14      CONTINUE
c           write to buffer
	   WRITE(CVAR,6,ERR=100)(IDEF(I),I=1,N)
c           write to screen
	   WRITE(6,*) ' [',(IDEF(I),I=1,N),']'
	   IF (JOU_TASK .EQ. 2) THEN
c we're building a journal file, write out the default
	      WRITE(27,6)(IDEF(I),I=1,N)
	   END IF
	END IF
	RETURN
	END SUBROUTINE IACCEPT

	SUBROUTINE IACCEPT_S(IVAR,IDEF,N,TEXT)
	INTEGER N
        CHARACTER*(*) TEXT
        INTEGER IVAR,IDEF, ivar_arr(1), idef_arr(1)
	ivar_arr(1) = IVAR; idef_arr(1) = IDEF
	CALL IACCEPT(ivar_arr,idef_arr,1,TEXT)
	IVAR = ivar_arr(1)
	END SUBROUTINE IACCEPT_S

#include "iosub.inc"
	End Module smacosio_mod

C***********************************************************************
C****** END file scompio.F *********************************************
C***********************************************************************
