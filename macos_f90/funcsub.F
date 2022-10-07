C***********************************************************************
C	Begin file funcsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2008, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	SUBROUTINE CPERTURB(iElt)
	use param_mod
	use src_mod
	use elt_mod
	use math_mod
	use loh_I1, only : SegXgrid_FLG

#ifdef CSMACOS
	use smacosio_mod
#else
	use macosio_mod
#endif
	use traceutil_mod, only : CRIncidPosNom,ShowMetData_FLG

	IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: ANS

	INTEGER :: i,j,k,iElt
	LOGICAL :: LLocCoord,LCMP

	REAL*8 :: th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9),DDEF(7),PV(6)

 400	FORMAT(' Computing new perturbed system parameters')

	CALL DZERO(DDEF,7)
	LLocCoord=.true.  ! default

C  1. Compute perturbed entrance pupil data

	IF (iElt.EQ.0) THEN  ! source is perturbed
	  CALL DZERO(D1,6)
	  CALL DACCEPT(th,D1,3,
     &       'Enter rotational perturbation vector (x,y,z):')
	  CALL DACCEPT(del,D1,3,
     &       'Enter translational perturbation vector (x,y,z):')

c	  WRITE(*,400)

          If (SrcLF_FLG) Then
            ! Source local frames defined, so convert source perturbations 
            ! from local frames to global frames
            If (CalcSrcZFrame_FLG) Then
              ! Compute zFrame only when it's not defined in Rx
              Call DXPROD(TElt(1:3,3,0),TElt(1:3,1,0),TElt(1:3,2,0))
              Call DUNITIZE(TElt(1:3,3,0))
              Call DXPROD(TElt(1:3,2,0),TElt(1:3,3,0),TElt(1:3,1,0))
              Call DUNITIZE(TElt(1:3,2,0))
              Call DXPROD(TElt(1:3,1,0),TElt(1:3,2,0),TElt(1:3,3,0))
              print*,'****** zFrame =', TElt(1:3,3,0)
              CalcSrcZFrame_FLG=.FALSE.
            End If
            TElt(4:6,4:6,0)=TElt(1:3,1:3,0)
            D1(1:3)=th; D1(4:6)=del
            k=6
            CALL DMPROD(PV,TElt(1,1,0),D1,6,k,1) ! convert to global
            Do i=1,3
              th(i)=PV(i); del(i)=PV(i+3)
            End Do
          End If
          print*,'*** CPERTURB src: th = ',th(1:3),', del =',del(1:3) 

	  CALL Qform(Q,dQ,th)
	  CALL DMPROD(D1,Q,ChfRayDir,3,3,1)
	  CALL DEQUATE(ChfRayDir,D1,3)
	  CALL DUNITIZE(ChfRayDir)

	  CALL DMPROD(D1,Q,xGrid,3,3,1)
	  CALL DEQUATE(xGrid,D1,3)
	  CALL DMPROD(D1,Q,yGrid,3,3,1)
	  CALL DEQUATE(yGrid,D1,3)

	  CALL DSUB(D1,ChfRayPos,StopPos,3)
	  CALL DMPROD(D2,dQ,D1,3,3,1)
	  CALL DADD(D1,ChfRayPos,D2,3)
	  CALL DADD(ChfRayPos,D1,del,3)

	  If (SegXgrid_FLG) Then
            CALL DMPROD(D1,Q,SegXGrid,3,3,1)
            CALL DEQUATE(SegXGrid,D1,3)
            CALL DUNITIZE(SegXGrid)
          End If

C  2. Compute perturbed element data

	ELSE  ! other element (than source) perturbed
          !print*,'******funcsub: iElt, nECoord(iElt) =',iElt,nECoord(iElt) 

	  CALL CACCEPT(ANS,'GLOBAL',
     &      'Enter coordinate system for perturbation '//
     &	    '(ELEMENT or GLOBAL): ')
	  IF (LCMP(ANS,'Y',1).OR.LCMP(ANS,'E',1)) THEN
c sws new version
	      k= nECoord(iElt)
	      IF (k .EQ. 0) THEN
c k = 0, no dof
	         WRITE (*,*) 'No local DOF available!'
	         CALL DZERO(PV,6)
	      ELSE
	         IF (k .LT. 0) THEN
c just global coordinates, continue anyway
	            k = -k
	            LLocCoord=.false.
	         END IF
	         CALL DACCEPT(D1,DDEF,k,
     &	            'Enter perturbation vector in element coordinates:')

	        If (nECoord(iElt)<0) Then
                  PV(1:6)=D1(1:6)  ! assuming global frame
	        Else 
  	          CALL DMPROD(PV,TElt(1,1,iElt),D1,6,k,1) ! convert to global 
	        End If
	      END IF

	  ELSE
	    ! global coord used
	    LLocCoord=.false.
	    CALL DACCEPT(PV(1:3),DDEF,3,
     &        'Enter rotational perturbation vector (x,y,z):')
	    CALL DACCEPT(PV(4:6),DDEF,3,
     &        'Enter translational perturbation vector (x,y,z):')
C	    print*,'***CPERTURB PV(1:6) = ',
C    &         PV(1),PV(2),PV(3),PV(4),PV(5),PV(6)
	  END IF

c	  WRITE(*,400)

	  Do i=1,3
	    j=i+3
	    th(i)=PV(i)
	    del(i)=PV(j)
	  End Do

#if 0
	  if (iElt==3) then
            Open(Unit=13,File='smacos_d1.txt',Status='replace')
            Open(Unit=12,File='smacos_pv.txt',Status='replace')
          end if
          Write(13,*) D1(1:6)
          Write(12,*) PV(1:6)
          if (iElt==3+162-1) then
            Close(12); Close(13)
          end if
#endif

	  CALL Qform(Q,dQ,th)
	  CALL DMPROD(D1,Q,psiElt(1,iElt),3,3,1)
	  CALL DEQUATE(psiElt(1,iElt),D1,3)

	  CALL DSUB(D1,VptElt(1,iElt),RptElt(1,iElt),3)
	  CALL DMPROD(D2,dQ,D1,3,3,1)
	  CALL DADD(D1,VptElt(1,iElt),D2,3) ! = RptElt + Q*(VptElt-RptElt) -jzlou
	  CALL DADD(VptElt(1,iElt),D1,del,3) ! add translation part -jzlou

	  If (iElt==3 .OR. iElt==4) Then
            !print*,'**cpert: RptElt =',RptElt(1:3,iElt)
            !print*,'**cpert: psiElt aft =',psiElt(1:3,iElt)
            !print*,'**cpert: th =',th(1:3)
            !print*,'**cpert: del =',del(1:3)
          End If

	  ! Since RptElt is not affected by the rotation (about itself),
          ! its update is only affected by the translation part of
          ! perturbation.
	  CALL DADD(D1,RptElt(1,iElt),del,3)
	  CALL DEQUATE(RptElt(1,iElt),D1,3)

	  ! Update CRIncidPosNom, added by jzlou
          CALL DSUB(D1,CRIncidPosNom(1,iElt),RptElt(1,iElt),3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,CRIncidPosNom(1,iElt),D2,3)
          CALL DADD(CRIncidPosNom(1,iElt),D1,del,3)

	  If (nMetPos(iElt)>0) Then
	    ! Update all metrology points on the surface
	    Do j=1,nMetPos(iElt)
	      CALL DSUB(D1,SrfMetPos(1,j,iEltToMetSrf(iElt)),
     &                  RptElt(1,iElt),3)
              CALL DMPROD(D2,dQ,D1,3,3,1)
              CALL DADD(D1,SrfMetPos(1,j,iEltToMetSrf(iElt)),D2,3)
              CALL DADD(SrfMetPos(1,j,iEltToMetSrf(iElt)),D1,del,3)
	      If (ShowMetData_FLG) Then
	        print*,'At Elt',iElt,': Updated met point =',
     &            SrfMetPos(1:3,j,iEltToMetSrf(iElt))
	      End If
	    End Do
	  End If


C fixed by Gopal Vasudevan:
	  !IF (nObs(iElt).NE.0) THEN
          ! This IF test is removed by jzlou, since an aperture needs
          ! to be perturbed with the element.
	      CALL DMPROD(D1,Q,xObs(1,iElt),3,3,1)
	      CALL DEQUATE(xObs(1,iElt),D1,3)
	  !END IF

	  IF ((EltID(iElt).EQ.4).OR.(EltID(iElt).EQ.5).OR.
     &	       (EltID(iElt).EQ.13)) THEN
	    CALL DSUB(D1,h1HOE(1,iElt),RptElt(1,iElt),3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
	    CALL DADD(D1,h1HOE(1,iElt),D2,3)
	    CALL DADD(h1HOE(1,iElt),D1,del,3)
	    CALL DSUB(D1,h2HOE(1,iElt),RptElt(1,iElt),3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
	    CALL DADD(D1,h2HOE(1,iElt),D2,3)
	    CALL DADD(h2HOE(1,iElt),D1,del,3)
	  END IF

	  IF (((SrfType(iElt).GE.4).AND.(SrfType(iElt).LE.9))
     &        .OR. SrfType(iElt).EQ.11) THEN  ! added type 11 for AsGrData, 07/2019
	    CALL DSUB(D1,pMon(1,iElt),RptElt(1,iElt),3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
	    CALL DADD(D1,pMon(1,iElt),D2,3)
	    CALL DADD(pMon(1,iElt),D1,del,3)
	    CALL DMPROD(D1,Q,xMon(1,iElt),3,3,1)
	    CALL DEQUATE(xMon(1,iElt),D1,3)
	    CALL DMPROD(D1,Q,yMon(1,iElt),3,3,1)
	    CALL DEQUATE(yMon(1,iElt),D1,3)
	    CALL DMPROD(D1,Q,zMon(1,iElt),3,3,1)
	    CALL DEQUATE(zMon(1,iElt),D1,3)
	  END IF
	END IF

	If (LLocCoord .and. LUpdateTElt_FLG(iElt)) Then
          ! TElt(1:6,1:6,iElt) is updated here
          Do j=1,3
            CALL DMPROD(D1,Q,TElt(1:3,j,iElt),3,3,1)
            CALL DEQUATE(TElt(1:3,j,iElt),D1,3)
            CALL DMPROD(D1,Q,TElt(4:6,3+j,iElt),3,3,1)
            CALL DEQUATE(TElt(4:6,3+j,iElt),D1,3)
          End Do
        End If

#if 1
        ! Perturb linked element, if any
	Do i=1, nElt
          if (LnkElt(i)==iElt) then
	    print*,' **CPERTURB: perturbing linked element: ', i
	    Call LnkEltCPERTURB(i,th,del,LLocCoord,LUpdateTElt_FLG(iElt))
	    exit  ! only one linked element allowed
	  end if
	End Do
#endif
	RETURN
	END SUBROUTINE CPERTURB 


C***********************************************************************
#include "lnk_pert.inc"
C***********************************************************************

	SUBROUTINE CPRead(PertVec,ifPertRead,iPertElt)

C       A routine to read input perturbation data.  The data is
C       assumed to be in the form of transposed u-vectors, starting
C       with element 1 and proceeding to the last element.  Deformable
C       elements are loaded completely (all rays); non-deformable
C       elements are loaded for ray 1 only.

	use param_mod
	use cfiles_mod
	use elt_mod
	use math_mod

#ifdef CSMACOS
	use smacosio_mod
#else
	use macosio_mod
#endif

	IMPLICIT NONE


	LOGICAL ifPertRead,exist
	INTEGER i,j,iaxis,iElt,iPertElt,ICLEN
	REAL*8 invec6(6),th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9)
	SREAL :: PertVec(6,mElt)

 100	FORMAT(' Input file ',A,' not compatible with current system:'/
     &	'  iPertElt=',i3,'  iElt=',i3)

C  Read .pert file

	CALL CACCEPT(pertfil,pertfil,'Enter new .pert-file name:')
 1	CONTINUE
	  INQUIRE (file=pertfil,exist=exist)
	  IF (.NOT.exist) THEN
	    WRITE(*,*)
     &	   ' Input file '//pertfil(1:ICLEN(pertfil))//' does not exist.'
	    CALL CACCEPT(pertfil,'q','Enter new .pert-file name:')
	    IF ((ICLEN(pertfil).EQ.1).AND.(pertfil(1:1).EQ.'q')) THEN
	      ifPertRead=.FALSE.
	      RETURN
	    END IF
	    GO TO 1
	  ELSE
	    OPEN (7,FILE=PERTFIL,STATUS='OLD')
	    DO 4 iElt=1,iPertElt
	      READ(7,*,ERR=7,END=7)invec6
	      DO 3 iaxis=1,6
	        PertVec(iaxis,iElt)=invec6(iaxis)
 3            CONTINUE
 4          CONTINUE
	  END IF
	  ifPertRead=.TRUE.
	  WRITE(*,*)' File '//pertfil(1:ICLEN(pertfil))//' read.'
	  CLOSE(7)

C  Apply perturbations to each element

	DO 6 iElt=1,iPertElt
	  DO 5 i=1,3
	    j=i+3
	    th(i)=PertVec(i,iElt)
	    del(i)=PertVec(j,iElt)
  5	  CONTINUE
	  CALL Qform(Q,dQ,th)
	  CALL DMPROD(D1,Q,psiElt(1,iElt),3,3,1)
	  CALL DEQUATE(psiElt(1,iElt),D1,3)

	  CALL DSUB(D1,VptElt(1,iElt),RptElt(1,iElt),3)
	  CALL DMPROD(D2,dQ,D1,3,3,1)
	  CALL DADD(D1,VptElt(1,iElt),D2,3)
	  CALL DADD(VptElt(1,iElt),D1,del,3)

	  CALL DADD(D1,RptElt(1,iElt),del,3)
	  CALL DEQUATE(RptElt(1,iElt),D1,3)

	  IF (EltID(iElt).EQ.4) THEN
	    CALL DSUB(D1,h1HOE(1,iElt),RptElt(1,iElt),3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
	    CALL DADD(D1,h1HOE(1,iElt),D2,3)
	    CALL DADD(h1HOE(1,iElt),D1,del,3)
	    CALL DSUB(D1,h2HOE(1,iElt),RptElt(1,iElt),3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
	    CALL DADD(D1,h2HOE(1,iElt),D2,3)
	    CALL DADD(h2HOE(1,iElt),D1,del,3)
	  END IF

	  IF ((SrfType(iElt).GE.4).AND.(SrfType(iElt).LE.9)) THEN
	    CALL DSUB(D1,pMon(1,iElt),RptElt(1,iElt),3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
	    CALL DADD(D1,pMon(1,iElt),D2,3)
	    CALL DADD(pMon(1,iElt),D1,del,3)
	    CALL DMPROD(D1,Q,xMon(1,iElt),3,3,1)
	    CALL DEQUATE(xMon(1,iElt),D1,3)
	    CALL DMPROD(D1,Q,yMon(1,iElt),3,3,1)
	    CALL DEQUATE(yMon(1,iElt),D1,3)
	    CALL DMPROD(D1,Q,zMon(1,iElt),3,3,1)
	    CALL DEQUATE(zMon(1,iElt),D1,3)
	  END IF
  6	CONTINUE
	RETURN

C  Here if read failed

 7	CONTINUE
	WRITE(*,100)pertfil(1:ICLEN(pertfil)),iPertElt,iElt
	ifPertRead=.FALSE.
	CLOSE(7)
	RETURN

	END SUBROUTINE CPRead

C***********************************************************************
C***********************************************************************

	SUBROUTINE CEXPORT(ifBuild,npts,GridType,RayMap,SegMap,RayIndex,
     &		  Cmatrix,iPertElt,iCurRayElt,ifDMBuild,iObsOpt)
	use param_mod
	use cfiles_mod
	use elt_mod
	use math_mod

#ifdef CSMACOS
	use smacosio_mod
#else
	use macosio_mod
#endif

	IMPLICIT NONE


	INTEGER i,j,iElt,iNode,iDOF,npts,RayMap(mpts,mpts),
     &	 iRay,jRay,RayIndex(2,mRay),SegMap(mpts,mpts),NCROWS,NCCOLS,
     &	 GridType,nFirstRay,nLastRay,nExpRay,nCMCoords(0:mElt),iObsOpt,
     &	 irow,icol,iPertElt,LENGTH,ICLEN,iExMode,iCurRayElt,eRay,
     &   iIntSrf,totRow,totCol

	SREAL :: Cmatrix(7,mCm,bRay),Celt(49),Cout(49),Tcompact(49),
     &	         TEsp(6,6,mElt),CoutDM(21)

	LOGICAL IfEltCoords,IfInPartials,ifBuild,exist,LCMP,ifDMBuild

	CHARACTER(len=MacosCharLen) :: COMMAND,ANS

C       File setup and initial conditions read

 100	FORMAT('% RayPos 3-vectors for',i8,' rays at Element',i6)
 101    FORMAT(i7,1X,3D23.15)
 102	FORMAT(' RayPos(1:3,',i7,')=[',3D23.15,3H]';)
 103	FORMAT('% RayDir 3-vectors for',i8,' rays at Element',i6)
 104	FORMAT(' RayDir(1:3,',i7,')=[',3D23.15,3H]';)
 105	FORMAT('% CumRayL scalars for',i8,' rays at Element',i6)
 106    FORMAT(i7,1X,3D23.15)
 107	FORMAT(' CumRayL(',i7,')=[',D23.15,3H]';)
 108	FORMAT('% RayIndex 2-vectors for',i8,' rays at Element',i6)
 109    FORMAT(i7,1X,i7,1x,i7)
 110    FORMAT(i7,1X,i7,1x,i7)

 200    FORMAT(/'COMP Export Function'/
     &  '  Prompts for EXPORT mode, then for the data to be exported.',
     &  ' Typing'/'  Quit returns you to the COMP executive.'/
     &	'     Export modes are:'/
     &  '         BINary: data is saved to a Fortran binary file'/
     &  '         TEXt: data is saved to an ASCII text file'/
     &  '         MFILE: data is saved to a Matlab .m-file'/
     &  '         MAT: data is saved to a Matlab .mat-file'/
     &	'     Data available for EXPORT:'/
     &  '         CMAtrix: optical sensitivity matrix in global or',
     &	' local coordinates'/
     &  '         RayIndex: pupil position 2-vectors for each ray'/
     &  '         RayPos: current position 3-vectors for each ray'/
     &  '         RayDir: current ray direction 3-vectors for each ray'/
     &  '         CumRayL: current optical pathlengths for each ray'//
     &  '  Note: Exported data is sensitive to obscuration option.')
 202    FORMAT('Command not understood.  Please reenter...')

C  Set EXPORT mode

	mfil=filnam
	LENGTH=ICLEN(filnam)
  1      CONTINUE
	  CALL CACCEPT(command,'BINARY',
     &	  ' Enter EXPORT mode (BINary, TEXt, MFIle, NAStran,'//
     &    ' H or Quit): ')

	  IF (LCMP(command,'QUIT',1)) THEN
	    RETURN

	  ELSE IF (LCMP(command,'BINary',3)) THEN
	    iExMode=1
	    mfil=mfil(1:LENGTH)//'.bin'
  2	    CONTINUE
	      INQUIRE (FILE=mfil,EXIST=exist)
	      IF (exist) THEN
	        LENGTH=ICLEN(mfil)
	        CALL CACCEPT(ANS,'YES',
     &	        'File '//mfil(1:LENGTH)//' exists.  Overwrite?')
	        IF (LCMP(ans,'NO',1)) THEN
	          CALL CACCEPT(mfil,mfil,'Enter new .bin-file name:')
	          IF ((ICLEN(mfil).EQ.1).AND.(mfil(1:1).EQ.'q')) RETURN
	          LENGTH=ICLEN(mfil)
	          mfil=mfil(1:LENGTH)//'.bin'
	          GO TO 2
	        ELSE
	          OPEN (UNIT=8,FILE=mfil,STATUS='OLD')
	          CLOSE (UNIT=8,STATUS='DELETE')
	        END IF
	      END IF
	    OPEN (8,FILE=mfil,STATUS='NEW',FORM='UNFORMATTED')

	  ELSE IF (LCMP(command,'TEXT',3)) THEN
	    iExMode=2
	    mfil=mfil(1:LENGTH)//'.txt'
  3	    CONTINUE
	      INQUIRE (FILE=mfil,EXIST=exist)
	      IF (exist) THEN
	        LENGTH=ICLEN(mfil)
	        CALL CACCEPT(ANS,'YES',
     &	        'File '//mfil(1:LENGTH)//' exists.  Overwrite?')
	        IF (LCMP(ans,'NO',1)) THEN
	          CALL CACCEPT(mfil,mfil,'Enter new .txt-file name:')
	          IF ((ICLEN(mfil).EQ.1).AND.(mfil(1:1).EQ.'q')) RETURN
	          LENGTH=ICLEN(mfil)
	          mfil=mfil(1:LENGTH)//'.txt'
	          GO TO 3
	        ELSE
	          OPEN (UNIT=8,FILE=mfil,STATUS='OLD')
	          CLOSE (UNIT=8,STATUS='DELETE')
	        END IF
	      END IF
	    OPEN (8,FILE=mfil,STATUS='NEW',FORM='FORMATTED')

	  ELSE IF (LCMP(command,'MFILE',3)) THEN
	    iExMode=3
	    mfil=mfil(1:LENGTH)//'.m'
  4	    CONTINUE
	      INQUIRE (FILE=mfil,EXIST=exist)
	      IF (exist) THEN
	        LENGTH=ICLEN(mfil)
	        CALL CACCEPT(ANS,'YES',
     &	        'File '//mfil(1:LENGTH)//' exists.  Overwrite?')
	        IF (LCMP(ans,'NO',1)) THEN
	          CALL CACCEPT(mfil,mfil,'Enter new .m-file name:')
	          IF ((ICLEN(mfil).EQ.1).AND.(mfil(1:1).EQ.'q')) RETURN
	          LENGTH=ICLEN(mfil)
	          mfil=mfil(1:LENGTH)//'.m'
	          GO TO 4
	        ELSE
	          OPEN (UNIT=8,FILE=mfil,STATUS='OLD')
	          CLOSE (UNIT=8,STATUS='DELETE')
	        END IF
	      END IF
	    OPEN (8,FILE=mfil,STATUS='NEW',FORM='FORMATTED')

	  ELSE IF (LCMP(command,'NAStran',3)) THEN
	    iExMode=4
	    mfil=mfil(1:LENGTH)//'.dat'
  5         CONTINUE
	    INQUIRE (FILE=mfil,EXIST=exist)
            IF (exist) THEN
              LENGTH=ICLEN(mfil)
              CALL CACCEPT(ANS,'YES',
     &        'File '//mfil(1:LENGTH)//' exists.  Overwrite?')
	      IF (LCMP(ans,'NO',1)) THEN
                CALL CACCEPT(mfil,mfil,'Enter new .dat-file name:')	  	    
	        IF ((ICLEN(mfil).EQ.1).AND.(mfil(1:1).EQ.'q')) RETURN
	        LENGTH=ICLEN(mfil)
                mfil=mfil(1:LENGTH)//'.dat'
	        GO TO 5
	      ELSE
                OPEN (UNIT=8,FILE=mfil,STATUS='OLD')
                CLOSE (UNIT=8,STATUS='DELETE')
              END IF
	    END IF
	    OPEN (8,FILE=mfil,STATUS='NEW',FORM='FORMATTED')

	  ELSE IF (LCMP(command,'HELP',1)) THEN
	    WRITE(*,200)

	  ELSE
	    GO TO 1
	  END IF

	  WRITE(*,*)' Exporting data to '//mfil

	  CALL CACCEPT(ANS,'CHFRay',
     &	   'Export data for CHFRay only, ALLRays, or NOChfray? ')
	  IF (LCMP(ANS,'CHFRAY',2)) THEN
	    nFirstRay=1
	    nLastRay=1
	  ELSE IF (LCMP(ANS,'NOCHFRAY',2)) THEN
	    nFirstRay=2
	    nLastRay=nRay
	  ELSE
	    nFirstRay=1
	    nLastRay=nRay
	  END IF

	  ! Currently for test, when iExMode=4 (Nastran), 
	  ! "CHFRay only" is the only option -jzlou
	  if (iExMode==4) then
	    nFirstRay=1; nLastRay=1
	  end if

	  eRay=0
	  DO iRay=nFirstRay,nLastRay
	    IF ((iObsOpt.EQ.0)
     &	    .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	    .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	      eRay=eRay+1
	    END IF
 	  END DO

C  Write data to file

 14	CONTINUE
	  CALL CACCEPT(command,'QUIT',
     &	  ' EXPORT (CMAtrix, RayIndex, RAYPos, RAYDir,'//
     &	  ' RAYL, H or Q): ')

	  IF (LCMP(command,'QUIT',1)) THEN
	    CLOSE(8)
	    RETURN

	  ELSE IF (LCMP(command,'RAYINDEX',4)) THEN
	    IF (iExMode.EQ.1) THEN
	      WRITE(*,100)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8)(RayIndex(i,iRay),i=1,2)
	        END IF
 	      END DO
	    ELSE IF (iExMode.EQ.2) THEN
	      WRITE(*,108)eRay,iCurRayElt
	      WRITE(8,108)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,109)(RayIndex(i,iRay),i=1,2)
	        END IF
 	      END DO
	    ELSE IF (iExMode.EQ.3) THEN
	      WRITE(*,100)eRay,iCurRayElt
	      WRITE(8,100)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,110)jRay,(RayIndex(i,iRay),i=1,2)
	        END IF
 	      END DO
	    END IF

	  ELSE IF (LCMP(command,'RAYPOS',4)) THEN
	    IF (iExMode.EQ.1) THEN
	      WRITE(*,100)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8)jRay,(RayPos(i,iRay),i=1,3)
	        END IF
 	      END DO
	    ELSE IF (iExMode.EQ.2) THEN
	      WRITE(*,100)eRay,iCurRayElt
	      WRITE(8,100)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,101)jRay,(RayPos(i,iRay),i=1,3)
	        END IF
 	      END DO
	    ELSE IF (iExMode.EQ.3) THEN
	      WRITE(*,100)eRay,iCurRayElt
	      WRITE(8,100)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,102)jRay,(RayPos(i,iRay),i=1,3)
	        END IF
 	      END DO
	    END IF

	  ELSE IF (LCMP(command,'RAYDIR',4)) THEN
	    IF (iExMode.EQ.1) THEN
	      WRITE(*,103)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0).OR.
     &	        ((iObsOpt.EQ.1).AND.LRayPass(iRay)).OR.
     &	        ((iObsOpt.EQ.1).AND.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8)jRay,(RayDir(i,iRay),i=1,3)
	        END IF
 	      END DO
	    ELSE IF (iExMode.EQ.2) THEN
	      WRITE(*,103)eRay,iCurRayElt
	      WRITE(8,103)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,101)jRay,(RayDir(i,iRay),i=1,3)
	        END IF
 	      END DO
	    ELSE IF (iExMode.EQ.3) THEN
	      WRITE(*,103)eRay,iCurRayElt
	      WRITE(8,103)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,104)jRay,(RayDir(i,iRay),i=1,3)
	        END IF
	      END DO
	    END IF

	  ELSE IF (LCMP(command,'RAYL',4)) THEN
	    IF (iExMode.EQ.1) THEN
	      WRITE(*,105)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8)jRay,CumRayL(iRay)
	        END IF
	      END DO
	    ELSE IF (iExMode.EQ.2) THEN
	      WRITE(*,105)eRay,iCurRayElt
	      WRITE(8,105)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,106)jRay,CumRayL(iRay)
	        END IF
	      END DO
	    ELSE IF (iExMode.EQ.3) THEN
	      WRITE(*,105)eRay,iCurRayElt
	      WRITE(8,105)eRay,iCurRayElt
	      jRay=0
	      DO iRay=nFirstRay,nLastRay
	        IF ((iObsOpt.EQ.0)
     &	        .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	        .OR.((iObsOpt.EQ.2).AND..NOT.LRayPass(iRay))) THEN
	          jRay=jRay+1
	          WRITE(8,107)jRay,CumRayL(iRay)
	        END IF
	      END DO
	    END IF

	  ELSE IF (LCMP(command,'CMATRIX',3)) THEN
	    IF (.NOT.ifBuild) THEN
	      WRITE (*,*)' C-matrix has not been BUILT.'
	      GO TO 14
	    END IF

	    CALL CACCEPT(ANS,'LOCAL',
     &	      'Export Cmatrix in GLObal or LOCal coordinates? ')
	    NCCOLS=0
	    IF (LCMP(ANS,'LOCAL',3)) THEN
	      IfEltCoords=.TRUE.
	      DO i=1,iPertElt
	        j=IABS(nECoord(i))
		nCMCoords(i)=j
	        NCCOLS=NCCOLS+j
	        IF (j.NE.0) CALL DTSP2(TEsp,TElt,6,6,i,j)
	      END DO
	      CALL COMPACT(Tcompact,Tout,nOutCord,7,7,7)
	      NCROWS=nOutCord*(nLastRay-nFirstRay+1)
	    ELSE
	      IfEltCoords=.FALSE.
	      DO i=1,iPertElt
		nCMCoords(i)=6
	        NCCOLS=NCCOLS+6
	      END DO
	      NCROWS=7*(nLastRay-nFirstRay+1)
	    END IF

	    CALL CACCEPT(ANS,'NO',
     &	      'Include input ray partials in Cmatrix? ')
	    IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	      IfInPartials=.TRUE.
	      nCMCoords(0)=7
	      NCCOLS=NCCOLS+7
	    ELSE
	      IfInPartials=.FALSE.
	      nCMCoords(0)=0
	    END IF

	    ! Write out header 
	    IF ((iExMode.EQ.2).OR.(iExMode.EQ.3)) THEN
      	      CALL MFil1(iPertElt,nCMCoords,nOutCord,
     &	                 nFirstRay,nLastRay,eRay) 
	    ELSE IF (iExMode.EQ.4) THEN
	      ! NAStran format
	      Call NasFil1(iPertElt,IfInPartials,
     &                     nOutCord,nECoord(1:iPertElt),
     &                     totRow,totCol)
	    END IF

	    irow=1; jRay=0
	    DO iRay=nFirstRay,nLastRay
	      IF ((iObsOpt.EQ.0).OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))
     &	           .OR.((iObsOpt.EQ.1).AND.LRayPass(iRay))) THEN
	        jRay=jRay+1
	        iElt=0
	        icol=1
	        IF (IfEltCoords.AND.IfInPartials) THEN
	          ! write out Cmatrix wrt src in local coord
		  CALL RMPRDB(Cout,Tcompact,Cmatrix(1,1,iRay),
     &	          	      nOutCord,7,7)
	          if (iExMode/=4) then
	            CALL MFil2(Cout,nOutCord,7,irow,icol,iRay,iElt,
     &	          	     iExMode,jRay)
	          else 
	            Call NasFil2(Cout,nOutCord,7,irow,icol)
	          end if
	          icol=8
	        ELSE IF (IfInPartials) THEN
	          ! write out Cmatrix wrt src in global coord
	          if (iExMode/=4) then
	            CALL MFil2(Cmatrix(1,1,iRay),7,7,irow,icol,iRay,iElt,
     &	          	     iExMode,jRay)
	          else
	            Call NasFil2(Cmatrix(1,1,iRay),7,7,irow,icol)
	          end if
	          icol=8
	        END IF

	        i=8
	        DO iElt=1,iPertElt
	          IF ((IfEltCoords).AND.(nECoord(iElt).NE.0)) THEN
	            ! local coord
	            CALL RMPRDC(Celt,Cmatrix(1,i,iRay),TEsp(1,1,iElt),
     &	            	   7,6,IABS(nECoord(iElt)))
	            CALL RMPRDB(Cout,Tcompact,Celt,
     &	                   nOutCord,7,IABS(nECoord(iElt)))
	            if (iExMode/=4) then
	              CALL MFil2(Cout,nOutCord,IABS(nECoord(iElt)),
     &	            	     irow,icol,iRay,iElt,iExMode,jRay)
	            else 
	              Call NasFil2(Cout,nOutCord,IABS(nECoord(iElt)),
     &                       irow,icol)
	            end if
	            icol=icol+IABS(nECoord(iElt))
	          ELSE IF (nECoord(iElt).NE.0) THEN
	            ! global coord
	            if (iExMode/=4) then
	              CALL MFil2(Cmatrix(1,i,iRay),7,6,irow,
     &	            	       icol,iRay,iElt,iExMode,jRay)
	            else 
	              Call NasFil2(Cmatrix(1,i,iRay),7,6,irow,icol)
	            end if
	            icol=icol+6
	          END IF
	          i=i+6
                END DO
	        IF (IfEltCoords) THEN
	          irow=irow+nOutCord
	        ELSE
		  irow=irow+7
	        END IF
	      END IF
	    END DO  ! iRay loop

	  ELSE IF (LCMP(command,'HELP',1)) THEN
	    WRITE(*,200)

	  ELSE
	    WRITE(*,202)
	    print*,'CEXPORT(): command = ', command
	  END IF

	GO TO 14
	END SUBROUTINE CEXPORT

C***********************************************************************

	Subroutine CPERTURB_GRP_DVR(iElt)
        use param_mod
        use src_mod
        use elt_mod
	use math_mod

#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

        IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: ANS
        REAL*8 :: th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9),
     &            DDEF(7),PV(6),RptVec(3)
        INTEGER :: i,j,k,iElt
        LOGICAL :: LLocCoord, LCMP

        CALL DZERO(DDEF,7)
	LLocCoord=.true.

        CALL CACCEPT(ANS,'GLOBAL',
     &      'Enter coordinate system for perturbation '//
     &      '(ELEMENT or GLOBAL): ')

        IF (LCMP(ANS,'Y',1).OR.LCMP(ANS,'E',1)) THEN
          ! element coord system
          k=nECoord(iElt)
          IF (k .EQ. 0) THEN
             ! k = 0, no dof
             WRITE (*,*) 'No local DOF available!'
             CALL DZERO(PV,6)
          ELSE
             IF (k .LT. 0) THEN
               ! just global coordinates, continue anyway
               k = -k
	       LLocCoord=.false.
             END IF
             CALL DACCEPT(D1,DDEF,k,
     &         'Enter perturbation vector in element coordinates:')

             ! Compute PV=TElt*D1, converting perturbations specified
             ! in local coords into values in global coords. -jzlou
             CALL DMPROD(PV,TElt(1,1,iElt),D1,6,k,1)
          END IF
        ELSE
          ! global coord system
	  LLocCoord=.false.
          CALL DACCEPT(PV(1:3),DDEF,3,
     &      'Enter rotational perturbation vector (x,y,z):')
          CALL DACCEPT(PV(4:6),DDEF,3,
     &      'Enter translational perturbation vector (x,y,z):')
          ! Perturbations entered in global coords are used as they are. -jzlou
        END IF

!       print*,'* CPERTURB_GRP_DVR: iElt =',iElt
!       print*,'* CPERTURB_GRP_DVR: EltGrp(0,iElt) =',EltGrp(0,iElt)
!       print*,'* CPERTURB_GRP_DVR: EltGrp(1:EltGrp(0,iElt),iElt) =',
!    &              EltGrp(1:EltGrp(0,iElt),iElt)
!       STOP

	RptVec(1:3)=RptElt(1:3,iElt)
        Do i=1,EltGrp(0,iElt)
!         print*,'****B4 CPERTURB_GRP: EltGrp(i,iElt) =',
!    &           EltGrp(i,iElt)  
          Call CPERTURB_GRP(RptVec,EltGrp(i,iElt),th,del,Q,dQ,D1,D2,PV,
     &                      LLocCoord)
        End Do

        ! Update the group rotation point with group translations
        RptElt(1:3,iElt) = RptElt(1:3,iElt) + del ! jzlou, 01/15/2016
        End Subroutine CPERTURB_GRP_DVR

C***********************************************************************

        Subroutine CPERTURB_GRP(rptVec,jElt,th,del,Q,dQ,D1,D2,PV,lloc)
        use param_mod
        use src_mod
        use elt_mod
	use math_mod
	use loh_I1, only : SegXgrid_FLG

#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif
	use traceutil_mod, only : CRIncidPosNom

        IMPLICIT NONE

        INTEGER :: jElt,i,j,k
        LOGICAL :: lloc  ! global or local coord system
        REAL*8 :: rptVec(3),th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9),
     &            DDEF(7),PV(6)

        Do i=1,3
          j=i+3
          th(i)=PV(i); del(i)=PV(j)
        End Do

	If (jElt==0) Then
	  ! source is perturbed
          CALL Qform(Q,dQ,th)

          CALL DMPROD(D1,Q,ChfRayDir,3,3,1)
          CALL DEQUATE(ChfRayDir,D1,3)
          CALL DUNITIZE(ChfRayDir)

          CALL DMPROD(D1,Q,xGrid,3,3,1)
          CALL DEQUATE(xGrid,D1,3)
          CALL DUNITIZE(xGrid)
          CALL DMPROD(D1,Q,yGrid,3,3,1)
          CALL DEQUATE(yGrid,D1,3)
          CALL DUNITIZE(yGrid)

          print*,'*** Grp rotation of src: Rpt =',rptVec(1:3)  

          !CALL DSUB(D1,ChfRayPos,StopPos,3)
        ! When source is part of a group perturbation, the source is
        ! rotated about the group rotation point 'rptVec' instead of
        ! about the STOP, which makes sense. -jzlou 08/30/2014
          CALL DSUB(D1,ChfRayPos,rptVec,3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,ChfRayPos,D2,3)
          CALL DADD(ChfRayPos,D1,del,3)

	  If (SegXgrid_FLG) Then
            CALL DMPROD(D1,Q,SegXGrid,3,3,1)
            CALL DEQUATE(SegXGrid,D1,3)
            CALL DUNITIZE(SegXGrid)
          End If
	  return
	End If

        ! Perturb a non-source element below
        CALL Qform(Q,dQ,th)

        ! update psiElt(*,jElt)
        CALL DMPROD(D1,Q,psiElt(1,jElt),3,3,1)
        CALL DEQUATE(psiElt(1,jElt),D1,3)

        ! update VptElt(*,jElt)
        CALL DSUB(D1,VptElt(1,jElt),rptVec,3)
        CALL DMPROD(D2,dQ,D1,3,3,1)
        CALL DADD(D1,VptElt(1,jElt),D2,3)
        CALL DADD(VptElt(1,jElt),D1,del,3)

        ! update RptElt(*,jElt)
        CALL DSUB(D1,RptElt(1,jElt),rptVec,3)
        CALL DMPROD(D2,dQ,D1,3,3,1)
        CALL DADD(D1,RptElt(1,jElt),D2,3)
        CALL DADD(RptElt(1,jElt),D1,del,3)

        !CALL DADD(D1,rptVec,del,3)
        !CALL DEQUATE(RptElt(1,jElt),D1,3)

	! update CRIncidPosNom, -jzlou
        CALL DSUB(D1,CRIncidPosNom(1,jElt),rptVec,3)
        CALL DMPROD(D2,dQ,D1,3,3,1)
        CALL DADD(D1,CRIncidPosNom(1,jElt),D2,3)
        CALL DADD(CRIncidPosNom(1,jElt),D1,del,3)

	! update metrology points on surface
	If (nMetPos(jElt)>0) Then
	  Do j=1,nMetPos(jElt)
	    CALL DSUB(D1,SrfMetPos(1,j,iEltToMetSrf(jElt)),
     &                rptVec,3)
	    CALL DMPROD(D2,dQ,D1,3,3,1)
            CALL DADD(D1,SrfMetPos(1,j,iEltToMetSrf(jElt)),D2,3)
            CALL DADD(SrfMetPos(1,j,iEltToMetSrf(jElt)),D1,del,3)
	  End Do
	End If

C fixed by Gopal Vasudevan:
        !IF (nObs(jElt).NE.0) THEN
        ! This IF test is removed by jzlou, since an aperture needs
        ! to be perturbed with the element.
        ! update xObs(*,jElt)
          CALL DMPROD(D1,Q,xObs(1,jElt),3,3,1)
          CALL DEQUATE(xObs(1,jElt),D1,3)
        !END IF

        IF ((EltID(jElt).EQ.4).OR.(EltID(jElt).EQ.5).OR.
     &      (EltID(jElt).EQ.13)) THEN
          ! update h1HOE(*,jElt)
          CALL DSUB(D1,h1HOE(1,jElt),rptVec,3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,h1HOE(1,jElt),D2,3)
          CALL DADD(h1HOE(1,jElt),D1,del,3)

          ! update h2HOE(*,jElt)
          CALL DSUB(D1,h2HOE(1,jElt),rptVec,3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,h2HOE(1,jElt),D2,3)
          CALL DADD(h2HOE(1,jElt),D1,del,3)
        END IF

        IF ((SrfType(jElt).GE.4).AND.(SrfType(jElt).LE.9)) THEN
          ! update pMon(*,jElt)
          CALL DSUB(D1,pMon(1,jElt),rptVec,3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,pMon(1,jElt),D2,3)
          CALL DADD(pMon(1,jElt),D1,del,3)

          ! update xMon(*,jElt)
          CALL DMPROD(D1,Q,xMon(1,jElt),3,3,1)
          CALL DEQUATE(xMon(1,jElt),D1,3)

          ! update yMon(*,jElt)
          CALL DMPROD(D1,Q,yMon(1,jElt),3,3,1)
          CALL DEQUATE(yMon(1,jElt),D1,3)

          ! update zMon(*,jElt)
          CALL DMPROD(D1,Q,zMon(1,jElt),3,3,1)
          CALL DEQUATE(zMon(1,jElt),D1,3)
        END IF

        If (lloc .and. LUpdateTElt_FLG(jElt)) Then
          ! TElt(1:6,1:6,jElt) is updated here
          Do j=1,3
            CALL DMPROD(D1,Q,TElt(1:3,j,jElt),3,3,1)
            CALL DEQUATE(TElt(1:3,j,jElt),D1,3)
            CALL DMPROD(D1,Q,TElt(4:6,3+j,jElt),3,3,1)
            CALL DEQUATE(TElt(4:6,3+j,jElt),D1,3)
          End Do
        End If
        End Subroutine CPERTURB_GRP
	
C***********************************************************************
C****** End file funcsub.F *********************************************
C***********************************************************************

