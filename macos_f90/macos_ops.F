C***********************************************************************
C       Begin file macos_ops.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module macos_ops_mod
	use math_mod

	Implicit none
	Save

	Integer :: nSpotRays

	Contains

	Subroutine MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                       OPDMap,SpotMap,WFErms,PixMap,
     &                       zernIDArr,zernCoefArr,
     &                       asphIDArr,asphCoefArr)
	use param_mod
        use src_mod
        use elt_mod
	use macos_mod

#ifdef CMACOS
	use macos_vars_mod
#else
  	use smacos_vars_mod
#endif
	use sourcsub_mod
	use tracesub_mod

	CHARACTER(len=MacosCharLen) :: cmd,CARG(9)
        REAL*8 :: DARG(9)
        LOGICAL :: LARG
        INTEGER :: IARG(9)
        SREAL :: RARG(9)
        REAL*8, dimension(:,:),intent(out) :: OPDMap,SpotMap
	SREAL, dimension(:,:),intent(out) :: PixMap
	REAL*8 :: WFErms,WFEPV_loc
	INTEGER, dimension(:), optional :: zernIDArr,asphIDArr
        REAL*8, dimension(:), optional :: zernCoefArr,asphCoefArr
	LOGICAL :: LCMP

C  FORMAT statements
#include "macos_formats.inc"

	PixMap=0d0  ! to avoid compiler warning; it's not used later
	iElt=IARG(1)

	If (LCMP(cmd,'STOP',4)) Then
	  StopPos(1:3)=DARG(1:3)
	  ifStopSet=.TRUE.

	Else If (LCMP(cmd,'FEX',3)) Then
          !print*,'In OPT FEX ...'
	  !stop
	  ifGrid=.FALSE.
          ifPropagate=.FALSE.
          ifTrace=.FALSE.
          ifBuild=.FALSE.
          ifOPDCalc=.FALSE.

	  CALL FEX(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &      ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &      iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &      nBadRays,StartIndRef,StartExtinc,
     &      Cmatrix,DrawRayVec,DrawEltVec,
     &      iElt,zp,psip,CrossPt,ifCentroid,
     &      CRSpot,CentroidSpot,dCRCentroid)

	    eElt(iElt)=0d0
            fElt(iElt)=ABS(zp)
            KcElt(iElt)=0d0
            KrElt(iElt)=-fElt(iElt)
            zElt(iElt)=zp
            CALL DEQUATE(psiElt(1,iElt),psip,3)
            CALL DEQUATE(VptElt(1,iElt),CrossPt,3)
            CALL DEQUATE(RptElt(1,iElt),CrossPt,3)

	Else If (LCMP(cmd,'PERT',4)) Then
	  IF ((iElt.EQ.0).AND.(.NOT.ifStopSet)) THEN
!               WRITE(*,*) '** Error (MACOS_OPS): Must set stop'//
!    &      'before perturbing source. Use STOP command.'
!            stop
          END IF
	  !
	  Call CPERTURB_2(iElt,CARG,DARG)
	  ifTrace=.FALSE.
          ifBuild=.FALSE.
          ifPropagate=.FALSE.
          ifPerturb=.FALSE.
          ifGrid=.FALSE.
          iCurRayElt=0
	  !
	  ! In double-pass system, also need to perturb its
	  ! return-pass counterpart if exists.
	  Do i=1, nElt
	    if (LnkElt(i)==iElt) then
	      Call CPERTURB_2(i,CARG,DARG)
              ifTrace=.FALSE.
              ifBuild=.FALSE.
              ifPropagate=.FALSE.
              ifPerturb=.FALSE.
              ifGrid=.FALSE.
              iCurRayElt=0
	    end if
	  End Do

	Else If (LCMP(cmd,'GPER',4)) Then
          Call CPERTURB_GRP_DVR_2(iElt,CARG,DARG)
          ifTrace=.FALSE.
          ifBuild=.FALSE.
          ifPropagate=.FALSE.
          ifPerturb=.FALSE.
          ifGrid=.FALSE.
          iCurRayElt=0

	Else If (LCMP(cmd,'SPOT',4)) Then
          IF (ifLoad) THEN
            ! use 'BEAM' option in 'SPOT' command
            ifGrid=.FALSE.
            ifDraw=.FALSE.
            ifPropagate=.FALSE.
            ifTrace=.FALSE.
            ifBuild=.FALSE.
            ifOPDCalc=.FALSE.
            iLocal=iElt
            !spcOption=1  ! so that spot diagram centered around element
            spcOption=2  ! so that spot diagram centered around chief ray

            CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &        ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &        ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,
     &        dxSource,nBadRays,StartIndRef,StartExtinc,
     &        Cmatrix,DrawRayVec,DrawEltVec,
     &        iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)

            IF (ifLocalCoord) THEN
              !WRITE(*,311)xLocal,yLocal,zLocal
              CALL DZERO(Tspot(1,1),49)
              DO iCoord=1,3
                Tspot(1,iCoord)=xLocal(iCoord)
                Tspot(2,iCoord)=yLocal(iCoord)
                Tspot(3,3+iCoord)=xLocal(iCoord)
                Tspot(4,3+iCoord)=yLocal(iCoord)
              END DO
              Tspot(5,7)=1
              ! Now Tspot can be viewed as a 5x7 matrix, as
              ! xL(1) xL(2) xL(3) 0 0 0 0 0
              ! yL(1) yL(2) yL(3) 0 0 0 0 0
              ! 0     0     0     xL(1) xL(2) xL(3) 0
              ! 0     0     0     yL(1) yL(2) yL(3) 0
              ! 0     0     0     0     0     0     1
              CALL COMPACT(Tcompact,Tspot,5,7,7,7)
              CALL DCOMPACT(DTcompact,Tspot,5,7,7,7)
            ELSE
              WRITE(*,*)'** macos_ops: SPOT command error!'
              STOP
            END IF
            RefSpot(1)=DDOTC(xLocal,VptElt(1,iEndElt))
            RefSpot(2)=DDOTC(yLocal,VptElt(1,iEndElt))

            ! Trace rays
            iEndElt=iElt; iStartElt=0
            IF ((iStartElt.EQ.0).OR.(iStartElt.LT.iEndElt)) THEN
              ifRayPrint=.FALSE.
              IF (iStartElt.EQ.0) THEN
                IF (ABS(zSource).GE.1d10) THEN
                  CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &            ifPropagate)
                ELSE
                  CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &            ifPropagate)
                END IF
                ifGrid=.TRUE.
              END IF
              WRITE(*,*)' Tracing ',nRay,' rays...'
              iStartRay=1
              iEndRay=nRay
              CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &                    ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &                    iStartRay,iEndRay,RayMap,SegMap,
     &                    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &                    StartIndRef,StartExtinc,DrawRayVec,
     &                    DrawEltVec,nDrawElt(1),xDraw,yDraw)
              iCurRayElt=iEndElt
              ifTrace=.TRUE.
              ifBuild=.FALSE.
              ifPropagate=.FALSE.
            END IF

	    ! Compute spot diagram data from ray trace
            iSpot=0
            nGridPts=npts+1
            mpt2=mpts*mpts
            xcent=0d0; ycent=0d0
            DO iRay=1,nRay
              IF (LRayOK(iRay).AND.((iRay.EQ.1).OR.(iObsOpt.EQ.0).OR.
     &          (LRayPass(iRay).AND.(iObsOpt.EQ.1)).OR.
     &          (.NOT.LRayPass(iRay).AND.(iObsOpt.EQ.2))))THEN
                DO iaxis=1,3
                  j=iaxis+3
                  xdum(iaxis)=RayDir(iaxis,iRay)
                  xdum(j)=RayPos(iaxis,iRay)
                END DO
                xdum(7)=0e0
                CALL DMPRDB(xnom,DTcompact,xdum,5,7,1)
                IF (iRay.EQ.1) THEN
                  CntrSpot(1)=xnom(3)
                  CntrSpot(2)=xnom(4)
                END IF
	        IF ((iRay.NE.1).OR.
     &            ((iRay.EQ.1).AND.(nGridpts.EQ.1))) THEN
                  iSpot=iSpot+1
                  IF (spcOption.EQ.1) THEN
                    SpotMap(iSpot,1)=xnom(3)-RefSpot(1)
                    SpotMap(iSpot,2)=xnom(4)-RefSpot(2)
                  ELSE
                    SpotMap(iSpot,1)=xnom(3)-CntrSpot(1)
                    SpotMap(iSpot,2)=xnom(4)-CntrSpot(2)
                  END IF
                  xcent=xcent+xnom(3)
                  ycent=ycent+xnom(4)
                END IF
              END IF
            END DO
            nSpotRays=iSpot
            NBEAMS=1
            IF (iSpot.GT.0) THEN
              xcent=xcent/DBLE(iSpot)
              ycent=ycent/DBLE(iSpot)
              WRITE(*,408)(CntrSpot(1)-RefSpot(1)),
     &            (CntrSpot(2)-RefSpot(2)),xcent,ycent
            ELSE
              WRITE(*,*)' All rays were obscured or lost'
            END IF
          ELSE
            WRITE(*,*) '** Error (MACOS_OPS): load Rx first!'
            STOP
          END IF

	Else If (LCMP(cmd,'OPD',3)) Then
	  IF (ifLoad) THEN
	    isDOP=.true.
            iEndElt=iElt
            CALL OPD(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &       ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &       iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &       nBadRays,StartIndRef,StartExtinc,OPDMap,WFErms,WFEPV_loc,
     &       DAvgL,L1,Cmatrix,DrawRayVec,DrawEltVec)
	    Avgl = DAvgl
	    ifOPDCalc=.TRUE.
            CALL DXCALC(iElt,iEm1,dxMin(iElt),dxMax(iElt),
     &               dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
	    !WRITE(*,442) WFErms,DAvgl
	    isDOP=.false.

	  Else
	    WRITE(*,*) '** Error (MACOS_OPS): load Rx first!'
	    stop
	  End If

	Else If (LCMP(cmd,'CRT',3)) Then
	  IF (ifLoad) THEN
            isDOP=.true.
            iEndElt=iElt
	    CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &        ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &        ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &        iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &        nBadRays,StartIndRef,StartExtinc,
     &        Cmatrix,DrawRayVec,DrawEltVec)
            isDOP=.false.
          Else
            WRITE(*,*) '** Error (MACOS_OPS): load Rx first!'
            stop
          End If

	Else If (LCMP(cmd,'GBS',3)) Then
	  isDOP=.true.
	  iEndElt=iElt
          Call GetBeamSizeCmd(bmRad,npts,ifPol,ifGrid,ifPropagate,
     &            ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &            ifIntSrfSet,ifRayPrint,iObsOpt,iEndElt,iCurWFElt,
     &            iCurRayElt,Dicr,Djcr,dxSource,nBadRays,StartIndRef,
     &            StartExtinc,OPDMap,WFE,DAvgL,L1,Cmatrix,DrawRayVec,
     &            DrawEltVec)
          isDOP=.false.
	  BeamSizeAtElt=2d0*bmRad

	Else If (LCMP(cmd,'RefRay',6)) Then
	  isDOP=.true.
          iEndElt=iElt
	  !print*,' ** macos_ops: refRay =',refRay(1:5)
!	  Call TraceRefRays(refRay,nRefRay,iEndElt,
!    &                      RefRayDir(1:3,nRefRay))
	  Call TraceRefRays(iEndElt)
	  isDOP=.false.

	Else If (LCMP(cmd,'ROC_PERT',8)) Then
	  KrElt(iElt)=KrElt(iElt)+DARG(1)
	  Do i=1, nElt
            if (LnkElt(i)==iElt)
     &	      KrElt(i)=KrElt(i)+DARG(1)
	  End Do

	Else If (LCMP(cmd,'CONIC_PERT',10)) Then
	  KcElt(iElt)=KcElt(iElt)+DARG(1)
	  Do i=1, nElt
            if (LnkElt(i)==iElt)
     &        KcElt(i)=KcElt(i)+DARG(1)
          End Do

	Else If (LCMP(cmd,'ZERN_PERT',9)) Then
	  ! IARG(2) = no. of Zern terms
	  Do i=1,IARG(2)
	    j=zernIDArr(i)  ! Zern term id
	    ZernCoef(j,iElt)=ZernCoef(j,iElt)+zernCoefArr(i)
	  End Do

	  Do i=1,nElt  ! for double-pass
            if (LnkElt(i)==iElt) then
	      Do k=1,IARG(2)
	        j=zernIDArr(k)  ! Zern term id
	        ZernCoef(j,i)=ZernCoef(j,i)+zernCoefArr(k)
	      End Do
	    end if
          End Do

	Else If (LCMP(cmd,'ASPH_PERT',9)) Then
          ! IARG(2) = no. of Asph coeffs
	  Do i=1,IARG(2)
            j=asphIDArr(i)  ! Asph term id
            AsphCoef(j,iElt)=AsphCoef(j,iElt)+asphCoefArr(i)
          End Do

	  Do i=1,nElt  ! for double-pass
            if (LnkElt(i)==iElt) then
              Do k=1,IARG(2)
                j=asphIDArr(k)  ! Asph term id
                AsphCoef(j,i)=AsphCoef(j,i)+asphCoefArr(k)
              End Do
            end if
          End Do

	Else If (LCMP(cmd,'RESET',5)) Then
	    ifTrace=.FALSE.
            ifBuild=.FALSE.
            ifPropagate=.FALSE.
            ifCmt=.FALSE.
            ifPerturb=.FALSE.
            DO iElt=1,mElt
              ifRegrid(iElt)=.FALSE.
            END DO
            ifPol=.FALSE.
            ifVecDif3=.FALSE.
            ifOld=.FALSE.
            ifRayPrint=.FALSE.
            ifGrid=.FALSE.
            ifPixLoc=.FALSE.
#ifdef MACOS_CMD
            ifAdd=.FALSE.
            ifPixArray=.FALSE.
#endif
            ifDMBuild=.FALSE.
            ifCentroid=.TRUE.
            ifFilterData=.FALSE.
            ifLNegOK=.TRUE.
            ifUDSinit=.FALSE.
            ifStopSet=.FALSE.
            i=mElt
            CALL LZERO(ifIntSrfSet,i)
            iCurRayElt=0
            iCurWFElt=0
            iObsOpt=1
            StretchType=1
            spcOption=1
            BeamType=1
            sizPix=0d0
            EltPix(1)=0d0
            EltPix(2)=0d0
            WinCenPix(1)=0d0
            WinCenPix(2)=0d0
	End If
	End Subroutine MACOS_OPS

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Adapted from CPERTURB routine in funcsub.F

	Subroutine CPERTURB_2(iElt,CARG,DARG)
        use param_mod
        use src_mod
        use elt_mod
	use dopt_mod, only : OptSrcRpt_FLG,OptSrcRpt

#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

        IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: CARG(9),ANS
        REAL*8 :: DARG(9)
        INTEGER i,j,k, iElt
        LOGICAL LCMP

        REAL*8 th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9),
     &         DDEF(7),PV(6),del2(3)

 400    FORMAT(' Computing new perturbed system parameters')

        CALL DZERO(DDEF,7)

C  1. Compute perturbed entrance pupil data

        IF (iElt.EQ.0) THEN  ! source is perturbed
c         CALL DZERO(D1,6)
c         CALL DACCEPT(th,D1,3,
c    &       'Enter rotational perturbation vector (x,y,z):')
c         CALL DACCEPT(del,D1,3,
c    &       'Enter translational perturbation vector (x,y,z):')

	  th=DARG(1:3)
	  del=DARG(4:6)

          CALL Qform(Q,dQ,th)
          CALL DMPROD(D1,Q,ChfRayDir,3,3,1)
          CALL DEQUATE(ChfRayDir,D1,3)
          CALL DUNITIZE(ChfRayDir)

          CALL DMPROD(D1,Q,xGrid,3,3,1)
          CALL DEQUATE(xGrid,D1,3)
          CALL DMPROD(D1,Q,yGrid,3,3,1)
          CALL DEQUATE(yGrid,D1,3)

	  If (.not.OptSrcRpt_FLG) Then
            CALL DSUB(D1,ChfRayPos,StopPos,3)
	  Else
	    CALL DSUB(D1,ChfRayPos,OptSrcRpt,3)
	  End If

          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,ChfRayPos,D2,3)
          !CALL DADD(ChfRayPos,D1,del,3)  % old
	  ! Assume del is specified in (xG,yG,ChfRayDir) frame
	  del2=del(1)*xGrid+del(2)*yGrid+del(3)*ChfRayDir
	  CALL DADD(ChfRayPos,D1,del2,3)

	  !print*,'macos_ops: StopPos =',StopPos
	  !print*,'macos_ops: Purturbed ChfRayPos =',ChfRayPos

C  2. Compute perturbed element data

        ELSE
	  ! not source element
c	  CALL CACCEPT(ANS,'GLOBAL',
c    &      'Enter coordinate system for perturbation '//
c    &      '(ELEMENT or GLOBAL): ')
c         IF (LCMP(ANS,'Y',1).OR.LCMP(ANS,'E',1)) THEN
	  IF ((LCMP(CARG(1),'ELE',3)).OR.(LCMP(CARG(2),'ELE',3)).OR.
     &        (LCMP(CARG(3),'ELE',3))) THEN
            ! element coord system
            k= nECoord(iElt)
            IF (k .EQ. 0) THEN
               ! k = 0, no dof
               WRITE (*,*) 'No local DOF available!'
               CALL DZERO(PV,6)
            ELSE
               IF (k .LT. 0) THEN
                 ! just global coordinates, continue anyway
                 k = -k
               END IF
c              CALL DACCEPT(D1,DDEF,k,
c    &           'Enter perturbation vector in element coordinates:')
	       D1(1:6)=DARG(1:6)

C              print*,'***CPERTURB DDEF = ', DDEF,', D1(1:6) = ',
C    &             D1(1),D1(2),D1(3),D1(4),D1(5),D1(6)

               ! Compute PV=TElt*D1, converting perturbation specified
               ! local coords into global coords. -jzlou
               CALL DMPROD(PV,TElt(1,1,iElt),D1,6,k,1)
            END IF
          ELSE
            ! global corrd system
c           CALL DACCEPT(PV(1:3),DDEF,3,
c    &        'Enter rotational perturbation vector (x,y,z):')
c           CALL DACCEPT(PV(4:6),DDEF,3,
c    &        'Enter translational perturbation vector (x,y,z):')
	    PV(1:6)=DARG(1:6)
            ! Perturbations entered in global coords are used as they are. -jzlou

C           print*,'***CPERTURB PV(1:6) = ',
C    &         PV(1),PV(2),PV(3),PV(4),PV(5),PV(6)
          END IF

c         WRITE(*,400)

          DO 8 i=1,3
            j=i+3
            th(i)=PV(i)
            del(i)=PV(j)
  8       CONTINUE
          CALL Qform(Q,dQ,th)
          CALL DMPROD(D1,Q,psiElt(1,iElt),3,3,1)
          CALL DEQUATE(psiElt(1,iElt),D1,3)

          CALL DSUB(D1,VptElt(1,iElt),RptElt(1,iElt),3)
          CALL DMPROD(D2,dQ,D1,3,3,1)
          CALL DADD(D1,VptElt(1,iElt),D2,3)
          CALL DADD(VptElt(1,iElt),D1,del,3)

          CALL DADD(D1,RptElt(1,iElt),del,3)
          CALL DEQUATE(RptElt(1,iElt),D1,3)
C fixed by Gopal Vasudevan:
          !IF (nObs(iElt).NE.0) THEN
          ! This IF test is removed by jzlou, since an aperture needs
          ! to be perturbed with the element.
             CALL DMPROD(D1,Q,xObs(1,iElt),3,3,1)
             CALL DEQUATE(xObs(1,iElt),D1,3)
          !END IF

          IF ((EltID(iElt).EQ.4).OR.(EltID(iElt).EQ.5).OR.
     &        (EltID(iElt).EQ.13)) THEN
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
        END IF

        RETURN
        END Subroutine CPERTURB_2

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Subroutine CPERTURB_GRP_DVR_2(iElt,CARG,DARG)
        use param_mod
        use src_mod
        use elt_mod

        Implicit none

	CHARACTER(len=MacosCharLen) :: CARG(9)
        REAL*8 :: DARG(9)
        INTEGER :: i,j,k,iElt
        LOGICAL :: LLocCoord, LCMP
        REAL*8 :: th(3),del(3),Q(3,3),dQ(3,3),D1(9),D2(9),
     &            DDEF(7),PV(6),RptVec(3)


        CALL DZERO(DDEF,7)

        if (nECoord(iElt)<0) then
          LLocCoord=.false.  ! use global coords
          !print*,' **CPERTURB_GRP_DVR_2: use global coords'
        else
          LLocCoord=.true.   ! use local coords
          !print*,' **CPERTURB_GRP_DVR_2: use local coords'
        end if

	! Not for source element
        IF (LCMP(CARG(1),'ELE',3)) THEN
          ! element coord system
          k= nECoord(iElt)

          IF (k .EQ. 0) THEN
            ! k = 0, no dof
            WRITE (*,*) 'No local DOF available!'
            CALL DZERO(PV,6)
          ELSE
            IF (k .LT. 0) THEN
              ! just global coordinates, continue anyway
              k = -k
            END IF
            D1(1:6)=DARG(1:6)

            ! Compute PV=TElt*D1, converting perturbation specified
            ! local coords into global coords.
            ! TElt was initialized to identity matrix in the case of global
            ! coord. -jzlou
            CALL DMPROD(PV,TElt(1,1,iElt),D1,6,k,1)
          END IF
        ELSE
          ! global corrd system
          PV(1:6)=DARG(1:6)
          ! Perturbations entered in global coords are used as they are. -jzlou
        END IF

        if (.false.) then
         print*,'macos_ops(): EltGrp(0,iElt) =',EltGrp(0,iElt)
         print*,'macos_ops(): EltGrp(1:,iElt) =',
     &         EltGrp(1:EltGrp(0,iElt),iElt)
         print*,'macos_ops(): RptElt(1:3,iElt) =',RptElt(1:3,iElt)
         stop
        end if

        RptVec(1:3)=RptElt(1:3,iElt)
        Do i=1,EltGrp(0,iElt)
          Call CPERTURB_GRP(RptVec,EltGrp(i,iElt),th,del,Q,dQ,D1,D2,PV,
     &                      LLocCoord)
        End Do
        End Subroutine CPERTURB_GRP_DVR_2

        End Module macos_ops_mod
