C***********************************************************************
C       Begin file stop_set.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************
C
C  Use MACOS nonlinear optimizer to enforce system STOP set at an element
C  for systems with a point source.
C

#include "realtype.h"

	Module  stop_set_mod

	Contains

	SUBROUTINE SetStop(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &    ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &    ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &    Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec,iLocal,xLocal,yLocal,
     &    xTarg,sqErr,ifBeamCentered)

        use param_mod
        use src_mod
        use cfiles_mod
        use elt_mod
	use tracesub_mod
        use dopt_mod, only : beamOpt,mVarDOF,old_wfe,new_wfe,rtn_flag,
     &   		     ifOptPrint,OptAlg,NonLin
        use design_optim_mod

        IMPLICIT NONE

	LOGICAL :: ifPol,ifGrid,ifPropagate,ifRayPrint,
     &    ifLNegOK,ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &    ifDMBuild,ifDraw,ifBeamCentered

        INTEGER :: npts,iCurWFElt,iCurRayElt,
     &          nBadRays,iEndElt,iObsOpt,
     &          DrawEltVec(mDrawElt,mDrawRay),iLocal

        SREAL :: Cmatrix(7,mCm,bRay),
     &           DrawRayVec(2,mDrawElt,mDrawRay)

        REAL*8 :: Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &   	xTarg(3),th(3),
     &   	oldChfRayDir(3),ChfRayDirSv(3),ChfRayPosSv(3),CrsDist,
     &   	xLocal(3),yLocal(3),sqErr,errTol,zp

	Logical :: LOptIfFEX_loc,tgtWF_flg_loc
        Integer :: OptTarget_loc,nVarElt_loc,nOptFov_loc,
     &             OptTgtElt_loc,nOptWavelen_loc,OptAlg_sv,
     &             nitrs_dopt_loc,nWFZern_loc,npts_sv,opt_npts_sv
        Integer :: varElts_loc(1),varDOFs_loc(mVarDOF),
     &             zernElt_loc(1),asphElt_loc(1),
     &             WFZernMode_loc(1),varZerns_loc(1),varAsphs_loc(1)
        Character(len=32) :: tgtWF_fn_loc, tgtWF_fmt_loc
        REAL*8 :: opt_wavelen_loc(1),ddx_loc(3),ddy_loc(3),dopt_tol_loc,
     &            opt_fov_loc(3,2,12),fov_wt_loc(1)
        Type(OptBeamStr) :: beamOpt_sv

	! Compute STOP position in 3D global frame
        CALL DSMPROD(ddx_loc,xLocal,xTarg(1),3)
        CALL DSMPROD(ddy_loc,yLocal,xTarg(2),3)
        !print*,'** xTarg(1:2) =',xTarg(1:2)
        beamOpt%nomBeamPos= VptElt(1:3,iLocal)+ddx_loc+ddy_loc

#if 0
    ! Test only: compute initial STOP position error
        iEndElt=iLocal  ! = id of STOP element
        ! Trace chief ray from source to STOP element (=iEndElt)
        CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)

        del=RayPos(1:3,1)-beamOpt%nomBeamPos(1:3)

        print*,'beamOpt%nomBeamPos(1:3) =',beamOpt%nomBeamPos(1:3)
        print*,'RayPos(1:3,1) =',RayPos(1:3,1)
        sqErr=DDOTC(del,del)
        WRITE(*,*)'Initial ChfRayDir =',ChfRayDir(1:3)
        WRITE(*,*)'Initial ChfRayPos =',ChfRayPos(1:3)
        WRITE(*,*)'** Initial ChfRay STOP error =',
     &            sqrt(sqErr)*1d06,' nm'
#endif 

        ! Use MACOS nonlinear optimizer to enforce element STOP for pt source
	OptTarget_loc=BEAM_TARGET
	nVarElt_loc=1
	varElts_loc(1)=0  ! varElt is source
	varDOFs_loc(1:mVarDOF)=0; varDOFs_loc(1:2)=1  ! DOF = tip/tilt of source  
	zernElt_loc(1)=0
	asphElt_loc(1)=0
	nOptFov_loc=1
	opt_fov_loc(1:3,1,1)=ChfRayDir(1:3)
	opt_fov_loc(1:3,2,1)=ChfRayPos(1:3)
	fov_wt_loc(1)=1d0
	nOptWavelen_loc=1
	opt_wavelen_loc(1)=Wavelen
	tgtWF_flg_loc=.false.
	nWFZern_loc=0
	OptTgtElt_loc=iLocal  ! STOP element 
	LOptIfFEX_loc=.false.
	nitrs_dopt_loc=12  ! max iterations
	errTol=1d-22
	dopt_tol_loc=0.5*sqrt(errTol)
	oldChfRayDir(1:3)=ChfRayDir(1:3)

	beamOpt_sv=beamOpt  ! save global beamOpt

	beamOpt%ifOptBeamDir=.false.
        beamOpt%ifOptBeamPos=.true.
	beamOpt%beamPosElt=iLocal

	!StopPos=ChfRayPos(1:3)  ! rotation point for chief ray, temporary
 	StopPos=ChfRayPos(1:3)+zSource*ChfRayDir(1:3) ! correct way, keep real source point
						      ! unchanged.
	ifStopSet=.TRUE.  ! tentatively

	OptAlg_sv=OptAlg; OptAlg=NonLin
	npts_sv=npts; npts=2
        opt_npts_sv=opt_npts; opt_npts=npts
	
	ifOptPrint=.false.

        ! Unconstrained system optimization
        Call nls_optim_dvr(.true.,OptTarget_loc,
     &                     nVarElt_loc,varElts_loc,varDOFs_loc,
     &                     nWFZern_loc,WFZernMode_loc,
     &                     zernElt_loc(1:nVarElt_loc),varZerns_loc,
     &                     asphElt_loc(1:nVarElt_loc),varAsphs_loc,
     &                     opt_fov_loc,nOptFov_loc,fov_wt_loc,
     &                     opt_wavelen_loc,nOptWavelen_loc,
     &                     dopt_tol_loc,nitrs_dopt_loc,
     &                     old_wfe,new_wfe,rtn_flag,
     &                     tgtWF_flg_loc,tgtWF_fn_loc,tgtWF_fmt_loc,
     &                     OptTgtElt_loc,LOptIfFEX_loc,beamOpt)

#if 0
	! Check convergence
        iEndElt=iLocal  ! = id of STOP element
        ! Trace chief ray from source to STOP element (=iEndElt)
        CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)

        del=RayPos(1:3,1)-beamOpt%nomBeamPos(1:3)
	sqErr=DDOTC(del,del)
	WRITE(*,*)'Final ChfRayDir =',ChfRayDir(1:3)
	WRITE(*,*)'Final ChfRayPos =',ChfRayPos(1:3)
	WRITE(*,*)'** Final ChfRay STOP error =',
     &            sqrt(sqErr)*1d06,' nm'	
#endif
	
	WRITE(*,*)' ' 
	If (sqErr.LE.errTol) Then
          WRITE(*,*)' STOP set successfully at element',iLocal 
	  ifBeamCentered=.TRUE.
	Else
	  WRITE(*,*)' WARNING: Stop calculation failed to converge!'
          CALL DEQUATE(ChfRayDir,oldChfRayDir,3)
          ifBeamCentered=.FALSE.
	  ifStopSet=.FALSE.
	End If


	! ** Now compute entrance pupil position **

	ChfRayDirSv(1:3)=ChfRayDir(1:3)
	ChfRayPosSv(1:3)=ChfRayPos(1:3)

C     Compute new differential chief ray, displaced in -yGrid:
        CALL DSMPROD(th,-yGrid,1d-5,3)
	CALL DADD(ChfRayPos,ChfRayPosSv,th,3)

#if 0
    ! Test only: compute initial differential STOP position error
        iEndElt=iLocal  ! = id of STOP element
        ! Trace differential chief ray from source to STOP element (=iEndElt)
        CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)

        del=RayPos(1:3,1)-beamOpt%nomBeamPos(1:3)

        print*,'beamOpt%nomBeamPos(1:3) =',beamOpt%nomBeamPos(1:3)
        print*,'RayPos(1:3,1) =',RayPos(1:3,1)
        sqErr=DDOTC(del,del)
        WRITE(*,*)'Initial ChfRayDir =',ChfRayDir(1:3)
        WRITE(*,*)'Initial ChfRayPos =',ChfRayPos(1:3)
        WRITE(*,*)'** Initial differential ChfRay STOP error =',
     &            sqrt(sqErr)*1d06,' nm'
#endif

	! Enforce STOP for differential chief ray 

	opt_fov_loc(1:3,1,1)=ChfRayDir(1:3)
        opt_fov_loc(1:3,2,1)=ChfRayPos(1:3)

	xGrid_sv=xGrid; yGrid_sv=yGrid
        Call nls_optim_dvr(.true.,OptTarget_loc,
     &                     nVarElt_loc,varElts_loc,varDOFs_loc,
     &                     nWFZern_loc,WFZernMode_loc,
     &                     zernElt_loc(1:nVarElt_loc),varZerns_loc,
     &                     asphElt_loc(1:nVarElt_loc),varAsphs_loc,
     &                     opt_fov_loc,nOptFov_loc,fov_wt_loc,
     &                     opt_wavelen_loc,nOptWavelen_loc,
     &                     dopt_tol_loc,nitrs_dopt_loc,
     &                     old_wfe,new_wfe,rtn_flag,
     &                     tgtWF_flg_loc,tgtWF_fn_loc,tgtWF_fmt_loc,
     &                     OptTgtElt_loc,LOptIfFEX_loc,beamOpt)
	xGrid=xGrid_sv; yGrid=yGrid_sv

#if 0		
	! Check STOP convergence for differential chief ray
        iEndElt=iLocal  ! = id of STOP element
        ! Trace chief ray from source to STOP element (=iEndElt)
        CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &    ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &    ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &    iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,
     &    Cmatrix,DrawRayVec,DrawEltVec)

        del=RayPos(1:3,1)-beamOpt%nomBeamPos(1:3)
        sqErr=DDOTC(del,del)
        WRITE(*,*)'Final ChfRayDir =',ChfRayDir(1:3)
        WRITE(*,*)'Final ChfRayPos =',ChfRayPos(1:3)
        WRITE(*,*)'** Final differential ChfRay STOP error =',
     &            sqrt(sqErr)*1d06,' nm'
#endif

	! Find crossing point of the two chief rays:
	CALL FindCrossPt(ChfRayPosSv,ChfRayDirSv,
     &                ChfRayPos,ChfRayDir,zp,CrsDist,StopPos)
	WRITE(*,*)' Computed StopPos =',StopPos(1:3)
	WRITE(*,*)' '

	ChfRayPos=ChfRayPosSv; ChfRayDir=ChfRayDirSv

	! Restore parameters to pre-StopSet values 
        OptAlg=OptAlg_sv
        npts=npts_sv; opt_npts=opt_npts_sv
	nGridPts=npts+1
        beamOpt=beamOpt_sv
	ifOptPrint=.true.
	END SUBROUTINE SetStop

	End Module stop_set_mod
 
