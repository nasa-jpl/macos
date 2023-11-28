C***********************************************************************
C       Begin file dopt_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2009, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	MODULE dopt_mod
          USE param_mod
          USE elt_mod

          IMPLICIT NONE
          SAVE

	! Design optimization related variables

	Integer, parameter :: WFE_TARGET       = 1,
     &                        WFE_ZMODE_TARGET = 2,
     &                        BEAM_TARGET      = 3,
     &                        SPOT_TARGET      = 4,
     &                        OPL_TARGET       = 5

	Integer, parameter :: mVarDOF    =  8,
     &                        mOptZern   = 45,   ! ToDo -- update to new values
     &                        mOptAsph   =  9,
     &                        max_fov    = 12,
     &                        max_wl     =  6

        INTEGER, PARAMETER :: Linear = 1,
     &                        NonLin = 2

        CHARACTER(len=*), PARAMETER :: DOPT_TARGET(5)=
     &                    (/'RMS Wavefront Error',
     &                      'WFE Zernike Modes',
     &                      'Beam Property',
     &                      'Spot Size',
     &                      'OPL'/)

        CHARACTER(len=*), PARAMETER :: DOF_NameList(8)=
     &                    (/'TIP', 'TILT', 'CLOCK',
     &                       'DX', 'DY', 'PIST','ROC', 'CONIC'/)

	! total DOFs = 6 rb + ROC + Conic + mOptZern
	Character(len=6), dimension(:,:), allocatable :: DOF_VarElt
	Character(len=7) :: valbuf
        Character(len=32) :: tgtWF_fn,tgtWF_fmt
        Character(len=MacosValLen) :: PinvFile
	Logical, dimension(:), allocatable :: isVarElt,isVarElt_Rx
	Integer, dimension(:,:),allocatable :: varEltDOF,
     &           OptEltZernTerm,OptEltAsphTerm
	Integer :: OptTarget,nVarElt,nOptWavelen,nVarElt_Rx,
     &             nOptFov,nOptChfRayDir,nOptChfRayPos,off
	Integer, dimension(:), allocatable :: varElts,varDOFs,varZerns,
     &            nDOF_VarElt,nOptEltZern,zernElt,WFZernMode,
     &            nOptEltAsph,asphElt,varAsphs
	Integer :: varEltID(2),zernRange(2),nID,slen,ifov,jfov,iwl,
     &             nWFZern,m_err_dopt
	Real*8 :: fov_ChfRayDir(3),fov_ChfRayPos(3),
     &		  ChfRayDir_Rx(3),ChfRayPos_Rx(3),OptSrcRpt(3),
     &            opt_fov(3,2,max_fov), ! max. FOV to optimize
     &            fov_weight(max_fov),
     &            opt_wavelen(max_wl)   ! max wavelens to optimize

        Integer :: iVarElt,nitrs_dopt,rtn_flag,zern_off,asph_off
        Real*8, dimension(:,:), allocatable :: psiElt_old,vptElt_old,
     &                          ZernCoef_old,AsphCoef_old,varEltCons
	Real*8, dimension(:,:,:), allocatable :: EltCons
        Real*8, dimension(:), allocatable :: KrElt_old, KcElt_old
        Real*8 :: old_wfe(max_fov,max_wl),new_wfe(max_fov,max_wl),
     &            dopt_tol,eltRotations(3,3),consArr(2*mVarDOF),
     &            OptSpotSize,SvdSvCut
	Integer :: OptTgtElt,OptAlg,aparam_off
	Logical :: LOptIfFEX,isSrcOpt,LOptCons,tgtWF_flg,ifOptPrint,
     &             LOptNomSens,OptSrcRpt_FLG,LCmdCalibElts,
     &             save_pinv, ! whether to save linear optimizer pseduo inverse
     &             use_saved_pinv

	Type(OptBeamStr) :: beamOpt


	Contains

        ! --------------------------------------------------------------
        SUBROUTINE dopt_init

          IF (ALLOCATED(varElts)) THEN
            DEALLOCATE(varElts, varDOFs, isVarElt,varEltDOF,
     &                 psiElt_old, vptElt_old, KrElt_old, KcElt_old,
     &                 ZernCoef_old,AsphCoef_old,
     &                 nDOF_VarElt,DOF_VarElt,
     &                 nOptEltZern,zernElt,nOptEltAsph,asphElt,
     &                 varZerns,varAsphs,isVarElt_Rx,
     &                 OptEltZernTerm,OptEltAsphTerm,WFZernMode,
     &                 EltCons,varEltCons,stat=m_err_dopt)
            IF (m_err_dopt /= 0)
     &        CALL macos_memory_failure
     &                    ('dopt_mod: deallocate failed!')
	  END IF
          ALLOCATE(varElts(mElt),varDOFs(mVarDOF*mElt),
     &             isVarElt(0:mElt),varEltDOF(mVarDOF,0:mElt),
     &             psiElt_old(3,mElt),vptElt_old(3,mElt),
     &             KrElt_old(mElt), KcElt_old(mElt),
     &             ZernCoef_old(mOptZern,mElt),
     &             AsphCoef_old(mOptAsph,mElt),
     &             nDOF_VarElt(mElt),DOF_VarElt(mVarDOF,mElt),
     &             nOptEltZern(mElt),zernElt(mElt),
     &             nOptEltAsph(mElt),asphElt(mElt),
     &             varZerns(mOptZern*mElt),
     &             varAsphs(mOptAsph*mElt),
     &             isVarElt_Rx(0:mElt),
     &             OptEltZernTerm(0:mOptZern,mElt),
     &             OptEltAsphTerm(0:mOptAsph,mElt),
     &             EltCons(2,mVarDOF,mElt),varEltCons(2,mVarDOF*mElt),
     &	           WFZernMode(mOptZern),stat=m_err_dopt)
          IF (m_err_dopt /= 0)
     &      CALL macos_memory_failure
     &                  ('dopt_mod: allocate failed!')

          ! Initialisation
          CALL dopt_init_vars()

	END SUBROUTINE dopt_init

        ! --------------------------------------------------------------
        SUBROUTINE dopt_init_vars()
          USE constants, ONLY : EPS
          IMPLICIT NONE

          varElts(:)          = 0
          varDOFs(:)          = 0
          isVarElt(:)         = .FALSE.
          varEltDOF(:,:)      =   0
          psiElt_old(:,:)     = 0d0
          vptElt_old(:,:)     = 0d0
          KrElt_old(:)        = 0d0
          KcElt_old(:)        = 0d0
          ZernCoef_old(:,:)   = 0d0
          AsphCoef_old(:,:)   = 0d0
          nDOF_VarElt(:)      =   0
          DOF_VarElt(:,:)     =  ''
          nOptEltZern(:)      =   0
          zernElt(:)          =   0
          nOptEltAsph(:)      =   0
          asphElt(:)          =   0
          varZerns(:)         =   0
          varAsphs(:)         =   0
          isVarElt_Rx(:)      = .FALSE.
          OptEltZernTerm(:,:) =   0
          OptEltAsphTerm(:,:) =   0
          EltCons(:,:,:)      = 0d0
          varEltCons(:,:)     = 0d0
          WFZernMode(:)       =   0

          ! Design optimization related initializations
          OptTarget           = WFE_TARGET

          nitrs_dopt          = 10
          dopt_tol            = 1d-12
          nOptChfRayDir       = 0        ! in msmacosio.inc when reading ChfRayDir: nOptChfRayDir=nOptChfRayDir+1
          nOptChfRayPos       = 0        ! also: in macos_cmd_loop.inc: values are set to 0
          ifOptPrint          = .TRUE.
          SvdSvCut            = EPS      ! defined in lsq_mod
          save_pinv           = .FALSE.  ! whether to save linear optimizer pseduo inverse
          use_saved_pinv      = .FALSE.
          PinvFile            = 'dady0_saved.txt'

          valbuf              = ''
          tgtWF_fn            = ''
          tgtWF_fmt           = ''

          nVarElt             = 0
          nOptWavelen         = 0
          nVarElt_Rx          = 0
          nOptFov             = 0
          off                 = 0

          varEltID(:)         = 0
          zernRange(:)        = 0
          nID                 = 0
          slen                = 0
          ifov                = 0
          jfov                = 0
          iwl                 = 0
          nWFZern             = 0
          m_err_dopt          = 0

          fov_ChfRayDir(:)    = 0
          fov_ChfRayPos(:)    = 0
       	  ChfRayDir_Rx(:)     = 0
          ChfRayPos_Rx(:)     = 0
          OptSrcRpt(:)        = 0
          opt_fov(:,:,:)      = 0 ! max. FOV to optimize
          fov_weight(:)       = 0
          opt_wavelen(:)      = 0 ! max wavelens to optimize

          iVarElt             = 0
          nitrs_dopt          = 0
          rtn_flag            = 0
          zern_off            = 0
          asph_off            = 0

          old_wfe(:,:)        = 0d0
          new_wfe(:,:)        = 0d0
          dopt_tol            = 0d0
          eltRotations(:,:)   = 0d0
          consArr(:)          = 0d0
          OptSpotSize         = 0d0
          SvdSvCut            = 0d0

          OptTgtElt           = 0
          OptAlg              = 0
          aparam_off          = 0

          LOptIfFEX           = .FALSE.
          isSrcOpt            = .FALSE.
          LOptCons            = .FALSE.
          tgtWF_flg           = .FALSE.
          LOptNomSens         = .FALSE.
          OptSrcRpt_FLG       = .FALSE.
          LCmdCalibElts       = .FALSE.

        END SUBROUTINE dopt_init_vars

        ! --------------------------------------------------------------
	! Set one of design optimization parameters in SMACOS
        SUBROUTINE SMACOS_SetOptParam(pname,carr,iarr,darr)

          IMPLICIT NONE
          Character(len=*) :: pname,carr(2)
          Integer, dimension(:) :: iarr
          Real*8, dimension(:) :: darr
	  Logical :: LCMP

	  If (LCMP(pname,'OptTarget',9)) Then

	    if (LCMP(carr(1),'WFE',3)) then
	      OptTarget=WFE_TARGET
	    else if (LCMP(carr(1),'ZWF',3)) then
              OptTarget=WFE_ZMODE_TARGET
	    else if (LCMP(carr(1),'Beam',4)) then
	      OptTarget=BEAM_TARGET
	    end if

	  Else If (LCMP(pname,'NFov',4)) Then
	    nOptFov=iarr(1)
	    fov_weight(1:nOptFov)=darr(1:nOptFov)
	  Else If (LCMP(pname,'NItr',4)) Then
            nitrs_dopt=iarr(1)
	  Else If (LCMP(pname,'ErrTol',6)) Then
	    dopt_tol=darr(1)
	  Else If (LCMP(pname,'TgtWFFile',9)) Then
	     tgtWF_fmt=carr(1)
	     tgtWF_fn=carr(2)
	  End If
        END SUBROUTINE SMACOS_SetOptParam

        END MODULE dopt_mod
