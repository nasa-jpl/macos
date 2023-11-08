C***********************************************************************
C       Begin file dopt_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2009, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module dopt_mod
	use param_mod
	use nn_type, only : eps
	use elt_mod

	Implicit none
	Save

	! Design optimization related variables

	Integer, parameter :: WFE_TARGET=1,WFE_ZMODE_TARGET=2,
     &                        BEAM_TARGET=3,SPOT_TARGET=4,OPL_TARGET=5

	Integer, parameter :: mVarDOF=8, mOptZern=45,mOptAsph=9,
     &           max_fov=12, max_wl=6 
	! total DOFs = 6 rb + ROC + Conic + mOptZern
	Character(len=6), dimension(8) :: DOF_NameList
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

        Character(len=40), dimension(5) :: DOPT_TARGET
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
	Integer, parameter :: Linear=1,NonLin=2
	

	Contains

	Subroutine dopt_init
	If (allocated(varElts)) Then
	  deallocate(varElts,varDOFs,isVarElt,varEltDOF,
     &               psiElt_old,vptElt_old,KrElt_old,KcElt_old,
     &               ZernCoef_old,AsphCoef_old,
     &               nDOF_VarElt,DOF_VarElt,
     &               nOptEltZern,zernElt,nOptEltAsph,asphElt,
     &               varZerns,varAsphs,isVarElt_Rx,
     &               OptEltZernTerm,OptEltAsphTerm,WFZernMode,
     &               EltCons,varEltCons,stat=m_err_dopt) 
	  if (m_err_dopt /= 0) then
              call macos_memory_failure
     &                  ('dopt_mod: deallocate failed!')
          end if
	End If
        allocate(varElts(mElt),varDOFs(mVarDOF*mElt),
     &           isVarElt(0:mElt),varEltDOF(mVarDOF,0:mElt),
     &           psiElt_old(3,mElt),vptElt_old(3,mElt),
     &           KrElt_old(mElt), KcElt_old(mElt),
     &           ZernCoef_old(mOptZern,mElt),
     &           AsphCoef_old(mOptAsph,mElt),
     &           nDOF_VarElt(mElt),DOF_VarElt(mVarDOF,mElt),
     &           nOptEltZern(mElt),zernElt(mElt),
     &           nOptEltAsph(mElt),asphElt(mElt),
     &           varZerns(mOptZern*mElt),
     &           varAsphs(mOptAsph*mElt),
     &           isVarElt_Rx(0:mElt),
     &           OptEltZernTerm(0:mOptZern,mElt),
     &           OptEltAsphTerm(0:mOptAsph,mElt), 
     &           EltCons(2,mVarDOF,mElt),varEltCons(2,mVarDOF*mElt),
     &	         WFZernMode(mOptZern),stat=m_err_dopt)
        if (m_err_dopt /= 0) then
           call macos_memory_failure
     &                ('dopt_mod: allocate failed!')
        end if
        isVarElt=.false.; varEltDOF=0
        varElts=0; varDOFs=0
        psiElt_old=0d0; vptElt_old=0d0

	! Design optimization related initializations
	OptTarget=WFE_TARGET  ! default
        DOPT_TARGET(1)='RMS Wavefront Error'
        DOPT_TARGET(2)='WFE Zernike Modes'
        DOPT_TARGET(3)='Beam Property'
        DOPT_TARGET(4)='Spot Size'
	DOPT_TARGET(5)='OPL'
        DOF_NameList(1)='TIP';DOF_NameList(2)='TILT'
        DOF_NameList(3)='CLOCK'; DOF_NameList(4)='DX'
        DOF_NameList(5)='DY'; DOF_NameList(6)='PIST'
        DOF_NameList(7)='ROC'; DOF_NameList(8)='CONIC'
        !nitrs_dopt=10; dopt_tol=1d-08  ! defaults
        nitrs_dopt=10; dopt_tol=1d-12  ! defaults
        nOptChfRayDir=1; nOptChfRayPos=1
	ifOptPrint=.TRUE.
	SvdSvCut=eps  ! defined in lsq_mod
        save_pinv=.false.; use_saved_pinv=.false. ! default
        PinvFile='dady0_saved.txt'  ! default
	End Subroutine dopt_init

C***********************************************************************

	! Set one of design optimization parameters in SMACOS
        Subroutine SMACOS_SetOptParam(pname,carr,iarr,darr)
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
        End Subroutine SMACOS_SetOptParam

	End Module dopt_mod
