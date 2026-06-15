C***********************************************************************
C  File macos.F
C  MACOS main program
C
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2026, California Institute of Technology.  |
C    |  ALL RIGHTS RESERVED. United State Government Sponsorship Is   |
C    |  Acknowledged. Any commercial use must be negotiated with the  |
C    !  Office of Technology Transfer at the California Institute of  |
C    |  Technology.                                                   |
C    |  							      |
C    | This software may be subject to U.S. export control laws. By   |
C    ! accepting this software, the user agrees to comply with all    |
C    | applicable U.S. export laws and regulations. User has the      |
C    | responsibility to obtain export licenses, or other export      |
C    | authority as may be required before exporting such information |
C    | to foreign countries or providing access to foreign persons.   |
C    +----------------------------------------------------------------+
C
C***********************************************************************

C***********************************************************************
C       Modeling and Analysis for Controlled Optical Systems (MACOS)
C	MACOS Version 4.00		       Work in progress
C***********************************************************************
C
C       MACOS Version 4.00 	               Last updated: 4/2026
C       Contacts:    Scott Basinger, Scott.A.Basinger@jpl.nasa.gov
C                    John Lou, John.Z.Lou@jpl.nasa.gov
C                    Luis Marchen, Luis.F.Marchen@jpl.nasa.gov
C                    Carl Nissly, Carl.R.Nissly@jpl.nasa.gov
C                    David C. Redding, david.c.redding@gmail.com
C                    Norbert Sigrist, Norbert.Sigrist@jpl.nasa.gov
C
C       Based on MACOS Versions 3.31, 3.2 and 2.86, and COMP, the
C       Controlled Optics Modeling Program
C       Authors:     David C. Redding, Laura Needels, Sam Sirlin,
C                    Scott Basinger, John Lou, Luis Marchen, Carl 
C                    Nissly, and Norbert Sigrist
C
C***********************************************************************

#include "realtype.h"

	Program MACOS
	use param_mod
	use elt_mod
	use src_mod
	use cfiles_mod
	use macos_mod
	use sysprop_mod, only: DeriveSysProps, PrintSysProps
	use plot_mode_mod
        use surfsub

	use macosio_mod
        use sourcsub_mod
        use tracesub_mod
	use traceutil_mod
        use propsub_mod
        use lohpars_mod
	use pgplot_mod
	use math_mod

	use macos_vars_mod  ! variables specific for macos.F
	use dopt_mod
	use zern_wf_mod
	use design_optim_mod
#ifdef USE_NPSOL
	use design_cons_optim_mod
#endif
	use design_slsqp_optim_mod
        use stop_set_mod
	use loh_I1, only : BaseUnits_FLG
        use validate_prescription_mod, only: ValidatePrescription

	Implicit None

	Interface
          Subroutine macos_init_all(psize)
          implicit none
          integer, optional :: psize
          End Subroutine macos_init_all
          Subroutine install_sigint_handler() bind(C,name='install_sigint_handler')
          use iso_c_binding
          End Subroutine install_sigint_handler
        End Interface

C  FORMAT statements
#include "macos_formats.inc"

	! External function declearations
	Logical :: LCMP,inCmdHist
	INTEGER :: ICLEN,StrLen
        !INTEGER :: IDFLOOR,STRLEN
	!Real*8 :: DDOTC,GAUSS,POISSON
	!SREAL :: POIDEV

	! Threading related
        !Integer :: macos_threads,new_macos_threads,iam

C  Execution starts

        ! Scene pthread defined in C
	!Call startpthread

	Call install_sigint_handler()
	Call set_macos_problem_size(macos_model_size)

	mVarDOF_nls=mVarDOF    ! mVarDOF_nls, mOptZern_nls and mOptAsph_nls
	mOptZern_nls=mOptZern  ! are declared in design_optim_mod
	mOptAsph_nls=mOptAsph
#ifdef USE_NPSOL
	mVarDOF_np=mVarDOF     ! mVarDOF_np and mOptZern_np
        mOptZern_np=mOptZern   ! are declared in design_cons_optim_mod
#endif
	mVarDOF_sls=mVarDOF    ! design_slsqp_optim_mod copies
	mOptZern_sls=mOptZern

        ! Initialize MACOS
        Call macos_init_all(macos_model_size)
	macos_init=.true.

	! Local dynamic arrays
        allocate(L1(md2), DrawEltVec(mDrawElt,mDrawRay),
     &           R1(mdttl,mdttl), R2(mdttl,mdttl),
     &           PertVec(6,mElt),
     &           DrawRayVec(2,mDrawElt,mDrawRay),
     &           D2(mdttl,mdttl),DV1(md2), DV2(md2),
     &	         CRWalk(3,mElt), CRWalkProjLoc(3,mElt),
     &           CD1(mdttl,mdttl), CD2(mdttl,mdttl),
     &           DWF(mdttl,mdttl), stat=m_err_macos)
	If (m_err_macos /= 0)
     &	  call macos_memory_failure('MACOS: allocate failed!')
        L1 = .false.; DrawEltVec = 0
        R1 = 0.0; R2 = 0.0; PertVec = 0.0
        DrawRayVec = 0.0; D2 = dzero_
        DV1 = dzero_; DV2 = dzero_
        CD1 = (0d0,0d0); CD2 = (0d0,0d0)
        DWF = czero_

	Call dopt_init

        WRITE(*,303) mElt,mEm1,mRay,bRay,mWF,mIntSrf,mDP,mPix,mPix,
     &		     mdttl,mdttl


C  Load glass table
#include "macos_glass.inc"

C  Variable initializations
#include "macos_init.inc"


C  *** Begin MACOS top-level command loop ***

#undef  SMACOS_CMD
#define MACOS_CMD
#define DESIGN_OPTIM

#include "macos_cmd_loop.inc"

#ifdef VIS3D_DEMO
        Call killdemo
#endif
C	STOP
        END Program MACOS

