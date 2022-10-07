C***********************************************************************
C  File macos.F
C  MACOS main program
C
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2012, California Institute of Technology.  |
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
C       MACOS Version 3.50 	               Last updated: 07/2019
C       (Based on MACOS Versions 3.34, 3.31, 3.2 and 2.86)
C       MACOS Version 3.34 	               Last updated: 11/2011
C       (Based on MACOS Versions 3.31, 3.2 and 2.86)
C
C       Contact:     John Z. Lou
C       Address:     Jet Propulsion Laboratory 
C                    MS 306-451 
C                    4800 Oak Grove Dr., Pasadena CA 91109
C       Email:       John.Z.Lou@jpl.nasa.gov
C        Cell:       (818)653-8024
C      Office:       (818)354-4870
C***********************************************************************

C***********************************************************************
C       ** History **
C       MACOS Version 2.86		           08/2001
C
C       Contact:     David C. Redding
C       Address:     Jet Propulsion Laboratory
C                    MS 306-451
C                    4800 Oak Grove Dr., Pasadena CA 91109
C       Email:       dcr@huey.jpl.nasa.gov
C       Phone:       (818)354-3696
C***********************************************************************

#include "realtype.h"

	Program MACOS
	use param_mod
	use elt_mod
	use src_mod
	use cfiles_mod
	use macos_mod

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
	use design_cons_optim_mod
        use stop_set_mod
	use loh_I1, only : BaseUnits_FLG

	Implicit None

	Interface
          Subroutine macos_init_all(psize)
          implicit none
          integer, optional :: psize
          End Subroutine macos_init_all
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

	Call set_macos_problem_size(macos_model_size) 

	mVarDOF_nls=mVarDOF    ! mVarDOF_nls, mOptZern_nls and mOptAsph_nls 
	mOptZern_nls=mOptZern  ! are declared in design_optim_mod
	mOptAsph_nls=mOptAsph
	mVarDOF_np=mVarDOF     ! mVarDOF_np and mOptZern_np
        mOptZern_np=mOptZern   ! are declared in design_cons_optim_mod 

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

