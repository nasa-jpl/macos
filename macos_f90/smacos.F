C***********************************************************************
C	File smacos.F
C
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  ALL RIGHTS RESERVED. United State Government Sponsorship Is   |
C    |  Acknowledged. Any commercial use must be negotiated with the  |
C    !  Office of Technology Transfer at the California Institute of  |
C    |  Technology.                                                   |
C    |                                                                |
C    | This software may be subject to U.S. export control laws. By   |
C    ! accepting this software, the user agrees to comply with all    |
C    | applicable U.S. export laws and regulations. User has the      |
C    | responsibility to obtain export licenses, or other export      |
C    | authority as may be required before exporting such information |
C    | to foreign countries or providing access to foreign persons."  |
C    +----------------------------------------------------------------+
C
C***********************************************************************

C***********************************************************************
C       Subroutine Modeling and Analysis for Controlled Optical Systems  
C       SMACOS Version 3.31            	     05/2007
C       (Based On SMACOS Versions 2.86 and 3.2)
C
C       Contact:     John Z. Lou
C       Address:     Jet Propulsion Laboratory
C                    MS 306-451
C                    4800 Oak Grove Dr., Pasadena CA 91109
C       Email:       John.Z.Lou@jpl.nasa.gov
C       Phone:       (818)354-4870
C***********************************************************************

C***********************************************************************
C       ** History **
C       SMACOS Version 2.86                  08/16/2001
C
C       Contact:     David Redding
C       Address:     Jet Propulsion Laboratory 
C                    MS 306-451 
C                    4800 Oak Grove Dr., Pasadena CA 91109
C       Email:       dcr@huey.jpl.nasa.gov
C       Phone:       (818)354-3696
C***********************************************************************

#include "realtype.h"

	Module smacos_mod
	Implicit none

	Integer, parameter :: max_cmd_len=MacosCharLen

	Contains

	SUBROUTINE SMACOS(commandARG,CARG,DARG,IARG,LARG,RARG,
     &	                  OPDMat,RaySpot,RMSWFE,PixArray)

	use param_mod
	use elt_mod
	use src_mod
	use cfiles_mod
	use macos_mod

	use math_mod
	use smacosio_mod
	use sourcsub_mod
	use smacosutil_mod
	use tracesub_mod
	use traceutil_mod
        use propsub_mod
	use lohpars_mod

	use dopt_mod
	use smacos_vars_mod
 	use design_optim_mod
	use design_cons_optim_mod 
	use stop_set_mod
	use loh_I1, only : BaseUnits_FLG

	IMPLICIT NONE

	! Subroutine arguments
	Character(len=max_cmd_len) :: commandARG 
        Character(len=max_cmd_len), dimension(:), intent(in) :: CARG(9)
	Logical, intent(in) :: LARG
	Integer, dimension(:), intent(in) :: IARG
 	SREAL, dimension(:), intent(in) :: RARG
	Real*8, dimension(:), intent(in) :: DARG
	
	!Integer, intent(in) :: IARG(9)
	!SREAL, intent(in) :: RARG(9)  
	!Real*8, intent(in) :: DARG(9)

        SREAL, intent(out) :: RMSWFE
	SREAL, intent(out), dimension(:,:) :: PixArray
        Real*8, intent(out), dimension(:,:) :: OPDMat,RaySpot

	INTEGER :: ICLEN,STRLEN
	LOGICAL :: LCMP


c***loh	Declares for journal and batch files
	INTEGER row_sz,col_sz
cc***loh
C
C	DATA ifInit/.FALSE./
C	DATA ifEcho/.FALSE./
C

	! SMACOS array storage management
	If (first_enter .or. macos_realloc) then
	  if (allocated(L1)) then
	    deallocate(L1, DrawEltVec,R1,R2,PertVec,
     &                 DrawRayVec,D2,DV1,DV2, 
     &                 CD1,CD2,DWF,
     &		       stat=m_err_smacos)
	    if (m_err_smacos /= 0) call macos_memory_failure
     &               ('SMACOS: deallocate failed!')
	  end if
	  !
	  allocate(L1(md2), DrawEltVec(mDrawElt,mDrawRay),
     &	           R1(mdttl,mdttl), R2(mdttl,mdttl),
     &             PertVec(6,mElt), 
     &             DrawRayVec(2,mDrawElt,mDrawRay),
     &	  	   D2(mdttl,mdttl),
     &             DV1(md2), DV2(md2),
     &		   CD1(mdttl,mdttl), CD2(mdttl,mdttl),
     &		   DWF(mdttl,mdttl), stat=m_err_smacos)
	  if (m_err_smacos /= 0) call macos_memory_failure
     &               ('SMACOS: allocate failed!') 
	  L1 = .false.; DrawEltVec = 0
	  R1 = 0.0; R2 = 0.0; PertVec = 0.0
	  DrawRayVec = 0.0; D2 = dzero_
	  DV1 = dzero_; DV2 = dzero_
	  CD1 = (0d0,0d0); CD2 = (0d0,0d0)
	  DWF = czero_ 
	  ifInit = .false.  ! for glass-table etc.
	  ifEcho = .false.
	  macos_realloc = .false.; first_enter = .false.

	  mVarDOF_nls=mVarDOF    ! mVarDOF_nls and mOptZern_nls
	  mOptZern_nls=mOptZern  ! are defined in design_optim_mod

	  Call dopt_init
	End If


C  FORMAT statements
#include "macos_formats.inc"

C       Glass table
C       INCLUDE 'glass.inc'


C  3. Main program

C	*** Initialize stuff first time around ***

        IF (.NOT.ifInit) THEN

#include "smacos_glass.inc"

          bat_task=1
          jou_task=1
c sws
          pstack = 1
          read_len(pstack) = 0
c start of file unit numbers (could be n+10)
          ustack= 50
c sws
        END IF

C  ***Load command stack***

        command=commandARG

        ! -jzlou
        ! Added for converting from ZMAX input to MACOS input.
        ! If the filename (CARG(1)) has the extension 'ZMX' or 'zmx',
        ! the input file will be converted from ZMX format to a 
        ! MACOS prescription file with the extension '.in'.
        ! The returned CARC(1) will be contain the same filename 
        ! but without an extension. 
        !
        !if (LCMP(command,'LOAD',3).OR.LCMP(command,'OLD',2)) then
        !  Call MACOSZ2MDVR(CARG(1))
        !end if

        CALL LoadStack(command,CARG,DARG,IARG,LARG,RARG)

        !print*,'==========> smacos.F: command =', command


C  ***Begin SMACOS top-level command loop***

#undef  MACOS_CMD
#define SMACOS_CMD
#define DESIGN_OPTIM

#include "macos_cmd_loop.inc"


C       ***Here to quit***

c  99   CONTINUE
c       STOP

        RETURN
        END SUBROUTINE SMACOS

	End module smacos_mod

