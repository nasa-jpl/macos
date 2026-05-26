C
C smacos_dvr.F
C Purpose: for testing SMACOS code (e.g. when crashing in GMI)
C John Z. Lou, JPL
C Last updated: 05-14-2009
C

#include "realtype.h"

	Program smacos_dvr
	use param_mod
	use elt_mod
	use src_mod
	use cfiles_mod
	use macos_mod
        use surfsub

	use macosio_mod
        use sourcsub_mod
        use tracesub_mod
	use traceutil_mod
        use propsub_mod
        use lohpars_mod
	use pgplot_mod
	use math_mod

	use macos_vars_mod  ! variables specific for macos.F (incl. nPgPanel)
	use dopt_mod
	use zern_wf_mod
	use design_optim_mod
#ifdef USE_NPSOL
	use design_cons_optim_mod
#endif
        use stop_set_mod
	use loh_I1, only : BaseUnits_FLG
	use smacos_mod

	Implicit none

	CHARACTER(len=MacosCharLen) :: CARG(9),RxName
        REAL*8 :: DARG(9),tmpDARG(9)
        INTEGER :: IARG(9),modelSize
        LOGICAL :: LARG
        SREAL :: RARG(9)

C  Execution starts
	modelSize=512; nGridPts=255
	Call macos_init_all(modelSize)

C  Open graphics device once for OPD/SPOT plot commands below.
C  Interactive macos.F does this via GRAINI calls scattered through
C  macos_cmd_loop.inc; smacos_dvr must open the device itself.
C  nPgPanel selects 1x1/1x2/1x3/2x2 layout (see pgplotsub.F:GRAINI);
C  default 1 = single panel.
        nPgPanel = 1
        CALL GRAINI

	RxName='SegDemo3'		! 'most3frntdbl'

	command='OLD'
        CARG(1)=RxName
	CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RmsWfe,PixArray)
	 
!	command='SPOT'
!	!command='OPD'
!        CARG='BEA'
!        IARG(1)=50
!        CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
!     &              OPDMat,RaySpot,RmsWfe,PixArray)
	command='OPD'
        IARG(1)=9
        CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RmsWfe,PixArray)
	command='SPOT'
        CARG='TOUT'
        IARG(1)=10
        CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RmsWfe,PixArray)
	Stop
	End Program smacos_dvr
