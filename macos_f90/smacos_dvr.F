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

	use macos_vars_mod  ! variables specific for macos.F
	use dopt_mod
	use zern_wf_mod
	use design_optim_mod
	use design_cons_optim_mod
        use stop_set_mod
	use loh_I1, only : BaseUnits_FLG

	Implicit none

	CHARACTER(len=MacosCharLen) :: command,CARG(9),RxName
        REAL*8 :: DARG(9),tmpDARG(9)
        INTEGER :: IARG(9),modelSize
        LOGICAL :: LARG
        SREAL :: RARG(9),RmsWfe
	SREAL, dimension(:,:), allocatable, save :: PixArray
        REAL*8, dimension(:,:), allocatable, save ::OPDMat,RaySpot


C  Execution starts
	modelSize=512; nGridPts=255
	Call macos_init_all(modelSize)

	allocate(PixArray(mPix,mPix),
     &           OPDMat(mpts,mpts),RaySpot(mRay,2));

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
