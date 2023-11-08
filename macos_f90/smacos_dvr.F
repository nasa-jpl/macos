C
C smacos_dvr.F
C Purpose: for testing SMACOS code (e.g. when crashing in GMI)
C John Z. Lou, JPL
C Last updated: 05-14-2009
C

#include "realtype.h"

	Program smacos_dvr
	use param_mod
        use src_mod
        use elt_mod
        use macos_mod
        use smacos_mod

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

	RxName='most3frntdbl'

	command='OLD'
        CARG(1)=RxName
	CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RmsWfe,PixArray)
	 
	command='SPOT'
	!command='OPD'
        CARG='BEA'
        IARG(1)=50
        CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RmsWfe,PixArray)
	Stop
	End Program smacos_dvr
