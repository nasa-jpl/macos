C
C A Test Program for MACOS
C

#include "realtype.h"
 
	Program dopt_main
	use smacos_mod
	use design_optim_smacos_mod

	Implicit none

C  These are the SMACOS call-line variables:	
	CHARACTER(len=MacosCharLen) :: command,CARG(9)
        REAL*8 DARG(9)
        INTEGER IARG(9)
        LOGICAL LARG
        SREAL :: RARG(9),RMSWFE

C Dynamic local arrays
	SREAL, dimension(:,:), allocatable, save :: PixArray
        REAL*8, dimension(:,:), allocatable, save :: OPDMat, RaySpot

	Integer :: nVarElts, elts_arr(1), dof_arr(7),
     &             total_itrs, return_flag
	Integer, parameter :: model_size=128

	call macos_init_all(model_size)

	allocate(PixArray(mPix,mPix),OPDMat(mpts, mpts), 
     &           RaySpot(mRay,2))

	command='OLD'
        CARG(1)='/home/jzlou/Pub/j18sd_s'
        Call SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RMSWFE,PixArray)

	command='STOP'
        CARG(1)='OBJECT'
        CARG(2)='YES'	
	DARG(1:3)=0.0
	Call SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RMSWFE,PixArray)

	command='PERTURB'	
	CARG(1)='ELEMENT'
	IARG(1)=23
        DARG(1:2)=1d-04; DARG(3:6)=0d0
	Call SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     &              OPDMat,RaySpot,RMSWFE,PixArray)

	nVarElts=1
        elts_arr(1)=23
        dof_arr(1:2)=1; dof_arr(3:7)=0
	total_itrs=20
	Call nls_optim_dvr(nVarElts,elts_arr,dof_arr,
     &                 WFE_TARGET,total_itrs,return_flag)

	stop
	End Program dopt_main
