C***********************************************************************
C       Begin file macos_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C
C This module defines variables used by smacos.F, macos.F and other
C application code, such as dda.f and nng.f, etc.
C

#include "realtype.h"
 	
 	Module macos_mod
	use traceutil_mod, only : RayStateStr
	Implicit none
	Save

	Logical, dimension(:), allocatable :: ifRegrid,ifIntSrfSet,
     &			                      zernUsrOpt
        Integer, dimension(:), allocatable :: nDrawElt
        SREAL, dimension(:,:,:), allocatable :: Cmatrix,TEsp
        Real*8, dimension(:), allocatable :: dxMin,dxMax, 
     &	  		      DXavg,dxDev
	Real*8, parameter :: smallestReal=1.0d-15
	Integer, parameter :: CohImg=-1,InCohImg=1
	Integer :: ImgMode
	Real*8 :: bmRad  ! beam radius
	Character(len=2) :: cUnit

#ifdef CMACOS
        SREAL, dimension(:,:), allocatable :: PixArray
        Real*8, dimension(:,:), allocatable :: RaySpot,OPDMat,
     &                   OPDSave,OPDbuf,RaySpotSave,RaySpotBuf,
     &                   RefOPD
	Logical :: ifPGColor
	
        Integer :: OPDElt,OPDEltSave,SpotElt,SpotEltSave
#endif

	LOGICAL ifLoad,ifTrace,ifCmt,ifPerturb,ifPol,ifVecDif3,
     &    ifPlot,ifBuild,ifPropagate,ifOld,ifRayPrint,ifEFLCalc,
     &    ifGrid,ifOPDCalc,ifPixLoc,ifAdd,ifPertRead,ifAddCA,
     &    ifDMBuild,ifDraw,ifCentroid,ifFilterData,ifPixArray,ifLNegOK,
     &    ifStopSet,ifBeamCentered,ifLocalCoord,ifGraphics,ifMacro,
     &    ifUDSinit,ifPlot_sv,ifObjImgSet,LRtn_FLG,
     &    firstOPDCalc,multiOPDCalc,firstSPOTCalc,multiSPOTCalc,
     &    firstIntCalc

	Character(len=MacosValLen) :: ObjIntFn
	LREAL, dimension(:,:), allocatable :: ObjInt
        Character(len=80) :: FitsFilename

	Integer, parameter :: ZernCircular=1,ZernAnnular=2,ZernHex=3
	Integer :: zcoType,quarter,m_err_macos
	Logical :: macos_realloc = .true.

	Type(RayStateStr), pointer :: rayState


	Contains

	Subroutine macos_mod_init
        use param_mod

	! Dealloocate all dynamic arrays if already allocated.
        If (allocated(ifRegrid)) Then
	  deallocate(ifRegrid, ifIntSrfSet,
     &		     nDrawElt,zernUsrOpt,
     &	             Cmatrix, TEsp,   
     &		     dxMin, dxMax, DXavg, dxDev,
     &		     stat=m_err_macos)
          if (m_err_macos /= 0)
     &      call macos_memory_failure(
     &          'macos_mod_init: deallocate failed!')

#ifdef CMACOS
	  deallocate(PixArray,RaySpot,OPDMat,OPDSave,
     &               OPDbuf,RaySpotSave,RaySpotBuf,
     &               RefOPD,stat=m_err_macos)
          if (m_err_macos /= 0)
     &      call macos_memory_failure(
     &          'macos_mod_init: deallocate failed!') 		    
#endif
	End If

	! Logical arrays
	allocate(ifRegrid(mElt), ifIntSrfSet(mElt), 
     &           zernUsrOpt(mElt),stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &		'macos_mod_init: allocate failed!')

	! Integer arrays
	allocate(nDrawElt(mDrawRay), 
     &		 stat=m_err_macos)
	if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')
	nDrawElt=0  ! very important -jzlou

	! SREAL arrays
	allocate(Cmatrix(7,mCm,bRay), TEsp(6,6,mElt),
     &	         stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')
	Cmatrix=0.0; TEsp=0.0

	! Real*8 arrays
	allocate(dxMin(0:mElt),dxMax(0:mElt),
     &           dxAvg(0:mElt), dxDev(0:mElt),
     &		 stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')
	dxMin=0d0; dxMax=0d0; dxAvg=0d0; dxDev=0d0

#ifdef CMACOS
	! PixArray = SREAL, RaySpot and OPDMat = Real*8
	allocate(PixArray(mPix,mPix), RaySpot(mRay,2), 
     &           OPDMat(mpts,mpts),OPDSave(mpts,mpts),
     &           OPDbuf(mpts,mpts),RaySpotSave(mRay,2),
     &           RaySpotBuf(mRay,2),RefOPD(mpts,mpts),
     &           stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')
	PixArray=0.0; RaySpot=0d0; OPDMat=0d0; OPDSave=0d0
        RaySpotSave=0d0; RefOPD=0d0

	ifPGColor=.false. ! default 
#endif
	End Subroutine macos_mod_init	

#ifdef CMACOS
	Subroutine set_macos_problem_size(problem_size)
        use param_mod

	Integer :: problem_size
	Integer, parameter :: default_problem_size = 128
	Character(len=MacosCharLen) :: input_buf
	Logical :: use_default = .false.

	print*,' '
	print*,
     &  '******************'//
     &  '******************************************************'
	print*,
     &  '*                 '//
     &  '                                                     *'
	print*,
     &  '* Modeling and Analysis for Controlled Optical Systems '//
     &  'Program (MACOS) *'
        print*,
     &  '* MACOS Version 3.34                                '//
     &  '                   *' 
        print*,
     &  '* Copyright (C) 1995-2019, California Institute of T'//
     &  'echnology          *' 
        print*,
     &  '* U.S. Government Sponsorship Is Acknowledged.      '//
     &  '                   *' 
        print*, 
     &  '*                 '//
     &  '                                                     *'
	print*,
     &  '******************'//
     &  '******************************************************'

	!problem_size = default_problem_size
	problem_size = -1
  15    print*,' '
        print*,'Please select a model size to run MACOS.'
        print*,
     &  'Valid model sizes are 128, 256, 512, 1024, 2048, 4096, 8192.'
!       print*,'Hit Return key to use default model size ',
!    &    default_problem_size, '.'
        write(*,19) ' Hit Return key to use default model size ',
     &              default_problem_size, '.' 
  19    format(A,I3,A)
        print*,' '
        write(*,14) ' MACOS model size: '
  14    format(A,$)
        read(*,16,end=17,err=17) input_buf
  16    format(A)
  17    continue
        read(input_buf,*,end=18) problem_size
  18    if (problem_size==-1) then
	  problem_size = default_problem_size
          use_default = .true.
        end if
        if (problem_size /= 128 .and. problem_size /= 256
     &      .and. problem_size /= 512 .and. 
     &      problem_size /= 1024 .and. 
     &      problem_size /= 2048 .and. 
     &      problem_size /= 4096 .and.
     &      problem_size /= 8192) then
          print*,'** Invalid problem size ', problem_size
          go to 15
        end if

        if (use_default) then
          print*,' Using default model size ', problem_size
        else
          print*,' Using specified model size ', problem_size
        end if
	End Subroutine set_macos_problem_size
#endif


#ifdef CSMACOS
	Subroutine smacos_set_model_size(model_size)
        use param_mod

	interface
	  Subroutine macos_init_all(psize)
	  implicit none
	  integer, optional :: psize
	  End Subroutine macos_init_all
	end interface

	Integer :: model_size
	if (model_size /= 128 .and. model_size /= 256 .and.
     &	    model_size /= 512 .and.
     &      model_size /= 1024 .and. model_size /= 2048)
     &  then
	  write(*,*) '**SMACOS Error: Invalid model size: ',
     &      model_size 
	  stop
	end if
	!
	call macos_init_all(model_size)
	macos_realloc = .true.  ! so smacos will reallocate
	End Subroutine smacos_set_model_size
#endif


	Real Function MACOS_Timer()
	real dtime
	real, save :: time(2)
	MACOS_Timer = dtime(time)
	End Function MACOS_Timer	

 	End Module macos_mod
