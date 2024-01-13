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
#include "macos_definitions.h"

        Module macos_mod
        use traceutil_mod, only: RayStateStr
        use Constants, only: smallestReal => EPS   ! used in macos_cmd_loop.inc

        Implicit none
        Save

        Logical, dimension(:), allocatable :: ifRegrid,ifIntSrfSet,
     &                                        zernUsrOpt
        Integer, dimension(:), allocatable :: nDrawElt
        SREAL, dimension(:,:,:), allocatable :: Cmatrix,TEsp
        Real*8, dimension(:), allocatable :: dxMin,dxMax,DXavg,dxDev
        ! Real*8, parameter :: smallestReal=1.0d-15
        Integer, parameter :: CohImg=-1,InCohImg=1
        Integer :: ImgMode
        Real*8 :: bmRad  ! beam radius
        Character(len=2) :: cUnit

#ifdef CMACOS
        SREAL, dimension(:,:), allocatable :: PixArray
        Real*8, dimension(:,:), allocatable :: RaySpot,OPDMat,
     &                          OPDSave,OPDbuf
	Logical :: ifPGColor

        Integer :: OPDElt,OPDEltSave
#endif

	LOGICAL ifLoad,ifTrace,ifCmt,ifPerturb,ifPol,ifVecDif3,
     &    ifPlot,ifBuild,ifPropagate,ifOld,ifRayPrint,ifEFLCalc,
     &    ifGrid,ifOPDCalc,ifPixLoc,ifAdd,ifPertRead,ifAddCA,
     &    ifDMBuild,ifDraw,ifCentroid,ifFilterData,ifPixArray,ifLNegOK,
     &    ifStopSet,ifBeamCentered,ifLocalCoord,ifGraphics,ifMacro,
     &    ifUDSinit,ifPlot_sv,ifObjImgSet,LRtn_FLG,
     &    firstOPDCalc,multiOPDCalc,firstIntCalc

	Character(len=MacosValLen) :: ObjIntFn
	LREAL, dimension(:,:), allocatable :: ObjInt
        Character(len=80) :: FitsFilename

	Integer, parameter :: ZernCircular=1,ZernAnnular=2,ZernHex=3
	Integer :: zcoType,quarter,m_err_macos
	Logical :: macos_realloc = .true.

	Type(RayStateStr), pointer :: rayState


	Contains

	! --------------------------------------------------------------
        Subroutine macos_mod_init
        use param_mod

	! Dealloocate all dynamic arrays if already allocated.
        IF (allocated(ifRegrid)) THEN
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
     &               OPDbuf,stat=m_err_macos)
          if (m_err_macos /= 0)
     &      call macos_memory_failure(
     &          'macos_mod_init: deallocate failed!')
#endif
	END IF

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

	! SREAL arrays
	allocate(Cmatrix(7,mCm,bRay), TEsp(6,6,mElt),
     &	         stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')

        ! Real*8 arrays
        allocate(dxMin(0:mElt),dxMax(0:mElt),
     &           dxAvg(0:mElt), dxDev(0:mElt),
     &		 stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')

#ifdef CMACOS
        ! PixArray = SREAL, RaySpot and OPDMat = Real*8
        allocate(PixArray(mPix,mPix), RaySpot(mRay,2),
     &           OPDMat(mpts,mpts),OPDSave(mpts,mpts),
     &           OPDbuf(mpts,mpts),stat=m_err_macos)
        if (m_err_macos /= 0)
     &    call macos_memory_failure(
     &          'macos_mod_init: allocate failed!')
#endif

          ! Variable initialisation
          call macos_mod_init_vars()

	End Subroutine macos_mod_init

	! --------------------------------------------------------------

        subroutine macos_mod_init_vars()
          use Kinds
          implicit none

          cUnit           = ''
          FitsFilename    = ''

          ObjIntFn        = ''
          ! ObjInt(:,:)          ! allocated in utilsub.F -- ToDo: Why declare here???

          ifRegrid(:)     = .FALSE.
          ifIntSrfSet(:)  = .FALSE.
          zernUsrOpt(:)   = .FALSE.

          nDrawElt(:)     = 0
          ImgMode         = 0

          zcoType         = 0
          quarter         = 0
          m_err_macos     = 0

          Cmatrix(:,:,:)  = 0e0_pr   ! SREAL
          TEsp(:,:,:)     = 0e0_pr   ! SREAL

          dxMin(:)        = 0e0_pr
          dxMax(:)        = 0e0_pr
          dxAvg(:)        = 0e0_pr
          dxDev(:)        = 0e0_pr

          bmRad           = 0e0_pr   ! beam radius

          ! for SMACOS -- ToDo: need to combine SMACOS with MACOS init
#ifdef CMACOS
          PixArray(:,:)   = 0e0_pr
          RaySpot(:,:)    = 0e0_pr
          OPDMat(:,:)     = 0e0_pr
          OPDSave(:,:)    = 0e0_pr
          OPDbuf(:,:)     = 0e0_pr

          ifPGColor       = .FALSE.
          OPDElt          =  0
          OPDEltSave      =  0
#endif

          ifLoad          = .FALSE.
          ifTrace         = .FALSE.
          ifCmt           = .FALSE.
          ifPerturb       = .FALSE.
          ifPol           = .FALSE.
          ifVecDif3       = .FALSE.
          ifBuild         = .FALSE.
          ifPropagate     = .FALSE.
          ifPlot          = .FALSE.
          ifRayPrint      = .FALSE.
          ifPixLoc        = .FALSE.
          ifEFLCalc       = .FALSE.
          ifGrid          = .FALSE.
          ifAdd           = .FALSE.
          ifPixArray      = .FALSE.
          ifDMBuild       = .FALSE.
          ifCentroid      = .TRUE.   ! default: Centroid and not ChiefRay
          ifLNegOK        = .TRUE.   ! default: incl. Neg. Ray Lengths
          ifFilterData    = .FALSE.
          ifUDSinit       = .FALSE.
          ifOPDCalc       = .FALSE.
          ifPertRead      = .FALSE.
          ifAddCA         = .FALSE.
          ifObjImgSet     = .FALSE.
          LRtn_FLG        = .FALSE.

          ifStopSet       = .FALSE.
          ifBeamCentered  = .FALSE.
          ifLocalCoord    = .FALSE.

          ifMacro         = .FALSE.

          ifPlot_sv       = .FALSE.
          ifDraw          = .FALSE.
          ifGraphics      = .FALSE.

          firstOPDCalc    = .FALSE.
          multiOPDCalc    = .FALSE.

          ifEFLCalc       = .FALSE.
          firstIntCalc    = .FALSE.

          ! ifOld  =.FALSE.            ! note: it set to FALSE, SMACOS cannot load an Rx (ToDo)

        end subroutine macos_mod_init_vars

	! --------------------------------------------------------------

#ifdef CMACOS
	Subroutine set_macos_problem_size(problem_size)
        use param_mod

	Integer :: problem_size
	Integer, parameter :: default_problem_size = 128
	Character(len=MacosCharLen) :: input_buf
	Logical :: use_default = .false.

        print*,' '
        write(*, "(' ',A)"), repeat('*',72)
        write(*, "(' *',70x,'*')")
        write(*, "(' * ',A,' *')") MACOS_ACRONYM
        write(*, "(' * Copyright (C) 1995-',I4,', California Institute of Technology',10X,'*')"), MACOS_VERSION_YEAR
        write(*, "(' * U.S. Government Sponsorship is Acknowledged.',25x,'*')")
        write(*, "(' *',70x,'*')")
        write(*, "(' * ',20x,'MACOS Version: ',A7,27X,'*')"), MACOS_VERSION
        write(*, "(' *',70x,'*')")
        write(*, "(' ',A)"), repeat('*',72)

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
