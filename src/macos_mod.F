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

!-->  The five CLI session modes (POLARIZATION, VECTOR/SCALAR, CENTROID/
!     CHIEFRAY, LNEG/NOLNEG, CIR/GIR) are set HERE, once, and survive
!     every load until changed (Dave 2026-10-11, PLAN_CONSOLIDATION 4
!     item 2; rule 6).  A LOAD names the non-default ones in force.
          ifPol           = .FALSE.
          ifVecDif3       = .FALSE.
          ifCentroid      = .FALSE.  ! default: CHIEF RAY (Dave 2026-08-27);
                                     ! centroid = opt-in (CENTROID cmd / api mode 0)
          ifLNegOK        = .TRUE.   ! default: incl. Neg. Ray Lengths
#ifdef CMACOS
          ifPGColor       = .FALSE.
#endif
!<--

	End Subroutine macos_mod_init

	! --------------------------------------------------------------

        subroutine macos_mod_init_vars()
          use Kinds
          implicit none

          cUnit           = ''
          FitsFilename    = ''

          ObjIntFn        = ''
          ! ObjInt(:,:)          ! allocated in utilsub.F (notes/hygiene.md, the census)

          ifRegrid(:)     = .FALSE.
          ifIntSrfSet(:)  = .FALSE.
          zernUsrOpt(:)   = .FALSE.

          nDrawElt(:)     = 0
          ImgMode         = 0

!-->      the ZCO / Zernike-fit basis: CIRCULAR unless the deck says
!         ZCOZernType= (2026-10-10).  It was reset to 0 here after the
!         LOAD handler set ZernCircular, so every fit after a load took
!         zern_fit's non-circular branch with type 0 and read
!         ZernAnnuRatio(ie_m=0) out of bounds -- the CALIB ZWF target
!         (WFE_ZMODE) then STOPped the host in nerror.  The rule-6 trap
!         again: the per-load reset must reset TO the default.
          zcoType         = ZernCircular
!<--
          quarter         = 0
          m_err_macos     = 0

          Cmatrix(:,:,:)  = 0e0_pr   ! SREAL
          TEsp(:,:,:)     = 0e0_pr   ! SREAL

          dxMin(:)        = 0e0_pr
          dxMax(:)        = 0e0_pr
          dxAvg(:)        = 0e0_pr
          dxDev(:)        = 0e0_pr

          bmRad           = 0e0_pr   ! beam radius

          ! for SMACOS -- combining with SMACOS's init is open (notes/hygiene.md)
#ifdef CMACOS
          PixArray(:,:)   = 0e0_pr
          RaySpot(:,:)    = 0e0_pr
          OPDMat(:,:)     = 0e0_pr
          OPDSave(:,:)    = 0e0_pr
          OPDbuf(:,:)     = 0e0_pr

          OPDElt          =  0
          OPDEltSave      =  0
#endif

          ifLoad          = .FALSE.
          ifTrace         = .FALSE.
          ifCmt           = .FALSE.
          ifPerturb       = .FALSE.
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
!         ifPol, ifVecDif3, ifCentroid, ifLNegOK, ifPGColor: SESSION
!         modes, set at allocation (macos_mod_init), NOT per load (Dave
!         2026-10-11, PLAN_CONSOLIDATION 4 item 2).
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

          ! ifOld  =.FALSE.            ! note: it set to FALSE, SMACOS cannot load an Rx -- so it stays commented out

        end subroutine macos_mod_init_vars

	! --------------------------------------------------------------

#ifdef CMACOS
	Subroutine set_macos_problem_size(problem_size)
        use param_mod

	Integer :: problem_size
	Integer, parameter :: default_problem_size = 256
	Character(len=MacosCharLen) :: input_buf
	Logical :: use_default
!-->  PLAN_CONSOLIDATION 1.2 (yy), Dave 2026-10-20: an answer that is
!     not one of the sizes is refused by name and asked again, at most
!     mTries times (a journal or a pipe can never spin here).
	Integer, parameter :: mTries = 5
	Integer :: iTry, ios, k, n
	Logical :: isDigits
!<--

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
        iTry = 0
  15    problem_size = -1
        use_default = .false.
        input_buf = ' '
        iTry = iTry + 1
        print*,' '
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
        read(*,16,iostat=ios) input_buf
  16    format(A)
!-->  (yy): end of input is no answer -- say so and stop (it was a
!     list-directed read of a blank or stale buffer: a forrtl crash, or
!     the last refused size re-read forever).  A blank line is the
!     default.  Anything else must be the digits of one of the sizes:
!     a word, '256abc', '256.5', '-1', '512 x' or a number too long for
!     an integer was a crash or a silent guess (ifx read '256.5' as 256,
!     '-1' as the default, '512 x' as 512); each is now refused.
        if (ios /= 0) then
          print*,' '
          print*,'** MACOS model size: no answer (end of input); ',
     &      'stopping.'
          stop 1
        end if
        do k = 1, LEN(input_buf)
          if (input_buf(k:k) == CHAR(9)) input_buf(k:k) = ' '
        end do
        input_buf = ADJUSTL(input_buf)
        n = LEN_TRIM(input_buf)
        if (n == 0) then
	  problem_size = default_problem_size
          use_default = .true.
        else
          isDigits = n <= 5
          do k = 1, n
            if (input_buf(k:k) < '0' .or. input_buf(k:k) > '9')
     &        isDigits = .false.
          end do
          if (isDigits) read(input_buf(1:n),'(I5)') problem_size
        end if
        if (problem_size /= 128 .and. problem_size /= 256
     &      .and. problem_size /= 512 .and.
     &      problem_size /= 1024 .and.
     &      problem_size /= 2048 .and.
     &      problem_size /= 4096 .and.
     &      problem_size /= 8192) then
          write(*,'(A,A,A,I1,A,I1,A)') ' ** Model size "',
     &      input_buf(1:MIN(n,40)), '" is not one of the sizes'//
     &      ' 128 ... 8192 (attempt ', iTry, ' of ', mTries, ')'
          if (iTry >= mTries) then
            write(*,'(A,I1,A)') ' ** MACOS model size: no valid '//
     &        'answer in ', mTries, ' attempts; stopping.'
            stop 1
          end if
          go to 15
        end if
!<--

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
!-->  One line at every LOAD naming the session modes NOT at their defaults
!     (Dave 2026-10-11, PLAN_CONSOLIDATION 4 item 2): they survive loads,
!     so a deck loaded after a CENTROID or a POLARIZATION runs with it --
!     say so.  Nothing printed when every mode is at its default.
        Subroutine ModesGet(pol, vecdif, centroid, lneg, cir)
!         The five session modes as the engine holds them -- the ONE source
!         for the LOAD notice and the bindings' modes_get (cir is .FALSE.
!         outside the CMACOS build, where ifPGColor does not exist).
          Implicit None
          Logical, Intent(Out) :: pol, vecdif, centroid, lneg, cir
          pol = ifPol; vecdif = ifVecDif3; centroid = ifCentroid
          lneg = ifLNegOK; cir = .FALSE.
#ifdef CMACOS
          cir = ifPGColor
#endif
        End Subroutine ModesGet

        Subroutine ModesInForceNote()
!         The five session modes only: OBS is NOT one (Dave 2026-10-28, a load
!         resets it to POSITIVE; macos_cmd_loop.inc LOAD says so).
          Implicit None
          Character(len=96) :: line
          Integer :: n
          Logical :: pol, vecdif, centroid, lneg, cir
          Call ModesGet(pol, vecdif, centroid, lneg, cir)
          line = ' '; n = 0
          If (pol) Then
            line = trim(line)//' POLARIZATION'; n = n + 1
            If (vecdif) Then
              line = trim(line)//', VECTOR'
            Else
              line = trim(line)//', SCALAR'
            End If
          End If
          If (centroid) Then
            If (n > 0) line = trim(line)//','
            line = trim(line)//' CENTROID'; n = n + 1
          End If
          If (.not. lneg) Then
            If (n > 0) line = trim(line)//','
            line = trim(line)//' NOLNEG'; n = n + 1
          End If
          If (cir) Then
            If (n > 0) line = trim(line)//','
            line = trim(line)//' CIR'; n = n + 1
          End If
          If (n > 0) WRITE(*,'(2A)') '  LOAD: session modes in force:',
     &        trim(line)
        End Subroutine ModesInForceNote
!<--


 	End Module macos_mod
