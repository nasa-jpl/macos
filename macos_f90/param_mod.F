C***********************************************************************
C    Begin file parm_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship is acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C
C mElt = number of surfaces 
C mEm1 = number of segments
C mRay = max numberof rays
C bRay = max numberof rays for sensitivity calc.
C mWF =  max number diffraction planes
C mIntSrf = max number of interploation surfaces
C mDP = number of data points on each interploation surface
C mdttl = diffraction grid size
C 

#include "realtype.h"

	module param_mod
	Implicit none
	Save

	! Application Problem Size
        Integer :: macos_model_size

	Integer, parameter :: SellmeierGlass=1,SchottGlass=2,
     &                        HartmannGlass=3

	INTEGER :: mElt, mEm1, mRay, bRay, mCm, mpts, mdttl, md2,
     &    mWF, mPix, mIntSrf, mDP, LIWORK, mDMdof, mDrawRay, 
     &    mDrawElt, mObs, mWavelen, mStarType, mSegMirs, mUDS, 
     &    mUDScoef, mGridSrf, mGridMat, mGlass, mCplxObs, mObsMat,
     &    mAmplMat, mAmplSrf

	SREAL, parameter :: fzero_ = 0.0
        Real*8, parameter :: dzero_ = 0.0d0
        Complex*16, parameter :: czero_ = (0d0,0d0)
!	Character(len=MacosCharLen) :: macos_env_var

!       Set parameter strings based on preprocessor values for
!       the build location of MACOS and current subversion revision.
        CHARACTER(*), PARAMETER :: build_loc = 'BUILD_LOC'
        CHARACTER(*), PARAMETER :: svn_rev   = 'SVN_REV'
	
	Contains

	Subroutine param_mod_init(psize)
        implicit none

        Integer, optional :: psize  
        Integer :: ICLEN
!	character(len=64), parameter :: macos_param_file =
!    &                     'macos_param.txt' 
	character(len=MacosCharLen) :: macos_param_file 
	Integer :: macos_param_unit, ios
	Logical :: EXIST

	Namelist /macos_param128/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef
	Namelist /macos_param256/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef
	Namelist /macos_param512/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef
	Namelist /macos_param1024/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef
	Namelist /macos_param2048/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef
	Namelist /macos_param4096/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef
	Namelist /macos_param8192/ mElt, mpts, bRay, mWF,
     &         mIntSrf, mDP, mDrawRay, mObs,
     &         mWavelen, mStarType, mSegMirs, mGridSrf,
     &         mGridMat, mGlass, mCplxObs, mObsMat, mUDScoef


        If (present(psize)) then
           if (psize /= 128 .and. psize /= 256 .and.
     &         psize /= 512 .and. psize /= 1024 .and.
     &         psize /= 2048 .and. psize /= 4096 .and.
     &         psize /= 8192) then
              print*,'** SMACOS Error: Invalid problem size: ',
     &          psize
              stop
           end if
           macos_model_size = psize
        Else
           macos_model_size = 128 ! default problem size
        End If


        ! Problem-size dependent parameters
	! First read 'independent' parameters from namelist file

        macos_param_file='macos_param.txt'
        CALL find_macos_file('MACOS_HOME', macos_param_file)
#if 0
#ifdef SUNOS
        CALL GETENV('MACOS_HOME',macos_env_var)
        if (macos_env_var(1:2).EQ.'  ') then
          macos_param_file = 'macos_param.txt'
        else
	  print*,' '
          print*,'macos_env_var = ', 
     &      macos_env_var(1:ICLEN(macos_env_var))
          macos_param_file =
     &      macos_env_var(1:ICLEN(macos_env_var))
     &      //'/macos_param.txt'
          print*,'macos_param_file = ', 
     &      macos_param_file(1:ICLEN(macos_param_file))
        end if
#else
# ifdef LNXOS
        ! On Linux/lf95, GETENV doesn't work, and
        ! calling it crashes code. -jzlou
        macos_param_file = 'macos_param.txt'
# endif
# ifdef MSWIN
        ! On Windows/intel f95 compiler, GETENV doesn't work.
        macos_param_file = 'macos_param.txt'
# endif
#endif

	INQUIRE (FILE=macos_param_file,EXIST=EXIST)	
	if (.not. EXIST) then
	  ! The macos parameter file must be either in the
          ! MACOS_HOME directory or in the current execution
          ! directory. If not, it's an error, so quit.  
	  call macos_io_failure
     &           ('** macos parameter file '//
     &            macos_param_file//'does not exist, bye!') 
	end if
#endif

	if (ICLEN(macos_param_file).EQ.0) then
          call macos_io_failure
     &           ('** macos parameter file '//
     &            'macos_param.txt could not be found, bye!')
        else
           print *, ' Using parameter file "',
     &              macos_param_file(1:ICLEN(macos_param_file)), '"'
        end if

	macos_param_unit = 10
	Open(unit=macos_param_unit,file=macos_param_file,
     &       iostat=ios)
	if (ios /= 0) call macos_io_failure(
     &    '** Error open file '//macos_param_file)

	if (macos_model_size == 128) then
	  read(unit=macos_param_unit,nml=macos_param128,
     &         iostat=ios)
	else if (macos_model_size == 256) then
          read(unit=macos_param_unit,nml=macos_param256,
     &         iostat=ios)
	else if (macos_model_size == 512) then
          read(unit=macos_param_unit,nml=macos_param512,
     &         iostat=ios)
	else if (macos_model_size == 1024) then
          read(unit=macos_param_unit,nml=macos_param1024,
     &         iostat=ios)
	else if (macos_model_size == 2048) then
          read(unit=macos_param_unit,nml=macos_param2048,
     &         iostat=ios)
	else if (macos_model_size == 4096) then
          read(unit=macos_param_unit,nml=macos_param4096,
     &         iostat=ios)
	else if (macos_model_size == 8192) then
          read(unit=macos_param_unit,nml=macos_param8192,
     &         iostat=ios)
	end if
        if (ios /= 0) call macos_io_failure(
     &    '** Error read file '//macos_param_file)

	if (ios /= 0) call macos_io_failure(
     &     '** Error read file '//macos_param_file)
	
	mdttl = macos_model_size
        mPix = macos_model_size
        mEm1 = mElt-1; mRay=mpts**2+1; mCm=7+6*mElt
        md2 = mdttl**2
        LIWORK = (36*mDP)+1
        mDMdof = 3*mIntSrf*mDP; mDrawElt = 2*mElt
	!
        Close(macos_param_unit)
	End Subroutine param_mod_init


	SUBROUTINE find_macos_file(env_var_name, file_name)
        IMPLICIT NONE

        CHARACTER(*), INTENT(in)                   :: env_var_name
        CHARACTER(len=MacosCharLen), INTENT(inout) :: file_name

        INTEGER :: ICLEN

        CHARACTER(len=MacosCharLen)         :: env_var_value
        CHARACTER(len=MacosCharLen)         :: full_file_name
        LOGICAL   :: doesExist


#ifdef DEBUG
        print *, '==> find_macos_file(', env_var_name, ',', 
     &    file_name(1:ICLEN(file_name)), ')'
#endif

! 1) Look in directory specified by environment variable
#ifndef MSWIN
        ! On Windows/intel f95 compiler, GETENV doesn't work.
        IF (LEN(env_var_name).GT.0) then
           CALL GETENV(env_var_name,env_var_value)
           IF (ICLEN(env_var_value).GT.0) then
              ! See if the file exists
              full_file_name = 
     &            env_var_value(1:ICLEN(env_var_value)) // '/'
     &            // file_name(1:ICLEN(file_name))
#ifdef DEBUG
              print *, '  Looking for "',
     &                  full_file_name(1:ICLEN(full_file_name)), '"'
#endif
              INQUIRE (FILE=full_file_name,EXIST=doesExist)
              IF (doesExist) THEN
		! Found it and we're done.
                 file_name = full_file_name
                 RETURN
              ELSE
                 ! Warn user that it doesn't exist.
                 print *,'Warning: $', 
     &                   env_var_name(1:LEN(env_var_name)),
     &                   ' is defined to be ',
     &                    env_var_value(1:ICLEN(env_var_value)), ', but'
                 print *, '         failed to find file: ',
     &                    '$', env_var_name(1:LEN(env_var_name)), '/',
     &                    file_name(1:ICLEN(file_name))
              END IF
           END IF
        END IF
#endif

! 2) Look in the current directory
#ifdef DEBUG
        print *, '  Looking for "', file_name(1:ICLEN(file_name)), '"'
#endif
        INQUIRE (FILE=file_name,EXIST=doesExist)
        IF (doesExist) THEN
           ! Found it and we're done.
           RETURN
        END IF
        ! If it's not found in the current directory, that's normal so
        ! don't warn the user.

! 3) Look in the BUILD directory
        full_file_name = 
     &       build_loc // '/' // file_name(1:ICLEN(file_name))
#ifdef DEBUG
              print *, '  Looking for "',
     &                  full_file_name(1:ICLEN(full_file_name)), '"'
#endif
        INQUIRE (FILE=full_file_name,EXIST=doesExist)
        IF (doesExist) THEN
         ! Found it and we're done.
         file_name = full_file_name
         RETURN
        ELSE
         ! If it's not found in the build directory, then we're in trouble.
         ! Warn the user we couldn't find it anywhere.
         print *, 
     &   'Warning: Failed to find ', file_name(1:ICLEN(file_name))
         print *, 
     &   '  in the current directory or in the build directory:'
         print *, '        ', build_loc
         ! Return an empty string.
         file_name = ''
        END IF
        RETURN
        END SUBROUTINE find_macos_file

	End module param_mod
