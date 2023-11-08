C***********************************************************************
C       Begin file realtype.h
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#define LREAL   REAL*8
! LREAL replaced RealType in source files  -- 09-2007, -jzlou

#define SREAL   REAL*8
! SREAL replaced REAL*4 in source files    -- 09-2007, -jzlou

#define realRecLen  8
  ! 'length' of unit of each real value, used in open Fortran direct access file.

#define MacosCharLen   256
#define MacosVarNamLen 256
#define MacosValLen    220

#ifdef SUNOS
#define  NoConOpt
#endif

#ifdef LNXOS32
#define LNXOS
#endif

#ifdef LNXOS64
#define LNXOS
#endif

#ifndef BUILD_LOC
#define BUILD_LOC ''
#endif
#ifndef SVN_REV
#define SVN_REV ''
#endif
