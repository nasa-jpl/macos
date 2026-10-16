C***********************************************************************
C       Begin file cfiles_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module cfiles_mod
	use param_mod

	Implicit none
	Save

	CHARACTER(len=MacosCharLen) ::  filnam,infil,outfil,cmatfil,
     &                plotfil, mfil,pertfil,nasfil,imgfil,hdrfil,
     &                srffil,filtfil,chBuf,apodType

	CHARACTER(len=MacosCharLen),dimension(:),allocatable :: UDSrfFile
!-->  Was CHARACTER*24: a GridFile=/AmplFile= name longer than 24 chars
!     was TRUNCATED on assignment, so the open failed ("does not exist")
!     and the surface went silently flat -- the ifo-l2 case was a 25-char
!     name, one over.  Use the same length as UDSrfFile above.  SAVE
!     writes these with TRIM + "2x,A", so the emitted file is unchanged.
        CHARACTER(len=MacosCharLen), dimension(:), allocatable ::
     &                                             GridFile,AmplFile
!<--

	Integer, parameter :: stepApod  = 1,
     &                        linApod   = 2,
     &                        gaussApod = 3

        Integer :: iApodType,m_err_cfile

        Contains

        ! --------------------------------------------------------------
        Subroutine cfiles_mod_init

        if (allocated(UDSrfFile)) then
          deallocate(UDSrfFile,GridFile,AmplFile, stat=m_err_cfile)
          if (m_err_cfile /= 0)
     &      call macos_memory_failure
     &               ('cfile_mod_init: deallocate failed!')
        end if

        ! Character arrays
        allocate(UDSrfFile(mElt),GridFile(mElt),AmplFile(mElt),
     &                     stat=m_err_cfile)
        if (m_err_cfile /= 0)
     &      call macos_memory_failure
     &               ('cfile_mod_init: allocate failed!')

          ! Initalise
          CALL cfiles_mod_init_vars()

        End Subroutine cfiles_mod_init


        ! --------------------------------------------------------------
        SUBROUTINE cfiles_mod_init_vars()

          filnam   = ''
          infil    = ''
          outfil   = ''
          cmatfil  = ''
          plotfil  = ''
          mfil     = ''
          pertfil  = ''
          nasfil   = ''
          imgfil   = ''
          hdrfil   = ''
          srffil   = ''
          filtfil  = ''
          chBuf    = ''
          apodType = ''

          UDSrfFile(:) = ''
          GridFile(:)  = ''
          AmplFile(:)  = ''

        END SUBROUTINE cfiles_mod_init_vars


	End Module cfiles_mod
