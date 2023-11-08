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
        CHARACTER*24, dimension(:), allocatable :: GridFile,AmplFile

	Integer, parameter :: stepApod=1,linApod=2,gaussApod=3
        Integer :: iApodType,m_err_cfile

        Contains

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

        End Subroutine cfiles_mod_init

	End Module cfiles_mod
