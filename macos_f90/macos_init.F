C***********************************************************************
C       Begin file macos_init.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C
C This subroutine reallocates most dynamic arrays in MACOS.
C
	Subroutine macos_init_all(psize)
        use param_mod
        use src_mod
	use cfiles_mod
	use elt_mod
	use macos_mod
	use traceutil_mod

        Implicit none

        integer, optional :: psize

        if (present(psize)) then
          call param_mod_init(psize)
        else
          call param_mod_init
        end if

        call src_mod_init
	call cfiles_mod_init
	call elt_mod_init
	call macos_mod_init
	call ray_mod_init(mElt,mRay)
        End Subroutine macos_init_all

  
        Subroutine macos_io_failure(err_msg)
        character(len=*) err_msg
        print*, err_msg
        stop
        End Subroutine macos_io_failure


        Subroutine macos_memory_failure(err_msg)
        character(len=*) err_msg
        print*, err_msg
        End Subroutine macos_memory_failure
