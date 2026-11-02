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

      subroutine reinitialise_variables()
        use src_mod
        use cfiles_mod
        use elt_mod
        use macos_mod
        use traceutil_mod
        use dopt_mod

        implicit none
        ! ---------------------------
        ! call src_mod_init
        ! call cfiles_mod_init
        ! call elt_mod_init
        ! call macos_mod_init
        ! call ray_mod_init(mElt,mRay)
        ! call dopt_init

!-->    PLAN_CONSOLIDATION 1.2 (ooo), 2026-10-24: POLARIZATION is a session
!       mode and survives a load (PLAN 4 item 2, 2026-10-11), but the source
!       Jones vector Ex0/Ey0 it was switched on with lives in src_mod and was
!       zeroed by src_mod_init_vars on every load -- the mode survived with NO
!       field (an INTENSITY after polarization-then-load summed to zero; found
!       by macos-fd on the tg_psi_dm templates).  Keep Ex0/Ey0 across the load
!       when the mode is in force; a deck's own Ex0Ey0= line, parsed after
!       this, still overrides.
        block
          complex(kind(Ex0)) :: svEx0, svEy0
          logical :: keepPol
          keepPol = ifPol
          svEx0 = Ex0; svEy0 = Ey0
          call src_mod_init_vars()
          call cfiles_mod_init_vars()
          call elt_mod_init_vars()
          call macos_mod_init_vars()
          call dopt_init_vars()
          call ray_mod_init(mElt,mRay)
          if (keepPol) then
            Ex0 = svEx0; Ey0 = svEy0
          end if
        end block
!<--

      end subroutine reinitialise_variables


      subroutine initialise_all()
        use src_mod
        use cfiles_mod
        use elt_mod
        use macos_mod
        use traceutil_mod
        use dopt_mod

        implicit none
        ! ---------------------------
        call src_mod_init
        call cfiles_mod_init
        call elt_mod_init
        call macos_mod_init
        call opd_ref_init_opts  ! session option, before ray_mod_init
        call ray_mod_init(mElt,mRay)
        call dopt_init

      end subroutine initialise_all


      Subroutine macos_init_all(psize)
        use param_mod
        use src_mod
	use cfiles_mod
	use elt_mod
	use macos_mod
	use traceutil_mod

        implicit none

        integer, optional :: psize

        ! initialise via "macos_param.txt"
        if (present(psize)) then
          call param_mod_init(psize)
        else
          call param_mod_init
        end if

        call initialise_all

        ! call src_mod_init
	! call cfiles_mod_init
	! call elt_mod_init
	! call macos_mod_init
	! call ray_mod_init(mElt,mRay)
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
