C***********************************************************************
C       Begin file lohpars_mod.F 
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C***********************************************************************
C This module replaces lohpars.cmn
C***********************************************************************

#include "realtype.h"

	Module lohpars_mod

C	Begin Source Definition Include File - 'lohpars.cmn' - comp v1.4
c sws 9/21/93

c bat_task:
c        1 read from console
c        2 read from .jou file

c jou_task:
c        1 journal mode off
c        2 write commands to .jou file

        INTEGER bat_task,jou_task

c exec (LQW) stack
	INTEGER, parameter :: LQWSTACK=10

c must have a common line size CVAR, var_store
c collection of command words, 1/line
	INTEGER, parameter :: cvarlen=MacosCharLen 
!-->  Tokens per command line: the capacity of the MacosCharLen line (one
!     character and a delimiter each).  Was 10, written unbounded -- an
!     11th token overran var_store (PLAN_CONSOLIDATION 1.2 (aaa)).  A line
!     that still counts more (empty ';' tokens) is refused by the parsers.
	INTEGER, parameter :: mLohTok=MacosCharLen/2
        CHARACTER(len=MacosCharLen) :: var_store(mLohTok,LQWSTACK)
!<--

c dummy not passed in commons
        CHARACTER(len=MacosCharLen) :: CVAR, 
     &	         cmdArg,  ! to hold one command arg read in PARSE_LOH
     &  	 cmdline  ! to hold entire command line

c read_len  number of commands in var_store
c read_cur  current position in var_store
c pstack    pointer to current stack 
c ustack    current unit number for current file
        INTEGER read_len(LQWSTACK),read_cur(LQWSTACK),pstack,ustack
!-->  The token READ_LOH consumed last, so UNREAD_LOH can hand it back
!     (PLAN 3.3, CALIB's yes/no prompts; macosio.F).  last_stk=0: none.
        INTEGER :: last_cur=0, last_len=0, last_stk=0
!<--
	Logical, save :: RX_Flg=.FALSE.

	End Module lohpars_mod

