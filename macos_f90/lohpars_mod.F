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
        CHARACTER(len=MacosCharLen) :: var_store(10,LQWSTACK)

c dummy not passed in commons
        CHARACTER(len=MacosCharLen) :: CVAR, 
     &	         cmdArg,  ! to hold one command arg read in PARSE_LOH
     &  	 cmdline  ! to hold entire command line

c read_len  number of commands in var_store
c read_cur  current position in var_store
c pstack    pointer to current stack 
c ustack    current unit number for current file
        INTEGER read_len(LQWSTACK),read_cur(LQWSTACK),pstack,ustack
	Logical, save :: RX_Flg=.FALSE.

	End Module lohpars_mod

