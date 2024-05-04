C***********************************************************************
C	Begin file smacosutil.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	module smacosutil_mod

	Contains

	SUBROUTINE LoadStack(command,CARG,DARG,IARG,LARG,RARG)

C       These are SCOMP-specific variables

	IMPLICIT NONE
	LOGICAL LARG,LCMP
	CHARACTER(len=*) :: command
	CHARACTER(len=*), dimension(:) :: CARG
	CHARACTER*1000 STACK
	INTEGER, dimension(:) :: IARG
	INTEGER ls,K
	SREAL, dimension(:) :: RARG
	REAL*8, dimension(:) :: DARG
	COMMON /SCCOM/ls,STACK

	ls=0
	STACK(1:1000) = ' '

C  All commands load the command name
	CALL C2CW(command,ls,STACK,' ')

C  Load CARG(1) for some commands
	IF (LCMP(command,'LOAD',3).OR.LCMP(command,'OLD',2).OR.
     &	 LCMP(command,'EXECUTE',3).OR.LCMP(command,'LQW',3).OR.
     &	 LCMP(command,'SPCENTER',3).OR.LCMP(command,'SAO',3).OR.
     &	 LCMP(command,'PREAD',3).OR.LCMP(command,'LPREAD',3))
     &   THEN
      	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load CARG(1 and 2) for some commands
	IF (LCMP(command,'SAVE',3).OR.LCMP(command,'JOURNAL',3))
     &	THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  CALL C2CW(CARG(2),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load IARG(1) for other commands
	IF (LCMP(command,'SHOW',3).OR.LCMP(command,'RAY',3).OR.
     &	  LCMP(command,'PARTIALS',3).OR.LCMP(command,'OPD',3).OR.
     &    LCMP(command,'GBS',3).OR.
     &	  LCMP(command,'ROW',3).OR.LCMP(command,'COL',3).OR.
     &	  LCMP(command,'PARTIALS',3).OR.LCMP(command,'OPD',3).OR.
     &	  LCMP(command,'AM',2).OR.LCMP(command,'REAL',3).OR.
     &	  LCMP(command,'DMBUILD',3).OR.LCMP(command,'COORD',3).OR.
     &	  LCMP(command,'BUILD',3).OR.LCMP(command,'INTENSITY',3).OR.
     &	  LCMP(command,'GAIN',3).OR.LCMP(command,'LOG',3).OR.
     &    LCMP(command,'DVAR',4)) THEN
      	  CALL I2C(IARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load IARG(1) and stick a zero in after for RAY
C	IF (LCMP(command,'RAY',3)) ! these two lines were commented out by jzlou
C     	  CALL I2C(0,ls,STACK,' ') ! since it caused problem in *loop*.inc

C  Load IARG(1) and IARG(2) for COL
	IF (LCMP(command,'COL',3).OR.LCMP(command,'GBS',3)) THEN
      	  CALL I2C(0,ls,STACK,' ')
	  GO TO 110
	END IF

C  Load seed for SEED
	IF (LCMP(command,'SEED',3)) THEN
      	  CALL I2C(IARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load r0, wavelength and tilt participation for ATMOS
	IF (LCMP(command,'ATMOS',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL R2C(RARG(1),ls,STACK,' ')
	  CALL R2C(RARG(2),ls,STACK,' ')
	  CALL R2C(RARG(3),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  IF (.NOT.LCMP(CARG(1),'RAY',1)) THEN
	    CALL D2C(DARG(1),ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load IARG(1) and UDS file name for UDSINIT
	IF (LCMP(command,'UDSINIT',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load IARG(1) and coordinate type for SPOT
	IF (LCMP(command,'SPOT',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')   ! iElt
	  CALL C2CW(CARG(1),ls,STACK,' ')  ! BEAM, TOUT, ENTER
	  GO TO 110
	END IF

C  Load FDP command arguments
	IF (LCMP(command,'FDP',3)) THEN
          CALL I2C(IARG(1),ls,STACK,' ')
          CALL C2CW(CARG(1),ls,STACK,' ')
          GO TO 110
        END IF

C  Load 1st & last elt, plus plane, for DRAW
	IF (LCMP(command,'DRAW',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL I2C(IARG(2),ls,STACK,' ')
	  CALL I2C(IARG(3),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load CARG(1) with obscuration option for OBS
	IF (LCMP(command,'OBS',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load CARG(1) with coordinate type for PLOCate / WINdow
	IF (LCMP(command,'PLOCATE',4).OR.
     &	    LCMP(command,'WINDOW',3)) THEN
	  CALL C2C('TOUT',ls,STACK,' ')
	  CALL D2C(DARG(1),ls,STACK,' ')
	  CALL D2C(DARG(2),ls,STACK,',')
	  CALL D2C(DARG(3),ls,STACK,' ')
	  CALL D2C(DARG(4),ls,STACK,',')
	  CALL D2C(DARG(5),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load IARG(1), zero and coordinate type for PARTIALS
	IF (LCMP(command,'PARTIALS',3)) THEN
	  CALL I2C(0,ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load POL
	IF (LCMP(command,'POLARIZATION',3)) THEN
	  CALL D2C(DARG(1),ls,STACK,',')
	  CALL D2C(DARG(2),ls,STACK,' ')
	  CALL D2C(DARG(3),ls,STACK,',')
	  CALL D2C(DARG(4),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load ORS
	IF (LCMP(command,'ORS',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  IF (LCMP(CARG(1),'NO',1)) THEN
	    CALL C2C('NO',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load SRS
	IF (LCMP(command,'SRS',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL I2C(IARG(2),ls,STACK,' ')
	  IF (LCMP(CARG(1),'NO',1)) THEN
	    CALL C2C('NO',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C Added by GV, 26 Oct 1999
C  Load ray #, 1st & last elt, for SRAY
	IF (LCMP(command,'SRAY',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL I2C(IARG(2),ls,STACK,' ')
	  CALL I2C(IARG(3),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load FEXIT
	IF (LCMP(command,'FEXIT',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
c	  CALL D2C(DARG(1),ls,STACK,',')
c	  CALL D2C(DARG(2),ls,STACK,',')
c	  CALL D2C(DARG(3),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load STOP
	IF (LCMP(command,'STOP',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  IF (LCMP(CARG(1),'ELT',1)) THEN
	    CALL I2C(IARG(1),ls,STACK,' ')
	    CALL D2C(DARG(1),ls,STACK,',')
	    CALL D2C(DARG(2),ls,STACK,' ')
	  ELSE
	    CALL D2C(DARG(1),ls,STACK,',')
	    CALL D2C(DARG(2),ls,STACK,',')
	    CALL D2C(DARG(3),ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load PFP
	IF (LCMP(command,'PFP',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL D2C(DARG(1),ls,STACK,' ')
	  CALL D2C(DARG(2),ls,STACK,',')
	  CALL D2C(DARG(3),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load FFP
	IF (LCMP(command,'FFP',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL D2C(DARG(1),ls,STACK,',')
	  CALL D2C(DARG(2),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load CENTER
	IF (LCMP(command,'CEN',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL D2C(DARG(1),ls,STACK,',')
	  CALL D2C(DARG(2),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load MULTISPEC
	IF (LCMP(command,'MULTISPEC',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL I2C(IARG(2),ls,STACK,' ')
	  CALL I2C(IARG(3),ls,STACK,' ')
	  CALL R2C(RARG(1),ls,STACK,' ')
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load RFILT
	IF (LCMP(command,'RFILT',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load SFILT
	IF (LCMP(command,'SFILT',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load BEAM
	IF (LCMP(command,'BEAM',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  IF (.NOT.(LCMP(CARG(1),'U',1).OR.LCMP(CARG(1),'1',1))) THEN
	    CALL D2C(DARG(1),ls,STACK,' ')
	    CALL D2C(DARG(2),ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load PIXILATE
	IF (LCMP(command,'PIXILATE',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL I2C(IARG(2),ls,STACK,' ')
	  CALL R2C(RARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load LPIXILATE
	IF (LCMP(command,'LPIXILATE',3)) THEN
	  CALL I2C(IARG(2),ls,STACK,' ')
	  CALL R2C(RARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load COMPOSE
	IF (LCMP(command,'COMPOSE',3)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  CALL I2C(IARG(2),ls,STACK,' ')
	  CALL R2C(RARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load ADD
	IF (LCMP(command,'ADD',3)) THEN
	  IF (LCMP(CARG(1),'NO',1)) THEN
	    CALL C2C('NOPLOT',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load DADD
	IF (LCMP(command,'DADD',3)) THEN
	  IF (LCMP(CARG(1),'NO',1)) THEN
	    CALL C2C('NOPLOT',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load NOISE
	IF (LCMP(command,'NOISE',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  IF (LCMP(CARG(1),'GAUSS',2)) THEN
	    CALL D2C(DARG(1),ls,STACK,' ')
	  ELSE IF (LCMP(CARG(1),'READ',2)) THEN
	    CALL D2C(DARG(1),ls,STACK,' ')
	    CALL D2C(DARG(2),ls,STACK,' ')
	  ELSE IF (LCMP(CARG(1),'POISSON',2)) THEN
	    !CALL D2C(RARG(1),ls,STACK,' ')
	    CALL R2C(RARG(1),ls,STACK,' ')
	    CALL D2C(DARG(1),ls,STACK,' ')
	  ELSE IF (LCMP(CARG(1),'BACKGROUND',2)) THEN
	    !CALL D2C(RARG(1),ls,STACK,' ')
	    CALL R2C(RARG(1),ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load BLUR
	IF (LCMP(command,'BLUR',3)) THEN
	  CALL R2C(RARG(1),ls,STACK,' ')
	  IF (LCMP(CARG(1),'NO',1)) THEN
	    CALL C2C('NOPLOT',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load GBLUR
	IF (LCMP(command,'GBLUR',3)) THEN
	  CALL R2C(RARG(1),ls,STACK,',')
	  CALL R2C(RARG(2),ls,STACK,' ')
	  IF (LCMP(CARG(1),'NO',1)) THEN
	    CALL C2C('NOPLOT',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load STRETCH
	IF (LCMP(command,'STRETCH',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  GO TO 110
	END IF

C  Load PERTURB and LPERTURB
	IF (LCMP(command,'PERT',3).OR.LCMP(command,'LPERT',3)
     &      .OR.LCMP(command,'GPER',4)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ')
	  IF (IARG(1).NE.0) THEN
	   IF (LCMP(CARG(1),'NO',1).OR.
     &	       LCMP(CARG(1),'GLO',2)) THEN
	    CALL C2C('NO',ls,STACK,' ')
	    CALL D2C(DARG(1),ls,STACK,',')
	    CALL D2C(DARG(2),ls,STACK,',')
	    CALL D2C(DARG(3),ls,STACK,' ')
	    CALL D2C(DARG(4),ls,STACK,',')
	    CALL D2C(DARG(5),ls,STACK,',')
	    CALL D2C(DARG(6),ls,STACK,' ')
	   ELSE
	    CALL C2C('YES',ls,STACK,' ')
	    CALL D2C(DARG(1),ls,STACK,',')
	    CALL D2C(DARG(2),ls,STACK,',')
	    CALL D2C(DARG(3),ls,STACK,',')
	    CALL D2C(DARG(4),ls,STACK,',')
	    CALL D2C(DARG(5),ls,STACK,',')
	    CALL D2C(DARG(6),ls,STACK,' ')
	   END IF
	  END IF
	  IF (IARG(1).EQ.0) THEN
	    CALL D2C(DARG(1),ls,STACK,',')
	    CALL D2C(DARG(2),ls,STACK,',')
	    CALL D2C(DARG(3),ls,STACK,' ')
	    CALL D2C(DARG(4),ls,STACK,',')
	    CALL D2C(DARG(5),ls,STACK,',')
	    CALL D2C(DARG(6),ls,STACK,' ')
	  END IF
	  GO TO 110
	END IF

C  Load EXPORT
	IF (LCMP(command,'EXPORT',3)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  IF (LCMP(CARG(2),'NO',1)) THEN
	    CALL C2C('NO',ls,STACK,' ')
	  ELSE
	    CALL C2C('YES',ls,STACK,' ')
	  END IF
	  CALL C2CW(CARG(3),ls,STACK,' ')
	  CALL C2CW(CARG(4),ls,STACK,' ')
	  CALL C2CW(CARG(5),ls,STACK,' ')
	  CALL C2C('quit',ls,STACK,' ')
	  GO TO 110
	END IF

	IF (LCMP(command,'AVAR',4)) THEN
	  CALL I2C(IARG(1),ls,STACK,' ') ! Elt ID
	  CALL I2C(IARG(2),ls,STACK,',') ! 8 DOF (6 DOF+ROC+CC)
	  CALL I2C(IARG(3),ls,STACK,',')
	  CALL I2C(IARG(4),ls,STACK,',')
	  CALL I2C(IARG(5),ls,STACK,',')
	  CALL I2C(IARG(6),ls,STACK,',')
	  CALL I2C(IARG(7),ls,STACK,',')
	  CALL I2C(IARG(8),ls,STACK,',')
	  CALL I2C(IARG(9),ls,STACK,' ')
	  !CALL I2C(IARG(10),ls,STACK,',') ! Zernike coef range
	  !CALL I2C(IARG(11),ls,STACK,' ')
	  GO TO 110
	END IF

	IF (LCMP(command,'SETC',4)) THEN
	  CALL C2CW(CARG(1),ls,STACK,' ')
	  IF (LCMP(CARG(1),'OptTgt',6)) THEN
	    CALL C2CW(CARG(2),ls,STACK,' ')
	    IF (LCMP(CARG(2),'ZWF',3)) THEN
	      CALL I2C(IARG(1),ls,STACK,' ')
	      DO K=1,IARG(1)-1
	        CALL I2C(IARG(K+1),ls,STACK,',')
	      END DO
	      CALL I2C(IARG(IARG(1)+1),ls,STACK,' ')
	    END IF
	  ELSE
	    ! Other parameter types ...

	  END IF
	  GO TO 110
	END IF

C  Load SRT
        IF (LCMP(command,'SRT',3)) THEN
          CALL I2C(IARG(1),ls,STACK,' ')
          CALL D2C(DARG(1),ls,STACK,',')
          CALL D2C(DARG(2),ls,STACK,' ')
	  GO TO 110
        END IF

C  Load RETURN at end of each string
 110	CALL C2C('RETURN',ls,STACK,' ')
	END SUBROUTINE LoadStack

C***********************************************************************

	SUBROUTINE I2C(Int,ls,STACK,ctag)
	IMPLICIT NONE
	CHARACTER*8 C1
	CHARACTER*1 ctag
	CHARACTER*1000 STACK
	INTEGER i,nend,nstart,len,Int,len1,ls
 1	FORMAT(i8)
	WRITE(C1,1)Int
	i=0
 2      CONTINUE
	  i=i+1
	  IF ((C1(i:i).EQ.' ').AND.(i.LT.8)) GO TO 2
	nstart=i
	nend=8
	len=nend-nstart+1
	len1=ls+1
	ls=ls+len+1
	STACK(len1:ls)=c1(nstart:nend)//ctag
	RETURN
	END SUBROUTINE I2C

C***********************************************************************

	SUBROUTINE C2C(CARG,ls,STACK,ctag)
	IMPLICIT NONE
	CHARACTER*(*) CARG
	CHARACTER*1 ctag
	CHARACTER(len=MacosCharLen) :: C1
	CHARACTER*1000 STACK
	INTEGER i,nend,nstart,len1,ls
	INTEGER ICLEN

 1	FORMAT(i8)

csws	C1=CARG
c        i = 132
c        i = min( i, 1+iclen(carg) )
c        C1(1:i)=CARG(1:i)
c sws
c assume carg is just characters
	c1 = carg // ' '

	i=0
 2      CONTINUE
	  i=i+1
	  IF ((C1(i:i).EQ.' ').AND.(i.LT.MacosCharLen)) GO TO 2
	nstart=i
 3      CONTINUE
	  i=i+1
	  iF ((C1(i:i).NE.' ').AND.(i.LE.MacosCharLen)) GO TO 3
	nend=i-1
	len1=ls+1
	ls=ls+nend-nstart+2
	STACK(len1:ls)=c1(nstart:nend)//ctag
	RETURN
	END SUBROUTINE C2C

C***********************************************************************
C***********************************************************************

c sws new function, for words (terminated by ' ') embedded in long
c     strings

	SUBROUTINE C2CW(CARG,ls,STACK,ctag)
	IMPLICIT NONE
	CHARACTER*(*) CARG
	CHARACTER*1 ctag
	CHARACTER(len=MacosCharLen) :: C1
	CHARACTER*1000 STACK
	INTEGER i,nend,nstart,len1,ls
	INTEGER ICLEN


 1	FORMAT(i8)

csws	C1=CARG
        i = MacosCharLen-1
        i = min( i, iclen(carg) )
        C1(1:i+1)=CARG(1:i) // ' '

	i=0
 2      CONTINUE
	  i=i+1
	  IF ((C1(i:i).EQ.' ').AND.(i.LT.MacosCharLen)) GO TO 2
	nstart=i
 3      CONTINUE
	  i=i+1
	  iF ((C1(i:i).NE.' ').AND.(i.LE.MacosCharLen)) GO TO 3
	nend=i-1
	len1=ls+1
	ls=ls+nend-nstart+2
	STACK(len1:ls)=c1(nstart:nend)//ctag
	RETURN
	END SUBROUTINE C2Cw

C***********************************************************************
C***********************************************************************

	SUBROUTINE R2C(R1,ls,STACK,ctag)
	IMPLICIT NONE
	CHARACTER*1 ctag
	CHARACTER*14 C1,Char
	CHARACTER*1000 STACK
	INTEGER i,nstart,len1,ls
	SREAL R1
 1	FORMAT(e14.7)
	WRITE(C1,1)R1
	i=0
 2      CONTINUE
	  i=i+1
	  IF ((C1(i:i).EQ.' ').AND.(i.LT.14)) GO TO 2
	nstart=i
	len1=ls+1
	ls=ls-nstart+16
	STACK(len1:ls)=c1(nstart:14)//ctag
	RETURN
	END SUBROUTINE R2C

C***********************************************************************
C***********************************************************************

	SUBROUTINE D2C(D1,ls,STACK,ctag)
	IMPLICIT NONE
	CHARACTER*1 ctag
	CHARACTER*23 C1
	CHARACTER*1000 STACK
	INTEGER i,nstart,len1,ls
	REAL*8 D1
 1	FORMAT(d23.16)
	WRITE(C1,1)D1
	i=0
 2      CONTINUE
	  i=i+1
	  IF ((C1(i:i).EQ.' ').AND.(i.LT.23)) GO TO 2
	nstart=i
	len1=ls+1
	ls=ls-nstart+25
	STACK(len1:ls)=c1(nstart:23)//ctag
	RETURN
	END SUBROUTINE D2C

C***********************************************************************

	SUBROUTINE DQUERY(DVAR,DDEF,N,TEXT)
	CHARACTER*(*) TEXT
	CHARACTER*80 ASTR
	REAL*8 DVAR(N),DDEF(N)
  1	FORMAT(' ',A,$)
  2	FORMAT(A)
  3	WRITE(*,1)TEXT
	ASTR=' '
	DO 4 I=1,N
	  DVAR(I)=DDEF(I)
  4	CONTINUE
	READ(*,2,ERR=5) ASTR
	IF(ASTR.NE.' ')THEN
	  READ(ASTR,*,ERR=3)(DVAR(I),I=1,N)
	ENDIF
  5	CONTINUE
	RETURN
	END SUBROUTINE DQUERY

C***********************************************************************

	SUBROUTINE RQUERY(RVAR,RDEF,N,TEXT)
	CHARACTER*(*) TEXT
	SREAL RVAR(N),RDEF(N)
	CHARACTER*80 ASTR
  1	FORMAT(' ',A,$)
  2	FORMAT(A)
  3	WRITE(*,1)TEXT
	ASTR=' '
	DO 4 I=1,N
	  RVAR(I)=RDEF(I)
  4	CONTINUE
	READ(*,2,ERR=5) ASTR
	IF(ASTR.NE.' ')THEN
	  READ(ASTR,*,ERR=3)(RVAR(I),I=1,N)
	ENDIF
  5	CONTINUE
	RETURN
	END SUBROUTINE RQUERY

C***********************************************************************

	SUBROUTINE IQUERY(IVAR,IDEF,TEXT)
C	INTEGER INPUT ROUTINE IQUERY
C	TYPES OUT TEXT
C	IF INPUT IS CR, IVAR=IDEF
C	OTHERWISE, IVAR=INPUT

	CHARACTER*(*) TEXT
	INTEGER IVAR,IDEF

	CHARACTER*31 ASTR

  2	FORMAT(A)
  3	FORMAT(I31)
  4     FORMAT(' ',A,$)

100	WRITE(*,4)TEXT
	ASTR=' '
	READ(*,2,ERR=5) ASTR
	IF(ASTR.EQ.' ')THEN
	  IVAR=IDEF
	ELSE
	  READ(ASTR,3,ERR=100)IVAR
	END IF
  5	CONTINUE
	RETURN
	END SUBROUTINE IQUERY

C***********************************************************************

	SUBROUTINE CQUERY(CVAR,CDEF,TEXT)
	CHARACTER*(*) TEXT,CDEF
	CHARACTER(len=*) :: CVAR
  1	FORMAT(' ',A,' (',A,'): ',$)
  2	FORMAT(' ',A,$)
  3	FORMAT(A)
  4	CONTINUE
	CVAR=' '
	IF (CDEF.EQ.' ')THEN
	  WRITE(*,2)TEXT
	  READ(*,3) CVAR
	  IF(CVAR.EQ.' ')GO TO 4
	ELSE
	  WRITE(*,1)TEXT,CDEF
	  READ(*,3) CVAR
	  IF(CVAR.EQ.' ')CVAR=CDEF
	END IF
	RETURN
	END SUBROUTINE CQUERY

C***********************************************************************

	SUBROUTINE QUERY(CVAR,TEXT)
	CHARACTER*(*) TEXT
	CHARACTER(len=*) :: CVAR
  5     FORMAT(' ',A,$)
  6     FORMAT(A)
100     CONTINUE
	CVAR=' '
	WRITE(*,5)TEXT
	READ(*,6) CVAR
	IF(CVAR.EQ.' ')GO TO 100
	RETURN
	END SUBROUTINE QUERY

	End module smacosutil_mod

C***********************************************************************
C*******End File scomputil.f********************************************
C***********************************************************************
