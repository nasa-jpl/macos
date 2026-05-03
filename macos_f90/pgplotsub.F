C***********************************************************************
C       Begin file pgplotsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************


#include "realtype.h"

	Module pgplot_mod
	!include 'pgplot.inc'

	Contains

C***********************************************************************

        SUBROUTINE GRAEND
C       Use this routine to terminate graphic output

        CALL PGEND

        RETURN
        END SUBROUTINE GRAEND

C***********************************************************************

C  Use this routine to initialize graphic output - such as '/xw' or '/gif'
        SUBROUTINE GRAINI
	use macos_vars_mod, only: nPgPanel

	if (nPgPanel==1) then
          CALL PGBEGIN(0,'?',1,1)
	else if (nPgPanel==2) then
          CALL PGBEGIN(0,'?',2,1) ! 1x2 horizontal panels
	  CALL PGPAP(9.0,0.4)
	else if (nPgPanel==3) then
          CALL PGBEGIN(0,'?',3,1) ! 1x3 horizontal panels
          CALL PGPAP(13.0,0.4)
	else if (nPgPanel==4) then
          CALL PGBEGIN(0,'?',2,2) ! 2x2 panels
	  CALL PGPAP(7.5,1.0)
	end if
	!print*,'*** nPgPanel =',nPgPanel

        !CALL PGSVP(0.0,2.0,0.0,1.0)
	!CALL PGSWIN (0.0,2.0,0.0,1.0)
	!CALL PGWNAD(0.0,1.0,0.0,1.0)
        !CALL PGENV(0.0,1.0,0.0,1.0,0,0)
        RETURN
        END SUBROUTINE GRAINI


C***********************************************************************

        SUBROUTINE CONTOUR(X_in,Y_in,Z_in,NROW,NCOL,NDROW,NDCOL,
     &  		   NCONTOURS,CXLABL,CZLABL,CTITLE)

        IMPLICIT NONE

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,NCONTOURS
	SREAL :: X_in(NROW),Y_in(NCOL),Z_in(NDROW,NDCOL)
        REAL*4 X(NROW),Y(NCOL),Z(NDROW,NDCOL),TRCONT(6),
     &         DX,DY,ZMIN,ZMAX,CONTOURS(16)

	X=X_in; Y=Y_in; Z=Z_in

        CALL PGBBUF

        DX=(X(NROW)-X(1))/FLOAT(NROW)
        DY=(Y(NCOL)-Y(1))/FLOAT(NCOL)
        TRCONT(1)=X(1)-DX
        TRCONT(2)=DX
        TRCONT(3)=0e0
        TRCONT(4)=Y(1)-DY
        TRCONT(5)=0e0
        TRCONT(6)=DY

C       SET PLOT DIMENSIONS

        ZMIN=1E31
        ZMAX=-1E31
        DO 2 I=1,NROW
          DO 1 J=1,NCOL
            IF (Z(I,J).LT.ZMIN)ZMIN=Z(I,J)
            IF (Z(I,J).GT.ZMAX)ZMAX=Z(I,J)
   1      CONTINUE
   2    CONTINUE

        IF (NCONTOURS.GT.16)NCONTOURS=16
        DO 3 I=1,NCONTOURS
          CONTOURS(I)=ZMIN+FLOAT(I-1)*(ZMAX-ZMIN)/FLOAT(NCONTOURS)
  3     CONTINUE

C       CALL PLOT ROUTINES

        CALL PGENV(X(1),X(NROW),Y(1),Y(NCOL),1,0)
        CALL PGLABEL(CXLABL,CZLABL,CTITLE)
        CALL PGCONT(Z,NDROW,NDCOL,1,NROW,1,NCOL,
     &  CONTOURS,NCONTOURS,TRCONT)

        CALL PGEBUF

        RETURN
        END SUBROUTINE CONTOUR


C***********************************************************************


        SUBROUTINE SPOTDIAG(DXYMAT,NROW,NCOL,NDROW,NDCOL,NBEAMS,
     &  		    CXLABL,CZLABL,CTITLE)

        IMPLICIT NONE

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,K,NBEAMS,ISYM(4),CI
        REAL*8 DXYMAT(NDROW,NDCOL)
        REAL*4 XYMAT(NDROW,NDCOL),XMIN,XMAX,YMIN,YMAX,PAD,DX,DY
        DATA ISYM/20,5,8,16/

        CALL PGBBUF

C       XYMAT = XYMAT(n_pts,2), storing corrds of points in 
C       the spot diagram. -jzlou	
	XYMAT(1:NROW,1:NCOL) = DXYMAT(1:NROW,1:NCOL)


C       SET PLOT DIMENSIONS

        XMIN=1E30; XMAX=-1E30
        YMIN=1E30; YMAX=-1E30

        IF (NBEAMS.EQ.1) THEN
          J=1
        ELSE
          J=3
        END IF

        DO I=1,NROW
          IF (XYMAT(I,J).LT.XMIN) XMIN=XYMAT(I,J)
          IF (XYMAT(I,J).GT.XMAX) XMAX=XYMAT(I,J)
        END DO

        IF (NBEAMS.EQ.1) THEN
          J=2
        ELSE
          J=4
        END IF

        DO I=1,NROW
          IF (XYMAT(I,J).LT.YMIN) YMIN=XYMAT(I,J)
          IF (XYMAT(I,J).GT.YMAX) YMAX=XYMAT(I,J)
        END DO

        DX=0.05e0*(XMAX-XMIN) 
        DY=0.05e0*(YMAX-YMIN)
        DX=MAX(DX,DY)
        XMIN=XMIN-DX; XMAX=XMAX+DX
        YMIN=YMIN-DX; YMAX=YMAX+DX

C       CALL PLOT ROUTINES

        CALL PGENV(XMIN,XMAX,YMIN,YMAX,1,0)
        CALL PGLABEL(CXLABL,CZLABL,CTITLE)
        CALL PGQCI(CI)

C       When nbeams=1, j=1 and k=2; when nbeams == 2, j=3 and k=4 
        I=NBEAMS; K=2*I; J=K-1
        CALL PGPOINT(NROW,XYMAT(1,J),XYMAT(1,K),ISYM(I))

        CALL PGEBUF

        RETURN
        END SUBROUTINE SPOTDIAG


C***********************************************************************
C	Not used in macos -jzlou
C***********************************************************************

        SUBROUTINE LINSPOTDIAG(XYMAT_in,NROW,NCOL,NDROW,NDCOL,NBEAMS,
     &  		       CXLABL,CZLABL,CTITLE)

        IMPLICIT NONE

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,K,NBEAMS,ISYM(4),CI
	SREAL :: XYMAT_in(NDROW,NDCOL)
        REAL*4 XYMAT(NDROW,NDCOL),XMIN,XMAX,YMIN,YMAX,PAD,DX,DY
        DATA ISYM/20,5,8,16/

	XYMAT = XYMAT_in

        CALL PGBBUF

C       SET PLOT DIMENSIONS

        XMIN=1E30
        XMAX=-1E30
        YMIN=1E30
        YMAX=-1E30
        DO 2 J=1,NCOL,2
          DO 1 I=1,NROW
            IF (XYMAT(I,J).LT.XMIN)XMIN=XYMAT(I,J)
            IF (XYMAT(I,J).GT.XMAX)XMAX=XYMAT(I,J)
   1      CONTINUE
   2    CONTINUE
        DO 4 J=2,NCOL,2
          DO 3 I=1,NROW
            IF (XYMAT(I,J).LT.YMIN)YMIN=XYMAT(I,J)
            IF (XYMAT(I,J).GT.YMAX)YMAX=XYMAT(I,J)
   3      CONTINUE
   4    CONTINUE

        DX=0.05e0*(XMAX-XMIN)
        DY=0.05e0*(YMAX-YMIN)
        XMIN=XMIN-DX
        XMAX=XMAX+DX
        YMIN=YMIN-DY
        YMAX=YMAX+DY

C       CALL PLOT ROUTINES

        !CALL PGENV(XMIN,XMAX,YMIN,YMAX,1,0)
        CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
        CALL PGLABEL(CXLABL,CZLABL,CTITLE)
        CALL PGQCI(CI)
        DO 5 I=1,NBEAMS
          K=2*I
          J=K-1
C         CALL PGSCI((I+1))
          CALL PGPOINT(NROW,XYMAT(1,J),XYMAT(1,K),ISYM(I))
  5     CONTINUE
C       CALL PGSCI(CI)

        CALL PGEBUF

        RETURN
        END SUBROUTINE LINSPOTDIAG


C***********************************************************************
C	SLICE displays array data as a 3D surface
C***********************************************************************

        SUBROUTINE SLICE(X,Y,Z,NROW,NCOL,NDROW,NDCOL,
     &  		 PX,PZ,CXLABL,CYLABL,CTITLE)

        IMPLICIT NONE

        CHARACTER*(*) CXLABL,CYLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J
 	SREAL :: X(NROW),Y(NCOL),Z(NDROW,NDCOL),
     &           PX(NROW),PZ(NROW),DX,DZ,ZMIN
        REAL*4 ZMAX,XSKEW,ZSKEW,XBMIN,XBMAX,ZBMIN,ZBMAX

C       SET PLOT DIMENSIONS

        XSKEW=0.35E0
        DX=(X(2)-X(1))*XSKEW
        XBMIN=X(1)
        XBMAX=X(NROW)+FLOAT(NCOL)*DX

        ZSKEW=0.45E0
        ZMIN=1E31
        ZMAX=-1E31
        DO 2 I=1,NROW
          DO 1 J=1,NCOL
            IF (Z(I,J).LT.ZMIN)ZMIN=Z(I,J)
            IF (Z(I,J).GT.ZMAX)ZMAX=Z(I,J)
   1      CONTINUE
   2    CONTINUE
        DZ=ZSKEW*(ZMAX-ZMIN)/FLOAT(NCOL)
        ZBMIN=ZMIN
        ZBMAX=ZMAX+FLOAT(NCOL)*DZ

C       CALL PLOT ROUTINES

        CALL PGENV(XBMIN,XBMAX,ZBMIN,ZBMAX,0,0)
        CALL PGLABEL(CXLABL,CYLABL,CTITLE)
        CALL PGSFS(1)
        CALL PGSCI(1)
        CALL PGSLICE(X,Y,Z,PX,PZ,NROW,NCOL,NDROW,NDCOL,
     &               DX,DZ,ZMIN)

        RETURN
        END SUBROUTINE SLICE

C***********************************************************************

        SUBROUTINE PGSLICE(X_in,Y_in,Z_in,PX_in,PZ_in,
     &                     NROW,NCOL,NDROW,NDCOL,DX_in,DZ_in,ZMIN_in)

        IMPLICIT NONE

        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,IP1,CI
	SREAL :: X_in(NROW),Y_in(NROW),Z_in(NDROW,NDCOL),
     &           PX_in(NROW),PZ_in(NROW),DX_in,DZ_in,ZMIN_in
        REAL*4 X(NROW),Y(NROW),Z(NDROW,NDCOL),PX(NROW),PZ(NROW),DX,DZ,
     &         XSTEP,ZSTEP,ZMIN,POLYX(4),POLYZ(4)

	X=X_in; Y=Y_in; Z=Z_in; PX=PX_in; PZ=PZ_in
	DX=DX_in; DZ=DZ_in; ZMIN=ZMIN_in

        CALL PGBBUF
        CALL PGQCI(CI)

        XSTEP=FLOAT(NCOL-1)*DX
        ZSTEP=FLOAT(NCOL)*DZ
        DO 5 J=NCOL,1,-1
          DO 3 I=1,NROW
            PX(I)=X(I)+XSTEP
            PZ(I)=Z(I,J)+ZSTEP
 3        CONTINUE
          CALL PGSCI(0)
          DO 4 I=1,(NROW-1)
            IP1=I+1
            POLYX(1)=PX(I)
            POLYX(2)=PX(I)
            POLYX(3)=PX(IP1)
            POLYX(4)=PX(IP1)
            POLYZ(1)=ZSTEP+ZMIN
            POLYZ(2)=PZ(I)
            POLYZ(3)=PZ(IP1)
            POLYZ(4)=ZSTEP+ZMIN
            CALL PGPOLY(4,POLYX,POLYZ)
 4        CONTINUE
          CALL PGSCI(CI)
          CALL PGLINE(NROW,PX,PZ)
          XSTEP=XSTEP-DX
          ZSTEP=ZSTEP-DZ
5       CONTINUE

        CALL PGEBUF

        RETURN
        END SUBROUTINE PGSLICE


C***********************************************************************
C   	Rendering a gray-scale or color image
C       Color representations were similar to those in Matlab		
C       -jzlou, 03/2008
C***********************************************************************

        SUBROUTINE GRAY(X_in,Y_in,Z_in,NROW,NCOL,NDROW,NDCOL,
     &  	        CXLABL,CZLABL,CTITLE)
	use macos_mod, only : ifPGColor
        IMPLICIT NONE

        INTEGER :: NROW,NCOL,NDROW,NDCOL,I,J
        !CHARACTER*(*) CXLABL,CZLABL,CTITLE
	REAL*8 :: X_in(NDROW),Y_in(NDCOL),Z_in(NDROW,NDCOL)
        CHARACTER(len=*) :: CXLABL,CZLABL,CTITLE
        !REAL*4 :: X(NDROW),Y(NDCOL)

        Real*4, dimension(:), allocatable, save :: X,Y
        Real*4, dimension(:,:), allocatable, save :: Z
        REAL*4 :: TRCONT(6),DX,DY,ZMIN,ZMAX

        Integer :: m_err_gray
#if 0
        REAL*4 :: X(NDROW),Y(NDCOL),Z(NDROW,NDCOL),TRCONT(6),
     &            DX,DY,ZMIN,ZMAX
#endif
	INTEGER, PARAMETER :: mColor=64
	Integer, Save :: ICILO,ICIHI,NDROW_sv
        REAL, Save :: CInd(3,mColor),CR(512),CG(512),CB(512)
	LOGICAL, Save :: first_entry=.true.

!	print*,'*** In GRAY: NROW,NCOL,NDROW,NDCOL =',
!    &    NROW,NCOL,NDROW,NDCOL

	If (first_entry .OR. (NDROW_sv .NE. NDROW)) Then
          If (allocated(X)) deallocate(X,Y,Z)
          allocate(X(NDROW),Y(NDCOL),Z(NDROW,NDCOL),
     &             stat=m_err_gray)
          If (m_err_gray /= 0) Then
            call macos_memory_failure('GRAY: allocate failed')
          End If
          NDROW_sv = NDROW
        End If

        If (first_entry) Then
	  ! Save default color representations for gray image
	  CALL PGQCIR(ICILO,ICIHI)
          Do J=ICILO,ICIHI
            CALL PGQCR(J,CR(J),CG(J),CB(J))
          End Do

	  ! Initialize gray and color indices
#include "mcolors.inc"
	  first_entry=.false.
	End If

	If (ifPGColor) Then
	  ! Set color representation
	  Do J=1,mColor
	    Call PGSCR(J-1+16,CInd(1,J),CInd(2,J),CInd(3,J))
	  End Do
	  ! Set color range
	  !CALL PGSCIR(0,mColor-1)
	  CALL PGSCIR(16,(mColor-1)+16)
	Else
	  ! Set default color representation for gray image
	  Do J=ICILO,ICIHI
            Call PGSCR(J,CR(J),CG(J),CB(J))
          End Do
	  ! Set color range
          CALL PGSCIR(ICILO,ICIHI) 
	End If	

	X(1:NROW)=X_in(1:NROW); Y(1:NCOL)=Y_in(1:NCOL)
        Z(1:NROW,1:NCOL)=Z_in(1:NROW,1:NCOL)

#if 0
	Open(unit=12,file='opd_chk.txt',status='replace')
	Open(unit=13,file='X_chk.txt',status='replace')
	Open(unit=14,file='X_in_chk.txt',status='replace')
        Do J=1,NCOL
          WRITE(12,*) Z(1:NROW,J)
	End Do
	Close(12)
        Do i=1,NROW
          WRITE(13,*) X(i)
        End Do
        Do i=1,NROW
          WRITE(14,*) X_in(i)
        End Do
	Close(13); Close(14)
#endif

        CALL PGBBUF

        DX=(X(NROW)-X(1))/FLOAT(NROW)
        DY=(Y(NCOL)-Y(1))/FLOAT(NCOL)
        TRCONT(1)=X(1)-0.5*DX
        TRCONT(2)=DX
        TRCONT(3)=0e0
        TRCONT(4)=Y(1)-0.5*DY
        TRCONT(5)=0e0
        TRCONT(6)=DY

!       print*,'***** X(NROW),X(1),FLOAT(NROW) =',
!    &          X(NROW),X(1),FLOAT(NROW)

C       SET PLOT DIMENSIONS

        ZMIN=1E31; ZMAX=-1E31
        DO I=1,NROW
          DO J=1,NCOL
            IF (Z(I,J).LT.ZMIN)ZMIN=Z(I,J)
            IF (Z(I,J).GT.ZMAX)ZMAX=Z(I,J)
          END DO
   	END DO

C       CALL PLOT ROUTINES

!        print*,'*** In GRAY: NROW,NCOL,NDROW,NDCOL =',
!    &    NROW,NCOL,NDROW,NDCOL
!        print*,'TRCONT =',TRCONT

	If (ifPGColor) Then
	  CALL PGSCI(1)
	Else
	  CALL PGSCI(ICILO)
	End If
	
        CALL PGENV(X(1),X(NROW),Y(1),Y(NCOL),1,0)
        CALL PGLABEL(CXLABL,CZLABL,CTITLE)

	If (ifPGColor) Then
 	  CALL PGIMAG(Z,NDROW,NDCOL,1,NROW,1,NCOL,
     &                ZMIN,ZMAX,TRCONT)
	  CALL PGWEDG('RI',2.2,5.0,ZMIN,ZMAX,'pixel value')
	Else
          CALL PGGRAY(Z,NDROW,NDCOL,1,NROW,1,NCOL,
     &                ZMIN,ZMAX,TRCONT)
	  CALL PGWEDG('RG',2.2,5.0,ZMIN,ZMAX,'pixel value')
	End If
	CALL PGSCH(1.0)

        CALL PGEBUF
!#endif
        RETURN
        END SUBROUTINE GRAY


C***********************************************************************
C	Called by other routine in this file. -jzlou
C***********************************************************************

        SUBROUTINE FREDDY(ARRAY,KX,NY,SIZE,ANGLE)
        INTEGER KX, NY
        REAL ARRAY(KX,NY), SIZE, ANGLE
C
C Draws isometric plot of array
C
        REAL FMAX,FMIN,DELTAX,DELTAY,DELTAV,SINE,PEAK,X,DX,HEIGHT
        INTEGER I,J,KI,KJ,NX,MX,MY,STEP,LEFT,RIGHT,IT,MN,INCX
        LOGICAL VISBLE
        COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
        MN = KX*NY
        NX = KX
C     Check array size:
        IF(NX.LT.2 .OR. NY.LT.2) RETURN
        FMAX = ARRAY(1,1)
        FMIN = FMAX
        DO 20 J=1,NY
            DO 10 I=1,NX
                FMIN = AMIN1(ARRAY(I,J),FMIN)
                FMAX = AMAX1(ARRAY(I,J),FMAX)
   10     CONTINUE
   20 CONTINUE
        DELTAX = SIZE/(NX+NY)
        SINE = SIN(ANGLE/58.)
        DELTAY = DELTAX*SINE
        HEIGHT = SIZE*(1.-ABS(SINE))
        DELTAV = HEIGHT
        FMAX = FMAX-FMIN
        IF(FMAX.LT.0.0001) FMAX = DELTAV
        DELTAV = DELTAV/FMAX
        MX = NX+1
        MY = NY+1
        STEP = MX
C
C Start PGPLOT buffering.
C
        CALL PGBBUF
C
C Work our way down the Y axis, THEN up the X axis,
C calculating the Y plotter coordinates for each
C column of the plot, doing the hidden-line suppression
C at the same time.
C
        DO 50 J=1,NY
            KJ = MY-J
            KI = 1
C               ( KI,KJ are coordinates of bottom of column)
            ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
   30     PEAK = ARRAY(KI,KJ)
   40     KI = KI+1
            KJ = KJ+1
            IF(KI.GT.NX .OR. KJ.GT.NY) GO TO 50
            ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
            IF(ARRAY(KI,KJ).GT.PEAK) GO TO 30
            IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
            GO TO 40
   50 CONTINUE
C
C Now to work our way up the X axis
C
        DO 80 I=2,NX
            KI = I
            KJ = 1
            ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
   60     PEAK = ARRAY(KI,KJ)
   70     KI = KI+1
            KJ = KJ+1
            IF(KI.GT.NX .OR. KJ.GT.NY) GO TO 80
            ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
            IF(ARRAY(KI,KJ).GT.PEAK) GO TO 60
            IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
            GO TO 70
   80 CONTINUE
C
C Draw a line along the bottom of the vertical faces
C
        CALL PGMOVE(DELTAX*(NX+NY-2), DELTAY*(MX))
        CALL PGDRAW(DELTAX*(NY-1),    DELTAY*2)
        CALL PGDRAW(0.0,              DELTAY*MY)
C
C Array is now ready for plotting.  IF a point is
C positive, THEN it is to be plotted at that Y
C coordinate; IF it is negative, THEN it is
C invisible, but at minus that Y coordinate (the point
C where the line heading towards it disappears has to
C be determined by finding the intersection of it and
C the cresting line).
C
C Plot rows:
C
        DO 110 J=1,NY,2
            KJ = MY-J
            DX = DELTAX*(J-2)
            X = DX+DELTAX
            CALL PGMOVE(X,DELTAY*(KJ+1))
            CALL PGDRAW(X,ARRAY(1,KJ))
            VISBLE = .TRUE.
            DO 90 I=2,NX
                RIGHT = I+NX*(KJ-1)
                LEFT = RIGHT-1
                IT = RIGHT
                X = DX+DELTAX*I
                CALL FREDGO(ARRAY,MN)
   90     CONTINUE
C
C Now at far END of row so come back
C
            KJ = KJ-1
            IF(KJ.LE.0) GO TO 170
            VISBLE = ARRAY(NX,KJ).GE.0.0
            DX = DELTAX*(NX+J)
            IF(VISBLE) CALL PGMOVE(DX-DELTAX,ARRAY(NX,KJ))
            DELTAX = -DELTAX
            DO 100 I=2,NX
                KI = MX-I
                LEFT = KI+NX*(KJ-1)
                RIGHT = LEFT+1
                IT = LEFT
                X = DX+DELTAX*I
                CALL FREDGO(ARRAY,MN)
  100     CONTINUE
C
            X = DX+DELTAX*NX
            IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(1,KJ))
            CALL PGDRAW(X,DELTAY*(KJ+1))
C               (set DELTAX positive for RETURN trip)
            DELTAX = -DELTAX
  110 CONTINUE
C
C Now DO the columns:
C as we fell out of the last DO-loop we DO the
C columns in ascending-X order
C
        INCX = 1
        KI = 1
C               (set DELTAX -ve since scanning R to L)
  120 DX = DELTAX*(KI+NY-1)
        DELTAX = -DELTAX
        X = DX+DELTAX
        CALL PGMOVE(X,ARRAY(1,1))
  130 VISBLE = .TRUE.
        DO 140 J=2,NY
            LEFT = KI+NX*(J-1)
            RIGHT = LEFT-NX
            IT = LEFT
            X = DX+DELTAX*J
            CALL FREDGO(ARRAY,MN)
  140 CONTINUE
C
C At far END, increment X and check still inside array
C
        KI = KI+INCX
        IF(KI.LE.0 .OR. KI.GT.NX) GO TO 180
        VISBLE = ARRAY(KI,NY).GE.0.0
        DELTAX = -DELTAX
        DX = DELTAX*(KI-2)
        X = DX+DELTAX
        IF(VISBLE) CALL PGMOVE(X,ARRAY(KI,NY))
        DO 150 J=2,NY
            KJ = MY-J
            RIGHT = KI+NX*(KJ-1)
            LEFT = RIGHT+NX
            IT = RIGHT
            X = DX+DELTAX*J
            CALL FREDGO(ARRAY,MN)
  150 CONTINUE
C
        X = DX+DELTAX*NY
        IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(KI,1))
        IF(KI.EQ.1) GO TO 180
        CALL PGDRAW(X,DELTAY*(KI+1))
        KI = KI+INCX
        IF(KI.GT.NX) GO TO 180
        IF(KI.EQ.1) GO TO 120
  160 DELTAX = -DELTAX
        DX = DELTAX*(1-KI-NY)
        X = DX+DELTAX
        CALL PGMOVE(X,DELTAY*(KI+1))
        CALL PGDRAW(X,ARRAY(KI,1))
        GO TO 130
C
C DO columns backwards because ended rows at far END of X
C
  170 KI = NX
        INCX = -1
        DX = DELTAX*(KI+NY)
        GO TO 160
C
C
  180 CALL PGEBUF
        END SUBROUTINE FREDDY


C***********************************************************************
C	Called by other routine in this file
C***********************************************************************

        SUBROUTINE FREDGO(ARRAY,MN)
        INTEGER MN
        REAL ARRAY(MN)
C
        INTEGER STEP,LEFT,RIGHT,IT,NX
        LOGICAL VISBLE
        REAL AL,AR,BL,EM,XX,X,Y,DELTAX
        COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
C Test visibility
C
        IF(ARRAY(IT).LT.0.0) GO TO 80
C
C This point is visible - was last?
C
        IF(VISBLE) GO TO 50
C
C No: calculate point where this line vanishes
C
   10 IF(LEFT.LE.NX .OR. MOD(LEFT-1,NX).EQ.0 .OR.
     &      RIGHT.LE.NX .OR. MOD(RIGHT-1,NX).EQ.0) GO TO 100
        AL = ABS(ARRAY(LEFT))
        AR = ABS(ARRAY(RIGHT))
        IF(ARRAY(LEFT).LT.0.0) GO TO 70
C               Right-hand point is crested
   20 RIGHT = RIGHT-STEP
        IF(ARRAY(RIGHT).LT.0.0) GO TO 20
C               Left-hand END of cresting line is either
C               RIGHT+NX or RIGHT-1
        LEFT = RIGHT+NX
        IF(ARRAY(LEFT).LT.0.0) LEFT = RIGHT-1
C
C               RIGHT and LEFT index into the endpoints of the
C               cresting line
   30 BL = ABS(ARRAY(LEFT))
        EM = ABS(ARRAY(RIGHT))-BL
        XX = EM-AR+AL
        IF(ABS(XX).LT.0.0001) GO TO 60
        XX = (AL-BL)/XX
   40 Y = EM*XX+BL
        IF(DELTAX.GT.0.0) XX = 1.0-XX
        XX = X-XX*DELTAX
        IF(VISBLE) GO TO 90
C               Drawing a line from an invisible point
C               to a visible one
        CALL PGMOVE(XX,Y)
        VISBLE = .TRUE.
   50 CALL PGDRAW(X,ARRAY(IT))
        RETURN
C
   60 XX = 0.5
        GO TO 40
C
C Left-hand point crested
C
   70 LEFT = LEFT-STEP
        IF(ARRAY(LEFT).LT.0.0) GO TO 70
C
C Right-hand END of cresting line is either LEFT+1 or LEFT-NX
C
        RIGHT = LEFT+1
        IF(ARRAY(RIGHT).LT.0.0) RIGHT = LEFT-NX
        GO TO 30
C
C This point is invisible; IF last one was too, THEN forget it;
C ELSE draw a line towards it
C
   80 IF(.NOT.VISBLE) RETURN
        GO TO 10
C
   90 CALL PGDRAW(XX,Y)
  100 VISBLE = .FALSE.
        RETURN
        END SUBROUTINE FREDGO


C***********************************************************************
C	Set array data plot to 3D surface type
C***********************************************************************

        SUBROUTINE WIRE(X_in,Y_in,Z_in,ZVEC_in,NROW,NCOL,NDROW,NDCOL,
     &                  ND2,CXLABL,CZLABL,CTITLE)

        IMPLICIT NONE

        CHARACTER*(*) CXLABL,CZLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,ND2,I,J,K
	SREAL :: X_in(NROW),Y_in(NCOL),Z_in(NDROW,NDCOL),
     &           ZVEC_in(ND2)
        REAL*4 X(NROW),Y(NCOL),Z(NDROW,NDCOL),ZVEC(ND2),DX,DZ,SIZE,
     &         ZMIN,ZMAX

	X=X_in; Y=Y_in; Z=Z_in; ZVEC=ZVEC_in

C       SET PLOT DIMENSIONS

        ZMIN=1E31
        ZMAX=-1E31
        DO 2 I=1,NROW
          DO 1 J=1,NCOL
            IF (Z(I,J).LT.ZMIN)ZMIN=Z(I,J)
            IF (Z(I,J).GT.ZMAX)ZMAX=Z(I,J)
   1      CONTINUE
   2    CONTINUE

        DZ=ZMAX-ZMIN
        DX=X(NROW)-X(1)
        K=0
        DO 4 J=1,NCOL
          DO 3 I=1,NROW
            K=K+1
            ZVEC(K)=Z(I,J)*DX/DZ
   3      CONTINUE
   4    CONTINUE


C       CALL PLOT ROUTINES

        SIZE=X(NROW)
        CALL PGENV(0., SIZE, 0., SIZE, 1, -2)
        CALL FREDDY(ZVEC,NROW,NROW,SIZE,25.0)

        RETURN
        END SUBROUTINE WIRE


C***********************************************************************

        SUBROUTINE PLOTCOL(X_in,Y_in,Z_in,NROW,NCOL,NDROW,NDCOL,
     &  		   CXLABL,CYLABL,CTITLE)

        IMPLICIT NONE

        CHARACTER*(*) CXLABL,CYLABL,CTITLE
        INTEGER NROW,NCOL,NDROW,NDCOL,I,J,IMAX,JMAX,CI
	SREAL :: X_in(NROW),Y_in(NCOL),Z_in(NDROW)
        REAL*4 X(NROW),Y(NCOL),Z(NDROW),ZMIN,ZMAX

	X=X_in; Y=Y_in; Z=Z_in

        CALL PGBBUF
        CALL PGQCI(CI)

C       FIND PLOT DIMENSIONS

        ZMIN=1E31
        ZMAX=-1E31
        DO 1 J=1,NCOL
          IF (Z(J).LT.ZMIN) ZMIN=Z(J)
          IF (Z(J).GT.ZMAX) ZMAX=Z(J)
   1    CONTINUE

C       CALL PLOT ROUTINES

        CALL PGENV(X(1),X(NROW),ZMIN,ZMAX,0,0)
        CALL PGLABEL(CXLABL,CYLABL,CTITLE)
        CALL PGSCI(1)
        CALL PGLINE(NROW,X,Z)
        CALL PGEBUF

        RETURN
        END SUBROUTINE PLOTCOL


C***********************************************************************

        SUBROUTINE DRAW(CTITLE,nDrawElt,nDrawRay,DrawEltVec,
     &  		DrawRayVec_in,XLabel,YLabel)
        use param_mod
        use elt_mod
	use traceutil_mod, only : 
     &        IsGrinRefElt,gradInd_kmax,gradIndRayPath,
     &        pgplotDrawPlane,iElt2iGradInd,IsOverDraw

        IMPLICIT NONE

        CHARACTER*(*) CTITLE
        CHARACTER*8 XLabel,YLabel
        INTEGER i,j,k,ii,nDrawRay,iEltOccur,iElt,iVec(2),
     &  	nDrawElt(mDrawRay),DrawEltVec(mDrawElt,mDrawRay),
     &	        iRayInc,nEltOccur(0:mElt),RayColor,EltColor,
     &          RefSrfColor,ObsSrfColor,FrameColor,midRay,
     &          GradIndRayColor
        Integer, save :: odRayColor=8,odGradIndRayColor=10
	SREAL :: DrawRayVec_in(2,mDrawElt,mDrawRay)
        REAL*4 DrawRayVec(2,mDrawElt,mDrawRay),XMAX,XMIN,YMAX,YMIN,
     &         x,y,xVec(mDrawRay,2),yVec(mDrawRay,2),DX,DY,
     &         xxVec(512),yyVec(512)
	Logical, save :: LfastDraw

	DrawRayVec = DrawRayVec_in

	! Draw color settings
        EltColor=3; RefSrfColor=4; ObsSrfColor=5
	If (.not. IsOverDraw) Then
          RayColor=2
	  GradIndRayColor=6
	Else
	  RayColor=odRayColor; odRayColor=odRayColor+1
	  GradIndRayColor=odGradIndRayColor  
	  odGradIndRayColor=odGradIndRayColor+1
	End If

C  Find plot dimensions

        XMIN=1E31; XMAX=-1E31
        YMIN=1E31; YMAX=-1E31

        DO j=1,nDrawRay
          DO i=1,nDrawElt(j)
            x=DrawRayVec(1,i,j)
            IF (x.LT.XMIN) XMIN=x
            IF (x.GT.XMAX) XMAX=x
            y=DrawRayVec(2,i,j)
            IF (y.LT.YMIN) YMIN=y
            IF (y.GT.YMAX) YMAX=y
          END DO
        END DO

        DX=0.05e0*(XMAX-XMIN)
        DY=0.05e0*(YMAX-YMIN)
        DX=MAX(DX,DY)
        XMIN=XMIN-DX
        XMAX=XMAX+DX
        YMIN=YMIN-DX
        YMAX=YMAX+DX

C  CALL plot routines

        if (.not. IsOverDraw) then
          CALL PGBBUF
          CALL PGQCI(FrameColor)
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,1,0)
	end if

        CALL PGLABEL(XLabel,YLabel,CTITLE)

C  Draw rays

        CALL PGSCI(RayColor)
                
        IF (nDrawRay.LT.99) THEN
	  ! draw one of every eight rays -jzlou
          iRayInc=(nDrawRay+1)/8 
        ELSE
          iRayInc=10
        END IF
        midRay=nDrawRay/2+1

	Write(*,*)' DRAW: nDrawRay,midRay,iRayInc =',nDrawRay,midRay,iRayInc

	LfastDraw=.true.
	Do j=1,nDrawRay
          Do i=1,nDrawElt(j)
	    If (IsGrinRefElt(i)) Then
	      ! Has graded-index lens, must draw rays segment to segment
	      LfastDraw=.false.  
	      Go to 105
	    End If
	  End Do 	  
	End Do
 105    Continue

#if 0
	print*,'** nDrawElt(33),midRay =',nDrawElt(33),midRay
	DO j=1,midRay,iRayInc
	  print*,'** nDrawElt(j) =',nDrawElt(j)
	END DO
#endif

C DrawRayVec(1:2,i,j): position coords of j-th ray hitting i-th element

	IF (LfastDraw) THEN
	  ! Efficient way, when no graded index lens in system
          DO j=1,midRay,iRayInc
	    !print*,'** nDrawElt(j) =',nDrawElt(j)
            DO i=1,nDrawElt(j)
              xVec(i,1)=DrawRayVec(1,i,j)
              yVec(i,1)=DrawRayVec(2,i,j)
    	      if (.false.) then
                print*,'** Draw: i =',i,',  xVec,yVec=',
     &            xVec(i,1),yVec(i,1)
	      end if
            END DO
            CALL PGLINE(nDrawElt(j),xVec(1,1),yVec(1,1))
          END DO
	ELSE
	  ! Draw one ray segment (between two surfaces) at a time
	  DO j=1,midRay,iRayInc
	    xVec(1,1)=DrawRayVec(1,1,j)
            yVec(1,1)=DrawRayVec(2,1,j)
	    DO i=2,nDrawElt(j)
!	      print*,'** i,nDrawElt(j),IsGrinRefElt(i-1) =',
!    &          i,nDrawElt(j),IsGrinRefElt(i-1)
	      xVec(i,1)=DrawRayVec(1,i,j)
              yVec(i,1)=DrawRayVec(2,i,j)
	      If (.not.IsGrinRefElt(i-1)) Then
	        ! draw straight line segment to next surface
	        CALL PGLINE(2,xVec(i-1,1),yVec(i-1,1)) 
              Else
	        ! Is front surface of graded-index lens, draw curved ray path within
	        ! the lens instead.  
	        ii=iElt2iGradInd(i-1)
	        !print*,'**1: j,nDrawElt(j) =',j,nDrawElt(j)
	        !print*,'gradInd_kmax(j,ii) =',gradInd_kmax(j,ii)
	        If (pgplotDrawPlane==1) Then ! 'YZ'
	          xxVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(3,1:gradInd_kmax(j,ii),j,ii)
	          yyVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(2,1:gradInd_kmax(j,ii),j,ii)
	        Else If (pgplotDrawPlane==3) Then ! 'XY'
	          xxVec(1:gradInd_kmax(j,ii)) 
     &              =gradIndRayPath(1,1:gradInd_kmax(j,ii),j,ii)
                  yyVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(2,1:gradInd_kmax(j,ii),j,ii)
	        Else ! default to 'XZ'
	          xxVec(1:gradInd_kmax(j,ii)) 
     &              =gradIndRayPath(3,1:gradInd_kmax(j,ii),j,ii)
                  yyVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(1,1:gradInd_kmax(j,ii),j,ii)
	        End If
	        CALL PGSCI(GradIndRayColor)
	        Do k=1,gradInd_kmax(j,ii)-1,3
	          CALL PGLINE(2,xxVec(k:k+1),yyVec(k:k+1)) 
	        End Do
	        CALL PGSCI(RayColor)
	      End If
	    END DO
	  END DO
	END IF

	!print*,'DrawRayVec(1,1:2,1) =',DrawRayVec(1,1:2,1)
	!print*,'DrawRayVec(2,1:2,1) =',DrawRayVec(2,1:2,1)

	IF (LfastDraw) THEN
          DO j=nDrawRay,midRay,-iRayInc
            DO i=1,nDrawElt(j)
              xVec(i,1)=DrawRayVec(1,i,j)
              yVec(i,1)=DrawRayVec(2,i,j)
            END DO
            CALL PGLINE(nDrawElt(j),xVec(1,1),yVec(1,1))
          END DO
	ELSE
	  ! Draw one ray segment (between two surfaces) at a time
	  DO j=nDrawRay,midRay,-iRayInc
            xVec(1,1)=DrawRayVec(1,1,j)
            yVec(1,1)=DrawRayVec(2,1,j)
            DO i=2,nDrawElt(j)
	      xVec(i,1)=DrawRayVec(1,i,j)
              yVec(i,1)=DrawRayVec(2,i,j)
              If (.not.IsGrinRefElt(i-1)) Then
                ! draw straight line segment to next surface
                CALL PGLINE(2,xVec(i-1,1),yVec(i-1,1))
              Else
                ! Is front surface of graded-index lens, draw curved ray path within
                ! the lens instead.
	        ii=iElt2iGradInd(i-1)
	        !print*,'**2: j,nDrawElt(j) =',j,nDrawElt(j)
                !print*,'gradInd_kmax(j,ii) =',gradInd_kmax(j,ii)
	        If (pgplotDrawPlane==1) Then ! 'YZ'
                  xxVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(3,1:gradInd_kmax(j,ii),j,ii)
                  yyVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(2,1:gradInd_kmax(j,ii),j,ii)
                Else If (pgplotDrawPlane==3) Then ! 'XY'
                  xxVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(1,1:gradInd_kmax(j,ii),j,ii)
                  yyVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(2,1:gradInd_kmax(j,ii),j,ii)
                Else ! default to 'XZ'
                  xxVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(3,1:gradInd_kmax(j,ii),j,ii)
                  yyVec(1:gradInd_kmax(j,ii))
     &              =gradIndRayPath(1,1:gradInd_kmax(j,ii),j,ii)
                End If
	        CALL PGSCI(GradIndRayColor)
                Do k=1,gradInd_kmax(j,ii)-1,3
                  CALL PGLINE(2,xxVec(k:k+1),yyVec(k:k+1))
                End Do
	        CALL PGSCI(RayColor)
              End If
	    END DO
	  END DO
	END IF


C  Draw surfaces
C  DrawEltVec = [0 1 2 3 3 2 1 6 7 8]

C nDrawRay: total rays to draw
C nDrawElt(j): draw thru nDrawElt(j) elements for j-th ray
C DrawRayVec(1:2,i,j): position coords of j-th ray hitting i-th element 
C DrawEltVec(i,j): given j (j-th ray) and i (i-th element j-th ray will
C                  hit), DrawEltVec(i,j) gives a valid element id.
C -jzlou

        DO iElt=0,mElt
          iVec(1:2)=0
          nEltOccur(iElt)=0

          DO j=1,nDrawRay
            iEltOccur=0
            DO i=1,nDrawElt(j)
              IF (DrawEltVec(i,j).EQ.iElt) THEN
                iEltOccur=iEltOccur+1
                IF (iEltOccur.LE.2) THEN
                  iVec(iEltOccur)=iVec(iEltOccur)+1
                  xVec(iVec(iEltOccur),iEltOccur)=DrawRayVec(1,i,j)
                  yVec(iVec(iEltOccur),iEltOccur)=DrawRayVec(2,i,j)
                  IF (nEltOccur(iElt).LT.iEltOccur)
     &              nEltOccur(iElt)=iEltOccur
                END IF
              END IF
            END DO
          END DO

          DO iEltOccur=1,nEltOccur(iElt)
            IF ((EltID(iElt).EQ.3).OR.(EltID(iElt).EQ.8)) THEN
              CALL PGSCI(RefSrfColor)  ! use ref surf color
            ELSE IF (EltID(iElt).EQ.7) THEN
              CALL PGSCI(ObsSrfColor)  ! use obscuration color
            ELSE
              CALL PGSCI(EltColor)     ! use element color
            END IF
            CALL PGLINE(iVec(iEltOccur),xVec(1,iEltOccur),
     &      	        yVec(1,iEltOccur))
          END DO
        END DO

	if (.not. IsOverDraw) then
          CALL PGSCI(FrameColor)
          CALL PGEBUF
	end if

        RETURN
        END SUBROUTINE DRAW

C***********************************************************************

	End Module pgplot_mod
