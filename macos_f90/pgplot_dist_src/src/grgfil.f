C*GRGFIL -- find data file
C+
      SUBROUTINE GRGFIL(TYPE, NAME)
      CHARACTER*(*) TYPE, NAME
C
C This routine encsapsulates the algorithm for finding the PGPLOT
C run-time data files.
C
C 1. The binary font file: try the following in order:
C     file specified by PGPLOT_FONT
C     file "grfont.dat" in directory specified by PGPLOT_DIR
C                       (with or without '/' appended)
C     file "grfont.dat" in directory /usr/local/pgplot/
C
C 2. The color-name database: try the following in order:
C     file specified by PGPLOT_RGB
C     file "rgb.txt" in directory specified by PGPLOT_DIR
C                       (with or without '/' appended)
C     file "rgb.txt" in directory /usr/local/pgplot/
C
C Arguments:
C  TYPE (input)  : either 'FONT' or 'RGB' to request the corresponding
C                  file.
C  NAME (output) : receives the file name.
C--
C  2-Dec-1994 - new routine [TJP].
C-----------------------------------------------------------------------
      !CHARACTER*(*) DEFDIR, DEFFNT, DEFRGB
      CHARACTER*255 DEFDIR, DEFFNT, DEFRGB
      !PARAMETER  (DEFDIR='/usr/local/pgplot/')

!#ifdef LNX32
!     PARAMETER (DEFDIR='/proj/jwst/pgplot/lnx/32/')
!#endif

!#ifdef LNX64
!     PARAMETER 
!    &(DEFDIR='/proj/jwst/jzlou/tst_malfoy/MACOS_BUILD/pgplot_lnx64/')
!#endif

      PARAMETER (DEFFNT='grfont.dat') 
!    &(DEFFNT='/mnt/BATBRAIN/USERS/jzlou/MACOS/pgplot_lnx64/grfont.dat')
      PARAMETER  (DEFRGB='rgb.txt')
      CHARACTER*255 FF
      CHARACTER*16 DEFLT
      INTEGER I, L, LD, PGSTRLEN
      LOGICAL TEST, DEBUG
C
C Is debug output requested?
C
      CALL GRGENV('DEBUG', FF, L)
      DEBUG = L.GT.0
C
C Which file?
C
      IF (TYPE.EQ.'FONT') THEN
         DEFLT = DEFFNT
         LD = LEN(DEFFNT)
      ELSE IF (TYPE.EQ.'RGB') THEN
         DEFLT = DEFRGB
         LD = LEN(DEFRGB)
      ELSE
         CALL GRWARN('Internal error in routine GRGFIL')
      END IF
C
C Try each possibility in turn.
C
      DO 10 I=1,4
         IF (I.EQ.1) THEN
            CALL GRGENV(TYPE, FF, L)
         ELSE IF (I.EQ.2) THEN
            CALL GRGENV('DIR', FF, L)
            IF (L.GT.0) THEN
               FF(L+1:) = DEFLT
               L = L+LD
            END IF
         ELSE IF (I.EQ.3) THEN
            CALL GRGENV('DIR', FF, L)
            IF (L.GT.0) THEN
               FF(L+1:L+1) = '/'
               FF(L+2:) = DEFLT
               L = L+1+LD
            END IF
         ELSE IF (I.EQ.4) THEN
	    LD = PGSTRLEN(DEFLT,16)
	    L = PGSTRLEN(DEFDIR,255)
            FF = DEFDIR(1:L)//DEFLT
            L = L+LD
         END IF

         IF (L.GT.0) THEN
            IF (DEBUG) THEN
               CALL GRWARN('Looking for '//FF(:L))
            END IF
            INQUIRE (FILE=FF(:L), EXIST=TEST)
            IF (TEST) THEN
               NAME = FF(:L)
               RETURN
            ELSE IF (DEBUG) THEN
               CALL GRWARN('WARNING: file not found')
            END IF
         END IF
 10   CONTINUE

      
      print*, '***grgfil: FF =',FF(1:L)
      print*, '***grgfil: NAME =',NAME
C
C Failed to find the file.
C
      NAME = DEFLT
C-----------------------------------------------------------------------
      END

C-----------------------------------------------------------------------

      FUNCTION PGSTRLEN(STR,N)
      Integer CLEN,PGSTRLEN,I
      Character(Len=*) :: STR

      CLEN=0
      DO I=1,N
	If (STR(I:I) /= ' ') THEN
          CLEN=CLEN+1
	Else
	  PGSTRLEN=CLEN
	  RETURN
	End If 
      END DO
      END FUNCTION PGSTRLEN 

