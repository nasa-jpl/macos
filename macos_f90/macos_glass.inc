C  ***Load glass tables***
C  Resolution order (via find_macos_file in param_mod):
C    1) CWD
C    2) directory of the running macos executable
C    3) $MACOS_HOME
C    4) compile-time build directory (BUILD_LOC)
C  If none of those have it, fall back to the hardcoded
C  Air / BK7 / LAK9 set below.

        glassFile = 'macos_glass_list.txt'
        CALL find_macos_file('MACOS_HOME', glassFile)
        IF (ICLEN(glassFile) .GT. 0) THEN
          INQUIRE (FILE=glassFile, EXIST=EXIST)
        ELSE
          EXIST = .FALSE.
        END IF

	!print*,'**INQUIRE: EXIST =',EXIST
	!STOP ! test only

        GlassName(1)='Air'
        CALL DZERO(GlassTable(1,1),6)
        IF (EXIST) THEN
          OPEN (UNIT=32,FILE=glassFile,STATUS='OLD',FORM='FORMATTED')
          k=0
          DO i=2,mGlass
            READ(32,*,END=476) glassName(i),(glassTable(j,i),j=1,6)
            k=k+1
          END DO
 476      CONTINUE
          ! write out how many glasses read in from the glass_list file
          WRITE(*,312) glassFile(1:ICLEN(glassFile)),k
	  WRITE(*,*)' ' 
        ELSE
	  ! When file 'macos_glass_list.txt' doest not exist, define
          ! following glasses only
          WRITE(*,314) ! write out 'No glass table were found'
	  WRITE(*,*)' '
          GlassName(2)='BK7'
          GlassTable(1,2)=1.03961212
          GlassTable(2,2)=2.31792344e-1
          GlassTable(3,2)=1.01046945
          GlassTable(4,2)=6.00069867e-3
          GlassTable(5,2)=2.00179144e-2
          GlassTable(6,2)=1.0350653e+2
          GlassName(3)='LAK9'
          GlassTable(1,3)=1.46231905d0
          GlassTable(2,3)=3.44399589d-1
          GlassTable(3,3)=1.15508372d0
          GlassTable(4,3)=7.24270156d-3
          GlassTable(5,3)=2.43353131d-2
          GlassTable(6,3)=8.54686868d1
	  print*,' **Glass table file macos_glass_list.txt not found;'
          print*,
     &  '   only Air, and BK7 and LAK9 glass types have been defined'
        END IF


