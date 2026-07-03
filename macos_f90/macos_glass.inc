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
!--> Built-in catalog fallback (PLAN section 0: no separate runtime
!    glass file needed).  glass_builtin.f90 is GENERATED from
!    macos_glass_list.txt by tools/gen_glass_builtin.py; a glass
!    file found at runtime still overrides (branch above).
          CALL LoadBuiltinGlass(GlassName, GlassTable, mGlass, k)
          WRITE(*,*)' No glass table file found;',
     &              ' using built-in catalog:',k,' glasses.'
	  WRITE(*,*)' '
!<--
        END IF


