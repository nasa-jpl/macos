C  ***For Re-Load glass tables***

        CALL GETENV('MACOS_HOME',macosHome)
        IF (.false. .and. macosHome(1:2).EQ.'  ') THEN
          WRITE(*,313)
          macosHome='~/macosData/'
          glassFile='~/macosData/macos.glass'
          INQUIRE (FILE=glassFile,EXIST=EXIST)
          IF (.NOT.EXIST) THEN
            macosHome='/home/comp/macosData/'
            glassFile='/home/comp/macosData/macos.glass'
            INQUIRE (FILE=glassFile,EXIST=EXIST)
            IF (.NOT.EXIST) THEN
              macosHome=''
              glassFile='macos.glass'
              INQUIRE (FILE=glassFile,EXIST=EXIST)
            END IF
          END IF
        ELSE
          !glassFile=macosHome(1:ICLEN(macosHome))//'macos.glass'
          glassFile='macos.glass'
          INQUIRE (FILE=glassFile,EXIST=EXIST)
        END IF

        GlassName(1)='Air'
        CALL DZERO(GlassTable(1,1),6)
        IF (EXIST) THEN
          OPEN (UNIT=32,FILE=glassFile,STATUS='OLD',FORM='FORMATTED')
          k=0
          DO i=2,mGlass
            READ (32,*,END=477) glassName(i),(glassTable(j,i),j=1,6)
            k=k+1
          END DO
 477      CONTINUE
          WRITE(*,312)glassFile(1:ICLEN(glassFile)),k
        ELSE
!--> Built-in catalog fallback (PLAN section 0) -- generated
!    glass_builtin.f90; a found glass file still overrides.
          CALL LoadBuiltinGlass(GlassName, GlassTable, mGlass, k)
          WRITE(*,*)' No glass table file found;',
     &              ' using built-in catalog:',k,' glasses.'
!<--
        END IF


