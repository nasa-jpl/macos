C
C smacos_glass.inc
C

C  ***Load glass tables***

#ifdef SUNOS
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
            glassFile='macos_glass_list.txt'
            INQUIRE (FILE=glassFile,EXIST=EXIST)
          END IF
#else
	  ! Don't GETENV MacosHome on non-UNIX machines,
	  ! assume glass_list file in local folder -jzlou
          !glassFile=macosHome(1:ICLEN(macosHome))//'macos.glass'
          glassFile='macos_glass_list.txt'
          INQUIRE (FILE=glassFile,EXIST=EXIST)    
#endif

          GlassName(1)='Air'
          CALL dZERO(GlassTable,(6*mGlass))

          IF (EXIST) THEN
            OPEN (UNIT=32,FILE=glassFile,STATUS='OLD',FORM='FORMATTED')
            k=0
            DO i=2,mGlass
              READ (32,*,ERR=476,END=476) 
     &		glassName(i),(glassTable(j,i),j=1,6)
              k=k+1
            END DO
 476        CONTINUE
	    ! Write out how many glasses read in from the glass_list file
            WRITE(*,312)glassFile(1:ICLEN(glassFile)),k
          ELSE
	    ! When file 'macos_glass_list.txt' doest not exist, define
	    ! following glasses only
            WRITE(*,314)
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

#include "prop_defs.inc"

	  ifInit=.TRUE.
          ifLoad=.FALSE.
          ifTrace=.FALSE.
          ifBuild=.FALSE.
          ifPropagate=.FALSE.
          ifCmt=.FALSE.
          ifPerturb=.FALSE.
          ifPlot=.FALSE.
          ifPol=.FALSE.
          ifVecDif3=.FALSE.
          ifRegrid(1:mElt)=.FALSE.
          ifOld=.FALSE.
          ifRayPrint=.FALSE.
          ifGrid=.FALSE.
          ifPixLoc=.FALSE.
          ifAdd=.FALSE.
          ifPixArray=.FALSE.
          ifDMBuild=.FALSE.
          ifCentroid=.TRUE.
          ifFilterData=.FALSE.
          ifLNegOK=.TRUE.
          ifUDSinit=.FALSE.
          ifStopSet=.FALSE.
          i=mElt
          CALL LZERO(ifIntSrfSet,i)
          iCurRayElt=0
          iCurWFElt=0
          CALL ArchInit(iPloTyp,ifGraphics)
          iObsOpt=1
          BeamType=1
          StretchType=1
          spcOption=1
          sizPix=0d0
          EltPix(1)=0d0
          EltPix(2)=0d0
          WinCenPix(1)=0d0
          WinCenPix(2)=0d0
          CALL RTIME(StartTime)
          i=ABS(MOD(StartTime,1000.0))
          NoiseSeed(1)=29*i
          NoiseSeed(2)=32000-NoiseSeed(1)
