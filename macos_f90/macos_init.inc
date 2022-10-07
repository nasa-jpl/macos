C  *** Initialize stuff ***
#include "prop_defs.inc"

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

        bat_task=1; jou_task=1
c sws
        pstack=1
        read_len(pstack)=0
c start of file unit numbers (could be n+10)
        ustack= 50

	nPgPanel=1
