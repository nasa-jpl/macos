C
C File save_nom.inc  
C to be included in macos_cmd_loop.inc
C

#if 1
	! Save Rx state
	ChfRayPosRx(1:3)=ChfRayPos(1:3)
        ChfRayDirRx(1:3)=ChfRayDir(1:3)
        xGridRx(1:3)=xGrid(1:3)
        yGridRx(1:3)=yGrid(1:3)
        zGridRx(1:3)=zGrid(1:3)
        WavelenRx=Wavelen

        Do id=1,nElt
          KrRx(id)=KrElt(id)
          KcRx(id)=KcElt(id)
          nObsRx(id)=nObs(id)
          ObsTypeRx(1,id)=ObsType(1,id)
          psiRx(1:3,id)=psiElt(1:3,id)
          vptRx(1:3,id)=vptElt(1:3,id)
          rptRx(1:3,id)=rptElt(1:3,id)
          pMonRx(1:3,id)=pMon(1:3,id)
          xMonRx(1:3,id)=xMon(1:3,id)
          yMonRx(1:3,id)=yMon(1:3,id)
          zMonRx(1:3,id)=zMon(1:3,id)
          ObsVecRx(1:3,1,id)=ObsVec(1:3,1,id)
          ZernCoefRx(1:45,id)=ZernCoef(1:45,id)
          CALL DEQUATE(TEltRx(1,1,id),TElt(1,1,id),36)
        End Do
        CALL DEQUATE(ToutRx,Tout,49)
#endif

        ! Now compute and save nominal OPD (in OPDMatRx)
	iEndElt=nElt-1
	CALL OPD(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &        ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &        iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &        nBadRays,StartIndRef,StartExtinc,OPDMatRx,WFE,DAvgL,
     &        L1,Cmatrix,DrawRayVec,DrawEltVec)

	!print*,'*** Rx RMS WFE =',WFE

