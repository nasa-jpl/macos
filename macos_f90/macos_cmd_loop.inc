C***********************************************************************
C       Begin file macos_cmd_loop_.inc
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C
C  MACOS top-level command loop
C

C  *** Begin MACOS top-level command loop ***

 1      CONTINUE
        command=' '


#ifdef MACOS_CMD

	  CALL PROMPT(command,'MACOS> ')

	  IF (LCMP(command,'QUIT',1) .OR.
     &        LCMP(command,'EXIT',4) .OR.
     &        LCMP(command,'BYE',3)) THEN
	    CALL GRAEND
	    CLOSE(1)
	    CLOSE(2)
	    CLOSE(3)
c***loh
	    IF (jou_task .eq. 2)  THEN
	      CLOSE (unit=27,STATUS='KEEP')
	    END IF
c***loh
#ifdef VIS3D_DEMO
            Call killdemo
#endif
	    STOP

#else
! for SMACOS
	  CALL PROMPT(command,'S-MACOS>')

          IF (LCMP(command,'RETURN',4)) THEN
            RETURN

          ELSE IF (LCMP(command,'QUIT',1)) THEN
            CALL GRAEND
            CLOSE(1)
            CLOSE(2)
            CLOSE(3)
c***loh
            IF (jou_task .eq. 2)  THEN
              CLOSE (unit=27,STATUS='KEEP')
            END IF
c***loh
            RETURN
#endif


C         *** print help message***

	  ELSE IF (LCMP(command,'HELP',3)) THEN
	    WRITE(*,300)
	    WRITE(*,307)
	    WRITE(*,304)
	    CALL IACCEPT_S(i,0,1,' Hit return to continue...')
	    WRITE(*,305)
	    WRITE(*,306)
	    CALL IACCEPT_S(i,0,1,' Hit return to continue...')
	    WRITE(*,321)


! For command-line MACOS only
#ifdef MACOS_CMD
C         *** reset model size ***

	  ELSE IF (LCMP(command,'MRESET',4)) THEN

 123	    CALL IACCEPT_S(macos_model_size,128,1,
     &       'Enter new size '//
     &       '(128, 256, 512, 1024, 2048, 4096, 8192):')

	    print*,'macos_model_size = ',
     &        macos_model_size

	    if ((macos_model_size == 128) .or.
     &	        (macos_model_size == 256) .or.
     &          (macos_model_size == 512) .or.
     &		(macos_model_size == 1024) .or.
     &          (macos_model_size == 2048) .or.
     &          (macos_model_size == 4096) .or.
     &          (macos_model_size == 8192)) then

                 write(*,*) 'Resetting MACOS model size...'

	         ! Reinitialize MACOS dynamic arrays
                 call macos_init_all(macos_model_size)
                 call dopt_init
                 macos_init=.true.

	         ! Local dynamic arrays
	         if (allocated(L1)) then
		   deallocate(L1,DrawEltVec,
     &		     R1,R2,PertVec,DrawRayVec,
     &		     D2,DV1,DV2,CD1,CD2,DWF,stat=m_err_macos)
		   if (m_err_macos /= 0)
     &               call macos_memory_failure(
     &                 'macos_mod_init: deallocate failed!')
		 end if
	         !
        	 allocate(L1(md2), DrawEltVec(mDrawElt,mDrawRay),
     &             R1(mdttl,mdttl), R2(mdttl,mdttl),
     &             PertVec(6,mElt),
     &             DrawRayVec(2,mDrawElt,mDrawRay),
     &             D2(mdttl,mdttl),
     &             DV1(md2), DV2(md2),
     &             CRWalk(3,mElt), CRWalkProjLoc(3,mElt),
     &             CD1(mdttl,mdttl), CD2(mdttl,mdttl),
     &             DWF(mdttl,mdttl), stat=m_err_macos)
                 if (m_err_macos /= 0) then
                   call macos_memory_failure('MACOS: allocate failed!')
	         else
          	   L1 = .false.; DrawEltVec = 0
          	   R1 = 0.0; R2 = 0.0; PertVec = 0.0
          	   DrawRayVec = 0.0; D2 = dzero_
          	   DV1 = dzero_; DV2 = dzero_
                   CRWalk = dzero_; CRWalkProjLoc = dzero_
          	   CD1 = (0d0,0d0); CD2 = (0d0,0d0)
                   DWF = czero_
                 end if

		 WRITE(*,303) mElt, mEm1, mRay, bRay, mWF,
     &               mIntSrf, mDP, mPix, mPix, mdttl, mdttl

C Reload glass table
#include "rl_macos_glass.inc"

C  Variable re-initializations (including ifLoad)
                 ifPlot_sv = ifPlot  ! save old state for some variables ...
#include "macos_init.inc"
	      ifPlot = ifPlot_sv  ! restore old state

	    else
	      write(*,*) '**Error: Invalid model size: ',
     &          macos_model_size
	      goto 123
	    end if
#endif
	  ELSE IF (LCMP(command,'JWST_V3D',4) .OR.
     &             LCMP(command,'Vis3D',4)) THEN

#ifdef VIS3D_DEMO
	    !Call startpthread
	    Call startchild
#else
	    WRITE(*,*)'  **JWST 3D rendering demo not included in MACOS!'
#endif


#ifdef MACOS_CMD
C 	Select PGPLOT Panels among 1, 2 (1x2), 3 (1x3) and 4 (2x2)
	  ELSE IF (LCMP(command,'PGP',3)) THEN
	    CALL IACCEPT_S(nPgPanel,1,1,
     &             'Enter PGPLOT display panel option (1,2,3 or 4):')
	    If (nPgPanel/=1.and.nPgPanel/=2.and.
     &          nPgPanel/=3.and.nPgPanel/=4) Then
	      nPgPanel=1
	      WRITE(*,*)'  Invalid input! Use 1 PGPLOT panel'
	    Else
	      if (nPgPanel==1) then
	        WRITE(*,*)'  Use 1 PGPLOT panel'
	      else if (nPgPanel==2) then
	        WRITE(*,*)'  Use 1 x 2 PGPLOT panels'
	      else if (nPgPanel==3) then
	        WRITE(*,*)'  Use 1 x 3 PGPLOT panels'
	      else if (nPgPanel==4) then
	        WRITE(*,*)'  Use 2 x 2 PGPLOT panels'
	      end if
	    End If
	    ifPlot=.FALSE.  ! enable reinitialize PGPLOT

          ELSE IF (LCMP(command,'CIR',3)) THEN
            ifPGColor=.TRUE.
            WRITE(*,*)' Color image rendering is on'
          ELSE IF (LCMP(command,'GIR',3)) THEN
            ifPGColor=.FALSE.
            WRITE(*,*)' Gray image rendering is on'
#endif


#ifdef MACOS_CMD
C         *** A few OS shell utility commands ***

	  ELSE IF (LCMP(command,'!',1)) THEN
	    Call SHELLCMD(cmdline,130)
	    read_cur(pstack)=0; read_len(pstack)=0 ! clean up input deck

	  ELSE IF (LCMP(command,'PWD',3)) THEN
            Call PWD

	  ELSE IF (LCMP(command,'CD',2)) THEN
	    CALL CACCEPT(ANS,'no default', 'Enter folder path ')
	    Call ChgDir(ANS,ICLEN(ANS))

	  ELSE IF (LCMP(command,'RX',2) .OR.
     &             LCMP(command,'LS',2) .OR. LCMP(command,'DIr',2)
     &             .OR. LCMP(command,'LL',2)) THEN
	    IF (read_len(pstack) .NE. 0) THEN
	      cmdArg=var_store(read_cur(pstack),pstack)
	      read_cur(pstack)=0; read_len(pstack)=0
	    ELSE
              cmdArg=' '
	    END IF

	    If (LCMP(command,'RX',2)) Then
	      Call LSFILES(0,cmdArg,0)
	      RX_Flg=.TRUE.
	      WRITE(*,*)' '
	      WRITE(*,*)' The number next to a Rx filename is file_id.'
	      WRITE(*,*)' Use FID command to load a Rx file.'//
     &                  ' Syntax: FID file_id'; WRITE(*,*)' '
	    Else If (LCMP(command,'LL',2)) Then
              Call LSFILES(2,cmdArg,ICLEN(cmdArg))
            Else
              Call LSFILES(1,cmdArg,ICLEN(cmdArg))
            End If

	  ELSE IF (LCMP(command,'VI',2).OR.LCMP(command,'EMAcs',3)) THEN
	    CALL CACCEPT(ANS,'no default', 'Enter filename ')
            Call StartEditor(command,ICLEN(command),ANS,ICLEN(ANS))
	    WRITE(*,*)' '
#endif


C         ***next command: Load optical description data***

	  ELSE IF (LCMP(command,'LOAD',3).OR.
     &             LCMP(command,'OLD',2).OR.LCMP(command,'NEW',2)
     &             .OR. LCMP(command,'FID',3)) THEN

	    IF (LCMP(command,'NEW',2)) THEN
	      ifOld=.FALSE.
	    ELSE
	      ifOld=.TRUE.
	    END IF

#ifdef MACOS_CMD
	    IF (LCMP(command,'FID',3)) THEN
	      ifFID=.TRUE.
	    ELSE
	      ifFID=.FALSE.
	    END IF
#endif

	    nOptChfRayDir=0; nOptChfRayPos=0
	    nPrvGlass=0; BaseUnits_FLG=.FALSE.
	    zcoType=ZernCircular

#ifdef DESIGN_OPTIM
	    ! Design optimization related
	    isVarElt=.false.; varEltDOF=0
            varElts=0; varDOFs=0
            psiElt_old=0d0; vptElt_old=0d0
	    nVarElt=0
            nOptFov=0
	    tgtWF_flg=.false.
	    fov_weight=1.0  ! uniform weighting by default
	    nOptWavelen=0

	    ! Zernike and/or AsphCoef optimization
	    nOptEltZern=0; OptEltZernTerm=-1
	    nOptEltAsph=0; OptEltAsphTerm=-1

	    isSrcOpt=.FALSE.   ! source optimization, default to false
	    LOptIfFEX=.TRUE.   ! default
	    OptTgtElt=-1  ! init value
	    OptAlg=NonLin ! init value
	    LOptNomSens=.FALSE. ! init value
	    LOptCons=.FALSE.  ! default to non-constrained sytem optimization
	    LUseChfRayIfOK=.TRUE.

	    beamOpt%ifOptBeamDir=.FALSE.; beamOpt%ifOptBeamPos=.FALSE.
	    beamOpt%ifOptBeamSize=.FALSE.
	    beamOpt%ifOptBeamRefRayDir=.FALSE.

#endif
	    LUpdateTElt_FLG=.FALSE.
            saveOPD_Flg=.FALSE.  ! default
            SaveRayPosHist_Flg=.FALSE. ! default
            DumpRayPosHist_Flg=.FALSE. ! default
            SaveVis3d_dat_Flg=.FALSE.  ! default
	    EltGrp(0,0:mElt)=0
	    RxStopSet=.FALSE.; ifStopSet=.FALSE.
            RxNoStopSet_Flg=.FALSE.
	    IsVarIndRefElt=.FALSE.; IsGrinRefElt=.FALSE.
	    TotalLensletsArray=0
	    Poly3DApVecL=.FALSE.; ifSave3DApVec=.FALSE.
	    zernUsrOpt=.FALSE.; LxObsDef=.FALSE.;
	    UDBeam_FLG=.FALSE.; ifObjImgSet=.FALSE.
	    ifRegrid(1:mElt)=.FALSE.; ifOPDModGrid=.FALSE.
	    LRayTgtElt=.FALSE. ! default
            refRay = -1d2   ! default initial value, indicating state RRAY before running
	                    ! RRAY command  01/27/2017
            GridSrfOrder=1  ! default GridSrf interpolation to linear order
	    iElt2iGradInd=0 ! for graded index lenses, mapping from iElt to
	                    ! graded index lens id, initialized in msmacosio.inc
	    numGradInd=0; nDOE=0
	    iEltToMetSrf=0; tMetSrf=0; metBeamFlg=0 ! clean up metrology setting
            OPDRefRayLen_FLG=.FALSE.; ShowMetData_FLG=.FALSE.
!	    Do i=1,mElt
!	      LensletStrArray(1:mLenslet,i)%IndRef=IndRef(i)
!	    End Do

	    ! For source local frames (will set to true if defined in Rx)
            SrcLF_FLG=.FALSE.  ! default
            CalcSrcZFrame_FLG=.TRUE.  ! default

#ifdef MACOS_CMD
	    PrtPrinRots=.FALSE.; LudDrawGrid=.FALSE.
	    CALL MBFile6(npts,ifOld,ifFID,LRtn_FLG)
	    if (.not.LRtn_FLG) GO TO 1
#else
	    ! SMACOS
	    CALL MBFile6(npts,ifOld,LRtn_FLG)
#endif

#ifdef DESIGN_OPTIM
	    if (OptTgtElt < 0) then
	      OptTgtElt=nElt-1      ! default
	    end if
#endif

	    ifLoad=.TRUE.; ifBuild=.FALSE.
	    ifPropagate=.FALSE.; ifPerturb=.FALSE.
	    ifPol=.FALSE.; ifPixLoc=.FALSE.
	    ifVecDif3=.FALSE.; ifTrace=.FALSE.
	    ifGrid=.FALSE.; ifAdd=.FALSE.;
	    ifAddCA=.FALSE.; ifPixArray=.FALSE.
	    ifCmt=.FALSE.; ifEFLCalc=.FALSE.
	    ifCentroid=.TRUE.; ifFilterData=.FALSE.
	    ifDMBuild=.FALSE.; ifLNegOK=.TRUE.; i=mElt
	    CALL LZERO(ifIntSrfSet,i)
	    iCurRayElt=0; iCurWFElt=0; iObsOpt=1
	    If (UDBeam_FLG) Then
	      BeamType=5
	    Else
              BeamType=1
	    End If
	    nGridPts=npts+1; sizPix=0d0
	    nTestPts=MIN(nGridPts,21)
	    nMinPts=3; nSpotPts=nGridPts
	    EltPix(1:2)=0d0; WinCenPix(1:2)=0d0
	    BWK_FLG=.FALSE.
	    LCalcWFELT=.FALSE.  ! default
#ifdef MACOS_CMD
	    curFov=1; curwl=1; OPDSave=0d0
            firstOPDCalc=.FALSE.; multiOPDCalc=.FALSE.
	    firstIntCalc=.FALSE.
#endif

#if 0
	    Do iElt=1,nElt
	      If (IsVarIndRefElt(iElt)) Then
	        Call InitLensletIndex(LensArrayType(iElt),
     &                          TotalLensletsArray(iElt),iElt)
	      End If
	    End Do
#endif

#ifdef DESIGN_OPTIM
	    if (nOptFov==0) then
	      ! Rx FOV, by default, is used for design optimization
	      opt_fov(1:3,1,1)=ChfRayDir(1:3)
              opt_fov(1:3,2,1)=ChfRayPos(1:3)
	      nOptFov=1
	    end if
	    if (nOptWavelen==0) then
	      ! Rx Wavelength, by default, is used for design optimization
	      opt_wavelen(1)=Wavelen
	      nOptWavelen=1
	    end if

	    OptSrcRpt_FLG=.FALSE.

            ! ** Support Linear WFC, calc system sensitivity around nominal
            !    state. Sensitivity will be kept in system optimizer **
	    If (LOptNomSens) Then
	      if (nVarElt==0) then
                WRITE(*,*)' '
                WRITE(*,*)' No Variable Element for Optimization!'
                goto 108
              else if (nOptFov==0) then
                WRITE(*,*)' '
                WRITE(*,*)' No Field of View for Optimization!'
                goto 108
              end if

	      if (LOptIfFEX .and. .not.ifStopSet) then
                WRITE(*,*)' Must set STOP before optimization with FEX!'
                WRITE(*,*)' '
                goto 108
              end if

	      ! Process input arguments to system optimer drivers
              varZerns=0; zern_off=1
	      varAsphs=0; asph_off=1
              k=1
              iVarElt=0

	      ! Process source first
              If (isVarElt(0)) Then
                iVarElt=iVarElt+1
                varElts(iVarElt)=0
                varDOFs((iVarElt-1)*mVarDOF+1:iVarElt*mVarDOF)
     &               =varEltDOF(1:mVarDOF,0)
                If (.not.ifStopSet .and. .not.OptSrcRpt_FLG) Then
                  if (varEltDOF(1,0)/=0 .or. varEltDOF(2,0)/=0
     &                .or. varEltDOF(3,0)/=0) then
                    WRITE(*,*)' Must set STOP before optimizing source!'
                    go to 108
                  end if
                End If
                EltName(0)='Source'
                nOptFov=1
                opt_fov(1:3,1,1)=ChfRayDir; opt_fov(1:3,2,1)=ChfRayPos
                zernElt(iVarElt)=0; asphElt(iVarElt)=0
              End If

	      Do iElt=1,nElt
                If (isVarElt(iElt)) Then
                  iVarElt=iVarElt+1
                  varElts(iVarElt)=iElt
                  varDOFs((iVarElt-1)*mVarDOF+1:iVarElt*mVarDOF)
     &               =varEltDOF(1:mVarDOF,iElt)

                  zernElt(iVarElt)=0  ! range 1:nVarElt
                  if (nOptEltZern(iElt)>0) then
                    zernElt(iVarElt)=nOptEltZern(iElt)
                    varZerns(zern_off:zern_off+nOptEltZern(iElt)-1)
     &                = OptEltZernTerm(1:nOptEltZern(iElt),iElt)
                    zern_off=zern_off+nOptEltZern(iElt)
                  end if

	          asphElt(iVarElt)=0  ! range 1:nVarElt
                  if (nOptEltAsph(iElt)>0) then
                    asphElt(iVarElt)=nOptEltAsph(iElt)
                    varAsphs(asph_off:asph_off+nOptEltAsph(iElt)-1)
     &                = OptEltAsphTerm(1:nOptEltAsph(iElt),iElt)
                    asph_off=asph_off+nOptEltAsph(iElt)
                  end if

                  if (LOptCons) then
                    ! process variable constraints
                    do i=1,mVarDOF
                      if (varEltDOF(i,iElt)/=0) then
                        varEltCons(1:2,k)=EltCons(1:2,i,iElt)
                        k=k+1
                      end if
                    end do
                  end if
	        End If
              End Do

	      ! Save Rx npts, use opt_npts for design optimization
              npts_Rx=npts; npts=opt_npts

#if 1
	      WRITE(*,*)' '
	      NomSensOnly=.TRUE.; LinNomCtrl=.TRUE.
	      Call nls_optim_dvr(.true.,OptTarget,
     &                   nVarElt,varElts,varDOFs,
     &                   nWFZern,WFZernMode,
     &                   zernElt(1:nVarElt),varZerns,
     &                   asphElt(1:nVarElt),varAsphs,
     &                   opt_fov,nOptFov,fov_weight,
     &                   opt_wavelen,nOptWavelen,
     &                   dopt_tol,1,
     &                   old_wfe,new_wfe,rtn_flag,
     &                   tgtWF_flg,tgtWF_fn,tgtWF_fmt,
     &                   OptTgtElt,LOptIfFEX,beamOpt)
	     NomSensOnly=.FALSE.
	     WRITE(*,*)' '

	     If (rtn_flag /= 0) goto 108

	     ! Restore Rx state
	     ! Restore pre-optimization optical state by calling
             ! this routine with start_optim=.false. and
             ! nitr_dopt < 0.
             Call nls_optim_dvr(.false.,OptTarget,
     &                   nVarElt,varElts,varDOFs,
     &                   nWFZern,WFZernMode,
     &                   zernElt(1:nVarElt),varZerns,
     &                   asphElt(1:nVarElt),varAsphs,
     &                   opt_fov,nOptFov,fov_weight,
     &                   opt_wavelen,nOptWavelen,
     &                   dopt_tol,-1,
     &                   old_wfe,new_wfe,rtn_flag,
     &                   tgtWF_flg,tgtWF_fn,tgtWF_fmt,
     &                   OptTgtElt,LOptIfFEX,beamOpt)
#endif
	     npts=npts_Rx; nGridPts=npts+1

 108          If (rtn_flag /= 0) Then
                WRITE(*,*) ' ** Nominal sensitivity  calc aborted!'
              End If
	    End If  ! If (LOptNomSens)

#endif

#ifdef MACOS_CMD
           ! Preprocessing for 'BWK' command -jzlou
  	   If (BWK_FLG) Then
             iRay=1
             IF (ABS(zSource).GE.1d10) THEN
               CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate)
             ELSE
               CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                       ifPropagate)
             END IF
             iStartRay=iRay
             iEndRay=iRay
             iStartElt=0
             iEndElt=nElt
             ifRayPrint=.TRUE.; LBWKnoPrint=.TRUE.
             CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &                   ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &                   iStartRay,iEndRay,RayMap,SegMap,
     &                   ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &                   StartIndRef,StartExtinc,
     &                   DrawRayVec,DrawEltVec,nDrawElt(1),
     &                   xDraw,yDraw)
             ifRayPrint=.FALSE.; LBWKnoPrint=.FALSE.
             Do j=1,nElt
               CRIncidPosNom(1:3,j)=CRIncidPos(1:3,j)
             End Do
	  End If
#endif
	  ! Set default option for source adjustment
	  IF (ABS(zSource).LE.1d10) THEN
	    SAOpt=SADir  ! point source, ChfRayDir adjusted
	  ELSE
	    SAOpt=SAPos  ! collimated source, ChfRayPos adjusted
	  END IF

	  ! Set aperture stop if needed - must be last step in 'LOAD' cmd,
          ! since it jumps to label 210 -jzlou
          If (RxStopSet) Then
            if (ifStopSet) then
              WRITE(*,*)' '
              WRITE(*,120)
     &          ' Aperture stop set at (global):',StopPos(1:3)
              WRITE(*,*)' '
              GO TO 205
            else
              if (EltStopSet) then
                WRITE(*,*)' '
                WRITE(*,315)StopElt,StopOffset
                GO TO 210
              else
                WRITE(*,*)' *** Setting aperture stop failed!'
              WRITE(*,*)' '
              end if
            end if
          End If

#ifdef MACOS_CMD
C         *** Select a new FOV
          ELSE IF (LCMP(command,'SFOV',4)) THEN
            IF (ifLoad) THEN
	      WRITE(*,*)'Total FOV available: ',nOptChfRayPos
	      CALL IACCEPT_S(ifov,curFov,1,'Enter ID of FOV to use:')
	      If (ifov<1 .OR. ifov>nOptChfRayPos) Then
	        WRITE(*,*)' ** Invalid ID for FOV selection'
	      Else
	        curFov=ifov
	        ChfRayDir(1:3)=opt_fov(1:3,1,curFov)
	        ChfRayPos(1:3)=opt_fov(1:3,2,curFov)
	        WRITE(*,*)'   Set ChfRayDir =',opt_fov(1:3,1,curFov)
	        WRITE(*,*)'   Set ChfRayPos =',opt_fov(1:3,2,curFov)
	        WRITE(*,*)' '
	        ! Reinitialize some variables
                ifBuild=.FALSE.; ifPropagate=.FALSE.
                ifPerturb=.FALSE.; ifTrace=.FALSE.
                ifGrid=.FALSE.; ifDMBuild=.FALSE.
                i=mElt; iCurRayElt=0; iCurWFElt=0
	      End If
	    ELSE
	      WRITE(*,150)
	    END IF

C	  *** Select a new wavelength
          ELSE IF (LCMP(command,'SWL',3)) THEN
	    IF (ifLoad) THEN
	      WRITE(*,*)nOptWavelen,'wavelengths available: '
	      If (nOptWavelen>0) Then
		WRITE(*,*)' ',opt_wavelen(1:nOptWavelen)
	        WRITE(*,*)' '
	      End If
	      CALL IACCEPT_S(iwl,curwl,1,'Enter ID of wavelength to use')
              If (iwl<1 .OR. iwl>nOptWavelen) Then
                WRITE(*,*)' ** Invalid wavelength ID ',iwl
	      Else
	        curwl=iwl
	        Wavelen=opt_wavelen(curwl)
	        WRITE(*,*)'  Wavelength set to ',Wavelen
	        WRITE(*,*)' '
	        Call SetIndRef4Refractors(curwl)

	        ! Reinitialize some variables
                ifBuild=.FALSE.; ifPropagate=.FALSE.
                ifPerturb=.FALSE.; ifTrace=.FALSE.
                ifGrid=.FALSE.; ifDMBuild=.FALSE.
                i=mElt; iCurRayElt=0; iCurWFElt=0
	      End If
	    ELSE
              WRITE(*,150)
            END IF
#endif


C         ***next command: Find a source ray for a given position in source
C            aperture, -jzlou***

	  ELSE IF (LCMP(command,'FSR',3)) THEN
            IF (ifLoad) THEN
	      dsrp(1:2)=0d0
	      Call DACCEPT(srp,dsrp,2,
     &          'Enter ray position in source aperture:')

	      Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),srid,.TRUE.)
	      write(*,*) 'Ray id: ', srid
	    ELSE
              WRITE(*,150)
            END IF

C         *** Next command: Find reference ('marginal') rays, fairly
C           * robust now; works for any shape of source aperture and
C           * sparse pupil aperture. -jzlou

	  ELSE IF (LCMP(command,'RRAY',3)) THEN
	    IF (ifLoad) THEN
	      WRITE(*,163) '  Source aperture diameter = ',
     &                     Aperture,BaseUnits
	      ! Try to trace chief ray
              iStartElt=0; iEndElt=nElt
              ifRayPrint=.FALSE.
              CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &         ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &         ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,iCurRayElt,
     &         iEndElt,Dicr,Djcr,dxSource,nBadRays,StartIndRef,
     &         StartExtinc,Cmatrix,DrawRayVec,DrawEltVec)
              !
	      If (LRayOK(1)) Then
	        refRay(1)=1;
	      Else
		refRay(1)=-1;
	      End If
	      WRITE(*,*)' At source aperture: '
	      srp(1)=1.0; srp(2)=0.0
	      Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(2),.TRUE.)
	      srp(1)=0.0; srp(2)=1.0
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(3),.TRUE.)
	      srp(1)=-1.0; srp(2)=0.0
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(4),.TRUE.)
	      srp(1)=0.0; srp(2)=-1.0
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(5),.TRUE.)
	      If (refRay(1)==1) Then
	        WRITE(*,*)
     &             ' Reference ray indices (chief + 4 marginal) = ',
     &             refRay(1:5)
	      Else
	        WRITE(*,*)
     &             ' Reference ray indices (4 marginal) = ',
     &             refRay(2:5)
	        WRITE(*,*) ' Chief ray failed to pass through system'
	      End If
	    ELSE
              WRITE(*,150)
            END IF

C
C    Find "pupil rays" indexes, eight rays total
C
	  ELSE IF (LCMP(command,'PRAY',3)) THEN
            IF (ifLoad) THEN
              WRITE(*,163) '  Source aperture diameter = ',
     &                     Aperture,BaseUnits
              ! Try to trace chief ray
              iStartElt=0; iEndElt=nElt
              ifRayPrint=.FALSE.
              CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &         ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &         ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,iCurRayElt,
     &         iEndElt,Dicr,Djcr,dxSource,nBadRays,StartIndRef,
     &         StartExtinc,Cmatrix,DrawRayVec,DrawEltVec)
              !
              If (LRayOK(1)) Then
                refRay(1)=1;
              Else
                refRay(1)=-1;
              End If
              WRITE(*,*)' At source aperture: '
              srp(1)=9.0; srp(2)=0.0
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(2),.TRUE.)
              srp(1)=0.6364; srp(2)=0.6364
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(3),.TRUE.)
              srp(1)=0.0; srp(2)=0.9
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(4),.TRUE.)
              srp(1)=-0.6364; srp(2)=0.6364
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(5),.TRUE.)
              srp(1)=-0.9; srp(2)=0.0
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(6),.TRUE.)
              srp(1)=-0.6364; srp(2)=-0.6364
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(7),.TRUE.)
              srp(1)=0.0; srp(2)=-0.9
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(8),.TRUE.)
              srp(1)=0.6364; srp(2)=-0.6364
              Call FSR(npts,nElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &                 srp(1),srp(2),refRay(9),.TRUE.)
	    ELSE
              WRITE(*,150)
	    End IF

C         *** Trace pupil rays, need to run 'PRAY' command first ***
C         *** Pupil rays are chief ray + eight marginal rays ***
C         *** TPR = Trace Pupil Rays (and print ray footprint at the final surface)

	  ELSE IF (LCMP(command,'TPR',3)) THEN
	    IF (refRay(1) < -2) THEN
              WRITE(*,163)
     &     '   ** Need to run RRAY command before executing TRR command'
            ELSE
	      ! Trace all reference rays
	      IF (ABS(zSource).GE.1d10) THEN
               CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
              ELSE
               CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
	      END IF
              Do iRay=1,mPupRay
                If (refRay(iRay) < 0) Then
	          Exit  ! Cycle would do the same
	        End If
                ! Trace a ray
                iStartElt=0
	        If (LRayTgtElt==.TRUE.) Then
                  iEndElt=RayTgtElt
                Else
                  iEndElt=nElt
                End If
                ifRayPrint=.FALSE.; ifBuild=.FALSE.;
                ifTrace=.FALSE.; ifPropagate=.FALSE.
                iStartRay=refRay(iRay)  ! refRay(iRay) is one of reference rays
                iEndRay=iStartRay
                iStartElt=0
                If (.not. LRayTgtElt) Then
                  iEndElt=nElt
                Else
                  iEndElt=RayTgtElt
                End If
                CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &                 ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &                 iStartRay,iEndRay,RayMap,SegMap,
     &                 ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &                 StartIndRef,StartExtinc,
     &                 DrawRayVec,DrawEltVec,nDrawElt(1),
     &                 xDraw,yDraw)
                print*,' **CTRACE iRay =',refRay(iRay),',
     &                 iEndElt =',iEndElt
                print*,'     RayPosHist(1:3,iray,iEndElt) =',
     &                 RayPosHist(1:3,refRay(iRay),iEndElt)

	        ! Pupil rays are indexed 2 to 9
	        RRayPosHis(1:3,iRay) =
     &              RayPosHist(1:3,refRay(iRay),iEndElt)
	      End Do ! Do iRay
	    END IF


C         ***Set up interpolated surface calculations***

	  ELSE IF (LCMP(command,'SINT',3)) THEN
	    IF (ifLoad) THEN
	      WRITE(*,*)' Computing interpolated surface data...'

	      CALL RTIME(StartTime)
	      CALL SetIntSrf(ifIntSrfSet)
	      CALL RTIME(EndTime)
	      DTime=EndTime-StartTime

	      WRITE(*,116)DTime
	      ifBuild=.FALSE.
	      ifPropagate=.FALSE.
	      ifPerturb=.FALSE.
	      ifTrace=.FALSE.
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Set up interpolated surface calculations***

	  ELSE IF (LCMP(command,'UDSINIT',3)) THEN
	    IF (ifLoad) THEN
	      CALL IACCEPT_S(iElt,0,1,
     &        'Enter number of element (0=aperture):')
	      IF (iElt.GT.nElt) iElt=nElt
	      IF (iElt.LT.0) iElt=0
	      CALL CACCEPT(UDSrfFile(iElt),'none',
     &	     ' Enter user-defined surface file name ')
	      CALL UDSinit(ifUDSinit,iElt)
	      IF (ifUDSinit) THEN
	        ifBuild=.FALSE.
	        ifPropagate=.FALSE.
	        ifPerturb=.FALSE.
	        ifTrace=.FALSE.
	      END IF
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Modify optical description data***

	  ELSE IF (LCMP(command,'MODIFY',3)) THEN
#ifdef MACOS_CMD
	    IF (ifLoad) THEN
	      !print*,'b4 CALL Mod_Loh'
	      CALL MOD_LOH(npts)
#else
   ! for SMACOS
	      WRITE(*,*)' SMACOS MODIFY reset performed.'
#endif
	      ifBuild=.FALSE.
	      ifPropagate=.FALSE.
	      ifPerturb=.FALSE.
	      ifTrace=.FALSE.
	      ifGrid=.FALSE.
	      ifDMBuild=.FALSE.
	      i=mElt
	      iCurRayElt=0
	      iCurWFElt=0
#ifdef MACOS_CMD
	    ELSE
	      WRITE(*,150)
	    END IF
#endif

C         ***Summarize optical description data***

	  ELSE IF (LCMP(command,'SUMMARIZE',3)) THEN
	    IF (ifLoad) THEN
	      CALL SumData(GridType,npts)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Reset defaults, options, etc.***

	  ELSE IF (LCMP(command,'RESET',3)) THEN
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifPropagate=.FALSE.
	    ifCmt=.FALSE.
	    ifPerturb=.FALSE.
	    DO iElt=1,mElt
	      ifRegrid(iElt)=.FALSE.
	    END DO
	    ifPol=.FALSE.
	    ifVecDif3=.FALSE.
	    ifOld=.FALSE.
	    ifRayPrint=.FALSE.
	    ifGrid=.FALSE.
	    ifPixLoc=.FALSE.
#ifdef MACOS_CMD
	    ifAdd=.FALSE.; ifAddCA=.FALSE.
	    ifPixArray=.FALSE.
#endif
	    ifDMBuild=.FALSE.
	    ifCentroid=.TRUE.
	    ifFilterData=.FALSE.
	    ifLNegOK=.TRUE.
	    ifUDSinit=.FALSE.
            If (.not. RxStopSet) Then
	      ! don't want to reset STOP again when
              ! STOP is pre-set in Rx. 12/2015
 	      ifStopSet=.FALSE.
	    End If
	    i=mElt
	    CALL LZERO(ifIntSrfSet,i)
	    iCurRayElt=0
	    iCurWFElt=0
	    iObsOpt=1
	    StretchType=1
	    spcOption=1
	    If (UDBeam_FLG) Then
              BeamType=5
            Else
              BeamType=1
            End If
	    !BeamType=1
	    sizPix=0d0
	    EltPix(1)=0d0
	    EltPix(2)=0d0
	    WinCenPix(1)=0d0
	    WinCenPix(2)=0d0

C         ***Status check***

	  ELSE IF (LCMP(command,'STATUS',3)) THEN
	    IF (ifLoad) THEN
	      IF (iPloTyp.EQ.1) THEN
	        CTITLE='SLICE'
	      ELSE IF (iPloTyp.EQ.2) THEN
	        CTITLE='GRAY'
	      ELSE IF (iPloTyp.EQ.3) THEN
	        CTITLE='CONTOUR'
	      ELSE IF (iPloTyp.EQ.4) THEN
	        CTITLE='WIRE'
	      ELSE IF (iPloTyp.EQ.5) THEN
	        CTITLE='TEXT'
	      ELSE IF (iPloTyp.EQ.6) THEN
	        CTITLE='BINARY'
	      ELSE IF (iPloTyp.EQ.7) THEN
	        CTITLE='LINE'
	      ELSE IF (iPloTyp.EQ.8) THEN
	        CTITLE='FITS'
	      ELSE IF (iPloTyp.EQ.9) THEN
	        CTITLE='MAT'
	      ELSE IF (iPloTyp.EQ.0) THEN
	        CTITLE='NULL'
	      END IF
	      WRITE(*,207)filnam(1:ICLEN(filnam)),iCurRayElt,iCurWFElt,
     &	       iObsOpt,CTITLE(1:ICLEN(CTITLE)),ifAdd,ifPixLoc,ifPol
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Save new optical description data***

	  ELSE IF (LCMP(command,'SAVE',3)) THEN
	    IF (ifLoad) THEN
	      CALL CSAVE(npts)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Show element data***

	  ELSE IF (LCMP(command,'SHOW',3)) THEN
	    IF (ifLoad) THEN
	      CALL IACCEPT_S(iElt,0,1,
     &        'Enter number of element (0=aperture):')
	      IF (iElt.GT.nElt) iElt=nElt
	      IF (iElt.LT.0) iElt=0
	      CALL ShowData(npts,iElt)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Set obscuration option***

	  ELSE IF (LCMP(command,'OBS',3)) THEN
	    CALL CACCEPT(ANS,'POSITIVE',
     &	    ' Enter ray-trace obscuration option '//
     &	    '(ALL, POSITIVE, NEGATIVE):')
	    IF (LCMP(ANS,'ALL',1).OR.LCMP(ANS,'0',1)) THEN
	      iObsOpt=0
	      WRITE(*,*)' All rays plotted on spot diagrams'
	    ELSE IF (LCMP(ANS,'POS',1).OR.LCMP(ANS,'1',1)) THEN
	      iObsOpt=1
	      WRITE(*,*)' Unobscured rays only plotted on spot diagrams'
	    ELSE IF (LCMP(ANS,'NEG',1).OR.LCMP(ANS,'2',1)) THEN
	      iObsOpt=2
	      WRITE(*,*)' Obscured rays only plotted on spot diagrams'
	    END IF

C         ***Turn on polarization option***

	  ELSE IF (LCMP(command,'POLARIZATION',3)) THEN
	    ifPol=.TRUE.
    	    ifPropagate=.FALSE.
	    WRITE(*,*)' Polarization turned on'
	    IF (mWF.GE.3) THEN
	      WRITE(*,*)' Vector diffraction enabled'
	      ifVecDif3=.TRUE.
	    ELSE
	      ifVecDif3=.FALSE.
	    END IF
	    DV1(1)=0d0
	    DV1(2)=0d0
	    CALL DACCEPT(DV2,DV1,2,'Enter Ex0 (real,imag):')
	    Ex0=DCMPLX(DV2(1),DV2(2))
	    CALL DACCEPT(DV2,DV1,2,'Enter Ey0 (real,imag):')
	    Ey0=DCMPLX(DV2(1),DV2(2))

C         ***Turn off polarization option***

	  ELSE IF (LCMP(command,'NOPOLARIZATION',4)) THEN
	    IF (ifPol) ifPropagate=.FALSE.
	    ifPol=.FALSE.
      	    ifVecDif3=.FALSE.
	    WRITE(*,*)' Polarization turned off'

C         ***Turn on scalar diffraction***

	  ELSE IF (LCMP(command,'SCALAR',3)) THEN
	    ifVecDif3=.FALSE.
    	    ifPropagate=.FALSE.
	    WRITE(*,*)' Scalar diffraction enabled'

C         ***Turn on vector diffraction***

	  ELSE IF (LCMP(command,'VECTOR',3)) THEN
	    IF (ifPol) THEN
	      WRITE(*,*)' Vector diffraction enabled'
	      ifVecDif3=.TRUE.
    	      ifPropagate=.FALSE.
	    ELSE
	      WRITE(*,*)' Must turn on POLarization first'
	      ifVecDif3=.FALSE.
	    END IF

C         ***Set centroid option***

	  ELSE IF (LCMP(command,'CENTROID',5)) THEN
	    ifCentroid=.TRUE.
	    WRITE(*,*)' Image referenced to ray centroid'

C         ***Set chief ray option***

	  ELSE IF (LCMP(command,'CHIEFRAY',3)) THEN
	    ifCentroid=.FALSE.
	    WRITE(*,*)' Image referenced to chief ray'

C         ***Set regridding on for specified element***

	  ELSE IF (LCMP(command,'REGRID',3)) THEN
 	    CALL IACCEPT_S(iElt,0,1,'Enter element number:')
	    ifRegrid(iElt)=.TRUE.
	    ifPropagate=.FALSE.
	    !CALL DZERO(WFbuff,2*md2)
	    WFbuff=(0d0,0d0)
	    WRITE(*,*)' Rays will be regridded at element ',iElt

C         ***Turn regridding off for specified element***

	  ELSE IF (LCMP(command,'NOREGRID',5)) THEN
 	    CALL IACCEPT_S(iElt,0,1,'Enter element number:')
	    ifRegrid(iElt)=.FALSE.
	    WRITE(*,*)' Rays will not be regridded at element ',iElt

C         ***Allow negative ray lengths at reference surfaces***

	  ELSE IF (LCMP(command,'LNEG',4)) THEN
	    ifLNegOK=.TRUE.
	    WRITE(*,*)' Negative-length rays OK at reference surfaces'

C         ***Prohibit negative ray lengths***

	  ELSE IF (LCMP(command,'NOLNEG',4)) THEN
	    ifLNegOK=.FALSE.
	    WRITE(*,*)' No negative-length rays allowed'

C         ***Turn off pixel array location option***

	  ELSE IF (LCMP(command,'NOPLOC',4)) THEN
	    ! will remove the effect of 'WINDOW' command. -jzlou
	    ifPixLoc=.FALSE.
	    WRITE(*,*)' Images will be placed at center of pixel array'

C         ***Turn on pixel array location option***

	  ELSE IF (LCMP(command,'PLOCATE',4).OR.
     &	           LCMP(command,'WINDOW',3)) THEN
	    ifPixLoc=.TRUE.

C	  Get local coordinates:
	    IF (nOutCord.NE.5) THEN
	        CALL CACCEPT(ANS,'Beam',
     &	         'Enter output coordinate option (Enter or Beam):')
	    ELSE
	        CALL CACCEPT(ANS,'Tout',
     &	        'Enter output coordinate option (Tout, Enter or Beam):')
	    END IF

	    IF (LCMP(ANS,'YES',1).OR.LCMP(ANS,'TOUT',3)) THEN
	        ! Use TOUT option
	        CALL COMPACT(Tcompact,Tout,5,7,7,7)
	        CALL DCOMPACT(DTcompact,Tout,5,7,7,7)
	          ! COMPACT assigns 5x7 data from Tout to Tcompact,
	          ! Tout is 7x7 array and Tcompact is 5x7 array. -jzlou
	    ELSE IF (LCMP(ANS,'ENTER',3)) THEN
	        ! Use ENTER option
	        invec7(1)=-psiElt(3,iElt)
	        invec7(2)=-psiElt(1,iElt)
	        invec7(3)=-psiElt(2,iElt)
	        CALL DACCEPT(invec7,invec7,3,'Enter xOut (x,y,z) ')
	        CALL DUNITIZE(invec7)
	        CALL DXPROD(DDEF,psiElt,invec7)
	        CALL DNEGATE(DDEF,DDEF,3)
	        CALL DACCEPT(DDEF,DDEF,3,'Enter yOut (x,y,z) ')
	        CALL DUNITIZE(DDEF)
	        CALL DZERO(Tspot(1,1),49)
	        DO iCoord=1,3
	          Tspot(1,iCoord)=invec7(iCoord)
	          Tspot(3,3+iCoord)=invec7(iCoord)
	          Tspot(2,iCoord)=DDEF(iCoord)
	          Tspot(4,3+iCoord)=DDEF(iCoord)
	        END DO
	        Tspot(5,7)=1
	        CALL COMPACT(Tcompact,Tspot,5,7,7,7)
	        CALL DCOMPACT(DTcompact,Tspot,5,7,7,7)
	    ELSE
	        ! Use Beam option
	        ifGrid=.FALSE.; ifDraw=.FALSE.
	        ifPropagate=.FALSE.;  ifTrace=.FALSE.
	        ifBuild=.FALSE.;  ifOPDCalc=.FALSE.
	        !iLocal=iElt  ! original from macos-2.86
	        iLocal=nElt   ! change by jzlou, 08/22/05
	        CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	         ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &	         ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &           Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &	         Cmatrix,DrawRayVec,DrawEltVec,
     &	         iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)
	        IF (ifLocalCoord) THEN
	          WRITE(*,311)xLocal,yLocal,zLocal
	          CALL DZERO(Tspot(1,1),49)
	          DO iCoord=1,3
	            Tspot(1,iCoord)=xLocal(iCoord)
	            Tspot(3,3+iCoord)=xLocal(iCoord)
	            Tspot(2,iCoord)=yLocal(iCoord)
	            Tspot(4,3+iCoord)=yLocal(iCoord)
	          END DO
	          Tspot(5,7)=1
	          CALL COMPACT(Tcompact,Tspot,5,7,7,7)
	          CALL DCOMPACT(DTcompact,Tspot,5,7,7,7)
	        ELSE
	          GO TO 1
	        END IF
	    END IF  ! Tout, ENTER or BEAM option

	    CALL DACCEPT_S(sizPix,DDEF(1), 1,
     &	    	   'Enter pixel size for placing window:')
	    IF (sizPix.GT.0d0) THEN
	      CALL DACCEPT(EltPix,DDEF,2,
     &	        'Enter pixel coords of element vertex (x,y):')
	      CALL DACCEPT(WinCenPix,DDEF,2,
     &	        'Enter window location in pixel coords (x,y):')
	    ELSE
	      sizPix=0d0
	      EltPix(1)=0d0
	      EltPix(2)=0d0
	      WinCenPix(1)=0d0
	      WinCenPix(2)=0d0
	    END IF

C         ***Print ray maps***

	  ELSE IF (LCMP(command,'MAP',3)) THEN
	    IF (ifLoad) THEN
	      IF (.NOT.ifGrid) THEN
	        IF (ABS(zSource).GE.1d10) THEN
	          CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	          ifPropagate)
	        ELSE
	          CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	          ifPropagate)
	        END IF
	        ifGrid=.TRUE.
	      END IF
	      WRITE(*,217)
	      iLognRay=LOG10(FLOAT(nRay))
	      IF (iLognRay.LE.0) iLognRay=1
	      nAcross(1)=26
	      nAcross(2)=20
	      nAcross(3)=16
	      nAcross(4)=13
	      nAcross(5)=11
	      nAcross(6)=9
	      k=npts+1
	      iSkip=1+k/nAcross(iLognRay)
#ifdef MACOS_CMD
	      DO j=k,1,-iSkip
	        IF (iLognRay.EQ.1) THEN
	          WRITE(*,221)(RayMap(i,j),i=1,k,iSkip)
	        ELSE IF (iLognRay.EQ.2) THEN
	          WRITE(*,222)(RayMap(i,j),i=1,k,iSkip)
	        ELSE IF (iLognRay.EQ.3) THEN
	          WRITE(*,223)(RayMap(i,j),i=1,k,iSkip)
	        ELSE IF (iLognRay.EQ.4) THEN
	          WRITE(*,224)(RayMap(i,j),i=1,k,iSkip)
	        ELSE IF (iLognRay.GE.5) THEN
	          WRITE(*,225)(RayMap(i,j),i=1,k,iSkip)
	        END IF
              END DO
#else
   ! SMACOS
	      if (.not.allocated(srcMap%rayMap)) allocate(srcMap%rayMap(k,k))
              j1=1
	      DO j=k,1,-iSkip
                i1=1
	        !DO i=k,1,-iSkip
	        DO i=1,k,iSkip
                  srcMap%rayMap(i1,j1)=RayMap(i,j); i1=i1+1
	        END DO
	        j1=j1+1
	      END DO
	      srcMap%size=j1-1
#endif
	      IF (GridType.GE.3) THEN
#ifdef MACOS_CMD
	        CALL IACCEPT_S(I,0,1,'Hit return to continue ')
	        WRITE(*,219)
	        DO j=k,1,-iSkip
	          IF (iLognRay.EQ.1) THEN
	            WRITE(*,221)(SegMap(i,j),i=1,k,iSkip)
	          ELSE IF (iLognRay.EQ.2) THEN
	            WRITE(*,222)(SegMap(i,j),i=1,k,iSkip)
	          ELSE IF (iLognRay.EQ.3) THEN
	            WRITE(*,223)(SegMap(i,j),i=1,k,iSkip)
	          ELSE IF (iLognRay.EQ.4) THEN
	            WRITE(*,224)(SegMap(i,j),i=1,k,iSkip)
	          ELSE IF (iLognRay.GE.5) THEN
	            WRITE(*,225)(SegMap(i,j),i=1,k,iSkip)
	          END IF
                END DO
#else
   ! SMACOS
	        if (.not.allocated(srcMap%segMap)) allocate(srcMap%segMap(k,k))
                j1=1
	        DO j=k,1,-iSkip
	          i1=1
	          DO i=k,1,-iSkip
                    srcMap.segMap(i1,j1)=SegMap(i,j); i1=i1+1
		  END DO
	          j1=j1+1
	        END DO
#endif
	      END IF ! if GridType.GE.3
	    ELSE
	      WRITE(*,150)
	    END IF

C         *** 'ray' command ***

	  ELSE IF (LCMP(command,'RAY',2)) THEN

	    IF (ifLoad) THEN
	      ifRayPrint=.TRUE.
	      ifBuild=.FALSE.
	      ifTrace=.FALSE.
	      ifPropagate=.FALSE.
  4	      CONTINUE
	        CALL IACCEPT_S(iRay,0,1,
     &          	'Enter number of ray (1=chief ray, 0=quit):')
	        IF (iRay.EQ.0) GO TO 5
	        IF (ABS(zSource).GE.1d10) THEN
	          CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	          		 ifPropagate)
	        ELSE
	          CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	          	        ifPropagate)
	        END IF
	        IF ((iRay.LT.1).OR.(iRay.GT.nRay)) THEN
		  WRITE(*,*)' Invalid ray number. '
		  GO TO 4
		END IF
	        iStartRay=iRay
	        iEndRay=iRay
	        iStartElt=0
	        If (.not. LRayTgtElt) Then
 	          iEndElt=nElt
	        Else
	          iEndElt=RayTgtElt
	        End If
	        ifRayPrint=.TRUE.
	        CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	               ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &	               iStartRay,iEndRay,RayMap,SegMap,
     &	               ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &                 StartIndRef,StartExtinc,
     &                 DrawRayVec,DrawEltVec,nDrawElt(1),
     &		       xDraw,yDraw)
C	        print*,'**CTRACE cmd: RayPosHist(1:3,iray,iEndElt) =',
C    &                 RayPosHist(1:3,iRay,iEndElt)
C	        GO TO 4
  5	      CONTINUE
	      ifRayPrint=.FALSE.
	    ELSE
	      WRITE(*,150)
	    END IF

#if 0
! for McMichael's alignment work only
	    Do iElt=1,nElt
	      print*,'iElt=',iElt,', CRIncidDir=',
     &          CRIncidDir(1:3,iElt)
	      print*,'iElt=',iElt,', CRIncidPos=',
     &          CRIncidPos(1:3,iElt)
              print*,'iElt=',iElt,', CRNhat=',
     &          CRNhat(1:3,iElt)
	    End Do
#endif


#ifdef MACOS_CMD
C         *** 'BWK' command ***

          ELSE IF (LCMP(command,'BWK',3)) THEN
            IF (ifLoad) THEN
	      if (.not. BWK_FLG) then
	        WRITE(*,*) '** Need set CalBWk= Y in prescription'
                GO TO 1
	      end if

              ! First trace chief ray
              iRay=1
              IF (ABS(zSource).GE.1d10) THEN
                CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                         ifPropagate)
              ELSE
                CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate)
              END IF
              iStartRay=iRay
              iEndRay=iRay
              iStartElt=0
              iEndElt=nElt
              ifRayPrint=.TRUE.; LBWKnoPrint=.TRUE.
              CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &               ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &               iStartRay,iEndRay,RayMap,SegMap,
     &               ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &               StartIndRef,StartExtinc,
     &               DrawRayVec,DrawEltVec,nDrawElt(1),xDraw,yDraw)
              ifRayPrint=.FALSE.; LBWKnoPrint=.FALSE.

	      ! Calc beam walk at perturbed state
              Do j=1,nElt
                CRIncidSurfN(1:3,j)=CRNhat(1:3,j)
                CALL DUNITIZE(CRIncidSurfN(1:3,j))
                CALL DXPROD(CRYLoc(1:3,j),
     &                  -CRIncidDir(1:3,j),-CRIncidSurfN(1:3,j))
                CALL DUNITIZE(CRYLoc(1:3,j))
                CALL DXPROD(CRXLoc(1:3,j),
     &                      CRYLoc(1:3,j),CRIncidSurfN(1:3,j))
                CALL DUNITIZE(CRXLoc(1:3,j))

                ! Chief ray walk in global coords
                CRWalk(1:3,j)=CRIncidPos(1:3,j)-CRIncidPosNom(1:3,j)

                ! Chief ray walk projected onto element-local coords
              CRWalkProjLoc(1,j)=DDOTC(CRWalk(1:3,j),CRXLoc(1:3,j))
              CRWalkProjLoc(2,j)=DDOTC(CRWalk(1:3,j),CRYLoc(1:3,j))
              CRWalkProjLoc(3,j)
     &          =DDOTC(CRWalk(1:3,j),CRIncidSurfN(1:3,j))

             WRITE(*,162)' iElt=',j,
     &         '   ChfRayWalk =',CRWalkProjLoc(1:3,j)
             !WRITE(*,*)'iElt=',j,', ChfRay walk = ',CRWalkProjLoc(1:3,j)
              End Do
            ELSE
              WRITE(*,150)
            END IF
#endif


C	  ***next command: Calculate system effective focal length***

	  ELSE IF (LCMP(command,'EFL',3)) THEN

! EFL command code contained in efl.inc
#include "efl.inc"
	    If (ifLoad) Then
#ifdef MACOS_CMD
	      If (BaseUnits_FLG) Then
                WRITE(*,*)'  Estimated system effective focal length =',
     &                     efl,BaseUnits(1:ICLEN(BaseUnits))
	      Else
	        WRITE(*,*)'  Estimated system effective focal length =',
     &                    efl
	      End If
#endif
	      If (abs(Aperture) > smallestReal) Then
                focalRatio=efl/Aperture
	        If (LCMP(BaseUnits,'mm',2) .OR.
     &	            LCMP(BaseUnits,'none',4)) Then
	          plateScale=206265d0/(Aperture*focalRatio)
	        Else If (LCMP(BaseUnits,'m',1)) Then
	          plateScale=206265d0/(Aperture*1d3*focalRatio)
	        End If
#ifdef MACOS_CMD
      	        WRITE(*,*) '  Estimated system focal ratio =',
     &                     focalRatio
	        WRITE(*,*) '  Estimated system plate scale =',
     &                     plateScale,'arcsec/mm'
	        If (LCMP(BaseUnits,'none',4)) Then
	          WRITE(*,*) '  Plate scale calculated assuming base unit'
	          WRITE(*,*) '  is mm'
	        End If
#endif
	      End If
C               WRITE(*,*)'(** Warning: for decentered system, '//
C    &                    'the computed EFL could be inaccurate)'
	      ifEFLCalc=.true.
	    End If

C         ***next command: Print rays surface - to - surface***

	  ELSE IF (LCMP(command,'SRAY',3)) THEN
	    !
	    ! This command is buggy, doesn't work  - jzlou
	    !
	    IF (ifLoad) THEN
	      ifRayPrint=.TRUE.
	      ifBuild=.FALSE.
	      ifTrace=.FALSE.
	      ifPropagate=.FALSE.
 	      CALL IACCEPT_S(iRay,0,1,
     &               'Enter number of ray (1=chief ray, 0=quit):')
	      CALL IACCEPT_S(iElt,iCurRayElt,1,
     &	                     'Enter starting element:')
 	      CALL IACCEPT_S(iEndElt,iCurRayElt,1,
     &	                     'Enter end element:')

	      IF ((iElt.GT.nElt).OR.(iEndElt.GT.nElt).OR.
     &	      (iElt.GT.iEndElt).OR.(iElt.LT.0).OR.
     &	      (iEndElt.LT.0)) THEN
	        WRITE (*,*)' Invalid element sequence.'
	        GO TO 1
	      END IF

	      print*,'iElt,iCurRayElt =',iElt,iCurRayElt

	      IF (iElt.GT.iCurRayElt) THEN
	        iStartElt=0
	        IF (ABS(zSource).GE.1d10) THEN
	          CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	          ifPropagate)
	        ELSE
	          CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	          ifPropagate)
	        END IF
	      ELSE
	        iStartElt=iCurRayElt
	      END IF

	      iStartRay=iRay
	      iEndRay=iRay
	      IF (iStartElt.LT.iElt) THEN
	        ifRayPrint=.FALSE.
	        CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	         ifLNegOK,Cmatrix,iStartElt,iElt,iStartRay,iEndRay,
     &	   RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &     StartIndRef,StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &		 xDraw,yDraw)
	      END IF
	      ifRayPrint=.TRUE.
	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	       ifLNegOK,Cmatrix,iElt,iEndElt,iStartRay,iEndRay,
     &    RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &    StartIndRef,StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &	       xDraw,yDraw)
	      ifRayPrint=.FALSE.
	    ELSE
	      WRITE(*,150)
	    END IF

#ifdef MACOS_CMD
  ! Trace a single ray with a designated position in source aperture
	  ELSE IF (LCMP(command,'SRT',3)) THEN
	    IF (ifLoad) THEN
	      ifSRT=.TRUE.
	      ! Get last element to trace to
	      CALL IACCEPT_S(iEndElt,nElt,1,
     &                       'Enter last element to trace to:')
	      ! Get ray position in source aperture
	      Call DACCEPT(srp,dsrp,2,
     &          'Enter ray position in normalized source aperture:')
	      If (abs(srp(1))>0.5d0 .OR.abs(srp(2))>0.5d0) Then
	        WRITE(*,*)' **Ray position is outside source aperture'
	        GO TO 1
	      End If

	      ifRayPrint=.TRUE.
	      Call SRTrace(srp(1),srp(2),ifPol,ifGrid,ifPropagate,
     &         ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &         ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &         iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &         nBadRays,StartIndRef,StartExtinc,
     &         Cmatrix,DrawRayVec,DrawEltVec,rayState)
	      ifRayPrint=.FALSE.; ifSRT=.FALSE.; ifTrace=.FALSE.
	    ELSE
              WRITE(*,150)
            END IF
#endif

#ifdef MACOS_CMD
C	  ***Reset pgplot graphics device***
	  ELSE IF (LCMP(command,'PGD',3)) THEN
	    ifPlot=.true.
	    CALL GRAINI
#endif

C         ***Draw system layout and rays***

	  ELSE IF (LCMP(command,'DRAW',3) .OR. LCMP(command,'ODRAW',3)) THEN
	    DrawMode=.TRUE.  ! turn on DRAW mode, defined in src_mod
	    IF (LCMP(command,'ODRAW',3)) Then
              IsOverDraw=.TRUE.
	    ELSE
	      IsOverDraw=.FALSE.
	    END IF
	    IF (ifLoad) THEN
	      IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.5).AND.(iPloTyp.NE.6)
     &	          .AND.(iPloTyp.NE.8).AND.(iPloTyp.NE.9)
     &	          .AND.(.NOT.ifPlot)) THEN
	        ifPlot=.TRUE.
	        If (.not.IsOverDraw) CALL GRAINI
	      END IF
	      ifRayPrint=.FALSE.
	      ifBuild=.FALSE.
	      ifTrace=.FALSE.
	      ifPropagate=.FALSE.
	      nDrawRay=MIN(mRay,mDrawRay)
	      IF (MOD(nDrawRay,2).EQ.0) nDrawRay=nDrawRay-1
	      k=1

	      IF (ABS(zSource).GE.1d10) THEN
	        CALL ColSource(nDrawRay,Dicr,Djcr,dxSource,ifPol,
     &	                       ifPropagate)
	      ELSE
	        CALL PtSource(nDrawRay,Dicr,Djcr,dxSource,ifPol,
     &	                      ifPropagate)
	      END IF

	      iStartElt=0
	      iEndElt=nElt
	      CALL IACCEPT_S(iElt,0,1,
     &	                     'Enter first element to include:')
	      IF (iStartElt.LT.0) iStartElt=0
	      IF (iStartElt.GT.nElt) iStartElt=nElt
	      CALL IACCEPT_S(iEndElt,nElt,1,
     &	                     'Enter last element to include:')
	      IF (iEndElt.LT.iStartElt) iEndElt=iStartElt
	      IF (iEndElt.GT.nElt) iEndElt=nElt

	      CALL CACCEPT(ANS,'XZ',
     &         'Enter drawing plane (XZ, YZ or XY in source coords): ')
	      IF (LCMP(ANS,'YZ',2)) THEN
	        k=1; ANS='YZ'
	        If (.not.LudDrawGrid) Then
	          ! Original: draw frames dictated by zGrid,yGrid recomputed in source routine
	          CALL DEQUATE(xDraw,zGrid,3)
	          CALL DEQUATE(yDraw,yGrid,3)
	        Else
	          Call DXPROD(UDzDrawGrid,UDxDrawGrid,UDyDrawGrid)
	          xDraw(1:3)=UDzDrawGrid(1:3); yDraw(1:3)=UDyDrawGrid(1:3)
	        End If
	        XLabel='Z Axis'; YLabel='Y Axis'
	      ELSE IF (LCMP(ANS,'XY',2)) THEN
	        k=3; ANS='XY'
	        If (.not.LudDrawGrid) Then
 	          CALL DEQUATE(xDraw,xGrid,3)
	          CALL DEQUATE(yDraw,yGrid,3)
	        Else
	          xDraw(1:3)=UDxDrawGrid(1:3); yDraw(1:3)=UDyDrawGrid(1:3)
		End If
	        XLabel='X Axis'; YLabel='Y Axis'
	      ELSE
	        k=2; ANS='XZ'
	        If (.not.LudDrawGrid) Then
	          CALL DEQUATE(xDraw,zGrid,3)
	          CALL DEQUATE(yDraw,xGrid,3)
	        Else
	          Call DXPROD(UDzDrawGrid,UDxDrawGrid,UDyDrawGrid)
            xDraw(1:3)=UDzDrawGrid(1:3); yDraw(1:3)=UDxDrawGrid(1:3)
		End If
	        XLabel='Z Axis'; YLabel='X Axis'
	      END IF
	      pgplotDrawPlane=k  ! used in DRAW routine

C    Trace up to first element to be drawn
	      i=nDrawRay/2+1
	      IF (iElt.NE.0) THEN
	        ifDraw=.FALSE.
	        DO iDrawRay=1,nDrawRay
	          IF (k.EQ.1) THEN  ! 'YZ' option
	            iRay=RayMap(i,iDrawRay)
	          ELSE
	            iRay=RayMap(iDrawRay,i)
	          END IF
	          IF (iRay.GT.0) THEN
	            iStartRay=iRay
	            iEndRay=iRay
	            iDrawRay_global=iDrawRay
	            CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	              ifLNegOK,Cmatrix,
     &	              iStartElt,iElt,iStartRay,iEndRay,RayMap,SegMap,
     &	              ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &	              StartIndRef,StartExtinc,DrawRayVec(1,1,iDrawRay),
     &	              DrawEltVec(1,iDrawRay),nDrawElt(iDrawRay),
     &	              xDraw,yDraw)
	            iDrawRay_global=0
	            if (.false. .and. iDrawRay==1) then
	              print*,'DrawRayVec(1,1:2,iDrawRay)=',
     &		  	  DrawRayVec(1,1:2,iDrawRay)
		    end if
	          END IF
 	        END DO
	      END IF

	      if (.false.) then
	        print*,'==> DRAW: nDrawRay =',nDrawRay
	        print*,'==> (mid ray index) i =',i
	        print*,'==> iPloTyp =',iPloTyp
	        print*,'==> RayMap(51:55,1) =',RayMap(51:55,1)
	      end if

	      ifDraw=.TRUE.
	      DO iDrawRay=1,nDrawRay
	        IF (k.EQ.1) THEN  ! k value indicates draw plane, 1 => 'YZ'
	          iRay=RayMap(i,iDrawRay)
	        ELSE
	          iRay=RayMap(iDrawRay,i)
	        END IF

	        IF (iRay.GT.0) THEN
	          iStartRay=iRay
	          iEndRay=iRay
	          iDrawRay_global=iDrawRay
	          CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	            ifLNegOK,Cmatrix,iElt,iEndElt,
     &	            iStartRay,iEndRay,RayMap,SegMap,
     &	            ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &		    StartIndRef,StartExtinc,DrawRayVec(1,1,iDrawRay),
     &	            DrawEltVec(1,iDrawRay),nDrawElt(iDrawRay),
     &              xDraw,yDraw)
	          iDrawRay_global=0
	          if (.false. .and. iDrawRay==1) then
                    print*,'****** DrawRayVec(1:2,1,iDrawRay)=',
     &                DrawRayVec(1:2,1,iDrawRay)
                  end if
	        END IF
	        if (.false. .and. iDrawRay==33) then
                  print*,'****** iRay,nDrawElt(iDrawRay=33) =',
     &                   iRay,nDrawElt(iDrawRay)
                end if
 	      END DO

	      j=ICLEN(filnam)
	      IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &            .OR.(iPloTyp.EQ.9)) THEN

C       iPloTyp=1: Slice
C	iPloTyp=2: Screen output
C	iPloTyp=5: Text format file
C       iPloTyp=6: Binary format file
C	iPloTyp=8: Fits format file
C   	iPloTyp=9: Matlab format file
C	-jzlou

	        imgfil=filnam(1:j)//'.draw'
	        j=ICLEN(filnam)
	        CTITLE='Layout, '//ANS(1:2)//' Plane, File='//
     &	        filnam(1:j)
	        CALL DRAWOUT(iPloTyp,CTITLE,nDrawElt,nDrawRay,
     &	        	     DrawEltVec,DrawRayVec,imgfil,R1)

	      ELSE IF (iPloTyp.NE.0) THEN
C 	Produce PGPLOT output  -jzlou
	        CTITLE='Layout, '//ANS(1:2)//' Plane, File='//
     &	               filnam(1:j)
	        CALL DRAW(CTITLE,nDrawElt,nDrawRay,
     &	                  DrawEltVec,DrawRayVec,XLabel,YLabel)
	      END IF

	      ifBuild=.FALSE.
	      ifTrace=.FALSE.
	      ifPropagate=.FALSE.
	      ifDraw=.FALSE.
	    ELSE
	      WRITE(*,150)
	    END IF
	    DrawMode=.FALSE.  ! reset to default

C         ***Add atmospheric effects***

	  ELSE IF (LCMP(command,'ATMOSPHERE',3)) THEN
	    IF (.NOT.ifLoad) THEN
	      WRITE(*,150)
	      GO TO 1
	    END IF

C         Get element number:
	    IF (ifTrace.OR.ifPropagate) THEN
	      iElt=iCurRayElt
	    ELSE
	      iElt=1
	    END IF
	    CALL IACCEPT_S(iElt,iCurRayElt,1,
     &	           'Enter number of atmosphere phase screen element:')
	    IF ((iElt.LE.0).OR.(iElt.GT.nElt)) THEN
	      WRITE (*,*)'ATMOS: Invalid element number: ',iElt
	      GO TO 1
	    ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.11).OR.
     &	             (EltID(iElt).EQ.12)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iElt)
	      GO TO 1
	    END IF
	    iEndElt=iElt
	    iLocal=iElt

C	  Get atmospheric parameters
	    r0=Aperture/16d0
	    CALL RACCEPT_S(r0,r0,1,'Enter atmosphere r0:')
	    atmWavelen=Wavelen
	    CALL RACCEPT_S(atmWavelen,atmWavelen,1,
     &	                   'Enter atmosphere wavelength:')
	    pTilt=1d0
	    CALL RACCEPT_S(pTilt,pTilt,1,
     &	                   'Enter atmosphere tilt participation (0-1):')

C	  Get grid option:
	    CALL CACCEPT(ANS,'RAY',
     &                   'Use RAY grid or separate ATMospheric grid? ')

C	  Get local beam coordinates:
	    ifGrid=.FALSE.
	    ifDraw=.FALSE.
	    ifPropagate=.FALSE.
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifOPDCalc=.FALSE.
	    IF (.NOT.LCMP(ANS,'RAY',1)) THEN
	      iLocal=iElt
	      CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,
     &	      	ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,
     &	      	ifRegrid,ifIntSrfSet,ifRayPrint,iObsOpt,iEndElt,
     &	      	iCurWFElt,iCurRayElt,Dicr,Djcr,
     &	        dxSource,nBadRays,StartIndRef,StartExtinc,
     &	        Cmatrix,DrawRayVec,DrawEltVec,
     &	        iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)
	      IF (ifLocalCoord) THEN
	        WRITE(*,311)xLocal,yLocal,zLocal
	      ELSE
	        WRITE(*,*)' Unable to compute local coordinates.'
	        GO TO 1
	      END IF
	    END IF

C         Trace beam up to atmosphere surface:
	    CALL OPD(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	      ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	      iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	      nBadRays,StartIndRef,StartExtinc,OPDMat,WFE,WFEPV,DAvgL,
     &	      L1,Cmatrix,DrawRayVec,DrawEltVec)
	      ifOPDCalc=.TRUE.
	    CALL DXCALC(iElt,iEm1,dxMin(iElt),dxMax(iElt),
     &	      dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)

C         Add atmospheric effects to each ray
	    IF (LCMP(ANS,'RAY',1)) THEN
	      !NoiseSeed=0d0  ! test only
	      CALL AddAtmos(iElt,npts,iWF,R1,CD2,r0,atmWavelen,pTilt,
     &	        NoiseSeed,ifPol,ifGrid,ifPropagate,Dicr,Djcr,
     &	        dxSource,xGrid,yGrid,zGrid)

C	  Otherwise add atmospheric effects using separate grid
	    ELSE
	      rAtm=dxElt(iElt)
	      CALL DACCEPT_S(rAtm,rAtm,1,'Enter grid spacing:')
	      IF (rAtm.LT.dxElt(iElt)) THEN
	        dx=dxElt(iElt)
	        rAtm=dxElt(iElt)
	      ELSE
	        dx=rAtm
	      END IF
	      nAtmPts=INT(dx*npts/rAtm)+3
	      nAtmPts=MAX(mpts,nAtmPts)
	      j=(nAtmPts+1)/2
	      CALL AtmosPhase(CD2,R1,dx,r0,atmWavelen,nAtmPts,pTilt,
     &	                      NoiseSeed,mdttl)

	      DO i=1,3
	        crLocalPos(i)=RayPos(i,1)
	      END DO

	      C1=1d0/rAtm
	      DO iRay=2,nRay
	        IF (LRayOK(iRay)) THEN
	          DO i=1,3
	            CRSpot(i)=RayPos(i,iRay)-crLocalPos(i)
	          END DO
	          x=DDOTC(xLocal,CRSpot)/rAtm
	          y=DDOTC(yLocal,CRSpot)/rAtm
	          i0=IDFLOOR(x)+j
	          i1=i0+1
	          j0=IDFLOOR(y)+j
	          j1=j0+1
	          x=x-DBLE(i0-j)
	          y=y-DBLE(j0-j)
	          u00=R1(i0,j0)
	          u10=R1(i1,j0)
	          u01=R1(i0,j1)
	          u11=R1(i1,j1)

	          CALL INTNORM(x,y,u00,u01,u10,u11,dL,tiltx,tilty)

	          CumRayL(iRay)=CumRayL(iRay)+dL
	          DO i=1,3
	            RayDir(i,iRay)=RayDir(i,iRay)+
     &	            C1*(tiltx*xLocal(i)+tilty*yLocal(i))
	          END DO
	          CALL DUNITIZE(RayDir(1,iRay))

	        END IF
	      END DO
	    END IF

C         ***Compute linear models***

	  ELSE IF (LCMP(command,'BUILD',2).OR.LCMP(command,'DMBUILD',3))
     &	  THEN
	    IF (ifLoad) THEN
 26	      CALL IACCEPT_S(iElt,1,1,
     &	      'Enter terminal element number:')
	      IF (iElt.LE.0) THEN
	        GO TO 1
	      ELSE IF (iElt.GT.nElt) THEN
	        WRITE (*,*)'BUILD: Invalid element number: ',iElt
	        GO TO 26
	      ELSE IF ((EltID(iElt).EQ.9).OR.
     &	      (EltID(iElt).EQ.11).OR.(EltID(iElt).EQ.12)) THEN
	        WRITE(*,*)'Invalid element type =',EltID(iElt)
	        GO TO 26
	      END IF
	      iStartElt=0
	      iEndElt=iElt
	      iPertElt=iElt
	      iCurRayElt=iEndElt

	      CALL RTIME(StartTime)
	      IF (ABS(zSource).GE.1d10) THEN
	        CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	        ifPropagate)
	      ELSE
	        CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	        ifPropagate)
	      END IF
	      ifGrid=.TRUE.
              IF (nRay.LE.bRay) THEN
	        WRITE(*,118)nRay
	      ELSE
	        WRITE(*,119)nRay,bRay
	        GO TO 1
	      END IF
	      ifPropagate=.FALSE.
	      ifTrace=.TRUE.
	      ifBuild=.TRUE.
	      IF (LCMP(command,'DMBUILD',3)) ifDMBuild=.TRUE.
	      ifOPDCalc=.TRUE.
	      iCurRayElt=iEndElt
	      iLinElt=iEndElt
	      iStartRay=1
	      iEndRay=nRay

	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	      	ifLNegOK,Cmatrix,
     &	      	iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	      	ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	      	StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(nRay),
     &        	xDraw,yDraw)

	      CALL RTIME(EndTime)
	      DTime=EndTime-StartTime
	      WRITE(*,116)DTime
	      CALL DXCALC(iLinElt,iLinElt,dxMin(iElt),dxMax(iElt),
     &	      dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
	      IF (nBadRays.GT.0) CALL WARN(LEltOK,iEndElt,nBadRays)
	      DO 75 i=1,6
	        DO 74 j=1,mElt
		  PertVec(i,j)=0d0
 74	        CONTINUE
 75	      CONTINUE
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Print ray partials***

	  ELSE IF (LCMP(command,'PARTIALS',2)) THEN
	    IF (ifBuild) THEN
c 71	      CONTINUE
	      CALL IACCEPT_S(iRay,0,1,
     &        'Enter ray number for partials (0=done, 1=chief ray):')
	      IF (iRay.EQ.0) GO TO 73
	      IF ((iRay.LT.1).OR.(iRay.GT.nRay)
     &        .OR.(.NOT.LRayOK(iRay))) THEN
		WRITE(*,*)' Invalid ray number. '
		GO TO 73
	      END IF
	      CALL CACCEPT(ANS,'YES',
     &        'Print partials in actuator/sensor coordinates? ')
	      IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
		CALL COMPACT(Tcompact,Tout,nOutCord,7,7,7)
		WRITE(*,110)iRay,iPertElt,EltName(iPertElt),0,EltName(0)
		CALL RMPRDB(Cout,Tcompact,Cmatrix(1,1,iRay),
     &          nOutCord,7,7)
		CALL RMOUT(Cout,nOutCord,7)
		DO 78 iElt=1,iPertElt
		  IF (LRayToElt(iElt,iRay).AND.(nECoord(iElt).NE.0))
     &            THEN
		    i=8+(iElt-1)*6
		    WRITE(*,111)iRay,iPertElt,EltName(iPertElt),iElt,
     &              EltName(iElt)
		    CALL DTSP2(TEsp,Telt,6,6,iElt,IABS(nECoord(iElt)))
		    CALL RMPRDC(Celt,Cmatrix(1,i,iRay),TEsp(1,1,iElt),
     &              7,6,IABS(nECoord(iElt)))
		    CALL RMPRDB(Cout,Tcompact,Celt,
     &              nOutCord,7,IABS(nECoord(iElt)))
		    CALL RMOUT(Cout,nOutCord,IABS(nECoord(iElt)))
		  END IF
 78	        CONTINUE
	      ELSE
	        !print*,'*** printing partials in global coord'
		WRITE(*,110)iRay,iPertElt,EltName(iPertElt),0,EltName(0)
		CALL RMOUT(Cmatrix(1,1,iRay),7,7)
		DO 79 iElt=1,iPertElt
		  IF (LRayToElt(iElt,iRay)) THEN
		    i=8+(iElt-1)*6
		    WRITE(*,111)iRay,iPertElt,EltName(iPertElt),iElt,
     &              EltName(iElt)
		    CALL RMOUT(Cmatrix(1,i,iRay),7,6)
		  END IF
 79	        CONTINUE
	      END IF
c	      GO TO 71
 73           CONTINUE
	    ELSE
	      WRITE (*,152)
	    END IF

C         ***Set beam illumination type***

	  ELSE IF (LCMP(command,'BEAM',3)) THEN
	    IF (ifLoad) THEN
	      CALL CACCEPT(ANS,'UNIFORM',
     &       'Enter beam type (UNIFORM, GAUSSIAN, COS**POWER, DIPOLE):')
	      IF (LCMP(ANS,'UNIFORM',1).OR.LCMP(ANS,'1',1)) THEN
	        BeamType=1
	      ELSE IF (LCMP(ANS,'GAUSS',1).OR.LCMP(ANS,'2',1)) THEN
	        BeamType=2
	        S1=Aperture/4d0
	        CALL DACCEPT_S(rxBeam,S1,1,'Enter x beam waist radius:')
	        CALL DACCEPT_S(ryBeam,S1,1,'Enter y beam waist radius:')
	      ELSE IF (LCMP(ANS,'COS',1).OR.LCMP(ANS,'3',1)) THEN
	        BeamType=3
	        S1=Aperture/2d0
	        CALL DACCEPT_S(rxBeam,S1,1,'Enter cosine beam radius:')
	        ryBeam=rxBeam
	        S1=1d0
	        CALL DACCEPT_S(CosPower,S1,1,'Enter cosine exponent:')
	      ELSE IF (LCMP(ANS,'DIPOLE',1).OR.LCMP(ANS,'4',1)) THEN
	        BeamType=4
	      ELSE
	        WRITE(*,*) ' Unknown beam type.'
	      END IF
	      iCurWFElt=0
	      ifPropagate=.FALSE.
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Optimize reference surface***

	  ELSE IF (LCMP(command,'ORS',3)) THEN
 29	    CALL IACCEPT_S(iElt,1,1,
     &	           'Enter number of element to be optimized:')
	    IF (iElt.LE.0) THEN
	      GO TO 1
	    ELSE IF (iElt.GT.nElt) THEN
	      WRITE (*,*)'ORS: Invalid element number: ',iElt
	      GO TO 29
	    ELSE IF ((EltID(iElt).NE.3).AND.(EltID(iElt).NE.7)
     &	             .AND.(EltID(iElt).NE.8)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iElt)
	      GO TO 29
	    END IF
	    iStartElt=0
	    iEndElt=iElt-1

	    IF (ABS(zSource).GE.1d10) THEN
	      CALL ColSource(nTestPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    ELSE
	      CALL PtSource(nTestPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    END IF
	    ifGrid=.FALSE.
	    ifPropagate=.FALSE.
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifOPDCalc=.FALSE.
	    iStartRay=1
	    iEndRay=nRay
	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	          ifLNegOK,Cmatrix,
     &	          iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	          ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &            StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(nRay),
     &            xDraw,yDraw)

cdebug	    CALL RSTRACE(fElt(iElt),iElt,S1,i,VptElt(1,iElt),
c     &	    RptElt(1,iElt),psiElt(1,iElt),S2)
cdebug	    WRITE(*,*) ' Initial OPD=',S1

            CALL CRSOPTIMIZE(iElt)

C         ***Slave reference surface***

	  ELSE IF (LCMP(command,'SRS',3)) THEN
 57	    CALL IACCEPT_S(iSlv1,1,1,
     &	    'Enter number of element to be slaved:')
	    IF (iSlv1.LE.0) THEN
	      GO TO 1
	    ELSE IF (iSlv1.GT.nElt) THEN
	      WRITE (*,*)'SRS: Invalid element number: ',iElt
	      GO TO 57
	    ELSE IF ((EltID(iSlv1).NE.2).AND.(EltID(iSlv1).NE.3)
     &	    .AND.(EltID(iSlv1).NE.7).AND.(EltID(iSlv1).NE.8)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iElt)
	      GO TO 57
	    END IF

 58	    CALL IACCEPT_S(iSlv2,1,1,
     &	    'Enter number of element to slave to:')
	    IF (iSlv2.LE.0) THEN
	      GO TO 1
	    ELSE IF ((iSlv2.GT.nElt).OR.(iSlv2.EQ.iSlv1)) THEN
	      WRITE (*,*)'SRS-2: Invalid element number: ', iElt
	      GO TO 58
	    ELSE IF ((EltID(iSlv2).NE.2).AND.(EltID(iSlv2).NE.3)
     &	    .AND.(EltID(iSlv2).NE.7).AND.(EltID(iSlv2).NE.8)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iElt)
	      GO TO 58
	    END IF

	    ifGrid=.FALSE.
	    ifPropagate=.FALSE.
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifOPDCalc=.FALSE.
	    iStartRay=1
	    iEndRay=1

	    iStartElt=0
	    iEndElt=iSlv1
	    IF (ABS(zSource).GE.1d10) THEN
	      CALL ColSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    ELSE
	      CALL PtSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    END IF
	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	     ifLNegOK,Cmatrix,iStartElt,iEndElt,iStartRay,iEndRay,
     &       RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &       StartIndRef,StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &       xDraw,yDraw)

	    LSlv1=CumRayL(1)
	    signpsi=SIGN(1D0,DDOTC(RayDir(1,1),psiElt(1,iSlv1)))
	    CALL DSMPROD(psip,RayDir,signpsi,3)
	    CALL DEQUATE(Vptp,RayPos,3)

	    iEndElt=iSlv2
	    IF (ABS(zSource).GE.1d10) THEN
	      CALL ColSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    ELSE
	      CALL PtSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    END IF
	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	     ifLNegOK,Cmatrix,iStartElt,iEndElt,iStartRay,iEndRay,
     &       RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &       StartIndRef,StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &	     xDraw,yDraw)

	    ifReturn1=.FALSE.
	    DO j=2,(iSlv1-1)
	      IF (EltID(j).EQ.8) ifReturn1=.NOT.ifReturn1
	    END DO
	    ifReturn2=.FALSE.
	    DO j=2,(iSlv2-1)
	      IF (EltID(j).EQ.8) ifReturn2=.NOT.ifReturn2
	    END DO

	    IF (.NOT.(ifReturn1.OR.ifReturn2)) THEN
	      IF (fElt(iSlv1).GE.BIG) THEN
	        IF (fElt(iSlv2).LT.BIG) THEN
	          DO i=1,3
	            Vptp(i)=RayPos(i,1)+zElt(iSlv2)*RayDir(i,1)
	          END DO
	          zp=zElt(iSlv2)+CumRayL(1)-LSlv1
	        ELSE
	          IF ((zElt(iSlv1).GE.BIG).OR.(zElt(iSlv2).GE.BIG)) THEN
	            WRITE(*,*) 'WARNING: Plane-to-plane propagation',
     &	            ' should be set up with zElt < INF'
	          END IF
	          IF (CumRayL(1).GT.LSlv1) THEN
	            zp=zElt(iSlv2)+(CumRayL(1)-LSlv1)
	          ELSE
	            zp=zElt(iSlv2)-(CumRayL(1)-LSlv1)
	          END IF
	          DO i=1,3
	            Vptp(i)=RayPos(i,1)+zElt(iSlv2)*RayDir(i,1)
	          END DO
	          fp=INF
	        END IF
	      ELSE
	        IF (fElt(iSlv2).LT.BIG) THEN
	          zp=zElt(iSlv2)+CumRayL(1)-LSlv1
	          fp=ABS(zp)
	        ELSE
	          zp=CumRayL(1)-LSlv1
	          fp=INF
	        END IF
	      END IF

	    ELSE IF (ifReturn1.AND.ifReturn2) THEN
	      IF (fElt(iSlv1).GE.BIG) THEN
	        zp=INF
	        IF (fElt(iSlv2).LT.BIG) THEN
	          DO i=1,3
	            Vptp(i)=RayPos(i,1)+zElt(iSlv2)*RayDir(i,1)
	          END DO
	        END IF
	      ELSE
	        IF (fElt(iSlv2).LT.BIG) THEN
	          zp=zElt(iSlv2)-CumRayL(1)+LSlv1
	        ELSE
	          zp=CumRayL(1)-LSlv1
	        END IF
	      END IF
	      fp=ABS(zp)

	    ELSE
	      WRITE(*,*)' WARNING: SRS is confused about return status'
	      !fp=fElt(iSlv1)
	      !zp=zElt(iSlv1)
	    END IF

	    WRITE(*,109)fElt(iSlv1),fp,zElt(iSlv1),zp,psiElt(1,iSlv1),
     &	    psip(1),psiElt(2,iSlv1),psip(2),psiElt(3,iSlv1),psip(3),
     &	    VptElt(1,iSlv1),Vptp(1),VptElt(2,iSlv1),Vptp(2),
     &	    VptElt(3,iSlv1),Vptp(3)
	    CALL CACCEPT(ANS,'YES',
     &      'Accept the new element? ')
	    IF ((ANS(1:1).NE.'n').OR.(ANS(1:1).NE.'N')) THEN
	      eElt(iSlv1)=0d0
	      fElt(iSlv1)=fp
	      KcElt(iSlv1)=0d0
	      KrElt(iSlv1)=-fp
	      zElt(iSlv1)=zp
	      DO 59 i=1,3
	        psiElt(i,iSlv1)=psip(i)
	        VptElt(i,iSlv1)=Vptp(i)
	        RptElt(i,iSlv1)=Vptp(i)
  59	      CONTINUE
	    END IF


C         *** Find Detector (Focal) Plane Position and Orientation ***

          ELSE IF (LCMP(command,'FDP',3)) THEN
            IF (ifLoad) THEN
           CALL IACCEPT_S(iElt,nElt,1,'Enter detector element number:')
           Call FDPN(iElt,npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &       ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &       ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &       Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &       Cmatrix,DrawRayVec,DrawEltVec)
            ELSE
          WRITE(*,*)' Must load optical data first: use NEW or OLD.'
              GO TO 1
            END IF


C	  *** Compute system Line of Sight (LOS) ***
	  ! Based on V.N. Mahajan's paper: 'Line of Sight of an Aberrated Optical System',
	  ! using its pupil function formulation

	  ELSE IF (LCMP(command,'LOS',3)) THEN
	    IF (ifLoad) THEN
	      CALL IACCEPT_S(iElt,nElt-1,1,
     &          'Enter element where LOS is to be evaluated:')
	      If (iElt<0 .OR. iElt>nElt) Then
	        iElt=nElt-1  ! default to XP element
	      End If
	      iEndElt=iElt;
	      LCalcWFElt=.TRUE.  ! enabling computing WFElt in OPD()

	       CALL RunProp(iElt,npts,ifPol,ifGrid,ifPropagate,
     &     ifBuild,ifTrace,ifRegrid,ifVecDif3,ifIntSrfSet,ifLNegOK,
     &     iEndElt,iCurWFElt,iCurRayElt,iWF,Dicr,Djcr,dxSource,
     &     iObsOpt,nBadRays,StartIndRef,StartExtinc,DWF,OPDMat)

	       CALL DXCALC(iElt,iEm1,dxMin(iElt),dxMax(iElt),
     &            dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
	       ! dxElt(iElt) = dx of XP OPD grid for LOS command

	       ! Calc LOS, both linear and angular
	       iWF=1  ! assuming non-polarized
	       Call LOSCALC(zElt(iElt),dxElt(iElt),nGridpts,
     &                      OPDMat(1:mpts,1:mpts),
     &    		    mdttl,mpts,WFElt(1:mdttl,1:mdttl,iWF),
     &                      sysLOS,sysLOSAng)

	    ELSE
	 WRITE(*,*)' Must load optical data first: use NEW or OLD.'
              GO TO 1
	    END IF


C         *** Surface metrology calculation ***
	  ELSE IF (LCMP(command,'METcalc',3)) THEN
            IF (ifLoad) THEN
	      Call SrfMetCalc(nElt,iEltToMetSrf,nMetPos,tMetSrf,
     &                        SrfMetPos,SrfMetMea)
	    ELSE
	      WRITE(*,*)' Must load optical data first: use NEW or OLD.'
              GO TO 1
	    END IF

C         *** Find exit pupil (with option to minimize tilt) ***

	  ELSE IF (LCMP(command,'FEXIT',3)) THEN
            IF (ifLoad) THEN
 54	      CALL IACCEPT_S(iElt,(nElt-1),1,
     &	        'Enter number of exit pupil return surface:')
	      iEm1=iElt-1
	      IF (iElt.LE.0) THEN
	        GO TO 1
	      ELSE IF (iElt.GT.nElt) THEN
	        WRITE (*,*)'FEXIT: Invalid element number: ',iElt
	        GO TO 54
	      ELSE IF ((EltID(iElt).NE.8).AND.(EltID(iElt).NE.3)) THEN
	        WRITE(*,*)'Invalid element type =',EltID(iElt) ! must be Return
                                                               ! or Reference surface
	        GO TO 54
	      END IF
	    ELSE
	      WRITE(*,*)' Must load optical data first: use NEW or OLD.'
	      GO TO 1
	    END IF

	    IF (.NOT.ifStopSet) THEN
	      WRITE(*,*)
     &	      ' Must set stop before finding pupil. Use STOP command.'
	      GO TO 1
	    END IF

	    ifGrid=.FALSE.
	    ifPropagate=.FALSE.
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifOPDCalc=.FALSE.

	    CALL FEX(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	      ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	      iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	      nBadRays,StartIndRef,StartExtinc,
     &	      Cmatrix,DrawRayVec,DrawEltVec,
     &	      iElt,zp,psip,CrossPt,ifCentroid,
     &	      CRSpot,CentroidSpot,dCRCentroid)

	    WRITE(*,407)CRSpot,CentroidSpot,dCRCentroid

	    WRITE(*,107)fElt(iElt),ABS(zp),zElt(iElt),zp,psiElt(1,iElt),
     &	      psip(1),psiElt(2,iElt),psip(2),psiElt(3,iElt),psip(3),
     &	      VptElt(1,iElt),CrossPt(1),VptElt(2,iElt),CrossPt(2),
     &	      VptElt(3,iElt),CrossPt(3)
	    CALL CACCEPT(ANS,'YES',
     &                   'Accept the new element? ')
	    IF ((ANS(1:1).NE.'n').AND.(ANS(1:1).NE.'N')) THEN
	      eElt(iElt)=0d0
	      fElt(iElt)=ABS(zp)
	      KcElt(iElt)=0d0
	      KrElt(iElt)=-fElt(iElt)
	      zElt(iElt)=zp
	      CALL DEQUATE(psiElt(1,iElt),psip,3)
	      CALL DEQUATE(VptElt(1,iElt),CrossPt,3)
	      CALL DEQUATE(RptElt(1,iElt),CrossPt,3)
	    END IF

C         ***Find Field Point using global or pixel coordinates***

	  ELSE IF (LCMP(command,'FFP',3)
     &             .OR.LCMP(command,'PFP',3)) THEN
	    CALL IACCEPT_S(iSlv1,nElt,1,
     &	      'Enter element number:')
	    IF ((EltID(iSlv1).EQ.9).OR.(EltID(iSlv1).EQ.11)
     &	      .OR.(EltID(iSlv1).EQ.12)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iSlv1)
	      GO TO 384
	    ELSE IF ((iSlv1.GT.nElt).OR.(iSlv1.LE.0)) THEN
	      WRITE(*,*)'FFP: Invalid element number: ',iElt
	      GO TO 384
	    END IF
	    CALL DZERO(DDEF,3)

	    IF (.NOT.ifStopSet) THEN
	      WRITE(*,*)
     &	      ' Must set stop before FFP or PFP. Use STOP command.'
	      GO TO 1
	    END IF

	    IF (LCMP(command,'FFP',3)) THEN
	      CALL DACCEPT(xTarg,DDEF,2,
     &	        'Enter offset vector in global units (dx,dy):')
	    ELSE
	      CALL DACCEPT_S(S1,SizPix,1,
     &	      'Enter pixel size for placing image on pixel array:')
	      CALL DACCEPT(xTarg,DDEF,2,
     &	      'Enter image position in pixel units (dx,dy):')
	      xTarg(1)=(xTarg(1)-EltPix(1))*S1
	      xTarg(2)=(xTarg(2)-EltPix(2))*S1
	    END IF
	    dthPt=1d-9
	    ifGrid=.FALSE.
	    ifPropagate=.FALSE.
	    ifTrace=.FALSE.
	    ifBuild=.FALSE.
	    ifOPDCalc=.FALSE.
	    iStartRay=1
	    iEndRay=1
	    iStartElt=0
	    iEndElt=iSlv1
c--Nominal incidence point
	    CALL DEQUATE(cr1indir,ChfRayDir,3)
	    CALL DEQUATE(cr1inpos,ChfRayPos,3)
	    IF (ABS(zSource).GE.1d10) THEN
	      CALL ColSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    ELSE
	      CALL PtSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	      ifPropagate)
	    END IF
	    CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	      ifLNegOK,Cmatrix,
     &	    iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	      StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &	          xDraw,yDraw)
	    CALL DEQUATE(crPos,RayPos,3)
	    IF (nBadRays.NE.0) THEN
	      WRITE(*,*)'Nominal chief ray becomes undefined'
	      GO TO 385
	    END IF
c--Define elt coords
	    IF ((iSlv1.EQ.nElt).AND.(nOutCord.EQ.5)) THEN
	      DO i=1,3
	        xPt(i)=Tout(1,i)
	        yPt(i)=Tout(2,i)
 	      END DO
	    ELSE
	      CALL DSMPROD(th,xGrid,dthPt,3)
	      CALL DEQUATE(cr2indir,ChfRayDir,3)
	      CALL DEQUATE(cr2inpos,ChfRayPos,3)
	      CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
	      CALL DEQUATE(RayDir,cr2indir,3)
	      CALL DEQUATE(RayPos,cr2inpos,3)
	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	       ifLNegOK,Cmatrix,
     &	       iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	       ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	       StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &         xDraw,yDraw)
	      CALL DSUB(xPt,RayPos,crPos,3)  ! should be yPt instead? -jzlou
	      !
	      CALL DSMPROD(th,yGrid,dthPt,3)
	      CALL DEQUATE(cr2indir,ChfRayDir,3)
	      CALL DEQUATE(cr2inpos,ChfRayPos,3)
	      CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
	      CALL DEQUATE(RayDir,cr2indir,3)
	      CALL DEQUATE(RayPos,cr2inpos,3)
	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	       ifLNegOK,Cmatrix,
     &	       iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	       ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	       StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &         xDraw,yDraw)
	      CALL DSUB(yPt,RayPos,crPos,3) ! should be xPt instead? -jzlou
	    END IF
	    CALL DUNITIZE(xPt)
	    CALL DUNITIZE(yPt)
c--Incidence pt in local coords
	    CALL DSUB(DV1,crPos,VptElt(1,iSlv1),3)
	    xInt(1)=DDOTC(xPt,DV1)
	    xInt(2)=DDOTC(yPt,DV1)
	    xInt0(1)=xInt(1)
	    xInt0(2)=xInt(2)
c--Loop until desired pt is matched (max 100)
	    !DO k=1,100
            Do k=1,100
c----Check for convergence
	      CALL DSUB(dxPt,xInt,xTarg,2)
	      sqErr=dxPt(1)*dxPt(1)+dxPt(2)*dxPt(2)
	      IF (sqErr.LE.1d-18) GO TO 384
c----Compute x-gradient
	      CALL DSMPROD(th,xGrid,dthPt,3)
	      CALL DEQUATE(cr2indir,cr1indir,3)
	      CALL DEQUATE(cr2inpos,cr1inpos,3)
	      CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
	      CALL DEQUATE(RayDir,cr2indir,3)
	      CALL DEQUATE(RayPos,cr2inpos,3)
	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	       ifLNegOK,Cmatrix,
     &	       iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	       ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	       StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &	       xDraw,yDraw)
	      CALL DSUB(DV1,RayPos,crPos,3)
	      ddxPt(1)=DDOTC(DV1,xPt)
	      ddxPt(2)=DDOTC(DV1,yPt)
	      GPt(1,1)=ddxPt(1)/dthPt
	      GPt(2,1)=ddxPt(2)/dthPt
c----Compute y-gradient
	      CALL DSMPROD(th,yGrid,dthPt,3)
	      CALL DEQUATE(cr2indir,cr1indir,3)
	      CALL DEQUATE(cr2inpos,cr1inpos,3)
	      CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
	      CALL DEQUATE(RayDir,cr2indir,3)
	      CALL DEQUATE(RayPos,cr2inpos,3)
	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	       ifLNegOK,Cmatrix,
     &	       iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	       ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	       StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &         xDraw,yDraw)
	      CALL DSUB(DV1,RayPos,crPos,3)
	      ddxPt(1)=DDOTC(DV1,xPt)
	      ddxPt(2)=DDOTC(DV1,yPt)
	      GPt(1,2)=ddxPt(1)/dthPt
	      GPt(2,2)=ddxPt(2)/dthPt
c----Compute next step
	      S1=1d0/(GPt(1,1)*GPt(2,2)-GPt(1,2)*GPt(2,1))
	      GPtInv(1,1)=S1*GPt(2,2)
	      GPtInv(1,2)=-S1*GPt(1,2)
	      GPtInv(2,1)=-S1*GPt(2,1)
	      GPtInv(2,2)=S1*GPt(1,1)
	      thx=-(GPtInv(1,1)*dxPt(1)+GPtInv(1,2)*dxPt(2))
	      thy=-(GPtInv(2,1)*dxPt(1)+GPtInv(2,2)*dxPt(2))
	      CALL DSMPROD(DV1,xGrid,thx,3)
	      CALL DSMPROD(DV2,yGrid,thy,3)
	      CALL DADD(th,DV1,DV2,3)
	      CALL DEQUATE(cr2indir,cr1indir,3)
	      CALL DEQUATE(cr2inpos,cr1inpos,3)
	      CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
	      CALL DEQUATE(RayDir,cr2indir,3)
	      CALL DEQUATE(RayPos,cr2inpos,3)
	      CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	       ifLNegOK,Cmatrix,
     &	       iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	       ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	       StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &	       xDraw,yDraw)
	      IF (nBadRays.NE.0) THEN
	        WRITE(*,*)'Chief ray becomes undefined during search'
	        GO TO 385
	      END IF
	      CALL DEQUATE(crPos,RayPos,3)
	      CALL DEQUATE(cr1indir,cr2indir,3)
	      CALL DEQUATE(cr1inpos,cr2inpos,3)
	      CALL DSUB(DV1,RayPos,VptElt(1,iSlv1),3)
	      xInt(1)=DDOTC(xPt,DV1)
	      xInt(2)=DDOTC(yPt,DV1)
	      S1=SQRT(thx*thx+thy*thy)
	      IF ((dthPt.GT.S1).AND.(dthPt.GT.5d-12)) dthPt=0.2*dthPt
    	    End Do
 384	    CONTINUE
	    k=k-1
	    IF (k.EQ.0) THEN
	      WRITE(*,*)' No change in source required'
	    ELSE
	      WRITE(*,308)k,DSQRT(sqErr),
     &	      xInt0(1),xInt(1),xInt0(2),xInt(2),
     &	      ChfRayDir(1),cr2indir(1),ChfRayDir(2),cr2indir(2),
     &	      ChfRayDir(3),cr2indir(3),
     &	      ChfRayPos(1),cr2inpos(1),ChfRayPos(2),cr2inpos(2),
     &	      ChfRayPos(3),cr2inpos(3)
	      CALL CACCEPT(ANS,'YES',
     &        'Accept the new chief ray? ')
	      IF (LCMP(ANS,'YES',1)) THEN
	        CALL DEQUATE(ChfRayPos,cr2inpos,3)
	        CALL DEQUATE(ChfRayDir,cr2indir,3)
	      END IF
	    END IF
 385	    CONTINUE

C         ***Center Beam in System Stop***

	  ELSE IF (LCMP(command,'CENTER',5)) THEN
	    CALL IACCEPT_S(iStop,1,1,
     &	      'Enter system stop element number:')
	    IF ((EltID(iStop).EQ.9).OR.(EltID(iStop).EQ.11)
     &	      .OR.(EltID(iStop).EQ.12)) THEN
	      WRITE(*,*)'Invalid element type =', EltID(iStop)
	    ELSE IF ((iStop.GT.nElt).OR.(iStop.LE.0)) THEN
	      WRITE(*,*)'CENTER: Invalid element number: ', iElt
	    ELSE
	      CALL DEQUATE(DDEF,VptElt,3)
	      CALL DACCEPT(xTarg,DDEF,2,'Enter offset vector (dx,dy):')
	      ifGrid=.FALSE.
	      ifDraw=.FALSE.
	      ifPropagate=.FALSE.
	      ifTrace=.FALSE.
	      ifBuild=.FALSE.
	      ifOPDCalc=.FALSE.

	      CALL STOP(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	        ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	        ifIntSrfSet,ifRayPrint,iObsOpt,iEndElt,iCurWFElt,
     &	        iCurRayElt,Dicr,Djcr,dxSource,
     &	        nBadRays,StartIndRef,StartExtinc,
     &	        Cmatrix,DrawRayVec,DrawEltVec,
     &	        iStop,xTarg,k,sqErr,ifBeamCentered)

	      IF (ifBeamCentered) THEN
	        WRITE(*,309)k,DSQRT(sqErr)
	      END IF

	    END IF

C         ***Compute beam coordinates at an element***

	  ELSE IF (LCMP(command,'COORD',3)) THEN
	    CALL IACCEPT_S(iLocal,1,1,
     &	    'Enter element number:')
	    IF ((EltID(iLocal).EQ.9).OR.(EltID(iLocal).EQ.11)
     &	    .OR.(EltID(iLocal).EQ.12)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iLocal)
	    ELSE IF ((iLocal.GT.nElt).OR.(iLocal.LE.0)) THEN
	      WRITE(*,*)'COORD: Invalid element number: ', iLocal
	    ELSE
	      ifGrid=.FALSE.
	      ifDraw=.FALSE.
	      ifPropagate=.FALSE.
	      ifTrace=.FALSE.
	      ifBuild=.FALSE.
	      ifOPDCalc=.FALSE.

	      CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	      ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &	      ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,
     &	      dxSource,nBadRays,StartIndRef,StartExtinc,
     &	      Cmatrix,DrawRayVec,DrawEltVec,
     &	      iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)
	      IF (ifLocalCoord) THEN
	        WRITE(*,311)xLocal,yLocal,zLocal
	      END IF
	    END IF

	  ELSE IF (LCMP(command,'ACOOR',3)) THEN
	    IF (.not.ifLoad) THEN
              WRITE(*,150)
              GO TO 1
            END IF
            CALL IACCEPT_S(iLocal,1,1,
     &                     'Enter first element number:')
            CALL IACCEPT_S(jLocal,1,1,
     &                     'Enter last element number:')
            ! Check input validity
            IF ((iLocal.GT.nElt).OR.(iLocal.LE.0)) THEN
              WRITE(*,*)'ACOOR: Invalid element number: ', iLocal
              GO TO 1
            ELSE IF ((jLocal.GT.nElt).OR.(jLocal.LE.0)) THEN
              WRITE(*,*)'ACOOR: Invalid element number: ', jLocal
               GO TO 1
            ELSE IF (iLocal>jLocal) THEN
       WRITE(*,*)'ACOOR: last element must be greater than first one'
               GO TO 1
            END IF

	    ! Compute and save
            Open(unit=12,file='xyzLocal.txt',status='replace')
	    Write(*,*) ' Result is saved in file xyzLocal.txt'
            Do ie=iLocal,jLocal
              If ((EltID(ie)/=9).AND.(EltID(ie)/=11)
     &            .AND.(EltID(ie)/=12)) Then
                ifGrid=.FALSE.; ifDraw=.FALSE.
                ifPropagate=.FALSE.; ifTrace=.FALSE.
                ifBuild=.FALSE.; ifOPDCalc=.FALSE.
                !
                CALL LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &           ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &           ifIntSrfSet,
     &           ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &           Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &           Cmatrix,DrawRayVec,DrawEltVec,
     &           ie,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)
                If (.not.ifLocalCoord) Then
          WRITE(*,*)'ACOOR: Compute local coord failed on elmt ',ie
                Else
                  WRITE(12,*) 'iElt =',ie
                  WRITE(12,316) xLocal,yLocal,zLocal
                End If
	      Else
	        WRITE(*,*)'ACOOR: Invalid element type =',EltID(ie)
              End If
            End Do
            Close(12)

C         *** Set option for source adjustment, e.g. in STOP command.
C         *** If this command is not run, the default option will be used.
C         *** Default option:
C         ***   For collimated source, ChfRayDir is fixed.
C         ***   For point source, ChfRayPos is fixed.
	  ELSE IF (LCMP(command,'SAOpt',3)) THEN
	    CALL CACCEPT(ANS,'DIr',
     &              'Choose a source adjustment option, DIr or POs? ')
	    IF (LCMP(ANS,'DI',2) .AND. ABS(zSource).LE.1d10) THEN
       WRITE(*,*)'- Source direction changes when ajusting source'
	      SAOpt=SADir
	    ELSE IF (LCMP(ANS,'PO',2)) THEN
       WRITE(*,*)'- Source position changes when ajusting source'
	      SAOpt=SAPos
	    ELSE
	      If (LCMP(ANS,'DI',2).AND.(ABS(zSource).GT.1d10)) THEN
	       WRITE(*,*)
     &   '** Can not change source direction for collimated source'
	      End If
	      If (SAOpt==SADir) Then
	 WRITE(*,*)'- Source direction changes when ajusting source'
	      Else
	 WRITE(*,*)'- Source position changes when ajusting source'
	      End If
	    END IF

	  ELSE IF (LCMP(command,'SZCO',3)) THEN
            IF (ifLoad) THEN
              CALL CACCEPT(ANS,'Cir',
     &        'Choose Zernike type for ZCO command (Cir,Ann,Hex): ')
              If (LCMP(ANS,'Cir',3)) Then
                zcoType=ZernCircular
                WRITE(*,*)' -- Set to using NormNoll'
              Else If (LCMP(ANS,'Ann',3)) Then
                zcoType=ZernAnnular
                WRITE(*,*)' -- Set to using NormAnnular'
              Else If (LCMP(ANS,'Hex',3)) Then
                zcoType=ZernHex
                WRITE(*,*)' -- Set to using NormHex'
              Else
                WRITE(*,*)' ** Unknown Zernike type'
              End If
            ELSE
              WRITE(*,150)
            END IF

C         *** Set System Stop ***

	  ELSE IF (LCMP(command,'STOP',3)) THEN
	    CALL CACCEPT(ANS,'OBJ','Stop at ELT or OBJect point? ')
	    IF (LCMP(ANS,'OBJ',3)) THEN
	      CALL DZERO(DDEF,7)
	      CALL DACCEPT(StopPos,DDEF,3,
     &	 'Enter stop (entrance pupil) position in object space (x,y,z):')

	      ! Chief ray aiming
	      ! 'Clean up' xGrid, yGrid first, as in source routines
 205	      CALL DXPROD(zGrid,xGrid,yGrid)
              IF (DDOTC(zGrid,ChfRayDir).GT.0d0) THEN
                zGrid = ChfRayDir
              ELSE
                zGrid = -ChfRayDir
              END IF
              CALL DUNITIZE(zGrid)
              CALL DXPROD(yGrid,zGrid,xGrid)
              CALL DUNITIZE(yGrid)
              CALL DXPROD(xGrid,yGrid,zGrid)

	      IF (ABS(zSource).LE.1d10) THEN
	        ! Pt source
	        If (SAOpt==SADir) Then
	          cr1indir(1:3)=ChfRayDir(1:3)
	          !ChfRayDir(1:3)=StopPos(1:3)-ChfRayPos(1:3)
	          srcPos(1:3)=ChfRayPos(1:3)+zSource*ChfRayDir
	          ChfRayDir(1:3)=StopPos(1:3)-srcPos
	          CALL DUNITIZE(ChfRayDir)
	          c1=DDOTC(ChfRayDir,ChfRayDir)
	          If (c1.LE.1d-22) Then
                    WRITE(*,*)' Error in setting STOP, abort ...'
	            ChfRayDir(1:3)=cr1indir(1:3)
	          Else
	            ! Need update ChfRayPos as well
	            ChfRayPos(1:3)=srcPos(1:3)-zSource*ChfRayDir(1:3)
		  End If
	        Else If (SAOpt==SAPos) Then
	          ! Translate chief ray in (xGrid,yGrid) plane
	          x=DDOTC(StopPos(1:3)-ChfRayPos(1:3),xGrid)
                  y=DDOTC(StopPos(1:3)-ChfRayPos(1:3),yGrid)
            ChfRayPos(1:3)=ChfRayPos(1:3)+x*xGrid(1:3)+y*yGrid(1:3)
	        End If
	      ELSE
	        ! Collimated source, translate chief ray in (xGrid,yGrid) plane
	        x=DDOTC(StopPos(1:3)-ChfRayPos(1:3),xGrid)
	        y=DDOTC(StopPos(1:3)-ChfRayPos(1:3),yGrid)
           ChfRayPos(1:3)=ChfRayPos(1:3)+x*xGrid(1:3)+y*yGrid(1:3)
	      END IF
	      ifStopSet=.TRUE.

	    ELSE
	      ! Stop specified at ELT point
	      CALL IACCEPT_S(iStop,1,1,
     &	                     'Enter system stop element number:')
	      IF ((EltID(iStop).EQ.9).OR.(EltID(iStop).EQ.11)
     &	          .OR.(EltID(iStop).EQ.12)) THEN
	        WRITE(*,*)'Invalid element type =',EltID(iStop)
	      ELSE IF ((iStop.GT.nElt).OR.(iStop.LE.0)) THEN
	        WRITE(*,*)'STOP: Invalid element number: ',iStop
	      ELSE
	        CALL DEQUATE(DDEF,VptElt(1:2,iStop),3)
	        CALL DACCEPT(xTarg,DDEF,2,
     &	        	     'Enter offset vector (dx,dy):')

 210	        If (RxStopSet .and. EltStopSet) Then
	          iStop=StopElt
	          xTarg(1:2)=StopOffset(1:2)
	          RxStopSet=.false.
	        End If

	        ifGrid=.FALSE.
	        ifDraw=.FALSE.
	        ifPropagate=.FALSE.
	        ifTrace=.FALSE.
	        ifBuild=.FALSE.
	        ifOPDCalc=.FALSE.

                IF (ABS(zSource).LE.1d10 .AND. SAOpt==SADir) THEN
		  ! With point source (source at finite distance) and
		  ! when SAOpt is set to SAPos, enforce ChfRay
                  ! thru STOP position at element by tip/tilt
                  ! chief ray in source plane. -jzlou
	          CALL SetStop(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &              ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &              ifIntSrfSet,ifRayPrint,iObsOpt,iEndElt,iCurWFElt,
     &              iCurRayElt,Dicr,Djcr,dxSource,
     &              nBadRays,StartIndRef,StartExtinc,
     &              Cmatrix,DrawRayVec,DrawEltVec,
     &              iStop,TElt(1:3,1,iStop),TElt(1:3,2,iStop),
     &              xTarg,sqErr,ifBeamCentered)

	          !WRITE(*,310) StopPos
                   ! The "Computed StopPos" is actually system entrance pupil
                  ifStopSet=.TRUE.
	          GO TO 1
	        END IF

	        ! Otherwise, for source at infinity OR finite source with
                ! SAOpt set to SAPos,
                ! Center beam in stop to get new chief ray:
	        CALL STOP(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	            ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	            ifIntSrfSet,ifRayPrint,iObsOpt,iEndElt,iCurWFElt,
     &	            iCurRayElt,Dicr,Djcr,dxSource,
     &	            nBadRays,StartIndRef,StartExtinc,
     &	            Cmatrix,DrawRayVec,DrawEltVec,
     &	            iStop,xTarg,k,sqErr,ifBeamCentered)

	        CALL DEQUATE(cr1indir,ChfRayDir,3)
	        CALL DEQUATE(cr1inpos,ChfRayPos,3)

                ! Starting guess for StopPos is ray at first surface:
	        iStartRay=1
	        iEndRay=1
	        iStartElt=0
	        iEndElt=1
	        nMinPts=3
	        IF (ABS(zSource).GE.1d10) THEN
	          CALL ColSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	          ifPropagate)
	        ELSE
	          CALL PtSource(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &	          ifPropagate)
	        END IF
	        CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	          ifLNegOK,Cmatrix,
     &	          iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &	          ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &	          StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &            xDraw,yDraw)
	        CALL DEQUATE(StopPos,RayPos,3)
	          ! After ray 1 hit first surface, that position is designated
	          ! as first guess of StopPos. -jzlou

C         Compute new differential chief ray:
c	        th(1)=3.5355d-6
c	        th(2)=3.5355d-6
c	        th(3)=3.5355d-6
	        CALL DSMPROD(th,xGrid,5d-6,3)
	        CALL RotChfRay(th,ChfRayDir,ChfRayPos,StopPos)

	        !print*,' **B4 2nd STOP , yGrid =',yGrid

	        ! When differential chief ray is used to set up source
	        ! xGrid and ygrid will be recomputed, so save the existing
                ! ones for 'main chief ray'. Otherwise it would be at least
	        ! confusing to the user when running 'show 0' command
	        xGrid_sv=xGrid; yGrid_sv=yGrid

	        CALL STOP(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	          ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &	          ifIntSrfSet,ifRayPrint,iObsOpt,iEndElt,iCurWFElt,
     &	          iCurRayElt,Dicr,Djcr,dxSource,
     &	          nBadRays,StartIndRef,StartExtinc,
     &	          Cmatrix,DrawRayVec,DrawEltVec,
     &	          iStop,xTarg,k,sqErr,ifBeamCentered)
	        xGrid=xGrid_sv; yGrid=yGrid_sv

	        !print*,' **Aft 2nd STOP , yGrid =',yGrid

C         Find crossing point of the 2 chief rays:

	        CALL FindCrossPt(cr1inpos,cr1indir,
     &	               ChfRayPos,ChfRayDir,zp,CrossDist,StopPos)

	        WRITE(*,310) StopPos
	        WRITE(*,*)' '
	          ! The "Computed" StopPos is actually the entrance pupil,
                  ! which is different from the user specified position
	          ! on the element where system STOP is set. -jzlou
	        ifStopSet=.TRUE.
	        CALL DEQUATE(ChfRayDir,cr1indir,3)
	        CALL DEQUATE(ChfRayPos,cr1inpos,3)
	      END IF
	    END IF


#ifdef DESIGN_OPTIM

C	  *** System Optimization Commands ***

	  ! This command lists all optical elements in Rx, which is not
	  ! optimization specific.
	  ELSE IF (LCMP(command,'ELTS',4)) THEN
	    IF (ifLoad) THEN
	       WRITE(*,*)' '
	       WRITE(*,*)' Elements in optical prescription:'
	       Do iElt=1,nElt
	         WRITE(*,*)' Element name = ', EltName(iElt),
     &             'ID =', iElt
	       End Do
	    ELSE
              WRITE(*,150)
            END IF

C         *** Add a variable element, supported in both MACOS and SMACOS

	  ELSE IF (LCMP(command,'AVAR',4)) THEN
	    IF (ifLoad) THEN
 21	      CALL IACCEPT_S(iElt,iCurWFElt,1,
     &               'Enter variable element id for optimization:')
	      IF ((iElt.LE.0).OR.(iElt.GT.nElt)) THEN
                WRITE (*,*)'AVAR: Invalid element number: ',iElt
#ifdef MACOS_CMD
                GO TO 21
#else
	        WRITE(*,*)'SMACOS quit ...'
	        stop
#endif

C	      ELSE IF (EltId(iElt)==3) THEN
C	        WRITE(*,*)' '
C	        WRITE (*,*)
C    &           'Variable element can not be a reference surface!'
C	        GO TO 24
	      ELSE
	        If (.not. isVarElt(iElt)) then
	          isVarElt(iElt)=.true.
	          nVarElt=nVarElt+1
	          varElts(nVarElt)=iElt
	          CALL IACCEPT(varEltDOF(1:mVarDOF,iElt),
     &			   varEltDOF(1:mVarDOF,iElt),
     &                     mVarDOF,'Enter variable degree of freedom:')
	          nDOF_VarElt(nVarElt)=0
                  Do i=1,mVarDOF
	            if (varEltDOF(i,iElt)/=0) then
                      nDOF_VarElt(nVarElt)=nDOF_VarElt(nVarElt)+1
                      DOF_VarElt(nDOF_VarElt(nVarElt),nVarElt)
     &                  = DOF_NameList(i)
	            end if
                  End Do

	          ! Add Zernike terms (range) for optimization
	          zernRange(1)=1; zernRange(2)=45
 156	          CALL IACCEPT(zernRange,zernRange,2,
     &			       'Enter range of Zernike terms:')
	          if (zernRange(1)>0 .and. zernRange(1)<46 .and.
     &                zernRange(2)>0 .and. zernRange(2)<46 .and.
     &                zernRange(1)<zernRange(2)) then
	            nOptEltZern(iElt)=zernRange(2)-zernRange(1)+1
	            Do i=zernRange(1),zernRange(2)
	              OptEltZernTerm(i-zernRange(1)+1,iElt)=i
	            End Do
	          else
#ifdef MACOS_CMD
	            if (zernRange(1)==0 .and. zernRange(2)==0) then
	              nOptEltZern(iElt)=0  ! no Zernike optimization
	            else
	              WRITE(*,*)' Invalid range for Zernike terms!'
	              zernRange(1)=1; zernRange(2)=45
	              go to 156
	            end if
#else
	            ! SMACOS, and no Zernike optimization
	            nOptEltZern(iElt)=0
#endif
	          end if
	        Else
	          WRITE(*,*)' Variable element',iElt,'already defined!'
	        End If
              END IF
	    ELSE
              WRITE(*,150)
            END IF
 24	    CONTINUE

C	  *** Modify a defined variable element
	  ELSE IF (LCMP(command,'MVAR',4)) THEN
	    IF (ifLoad) THEN
 25           CALL IACCEPT_S(iElt,iCurWFElt,-1,
     &               'Enter variable element ID for optimization:')
              IF ((iElt.LT.0).OR.(iElt.GT.nElt)) THEN
                WRITE (*,*)' MVAR: Invalid element ID: ',iElt
                GO TO 25
	      ELSE
	        ! Search
	        if (.not.isVarElt(iElt)) then
	      	  WRITE(*,*)' Element',iElt,'is not a variable'
	          GO TO 27
	        end if
	        Do iVarElt=1,nVarElt
	          if (varElts(iVarElt)==iElt) exit
	        End Do
	        CALL IACCEPT(varEltDOF(1:mVarDOF,iElt),
     &                     varEltDOF(1:mVarDOF,iElt),
     &                     mVarDOF,'Enter variable degree of freedom:')
	        nDOF_VarElt(iVarElt)=0
                Do i=1,mVarDOF
                  if (varEltDOF(i,iElt)/=0) then
                    nDOF_VarElt(iVarElt)=nDOF_VarElt(iVarElt)+1
                    DOF_VarElt(nDOF_VarElt(iVarElt),iVarElt)
     &                = DOF_NameList(i)
                  end if
                End Do
	      END IF
	    ELSE
              WRITE(*,150)
            END IF
 27	    CONTINUE

C         *** Delete a variable element, supported in both MACOS and
C             SMACOS
	  ELSE IF (LCMP(command,'DVAR',4)) THEN
	    IF (ifLoad) THEN
 23	      CALL IACCEPT_S(iElt,iCurWFElt,-1,
     &              'Enter element ID to remove from optimization:')
#ifdef MACOS_CMD
              IF ((iElt.LT.0).OR.(iElt.GT.nElt)) THEN
	        WRITE (*,*)'DVAR: Invalid element ID: ',iElt
                GO TO 23
	      END IF
#endif
	      IF ((iElt.GE.0) .AND. (iElt.LE.nElt)) THEN
	        If (isVarElt(iElt)) then
	          isVarElt(iElt)=.false.
                  nVarElt=nVarElt-1
	          varEltDOF(1:mVarDOF,iElt)=0
	          nOptEltZern(iElt)=0
	        Else
	          WRITE(*,*)'Element ',iElt,
     &                      ' is not a variable element!'
	        End If
	      END IF

	    ELSE
              WRITE(*,150)
            END IF

	  ELSE IF (LCMP(command,'AFOV',4)) THEN
	    IF (ifLoad) THEN
	      WRITE(*,*)'Add FOV for design optimization ...'
	      CALL DACCEPT(fov_ChfRayDir,
     &               fov_ChfRayDir,-3,
     &               'Enter Chief Ray Direction')
	      CALL DACCEPT(fov_ChfRayPos,
     &               fov_ChfRayDir,-3,
     &               'Enter Chief Ray Position')
	      nOptFov=nOptFov+1
	      opt_fov(1:3,1,nOptFov)=fov_ChfRayDir
	      opt_fov(1:3,2,nOptFov)=fov_ChfRayPos
	    ELSE
              WRITE(*,150)
            END IF

	  ELSE IF (LCMP(command,'DFOV',4)) THEN
            IF (ifLoad) THEN
	      CALL IACCEPT_S(jfov,
     &               jfov,-1,
     &               'Enter ID of FOV to delete')
	      If (jfov<1 .or. jfov>nOptFov) Then
	        WRITE(*,*)'DFOV: Invalid FOV ID'
	      Else
	        Do ifov=1,nOptFov
	          if (ifov==jfov) then
                    if (jfov==nOptFov) then
	              nOptFov=nOptFov-1
	            else
                      Do i=ifov,nOptFov-1
	                opt_fov(1:3,1:2,i)=opt_fov(1:3,1:2,i+1)
	              End Do
	              nOptFov=nOptFov-1
	            end if
	          end if
	        End Do
	      End If
	    ELSE
              WRITE(*,150)
            END IF

	  ELSE IF (LCMP(command,'FOVS',4)) THEN
	    ! List all FOVs to be used for design optimization
	    IF (ifLoad) THEN
	      WRITE(*,*)' '
	      WRITE(*,*)
     &        ' Number of field of view used for design optimization =',
     &        nOptFov
	      Do ifov=1,nOptFov
	        WRITE(*,*)' FOV ID =',ifov,', Weight =',fov_weight(ifov)
	        WRITE(*,*)'    ChfRayPos =',
     &            '[',opt_fov(1:3,1,ifov),']'
	        WRITE(*,*)'    ChfRayDir =',
     &            '[',opt_fov(1:3,2,ifov),']'
	      End Do
	      WRITE(*,*)' '
	    ELSE
              WRITE(*,150)
            END IF

	  ELSE IF (LCMP(command,'WLENS',5)) THEN
	    ! List all wavelengths to be used for design optimization
            IF (ifLoad) THEN
              WRITE(*,*)' '
              WRITE(*,*)
     &        ' Number of wavelengths used for design optimization =',
     &        nOptWavelen
	      Do iwl=1,nOptWavelen
	        WRITE(*,151)' Wavelength',iwl,' =',opt_wavelen(iwl)
	      End Do
	      WRITE(*,*)' '
	    ELSE
              WRITE(*,150)
            END IF

	  ELSE IF (LCMP(command,'VARS',4)) THEN
	    IF (ifLoad) THEN
	      if (nVarElt>0) then
	        WRITE(*,*)' '
	        WRITE(*,202)
     &            'Variable elements defined for optimization: '
	        Do iElt=1,nVarElt
	          nID=1
                  varEltID(nID)=varElts(iElt)
                  Do i=1, nElt
                    if (LnkElt(i)==varElts(iElt)) then
                      nID=nID+1
                      varEltID(nID)=i
                    end if
                  End Do
                  slen=StrLen(EltName(varElts(iElt)))
                  WRITE(*,371)
     &              ' Elt Name= ', EltName(varElts(iElt))(1:slen)
	          WRITE(*,372) ', Elt ID='
	          Do i=1,nID
                    valbuf(1:len(valbuf))=' '
                    WRITE(valbuf,373) varEltID(i)
                    valbuf=ADJUSTL(valbuf)
                    if (i<nID) then
                      WRITE(*,372) valbuf(1:strlen(valbuf)+1)
                    else
                      WRITE(*,372) valbuf(1:strlen(valbuf))
                    end if
                  End Do
	          if (nDOF_VarElt(iElt)>0) then
	            WRITE(*,374)', Variable DOF= '
                    Do i=1,nDOF_VarElt(iElt)
                      WRITE(*,372)
     &                 DOF_VarElt(i,iElt)(1:strlen(DOF_VarElt(i,iElt)))
                      WRITE(*,372) ' '
                    End Do
	          end if
	          if (nOptEltZern(varElts(iElt))>0) then
	            WRITE(*,374)'  Zern Terms='
	            Do i=1,nOptEltZern(varElts(iElt))
	              WRITE(*,375) OptEltZernTerm(i,varElts(iElt))
	            End Do
	          end if
	          WRITE(*,*)' '
	        End Do
	        WRITE(*,*)' '
	      else
	        WRITE(*,*)' '
                WRITE(*,*)' No Variable Element for Optimization!'
	        WRITE(*,*)' '
	      end if
	    ELSE
              WRITE(*,150)
            END IF
 375        FORMAT(I3,$)


C       *** Set parameters for design optimization

	  ELSE IF (LCMP(command,'SETC',4).OR.
     &             LCMP(command,'SETP',4)) THEN
	    IF (ifLoad) THEN
#ifdef MACOS_CMD
	      WRITE(*,*)' '
	      WRITE(*,*)
     &          'Setting options for MACOS design optimization ...'
	      WRITE(*,*)' '

	      CALL CACCEPT(ANS,'WFE',
     &             'Set target of design optimization (WFE,ZWF,Beam):')
	      If (LCMP(ANS,'WFE',3)) Then
	        OptTarget=WFE_TARGET
	      Else If (LCMP(ANS,'ZWF',3)) Then
	        OptTarget=WFE_ZMODE_TARGET
	        CALL IACCEPT_S(nWFZern,0,-1,
     &                         'Enter number of Zernike modes:')
	        WFZernMode=0
	        CALL IACCEPT(WFZernMode,WFZernMode,-nWFZern,
     &                      'Enter Zernike modes')
	      Else If (LCMP(ANS,'Beam',4)) Then
	        OptTarget=BEAM_TARGET
	      Else
	        OptTarget=WFE_TARGET  ! default
	      End If

	      CALL IACCEPT_S(nitrs_dopt,nitrs_dopt,1,
     &          'Set maximum iterations for optimization:')
	      CALL DACCEPT_S(dopt_tol,dopt_tol,1,
     &          'Set convergence tolerance for optimization:')

	      CALL IACCEPT_S(i,nOptFov,1,
     &          'Set number of field of view:')
              if (i .LE. nOptFov) then
                nOptFov=i
              else
                WRITE(*,*)' Set number of field of view to',nOptFov
              end if

	      CALL DACCEPT(fov_weight(1:nOptFov),
     &                     fov_weight(1:nOptFov),nOptFov,
     &          	   'Set field of view weights:')

	      CALL CACCEPT(ANS,'no',
     &         'Set target wavefront file for optimization (no/yes):')
	      if (LCMP(ANS,'Y',1)) then
	        CALL CACCEPT(tgtWF_fmt,'bin',
     &         'Enter file format for target wavefront (bin/ascii):')
	        CALL CACCEPT(tgtWF_fn,' ',
     &               'Enter file name for target wavefront: ')
	        tgtWF_flg=.true.
	      else
	        tgtWF_flg=.false.
	      end if

#else
	  ! SMACOS portion
	  ! Call SMACOS_SetOptParam(pname,carr,iarr,darr) to set
	  ! individual parameters - see dopt_mod.F
#endif

	    ELSE
	      ! not ifLoad
              WRITE(*,150)
            END IF


C         *** System optimization command, supported in both MACOS and SMACOS ***

	  ELSE IF (LCMP(command,'CALIB',3)) THEN

	    IF (ifLoad) THEN

#ifdef MACOS_CMD
	! This part is to increase 'calib' command's flexibility with
        ! optional arguments - for command-line macos only
	      LCmdCalibElts=.FALSE.
              i=read_len(pstack)-1  ! # arguments to calib command
              If (i.GE.1) Then
                LCmdCalibElts=.True.
                isVarElt_Rx=isVarElt; nVarElt_Rx=nVarElt
                isVarElt=.FALSE.; nVarElt=0
                Do j=1,i
                  CALL IACCEPT_S(iElt,iElt,-1,
     &                  'Enter element to be calibrated:')
                  If (isVarElt_Rx(iElt)) Then
                    isVarElt(iElt)=.TRUE.; nVarElt=nVarElt+1
                  Else
                    WRITE(*,*)' ** Calib error: Elt',iElt,
     &               ' is not variable; optimization aborted!'
                    GO TO 211
                  End If
                End Do
              End If
#endif

	      if (nVarElt==0) then
	        WRITE(*,*)' '
                WRITE(*,*)' No Variable Element for Optimization!'
                goto 105
	      else if (nOptFov==0) then
	        WRITE(*,*)' '
                WRITE(*,*)' No Field of View for Optimization!'
                goto 105
              end if

	      if (LOptIfFEX .and. .not.ifStopSet) then
	        WRITE(*,*)' Must set STOP before optimization with FEX!'
	        WRITE(*,*)' '
	        goto 105
	      end if

	      ! Process input arguments to system optimzer drivers
	      varZerns=0; zern_off=1
	      varAsphs=0; asph_off=1
              k=1
              iVarElt=0

	      ! Process source first
	      If (isVarElt(0)) Then
	        iVarElt=iVarElt+1
                varElts(iVarElt)=0
	        varDOFs((iVarElt-1)*mVarDOF+1:iVarElt*mVarDOF)
     &               =varEltDOF(1:mVarDOF,0)
	        If (.not.ifStopSet .and. .not.OptSrcRpt_FLG) Then
	          if (varEltDOF(1,0)/=0 .or. varEltDOF(2,0)/=0
     &                .or. varEltDOF(3,0)/=0) then
	            WRITE(*,*)' Must set STOP before optimizing source!'
	            go to 105
	          end if
	        End If
	        EltName(0)='Source'
	        nOptFov=1
	        opt_fov(1:3,1,1)=ChfRayDir; opt_fov(1:3,2,1)=ChfRayPos
	        zernElt(iVarElt)=0; asphElt(iVarElt)=0
	      End If

	      Do iElt=1,nElt
		if (isVarElt(iElt)) then
	          iVarElt=iVarElt+1
	          varElts(iVarElt)=iElt
	          varDOFs((iVarElt-1)*mVarDOF+1:iVarElt*mVarDOF)
     &               =varEltDOF(1:mVarDOF,iElt)

	          zernElt(iVarElt)=0  ! range 1:nVarElt
	          if (nOptEltZern(iElt)>0) then
	            zernElt(iVarElt)=nOptEltZern(iElt)
	            varZerns(zern_off:zern_off+nOptEltZern(iElt)-1)
     &	              = OptEltZernTerm(1:nOptEltZern(iElt),iElt)
	            zern_off=zern_off+nOptEltZern(iElt)
	          end if

		  asphElt(iVarElt)=0  ! range 1:nVarElt
                  if (nOptEltAsph(iElt)>0) then
                    asphElt(iVarElt)=nOptEltAsph(iElt)
                    varAsphs(asph_off:asph_off+nOptEltAsph(iElt)-1)
     &                = OptEltAsphTerm(1:nOptEltAsph(iElt),iElt)
                    asph_off=asph_off+nOptEltAsph(iElt)
                  end if

	          If (LOptCons) Then
                    ! process variable constraints
                    Do i=1,mVarDOF
                      if (varEltDOF(i,iElt)/=0) then
                        varEltCons(1:2,k)=EltCons(1:2,i,iElt)
                        k=k+1
                      end if
                    End Do
                  End If

	          ! Save state of variable element
	          psiElt_old(1:3,iVarElt)=psiElt(1:3,iElt)
	          vptElt_old(1:3,iVarElt)=vptElt(1:3,iElt)
	          KrElt_old(iVarElt)=KrElt(iElt)
	          KcElt_old(iVarElt)=KcElt(iElt)

#if 0
	           nDOF_VarElt(iVarElt)=0
	           Do i=1,7
	             if (varEltDOF(i,iElt)/=0) then
	 	       nDOF_VarElt(iVarElt)=nDOF_VarElt(iVarElt)+1
	               DOF_VarElt(nDOF_VarElt(iVarElt),iVarElt)
     &                   = DOF_NameList(i)
	             end if
	           End Do
#endif
	        end if
	      End Do

	      ! Save Rx Chief Ray
	      ChfRayDir_Rx=ChfRayDir
	      ChfRayPos_Rx=ChfRayPos

	      ! Save Rx npts, use opt_npts for design optimization
	      npts_Rx=npts; npts=opt_npts

	      !print*,'B4 calling nls_optim_dvr ...'
	      !print*,'nOptEltZern= ',nOptEltZern(1:nElt)
	      !print*,'varZerns= ',varZerns

	      If (OptTarget==BEAM_TARGET .and.
     &          beamOpt%ifOptBeamRefRayDir) Then
	        ! Find number of ref rays and their ray ids
	        iEndElt=beamOpt%beamRefRayElt
	        Call getRefRays(nRefRay,refRay,npts,ifPol,ifGrid,
     &                 ifPropagate,ifDraw,ifLNegOK,ifBuild,ifDMBuild,
     &                 ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,iObsOpt,
     &                 iCurWFElt,iCurRayElt,iEndElt,Dicr,Djcr,
     &                 dxSource,nBadRays,StartIndRef,StartExtinc,
     &                 Cmatrix,DrawRayVec,DrawEltVec)
	      End If

	      If (.not. LOptCons) Then
	        print*,'*** Before nls_driver: ifStopSet =',ifStopSet
                ! Unconstrained system optimization
	        opt_fov(1:3,1,1)=ChfRayDir; opt_fov(1:3,2,1)=ChfRayPos
                Call nls_optim_dvr(.true.,OptTarget,
     &                             nVarElt,varElts,varDOFs,
     &                             nWFZern,WFZernMode,
     &                             zernElt(1:nVarElt),varZerns,
     &                             asphElt(1:nVarElt),varAsphs,
     &                             opt_fov,nOptFov,fov_weight,
     &                             opt_wavelen,nOptWavelen,
     &                             dopt_tol,nitrs_dopt,
     &                             old_wfe,new_wfe,rtn_flag,
     &                             tgtWF_flg,tgtWF_fn,tgtWF_fmt,
     &                             OptTgtElt,LOptIfFEX,beamOpt)
	        print*,'*** After nls_driver: ifStopSet =',ifStopSet
	      Else
                ! Constrained system optimization
                Call np_optim_dvr(.true.,OptTarget,
     &                             nVarElt,varElts,varDOFs,
     &                             varEltCons,
     &                             nWFZern,WFZernMode,
     &                             zernElt(1:nVarElt),varZerns,
     &                             opt_fov,nOptFov,fov_weight,
     &                             opt_wavelen,nOptWavelen,
     &                             dopt_tol,nitrs_dopt,
     &                             old_wfe,new_wfe,rtn_flag,
     &                             tgtWF_flg,tgtWF_fn,tgtWF_fmt,
     &                             OptTgtElt,LOptIfFEX,beamOpt)
              End If

	      npts=npts_Rx; nGridPts=npts+1
	      If (rtn_flag /= 0) Then
	        goto 105
	      Else
	        If (isVarElt(0)) Then
		  ! set STOP
!#include "srcaim.inc"
	        End If
              End If

#ifdef MACOS_CMD
	      CALL CACCEPT(ANS,'yes',
     &               'Display optimization summary (yes,no):')
#else
	      ANS(1:3)='yes'  ! default to 'yes' for SMACOS
#endif

              IF (LCMP(ANS,'yes',1)) THEN
	        WRITE(*,*)'  '
	   WRITE(*,*)'    ==============================================='
	   WRITE(*,*)'    ***    MACOS System Optimization Summary    ***'
	   WRITE(*,*)'    ==============================================='
	        WRITE(*,*)'  '
	        If (.not.LOptCons) Then
                  if (OptAlg==1) then
               WRITE(*,*)' Optimization Scheme: Linear Unconstrained'
                  else
               WRITE(*,*)' Optimization Scheme: Nonlinear Unconstrained'
                  end if
                Else
               WRITE(*,*)' Optimization Scheme: Nonlinear Constrained'
                End If

	        If (OptTarget==WFE_TARGET) Then
	          WRITE(*,*)' Optimization Target = ',
     &		          DOPT_TARGET(OptTarget)(1:19),' at Element ',
     &                    OptTgtElt
	        Else If (OptTarget==WFE_ZMODE_TARGET) Then
	          WRITE(*,*)' Optimization Target = ',
     &                    DOPT_TARGET(OptTarget)(1:17),' at Element ',
     &                    OptTgtElt
	        Else If (OptTarget==SPOT_TARGET) Then
                  WRITE(*,*)' Optimization Target: ',
     &                    DOPT_TARGET(OptTarget)(1:14),'at Element ',
     &                    OptTgtElt
	        Else If (OptTarget==Beam_TARGET) Then
                  WRITE(*,*)' Optimization Target = ',
     &                    DOPT_TARGET(OptTarget)(1:17),' at Element ',
     &                    OptTgtElt
	        End If
	        WRITE(*,*)' Number of Rays Used          = ', nRay
	        WRITE(*,*)' Number of Field of View Used = ', nOptFov
	        WRITE(*,*)' Number of Wavelength Used    = ', nOptWavelen
                Do iwl=1,nOptWavelen
	          WRITE(*,198)'  Wavelength = ', opt_wavelen(iwl),
     &            '  ',WaveUnits
 198    	  Format(1P,A15,D14.6,A2,A6)
	          if (OptTarget==WFE_TARGET) then
	            WRITE(*,*)
     &    	     ' * RMS WFE at FOVs BEFORE Optimization:',
     &               old_wfe(1:nOptFov,iwl)
	            WRITE(*,*)
     &               ' * RMS WFE at FOVs AFTER  Optimization:',
     &               new_wfe(1:nOptFov,iwl)
	          else if (OptTarget==SPOT_TARGET) then
	            WRITE(*,*)
     &               ' * Spot size at FOVs BEFORE Optimization:',
     &               SpotSize_old
	            WRITE(*,*)
     &               ' * Spot size at FOVs AFTER Optimization:',
     &               SpotSize_new
	          end if
	        End Do

	        WRITE(*,202)'Variable Elements for This Optimization: '
		aparam_off=0
	        Do iElt=1,nVarElt
	          nID=1
	          varEltID(nID)=varElts(iElt)
	          Do i=1,nElt
                    if (LnkElt(i)==varElts(iElt)) then
     	              nID=nID+1
	              varEltID(nID)=i
	            end if
                  End Do
	          slen=StrLen(EltName(varElts(iElt)))

C	          WRITE(*,*)' Elt Name= ', EltName(varElts(iElt))(1:slen),
C    &                      ', Elt ID=', varEltID(1:nID),
C    &              ',  Variable DOF= ',
C    &              DOF_VarElt(1:nDOF_VarElt(iElt),iElt)

	          ! Very detailed adjustment of output format
	          WRITE(*,371)' Elt Name = ',EltName(varElts(iElt))(1:slen)
 371              FORMAT(3x,2A,$)
                  WRITE(*,372) ', Elt ID = '
 372              FORMAT(A,$)
	          Do i=1,nID
	            valbuf(1:len(valbuf))=' '
	            WRITE(valbuf,373) varEltID(i)
	            valbuf=ADJUSTL(valbuf)
	            if (i<nID) then
	              WRITE(*,372) valbuf(1:strlen(valbuf)+1)
	            else
	              WRITE(*,372) valbuf(1:strlen(valbuf))
	            end if
	            !WRITE(*,373) varEltID(i)
	          End Do
 373              FORMAT(I4)

	          If (nDOF_VarElt(iElt)>0 .OR. zernElt(iElt)>0
     &                .OR. asphElt(iElt)>0) Then
                    WRITE(*,374)', Variable DOF = '
                    Do i=1,nDOF_VarElt(iElt)
                      WRITE(*,372)
     &                DOF_VarElt(i,iElt)(1:strlen(DOF_VarElt(i,iElt)))
                      if (i<nDOF_VarElt(iElt)) WRITE(*,372) ' '
                    End Do
                    If (zernElt(iElt)>0) Then
                      If (nDOF_VarElt(iElt)>0) WRITE(*,374)','
                      WRITE(*,374)' Surface Zernike modes'
                    End If
	            If (asphElt(iElt)>0) Then
                      If (nDOF_VarElt(iElt)>0) WRITE(*,374)','
                      WRITE(*,374)' Surface aspheric coeffs'
                    End If
                    WRITE(*,*)' '
                  Else
                    WRITE(*,*)', Variable DOF = None'
                  End If
 374              FORMAT(A,$)

	          If (varElts(iElt)==0) Then ! source
                    WRITE(*,201)'   Old ChfRayDir =',ChfRayDir_Rx
                    WRITE(*,201)'   Old ChfRayPos =',ChfRayPos_Rx
                    WRITE(*,201)'   New ChfRayDir =',ChfRayDir
                    WRITE(*,201)'   New ChfRayPos =',ChfRayPos
	          Else
	            WRITE(*,201)'   Old PsiElt =',psiElt_old(1:3,iElt)
	            WRITE(*,201)'   Old VptElt =',vptElt_old(1:3,iElt)
	            WRITE(*,201)'   Old KrElt  =', KrElt_old(iElt)
	            WRITE(*,201)'   Old KcElt  =', KcElt_old(iElt)
		    If (zernElt(iElt)>0) Then
                      WRITE(*,*)'     Old ZernCoef ='
                      WRITE(*,204) ZernCoef_old(1:6,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(7:12,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(13:18,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(19:24,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(25:30,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(31:36,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(37:42,varElts(iElt))
                      WRITE(*,204) ZernCoef_old(43:45,varElts(iElt))
                    End If
                    If (asphElt(iElt)>0) Then
                      WRITE(*,376)'     Old AsphCoef ='
                      WRITE(*,*) AsphCoef_old(1:3,varElts(iElt))
                      if (maxAsphCoef>3) then
                        WRITE(*,*)'       ',
     &                   AsphCoef_old(4:maxAsphCoef,varElts(iElt))
                      end if
                    End If
 376                FORMAT(A,$)

                WRITE(*,201)'   New PsiElt =',psiElt(1:3,varElts(iElt))
                WRITE(*,201)'   New VptElt =',vptElt(1:3,varElts(iElt))
                WRITE(*,201)'   New KrElt  =', KrElt(varElts(iElt))
	        WRITE(*,201)'   New KcElt  =', KcElt(varElts(iElt))
	            j=0
	            Do i=1,mVarDOF
	              if (varEltDOF(i,varElts(iElt))>0) then
	                j=1; exit
	              end if
	            End DO
	            If (j>0) Then
                      ! there were DOF to adjust in VarElt
	              If (nECoord(varElts(iElt))>0) Then
	                WRITE(*,164)
     &                   '   Variable DOF Adjustments in Local Frame: '
		      Else
	                WRITE(*,164)
     &                   '   Variable DOF Adjustments in Global Frame: '
	              End If
	            End If
		    If (.not. LOptCons .and. j>0) Then
	              j=0; WRITE(*,165)'     '
		      Do i=1,mVarDOF
		        if (varEltDOF(i,varElts(iElt))>0) then
		 	  j=j+1
                          WRITE(*,166) aparams_nls_m(aparam_off+j)
	                end if
		      End Do
		      WRITE(*,'(x)')
		      aparam_off=aparam_off+j
		    Else
	              ! for constrained optimizer ...
		    End If

	            If (zernElt(iElt)>0) Then
                      WRITE(*,*)'    New Zernike coef='
                      WRITE(*,204) ZernCoef(1:6,varElts(iElt))
                      WRITE(*,204) ZernCoef(7:12,varElts(iElt))
                      WRITE(*,204) ZernCoef(13:18,varElts(iElt))
                      WRITE(*,204) ZernCoef(19:24,varElts(iElt))
                      WRITE(*,204) ZernCoef(25:30,varElts(iElt))
                      WRITE(*,204) ZernCoef(31:36,varElts(iElt))
                      WRITE(*,204) ZernCoef(37:42,varElts(iElt))
                      WRITE(*,204) ZernCoef(43:45,varElts(iElt))
                    End If
	            !
	            If (asphElt(iElt)>0) Then
                      WRITE(*,376)'     New AsphCoef ='
                      WRITE(*,*) AsphCoef(1:3,varElts(iElt))
                      If (maxAsphCoef>3) Then
                        WRITE(*,*)'       ',
     &                    AsphCoef(4:maxAsphCoef,varElts(iElt))
                      End If
                    End If
	          End If ! not source
	        End Do
	        WRITE(*,*)'  '
	      END IF

#ifdef MACOS_CMD
	      ! Should output 'aparam' vector from nls_optim_dvr routine,
	      ! which provides total adjustments of calibrated parameters
	      ! from initial state, and thus provides rotations and
	      ! translations of each calibrated elements. So the following
	      ! PrtPrinRots calculation is not necessary. -jzlou, 01/2008

	      If (PrtPrinRots) Then
	        ! Calculate three PRINCIPAL rotations needed to go from initial
	        ! element orientations to optimized element orientations
	        WRITE(*,*)
     &           ' Calculated Element Principal Rotations (Z-Y-X order)'
	        Do i=1,nVarElt
	          iElt=varElts(i)
	          Call Calc_EltPrRotations(psiElt_old(1:3,i),
     &                      psiElt(1:3,iElt),nECoord(iElt),
     &                      TElt(1:3,1:3,iElt),eltRotations,rtnOK)
	          if (nECoord(iElt) .lt. 3) then
	            WRITE(*,208) '  Element',iElt,'  (global coord frame):'
	          else
	            WRITE(*,208) '  Element',iElt,'  (element coord frame):'
	          end if
              WRITE(*,209) '   Z-rotation vector: ',eltRotations(1:3,3)
              WRITE(*,209) '   Y-rotation vector: ',eltRotations(1:3,2)
              WRITE(*,209) '   X-rotation vector: ',eltRotations(1:3,1)
	          print*,' '
 208              FORMAT(A,I3,A)
 209              FORMAT(A,3D17.9)
	        End Do
	      End If
#endif

#ifdef MACOS_CMD
	      CALL CACCEPT(ANS,'yes',
     &          'Save optimized optical state (yes,no):')
#else
	      ANS(1:3)='yes'  ! default to 'yes' for SMACOS
#endif

              IF (LCMP(ANS,'no',1)) THEN
	        ! Save Rx npts, use opt_npts for system optimization
                npts_Rx=npts; npts=opt_npts

	        ! Restore pre-optimization optical state by calling
	        ! this routine with start_optim=.false. and
	        ! nitr_dopt < 0.
	        If (.not. LOptCons) Then
                  ! Unconstrained version
	          Call nls_optim_dvr(.false.,OptTarget,
     &			 nVarElt,varElts,varDOFs,
     &                   nWFZern,WFZernMode,
     & 			 zernElt(1:nVarElt),varZerns,
     &                   asphElt(1:nVarElt),varAsphs,
     &                   opt_fov,nOptFov,fov_weight,
     &                   opt_wavelen,nOptWavelen,
     &                   dopt_tol,-1,
     &                   old_wfe,new_wfe,rtn_flag,
     &			 tgtWF_flg,tgtWF_fn,tgtWF_fmt,
     &                   OptTgtElt,LOptIfFEX,beamOpt)
	        Else
                  ! Constrained version
                  Call np_optim_dvr(.false.,OptTarget,
     &                              nVarElt,varElts,varDOFs,
     &                              varEltCons,
     &                              nWFZern,WFZernMode,
     &                              zernElt(1:nVarElt),varZerns,
     &                              opt_fov,nOptFov,fov_weight,
     &                              opt_wavelen,nOptWavelen,
     &                              dopt_tol,nitrs_dopt,
     &                              old_wfe,new_wfe,rtn_flag,
     &                              tgtWF_flg,tgtWF_fn,tgtWF_fmt,
     &                              OptTgtElt,LOptIfFEX,beamOpt)
                End If
	        npts=npts_Rx
                WRITE(*,*)'* Optimized optical state NOT saved!'
	      ELSE
                WRITE(*,*)'* Optimized optical state has been saved'
	      END IF
#ifdef MACOS_CMD
	      WRITE(*,*)' '
#endif
            ELSE ! else of ifLoad
              WRITE(*,150)
	      rtn_flag=1
            END IF

 105	    If (rtn_flag /= 0) Then
              WRITE(*,*) ' ** Optimzation aborted!'
            End If

	    If (opt_npts /= npts) Then
	      ! trace rays again when doing OPD
	      ifOPDCalc=.FALSE.; iCurRayElt=0
	    End If

#ifdef MACOS_CMD
  211       If (LCmdCalibElts) Then
              isVarElt=isVarElt_Rx; nVarElt=nVarElt_Rx
            End If
#endif


#endif
! End of system optimization commands


C         ***Plot OPD map***

	  ELSE IF (LCMP(command,'OPD',3)) THEN
	    IF (ifLoad) THEN
	      IF (ifTrace.OR.ifPropagate) THEN
	        i=iCurRayElt
	      ELSE
	        i=0
	      END IF
  6	      CALL IACCEPT_S(iElt,iCurRayElt,1,
     &	        'Enter element where OPD is to be evaluated:')
	      IF (iElt.LE.0) THEN
	        GO TO 1
	      ELSE IF (iElt.GT.nElt) THEN
	        WRITE (*,*)'OPD: Invalid element number: ', iElt
	        GO TO 6
	      ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.11).OR.
     &	      (EltID(iElt).EQ.12)) THEN
	        WRITE(*,*)'Invalid element type =',EltID(iElt)
	        GO TO 6
	      END IF
	      iEndElt=iElt

#ifdef MACOS_CMD
	      if (firstOPDCalc) then
	        OPDSave=OPDMat; OPDEltSave=OPDElt
	      end if
#endif

	      CALL OPD(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &	       ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &	       iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &	       nBadRays,StartIndRef,StartExtinc,OPDMat,WFE,WFEPV,DAvgL,
     &	       L1,Cmatrix,DrawRayVec,DrawEltVec)
	      Avgl=DAvgl; ifOPDCalc=.TRUE.
	      If (.not. firstOPDCalc) Then
                firstOPDCalc=.TRUE.
	      Else
	        multiOPDCalc=.TRUE.
	      End if
#ifdef MACOS_CMD
              OPDElt=iElt
#endif

              ! Dump OPD to a file
	      If (saveOPD_Flg) Then
	        print*,'** Dumping OPD: npts =',npts
	        Open(Unit=15,File='Opd_macos.txt')
                Do j=1,npts+1
                   Write(15,*) OPDMat(1:npts+1,j)
                End Do
	        Close(15)
	      End If

 ! Dump ray positions on all elements for 3D system rendering  09/27/2014
              If (SaveRayPosHist_Flg .AND. DumpRayPosHist_Flg) Then
                print*,'** Dumping RayPosHist array ...'
                Open(Unit=15,File='RayPosHist.txt')
	        Write(*,*)'** Writing RayPosHist: nElt,nRay =',nElt,nRay
	        i=0  ! line counter
	        Do j=0,nElt
	          Do iRay=1,nRay
	            If (LRayOKHist(iRay,j)) Then
                      i=i+1
	              !Write(15,*) j, iRay, RayPosHist(1:3,iRay,j)
	              Write(15,320) j, iRay, RayPosHist(1:3,iRay,j)
	            End If
	          End Do
	        End Do
	        Write(*,*) '** Dumping RayPosHist: Total lines =',i
	        Close(15)
	        !
	        Open(Unit=16,File='ChfRayPosHist.txt')
	        Write(*,*)'** Writing ChfRayPosHist: nElt =',nElt
	        Write(16,*)'** Writing ChfRayPosHist: nElt =',nElt
	        Do j=1,nElt
                  If (LRayOKHist(1,j)) Then
	            Write(16,*) j, vptElt(1:3,j), RayPosHist(1:3,1,j)
		  End If
                End Do
		Close(16)
	      End If

#if 0
	      Open(Unit=15,File='RBOpd.txt')
              Do iRay=2,nRay
	        If (opdRayMask(RayIndex(1,iRay),RayIndex(2,iRay))==1)
     &          Then
	          Write(15,*) iRay,RayPos(1:3,iRay)
	        End If
	      End Do
	      Close(15)
#endif


	      iEm1=iElt-1
	      !print*,'*******OPD: iElt,iEm1 =',iElt,iEm1

	      If (ifCtraceCalled) Then
	        CALL DXCALC(iElt,iEm1,dxMin(iElt),dxMax(iElt),
     &	               dxAvg(iElt),dxDev(iElt),ifRegrid,ifOPDCalc,npts)
	      End If
	      RMSWFE=WFE
	      If (BaseUnits_FLG) Then
	        WRITE(*,443) RMSWFE,cUnit,WFEPV,cUnit,DAvgl,cUnit
	      Else
	        WRITE(*,442) RMSWFE,WFEPV,DAvgl
	      End If

	      !iEm1=iElt-1

	      IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.5).AND.(iPloTyp.NE.6)
     &             .AND.(iPloTyp.NE.8).AND.(iPloTyp.NE.9)
     &	           .AND.(.NOT.ifPlot)) THEN
	        ifPlot=.TRUE.
	        CALL GRAINI
	      END IF


	      DXO=dxElt(iElt)
	      DYO=DXO
	      XMIN=0d0-DXO*FLOAT(nGridPts/2)
	      YMIN=0d0-DYO*FLOAT(nGridPts/2)-DYO
	      DO i=1,nGridPts
	        XFP(i)=XMIN
	        YFP(i)=YMIN
	        XMIN=XMIN+DXO
	        YMIN=YMIN+DYO
   	      END DO

#if 0
	      print*,'******OPD: dxElt(iElt),dxAvg(iElt) =', dxElt(iElt),
     &               dxAvg(iElt)
	      print*,'******OPD: DXO,XMIN,YMIN =', DXO,XMIN,YMIN
#endif

	      CALL IntToChar(Cinteger,iElt,i)
	      IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &            .OR.(iPloTyp.EQ.9)) THEN
	        j=ICLEN(filnam)
	        imgfil=filnam(1:j)//'.OPD'//Cinteger(1:i)
	        CTITLE='OPD Map, Elt='//Cinteger(1:i)
	        CALL DBINOUT(iPloTyp,XFP,YFP,OPDMat,nGridPts,nGridPts,
     &	        	mpts,mpts,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                  dxElt(iElt),dxElt(iElt),iElt,R1)
	      ELSE IF (iPloTyp.NE.0) THEN
	        ! Render in pgplot window, XFP and YFP are grid coord vectors
	        CTITLE='OPD Map, Elt='//Cinteger(1:i)//', File='//filnam
	        CALL DSRFPLT(iPloTyp,XFP,YFP,OPDMat,nGridPts,nGridPts,
     &	               mpts,mpts,'X-Axis','Y-Axis',CTITLE,R1,R2,1)
	      END IF
	      IF (nBadRays.GT.0) CALL WARN(LEltOK,iEndElt,nBadRays)
	    ELSE
	      WRITE(*,150)
	    END IF

#ifdef MACOS_CMD

C  Plot differential OPD
          ELSE IF (LCMP(command,'DOPD',3)) THEN
            If (ifLoad .and. multiOPDCalc) Then
              If (OPDEltSave /= OPDElt) Then
                WRITE(*,*)' **Warning: OPDs are not on same element!'
              End If
              OPDbuf=OPDMat-OPDSave
              WRITE(*,*)' RMS dOPD =',RMS2(OPDbuf,nGridPts,nGridPts)
              Call FindPV(OPDbuf,nGridPts,nGridPts,c1,c2)
              WRITE(*,*)' P-V dOPD =',c1,c2

              WRITE(*,*)' Plotting differential OPD map ...'
              IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &            .OR.(iPloTyp.EQ.9)) THEN
                j=ICLEN(filnam)
                imgfil=filnam(1:j)//'.OPD'//Cinteger(1:i)
                CTITLE='OPD Map, Elt='//Cinteger(1:i)
                CALL DBINOUT(iPloTyp,XFP,YFP,OPDbuf,nGridPts,nGridPts,
     &                  mpts,mpts,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                  dxElt(iElt),dxElt(iElt),iElt,R1)
              ELSE IF (iPloTyp.NE.0) THEN
                ! Render in pgplot window, XFP and YFP are grid coord vectors
           CTITLE='DOPD Map, Elt='//Cinteger(1:i)//', File='//filnam
           CALL DSRFPLT(iPloTyp,XFP,YFP,OPDbuf,nGridPts,nGridPts,
     &                 mpts,mpts,'X-Axis','Y-Axis',CTITLE,R1,R2,1)
              END IF
            Else
	      WRITE(*,*)' Before run DOPD command, load Rx and'
              WRITE(*,*)' run OPD command at least twice!'
            End If

C   Remove Zernike components from OPD. Intended for command-line MACOS only -jzlou
          ELSE IF (LCMP(command,'ZRM',3)) THEN
            IF (.not.ifLoad) THEN
              WRITE(*,150)
              GO TO 1
            ELSE IF (.not.ifOPDCalc) THEN
              WRITE(*,158)
              GO TO 1
            END IF
	    WRITE(*,160)

	    Call zern_interface(OPDElt,OPDMat,OPDbuf,WFE,
     &                          zernUsrOpt(OPDElt),ZRMcmd)
              ! new OPD = OPDbuf, WFE = WFE of OPDbuf
	    WRITE(*,*)'  RMS WFE (after ZRM) =',WFE
            print*,' '

#if 0
 ! Dump OPD for debugging
              Open(Unit=15,File='OpdZrm.txt')
              Do j=1,npts+1
                 Write(15,*) OPDbuf(1:npts+1,j)
              End Do
              Close(15)
#endif
            ! Plot Zernike-removed OPD map
            CALL IntToChar(Cinteger,iElt,i)
            IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &          .OR.(iPloTyp.EQ.9)) THEN
              j=ICLEN(filnam)
              imgfil=filnam(1:j)//'.OPD'//Cinteger(1:i)
              CTITLE='OPD Map, Elt='//Cinteger(1:i)
              CALL DBINOUT(iPloTyp,XFP,YFP,OPDbuf,nGridPts,nGridPts,
     &                mpts,mpts,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                dxElt(iElt),dxElt(iElt),iElt,R1)
            ELSE IF (iPloTyp.NE.0) THEN
              CTITLE='OPD Map, Elt='//Cinteger(1:i)//', File='//filnam
              CALL DSRFPLT(iPloTyp,XFP,YFP,OPDbuf,nGridPts,nGridPts,
     &               mpts,mpts,'X-Axis','Y-Axis',CTITLE,R1,R2,1)
            END IF


            CALL CACCEPT(ANS,'no','Save new OPD (no,yes):')
            IF (LCMP(ANS,'yes',1)) THEN
              OPDMat=OPDbuf
	      ifTrace=.FALSE.
              ifBuild=.FALSE.
              ifPropagate=.FALSE.
              ifPerturb=.FALSE.
              ifGrid=.FALSE.
              iCurRayElt=0
            END IF

C         *** Compute and display Zernike mode coefficients
          ELSE IF (LCMP(command,'ZCOEF',3)) THEN
            IF (.not.ifLoad) THEN
              WRITE(*,150)
              GO TO 1
            ELSE IF (.not.ifOPDCalc) THEN
              WRITE(*,158)
              GO TO 1
            END IF
	    WRITE(*,159)
	    Call zern_interface(OPDElt,OPDMat,OPDbuf,WFE,
     &                          zernUsrOpt(OPDElt),ZCOcmd)

	  ELSE IF (LCMP(command,'ZABerr',3)) THEN
	    IF (.not.ifLoad) THEN
              WRITE(*,150)
              GO TO 1
            ELSE IF (.not.ifOPDCalc) THEN
              WRITE(*,158)
              GO TO 1
            END IF
            WRITE(*,161)
            Call zern_interface(OPDElt,OPDMat,OPDbuf,WFE,
     &                          zernUsrOpt(OPDElt),ZABcmd)
	    WRITE(*,*)'  RMS WFE (after ZAB) =',WFE
            print*,' '

	    ! Plot OPD map with selected Zernike modes
            CALL IntToChar(Cinteger,iElt,i)
            IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &          .OR.(iPloTyp.EQ.9)) THEN
              j=ICLEN(filnam)
              imgfil=filnam(1:j)//'.OPD'//Cinteger(1:i)
              CTITLE='OPD Map with Selected Zernike modes, Elt='
     &               //Cinteger(1:i)
              CALL DBINOUT(iPloTyp,XFP,YFP,OPDbuf,nGridPts,nGridPts,
     &                mpts,mpts,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                dxElt(iElt),dxElt(iElt),iElt,R1)
            ELSE IF (iPloTyp.NE.0) THEN
              CTITLE='OPD Map, Elt='//Cinteger(1:i)//', File='//filnam
              CALL DSRFPLT(iPloTyp,XFP,YFP,OPDbuf,nGridPts,nGridPts,
     &               mpts,mpts,'X-Axis','Y-Axis',CTITLE,R1,R2,1)
            END IF
            CALL CACCEPT(ANS,'no','Save new OPD (no,yes):')
            IF (LCMP(ANS,'yes',1)) THEN
              OPDMat=OPDbuf
              ifTrace=.FALSE.
              ifBuild=.FALSE.
              ifPropagate=.FALSE.
              ifPerturb=.FALSE.
              ifGrid=.FALSE.
              iCurRayElt=0
            END IF
#endif

C  	  ELSE IF (LCMP(command,'DOPL',4)) THEN
C           IF (.not.ifLoad) THEN
C	      CALL CACCEPT(Rx1,' ',
C    &         'Enter first prescription file name: ')
C 	    END IF
C	    CALL CACCEPT(Rx2,' ',
C    &        'Enter second prescription file name: ')
C  	    CALL IACCEPT_S(iElt,-1,1,
C    &         'Enter last element number:')
C 	    CALL CACCEPT(command,'CHFRay',
C    &                   'Enter dOPL mode (CHFRay,BEAM): ')
C
C	    CALL dOPL(ifLoad,Rx1,Rx2,iElt,command)  ! defined in utilsub.F

C#if 1
C          *** Get optical beam size at element surface(s)
	   ELSE IF (LCMP(command,'GBS',3)) THEN
	     IF (ifLoad) THEN
               i=1
               CALL IACCEPT_S(iLocal,i,1,'Enter first element number:')
               CALL IACCEPT_S(jLocal,i,1,'Enter last element number:')
               If (iLocal < 0 .OR. jLocal < 0) Then
                 Write(*,*) 'GBS: Invalid element number input!'
                 GO TO 1
               Else If (iLocal > jLocal) Then
                 Write(*,*) 'first element must not be larger than last'
                 GO TO 1
               End If
	       LShowDtime=.FALSE.
               If (iLocal == jLocal) Then
                 If (iLocal==0) Then
                   ! source element
                   WRITE(*,167)' Aperture at source =',Aperture,cUnit
                 Else
                   iElt=iLocal
	           ifRayPrint=.FALSE. ! don't print raytrace output in 'GBS' cmd -jzlou
	           !print*,'** B4 GetBeamSizeCmd(), ifRayPrint =',ifRayPrint
               Call
     &            GetBeamSizeCmd(bmRad,npts,ifPol,ifGrid,ifPropagate,
     &            ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &            ifIntSrfSet,ifRayPrint,iObsOpt,iLocal,iCurWFElt,
     &            iCurRayElt,Dicr,Djcr,dxSource,nBadRays,StartIndRef,
     &            StartExtinc,OPDMat,WFE,DAvgL,L1,Cmatrix,DrawRayVec,
     &            DrawEltVec)
              WRITE(*,168)' At element ',iElt,' (',
     &                    EltTypeName(EltID(iElt)),')',
     &                    ', Optical beam diameter =',
     &                    2d0*bmRad,cUnit
                 End If
               Else
                 ! Save result to a file
                 Open(unit=12,file='beamSize.txt',status='replace')
                 Do iElt=iLocal,jLocal
                   If (iElt==0) Then
                     ! source element
	             If (SaveVis3d_dat_Flg) Then
		       WRITE(12,318) 'ChfRayDir =',ChfRayDir(1:3)
		       WRITE(12,318) 'ChfRayPos =',ChfRayPos(1:3)
                     End If
                     WRITE(12,*) 'iElt =  0, Source aperture =',
     &                           Aperture
                     WRITE(*,*)  'iElt =  0, Source aperture =',
     &                           Aperture
       WRITE(12,*)'-------------------------------------------------'
                   Else
                     Call
     &            GetBeamSizeCmd(bmRad,npts,ifPol,ifGrid,ifPropagate,
     &             ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &             ifIntSrfSet,ifRayPrint,iObsOpt,iElt,iCurWFElt,
     &             iCurRayElt,Dicr,Djcr,dxSource,nBadRays,StartIndRef,
     &             StartExtinc,OPDMat,WFE,DAvgL,L1,Cmatrix,DrawRayVec,
     &             DrawEltVec)

!                    Write(12,*) 'iElt =',iElt,
!    &                           ' (', EltTypeName(EltID(iElt)),')'
!                    Write(12,*) '  Beam diameter =',2d0*bmRad

                     Write(12,319)'iElt =',iElt,
     &                            ' (',EltName(iElt),', ',
     &                            EltTypeName(EltID(iElt)),')',
     &                            '  Beam diameter =',2d0*bmRad
                     Write(*,169)  'iElt =',iElt,' (',
     &                 EltTypeName(EltID(iElt)),')',
     &                 '  Optical beam diameter ='
                     WRITE(*,*) 2d0*bmRad
		     If (SaveVis3d_dat_Flg) Then
	               ! Kr,Kc,Psi,Vpt
                       Write(12,*) '  KrElt =',KrElt(iElt)
                       Write(12,*) '  KcElt =',KcElt(iElt)
                       !Write(12,*) '  PsiElt =',PsiElt(1:3,iElt)
                       Write(12,317)'  PsiElt =',PsiElt(1:3,iElt)
                       Write(12,317)'  VptElt =',VptElt(1:3,iElt)
	               If (LRayOKHist(1,iElt)) Then
                     Write(12,317) '  CRPos = ',RayPosHist(1:3,1,iElt)
                       End If
                     End If
                   End If
                 End Do
	         LShowDtime=.TRUE.
                 Close(12)
                 Write(*,*) '==> Result is saved in file beamSize.txt'
               End If
             ELSE
          WRITE(*,*)' Must load optical data first: use NEW or OLD.'
               GO TO 1
             END IF
C#endif

C         ***Set spot diagram centering option***

	  ELSE IF (LCMP(command,'SPCENTER',3)) THEN
 	    CALL CACCEPT(ANS,'elt',
     &	  'Enter spot diagram centering option (elt,chfray):')
	    IF (LCMP(ANS,'ELT',1)) THEN
	      spcOption=1
	    ELSE
	      spcOption=2
	    END IF

C         ***Plot spot diagram***

	  ELSE IF (LCMP(command,'SPOT',3)) THEN
	    IF (ifLoad) THEN
	      IF (ifTrace) THEN
	        i=iCurRayElt
	      ELSE IF (ifPropagate) THEN
	        i=iCurWFElt
	      ELSE
	        i=0
	      END IF
 10	      CALL IACCEPT_S(iElt,i,1,
     &	  'Enter element where spot diagram is to be computed:')
	      IF (iElt.LE.0) THEN
	        GO TO 1
	      ELSE IF (iElt.GT.nElt) THEN
	        WRITE (*,*)'SPOT: Invalid element number: ',iElt
	        GO TO 10
	      ELSE IF
     &        ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.11).OR.
     &	      (EltID(iElt).EQ.12)) THEN
	        WRITE(*,*)'Invalid element type =',EltID(iElt)
	        GO TO 10
	      END IF

C	  Get local coordinates:
	      IF (nOutCord.NE.5) THEN
	        CALL CACCEPT(ANS,'Beam',
     &	  'Enter output coordinate option (Enter or Beam):')
	      ELSE
	        CALL CACCEPT(ANS,'Tout',
     &	  'Enter output coordinate option (Tout, Enter or Beam):')
	      END IF
	      IF (LCMP(ANS,'YES',1).OR.LCMP(ANS,'TOUT',3)) THEN
	        CALL COMPACT(Tcompact,Tout,5,7,7,7)
	        CALL DCOMPACT(DTcompact,Tout,5,7,7,7)
	        iEndElt=iElt
	        DO i=1,3
	          xLocal(i)=Tout(1,i)
	          yLocal(i)=Tout(2,i)
	        END DO
	      ELSE IF (LCMP(ANS,'ENTER',3)) THEN
	        invec7(1)=-psiElt(3,iElt)
	        invec7(2)=-psiElt(1,iElt)
	        invec7(3)=-psiElt(2,iElt)
	    CALL DACCEPT(invec7,invec7,3,'Enter xOut (x,y,z) ')
	    CALL DUNITIZE(invec7)
	    CALL DXPROD(DDEF,psiElt,invec7)
	    CALL DNEGATE(DDEF,DDEF,3)
	    CALL DACCEPT(DDEF,DDEF,3,'Enter yOut (x,y,z) ')
	    CALL DUNITIZE(DDEF)
	    CALL DZERO(Tspot(1,1),49)
	        DO iCoord=1,3
	          Tspot(1,iCoord)=invec7(iCoord)
	          Tspot(3,3+iCoord)=invec7(iCoord)
	          Tspot(2,iCoord)=DDEF(iCoord)
	          Tspot(4,3+iCoord)=DDEF(iCoord)
	        END DO
	        Tspot(5,7)=1
	        CALL COMPACT(Tcompact,Tspot,5,7,7,7)
	        CALL DCOMPACT(DTcompact,Tspot,5,7,7,7)
	        iEndElt=iElt
	        DO i=1,3
	          xLocal(i)=Tspot(1,i)
	          yLocal(i)=Tspot(2,i)
	        END DO
	      ELSE
	        ifGrid=.FALSE.
	        ifDraw=.FALSE.
	        ifPropagate=.FALSE.
	        ifTrace=.FALSE.
	        ifBuild=.FALSE.
	        ifOPDCalc=.FALSE.
	        iLocal=iElt
	        CALL
     &      LocalCoord(npts,ifPol,ifGrid,ifPropagate,ifDraw,
     &	    ifLNegOK,ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,
     &	    ifRayPrint,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &	    Dicr,Djcr,dxSource,nBadRays,StartIndRef,StartExtinc,
     &	    Cmatrix,DrawRayVec,DrawEltVec,
     &	    iLocal,xLocal,yLocal,zLocal,crLocalPos,ifLocalCoord)
	        IF (ifLocalCoord) THEN
	          WRITE(*,311)xLocal,yLocal,zLocal
	          CALL DZERO(Tspot(1,1),49)
	          DO iCoord=1,3
	            Tspot(1,iCoord)=xLocal(iCoord)
	            Tspot(3,3+iCoord)=xLocal(iCoord)
	            Tspot(2,iCoord)=yLocal(iCoord)
	            Tspot(4,3+iCoord)=yLocal(iCoord)
	          END DO
	          Tspot(5,7)=1
	          CALL COMPACT(Tcompact,Tspot,5,7,7,7)
	          CALL DCOMPACT(DTcompact,Tspot,5,7,7,7)
	        ELSE
	          GO TO 10
	        END IF
	      END IF
	      RefSpot(1)=DDOTC(xLocal,VptElt(1,iEndElt))
	      RefSpot(2)=DDOTC(yLocal,VptElt(1,iEndElt))

C	  Trace rays:
	      iEndElt=iElt
	      IF (ifTrace.AND.(iCurRayElt.LE.iEndElt)) THEN
	        iStartElt=iCurRayElt
	      ELSE IF
     &          (ifPropagate.AND.(iCurWFElt.EQ.iEndElt)) THEN
	        iStartElt=iCurWFElt
	      ELSE
	        iStartElt=0
	      END IF

	      IF ((iStartElt.EQ.0)
     &            .OR.(iStartElt.LT.iEndElt)) THEN
	        ifRayPrint=.FALSE.
	        CALL RTIME(StartTime)
	        IF (iStartElt.EQ.0) THEN
	          IF (ABS(zSource).GE.1d10) THEN
	            CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	            ifPropagate)
	          ELSE
	            CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	            ifPropagate)
	          END IF
	          ifGrid=.TRUE.
	        END IF
	        WRITE(*,*)' STOP cmd: Tracing ',nRay,' rays...'
	        iStartRay=1
	        iEndRay=nRay
	        CALL
     &     CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &	   ifLNegOK,Cmatrix,iStartElt,iEndElt,iStartRay,iEndRay,
     &     RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &           StartIndRef,StartExtinc,DrawRayVec,DrawEltVec,
     & 		 nDrawElt(1),xDraw,yDraw)
	        iCurRayElt=iEndElt
	        ifTrace=.TRUE.
	        ifBuild=.FALSE.
	        ifPropagate=.FALSE.
	        CALL RTIME(EndTime)
	        DTime=EndTime-StartTime
	        WRITE(*,816)DTime
	      END IF

	      !WRITE(*,*)' STOP cmd: iPloTyp =', iPloTyp

C	      Initialize graphics:
	      ! change by jzlou, 01/07/05
	      IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.6).AND.
     &              (.NOT.ifPlot)) THEN
	        ifPlot=.TRUE.
	        CALL GRAINI
	      END IF

C	      Compute spot diagram data from ray trace and plot:
	      iSpot=0
	      nGridPts=npts+1
	      mpt2=mpts*mpts
	      xcent=0d0
	      ycent=0d0
	      DO 823 iRay=1,nRay
	   IF (LRayOK(iRay).AND.((iRay.EQ.1).OR.(iObsOpt.EQ.0).OR.
     &	         (LRayPass(iRay).AND.(iObsOpt.EQ.1)).OR.
     &	         (.NOT.LRayPass(iRay).AND.(iObsOpt.EQ.2))))THEN
	          DO 22 iaxis=1,3
	            j=iaxis+3
	            xdum(iaxis)=RayDir(iaxis,iRay)
	            xdum(j)=RayPos(iaxis,iRay)
 22	          CONTINUE
	          xdum(7)=0e0
	          CALL DMPRDB(xnom,DTcompact,xdum,5,7,1)
	          IF (iRay.EQ.1) THEN
	            CntrSpot(1)=xnom(3)
	            CntrSpot(2)=xnom(4)
	          END IF
	          IF ((iRay.NE.1).OR.
     &	            ((iRay.EQ.1).AND.(nGridpts.EQ.1))) THEN
	            iSpot=iSpot+1
	            IF (spcOption.EQ.1) THEN
	              RaySpot(iSpot,1)=xnom(3)-RefSpot(1)
	              RaySpot(iSpot,2)=xnom(4)-RefSpot(2)
	            ELSE
	              RaySpot(iSpot,1)=xnom(3)-CntrSpot(1)
	              RaySpot(iSpot,2)=xnom(4)-CntrSpot(2)
	            END IF
	            xcent=xcent+xnom(3)
	            ycent=ycent+xnom(4)
	          END IF
	        END IF
 823	      CONTINUE

	      NBEAMS=1
	      IF (iSpot.GT.0) THEN
	        xcent=xcent/DBLE(iSpot)
	        ycent=ycent/DBLE(iSpot)
	        WRITE(*,408)(CntrSpot(1)-RefSpot(1)),
     &	          (CntrSpot(2)-RefSpot(2)),xcent,ycent
	      ELSE
	        WRITE(*,*)' All rays were obscured or lost'
	      END IF

#if 0
 ! dump ray spot
               Open(Unit=15,File='sspot.txt')
               Do j=1,iSpot
                 Write(15,*) RaySpot(j,1:2)
               End Do
               Close(15)
#endif

	      CALL IntToChar(Cinteger,iElt,i)
	  IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &            .OR.(iPloTyp.EQ.9)) THEN
	        j=ICLEN(filnam)
	        imgfil=filnam(1:j)//'.spot'//Cinteger(1:i)
	        CTITLE='Spot Diagram, Elt='//Cinteger(1:i)

	        ! write to a file
	  CALL SPOTOUT(RaySpot,iSpot,2,mRay,2,1,
     &	         'X-Axis','Y-Axis',CTITLE,iPloTyp,imgfil,iElt)

	      ELSE IF (iPloTyp.NE.0) THEN
	        CTITLE='Spot Diagram, Elt='//Cinteger(1:i)//
     &	        ', File='//filnam

	        ! go to pgplot
	        CALL SpotDiag(RaySpot,iSpot,2,mRay,2,1,
     &	                      'X-Axis','Y-Axis',CTITLE)
	      END IF
	 IF (nBadRays.GT.0) CALL WARN(LEltOK,iEndElt,nBadRays)

	    ELSE
	 WRITE(*,*)' Must load optical data first: use NEW or OLD.'
	      GO TO 1
	    END IF

C         ***Export data***

	  ELSE IF (LCMP(command,'EXPORT',3)) THEN
	    IF (ifBuild.OR.ifGrid.OR.ifTrace.OR.ifPropagate) THEN
	 WRITE(*,*)' Current element is ',iCurRayElt
	 CALL CEXPORT(ifBuild,npts,GridType,RayMap,SegMap,RayIndex,
     &        Cmatrix,iPertElt,iCurRayElt,ifDMBuild,iObsOpt)
	    ELSE
	      WRITE(*,*)' No data has been generated yet.'
	    END IF

C         ***next command: Plot wavefront stuff***
C                          -- original comment, doesn't make sense -jzlou
C         *** Should perform ray trace and diffraction here -- jzlou

	  ELSE IF (LCMP(command,'IN',2).OR.LCMP(command,'PH',2).OR.
     &      LCMP(command,'AMPL',3).OR.LCMP(command,'REAL',3).OR.
     &      LCMP(command,'LOG',3).OR.LCMP(command,'PI',2).OR.
     &      LCMP(command,'MTF',3).OR. LCMP(command,'IMG',3).OR.
     &      LCMP(command,'DADD',3).OR.LCMP(command,'BLUR',3).OR.
     &      LCMP(command,'ADD',2).OR.LCMP(command,'GAIN',3).OR.
     &      LCMP(command,'CADD',3).OR.LCMP(command,'GBLUR',3)) THEN

	    ! 'INtensity' and 'PIxel' commands come here -jzlou

	    IF (nGridPts.EQ.1) THEN
	      WRITE(*,*)' Too few rays, reset nGridPts.'
	      GO TO 1
	    END IF
	    IF (.NOT.ifLoad) THEN
	      WRITE(*,150)
	      GO TO 1
	    END IF

	    IF (LCMP(command,'IMG',3)) THEN
	      ! Abort if in vector mode or input object intensity
	      ! file not available ...
	      If (.not.ifObjImgSet) Then
	        WRITE(*,*)' * Object intensity input not available'
	      Else If (ifVecDif3) Then
	        WRITE(*,*)
     &         ' * Imaging not available in vector diffraction mode'
	      End If
	      If (ifVecDif3 .OR. (.not.ifObjImgSet)) GO TO 1
	    End If

!  'ADD', 'DADD' and 'CADD': should be run after executing 'COMPOSE' cmd;
!  each of these cmds just sets 'iElt' to 'iComElt' that is set by COMPOSE,
!  so the subsequent RunProp() will propagate to the proper iElt. -jzlou
	    IF (LCMP(command,'ADD',2)) THEN
	      IF (.NOT.ifAdd) THEN
	        WRITE(*,*)' Must COMpose before ADDing images'
	        GO TO 1
	      END IF
	      iElt=iComElt
	    ELSE IF (LCMP(command,'DADD',3)) THEN
	      IF (.NOT.ifPixArray) THEN
	        WRITE(*,*)' No data to display'
	        GO TO 1
	      END IF
	      iElt=iComElt
#if 0
	    ELSE IF (LCMP(command,'CADD',3)) THEN
	      IF (.NOT.ifRegrid(iElt)) THEN
	        WRITE(*,*)' Must REGRID before adding complex WFs'
	        GO TO 1
	      END IF
	      iElt=iComElt
#endif
	    ELSE IF (LCMP(command,'BLUR',3)) THEN
	      IF (.NOT.ifPixArray) THEN
	        WRITE(*,*)' No data to display'
	        GO TO 1
	      END IF
	      iElt=iComElt
	    ELSE IF (LCMP(command,'GBLUR',3)) THEN
	      IF (.NOT.ifAdd) THEN
	        WRITE(*,*)' Must COMpose before GBLURing images'
	        GO TO 1
	      END IF
	      iElt=iComElt
	    ELSE
	      IF ((iCurWFElt.LE.0).OR.(iCurWFElt.GT.nElt)) THEN
	        i=nElt
	      ELSE
	        i=iCurWFElt
	      END IF
19	      CALL IACCEPT_S(iElt,i,1,
     &	        'Enter number of element where data is to be generated:')
	      IF ((iElt.LE.0).OR.(iElt.GT.nElt)) THEN
	        WRITE (*,*)' WF-stuff: Invalid element number: ',iElt
	        GO TO 1
	      ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.11).OR.
     &	      (EltID(iElt).EQ.12)) THEN
	        WRITE(*,*)'Invalid element type =',EltID(iElt)
	        GO TO 1
	      END IF
	    END IF

	    iEm1=iElt-1

	    ! iEndElt=0 here, normally; iEndElt is actually un-initialized
            ! at this moment, and it is initialized in RunProp();
            ! iElt at this moment is supposed to be the 'End' element the
            ! propagation finishes at. -jzlou

	    CALL RunProp(iElt,npts,ifPol,ifGrid,ifPropagate,
     &	     ifBuild,ifTrace,ifRegrid,ifVecDif3,ifIntSrfSet,ifLNegOK,
     &	     iEndElt,iCurWFElt,iCurRayElt,iWF,Dicr,Djcr,dxSource,
     &	     iObsOpt,nBadRays,StartIndRef,StartExtinc,DWF,OPDMat)

	    IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.5).AND.(iPloTyp.NE.6)
     &      	.AND.(iPloTyp.NE.8).AND.(iPloTyp.NE.9)
     &	    	.AND.(.NOT.ifPlot)) THEN
	      ifPlot=.TRUE.
	      CALL GRAINI
	    END IF

	    du=dxElt(iEm1)
	    udia=du*FLOAT(mdttl-1)
	    dx=dxElt(iElt)
	    xdia=FLOAT(mdttl-1)*dx

	    ! du is grid spacing of complex optical field at starting element;
            ! dx is grid spacing of complex optical field at 'End'
            ! element. -jzlou

 206    FORMAT(1P,' Wavefront Propagation Data Summary:'/
     &    ' Wavelength=',d14.7,';  Source Flux=',d14.7/
     &    ' u-v plane diam=',D14.7,'  du=',D14.7/
     &    ' x-y plane diam=',D14.7,'  dx=',D14.7)
 418	FORMAT(1P,' Peak intensity=',D14.7,
     &    '; Peak occurs at i=',i5,', j=',i5/
     &	  ' Sum of intensity=',D14.7)
 419	FORMAT(1P,' (2)Peak intensity=',D14.7,
     &    '; Peak occurs at i=',i5,', j=',i5/
     &	  ' Sum of intensity=',D14.7/' Maximum gain=',D14.7)
 420	FORMAT(1P,' Peak amplitude=',D14.7,'; Peak occurs at i=',i5,
     &    ', j=',i5)

	    WRITE(*,206) Wavelen,Flux,udia,du,xdia,dx

	    ! Calc physical (image plane) coords for intensity
	    X0PFP=0d0-dx*FLOAT((mdttl)/2)
	    DO i=1,mdttl
	      XFP(i)=X0PFP; YFP(i)=X0PFP
	      X0PFP=X0PFP+dx
   	    END DO ! XFP and YFP are used in plotting routines below

	    IF (LCMP(command,'INTENSITY',2)) THEN

	      CALL Ca2Int(ifVecDif3,Power,MaxInt,imax,jmax)
	        ! In Ca2Int, the computed intensity is stored
	        ! in MWFFT array
	      WRITE(*,418) MaxInt,imax,jmax,Power
	      IF (StretchType.EQ.2) THEN
	        ANS='LOG10 Intensity#%'
	      ELSE IF (StretchType.EQ.3) THEN
	        ANS='SQRT Intensity#%'
	      ELSE
	        ANS='Intensity#%'
	      END IF

	      CALL SrfOut(iElt,iPloTyp,ANS,
     &	        'X-Axis','Y-Axis','.int',mdttl,mdttl,R1,R2,StretchType)
	      !print*,'****Int: after SrfOut, mdttl =',mdttl
	      firstIntCalc=.true.

#if 0
   ! Dump intensity for debugging
               Open(Unit=15,File='Int.txt')
               Do j=1,mdttl
                 Write(15,*) MWFFT(1:mdttl,j)
               End Do
               Close(15)
	       print*,'******* Int.txt created!'
#endif

	    ELSE IF (LCMP(command,'MTF',3)) THEN
	      CALL Ca2Int(ifVecDif3,Power,MaxInt,imax,jmax)
	      CALL Int2MTF(MWFFT,mdttl,dxElt(iElt),XFP,YFP)

#ifdef MACOS_CMD
	      if (.not.ifEFLCalc) then
#include "efl.inc"
	        if (abs(Aperture) > smallestReal)
     &            focalRatio=efl/Aperture
	      end if

	      saviPloTyp=iPloTyp
	      if (iPloTyp==2) then
	         iPloTyp=1 ! set to slice (no gray image for MTF)
	      end if
	      ANS='MTF Plot#%'
	      CALL SrfOut(iElt,iPloTyp,ANS,
     &     'Cycle/'//BaseUnits(1:ICLEN(BaseUnits)),'Normalized MTF',
     &     '.mtf',mdttl,mdttl,R1,R2,StretchType)
	      if (iPloTyp==1) then
	        print*,' '
        print*,'Use ROW or COL command to plot 1D MTF in each dimension'
	        print*,' '
	      end if
	      !print*,'Image plane sampling dx =', 1d0/(mdttl*dxElt(iElt))
	      !print*,'Nyquist sampling dx =', (WaveLen*focalRatio/2d0)
	      if ((WaveLen*focalRatio/2d0)<(1d0/(mdttl*dxElt(iElt)))) then
	        !print*,'**Warning: There could be aliasing in MTF,'
	        !print*,'**         use larger ray-tracing grid and macos model'
	      end if
	      iPloTyp=saviPloTyp
#endif

	    ELSE IF (LCMP(command,'AMPLITUDE',2)) THEN
	      IF (ifVecDif3) THEN
	        CALL Ca2Amp(1,power,MaxInt,imax,jmax)
	        WRITE(*,420)MaxInt,imax,jmax
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Amplitude X#%',
     &	        'X-Axis','Y-Axis','.ampX',mdttl,mdttl,R1,R2,1)
	        CALL Ca2Amp(2,power,MaxInt,imax,jmax)
	        WRITE(*,420)MaxInt,imax,jmax
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Amplitude Y#%',
     &	        'X-Axis','Y-Axis','.ampY',mdttl,mdttl,R1,R2,1)
	        CALL Ca2Amp(3,power,MaxInt,imax,jmax)
	        WRITE(*,420)MaxInt,imax,jmax
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Amplitude Z#%',
     &	        'X-Axis','Y-Axis','.ampZ',mdttl,mdttl,R1,R2,1)
	      ELSE
	        CALL Ca2Amp(iWF,power,MaxInt,imax,jmax)
	        WRITE(*,420)MaxInt,imax,jmax
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Amplitude#%',
     &	        'X-Axis','Y-Axis','.amp',mdttl,mdttl,R1,R2,1)
	      END IF

#if 0
 ! Dump Amplitude for debugging
              Open(Unit=15,File='ampl.txt')
              Do j=1,mdttl
                 Write(15,*) MWFFT(1:mdttl,j)
              End Do
              Close(15)
#endif

	    ELSE IF (LCMP(command,'PHASE',3)) THEN
	      IF (ifVecDif3) THEN
	        CALL Ca2Pha(1)
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Phase X#%',
     &	          'X-Axis','Y-Axis','.phaseX',mdttl,mdttl,
     &	          R1,R2,1)
	        CALL Ca2Pha(2)
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Phase Y#%',
     &	          'X-Axis','Y-Axis','.phaseY',mdttl,mdttl,
     &	          R1,R2,1)
	        CALL Ca2Pha(3)
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Phase Z#%',
     &	          'X-Axis','Y-Axis','.phaseZ',mdttl,mdttl,
     &	        R1,R2,1)
	      ELSE
	        CALL Ca2Pha(iWF)
	        CALL SrfOut(iElt,iPloTyp,'Wavefront Phase#%',
     &	          'X-Axis','Y-Axis','.phase',mdttl,mdttl,
     &	          R1,R2,1)
#if 0
 ! Dump PHASE for debugging
              Open(Unit=15,File='phase.txt')
              Do j=1,mdttl
                 Write(15,*) MWFFT(1:mdttl,j)
              End Do
              Close(15)
#endif
	      END IF

	    ELSE IF (LCMP(command,'IMG',3)) THEN
	      If (iElt==nElt) Then
	        WRITE(*,*)
     &		' Simulating incoherent or coherent imaging ...'
	      Else
	        WRITE(*,*) ' Imaging must be on focal plane!'
	        GO TO 1
	      End If

	      ImgMode=InCohImg  ! incoherent is default
	      ! Simulate coherent or incoherent imaging with computed INT
	      ! or CA; focal plane psf has been computed.
	      CALL CACCEPT(ANS,'INC','INCoherent or COHerent imaging? ')
	      If (.not.LCMP(ANS,'INC',3)) Then
	        ! coherent
	        ImgMode=CohImg
	      End If

	      ! defined in sunsub.F
 	      Call SysImging(ImgMode,WFElt(:,:,1),ObjInt,
     &                       MWFFT,mdttl)

#if 0
 ! for test only
	      Open(unit=13,file='caReal.txt',status='replace')
	      Open(unit=14,file='caImag.txt',status='replace')
	      Open(unit=15,file='mwfft.txt',status='replace')
	      !MWFFT=WFElt(:,:,1)*DCONJG(WFElt(:,:,1))
	      Do j=1,mdttl
	        Write(13,*) DREAL(WFElt(1:mdttl,j,1))
	        !Write(13,*) DBLE(WFElt(1:mdttl,j,1))
	        Write(14,*) DIMAG(WFElt(1:mdttl,j,1))
	        !Write(15,*) MWFFT(1:mdttl,j)
	      End Do
	      Close(13)
	      Close(14)
	      Close(15)
#endif

	      If (ImgMode==InCohImg) Then
	       CALL SrfOut(iElt,iPloTyp,'Incoherent Imaging#%',
     &         'X-Axis','Y-Axis','.intensity',mdttl,mdttl,R1,R2,1)
	      Else
	       CALL SrfOut(iElt,iPloTyp,'Coherent Imaging#%',
     &         'X-Axis','Y-Axis','.intensity',mdttl,mdttl,R1,R2,1)
	      End If

	    ELSE IF (LCMP(command,'GAIN',3)) THEN
	      CALL Ca2Gain(ifVecDif3,iElt,Power,MaxInt,imax,jmax)
	      WRITE(*,419)MaxInt,imax,jmax,Power,MWFFT(imax,jmax)
	      CALL SrfOut(iElt,iPloTyp,'Gain#%',
     &	      'X-Axis','Y-Axis','.gain',mdttl,mdttl,R1,R2,1)

	    ELSE IF (LCMP(command,'LOG',3)) THEN
	      CALL Ca2Log(ifVecDif3,Power,MaxInt,imax,jmax)
	      WRITE(*,418)MaxInt,imax,jmax,Power

	      print*,'*****B4 SrfOut: mdttl =',mdttl

	      CALL SrfOut(iElt,iPloTyp,'Log10 Wavefront Intensity#%',
     &	                  'X-Axis','Y-Axis','.log',mdttl,mdttl,R1,R2,1)

#if 1
   ! Dump intensity for debugging
               Open(Unit=15,File='IntLog.txt')
               Do j=1,mdttl
                 Write(15,*) MWFFT(1:mdttl,j)
               End Do
               Close(15)
               print*,'******* IntLog.txt created!'
#endif


	    ELSE IF (LCMP(command,'REAL',3)) THEN
	      DO j=1,mdttl
	        DO i=1,mdttl
	  	    MWFFT(i,j)=SNGL(DREAL(WFElt(i,j,iWF))) ! sab update
	  	    !MWFFT(i,j)=REAL(WFElt(i,j,iWF)) ! original
   	        END DO
   	      END DO
	      CALL SrfOut(iElt,iPloTyp,'Real Part of the Wavefront#%',
     &	      	'X-Axis','Y-Axis','.real',mdttl,mdttl,R1,R2,1)
	      DO j=1,mdttl
	        DO i=1,mdttl
	  	    MWFFT(i,j)=DIMAG(WFElt(i,j,iWF))
   	        END DO
   	      END DO
	      CALL SrfOut(iElt,iPloTyp,
     &	       	'Imaginary Part of the Wavefront#%','X-Axis','Y-Axis',
     &	        '.imag',mdttl,mdttl,R1,R2,1)

	    ELSE IF (LCMP(command,'PIXEL',2)) THEN
	      ! This 'ELSE IF' is under the main 'ELSE IF' command structure
	      ! for all diffraction-image commands.
	      ! By now, RunProp() has been run; ray trace and diffraction all
	      ! have been performed. -jzlou
	      i=mPix
	      iComElt=iElt  ! typically iElt=nElt at this moment
	      ifPixArray=.TRUE.
	      IF (ifAdd) THEN
	        WRITE(*,*)' Overwriting previous COMPOSED image...'
	        ifAdd=.FALSE.
	      END IF
	      CALL IACCEPT_S(npix,i,1,'Enter number of pixels per side:')
	      IF (npix.GT.mPix) THEN
	        npix=mPix
	        WRITE(*,*)' WARNING: NUMBER OF PIXELS RESET TO', npix
	      END IF
	      nPixAdd=nPix
	      dxWF=DABS(dxElt(iElt))  ! dxWF is default pixel size, -jzlou
	      CALL RACCEPT_S(dxpix,dxWF,1,'Enter size of pixel:')
	      S1=dxWF/2d0
	      IF (dxpix.LT.S1) WRITE(*,*)' WARNING: PIXELS TOO SMALL'

	      CALL Ca2Int(ifVecDif3,Power,MaxInt,imax,jmax)

	      WRITE(*,418) MaxInt,imax,jmax,Power

	      ! ifPixLoc will be 'true' after executing 'WINDOW' command,
	      ! -jzlou
	      IF (ifPixLoc) THEN
	        CALL DSUB(DV1,VptElt(1,iElt),RayPos(1,1),3) ! DV1=VptElt-ChfRayPos
	        drx=DTcompact(3,4)*DV1(1)+DTcompact(3,5)*DV1(2)+
     &	          DTcompact(3,6)*DV1(3)+(WinCenPix(1)-EltPix(1))*sizPix
	        dry=DTcompact(4,4)*DV1(1)+DTcompact(4,5)*DV1(2)+
     &	          DTcompact(4,6)*DV1(3)+(WinCenPix(2)-EltPix(2))*sizPix
	      ELSE
	        drx=0.0; dry=0.0
	      END IF

	      !print*,'PIX: DV1(1:3) =', DV1(1:3)
	      !print*,'PIX: iCurWFElt =', iCurWFElt
	      !print*,'PIX: PropType(iCurWFElt-1) =', PropType(iCurWFElt-1)

	      IF (PropType((iCurWFElt-1)).NE.1) THEN
	        ! =1 for Geometric, other values for diffraction types. -jzlou
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	               npix,mPix,mdttl,dxWF,dxpix,iPloTyp,ifPixLoc,ifAdd)
	      ELSE
	        ! Geometric propagation
C	        CALL CPIX_J(PixArray,MWFFT,XFP,YFP,drx,dry,
C     &	          npix,mpix,dxpix,iPloTyp,RayPos,RayID,
C     &	          psiElt(1,iElt),VptElt(1,iElt),LRayOK,ifPixLoc,ifAdd)
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	               npix,mPix,mdttl,dxWF,dxpix,iPloTyp,ifPixLoc,ifAdd)
	      END IF
	      IF (StretchType.EQ.2) THEN
	        ANS='LOG10 Pixilated Intensity#%'
	      ELSE IF (StretchType.EQ.3) THEN
	        ANS='SQRT Pixilated Intensity#%'
	      ELSE
	        ANS='Pixilated Intensity#%'
	      END IF

	      if (.false.) then ! debugging
	        pixMax=maxval(PixArray); pixMin=minval(PixArray)
	        print*,'PIX: pixMax, pixMin=', pixMax,pixMin
	      end if

	      CALL PixOut(PixArray,iElt,iPloTyp,ANS,
     &	        'X-Axis','Y-Axis','.pix',npix,mpix,R1,R2,StretchType)

#if 0
               ! Dump PixArray for debugging
               Open(Unit=15,File='PixArray.txt')
               Do j=1,npix
                 Write(15,*) PixArray(1:npix,j)
               End Do
               Close(15)
#endif

	    ELSE IF (LCMP(command,'ADD',2)) THEN
	      i=mdttl/2
	      dxWF=DABS(dxElt(iElt))
	      S1=dxWF/2d0
	      IF (dxPixAdd.LT.S1) WRITE(*,*)' WARNING: PIXELS TOO SMALL'
	      !print*,'** dxPixAdd =',dxpixAdd
	      !print*,'** iElt =',iElt,', dxWF =',dxWF
	      CALL Ca2Int(ifVecDif3,Power,MaxInt,imax,jmax)
	      WRITE(*,418)MaxInt,imax,jmax,Power
	      IF (ifPixLoc) THEN
	        CALL DSUB(DV1,VptElt(1,iElt),RayPos(1,1),3)
	        drx=Tcompact(3,4)*DV1(1)+Tout(3,5)*DV1(2)+
     &	        Tout(3,6)*DV1(3)+(WinCenPix(1)-EltPix(1))*sizPix
	        dry=Tcompact(4,4)*DV1(1)+Tout(4,5)*DV1(2)+
     &	        Tout(4,6)*DV1(3)+(WinCenPix(2)-EltPix(2))*sizPix
	      ELSE
	        drx=0.0
	        dry=0.0
	      END IF
	      IF (PropType((iCurWFElt-1)).NE.1) THEN
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	        nPixAdd,mPix,mdttl,dxWF,dxPixAdd,iPloTyp,ifPixLoc,ifAdd)
	      ELSE
c	        CALL CPIX_J(PixArray,MWFFT,XFP,YFP,drx,dry,
c     &	        nPixAdd,mdttl,dxpix,iPloTyp,RayPos,RayID,
c     &	        psiElt(1,iElt),VptElt(1,iElt),LRayOK,ifPixLoc,ifAdd)
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	        nPixAdd,mPix,mdttl,dxWF,dxPixAdd,iPloTyp,ifPixLoc,ifAdd)
	      END IF
	      CALL IntToChar(Cinteger,iElt,i)

	      ifPixArray=.TRUE.
	      CALL CACCEPT(ANS,'YES',
     &          'Image added. Plot composed image? ')
	      IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	        IF (StretchType.EQ.2) THEN
	          ANS='LOG10 Composed Image#%'
	        ELSE IF (StretchType.EQ.3) THEN
	          ANS='SQRT Composed Image#%'
	        ELSE
	          ANS='Composed Image#%'
	        END IF
	        CALL PixOut(PixArray,iElt,iPloTyp,ANS,
     &	        'X-Axis','Y-Axis','.add',nPixAdd,mpix,R1,R2,StretchType)
	      END IF

	    ELSE IF (LCMP(command,'DADD',2)) THEN
	      ! display composed (added) image
	      IF (StretchType.EQ.2) THEN
	        ANS='LOG10 Composed Image#%'
	      ELSE IF (StretchType.EQ.3) THEN
	        ANS='SQRT Composed Image#%'
	      ELSE
	        ANS='Composed Image#%'
	      END IF
	      CALL PixOut(PixArray,iElt,iPloTyp,ANS,
     &	      'X-Axis','Y-Axis','.add',nPixAdd,mpix,R1,R2,StretchType)

	    ELSE IF (LCMP(command,'CADD',3)) THEN
	      If (.not.ifAddCA) Then
	        ifAddCA=.true.
	        iComElt=iElt  ! target element for combining complex WF

	        ! Resample CA to a rectangular grid at the image plane.
	        ! The Grid has dims mdttl by mdttl, centered around the
	        ! chief ray of first CA.
	        !Call CAReGrid(...)
	      End If

	      ! Resample and add complex WFs, need implementatioin ...
              ! CaArray stores resampled and accumulated CA fields.
	      If (ifAddCA) Then
	        !Call ResampAddCa(CaArray,...)
	      End If

	      CALL CACCEPT(ANS,'YES',
     &          'Complex amplitude added. Plot composed image? ')
              IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	        ! need implementation
	        ! print out the intensity of CaArray

	        ifAddCA=.false.  ! end compose complex amplitude
	      END If

	    ELSE IF (LCMP(command,'BLUR',2)) THEN
	      IF (StretchType.EQ.2) THEN
	        ANS='LOG10 Composed Image#%'
	      ELSE IF (StretchType.EQ.3) THEN
	        ANS='SQRT Composed Image#%'
	      ELSE
	        ANS='Composed Image#%'
	      END IF
	      CALL RACCEPT_S(kval,1d0,1,'Enter blur kernel peak (1e0):')
	      CALL CBLUR(PixArray,nPixAdd,mPix,kval,R1)
	      CALL PixOut(PixArray,iElt,iPloTyp,ANS,
     &	      'X-Axis','Y-Axis','.add',nPixAdd,mpix,R1,R2,StretchType)

	    ELSE IF (LCMP(command,'GBLUR',2)) THEN
c  set plot label
	      IF (StretchType.EQ.2) THEN
	        ANS='LOG10 Composed Image#%'
	      ELSE IF (StretchType.EQ.3) THEN
	        ANS='SQRT Composed Image#%'
	      ELSE
	        ANS='Composed Image#%'
	      END IF
c  convert WF amplitude to intensity
	      i=mdttl/2
	      dxWF=DABS(dxElt(iElt))
	      S1=dxWF/2d0
	      IF (dxPixAdd.LT.S1) WRITE(*,*)' WARNING: PIXELS TOO SMALL'
	      CALL Ca2Int(ifVecDif3,Power,MaxInt,imax,jmax)
	      WRITE(*,418)MaxInt,imax,jmax,Power
	      IF (ifPixLoc) THEN
	        CALL DSUB(DV1,VptElt(1,iElt),RayPos(1,1),3)
	        drx=Tcompact(3,4)*DV1(1)+Tout(3,5)*DV1(2)+
     &	        Tout(3,6)*DV1(3)+(WinCenPix(1)-EltPix(1))*sizPix
	        dry=Tcompact(4,4)*DV1(1)+Tout(4,5)*DV1(2)+
     &	        Tout(4,6)*DV1(3)+(WinCenPix(2)-EltPix(2))*sizPix
	      ELSE
	        drx=0.0; dry=0.0
	      END IF

c  kl--blur the intensity
C add blur function  3/11/99
c  get blur parameters
	      CALL DZERO(DDEF,9)

 	      CALL DACCEPT(DDEF,DDEF,2,
     &	        'Enter Gaussian blur kernel width in base units (0,0):')
	      sigx=DDEF(1); sigy=DDEF(2)

C  This part calculates the FT of the Gaussian
	      IF ((sigx.NE.0d0).AND.(sigy.NE.0d0)) THEN
	        NN(1)=mdttl
	        NN(2)=mdttl
	        sigx=sigx/dxWF
	        sigy=sigy/dxWF
	        c1=sigx*sigy*TWOPI
	        DO j=1,mdttl
	          j0 = j - mdttl/2 - 1
	          c3=DEXP(-DBLE(j0*j0)/(2.d0*sigy*sigy))
	          DO i=1,mdttl
	            i0 = i - mdttl/2 - 1
	            c2=DEXP(-DBLE(i0*i0)/(2.d0*sigx*sigx))
	            DWF(i,j)=c2*c3/c1
	          END DO
	        END DO
	        CALL DFOURN(DWF,NN,2,1)

C  This part takes the FT of the field
	        DO j=1,mdttl
	          DO i=1,mdttl
	            WFbuff(i,j)=DBLE(MWFFT(i,j))
	          END DO
	        END DO
	        CALL DFOURN(WFbuff,NN,2,1)

C  Now multiply the two together and take the Inverse FT
	        DO j=1,mdttl
	          DO i=1,mdttl
	            WFbuff(i,j)=WFbuff(i,j)*DWF(i,j)
	          END DO
	        END DO
	        CALL DFOURN(WFbuff,NN,2,-1)

C  Shift and divide by mdttl squared
	        CALL CRSWAP2(MWFFT,WFbuff,mdttl)
	      END IF

C end blur function  3/11/99

c  pixilate intensity
	      IF (PropType((iCurWFElt-1)).NE.1) THEN
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	        nPixAdd,mPix,mdttl,dxWF,dxPixAdd,iPloTyp,ifPixLoc,ifAdd)
	      ELSE
c	        CALL CPIX_J(PixArray,MWFFT,XFP,YFP,drx,dry,
c     &	        nPixAdd,mdttl,dxpix,iPloTyp,RayPos,RayID,
c     &	        psiElt(1,iElt),VptElt(1,iElt),LRayOK,ifPixLoc,ifAdd)
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	        nPixAdd,mPix,mdttl,dxWF,dxPixAdd,iPloTyp,ifPixLoc,ifAdd)
	      END IF
	      CALL IntToChar(Cinteger,iElt,i)

	      ifPixArray=.TRUE.
	      CALL CACCEPT(ANS,'YES',
     &        	           'Image added. Plot composed image? ')

	      IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	      IF (StretchType.EQ.2) THEN
	        ANS='LOG10 Composed Image#%'
	      ELSE IF (StretchType.EQ.3) THEN
	        ANS='SQRT Composed Image#%'
	      ELSE
	        ANS='Composed Image#%'
	      END IF
	        CALL PixOut(PixArray,iElt,iPloTyp,ANS,
     &	        'X-Axis','Y-Axis','.add',nPixAdd,mpix,R1,R2,StretchType)
	      END IF
	    END IF

C	  ***next command: Calculate system modular transfer function (MTF) ***
#if 1
	  ELSE IF (LCMP(command,'CMTF',4)) THEN
	     IF (.NOT.ifLoad) then
               WRITE(*,150); GO TO 1
	     END IF

	     ifInMTFCalc=.true.  ! to be used in CPROPAGATE in propsub.F

	     if (.not. ifEFLCalc) then
#include "efl.inc"
	       ifEFLCalc=.true.
	     end if
	     print*,' Computing 2D MTF ...'

	     ! Propagate ray beam through system
	     iElt=nElt
	     CALL RunProp(iElt,npts,ifPol,ifGrid,ifPropagate,
     &         ifBuild,ifTrace,ifRegrid,ifVecDif3,ifIntSrfSet,ifLNegOK,
     &         iEndElt,iCurWFElt,iCurRayElt,iWF,Dicr,Djcr,dxSource,
     &         iObsOpt,nBadRays,StartIndRef,StartExtinc,DWF,OPDMat)
	     !print*,'**After RunProp: WFElt(64,64,iWF)=',WFElt(64,64,iWF)

	     ! Compute MTF by autocorrelating WFElt (complex field)
	     Call CalcMTF(WFElt(1:mdttl,1:mdttl,iWF),
     &                    MTFMat2(1:2*mdttl,1:2*mdttl),mdttl)
	     !print*,'**After calling CalcMTF()'

	     if (.false.) then
	       ! Test only
	       Open(Unit=15,File='MTFData.txt')
	       Do j=1,2*mdttl
	         Write(15,*) MTFMat2(1:2*mdttl,j)
	       End Do
	       Close(15)
	     end if

	     ! Important to set it to 'false' at the end of 'MTF' command
	     ifInMTFCalc=.false.
#endif

C         ***next command: Set plot stretch***

	  ELSE IF (LCMP(command,'STRETCH',3)) THEN
	    CALL CACCEPT(ANS,'LINEAR',
     &	    'Enter image stretch type (LINEAR, LOG10, SQRT):')
	    IF (LCMP(ANS,'LOG10',3)) THEN
	      StretchType=2
	    ELSE IF (LCMP(ANS,'SQRT',3)) THEN
	      StretchType=3
	    ELSE
	      StretchType=1
	    END IF

C         ***Set noise seed***

	  ELSE IF (LCMP(command,'SEED',3)) THEN
	    CALL IACCEPT_S(k,12345,1,
     &	      'Enter noise seed (any integer 1 to 32,700):')
	    CALL RNGIN(k,NoiseSeed)

C         ***Add noise to composed image***

	  ELSE IF (LCMP(command,'NOISE',3)) THEN
	    IF (.NOT.ifAdd) THEN
	      WRITE(*,*)' Must COMpose before adding noise.'
	      GO TO 1
	    END IF
	    CALL CACCEPT(ANS,'PHOTON',
     &	      'Enter noise type (PHOton, BACkground, REAd, GAUssian):')
	    IF (LCMP(ANS,'PHOTON',3)) THEN
	      DO  j=1,nPixAdd
	        DO  i=1,nPixAdd
                  PixArray(i,j)=POIDEV(PixArray(i,j),NoiseSeed(1))
	        END DO
	      END DO
	    ELSE IF (LCMP(ANS,'BACKGROUND',3)) THEN
	      CALL RACCEPT_S(expval,0d0,1,
     &	      'Enter expected number of background noise photons'//
     &	      ' per pixel:')
	      DO  j=1,nPixAdd
	        DO  i=1,nPixAdd
	          S2=POIDEV(expval,NoiseSeed(1))
                  PixArray(i,j)=PixArray(i,j)+S2
	        END DO
	      END DO
	    ELSE IF (LCMP(ANS,'POISSON',3)) THEN
	      CALL RACCEPT_S(expval,0d0,1,
     &	      'Enter expected number of noise photo-electrons:')
	      CALL DACCEPT_S(S2,1d0,1,
     &	      'Enter detector quantum efficiency:')
	      DO  j=1,nPixAdd
	        DO  i=1,nPixAdd
                  PixArray(i,j)=PixArray(i,j)+
     &	          POIDEV(expval,NoiseSeed(1))/S2
	        END DO
	      END DO
	    ELSE IF (LCMP(ANS,'READ',3)) THEN
	      CALL DACCEPT_S(S1,0d0,1,
     &	      'Enter standard deviation of noise (photo-electrons):')
	      CALL DACCEPT_S(S2,1d0,1,
     &	      'Enter detector quantum efficiency:')
	      DO  j=1,nPixAdd
	        DO  i=1,nPixAdd
                  PixArray(i,j)=PixArray(i,j)+S1*GAUSS(NoiseSeed)/S2
	        END DO
	      END DO
	    ELSE IF (LCMP(ANS,'GAUSSIAN',3)) THEN
	      CALL DACCEPT_S(S1,0d0,1,
     &	      'Enter flux standard deviation:')
	      DO  j=1,nPixAdd
	        DO  i=1,nPixAdd
                  PixArray(i,j)=PixArray(i,j)+S1*GAUSS(NoiseSeed)
	        END DO
	      END DO
	    ELSE
	      WRITE(*,*)' Unknown noise type.'
	    END IF

C         ***Set parameters for composed images***

	  ELSE IF (LCMP(command,'COMPOSE',3)) THEN
	    CALL IACCEPT_S(iComElt,iCurWFElt,1,
     &	      'Enter element where image is to be composed:')
	    iElt=iComElt
	    IF ((iComElt.LE.0).OR.(iComElt.GT.nElt)) THEN
	      WRITE (*,*)' COMPOSE: Invalid element number: ',iElt
	      GO TO 1
	    ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.11).OR.
     &	    	     (EltID(iElt).EQ.12)) THEN
	      WRITE(*,*)' COMPOSE: Invalid element type =',EltID(iElt)
	      GO TO 1
	    END IF
	    IF (ifAdd) THEN
	      WRITE(*,*)' Overwriting previous COMPOSED image...'
	    END IF

	    ifAdd=.TRUE.
	    ifPixArray=.FALSE.
	    i=mPix
	    CALL IACCEPT_S(nPixAdd,i,1,'Enter number of pixels per side:')
	    IF (nPixAdd.GT.mPix) THEN
	      nPixAdd=mPix
	      WRITE(*,*)' WARNING: Number of pixels reset to', nPixAdd
	    END IF
	    dxWF=1d-3
	    CALL RACCEPT_S(dxPixAdd,dxWF,1,'Enter size of pixel:')
	    DO 52 j=1,mPix
	      DO 51 i=1,mPix
	  	PixArray(i,j)=0e0
 51	      CONTINUE
 52	    CONTINUE

C         ***Read filter data file***

	  ELSE IF (LCMP(command,'RFILT',3)) THEN
	    IF (ifLoad) THEN
	      CALL CFiltRead(ifFilterData)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Load new filter data***

	  ELSE IF (LCMP(command,'NFILT',3)) THEN
	    IF (ifLoad) THEN
	      CALL CFiltLoad(ifFilterData)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Save filter data***

	  ELSE IF (LCMP(command,'SFILT',3)) THEN
	    IF (ifLoad) THEN
	      CALL CFiltSave(ifFilterData)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Plot multi-spectral image***

	  ELSE IF (LCMP(command,'MULTISPEC',3)) THEN

	    IF (nGridPts.EQ.1) THEN
	      WRITE(*,*)' Too few rays, reset nGridPts.'
	      GO TO 1
	    END IF
	    IF (.NOT.ifLoad) THEN
	      WRITE(*,150)
	      GO TO 1
	    END IF
	    IF (.NOT.ifFilterData) THEN
	      WRITE(*,*)' Load filter data first.'
	      GO TO 1
	    END IF

	    IF ((iCurWFElt.LE.0).OR.(iCurWFElt.GT.nElt)) THEN
	      i=nElt
	    ELSE
	      i=iCurWFElt
	    END IF
	    CALL IACCEPT_S(iElt,i,1,
     &	    'Enter number of element where data is to be generated:')
	    IF ((iElt.LE.0).OR.(iElt.GT.nElt)) THEN
	      WRITE (*,*)'MULTISPEC: Invalid element number: ',iElt
	      GO TO 1
	    ELSE IF ((EltID(iElt).EQ.9).OR.(EltID(iElt).EQ.11).OR.
     &	    (EltID(iElt).EQ.12)) THEN
	      WRITE(*,*)'Invalid element type =',EltID(iElt)
	      GO TO 1
	    END IF
	    ifPixArray=.TRUE.
	    iComElt=iElt

	    IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.5).AND.(iPloTyp.NE.6)
     &      	.AND.(iPloTyp.NE.8).AND.(iPloTyp.NE.9)
     &	    	.AND.(.NOT.ifPlot)) THEN
	      ifPlot=.TRUE.
	      CALL GRAINI
	    END IF

	    CALL IntToChar(Cinteger,nStarType,i)
	    CALL IACCEPT_S(iStarType,1,1,
     &	    'Enter StarType (1-'//Cinteger(1:i)//'): ')
	    IF (iStarType.GT.nStarType) THEN
	      WRITE(*,*)' Invalid star type.'
	      GO TO 1
	    END IF

	    IF (ifAdd) THEN
	      WRITE(*,*)' Overwriting previous COMPOSED image...'
	    END IF

	    ifAdd=.TRUE.
	    i=mPix
	    CALL IACCEPT_S(nPixAdd,i,1,'Enter number of pixels per side:')
	    IF (nPixAdd.GT.mPix) THEN
	      nPixAdd=mPix
	      WRITE(*,*)' WARNING: Number of pixels reset to', nPixAdd
	    END IF
	    dxWF=1d-3
	    CALL RACCEPT_S(dxPixAdd,dxWF,1,'Enter size of pixel:')
	    DO 852 j=1,mPix
	      DO 851 i=1,mPix
	  	PixArray(i,j)=0e0
 851	      CONTINUE
 852	    CONTINUE

	    iEm1=iElt-1
	    BaseFlux=Flux
	    BaseWavelen=Wavelen
	    DO 853 iWavelen=1,nWavelen
	      Flux=FiltFlux(iWavelen,iStarType)*BaseFlux
	      Wavelen=FiltWavelen(iWavelen)

	      ifPropagate=.FALSE.
	      CALL RunProp(iElt,npts,ifPol,ifGrid,ifPropagate,
     &	       ifBuild,ifTrace,ifRegrid,ifVecDif3,ifIntSrfSet,ifLNegOK,
     &	       iEndElt,iCurWFElt,iCurRayElt,iWF,Dicr,Djcr,dxSource,
     &	       iObsOpt,nBadRays,StartIndRef,StartExtinc,DWF,OPDMat)

	      du=dxElt(iEm1)
	      udia=du*FLOAT(mdttl-1)
	      dx=dxElt(iElt)
	      xdia=FLOAT(mdttl-1)*dx

	      WRITE(*,206)Wavelen,Flux,udia,du,xdia,dx

	      X0PFP=0d0-dx*FLOAT((mdttl)/2)
	      DO 834 i=1,mdttl
	        XFP(i)=X0PFP
	        YFP(i)=X0PFP
 834	      CONTINUE

	      i=mdttl/2
	      dxWF=DABS(dxElt(iElt))
	      CALL Ca2Int(ifVecDif3,Power,MaxInt,imax,jmax)
	      WRITE(*,418)MaxInt,imax,jmax,Power
	      S1=dxWF/2d0
	      IF (dxPixAdd.LT.S1) WRITE(*,*)' WARNING: PIXELS TOO SMALL'
	      IF (ifPixLoc) THEN
	        CALL DSUB(DV1,VptElt(1,iElt),RayPos(1,1),3)
	        drx=Tcompact(3,4)*DV1(1)+Tout(3,5)*DV1(2)+
     &	        Tout(3,6)*DV1(3)+(WinCenPix(1)-EltPix(1))*sizPix
	        dry=Tcompact(4,4)*DV1(1)+Tout(4,5)*DV1(2)+
     &	        Tout(4,6)*DV1(3)+(WinCenPix(2)-EltPix(2))*sizPix
	      ELSE
	        drx=0.0
	        dry=0.0
	      END IF
	      IF (PropType((iCurWFElt-1)).NE.1) THEN
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	        nPixAdd,mPix,mdttl,dxWF,dxPixAdd,iPloTyp,ifPixLoc,ifAdd)
	      ELSE
c	        CALL CPIX_J(PixArray,MWFFT,XFP,YFP,drx,dry,
c     &	        nPixAdd,mdttl,dxpix,iPloTyp,RayPos,RayID,
c     &	        psiElt(1,iElt),VptElt(1,iElt),LRayOK,ifPixLoc,ifAdd)
	        CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,drx,dry,
     &	        nPixAdd,mPix,mdttl,dxWF,dxPixAdd,iPloTyp,ifPixLoc,ifAdd)
	      END IF
	      CALL IntToChar(Cinteger,iElt,i)
 853	    CONTINUE

	    Flux=BaseFlux
	    Wavelen=BaseWavelen
	    IF (StretchType.EQ.2) THEN
	      ANS='LOG10 Image, Filter='//
     &	      filtfil(1:(ICLEN(filtfil)-4))//'#%'
	    ELSE IF (StretchType.EQ.3) THEN
	      ANS='SQRT Image, Filter='//
     &	      filtfil(1:(ICLEN(filtfil)-4))//'#%'
	    ELSE
	      ANS='Image, Filter='//filtfil(1:(ICLEN(filtfil)-4))//'#%'
	    END IF
	    CALL PixOut(PixArray,iElt,iPloTyp,
     &	     'Multispectral Image, Filter='//
     &	     filtfil(1:(ICLEN(filtfil)-4))//'#%',
     &	     'X-Axis','Y-Axis','.add',nPixAdd,mpix,R1,R2,StretchType)

C         ***Compute perturbed system parameters from entered data***

	  ELSE IF (LCMP(command,'PERTURB',3)) THEN
	    IF (ifLoad) THEN
61	      CALL IACCEPT_S(iElt,iCurWFElt,1,
     &	                     'Enter element to be perturbed:')
	      IF ((iElt.LT.0).OR.(iElt.GT.nElt)) THEN
	        WRITE (*,*)'PERTURB: Invalid element number: ',iElt
	        !GO TO 61
	      END IF

	      IF ((iElt.EQ.0).AND.(.NOT.ifStopSet)
     &            .AND.(.NOT.RxNoStopSet_Flg)) THEN
	        WRITE(*,*)
     &	    ' Must set STOP before perturbing source. Use STOP command.'
	        GO TO 1
	      END IF

	      CALL CPERTURB(iElt)

	      ifTrace=.FALSE.
	      ifBuild=.FALSE.
	      ifPropagate=.FALSE.
	      ifPerturb=.FALSE.
	      ifGrid=.FALSE.
	      iCurRayElt=0

	      If (EltStopSet .AND. (iElt<=StopElt)
     &            .AND. LStopAtEltAutoSet) Then
                ! Set aperture stop at element after perturb
                WRITE(*,*)' '
                WRITE(*,315)StopElt,StopOffset
                GO TO 210
              End If
	    ELSE
	      WRITE(*,150)
	    END IF


C        *** Group PERTURB, supported in both MACOS and SMACOS

	  ELSE IF (LCMP(command,'GPERTURB',4)) THEN
	    IF (ifLoad) THEN
 63	      CALL IACCEPT_S(iElt,iCurWFElt,1,
     &                       'Enter element to be perturbed:')
              IF ((iElt.LE.0).OR.(iElt.GT.nElt)) THEN
                WRITE (*,*)'GPERTURB: Invalid element number: ',iElt
                GO TO 63
              END IF

	      If (EltGrp(0,iElt)>0) Then
	        Call CPERTURB_GRP_DVR(iElt)
	      Else
	        WRITE(*,*)' No element group defined at element ',iElt
	      End If

	      ifTrace=.FALSE.
              ifBuild=.FALSE.
              ifPropagate=.FALSE.
              ifPerturb=.FALSE.
              ifGrid=.FALSE.
              iCurRayElt=0
	    ELSE
              WRITE(*,150)
            END IF

C         ***Compute perturbed system parameters from file data***

	  ELSE IF (LCMP(command,'PREAD',3)) THEN
	    IF (ifLoad) THEN
	      ifPertRead=.FALSE.
	      iPertElt=nElt
	      CALL CPRead(PertVec,ifPertRead,iPertElt)
	    ELSE
	      WRITE(*,150)
	    END IF

C         ***Reset linear model perturbation vector***

	  ELSE IF (LCMP(command,'LRESET',3)) THEN
	    DO 85 i=1,6
	      DO 84 j=1,nElt
	          PertVec(i,j)=0d0
 84	      CONTINUE
 85	    CONTINUE

C         ***Perturb linear model using file input***

	  ELSE IF (LCMP(command,'LPREAD',3)) THEN
	    IF (ifBuild.OR.ifCmt) THEN
	      ifPertRead=.FALSE.
	      CALL CLPRead(PertVec,ifPertRead,iPertElt)
	      IF (ifPertRead) ifPerturb=.TRUE.
	    ELSE
	      WRITE(*,152)
	    END IF

C         ***Perturb linear model***

	  ELSE IF (LCMP(command,'LPERTURB',3)) THEN
	    IF (ifBuild.OR.ifCmt) THEN
62	      CALL IACCEPT_S(iElt,iCurWFElt,1,
     &	      'Enter element to be perturbed:')
	      IF ((iElt.LE.0).OR.(iElt.GT.nElt)) THEN
	        WRITE (*,*)'LPERTURB: Invalid element number.'
	        GO TO 62
	      END IF

	      CALL CLPERTURB(PertVec,RaySpot,
     &        OPDMat,AvgL,Cmatrix,npts,RayIndex,iPertElt,iElt)
	      ifPerturb=.TRUE.
	    ELSE
	      WRITE(*,152)
	    END IF

C         ***Spot diagram from linear model***

	  ELSE IF (LCMP(command,'LSPOT',3)) THEN
	    IF (ifPerturb) THEN
	      IF ((iPloTyp.NE.0).AND.(.NOT.ifPlot)) THEN
	        ifPlot=.TRUE.
	        CALL GRAINI
	      END IF
	      CALL CLTRACE(PertVec,RaySpot,
     &        OPDMat,AvgL,Cmatrix,npts,RayIndex,iPertElt,iElt)
	      NBEAMS=2
	      CALL SpotDiag(RaySpot,nRay,2,mRay,2,1,
     &	                    'X-Axis','Y-Axis','Spot Diagram')
	    ELSE
	      WRITE(*,157)
	    END IF

C         ***OPD map from linear model***

	  ELSE IF (LCMP(command,'LOPD',3)) THEN
	    IF (ifPerturb) THEN
	      IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.5).AND.(iPloTyp.NE.6)
     &        .AND.(iPloTyp.NE.8).AND.(iPloTyp.NE.9)
     &        .AND.(.NOT.ifPlot)) THEN
	        ifPlot=.TRUE.
	        CALL GRAINI
	      END IF
	      CALL CLTRACE(PertVec,RaySpot,
     &        OPDMat,AvgL,Cmatrix,npts,RayIndex,iPertElt,iElt)
	      dx=dxElt(iLinElt)
	      X0PFP=RaySpot(1,1)-dx*FLOAT(nGridPts/2)
	      Y0PFP=RaySpot(1,2)-dx*FLOAT(nGridPts/2)
	      DO 13 i=1,nGridPts
	        XFP(i)=X0PFP
	        YFP(i)=Y0PFP
	        X0PFP=X0PFP+dx
	        Y0PFP=Y0PFP+dx
 13	      CONTINUE
	      CALL DSRFPLT(iPloTyp,XFP,YFP,OPDMat,nGridPts,nGridPts,
     &	      mpts,mpts,'X-Axis','Y-Axis','LinMod OPD Map',
     &	      R1,R2,1)
	    ELSE
	      WRITE(*,157)
	    END IF

C         ***Linear WF plots***

	  ELSE IF (LCMP(command,'LIN',3).OR.LCMP(command,'LPI',3)
     &	  .OR.LCMP(command,'LLO',3)) THEN
	    IF (ifPerturb) THEN
	      IF ((iPloTyp.NE.0).AND.(iPloTyp.NE.5).AND.(iPloTyp.NE.6)
     &        .AND.(iPloTyp.NE.8).AND.(iPloTyp.NE.9)
     &        .AND.(.NOT.ifPlot)) THEN
	        ifPlot=.TRUE.
	        CALL GRAINI
	      END IF
	      CALL CLTRACE(PertVec,RaySpot,
     &        OPDMat,AvgL,Cmatrix,npts,RayIndex,iPertElt,iElt)
	      CALL CLPLOT(command,npts,RaySpot,
     &        OPDMat,AvgL,RayIndex,RayMap,
     &	      iPloTyp,iPertElt,CD1,CD2,PixArray,StretchType)
	    ELSE
	      WRITE(*,157)
	    END IF

C         ***Set surface plot type***

	  ELSE IF (LCMP(command,'NULL',3).OR.LCMP(command,'NONE',3))
     &	  THEN
	    iPloTyp=0
	  ELSE IF (LCMP(command,'SLICE',3)) THEN
	    iPloTyp=1
	    WRITE(*,*)' Plot type set to SLICE'
	  ELSE IF (LCMP(command,'GRAY',2)) THEN
	    iPloTyp=2
	    WRITE(*,*)' Plot type set to GRAY'
	  ELSE IF (LCMP(command,'CONTOUR',3)) THEN
	    iPloTyp=3
	    WRITE(*,*)' Plot type set to CONTOUR'
	  ELSE IF (LCMP(command,'WIRE',3)) THEN
	    iPloTyp=4
	    WRITE(*,*)' Plot type set to WIRE'
	  ELSE IF (LCMP(command,'TEXT',3)) THEN
	    iPloTyp=5
	    WRITE(*,*)' Plot type set to TEXT'
	  ELSE IF (LCMP(command,'BINARY',3)) THEN
	    iPloTyp=6
	    WRITE(*,*)' Plot type set to BINARY'
	  ELSE IF (LCMP(command,'COLUMN',3)) THEN
	    iPloTyp=7
	    WRITE(*,*)' Plot type set to COLUMN'
	  ELSE IF (LCMP(command,'FITS',3)) THEN
	    iPloTyp=8
	    WRITE(*,*)' Plot type set to FITS'
	  ELSE IF (LCMP(command,'MAT',3)) THEN
	    iPloTyp=9
	    WRITE(*,*)' Plot type set to MAT'
	  ELSE IF (LCMP(command,'ROW',3)) THEN
	    iPloTyp=10
	    WRITE(*,*)' Plot type set to ROW'
c***loh
C
C         ***Set batch mode***

	  ELSE IF (LCMP(command,'LQW',2).OR.LCMP(command,'EXECUTE',3))
     &	  THEN
c sws
c used to be toggle, now just .jou file read
C            Go into batch mode
	     CALL MB_LOH(1)
c push the stack
             pstack= pstack+1
             read_len(pstack) = 0
             read_cur(pstack) = 0
c make sure we're in batch
             bat_task=2
	     GO TO 1
c sws

C         ***Toggle writing journal file***

          ELSE IF (LCMP(command,'JOURNAL',3)) THEN
            IF (jou_task .EQ. 1) THEN

C             Go into batch mode

              CALL JOU_LOH(1)
              jou_task=2
            ELSE IF (jou_task .EQ. 2) THEN

C             Go into interactive mode

              jou_task=1
              CALL JOU_LOH(2)
            END IF

#ifdef MACOS_CMD
          ELSE IF (LCMP(command,'WFITS',5)) THEN
	    CALL CACCEPT(ANS,'OPD','Save OPD or Intensity (OPD,INT):')
	    If (LCMP(ANS,'OPD',3)) THEN
              If (firstOPDCalc .OR. multiOPDCalc) Then
	        Write(*,*)'  Dumping current OPD into FITS file Opd.fits'
	        FitsFilename='Opd.fits'
	        Call fitswrite(FitsFilename,OPDMat,npts+1)
              Else
	        Write(*,*)'  OPD data is not available, run OPD command'
	      End If
	    Else If (LCMP(ANS,'INT',3)) THEN
	      If (firstIntCalc) Then
	        Write(*,*)'  Dumping current intensity into FITS file Int.fits'
	        FitsFilename='Int.fits'
	        Call fitswrite(FitsFilename,MWFFT,mdttl);
              Else
	        Write(*,*)'  Intensity data is not available, run INt command'
	      End If
	    Else
	      Write(*,*)'  Invalid answer: ',ANS
	    End If

	  ELSE IF (LCMP(command,'PWIN',4)) THEN
	    ! Start a Python-Tcl/Tk GUI window
            ANS='!he.py&'
  	    Call SHELLCMD(ANS,7)
#endif


c***loh

#ifdef SMACOS_CMD
C       ***Get matrix from matlab

          ELSE IF (LCMP(command,'GETMATLABMATRIX',4)) THEN

cwr         CALL RBININ(iPloTyp,R1,row_sz,col_sz)
#endif

C *************************************************
C  For testing srtrace.F only, -jzlou
	  ELSE IF (LCMP(command,'SRTRACE',6)) THEN
	    IF (ifLoad) THEN
	      print*,'Calling SRTRACE_Test ...'
	      CALL SRTRACE_Test(0,nElt,ChfRayPos,ChfRayDir)
	      print*,'After Calling SRTRACE_Test ...'
	    ELSE
	      WRITE(*,150)
	    END IF
C *************************************************

C         ***Here if no match***

	  ELSE

#ifdef MACOS_CMD
	    IF (bat_task.NE.2) WRITE(*,302)
#else
	    IF (ifEcho) WRITE(*,302)
#endif
	    print*,'** Unknown command = ', command
	  END IF
	GO TO 1



!
! Defined in prop_defs.inc, corresponding to EltID array. -jzlou
!
! 	EltTypeName(1)='Reflector'
!       EltTypeName(2)='FocalPlane'
!       EltTypeName(3)='Reference'
!       EltTypeName(4)='HOE'
!       EltTypeName(5)='Grating'
!       EltTypeName(6)='Refractor'
!       EltTypeName(7)='Obscuring'
!       EltTypeName(8)='Return'
!       EltTypeName(9)='NSRefractor'
!       EltTypeName(10)='LensArray'
!       EltTypeName(11)='Segment'
!       EltTypeName(12)='NSReflector'
!       EltTypeName(13)='TrGrating'
!       EltTypeName(14)='RfPolarizer'
!       EltTypeName(15)='TrPolarizer'
