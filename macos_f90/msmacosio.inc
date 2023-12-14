C Shared input deck code between macosio.F and smacosio.F
C -jzlou, 03/2006

C       Set entrance condition flags to .FALSE.

          CALL LOHIN1()

  50      CONTINUE   ! This section is for header part of Rx

          READ(RX_FID, 123, END=99) STRING
 354      CONTINUE

          CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,RX_FID)

          ! --------------------------------------------------------------------------------
	  IF (LCMP(VAR_NAM,'JZLouChfRayDir',9)) THEN  ! disabled
            READ(VALUE,*)(ChfRayDir(i),i=1,3)
            IF (DMAG(ChfRayDir(1),3) .EQ. 0.0D0) THEN
              WRITE(*,124)' ChfRayDir'
            ELSE
              CALL DUNITIZE(ChfRayDir)
              ChfRayDir_FLG=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'JZLouChfRayPos',9)) THEN  ! disabled
            READ(VALUE,*)(ChfRayPos(i),i=1,3)
            ChfRayPos_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'xFrame',6)) THEN
            READ(VALUE,*) TElt(1:3,1,0)
            SrcLF_FLG=.TRUE.;
	    print*,' **** msmacosio: src xFrame read'
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'yFrame',6)) THEN
            READ(VALUE,*) TElt(1:3,2,0)
            SrcLF_FLG=.TRUE.
	    print*,' **** msmacosio: src yFrame read'
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'zFrame',6)) THEN
            READ(VALUE,*) TElt(1:3,3,0)
            SrcLF_FLG=.TRUE.
            CalcSrcZFrame_FLG=.FALSE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OPDRefRayLen',12)) THEN
            READ(VALUE,*) OPDRefRayLen
            OPDRefRayLen_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'FEXCentroid',11)) THEN
            IF (.NOT.LCMP(VALUE,'Y',1)) Rx_FEXCentrFlg=.False.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'SaveVis3dDat',12)) THEN
            IF (LCMP(VALUE,'Y',1)) SaveVis3d_dat_Flg=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'ShowMetData',11)) THEN
	    IF (LCMP(VALUE,'Y',1)) ShowMetData_FLG=.TRUE.
	    GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'ObjFile',7)) THEN
	    ObjIntFn=' '
	    READ(VALUE,*) ObjIntFn
	    !print*,'ObjIntFn = ',ObjIntFn(1:StrLen(ObjIntFn))
	    Call ReadDblArray2D(ObjIntFn,mdttl)
	    ifObjImgSet=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'ApStop',6)) THEN
            READ(VALUE,*)(StopPos(i),i=1,3)
            ifStopSet=.TRUE.; RxStopSet=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'RxNoStopSet',11)) THEN
            IF (LCMP(VALUE,'Y',1)) Then
              RxNoStopSet_Flg = .TRUE.
            Else
              RxNoStopSet_Flg = .FALSE.
            End If
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'RayTgtElt',9)) THEN
	    LRayTgtElt=.TRUE.
            READ(VALUE,*) RayTgtElt
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'Save3DApVec',9)) THEN
            IF (LCMP(VALUE,'Y',1)) THEN
              ifSave3dApVec=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'SaveOPDMap',10)) THEN
	    IF (LCMP(VALUE,'Y',1)) THEN
              saveOPD_Flg=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'SaveRayPosH',11)) THEN
            IF (LCMP(VALUE,'Y',1)) THEN
    	      SaveRayPosHist_Flg=.TRUE.
	      print*,'** msmacosio: SaveRayPosHist is true ...'
	    END IF
 	    GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'DumpRayPosH',11)) THEN
            IF (LCMP(VALUE,'Y',1)) THEN
              DumpRayPosHist_Flg=.TRUE.
              print*,'** msmacosio: DumpRayPosHist is true ...'
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'CalBWK',6)) THEN
	    IF (LCMP(VALUE,'Y',1)) THEN
	      BWK_FLG=.TRUE.
	    END IF
	    GO TO 50

#ifdef CMACOS
          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'PrtEltRots',6)) THEN
            IF (LCMP(VALUE,'Y',1)) THEN
              PrtPrinRots=.TRUE.
            END IF
            GO TO 50
#endif

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptTarget',7)) THEN
!    &             .OR. LCMP(VAR_NAM,'OptTgt',6)) THEN
	    IF (LCMP(VALUE,'WFE',3)) THEN
	      OptTarget=WFE_TARGET
	    ELSE IF (LCMP(VALUE,'SPOT',3)) THEN
              OptTarget=SPOT_TARGET
	    ELSE IF (LCMP(VALUE,'ZWF',3)) THEN
	      OptTarget=WFE_ZMODE_TARGET
	    ELSE IF (LCMP(VALUE,'Beam',4)) THEN
	      OptTarget=BEAM_TARGET
	    ELSE IF (LCMP(VALUE,'OPL',3)) THEN
	      OptTarget=OPL_TARGET
	    END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptSpotSize',11)) THEN
            READ(VALUE,*)OptSpotSize
            if (OptSpotSize<0) OptSpotSize=0d0
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptRayGrid',6)) THEN
	    opt_npts_flg=.true.
	    READ(VALUE,*) opt_npts
	    IF (opt_npts>mpts) THEN
              WRITE(*,*)
     &          ' Too many grid points for optimization, resetting to',mpts
	      opt_npts=mpts-1
	    ELSE
	      opt_npts=opt_npts-1
	    END IF
	    IF (opt_npts<1) opt_npts=15
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptZernModes',8)) THEN
	    READ(VALUE,*) nWFZern
	    READ(VALUE,*) nWFZern,WFZernMode(1:nWFZern)
	    GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'VarSrc',6)) THEN
	    isVarElt(0)=.TRUE.
	    nVarElt=nVarElt+1
	    varElts(nVarElt)=0
	    ! Read in three rotation and three translation
	    ! flags for source
	    READ(VALUE,*)(varEltDOF(j,0),j=1,mVarDOF-2)
	    nDOF_VarElt(nVarElt)=0
            DO i=1,mVarDOF
              IF (varEltDOF(i,0)/=0) THEN
                nDOF_VarElt(nVarElt)=nDOF_VarElt(nVarElt)+1
                DOF_VarElt(nDOF_VarElt(nVarElt),nVarElt)
     &            = DOF_NameList(i)
              END IF
            END DO
	    GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptSrcRpt',9)) THEN
	    READ(VALUE,*)(OptSrcRpt(i),i=1,3)
	    OptSrcRpt_FLG=.TRUE.
	    GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'OptChfRayDir',12) .OR.
     &             LCMP(VAR_NAM,'ChfRayDir',9)) THEN
            nOptChfRayDir=nOptChfRayDir+1
            READ(VALUE,*) (opt_fov(i,1,nOptChfRayDir),i=1,3)
	    If (nOptChfRayDir==1) Then
	       ChfRayDir(1:3)=opt_fov(1:3,1,nOptChfRayDir)
	       IF (DMAG(ChfRayDir(1),3) .EQ. 0.0D0) THEN
                 WRITE(*,124)' ChfRayDir'
                 GO TO 99
               ELSE
                 CALL DUNITIZE(ChfRayDir)
                 ChfRayDir_FLG=.TRUE.
               END IF
	    End If
            if (nOptChfRayDir==nOptChfRayPos)
     &        nOptFov=nOptChfRayDir
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'OptChfRayPos',12).OR.
     &             LCMP(VAR_NAM,'ChfRayPos',9)) THEN
            nOptChfRayPos=nOptChfRayPos+1
            READ(VALUE,*) (opt_fov(i,2,nOptChfRayPos),i=1,3)
	    If (nOptChfRayPos==1) Then
	      ChfRayPos(1:3)=opt_fov(1:3,2,nOptChfRayPos)
	      ChfRayPos_FLG=.TRUE.
	    End If
	    !print*,'ChfRayPos(1:3) =',ChfRayPos(1:3)
            if (nOptChfRayPos==nOptChfRayDir)
     &        nOptFov=nOptChfRayPos
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptWFElt',8).OR.
     &             LCMP(VAR_NAM,'OptTgtElt',9)) THEN
            READ(VALUE,*) OptTgtElt
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptTgtWF',8)) THEN
#if 0
	    !print*,' **OptTgtWF: VALUE =',VALUE(1:ICLEN(VALUE))
	    print*,'  **In OptTgtWF'
	    tgtWF_fn(1:ICLEN(VALUE))=VALUE(1:ICLEN(VALUE))
	    If (LCMP(tgtWF_fn(ICLEN(tgtWF_fn)-2:ICLEN(tgtWF_fn)),
     &              'bin',3)) Then
              tgtWF_fmt(1:1)='B'
	    Else
              tgtWF_fmt(1:1)='T'
	    End If
#endif
	    tgtWF_fn(1:ICLEN(VALUE))=VALUE(1:ICLEN(VALUE))
	    !print*,' ** tgtWF_fn =',tgtWF_fn
	    If (LCMP(tgtWF_fn(ICLEN(VALUE):ICLEN(VALUE)),
     &          'T',1)) Then
	      tgtWF_fmt(1:1)='T'
	   Else
	      tgtWF_fmt(1:1)='B'
	    End If
            tgtWF_flg=.true.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'OptFEX',6)) THEN
            If (LCMP(VALUE,'N',1)) LOptIfFEX=.FALSE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptAlg',6)) THEN
            IF (LCMP(VALUE,'NonL',4)) THEN
              OptAlg=NonLin
	    ELSE IF (LCMP(VALUE,'Lin',3)) THEN
	      OptAlg=Linear
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptNomSens',9)) THEN
	    if (LCMP(VALUE,'Y',1)) LOptNomSens=.TRUE.
	    GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptSavePinv',11)) THEN
            if (LCMP(VALUE,'Y',1)) save_pinv=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          !  if OptPinvFile is not specified in Rx, the default file 'dady0_saved.txt'
          !  will be used when 'OptUseSavedPinv' is true. 08/2009

          ELSE IF (LCMP(VAR_NAM,'OptPinvFile',11)) THEN
            PinvFile=VALUE(1:MacosValLen)
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'OptUseSavedPinv',15)) THEN
            if (LCMP(VALUE,'Y',1)) use_saved_pinv=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptSvdSvCut',9)) THEN
	    READ(VALUE,*) SvdSvCut
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptFOVWt',5)) THEN
            READ(VALUE,*)(fov_weight(i),i=1,nOptFov)
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'OptWaveLen',10).OR.
     &             LCMP(VAR_NAM,'WaveLen',7)) THEN
            nOptWavelen=nOptWavelen+1
            READ(VALUE,*) opt_wavelen(nOptWavelen)
	    If (nOptWavelen==1) Then
	      Wavelen=opt_wavelen(nOptWavelen)
	      Wavelen_FLG=.TRUE.
	    End If
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'ArrWaveLen',10)) THEN
	    If (.not.Wavelen_FLG) Then
	      Call Get_Values(VALUE,opt_wavelen(1:max_wl),nTok)
	    Else
	      Call Get_Values(VALUE,
     &                 opt_wavelen(nOptWavelen+1:max_wl),nTok)
	    End If
	    nOptWavelen=nOptWavelen+nTok
	    If (.not.Wavelen_FLG) Then
	      Wavelen=opt_wavelen(1); Wavelen_FLG=.TRUE.
	    End If
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'OptMaxItrs',10)) THEN
            READ(VALUE,*) nitrs_dopt
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'ZCOZernType',10)) THEN
            ! zcoType is default to ZernCircular
            ! which uses NormNoll
            IF (LCMP(VALUE,'NormAnnular',10)) THEN
              zcoType=ZernAnnular
            Else If(LCMP(VALUE,'NormHex',7)) THEN
              zcoType=ZernHex
            End If
            GO TO 50

#ifdef CMACOS
          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'PgplotImage',11)) THEN
	    IF (LCMP(VALUE,'Color',5)) THEN
	      ifPGColor=.TRUE.
	    END IF
            GO TO 50
#endif

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'UseChfRay4OPD',10)) THEN
            If (LCMP(VALUE,'N',1)) Then
              LUseChfRayIfOK=.FALSE.  ! default to TRUE
            End If
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'zSource',7)) THEN
            READ(VALUE,*)zSource
            IF (ABS(zSource) <= EPS) GOTO 99
            zSource_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'IndRef',6)) THEN
            READ(VALUE,*)IndRef(0)
            IndRef_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'Extinc',6)) THEN
            READ(VALUE,*)Extinc(0)
            Extinc_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'BaseUnits',7)) THEN
            BaseUnits=VALUE
            BaseUnits_FLG=.TRUE.
            IF (LCMP(BaseUnits,'cm',2)) THEN
              CBM=1d-2; cUnit='cm'
            ELSE IF (LCMP(BaseUnits,'mm',2)) THEN
              CBM=1d-3; cUnit='mm'
            ELSE IF (LCMP(BaseUnits,'m',1)) THEN
              CBM=1d0; cUnit='m '
            ELSE IF (LCMP(WaveUnits,'in',2)) THEN
              CBM=2.54d-2; cUnit='in'
            ELSE IF (LCMP(WaveUnits,'none',2)) THEN
              CBM=1d0; cUnit='  '
            ELSE
              WRITE(*,*)' Unknown units, unitless operation assumed'
              BaseUnits='none'
              WaveUnits='none'
              CBM=1d0
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'WaveUnits',7)) THEN
            WaveUnits=VALUE
            WaveUnits_FLG=.TRUE.
            IF (LCMP(WaveUnits,'cm',2)) THEN
              DDEF(1)=6.328d-5
              CWM=1d-2
            ELSE IF (LCMP(WaveUnits,'mm',2)) THEN
              DDEF(1)=6.328d-4
              CWM=1d-3
            ELSE IF (LCMP(WaveUnits,'um',2)) THEN
              CWM=1d-6
            ELSE IF (LCMP(WaveUnits,'nm',2)) THEN
              CWM=1d-9
            ELSE IF (LCMP(WaveUnits,'m',1)) THEN
              CWM=1d0
            ELSE IF (LCMP(WaveUnits,'A',1)) THEN
              CWM=1d-10
            ELSE IF (LCMP(WaveUnits,'in',2)) THEN
              CWM=2.54d-2
            ELSE IF (LCMP(WaveUnits,'none',2)) THEN
              CWM=1d0
            ELSE
              WRITE(*,*)' Unknown units, unitless operation assumed'
              WaveUnits='none'
              CWM=1d0
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'JZLouWavelen',7)) THEN
            READ(VALUE,*)Wavelen
            Wavelen_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'Flux',4)) THEN
            READ(VALUE,*)Flux
            Flux_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'GridType',8)) THEN
            IF (LCMP(VALUE,GridTypeName(1),2).OR.
     &      LCMP(VALUE,'1',1)) THEN
              GridType=1  ! Circular
            ELSE IF (LCMP(VALUE,GridTypeName(2),2).OR.
     &      LCMP(VALUE,'2',1)) THEN
              GridType=2  ! Square
            ELSE IF (LCMP(VALUE,GridTypeName(3),2).OR.
     &      LCMP(VALUE,'3',1)) THEN
              GridType=3  ! Hex
            ELSE IF (LCMP(VALUE,GridTypeName(4),2).OR.
     &      LCMP(VALUE,'4',1)) THEN
              GridType=4  ! Pie
            ELSE IF (LCMP(VALUE,GridTypeName(5),2).OR.
     &      LCMP(VALUE,'5',1)) THEN
              GridType=5  ! Flower
            END IF
            GridType_FLG=.TRUE.
            IF (GridType .LT. 3) THEN
              nSeg_FLG=.TRUE.
              gap_FLG=.TRUE.
              width_FLG=.TRUE.
              nPetals_FLG=.TRUE.
              radCtr_FLG=.TRUE.
              SegCoord_FLG=.TRUE.
              nSeg=0
            ELSE IF ((GridType .EQ. 3).OR.(GridType .EQ. 4)) THEN
              nPetals_FLG=.TRUE.
              radCtr_FLG=.TRUE.
            ELSE IF (GridType .EQ. 5) THEN
              gap_FLG=.TRUE.
              width_FLG=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'Aperture',8)) THEN
            READ(VALUE,*)Aperture
            Aperture_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'Obscratn',8)) THEN
            READ(VALUE,*)Obscratn
            Obscratn_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'nGridpts',8)) THEN
            READ(VALUE,*) nGridpts
            npts=nGridpts-1
            nGridpts_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'xGrid',5)) THEN
            READ(VALUE,*)(xGrid(i),i=1,3)
            IF (DMAG(xGrid(1),3) .EQ. 0.0D0) THEN
              WRITE(*,124)' xGrid'
            ELSE
              xGrid_FLG=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'yGrid',5)) THEN
            READ(VALUE,*)(yGrid(i),i=1,3)
            IF (DMAG(yGrid(1),3) .EQ. 0.0D0) THEN
              WRITE(*,124)' yGrid'
            ELSE
              yGrid_FLG=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'drawXGrid',6)) THEN
            READ(VALUE,*)(UDxDrawGrid(i),i=1,3)
            LudDrawGrid=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'drawYGrid',6)) THEN
            READ(VALUE,*)(UDyDrawGrid(i),i=1,3)
            LudDrawGrid=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'UDBeam',6)) THEN
	    If (LCMP(VALUE,'Y',1)) Then
	      UDBeam_FLG=.TRUE.
	    End If
	    GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'UDSrcProf',9)) THEN
	    UDSrcProf=VALUE(1:MacosValLen)
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'nElt',4)) THEN
            READ(VALUE,*)nElt
            nElt_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'nSeg',4)) THEN
            READ(VALUE,*) nSeg
            nSeg_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'width',5)) THEN
	    READ(VALUE,*)width
            IF (GridType .NE. 5) THEN
              READ(VALUE,*)width
            ELSE
              READ(VALUE,*)nPetals
            END IF
            width_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'gap',3)) THEN
!#if 0
!            IF (GridType .NE. 5) THEN
!              READ(VALUE,*)gap
!            ELSE
!              READ(VALUE,*)radCtr
!            END IF
!#endif
            READ(VALUE,*)gap
            gap_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'nRings',5)) THEN
            READ(VALUE,*) nRings
            nRings_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'nPetals',5)) THEN
	    If (nRings_FLG) Then
              READ(VALUE,*) nPetals(0:nRings)
              nPetals_FLG=.TRUE.
	    Else
              Write(*,*)
     &        ' **Rx error: nRings must be specified before nPetals'
	      Stop
	    End If
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'radCtr',4)) THEN
            READ(VALUE,*)radCtr
            radCtr_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'RingRads',8)) THEN
            READ(VALUE,*) RingRads(1:nRings)
	    RingRads(0)=1d0
            RingRad_FLG=.TRUE.
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'SegCoord',8)) THEN
            IF (GridType .NE. 5) THEN
              READ(VALUE,*)(SegCoord(i,1),i=1,3)
              IF (nSeg_FLG) THEN
                IF (nSeg .GT. 1)  THEN
                  DO iSeg=2,nSeg
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(SegCoord(i,iSeg),i=1,3)
                  END DO
                END IF
                SegCoord_FLG=.TRUE.
              ELSE
                WRITE(*,*) '  Must input nSeg before SegCoord.'
              END IF
            ELSE
              READ(VALUE,*)SegCoord(1:2,1)
              IF (nSeg_FLG) THEN
                IF (nSeg .GT. 1)  THEN
                  DO iSeg=2,nSeg
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)SegCoord(1:2,iSeg)
                  END DO
                END IF
                SegCoord_FLG=.TRUE.
              ELSE
                WRITE(*,*) '  Must input nSeg before SegCoord.'
                GO To 99    ! abort reading Rx
              END IF
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'SegXgrid',8)) THEN
             READ(VALUE,*)(SegXgrid(i),i=1,3)
            IF (DMAG(SegXgrid(1),3) .EQ. 0.0D0) THEN
              WRITE(*,124)' SegXGrid'
            ELSE
              SegXgrid_FLG=.TRUE.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
	  ELSE IF (LCMP(VAR_NAM,'SegApType',9)) THEN
            ! Added by jzlou
            IF (.NOT.SegCoord_FLG) THEN
              print*,
     &          'macosio: SegCoord must be defined before SegApType'
              GO TO 99
            END IF
            IF (LCMP(VALUE,ApTypeName(6),3).OR.
     &          LCMP(VALUE,'6',1)) THEN
              ! Hexagonal aperture for segments
              SegApType=6; SegApType_FLG=.true.
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'SegApVec',8)) THEN
            ! Added by jzlou
            IF (SegApType_FLG) THEN
              IF (SegApType==6) THEN
                READ(VALUE,*)(SegApVec(j),j=1,2) ! hex h/2 and orient
              END IF
            END IF
            GO TO 50

          ! --------------------------------------------------------------------------------
          ! Header Exit Criteria
          ! --------------------------------------------------------------------------------
          ELSE IF (LCMP(VAR_NAM,'iElt',4)) THEN
            CALL ChkDf1(npts,ENTR_OK)
            IF (ENTR_OK) THEN
              GO TO 60
            ELSE
              WRITE(*,*)'  Input file not properly loaded'
              GO TO 99
            END IF
          ELSE
            GO TO 50
          END IF


          ! --------------------------------------------------------------------------------
          !
          ! Start the Element Definition Section
          !
          ! --------------------------------------------------------------------------------
          !
          ! Ignore "nElt" definition in Rx => updated at the end
          !

   60     CONTINUE

	  CWB=CWM/CBM
          WaveBU=CWB*Wavelen

          LApTypeDefined=.false.; LApVecDefined=.false.
	  LTaperedElt=.false.; LApSunBlades=.false.
	  LStopAtEltAutoSet=.false.

	  iSeg=0
          jGridSrf=0; iEltToGridSrf=0
          jAmplSrf=0; iEltToAmplSrf=0; ifAmplSrf=.FALSE.;
          jMetSrf=0; iEltToMetSrf=0; nMetPos=0
          tMetSrf=-1; tMetSrfPos=0

          nAsphCoef     = nAsphCoef_Default   ! reset to default
          AsphCoef(:,:) = 0d0

          nZernCoef     = nZernModes_Default  ! Default # of Zernike Coefficients to read in
          ZernCoef(:,:) = 0d0

	  iElt=0
          DO   ! this is a "DO ... END DO" loop

            ! update -- nElt were read in thus far
            !           note that duplicate Elements can only use iElt < nElt
            nElt = iElt

            ! prep for next Element reading
            iElt       = iElt + 1
	    DupElt_FLG = .FALSE.

            CALL LOHIN2()  ! set all elt flags to false. -jzlou

   61       CONTINUE
            READ(RX_FID, 123, END=99) STRING

  901       CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,RX_FID)
            IF (LCMP(STRING,'EOF',3)) GOTO 99

            ! ------------------------------------------------------------------------------
            IF (LCMP(VAR_NAM,'EltName',7)) THEN
              EltName(iElt)=VALUE
              EltName_FLG=.TRUE.

              Call RecordElementID(iElt,VALUE,.false.)
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'MRP',3)) THEN
              EltMRP(iElt)=VALUE
              LEltMRP_FLG(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'Element',7)) THEN
              TypeString=VALUE

              !  Assign EltID & set flags that aren't needed
              DO i=1,17
                IF (LCMP(TypeString,EltTypeName(i),6)) THEN
                  EltID(iElt)=i
	          If (i==17) Then
		    ! is a DOE element
	            nDOE=nDOE+1; DoeId(iElt)=nDOE
	          End If
                  EltID_FLG=.TRUE.
                END IF
              END DO
              IF (.NOT.EltID_FLG) THEN
                WRITE(*,130) TypeString
                GO TO 61
              END IF

              ! Set flags that aren't needed

	      IF ((EltID(iElt).EQ.2).OR.(EltID(iElt).EQ.3)
     &          .OR.(EltID(iElt).EQ.7).OR.(EltID(iElt).EQ.9)) THEN
                IndRef_FLG=.TRUE.
                Extinc_FLG=.TRUE.
                GlassCoef_FLG=.TRUE.
              END IF

              IF ((EltID(iElt).NE.4).AND.(EltID(iElt).NE.5).AND.
     &        (EltID(iElt).NE.13)) THEN
                h1HOE_FLG = .TRUE.
                h2HOE_FLG = .TRUE.
                WaveHOE_FLG = .TRUE.
                RuleWidth_FLG = .TRUE.
                OrderHOE_FLG = .TRUE.

              ELSE IF ((EltID(iElt).EQ.5).OR.(EltID(iElt).EQ.13)) THEN
                h2HOE_FLG = .TRUE.
                WaveHOE_FLG = .TRUE.

              ELSE
                RuleWidth_FLG = .TRUE.
              END IF

              IF (EltID(iElt).NE.10) THEN
                LensArrayWidth_FLG=.TRUE.
                LensArrayType_FLG=.TRUE.
              END IF

              IF (EltID(iElt).EQ.11) THEN
                iSeg=iSeg+1
                EltToSegMap(iElt)=iSeg
                SegToEltMap(iSeg)=iElt
              ELSE
                EltToSegMap(iElt)=0
              END IF

              ! Check to make sure no Element/Surface incompatability exists
	      IF (SrfType_FLG) THEN
                CALL EltSurfCompat(iElt,EltID_FLG,SrfType_FLG)
                IF (EltID_FLG.AND.SrfType_FLG) THEN
                  EltType_FLG=.TRUE.
                ELSE
                  WRITE(*,131)EltTypeName(EltID(iElt)),
     &            SrfTypeName(EltID(iElt))
                END IF
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'EltGrp',6)) THEN
              READ(VALUE,*) EltGrp(0,iElt)
              j=EltGrp(0,iElt)  ! # members in group
	      IF (j>0) THEN
                READ(VALUE,*)(EltGrp(i,iElt),i=0,j)
              ELSE IF (j<0) THEN
                ! element group range specified at this element
                EltGrp(0,iElt)=-EltGrp(0,iElt)
                !READ(VALUE,*) j,eltGrpStrt,eltGrpEnd
	        READ(VALUE,*) j,eltGrpStrt
                eltGrpEnd=eltGrpStrt+EltGrp(0,iElt)-1
                IF (eltGrpStrt<0 .OR. eltGrpEnd<0 .OR.
     &              (eltGrpStrt>eltGrpEnd)) THEN
                  WRITE(*,*)' **EltGrp: invalid specification in Rx!'
                  GO TO 99
                END IF
                k=1
                Do i=eltGrpStrt,eltGrpEnd
                  EltGrp(k,iElt)=i; k=k+1
                End Do
              end if
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'MrEltGrp',8)) THEN
              ! Allow an EltGrp to be specified with multiple ranges
              ! Syntax: #groups strt1 end1 strt2 end2 ...
              READ(VALUE,*) j ! number of element groups
              IF (j<=0 .or. j>10) THEN
                WRITE(*,*)
     &          ' **MrEltGrp: elt groups must be >=1 and <=10 in Rx!'
                GO TO 99
              END IF
              READ(VALUE,*) j,grpRange(1:2*j)
              EltGrp(0,iElt)=0; k=1
              Do g=1,j
                EltGrp(0,iElt) = EltGrp(0,iElt) +
     &            grpRange(2*g)-grpRange(2*g-1)+1
                Do i=grpRange(2*g-1),grpRange(2*g)
                  EltGrp(k,iElt)=i; k=k+1
                End Do
              End Do
              !print*,'iElt =',iElt,', EltGrp(0,iElt)=',EltGrp(0,iElt)
              !print*,'EltGrp ids =',EltGrp(1:EltGrp(0,iElt),iElt)
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'UpdateTElt',10)) THEN
              IF (LCMP(VALUE,'Y',1)) THEN
                LUpdateTElt_FLG(iElt)=.TRUE.
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'NSCount',7)) THEN
              READ(VALUE,*) NSCnt(iElt)
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'Surface',7)) THEN
              TypeString=VALUE

              ! Assign SrfType
              DO i=1,mSrfType  ! mSrfType = 11  -jzlou
                IF (LCMP(TypeString,SrfTypeName(i),4)) THEN
                  SrfType(iElt)=i
                  SrfType_FLG=.TRUE.
	          If (EltID(iElt)==3 .AND. SrfType(iElt)==9) Then
	            ! a reference (EltID=3) and GridData surface, a surface grid
		    ! data array will be defined to modulate wavefront. -jzlou
                    ifOPDModGrid(iElt)=.TRUE.
		  End If
	          If (.false. .AND. EltID(iElt)==1
     &                .AND. SrfType(iElt)==11) Then
	            print*,'** Found a reflector with AsphGrid srf'
	            STOP
	          End If
                 END IF
              END DO
              IF (.NOT.SrfType_FLG) THEN
                WRITE(*,132) TypeString
                GO TO 61
              END IF

              ! Set flags that aren't needed
              IF (SrfType(iElt).EQ.SrfType_Flat) THEN
                KcElt(iElt)=0d0
                KrElt(iElt)=-1d22
                fElt(iElt)=1d22
                eElt(iElt)=0d0
                KcElt_FLG=.TRUE.
                KrElt_FLG=.TRUE.
                fElt_FLG=.TRUE.
                eElt_FLG=.TRUE.
              ELSE IF (SrfType(iElt).EQ.SrfType_Anamorphic) THEN
                KcElt_FLG=.TRUE.
                KrElt_FLG=.TRUE.
                fElt_FLG=.TRUE.
                eElt_FLG=.TRUE.
              END IF
              IF ((SrfType(iElt).NE.SrfType_Aspheric).AND.
     &            (SrfType(iElt).NE.SrfType_AsGrData)) THEN
                AsphCoef_FLG=.TRUE.
              END IF
              IF (SrfType(iElt).NE.SrfType_Monomial) THEN
                MonCoef_FLG=.TRUE.
              END IF
              IF ((SrfType(iElt).NE.SrfType_Anamorphic).AND.
     &            (SrfType(iElt).NE.SrfType_Toric)) THEN
                AnaCoef_FLG=.TRUE.
              END IF
              IF (SrfType(iElt).NE.SrfType_UserDefined) THEN
                UDSrfCoef_FLG=.TRUE.
                UDSrfType_FLG=.TRUE.
                UDSrfParam_FLG=.TRUE.
              END IF
              IF (SrfType(iElt).NE.SrfType_Zernike) THEN
                ZernCoef_FLG=.TRUE.
                ZernTypeL_FLG=.TRUE.
              END IF
              IF (SrfType(iElt).NE.SrfType_GridData) THEN
                nGridMat_FLG=.TRUE.
                GridMat_FLG=.TRUE.
                GridSrfdx_FLG=.TRUE.
              END IF
              IF ((SrfType(iElt).LT.SrfType_Monomial).AND.
     &            (EltID(iElt).NE.SrfType_Toric).AND.
     &            (EltID(iElt).NE.SrfType_AsGrData)) THEN
                pMon_FLG=.TRUE.
                xMon_FLG=.TRUE.
                yMon_FLG=.TRUE.
                zMon_FLG=.TRUE.
              END IF
              IF ((SrfType(iElt).NE.SrfType_Monomial).AND.
     &            (SrfType(iElt).LT.SrfType_Zernike)) THEN
                lMon_FLG=.TRUE.
              END IF

              ! Check to make sure no Element/Surface incompatibility exists

              IF (EltID_FLG) THEN
                CALL EltSurfCompat(iElt,EltID_FLG,SrfType_FLG)
                IF (EltID_FLG.AND.SrfType_FLG) THEN
                  EltType_FLG=.TRUE.
                ELSE
                  WRITE(*,131)EltTypeName(EltID(iElt)),
     &            SrfTypeName(EltID(iElt))
                END IF
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'EltType',7)) THEN
              READ(VALUE,*) EltType(iElt)
              EltType_FLG=.TRUE.

              CALL EltConvert(EltType(iElt),EltID(iElt),SrfType(iElt))
C             WRITE(*,*)' Element='//EltTypeName(EltID(iElt))(1:12)
C             WRITE(*,*)' Surface='//SrfTypeName(SrfType(iElt))(1:12)
              EltID_FLG=.TRUE.
              SrfType_FLG=.TRUE.

              IF ((EltType(iElt).EQ.5).OR.(EltType(iElt).EQ.7).OR.
     &          (EltType(iElt).EQ.15).OR.(EltType(iElt).EQ.17).OR.
     &          (EltType(iElt).EQ.19)) THEN
                iSeg=iSeg+1
                EltToSegMap(iElt)=iSeg
                SegToEltMap(iSeg)=iElt
              ELSE
                EltToSegMap(iElt)=0
              END IF
              IF ((EltType(iElt).NE.7).AND.(EltType(iElt).NE.12).AND.
     &          (EltType(iElt).NE.20)) THEN
                AsphCoef_FLG = .TRUE.
              END IF
              IF (EltType(iElt).NE.25) THEN
                AnaCoef_FLG = .TRUE.
              ELSE
                fElt_FLG=.TRUE.
                eElt_FLG=.TRUE.
                KrElt_FLG = .TRUE.
                KcElt_FLG = .TRUE.
              END IF
              IF ((EltType(iElt).NE.14).AND.(EltType(iElt).NE.15)
     &            .AND.(EltType(iElt).NE.21))THEN
                ZernCoef_FLG = .TRUE.
                ZernTypeL_FLG = .TRUE.
                lMon_FLG = .TRUE.
              END IF
              IF ((EltType(iElt).NE.15).AND.(EltType(iElt).NE.16))
     &        THEN
                MonCoef_FLG = .TRUE.
                lMon_FLG = .TRUE.
              END IF
              IF ((EltType(iElt).NE.14).AND.(EltType(iElt).NE.15)
     &         .AND.(EltType(iElt).NE.16).AND.(EltType(iElt).NE.17)
     &         .AND.(EltType(iElt).NE.18).AND.(EltType(iElt).NE.19)
     &         .AND.(EltType(iElt).NE.21).AND.(EltType(iElt).NE.22)
     &         .AND.(EltType(iElt).NE.23).AND.(EltType(iElt).NE.25))
     &        THEN
                pMon_FLG = .TRUE.
                xMon_FLG = .TRUE.
                yMon_FLG = .TRUE.
                zMon_FLG = .TRUE.
              END IF
	      IF (EltType(iElt).NE.26) THEN
                UDSrfCoef_FLG = .TRUE.
                UDSrfType_FLG = .TRUE.
                UDSrfParam_FLG = .TRUE.
              END IF
              IF ((EltType(iElt).NE.22).AND.(EltType(iElt).NE.23)) THEN
                LensArrayType_FLG = .TRUE.
                LensArrayWidth_FLG = .TRUE.
              END IF
              IF ((EltType(iElt).NE.6).AND.(EltType(iElt).NE.7)) THEN
                h1HOE_FLG = .TRUE.
                h2HOE_FLG = .TRUE.
                OrderHOE_FLG = .TRUE.
                WaveHOE_FLG = .TRUE.
                RuleWidth_FLG = .TRUE.
                OrderHOE_FLG = .TRUE.
              ELSE IF (EltType(iElt).EQ.7) THEN
                h2HOE_FLG = .TRUE.
                WaveHOE_FLG = .TRUE.
              ELSE
                RuleWidth_FLG = .TRUE.
              END IF
              IF (EltType(iElt) .NE. 9) THEN  ! EltType 9 is NSRefractor
                ApType_FLG = .TRUE.
                ApVec_FLG = .TRUE.
                nObs_FLG = .TRUE.
                nObs(iElt) = 0
                xObs_FLG = .TRUE.
                DO I=1,mObs
                  ObsType_FLG(I) = .TRUE.
                  ObsVec_FLG(I) = .TRUE.
                END DO
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
            ! For System Optimization Only
            ! ------------------------------------------------------------------------------

            ELSE IF (LCMP(VAR_NAM,'VarElt',6)) THEN
              ! 2nd test below makes sure it's not a double-pass elt
              If (.not. isVarElt(iElt) .and. LnkElt(iElt)<0) Then
                isVarElt(iElt)=.true.
                nVarElt=nVarElt+1
	      End If
              !varElts(nVarElt)=iElt
              READ(VALUE,*)(varEltDOF(j,iElt),j=1,mVarDOF)
              nDOF_VarElt(nVarElt)=0
              DO i=1,mVarDOF
                IF (varEltDOF(i,iElt)/=0) THEN
                  nDOF_VarElt(nVarElt)=nDOF_VarElt(nVarElt)+1
                  DOF_VarElt(nDOF_VarElt(nVarElt),nVarElt)
     &              = DOF_NameList(i)
                END IF
              END DO
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'VarCons',7)) THEN
              LOptCons=.TRUE.
              j=0
              Do i=1,mVarDOF
                If (varEltDOF(i,iElt)/=0) Then
                  j=j+1
                End If
              End Do
              If (j>0) Then
                READ(VALUE,*) consArr(1:2*j)
                off=1
                Do i=1,mVarDOF
                  If (varEltDOF(i,iElt)/=0) Then
                    EltCons(1:2,i,iElt)=consArr(off:off+1)
                    off=off+2
                  End If
                End Do
              End If
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'OptBeamDir',10) .AND.
     &               (OptTarget==BEAM_TARGET)) THEN
	      If (.not. beamOpt%ifOptBeamDir) Then
		READ(VALUE,*) (beamOpt%nomBeamDir(j),j=1,3)
	        beamOpt%beamDirElt=iElt
	        beamOpt%ifOptBeamDir=.TRUE.
	      End If
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'OptBeamPos',10) .AND.
     &               (OptTarget==BEAM_TARGET)) THEN
              If (.not. beamOpt%ifOptBeamPos) Then
                READ(VALUE,*) (beamOpt%nomBeamPos(j),j=1,3)
                beamOpt%beamPosElt=iElt
                beamOpt%ifOptBeamPos=.TRUE.
	        If (OptTgtElt /= beamOpt%beamPosElt) Then
	          Write(*,*)'*** msmacosio warning: '
	          Write(*,*)'      OptTgtElt /= beamOpt%beamPosElt'
                 Write(*,*)'      OptTgtElt reset to beamOpt%beamPosElt'
                  OptTgtElt = beamOpt%beamPosElt
	        End If
              End If
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'OptBeamSize',11)) THEN
	      If (.not. beamOpt%ifOptBeamSize) Then
	        READ(VALUE,*) beamOpt%nomBeamSize
	        beamOpt%beamSizeElt=iElt
                beamOpt%ifOptBeamSize=.TRUE.
	      End If
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'OptBeamRefRayDir',15)) THEN
	      If (.not. beamOpt%ifOptBeamRefRayDir) Then
	        READ(VALUE,*) (beamOpt%nomBeamRefRayDir(j),j=1,3)
	        beamOpt%beamRefRayElt=iElt
	        beamOpt%ifOptBeamRefRayDir=.TRUE.
	      End If
	      GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'OptZern',7) .OR.
     &               LCMP(VAR_NAM,'VarZern',7)) THEN
	      If (.not. isVarElt(iElt) .and. LnkElt(iElt)<0) Then
                isVarElt(iElt)=.true.
                nVarElt=nVarElt+1
              End If
              READ(VALUE,*) nOptEltZern(iElt)
              READ(VALUE,*)
     &         (OptEltZernTerm(j,iElt),j=0,nOptEltZern(iElt))

!             print*,'iElt=', iElt,', nOptElt=',nOptEltZern(iElt),
!    &          ',  OptEltZernTerm=',
!    &          OptEltZernTerm(0:nOptEltZern(iElt),iElt)

              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'OptAsph',7)) THEN
              If (isVarElt(iElt)) Then
                READ(VALUE,*) nOptEltAsph(iElt)
                READ(VALUE,*)
     &          (OptEltAsphTerm(j,iElt),j=0,nOptEltAsph(iElt))
              End If
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'DpElt',5).OR.
     &               LCMP(VAR_NAM,'LINK',4)) THEN
              READ(VALUE,*) LnkElt(iElt)
              If (LnkElt(iElt).LE.0) Then
                WRITE(*,*)
     &          ' Link Elt value can not be negative at Elt',iElt
                WRITE(*,*)' Input file not properly loaded'
                GO TO 99
              End If
              GO TO 61

            ! ------------------------------------------------------------------------------ fElt
	    ELSE IF (LCMP(VAR_NAM,'fElt',4)) THEN
              READ(VALUE,*) fElt(iElt)
              fElt_FLG    = .TRUE.
              KrElt_FLG   = .TRUE.
              KcElt_FLG   = .TRUE.
              KcElt(iElt) = 0D0;
              KrElt(iElt) = 0D0
              GO TO 61

            ! ------------------------------------------------------------------------------ eElt
            ELSE IF (LCMP(VAR_NAM,'eElt',4)) THEN
              READ(VALUE,*)eElt(iElt)
              eElt_FLG=.TRUE.
              KrElt_FLG = .TRUE.; KcElt_FLG = .TRUE.
              KcElt(iElt) = 0.0D0; KrElt(iElt) = 0.0D0
              GO TO 61

            ! ------------------------------------------------------------------------------
            !  ...
            !  nAsphCoef=  7
            !  AsphCoefs=  3.825055423250000d-009 -1.713067033850000d-014  8.407870699740000d-020 -4.925133142010000d-025
            !             -1.623672889230000d-028  4.446562373320000d-033 -7.802857348350000d-038
            !  ...
            !
            ! Group of 4 Coefficients -- if not formatted this way, it will crash
            ! ------------------------------------------------------------------------------ nAsphCoefs
            ELSE IF (LCMP(VAR_NAM,'nAsphCoefs',7)) THEN

              READ(VALUE,*) nAsphCoef
              IF (nAsphCoef<1 .OR. nAsphCoef>mAsphCoef) THEN
                WRITE(*,*)'**WARNING: nAsphCoefs is out of range -- failure'
                GO TO 99
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------ AsphCoefs
            ELSE IF (LCMP(VAR_NAM,'AsphCoef',8)) THEN

              BLOCK
                INTEGER, PARAMETER :: Grp = 4  ! defines # of Coefs. within a Grp per line
                INTEGER            :: iZ, nZ
                ! --------------------------

                AsphCoef(:,iElt) = 0d0
                IF (nAsphCoef<=Grp) THEN
                  READ(VALUE,*)(AsphCoef(j,iElt),j=1,nAsphCoef)
                ELSE

                  READ(VALUE,*)(AsphCoef(j,iElt),j=1,Grp)

                  ! complete line group ("Grp" Coefficients)
                  DO iZ=1, nAsphCoef/Grp-1
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(AsphCoef(j,iElt),j=iZ*Grp+1,(iZ+1)*Grp)
                  END DO

                  ! incomplete Grp
                  nZ = MOD(nAsphCoef,Grp)
                  iZ = nAsphCoef/Grp
                  IF (nZ>0) THEN
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(AsphCoef(j,iElt),j=iZ*Grp+1,iZ*Grp+nZ)
                  END IF
                END IF

                nAsphCoef = nAsphCoef_Default  ! reset to default
              END BLOCK

              AsphCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------ AnaCoef
            ELSE IF (LCMP(VAR_NAM,'AnaCoef',7)) THEN
              READ(VALUE,*)(AnaCoef(j,iElt),j=1,4)
              AnaCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------ ZernType
            ELSE IF (LCMP(VAR_NAM,'ZernType',8)) THEN
              TypeString=VALUE
	      If (LCMP(TypeString,'Noll',4))
     &  	TypeString='NormNoll'  ! backward compatibility

              ! Assign ZernTypeL
              If (LCMP(TypeString,'Norm',4)) Then
                ! normalized Zernike coeff
                DO i=4,9
                  IF (LCMP(TypeString,ZernTypeNameL(i),7)) THEN
                    ZernTypeL(iElt)=i
                    ZernTypeL_FLG=.TRUE.
                  END IF
                END DO
	        If (ZernTypeL(iElt)==9) Then
                  ! Annual
                  i=1
                  Do ! find first ' ' in VALUE
                    If (VALUE(i:i) == ' ') exit
                    i=i+1
                  End Do
                  Do ! skip all leading ' '
		    If (VALUE(i:i) /= ' '.OR.(i==MacosValLen)) exit
                    i=i+1
                  End Do
	          If (i<MacosValLen)
     &              READ(VALUE(i:),*) ZernAnnuRatio(iElt)
                End If
              Else
                ! un-normalized Zernike coeff
                DO i=1,3
                  IF (LCMP(TypeString,ZernTypeNameL(i),3)) THEN
                    ZernTypeL(iElt)=i
                    ZernTypeL_FLG=.TRUE.
                  END IF
                END DO
!#if 0
!	        IF (LCMP(TypeString,ZernTypeNameL(7),7)) THEN
!                  ! Hex zernike
!                  ZernTypeL(iElt)=7
!                  ZernTypeL_FLG=.TRUE.
!                END IF
!#endif
              End If
              IF (.NOT.ZernTypeL_FLG) THEN
                WRITE(*,133) TypeString(1:ICLEN(TypeString))
		WRITE(*,*)'Check and reload Rx!'
		LOAD_SUCCESS=.FALSE.; GO TO 70
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------ nZernCoef
            ! if not defined, default value is used
            ELSE IF (LCMP(VAR_NAM,'nZernCoef',9)) THEN
              READ(VALUE,*) nZernCoef
              if ((nZernCoef<1).or.(nZernCoef>mZernModes)) GOTO 99
              GO TO 61

            ! ------------------------------------------------------------------------------ ZernCoef
            ELSE IF (LCMP(VAR_NAM,'ZernCoef',8)) THEN

              ! QUICK FIX to make it backwards compatible
              ! ==> Flag "nZernCoef" can be used to define max modes to read in
              BLOCK
                INTEGER, PARAMETER :: Grp = 6  ! defines # of Coefs. within a Grp per line
                INTEGER            :: iZ, nZ
                ! --------------------------

                ZernCoef(:,iElt) = 0d0
                IF (nZernCoef<=Grp) THEN
                  READ(VALUE,*)(ZernCoef(j,iElt),j=1,nZernCoef)
                ELSE

                  READ(VALUE,*)(ZernCoef(j,iElt),j=1,Grp)

                  ! complete line group ("Grp" Coefficients)
                  DO iZ=1, nZernCoef/Grp-1
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(ZernCoef(j,iElt),j=iZ*Grp+1,(iZ+1)*Grp)
                  END DO

                  ! incomplete Grp
                  nZ = MOD(nZernCoef,Grp)
                  iZ = nZernCoef/Grp
                  IF (nZ>0) THEN
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(ZernCoef(j,iElt),j=iZ*Grp+1,iZ*Grp+nZ)
                  END IF
                END IF
                nZernCoef = nZernModes_Default  ! reset to default
              END BLOCK

	      ! Added by jzlou, 03/2005
	      If ( ZernTypeL(iElt)>3 .and. ZernTypeL(iElt)<7)
     &        Then ! normalized Zern type
                ! Set up normalizing coef for circular Zernike terms
                Call SetZernNormCoef(ZernTypeL(iElt))
	      End If

	      ZernCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'ZernAnnularRatio',12)) THEN
              READ(VALUE,*) ZernAnnuRatio(iElt)
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'ZernCenter',8)) THEN
              READ(VALUE,*) ZernCent(1:3,iElt)
              ZernCent_FLG=.TRUE.
              zernUsrOpt(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'ZernXDir',5)) THEN
              READ(VALUE,*) ZernXGrid(1:3,iElt)
              zernUsrOpt(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'ZernYDir',5)) THEN
              READ(VALUE,*) ZernYGrid(1:3,iElt)
              zernUsrOpt(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'ZernRad',7)) THEN
              READ(VALUE,*) ZernRad(iElt)
              zernUsrOpt(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'MonCoef',7)) THEN
              READ(VALUE,*)(MonCoef(j,iElt),j=1,6)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=7,12)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=13,18)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=19,24)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=25,30)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=31,36)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=37,42)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=43,48)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=49,54)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=55,60)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=61,66)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=67,72)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=73,78)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=79,84)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=85,90)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=91,96)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=97,102)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=103,108)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=109,114)
              READ(RX_FID, 123, END=99)STRING
              READ(STRING,*)(MonCoef(j,iElt),j=115,120)
              MonCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'UDSrfFile',6)) THEN
              UDSrfFile(iElt)=VALUE
              IF (.NOT.UDSrfParam_FLG) THEN
                WRITE(*,*)
     &      ' WARNING: Must specify UDSrfParam before loading UDSrfFile'
                EXIST=.FALSE.
              ELSE
                CALL UDSinit(EXIST,iElt)
              END IF
              IF (.NOT.EXIST) THEN
                WRITE(*,*)' UDSrfFile not loaded correctly'
              END IF
              UDSrfCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'UDSrfType',6)) THEN
              READ(VALUE,*) UDSrfType(iElt)
              UDSrfType_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'UDSrfParam',9)) THEN
              READ(VALUE,*) (UDSrfParam(j,iElt),j=1,6)
              UDSrfParam_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'nGridMat',3)) THEN
              IF (iEltToGridSrf(iElt).EQ.0) THEN
                jGridSrf=jGridSrf+1
                iEltToGridSrf(iElt)=jGridSrf
              END IF
              READ(VALUE,*) nGridMat(iElt)
              nGridMat_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ! The two fields below should be combined into one -jzlou
	    ELSE IF (LCMP(VAR_NAM,'nMetPos',7)) THEN
	      IF (iEltToMetSrf(iElt).EQ.0) THEN
	        jMetSrf=jMetSrf+1
                iEltToMetSrf(iElt)=jMetSrf
	      END IF
	      READ(VALUE,*) nMetPos(iElt)
	      Do i=1,nMetPos(iElt)
                READ(RX_FID, 123, END=99) STRING
                READ(STRING,*) SrfMetPos(1:3,i,iEltToMetSrf(iElt))
              End Do
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'tMetElt',7)) THEN
	      If (nMetPos(iElt).LE.0) Then
                Write(*,*)' ** MACOS Rx error: iElt,nMetPos =',
     &                    iElt,nMetPos
                Write(*,*)
     &          ' **   nMetPos field must appear before tMetElt field,'
                Write(*,*)
     &          ' **   and nMetPos must be positive integer'
                STOP
              End If
              READ(VALUE,*) tMetSrf(iElt),ntMetPos(iElt)
!	      print*,' **tMetElt: iElt,iEltToMetSrf(iElt) =',
!    &               iElt,iEltToMetSrf(iElt)
!             print*,' **nMetPos(iElt),tMetSrf(iElt),ntMetPos(iElt) =',
!    &               nMetPos(iElt),tMetSrf(iElt),ntMetPos(iElt)
              If (tMetSrf(iElt)>0) Then
	        Do i=1,nMetPos(iElt)
	          READ(RX_FID, 123, END=99) STRING
                  READ(STRING,*)
     &              metBeamFlg(1:ntMetPos(iElt),i,iEltToMetSrf(iElt))
	        End Do
!               READ(VALUE,*)tMetSrf(iElt),ntMetPos(iElt),
!    &                       tMetSrfPos(1:ntMetPos(iElt),iElt)
              End If
              print*,' ** End tMetSrf(iElt) =',tMetSrf(iElt)
              GO TO 61

	    ELSE IF (LCMP(VAR_NAM,'nAmplMat',3)) THEN
              IF (.false. .and. iEltToAmplSrf(iElt).EQ.0) THEN
                jAmplSrf=jAmplSrf+1
                iEltToAmplSrf(iElt)=jAmplSrf
              END IF
              READ(VALUE,*) nAmplMat(iElt)
              nAmplMat_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'GridFile',8)) THEN
              IF (iEltToGridSrf(iElt).EQ.0) THEN
                jGridSrf=jGridSrf+1
                iEltToGridSrf(iElt)=jGridSrf
              END IF
              GridFile(iElt)=VALUE

	      k=ICLEN(GridFile(iElt))

	      If (.not.LCMP(GridFile(iElt)(k-2:k),'bin',3)) Then
	        ! ascii
                CALL GridInit(GridFile(iElt),jGridSrf,nGridMat(iElt),
     &                        ifGridInit,.false.)
	      Else
	        ! binary
	        CALL GridInit(GridFile(iElt),jGridSrf,nGridMat(iElt),
     &                        ifGridInit,.true.)
	      End If

	      ! added 09/20/2016 -jzlou
	      ifGridDataDefined(iEltToGridSrf(iElt)) = ifGridInit

              IF (.NOT.EXIST) THEN
                WRITE(*,*)' GridFile not loaded correctly, iELT=',iElt
              END IF
              GridMat_FLG=.TRUE.

	      ! if a reference surface (EltID=3), the grid is defined
              ! for wavefront/OPD modulation -jzlou
	      !If (EltID(iElt)==3) ifOPDModGrid(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'GridSrfOrder',10)) THEN
              READ(VALUE,*) GridSrfOrder(iElt)
	      If (GridSrfOrder(iElt) /= 1 .and.
     &            GridSrfOrder(iElt) /= 3) Then
	        WRITE(*,*)' ** Invalid GridSrfOrder value for Elt ',iElt
                WRITE(*,*)' ** Reset to bi-linear surface interpolator'
	      End If
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'AmplFile',8)) THEN
	      IF (iEltToAmplSrf(iElt).EQ.0) THEN
	        jAmplSrf=jAmplSrf+1
	        iEltToAmplSrf(iElt)=jAmplSrf
	        ifAmplSrf(iElt)=.TRUE.
	        LTaperedElt(iElt)=.true.
	      END IF

              AmplFile(iElt)=VALUE
	      Call AmplInit(AmplFile(iElt),jAmplSrf,nAmplMat(iElt),
     &                      ifAmplInit)
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'GridSrfdx',9)) THEN
              IF (iEltToGridSrf(iElt).EQ.0) THEN
                jGridSrf=jGridSrf+1
                iEltToGridSrf(iElt)=jGridSrf
              END IF
              READ(VALUE,*)GridSrfdx(iElt)
              GridSrfdx_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'AmplSrfdx',9)) THEN
              IF (.false. .and. iEltToAmplSrf(iElt).EQ.0) THEN
                jAmplSrf=jAmplSrf+1
                iEltToAmplSrf(iElt)=jAmplSrf
              END IF
              READ(VALUE,*)AmplSrfdx(iElt)
              AmplSrfdx_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'KrElt',5)) THEN
              READ(VALUE,*)KrElt(iElt)
              KrElt_FLG=.TRUE.
              eElt_FLG=.TRUE.
              fElt_FLG=.TRUE.
              eElt(iElt) = 0.0D0
              fElt(iElt) = 0.0D0
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'KcElt',5)) THEN
              READ(VALUE,*)KcElt(iElt)
              KcElt_FLG=.TRUE.
              eElt_FLG=.TRUE.
              fElt_FLG=.TRUE.
              eElt(iElt) = 0.0D0
              fElt(iElt) = 0.0D0
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'pMon',4)) THEN
              READ(VALUE,*)(pMon(j,iElt),j=1,3)
              pMon_FLG=.TRUE.
	      if (.not.ZernCent_FLG)
     &          ZernCent(1:3,iElt)=pMon(1:3,iElt)
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'xMon',4)) THEN
              READ(VALUE,*)(xMon(j,iElt),j=1,3)
              xMon_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'yMon',4)) THEN
              READ(VALUE,*)(yMon(j,iElt),j=1,3)
              yMon_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'zMon',4)) THEN
              READ(VALUE,*)(zMon(j,iElt),j=1,3)
              zMon_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'lMon',4)) THEN
              READ(VALUE,*)lMon(iElt)
              lMon_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'LensArrayType',12)) THEN
              READ(VALUE,*)LensArrayType(iElt)
              LensArrayType_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'LensArrayWidth',12)) THEN
              READ(VALUE,*)LensArrayWidth(iElt)
              LensArrayWidth_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'LensArrayIndRef',12)) THEN
	      ! Warning:
	      !   Here a certain ordering of the lens array is
	      !   assumed, for both hex and rec arrays.
	      IsVarIndRefElt(iElt)=.true.
! --------------------------------
! ToDo: Fix this calling -- not agreeing with Interface Definition
!
!              Call InitLensletIndexAndCtr(LensArrayType(iElt),
!     &		     LensArrayWidth(iElt),
!     &               TotalLensletsArray(iElt),iElt)

              Call InitLensletIndexAndCtr(LensArrayType(iElt),
     &		     LensArrayWidth(iElt),iElt)
! ---------------------------------
	      IF (abs(LensArrayType(iElt))==1) THEN
                TotalLensletsArray(iElt)=nHexLenlets
              ELSE
                TotalLensletsArray(iElt)=nRecLenlets
              END IF

	      READ(VALUE,*)TotalLensletsArray(iElt),nValsPerRow
	      If (nValsPerRow>0) Then
	        k=TotalLensletsArray(iElt); i=1
	        Do
	         if (k<nValsPerRow) exit
	          ! read in 'nValsPerRow' IndRef values from Rx
	          READ(RX_FID, 123, END=99)STRING
	          READ(STRING,*)
     &	           LensletStrArray(i:i+nValsPerRow-1,iElt)%IndRef
	          k=k-nValsPerRow; i=i+nValsPerRow
	        End Do
	        If (k>0) Then
	          READ(RX_FID, 123, END=99)STRING
	          READ(STRING,*)
     &             LensletStrArray(i:i+k-1,iElt)%IndRef
	        End If
	      End If
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'XYIndRefFile',12)) THEN
	      IF (EltID(iElt)==10) THEN ! only valid on LensArray
	        IsVarIndRefElt(iElt)=.true.
	        LensletStrArray(1:mLensLet,iElt)%IndRef=IndRef(iElt) ! default
	        Call InitLensletIndexAndCtr(LensArrayType(iElt),
     &	    	                            LensArrayWidth(iElt),iElt)
	        IF (abs(LensArrayType(iElt))==1) THEN
                  TotalLensletsArray(iElt)=nHexLenlets ! total hex elements
                ELSE
                  TotalLensletsArray(iElt)=nRecLenlets ! total rec elements
                END IF
	        Call InitLensArrIndRefFrmXYFile(VALUE,iElt,
     &                                          LensArrayType(iElt))
	      END IF
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'IndRefGridFile',12)) THEN
	      If (EltID(iElt)==6) Then ! only valid on Refractor
	        IsVarIndRefElt(iElt)=.true.
	        ! IndRefGridArr(iElt)%width is read in from index grid file -jzlou
	        Call InitIndRefGridFrmFile(VALUE,iElt)
	        !print*,'Stop after reading index grid'
	        !stop
	      End If
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'h1HOE',5)) THEN
              READ(VALUE,*)(h1HOE(j,iElt),j=1,3)
              h1HOE_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'h2HOE',5)) THEN
              READ(VALUE,*)(h2HOE(j,iElt),j=1,3)
              h2HOE_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'OrderHOE',8)) THEN
              READ(VALUE,*)OrderHOE(iElt)
              OrderHOE_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'WaveHOE',7)) THEN
              READ(VALUE,*)WaveHOE(iElt)
              WaveHOE_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'RuleWidth',9)) THEN
              READ(VALUE,*)RuleWidth(iElt)
              RuleWidth_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'DoeWL',5)) THEN
              READ(VALUE,*) DoeWL(nDOE)
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'DoePhase',8)) THEN
	      READ(VALUE,*)(DoePhProf(j,nDOE),j=1,6)
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'DoeFl',5)) THEN
              READ(VALUE,*) DoeFl(nDOE)
              !DoeFl_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'psiElt',6)) THEN
              READ(VALUE,*)(psiElt(j,iElt),j=1,3)
              CALL DUNITIZE(psiElt(1,iElt))
              psiElt_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'VptElt',6)) THEN
              READ(VALUE,*)(VptElt(j,iElt),j=1,3)
              VptElt_FLG=.TRUE.
	      If ((.not.pMon_FLG).and.(.not.ZernCent_FLG))
     &          ZernCent(1:3,iElt)=VptElt(1:3,iElt)
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'RptElt',6)) THEN
              READ(VALUE,*)(RptElt(j,iElt),j=1,3)
              RptElt_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'IndRef',6)) THEN
              READ(VALUE,*)IndRef(iElt)
              IndRefE_FLG=.TRUE.
              GlassElt_FLG=.TRUE.
              GlassCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'ArrIndRef',9)) THEN
	      ArrIndRef(1:nOptWavelen,iElt)=1d0 ! default first
              Call Get_Values(VALUE,ArrIndRef(:,iElt),nTok)
                     ! nTok returns number of values read into ArrIndRef
	      If (nTok /= nOptWavelen) Then
                WRITE(*,*)' '
                WRITE(*,*)
     &           '** Load Rx warning at element',iElt,':'
                WRITE(*,*)
     &   '   Inconsistent number of values between '//
     &   '   wavelengths and '//
     &   '   refractive indexes'
	        Pause
	      End If
	      IndRef(iElt)=ArrIndRef(1,iElt)
              IndRefE_FLG=.TRUE.
              GlassElt_FLG=.TRUE.
              GlassCoef_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'GradInd',7)) THEN
              READ(VALUE,*)GrinRefArr(iElt)%modelType
	      IndRefE_FLG=.TRUE.
	      GlassElt_FLG=.TRUE.
              GlassCoef_FLG=.TRUE.
	      IsGrinRefElt(iElt)=.TRUE.
	      numGradInd=numGradInd+1
	      iElt2iGradInd(iElt)=numGradInd
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'GradCoef',8)) THEN
	      If (GrinRefArr(iElt)%modelType <= 2 .AND.
     &            GrinRefArr(iElt)%modelType >= 5) Then
	        READ(VALUE,*)GrinRefArr(iElt)%modelCoef(1:8)
	      Else If (GrinRefArr(iElt)%modelType==1 .OR.
     &                 GrinRefArr(iElt)%modelType==6) Then
                READ(VALUE,*)GrinRefArr(iElt)%modelCoef(1:3)
	      Else If (GrinRefArr(iElt)%modelType==7) Then
	        READ(VALUE,*)GrinRefArr(iElt)%modelCoef(1:5)
	      Else
	        ! If unknown graded-index modelType or undefine,
	        ! return from function with error message
	       Write(*,*)'** msmacosio.inc: '//
     &         'invalid or undefined graded index'
	       Write(*,*)
     &         '** model type. Define a valid model type in Rx first!'
	       return
	      End If
	      GrinRefArr(iElt)%dt=GrinRefArr(iElt)%modelCoef(1)
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'GradLensZ',9)) THEN
	      READ(VALUE,*)GrinRefArr(iElt)%lensDepth
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'Coating',6)) THEN
	      READ(VALUE,*)EltCoat(iElt)
	      If (EltCoat(iElt)>0) Then
	        ! Coating layers are specified in Rx, from outmost
	        ! layer to innermost layer, each layer having three values:
	        ! real_part_IndRef,imag_part_IndRef,layer_thickness,
	        ! and layer_thickiness is in terms of wavelength lambda
	        i=EltCoat(iElt)
	        Do k=1,EltCoat(iElt)
	          READ(RX_FID, 123, END=99) STRING
                  READ(STRING,*)IndRefArr(k,iElt),ExtincArr(k,iElt),
     &				EltCoatThk(k,iElt)
	          EltCoatThk(k,iElt)
     &              =EltCoatThk(k,iElt)*Wavelen/IndRefArr(k,iElt)
	        End Do
	        IndRefArr(0,iElt)=IndRef(iElt-1) ! IndRef before coating layers
	        ExtincArr(0,iElt)=Extinc(iElt-1)
	        IndRefArr(i+1,iElt)=IndRef(iElt) ! IndRef after coating layers
	        ExtincArr(i+1,iElt)=Extinc(iElt)
	      End If

	      IF (.FALSE.) THEN
	       print*,'iElt=', iElt,'  Coating layers=',EltCoat(iElt)
	       DO k=0,EltCoat(iElt)+1
	         print*,'IndRefArr(k),ExtincArr(k),EltCoatThk(k)=',
     &	 	   IndRefArr(k,iElt),ExtincArr(k,iElt),
     &		   EltCoatThk(k,iElt)
	       END DO
	      END IF

              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'GlassElt',6)) THEN
              GlassElt(iElt) = VALUE
              IndRefE_FLG    = .TRUE.
              GlassElt_FLG   = .TRUE.
              ! Assume Sellmeier glass by default
              GlassType(iElt)=SellmeierGlass
              LGlass(iElt)=.FALSE.
              i=1
              DO WHILE (((.NOT.LGlass(iElt))).AND.(i.LE.mGlass))
                IF (LCMP(GlassElt(iElt),GlassName(i),
     &                   ICLEN(GlassElt(iElt)))) THEN
                  LGlass(iElt)=.TRUE.
                  CALL DEQUATE(GlassCoef(1,iElt),GlassTable(1,i),6)
                END IF
                i=i+1
              END DO
              IF (LGlass(iElt)) THEN
                GlassCoef_FLG=.TRUE.
              ELSE
                GlassCoef_FLG=.FALSE.
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'GlassType',7)) THEN
              IF (LCMP(VALUE,'Scho',4)) THEN
                GlassType(iElt)=SchottGlass
                LGlass(iElt)=.TRUE.
                GlassCoef_FLG=.TRUE.
	      ELSE IF (LCMP(VALUE,'Sell',4)) THEN
                GlassType(iElt)=SellmeierGlass
                LGlass(iElt)=.TRUE.
                GlassCoef_FLG=.TRUE.
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'GlassCoef',6)) THEN
              READ(VALUE,*)(GlassCoef(j,iElt),j=1,6)
              GlassCoef_FLG=.TRUE.
              LGlass(iElt)=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'Extinc',6)) THEN
              READ(VALUE,*)Extinc(iElt)
              ExtincE_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'Regrid',6)) THEN ! added by jzlou
              If (LCMP(VALUE,'Y',1)) Then
	        ifRegrid(iElt)=.TRUE.
	      End If
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'zElt',4)) THEN
              READ(VALUE,*)zElt(iElt)
              zElt_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'PropType',8)) THEN
              IF (LCMP(VALUE,PropTypeName(10),5).OR.
     &        LCMP(VALUE,'10',2)) THEN
                PropType(iElt)=10
              ELSE IF (LCMP(VALUE,PropTypeName(11),7).OR.
     &        LCMP(VALUE,'11',2)) THEN
                PropType(iElt)=11
              ELSE IF (LCMP(VALUE,PropTypeName(12),7).OR.
     &        LCMP(VALUE,'12',2)) THEN
                PropType(iElt)=12
              ELSE IF (LCMP(VALUE,PropTypeName(13),7).OR.
     &        LCMP(VALUE,'13',2)) THEN
                PropType(iElt)=13
              ELSE IF (LCMP(VALUE,PropTypeName(14),7).OR.
     &        LCMP(VALUE,'14',2)) THEN
                PropType(iElt)=14
              ELSE IF (LCMP(VALUE,PropTypeName(15),7).OR.
     &        LCMP(VALUE,'15',2)) THEN
                PropType(iElt)=15
              ELSE IF (LCMP(VALUE,PropTypeName(1),7).OR.
     &        LCMP(VALUE,'1',1)) THEN
                PropType(iElt)=1
              ELSE IF (LCMP(VALUE,PropTypeName(2),7).OR.
     &        LCMP(VALUE,'2',1)) THEN
                PropType(iElt)=2
              ELSE IF (LCMP(VALUE,PropTypeName(3),7).OR.
     &        LCMP(VALUE,'3',1)) THEN
                PropType(iElt)=3
              ELSE IF (LCMP(VALUE,PropTypeName(4),7).OR.
     &        LCMP(VALUE,'4',1)) THEN
                PropType(iElt)=4
              ELSE IF (LCMP(VALUE,PropTypeName(5),7).OR.
     &        LCMP(VALUE,'5',1)) THEN
                PropType(iElt)=5
              ELSE IF (LCMP(VALUE,PropTypeName(6),7).OR.
     &        LCMP(VALUE,'6',1)) THEN
                PropType(iElt)=6
              ELSE IF (LCMP(VALUE,PropTypeName(7),7).OR.
     &        LCMP(VALUE,'7',1)) THEN
                PropType(iElt)=7
              ELSE IF (LCMP(VALUE,PropTypeName(8),7).OR.
     &        LCMP(VALUE,'8',1)) THEN
                PropType(iElt)=8
              ELSE IF (LCMP(VALUE,PropTypeName(9),7).OR.
     &        LCMP(VALUE,'9',1)) THEN
                PropType(iElt)=9
              END IF
              PropType_FLG=.TRUE.
	      GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'PinHole',6)) THEN
              READ(VALUE,*)PinHole(iElt)
              PinHole_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'nCoatElt',5)) THEN
              READ(VALUE,*)nCoatElt(iElt)
              nCoat_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'CoatIndxElt',8)) THEN
              READ(VALUE,*)(CoatIndxElt(j,iElt),j=1,nCoatElt(iElt))
              CoatIndx_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'CoatThkElt',7)) THEN
              READ(VALUE,*)(CoatThkElt(j,iElt),j=1,nCoatElt(iElt))
              CoatThk_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'nECoord',7)) THEN
              READ(VALUE,*)nECoord(iElt)
              nECoord_FLG=.TRUE.
              IF ((nECoord(iElt) .EQ. 0) .OR.
     &          (nECoord(iElt) .EQ. -6)) THEN
                TElt_FLG=.TRUE.
                CALL DZERO(TElt(1,1,iElt),36)
                DO 65 j=1,6
                  TElt(j,j,iElt)=1.0D0
   65           CONTINUE
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'ApStop',6)) THEN
              RxStopSet=.TRUE.
              IF (EltID(iElt)==0 .OR.
     &          EltID(iElt)==NSRefractorElt .OR.
     &          EltID(iElt)==NSReflectorElt .OR.
     &          EltID(iElt)==SegmentElt) THEN
                  ifStopSet=.FALSE.
                  EltStopSet=.FALSE.
              ELSE
                StopElt=iElt
                READ(VALUE,*) StopOffset(1:2)
                EltStopSet=.TRUE.
	        nTok = GetNumToken(VALUE)
                LStopAtEltAutoSet = (nTok > 2)
                !IF (nTok > 2) THEN
                !  !print*,' ** STOP at element will be auto-reset!'
	        !  LStopAtEltAutoSet=.TRUE.
                !END IF
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'ApType',6)) THEN

              LApTypeDefined(iElt)=.TRUE.

              IF (LCMP(VALUE,ApTypeName(0),2).OR.
     &            LCMP(VALUE,'0',1)) THEN
                ApType(iElt)=0
                LApTypeDefined(iElt)=.FALSE.
              ELSE IF (LCMP(VALUE,ApTypeName(1),2).OR.
     &                 LCMP(VALUE,'1',1)) THEN
                ApType(iElt)=1
              ELSE IF (LCMP(VALUE,ApTypeName(2),2).OR.
     &                 LCMP(VALUE,'2',1)) THEN
                ApType(iElt)=2
              ELSE IF (LCMP(VALUE,ApTypeName(3),10).OR.
     &                 LCMP(VALUE,'3',1)) THEN
                ApType(iElt)=3 ! ApType = 'Tapered_Circular'
              ELSE IF (LCMP(VALUE,ApTypeName(4),5).OR.
     &                 LCMP(VALUE,'4',1)) THEN
                ApType(iElt)=4 ! ApType = 'Elliptical'
              ELSE IF (LCMP(VALUE,ApTypeName(5),10).OR.
     &                 LCMP(VALUE,'5',1)) THEN
                ApType(iElt)=5 ! ApType = 'Tapered_Elliptical'
              ELSE IF (LCMP(VALUE,ApTypeName(6),3).OR.
     &                 LCMP(VALUE,'6',1)) THEN
                ApType(iElt)=6 ! ApType = 'Hexagonal'
              ELSE IF (LCMP(VALUE,ApTypeName(7),7).OR.
     &                 LCMP(VALUE,'7',1)) THEN
                ApType(iElt)=7 ! ApType = 'Polygonal'

	      ELSE IF (LCMP(VALUE,ApTypeName(8),12).OR.
     &                 LCMP(VALUE,'8',1)) THEN
                ApType(iElt)=8 ! ApType = 'Tapered_Polygonal'

                !print*,'Tapered_Polygonal: VALUE =',VALUE
                !READ(VALUE,*) chBuf,apodType
                 ! get apodization type, can be one of Linear, Gaussian, Step
                 ! In Rx, one can have, for example
                 !   ApType=  Tapered_Polygonal Gaussian
                 !   ApVec=   6.6175000000d+02  -1.1461846219d+03 6 1.2d02
                 !            2.0160000000d+00   7.5852762191d+02
                 !            ...
                 !  where the last two values of ApVec (first line) are apodization region
                 !  width (from polygon edge) and the standard deviation of the Gaussian
                 !  distribution respectively.
#if 1
	         ! Get apodization type
                 READ(VALUE,*) chBuf,apodType
                 If (LCMP(apodType,'Step',4)) Then
                   iApodType=stepApod
                 Else If (LCMP(apodType,'Lin',3)) Then
                   iApodType=linApod
                 Else If (LCMP(apodType,'Gaussian',4)) Then
                   iApodType=gaussApod
                 Else
                   WRITE(*,*) '** Unknown apodization type at Elt ',iElt
                   GO TO 99
                 End If
#endif
            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VALUE,ApTypeName(9),8).OR.
     &                 LCMP(VALUE,'9',1)) THEN
	        If (.not.LApSunBlades) Then
                  ApType(iElt)=9 ! ApType = 'SunBlades'
		  LApSunBlades=.TRUE.
	        Else
                  WRITE(*,*)
     &          ' **Load error: ApType SunBlades multiply defined!'
	          GO TO 61
	        End If
              END IF
              ApType_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'ApVec',5)) THEN
              LApVecDefined(iElt)=.TRUE.

              IF (ApType_FLG) THEN
                IF (ApType(iElt).EQ.1) THEN
                  READ(VALUE,*)(ApVec(j,iElt),j=1,3)
                ELSE IF (ApType(iElt).EQ.2) THEN
                  READ(VALUE,*)(ApVec(j,iElt),j=1,4)
                ELSE IF (ApType(iElt).EQ.3) THEN
                  ! Tapered_Circular aperture,
                  ! last value read in is size of tapered region
                  ! as a percentage of aperture radius -jzlou
                  READ(VALUE,*)(ApVec(j,iElt),j=1,4)
                  LTaperedElt(iElt)=.true.
                ELSE IF (ApType(iElt).EQ.4) THEN
                  ! Elliptical aperture, read in a,b,c_x,c_y
                  READ(VALUE,*)(ApVec(j,iElt),j=1,4)
                ELSE IF (ApType(iElt).EQ.5) THEN
                  ! Tapered elliptical aperture, read in a,b,c_x,c_y,
                  ! percent_of_taper
                  READ(VALUE,*)(ApVec(j,iElt),j=1,5)
                  LTaperedElt(iElt)=.true.
                ELSE IF (ApType(iElt).EQ.6) THEN
                  ! Hexagonal aperture, read in h,xc,yc,orient,
                  ! orient specified in degree for orientation of hex
                  READ(VALUE,*)(ApVec(j,iElt),j=1,4)
                  Call SetHexBound(ApVec(1,iElt),iElt)
                ELSE IF (ApType(iElt).EQ.7) THEN
                  READ(VALUE,*) (ApVec(j,iElt),j=1,3)
                    ! # of vertices of polygon
                  nVertex=ApVec(3,iElt)
                  If (nVertex < 3) Then
                    WRITE(*,*)
     &               '** Error in MACOS Rx: a polygon aperture'//
     &               ' must have at least three vertices!'
                    GOTO 99
                  End If
                  Do k=1,nVertex
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(PolyApVtx(j,k,iElt),j=1,2)
                  End Do
                  Call SetCvxPolyApBound(iElt)

	        ELSE IF (ApType(iElt).EQ.8) THEN
                  ! Tapered (apodized) polygonal aperture
	          !
	          ! Now support three types of apodization types: step, linear and Gaussian
                  If (iApodType==stepApod) Then
                    ! read in xc,yc,nVertex,width,step_percent
                    READ(VALUE,*) (ApVec(j,iElt),j=1,5)
                  Else If (iApodType==linApod) Then
                    ! read in xc,yc,nVertex,width
                    READ(VALUE,*) (ApVec(j,iElt),j=1,4)
                  Else If (iApodType==gaussApod) Then
                    ! read in xc,yc,nVertex,width,stdv
                    READ(VALUE,*) (ApVec(j,iElt),j=1,5)
                  End If

	          nVertex=ApVec(3,iElt)
                  ! ApVec(4,iElt) = width of tapered region
                  If (ApVec(4,iElt)<1d-10) Then
                    write(*,*)'** Error in MACOS Rx: tapered width
     &                for polygonal aperture is too small or zero!'
                    GOTO 99
                  End If
                  LTaperedElt(iElt)=.true.
	          print*,'** iElt = ',iElt,', is Tapered_Polygonal ...'
	          !
	          Do k=1,nVertex
                    READ(RX_FID, 123, END=99)STRING
                    READ(STRING,*)(PolyApVtx(j,k,iElt),j=1,2)
                  End Do
                  Call SetCvxPolyApBound(iElt)

	        ELSE IF (ApType(iElt).EQ.9) THEN
		  READ(VALUE,*) (ApVec(j,iElt),j=1,5)
	          !print*,' **SunBlades, ApVec =',ApVec(1:5,iElt)
                END IF
                ApVec_FLG=.TRUE.
              ELSE
                WRITE(*,*)'  Must input ApType before ApVec'
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'PolyApVec',9).OR.
     &               LCMP(VAR_NAM,'Poly3DApVec',11)) THEN
              LApVecDefined(iElt)=.TRUE.
              IF (.not.ApType_FLG) THEN
               WRITE(*,*)
     &         ' Must specify ApType before PolyApVec, iElt =',
     &                   iElt
	        LOAD_SUCCESS=.false.; GO TO 70
              ELSE IF (.not.xObs_FLG) THEN
               WRITE(*,*)' Must specify xObs before PolyApVec, iElt =',
     &                   iElt
	        LOAD_SUCCESS=.false.; GO TO 70
              ELSE
                !READ(VALUE,*) (ApVec(j,iElt),j=1,4)
                !nVertex=ApVec(4,iElt) ! # of vertices of polygon
                READ(VALUE,*) ApVec(1,iElt)
                nVertex=ApVec(1,iElt) ! # of vertices of polygon
                If (nVertex < 3) Then
                  write(*,*)
     &              '** Error in MACOS Rx: a polygon aperture'//
     &              ' must have at least three vertices!'
                  GO TO 99
                End If
                Do k=1,nVertex
                  READ(RX_FID, 123, END=99)STRING
                  READ(STRING,*)(PolyVtx3D(j,k),j=1,3)
                End Do
                Call SetCvxPolyApVtx(iElt,nVertex,
     &                               PolyVtx3D(:,1:nVertex))
                Call SetCvxPolyApBound(iElt)
              END IF
              ApVec_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'nObs',4)) THEN
              READ(VALUE,*) nObs(iElt)
              nObs_FLG=.TRUE.
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'ObsType',7)) THEN
              IF (nObs_FLG) THEN
                DO iObs=1,nObs(iElt)
                  IF (iObs.GT.1) THEN
                    READ(RX_FID, 123, END=99) STRING
                    CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,RX_FID)
                  END IF
                  TypeString=VALUE
	          IF (LCMP(VALUE,ObsTypeName(-8),7).OR.
     &                     LCMP(VALUE,'-8',2)) THEN
                    ! Negtive convex polygon obscutation -jzlou
                    ObsType(iObs,iElt)=-8
	            ObsType_FLG(iObs)=.TRUE.
	          ELSE IF (LCMP(VALUE,ObsTypeName(-6),4).OR.
     &                LCMP(VALUE,'-6',2)) THEN
	            ObsType(iObs,iElt)=-6
	            ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-5),4).OR.
     &                LCMP(VALUE,'-5',2)) THEN
                    ObsType(iObs,iElt)=-5
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-4),4).OR.
     &                     LCMP(VALUE,'-4',2)) THEN
                    ObsType(iObs,iElt)=-4
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-3),4).OR.
     &            LCMP(VALUE,'-3',2)) THEN
                    ObsType(iObs,iElt)=-3
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-2),4).OR.
     &            LCMP(VALUE,'-2',2)) THEN
                    ObsType(iObs,iElt)=1
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-1),4).OR.
     &                     LCMP(VALUE,'-1',2)) THEN
                    ObsType(iObs,iElt)=-1
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(1),4).OR.
     &                     LCMP(VALUE,'1',1)) THEN
                    ObsType(iObs,iElt)=1
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(2),4).OR.
     &                     LCMP(VALUE,'2',1)) THEN
                    ObsType(iObs,iElt)=2
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(3),4).OR.
     &                     LCMP(VALUE,'3',1)) THEN
                    ObsType(iObs,iElt)=3
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(4),4).OR.
     &                     LCMP(VALUE,'4',1)) THEN
                    ObsType(iObs,iElt)=4
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(5),4).OR.
     &                     LCMP(VALUE,'5',1)) THEN
                    ObsType(iObs,iElt)=5
                    ObsType_FLG(iObs)=.TRUE.
	          ELSE IF (LCMP(VALUE,ObsTypeName(6),4).OR.
     &                     LCMP(VALUE,'6',1)) THEN
	   	    ObsType(iObs,iElt)=6
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(7),4).OR.
     &                     LCMP(VALUE,'7',1)) THEN
                    ! Tapered elliptical obscutation -jzlou
                    ObsType(iObs,iElt)=7
                    ObsType_FLG(iObs)=.TRUE.
                    LTaperedElt(iElt)=.TRUE.
	          ELSE IF (LCMP(VALUE,ObsTypeName(8),7).OR.
     &                     LCMP(VALUE,'8',1)) THEN
                    ! Convex polygon obscutation -jzlou
                    ObsType(iObs,iElt)=8
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(9),7).OR.
     &                     LCMP(VALUE,'9',1)) THEN
                    ! RotRectangle obscuration
                    ObsType(iObs,iElt)=9
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-9),9).OR.
     &                     LCMP(VALUE,'-9',2)) THEN
                    ! RotNegRectangle obscuration
                    ObsType(iObs,iElt)=-9
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(10),7).OR.
     &                     LCMP(VALUE,'10',2)) THEN
                    ! RotEllipse obscuration
                    ObsType(iObs,iElt)=10
                    ObsType_FLG(iObs)=.TRUE.
                  ELSE IF (LCMP(VALUE,ObsTypeName(-10),9).OR.
     &                     LCMP(VALUE,'-10',3)) THEN
                    ! RotNegEllipse obscuration
                    ObsType(iObs,iElt)=-10
                    ObsType_FLG(iObs)=.TRUE.
                  END IF

	          ! Read in ObsVec line
                  READ(RX_FID, 123, END=99) STRING
                  CALL GET_EQ(STRING,VAR_NAM,VALUE,IND,ARR_FLG,RX_FID)

	          ! Check whether ObsVec are 3D vertices
	          If (LCMP(VAR_NAM,'PolyObsVec',10).OR.
     &                LCMP(VAR_NAM,'Poly3DObsVec',12)) THEN
                    ! 3D vertices in global frame
                    IF (.not.xObs_FLG) THEN
                     WRITE(*,*)
     &               ' Must specify xObs before PolyObs3DVec, iElt=',
     &                         iElt
	              LOAD_SUCCESS=.false.; GO TO 70
	            ELSE IF (.not.ObsType_FLG(iObs)) THEN
                      WRITE(*,*)
     &                ' Must specify ObsType before PolyObsVec'
                      LOAD_SUCCESS=.false.; GO TO 70
                    END IF

                    IF (ABS(ObsType(iObs,iElt)).EQ.8) THEN
                      !READ(VALUE,*)(ObsVec(j,iObs,iElt),j=1,4)
                      !nVertex=ObsVec(4,iObs,iElt) ! # of vertices of polygon
                      READ(VALUE,*) ObsVec(1,iObs,iElt)
                      nVertex=ObsVec(1,iObs,iElt) ! # of vertices of polygon
                      If (nVertex < 3) Then
                        write(*,*)
     &                   '** Error in MACOS Rx: a polygon obscuration'//
     &                   ' must have at least three vertices!'
                        GO TO 99
                      End If
                      Do k=1,nVertex
                        READ(RX_FID, 123, END=99) STRING
                        READ(STRING,*)(PolyVtx3D(j,k),j=1,3)
                      End Do
                      Call SetCvxPolyObsVtx(iElt,nVertex,
     &                                      PolyVtx3D(:,1:nVertex),iObs)
                      Call SetCvxPolyObsBound(iObs,iElt)
	              !print*,' *** msmacosio: iObs,nVertex =',iObs,nVertex
	              !do i=1,nVertex
	              ! print*,' ** PolyVtx3D(1:3,i) =',PolyVtx3D(1:3,i)
	              !end do
                    ELSE
                      WRITE(*,*)' **Incorrect ObsType for PolyObs3DVec'
                      GO TO 99
                    END IF
                    ObsVec_FLG(iObs)=.TRUE.
                  Else
	            ! not 3D vertices, so as before
                    IF (ABS(ObsType(iObs,iElt)).EQ.1) THEN
                      ! Circle or NegCircle
                      READ(VALUE,*)(ObsVec(j,iObs,iElt),j=1,3)
                    ELSE IF (ABS(ObsType(iObs,iElt)).EQ.4) THEN
                      ! Triangle
                      READ(VALUE,*)(ObsVec(j,iObs,iElt),j=1,6)
                    ELSE IF (ABS(ObsType(iObs,iElt)).EQ.7) THEN
                      ! Tapered ellipse, get: a, b, xc, yc, dr
                      ! where a, b are half axises of ellipse,
                      ! xc, yc are center of ellipse, and dr
                      ! is percentage of 'a' and 'b' to be tapered
                      ! from boundary of ellipse. -jzlou
                      READ(VALUE,*)(ObsVec(j,iObs,iElt),j=1,5)
                    ELSE IF (ABS(ObsType(iObs,iElt)).EQ.8) THEN
                      ! Polygon type
                      READ(VALUE,*) (ObsVec(j,iObs,iElt),j=1,3)
                         ! # of vertices of polygon
                      nVertex=ObsVec(3,iObs,iElt)
                      IF (nVertex < 3) THEN
                        write(*,*)
     &                   '** Error in MACOS Rx: a polygon obscuration'//
     &                   ' must have at least three vertices!'
                       GO TO 99
                      END IF
                      Do k=1,nVertex
                        READ(RX_FID, 123, END=99)STRING
                        READ(STRING,*)(PolyObsVtx(j,k,iObs,iElt),j=1,2)
                      End Do
                      Call SetCvxPolyObsBound(iObs,iElt)
                    ELSE IF (ABS(ObsType(iObs,iElt)).EQ.9) THEN
                      ! RotRectangle type
                      ! Read in W,H,dx,dy,rot_angle
                      READ(VALUE,*) (ObsVec(j,iObs,iElt),j=1,5)
                    ELSE IF (ABS(ObsType(iObs,iElt)).EQ.10) THEN
                      ! RotEllipse type
                      ! Read in a,b,dx,dy,rot_angle
                      READ(VALUE,*) (ObsVec(j,iObs,iElt),j=1,5)
                    ELSE
		      ! For other obs shapes
                      READ(VALUE,*)(ObsVec(j,iObs,iElt),j=1,4)
                    END IF
                    ObsVec_FLG(iObs)=.TRUE.
                  End If
                END DO
              ELSE
                WRITE(*,*)' Must input nObs before ObsType and ObsVec'
                GO TO 99
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'xObs',4)) THEN
              READ(VALUE,*)(xObs(j,iElt),j=1,3)
              IF (DMAG(xObs(1,iElt),3) .EQ. 0.0D0) THEN
                WRITE(*,124)' xObs'
              ELSE
                xObs_FLG = .TRUE.; LxObsDef(iElt)=.TRUE.
              END IF
              GO TO 61

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'EdgeSensors',10)) THEN
              READ(VALUE,*) EdgeSensor(iElt) ! number of edge-sensors for iElt
              If (EdgeSensor(iElt)<1) Then
                print*,
     &      '**msmacosio error: EdgeSensor must be positive integer!'
                stop
              End If
              ! Read in position and measuring direction of all edge-sensors,
              ! in global frame
              Do k=1,EdgeSensor(iElt)
                !READ(RX_FID, 123, END=99)STRING
                READ(RX_FID, *)EdgeSensVec(1:9,k,iElt) ! (position + two directions)
                If (iElt==12) Then
                  print*,'EdgeSensVec(1:9,k,12) =',
     &                   EdgeSensVec(1:10,k,12)
	        End If
	      End Do
	      GO TO 61

            ELSE IF (LCMP(VAR_NAM,'TElt',4)) THEN
              IF (nECoord_FLG) THEN
                READ(VALUE,*)(TElt(1,j,iElt),j=1,nECoord(iElt))
                DO 63 k=2,6
                  READ(RX_FID, 123, END=99)STRING
                  READ(STRING,*)(TElt(k,j,iElt),j=1,nECoord(iElt))
 63             CONTINUE
                TElt_FLG=.TRUE.
                GO TO 61
              ELSE
                WRITE(*,*)'  must input nECoord before TElt'
                GO TO 99
              END IF

            ! ------------------------------------------------------------------------------
	    ELSE IF (LCMP(VAR_NAM,'CopyElt',7) .OR.
     &               LCMP(VAR_NAM,'DuplElt',7)) THEN
              ! Make a duplicate element at given element id, for douple-pass
              READ(VALUE,*) dupElt  ! duplicate Element id
              If (dupElt>0 .AND. dupElt<nElt .AND. iElt>dupElt) Then
                Call MakeDupElt(iElt,dupElt,jGridSrf)
                DupElt_FLG=.TRUE.
              Else
                WRITE(*,*)'  ** At Elt ',iElt,
     &                    ', duplicate Elt id ', dupElt,' is invalid'
              End If
              GO TO 61

            ! ------------------------------------------------------------------------------
            ELSE IF (LCMP(VAR_NAM,'iElt',4)) THEN
	      IF (DupElt_FLG) THEN
                ! If the element previsouly read in is a duplicate element,
                ! don't call ChkDf2 to check input consistency for that element.
                ELEM_OK = .TRUE.
              ELSE
                CALL ChkDf2(iElt,ELEM_OK)
              END IF
              IF (ELEM_OK) THEN
                GO TO 62
              ELSE
                WRITE(*,*)'  Input file not properly loaded'
                LOAD_SUCCESS = .FALSE.
                GO TO 99
              END IF

            ! ------------------------------------------------------------------------------
            ! Exit criteria for reading element definitions
            ! ------------------------------------------------------------------------------

            ELSE IF (LCMP(VAR_NAM,'nOutCord',8)) THEN

	      nElt = iElt

              CALL ChkDf2(iElt, LOAD_SUCCESS)
              IF (LOAD_SUCCESS) GO TO 70

              ! Failure case
              GO TO 99

            ! ------------------------------------------------------------------------------
            ELSE
              WRITE(*,*)' MBFile6: Unidentified string =',VAR_NAM,
     &                  ', iElt=',iElt
              GO TO 61
            END IF

 62	    CONTINUE
            ! ------------------------------------------------------------------------------
            ! Prep for next Element definition
            !    since Element Exit criteria was not yet triggered
            ! ------------------------------------------------------------------------------

            ! check array size
            IF ((iElt+1).GT.mElt) THEN
              nElt         = 0
              LOAD_SUCCESS = .FALSE.

              WRITE(*,*) ' Too many elements defined => increase allocation.'
              GO TO 99
            END IF

	  END DO  ! End of "DO ... END DO" loop for reading elements

