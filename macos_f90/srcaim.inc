
	      print*,'B4 chief ray aiming: StopPos =',StopPos(1:3)

	      ! Chief ray aiming
              ! 'Clean up' xGrid, yGrid first, as in source routines
              CALL DXPROD(zGrid,xGrid,yGrid)
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

	      print*,'After chief ray aim: ChfRayDir =',ChfRayDir(1:3)
	      print*,'After chief ray aim: ChfRayPos =',ChfRayPos(1:3)

