C  code for EFL command

            IF (ifLoad) THEN
              ifRayPrint=.FALSE.
              ifBuild=.FALSE.
              ifTrace=.FALSE.
              ifPropagate=.FALSE.

              IF (isPointSource()) THEN
                WRITE(*,*)' **Need source at infinity'//
     &                    ' to run EFL command!'
                GO TO 1
              END IF
              CALL SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
     &                              ifPropagate)

	      ! Trace chief ray
              iStartRay=1; iEndRay=1
              iStartElt=0; iEndElt=nElt
              CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &                    ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &                    iStartRay,iEndRay,RayMap,SegMap,
     &                    ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &                    StartIndRef,StartExtinc,
     &                    DrawRayVec,DrawEltVec,nDrawElt(1),
     &                    xDraw,yDraw)
              if (.not.LRayOK(iStartRay)) then
                WRITE(*,*)' **Tracing chief ray failed, abort command!'
                GO TO 1
              end if
              ! Save chief ray angle/direction to focal plane
              RayDirBuff(1:3,1)=RayDir(1:3,iStartRay)

              ! Run FindMarginalRay routine to identify an off-axis ray near
              ! aperture edge; must be able to trace this ray to the focal plane.
              off_axis_rayOK=.false.
              Do itry=1,100
                !Call random_number(srp(1:2))  ! [0,1]x[0,1] uniform distr
                !srp=0.4*srp+0.5   ! fall in [0.5,0.9]x[0.5,0.9]
                Call FindMarginalRay(npts,nElt,Dicr,Djcr,dxSource,ifPol,
     &                   ifPropagate,srp(1),srp(2),srid,.FALSE.)
                ! Ray coords are returned in srp. -jzlou
                rayHt= dsqrt(srp(1)*srp(1)+srp(2)*srp(2))

                ! Trace a parallel off-axis ray returned by FSR routine
                !print*,' ** srid =',srid, ', nElt =',nElt; pause

                ! source must be reset here, otherwise when the same ray is traced
                ! twice, once in FSR above and once below, the result will be incorrect.  (ToDo)
                ! Need take a closer look at CTRACE ... -jzlou
                CALL SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
     &                                ifPropagate)

                iStartRay=srid; iEndRay=srid
                iStartElt=0; iEndElt=nElt
                ifTrace=.FALSE.; ifRayPrint=.FALSE.
                CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &                 ifLNegOK,Cmatrix,iStartElt,iEndElt,
     &                 iStartRay,iEndRay,RayMap,SegMap,
     &                 ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &                 StartIndRef,StartExtinc,
     &                 DrawRayVec,DrawEltVec,nDrawElt(1),
     &                 xDraw,yDraw)
                if (LRayOK(iStartRay)) then
                  off_axis_rayOK=.true.
                  exit
                end if
              End Do  ! 100 trys

              if (.not.off_axis_rayOK) then
                WRITE(*,*)
     &            ' **Failed to find a valid off-aix ray to trace,'//
     &            'abort command!'
                GO TO 1
              end if
              ! Save off-axis ray angle/direction to focal plane
              RayDirBuff(1:3,2)=RayDir(1:3,iStartRay)

              ! Calculate EFL
              costh=RayDirBuff(1,1)*RayDirBuff(1,2)+
     &              RayDirBuff(2,1)*RayDirBuff(2,2)+
     &              RayDirBuff(3,1)*RayDirBuff(3,2)
              if (costh > 1d0) costh = 1d0;

              If (abs(costh-1d0)<1d-10) Then
                print*,'**Warning: EFL = infinity'
                efl=1d22
              Else
                efl=rayHt/dsin(dacos(costh))
              End If
           ELSE
             WRITE(*,150)  ! 'load Rx first' warning
           END IF

