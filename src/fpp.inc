C
C Find a pupil position near a given element surface in the optical system
C

        SUBROUTINE FPP(npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &    ifBuild,ifDMBuild,ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,
     &    iObsOpt,iEndElt,iCurWFElt,iCurRayElt,Dicr,Djcr,dxSource,
     &    nBadRays,StartIndRef,StartExtinc,Cmatrix,
     &    DrawRayVec,DrawEltVec,iElt,zp,psip,CrossPt,
     &    ifCentroid,CRSpot,CentroidSpot,dCRCentroid)

       use param_mod
       use src_mod
       use cfiles_mod
       use elt_mod
       use sourcsub_mod
       use traceutil_mod, only : Rx_FEXCentrFlg

       IMPLICIT NONE

       LOGICAL :: ifPol,ifGrid,ifPropagate,ifRayPrint,ifLNegOK,
     &     ifBuild,ifTrace,ifRegrid(mElt),ifIntSrfSet(mElt),
     &     ifDMBuild,ifDraw,ifCentroid

        INTEGER i,iSpot,iRay,mpt2,npts,iCurWFElt,iCurRayElt,
     &    nBadRays,iStartElt,iEndElt,iObsOpt,iStartRay,
     &    iEndRay,nDrawElt,DrawEltVec(mDrawElt,mDrawRay),
     &    iElt,iEm1

       SREAL Cmatrix(7,mCm,bRay),DrawRayVec(2,mDrawElt,mDrawRay)

       REAL*8 Dicr,Djcr,dxSource,StartIndRef,StartExtinc,
     &    xDraw(3),yDraw(3),CentroidSpot(3),dCRCentroid(3),CRSpot(3),
     &    DV1(9),S1,th(3),
     &    oldChfRayDir(3),oldChfRayPos(3),
     &    zp,psip(3),CrossPt(3),cr1dir(3),cr1pos(3),cr2dir(3),cr2pos(3),
     &    cr2indir(3),cr2inpos(3),CrossDist


C      *** Find a pupil location after tracing rays to a given surface ***

        CALL SetSourceRayGrid(nMinPts,Dicr,Djcr,dxSource,ifPol,
     &                        ifPropagate)

       iStartRay = 1
       iEndRay = 1
       iStartElt = 0
       iEndElt = iElt

       CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &          ifLNegOK,Cmatrix,
     &          iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &          ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &          StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)

       CALL DEQUATE(cr1dir,RayDir(1,1),3)
       CALL DEQUATE(cr1pos,RayPos(1,1),3)
       !print*,' ** FPP at destination surface: cr1pos =', cr1pos
       !print*,' ** FPP at destination surface: cr1dir =', cr1dir

       ! Rotate chief ray about STOP by a small angle
       CALL DSMPROD(th,xGrid,5d-6,3)
       CALL DEQUATE(cr2indir,ChfRayDir,3)
       CALL DEQUATE(cr2inpos,ChfRayPos,3)
       CALL RotChfRay(th,cr2indir,cr2inpos,StopPos)
       CALL DEQUATE(RayDir(1,1),cr2indir,3)
       CALL DEQUATE(RayPos(1,1),cr2inpos,3)

       CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &           ifLNegOK,Cmatrix,
     &           iStartElt,iEndElt,iStartRay,iEndRay,RayMap,SegMap,
     &           ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &           StartExtinc,DrawRayVec,DrawEltVec,nDrawElt,xDraw,yDraw)
       CALL DEQUATE(cr2dir,RayDir(1,1),3)
       CALL DEQUATE(cr2pos,RayPos(1,1),3)
       !print*,' ** FPP at destination surface: cr2pos =', cr2pos
       !print*,' ** FPP at destination surface: cr2dir =', cr2dir

       CALL FindCrossPt(cr1pos,cr1dir,cr2pos,cr2dir,zp,
     &                  CrossDist,CrossPt)

       WRITE(*,101) CrossPt
 101   FORMAT(1P,'   *** FPP: Pupil Location = ',3D17.9)
       END SUBROUTINE FPP

