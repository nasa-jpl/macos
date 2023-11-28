C
C Code for generation of reference rays direction - used for collimating a beam
C

	Subroutine getRefRays(nRefRay,refRay,npts,ifPol,ifGrid,
     &                ifPropagate,ifDraw,ifLNegOK,ifBuild,ifDMBuild,
     &                ifTrace,ifRegrid,ifIntSrfSet,ifRayPrint,iObsOpt,
     &                iCurWFElt,iCurRayElt,iEndElt,Dicr,Djcr,
     &                dxSource,nBadRays,StartIndRef,StartExtinc,
     &                Cmatrix,DrawRayVec,DrawEltVec)
	use param_mod
	use elt_mod, only : LRayOK,mRefRay
	use tracesub_mod, only : CRTrace

        Implicit none

        Integer :: nRefRay,refRay(mRefRay),iEndElt,npts,iObsOpt,
     &             iCurWFElt,
     &             iCurRayElt,nBadRays,iStartElt,
     &             DrawEltVec(mDrawElt,mDrawRay)
        Real*8 ::  Dicr,Djcr,dxSource,StartIndRef,StartExtinc,srp(2)
	SREAL ::   Cmatrix(7,mCm,bRay),DrawRayVec(2,mDrawElt,mDrawRay)
	Logical :: ifPol,ifGrid,ifPropagate,
     &             ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,
     &             ifRegrid(mElt),ifIntSrfSet(mElt),ifRayPrint

	nRefRay = 0

        ! Trace chief ray
        iStartElt=0; ifRayPrint=.FALSE.
        CALL CRTrace(npts,ifPol,ifGrid,ifPropagate,
     &         ifDraw,ifLNegOK,ifBuild,ifDMBuild,ifTrace,
     &         ifRegrid,
     &         ifIntSrfSet,ifRayPrint,iObsOpt,iCurWFElt,
     &         iCurRayElt,iEndElt,Dicr,Djcr,dxSource,
     &         nBadRays,StartIndRef,StartExtinc,
     &         Cmatrix,DrawRayVec,DrawEltVec)
	!print*,' **Aft CRTrace, npts =',npts
	ifRayPrint=.TRUE.
        If (LRayOK(1)) Then
	  nRefRay = nRefRay+1; refRay(nRefRay)=1
        End If

        srp(1)=1.0; srp(2)=0.0
        Call FSR(npts,iEndElt,Dicr,Djcr,dxSource,ifPol,
     &           ifPropagate,srp(1),srp(2),refRay(2),.TRUE.)
 	!
        srp(1)=0.0; srp(2)=1.0
        Call FSR(npts,iEndElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &           srp(1),srp(2),refRay(3),.TRUE.)
	!
        srp(1)=-1.0; srp(2)=0.0
        Call FSR(npts,iEndElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &           srp(1),srp(2),refRay(4),.TRUE.)
	!
        srp(1)=0.0; srp(2)=-1.0
        Call FSR(npts,iEndElt,Dicr,Djcr,dxSource,ifPol,ifPropagate,
     &           srp(1),srp(2),refRay(5),.TRUE.)
	nRefRay = nRefRay+4
	print*,' ** getRefRays: refRay(1:5) =', refRay(1:5)
        End Subroutine getRefRays

C  **********************************************************************

#if 1
	!Subroutine TraceRefRays(refRay,nRefRay,iLastElt,refRayDir)
	Subroutine TraceRefRays(iLastElt)
	use param_mod
	!use src_mod, only : SegMap,RayMap,RayIndex
	!use elt_mod, only : zSource,LRayPass,CurRefRayDir,mRefRay
	use src_mod
	use elt_mod
	use macos_mod
#ifdef CMACOS
        use macos_vars_mod
#else
        use smacos_vars_mod
#endif
	use sourcsub_mod
	use tracesub_mod, only : ctrace

	Implicit none

	Integer :: iLastElt
	!Integer :: iLastElt,nRefRay,refRay(mRefRay) ! chief plus four marginal rays
        !Real*8 ::  refRayDir(3,nRefRay)

	DO iRay=1,nRefRay
	  ifBuild=.FALSE.
          ifTrace=.FALSE.
          ifPropagate=.FALSE.
          CALL SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
     &                              ifPropagate)

	  iStartElt = 0
          iStartRay = refRay(iRay)
          iEndRay   = iStartRay

          CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &           ifLNegOK,Cmatrix,iStartElt,iLastElt,iStartRay,iEndRay,
     &           RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,
     &           nBadRays,StartIndRef,StartExtinc,DrawRayVec,
     &           DrawEltVec,nDrawElt(1),xDraw,yDraw)
	  refRayDir(1:3,iRay)=CurRefRayDir(1:3)

	END DO

	End Subroutine TraceRefRays
#endif

