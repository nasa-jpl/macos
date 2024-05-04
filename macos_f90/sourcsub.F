C***********************************************************************
C	Begin file sourcesub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module sourcsub_mod
          use Kinds
	  use math_mod

	  Integer, parameter :: SADir=-1, SAPos=-2
	  Integer :: SAOpt

          private :: FSEG,HSEG,PSEG,SetColSrcRay,SetPtSrcRay,ColSource,PtSource

          real(pr), parameter :: MaxPntSrcDist = 1e10_pr    ! if |zSource| <= MaxPntSrcDist ==> finite Obj.


        CONTAINS

        ! --------------------------------------------------------------
        !
        ! Returns "effective" Src. position as fnct. of ChfRay & zSource
        !
        FUNCTION SourcePos() RESULT(SrcPos)
          USE src_mod, ONLY: ChfRayPos, ChfRayDir
          USE elt_mod, ONLY: zSource

          IMPLICIT NONE
          REAL(pr) :: SrcPos(size(ChfRayPos))
          ! - - - - - - - - - - - - - - - - - - - - - - - - -

          IF (ABS(zSource) <= MaxPntSrcDist) THEN
            SrcPos(:) = ChfRayPos(:) + zSource*ChfRayDir(:)
          ELSE
            SrcPos(:) = ChfRayPos(:)
          END IF

        END FUNCTION SourcePos
        !
        ! --------------------------------------------------------------
        !
        logical function isCollimatedSource()
          use elt_mod, only: zSource
          implicit none
          ! - - - - - - - - - - - - - - - - - - - - - - - - -

          isCollimatedSource = (ABS(zSource) > MaxPntSrcDist)

        end function isCollimatedSource
        !
        ! --------------------------------------------------------------
        !
        logical function isPointSource()
          use elt_mod, only: zSource
          implicit none
          ! - - - - - - - - - - - - - - - - - - - - - - - - -

          isPointSource = (ABS(zSource) <= MaxPntSrcDist)

        end function isPointSource
        !
        ! --------------------------------------------------------------
        !
        ! Determine ray grid at source (initial)  => "ifGrid" is updated
        !
        ! Note: in SMACOS usage "ifGrid" was not always updated; hence,
        !       defined it as optional
        ! ToDo:
        !   [ ] in those instances where it was not updated, is it a bug?
        !
        SUBROUTINE SetSourceRayGrid(npts,Dicr,Djcr,dxSource,ifPol,
     &                              ifPropagate,ifGrid)
          USE elt_mod,   ONLY: zSource
        !   use macos_mod, only: ifGrid

          IMPLICIT NONE
          integer,           intent(in)   :: npts
          logical,           intent(in)   :: ifPol, ifPropagate
          real(pr),          intent(out)  :: Dicr, Djcr, dxSource
          logical, optional, intent(inout):: ifGrid
          ! - - - - - - - - - - - - - - - - - - - - - - - - -

          IF (ABS(zSource) .GE. MaxPntSrcDist) THEN
            CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
          ELSE
            CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
          END IF

          IF (PRESENT(ifGrid)) ifGrid = .TRUE.

        END SUBROUTINE SetSourceRayGrid
        !
        ! --------------------------------------------------------------
        !
        SUBROUTINE SetSourceRay(npts,Dicr,Djcr,dxSource,
     &                          xo,yo,ifPol,ifPropagate)
          USE elt_mod, ONLY: zSource

          implicit none
          integer,  intent(in)   :: npts
          real(pr), intent(out)  :: Dicr, Djcr, dxSource
          real(pr), intent(in)   :: xo, yo
          logical,  intent(in)   :: ifPol, ifPropagate
          ! - - - - - - - - - - - - - - - - - - - - - - - - -

          IF (abs(zSource) .GE. MaxPntSrcDist) THEN
            CALL SetColSrcRay(npts,Dicr,Djcr,dxSource,
     &                        xo,yo,ifPol,ifPropagate)
          ELSE
            CALL SetPtSrcRay(npts,Dicr,Djcr,dxSource,
     &                       xo,yo,ifPol,ifPropagate)
          END IF

        END SUBROUTINE SetSourceRay
        !
        ! --------------------------------------------------------------
        !

! Argument npts is declared in macos.F, and read in from Rx
! for command-line macos.
! Arguments Dicr,Djcr,dxSource are declared in macos.F and
! defined in ColSource
! Argument ifPol and ifPropagate are input arguments only.

	SUBROUTINE ColSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
        use Kinds
	use Constants, only: pio2=>HALFPI, tan30=>Tan30
	use param_mod
        use src_mod
        use elt_mod

	use traceutil_mod

	IMPLICIT NONE

	! When ColSource (and PtSource) is called by 'DRAW' command,
	! the 'npts' argument = nDrawRay -jzlou

	LOGICAL ifPol,ifPropagate,ifRay,inFSEG

	INTEGER npts,iRay,i,j,k,iSeg,iWFA,iwf,jwf,irayid,nDrawRay,iRing

	REAL*8 A,O2,A2,ds,x,y,r,r2,Tbg(3,2),dpr(3),ddr(3),
     &	       Dicr,Djcr,dxSource,w2,ciy,rx2,ry2,dRad(2),
     &	       Sy,yg,Sx,xg,Sg,S,S1,SegX2(3),D1(9),dinset

	COMPLEX*16 :: SC,SEx0,SEy0

	nDrawRay=MIN(mRay,mDrawRay)

	w2=0.5d0*width
	ciy=1.5d0*tan30*width

	nRayToSeg(1:mRay)=0

	CALL DXPROD(zGrid,xGrid,yGrid)
	IF (DDOTC(zGrid,ChfRayDir).GT.0d0) THEN
	  zGrid = ChfRayDir
	ELSE
	  zGrid = -ChfRayDir
	END IF
	CALL DUNITIZE(zGrid)  ! nomalize z
	CALL DXPROD(yGrid,zGrid,xGrid)  ! cross-product z and x to get y
	CALL DUNITIZE(yGrid)
	CALL DXPROD(xGrid,yGrid,zGrid)
	Tbg(1:3,1)=xGrid(1:3)
	Tbg(1:3,2)=yGrid(1:3)


	IF ((GridType.EQ.GridType_Hex).OR.
     &      (GridType.EQ.GridType_Pie).OR.
     &      (GridType.EQ.GridType_Flower)) THEN
	  ! Hex, Pie or Flower source grid type

	  ! D1 is cross-product of ChfRay and SegXgrid -jzlou
	  CALL DXPROD(D1,ChfRayDir,SegXgrid)
	  CALL DXPROD(SegXgrid,D1,ChfRayDir)
	  CALL DUNITIZE(SegXgrid)

	  ! xGrid,yGrid and SegXgrid are all unit vectors,
	  ! so computed SegX2(1) and SegX2(2) are cosines
          ! of angles of SegXgrid with xGrid and yGrid. -jzlou
	  SegX2(1)=DDOTC(xGrid,SegXgrid)
	  SegX2(2)=DDOTC(yGrid,SegXgrid)
	  SegX2(3)=0d0
	END IF


C  Find the segment for the chief ray
	nGridPts=npts+1  ! no. of grid pts in X (or Y) of source aperture
	iRay=1  ! chief ray

	IF (Obscratn.NE.0d0) THEN
	  O2=(Obscratn/2d0)**2
	ELSE
	  O2=0d0
	END IF

!ToDo: defining aperture for the sake of plotting is wrong => what is the correct physics

	!A2=(Aperture/2d0)**2  ! original
	A2= (Aperture/2d0 + Aperture/DBLE(npts))**2  ! for plotting a 'bottom' ray
	!A2=(Aperture*0.5001)**2  ! original

! nRayToSeg(iRay) - given iRay, (nRayToSeg(iRay)-1) gives the number of
! rays that have hit iSeg before iRay does. -jzlou

        ! GridType from 3 to 5 are segmented
        ! -------------------------------------
        ! GridType: Hex
        ! -------------------------------------
        IF (GridType.EQ.GridType_Hex) THEN
          DO iSeg=1,nSeg   ! nSeg is read from Rx -jzlou
            IF (HSEG(0d0,0d0,width,gap,SegCoord(1,iSeg),SegX2)) THEN
              nRayToSeg(iRay)=nRayToSeg(iRay)+1
              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
            END IF
          END DO
          IF (nRayToSeg(iRay).EQ.0) THEN
            ! if chief ray didn't hit any segment, assign it to seg 1. -jzlou
            nRayToSeg(iRay)=1
            RayToSegMap(iRay,nRayToSeg(iRay))=1
          END IF
        ! ---------------------------------------------
        ! GridType: Pie
        ! ---------------------------------------------
        ELSE IF (GridType.EQ.GridType_Pie) THEN
          DO iSeg=1,nSeg
            IF (PSEG(0d0,0d0,width,gap,SegCoord(1,iSeg),SegX2)) THEN
              nRayToSeg(iRay)=nRayToSeg(iRay)+1
              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
            END IF
          END DO
          IF (nRayToSeg(iRay).EQ.0) THEN
            nRayToSeg(iRay)=1
            RayToSegMap(iRay,nRayToSeg(iRay))=1
          END IF
        ! ---------------------------------------------
        ! GridType: Flower
        ! ---------------------------------------------
        ELSE IF (GridType.EQ.GridType_Flower) THEN
          DO iSeg=1,nSeg
            iRing=SegCoord(1,iSeg)
            If (iRing>0) Then
              dRad(1)=RingRads(iRing-1)*radCtr
              dRad(2)=RingRads(iRing)*radCtr
            End If
            IF (FSEG(0d0,0d0,0d0,radCtr,dRad,nPetals(iRing),gap,
     &          SegCoord(:3,iSeg),SegX2)) THEN
              nRayToSeg(iRay)=1
              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
            END IF
          END DO
          IF (nRayToSeg(iRay).EQ.0) THEN
            nRayToSeg(iRay)=nRayToSeg(iRay)+1
            RayToSegMap(iRay,nRayToSeg(iRay))=1
          END IF
        END IF

	! Statements below seem useless -jzlou
	ddr(1:2)=0d0
	CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	RayDir(1:3,iRay)=ChfRayDir(1:3)
        RayPos(1:3,iRay)=ChfRayPos(1:3)+dpr(1:3)

        SegMap(1:mpts,1:mpts) = 0
        RayID(1:mpts,1:mpts)  = 0

        ! In the code section below, jwf,iwf are indexes to put source ray grid
        ! (nGridPts x nGridPts) at the center of WF grid (mdttl x mdttl). -jzlou

	IF (nGridPts.GT.1) THEN
	  A        = Aperture/2d0
	 !ds       = Aperture/DBLE(npts)  ! original
	  ds       = Aperture/DBLE(npts-1)
          nGridPts = npts
	  dinset   = ds/2d0
          !ds=Aperture/DBLE(nGridPts-1-mod(nGridPts+1,2)) ! nGridPts even or odd
	  !ds=Aperture/DBLE(nGridPts-1)  ! nGridPts odd
	  if (DrawMode) then
	    y=-Aperture/2d0-ds+dinset  !1d0/DBLE(npts)
	  else
	    y=-Aperture/2d0-ds
	  end if

	  DO 10 j=1,nGridPts
	    !jwf=j+(mdttl-nGridPts)/2+1  ! original
	    !jwf=j+(mdttl-nGridPts+1)/2
	    jwf = j+(mdttl-nGridPts+mod(nGridPts,2))/2
	    If (DrawMode .AND. j==nGridPts) Then
	      y=y+ds-2d0*dinset  !2d0/DBLE(npts)
            Else
              y=y+ds
            End If
	    if (DrawMode) then
	      x=-Aperture/2d0-ds+dinset  !1d0/DBLE(npts)
	    else
	      x=-Aperture/2d0-ds
	    end if

	    DO 9 i=1,nGridPts
	      !iwf=i+(mdttl-nGridPts+1)/2  ! original
	      !iwf=i+(mdttl-nGridPts+1)/2
	      iwf=i+(mdttl-nGridPts+mod(nGridPts,2))/2
	      If (DrawMode .AND. i==nGridPts) Then
                x=x+ds-2d0*dinset  !2d0/DBLE(npts)
              Else
                x=x+ds
              End If

              ! ---------------------------------------
              ! GridType: Circular
              ! ---------------------------------------
              ! where O2 = (Obscratn/2d0)**2

	      IF (GridType.EQ.GridType_Circular) THEN
	        r2=y**2+x**2
	        !IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	        IF (((r2.GE.O2).AND.(r2.LE.(A2+1d-9))).AND.(iRay.LT.mRay)) THEN
	          iRay=iRay+1
	          RayIndex(1,iRay) = i
	          RayIndex(2,iRay) = j
	          RayMap(i,j)      = iRay
		  RayID(iwf,jwf)   = iRay
	          RayWFMap(1,iRay) = iwf  ! index for WF grid array 'WFElt'
	          RayWFMap(2,iRay) = jwf
	          ddr(:2)          = (/x, y/)

	          CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	          RayDir(:,iRay) = ChfRayDir(:)        ! collimated source
	          RayPos(:,iRay) = ChfRayPos(:)+dpr(:)
	        END IF
              ! ---------------------------------------
              ! GridType: Square
              ! ---------------------------------------
	      ELSE IF (GridType.EQ.GridType_Square) THEN
	        IF (((ABS(x).LE.A).AND.(ABS(y).LE.A))
     &	            .AND.(iRay.LT.mRay)) THEN
                  iRay             = iRay+1
                  RayIndex(1,iRay) = i
                  RayIndex(2,iRay) = j
                  RayMap(i,j)      = iRay
                  RayID(iwf,jwf)   = iRay
                  RayWFMap(1,iRay) = iwf
                  RayWFMap(2,iRay) = jwf
                  ddr(:2)          = (/x, y/)

                  CALL DMPROD(dpr,Tbg,ddr,3,2,1)
                  RayDir(:,iRay) = ChfRayDir(:)
                  RayPos(:,iRay) = ChfRayPos(:)+dpr(:)
                END IF
              ! ---------------------------------------
              ! GridType: Hex
              ! ---------------------------------------
              ELSE IF (GridType.EQ.GridType_Hex) THEN
                ifRay=.FALSE.
                r2=y**2+x**2

                IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
                  DO iSeg=1,nSeg
                    IF (HSEG(x,y,width,gap,SegCoord(1,iSeg),SegX2)) THEN
!                   IF (HSEG(x,y,width,gap,SegCoord(1,iSeg),SegX2) .OR.
!    &                  (i==1 .and. j==(nDrawRay/2+1)) .OR.
!    &                  (i==(nDrawRay/2+1) .and. j==1)) THEN
                      IF (.NOT.ifRay) THEN
                        iRay             = iRay+1
                        ifRay            = .TRUE.
                        RayIndex(1,iRay) = i
                        RayIndex(2,iRay) = j
                        RayMap(i,j)      = iRay

                        SegMap(i,j)      = iSeg
                        RayID(iwf,jwf)   = iRay
                        RayWFMap(1,iRay) = iwf
                        RayWFMap(2,iRay) = jwf
                        ddr(:2)          = (/x, y/)

                        CALL DMPROD(dpr,Tbg,ddr,3,2,1)
                        RayDir(:,iRay) = ChfRayDir(:)
                        RayPos(:,iRay) = ChfRayPos(:)+dpr(:)
                      END IF

                      nRayToSeg(iRay) = nRayToSeg(iRay)+1
                      RayToSegMap(iRay,nRayToSeg(iRay)) = iSeg
                    END IF
                  END DO
                END IF

              ! ---------------------------------------
              ! GridType: Pie
              ! ---------------------------------------
              ELSE IF (GridType.EQ.GridType_Pie) THEN
                ifRay=.FALSE.
                r2=y**2+x**2
                IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
                  DO iSeg=1,nSeg
                    IF (PSEG(x,y,width,gap,SegCoord(1,iSeg),SegX2)) THEN
                      IF (.NOT.ifRay) THEN
                        iRay             = iRay+1
                        ifRay            = .TRUE.
                        RayIndex(1,iRay) = i
                        RayIndex(2,iRay) = j
                        RayMap(i,j)      = iRay
                        SegMap(i,j)      = iSeg
                        RayID(iwf,jwf)   = iRay
                        RayWFMap(1,iRay) = iwf
                        RayWFMap(2,iRay) = jwf
                        ddr(:2)          = (/x, y/)

                        CALL DMPROD(dpr,Tbg,ddr,3,2,1)
                        RayDir(:,iRay) = ChfRayDir(:)
                        RayPos(:,iRay) = ChfRayPos(:)+dpr(:)
                      END IF
                      nRayToSeg(iRay) = nRayToSeg(iRay)+1
                      RayToSegMap(iRay,nRayToSeg(iRay)) = iSeg
                    END IF
                  END DO
                END IF

              ! ---------------------------------------
              ! GridType: Flower
              ! ---------------------------------------
              ELSE IF (GridType.EQ.GridType_Flower) THEN
                ifRay=.FALSE.
                r2=y**2+x**2

                IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
                  r2 = SQRT(r2)
                  DO iSeg=1,nSeg
                    iRing=SegCoord(1,iSeg)
                    If (iRing>0) Then
                      dRad(1)=RingRads(iRing-1)*radCtr
                      dRad(2)=RingRads(iRing)*radCtr
                    End If
                    IF (FSEG(x,y,r2,radCtr,dRad,nPetals(iRing),gap,
     &                  SegCoord(:3,iSeg),SegX2)) THEN
                      IF (.NOT.ifRay) THEN
                        iRay             = iRay+1
                        ifRay            = .TRUE.
                        RayIndex(1,iRay) = i
                        RayIndex(2,iRay) = j
                        RayMap(i,j)      = iRay
                        SegMap(i,j)      = iSeg
                        RayID(iwf,jwf)   = iRay
                        RayWFMap(1,iRay) = iwf
                        RayWFMap(2,iRay) = jwf
                        ddr(:2)          = (/x, y/)

                        CALL DMPROD(dpr,Tbg,ddr,3,2,1)
                        RayDir(:,iRay) = ChfRayDir(:)
                        RayPos(:,iRay) = ChfRayPos(:)+dpr(:)
                      END IF
                      nRayToSeg(iRay) = nRayToSeg(iRay)+1
                      RayToSegMap(iRay,nRayToSeg(iRay)) = iSeg
                    END IF
                  END DO
                END IF
              END IF
    9	    CONTINUE
   10	  CONTINUE
	  nRay=iRay

        ! ------------------------------
        ! nGridPts == 1
        ! ------------------------------
	ELSE IF ((GridType.EQ.GridType_Hex).OR.
     &           (GridType.EQ.GridType_Pie)) THEN
	  nRay=nSeg
	  DO iRay=1,nRay
	    iSeg=iRay
	    ddr(1)=DBLE(SegCoord(1,iSeg))*w2
	    ddr(2)=ciy*DBLE((SegCoord(2,iSeg)+SegCoord(3,iSeg)))/3D0
	    nRayToSeg(iRay)=nRayToSeg(iRay)+1
	    RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
c	    RayToSegMap(iRay)=iSeg
	    i=SegCoord(1,iSeg)+mpts/2
	    j=((SegCoord(2,iSeg)+SegCoord(3,iSeg))/3)+(mpts/2)
	    iwf=SegCoord(1,iSeg)+mdttl/2
	    jwf=((SegCoord(2,iSeg)+SegCoord(3,iSeg))/3)+(mdttl/2)

	    RayIndex(1,iRay) = i
	    RayIndex(2,iRay) = j
	    RayMap(i,j)      = iRay
	    RayID(iwf,jwf)   = iRay
	    RayWFMap(1,iRay) = iwf
	    RayWFMap(2,iRay) = jwf

	    CALL DMPROD(dpr,Tbg,ddr,3,2,1)
            RayDir(:,iRay) = ChfRayDir(:)
            RayPos(:,iRay) = ChfRayPos(:)+dpr(:)
          END DO

	ELSE
	  nRay=1
	  RayIndex(1,iRay)=1
	  RayIndex(2,iRay)=1
	  RayMap(1,1)=1
	END IF
	RayIndex(1,1)=0
	RayIndex(2,1)=0

c	Dicr=(DBLE(nGridPts)+1d0)/2d0+DBLE((mdttl-nGridPts)/2)
c	Djcr=(DBLE(nGridPts)+1d0)/2d0+DBLE((mdttl-nGridPts)/2+1)
	IF (MOD(nGridPts,2).EQ.0) THEN
	  Dicr=DBLE(mdttl/2)+0.5d0
	  Djcr=DBLE(mdttl/2)+1.5d0
	ELSE
	  Dicr=mdttl/2+1
	  Djcr=Dicr
	END IF
	dxSource=ds


C  Set source amplitude if propagating for polarized light

! RayE(1:3,iRay) is the source electric field for ith ray
! WFElt(i,j,1) = first wavefront field
! jzlou

! BeamTypes: (1) UNIFORM, (2) GAUSSIAN, (3) COS**POWER, (4) DIPOLE

	IF (ifPropagate.OR.ifPol) THEN
	  iWFA=1
	  IF (BeamType.EQ.2) THEN
C  Gaussian beam...
	    S=SQRT(Flux)
	    rx2=rxBeam*rxBeam
	    ry2=ryBeam*ryBeam
	    IF (ifPol) THEN
	      DO 35 k=1,3
	        RayE(k,1)=Ex0*DCMPLX((S*xGrid(k)),0d0)+
     &	          Ey0*DCMPLX((S*yGrid(k)),0d0)
   35	      CONTINUE
	    END IF
	    DO 23 j=1,mdttl
	      Sy=dxSource*(DBLE(j)-Djcr)/ryBeam
	      yg=Sy*Sy
	      DO 22 i=1,mdttl
	        irayid=RayID(i,j)
	        IF (irayid.EQ.0) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          Sx=dxSource*(DBLE(i)-Dicr)/rxBeam
	          xg=Sx*Sx
	          Sg=S*EXP(-(xg+yg))
	          WFElt(i,j,iWFA)=DCMPLX(Sg,0D0)
	          IF (ifPol) THEN
	            DO 31 k=1,3
	              RayE(k,irayid)=Ex0*DCMPLX((Sg*xGrid(k)),0d0)+
     &	                Ey0*DCMPLX((Sg*yGrid(k)),0d0)
   31	            CONTINUE
	          END IF
	        END IF
   22	      CONTINUE
   23	    CONTINUE
	  ELSE IF (BeamType.EQ.3) THEN
C  Cosine beam...
	    rx2=pio2/rxBeam
	    ry2=rx2
	    S=SQRT(Flux)
	    IF (ifPol) THEN
	      DO 34 k=1,3
	        RayE(k,1)=Ex0*DCMPLX((S*xGrid(k)),0d0)+
     &	        Ey0*DCMPLX((S*yGrid(k)),0d0)
   34	      CONTINUE
	    END IF
	    DO 29 j=1,mdttl
	      Sy=dxSource*(DBLE(j)-Djcr)*ry2
	      yg=Sy*Sy
	      DO 28 i=1,mdttl
	        irayid=RayID(i,j)
	        IF ((irayid.EQ.0).OR.(Sy.GT.pio2)) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          Sx=dxSource*(DBLE(i)-Dicr)*rx2
	          IF (Sx.GT.pio2) THEN
	            WFElt(i,j,iWFA)=0D0
	          ELSE
	            xg=Sx*Sx
	            Sg=SQRT(xg+yg)
	            S1=COS(Sg)
	            IF (S1.GE.0d0) THEN
	              S=SQRT(Flux*S1)**CosPower
	              WFElt(i,j,iWFA)=DCMPLX(S,0D0)
	              IF (ifPol) THEN
	                DO 32 k=1,3
	                  RayE(k,irayid)=Ex0*DCMPLX((S*xGrid(k)),0d0)+
     &	                  Ey0*DCMPLX((S*yGrid(k)),0d0)
   32	                CONTINUE
	              END IF
	            ELSE
	              WFElt(i,j,iWFA)=0D0
	            END IF
	          END IF
	        END IF
   28	      CONTINUE
   29	    CONTINUE
	  ELSE
C  Uniform or Dipole beam...
	    S=SQRT(Flux/DBLE(nRay-1))
	    SEx0=S*Ex0
	    SEy0=S*Ey0
	    SC=DCMPLX(S,0D0)
	    IF (ifPol) THEN
	      DO 36 k=1,3
	        RayE(k,1)=SEx0*DCMPLX(xGrid(k),0d0)+
     &	        SEy0*DCMPLX(yGrid(k),0d0)
   36	      CONTINUE
	    END IF
	    DO 25 j=1,mdttl
	      DO 24 i=1,mdttl
	        irayid=RayID(i,j)
	        IF (irayid.EQ.0) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          WFElt(i,j,iWFA)=SC
	          IF (ifPol) THEN
	            DO 33 k=1,3
	              RayE(k,irayid)=SEx0*DCMPLX(xGrid(k),0d0)+
     &	              SEy0*DCMPLX(yGrid(k),0d0)
   33	            CONTINUE
	          END IF
	        END IF
   24	      CONTINUE
   25	    CONTINUE
	  END IF
	END IF

#if 0
! testing code
	open(unit=9,file="src_ap_rays.txt")
	! Test only
	Do iRay=1,nRay
          dRay(1:3)=RayPos(1:3,iRay)-RayPos(1:3,1)
          rx=DDOTC(dRay(1:3),xGrid(1:3))
          ry=DDOTC(dRay(1:3),yGrid(1:3))
	  rx=dRay(1); ry=dRay(2)
	  write(9,*) rx,ry
          !if (abs(rx) < 350.d0 .and. abs(ry) < 20.*dxSource) then
            !print*,'ray id, rx,ry =', iRay, rx,ry
          if (abs(rx) < 350.d0 .and. abs(ry) < 15.*dxSource) then
            print*,'ray id, ry =', iRay, ry
          end if
        End Do
	close(9)
	stop
#endif
	RETURN
	END SUBROUTINE ColSource

C***********************************************************************

	SUBROUTINE PtSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
        use Kinds
	use Constants, only: pio2=>HALFPI
	use param_mod
        use src_mod
        use elt_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifPropagate,ifRay

	INTEGER npts,iRay,i,j,k,iSeg,iWFA,iwf,jwf,irayid,iRing,
     &          fu,m_err_ptsrc

	REAL*8 A,O2,A2,ds,thx,thy,r2,Tbg(3,2),
     &	       Dicr,Djcr,dxSource,th(3),RotPt(3),
     &	       rx2,ry2,Sy,yg,Sx,xg,Sg,S,xray(3),yray(3),zray(3),
     &	       S1,SegX2(3),D1(9)
        Real*8 :: du_loc,dx_loc,zElt_loc,dRad(2)

	COMPLEX*16 SC,SEx0,SEy0
        COMPLEX*16, dimension(:,:), allocatable :: dsc,DWF_loc

	nRayToSeg(1:nRay)=0

        ! zGrid = cross-product of xGrid and yGrid - jzlou
	CALL DXPROD(zGrid,xGrid,yGrid)

	IF (DDOTC(zGrid,ChfRayDir).GT.0d0) THEN
	  CALL DEQUATE(zGrid,ChfRayDir,3)
	ELSE
	  CALL DNEGATE(zGrid,ChfRayDir,3)
	END IF

	CALL DUNITIZE(zgrid)
	CALL DXPROD(yGrid,zGrid,xGrid)
	CALL DUNITIZE(yGrid)
	CALL DXPROD(xGrid,yGrid,zgrid)
	Tbg(1:3,1)=xgrid(1:3)
	Tbg(1:3,2)=ygrid(1:3)

	IF ((GridType.GE.3).AND.(GridType.LE.5)) THEN
          ! 'Hex' or 'Flower' types
	  CALL DXPROD(D1,ChfRayDir,SegXgrid)
	  CALL DXPROD(SegXgrid,D1,ChfRayDir)
	  CALL DUNITIZE(SegXgrid)
	  SegX2(1)=DDOTC(xGrid,SegXgrid)
	  SegX2(2)=DDOTC(yGrid,SegXgrid)
	  SegX2(3)=0d0
	END IF

C  Find the segment for the chief ray
	nGridPts=npts+1
	iRay=1
	IF (Obscratn.NE.0d0) THEN
	  O2=(Obscratn*0.49999d0)**2
	ELSE
	  O2=0d0
	END IF

	A2=(Aperture*0.50001)**2
	!A2= (Aperture/2d0 + Aperture/DBLE(npts))**2  ! for plotting a 'bottom' ray

        ! for segmented apetures (Hex,Pie,Flower) only
        ! -------------------------------------
        ! GridType: Hex
        ! -------------------------------------
	IF (GridType.EQ.GridType_Hex) THEN
	  DO iSeg=1,nSeg
	    IF (HSEG(0d0,0d0,width,gap,SegCoord(1,iSeg),SegX2)) THEN
	      nRayToSeg(iRay)=nRayToSeg(iRay)+1
	      RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	    END IF
 	  END DO
	  IF (nRayToSeg(iRay).EQ.0) THEN
	    nRayToSeg(iRay)=1
	    RayToSegMap(iRay,nRayToSeg(iRay))=1
	  END IF
        ! ---------------------------------------------
        ! GridType: Pie
        ! ---------------------------------------------
        ELSE IF (GridType.EQ.GridType_Pie) THEN
          ! 'Pie' type
          DO iSeg=1,nSeg
            IF (PSEG(0d0,0d0,width,gap,SegCoord(1,iSeg),SegX2)) THEN
              nRayToSeg(iRay)=nRayToSeg(iRay)+1
              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
            END IF
          END DO
          IF (nRayToSeg(iRay).EQ.0) THEN
            nRayToSeg(iRay)=1
            RayToSegMap(iRay,nRayToSeg(iRay))=1
          END IF
        ! ---------------------------------------------
        ! GridType: Flower
        ! ---------------------------------------------
        ELSE IF (GridType.EQ.GridType_Flower) THEN
          DO iSeg=1,nSeg
            iRing=SegCoord(1,iSeg)
            If (iRing>0) Then
              dRad(1)=RingRads(iRing-1)*radCtr
              dRad(2)=RingRads(iRing)*radCtr
            End If
            IF (FSEG(0d0,0d0,0d0,radCtr,dRad,nPetals(iRing),gap,
     &          SegCoord(:3,iSeg),SegX2)) THEN
              nRayToSeg(iRay)=nRayToSeg(iRay)+1
              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
            END IF
          END DO
          IF (nRayToSeg(iRay).EQ.0) THEN
            nRayToSeg(iRay)=1
            RayToSegMap(iRay,nRayToSeg(iRay))=1
          END IF
        END IF


        ! Define user-provided source ray beam - added by jzlou

        ! ---------------------------------------------
        ! GridType: UDG
        ! ---------------------------------------------
        If (GridType.EQ.GridType_UDG) Then
          RayMap(1:mpts,1:mpts)=0
          SegMap(1:mpts,1:mpts)=0  ! not needed, user-provided source rays is non-segmented
          RayID(1:mdttl,1:mdttl)=0

          ! Read in user-provided source beam data; source rays are defined 'zSource'
          ! distance from the point source location, on a spherical surface -jzlou
          fu=12
	  Open(Unit=fu,File='UDSrcBeam.txt',Status='OLD')
	  Read(fu,*) nGridPts,nRay  ! nRay= total number of valid rays
	  Do iRay=1,nRay
            Read(fu,*) i,j,RayPos(1:3,iRay),RayDir(1:3,iRay)
	    if (iRay==1) then
	      ! make sure chief ray is consistent with user data
	      ChfRayPos(1:3)=RayPos(1:3,iRay)
	      ChfRayDir(1:3)=RayDir(1:3,iRay)
	    end if
            jwf=j+(mdttl-nGridPts)/2+1
	    iwf=i+(mdttl-nGridPts+1)/2
            RayIndex(1,iRay)=i; RayIndex(2,iRay)=j;
            RayMap(i,j)=iRay
            RayID(iwf,jwf)=iRay
            RayWFMap(1,iRay)=iwf
            RayWFMap(2,iRay)=jwf
	  End Do
	  Close(fu)
	  GOTO 16  ! no need to go thru other GridType
	End If

        ! ---------------------------------------------

        ! RotPt defines the spatial position of the physical light source,
        ! see MACOS manual.
        RotPt = SourcePos()   ! ChfRayPos+zSource*ChfRayDir

	! Define RayDir and RayPos for chief ray
	th(1:3)=0d0
	CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,RayDir(1,iRay),
     &		   RayPos(1,iRay))

	! Define RayDir and RayPos for all other valid rays
	IF (nGridPts.GT.1) THEN
	  RayMap(1:mpts,1:mpts)=0
	  SegMap(1:mpts,1:mpts)=0
	  RayID(1:mdttl,1:mdttl)=0

	  A=Aperture*0.50001d0
	  ds=Aperture/DBLE(npts)
	  thy=-Aperture/2d0-ds

          DO 10 j=1,nGridPts
            jwf=j+(mdttl-nGridPts)/2+1
            thy=thy+ds
            thx=-Aperture/2d0-ds

            DO 9 i=1,nGridPts
              iwf=i+(mdttl-nGridPts+1)/2
              thx=thx+ds
              ! ---------------------------------------
              ! GridType: Circular
              ! ---------------------------------------
              IF (GridType.EQ.GridType_Circular) THEN
                r2=thy**2+thx**2
                IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
                  iRay             = iRay+1
                  RayIndex(1,iRay) = i
                  RayIndex(2,iRay) = j
                  RayMap(i,j)      = iRay
                  RayID(iwf,jwf)   = iRay
                  RayWFMap(1,iRay) = iwf
                  RayWFMap(2,iRay) = jwf

                  ! Euler vector
                  th(:) = -xGrid(:)*thy + yGrid(:)*thx

                  ! this routine sets ray direction and position -jzlou
                  CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,
     &	                     RayDir(1,iRay),RayPos(1,iRay))
                END IF
              ! ---------------------------------------
              ! GridType: Square
              ! ---------------------------------------
              ELSE IF (GridType.EQ.GridType_Square) THEN
                IF (((ABS(thx).LE.A).AND.(ABS(thy).LE.A))
     &	            .AND.(iRay.LT.mRay)) THEN
	          iRay             = iRay+1
	          RayIndex(1,iRay) = i
	          RayIndex(2,iRay) = j
	          RayMap(i,j)      = iRay
		  RayID(iwf,jwf)   = iRay
	          RayWFMap(1,iRay) = iwf
	          RayWFMap(2,iRay) = jwf

                  th(:) = -xGrid(:)*thy + yGrid(:)*thx
	          CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,
     &	                     RayDir(1,iRay),RayPos(1,iRay))
	        END IF
              ! ---------------------------------------
              ! GridType: Hex
              ! ---------------------------------------
              ELSE IF (GridType.EQ.GridType_Hex) THEN
	        ifRay=.FALSE.
	        r2=thy**2+thx**2
	        IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	          DO 7 iSeg=1,nSeg
		    IF (HSEG(thx,thy,width,gap,SegCoord(1,iSeg),
     &		             SegX2)) THEN
	              IF (.NOT.ifRay) THEN
	                iRay=iRay+1
	                ifRay=.TRUE.
	                RayIndex(1,iRay)=i
	                RayIndex(2,iRay)=j
	                RayMap(i,j)=iRay
	                SegMap(i,j)=iSeg
	                RayID(iwf,jwf)=iRay
	                RayWFMap(1,iRay)=iwf
	                RayWFMap(2,iRay)=jwf
	                DO k=1,3
	                  th(k)=-xGrid(k)*thy+yGrid(k)*thx
	                END DO
	                CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,
     &	                RayDir(1,iRay),RayPos(1,iRay))
	              END IF
	              nRayToSeg(iRay)=nRayToSeg(iRay)+1
	              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	            END IF
  7	  	  CONTINUE
	        END IF
              ! ---------------------------------------
              ! GridType: Pie
              ! ---------------------------------------
              ELSE IF (GridType.EQ.GridType_Pie) THEN
	        ifRay=.FALSE.
	        r2=thy**2+thx**2
	        IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	          DO 8 iSeg=1,nSeg
		    IF (PSEG(thx,thy,width,gap,SegCoord(1,iSeg),
     &		    SegX2)) THEN
	              IF (.NOT.ifRay) THEN
	                iRay=iRay+1
	                ifRay=.TRUE.
	                RayIndex(1,iRay)=i
	                RayIndex(2,iRay)=j
	                RayMap(i,j)=iRay
	                SegMap(i,j)=iSeg
	                RayID(iwf,jwf)=iRay
	                RayWFMap(1,iRay)=iwf
	                RayWFMap(2,iRay)=jwf
	                DO k=1,3
	                  th(k)=-xGrid(k)*thy+yGrid(k)*thx
	                END DO
	                CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,
     &	                           RayDir(1,iRay),RayPos(1,iRay))
	              END IF
	              nRayToSeg(iRay)=nRayToSeg(iRay)+1
	              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	            END IF
    8		  CONTINUE
	        END IF
	      END IF
    9	    CONTINUE
   10	  CONTINUE
	  nRay=iRay
	ELSE
	  nRay=1
	  RayIndex(1,iRay)=1
	  RayIndex(2,iRay)=1
	  RayMap(1,1)=1
	END IF
	RayIndex(1,1)=0
	RayIndex(2,1)=0

c	Dicr=(DBLE(nGridPts)+1d0)/2d0+DBLE((mdttl-nGridPts)/2)
c	Djcr=(DBLE(nGridPts)+1d0)/2d0+DBLE((mdttl-nGridPts)/2+1)
	IF (MOD(nGridPts,2).EQ.0) THEN
	  Dicr=DBLE(mdttl/2)+0.5d0
	  Djcr=DBLE(mdttl/2)+1.5d0
	ELSE
	  Dicr=mdttl/2+1
	  Djcr=Dicr
	END IF
	dxSource=ds

16      CONTINUE  ! end of defining all source rays


C  Set source amplitude if propagating or for polarized light.
C  ifPropagate is true when performing diffraction commands -jzlou

	IF (ifPropagate.OR.ifPol) THEN
	  iWFA=1
	  IF (BeamType.EQ.2) THEN
C  Gaussian beam...
	    S=SQRT(Flux)
	    rx2=rxBeam*rxBeam  ! rxBeam = beam waist in X axis -jzlou
	    ry2=ryBeam*ryBeam  ! ryBeam = ...
	    IF (ifPol) THEN
	      DO 35 k=1,3
	        RayE(k,1)=Ex0*DCMPLX((S*xGrid(k)),0d0)+
     &	        	  Ey0*DCMPLX((S*yGrid(k)),0d0)
   35	      CONTINUE
	    END IF
	    DO 23 j=1,mdttl
	      Sy=dxSource*(DBLE(j)-Djcr)/ryBeam
	      yg=Sy*Sy
	      DO 22 i=1,mdttl
	        irayid=RayID(i,j)
	        IF (irayid.EQ.0) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          Sx=dxSource*(DBLE(i)-Dicr)/rxBeam
	          xg=Sx*Sx
	          Sg=S*EXP(-(xg+yg))
	          WFElt(i,j,iWFA)=DCMPLX(Sg,0D0)
	          IF (ifPol) THEN
	            CALL DEQUATE(zray,RayDir(1,irayid),3)
	            CALL DXPROD(yray,zray,xGrid)
	            CALL DUNITIZE(yray)
	            CALL DXPROD(xray,yray,zray)
	            DO 31 k=1,3
	              RayE(k,irayid)=Ex0*DCMPLX((Sg*xray(k)),0d0)+
     &	                             Ey0*DCMPLX((Sg*yray(k)),0d0)
   31	            CONTINUE
	          END IF
	        END IF
   22	      CONTINUE
   23	    CONTINUE
	  ELSE IF (BeamType.EQ.3) THEN
C  Cosine beam...
	    rx2=pio2/rxBeam
	    ry2=rx2
	    S=SQRT(Flux)
	    IF (ifPol) THEN
	      DO 34 k=1,3
	        RayE(k,1)=Ex0*DCMPLX((S*xGrid(k)),0d0)+
     &	        Ey0*DCMPLX((S*yGrid(k)),0d0)
   34	      CONTINUE
	    END IF
	    DO 29 j=1,mdttl
	      Sy=dxSource*(DBLE(j)-Djcr)*ry2
	      yg=Sy*Sy
	      DO 28 i=1,mdttl
	        irayid=RayID(i,j)
	        IF ((irayid.EQ.0).OR.(Sy.GT.pio2)) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          Sx=dxSource*(DBLE(i)-Dicr)*rx2
	          IF (Sx.GT.pio2) THEN
	            WFElt(i,j,iWFA)=0D0
	          ELSE
	            xg=Sx*Sx
	            Sg=SQRT(xg+yg)
	            S1=COS(Sg)
	            IF (S1.GE.0d0) THEN
	              S=SQRT(Flux*S1)**CosPower
	              WFElt(i,j,iWFA)=DCMPLX(S,0D0)
	              IF (ifPol) THEN
	                CALL DEQUATE(zray,RayDir(1,irayid),3)
	                CALL DXPROD(yray,zray,xGrid)
	                CALL DUNITIZE(yray)
	                CALL DXPROD(xray,yray,zray)
	                DO 32 k=1,3
	                  RayE(k,irayid)=Ex0*DCMPLX((S*xray(k)),0d0)+
     &	                  Ey0*DCMPLX((S*yray(k)),0d0)
   32	                CONTINUE
	              END IF
	            ELSE
	              WFElt(i,j,iWFA)=0D0
	            END IF
	          END IF
	        END IF
   28	      CONTINUE
   29	    CONTINUE
	  ELSE IF (BeamType.EQ.4) THEN
C  Dipole beam...
	    S=SQRT(Flux)
	    IF (ifPol) THEN
	      DO 37 k=1,3
	        RayE(k,1)=S*xGrid(k)
   37	      CONTINUE
	    END IF
	    DO 40 j=1,mdttl
	      DO 39 i=1,mdttl
	        irayid=RayID(i,j)
	        IF (irayid.EQ.0) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          jwf=RayWFMap(1,irayid)
	          thx=DBLE(jwf-(mdttl/2+1))*dxsource
	          S1=COS(thx)
	          IF (S1.LT.0d0) THEN
	            S=0d0
	          ELSE
	            S=Flux*S1
	          END IF
	          WFElt(i,j,iWFA)=DCMPLX(S,0D0)
	          IF (ifPol) THEN
	            CALL DEQUATE(zray,RayDir(1,irayid),3)
	            CALL DXPROD(yray,zray,xGrid)
	            CALL DUNITIZE(yray)
	            CALL DXPROD(xray,yray,zray)
	            DO 38 k=1,3
	              RayE(k,irayid)=DCMPLX((S*xray(k)),0d0)
   38	            CONTINUE
	          END IF
	        END IF
   39	      CONTINUE
   40	    CONTINUE


	  ELSE IF (BeamType.EQ.5) THEN
	    If (.not. allocated(ds1)) Then
	      allocate(ds1(mpts,mpts),ds2(mpts,mpts),
     &               dsc(mdttl,mdttl),DWF_loc(mdttl,mdttl),
     &               stat=m_err_ptsrc)
	      If (m_err_ptsrc /= 0)
     &          call macos_memory_failure
     &              ('sourcsub (BeamType 5): allocate failed!')
	    End If

	    ! Input user defined beam profile ...
	    Open(Unit=12,File=UDSrcProf,Status='old')
	    Do j=1,nGridPts
	      Read(12,*) ds1(1:nGridPts,j)
	    End Do
	    Do j=1,nGridPts
	      Read(12,*) ds2(1:nGridPts,j)
	    End Do
	    Close(12)

!            If (.false.) Then
!              ! Form complex E field of input beam
!              dsc=DCMPLX(0d0,0d0)
!              Do j=1,nGridPts
!                jwf=j+(mdttl-nGridPts+mod(nGridPts,2))/2
!                Do i=1,nGridPts
!                  iwf=i+(mdttl-nGridPts+mod(nGridPts,2))/2
!                  dsc(iwf,jwf)=DCMPLX(ds1(i,j),ds2(i,j))
!                End Do
!              End Do
!
!              ! The complex E field need be 'back-propagated' from
!              ! image plane to pupil plane, since we eventually want
!              ! the input E field to end up at the detector plane.
!              ! -jzlou
!              !
!              ! These system-specific values are from a test run -jzlou
!              du_loc=2.8828222D-02; dx_loc=-6.0386600D-03
!              zElt_loc=1.408516486242871D+02
!              ! Back diffraction propagation
!              Call BKFFPROP(dsc,dx_loc,dsc,du_loc,zElt_loc,WaveBU,
!     &                      mdttl,DWF_loc)
!
!              WFElt(1:mdttl,1:mdttl,iWFA)=dsc
!
!              deallocate(ds1,ds2,dsc,DWF_loc,stat=m_err_ptsrc)
!              If (m_err_ptsrc /= 0)
!     &          call macos_memory_failure
!     &                ('sourcsub: deallocate failed!')
!
!            Else
              ! Assume complex E field is for pupil plane
              WFElt=DCMPLX(0d0,0d0)
              Do j=1,nGridPts
                jwf=j+(mdttl-nGridPts+mod(nGridPts,2))/2
                Do i=1,nGridPts
                  iwf=i+(mdttl-nGridPts+mod(nGridPts,2))/2
                  WFElt(iwf,jwf,iWFA)=DCMPLX(ds1(i,j),ds2(i,j))
                End Do
              End Do
!            End If


	  ELSE
C  Uniform beam...
	    S=SQRT(Flux/DBLE(nRay-1))  ! flux averaged to each ray -jzlou
	    SEx0=S*Ex0
	    SEy0=S*Ey0
	    SC=DCMPLX(S,0D0)
	    IF (ifPol) THEN
	      DO 36 k=1,3
	        RayE(k,1)=SEx0*DCMPLX(xGrid(k),0d0)+
     &	        SEy0*DCMPLX(yGrid(k),0d0)
   36	      CONTINUE
	    END IF
	    DO 25 j=1,mdttl
	      DO 24 i=1,mdttl
	        irayid=RayID(i,j)
	        IF (irayid.EQ.0) THEN
	          WFElt(i,j,iWFA)=0D0
	        ELSE
	          WFElt(i,j,iWFA)=SC
	          IF (ifPol) THEN
	            CALL DEQUATE(zray,RayDir(1,irayid),3)
	            CALL DXPROD(yray,zray,xGrid)
	            CALL DUNITIZE(yray)
	            CALL DXPROD(xray,yray,zray)
	            DO 33 k=1,3
	              RayE(k,irayid)=SEx0*DCMPLX(xray(k),0d0)+
     &	              SEy0*DCMPLX(yray(k),0d0)
	              S1=RayE(k,irayid)
   33	            CONTINUE
	          END IF
	        END IF
   24	      CONTINUE
   25	    CONTINUE
	  END IF
	END IF  ! ifPropagate

	RETURN
	END SUBROUTINE PtSource

C***********************************************************************

C HSEG is for hex segment
C 'S' argument below is SegCoord -jzlou

	FUNCTION HSEG(x,y,width,gap,S,SegXgrid)
	IMPLICIT NONE

	LOGICAL HSEG
	INTEGER S(3)
	REAL*8 x,y,width,gap,w2,L2,xBctr,dxB,slopex,bL,bLctr,
     &	       bR,bRctr,db,SegXgrid(3),xt,yt

C  **** rotate segment frame of reference
c	IF (SegXgrid(1) .EQ. 0) THEN
c	  xt=y
c	  yt=-x
c	ELSE
	  xt=x*SegXgrid(1)+y*SegXgrid(2)
	  yt=-x*SegXgrid(2)+y*SegXgrid(1)
c	END IF
C  ***  end rotate segemnt frame of reference

! for whatever reason, this part of calculation is artificially made
! more confusing, but the result is correct nevertheless. -jzlou, 10/2009
	w2=width/2d0
	L2=w2*1.154700538d0  ! 1.154700538d0 = 2d0/sqrt(3d0)
	xBctr=DBLE(S(1))*w2  ! X coord of seg ctr -jzlou
	dxB=w2-gap/2d0
	slopex=xt*0.5773502692d0  ! 0.5773502692d0 = 1d0/sqrt(3d0)
	bLctr=DBLE(S(2))*L2  ! L coord of seg ctr (actually increased by a factor 2d0/sqrt(3d0)) -jzlou
	bL=yt+slopex         ! L coord of ray at (xt,yt) (increased by a factor 2d0/sqrt(3d0))
	bRctr=DBLE(S(3))*L2  ! R coord of seg ctr -jzlou (increased by a factor 2d0/sqrt(3d0))
	bR=yt-slopex
	db=L2-gap*0.5773502692d0  ! = L2 - gap/2d0*(2/sqrt(3d0))
	IF ((xt.GE.(xBctr-dxB)).AND.(xt.LT.(xBctr+dxB)).AND.
     &	  (bL.GE.(bLctr-db)).AND.(bL.LT.(bLctr+db)).AND.
     &	  (bR.GE.(bRctr-db)).AND.(bR.LT.(bRctr+db))) THEN
	  HSEG=.TRUE.
	ELSE
	  HSEG=.FALSE.
	END IF

	RETURN
	END FUNCTION HSEG

C***********************************************************************

! Pie segment
	FUNCTION PSEG(x,y,width,gap,S,SegXgrid)
	IMPLICIT NONE

	LOGICAL PSEG
	INTEGER S(3)
	REAL*8 x,y,width,gap,w2,L2,xBctr,dxB,slopex,bL,bLctr,
     &	       bR,bRctr,db,SegXgrid(3),xt,yt

C  **** rotate segment frame of reference
c	IF (SegXgrid(1) .EQ. 0) THEN
c	  xt=y
c	  yt=-x
c	ELSE
	  xt=x*SegXgrid(1)+y*SegXgrid(2)
	  yt=-x*SegXgrid(2)+y*SegXgrid(1)
c	END IF
C  ***  end rotate segemnt frame of reference

	w2=width/2d0
	L2=w2*1.154700538d0
	xBctr=DBLE(S(2)-S(3))*w2
	dxB=w2-gap/2d0
	slopex=x*0.5773502692d0
	bLctr=DBLE(S(2))*L2
	bL=y+slopex
	bRctr=DBLE(S(3))*L2
	bR=y-slopex
	db=L2-gap*0.5773502692d0
	PSEG=.TRUE.
	IF (S(1).GE.0) PSEG=PSEG.AND.(x.GE.(xBctr-dxB))
	IF (S(1).LE.0) PSEG=PSEG.AND.(x.LT.(xBctr+dxB))
	IF (S(2).GE.0) PSEG=PSEG.AND.(bL.GE.(bLctr-db))
	IF (S(2).LE.0) PSEG=PSEG.AND.(bL.LT.(bLctr+db))
	IF (S(3).GE.0) PSEG=PSEG.AND.(bR.GE.(bRctr-db))
	IF (S(3).LE.0) PSEG=PSEG.AND.(bR.LT.(bRctr+db))

	RETURN
	END FUNCTION PSEG

C***********************************************************************
! radCtr=   Radius of the central region
! nRings=   Number of segment rings not counting the central region
! nPetals=   Number of segments in each ring. The first number is repeated for some reason...
! RingRads=   Ratio of ring radii to central region ring radius radCtr
!   nSeg=   Number of segments (1 for central region plus number in first and second ring)
!    gap=   Width of gap between segments
!

        PURE LOGICAL FUNCTION FSEG(x,y,r,radCtr,dRad,nPetal,gap,S,SegXgrid)
          use src_mod, only : FSEG_CCAT_Flg
          USE Constants, ONLY: TWOPI, EPS
          IMPLICIT NONE
          REAL*8,  INTENT(IN) :: x,           ! ray: x-pos in Glb. Source CSYS
     &                           y,           ! ray: y-pos in Glb. Source CSYS
     &                           r,           ! ray: radius from centre
     &                           radCtr,      ! Center Segment Radius
     &                           dRad(2)      ! (min, max) radius of ring (gap not included) scaled by radCtr
                                              !  min = RingRads(iRing-1)*radCtr
                                              !  max = RingRads(iRing)*radCtr
          INTEGER, INTENT(IN) :: nPetal       ! Number of segments in each ring. (was "number of petals in flower")
          REAL*8,  INTENT(IN) :: gap          ! Width of gap between rings and petals
          INTEGER, INTENT(IN) :: S(3)         ! SegCoord(1:3,mElt)
                                              ! S(1) == 0: central region (iRing=SegCoord(1,iSeg))
                                              ! S(2) == iPetal for Ring
          REAL*8,  INTENT(IN) :: SegXgrid(3)  ! Source x-axis

C In the flower pattern, prescription is defined by:
C   width=   8d0		% nPetal
C   gap=   1.61679d0		% radCtr

C Order of petals is fixed:
C   1. Center segment
C   2. First petal in SegXgrid coords (rays occupy wedge with
C      theta between 0 and 2*pi/nPetal)
C   3. Next petal.
C   4-nPetals. Etc.

C S(1) = ring ID, S(2) = petal id on that ring

C  radCtr = radius of the center segment. Actual segment boundaries
C           are curcumscribed polygon tangent to circle at outer
C           segment midpoints
C  nPetals = number of petals in flower

	! LOGICAL FSEG
	INTEGER :: iRing, ith
	REAL*8 xt,yt,cth,sth,th,
     &	       dnp,Q,thc,rv,hspan,dd,radmin,radmax,S2,dth

C  **** rotate segment frame of reference
c	IF (SegXgrid(1) .EQ. 0) THEN
c	  xt=y
c	  yt=-x
c	ELSE
	  xt=x*SegXgrid(1)+y*SegXgrid(2)
	  yt=-x*SegXgrid(2)+y*SegXgrid(1)
c	END IF
C  ***  end rotate segment frame of reference

        FSEG = .FALSE.    ! default
        ! ----------------------------------------- Central region (no sub-division)
        IF (S(1).EQ.0) THEN
          ! for central segment
          IF (r.LE.(radCtr-gap/2d0)) THEN
            FSEG=.TRUE.  ! within central circular region
          ELSE
            ! since radCtr is radius of max inscribed circle, a ray outside the circle
            ! could still belong to the central segment. -jzlou
            dnp = nPetal
            Q   = TWOPI/dnp              ! wedge angle value
            cth = xt/r                   ! cos(theta) of data point
            sth = yt/r                   ! sin(theta) of data point
            th  = ATAN2(sth,cth)         ! theta
            IF (th.LT.0d0) th=th+TWOPI   ! ensure 0 <= theta < 2*PI
            ith = INT(th/Q)+1
            thc = (DBLE(ith)-0.5d0)*Q
            rv  = xt*COS(thc)+yt*SIN(thc)
            !IF (rv.LE.(radCtr-gap/2d0)) THEN
            FSEG = (r.LE.(radCtr-gap/2d0))
          END IF

        ! ----------------------------------------- Outer Rings
        ELSE IF (r.GT.EPS) THEN
          ! for other segment
          iRing= S(1)  ! get ring id
          !radmin=radCtr+dRad*DBLE(iRing-1); radmax=radmin+dRad
          radmin=dRad(1); radmax=dRad(2)
          dnp = nPetal
          Q   = twopi/dnp
          cth = xt/r
          sth = yt/r
          th  = ATAN2(sth,cth)
          IF (th.LT.0d0) th=th+twopi
          ith = th/Q+1d0  ! in 'ith' petal
          if (th>(twopi-0.5d0*Q)) th=th-twopi  ! now 'th' for CCAT case only
          dth = th/Q
          S2  = S(2)-1d0
          ! This section was updated by jzlou to accommodate both CCAT and non-CCAT
          ! cases of FSEG segmentation, 08/2020
          thc = -1d18
          IF (.not.FSEG_CCAT_Flg .AND.
     &        ith.EQ.S(2)) THEN
!    &        (dth.GE.(S2) .AND. dth.LE.(S2+1))) THEN
            thc=(DBLE(ith)-0.5d0)*Q
          END IF
          IF (FSEG_CCAT_Flg .AND.
     &        (dth.GE.(S2-0.5d0) .AND. dth.LE.(S2+0.5d0))) THEN
            thc=S2*Q
          END IF
          !
          IF (thc > -1d10) THEN
              rv=xt*COS(thc)+yt*SIN(thc)
             !IF ((rv.GT.(radmin+gap/2d0)).AND.(rv.LT.(radmax-gap/2d0)))THEN
              IF ((r.GT.(radmin+gap/2d0)).AND.(r.LT.(radmax-gap/2d0)))THEN
                hspan = rv*TAN(0.5d0*Q)
                dd    = SQRT(abs(r*r-rv*rv))
              IF (dd.LE.(hspan-gap/2d0)) THEN
                FSEG=.TRUE.
              ELSE
                FSEG=.FALSE.
              END IF
            END IF
          END IF
        END IF

        END FUNCTION FSEG

C *************************************************************************

C**************************************************************************
!       ReGrid is called after near-field propagations, when ifRegrid(iElt)
!       is set to TRUE (default to false, use Regrid command to reset)
!       -jzlou
C**************************************************************************

	SUBROUTINE ReGrid(iElt)
	use param_mod
	! use src_mod, only : nGridPts   ! Resample OPDMat to new grid
        use elt_mod

	IMPLICIT NONE

	INTEGER i,iElt,iEm1,iWF,jWF,iCR,jCR,iRayCR,iRay

	REAL*8 th(3),CRDir(3),CRPos(3),RotPt(3),delx,dely,ddel,
     &	       xGrid(3),yGrid(3),zGrid(3),D1(9),D2(9),dth,thx,thy

C	Compute "chief ray" indices - these coincide with the real
C	chief ray only if an odd number of grid pts is used
C
	iRayCR=nRay/2+1
	iCR=RayWFMap(1,iRayCR)
	jCR=RayWFMap(2,iRayCR)

C	Compute rotation point, chief ray position and direction,
C	and local grid coordinates

	DO i=1,3
	  CRDir(i)=RayDir(i,1)
	  CRPos(i)=RayPos(i,1)
	  RotPt(i)=CRPos(i)+zElt(iElt)*CRDir(i)
	  D1(i)=RayPos(i,RayID((iCR+1),jCR))-
     &	        RayPos(i,RayID((iCR-1),jCR))  ! candidate for xGrid
	  D2(i)=RayDir(i,RayID((iCR-1),jCR))-
     &	        RayDir(i,RayID((iCR+1),jCR))  ! candidate for xGrid
    	END DO

	CALL DNEGATE(zGrid,RayDir(1,1),3)
	CALL DUNITIZE(zGrid)
	IF (DDOTC(D1,D1).GT.DDOTC(D2,D2)) THEN
	  CALL DXPROD(yGrid,zGrid,D1)
	ELSE
	  CALL DXPROD(yGrid,zGrid,D2)
	END IF
	CALL DUNITIZE(yGrid)
	CALL DXPROD(xGrid,yGrid,zGrid)

C	Compute new ray grid

	iEm1=iElt-1
	IF ((PropType(iEm1).EQ.2).OR.
     &	  (PropType(iEm1).EQ.5).OR.(PropType(iEm1).EQ.7).OR.
     &	  (PropType(iEm1).EQ.8).OR.(PropType(iEm1).EQ.11)) THEN
	  ! 2='NFS1surf', 5='NFSpherical', 7='SF1surf',
	  ! 8='SpatialFilter', 11='NF2'
C Regrid on a sphere
	  dth=ABS(dxElt(iElt)/zElt(iElt))
	  DO iRay=2,nRay
	    iWF=RayWFMap(1,iRay)
	    jWF=RayWFMap(2,iRay)
	    thx=DBLE(iWF-iCR)*dth
	    thy=DBLE(jWF-jCR)*dth
	    DO i=1,3
	      th(i)=-xGrid(i)*thy+yGrid(i)*thx
    	    END DO
	    CALL ReRay(th,CRDir,CRPos,RotPt,RayDir(1,iRay),
     &	    	       RayPos(1,iRay))
    	  END DO
	ELSE IF ((PropType(iEm1).EQ.3).OR.
     &	  (PropType(iEm1).EQ.4).OR.(PropType(iEm1).EQ.6).OR.
     &	  (PropType(iEm1).EQ.10)) THEN
	  ! 3='FarField', 4='NFP1surf', 6='NFPlane', 10='NF1'
C Regrid on a plane
	  ddel=ABS(dxElt(iElt))

#if 0
	  ! Resample OPDMat to new grid
	  i=1
	  Do iRay=1,nRay
	    oldRayPos(i:i+1)=RayPos(1:2,iRay)
	    i=i+2
	  End Do
	  i=1
	  Do j=1,nGridPts
            oldOPDMat(i:i+nGridPts-1)=OPDMat(1:nGridPts,j)
	    i=i+nGridPts
	  End Do
	  x0=CRPos(1)-ddel*nGridPts/2
	  y0=CRPos(2)-ddel*nGridPts/2
	  Call mOPDGridData(oldRayPos,nGridPts,x0,y0,ddel,oldOPDMat,
     &                      OPDMat)
#endif

	  DO iRay=2,nRay
	    iWF=RayWFMap(1,iRay)
	    jWF=RayWFMap(2,iRay)
	    delx=DBLE(iWF-iCR)*ddel
	    dely=DBLE(jWF-jCR)*ddel
	    DO i=1,3
	      RayDir(i,iRay)=CRDir(i)
	      RayPos(i,iRay)=CRPos(i)+xGrid(i)*delx+yGrid(i)*dely
    	    END DO
    	  END DO
	ELSE
	  WRITE(*,*)' Warning: Rays not regridded -- improper PropType'
	END IF

	RETURN
	END SUBROUTINE ReGrid

C***********************************************************************

	SUBROUTINE ReRay(th,CRDir,CRPos,RotPt,RayDir,RayPos)
	IMPLICIT NONE
	REAL*8 th(3),CRDir(3),CRPos(3),RotPt(3),Q(3,3),dQ(3,3),
     &	       D1(3),D2(3),RayDir(3),RayPos(3)

C       Q and dQ computed in Qform, and dQ = Q - I, -jzlou
	CALL Qform(Q,dQ,th)

C       CRDir = ChfRayDir, and from the routine
C       we get RayDir = Q*CRDir -jzlou
	CALL DMPROD(RayDir,Q,CRDir,3,3,1)

	D1(1:3) = CRPos(1:3) - RotPt(1:3)
	CALL DMPROD(D2,dQ,D1,3,3,1)
	RayPos(1:3) = CRPos(1:3)+D2(1:3) ! effectively RayPos(1:3)=Q*D1
	RETURN
	END SUBROUTINE ReRay

C***********************************************************************

! This routine is not called by 'btc' - jzlou

	SUBROUTINE AddAtmos(iElt,npts,iWFA,R1,CD2,r0,atmWavelen,pTilt,
     &	  NoiseSeed,ifPol,ifGrid,ifPropagate,Dicr,Djcr,
     &	  dxSource,xhat,yhat,zhat)
	use param_mod
        use src_mod
        use elt_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifGrid,ifPropagate,ifOPDCalc,
     &	 ifRegrid
	INTEGER i,j,iElt,npts,iWFA,iray,k,NoiseSeed(2),
     &	 im1,ip1,jm1,jp1,iwf,jwf
	SREAL R1(mdttl,mdttl),r0,rsamp,pTilt,atmWavelen
	 REAL*8 Dicr,Djcr,dxSource,xhat(3),yhat(3),zhat(3),ztdir,
     &	 dx,C1,dxMin,dxMax,dxAvg,dxDev,tiltx,tilty,
     &	 Rij,Rip1j,Rim1j,Rijp1,Rijm1,Rip1jp1,Rim1jp1,Rip1jm1,Rim1jm1
	COMPLEX*8 CD2(mdttl,mdttl)
	COMPLEX*16 cztdir

	ztdir = 0.0d0  ! -jzlou

C  Generate ray grid if not previously done
	IF ((iElt.EQ.0).AND..NOT.ifGrid) THEN
	  IF (ABS(zSource).GE.1d10) THEN
	    CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	    ifPropagate)
	  ELSE
	    CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &	    ifPropagate)
	  END IF
	  ifGrid=.TRUE.
	  dx=dxSource
	ELSE IF (iElt.EQ.0) THEN
	  dx=dxSource
	ELSE
	  ifOPDCalc=.FALSE.
	  ifRegrid=.FALSE.
	  CALL DXCALC(iElt,iElt,dxMin,dxMax,
     &	  dxAvg,dxDev,ifRegrid,ifOPDCalc,npts)
	  dx=dxElt(iElt)
	END IF

C  Generate phase screen: note that this should be at a pupil!
	nGridPts=npts+1
	rsamp=dx
	CALL AtmosPhase(CD2,R1,rsamp,r0,atmWavelen,nGridPts,pTilt,
     &	                NoiseSeed,mdttl)
c-- re-subtract piston start
c	AvgR=0e0
c	DO iRay=2,nRay
c	  i=RayIndex(1,iRay)
c	  j=RayIndex(2,iRay)
c	  AvgR=AvgR+R1(i,j)
c	END DO
c	AvgR=AvgR/REAL(nRay-1)
c	DO iRay=2,nRay
c	  i=RayIndex(1,iRay)
c	  j=RayIndex(2,iRay)
c	  R1(i,j)=R1(i,j)-AvgR
c	END DO
c-- re-subtract piston end

C  Compute ray tip and tilt (beamwalk = 0)

c	S1=DDOTC(zHat,RayDir(1,1))
c	IF (S1.GT.0d0) THEN
c	  ztdir=1d0
c	ELSE
c	  ztdir=-1d0
c	END IF
c	C1=ztdir*0.5d0/dx

	C1=0.5d0/dx

	DO iRay=2,nRay
	  iwf=RayWFMap(1,iRay)
	  jwf=RayWFMap(2,iRay)
	  im1=iwf-1
	  ip1=iwf+1
	  jm1=jwf-1
	  jp1=jwf+1
	  Rij=R1(iwf,jwf)
	  Rim1j=R1(im1,jwf)
	  Rip1j=R1(ip1,jwf)
	  Rijm1=R1(iwf,jm1)
	  Rijp1=R1(iwf,jp1)
	  Rip1jp1=R1(ip1,jp1)
	  Rim1jp1=R1(im1,jp1)
	  Rip1jm1=R1(ip1,jm1)
	  Rim1jm1=R1(im1,jm1)

C	  tiltx=C1*(Rim1j-Rip1j)
C	  tilty=C1*(Rijm1-Rijp1)
	  tiltx=C1*(Rim1j-Rip1j+Rim1jm1+Rim1jp1-Rip1jm1-Rip1jp1)
	  tilty=C1*(Rijm1-Rijp1+Rim1jm1+Rip1jm1-Rim1jp1-Rip1jp1)

C  Add atmospheric effects to ray states
c	  CumRayL(iRay)=CumRayL(iRay)-ztdir*Rij

	  CumRayL(iRay)=CumRayL(iRay)-Rij

c	  CALL Qform(Q,dQ,th)
c	  CALL DMPROD(D1,Q,RayDir(1,iRay),3,3,1)
c	  CALL DEQUATE(RayDir(1,iRay),D1,3)
	  DO k=1,3
	    RayDir(k,iRay)=RayDir(k,iRay)+tiltx*xhat(k)+tilty*yhat(k)
	  END DO
	  CALL DUNITIZE(RayDir(1,iRay))
	END DO

C  Update complex amplitude matrix if propagating
	cztdir=-(0d0,1d0)*ztdir
	IF (ifPropagate) THEN
	  DO j=1,mdttl
	    DO i=1,mdttl
	      iRay=RayID(i,j)
	      IF ((iRay.LE.nRay).AND.(iRay.GT.0).AND.LRayOK(iRay))
     &	        WFElt(i,j,iWFA)=WFElt(i,j,iWFA)*CDEXP(cztdir*R1(i,j))
	    END DO
	  END DO
	END IF

	RETURN
	END SUBROUTINE AddAtmos

C***********************************************************************
C***********************************************************************
C	AtmosPhase: subroutine for wavefront simulation
C	From routine "screenle" by Mark Colavita
C	Modified for COMP by Dave Redding 1/16/94
C
C	Subroutine AtmosPhase(x,rsamp,r0,nradius,itilt):
C	Generates a turbulent OPD array with piston errors
C	subtracted and tilt added back in (short exposure)
C	mmc 2/89
C
C	      x = npts*npts COMPLEX*8 array (used internally)
C	      f = npts*npts SREAL array containing the returned phase
C	     r0 = coherence diameter in cm
C	  rsamp = sampling interval in cm
C	nradius = aperture radius, i.e., radius in samples over which to
C		  subtract piston and tilt; nradius should be <= npts/4
C	  pTilt = percent of tilt added back in (ranges from 0 to 1)
C
C	This residual screen should have a phase variance within
C	a circle of radius nradius equal to 0.134(D/r0)^5/3,
C	D=2*nradius*rsamp.
C
C	Test results with oversample = 4 (SE=screen, LE=screenle),
C	100 Monte Carlo iterations
C	npts = 256, D/r0 = 9.45, SE var = 5.81 (th = 5.66)
C                                LE var = 55.2 (th = 42.0)
C	npts = 512, D/r0 = 18.9, SE var = 18.4 (th = 18.0)
C                                LE var = 119 (th = 134)
C	npts = 1024, D/r0 = 37.8, SE var = 54.5 (th = 57.1)
C                                LE var = 454 (th = 425)

	SUBROUTINE AtmosPhase(x,f,rsamp,r0,atmWavelen,
     &	                      nradius,pTilt,NoiseSeed,npts)
	IMPLICIT SREAL (A-H,O-Z)
	INTEGER npts,nn(2),noiseSeed(2)
	INTEGER*4 ii,jj,ii2,jj2
	COMPLEX*8 x(npts,npts),cs2o2
	SREAL f(npts,npts),sqrt2,s2o2,atmWavelen,Crad2d

	SAVE

	DATA rlast/-999d0/
	DATA s2o2/0.70710678118655d0/
	DATA cs2o2/(0.0d0,0.70710678118655d0)/
	DATA sqrt2/1.41421356237310d0/
	DATA Crad2d/0.15915494309190d0/

	r00=ABS(r0)
	nn(1)=npts
	nn(2)=npts

C  Initialization section
	IF (rlast.NE.r00) THEN
	  rlast=r00
	END IF

C  Generate filter for kolmogorov turbulence:
C  sqrt(psd)={0.151/r0**(5/6)}f**(-11/6)

	power=0.91666666666667d0
	coeff=(0.151d0/r00**0.83333333333333d0)*(1d0/(DBLE(npts)*rsamp))
	scale=(1d0/(DBLE(npts)*rsamp))**2
	rad2d=Crad2d*atmWavelen

	DO j=1,npts
	  DO i=1,npts
	    x(i,j)=0d0
	    f(i,j)=0d0
	  END DO
	END DO
	DO 50 i=1,npts
	  ii=i-1
	  IF (ii.GE.npts/2)ii=ii-npts
	  ii2=ii*ii
	  DO 50 j=1,npts
	    jj=j-1
	    IF (jj.GE.npts/2)jj=jj-npts
	    jj2=jj*jj
	    ff=scale*(ii2+jj2)
	    IF (ff.EQ.0.0)ff=scale
	    f(j,i)=coeff/ff**power
50	CONTINUE

C  Generate noise array with proper symmetry
C	  (using .707gas+j.707gas for all pts will yield two indep phase
C	  screens--one each in the REAL and imag parts of the transform.
C	  Using a symmetric array requires half as many random deviates
C	  and produces a single screen in the REAL part of the
C	  transform.)
	n2=npts/2
	n21=npts/2+1
	n22=npts/2+2
	np2=npts+2
	DO 100 i=1,n21
	  DO 100 j=1,n21
	    x(j,i)=s2o2*GAUSS(NoiseSeed)+cs2o2*GAUSS(NoiseSeed)
100	CONTINUE
	DO 110 j=n22,npts
	  x(j,1)=CONJG(x(np2-j,1))
	  x(j,n21)=CONJG(x(np2-j,n21))
110	CONTINUE
	DO 120 i=2,n2
	  DO 120 j=n22,npts
	    x(j,i)=s2o2*GAUSS(NoiseSeed)+cs2o2*GAUSS(NoiseSeed)
120	CONTINUE
	DO 130 i=n22,npts
	  np2i=np2-i
	  x(1,i)=CONJG(x(1,np2i))
	  DO 130 j=2,npts
	    x(j,i)=CONJG(x(np2-j,np2i))
130	CONTINUE
	x(1,1)=sqrt2*REAL(x(1,1))
	x(1,n21)=sqrt2*REAL(x(1,n21))
	x(n21,1)=sqrt2*REAL(x(n21,1))
	x(n21,n21)=sqrt2*REAL(x(n21,n21))

C    multiply by filter
	DO 200 i=1,npts
	  DO 200 j=1,npts
	    x(j,i)=x(j,i)*f(j,i)
200	CONTINUE

C    inverse transform
	CALL RFOURN(x,nn,2,-1)

C    calc piston and tilt
	IF ((nradius.GT.0).AND.(nradius.LE.npts)) THEN

	  bm=0d0
	  bx=0d0
	  by=0d0
	  cm=0d0
	  cx=0d0
	  cy=0d0
	  nr2=nradius*nradius
	  DO 300 i=-nradius,nradius
	    IF (i.LT.0) THEN
	      i1=i+npts+1
	    ELSE
	      i1=i+1
	    END IF
	  DO 300 j=-nradius,nradius
	    IF ((i*i+j*j).LE.nr2) THEN
	      IF (j.LT.0) THEN
	        j1=j+npts+1
	      ELSE
	        j1=j+1
	      END IF
	      foo=REAL(x(j1,i1))
	      bm=bm+foo
	      cm=cm+1
	      bx=bx+foo*i
	      cx=cx+i*i
	      by=by+foo*j
	      cy=cy+j*j
	    END IF
300	  CONTINUE
	  bm=bm/cm
	  bx=bx/cx
	  by=by/cy

C    IF selected, add back in approx. tilt
	  IF (pTilt .NE. 0.0) THEN
	    tiltamp = pTilt
     &	      *SQRT(1.8d0*(2d0*nradius*rsamp/r00)**1.66666666667)
	    bxt = tiltamp * GAUSS(NoiseSeed) / nradius
	    bx = bx + bxt
	    byt = tiltamp * GAUSS(NoiseSeed) / nradius
	    by = by + byt
	  END IF

C    subtract piston and tilt:
	  DO 400 i=-npts/2, npts/2-1
	    IF (i.LT.0) THEN
	      i1=i+npts+1
            ELSE
	      i1=i+1
	    END IF
	  DO 400 j=-npts/2, npts/2-1
	    IF (j.LT.0) THEN
	      j1=j+npts+1
	    ELSE
	      j1=j+1
	    END IF
	    x(j1,i1)=REAL(x(j1,i1))-bm-bx*i-by*j
400	  CONTINUE

	END IF

C  swap to put the center at pix npts/2+1
	CALL RSWAP3(x,npts)

C  copy to real array
	DO i=1,npts
	  DO j=1,npts
	    f(j,i)=rad2d*REAL(x(j,i))
	  END DO
	END DO

	RETURN
	END SUBROUTINE AtmosPhase

C***********************************************************************

#if 0
        SUBROUTINE BKFFPROP(WF2,DX2,WF1,DX1,DZ,LAMBDA,NFFT,DWF)
          use Constants, only: PI
          use dft_mod, only: REVAPPLYFAC2
          IMPLICIT NONE
          INTEGER I,J,NN(2),NFFT
          REAL*8 DX1,DX2,DZ,S,LAMBDA,C2
          COMPLEX*16 WF1(NFFT,NFFT),WF2(NFFT,NFFT),DWF(NFFT,NFFT),C1

        NN(1)=NFFT
        NN(2)=NFFT
        S=1D0/DBLE(NFFT)
        C1=DCMPLX(S)

        CALL DSWAP2(DWF,WF1,NFFT)
        CALL DFOURN(DWF,NN,2,-1)
        CALL DSWAP2(WF2,DWF,NFFT)

        DO J=1,NFFT
          DO I=1,NFFT
            !WF2(I,J)=C1*WF2(I,J)
            WF2(I,J)=WF2(I,J)  ! updated by jzlou, when DFOURN is NR based
          END DO
        END DO

        C1 = 1d0/DCMPLX(0d0,lambda*dz)
        C2 = PI/lambda/dz
        CALL RevApplyfac2(wf2,NFFT,NFFT,C1*dx1*dx1,C2,dx2)

        END SUBROUTINE BKFFPROP
#endif

C***********************************************************************
#include "ssrcray.inc"
C***********************************************************************
	End Module sourcsub_mod
