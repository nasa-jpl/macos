C
C ssrcset.inc
C

 	SUBROUTINE SetColSrcRay(npts,Dicr,Djcr,dxSource,
     &                          xo,yo,ifPol,ifPropagate)
	use param_mod
        use src_mod
        use elt_mod
	use traceutil_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifPropagate,ifRay
	INTEGER npts,iRay,i,j,k,iSeg,iWFA,iwf,jwf,irayid,nDrawRay,
     &          nGridPts_loc,iRing
	REAL*8 A,O2,A2,ds,x,y,r2,Tbg(3,2),dpr(3),ddr(3),
     &	       Dicr,Djcr,dxSource,w2,ciy,rx2,ry2,xo,yo,
     &	       Sy,yg,Sx,xg,Sg,S,S1,SegX2(3),D1(9),dRad(2)
	COMPLEX*16 SC,SEx0,SEy0
	Real*8, parameter :: tan30=5.773502691896257d-01,
     &			     pio2=1.57079632679490

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

! From "prop_defs.inc"
! GridTypeName(1)='Circular'
! GridTypeName(2)='Square'
! GridTypeName(3)='Hex'
! GridTypeName(4)='Pie'
! GridTypeName(5)='Flower'
 
	IF ((GridType.GE.3).AND.(GridType.LE.5)) THEN
	  ! Hex, Pie or Flower source grid type
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
	nGridPts_loc=npts+1  ! no. of grid pts in X (or Y) of source aperture
	iRay=1  ! chief ray

	IF (Obscratn.NE.0d0) THEN
	  O2=(Obscratn/2d0)**2
	ELSE
	  O2=0d0
	END IF

	A2=(Aperture/2d0)**2  ! original
	!A2= (Aperture/2d0 + Aperture/DBLE(npts))**2  ! for plotting a 'bottom' ray

! nRayToSeg(iRay) - given iRay, (nRayToSeg(iRay)-1) gives the number of
! rays that have hit iSeg before iRay does. -jzlou

        ! GridType from 3 to 5 are segmented 
        ! - jzlou
	IF (GridType.EQ.3) THEN
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
	ELSE IF (GridType.EQ.4) THEN
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
	ELSE IF (GridType.EQ.5) THEN
	  DO iSeg=1,nSeg
	    iRing=SegCoord(1,iSeg)
	    If (iRing>0) Then
              dRad(1)=RingRads(iRing-1)*radCtr                  
              dRad(2)=RingRads(iRing)*radCtr                    
            End If
	    IF (FSEG(0d0,0d0,0d0,radCtr,dRad,nPetals(iRing),gap,
     &	        SegCoord(1,iSeg),SegX2)) THEN
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

	SegMap(1:mpts,1:mpts)=0; RayID(1:mpts,1:mpts)=0

        ! In the code section below, jwf,iwf are indexes to put source ray grid 
        ! (nGridPts_loc x nGridPts_loc) at the center of WF grid (mdttl x mdttl). -jzlou

	IF (nGridPts_loc.GT.1) THEN
	  A=Aperture/2d0
	  ds=Aperture/DBLE(npts)
	  y=yo*(Aperture/2d0)
	  x=xo*(Aperture/2d0)
	  j=2; i=1; jwf=j; iwf=i 

          IF (GridType.EQ.1) THEN
	    r2=y**2+x**2
	    IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	      iRay=iRay+1
	      RayIndex(1,iRay)=i
	      RayIndex(2,iRay)=j
	      RayMap(i,j)=iRay
	      RayID(iwf,jwf)=iRay
	      RayWFMap(1,iRay)=iwf  ! index for WF grid array 'WFElt'
	      RayWFMap(2,iRay)=jwf
	      ddr(1)=x; ddr(2)=y
	      CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	      RayDir(1:3,iRay)=ChfRayDir(1:3)  ! collimated source
	      RayPos(1:3,iRay)=ChfRayPos(1:3)+dpr(1:3)
	    END IF
	  ELSE IF (GridType.EQ.2) THEN
	    IF (((ABS(x).LE.A).AND.(ABS(y).LE.A))
     &	       .AND.(iRay.LT.mRay)) THEN
	      iRay=iRay+1
	      RayIndex(1,iRay)=i
	      RayIndex(2,iRay)=j
	      RayMap(i,j)=iRay
	      RayID(iwf,jwf)=iRay
	      RayWFMap(1,iRay)=iwf
	      RayWFMap(2,iRay)=jwf
	      ddr(1)=x
	      ddr(2)=y
	      CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	      DO 12 k=1,3
	        RayDir(k,iRay)=ChfRayDir(k)
	        RayPos(k,iRay)=ChfRayPos(k)+dpr(k)
   12	      CONTINUE
	    END IF
	  ELSE IF (GridType.EQ.3) THEN
	    ifRay=.FALSE.
	    r2=y**2+x**2
	    IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	      DO 7 iSeg=1,nSeg
	        IF (HSEG(x,y,width,gap,SegCoord(1,iSeg),SegX2)) THEN
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
	                ddr(1)=x
	                ddr(2)=y
	                CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	                DO k=1,3
	                  RayDir(k,iRay)=ChfRayDir(k)
	                  RayPos(k,iRay)=ChfRayPos(k)+dpr(k)
 	                END DO
	          END IF
	          nRayToSeg(iRay)=nRayToSeg(iRay)+1
	          RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	        END IF
  7	      CONTINUE
	    END IF

	  ELSE IF (GridType.EQ.4) THEN
	    ifRay=.FALSE.
	    r2=y**2+x**2
	    IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	      DO 8 iSeg=1,nSeg
	        IF (PSEG(x,y,width,gap,SegCoord(1,iSeg),SegX2)) THEN
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
	                ddr(1)=x
	                ddr(2)=y
	                CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	                DO k=1,3
	                  RayDir(k,iRay)=ChfRayDir(k)
	                  RayPos(k,iRay)=ChfRayPos(k)+dpr(k)
 	                END DO
	              END IF
	              nRayToSeg(iRay)=nRayToSeg(iRay)+1
	              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	        END IF
    8	      CONTINUE
	    END IF
	  ELSE IF (GridType.EQ.5) THEN
            ifRay=.FALSE.
	    r2=y**2+x**2
	    IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	          r2=DSQRT(r2)
	          DO iSeg=1,nSeg
		    iRing=SegCoord(1,iSeg)
	            If (iRing>0) Then
              	      dRad(1)=RingRads(iRing-1)*radCtr        
                      dRad(2)=RingRads(iRing)*radCtr        
                    End If
		    IF (FSEG(x,y,r2,radCtr,dRad,nPetals(iRing),gap,
     &                  SegCoord(1,iSeg),SegX2)) THEN
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
	                ddr(1)=x
	                ddr(2)=y
	                CALL DMPROD(dpr,Tbg,ddr,3,2,1)
	                DO k=1,3
	                  RayDir(k,iRay)=ChfRayDir(k)
	                  RayPos(k,iRay)=ChfRayPos(k)+dpr(k)
 	                END DO
	              END IF
	              nRayToSeg(iRay)=nRayToSeg(iRay)+1
	              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	            END IF
   		  END DO
	    END IF
	  END IF
	  nRay=iRay
	ELSE
	  nRay=1
	  RayIndex(1,iRay)=1
	  RayIndex(2,iRay)=1
	  RayMap(1,1)=1
	END IF

	RayIndex(1,1)=0
	RayIndex(2,1)=0

	IF (MOD(nGridPts_loc,2).EQ.0) THEN
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
	RETURN
	END SUBROUTINE SetColSrcRay

C***********************************************************************

	SUBROUTINE SetPtSrcRay(npts,Dicr,Djcr,dxSource,
     &            	       xo,yo,ifPol,ifPropagate)
	use param_mod
        use src_mod
        use elt_mod

	IMPLICIT NONE

	LOGICAL ifPol,ifPropagate,ifRay
	INTEGER npts,iRay,i,j,k,iSeg,iWFA,iwf,jwf,irayid,
     &          nGridPts_loc,iRing
	REAL*8 A,O2,A2,ds,thx,thy,r2,Tbg(3,2),rx2,ry2,
     &	       Dicr,Djcr,dxSource,th(3),RotPt(3),xo,yo,
     &	       Sy,yg,Sx,xg,Sg,S,xray(3),yray(3),zray(3),
     &	       S1,SegX2(3),D1(9),dRad(2)
	COMPLEX*16 SC,SEx0,SEy0
	Real*8, parameter :: pio2=1.57079632679490

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
	  CALL DXPROD(D1,ChfRayDir,SegXgrid)
	  CALL DXPROD(SegXgrid,D1,ChfRayDir)
	  CALL DUNITIZE(SegXgrid)
	  SegX2(1)=DDOTC(xGrid,SegXgrid)
	  SegX2(2)=DDOTC(yGrid,SegXgrid)
	  SegX2(3)=0d0
	END IF

C  Find the segment for the chief ray 
	nGridPts_loc=npts+1
	iRay=1
	IF (Obscratn.NE.0d0) THEN
	  O2=(Obscratn*0.49999d0)**2
	ELSE
	  O2=0d0
	END IF

	A2=(Aperture*0.50001)**2
	!A2= (Aperture/2d0 + Aperture/DBLE(npts))**2  ! for plotting a 'bottom' ray


	IF (GridType.EQ.3) THEN
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
	ELSE IF (GridType.EQ.4) THEN
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
	ELSE IF (GridType.EQ.5) THEN
	  DO iSeg=1,nSeg
	    iRing=SegCoord(1,iSeg)
	    If (iRing>0) Then
              dRad(1)=RingRads(iRing-1)*radCtr        
              dRad(2)=RingRads(iRing)*radCtr        
            End If
	    IF (FSEG(0d0,0d0,0d0,radCtr,dRad,nPetals(iRing),gap,
     &          SegCoord(1,iSeg),SegX2)) THEN
	      nRayToSeg(iRay)=nRayToSeg(iRay)+1
	      RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	    END IF
	  END DO
	  IF (nRayToSeg(iRay).EQ.0) THEN
	    nRayToSeg(iRay)=1
	    RayToSegMap(iRay,nRayToSeg(iRay))=1
	  END IF
	END IF

C       RotPt defines the spatial position of the physical light source,
C       see MACOS manual. -jzlou		
	DO 16 i=1,3
	  RotPt(i)=ChfRayPos(i)+zSource*ChfRayDir(i)
   16	CONTINUE

	th(1:3)=0d0

	CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,RayDir(1,iRay),
     &		   RayPos(1,iRay))

	IF (nGridPts_loc.GT.1) THEN
	  RayMap(1:mpts,1:mpts)=0
	  SegMap(1:mpts,1:mpts)=0
	  RayID(1:mdttl,1:mdttl)=0

	  A=Aperture*0.50001d0
	  ds=Aperture/DBLE(npts)
	  thy=yo*Aperture/2d0
	  thx=xo*Aperture/2d0
	  j=2; i=1; jwf=j; iwf=i

	  IF (GridType.EQ.1) THEN
	    r2=thy**2+thx**2
	    IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	      iRay=iRay+1
	      RayIndex(1,iRay)=i
	      RayIndex(2,iRay)=j
	      RayMap(i,j)=iRay
	      RayID(iwf,jwf)=iRay
	      RayWFMap(1,iRay)=iwf
	      RayWFMap(2,iRay)=jwf

	      ! Eular vector
	      DO 12 k=1,3
	        th(k)=-xGrid(k)*thy+yGrid(k)*thx
 12	      CONTINUE

	      ! this routine sets ray direction and position -jzlou
	      CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,
     &	                 RayDir(1,iRay),RayPos(1,iRay))
	    END IF
	  ELSE IF (GridType.EQ.2) THEN
	    IF (((ABS(thx).LE.A).AND.(ABS(thy).LE.A))
     &	            .AND.(iRay.LT.mRay)) THEN
	          iRay=iRay+1
	          RayIndex(1,iRay)=i
	          RayIndex(2,iRay)=j
	          RayMap(i,j)=iRay
		  RayID(iwf,jwf)=iRay
	          RayWFMap(1,iRay)=iwf
	          RayWFMap(2,iRay)=jwf
	          DO 13 k=1,3
	            th(k)=-xGrid(k)*thy+yGrid(k)*thx
 13	          CONTINUE
	          CALL ReRay(th,ChfRayDir,ChfRayPos,RotPt,
     &	                     RayDir(1,iRay),RayPos(1,iRay))
	    END IF
	  ELSE IF (GridType.EQ.3) THEN
	    ifRay=.FALSE.
	    r2=thy**2+thx**2
	    IF (((r2.GE.O2).AND.(r2.LE.A2)).AND.(iRay.LT.mRay)) THEN
	      DO 7 iSeg=1,nSeg
	          IF (HSEG(thx,thy,width,gap,SegCoord(1,iSeg),
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
     &	                RayDir(1,iRay),RayPos(1,iRay))
	              END IF
	              nRayToSeg(iRay)=nRayToSeg(iRay)+1
	              RayToSegMap(iRay,nRayToSeg(iRay))=iSeg
	            END IF
  7	      CONTINUE
	    END IF
	  ELSE IF (GridType.EQ.4) THEN
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
    8	      CONTINUE
	    END IF
	  END IF
	  nRay=iRay
	ELSE
	  nRay=1
	  RayIndex(1,iRay)=1
	  RayIndex(2,iRay)=1
	  RayMap(1,1)=1
	END IF
	RayIndex(1,1)=0
	RayIndex(2,1)=0

	IF (MOD(nGridPts_loc,2).EQ.0) THEN
	  Dicr=DBLE(mdttl/2)+0.5d0
	  Djcr=DBLE(mdttl/2)+1.5d0
	ELSE
	  Dicr=mdttl/2+1
	  Djcr=Dicr
	END IF
	dxSource=ds

C  Set source amplitude if propagating or for polarized light

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
     &	              Ey0*DCMPLX((Sg*yray(k)),0d0)
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
	END IF
	RETURN
	END SUBROUTINE SetPtSrcRay

