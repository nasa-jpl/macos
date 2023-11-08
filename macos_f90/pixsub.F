C***********************************************************************
C       Begin file propsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	SUBROUTINE CPIXILATE(PixArray,WFArray,XDET,YDET,drx,dry,
     &		     npix,mPix,nWF,dxWF,dxpix,iPloTyp,ifPixLoc,ifAdd)
	use math_mod

C	This routine takes wavefront intensity data contained in WFArray
C	and redistributes it into detector pixels (PixArray).  WF grid
C	size is set by dxWF; grid is nWF by nWF.  Pixel size is dxpix;
C	the detector has npix by npix square pixels.

	IMPLICIT NONE

	LOGICAL ifPixLoc,ifAdd

	INTEGER i,j,iPloTyp,ipix,jpix,npix,nWF,ip1,im1,jp1,jm1,mPix

	SREAL WFArray(nWF,nWF),PixArray(mPix,mPix),xc,yc,xcpix,ycpix,
     &	 XDET(mPix),YDET(mPix),dxWF,dxpix,xij,yij,dxplot,xpix,ypix,
     &	 xWFmin,xpixmin,xpixmax,xWFmax,yWFmin,ypixmin,ypixmax,yWFmax,
     &	 A,A1,A2,A3,A4,A5,A6,A7,A8,A9,IWF,dx1,dx2,dx3,dy1,dy2,dy3,
     &	 drx,dry

C  1. Compute various parameters

	IF (ifPixLoc) THEN
	  ! xc = x_center coord (not pix coord) of WF grid
	  ! yc = y_center coord (not pix coord) of WF grid
	  xc=FLOAT(nWF+1)*dxWF/2e0+drx  
	  yc=FLOAT(nWF+1)*dxWF/2e0+dry
	  xcpix=FLOAT(npix)*dxpix/2e0
	  ycpix=FLOAT(npix)*dxpix/2e0
	ELSE
	  xc=FLOAT(nWF+1)*dxWF/2e0
	  yc=FLOAT(nWF+1)*dxWF/2e0
	  xcpix=FLOAT(npix)*dxpix/2e0
	  ycpix=FLOAT(npix)*dxpix/2e0
	END IF
	A=dxWF*dxWF

C  2. Initialize array

	dxplot=FLOAT(npix)/FLOAT(npix-1)
	IF (.NOT.ifAdd) THEN
	  DO 2 jpix=1,npix
	    DO 1 ipix=1,npix
	      PixArray(ipix,jpix)=0e0
  1	    CONTINUE
  2	  CONTINUE
	END IF

	DO 12 jpix=1,npix
          XDET(jpix)=(jpix-1)*dxplot
          YDET(jpix)=(jpix-1)*dxplot
 12	CONTINUE
	
C  3. Bin WF intensity values into pixels

	DO 4 j=1,nWF
	  yij=(FLOAT(j)-0.5e0)*dxWF-yc
	  jpix=(yij+ycpix)/dxpix+1e0
	  IF ((jpix.GE.0).AND.(jpix.LE.(npix+1))) THEN
	    ypix=(FLOAT(jpix)-0.5e0)*dxpix-ycpix
	    yWFmin=yij-dxWF/2e0
	    ypixmin=ypix-dxpix/2e0
	    ypixmax=ypix+dxpix/2e0
	    yWFmax=yij+dxWF/2e0

	    IF (yWFmin.LT.yPixmin) THEN
	      dy1=yPixmin-yWFmin
	      dy2=yPixmin
	    ELSE
	      dy1=0d0
	      dy2=yWFmin
	    END IF
	    IF (yWFmax.GT.yPixmax) THEN
	      dy3=yWFmax-yPixmax
	      dy2=yPixmax-dy2
	    ELSE
	      dy3=0d0
	      dy2=yWFmax-dy2
	    END IF

	    jp1=jpix+1; jm1=jpix-1

	    DO 3 i=1,nWF
	      xij=(FLOAT(i)-0.5e0)*dxWF-xc
	      ipix=(xij+xcpix)/dxpix+1e0

	      IF ((ipix.GE.0).AND.(ipix.LE.(npix+1))) THEN
	        xpix=(FLOAT(ipix)-0.5e0)*dxpix-xcpix
	        xWFmin=xij-dxWF/2e0
	        xpixmin=xpix-dxpix/2e0
	        xpixmax=xpix+dxpix/2e0
	        xWFmax=xij+dxWF/2e0

	        IF (xWFmin.LT.xPixmin) THEN
	          dx1=xPixmin-xWFmin
	          dx2=xPixmin
	        ELSE
	          dx1=0d0
	          dx2=xWFmin
	        END IF
	        IF (xWFmax.GT.xPixmax) THEN
	          dx3=xWFmax-xPixmax
	          dx2=xPixmax-dx2
	        ELSE
	          dx3=0d0
	          dx2=xWFmax-dx2
	        END IF

	        A1=dx1*dy3; A2=dx2*dy3; A3=dx3*dy3
	        A4=dx1*dy2; A5=dx2*dy2; A6=dx3*dy2
	        A7=dx1*dy1; A8=dx2*dy1; A9=dx3*dy1

	        ip1=ipix+1; im1=ipix-1
	        IWF=WFArray(i,j)/A

	        IF (jp1.LE.npix) THEN
	          IF (im1.GE.1) 
     &	            PixArray(im1,jp1)=PixArray(im1,jp1)+A1*IWF
	          IF ((ipix.GE.1).AND.(ipix.LE.npix)) 
     &	            PixArray(ipix,jp1)=PixArray(ipix,jp1)+A2*IWF
	          IF (ip1.LE.npix) 
     &	            PixArray(ip1,jp1)=PixArray(ip1,jp1)+A3*IWF
	        END IF
	        IF ((jpix.GE.1).AND.(jpix.LE.npix)) THEN
	          IF (im1.GE.1) 
     &	            PixArray(im1,jpix)=PixArray(im1,jpix)+A4*IWF
	          IF ((ipix.GE.1).AND.(ipix.LE.npix))
     &	            PixArray(ipix,jpix)=PixArray(ipix,jpix)+A5*IWF
	          IF (ip1.LE.npix) 
     &	            PixArray(ip1,jpix)=PixArray(ip1,jpix)+A6*IWF
	        END IF
	        IF (jm1.GE.1) THEN
	          IF (im1.GE.1) 
     &	            PixArray(im1,jm1)=PixArray(im1,jm1)+A7*IWF
	          IF ((ipix.GE.1).AND.(ipix.LE.npix)) 
     &	            PixArray(ipix,jm1)=PixArray(ipix,jm1)+A8*IWF
	          IF (ip1.LE.npix) 
     &	            PixArray(ip1,jm1)=PixArray(ip1,jm1)+A9*IWF
	        END IF
	      END IF
  3	    CONTINUE
	  END IF
  4	CONTINUE

	RETURN
	END SUBROUTINE CPIXILATE

C***********************************************************************

	SUBROUTINE CPIX_J(PixArray,WFArray,XDET,YDET,drx,dry,
     &	  npix,nWF,dxpix,iPloTyp,RayPos,RayID,
     &    psidet,vdet,LRayOk,ifPixLoc,ifAdd)

c This routine takes a wavefront intensity array contained in WFArray
c and the location of the rays in RayPos array at the last element 
c The location of the rays are mapped into detector coordinates using
c subroutine detmap (to be written).
c It is assumed that the axis of the detector coordinate system is near 
c parallel to the axis of the wavefront grid.  
c The dx between wavefront data is determined by the distance 
c between it and adjacent rays (given in AdjRay and RayID arrays) 
c The detector center is chosen to be at the chfray location if
c ifPixLoc is .TRUE., if not it is determined by array detcen.

	use param_mod
	use math_mod

	IMPLICIT NONE


        LOGICAL LRayOk(mRay),ifPixLoc,ifAdd
	INTEGER i,j,iPloTyp,ipix,jpix,npix,nWF,ip1,im1,jp1,jm1

	SREAL WFArray(nWF,nWF),PixArray(mPix,mPix),xcpix,ycpix,
     &	XDET(mPix),YDET(mPix),dxpix,xij,yij,dxplot,xpix,ypix,
     &	xWFmin,xpixmin,xpixmax,xWFmax,yWFmin,ypixmin,ypixmax,yWFmax,
     &	A,A1,A2,A3,A4,A5,A6,A7,A8,A9,IWF,dx1,dx2,dx3,dy1,dy2,dy3,
     &	drx,dry

c cpix declerations
        real*8 rayp(3), psidet(3), vdet(3),ddrx,ddry
        real*8 detcoord(3),detcoordmx(3),detcoordmy(3),detcoordpx(3)
     &           ,detcoordpy(3)
        integer RayID(mdttl,mdttl)
        real*8 RayPos(3,mRay)
        real*8 raychf(3),detcen(3),xdetvec(3),ydetvec(3)
        real*8 raypmx(3),raypmy(3),rayppx(3),rayppy(3)
        integer ii,iRayCR,iCR,jCR

c cpix_checker declerations

        real *8 dWFxx,dWFyy,dWFxxmax,dWFxxmin,dWFyymax,dWFyymin
        real *8 dxp,dxm,dyp,dym,dxerr,dyerr
        real *8 avgdxerr,avgdyerr,maxdxerr,maxdyerr
        integer errcnt

C  1. Compute various parameters
        print *, 'iamin cpix'
c	xc=(FLOAT(nWF)-1e0)*dxWF/2e0
c	yc=(FLOAT(nWF)-1e0)*dxWF/2e0
	xcpix=FLOAT(npix-1)*dxpix/2e0
	ycpix=FLOAT(npix-1)*dxpix/2e0
c	A=dxWF*dxWF

c sws ifpixloc needs a value
c        ifPixLoc = .false.
c Outdated?:
	IF (ifPixLoc) THEN
	  ddrx=drx
	  ddry=dry
	ELSE
	  ddrx=0d0
	  ddry=0d0
	END IF

c compute othervec as the vector from the WF center to adjacent pixel

        iCR=mdttl/2+1.0
        jCR=iCR
        print *, 'cr',iCR,jCR
        iRayCR=RayID(iCR,jCR)
        print *, 'iRayCR', iRayCR
        do 17 i=1,3
           xdetvec(i)=RayPos(i,iRayCR+1)-RayPos(i,iRayCR)
 17     continue
        call DUNITIZE(xdetvec)
        call DUNITIZE(psidet)
        
        call DXPROD(ydetvec,xdetvec,psidet)
        print *, 'xdetvec', xdetvec
        print *, 'psidet', psidet
        print *, 'ydetvec',ydetvec

c center detector to chief ray
        call dequate(detcen,RayPos(1,1),3)
        if (ifPixLoc) then
            call dequate(raychf,RayPos(1,1),3)
            call detmap(ddrx,ddry,raychf,xdetvec,ydetvec,vdet,detcen)
c            xc=raychf(1)
c            yc=raychf(2)
        endif


C  2. Initialize array

	IF (.NOT.ifAdd) THEN
	  DO 2 jpix=1,npix
	    DO 1 ipix=1,npix
	      PixArray(ipix,jpix)=0e0
  1	    CONTINUE
  2	  CONTINUE
	END IF
	dxplot=float(npix)/float(npix-1)
	DO 241 jpix=1,npix
          XDET(jpix)=(jpix-1)*dxplot
          YDET(jpix)=(jpix-1)*dxplot
  241	CONTINUE

c 72. Preprocessor and linearity checker


        dWFxxmin=1.e10
        dWFyymin=1.e10
        dWFxxmax=0.0
        dWFyymax=0.0
        avgdxerr=0.0
        avgdyerr=0.0
        errcnt=0
        maxdxerr=0.0
        maxdyerr=0.0

	DO 24 j=1,nWF
	  DO 23 i=1,nWF

C  3.a. Compute location and pixel indices of WF grid point on detector
            if (RayID(i,j).eq.0) go to 23

c get position of WF array in detector coordinates
           do  25 ii=1,3
            rayp(ii)=RayPos(ii,RayID(i,j))
 25        continue

            call detmap(ddrx,ddry,rayp,xdetvec,ydetvec,vdet,detcoord)

c get adjacent WF grid pts
c ignore edge rays

            if ((RayID(i-1,j).eq.0).OR.(RayID(i+1,j).eq.0).OR.
     &              (RayID(i,j-1).eq.0).OR.(RayID(i,j+1).eq.0))
     &              then
               go to 23
               else

            if (.NOT.((LRayOk(RayID(i-1,j))).AND.
     &         (LRayOK(RayID(i+1,j))).AND.
     &         (LRayOk(RayID(i,j-1))).AND.
     &         (LRayOk(RayID(i,j+1)))))
     &      go to 23
            endif

            call dequate(raypmx,RayPos(1,RayID(i-1,j)),3)
            call detmap(ddrx,ddry,raypmx,xdetvec,ydetvec,vdet,
     &      detcoordmx)
            call dequate(rayppx,RayPos(1,RayID(i+1,j)),3)
            call detmap(ddrx,ddry,rayppx,xdetvec,ydetvec,vdet,
     &      detcoordpx)
            call dequate(raypmy,RayPos(1,RayID(i,j-1)),3)
            call detmap(ddrx,ddry,raypmy,xdetvec,ydetvec,vdet,
     &      detcoordmy)
            call dequate(rayppy,RayPos(1,RayID(i,j+1)),3)
            call detmap(ddrx,ddry,rayppy,xdetvec,ydetvec,vdet,
     &      detcoordpy)


c compute WF coordinates to pixel coordinates
            xij=detcoord(1)-detcen(1)
            yij=detcoord(2)-detcen(2)


c	    ipix=ANINT((xij+xcpix)/dxpix+1E0)
c	    IF ((ipix.GE.0).AND.(ipix.LE.(npix+1))) THEN
c	      jpix=ANINT((yij+ycpix)/dxpix+1E0)
c	      IF ((jpix.GE.0).AND.(jpix.LE.(npix+1))) THEN

C  3.b. Compute grid and pixel dimension extremes

            xWFmin=xij-abs((detcoord(1)-detcoordmx(1)))/2e0
            xWFmax=xij+abs((detcoordpx(1)-detcoord(1)))/2e0
            yWFmin=yij-abs((detcoord(2)-detcoordmy(2)))/2e0
            yWFmax=yij+abs((detcoordpy(2)-detcoord(2)))/2e0

c compute maximum and minium dWFx and dWFy            

            dWFxx=abs(xWFmax-xWFmin)
            dWFyy=abs(yWFmax-yWFmin)

            if (dWFxx.lt.dWFxxmin) dWFxxmin=dWFxx
            if (dWFxx.gt.dWFxxmax) dWFxxmax=dWFxx
            if (dWFyy.lt.dWFyymin) dWFyymin=dWFyy
            if (dWFyy.gt.dWFyymax) dWFyymax=dWFyy

c compute amount of non-rectangularity


            if ((RayID(i-1,j-1).eq.0).OR.(RayID(i+1,j-1).eq.0).OR.
     &              (RayID(i-1,j+1).eq.0).OR.(RayID(i+1,j+1).eq.0))
     &              then
               go to 23
               else

            if (.NOT.((LRayOk(RayID(i-1,j-1))).AND.
     &         (LRayOK(RayID(i+1,j-1))).AND.
     &         (LRayOk(RayID(i-1,j+1))).AND.
     &         (LRayOk(RayID(i+1,j+1)))))
     &      go to 23
            endif


            call dequate(raypmx,RayPos(1,RayID(i-1,j-1)),3)
            call detmap(ddrx,ddry,raypmx,xdetvec,ydetvec,vdet,
     &      detcoordmx)
            call dequate(rayppx,RayPos(1,RayID(i+1,j-1)),3)
            call detmap(ddrx,ddry,rayppx,xdetvec,ydetvec,vdet,
     &      detcoordpx)
            call dequate(raypmy,RayPos(1,RayID(i-1,j+1)),3)
            call detmap(ddrx,ddry,raypmy,xdetvec,ydetvec,vdet,
     &      detcoordmy)
            call dequate(rayppy,RayPos(1,RayID(i+1,j+1)),3)
            call detmap(ddrx,ddry,rayppy,xdetvec,ydetvec,vdet,
     &      detcoordpy)
            
            dxp=abs(detcoordpy(1)-detcoordmy(1))
            dxm=abs(detcoordmx(1)-detcoordpx(1))
            dyp=abs(detcoordpy(2)-detcoordpx(2))
            dym=abs(detcoordmy(2)-detcoordmx(2))

            dxerr=((dxp+dxm)-4.0*(dWFxx))/(2*dWFxx)
            dyerr=((dyp+dym)-4.0*(dWFyy))/(2*dWFyy)

            avgdxerr=avgdxerr+dxerr
            avgdyerr=avgdyerr+dyerr
            errcnt=errcnt+1

            if (dxerr.gt.maxdxerr) maxdxerr=dxerr
            if (dyerr.gt.maxdyerr) maxdyerr=dyerr

c             print *, 'derr=',dxerr,dyerr
            
c            write (27,*) '-,-',detcoordmx
c            write (27,*) '-,+',detcoordmy
c            write (27,*) '+,-',detcoordpx
c            write (27,*) '+,+',detcoordpy
c            write (27,*) i,j, dxerr,dyerr 
c compute area based on new dx and dy
            A = abs((detcoordpx(1)-detcoordmx(1))*(detcoordpy(2)-
     &           detcoordmy(2))/4.0)

           
 23     continue
 24     continue

        avgdxerr=avgdxerr/errcnt
        avgdyerr=avgdyerr/errcnt

        print *, 'dxmax=',dWFxxmax,'dxmin=',dWFxxmin
        print *, 'dymax=',dWFyymax,'dymin=',dWFyymin
        print *, 'maxerror =',maxdxerr,maxdyerr
        print *, 'avgerr =', avgdxerr,avgdyerr
	
C  3. Bin WF intensity values into pixels

c	DO 4 j=nWF/2-1,nWF/2+1
c	  DO 3 i=nWF/2-1,nWF/2+1

	DO 4 j=1,nWF
	  DO 3 i=1,nWF

C  3.a. Compute location and pixel indices of WF grid point on detector
            if (RayID(i,j).eq.0) go to 3

c get position of WF array in detector coordinates
           do 15 ii=1,3
            rayp(ii)=RayPos(ii,RayID(i,j))
 15        continue

            call detmap(ddrx,ddry,rayp,xdetvec,ydetvec,vdet,detcoord)

c get adjacent WF grid pts
c ignore edge rays

            if ((RayID(i-1,j).eq.0).OR.(RayID(i+1,j).eq.0).OR.
     &              (RayID(i,j-1).eq.0).OR.(RayID(i,j+1).eq.0))
     &              then
               go to 3
            else

            if (.NOT.((LRayOk(RayID(i-1,j))).AND.
     &         (LRayOK(RayID(i+1,j))).AND.
     &         (LRayOk(RayID(i,j-1))).AND.
     &         (LRayOk(RayID(i,j+1)))))
     &      go to 3
            endif



            call dequate(raypmx,RayPos(1,RayID(i-1,j)),3)
            call detmap(ddrx,ddry,raypmx,xdetvec,ydetvec,vdet,
     &      detcoordmx)
            call dequate(rayppx,RayPos(1,RayID(i+1,j)),3)
            call detmap(ddrx,ddry,rayppx,xdetvec,ydetvec,vdet,
     &      detcoordpx)
            call dequate(raypmy,RayPos(1,RayID(i,j-1)),3)
            call detmap(ddrx,ddry,raypmy,xdetvec,ydetvec,vdet,
     &      detcoordmy)
            call dequate(rayppy,RayPos(1,RayID(i,j+1)),3)
            call detmap(ddrx,ddry,rayppy,xdetvec,ydetvec,vdet,
     &      detcoordpy)
           
c            print *, 'mx',detcoordmx
c            print *, 'my',detcoordmy
c            print *, 'px',detcoordpx
c            print *, 'py',detcoordpy


c	    xij=(FLOAT(i)-1e0)*dxWF-detcen(1)
c	    yij=(FLOAT(j)-1e0)*dxWF-detcen(2)



c compute WF coordinates to pixel coordinates
            xij=detcoord(1)-detcen(1)
            yij=detcoord(2)-detcen(2)


	    ipix=ANINT((xij+xcpix)/dxpix+1E0)
	    IF ((ipix.GE.0).AND.(ipix.LE.(npix+1))) THEN
	      jpix=ANINT((yij+ycpix)/dxpix+1E0)
	      IF ((jpix.GE.0).AND.(jpix.LE.(npix+1))) THEN

C  3.b. Compute grid and pixel dimension extremes

            xWFmin=xij-abs((detcoord(1)-detcoordmx(1)))/2e0
            xWFmax=xij+abs((detcoordpx(1)-detcoord(1)))/2e0
            yWFmin=yij-abs((detcoord(2)-detcoordmy(2)))/2e0
            yWFmax=yij+abs((detcoordpy(2)-detcoord(2)))/2e0



c compute area based on new dx and dy
            A = abs((detcoordpx(1)-detcoordmx(1))*(detcoordpy(2)-
     &           detcoordmy(2))/4.0)

	        xpix=(FLOAT(ipix)-1e0)*dxpix-xcpix
c	        xWFmin=xij-dxWF/2e0
	        xpixmin=xpix-dxpix/2e0
	        xpixmax=xpix+dxpix/2e0
c	        xWFmax=xij+dxWF/2e0
	        ypix=(FLOAT(jpix)-1e0)*dxpix-ycpix
c	        yWFmin=yij-dxWF/2e0
	        ypixmin=ypix-dxpix/2e0
	        ypixmax=ypix+dxpix/2e0
c	        yWFmax=yij+dxWF/2e0

C  3.c. Compute overlap conditions and areas

	        IF (xWFmin.LT.xPixmin) THEN
	          dx1=xPixmin-xWFmin
	          dx2=xPixmin
	        ELSE
	          dx1=0d0
	          dx2=xWFmin
	        END IF
	        IF (xWFmax.GT.xPixmax) THEN
	          dx3=xWFmax-xPixmax
	          dx2=xPixmax-dx2
	        ELSE
	          dx3=0d0
	          dx2=xWFmax-dx2
	        END IF

	        IF (yWFmin.LT.yPixmin) THEN
	          dy1=yPixmin-yWFmin
	          dy2=yPixmin
	        ELSE
	          dy1=0d0
	          dy2=yWFmin
	        END IF
	        IF (yWFmax.GT.yPixmax) THEN
	          dy3=yWFmax-yPixmax
	          dy2=yPixmax-dy2
	        ELSE
	          dy3=0d0
	          dy2=yWFmax-dy2
	        END IF

	        A1=dx1*dy3
	        A2=dx2*dy3
	        A3=dx3*dy3
	        A4=dx1*dy2
	        A5=dx2*dy2
	        A6=dx3*dy2
	        A7=dx1*dy1
	        A8=dx2*dy1
	        A9=dx3*dy1

C  3.d. Bin the pixel intensities

	        ip1=ipix+1
	        im1=ipix-1
	        jp1=jpix+1
	        jm1=jpix-1
	        IWF=WFArray(i,j)/A

	        IF (jp1.LE.npix) THEN
	          IF (im1.GE.1) 
     &	            PixArray(im1,jp1)=PixArray(im1,jp1)+A1*IWF
	          IF ((ipix.GE.1).AND.(ipix.LE.npix)) 
     &	            PixArray(ipix,jp1)=PixArray(ipix,jp1)+A2*IWF
	          IF (ip1.LE.npix) 
     &	            PixArray(ip1,jp1)=PixArray(ip1,jp1)+A3*IWF
	        END IF
	        IF ((jpix.GE.1).AND.(jpix.LE.npix)) THEN
	          IF (im1.GE.1) 
     &	            PixArray(im1,jpix)=PixArray(im1,jpix)+A4*IWF
	          IF ((ipix.GE.1).AND.(ipix.LE.npix))
     &	            PixArray(ipix,jpix)=PixArray(ipix,jpix)+A5*IWF
	          IF (ip1.LE.npix) 
     &	            PixArray(ip1,jpix)=PixArray(ip1,jpix)+A6*IWF
	        END IF
	        IF (jm1.GE.1) THEN
	          IF (im1.GE.1) 
     &	            PixArray(im1,jm1)=PixArray(im1,jm1)+A7*IWF
	          IF ((ipix.GE.1).AND.(ipix.LE.npix)) 
     &	            PixArray(ipix,jm1)=PixArray(ipix,jm1)+A8*IWF
	          IF (ip1.LE.npix) 
     &	            PixArray(ip1,jm1)=PixArray(ip1,jm1)+A9*IWF
	        END IF
	      END IF
	    END IF
  3	  CONTINUE
  4	CONTINUE

	RETURN
	END SUBROUTINE CPIX_J

C***********************************************************************


        subroutine detmap(ddrx,ddry,xin,xdetvec,ydetvec,vdet,xout)
        use math_mod
 
        real*8 xin(3),xdetvec(3),ydetvec(3),vdet(3),xout(3)
        real*8 delvec(3),zdetvec(3)
        real*8 ddrx,ddry
        CALL DSUB(delvec,xin,vdet,3)
        xout(1)=DDOTC(delvec,xdetvec)-ddrx
        xout(2)=DDOTC(delvec,ydetvec)-ddry
        CALL DXPROD(zdetvec,xdetvec,ydetvec)
        xout(3)=DDOTC(delvec,zdetvec)

c        print *, 'xout', xout
        return
        end subroutine detmap

C***********************************************************************

	SUBROUTINE CBLUR(PixArray,npix,mPix,kval,R1)
	use math_mod

	IMPLICIT NONE

	INTEGER npix,mPix,i,j,k,l,ik,jl

	SREAL PixArray(mPix,mPix),kval,R1(mPix,mPix),
     &	da,b,c,kernel(-1:1,-1:1),pv

	IF ((kval.GT.1e0).OR.(kval.LT.0e0)) RETURN

	da=1e0-kval
	b=.2*da
	c=.25*b

	kernel(-1,-1)=c
	kernel(-1,0)=b
	kernel(-1,1)=c
	kernel(0,-1)=b
	kernel(0,0)=kval
	kernel(0,1)=b
	kernel(1,-1)=c
	kernel(1,0)=b
	kernel(1,1)=c

	DO j=1,nPix
	  DO i=1,nPix
	    R1(i,j)=0e0
	  END DO
	END DO

	DO j=1,nPix
	  DO i=1,nPix
	    pv=PixArray(i,j)
	    DO l=-1,1
	      jl=j+l
	      DO k=-1,1
	        ik=k+i
	        IF ((ik.GE.1).AND.(ik.LE.nPix).AND.(jl.GE.1).AND.
     &	        (jl.LE.nPix)) THEN
	          R1(ik,jl)=R1(ik,jl)+pv*kernel(k,l)
	        END IF
	      END DO
	    END DO
	  END DO
	END DO

	DO j=1,nPix
	  DO i=1,nPix
	    PixArray(i,j)=R1(i,j)
	  END DO
	END DO	        

	RETURN
	END SUBROUTINE CBLUR

C***********************************************************************

	SUBROUTINE GBLUR(PixArray,npix,mPix,kernel,mk,krad,R1)
	use math_mod

	IMPLICIT NONE

	INTEGER npix,mPix,i,j,k,l,ik,jl,mk,jj

	SREAL PixArray(mPix,mPix),kval,R1(mPix,mPix),
     &	  da,b,c,kernel(-mk:mk,-mk:mk),pv,r2,krad,c1

	IF (krad.LE.0e0) RETURN

	c1=-1.17741002251547/2e0/krad

	DO j=0,mk
	  jj=j*j
	  DO i=0,mk
	    r2=i*i+jj
	    kernel(i,j)=0.39894228040143*EXP(r2*c1)
	    kernel(-i,j)=kernel(i,j)
	    kernel(i,-j)=kernel(i,j)
	    kernel(-i,-j)=kernel(i,j)
	  END DO
	END DO

	DO j=1,nPix
	  DO i=1,nPix
	  END DO
	END DO

	DO j=1,nPix
	  DO i=1,nPix
	    pv=PixArray(i,j)
	    DO l=-mk,mk
	      jl=j+l
	      DO k=-mk,mk
	        ik=k+i
	        IF ((ik.GE.1).AND.(ik.LE.nPix).AND.(jl.GE.1).AND.
     &	        (jl.LE.nPix)) THEN
	          R1(ik,jl)=R1(ik,jl)+pv*kernel(k,l)
	        END IF
	      END DO
	    END DO
	  END DO
	END DO

	DO j=1,nPix
	  DO i=1,nPix
	    PixArray(i,j)=R1(i,j)
	  END DO
	END DO	        

	RETURN
	END SUBROUTINE GBLUR

C***********************************************************************
C*******End file pixsub.F***********************************************
C***********************************************************************

