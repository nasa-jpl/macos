C***********************************************************************
C	Begin file linsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************


#include "realtype.h"

	SUBROUTINE CLPERTURB(PertVec,RaySpot,OPDMat,AvgL,Cmatrix,npts,
     & 			     RayIndex,iPertElt,iElt)
	use param_mod
	use elt_mod
	use math_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

	IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: ANS
	INTEGER iaxis,iElt,npts,RayIndex(2,mRay),iPertElt
	SREAL :: Cmatrix(7,mCm,bRay),PertVec(6,mElt)
	REAL*8 RaySpot(mRay,2),D1(9),PV_loc(6),DDEF(9),AvgL,
     &         OPDMat(mpts,mpts)

 100    FORMAT(1P,'   Rotational= ',3d14.6/'Translational= ',3d14.6)
 101    FORMAT(/)
 102    FORMAT(' Element ',i4,' Linear Perturbations:')
 400    FORMAT(' Tracing perturbed system')
 401    FORMAT(/' PERTURB expects 5 output coordinates in the order:'/
     &  ' 1-2: direction (dr); 3-4: beamwalk (gamma); 5: OPD (dL).'/
     &  ' Please MODIFY the last element to create the proper Tout.'/)
 402    FORMAT(1P,' RMS OPD including piston is ',e13.6)
 403    FORMAT(1P,' RMS WFE excluding piston is ',e13.6)
 404    FORMAT(1P,' Average OPD is ',e13.6)

	CALL DZERO(DDEF,9)

C  1. Check nOutCord

	IF (nOutCord.NE.5) THEN
	  WRITE(*,401)
	  RETURN
	END IF

C  2. Get perturbation data

	  CALL CACCEPT(ANS,'NO','Use element coordinates? ')
	  IF ((ANS(1:1).EQ.'y').OR.(ANS(1:1).EQ.'Y')) THEN
	    CALL DACCEPT(D1,DDEF,6,
     &      'Enter perturbation vector in element coordinates:')
	    CALL DMPROD(PV_loc,TElt(1,1,iElt),D1,
     &      6,nECoord(iElt),1)
	  ELSE
	    CALL DACCEPT(PV_loc(1:3),DDEF,3,
     &      'Enter rotational perturbation vector (x,y,z):')
	    CALL DACCEPT(PV_loc(4:6),DDEF,3,
     &      'Enter translational perturbation vector (x,y,z):')
	  END IF
	  DO 11 iaxis=1,6
	      PertVec(iaxis,iElt)=PV_loc(iaxis)
 11       CONTINUE

	WRITE(*,101)
	DO 5 iElt=1,nElt
	  WRITE(*,102) iElt
	  WRITE(*,100) (PertVec(iaxis,iElt),iaxis=1,6)
 5      CONTINUE

	RETURN
	END SUBROUTINE CLPERTURB

C***********************************************************************

	SUBROUTINE CLPLOT(command,npts,RaySpot,
     &  		  OPDMat,AvgL,RayIndex,RayMap,iPloTyp,
     &			  iLinElt,CD1,CD2,PixArray,StretchType)

	use param_mod
	use cfiles_mod
	use elt_mod
        use math_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

	IMPLICIT NONE

	CHARACTER(len=MacosCharLen) :: COMMAND
	LOGICAL LCMP

	INTEGER i,j,k,l,npts,iPloTyp,iElt,ICLEN,
     &  	RayIndex(2,mRay),nGridPts,iLinElt,
     &  	NN(2),pad,RayMap(mpts,mpts),npix,StretchType

	SREAL :: lambda,TransDis,Intensity,phi,udia,du,
     &           dx,xdia,X0PFP,Y0PFP,dxWF,dxpix,PixArray(mPix,mPix)

	REAL*8 RaySpot(mRay,2),OPDMat(mpts,mpts),AvgL,tpl

	COMPLEX*8 CD1(mdttl,mdttl),CD2(mdttl,mdttl)
	CHARACTER*8 CInteger
	CHARACTER*64 CTITLE

C       File setup and initial conditions read

 202	FORMAT(' Command not understood.  Please reenter...')
 206	FORMAT(1P,' FFT/Point Spread Function Data Summary:'/
     &	' Wavelength=',E17.10,';  Transmission Distance=',E17.10/
     &	' u-v plane diam=',E17.10,'  du=',E17.10/
     &	' x-y plane diam=',E17.10,'  dx=',E17.10)
 207    FORMAT(' Enter a no-blanks case description: ')

C  1. Command loop

	iElt=iLinElt
	nGridPts=npts+1

	lambda=Wavelen
	TransDis=zElt(iLinElt)
	Intensity=Flux/FLOAT(nRay)

	tpl=6.28318530717959D0/Wavelen
	du=dxElt(iLinElt)
	udia=du*FLOAT(mdttl-1)
	dx=lambda*TransDis/du/FLOAT(mdttl)
	xdia=FLOAT(mdttl-1)*dx
	WRITE(*,206)lambda,TransDis,udia,du,xdia,dx

	X0PFP=RaySpot(1,1)-dx*FLOAT(mdttl/2)
	Y0PFP=RaySpot(1,2)-dx*FLOAT(mdttl/2)

C  7. Pad OPD map for FFT

	pad=(mdttl-nGridPts+1)/2
	DO 2 j=1,mdttl
	  DO 1 i=1,mdttl
	    CD1(i,j)=0e0
 1	  CONTINUE
	  XFP(j)=X0PFP
	  YFP(j)=Y0PFP
	  X0PFP=X0PFP+dx
	  Y0PFP=Y0PFP+dx
 2	CONTINUE
	DO 4 j=1,mpts
	  l=j+pad
	  DO 3 i=1,mpts
	    IF (RayMap(i,j).NE.0) THEN
	      k=i+pad
	      phi=tpl*OPDMat(i,j)
	      CD1(k,l)=Intensity*CMPLX(cos(phi),sin(phi))
	    END IF
 3	  CONTINUE
 4	CONTINUE

C  8. Compute PSF

	NN(1)=mdttl
	NN(2)=mdttl
	CALL RSWAP2(CD2,CD1,mdttl)
	CALL RFOURN(CD2,NN,2,1)
	CALL RSWAP2(CD1,CD2,mdttl)
	k=0
	DO 6 j=1,mdttl
	  DO 5 i=1,mdttl
	    MWFFT(i,j)=(CD1(i,j)*CONJG(CD1(i,j)))
 5	  CONTINUE
 6	CONTINUE

	IF (LCMP(command,'LPI',2)) THEN
	  i=mdttl/2
	  CALL IACCEPT_S(npix,i,1,'Enter number of pixels per side:')
	  dxWF=dx
	  CALL RACCEPT_S(dxpix,dxWF,1,'Enter size of pixel:')
	  CALL CPIXILATE(PixArray,MWFFT,XFP,YFP,0d0,0d0,
     &	  npix,mPix,mdttl,dxWF,dxpix,iPloTyp,.FALSE.,.FALSE.)
	  CALL IntToChar(Cinteger,iLinElt,i)
	  IF (iPloTyp.GE.5) THEN
	    j=ICLEN(filnam)
	    imgfil=filnam(1:j)//'.lpx'//Cinteger(1:i)
	    CTITLE='Pixellated Intensity, Elt='
     &	    //Cinteger(1:i)
	    CALL BINOUT(iPloTyp,XFP,YFP,PixArray,npix,npix,
     &	                mPix,mPix,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                  dxElt(iElt), dxElt(iElt), iElt)
	  ELSE
	    IF (StretchType.EQ.2) THEN
	      CTITLE='LOG10 Linear Model Pixilated Intensity#%'
	    ELSE IF (StretchType.EQ.3) THEN
	      CTITLE='SQRT Linear Model Pixilated Intensity#%'
	    ELSE
	      CTITLE='Linear Model Pixilated Intensity#%'
	    END IF
	    CALL PixOut(PixArray,iElt,iPloTyp,CTITLE,
     &	    'X-Axis','Y-Axis','.pix',npix,mpix,REAL(CD1,KIND=8),
     &	    REAL(CD2,KIND=8),StretchType)
	  END IF

	ELSE IF (LCMP(command,'LLO',3)) THEN
	  CALL IntToChar(Cinteger,iLinElt,i)
	  IF (iPloTyp.GE.5) THEN
	    j=ICLEN(filnam)
	    imgfil=filnam(1:j)//'.llo'//Cinteger(1:i)
	    CTITLE='Linear Model Log10 Image Intensity, Elt='
     &	    //Cinteger(1:i)
	    CALL BINOUT(iPloTyp,XFP,YFP,MWFFT,mdttl,mdttl,
     &	                mdttl,mdttl,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                  dxElt(iElt), dxElt(iElt), iElt)
	  ELSE
	    CTITLE='Linear Model Log10 Image Intensity, Elt='
     &	    //Cinteger(1:i)//', File='//filnam
	    DO j=1,mdttl
	      DO i=1,mdttl
	        MWFFT(i,j)=LOG10(MWFFT(i,j))
	      END DO
	    END DO
	    CALL SRFPLT(iPloTyp,XFP,YFP,MWFFT,
     &	    mdttl,mdttl,mdttl,mdttl,'X-Axis','Y-Axis',CTITLE,
     &	    CD1,CD2,1)
	  END IF

	ELSE
	  CALL IntToChar(Cinteger,iLinElt,i)
	  IF (iPloTyp.GE.5) THEN
	    j=ICLEN(filnam)
	    imgfil=filnam(1:j)//'.lin'//Cinteger(1:i)
	    CTITLE='Linear Model Point Image Intensity, Elt='
     &	    //Cinteger(1:i)
	    CALL BINOUT(iPloTyp,XFP,YFP,MWFFT,mdttl,mdttl,
     &	                mdttl,mdttl,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                  dxElt(iElt), dxElt(iElt), iElt)
	  ELSE
	    IF (StretchType.EQ.2) THEN
	      CTITLE='LOG10 Linear Model Intensity#%'
	    ELSE IF (StretchType.EQ.3) THEN
	      CTITLE='SQRT Linear Model Intensity#%'
	    ELSE
	      CTITLE='Linear Model Intensity#%'
	    END IF
	    CALL SrfOut(iElt,iPloTyp,CTITLE,
     &	    'X-Axis','Y-Axis','.int',mdttl,mdttl,
     &	    REAL(CD1,KIND=8),REAL(CD2,KIND=8),StretchType)
	  END IF
	END IF

	RETURN
	END SUBROUTINE CLPLOT

C***********************************************************************

	SUBROUTINE CLTRACE(PertVec,RaySpot,OPDMat,AvgL,Cmatrix,
     &  		   npts,RayIndex,iPertElt,iElt)

	use param_mod
	use elt_mod
        use math_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

	IMPLICIT NONE

	INTEGER i,j,iaxis,iElt,iRay,npts,mpt2,
     &  	RayIndex(2,mRay),iPertElt

	SREAL :: Cmatrix(7,mCm,bRay),xout(7),PertVec(6,mElt),
     &         xdum(7),Tcompact(49),xnom(7),RMSWFE,SpotMag,RMSOPD

	REAL*8 AvgL,OPDMat(mpts,mpts),RaySpot(mRay,2)

 400    FORMAT(' Tracing perturbed system')
 401    FORMAT(/' PERTURB expects 5 output coordinates in the order:'/
     &  ' 1-2: direction (dr); 3-4: beamwalk (gamma); 5: OPD (dL).'/
     &  ' Please MODIFY the last element to create the proper Tout.'/)
 402    FORMAT(1P,' RMS OPD including piston is ',e13.6)
 403    FORMAT(1P,' RMS WFE excluding piston is ',e13.6)
 404    FORMAT(1P,' Average OPD is ',e13.6)

C  1. Check nOutCord

	IF (nOutCord.NE.5) THEN
	  WRITE(*,401)
	  RETURN
	END IF

	CALL COMPACT(Tcompact,Tout,5,7,7,7)

	SpotMag=1e0
c	CALL RACCEPT_S(SpotMag,1e0,1,
c     &  'Enter magnification factor for spot errors: ')

C  3. Compute perturbed rays

	WRITE(*,400)

	AvgL=0d0
	RMSWFE=0e0
	RMSOPD=0e0
	mpt2=mpts**2
	CALL DZERO(OPDMat,mpt2)
	DO 3 iRay=1,nRay
	  CALL RZERO(xdum,7)
	  DO 1 iElt=1,iPertElt
	    i=8+(iElt-1)*6
	    CALL RMPAC(xdum,Cmatrix(1,i,iRay),
     &      PertVec(1,iElt),7,6,1)
 1        CONTINUE
	  CALL RMPRDB(xout,Tcompact,xdum,5,7,1)
	  DO 2 iaxis=1,3
	    j=iaxis+3
	    xdum(iaxis)=RayDir(iaxis,iRay)
	    xdum(j)=RayPos(iaxis,iRay)
 2        CONTINUE
	  xdum(7)=0e0
	  CALL RMPRDB(xnom,Tcompact,xdum,5,7,1)
	  RaySpot(iRay,1)=xnom(3)+xout(3)*SpotMag
	  RaySpot(iRay,2)=xnom(4)+xout(4)*SpotMag
	  IF (iRay.GT.1) THEN
	    OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))=
     &      CumRayL(iRay)-CumRayL(1)+xout(5)
	    AvgL=AvgL+OPDMat(RayIndex(1,iRay),RayIndex(2,iRay))
	    RMSOPD=RMSOPD+
     &      (OPDMat(RayIndex(1,iRay),RayIndex(2,iRay)))**2
	  END IF
 3      CONTINUE
	AvgL=AvgL/DBLE(nRay-1)
	RMSOPD=SQRT(RMSOPD/FLOAT(nRay-1))
	WRITE(*,404) AvgL
	WRITE(*,402) RMSOPD

	RETURN
	END SUBROUTINE CLTRACE

C***********************************************************************

	SUBROUTINE CLPRead(PertVec,ifPertRead,iPertElt)

C       A routine to read input perturbation data.  The data is
C       assumed to be in the form of transposed u-vectors, starting
C       with element 1 and proceding to the last element.  Deformable
C       elements are loaded completely (all rays); non-deformable
C       elements are loaded for ray 1 only.

	use param_mod
	use cfiles_mod
	use elt_mod
#ifdef CSMACOS
        use smacosio_mod
#else
        use macosio_mod
#endif

	IMPLICIT NONE

	LOGICAL ifPertRead,exist
	INTEGER iaxis,iElt,iPertElt,ICLEN
	REAL*8 invec6(6)
	SREAL :: PertVec(6,mElt)

 100	FORMAT(' Input file ',A,' not compatible with current system:'/
     &	'  iPertElt=',i3,'  iElt=',i3)

	CALL CACCEPT(pertfil,pertfil(1:ICLEN(pertfil)),
     &	'Enter .pert-file name:')

 1	CONTINUE
	  INQUIRE (file=pertfil,exist=exist)
	  IF (.NOT.exist) THEN
	    WRITE(*,*)
     &     ' Input file '//pertfil(1:ICLEN(pertfil))//' does not exist.'
	    CALL CACCEPT(pertfil,'q','Enter new .pert-file name:')
	    IF ((ICLEN(pertfil).EQ.1).AND.(pertfil(1:1).EQ.'q')) THEN
	      ifPertRead=.FALSE.
	      RETURN
	    END IF
	    GO TO 1
	  ELSE
	    OPEN (7,FILE=PERTFIL,STATUS='OLD')
	    DO 4 iElt=1,iPertElt
	      READ(7,*,ERR=6,END=6)invec6
	      DO 3 iaxis=1,6
	        PertVec(iaxis,iElt)=invec6(iaxis)
 3            CONTINUE
 4          CONTINUE
	  END IF
	  ifPertRead=.TRUE.
	  WRITE(*,*)' File'//pertfil(1:ICLEN(pertfil))//' read.'
	  CLOSE(7)
	  RETURN

 6	CONTINUE
	WRITE(*,100)pertfil(1:ICLEN(pertfil)),iPertElt,iElt
	ifPertRead=.FALSE.
	CLOSE(7)
	RETURN

	END

C***********************************************************************
C*******End file linsub.f***********************************************
C***********************************************************************

