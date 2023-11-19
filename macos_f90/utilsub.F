C***********************************************************************
C	Begin file utilsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	SUBROUTINE IntToChar(Char,Int,len)
	IMPLICIT NONE
	CHARACTER*8 C1
	CHARACTER*8 Char
	INTEGER i,nend,nstart,len,Int
 1	FORMAT(i8)
	WRITE(C1,1)Int
	i=0
 2      CONTINUE
	  i=i+1
	  iF ((C1(i:i).EQ.' ').AND.(i.LT.8)) GO TO 2
	nstart=i
	nend=8
	Char=c1(nstart:nend)
	len=nend-nstart+1
	RETURN
	END SUBROUTINE IntToChar

C***********************************************************************

	FUNCTION LCMP(C1,C2,len)
	IMPLICIT NONE
	CHARACTER*(*)C1,C2
	LOGICAL LCMP
	INTEGER len,i,i1,i2

	LCMP=.TRUE.

	DO 1 i=1,len
	  IF (ICHAR(C1(i:i)).GE.97) THEN
	    i1=ICHAR(C1(i:i))-32
	  ELSE
	    i1=ICHAR(C1(i:i))
	  END IF
	  IF (ICHAR(C2(i:i)).GE.97) THEN
	    i2=ICHAR(C2(i:i))-32
	  ELSE
	    i2=ICHAR(C2(i:i))
	  END IF
	  IF (i1.NE.i2) THEN
	    LCMP=.FALSE.
	    RETURN
	  END IF
  1	CONTINUE
	RETURN
	END FUNCTION LCMP


! By jzlou
	FUNCTION LCMP2(C1,C2,clen)
	Implicit none
	Character(len=*) :: C1, C2
	Integer :: clen
	Character(len=clen) :: S1,S2
	Logical :: LCMP2
	S1=C1(1:clen); S2=C2(1:clen)
	LCMP2 = (.not. LGT(S1,S2)).and.(.not. LGT(S2,S1))
	END FUNCTION LCMP2


#ifdef not_defined
	FUNCTION STRLEN(S,MaxLen)
	CHARACTER*(*) :: S
	INTEGER :: MaxLen,i
	Len=0
 	Do i=1,MaxLen
	  if (S(i:i)/=' ') then
	    Len=Len+1
	  else
	    STRLEN=len
	    return
	  end if
	End Do
	END FUNCTION STRLEN
#endif

C***********************************************************************

       SUBROUTINE RMOUT(X,NA,NC)
       Implicit none

       SREAL X(NA,NC)
       Integer :: NA,NC,I,J
 102   FORMAT(1P,1X,8E13.4)
       DO 1 I=1,NA
	 WRITE(*,102)(X(I,J),J=1,NC)
 1     CONTINUE
       RETURN
       END SUBROUTINE RMOUT

C***********************************************************************

	FUNCTION ICLEN(CVAR)
	Implicit none
	CHARACTER*(*) CVAR
	INTEGER I,ICLEN
	I=0
 1      CONTINUE
	  I=I+1
	  IF (CVAR(I:I).NE.' ') GO TO 1
	ICLEN=I-1
	RETURN
	END FUNCTION ICLEN

C***********************************************************************

	FUNCTION ICLEN2(CVAR)
	Implicit none
	CHARACTER*(*) CVAR
	INTEGER I,ICLEN2
	I=0
 1      CONTINUE
	I=I+1
	IF ((CVAR(I:(I+1)).NE.'#%').AND.(I.LE.31)) GO TO 1
	ICLEN2=I-1
	RETURN
	END FUNCTION ICLEN2

C***********************************************************************

C***********************************************************************

	FUNCTION LINT(X)
	INTEGER LINT
	REAL*8 X,XP
	IF(X.GE.0)THEN
	  LINT=X
	ELSE
	  XP=X+1d-6
	  LINT=INT(XP)-1
	END IF
	RETURN
	END FUNCTION LINT

C***********************************************************************

	FUNCTION ClosestInt(X)
	Implicit None
	Real*8 :: X,R
	Integer :: ClosestInt
	R=X-Int(X);
        if (R<0.5) then
	  ClosestInt=Int(X)
	else
	  ClosestInt=Int(X)+1
	end if
	End FUNCTION ClosestInt

C***********************************************************************

	SUBROUTINE Ca2Int(ifVecDif3,power,MaxInt,imax,jmax)
	use param_mod
	use elt_mod

	Implicit None

	LOGICAL ifVecDif3
	INTEGER i,j,imax,jmax
	REAL*8 RS,MaxInt,power

	MaxInt=-1d33
	power=0e0
	IF (ifVecDif3) THEN
	  DO j=1,mdttl
	    DO i=1,mdttl
	      RS=WFElt(i,j,1)*DCONJG(WFElt(i,j,1))+
     &	         WFElt(i,j,2)*DCONJG(WFElt(i,j,2))+
     &	         WFElt(i,j,3)*DCONJG(WFElt(i,j,3))
	      MWFFT(i,j)=RS
	      power=power+RS
	      IF (RS.GT.MaxInt) THEN
	        imax=i; jmax=j
	        MaxInt=RS
	      END IF
  	    End Do
  	  End Do
	ELSE
	  Do j=1,mdttl
	    Do i=1,mdttl
	      RS=WFElt(i,j,1)*DCONJG(WFElt(i,j,1))
	      MWFFT(i,j)=RS
	      power=power+RS
	      IF (RS.GT.MaxInt) THEN
	        imax=i
	        jmax=j
	        MaxInt=RS
	      END IF
  	    End Do
  	  End Do
	END IF
	RETURN
	END SUBROUTINE Ca2Int

C***********************************************************************
	! Data: on input = Intensity data, on output = MTF data
	! dx (input): dx of intensity array
	! X,Y (output): spatial freq coords of MTF

	SUBROUTINE Int2MTF(Data,Dim,dx,X,Y)
	use param_mod
        use elt_mod

	Implicit None

	Integer :: Dim,nn(2),i,j
	SREAL :: Data(Dim,Dim),X(Dim),Y(Dim)
	Real*8 :: dx,df,S
	Complex*16 :: CData(Dim,Dim),DWF(Dim,Dim),W

        CData(1:Dim,1:Dim)=DCMPLX(Data(1:Dim,1:Dim))
	nn(1:2)=Dim
	CALL DSWAP2(DWF,CData,Dim)  ! this is needed because intensity data was DSWAP2,
				    ! so need to inverse that swap. -jzlou
	Call DFOURN(DWF,nn,2,1)
	CALL DSWAP2(CData,DWF,Dim)

	S=0d0;
	Do i=1,Dim
          Do j=1,Dim
	    Data(i,j) = sqrt(CData(i,j)*DCONJG(CData(i,j)))
            if (S<Data(i,j)) S=Data(i,j)
	  End Do
	End Do

	! Take norm of CData to get MTF
	Do i=1,Dim
	  Do j=1,Dim
	    W=CData(i,j)/S
	    Data(i,j) = W*DCONJG(W)
	    !Data(i,j) = sqrt(CData(i,j)*DCONJG(CData(i,j)))
	    !Data(i,j) = CData(i,j)*DCONJG(CData(i,j))
	  End Do
	End Do
	!Data=Data/(S*S)  ! normalize MTF

	df = 1d0/(Dim*dx)
	X(1)=-1d0/(2d0*dx); Y(1)=X(1)
	Do i=2,Dim
	  X(i)=X(i-1)+df; Y(i)=X(i)
	End Do
	END SUBROUTINE Int2MTF

C***********************************************************************

	SUBROUTINE Ca2Amp(iWF,power,MaxInt,imax,jmax)
	use param_mod
	use elt_mod

	Implicit None

	INTEGER i,j,iWF,imax,jmax
	REAL*8 RS,MaxInt,power
	MaxInt=-1d33
	power=0e0
	DO 2 j=1,mdttl
	  DO 1 i=1,mdttl
	    RS=ABS(WFElt(i,j,iWF))
	    MWFFT(i,j)=RS
	    power=power+RS
	    IF (RS.GT.MaxInt) THEN
	      imax=i
	      jmax=j
	      MaxInt=RS
	    END IF
 1	  CONTINUE
 2	CONTINUE

      if (.false.) then
! Test only -- write out WFElt array (pupil function data) -jzlou
	Open(unit=12,file='pupfun.txt')
	DO j=1,mdttl
	  WRITE(12,*) DREAL(WFElt(1:mdttl,j,iWF)),
     &                DIMAG(WFElt(1:mdttl,j,iWF))
	END DO
	Close(12)
      end if

	RETURN
	END SUBROUTINE Ca2Amp

C***********************************************************************

	SUBROUTINE Ca2Pha(iWF)
	use param_mod
	use elt_mod
	Implicit None

	INTEGER i,j,iWF
	REAL*8 RealPart,ImagPart,Mag,tmp,rtd,tpi
	DATA rtd/57.29577951308232d0/  ! radian to degree

	tpi = 4d0*atan(1d0)*2;  ! 2*pi, -jzlou

	DO 2 j=1,mdttl
	  DO 1 i=1,mdttl
	    RealPart=DBLE(WFElt(i,j,iWF))
	    ImagPart=DIMAG(WFElt(i,j,iWF))
	    Mag=SQRT(RealPart*RealPart+ImagPart*ImagPart)
	    IF (Mag.LT.1d-22) THEN
	      MWFFT(i,j)=0d0
	    ELSE
              tmp = ATAN2(ImagPart,RealPart)
	      if (tmp < 0d0) tmp = tmp + tpi;  ! atan2 has range [-pi,pi], but we need
				               ! the range [0,2*pi], which is 'input range';
	      !MWFFT(i,j)= tmp   ! in radian
	      MWFFT(i,j)=ATAN2(ImagPart,RealPart)*rtd  ! in degree
	      if (.false. .and. i==256 .and. j==130) then
	        print*,' *** Ca2Pha: ImagPart,RealPart =',ImagPart,RealPart
	        print*,' *** Ca2Pha: MWFFT(256,130) =',MWFFT(i,j)
	      end if
	    END IF
 1	  CONTINUE
 2	CONTINUE
	RETURN
	END SUBROUTINE Ca2Pha

C***********************************************************************

	SUBROUTINE Ca2Log(ifVecDif3,power,MaxInt,imax,jmax)
	use param_mod
	use elt_mod

	Implicit None

	LOGICAL ifVecDif3
	INTEGER i,j,imax,jmax
	REAL*8 RS,MaxInt,power
	MaxInt=-1d33
	IF (ifVecDif3) THEN
	  DO 4 j=1,mdttl
	    DO 3 i=1,mdttl
	      RS=WFElt(i,j,1)*DCONJG(WFElt(i,j,1))+
     &	      WFElt(i,j,2)*DCONJG(WFElt(i,j,2))+
     &	      WFElt(i,j,3)*DCONJG(WFElt(i,j,3))
	      IF (RS.LT.1E-24) THEN
	        MWFFT(i,j)=-24
	      ELSE
	        MWFFT(i,j)=LOG10(RS)
	      END IF
	      power=power+RS
	      IF (RS.GT.MaxInt) THEN
	        imax=i
	        jmax=j
	        MaxInt=RS
	      END IF
 3	    CONTINUE
 4	  CONTINUE
	ELSE
	  DO 6 j=1,mdttl
	    DO 5 i=1,mdttl
	      RS=WFElt(i,j,1)*DCONJG(WFElt(i,j,1))
	      IF (RS.LT.1E-24) THEN
	        MWFFT(i,j)=-24
	      ELSE
	        MWFFT(i,j)=LOG10(RS)
	      END IF
	      power=power+RS
	      IF (RS.GT.MaxInt) THEN
	        imax=i
	        jmax=j
	        MaxInt=RS
	      END IF
 5	    CONTINUE
 6	  CONTINUE
	END IF
	RETURN
	END SUBROUTINE Ca2Log

C***********************************************************************

	SUBROUTINE Ca2Gain(ifVecDif3,iElt,power,MaxInt,imax,jmax)
	use param_mod
	use elt_mod

	Implicit None

	LOGICAL ifVecDif3
	INTEGER i,j,iElt,iEm1,imax,jmax
	SREAL X0PFP,dangle
	REAL*8 DA,A,power,MaxInt,RS
	power=0e0
	MaxInt=-1d33
	IF (ifVecDif3) THEN
	  DO 4 j=1,mdttl
	    DO 3 i=1,mdttl
	      RS=WFElt(i,j,1)*DCONJG(WFElt(i,j,1))+
     &	      WFElt(i,j,2)*DCONJG(WFElt(i,j,2))+
     &	      WFElt(i,j,3)*DCONJG(WFElt(i,j,3))
	      MWFFT(i,j)=RS
	      power=power+RS
	      IF (RS.GT.MaxInt) THEN
	        imax=i
	        jmax=j
	        MaxInt=RS
	      END IF
 3	    CONTINUE
 4	  CONTINUE
	ELSE
	  DO 6 j=1,mdttl
	    DO 5 i=1,mdttl
	      RS=WFElt(i,j,1)*DCONJG(WFElt(i,j,1))
	      MWFFT(i,j)=RS
	      power=power+RS
	      IF (RS.GT.MaxInt) THEN
	        imax=i
	        jmax=j
	        MaxInt=RS
	      END IF
 5	    CONTINUE
 6	  CONTINUE
	END IF
	iEm1=iElt-1
	DA=dxElt(iElt)*dxElt(iElt)
	A=12.56637061435917d0*zElt(iEm1)*zElt(iEm1)
	RS=A/power/DA
	DO 2 j=1,mdttl
	  DO 1 i=1,mdttl
	    DA=MWFFT(i,j)
	    IF (DA.LT.1E-36) THEN
	      DA=1e-36
	    END IF
	    MWFFT(i,j)=10e0*LOG10(DA*RS)
 1	  CONTINUE
 2	CONTINUE
	dangle=57.29577951308232*dxElt(iElt)/zElt(iEm1)
	X0PFP=-dangle*FLOAT((mdttl)/2)
	DO 7 i=1,mdttl
	  XFP(i)=X0PFP
	  YFP(i)=X0PFP
	  X0PFP=X0PFP+dangle
 7	CONTINUE
	RETURN
	END SUBROUTINE Ca2Gain

C***********************************************************************

	SUBROUTINE SrfOut(iElt,iPloTyp,CTEXT,CXLABL,CYLABL,CTAG,np,na,
     &			  R1,R2,StretchType)
	use param_mod
	use cfiles_mod
	use elt_mod

	Implicit none

	INTEGER i,j,k,iPloTyp,ICLEN,np,na,StretchType,iElt,ICLEN2
	CHARACTER*(*) CXLABL,CYLABL,CTEXT,CTAG
	CHARACTER*72 CTITLE
	CHARACTER*(8) CEltNum
	SREAL R1(1),R2(1)
	CALL IntToChar(CEltNum,iElt,i)
	j=ICLEN(filnam)
	k=ICLEN2(CTEXT)
	IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &       .OR.(iPloTyp.EQ.9)) THEN
	  imgfil=filnam(1:j)//CTAG//CEltNum(1:i)
	  CTITLE=CTEXT(1:k)//', Elt='//CEltNum(1:i)
	  CALL BINOUT(iPloTyp,XFP,YFP,MWFFT,np,np,
     &	              na,na,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                dxElt(iElt),dxElt(iElt), iElt)
	ELSE IF (iPloTyp.NE.0) THEN
	  ! Screen output
	  CTITLE=CTEXT(1:k)//', Elt='//CEltNum(1:i)//', File='//filnam
	  CALL SRFPLT(iPloTyp,XFP,YFP,MWFFT,
     &	    np,np,na,na,CXLABL,CYLABL,CTITLE,R1,R2,StretchType)
C    &	    np,np,na,na,'X-Axis','Y-Axis',CTITLE,R1,R2,StretchType)
	END IF
	RETURN
	END SUBROUTINE SrfOut

C***********************************************************************

	SUBROUTINE PixOut(PixArray,iElt,iPloTyp,CTEXT,CXLABL,CYLABL,
     &	CTAG,np,na,R1,R2,StretchType)
	use param_mod
	use cfiles_mod
	use elt_mod

	Implicit None

	INTEGER i,j,k,iPloTyp,iElt,np,na,StretchType,ICLEN,ICLEN2
	CHARACTER*(*) CXLABL,CYLABL,CTEXT,CTAG
	CHARACTER*72 CTITLE
	CHARACTER*(8) CEltNum
	SREAL R1(1),R2(1),PixArray(na,na)

	CALL IntToChar(CEltNum,iElt,i)

	j=ICLEN(filnam)
	k=ICLEN2(CTEXT)
	IF ((iPloTyp.EQ.5).OR.(iPloTyp.EQ.6).OR.(iPloTyp.EQ.8)
     &       .OR.(iPloTyp.EQ.9)) THEN
	  imgfil=filnam(1:j)//CTAG//CEltNum(1:i)
	  CTITLE=CTEXT(1:k)//', Elt='//CEltNum(1:i)
	  CALL BINOUT(iPloTyp,XFP,YFP,PixArray,np,np,
     &	              na,na,'X-Axis','Y-Axis',CTITLE,imgfil,
     &                dxElt(iElt), dxElt(iElt), iElt)
	ELSE IF (iPloTyp.NE.0) THEN
	  CTITLE=CTEXT(1:k)//', Elt='//CEltNum(1:i)//', File='//filnam
	  CALL SRFPLT(iPloTyp,XFP,YFP,PixArray,
     &	  np,np,na,na,'X-Axis','Y-Axis',CTITLE,R1,R2,StretchType)
	END IF
	RETURN
	END SUBROUTINE PixOut


C***********************************************************************

        SUBROUTINE DXCALC(iElt,nEstart,dxMin,dxMax,dxAvg,dxDev,
     &                    ifRegrid,ifOPDCalc,npts)
        use param_mod
        use elt_mod
	use math_mod

        IMPLICIT NONE

C       These variables are the COMP global variables and reside in
C       common EltDble, EltReal, etc.

        REAL*8 dvec(3),dxMin,dx(mpts),dxAvg,dxDev,dxMax,dxRat,dxRef
        INTEGER i,j,k,l,iElt,iEm1,nEstart,nd2,npts,ii,iimax,lmax
        LOGICAL LFFProp,LSSProp,LFRProp,ifRegrid(mElt),ifOPDCalc

C       First see if dx is set by a far field or sphere-to-
C       sphere propagation

        iEm1=iElt-1
        LFFProp=.FALSE.
        LFRProp=.FALSE.
        LSSProp=.FALSE.
        IF (iElt.GT.0) THEN
          IF ((.NOT.ifOPDCalc).AND.((PropType(iEm1).EQ.3).OR.
     &    (PropType(iEm1).EQ.10).OR.(PropType(iEm1).EQ.15)))
     &      LFFProp=.TRUE.
          IF ((.NOT.ifOPDCalc).AND.ifRegrid(iElt).AND.
     &    ((PropType(iEm1).EQ.5).OR.(PropType(iEm1).EQ.14)))
     &      LSSProp=.TRUE.
          IF ((.NOT.ifOPDCalc).AND.((PropType(iEm1).EQ.12)))
     &      LFRProp=.TRUE.
        END IF

C       If far field, compute dx on that basis

        IF (LFFProp) THEN
          dxElt(iElt)=WaveBU*DBLE(zElt(nEstart)/dxElt(nEstart))/
     &    	      DBLE(mdttl)
          !print*,'***LFFProp: dxElt(nEstart) =',dxElt(nEstart)
          !print*,'***LFFProp: dxElt(iElt) =',dxElt(iElt)
	  ! For far field diffraction: how is grid-spacing on the image
	  ! plane computed as a function of grid-spacing on pupil plane?
	  !  Let grid-spacing on pupil plane dx1 = dxElt(nEstart),
          !  FFT variable on image plane f = u/(lamda*zElt(nEstart)),
          !  where u is the spatial variable on image plane.
	  !  According to FFT theory, ith freq f_i = i/(mdttl*dx1), so
	  !  df = 1/(mdttl*dx1) = 1/(mdttl*dxElt(nEstart)), and
          !  dxElt(iElt)= dx2 = du = df*(lamda*zElt(nEstart))
          !             = lamda*zElt(nEstart)/(mdttl*dxElt(nEstart))
	  !  The last equation is the code statement above.

	  ! Test only  -jzlou
	  !dxElt(iElt)=2d0*dxElt(iElt)

	  ! When dxElt(nEstart) is too large (and thus dx2 too small)
          ! compared to the 'needed spectrum band' on image plane,
          ! fft wrapping happens.

        ELSE IF (LSSProp) THEN
          dxElt(iElt)=ABS(dxElt(nEstart)*zElt(iElt)/zElt(nEstart))
          dxMin=0
          dxMax=0
          dxAvg=0
          dxDev=0

        ELSE IF (LFRProp) THEN
          dxElt(iElt)=WaveBU*
     &    	DBLE((zElt(nEstart)-zelt(ielt))/dxElt(nEstart)) /
     &    	DBLE(mdttl)

C        Otherwise, take a slice across the center of the pupil and
C        calculate the ray spacings along the slice.  Spacings are
C        linear rather than arc lengths!
        ELSE IF (nRay.NE.1) THEN
          nd2 = 1 + nRay/2 + npts/4
          j = nd2-npts/2
          IF (j.LT.2) j=2
          IF (nd2.GT.nRay) nd2=nRay
          dxMin=1d10; dxMax=0d0; dxAvg=0d0
          l=0
          DO i=j,(nd2-2)
            IF (LRayOK(i).AND.LRayOK(i+1)) THEN
              l=l+1
              Do k=1,3
                dvec(k)=RayPos(k,i)-RayPos(k,i+1)
              End Do
              dx(l)=DSQRT(dvec(1)*dvec(1)+dvec(2)*dvec(2)+
     &              dvec(3)*dvec(3))
              dxAvg=dxAvg+dx(l)
              IF (dx(l).LT.dxMin) dxMin=dx(l)
              IF (dx(l).GT.dxMax) dxMax=dx(l)
            END IF
          END DO

	  !print*,'Begin*** l =', l

	  !If (l==0) Then
	  If (l<4) Then
	    ! Got here because too few valid rays across center of pupil,
	    ! so try to search all line slices to find the slice with the
            ! most valid rays. -jzlou

	    lmax=0
	    Do ii=1,npts
              j=(ii-1)*npts+1; nd2=j+npts-1
	      l=0
	      DO i=j,nd2
	        IF (LRayOK(i).AND.LRayOK(i+1)) l=l+1
	      END DO
	      If (l>lmax) Then
	        lmax=l; iimax=ii
	      End If
	    END DO

	    j=(iimax-1)*npts+1; nd2=j+npts-1
            dxMin=1d10; dxMax=0d0; dxAvg=0d0
	    l=0
	    DO i=j,nd2
	      IF (LRayOK(i).AND.LRayOK(i+1)) THEN
                l=l+1
                Do k=1,3
                  dvec(k)=RayPos(k,i)-RayPos(k,i+1)
                End Do
	        dx(l)=DSQRT(dvec(1)*dvec(1)+dvec(2)*dvec(2)+
     &                  dvec(3)*dvec(3))
                dxAvg=dxAvg+dx(l)
                IF (dx(l).LT.dxMin) dxMin=dx(l)
                IF (dx(l).GT.dxMax) dxMax=dx(l)
	      END IF
	    End Do ! loop i
	    !if (l.GE.10) exit
	  End If

	  !print*,'** l =', l
	  !print*,'** dx(l) =',dx(l)

	  If (l==0) Then
            WRITE(*,*)
     &    '**Warning: Computing pupil diffraction grid spacing failed!'
            WRITE(*,*)'           Diffraction result is incorrect.'
          End If

          IF (l.GT.2) THEN
            !dxRef=(dxAvg-dxMin-dxMax)/DBLE(l-2)
	    dxRef = med(dx,l)
          ELSE
            dxRef=dx(l)
          END IF


C       Find the average and RMS deviation of the spacings; set dx to
C       the average value

          dxAvg=0d0
          j=0
          DO 3 i=1,l
            dxRat=dx(i)/dxRef  ! dx(i) divided by average dx
	    !print*,'** dx(i) =',dx(i)
            IF ((dxRat.GE.0.5).AND.(dxRat.LE.2)) THEN
	      ! above 'if' makes sure dx(i) used is not too far away
	      ! from average dx(i) -jzlou
              j=j+1
              dxAvg=dxAvg+dx(i)
            END IF
 3        CONTINUE

          IF (j.EQ.0) THEN
            dxElt(iElt)=1d10
	    !print*,'*** dxElt = 1d10'
	    !stop
            dxAvg=0d0; dxMax=0d0; dxDev=0d0
          ELSE
            dxAvg=dxAvg/DBLE(j)
            dxDev=0d0
            DO 5 i=1,l
              dxRat=dx(i)/dxRef
              IF ((dxRat.GE.0.5).AND.(dxRat.LE.2)) THEN
                dxDev=dxDev+(dx(i)-dxAvg)**2
              END IF
 5          CONTINUE
            dxDev=DSQRT(dxDev/DBLE(j))
            dxElt(iElt)=dxAvg
          END IF
        ELSE
          dxElt(iElt)=1d10
        END IF

        RETURN
        END SUBROUTINE DXCALC

C***********************************************************************

	Subroutine SetHexBound(h, thisElt)
	use param_mod
        use elt_mod

	Implicit none

        Real*8 :: h, a
	Integer :: thisElt
	Real*8, parameter :: eps=1.d-18

	if ((2*h)<=eps) then
	  print*,' '
          print*,'** Invalid hexagon aperture size, width=', 2*h,
     &           ' , quit ...'
	  stop
        end if

	a=h/sqrt(3.d0)

        if (.true.) then
	  ! sharp top hexagon
	  ApHexBound(1,1,thisElt)=h
	  ApHexBound(1,2,thisElt)=(2*a-a)/(0.0-h)
          ApHexBound(2,2,thisElt)=h; ApHexBound(3,2,thisElt)=a;
          ApHexBound(1,3,thisElt)=(a-2*a)/(-h-0.0)
          ApHexBound(2,3,thisElt)=0; ApHexBound(3,3,thisElt)=2*a
          ApHexBound(1,4,thisElt)=-h
          ApHexBound(1,5,thisElt)=(-2*a+a)/(0+h)
          ApHexBound(2,5,thisElt)=-h; ApHexBound(3,5,thisElt)=-a
          ApHexBound(1,6,thisElt)=(-a+2*a)/(h-0)
	  ApHexBound(2,6,thisElt)=0; ApHexBound(3,6,thisElt)=-2*a
        else
	  ! flat top hexagon
	  ApHexBound(1,1,thisElt)=(h-0.0)/(a-2*a)
	  ApHexBound(2,1,thisElt)=2*a; ApHexBound(3,1,thisElt)=0
	  ApHexBound(1,2,thisElt)=h
	  ApHexBound(1,3,thisElt)=(0.0-h)/(-2*a+a)
          ApHexBound(2,3,thisElt)=-a; ApHexBound(3,3,thisElt)=h
	  ApHexBound(1,4,thisElt)=(-h-0.0)/(-a+2*a)
          ApHexBound(2,4,thisElt)=-2*a; ApHexBound(3,4,thisElt)=0
          ApHexBound(1,5,thisElt)=-h
	  ApHexBound(1,6,thisElt)=(0.0+h)/(2*a-a)
	  ApHexBound(2,6,thisElt)=a; ApHexBound(3,6,thisElt)=-h
	end if
	End Subroutine SetHexBound


C***********************************************************************

        Subroutine SetCvxPolyApVtx(thisElt,nVertex,Vtx3D)
        use param_mod
        use elt_mod
        use math_mod

        Implicit none

        Integer :: thisElt,nVertex,j
        REAL*8 :: Vtx3D(3,nVertex)
        REAL*8 :: off(3),xObs_loc(3),yObs_loc(3),zObs_loc(3)

        CALL DEQUATE(zObs_loc,psiElt(1:3,thisElt),3)
        CALL DXPROD(yObs_loc,zObs_loc,xObs(1,thisElt))
        CALL DUNITIZE(yObs_loc)
        CALL DXPROD(xObs_loc,yObs_loc,zObs_loc)
        CALL DUNITIZE(xObs_loc)

	nVertex=ApVec(1,thisElt); ApVec(1:3,thisElt)=0d0
	Do j=1,nVertex
	  ApVec(1:3,thisElt)=ApVec(1:3,thisElt)+Vtx3D(1:3,j)
	End Do
	ApVec(1:3,thisElt)=ApVec(1:3,thisElt)/nVertex

        !nVertex=ApVec(4,thisElt)
        off=ApVec(1:3,thisElt)-VptElt(1:3,thisElt)
        ApVec(1,thisElt)=DDOTC(off,xObs_loc)
        ApVec(2,thisElt)=DDOTC(off,yObs_loc)
        ApVec(3,thisElt)=nVertex
        !
        Do j=1,nVertex
          off(1:3)=Vtx3D(1:3,j)-VptElt(1:3,thisElt)
          PolyApVtx(1,j,thisElt)=DDOTC(off,xObs_loc)
          PolyApVtx(2,j,thisElt)=DDOTC(off,yObs_loc)
          PolyApVtx(1,j,thisElt)=PolyApVtx(1,j,thisElt)-ApVec(1,thisElt)
          PolyApVtx(2,j,thisElt)=PolyApVtx(2,j,thisElt)-ApVec(2,thisElt)
        End Do
        End Subroutine SetCvxPolyApVtx

C***********************************************************************

        Subroutine SetCvxPolyObsVtx(thisElt,nVertex,Vtx3D,iObs)
        use param_mod
        use elt_mod
        use math_mod

        Implicit none

        Integer :: thisElt,nVertex,iObs,j
        REAL*8 :: Vtx3D(3,nVertex)
        REAL*8 :: off(3),xObs_loc(3),yObs_loc(3),zObs_loc(3)

        CALL DEQUATE(zObs_loc,psiElt(1:3,thisElt),3)
        CALL DXPROD(yObs_loc,zObs_loc,xObs(1,thisElt))
        CALL DUNITIZE(yObs_loc)
        CALL DXPROD(xObs_loc,yObs_loc,zObs_loc)
        CALL DUNITIZE(xObs_loc)

	nVertex=ObsVec(1,iObs,thisElt)
	ObsVec(1:3,iObs,thisElt)=0d0
	Do j=1,nVertex
	  ObsVec(1:3,iObs,thisElt)
     &       =ObsVec(1:3,iObs,thisElt)+Vtx3D(1:3,j)
	End Do
	ObsVec(1:3,iObs,thisElt)=ObsVec(1:3,iObs,thisElt)/nVertex

        !nVertex=ObsVec(4,iObs,thisElt)
        off=ObsVec(1:3,iObs,thisElt)-VptElt(1:3,thisElt)
        ObsVec(1,iObs,thisElt)=DDOTC(off,xObs_loc)
        ObsVec(2,iObs,thisElt)=DDOTC(off,yObs_loc)
        ObsVec(3,iObs,thisElt)=nVertex
        !
        Do j=1,nVertex
          off(1:3)=Vtx3D(1:3,j)-VptElt(1:3,thisElt)
          PolyObsVtx(1,j,iObs,thisElt)=DDOTC(off,xObs_loc)
          PolyObsVtx(2,j,iObs,thisElt)=DDOTC(off,yObs_loc)
          PolyObsVtx(1,j,iObs,thisElt)=PolyObsVtx(1,j,iObs,thisElt)
     &      -ObsVec(1,iObs,thisElt)
          PolyObsVtx(2,j,iObs,thisElt)=PolyObsVtx(2,j,iObs,thisElt)
     &      -ObsVec(2,iObs,thisElt)
        End Do
        End Subroutine SetCvxPolyObsVtx


C***********************************************************************


	Subroutine SetCvxPolyApBound(thisElt)
	use param_mod
	use elt_mod

	Implicit none

	Integer :: thisElt,nside,iside,idx1,idx2
	Real*8, parameter :: eps=1.d-8
	Real*8 :: pt1(2),pt2(2),det

	! Use PolyApVtx(:,:,:) to define PolyApEqn(3,MaxPolySide,mElt)
	! For each side, define a line equation: A*x+B*x+C=0

	nside= ApVec(3,thisElt) ! # vertices = # sides for N>=3

	Do iside=1,nside
	  idx1=iside
	  if (iside==nside) then
	    idx2=1
	  else
	    idx2=iside+1
	  end if
	  pt1(1:2)=PolyApVtx(1:2,idx1,thisElt)
	  pt2(1:2)=PolyApVtx(1:2,idx2,thisElt)

	  if (abs(pt2(2)-pt1(2))<eps) then
	    PolyApEqn(1,iside,thisElt)=0.d0
	    PolyApEqn(2,iside,thisElt)=1.d0
	    PolyApEqn(3,iside,thisElt)=-pt1(2)
	  else if (abs(pt2(1)-pt1(1))<eps) then
	    PolyApEqn(1,iside,thisElt)=1.d0
            PolyApEqn(2,iside,thisElt)=0.d0
            PolyApEqn(3,iside,thisElt)=-pt1(1)
	  else
	    det=(pt2(2)-pt1(2))/(pt2(1)-pt1(1))
	    PolyApEqn(1,iside,thisElt)=-det
            PolyApEqn(2,iside,thisElt)=1.d0
            PolyApEqn(3,iside,thisElt)=pt1(1)*det-pt1(2)
	    if (.false. .and. iside==2 .and. thisElt==7) then
              print*,'pt1(1:2) = ', pt1(1:2)
              print*,'pt2(1:2) = ', pt2(1:2)
              print*,'det = ', det
              print*,'PolyApEqn(1:3,iside,thisElt)=',
     &          PolyApEqn(1:3,iside,thisElt)
            end if
	  end if
	End Do

        if (.false. .and. thisElt==18) then
	  print*,'In SetCvxPolyBound: '
	  print*,'PolyApEqn(1:3,1,thisElt)=',PolyApEqn(1:3,1,thisElt)
	  print*,'PolyApEqn(1:3,2,thisElt)=',PolyApEqn(1:3,2,thisElt)
	  print*,'PolyApEqn(1:3,3,thisElt)=',PolyApEqn(1:3,3,thisElt)
	  print*,' '
	  stop
	end if
	End Subroutine SetCvxPolyApBound


C***********************************************************************

	Subroutine SetCvxPolyObsBound(iObs,thisElt)
	use param_mod
        use elt_mod

        Implicit none

        Integer :: iObs,thisElt,nside,iside,idx1,idx2
        Real*8, parameter :: eps=1d-8
        Real*8 :: pt1(2),pt2(2),det

        ! Use PolyObsVtx(:,:,:,:) to define PolyObsEqn(3,MaxPolySide,iObs,mElt)
        ! For each side, define a line equation: A*x+B*x+C=0

        nside= ObsVec(3,iObs,thisElt) ! # vertices = # sides for N>=3

	Do iside=1,nside
          idx1=iside
          if (iside==nside) then
            idx2=1
          else
            idx2=iside+1
          end if
          pt1(1:2)=PolyObsVtx(1:2,idx1,iObs,thisElt)
          pt2(1:2)=PolyObsVtx(1:2,idx2,iObs,thisElt)

          if (abs(pt2(2)-pt1(2))<eps) then
            PolyObsEqn(1,iside,iObs,thisElt)=0.d0
            PolyObsEqn(2,iside,iObs,thisElt)=1.d0
            PolyObsEqn(3,iside,iObs,thisElt)=-pt1(2)
          else if (abs(pt2(1)-pt1(1))<eps) then
            PolyObsEqn(1,iside,iObs,thisElt)=1.d0
            PolyObsEqn(2,iside,iObs,thisElt)=0.d0
            PolyObsEqn(3,iside,iObs,thisElt)=-pt1(1)
          else
            det=(pt2(2)-pt1(2))/(pt2(1)-pt1(1))
            PolyObsEqn(1,iside,iObs,thisElt)=-det
            PolyObsEqn(2,iside,iObs,thisElt)=1.d0
            PolyObsEqn(3,iside,iObs,thisElt)=pt1(1)*det-pt1(2)

	    if (.false. .and. iside==2) then
	      print*,'pt1(1:2) = ', pt1(1:2)
	      print*,'pt2(1:2) = ', pt2(1:2)
	      print*,'det = ', det
	      print*,'PolyObsEqn(1:3,iside,iObs)=',
     &          PolyObsEqn(1:3,iside,iObs,thisElt)
	    end if
          end if
	End Do

	if (.false. .and. (thisElt==18 .or. thisElt==26)) then
	 print*,'In SetCvxPolyObsBound: thisElt =',thisElt
         Do iside=1,nside
	   print*,'PolyObsEqn(1:3,iside,iObs)=',
     &	 	  PolyObsEqn(1:3,iside,iObs,thisElt)
	 End Do
         !stop
	end if
	End Subroutine SetCvxPolyObsBound


C********************************************************************

	Subroutine FSR(npts,lastElt,Dicr,Djcr,dxSource,ifPol,
     &                 ifPropagate,rxp,ryp,rid,ifPrt)
! rxp,ryp = input x,y position in source aperture (in terms of
!           percent of aperture); ray coord is defined by (xGrid,yGrid)
! rid = ray ID that matches most closely to the input rxp,ryp
        use param_mod
        use src_mod
        use elt_mod
        use sourcsub_mod
	use math_mod
 	use macos_mod, only : ifBuild,ifDMBuild,ifDraw,nDrawElt,
     &        	ifLNegOK,Cmatrix,ifRayPrint,ifIntSrfSet,ifTrace
#ifdef CMACOS
	use macos_vars_mod, only : xDraw,yDraw,DrawEltVec,
     &                      DrawRayVec,StartExtinc,StartIndRef,
     &                      nBadRays
#else
	use smacos_vars_mod, only : xDraw,yDraw,DrawEltVec,
     &                      DrawRayVec,StartExtinc,StartIndRef,
     &                      nBadRays
#endif
	use tracesub_mod

        Implicit none

        Integer :: npts,lastElt,rid
        Real*8 :: Dicr,Djcr,dxSource, rxp,ryp
        Logical :: ifPol,ifPropagate,ifPrt

        Integer :: iRay,iStartRay,iEndRay,iElt,iStartElt,iEndElt,ii
        Real*8 :: rx,ry,rx_max,rx_min,ry_max,ry_min,xD,yD,xray_p,yray_p,
     &            rd,tmp_d,dRay(3),dRaySave(3)
        Real*8 :: ray_d,ray_apd ! for test only

        ! Calc source ray grid
        IF (ABS(zSource).GE.1d10) THEN
          ! source at infinity - jzlou
          CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                   ifPropagate)
        ELSE
          CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,
     &                  ifPropagate)
        END IF

        CALL DUNITIZE(xGrid)
        CALL DUNITIZE(yGrid)

        ! Loop through nRay source rays to find source aperture diameter
	rx_max=-1d10; rx_min=1d10
        ry_max=-1d10; ry_min=1d10
	ray_apd=-1d0
        Do iRay=1,nRay
	  dRay(1:3)=RayPos(1:3,iRay)-RayPos(1:3,1)
          rx=DDOTC(dRay(1:3),xGrid(1:3))
          ry=DDOTC(dRay(1:3),yGrid(1:3))
          if (rx>rx_max) rx_max=rx
          if (rx<rx_min) rx_min=rx
          if (ry>ry_max) ry_max=ry
          if (ry<ry_min) ry_min=ry
	  ray_d=sqrt(rx**2+ry**2)
	  if (ray_d>ray_apd) ray_apd=ray_d
        End Do

	! X and y spans of source aperture
	xD=rx_max; yD=ry_max

        ! Ideal (target) source ray position
        xray_p=rxp*xD; yray_p=ryp*yD

        rid=-1

	! Initial distance of difference rays
        rd=1d14

        ! Loop through nRay to find a source ray that best fits (rxp,ryp)
        ifRayPrint=.FALSE.
	iStartElt=0; iElt=lastElt
        Do iRay=2,nRay
	  dRay(1:3)=RayPos(1:3,iRay)-RayPos(1:3,1)
          rx=DDOTC(dRay(1:3),xGrid(1:3))
          ry=DDOTC(dRay(1:3),yGrid(1:3))
	  !if (abs(rx) < 350.d0 .and. abs(ry) < 15.*dxSource) then
	  !  print*,'ray id, ry =', iRay, ry
	  !end if
          tmp_d=sqrt((rx-xray_p)**2+(ry-yray_p)**2)
          If (tmp_d<rd) Then
	    ! Make sure the ray can go thru system
	    iStartRay=iRay; iEndRay=iRay
	    ifTrace=.FALSE.; ifBuild=.FALSE.; ifPropagate=.FALSE.
            CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &        ifLNegOK,Cmatrix,iStartElt,iElt,iStartRay,iEndRay,
     &        RayMap,SegMap,ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,
     &        StartIndRef,StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &        xDraw,yDraw)
	    !print*,'** iRay,LRayOK(iRay) =',iRay,LRayOK(iRay); pause
	    if (LRayOK(iRay).AND.LRayPass(iRay)) then
              rd=tmp_d; rid=iRay; dRaySave=dRay
	    end if
          End If
        End Do
	ifRayPrint=.TRUE.

	if (rid>0) then
          rx=DDOTC(dRaySave(1:3),xGrid(1:3))
          ry=DDOTC(dRaySave(1:3),yGrid(1:3))
	  if (ifPrt) write(*,150) rx,ry,rid
	else
	  write(*,*) '** Failed to find current reference ray!'
	end if
 150    format('   Ref. ray XY offsets from chief ray:',f12.4,',',f12.4,
     &         ',   ID =',I8)
	rxp=rx; ryp=ry   ! return coord positions of ray, -jzlou
        End Subroutine FSR

C***********************************************************************

        Subroutine FindMarginalRay(npts,lastElt,Dicr,Djcr,dxSource,
     &                  ifPol,ifPropagate,rxp,ryp,rid,PrtFlg)
! rxp,ryp = ray position in source aperture of found marginal ray,
!           in source frame; return value
! rid = ray ID of found marginal ray, return value
        use param_mod
        use src_mod
        use elt_mod
        use sourcsub_mod
        use math_mod
        use macos_mod, only : ifBuild,ifDMBuild,ifDraw,nDrawElt,
     &          ifLNegOK,Cmatrix,ifRayPrint,ifIntSrfSet,ifTrace
#ifdef CMACOS
        use macos_vars_mod, only : xDraw,yDraw,DrawEltVec,
     &                      DrawRayVec,StartExtinc,StartIndRef,
     &                      nBadRays
#else
        use smacos_vars_mod, only : xDraw,yDraw,DrawEltVec,
     &                      DrawRayVec,StartExtinc,StartIndRef,
     &                      nBadRays
#endif
        use tracesub_mod

        Implicit none

	Integer :: npts,lastElt,rid
        Real*8 :: Dicr,Djcr,dxSource,rxp,ryp
        Logical :: ifPol,ifPropagate,PrtFlg

        Integer :: iRay,iStartRay,iEndRay,iElt,iStartElt
        Real*8 :: rx,ry,rx_max,rx_min,ry_max,ry_min,xD,yD,xray_p,yray_p,
     &            rd,tmp_d,dRay(3),dRaySave(3)
        Real*8 :: ray_d,ray_apd ! for test only

        ! Calc source ray grid
        IF (ABS(zSource).GE.1d10) THEN
          ! source at infinity - jzlou
          CALL ColSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
        ELSE
          CALL PtSource(npts,Dicr,Djcr,dxSource,ifPol,ifPropagate)
        END IF

        CALL DUNITIZE(xGrid)
        CALL DUNITIZE(yGrid)

	! Loop through nRay source rays to find source aperture diameter
        rx_max=-1d10; rx_min=1d10
        ry_max=-1d10; ry_min=1d10
        ray_apd=-1d0
        Do iRay=1,nRay
          dRay(1:3)=RayPos(1:3,iRay)-RayPos(1:3,1)
          rx=DDOTC(dRay(1:3),xGrid(1:3))
          ry=DDOTC(dRay(1:3),yGrid(1:3))
          if (rx>rx_max) rx_max=rx
          if (rx<rx_min) rx_min=rx
          if (ry>ry_max) ry_max=ry
          if (ry<ry_min) ry_min=ry
          ray_d=sqrt(rx**2+ry**2)
          if (ray_d>ray_apd) ray_apd=ray_d
        End Do

        ! X and y spans of source aperture
        xD=rx_max; yD=ry_max

        ! Central source ray position - most likely chief ray
        xray_p=0.5*xD; yray_p=0.5*yD

        rid=-1

        ! Initial distance of difference rays
        rd=0d0

	! Loop through nRay to find a source ray that maximize distance
	! from central ray
        ifRayPrint=.FALSE.
        iStartElt=0; iElt=nElt
        Do iRay=2,nRay
          dRay(1:3)=RayPos(1:3,iRay)-RayPos(1:3,1)
          rx=DDOTC(dRay(1:3),xGrid(1:3))
          ry=DDOTC(dRay(1:3),yGrid(1:3))
          tmp_d=sqrt((rx-xray_p)**2+(ry-yray_p)**2)
          If (tmp_d>rd) Then
            ! Make sure the the ray can go thru system
            iStartRay=iRay; iEndRay=iRay
            ifTrace=.FALSE.
            CALL CTRACE(RayIndex,ifBuild,ifDMBuild,ifPol,ifDraw,
     &        ifLNegOK,Cmatrix,
     &        iStartElt,iElt,iStartRay,iEndRay,RayMap,SegMap,
     &        ifRayPrint,ifIntSrfSet,LRayPass,nBadRays,StartIndRef,
     &        StartExtinc,DrawRayVec,DrawEltVec,nDrawElt(1),
     &        xDraw,yDraw)
            if (LRayOK(iRay).AND.LRayPass(iRay)) then
              rd=tmp_d; rid=iRay; dRaySave=dRay
            end if
          End If
        End Do
        ifRayPrint=.TRUE.

        if (rid>0) then
          !dRay(1:3)=RayPos(1:3,rid)-RayPos(1:3,1)
          rx=DDOTC(dRaySave(1:3),xGrid(1:3))
          ry=DDOTC(dRaySave(1:3),yGrid(1:3))
          if (PrtFlg) write(*,150) rx,ry,rid
        else
          write(*,*) '** Failed to find current referency ray!'
        end if
 150    format('   Ref. ray XY offsets from chief ray:',f12.4,',',f12.4,
     &         ',   ID =',I8)
        rxp=rx; ryp=ry   ! return coord positions of ray, -jzlou
	End Subroutine FindMarginalRay

C***********************************************************************

C Find system magnification factor M
	Function FindSysMag(npts,Dicr,Djcr,dxSource,
     & 		            ifPol,ifPropagate,ifPrt)
	use param_mod
        use src_mod
        use elt_mod
        use sourcsub_mod
        use math_mod
        use macos_mod, only : ifBuild,ifDMBuild,ifDraw,
     &       nDrawElt,
     &       ifLNegOK,Cmatrix,ifRayPrint,ifIntSrfSet,ifTrace
#ifdef CMACOS
        use macos_vars_mod, only : xDraw,yDraw,DrawEltVec,
     &                   DrawRayVec,StartExtinc,StartIndRef,
     &                   nBadRays
#else
        use smacos_vars_mod, only : xDraw,yDraw,DrawEltVec,
     &                   DrawRayVec,StartExtinc,StartIndRef,
     &                   nBadRays
#endif
        use tracesub_mod

	Implicit none

	Integer :: npts
	Real*8 :: Dicr,Djcr,dxSource,FindSysMag
	Logical :: ifPol,ifPropagate,ifPrt
        FindSysMag = 1.0  ! placeholder
	End Function FindSysMag

C***********************************************************************

	Subroutine CalcMTF(WF,MTF2,ndim)
	use math_mod

	Implicit none

	Complex*16, dimension(1:ndim,1:ndim) :: WF
	Real*8, dimension(1:2*ndim,1:2*ndim) :: MTF2
        Complex*16, dimension(1:ndim,1:ndim) :: tmp_carr
	Complex*16 :: csum
	Real*8 :: cv
	Integer :: ndim,dimi,dimj,i,j

	Do i=0,ndim-1
	  Do j=0,ndim-1
	     dimi=ndim-i; dimj=ndim-j
	     tmp_carr(1:dimi,1:dimj)
     &         =WF(i+1:ndim,j+1:ndim)*WF(1:ndim-i,1:ndim-j)
	     Call CMATSUM(tmp_carr(1:dimi,1:dimj),
     &                    csum,dimi,dimj)
             cv=csum*DCONJG(csum)
	     MTF2(ndim+i,ndim+j)=sqrt(cv)
	     !
             tmp_carr(1:dimi,1:dimj)
     &         =WF(1:ndim-i,j+1:ndim)*WF(i+1:ndim,1:ndim-j)
	     Call CMATSUM(tmp_carr(1:dimi,1:dimj),
     &                    csum,dimi,dimj)
             cv=csum*DCONJG(csum)
             MTF2(ndim-i,ndim+j)=sqrt(cv)
	     !
	     tmp_carr(1:dimi,1:dimj)
     &         =WF(i+1:ndim,1:ndim-j)*WF(1:ndim-i,j+1:ndim)
             Call CMATSUM(tmp_carr(1:dimi,1:dimj),
     &                    csum,dimi,dimj)
             cv=csum*DCONJG(csum)
             MTF2(ndim+i,ndim-j)=sqrt(cv)
	     !
	     dimi=ndim-i; dimj=ndim-j
             tmp_carr(1:dimi,1:dimj)
     &         =WF(1:ndim-i,1:ndim-j)*WF(i+1:ndim,j+1:ndim)
             Call CMATSUM(tmp_carr(1:dimi,1:dimj),
     &                    csum,dimi,dimj)
             cv=csum*DCONJG(csum)
             MTF2(ndim-i,ndim-j)=sqrt(cv)
	  End Do
        End Do

	MTF2=MTF2/MTF2(ndim,ndim) ! normalize with DC component
	End Subroutine CalcMTF


C***********************************************************************

	Subroutine Calc_EltPrRotations(old_psiElt,new_psiElt,
     &                               nECrd,TElt,elt_rots,Rtn_OK)
	use math_mod

	Implicit none

	Real*8 :: old_psiElt(3),new_psiElt(3),elt_rots(3,3),TElt(3,3),
     &            rotAxis(3),rotAng,C13,C23,C33,C12,C11,theta(3),pi,dd,
     &            tTElt(3,3),R(3)
	Logical :: Rtn_OK
	Integer :: nEcrd,i,j

	pi=4d0*atan(1d0)

	!print*,'Calc_EltRot: old_psiElt=',old_psiElt(1:3)
	!print*,'Calc_EltRot: new_psiElt=',new_psiElt(1:3)

	! Calc Eular rotation vector to go from old_psiElt to new_psiElt
	Call DXPROD(rotAxis,old_psiElt,new_psiElt)
	Call DUNITIZE(rotAxis)

	if (DDOTC(rotAxis,rotAxis)<1d-10) then
	  WRITE(*,*)
     &      '  Calib: calculated element total rotation = 0d0'
	  Rtn_OK=.false.
	  return
	else
	  Rtn_OK=.true.
	end if

	dd=DDOTC(old_psiElt,new_psiElt)
	dd=dd/(sqrt(DDOTC(old_psiElt,old_psiElt))*
     &         DDOTC(new_psiElt,new_psiElt))
	if (abs(dd)>1d0) then
	  dd=DSIGN(1d0,dd) ! make sure abs(dd) not exceeding one
	end if
	rotAxis=dacos(dd)*rotAxis

	If (nEcrd .ge. 3) Then
	  Do i=1,3
	    Do j=1,3
	      tTElt(i,j)=TElt(j,i)
	    End Do
	  End Do
	  ! Project rotAxis from global onto element coord system
	  CALL DMPROD(R,tTElt(1:3,1:3),rotAxis,3,3,1)
	  rotAxis=R
	End If

        dd=DDOTC(rotAxis,rotAxis)
	rotAng=sqrt(dd)
	Call DUNITIZE(rotAxis)

	!print*,'rotAng =',rotAng
	!print*,'rotAxis =',rotAxis

	C13=(1d0-dcos(rotAng))*rotAxis(1)*rotAxis(3)
     &      -rotAxis(2)*dsin(rotAng)
	C23=(1d0-dcos(rotAng))*rotAxis(2)*rotAxis(3)
     &      +rotAxis(1)*dsin(rotAng)
	C33=(1d0-dcos(rotAng))*rotAxis(3)**2+dcos(rotAng)
	C12=(1d0-dcos(rotAng))*rotAxis(1)*rotAxis(2)
     &      +rotAxis(3)*dsin(rotAng)
	C11=(1d0-dcos(rotAng))*rotAxis(1)**2+dcos(rotAng)

        theta(2)=-dasin(C13)
	if (abs(theta(2)-pi/2d0)<1d-5) then
	  print*,' Calc_EltRotations: singularity encountered!'
	  stop
	end if

	dd=dcos(theta(2))**2
	theta(1)=atan(C23*C33/dd)
	theta(3)=atan(C12*C11/dd)

	! Construct elt_rots
	elt_rots(1:3,1:3)=0d0
	elt_rots(1,1)=1d0; elt_rots(2,2)=1d0; elt_rots(3,3)=1d0;
	elt_rots(1:3,1)=theta(1)*elt_rots(1:3,1)
	elt_rots(1:3,2)=theta(2)*elt_rots(1:3,2)
	elt_rots(1:3,3)=theta(3)*elt_rots(1:3,3)
	End Subroutine Calc_EltPrRotations

C***********************************************************************

C Several routines supporting LensArray and Refractor with variable IndRef

	Subroutine InitLensletIndexAndCtr(type,width,iElt)
	use traceutil_mod, only : LensletStrArray,nHexLenlets,nRecLenlets,
     &                      HexLensletIdx,RecLensletIdx
	Implicit none
	Real*8, parameter :: tan30=5.773502691896257d-01
	Real*8 :: width,w2,xhat(2),yhat(2),xseg,yseg
	Integer :: type,iElt,totLensElts,i,k,ndim,nn,r,c,LINT,
     &             ILtri,IRtri,IXtri,IYtri,IXseg,IYseg
	Logical, save :: first_entry=.true.

	If (first_entry) Then
#include "lensarr_indexes.inc"
	  first_entry=.false.
	End If

	w2=0.5*width
	k=1
	If (abs(type)==1) Then
	  ! hex lenslets
	  totLensElts=nHexLenlets
	  Do i=1,totLensElts
	    LensletStrArray(i,iElt)%Xtri=HexLensletIdx(k)
    	    LensletStrArray(i,iElt)%Ltri=HexLensletIdx(k+1)
            LensletStrArray(i,iElt)%Rtri=HexLensletIdx(k+2)
	    k=k+3
	  End Do
	Else
	  totLensElts=nRecLenlets
	  ! rectangular lenslets
	  Do i=1,totLensElts
	    LensletStrArray(i,iElt)%RecXtri=RecLensletIdx(k)
            LensletStrArray(i,iElt)%RecYtri=RecLensletIdx(k+1)
	    k=k+2
	  End Do
	End If

	! Now computer XY centers of lenslets
        xhat(1)=1d0; xhat(2)=0d0
        yhat(1)=0d0; yhat(2)=1d0
	If (abs(type)==1) Then
	  ! Hex lenslet array
	  Do i=1,totLensElts
	    ILtri=LensletStrArray(i,iElt)%Ltri
	    IRtri=LensletStrArray(i,iElt)%Rtri
	    IXtri=LensletStrArray(i,iElt)%Xtri
	    IYtri=IRtri+ILtri

	    IYseg=LINT(DBLE(IYtri+2)/3d0)
            If (MOD(IYseg,2).EQ.0) Then
              IXseg=LINT(DBLE(IXtri+1)/2d0)
              xseg=IXseg*width
              yseg=1.5d0*DBLE(IYseg)*width*tan30 ! = ((sqrt(3)/2)*width)*IYseg
	    Else
              IXseg=LINT(DBLE(IXtri)/2d0)
              xseg=IXseg*width+w2
              yseg=1.5d0*DBLE(IYseg)*width*tan30
            End If
	    LensletStrArray(i,iElt)%Xctr=xseg
	    LensletStrArray(i,iElt)%Yctr=yseg
	  End Do
	Else
	  ! Rectangular lenslet array
	  Do i=1,totLensElts
	    IXtri=LensletStrArray(i,iElt)%RecXtri
	    IYtri=LensletStrArray(i,iElt)%RecYtri
	    IXseg=LINT(DBLE(IXtri+1)/2d0)
            LensletStrArray(i,iElt)%Xctr=IXseg*width
            IYseg=LINT(DBLE(IYtri+1)/2d0)
            LensletStrArray(i,iElt)%Yctr=IYseg*width
!	    print*,'RecXtri,RecYtri,Xctr,Yctr =',
!    &        IXtri,IYtri,IXseg*width,IYseg*width
	  End Do
	End If
	End Subroutine InitLensletIndexAndCtr


C   This routine is for use as a SMACOS interface to set IndRef values
C   of a lens array at a specific element, assuming counter-clockwise
C   ordering for both Hex and Rec arrays.
	Subroutine SetLensArrayIndRefVals(iElt,nval,vals)
	use traceutil_mod, only : LensletStrArray
	Integer :: iElt,nval
	Real*8 :: vals(nval)
 	LensletStrArray(1:nVal,iElt)%IndRef=vals(1:nVal)
	End Subroutine SetLensArrayIndRefVals

C

        Subroutine InitLensArrIndRefFrmXYFile(filename,iElt,type)
	use traceutil_mod, only : LensletStrArray,
     &                      mLenslet,nHexLenlets,nRecLenlets
	use math_mod

	Implicit none

	Character(len=*) :: filename
	Integer :: iElt,type,i,j,fu,nval
	Real*8 :: dd(2),dd2(2),ds,val(3)
	Real*8 :: xyzset(3,3*mLenslet)

	! Read in IndRef at XY points
	fu=10
	Open(unit=fu,file=filename)
	Read(fu,*) val(1:3)
	nval=val(1)
 	Do i=1,nval
	  Read(fu,*) xyzset(1:3,i)
	End Do
	close(fu)

	If (abs(type)==1) Then
	  ! Hex lens array
	  ! Nearest neighbor scheme
	  LensletStrArray(1:nHexLenlets,iElt)%NgbDist=1d22
	  Do j=1,nHexLenlets
	    dd(1)=LensletStrArray(j,iElt)%Xctr
	    dd(2)=LensletStrArray(j,iElt)%Yctr
	    Do i=1,nval
	      dd2(1:2)=xyzset(1:2,i)-dd(1:2)
	      ds=DDOTC2(dd2,dd2) ! diff vector from data pt to pmon
	      if (ds<LensletStrArray(j,iElt)%NgbDist) then
	        LensletStrArray(j,iElt)%NgbDist=ds
	        LensletStrArray(j,iElt)%IndRef=xyzset(3,i)
	      end if
	    End Do
	  End Do
	Else
	  ! Rec lens array
	  ! Nearest neighbor scheme
	  LensletStrArray(1:nRecLenlets,iElt)%NgbDist=1d22
          Do j=1,nRecLenlets
	    dd(1)=LensletStrArray(j,iElt)%Xctr
            dd(2)=LensletStrArray(j,iElt)%Yctr
	    Do i=1,nval
              dd2(1:2)=xyzset(1:2,i)-dd(1:2)
	      ds=DDOTC2(dd2,dd2) ! diff vector from data pt to pmon
              if (ds<LensletStrArray(j,iElt)%NgbDist) then
                LensletStrArray(j,iElt)%NgbDist=ds
                LensletStrArray(j,iElt)%IndRef=xyzset(3,i)
              end if
	    End Do
          End Do
	End If
        End Subroutine InitLensArrIndRefFrmXYFile

C

	! For matlab/mex interface
	Subroutine SetLensArrIndRefXYLoc(iElt,type,nval,xyzset)
	use elt_mod, only : IndRef
        use traceutil_mod
	use math_mod

        Implicit none

        Integer :: iElt,type,nval,i,j
        Real*8 :: xyzset(3,nval),dd(2),dd2(2),ds

	IsVarIndRefElt(iElt)=.true.
        LensletStrArray(1:mLensLet,iElt)%IndRef=IndRef(iElt) ! default

	If (abs(type)==1) Then
          ! Hex lens array
          ! Nearest neighbor scheme
          LensletStrArray(1:nHexLenlets,iElt)%NgbDist=1d22
          Do j=1,nHexLenlets
            dd(1)=LensletStrArray(j,iElt)%Xctr
            dd(2)=LensletStrArray(j,iElt)%Yctr
            Do i=1,nval
              dd2(1:2)=xyzset(1:2,i)-dd(1:2)
              ds=DDOTC2(dd2,dd2) ! diff vector from data pt to pmon
              if (ds<LensletStrArray(j,iElt)%NgbDist) then
                LensletStrArray(j,iElt)%NgbDist=ds
                LensletStrArray(j,iElt)%IndRef=xyzset(3,i)
              end if
            End Do
          End Do
        Else
          ! Rec lens array
          ! Nearest neighbor scheme
          LensletStrArray(1:nRecLenlets,iElt)%NgbDist=1d22
          Do j=1,nRecLenlets
            dd(1)=LensletStrArray(j,iElt)%Xctr
            dd(2)=LensletStrArray(j,iElt)%Yctr
            Do i=1,nval
              dd2(1:2)=xyzset(1:2,i)-dd(1:2)
              ds=DDOTC2(dd2,dd2) ! diff vector from data pt to pmon
              if (ds<LensletStrArray(j,iElt)%NgbDist) then
                LensletStrArray(j,iElt)%NgbDist=ds
                LensletStrArray(j,iElt)%IndRef=xyzset(3,i)
              end if
            End Do
          End Do
        End If
        End Subroutine SetLensArrIndRefXYLoc


C***********************************************************************

	Subroutine InitIndRefGridFrmFile(filename,iElt)
	use traceutil_mod, only : IndRefGridArr

        Implicit none

        Character(len=*) :: filename
        Integer :: iElt,fu,nx,ny,j
	Real*8 :: fnx,fny

        ! Read in IndRefGrid
        fu=12
        Open(unit=fu,file=filename)
        Read(fu,*) fnx,fny,IndRefGridArr(iElt)%width
	nx=fnx; ny=fny
	if (nx.LE.0 .OR. ny.LE.0) then
	  print*,'** InitIndRefGridFrmFile: invalid input'
	  stop
	end if
	if (allocated(IndRefGridArr(iElt)%IndRefGrid))
     &    deallocate(IndRefGridArr(iElt)%IndRefGrid)
	allocate(IndRefGridArr(iElt)%IndRefGrid(ny+1,nx+1))

	Write(*,*)'Reading in index gridfile ', filename
        Do j=1,ny
          Read(fu,*)IndRefGridArr(iElt)%IndRefGrid(j,1:nx)
        End Do
        close(fu)
	End Subroutine InitIndRefGridFrmFile


C***********************************************************************

	Subroutine CalcOrthoVec(vec,ovec)
	Implicit none
	Real*8 :: vec(3),ovec(3),vecl(3),dtmp
	Integer :: idx(3),i,j,itmp

	idx(1)=1; idx(2)=2; idx(3)=3
	vecl(1:3)=vec(1:3)
	Do i=1,2
	  Do j=i+1,3
            if (abs(vecl(i))<abs(vecl(j))) then
	      dtmp=vecl(i); vecl(i)=vecl(j); vecl(j)=dtmp
	      itmp=idx(i); idx(i)=idx(j); idx(j)=itmp
	    end if
	  End Do
	End Do

	! idx stores indexes of vec in decreasing absolute values
	ovec(idx(3))=0d0; ovec(idx(1))=vec(idx(2))
	ovec(idx(2))= -vec(idx(1))
	End Subroutine CalcOrthoVec

C***********************************************************************

	Function IsNSElt(iElt)
	Integer :: iElt
	Logical :: IsNSElt
	IsNSElt = (iElt==9).OR.(iElt==12)
	End Function IsNSElt

C***********************************************************************

	Subroutine ReadDblArray2D(fn,dim)
	use macos_mod, only : ObjInt
	Implicit none
	Character(len=MacosValLen) :: fn
	Integer :: dim,fu,i,m_err_rdarr2d

	If (.not.allocated(ObjInt)) Then
	  allocate(ObjInt(dim,dim),stat=m_err_rdarr2d)
	  if (m_err_rdarr2d /= 0)
     &      call macos_memory_failure
     &            ('ReadDblArray2D: allocate failed!')
	End If
	fu=12
	Open(unit=fu,file=fn,status='old')
	Do i=1,dim
	  Read(fu,*) ObjInt(1:dim,i)
	End Do
	Close(fu)
	End Subroutine ReadDblArray2D

C***********************************************************************

	Subroutine LOSCALC(zE,dx,ng,OPD,md,mp,CWF,los,losAng)
        use src_mod, only : RayMap
	use elt_mod, Only : LRayOK,LRayPass,BaseUnits

	Implicit none

	Integer :: ng,md,mp,off,i,j,iRay,jRay,ip1Ray,jp1Ray
	LREAL :: zE,dx,OPD(mp,mp),auv2ij,auv2ji,auv2,a2dwdx,a2dwdy,
     &           los(2),losAng(2)
	Complex*16 :: CWF(md,md)

	a2dwdx=0d0; auv2=0d0
	off=(md-ng)/2
	Do j=1,ng-1
	  Do i=1,ng-1
	    iRay=RayMap(i,j); ip1Ray=RayMap(i+1,j)
	    If (LRayOK(iRay) .and. LRayPass(iRay) .and.
     &          LRayOK(iRay) .and. LRayPass(ip1Ray)) Then
	      auv2ij=CWF(off+i,off+j)*DCONJG(CWF(off+i,off+j))
	      !a2dwdx=a2dwdx+auv2ij*(OPD(i+1,j)-OPD(i,j))/dx
	      a2dwdx=a2dwdx+(OPD(i+1,j)-OPD(i,j))/dx
	      !auv2=auv2+auv2ij
	      auv2=auv2+1d0
	    End If
	  End Do
	End Do
	losAng(1)=a2dwdx/auv2; los(1)=zE*losAng(1)

	a2dwdy=0d0; auv2=0d0
        Do i=1,ng-1
	  Do j=1,ng-1
	    jRay=RayMap(i,j); jp1Ray=RayMap(i,j+1)
	    If (LRayOK(jRay) .and. LRayPass(jRay) .and.
     &          LRayOK(jp1Ray) .and. LRayPass(jp1Ray)) Then
	      auv2ji=CWF(off+i,off+j)*DCONJG(CWF(off+i,off+j))
	      a2dwdy=a2dwdy+auv2ji*(OPD(i,j+1)-OPD(i,j))/dx
	      auv2=auv2+auv2ji
	    End If
          End Do
        End Do
	!print*,' ** LOSCALC: dx =',dx
	losAng(2)=a2dwdy/auv2; los(2)=zE*losAng(2)

        WRITE(*,100) losAng(2),',',losAng(1)
        WRITE(*,105) BaseUnits(1:2), ')) =',los(1),',',los(2)
 100    FORMAT(2P,'  System Line of Sight (X-tilt, Y-tilt (radian)) =',
     &         D16.7,A,D16.7)
 105    FORMAT(2P,'  System Line of Sight (PSF centroid at detector (',
     &         A,A,D16.7,A,D16.7)
	End Subroutine LOSCALC

C***********************************************************************

	Subroutine GetBeamSizeCmd(bmSize,npts_,
     &    ifPol_,ifGrid_,ifPropagate_,ifDraw_,
     &    ifLNegOK_,ifBuild_,ifDMBuild_,ifTrace_,ifRegrid_,
     &    ifIntSrfSet_,ifRayPrint_,iObsOpt_,iElt_,
     &    iCurWFElt_,iCurRayElt_,Dicr_,Djcr_,dxSource_,
     &    nBadRays_,StartIndRef_,StartExtinc_,OPDMat_,WFE_,
     &    DAvgL_,L1_,Cmatrix_,DrawRayVec_,DrawEltVec_)

        use elt_mod, only : EltID
        use tracesub_mod
        Implicit none

        Logical :: ifPol_,ifGrid_,ifPropagate_,ifDraw_,ifLNegOK_,
     &             ifBuild_,ifDMBuild_,ifTrace_,ifRayPrint_
        Logical :: ifRegrid_(mElt),ifIntSrfSet_(mElt),L1_(mRay),
     &             printDbg
        Integer :: npts_,iObsOpt_,iEndElt_,iCurWFElt_,iCurRayElt_,
     &             iElt_,nBadRays_
        Real*8 :: Dicr_,Djcr_,dxSource_,StartIndRef_,StartExtinc_,
     &            WFE_,DAvgL_,WFEPV_
        Real*8 :: OPDMat_(mpts,mpts)
        Integer :: DrawEltVec_(mDrawElt,mDrawRay)
        SREAL :: Cmatrix_(7,mCm,bRay),
     &           DrawRayVec_(2,mDrawElt,mDrawRay)

        Integer :: iRay,jRay,firstRayToLook
        Real*8 :: bmSize,g_ctr(3),dist(3),tmp_bs

        CALL OPD(npts_,ifPol_,ifGrid_,ifPropagate_,ifDraw_,
     &    ifLNegOK_,ifBuild_,ifDMBuild_,ifTrace_,ifRegrid_,
     &    ifIntSrfSet_,ifRayPrint_,iObsOpt_,iElt_,
     &    iCurWFElt_,iCurRayElt_,Dicr_,Djcr_,dxSource_,
     &    nBadRays_,StartIndRef_,StartExtinc_,OPDMat_,WFE_,WFEPV_,
     &    DAvgL_,L1_,Cmatrix_,DrawRayVec_,DrawEltVec_)

        ! Reset, so next OPD() will trace ray from source
        ifTrace_=.false.; iCurRayElt_=0

        ! Compute beam size at iElt_

        ! Note that command-line 'OPD' command refuses to compute OPD
        ! at a segment or a nonseq element, this routine however will call
        ! OPD (i.e. trace rays to) on such elements.
        ! Note when iElt_ is a segment or a non-sequential element, all rays
        ! will be traced to the group of non-seq elements, so the beam size
        ! estimated in such a case will be the *global* beam size for all
        ! non-seq elements in the group. -jzlou 08/22/2015

        firstRayToLook = 1  ! default to chief ray

        If (EltID(iElt_)==9 .OR. EltID(iElt_)==11 .OR.
     &      EltID(iElt_)==12) Then
          ! Segment or Nonseq Refractor or Reflector
          firstRayToLook = 2 ! don't include chief ray
        End If

        ! First estimate beam center position
        g_ctr(1:3)=0d0

        jRay=0; bmSize=0d0
        Do iRay=firstRayToLook,nRay
          If (LRayOK(iRay) .and. LRayPass(iRay)) Then
            g_ctr(1:3)=g_ctr(1:3)+RayPos(1:3,iRay)
            jRay=jRay+1
          End If
        End Do
        If (jRay>0) Then
          g_ctr=g_ctr/jRay
        Else
          Write(*,*)' ** All rays are lost at element',iElt_
          return
        End If
        !g_ctr(1:3)=RayPos(1:3,1)  ! use chief ray
        Do iRay=firstRayToLook,nRay
          If (LRayOK(iRay) .and. LRayPass(iRay)) Then
            dist(1:3)=RayPos(1:3,iRay)-g_ctr(1:3)
            tmp_bs=dist(1)**2+dist(2)**2+dist(3)**2
            if (tmp_bs>bmSize) bmSize=tmp_bs
          End If
        End Do
        bmSize=sqrt(bmSize)

        printDbg = .false. ! flag to control local debug print
        If (printDbg) Then
          print*,'** GetBeamSizeCmd: g_ctr =',g_ctr(1:3)
          print*,'** GetBeamSizeCmd: firstRayToLook =',firstRayToLook
          print*,'** GetBeamSizeCmd: nRay, iElt_ =',nRay,iElt_
        End If
        End Subroutine GetBeamSizeCmd

C***********************************************************************
C
C  Perform surface metrology measurement
C
	Subroutine SrfMetCalc(nElt,iEltToMetSrf,nMetPos,
     &                        tMetSrf,metPos,metMeas)
	use param_mod, only : mElt
	use elt_mod, only : mMetSrf,mMetPos,mMetBeam,ntMetPos,
     &                      metMeasBuf,metBeamFlg,nMetMeas
	use traceutil_mod, only: ShowMetData_FLG
	use math_mod, only : DDOTC

	Implicit none

	Integer :: nElt,iEltToMetSrf(mElt),nMetPos(mElt),tMetSrf(mElt),
     &             iElt,i,j,iMetBeam,jMetBeam
	Real*8 :: metPos(3,mMetPos,mMetSrf),
     &            metMeas(mMetBeam,mMetSrf),dd(3)

	jMetBeam=0;  metMeasBuf=0d0
	Do iElt=1,nElt
          If (iEltToMetSrf(iElt)>0 .AND. tMetSrf(iElt)>0) Then
            ! This surface has metrology point(s) and target metrology
	    ! element defined
            iMetBeam=0
	    Do j=1,nMetPos(iElt) ! loop thru metPos on this src srf
	      Do i=1,ntMetPos(iElt) ! loop thru target srf metPos
	        If (metBeamFlg(i,j,iEltToMetSrf(iElt))>0) Then
  	          dd(1:3)=metPos(1:3,j,iEltToMetSrf(iElt))-
     &                    metPos(1:3,i,iEltToMetSrf(tMetSrf(iElt)))
	          iMetBeam=iMetBeam+1
	          metMeas(iMetBeam,iEltToMetSrf(iElt))=sqrt(DDOTC(dd,dd))
#if 1
	          If (ShowMetData_FLG) Then
	            Write(*,*)' Met beam',iMetBeam,'length from Elt',
     &                iElt,'Pos',j,' to Elt',tMetSrf(iElt),'Pos',i,
     &                ' =',metMeas(iMetBeam,iEltToMetSrf(iElt))
	          End If
		  ! to output metrology measurement from GMI interface
                  jMetBeam=jMetBeam+1
                  metMeasBuf(jMetBeam)
     &              =metMeas(iMetBeam,iEltToMetSrf(iElt))
#endif
	        End If
	      End Do
	    End Do
	    nMetMeas=jMetBeam
            If (ShowMetData_FLG)
     &        Write(*,*)' --------------------------------------'
	  End If
	End Do
	If (.not.ShowMetData_FLG)
     &    Write(*,*)'** All metrology beam lengths updated'
	End Subroutine SrfMetCalc

#include "collim.inc"

C***********************************************************************

	Function GetNumToken(Val)
        Implicit none

        Character(len=MacosValLen) :: Val
        Integer :: GetNumToken, nbc,c
        Logical :: inTok

        Do nbc=1,MacosValLen
          if (Val(nbc:nbc) /= ' ') exit
        End Do
        If (nbc .GE. MacosValLen) Then
          GetNumToken=0
          return
        End If

        inTok=.false.; GetNumToken=0
        Do c=nbc,MacosValLen
          if (.not. inTok .and. Val(c:c) /= ' ') then
            GetNumToken = GetNumToken + 1; inTok=.true.
          else if (inTok .and. Val(c:c)==' ') then
            inTok=.false.
          end if
        End Do
        !print*,' ** GetNumToken =',GetNumToken
        End Function GetNumToken

C***********************************************************************

	Subroutine SetIndRef4Refractors(wl_id)
	use elt_mod, only : nElt,EltID,IndRef,ArrIndRef
	Implicit none
	Integer :: wl_id,ie

	Do ie=1,nElt
 	  If (EltID(ie).GT.3 .AND. EltID(ie).NE.7 .AND.
     &        EltID(ie).NE.8) Then
	    ! not Reflector,FocalPlane,Reference,Obscuring,Return surfaces
	    If (ArrIndRef(wl_id,ie) > 0) Then
	      IndRef(ie)=ArrIndRef(wl_id,ie)
	    Else
	      If (wl_id>1 .AND. IndRef(ie)/=1d0) Then
	        WRITE(*,*)' ** Invalid IndRef value in ArrIndRef for elt ',
     &          ie
	      End If
	    End If
	  End If
	End Do
	End Subroutine SetIndRef4Refractors

C***********************************************************************

	Function getIndex_grin5(crnt_pos,pv,zdir,modelCoef)
	use math_mod, only : DDOTC
	Implicit none
	Real*8 :: getIndex_grin5,crnt_pos(3),pv(3),zdir(3),modelCoef(10)
!	Real*8, parameter :: n0=1.62893, n2=-7.673417129133d-2,
!    &          n4=3.153771895506d-3
        Real*8 :: n0,n2,n4,rvec(3),r2,r4
	rvec=(crnt_pos-pv)-DDOTC(crnt_pos-pv,zdir)*zdir
	n0=modelCoef(2); n2=modelCoef(3); n4=modelCoef(4)
        r2=rvec(1)**2+rvec(2)**2+rvec(3)**2; r4=r2*r2
	getIndex_grin5=n0+n2*r2+n4*r4
	End Function getIndex_grin5


	Function getIndexGrad_grin5(crnt_pos,pv,zdir,modelCoef)
	use math_mod, only : DDOTC
	Implicit none
        Real*8 :: getIndexGrad_grin5(3),crnt_pos(3),pv(3),zdir(3),
     &            modelCoef(10)
	Real*8 :: n0,n2,n4,n,rvec(3),r2,r4,getIndex_grin5
!	Real*8, parameter :: n0=1.62893, n2=-7.673417129133d-2,
!    &          n4=3.153771895506d-3
        n0=modelCoef(2); n2=modelCoef(3); n4=modelCoef(4)
        n=getIndex_grin5(crnt_pos,pv,zdir,modelCoef)
        rvec=(crnt_pos-pv)-DDOTC(crnt_pos-pv,zdir)*zdir
	r2=rvec(1)**2+rvec(2)**2+rvec(3)**2; r4=r2*r2
        getIndexGrad_grin5(1)=2d0*rvec(1)*(n2+2*n4*r2)
	getIndexGrad_grin5(2)=2d0*rvec(2)*(n2+2*n4*r2)
	getIndexGrad_grin5(3)=2d0*rvec(3)*(n2+2*n4*r2)
        getIndexGrad_grin5(1:3)=n*getIndexGrad_grin5(1:3)
	End Function getIndexGrad_grin5

C***********************************************************************

	Function Calc_RayDielectricTC(wl)
	use elt_mod, only : czero_
	Implicit none
	Real*8 :: wl
	Complex*16 :: Calc_RayDielectricTC
	Calc_RayDielectricTC = czero_
	End Function Calc_RayDielectricTC

C***********************************************************************

	Function Calc_RayMetalTC(wl)
	use elt_mod, only : czero_
	Implicit none
	Real*8 :: wl
	Complex*16 :: Calc_RayMetalTC
	Calc_RayMetalTC = czero_
        End Function Calc_RayMetalTC

C***********************************************************************

        ! 'thsElt' is a duplicate element of 'dupElt' element
        Subroutine MakeDupElt(thsElt,dupElt,jGridSrf)
        use elt_mod
        use cfiles_mod, only : GridFile
        Implicit none
        Integer :: thsElt,dupElt,jGridSrf,iObs,nVertex,nside,iside,j

#if 1
        ! Make a copy of element 'dupElt' at element 'thsElt'
        EltName(thsElt)=EltName(dupElt)
        EltMRP(thsElt)=EltMRP(dupElt)
        EltID(thsElt)=EltID(dupElt)
        SrfType(thsElt)=SrfType(dupElt)
        LnkElt(thsElt)=dupElt
        NSCnt(thsElt)=NSCnt(dupElt)
        KrElt(thsElt)=KrElt(dupElt)
        KcElt(thsElt)=KcElt(dupElt)
        fElt(thsElt)=fElt(dupElt); eElt(thsElt)=eElt(dupElt)
          if (.false.) then
           print*,'dupElt, fElt(dupElt) =',dupElt,fElt(dupElt)
           print*,'dupElt, eElt(dupElt) =',dupElt,eElt(dupElt)
           print*,'thsElt, fElt(thsElt) =',thsElt,fElt(thsElt)
           !stop
          end if
        psiElt(1:3,thsElt)=psiElt(1:3,dupElt)
        VptElt(1:3,thsElt)=VptElt(1:3,dupElt)
        RptElt(1:3,thsElt)=RptElt(1:3,dupElt)
        xObs(1:3,thsElt)=xObs(1:3,dupElt)
        nObs(thsElt)=nObs(dupElt)
        If (nObs(thsElt)>0) then
          Do iObs=1,nObs(thsElt)
            ObsType(iObs,thsElt)=ObsType(iObs,dupElt)
            ObsVec(1:3,iObs,thsElt)=ObsVec(1:3,iObs,dupElt)
            If (ABS(ObsType(iObs,thsElt)).EQ.8) Then
              nVertex=ObsVec(3,iObs,thsElt)
              Do j=1,nVertex
                PolyObsVtx(1:2,j,iObs,thsElt)
     &              =PolyObsVtx(1:2,j,iObs,dupElt)
              End Do
              nside= ObsVec(3,iObs,dupElt)
              Do iside=1,nside
                PolyObsEqn(1:3,iside,iObs,thsElt)=
     &              PolyObsEqn(1:3,iside,iObs,dupElt)
              End Do
            End If
          End Do
        End If
        LApTypeDefined(thsElt)=LApTypeDefined(dupElt)
        LApVecDefined(thsElt)=LApVecDefined(dupElt)
        ApType(thsElt)=ApType(dupElt)
        ApVec(:,thsElt)=ApVec(:,dupElt)
        If (ApType(thsElt).EQ.6) Then
          Call SetHexBound(ApVec(1,thsElt),thsElt)
        End If
        IF (ApType(thsElt).EQ.7 .OR. ApType(thsElt).EQ.8) THEN
          nVertex=ApVec(3,thsElt)
          Do j=1,nVertex
            PolyApVtx(1:2,j,thsElt)=PolyApVtx(1:2,j,dupElt)
          End Do
          nside= ApVec(3,thsElt)
          Do iside=1,nside
            PolyApEqn(1:3,iside,thsElt)=
     &          PolyApEqn(1:3,iside,dupElt)
          End Do
        End If
        pMon(1:3,thsElt)= pMon(1:3,dupElt)
        xMon(1:3,thsElt)= xMon(1:3,dupElt)
        yMon(1:3,thsElt)= yMon(1:3,dupElt)
        zMon(1:3,thsElt)= zMon(1:3,dupElt)
        lMon(thsElt)= lMon(dupElt)
        nGridMat(thsElt)=nGridMat(dupElt)
        GridSrfdx(thsElt)=GridSrfdx(dupElt)
        GridFile(thsElt)=GridFile(dupElt)
        If (iEltToGridSrf(dupElt).NE.0) Then
          jGridSrf=jGridSrf+1
          iEltToGridSrf(thsElt)=jGridSrf
        End If
        ZernTypeL(thsElt)=ZernTypeL(dupElt)
        ZernCoef(:,thsElt)=ZernCoef(:,dupElt)
        IndRef(thsElt)=IndRef(dupElt)
        Extinc(thsElt)=Extinc(dupElt)
        LUpdateTElt_FLG(thsElt)=LUpdateTElt_FLG(dupElt)
        ifOPDModGrid(thsElt)=ifOPDModGrid(dupElt)
        zElt(thsElt)=zElt(dupElt)
        PropType(thsElt)=PropType(dupElt)
        nECoord(thsElt)=nECoord(dupElt)
        TElt(1:6,1:nECoord(thsElt),thsElt)=
     &    TElt(1:6,1:nECoord(dupElt),dupElt)
#endif
        End Subroutine MakeDupElt

#if 0
! Dump fits write routine
	Subroutine fitswrite(FitsFn,OPD,nsz)
	Character(len=80) :: FitsFn
	Real*8, dimension(:,:) :: OPD
	Integer :: nsz
	! ...
	End Subroutine  fitswrite
#endif

C***********************************************************************

#if 0
 	Subroutine dOPL(ifLoad,Rx1,Rx2,iElt,mode)
 	Logical :: ifLoad
 	Character(len=*) :: Rx1,Rx2,mode
 	Integer :: iElt
	! ...
 	End Subroutine dOPL
#endif

C***********************************************************************
C****** End file utilsub.F *********************************************
C***********************************************************************

