C***********************************************************************
C	Begin file elemsub.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C 11/16/03
C This file has been updated by J. Lou for adding SrfType 9 to Refractor
C routine, and other f90 stuff. It should be copied to MACOS on Unix.

C 09/21/05
C Converted the file into a f90 module, for the convenience of expanding
C AsphCoef array

C 06/16/06
C Trying to add support for multi-layer thin-film coating on optical
C surfaces during ray-tracing of polarized light.

C***********************************************************************

#include "realtype.h"

        MODULE elemsub_mod
        use math_mod
        use surfsub

        IMPLICIT NONE

        PUBLIC
        INTEGER, PRIVATE, PARAMETER :: max_coat_layers=10,mSB=10000
        REAL*8, PUBLIC              :: NhatP(3),Nsrfm(3)

        CONTAINS

        SUBROUTINE FindSrf(ifLNsrf,ifObs,
     &  LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,UDSrfCoef,UDSrfType,
     &  UDSrfParam,nGridMat,mGridMat,GridMat,GridSrfdx,
     &  psi,pv,prot,pin,ihat,
     &  rad,pmon,xhat,yhat,zhat,
     &  ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &  pout,L,SrfType,
     &  XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,nAsphCoef)

        use traceutil_mod
        use didesub, only: IntSrf
        use math_mod
        use elt_mod, only: mMonCoef,mAnaCoef


        ! Compute reflection point, normal, reflection and transition
        ! matrices for a reflecting surface

        IMPLICIT NONE

        INTEGER nGridMat,mGridMat,nAsphCoef
        REAL*8 GridMat(mGridMat,mGridMat),GridSrfdx

        INTEGER SrfType,NDP,IWORK(*),IERROR,UDSrfType,mDP,
     &          ApType,nObs,ObsType(nObs)

        LOGICAL LROK,ifBuild,ifLNsrf,ifObs,LRT,ifPol !,ChkRayTrans

        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &   MonCoef(mMonCoef),pr(3),L,Nhat(3),pout(3),rho(3),
     &   AnaCoef(mAnaCoef),UDSrfCoef(*),DWORK(*),xGrid(3),yGrid(3),
     &   zGrid(3),ApVec(6),ObsVec(6,nObs),UDSrfParam(6),
     &   AsphCoef(nAsphCoef)

        REAL*8 XYZSRF(3,NDP),DZSRF(2,NDP),Nvec(3),Nmag,dNdp(3,3),wVec(3)
        REAL*8 Kc,Kr,rad,pmon(3),xhat(3),yhat(3),zhat(3)

        ! Compute incidence point of ray

        SELECT CASE (SrfType)
          CASE (1)
          CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &                   pr,L,Nvec,Nmag,Nhat,dNdp)

          CASE (2)
            CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &                  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)

          CASE (3)
            CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &                   prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,
     &                   nAsphCoef)

          CASE (4,8)
            CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &                  prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &                  Nvec,Nmag,Nhat,dNdp)

          CASE (5)
            CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &                  XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,wVec,pr,
     &                  L,pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)

          CASE (6)
            CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &                  prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)

          CASE (7)
            CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &                 UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,
     &                 pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)

          CASE (9)
            CALL GridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,
     &                   ihat,pr,L,nGridMat,mGridMat,GridMat,GridSrfdx,
     &                   pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)

          CASE (10)
            CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &                    prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)

          CASE (11)
            CALL AsphGridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AnaCoef,psi,pv,
     &                       prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,
     &                       nAsphCoef,nGridMat,mGridMat,GridMat,
     &                       GridSrfdx,pmon,xhat,yhat,zhat)

          CASE DEFAULT
            LROK = .FALSE.
            WRITE(*,*)'FindSrf(): Unknown surface type = ',SrfType,
     &                ', Aborting trace.'
            RETURN

        END SELECT

        ! If not obscured, check if ray has negative length

        IF (LRT.AND.(L.LT.1d-22)) THEN

          ! Here if ray misses surface
          LROK=.FALSE.

        ELSE

          ! Check if ray is obscured
          IF (ifObs) THEN
            rho = pin - pv + L*ihat

            LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,psi,rho,ApType,
     &                      ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,
     &                      zGrid,taperFactor(crnt_traced_ray))
          ELSE
            LRT=.TRUE.
          END IF
        END IF

        END SUBROUTINE FindSrf

C***********************************************************************

	SUBROUTINE Reflector(
     &    ifBuild,ifDMBuild,ifPol,ifLNsrf,
     &    ifObs,LROK,LRT,Kc,
     &    Kr,AsphCoef,MonCoef,
     &    AnaCoef,UDSrfCoef,UDSrfType,
     &	  UDSrfParam,nGridMat,mGridMat,
     &    GridMat,GridSrfdx,
     &	  psi,pv,prot,
     &    pin,ihat,Evec, ! RayPos(1,iRay),RayDir(1,iRay),RayE(1,iRay)
     &    rad,pmon,
     &    xhat,yhat,zhat, ! xMon,yMon,zMon
     &	  na,nb,
     &    kxa,kxb,lambda,
     &	  ApType,ApVec,nObs,ObsType,
     &    ObsVec,xGrid,yGrid,zGrid, ! ObsVec,xObs,yObs,zObs
     &	  pout,rout,L,Eout, ! RayPos(1,iRay),RayDir(1,iRay),RayL(iRay),Evec
	                    ! Evec, Eout are both RayE in the calling routine. -jzlou
     &    dxidxim1,dxidui,SrfType,
     &	  XYZSRF,DZSRF,IERROR,
     &    IWORK,DWORK,NDP,mDP,nCoat,CoatW,nb_arr,kxb_arr)

!  pin = RayPos(1,iRay)
!  ihat = RayDir(1,iRay)
!  prot = RayE(1,iRay),

	USE Constants, ONLY: TWOPI
	use macos_debug
	use traceutil_mod
	use didesub, only: IntSrf
        use elt_mod, only: mMonCoef,mAnaCoef

C  ierror,iwork,ndp,wvec,anacoef,udsrftype,udsrfcoef + intsrf stuff

C  Compute reflection point, normal, reflection and transition
C  matrices for a reflecting surface

	IMPLICIT NONE

	INTEGER nGridMat,mGridMat
	REAL*8 GridMat(mGridMat,mGridMat),GridSrfdx

	INTEGER i,j,ip3,jp3,SrfType,NDP,IWORK(*),IERROR,
     &	  iDebug,UDSrfType,mDP,
     &	  ApType,nObs,ObsType(nObs),nCoat,ic

	LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &	        LRT

	REAL*8, dimension(:) :: CoatW
        REAL*8, dimension(0:) :: nb_arr,kxb_arr

	SREAL :: dxidxim1(7,7),dxidui(7,6)

        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &    MonCoef(mMonCoef),pr(3),rhat(3),L,Nhat(3),lambda,pout(3),rout(3),
     &	  na,nb,kxa,kxb,AnaCoef(mAnaCoef),UDSrfCoef(*),DWORK(*),
     &	  xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs),
     &    UDSrfParam(6),dd
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 XYZSRF(3,NDP),DZSRF(2,NDP),
     &	  wVec(3),rho(3),drdthp(3,3)

	REAL*8 Nvec(3),Nmag,dNdp(3,3),R(3,3),dpdd(3,3),dNdd(3,3),
     &    dNdi(3,3),Ndoti,Niinv,Nouti(3,3),drdN(3,3),drdd(3,3),
     &    drdth(3,3),
     &    drdi(3,3),drda(3,3),dgda(3,3),dgdi(3,3),dgdd(3,3),dldd(3),
     &    dprot(3),dddth(3,3),dgdth(3,3),dldth(3),Kc,Kr,
     &    NoidNdi(3,3),rad,pmon(3),xhat(3),yhat(3),zhat(3),
     &	  svec(3),shat(3),prhat(3),pihat(3),sfa

	REAL*8 D1(3,3),D2(3,3),D3(3,3),S1

	COMPLEX*16 Evec(3),Esrf(3),Eout(3),RP,RS,Epi,Esi,Epr,Esr,
     &	           C1,C2,NaCmplx,NbCmplx,ccfa,ccfb,cmu,cNdi,cone,
     &             ccfb_arr(0:max_coat_layers),RP1,RS1,TP,TS,TP1,TS1,
     &             DD1,DD2,RPnew,RSnew

	DATA idebug/0/

C  Compute incidence point of ray

	IF (SrfType.EQ.1) THEN
	  CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &	    pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.2) THEN
	  CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &	    prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.3) THEN
	  CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &	    prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
	ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
	  CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &	    prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &	    Nvec,Nmag,Nhat,dNdp)
cides
	ELSE IF (SrfType.EQ.5) THEN
	  CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &	    XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	    wVec,pr,L,pmon,xhat,yhat,zhat,
     &	    Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.6) THEN
	  CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	    prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.7) THEN
	  CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &	    UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,
     &	    pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.9) THEN
	  CALL GridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,
     &	    psi,pv,prot,pin,ihat,pr,L,
     &	    nGridMat,mGridMat,GridMat,GridSrfdx,
     &	    pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.10) THEN
	  CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	         prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	  !print*,'In Reflector(), Aft ToricSrf, LROK=',LROK
        ELSE IF (SrfType.EQ.11) THEN
          !print*,'***** In Reflector, B4 CALL AsphGridSrf'
	  !STOP
          CALL AsphGridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,
     &      psi,pv,prot,pin,
     &      ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef,
     &      nGridMat,mGridMat,GridMat,GridSrfdx,pmon,
     &      xhat,yhat,zhat)
	ELSE
	  WRITE(*,*)'Reflector: Unknown surface type. Aborting trace.'
	END IF


C  Return if ray missed surface

	IF (.NOT.LROK) RETURN
cdebug

	  Nhat(1)=Nvec(1)/Nmag
	  Nhat(2)=Nvec(2)/Nmag
	  Nhat(3)=Nvec(3)/Nmag
          Nsrfm=Nhat
cdebug

C  Compute reflected ray

	S1=2d0*(Nhat(1)*ihat(1)+Nhat(2)*ihat(2)+Nhat(3)*ihat(3))
	DO 1 i=1,3
	  rhat(i)=ihat(i)-S1*Nhat(i)
   1	CONTINUE

	if (.false.) then !if (iRay_debug == 2526 .and. iElt_debug == 90) then
	   print*,'**** rhat = ', rhat, ' Nhat(1),Nmag =',
     &        Nhat(1), Nmag, ' Nvec = ', Nvec,
     &        'SrfType = ', SrfType
	end if


C  Polarization calcs

        ! Propagation through isotropic medium to a surface,
	! applicable to both non-coated and coated surfaces
	If (ifPol) Then
          cone=1d0
          S1=-TWOPI*L/lambda
          C2=DCMPLX(0d0,S1)*NaCmplx
          C1=CDEXP(C2)
          Esrf(1:3)=C1*Evec(1:3) ! Esrf = electric field at the surface,
	                         ! C1 is a phase factor resulted from OPL,
	                         ! and Evec is "input" E field.
                                 ! -jzlou

	  !print*,'** E at surface: Esrf=',Esrf(1:3)

	  ! svec = incident S pol vector, normal to plane of incident ray
          ! and surf normal. -jzlou
          CALL DXPROD(svec,ihat,Nhat)
          sfa=DSQRT(svec(1)*svec(1)+svec(2)*svec(2)+svec(3)*svec(3))
          IF (DABS(sfa).LT.1d-10) THEN
            IF (ihat(2).NE.0d0) THEN
              svec(1)=-ihat(2)
              svec(2)=ihat(3)
              svec(3)=ihat(1)
            ELSE IF (ihat(3).NE.0d0) THEN
              svec(1)=ihat(2)
              svec(2)=-ihat(3)
              svec(3)=ihat(1)
            ELSE IF (ihat(1).NE.0d0) THEN
              svec(1)=ihat(2)
              svec(2)=ihat(3)
              svec(3)=-ihat(1)
            END IF
            CALL DXPROD(shat,ihat,svec)
            CALL DUNITIZE(shat)
          ELSE
            shat(1)=svec(1)/sfa  ! normalized incident S pol vector, -jzlou
            shat(2)=svec(2)/sfa
            shat(3)=svec(3)/sfa
          END IF
          CALL DXPROD(pihat,shat,ihat) ! incident P pol vector, -jzlou
          CALL DXPROD(prhat,shat,rhat) ! reflected P pol vector, -jzlou

	  Epi=DCMPLX(pihat(1))*Esrf(1)+DCMPLX(pihat(2))*Esrf(2)+
     &        DCMPLX(pihat(3))*Esrf(3) ! Incident E field at surface in P direction,-jzlou
          Esi=DCMPLX(shat(1))*Esrf(1)+DCMPLX(shat(2))*Esrf(2)+
     &        DCMPLX(shat(3))*Esrf(3) ! Incident E field at surface in S direction,-jzlou
	End if

	! Applicable to non-coated surface only
	IF (ifPol .and. nCoat==0) THEN
	  NaCmplx=DCMPLX(na,(-kxa)) ! complex index of refraction = na-j*extinc, -jzlou
	  NbCmplx=DCMPLX(nb,(-kxb)) ! extinc = kxa
	  cmu=NaCmplx/NbCmplx

	  ! Compute refracted ray cosine
	  cNdi=ihat(1)*Nhat(1)+ihat(2)*Nhat(2)+ihat(3)*Nhat(3)
	  C1=cone-cmu*cmu*(cone-cNdi*cNdi)
	  ccfb=CDSQRT(C1)  ! = cos(theta2), -jzlou

	  ! cosine of incidence ray with surface normal -jzlou
	  ccfa=ABS(cNdi) ! = cos(theta1), -jzlou

	  ! RP=Reflectance in P pol, RS=Reflectance in S pol, -jzlou
	  !RP=(NbCmplx*ccfa-NaCmplx*ccfb)/(NbCmplx*ccfa+NaCmplx*ccfb) ! dcr's original
	  RP=(NaCmplx*ccfb-NbCmplx*ccfa)/(NbCmplx*ccfa+NaCmplx*ccfb) ! fixed by jzlou
	  RS=(NaCmplx*ccfa-NbCmplx*ccfb)/(NaCmplx*ccfa+NbCmplx*ccfb)
	END IF ! ifPol & nCoat=0


	! Applicable to surface with thin-film coating layers
	IF (ifPol .and. nCoat/=0) THEN

	  ! First calc refraction angles for each layer, from outmost inward
          NaCmplx=DCMPLX(nb_arr(0),(-kxb_arr(0)))
          NbCmplx=DCMPLX(nb_arr(1),(-kxb_arr(1)))

	  !print*,'Reflector/coatings: nb_arr(3)=',nb_arr(3)
          !print*,'Reflector/coatings: kxb_arr(3)=',kxb_arr(3)

	  ccfb_arr(0)=DDOTC(ihat,Nhat) ! cosine of incident ray angle
          cmu=NaCmplx/NbCmplx
          cNdi=ihat(1)*Nhat(1)+ihat(2)*Nhat(2)+ihat(3)*Nhat(3)
          C1=cone-cmu*cmu*(cone-cNdi*cNdi)
          ccfb_arr(1)=CDSQRT(C1) ! cosine of refraction angle for outmost interface

	  !print*,'cosine: ccfb_arr(0) = ', ccfb_arr(0)
	  !print*,'nCoat = ', nCoat

          Do ic=1,nCoat
            NaCmplx=DCMPLX(nb_arr(ic),(-kxb_arr(ic)))
            NbCmplx=DCMPLX(nb_arr(ic+1),(-kxb_arr(ic+1)))
            cmu=NaCmplx/NbCmplx
            C1=cone-cmu*cmu*(cone-ccfb_arr(ic)*ccfb_arr(ic))
            ccfb_arr(ic+1)=CDSQRT(C1)
          End Do

	  !print*,'Reflector/coatings: nb_arr = ', nb_arr(1:nCoat+1)
	  !print*,'Reflector/coatings: kxb_arr = ', kxb_arr(1:nCoat+1)

	  ! Now calc P and S reflectances recursively for each interface,
          ! innermost interface first, then for each outer coating layer

          ! For innermost (substrate) interface, just use the Fresnel
          ! reflectance eqn
          NaCmplx=DCMPLX(nb_arr(nCoat),(-kxb_arr(nCoat)))
          NbCmplx=DCMPLX(nb_arr(nCoat+1),(-kxb_arr(nCoat+1)))
          RP=(NaCmplx*ccfb_arr(nCoat+1)-NbCmplx*ccfb_arr(nCoat))
     &       /(NbCmplx*ccfb_arr(nCoat)+NaCmplx*ccfb_arr(nCoat+1))
          RS=(NaCmplx*ccfb_arr(nCoat)-NbCmplx*ccfb_arr(nCoat+1))
     &       /(NaCmplx*ccfb_arr(nCoat)+NbCmplx*ccfb_arr(nCoat+1))

	  !print*,'Innermost RP = ', RP
	  !print*,'Innermost RS = ', RS

	  if (.false.) then
	   ! Add transmittance
           TS=2d0*NaCmplx*ccfb_arr(nCoat)
     &        /(NaCmplx*ccfb_arr(nCoat)+NbCmplx*ccfb_arr(nCoat+1))
           TP=2d0*NaCmplx*ccfb_arr(nCoat)
     &        /(NaCmplx*ccfb_arr(nCoat+1)+NbCmplx*ccfb_arr(nCoat))
	  end if

	  ! For each coating layer going outward, use thin film
          ! reflectance eqn,
          Do ic=nCoat,1,-1
            NaCmplx=DCMPLX(nb_arr(ic-1),(-kxb_arr(ic-1)))
            NbCmplx=DCMPLX(nb_arr(ic),(-kxb_arr(ic)))

            ! Calc phase delay due to pass thru coating
            !dd=CoatW(ic)/ccfb_arr(ic) ! ray pathlength in this layer
            !C2=DCMPLX(0d0,1d0)*TwoPi*2d0*dd*NbCmplx/lambda
            dd=CoatW(ic)*ccfb_arr(ic)
            C2=DCMPLX(0d0,1d0)*2d0*(TWOPI/lambda)*dd*NbCmplx
            C1=CDEXP(-C2)

            ! Reflectance
            RP1=(NaCmplx*ccfb_arr(ic)-NbCmplx*ccfb_arr(ic-1))
     &          /(NbCmplx*ccfb_arr(ic-1)+NaCmplx*ccfb_arr(ic))
            !RP1=-RP1 ! testing

C       print*,'layer: NaCmplx,NbCmplx,ccfb_arr(ic),ccfb_arr(ic-1)=',
C    &        NaCmplx,NbCmplx,ccfb_arr(ic),ccfb_arr(ic-1)
C       print*,'layer RP1 = ', RP1

            RPnew=(RP1+RP*C1)/(1d0+RP1*RP*C1)

            RS1=(NaCmplx*ccfb_arr(ic-1)-NbCmplx*ccfb_arr(ic))
     &          /(NaCmplx*ccfb_arr(ic-1)+NbCmplx*ccfb_arr(ic))
            RSnew=(RS1+RS*C1)/(1d0+RS1*RS*C1)

	    if (.false.) then
             ! Transmittance
             DD2=DCMPLX(0d0,1d0)*(TWOPI/lambda)*dd*NbCmplx
             DD1=CDEXP(-DD2)
             TP1= 2d0*NaCmplx*ccfb_arr(ic-1)
     &            /(NbCmplx*ccfb_arr(ic-1)+NaCmplx*ccfb_arr(ic))
             TP= (TP1*TP*DD1)/(1d0+RP1*RP*C1)
             TS1= 2d0*NaCmplx*ccfb_arr(ic-1)
     &            /(NaCmplx*ccfb_arr(ic-1)+NbCmplx*ccfb_arr(ic))
             TS= (TS1*TS*DD1)/(1d0+RS1*RS*C1)
	    end if

            RP=RPnew; RS=RSnew
          End Do
	END IF ! ifPol & nCoat/=0

	! Check result
	!if (.false. .and. ifPol) then
	!  print*,'*** In Reflector: RP = ', RP
	!  print*,'*** In Reflector: RS = ', RS
	!  !print*,'*** In Reflector: Epi = ', Epi
	!  !print*,'*** In Reflector: Esi = ', Esi
	!end if

	IF (ifPol) THEN
	  Epr=RP*Epi
	  Esr=RS*Esi
	  ! Note that Eout, the E field, is a 3D vector field.
	  Eout(1)=DCMPLX(prhat(1))*Epr+DCMPLX(shat(1))*Esr
	  Eout(2)=DCMPLX(prhat(2))*Epr+DCMPLX(shat(2))*Esr
	  Eout(3)=DCMPLX(prhat(3))*Epr+DCMPLX(shat(3))*Esr
	END IF


C  Compute partials

	IF (ifBuild) THEN

	  CALL DREFLECT(R,Nhat)

C    Partials of the normal

	  CALL DOUTER(Nouti,ihat,Nhat)
	  Ndoti=DDOTC(Nhat,ihat)
	  Niinv=1d0/Ndoti
	  CALL DSMPROD(NoidNdi,Nouti,Niinv,9)
	  CALL DEQUATE(dpdd,NoidNdi,9)
	  dpdd(1,1)=dpdd(1,1)-1d0
	  dpdd(2,2)=dpdd(2,2)-1d0
	  dpdd(3,3)=dpdd(3,3)-1d0
	  CALL DMPROD(dNdd,dNdp,dpdd,3,3,3)
	  S1=-L
	  CALL DSMPROD(dNdi,dNdd,S1,9)

C    Rotation-to-translation partials

	  CALL DSUB(dprot,prot,pr,3)
	  CALL DXMAT(dddth,dprot)

C    Ray direction partials

	  CALL DOUTER(D3,Nhat,ihat)
	  D3(1,1)=D3(1,1)+Ndoti
	  D3(2,2)=D3(2,2)+Ndoti
	  D3(3,3)=D3(3,3)+Ndoti
	  S1=-2d0
	  CALL DSMPROD(drdN,D3,S1,9)
	  CALL DMPROD(drdd,drdN,dNdd,3,3,3)
	  CALL DSMPROD(D2,rhat,S1,3)
	  CALL DXMAT(D1,D2)
	  CALL DPROJECT(D2,Nhat)
	  CALL DMPROD(drdthp,D1,D2,3,3,3)
	  CALL DMPROD(D2,drdd,dddth,3,3,3)
	  CALL DADD(drdth,D2,drdthp,9)
	  CALL DMPROD(D1,drdN,dNdi,3,3,3)
	  CALL DADD(drdi,R,D1,9)
	  CALL DNEGATE(drda,drdd,9)

C    Beamwalk partials

	  CALL DEQUATE(dgda,R,9)
	  CALL DSMPROD(dgdi,R,L,9)
	  CALL DPROJECT(D1,rhat)
	  CALL DMPROD(dgdd,D1,NoidNdi,3,3,3)
	  CALL DMPROD(dgdth,dgdd,dddth,3,3,3)

C    Pathlength partials

	  S1=((1d0-DDOTC(rhat,ihat))/Ndoti)
	  CALL DSMPROD(dldd,Nhat,S1,3)  ! dldd = S1*Nhat (scaling)

	  if (.false.) then
	   print*,'** Before computing pathlength partials: '
	   print*,'** rhat = ', rhat
	   print*,'** ihat = ', ihat
	   print*,'** Ndoti = ', Ndoti
	   print*,'** S1 = ', S1
	   print*,'** Nhat = ', Nhat
	   print*,'** dddth matrix = ', dddth
	   print*,'** dldd = ', dldd
	  end if

	  CALL DMPROD(dldth,dldd,dddth,1,3,3)

	  !print*,'** dldth = ', dldth

C    Combined partial matrices

	  DO 4 i=1,3
	    ip3=i+3
	    DO 3 j=1,3
	      jp3=j+3
	      dxidxim1(i,j)=drdi(i,j)
	      dxidxim1(i,jp3)=drda(i,j)
	      dxidxim1(ip3,j)=dgdi(i,j)
	      dxidxim1(ip3,jp3)=dgda(i,j)
	      dxidui(i,j)=drdth(i,j)
	      dxidui(i,jp3)=drdd(i,j)
	      dxidui(ip3,j)=dgdth(i,j)
	      dxidui(ip3,jp3)=dgdd(i,j)
    3	    CONTINUE
	    dxidxim1(i,7)=0e0
	    dxidxim1(ip3,7)=0e0
	    dxidxim1(7,i)=0E0
	    dxidxim1(7,ip3)=0E0
	    dxidui(7,i)=dldth(i)
	    dxidui(7,ip3)=dldd(i)
    4	  CONTINUE
	  dxidxim1(7,7)=1d0
	END IF


C  Interpolated surface partials
C	IF (ifDMBuild.AND.(SrfType.EQ.5)) THEN
C	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO

 	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	        psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &          xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))

! 	  LRT=ChkRayTrans_TestOnly(ifBuild,ifPol,ifLNsrf,
!    &          psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
!    &          xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray),pv)

	  IF (.NOT.LRT) Eout(1:3)=0d0
	ELSE
	  LRT=.TRUE.
	END IF

	IF ((SrfType.EQ.5).AND.(IERROR.NE.0)) LRT=.FALSE.

C  Set output values for ray position, direction
	rout(1:3)=rhat(1:3)
	pout(1:3)=pr(1:3)

	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE Reflector

C***********************************************************************

! Adding Graded Index Refractor ...
! Starting 02/18/2011

	SUBROUTINE Refractor(ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,
     &   UDSrfCoef,UDSrfType,UDSrfParam,
     &   nGridMat,mGridMat,GridMat,GridSrfdx,
     &   psi,pv,prot,pin,ihat,Evec,
     &	 rad,pmon,xhat,yhat,zhat,
     &	 na,nb,ifnc,nc,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &	 XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &   nCoat,CoatW,nb_arr,kxb_arr,
     &   IsVarIndRef,IndRefStr,IsGrinRef,thisGrinRef)

	 use Constants, only: PI, TWOPI
	 use traceutil_mod
	 use didesub, only: IntSrf
         use elt_mod, only: mMonCoef,mAnaCoef

	 IMPLICIT NONE

	 Interface
           Function getIndexGrad_grin5(pt,pv,zdir,mc)
           Implicit none
           Real*8, dimension(3) :: getIndexGrad_grin5,pt,pv,zdir
           Real*8, dimension(10) :: mc
	   End Function getIndexGrad_grin5
	 End Interface

C  Compute refraction point, normal and transition matrices
C  for a general lens surface


	INTEGER :: nGridMat,mGridMat
        REAL*8 :: GridMat(mGridMat,mGridMat),GridSrfdx,getIndex_grin5

	LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &	        LRT,IsVarIndRef,IsGrinRef,ifnc

	REAL*8, dimension(:) :: CoatW
        REAL*8, dimension(0:) :: nb_arr,kxb_arr
	Type(IndRefGridStr) :: IndRefStr
        Type(GrinRefStr) :: thisGrinRef

	INTEGER :: i,j,k,k1,ip3,jp3,SrfType,UDSrfType,
     &	  NDP,IWORK(*),IERROR,mDP,
     &	  ApType,nObs,ObsType(nObs),nCoat,ic,segGap,totStepsGuess,
     &    xoff,yoff,xseg,yseg,ClosestInt

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 :: pv(3),prot(3),pin(3),ihat(3),psi(3),na,nb,nc,
     &   MonCoef(mMonCoef),pr(3),rhat(3),L,Nhat(3),lambda,pout(3),rout(3),
     &	 AnaCoef(mAnaCoef),UDSrfCoef(*),XYZSRF(3,NDP),DZSRF(2,NDP),DWORK(*),
     &	 rho(3),wVec(3),xGrid(3),yGrid(3),zGrid(3),ApVec(6),
     &   ObsVec(6,nObs),UDSrfParam(6),dd,dpface(3),dt,thkZ,pcrnt(3),
     &   Tcrnt(3),
     &   A_(3),B_(3),C_(3),pcrnt_prev(3),zdir(3)
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 :: Nvec(3),Nmag,dNdp(3,3),dpdd(3,3),dNdd(3,3),dNdi(3,3),
     &   Ndoti,Niinv,Nouti(3,3),drdN(3,3),drdd(3,3),drdth(3,3),
     &   drdi(3,3),
     &   drda(3,3),dgda(3,3),dgdi(3,3),dgdd(3,3),dldd(3),n_loc,
     &   dprot(3),dddth(3,3),dgdth(3,3),dldth(3),Kc,Kr,kxa,kxb,
     &   NoidNdi(3,3),drdii(3,3),dNda(3,3),dlda(3),dldi(3)

	REAL*8 :: mu,mu2,Nx(3,3),NxNx(3,3),cvec(3),Nxix(3,3),
     &	 Perpr(3,3),D1(9),D2(9),D3(3,3),S1,S2,Ndi,
     &	 rad,pmon(3),xhat(3),yhat(3),zhat(3),muNdi,cfa,cfb,
     &	 svec(3),shat(3),prhat(3),pihat(3),sfa,
     &	 evpha(3),espha(3),eopha(3),dpha,dtr

	COMPLEX*16 :: Evec(3),Esrf(3),Eout(3),RP,RS,Epi,Esi,Epr,Esr,
     &		      C1,C2,NaCmplx,NbCmplx,ccfa,ccfb,cone,cmu,cNdi,
     &            ccfb_arr(0:max_coat_layers),RP1,RS1,TP,TS,TP1,TS1,
     &                DD1,DD2,RPnew,RSnew
        Logical, save :: saveRayPath=.TRUE.

C  Compute surface

	mu=na/nb; mu2=mu*mu
	ifnc=.false.

C  Compute incidence point of ray

	IF (SrfType.EQ.1) THEN
	  CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &	      	       pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.2) THEN
	  CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &	  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.3) THEN
	  CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
	ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
	  CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.5) THEN
	  CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &	  XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	  wVec,pr,L,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.6) THEN
	  CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.7) THEN
	  CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &	  UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,
     &	  pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)

        ! Added by J.Lou, 11-16-03
 	ELSE IF (SrfType.EQ.9) THEN
          CALL GridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,
     &      psi,pv,prot,pin,ihat,pr,L,
     &      nGridMat,mGridMat,GridMat,GridSrfdx,
     &      pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)

	ELSE IF (SrfType.EQ.10) THEN
	  CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	           prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)

        ELSE IF (SrfType.EQ.11) THEN
          CALL AsphGridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,
     &      psi,pv,prot,pin,
     &      ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef,
     &      nGridMat,mGridMat,GridMat,GridSrfdx,pmon,
     &      xhat,yhat,zhat)
	ELSE
	  WRITE(*,*)'2: Unknown surface type. Aborting trace.'
	  RETURN
	END IF

	! The surface routines above computes ray intersection point
        ! 'pr'; the code below computes the direction of the refracted
        ! ray. -jzlou

C  Return if ray missed surface

	IF (.NOT.LROK) RETURN

	Nsrfm=Nvec/Nmag

#if 1
	! Support inhomogeneous lens, but not a correct way to do it
        ! since refraction surface is no longer aligned with lens surface,
        ! feature added by jzlou, 02-15-07.
	If (IsVarIndRef) Then
	  !print*,'******* using variable index of refraction'
	  dpface=pr-pmon
	  xseg=ClosestInt(DDOTC(xhat,dpface)/IndRefStr%width)
	  yseg=ClosestInt(DDOTC(yhat,dpface)/IndRefStr%width)
	  xoff=(size(IndRefStr%IndRefGrid,1)-1)/2+1
	  yoff=(size(IndRefStr%IndRefGrid,2)-1)/2+1
          ! xoff & yoff give half length of surface size, xseg & yseg give
          ! offset from center of surface (+ or - value), that's why xseg must
          ! be smaller than xoff, and yseg must be smaller than yoff.
	  If (abs(xseg)<xoff.AND.abs(yseg)<yoff) Then
	    nc=IndRefStr%IndRefGrid(xoff+xseg,yoff+yseg)
	    ifnc=.true.
	    mu=na/nc
	    mu2=mu*mu
	  End If
	End If
#endif

#if 1
	If (IsGrinRef) Then
          ! Implement graded-index media (lens) here, using ray eqns for inhomogeneous
          ! media ! and based on the numerical procedure in A. Sharma's paper in Applied
          ! Physics / Vol. 21, No. 6, March 15, 1982, -jzlou
          ! Assumed input at the front surface of the lens:
          ! 1) Lens type, which dictates a specific refractive index model.
	  ! 2) Lens depth, which determines how far the integration of ray path
          !    proceeds. Lens depth is measured along the (negative) direction of
	  !    surface normal.
	  ! 3) dt, the size of numerical integration step.
          ! Ray path integration proceeds until the projection of ray position along
          ! front surface normal is within 'delta' distance to the end of lens depth.
          ! At that point, the final ray position and direction is returned from
          ! the defractor routine to go on ray tracing to the next optical surface.
          ! jzlou, 02/2011.
#if 0
	  ! Test only
	  print*,'Refractor: thisGrinRef%modelType =',thisGrinRef%modelType
          print*,'Refractor: thisGrinRef%modelCoef(1:8) =',
     &           thisGrinRef%modelCoef(1:8)
          print*,'Refractor: thisGrinRef%lensDepth =',
     &           thisGrinRef%lensDepth
	  pcrnt=pr
	  A_(1:3)=getIndexGrad_grin5(pcrnt,thisGrinRef%modelCoef)
	  print*,'return from grad func =',A_(1:3)
	  stop
#endif

          If (thisGrinRef%modelType==5) Then
	    ! one of the graded-index lens supported by Zemax
            L = L*na  ! ray segment OPL from previous surface to this surface
	    dt=thisGrinRef%dt ! integration path step size
	    thkZ=thisGrinRef%lensDepth

	    ! starting ray point is 'pr', computed by one of the above surface routines,
	    ! starting ray direction should be 'ihat'.
	    If (DDOTC(ihat,psi)<0) Then
	      zdir=-psi;
	    Else
	      zdir=psi;
	    End If
            pcrnt=pr
            n_loc=getIndex_grin5(pcrnt,pv,zdir,thisGrinRef%modelCoef)
	    totStepsGuess=INT(thkZ/(dt*n_loc)*1.2) ! estimated total integration steps w/ margin
            If (totStepsGuess > mGradIndRaySegs) Then
	      segGap=INT(totStepsGuess/mGradIndRaySegs)+1
	    Else
              segGap=1
	    End If
	    ! Tcrnt below is intial optical ray vector
            Tcrnt(1:3)=n_loc*ihat
	    k=0; k1=0  ! for saving ray path only
	    Do while (abs(DDOTC(pcrnt,psi))<=(thkZ-dt))
	    !Do while (abs(DDOTC(pcrnt,psi))<thkZ)
              A_(1:3)=
     &        dt*getIndexGrad_grin5(pcrnt,pv,zdir,thisGrinRef%modelCoef)
              B_(1:3)=
     &          dt*getIndexGrad_grin5(pcrnt+dt/2d0*Tcrnt+dt/8d0*A_,
     &          pv,zdir,thisGrinRef%modelCoef)
	      C_(1:3)=dt*getIndexGrad_grin5(pcrnt+dt*Tcrnt+dt/2d0*B_,
     &                                 pv,zdir,thisGrinRef%modelCoef)
	      pcrnt_prev(1:3)=pcrnt(1:3)
	      pcrnt(1:3)=pcrnt_prev(1:3)+dt*(Tcrnt+1d0/6d0*(A_+2*B_))
              Tcrnt(1:3)=Tcrnt(1:3)+1d0/6d0*(A_+4d0*B_+C_)
              if (saveRayPath .and. iDrawRay_global>0 .and.
     &            (mod(k1,segGap)==0)) then
	        k=k+1
           gradIndRayPath(1:3,k,iDrawRay_global,iElt2iGradInd(crnt_Elt))
     &            =pcrnt(1:3)
!               gradIndRayPath(1:3,k,iDrawRay_global,1)=pcrnt(1:3)
	      end if
	      ! need to update ray path length
	      n_loc=0.5*(getIndex_grin5(pcrnt_prev,pv,zdir,
     &          thisGrinRef%modelCoef)+
     &          getIndex_grin5(pcrnt,pv,zdir,thisGrinRef%modelCoef))
	      L=L+dt*n_loc*n_loc  ! since ds = dt*n_loc, optical ds = ds*n_loc
              k1=k1+1
            End Do ! R-K integration loop
	    gradInd_kmax(iDrawRay_global,iElt2iGradInd(crnt_Elt))=k
!	    print*,'crnt_traced_ray =',crnt_traced_ray,
!    &        ', gradInd_kmax(crnt_traced_ray,1) =',gradInd_kmax(crnt_traced_ray,1)
          Else
            ! other lens models -- to implement
          End If
          !
          ! Now the last point on integrated ray path becomes 'pr', and the ray
          ! direction at 'pr' becomes rhat.
          pr(1:3)=pcrnt(1:3)
          rhat(1:3)=Tcrnt(1:3)/getIndex_grin5(pcrnt,
     &              pv,zdir,thisGrinRef%modelCoef)
	  Call DUNITIZE(rhat)
          Go to 110
        End If
#endif


C  Compute refracted ray (with lens of homogeneous media)

	Ndi=ihat(1)*Nhat(1)+ihat(2)*Nhat(2)+ihat(3)*Nhat(3)
	muNdi=mu*Ndi
	DO i=1,3
	  cvec(i)=Ndi*Nhat(i)-ihat(i)
   	END DO

C  Check for total internal reflection;
C  possibly occurs only when from dense to loose media -jzlou
	S2=1d0-mu2+muNdi*muNdi
	IF (S2.LT.0d0) THEN
	  LROK=.FALSE.
	  RETURN
	END IF

C  Compute direction of refracted ray -jzlou
	cfb=DSQRT(S2)
	S1=0d0
	DO i=1,3
	  rhat(i)=-mu*cvec(i)-cfb*Nhat(i)
	  S1=S1+rhat(i)*rhat(i)
        END DO
	S2=DSQRT(S1)
	rhat(1:3)=rhat(1:3)/S2

 110    Continue  ! used to skip computing refracted ray in homogeneous media
	          ! as done above -jzlou, 02/2011


C  Polarization calcs - only then coating effects could be included

	IF (ifPol) THEN
	  NaCmplx=DCMPLX(na,(-kxa))
	  NbCmplx=DCMPLX(nb,(-kxb))
	  cmu=NaCmplx/NbCmplx
	  cone=1d0

C    Propagation through isotropic medium a to surface

c	  S1=MOD((TwoPi*L/lambda),TwoPi)
	  S1=-TWOPI*L/lambda
	  C2=DCMPLX(0d0,S1)*NaCmplx
	  C1=CDEXP(C2)
	  Esrf(1)=C1*Evec(1)
	  Esrf(2)=C1*Evec(2)
	  Esrf(3)=C1*Evec(3)
cdebug
	  dtr=TWOPI/360d0
	  dpha=datan2(DIMAG(c1),DBLE(c1))/dtr
	  do i=1,1
	    evpha(i)=datan2(DIMAG(Evec(i)),DBLE(Evec(i)))/dtr
	    espha(i)=datan2(DIMAG(Esrf(i)),DBLE(Esrf(i)))/dtr
	  end do
cdebug

C    Surface transmittance effects
	  cfa=ABS(Ndi)
	  ccfa=cfa
	  CALL DXPROD(svec,ihat,Nhat)
	  sfa=DSQRT(svec(1)*svec(1)+svec(2)*svec(2)+svec(3)*svec(3))
	  IF (DABS(sfa).LT.1d-10) THEN
	    IF (ihat(2).NE.0d0) THEN
	      svec(1)=-ihat(2)
	      svec(2)=ihat(3)
	      svec(3)=ihat(1)
	    ELSE IF (ihat(3).NE.0d0) THEN
	      svec(1)=ihat(2)
	      svec(2)=-ihat(3)
	      svec(3)=ihat(1)
	    ELSE IF (ihat(1).NE.0d0) THEN
	      svec(1)=ihat(2)
	      svec(2)=ihat(3)
	      svec(3)=-ihat(1)
	    END IF
	    CALL DXPROD(shat,ihat,svec)
	    CALL DUNITIZE(shat)
	  ELSE
	    shat(1)=svec(1)/sfa
	    shat(2)=svec(2)/sfa
	    shat(3)=svec(3)/sfa
	  END IF
	  CALL DXPROD(pihat,shat,ihat)
	  CALL DXPROD(prhat,shat,rhat)

	  Epi=DCMPLX(pihat(1),0d0)*Esrf(1)+DCMPLX(pihat(2),0d0)*Esrf(2)
     &	      +DCMPLX(pihat(3),0d0)*Esrf(3)
	  Esi=DCMPLX(shat(1),0d0)*Esrf(1)+DCMPLX(shat(2),0d0)*Esrf(2)
     &	      +DCMPLX(shat(3),0d0)*Esrf(3)
	End If

	! Applicable to non-coated surface only
	IF (ifPol .and. nCoat==0) THEN
C    Compute refracted ray cosine
          cNdi=Ndi
          C1=cone-cmu*cmu*(cone-cNdi*cNdi)
          ccfb=CDSQRT(C1)
          cfb=ccfb

	  S1=2d0*DSQRT(cfb/cfa/mu)
	  C1=DCMPLX(S1,0d0)*ccfa*NaCmplx
	  TP=C1/(NbCmplx*ccfa+NaCmplx*ccfb)
	  TS=C1/(NaCmplx*ccfa+NbCmplx*ccfb)
	END IF

	! Applicable to surface with thin-film coating layers
        IF (ifPol .and. nCoat/=0) THEN
	  ! First calc refraction angles for each layer, from outmost inward
          NaCmplx=DCMPLX(nb_arr(0),(-kxb_arr(0)))
          NbCmplx=DCMPLX(nb_arr(1),(-kxb_arr(1)))

	  ccfb_arr(0)=DDOTC(ihat,Nhat) ! cosine of incident ray angle
          cmu=NaCmplx/NbCmplx
          cNdi=ihat(1)*Nhat(1)+ihat(2)*Nhat(2)+ihat(3)*Nhat(3)
          C1=cone-cmu*cmu*(cone-cNdi*cNdi)
          ccfb_arr(1)=CDSQRT(C1) ! cosine of refraction angle for outmost interface
          !
          Do ic=1,nCoat
            NaCmplx=DCMPLX(nb_arr(ic),(-kxb_arr(ic)))
            NbCmplx=DCMPLX(nb_arr(ic+1),(-kxb_arr(ic+1)))
            cmu=NaCmplx/NbCmplx
            C1=cone-cmu*cmu*(cone-ccfb_arr(ic)*ccfb_arr(ic))
            ccfb_arr(ic+1)=CDSQRT(C1)
          End Do

	  ! Now calc P and S reflectances recursively for each interface,
          ! innermost interface first, then for each outer coating layer

          ! For innermost (substrate) interface, just use the Frasnel
          ! reflectance eqn
          NaCmplx=DCMPLX(nb_arr(nCoat),(-kxb_arr(nCoat)))
          NbCmplx=DCMPLX(nb_arr(nCoat+1),(-kxb_arr(nCoat+1)))
          RP=(NaCmplx*ccfb_arr(nCoat+1)-NbCmplx*ccfb_arr(nCoat))
     &       /(NbCmplx*ccfb_arr(nCoat)+NaCmplx*ccfb_arr(nCoat+1))
          RS=(NaCmplx*ccfb_arr(nCoat)-NbCmplx*ccfb_arr(nCoat+1))
     &       /(NaCmplx*ccfb_arr(nCoat)+NbCmplx*ccfb_arr(nCoat+1))

          ! Add transmittance
          TS=2d0*NaCmplx*ccfb_arr(nCoat)
     &       /(NaCmplx*ccfb_arr(nCoat)+NbCmplx*ccfb_arr(nCoat+1))
          TP=2d0*NaCmplx*ccfb_arr(nCoat)
     &       /(NaCmplx*ccfb_arr(nCoat+1)+NbCmplx*ccfb_arr(nCoat))

	  ! For each coating layer going outward, use thin film
          ! reflectance eqn,
          Do ic=nCoat,1,-1
            NaCmplx=DCMPLX(nb_arr(ic-1),(-kxb_arr(ic-1)))
            NbCmplx=DCMPLX(nb_arr(ic),(-kxb_arr(ic)))

            ! Calc phase delay due to pass thru coating
            !dd=CoatW(ic)/ccfb_arr(ic) ! ray pathlength in this layer
            !C2=DCMPLX(0d0,1d0)*TwoPi*2d0*dd*NbCmplx/lambda
            dd=CoatW(ic)*ccfb_arr(ic)
            C2=DCMPLX(0d0,1d0)*2d0*(TWOPI/lambda)*dd*NbCmplx
            C1=CDEXP(-C2)

            ! Reflectance
            RP1=(NaCmplx*ccfb_arr(ic)-NbCmplx*ccfb_arr(ic-1))
     &          /(NbCmplx*ccfb_arr(ic-1)+NaCmplx*ccfb_arr(ic))
            !RP1=-RP1 ! testing
            RPnew=(RP1+RP*C1)/(1d0+RP1*RP*C1)

            RS1=(NaCmplx*ccfb_arr(ic-1)-NbCmplx*ccfb_arr(ic))
     &          /(NaCmplx*ccfb_arr(ic-1)+NbCmplx*ccfb_arr(ic))
            RSnew=(RS1+RS*C1)/(1d0+RS1*RS*C1)

            ! Transmittance
            DD2=DCMPLX(0d0,1d0)*(TWOPI/lambda)*dd*NbCmplx
            DD1=CDEXP(-DD2)
            TP1= 2d0*NaCmplx*ccfb_arr(ic-1)
     &          /(NbCmplx*ccfb_arr(ic-1)+NaCmplx*ccfb_arr(ic))
            TP= (TP1*TP*DD1)/(1d0+RP1*RP*C1)
            TS1= 2d0*NaCmplx*ccfb_arr(ic-1)
     &          /(NaCmplx*ccfb_arr(ic-1)+NbCmplx*ccfb_arr(ic))
            TS= (TS1*TS*DD1)/(1d0+RS1*RS*C1)

            RP=RPnew; RS=RSnew
	  End Do
	END IF

	IF (ifPol) THEN
	  Epr=TP*Epi
	  Esr=TS*Esi
	  Eout(1)=DCMPLX(prhat(1),0d0)*Epr+DCMPLX(shat(1),0d0)*Esr
	  Eout(2)=DCMPLX(prhat(2),0d0)*Epr+DCMPLX(shat(2),0d0)*Esr
	  Eout(3)=DCMPLX(prhat(3),0d0)*Epr+DCMPLX(shat(3),0d0)*Esr
cdebug
	  DO i=1,1
	    eopha(i)=datan2(DIMAG(Eout(i)),DBLE(Eout(i)))/dtr
	  END DO
cdebug
	END IF

C  Partials

	IF (ifBuild) THEN
	  S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	  IF (S1.GT.0d0) THEN
	    CALL DNEGATE(dNdp,dNdp,9)
	  END IF

C    Lens partials

	  CALL DXMAT(Nx,Nhat)
	  CALL DMPROD(NxNx,Nx,Nx,3,3,3)

	  CALL DXMAT(D1,ihat)
	  CALL DMPROD(Nxix,Nx,D1,3,3,3)
	  CALL DMPROD(D3,ihat,Nxix,1,3,3)
	  S1=mu2/cfb
	  CALL DSMPROD(D2,Nhat,S1,3)
	  CALL DMPROD(D1,D2,D3,3,1,3)
	  D1(1)=D1(1)-cfb
	  D1(5)=D1(5)-cfb
	  D1(9)=D1(9)-cfb
	  CALL DMPROD(D2,Nx,ihat,3,3,1)
	  CALL DXMAT(D3,D2)
	  CALL DADD(D2,D3,Nxix,9)
	  CALL DSMPROD(D3,D2,mu,9)
	  CALL DADD(drdN,D1,D3,9)

	  CALL DMPROD(D1,ihat,NxNx,1,3,3)
	  CALL DSMPROD(D2,Nhat,(-S1),3)
	  CALL DMPROD(D3,D2,D1,3,1,3)
	  CALL DSMPROD(D1,NxNx,(-mu),9)
	  CALL DADD(drdii,D1,D3,9)

C    Partials of the normal

	  CALL DOUTER(Nouti,ihat,Nhat)
	  Ndoti=DDOTC(Nhat,ihat)
	  Niinv=1d0/Ndoti
	  CALL DSMPROD(NoidNdi,Nouti,Niinv,9)
	  CALL DEQUATE(dpdd,NoidNdi,9)
	  dpdd(1,1)=dpdd(1,1)-1d0
	  dpdd(2,2)=dpdd(2,2)-1d0
	  dpdd(3,3)=dpdd(3,3)-1d0
	  CALL DMPROD(dNdd,dNdp,dpdd,3,3,3)
	  S1=-L
	  CALL DSMPROD(dNdi,dNdd,S1,9)
	  CALL DNEGATE(dNda,dNdd,9)

C    Rotation-to-translation partials

	  CALL DSUB(dprot,prot,pr,3)
	  CALL DXMAT(dddth,dprot)

C    Ray direction partials

	  CALL DMPROD(drdd,drdN,dNdd,3,3,3)
	  CALL DMPROD(D1,drdN,dNdi,3,3,3)
	  CALL DADD(drdi,drdii,D1,9)
	  CALL DMPROD(D1,drdN,Nx,3,3,3)
	  CALL DMPROD(D2,drdd,dddth,3,3,3)
	  CALL DSUB(drdth,D2,D1,9)
	  CALL DMPROD(drda,drdN,dNda,3,3,3)

C    Beamwalk partials

	  CALL DPROJECT(Perpr,rhat)
	  CALL DMPROD(D1,Perpr,dpdd,3,3,3)
	  CALL DNEGATE(dgda,D1,9)
	  CALL DSMPROD(dgdi,dgda,L,9)
	  S1=1d0/Ndoti
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DMPROD(D2,Perpr,ihat,3,3,1)
	  CALL DMPROD(dgdd,D2,D1,3,1,3)
	  CALL DMPROD(dgdth,dgdd,dddth,3,3,3)

C    Pathlength partials

	  CALL DSUB(D3,rhat,ihat,3)
	  CALL DMPROD(dlda,D3,dpdd,1,3,3)
	  CALL DSMPROD(dldi,dlda,L,3)
	  S1=-DDOTC(D3,ihat)/Ndoti
	  CALL DSMPROD(dldd,Nhat,S1,3)
	  CALL DMPROD(dldth,dldd,dddth,1,3,3)

C    Combined partial matrices

	  DO 7 i=1,3
	    ip3=i+3
	    DO 6 j=1,3
	      jp3=j+3
	      dxidxim1(i,j)=drdi(i,j)
	      dxidxim1(i,jp3)=drda(i,j)
	      dxidxim1(ip3,j)=dgdi(i,j)
	      dxidxim1(ip3,jp3)=dgda(i,j)
	      dxidui(i,j)=drdth(i,j)
	      dxidui(i,jp3)=drdd(i,j)
	      dxidui(ip3,j)=dgdth(i,j)
	      dxidui(ip3,jp3)=dgdd(i,j)
    6	    CONTINUE
	    dxidxim1(i,7)=0e0
	    dxidxim1(ip3,7)=0e0
	    dxidxim1(7,i)=dldi(i)
	    dxidxim1(7,ip3)=dlda(i)
	    dxidui(7,i)=dldth(i)
	    dxidui(7,ip3)=dldd(i)
    7	  CONTINUE
	  dxidxim1(7,7)=1d0
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO

	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,psi,rho,
     &	         ApType,ApVec,nObs,ObsType,
     &           ObsVec,xGrid,yGrid,zGrid,
     &           taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) Eout(1:3)=0d0
	ELSE
	  LRT=.TRUE.
	END IF

	rout(1:3)=rhat(1:3)
	pout(1:3)=pr(1:3)

	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE Refractor

C***********************************************************************

	SUBROUTINE NSRefractor(ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,UDSrfCoef,UDSrfType,
     &	 UDSrfParam,
     &   nGridMat,mGridMat,GridMat,GridSrfdx,
     &   psi,pv,prot,pin,ihat,Evec,
     &	 rad,pmon,xhat,yhat,zhat,
     &	 na,nb,kxa,kxb,lambda,
     &	 nCoat,CoatIndx,CoatThk,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &	 XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

         use Constants, only: PI, TWOPI
	 use traceutil_mod
	 use didesub, only: IntSrf
         use elt_mod, only: mMonCoef,mAnaCoef

C  Compute refraction point, normal and transition matrices
C  for a general lens surface

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &	LRT

	INTEGER i,j,ip3,jp3,SrfType,nCoat,UDSrfType,
     &	  NDP,IWORK(*),IERROR,mDP,
     &	  ApType,nObs,ObsType(nObs)

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	INTEGER nGridMat,mGridMat
        REAL*8 GridMat(mGridMat,mGridMat),GridSrfdx

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),na,nb,
     &   MonCoef(mMonCoef),pr(3),rhat(3),L,Nhat(3),lambda,pout(3),rout(3),
     &	 AnaCoef(mAnaCoef),UDSrfCoef(*),XYZSRF(3,NDP),DZSRF(2,NDP),
     &	 rho(3),DWORK(*),wVec(3),xGrid(3),yGrid(3),zGrid(3),ApVec(6),
     &   ObsVec(6,nObs),UDSrfParam(6)
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 Nvec(3),Nmag,dNdp(3,3),dpdd(3,3),dNdd(3,3),dNdi(3,3),
     &  Ndoti,Niinv,Nouti(3,3),drdN(3,3),drdd(3,3),drdth(3,3),drdi(3,3),
     &  drda(3,3),dgda(3,3),dgdi(3,3),dgdd(3,3),dldd(3),
     &  dprot(3),dddth(3,3),dgdth(3,3),dldth(3),Kc,Kr,kxa,kxb,
     &  NoidNdi(3,3),drdii(3,3),dNda(3,3),dlda(3),dldi(3),
     &	CoatIndx(20),CoatThk(20)

	REAL*8 mu,mu2,Nx(3,3),NxNx(3,3),cvec(3),Nxix(3,3),
     &	Perpr(3,3),D1(9),D2(9),D3(3,3),S1,S2,Ndi,
     &	rad,pmon(3),xhat(3),yhat(3),zhat(3),muNdi,cfa,cfb,
     &	svec(3),shat(3),prhat(3),pihat(3),sfa

	COMPLEX*16 Evec(3),Esrf(3),Eout(3),TP,TS,Epi,Esi,Epr,Esr,
     &	C1,C2,NaCmplx,NbCmplx,ccfa,ccfb,cone,cmu,cNdi

C  Compute surface

	mu=na/nb
	mu2=mu*mu

C  Compute incidence point of ray

	IF (SrfType.EQ.1) THEN
	  CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &	  pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.2) THEN
	  CALL NSConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &	  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.3) THEN
	  CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
	ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
	  CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.5) THEN
	  CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &	  XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	  wVec,pr,L,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.6) THEN
	  CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.7) THEN
	  CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &	  UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,
     &	  pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.9) THEN
          CALL GridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,
     &    psi,pv,prot,pin,ihat,pr,L,
     &    nGridMat,mGridMat,GridMat,GridSrfdx,
     &    pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.10) THEN
	  CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.11) THEN
          CALL AsphGridSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,
     &      psi,pv,prot,pin,
     &      ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef,
     &      nGridMat,mGridMat,GridMat,GridSrfdx,pmon,
     &      xhat,yhat,zhat)
	ELSE
	  WRITE(*,*)'3: Unknown surface type. Aborting trace.'
	END IF

	Nsrfm=Nvec

C  Return if ray missed surface

	IF (.NOT.LROK) RETURN

	Nsrfm=Nvec/Nmag

C  Compute refracted ray

	Ndi=ihat(1)*Nhat(1)+ihat(2)*Nhat(2)+ihat(3)*Nhat(3)
	muNdi=mu*Ndi
	DO 1 i=1,3
	  cvec(i)=Ndi*Nhat(i)-ihat(i)
  1	CONTINUE

C  Check for total internal reflection

	S2=1d0-mu2+muNdi*muNdi
	IF (S2.LT.0d0) THEN
	  LROK=.FALSE.
	  RETURN
	END IF
	cfb=DSQRT(S2)
	S1=0d0
	DO 2 i=1,3
	  rhat(i)=-mu*cvec(i)-cfb*Nhat(i)
	  S1=S1+rhat(i)*rhat(i)
  2	CONTINUE
	S2=DSQRT(S1)
	DO 3 i=1,3
	  rhat(i)=rhat(i)/S2
  3	CONTINUE

C  Polarization calcs

	IF (ifPol) THEN
	  NaCmplx=DCMPLX(na,(-kxa))
	  NbCmplx=DCMPLX(nb,(-kxb))
	  cmu=NaCmplx/NbCmplx
	  cone=1d0

C    Propagation through isotropic medium a to surface

c	  S1=MOD((TwoPi*L/lambda),TwoPi)
	  S1=-TWOPI*L/lambda
	  C2=DCMPLX(0d0,S1)*NaCmplx
	  C1=CDEXP(C2)
	  Esrf(1)=C1*Evec(1)
	  Esrf(2)=C1*Evec(2)
	  Esrf(3)=C1*Evec(3)

C    Compute refracted ray cosine

	  cNdi=Ndi
	  C1=cone-cmu*cmu*(cone-cNdi*cNdi)
	  ccfb=CDSQRT(C1)
	  cfb=ccfb

C    Surface transmittance effects

	  cfa=ABS(Ndi)
	  ccfa=cfa
	  CALL DXPROD(svec,ihat,Nhat)
	  sfa=DSQRT(svec(1)*svec(1)+svec(2)*svec(2)+svec(3)*svec(3))
	  IF (DABS(sfa).LT.1d-10) THEN
	    IF (ihat(2).NE.0d0) THEN
	      svec(1)=-ihat(2)
	      svec(2)=ihat(3)
	      svec(3)=ihat(1)
	    ELSE IF (ihat(3).NE.0d0) THEN
	      svec(1)=ihat(2)
	      svec(2)=-ihat(3)
	      svec(3)=ihat(1)
	    ELSE IF (ihat(1).NE.0d0) THEN
	      svec(1)=ihat(2)
	      svec(2)=ihat(3)
	      svec(3)=-ihat(1)
	    END IF
	    CALL DXPROD(shat,ihat,svec)
	    CALL DUNITIZE(shat)
	  ELSE
	    shat(1)=svec(1)/sfa
	    shat(2)=svec(2)/sfa
	    shat(3)=svec(3)/sfa
	  END IF
	  CALL DXPROD(pihat,shat,ihat)
	  CALL DXPROD(prhat,shat,rhat)

	  Epi=DCMPLX(pihat(1),0d0)*Esrf(1)+DCMPLX(pihat(2),0d0)*Esrf(2)+
     &	  DCMPLX(pihat(3),0d0)*Esrf(3)
	  Esi=DCMPLX(shat(1),0d0)*Esrf(1)+DCMPLX(shat(2),0d0)*Esrf(2)+
     &	  DCMPLX(shat(3),0d0)*Esrf(3)

	  IF (nCoat.GT.0) THEN
	    CALL AirGap(TP,TS,cfa,sfa,cfb,na,nb,mu,lambda,nCoat,
     &	    CoatIndx,CoatThk)
	  ELSE
	    S1=2d0*DSQRT(cfb/cfa/mu)
	    C1=DCMPLX(S1,0d0)*ccfa*NaCmplx
	    TP=C1/(NbCmplx*ccfa+NaCmplx*ccfb)
	    TS=C1/(NaCmplx*ccfa+NbCmplx*ccfb)
	  END IF

	  Epr=TP*Epi
	  Esr=TS*Esi
	  Eout(1)=DCMPLX(prhat(1),0d0)*Epr+DCMPLX(shat(1),0d0)*Esr
	  Eout(2)=DCMPLX(prhat(2),0d0)*Epr+DCMPLX(shat(2),0d0)*Esr
	  Eout(3)=DCMPLX(prhat(3),0d0)*Epr+DCMPLX(shat(3),0d0)*Esr
	END IF

C  Partials

	IF (ifBuild) THEN
	  S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	  IF (S1.GT.0d0) THEN
	    CALL DNEGATE(dNdp,dNdp,9)
	  END IF

C    Lens partials

	  CALL DXMAT(Nx,Nhat)
	  CALL DMPROD(NxNx,Nx,Nx,3,3,3)

	  CALL DXMAT(D1,ihat)
	  CALL DMPROD(Nxix,Nx,D1,3,3,3)
	  CALL DMPROD(D3,ihat,Nxix,1,3,3)
	  S1=mu2/cfb
	  CALL DSMPROD(D2,Nhat,S1,3)
	  CALL DMPROD(D1,D2,D3,3,1,3)
	  D1(1)=D1(1)-cfb
	  D1(5)=D1(5)-cfb
	  D1(9)=D1(9)-cfb
	  CALL DMPROD(D2,Nx,ihat,3,3,1)
	  CALL DXMAT(D3,D2)
	  CALL DADD(D2,D3,Nxix,9)
	  CALL DSMPROD(D3,D2,mu,9)
	  CALL DADD(drdN,D1,D3,9)

	  CALL DMPROD(D1,ihat,NxNx,1,3,3)
	  CALL DSMPROD(D2,Nhat,(-S1),3)
	  CALL DMPROD(D3,D2,D1,3,1,3)
	  CALL DSMPROD(D1,NxNx,(-mu),9)
	  CALL DADD(drdii,D1,D3,9)

C    Partials of the normal

	  CALL DOUTER(Nouti,ihat,Nhat)
	  Ndoti=DDOTC(Nhat,ihat)
	  Niinv=1d0/Ndoti
	  CALL DSMPROD(NoidNdi,Nouti,Niinv,9)
	  CALL DEQUATE(dpdd,NoidNdi,9)
	  dpdd(1,1)=dpdd(1,1)-1d0
	  dpdd(2,2)=dpdd(2,2)-1d0
	  dpdd(3,3)=dpdd(3,3)-1d0
	  CALL DMPROD(dNdd,dNdp,dpdd,3,3,3)
	  S1=-L
	  CALL DSMPROD(dNdi,dNdd,S1,9)
	  CALL DNEGATE(dNda,dNdd,9)

C    Rotation-to-translation partials

	  CALL DSUB(dprot,prot,pr,3)
	  CALL DXMAT(dddth,dprot)

C    Ray direction partials

	  CALL DMPROD(drdd,drdN,dNdd,3,3,3)
	  CALL DMPROD(D1,drdN,dNdi,3,3,3)
	  CALL DADD(drdi,drdii,D1,9)
	  CALL DMPROD(D1,drdN,Nx,3,3,3)
	  CALL DMPROD(D2,drdd,dddth,3,3,3)
	  CALL DSUB(drdth,D2,D1,9)
	  CALL DMPROD(drda,drdN,dNda,3,3,3)

C    Beamwalk partials

	  CALL DPROJECT(Perpr,rhat)
	  CALL DMPROD(D1,Perpr,dpdd,3,3,3)
	  CALL DNEGATE(dgda,D1,9)
	  CALL DSMPROD(dgdi,dgda,L,9)
	  S1=1d0/Ndoti
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DMPROD(D2,Perpr,ihat,3,3,1)
	  CALL DMPROD(dgdd,D2,D1,3,1,3)
	  CALL DMPROD(dgdth,dgdd,dddth,3,3,3)

C    Pathlength partials

	  CALL DSUB(D3,rhat,ihat,3)
	  CALL DMPROD(dlda,D3,dpdd,1,3,3)
	  CALL DSMPROD(dldi,dlda,L,3)
	  S1=-DDOTC(D3,ihat)/Ndoti
	  CALL DSMPROD(dldd,Nhat,S1,3)
	  CALL DMPROD(dldth,dldd,dddth,1,3,3)

C    Combined partial matrices

	  DO 7 i=1,3
	    ip3=i+3
	    DO 6 j=1,3
	      jp3=j+3
	      dxidxim1(i,j)=drdi(i,j)
	      dxidxim1(i,jp3)=drda(i,j)
	      dxidxim1(ip3,j)=dgdi(i,j)
	      dxidxim1(ip3,jp3)=dgda(i,j)
	      dxidui(i,j)=drdth(i,j)
	      dxidui(i,jp3)=drdd(i,j)
	      dxidui(ip3,j)=dgdth(i,j)
	      dxidui(ip3,jp3)=dgdd(i,j)
    6	    CONTINUE
	    dxidxim1(i,7)=0e0
	    dxidxim1(ip3,7)=0e0
	    dxidxim1(7,i)=dldi(i)
	    dxidxim1(7,ip3)=dlda(i)
	    dxidui(7,i)=dldth(i)
	    dxidui(7,ip3)=dldd(i)
    7	  CONTINUE
	  dxidxim1(7,7)=1d0
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,psi,rho,
     &	         ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) Eout(1:3)=0d0
	ELSE
	  LRT=.TRUE.
	END IF

	DO 8 i=1,3
	  rout(i)=rhat(i)
	  pout(i)=pr(i)
   8	CONTINUE

	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE NSRefractor

C***********************************************************************

	SUBROUTINE AirGap(TP,TS,cfa,sfa,cfb,n1,n3,mu,Wavelen,nCoat,
     &	CoatIndx,CoatThk)
        USE Constants, ONLY: TWOPI
	IMPLICIT NONE
	INTEGER nCoat
	REAL*8 cfa,sfa,cfb,n1,n3,mu,Wavelen,CoatIndx(20),CoatThk(20),
     &	d,k,ky,k1x,k2x,k3x,ky2,k2x2,n2,sfac,RP,IP,A2,SA2,CA2,
     &	y1,sfbc
	COMPLEX*16 TP,TS,C1

	n2=1d0
	d=CoatThk(1)
	k=TWOPI/Wavelen
	ky=n1*k*sfa
	ky2=ky*ky
	k1x=n1*k*cfa
	k2x2=ABS(ky2-(n2*k)**2)
	k2x=DSQRT(k2x2)
	k3x=DSQRT(ABS((n3*k)**2-ky2))

	sfac=n2/n1
	sfbc=n3/n1
	A2=k2x*d
	IF (sfa.GT.sfac) THEN
	  IF (sfa.LE.sfbc) THEN
c	    y1=d*DTAN(DACOS(cfb))
	    y1=0d0
	    CA2=DCOSH(A2)
	    SA2=DSINH(A2)
	    C1=DCMPLX(0d0,DSQRT(n3*cfb/n1/cfa)*2d0*k1x*k2x)
     &	    *CDEXP(DCMPLX(0d0,ky*y1))
     &	    *CDEXP(DCMPLX(0d0,-n3*k*d/cfb))
	    RP=(k1x*k3x-k2x2)*SA2
	    IP=k2x*(k1x+k3x)*CA2
	    TP=C1/DCMPLX(RP,IP)
	    RP=(k1x*k3x*n2*n2/n1/n3-k2x2*n1*n3/n2/n2)*SA2
	    IP=k2x*(n3*k1x/n1+n1*k3x/n3)*CA2
	    TS=C1/DCMPLX(RP,IP)
	  ELSE
	    TP=0d0
	    TS=0d0
	  END IF
	ELSE
	  y1=d*DTAN(DASIN(sfa/sfac))
	  CA2=DCOS(A2)
	  SA2=DSIN(A2)
	  C1=DSQRT(n3*cfb/n1/cfa)*2d0*k1x*k2x
     &	  *CDEXP(DCMPLX(0d0,ky*y1))
     &	  *CDEXP(DCMPLX(0d0,-n3*k*d/cfb))

	  IP=-(k1x*k3x+k2x2)*SA2
	  RP=k2x*(k1x+k3x)*CA2
	  TP=C1/DCMPLX(RP,IP)

	  IP=-(k1x*k3x*n2*n2/n1/n3+k2x2*n1*n3/n2/n2)*SA2
	  RP=k2x*(n3*k1x/n1+n1*k3x/n3)*CA2
	  TS=C1/DCMPLX(RP,IP)
	END IF

	RETURN
	END SUBROUTINE AirGap

C***********************************************************************

	SUBROUTINE LensArray(ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,UDSrfCoef,UDSrfType,
     &	 UDSrfParam,psi,pv,prot,pin,ihat,Evec,
     &	 rad,pmon,xhat,yhat,zhat,
     &	 na,nb, ! na,nb = CurIndRef,IndRef(iElt) in calling routine
     &   ifnc,nc,
     &   kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &	 LensArrayType,width,
     &   IsVarIndRefLensArr,TotLenslets,thisElt,
     &	 XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

	!use traceutil_mod, only : LensletStrArray
	use traceutil_mod
        use elt_mod, only: mMonCoef,mAnaCoef

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifLNsrf,ifDMBuild,ifPol,ifObs,!ChkRayTrans,
     &          LRT,IsVarIndRefLensArr,ifnc

	INTEGER i,j,ip3,jp3,SrfType,LensArrayType,IXseg,IYseg,IXtri,
     &	  ILtri,IRtri,IYtri,LINT,UDSrfType,
     &	  NDP,IWORK(*),IERROR,mDP,
     &	  ApType,nObs,ObsType(nObs),TotLenslets,thisElt

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),na,nb,nc,
     &    MonCoef(mMonCoef),pr(3),rhat(3),L,Nhat(3),pout(3),rout(3),
     &	  kxa,kxb,lambda,rho(3),DWORK(*),AnaCoef(mAnaCoef),
     &	  UDSrfCoef(*),XYZSRF(3,NDP),DZSRF(2,NDP),xGrid(3),yGrid(3),
     &    zGrid(3),ApVec(6),ObsVec(6,nObs),UDSrfParam(6)
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 Nvec(3),Nmag,dNdp(3,3),dpdd(3,3),dNdd(3,3),dNdi(3,3),
     &    Ndoti,Niinv,Nouti(3,3),drdN(3,3),drdd(3,3),drdth(3,3),
     &    drdi(3,3),
     &    drda(3,3),dgda(3,3),dgdi(3,3),dgdd(3,3),dldd(3),
     &    dprot(3),dddth(3,3),dgdth(3,3),dldth(3),Kc,Kr,
     &    NoidNdi(3,3),drdii(3,3),dNda(3,3),dlda(3),dldi(3)

	REAL*8 mu,mu2,Nx(3,3),NxNx(3,3),cvec(3),cfb,Nxix(3,3),
     &	  Perpr(3,3),D1(9),D2(9),D3(3,3),S1,width,tan30,
     &    sin30,cos30,
     &	  rad,pmon(3),xhat(3),yhat(3),zhat(3),vptL(3),pface(3),
     &    dpface(3),xseg,yseg,Lhat(3),ARhat(3),w2,dd(2),ds

	COMPLEX*16 Evec(3),Eout(3)

	DATA
     &    tan30/5.773502691896257d-01/cos30/8.660254037844387d-01/
	DATA sin30/5d-01/

        ! nb,mu,mu2 will be updated below for LensArray
        ! with variable IndRef
	mu=na/nb
	mu2=mu*mu
	w2=0.5d0*width
        ifnc=.false.

C  Compute incidence point of ray on face plane

#if 0
	print*,'b4 L =',L
	print*,'b4: psi =', psi
	print*,'b4: Kc =', Kc
	print*,'b4: Kr =', Kr
	print*,'b4: pv =',  pv
	print*,'b4: prot =',prot
	print*,'b4: pin =',pin
	print*,'b4: ihat =',ihat
	print*,'b4: Nhat =',Nhat
	print*,'b4: Nvec =',Nvec
	print*,'b4 ifLNsrf =',ifLNsrf
	print*,'b4: Nmag =',Nmag
	print*,'b4: dNdp =',dNdp
	print*,'b4: LROK =',LROK
	print*,' '
#endif

	!ifLNsrf=.true.  ! added by jzlou
 	CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &		     pface,L,Nvec,Nmag,Nhat,dNdp)

#if 0
	! test only
	print*,'aft: psi =', psi
	print*,'aft: Kc =', Kc
	print*,'aft: Kr =', Kr
        print*,'aft: pv =',  pv
        print*,'aft: prot =',prot
        print*,'aft: pin =',pin
        print*,'aft: ihat =',ihat
        print*,'aft: Nhat =',Nhat
	print*,'aft: Nvec =',Nvec
	print*,'aft ifLNsrf =',ifLNsrf
	print*,'aft: Nmag =',Nmag
	print*,'aft: dNdp =',dNdp
	print*,'aft LROK =',LROK

	CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &              prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)

	print*,'LensArray: LROK =',LROK
        stop
#endif


C  Find vertex of appropriate lenslet

	DO i=1,3
	  !dpface(i)=pface(i)-pv(i)  ! original
	  dpface(i)=pface(i)-pmon(i) ! change by jzlou
  	END DO

#if 0
          if (IsVarIndRefLensArr) then
            print*,'iElt,TotLenslets =',thisElt,TotLenslets
            print*,'LensletStrArray%Xtri,Ltri,Rtri,IndRef ='
	    j=thisElt
	    Do i=1,TotLenslets
	      if (abs(LensArrayType).EQ.1) then
	   print*,LensletStrArray(i,j)%Xtri,LensletStrArray(i,j)%Ltri,
     &            LensletStrArray(i,j)%Rtri,LensletStrArray(i,j)%IndRef
	      else
           print*,LensletStrArray(i,j)%RecXtri,
     &            LensletStrArray(i,j)%RecYtri,
     &            LensletStrArray(i,j)%IndRef
	      end if
	    End Do
	    stop
	  end if
#endif

	IF (abs(LensArrayType).EQ.1) THEN
	  ! Hex lenslet array
	  DO i=1,3
	    Lhat(i)=sin30*xhat(i)+cos30*yhat(i)   ! L axis
	    ARhat(i)=-sin30*xhat(i)+cos30*yhat(i) ! R axis
  	  END DO

	  ILtri=LINT(DDOTC(Lhat,dpface)/w2)  ! coord value on L axis
	  IRtri=LINT(DDOTC(ARhat,dpface)/w2) ! coord value on R axis
	  IYtri=IRtri+ILtri
	  IXtri=LINT(DDOTC(xhat,dpface)/w2)  ! coord value on X axis

	  IYseg=LINT(DBLE(IYtri+2)/3d0)
	  IF (MOD(IYseg,2).EQ.0) THEN
	    IXseg=LINT(DBLE(IXtri+1)/2d0)
	    xseg=IXseg*width
	    yseg=1.5d0*DBLE(IYseg)*width*tan30 ! = ((sqrt(3)/2)*width)*IYseg
	  ELSE
	    IXseg=LINT(DBLE(IXtri)/2d0)
	    xseg=IXseg*width+w2
	    yseg=1.5d0*DBLE(IYseg)*width*tan30
	  END IF

	  ! added by jzlou
          If (IsVarIndRefLensArr) Then
	    !dsmin=1d22
	    Do i=1,TotLenslets
              dd(1)=LensletStrArray(i,thisElt)%Xctr
              dd(2)=LensletStrArray(i,thisElt)%Yctr
	      dd(1)=dd(1)-xseg
	      dd(2)=dd(2)-yseg
	      ds=sqrt(DDOTC2(dd,dd))
	      If (ds<w2/3) Then
	        nc=LensletStrArray(i,thisElt)%IndRef
                ifnc=.true.
                mu=na/nc; mu2=mu*mu
	        exit
              End If
	    End Do
	  End If

	ELSE
	  ! Rectangular lenslet array
	  IYtri=LINT(DDOTC(yhat,dpface)/w2)
	  IXtri=LINT(DDOTC(xhat,dpface)/w2)

	  IXseg=LINT(DBLE(IXtri+1)/2d0)
	  xseg=IXseg*width
	  IYseg=LINT(DBLE(IYtri+1)/2d0)
	  yseg=IYseg*width

	  ! added by jzlou
          If (IsVarIndRefLensArr) Then
            Do i=1,TotLenslets
	      dd(1)=LensletStrArray(i,thisElt)%Xctr
              dd(2)=LensletStrArray(i,thisElt)%Yctr
              dd(1)=dd(1)-xseg
              dd(2)=dd(2)-yseg
              ds=DDOTC2(dd,dd)
              If (ds<w2/3) Then
                nc=LensletStrArray(i,thisElt)%IndRef
	        ifnc=.true.
                mu=na/nc; mu2=mu*mu
              End If
	    End Do
	  End If
	END IF

	If (LensArrayType>0) Then
	  DO i=1,3
	    vptL(i)=xseg*xhat(i)+yseg*yhat(i)+pv(i)
  	  END DO
	Else
	  ! Use VptElt instead when LensArrayType < 0
	  vptL(1:3)=pv(1:3)
	End If

C  Compute incidence point of ray on lenslet

	IF (SrfType.EQ.2) THEN
	  CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,vptL,
     &	  	      prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.3) THEN
	  CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,vptL,
     &	  	     prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
	ELSE IF (SrfType.EQ.6) THEN
	  CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	         prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.7) THEN
	  CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &	         UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,
     &	         pmon,xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.10) THEN
	  CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	         prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE
	  WRITE(*,*)' Unknown lens array surface type. Aborting trace.'
	END IF


#if 0
	! debugging only
	If (.not. LROK) then
	  print*,'*** Ray missed surface ...'
	  print*,'ihat =',ihat
	  print*,'pin =',pin
	  print*,'vptL =',vptL
	  stop
	Else
	  print*,'=== Ray did not miss surface ...'
	End If
#endif


C  Return if ray missed surface

	IF (.NOT.LROK) RETURN

	Nsrfm=Nvec/Nmag

C  Compute refracted ray

	CALL DXMAT(Nx,Nhat)
	CALL DMPROD(NxNx,Nx,Nx,3,3,3)
	CALL DMPROD(cvec,NxNx,ihat,3,3,1)
	cfb=DSQRT(1d0+mu2*DDOTC(ihat,cvec))
	CALL DSMPROD(D1,cvec,(-mu),3)
	CALL DSMPROD(D2,Nhat,(-cfb),3)
	CALL DADD(rhat,D1,D2,3)

	IF (ifBuild) THEN

C    Lens partials

	  CALL DXMAT(Nx,Nhat)
	  CALL DMPROD(NxNx,Nx,Nx,3,3,3)

	  CALL DXMAT(D1,ihat)
	  CALL DMPROD(Nxix,Nx,D1,3,3,3)
	  CALL DMPROD(D3,ihat,Nxix,1,3,3)
	  S1=mu2/cfb
	  CALL DSMPROD(D2,Nhat,S1,3)
	  CALL DMPROD(D1,D2,D3,3,1,3)
	  D1(1)=D1(1)-cfb
	  D1(5)=D1(5)-cfb
	  D1(9)=D1(9)-cfb
	  CALL DMPROD(D2,Nx,ihat,3,3,1)
	  CALL DXMAT(D3,D2)
	  CALL DADD(D2,D3,Nxix,9)
	  CALL DSMPROD(D3,D2,mu,9)
	  CALL DADD(drdN,D1,D3,9)

	  CALL DMPROD(D1,ihat,NxNx,1,3,3)
	  CALL DSMPROD(D2,Nhat,(-S1),3)
	  CALL DMPROD(D3,D2,D1,3,1,3)
	  CALL DSMPROD(D1,NxNx,(-mu),9)
	  CALL DADD(drdii,D1,D3,9)

C    Partials of the normal

	  CALL DOUTER(Nouti,ihat,Nhat)
	  Ndoti=DDOTC(Nhat,ihat)
	  Niinv=1d0/Ndoti
	  CALL DSMPROD(NoidNdi,Nouti,Niinv,9)
	  CALL DEQUATE(dpdd,NoidNdi,9)
	  dpdd(1,1)=dpdd(1,1)-1d0
	  dpdd(2,2)=dpdd(2,2)-1d0
	  dpdd(3,3)=dpdd(3,3)-1d0
	  CALL DMPROD(dNdd,dNdp,dpdd,3,3,3)
	  S1=-L
	  CALL DSMPROD(dNdi,dNdd,S1,9)
	  CALL DNEGATE(dNda,dNdd,9)

C    Rotation-to-translation partials

	  CALL DSUB(dprot,prot,pr,3)
	  CALL DXMAT(dddth,dprot)

C    Ray direction partials

	  CALL DMPROD(drdd,drdN,dNdd,3,3,3)
	  CALL DMPROD(D1,drdN,dNdi,3,3,3)
	  CALL DADD(drdi,drdii,D1,9)
	  CALL DMPROD(D1,drdN,Nx,3,3,3)
	  CALL DMPROD(D2,drdd,dddth,3,3,3)
	  CALL DSUB(drdth,D2,D1,9)
	  CALL DMPROD(drda,drdN,dNda,3,3,3)

C    Beamwalk partials

	  CALL DPROJECT(Perpr,rhat)
	  CALL DMPROD(D1,Perpr,dpdd,3,3,3)
	  CALL DNEGATE(dgda,D1,9)
	  CALL DSMPROD(dgdi,dgda,L,9)
	  S1=1d0/Ndoti
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DMPROD(D2,Perpr,ihat,3,3,1)
	  CALL DMPROD(dgdd,D2,D1,3,1,3)
	  CALL DMPROD(dgdth,dgdd,dddth,3,3,3)

C    Pathlength partials

	  CALL DSUB(D3,rhat,ihat,3)
	  CALL DMPROD(dlda,D3,dpdd,1,3,3)
	  CALL DSMPROD(dldi,dlda,L,3)
	  S1=-DDOTC(D3,ihat)/Ndoti
	  CALL DSMPROD(dldd,Nhat,S1,3)
	  CALL DMPROD(dldth,dldd,dddth,1,3,3)

C    Combined partial matrices

	  DO 7 i=1,3
	    ip3=i+3
	    DO 6 j=1,3
	      jp3=j+3
	      dxidxim1(i,j)=drdi(i,j)
	      dxidxim1(i,jp3)=drda(i,j)
	      dxidxim1(ip3,j)=dgdi(i,j)
	      dxidxim1(ip3,jp3)=dgda(i,j)
	      dxidui(i,j)=drdth(i,j)
	      dxidui(i,jp3)=drdd(i,j)
	      dxidui(ip3,j)=dgdth(i,j)
	      dxidui(ip3,jp3)=dgdd(i,j)
    6	    CONTINUE
	    dxidxim1(i,7)=0e0
	    dxidxim1(ip3,7)=0e0
	    dxidxim1(7,i)=dldi(i)
	    dxidxim1(7,ip3)=dlda(i)
	    dxidui(7,i)=dldth(i)
	    dxidui(7,ip3)=dldd(i)
    7	  CONTINUE
	  dxidxim1(7,7)=1d0
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	        psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &		xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) THEN
	    Eout(1)=0d0
	    Eout(2)=0d0
	    Eout(3)=0d0
	  END IF
	ELSE
	  LRT=.TRUE.
	END IF

	DO 8 i=1,3
	  rout(i)=rhat(i)
	  pout(i)=pr(i)
   8	CONTINUE

	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE LensArray

C***********************************************************************

	SUBROUTINE ReturnSrf(ifBuild,ifPol,ifReturn,ifTransmit,ifLNsrf,
     &	 ifObs,LROK,LRT,Kc,Kr,
     &	 psi,pv,prot,pin,ihat,Evec,
     &	 na,nb,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui)

	USE Constants, ONLY: TWOPI
        USE traceutil_mod

C       Compute reflection point, normal, reflection and transition
C	matrices for a conic-of-revolution retroreflecting surface

	IMPLICIT NONE

	INTEGER i,j,ApType,nObs,ObsType(nObs)

	LOGICAL LROK,ifBuild,ifPol,ifReturn,ifTransmit,ifLNsrf,ifObs,LRT
!    &	        ChkRayTrans,LRT

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),S1,
     &    pr(3),rhat(3),L,pout(3),rout(3),na,nb,kxa,kxb,lambda,
     &	  xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs)

	REAL*8 prel(3),a,b,c,ksqrt,kp,km,mpr,rho(3),
     &         Kc,Kr,k2,itpsi,psitprel

	COMPLEX*16 Evec(3),Eout(3),C1,C2,NaCmplx

C	Compute mirror model
c		Kc=-e*e
c		Kr=-(1d0+e)*f

	!print*,' **pin(1:3) =', pin(1:3)
	!print*,' **pv(1:3) =', pv(1:3)

	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
   1	CONTINUE
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d10) THEN
	  mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
	  a=1d0+Kc*itpsi*itpsi
	  b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &	         ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
	  c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
	ELSE
	  a=0d0
	  b=itpsi
	  c=psitprel
	END IF

C	Compute intersection point and normal vector

	IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d10)) THEN
	  L=-c/b
	  IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98
	ELSE
	  k2=b*b-4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  IF (b.GE.0d0) THEN
	    ksqrt=-0.5d0*(b+DSQRT(k2))
	  ELSE
	    ksqrt=-0.5d0*(b-DSQRT(k2))
	  END IF
	  kp=ksqrt/a
	  km=c/ksqrt
	  IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	    IF (kp.GT.0D0) THEN
	      L=kp
	    ELSE
	      L=km
	    END IF
	  ELSE
	    IF (km.GT.0D0) THEN
	      L=km
	    ELSE
	      L=kp
	    END IF
	  END IF
	END IF

        !print*,' ** ReturnSrf: L =',L
        !print*,' ** ReturnSrf: ihat =',ihat

cLnegOK	-- Return surface will always have L positive;
c      *-- that would cause a confusion when settng XP which doesn't
c      *-- need a return surface. When the erroneous return surface is
c      *-- present, rays end up at a wrong 'XP' surface (a mirror image
c      *-- of the XP with respect to the focal plane) -jzlou, 02/28/2009
	IF (L.LT.0d0) GO TO 98

	IF (ifTransmit) THEN
	  DO 2 i=1,3
	    pr(i)=pin(i)+L*ihat(i)
	    rhat(i)=ihat(i)
   2	  CONTINUE
	ELSE
	  DO 11 i=1,3
	    pr(i)=pin(i)+L*ihat(i)
	    rhat(i)=-ihat(i)
  11	  CONTINUE
	END IF

C  Polarization calcs

	IF (ifPol) THEN
	  NaCmplx=DCMPLX(na,(-kxa))

C    Propagation through isotropic medium a to surface

c	  S1=MOD((TwoPi*L/lambda),TwoPi)
	  IF (ifReturn) THEN
	    S1=TWOPI*L/lambda
	  ELSE
	    S1=-TWOPI*L/lambda
	  END IF
	  C2=DCMPLX(0d0,S1)*NaCmplx
	  C1=CDEXP(C2)
	  Eout(1)=C1*Evec(1)
	  Eout(2)=C1*Evec(2)
	  Eout(3)=C1*Evec(3)
	END IF

	IF (ifBuild) THEN
	  DO 7 j=1,7
	    DO 6 i=1,7
	      dxidxim1(i,j)=0e0
    6	    CONTINUE
    7	  CONTINUE
	  DO 3 i=1,7
            dxidxim1(i,i)=1e0
    3	  CONTINUE
	  DO 9 j=1,6
	    DO 8 i=1,6
	      dxidui(i,j)=0e0
    8	    CONTINUE
    9	  CONTINUE
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	         psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) THEN
	    Eout(1)=0d0
	    Eout(2)=0d0
	    Eout(3)=0d0
	  END IF
	ELSE
	  LRT=.TRUE.
	END IF

	DO 10 i=1,3
	  rout(i)=rhat(i)
	  pout(i)=pr(i)
  10	CONTINUE
	RETURN

C	Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.

	NhatP(1:3)=psi(1:3)

	RETURN
	END SUBROUTINE ReturnSrf

C***********************************************************************

	SUBROUTINE FocalPln(ifBuild,ifPol,ifLNsrf,ifObs,
     &	LROK,LRT,
     &	psi,pv,prot,pin,ihat,Evec,
     &	na,nb,kxa,kxb,lambda,
     &	ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	pout,rout,L,Eout,s,dxidxim1,dxidui)

	USE Constants, ONLY: TWOPI
	use traceutil_mod

C       Compute reflection point, normal, reflection and transition
C	matrices for a general focal plane or reference surface

	IMPLICIT NONE

	INTEGER i,j,ip3,jp3,
     &	ApType,nObs,ObsType(nObs)

	LOGICAL LROK,ifBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &          LRT

	SREAL :: dxidxim1(7,7),dxidui(7,6),s(3)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),rhat(3),
     &  pr(3),L,Nhat(3),pout(3),rout(3),na,nb,kxa,kxb,lambda,
     &	xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs)

	REAL*8 prel(3),b,c,Ndoti,drdi(3,3),rho(3),
     &	dsda(3,3),dsdi(3,3),dsdd(3,3),dldd(3),
     &	dddth(3,3),dsdth(3,3),dldth(3),dlda(3),dldi(3),
     &	PN(3,3),dpvdth(3,3)

	REAL*8 D1(9),D2(9),S1,eye3(9)

	COMPLEX*16 Evec(3),Eout(3),C1,C2,NaCmplx

C	Compute mirror model

	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
	  Nhat(i)=psi(i)
   1	CONTINUE
	b=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	c=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)
	L=-c/b
cLnegOK -- focal plane always has L pos
	IF (L.LT.0d0) GO TO 98

	DO 2 i=1,3
	  S1=L*ihat(i)
	  pr(i)=pin(i)+S1
	  s(i)=prel(i)+S1
	  rhat(i)=ihat(i)
   2	CONTINUE

C  Polarization calcs

	IF (ifPol) THEN
	  NaCmplx=DCMPLX(na,(-kxa))

C    Propagation through isotropic medium a to surface

c	  S1=MOD((TwoPi*L/lambda),TwoPi)
	  S1=-TWOPI*L/lambda
	  C2=DCMPLX(0d0,S1)*NaCmplx
	  C1=CDEXP(C2)
	  Eout(1)=C1*Evec(1)
	  Eout(2)=C1*Evec(2)
	  Eout(3)=C1*Evec(3)
	END IF

C	Compute partials

	IF (ifBuild) THEN
	  CALL DZERO(eye3,9)
	  eye3(1)=1d0
	  eye3(5)=1d0
	  eye3(9)=1d0
	  CALL DPROJECT(PN,Nhat)
	  Ndoti=DDOTC(Nhat,ihat)

C	  Rotation-to-translation partials

	  CALL DSUB(D1,prot,pr,3)
	  CALL DXMAT(dddth,D1)
	  CALL DSUB(D1,pv,prot,3)
	  CALL DXMAT(dpvdth,D1)

C	  Ray direction partials

	  CALL DEQUATE(drdi,eye3,9)

C	  Beamwalk partials

	  S1=1d0/Ndoti
	  CALL DSMPROD(D1,ihat,S1,3)
	  CALL DOUTER(D2,D1,Nhat)
	  CALL DSUB(dsda,eye3,D2,9)
	  CALL DSMPROD(dsdi,dsda,L,9)
	  CALL DNEGATE(dsdd,dsda,9)
	  CALL DMPROD(D1,D2,dddth,3,3,3)
	  CALL DADD(dsdth,dpvdth,D1,9)

C	  Pathlength partials

	  CALL DSMPROD(dldd,Nhat,S1,3)
	  CALL DMPROD(dldth,dldd,dddth,1,3,3)
	  CALL DNEGATE(dlda,dldd,3)
	  CALL DSMPROD(dldi,dlda,L,3)

C	  Combined partial matrices

	  DO 4 i=1,3
	    ip3=i+3
	    DO 3 j=1,3
	      jp3=j+3
	      dxidxim1(i,j)=drdi(i,j)
	      dxidxim1(i,jp3)=0E0
	      dxidxim1(ip3,j)=dsdi(i,j)
	      dxidxim1(ip3,jp3)=dsda(i,j)
	      dxidui(i,j)=0E0
	      dxidui(i,jp3)=0E0
	      dxidui(ip3,j)=dsdth(i,j)
	      dxidui(ip3,jp3)=dsdd(i,j)
    3	    CONTINUE
	    dxidxim1(i,7)=0E0
	    dxidxim1(ip3,7)=0E0
	    dxidxim1(7,i)=dldi(i)
	    dxidxim1(7,ip3)=dlda(i)
	    dxidui(7,i)=dldth(i)
	    dxidui(7,ip3)=dldd(i)
    4	  CONTINUE
	  dxidxim1(7,7)=1d0
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	         psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &	         xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) THEN
	    Eout(1)=0d0
	    Eout(2)=0d0
	    Eout(3)=0d0
	  END IF
	ELSE
	  LRT=.TRUE.
	END IF

	DO 5 i=1,3
	  rout(i)=rhat(i)
	  pout(i)=pr(i)
   5	CONTINUE
        !print*,'** FocalPln: L =',L
	RETURN

C	Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.

	NhatP(1:3)=Nhat(1:3)
	RETURN
	END SUBROUTINE FocalPln

C***********************************************************************

	SUBROUTINE RefSrf(ifBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,
     &	 psi,pv,prot,pin,ihat,Evec,
     &	 na,nb,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui)
! Matched to calling routine in tracesub.F,
! pout = RayPos(1:3,iRay)
! rout = RayDir(1:3,iRay)
! L = RayL(iRay)

	USE Constants, ONLY: TWOPI
        use elt_mod, only : ifOPDModGrid,xMon,yMon,pMon,
     &                      GridSrfdx,nGridMat,GridMat,iEltToGridSrf
	use traceutil_mod

C       Compute reflection point, normal, reflection and transition
C	matrices for a general conic-of-revolution reference surface

	IMPLICIT NONE

	INTEGER i,j,ApType,nObs,ObsType(nObs),idx,jdx

	LOGICAL LROK,ifBuild,ifPol,ifLNsrf,ifObs,LRT
C    &     ChkRayTrans,LRT

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),pr(3),rhat(3),L,
     &   pout(3),rout(3),na,nb,kxa,kxb,lambda,rho(3),poff(3),
     &	 xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs),
     &   poffx,poffy,grid_dx,dL

	REAL*8 prel(3),a,b,c,ksqrt,kp,km,mpr,S1,
     &   Kc,Kr,k2,itpsi,psitprel

	COMPLEX*16 Evec(3),Eout(3),C1,C2,NaCmplx

C	Compute mirror model
c		Kc=-e*e
c		Kr=-(1d0+e)*f

	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
   1	CONTINUE
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d10) THEN
	  mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
	  a=1d0+Kc*itpsi*itpsi
	  b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &	    ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
	  c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
	ELSE
	  a=0d0
	  b=itpsi
	  c=psitprel
	END IF

C	Compute intersection point and normal vector

	IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d10)) THEN
	  L=-c/b
	  IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98
	ELSE
	  k2=b*b-4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  IF (b.GE.0d0) THEN
	    ksqrt=-0.5d0*(b+DSQRT(k2))
	  ELSE
	    ksqrt=-0.5d0*(b-DSQRT(k2))
	  END IF
	  kp=ksqrt/a
	  km=c/ksqrt
	  IF(ifLNsrf) THEN
	    IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	      L=kp
	    ELSE
	      L=km
	    END IF
	  ELSE
	    IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	      IF (kp.GT.0D0) THEN
	        L=kp
	      ELSE
	        L=km
	      END IF
	    ELSE
	      IF (km.GT.0D0) THEN
	        L=km
	      ELSE
	        L=kp
	      END IF
	    END IF
	    IF (L.LT.0d0) GO TO 98
	  END IF
	END IF

	DO 2 i=1,3
	  pr(i)=pin(i)+L*ihat(i)
	  rhat(i)=ihat(i)
   2	CONTINUE


C  Polarization calcs

	IF (ifPol) THEN
	  NaCmplx=DCMPLX(na,(-kxa))

C    Propagation through isotropic medium a to surface

c	  S1=MOD((TwoPi*L/lambda),TwoPi)
	  S1=-TWOPI*L/lambda
	  C2=DCMPLX(0d0,S1)*NaCmplx
	  C1=CDEXP(C2)
	  Eout(1)=C1*Evec(1)
	  Eout(2)=C1*Evec(2)
	  Eout(3)=C1*Evec(3)
	END IF

	IF (ifBuild) THEN

C	  Combined partial matrices

	  DO 7 j=1,7
	    DO 6 i=1,7
	      dxidxim1(i,j)=0e0
    6	    CONTINUE
	    dxidxim1(j,j)=1e0
    7	  CONTINUE
	  DO 9 j=1,6
	    DO 8 i=1,7
	      dxidui(i,j)=0e0
    8	    CONTINUE
    9	  CONTINUE
	END IF

	! This if (true/false) stmt was added by jzlou for debugging
	if (.true.) then

C  Check if ray is obscured
	IF (ifObs) THEN
	  DO i=1,3
	    rho(i)=pin(i)-pv(i)+L*ihat(i)
	  END DO
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	         psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) THEN
	    Eout(1)=0d0
	    Eout(2)=0d0
	    Eout(3)=0d0
	  END IF
	ELSE
	  LRT=.TRUE.
	END IF

	end if

	DO 10 i=1,3
	  rout(i)=rhat(i)
	  pout(i)=pr(i)
  10	CONTINUE

#if 1
	! If an OPD-modulation grid is defined at this reference surface,
        ! modify L (ray path length) accordingly. -jzlou
        If (ifOPDModGrid(crnt_Elt)) Then
	  poff=pout(1:3)-pMon(1:3,crnt_Elt)
	  poffx = DDOTC(poff,xMon(1:3,crnt_Elt))
	  poffy = DDOTC(poff,yMon(1:3,crnt_Elt))
          grid_dx = GridSrfdx(crnt_Elt)
	  jdx = poffx/grid_dx; idx = poffy/grid_dx
          jdx = jdx + (nGridMat(crnt_Elt)+1)/2
          idx = idx + (nGridMat(crnt_Elt)+1)/2
          if (idx<1 .OR. idx>nGridMat(crnt_Elt) .OR.
     &        jdx<1 .OR. jdx>nGridMat(crnt_Elt)) then
            GO TO 99  ! ray outside grid region
	  end if
	  dL=GridMat(idx,jdx,iEltToGridSrf(crnt_Elt))
	  L=L+dL
        End If
#endif
  99	RETURN

C	Here if ray misses surface

  98	CONTINUE
	LROK=.FALSE.

	NhatP(1:3)=psi(1:3) ! Ref surf is always flat -jzlou

	RETURN
	END SUBROUTINE RefSrf

C***********************************************************************

	SUBROUTINE ObsSrf(ifBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,
     &	 psi,pv,prot,pin,ihat,Evec,
     &	 na,nb,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui)

	USE Constants, ONLY: TWOPI
	use traceutil_mod

C  Compute incidence point for a ray on a conic-of-revolution
C  surface and determine if the ray is obscured.

	IMPLICIT NONE

	INTEGER i,j,ApType,nObs,ObsType(nObs)

	LOGICAL LRT,LROK,ifBuild,ifPol,ifLNsrf,!ChkRayTrans,
     &          ifObs

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &   pr(3),rhat(3),L,pout(3),rout(3),
     &	 na,nb,kxa,kxb,lambda

	REAL*8 prel(3),a,b,c,ksqrt,kp,km,mpr,
     &   Kc,Kr,rho(3),itpsi,psitprel,S1,k2,
     &	 xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs)

	COMPLEX*16 Evec(3),Eout(3),C1,C2,NaCmplx

	DO 1 i=1,3
	  prel(i)=pin(i)-pv(i)
   1	CONTINUE
	itpsi=ihat(1)*psi(1)+ihat(2)*psi(2)+ihat(3)*psi(3)
	psitprel=psi(1)*prel(1)+psi(2)*prel(2)+psi(3)*prel(3)

	IF (ABS(Kr).LT.1d10) THEN
	  mpr=prel(1)*prel(1)+prel(2)*prel(2)+prel(3)*prel(3)
	  a=1d0+Kc*itpsi*itpsi
	  b=2d0*(Kr*itpsi+Kc*itpsi*psitprel+
     &	  ihat(1)*prel(1)+ihat(2)*prel(2)+ihat(3)*prel(3))
	  c=(2d0*Kr+Kc*psitprel)*psitprel+mpr
	ELSE
	  a=0d0
	  b=itpsi
	  c=psitprel
	END IF

C	Compute intersection point and normal vector

	IF ((ABS(a).LE.1d-10).OR.(ABS(Kr).GE.1d10)) THEN
	  L=-c/b
	  IF (.NOT.ifLNsrf.AND.(L.LT.0d0)) GO TO 98
	ELSE
	  k2=b*b-4d0*a*c
	  IF (k2.LT.0d0) GO TO 98
	  IF (b.GE.0d0) THEN
	    ksqrt=-0.5d0*(b+DSQRT(k2))
	  ELSE
	    ksqrt=-0.5d0*(b-DSQRT(k2))
	  END IF
	  kp=ksqrt/a
	  km=c/ksqrt
	  IF(ifLNsrf) THEN
	    IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	      L=kp
	    ELSE
	      L=km
	    END IF
	  ELSE
	    IF (ABS((kp*kp)-mpr).LE.ABS((km*km)-mpr)) THEN
	      IF (kp.GT.0D0) THEN
	        L=kp
	      ELSE
	        L=km
	      END IF
	    ELSE
	      IF (km.GT.0D0) THEN
	        L=km
	      ELSE
	        L=kp
	      END IF
	    END IF
	    IF (L.LT.0d0) GO TO 98
	  END IF
	END IF

C  Compute point of incidence and output ray
	DO 2 i=1,3
	  S1=L*ihat(i)
	  pr(i)=pin(i)+S1
	  rho(i)=prel(i)+S1
	  rhat(i)=ihat(i)
   2	CONTINUE

C  Polarization calcs
	IF (ifPol) THEN
	  NaCmplx=DCMPLX(na,(-kxa))

C  Propagation through isotropic medium a to surface
	  S1=-TWOPI*L/lambda
	  C2=DCMPLX(0d0,S1)*NaCmplx
	  C1=CDEXP(C2)
	  Eout(1)=C1*Evec(1)
	  Eout(2)=C1*Evec(2)
	  Eout(3)=C1*Evec(3)
	END IF

C  Combined partial matrices
	IF (ifBuild) THEN
	  DO 7 j=1,7
	    DO 6 i=1,7
	      dxidxim1(i,j)=0e0
    6	    CONTINUE
	    dxidxim1(j,j)=1e0
    7	  CONTINUE
	  DO 9 j=1,6
	    DO 8 i=1,6
	      dxidui(i,j)=0e0
    8	    CONTINUE
    9	  CONTINUE
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	        psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &	        xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) Eout(1:3)=0d0
	ELSE
	  LRT=.TRUE.
	END IF

	DO 11 i=1,3
	  rout(i)=rhat(i)
	  pout(i)=pr(i)
  11	CONTINUE
	RETURN

C  Here if ray becomes undefined
  98	CONTINUE
	LROK=.FALSE.
	LRT=.FALSE.
	RETURN
	END SUBROUTINE ObsSrf


C***********************************************************************

	FUNCTION ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &		    psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &              xGrid,yGrid,zGrid,taperFac)
        use Constants, only: PI
        use elt_mod, only : ApHexBound,PolyApEqn,VptElt,
     &                      PolyObsEqn,ifAmplSrf,nAmplMat,AmplMat,
     &                      AmplSrfdx,iEltToAmplSrf,xMon,yMon,pMon
        use cfiles_mod, only : stepApod,linApod,gaussApod,iApodType
	use traceutil_mod

! Check wether a given ray passes an aperture or hits an obscuration
! -jzlou

! A 'taperFac' is defined for this ray on this element if the ray falls
! into the taperred region of a tapering aperture or obscuration. -jzlou

	IMPLICIT NONE

	!LOGICAL ChkRayTrans,ChkLineSide
	LOGICAL ChkRayTrans
	INTEGER ApType,nObs,ObsType(nObs),iObs
	LOGICAL ifBuild,ifPol,ifHitsPosObs,ifHitsNegObs,ifLNsrf

	REAL*8 psi(3),rho(3),xGrid(3),yGrid(3),zGrid(3),ApVec(6),
     &	  ObsVec(6,nObs),r2,A2,A3,px,py,
     &	  ax,ay,bx,by,cx,cy,rc,xcmin,xcmax,ycmin,ycmax,yt,
     &	  sqrt_A1, sqrt_A2, sqrt_r2, rd1,rd2,r1,taperFac,
     &    px2,py2,coss,sins,ApRot,ObsRot,tmpx,tmpy,
     &	  inPt(2),inPt_val,ray_val,dd,tf,sig,pOff(3),pOffx,pOffy,
     &    xsb(3),ysb(3),th1,th2,dth
	Real*8, save :: Asb(3,mSB),Bsb(3,mSB),Csb(3,mSB),
     &                  xcm(mSB),ycm(mSB)
	Integer :: nSides,iside,iamp,jamp,jAmpSrf,nBlades,
     &             ib,j,i1,i2
	Logical :: inHex,inPoly,inSB
	Logical, save :: firstEntry_sunBlades=.true.


	ifHitsPosObs=.FALSE.
	ifHitsNegObs=.FALSE.

C  Compute aperture and obscurations

	! rho(3), an argument, was computed as
        ! rho(i)=pin(i)-pv(i)+L*ihat(i)
        ! in the calling routine, where pin+L*ihat gives the intersection
        ! point on the element surface (in global coords) and pv(i) is
        ! element vertex point. So rho is a vector from element vertex to
        ! the ray intersection point on the surface. -jzlou

	! Project ray position vector onto x and y axes in aperture
	! coord frame -jzlou
	px=DDOTC(xGrid,rho); py=DDOTC(yGrid,rho)

	IF (ApType.EQ.1) THEN  ! circular aperture -jzlou
	  r2=(ApVec(2)-px)**2+(ApVec(3)-py)**2
	  A2=ApVec(1)**2
	  IF (r2.GT.A2) GO TO 99
	ELSE IF (ApType.EQ.2) THEN ! rectangular aperture -jzlou
	  IF ((px.LT.ApVec(1)).OR.(px.GT.ApVec(2)).OR.(py.LT.ApVec(3))
     &	      .OR.(py.GT.ApVec(4))) GO TO 99
	ELSE IF (ApType.EQ.3) THEN ! tapered circulcar aperture -jzlou
	  r2=(ApVec(2)-px)**2+(ApVec(3)-py)**2
	  A2=ApVec(1)**2
	  IF (r2.GT.A2) then
	    GO TO 99
	  ELSE
            ! As ray passes thru aperture, taper its amplitude
            ! linearly if ray falls into tapered region
	    sqrt_A1=ApVec(1)*(1.0-ApVec(4)); sqrt_A2=ApVec(1)
	    sqrt_r2=dsqrt(r2)
	    if (sqrt_r2>sqrt_A1) then ! fall into tapered region
	      if ((sqrt_A2-sqrt_A1)<1d-20) then
	        print*,'ChkRayTrans warning: tapered size = 0!'
	        GO TO 99
	      end if
	      ! tapering factor for this ray on this element only
	      taperFac = (sqrt_A2-sqrt_r2)/(sqrt_A2-sqrt_A1)
	      !print*,'ChkRay(): taperFac = ...'
            else
	      taperFac = 1.0 ! ray is outside tapered region, but inside aperture
	    end if
	  END IF
        ELSE IF (ApType.EQ.4) THEN ! elliptical aperture -jzlou
#if 0
          Test rotation of aperture
          px2=px*0.866+py*0.5; py2=-px*0.5+py*0.866
          Cx=ApVec(3)*0.866+ApVec(4)*0.5
          Cy=-ApVec(3)*0.5+ApVec(4)*0.866
 	  r2=((Cx-px2)/ApVec(1))**2
     &       +((Cy-py2)/ApVec(2))**2
#endif
          r2=((ApVec(3)-px)/ApVec(1))**2
     &       +((ApVec(4)-py)/ApVec(2))**2
          if (r2.GT.1.d0) GO TO 99
        ELSE IF (ApType.EQ.5) THEN ! tapered elliptical aperture -jzlou
          taperFac = 1.0  ! default
          rd1 = 1.d0-ApVec(5)  ! percent of non-tapered
          rd2 = rd1**2
          r2=((px-ApVec(3))/ApVec(1))**2
     &         +((py-ApVec(4))/ApVec(2))**2
          if (r2 > 1.d0) then
            GO TO 99 ! Ray lies outside aperture
          else if (r2 > rd2) then
            ! Ray inside aperture, and fall in tapered region
            r1=sqrt(r2)
            taperFac = (1.d0-r1)/(1.0-rd1)
          end if

	ELSE IF (ApType.EQ.6) THEN ! Hexagonal aperture

	   inHex=.false.

	   ! Get aperture rotation specified in degree, convert
	   ! it to radians
	   !ApRot = ApVec(4)*pi/180.d0  ! input in degree

	   ApRot = ApVec(4)   ! input in radian
	   coss=cos(ApRot); sins=sin(ApRot)

!	   px2=px*coss+py*sins; py2=-px*sins+py*coss
!          Cx=ApVec(2)*coss+ApVec(3)*sins
!          Cy=-ApVec(2)*sins+ApVec(3)*coss
!	   px=px2-Cx; py=py2-Cy

	   px2=px*coss+py*sins; py2=-px*sins+py*coss
    	   px2=px2-ApVec(2); py2=py2-ApVec(3)

	   if (.true.) then
	     ! sharp top Hex
	     if (px2 < ApHexBound(1,1,crnt_Elt).and.
     &     py2<(ApHexBound(1,2,crnt_Elt)*(px2-ApHexBound(2,2,crnt_Elt))
     &     +ApHexBound(3,2,crnt_Elt)).and.
     &     py2<(ApHexBound(1,3,crnt_Elt)*(px2-ApHexBound(2,3,crnt_Elt))
     &     +ApHexBound(3,3,crnt_Elt)).and.
     &     px2>ApHexBound(1,4,crnt_Elt) .and.
     &     py2>(ApHexBound(1,5,crnt_Elt)*(px2-ApHexBound(2,5,crnt_Elt))
     &     +ApHexBound(3,5,crnt_Elt)).and.
     &     py2>(ApHexBound(1,6,crnt_Elt)*(px2-ApHexBound(2,6,crnt_Elt))
     &         +ApHexBound(3,6,crnt_Elt))) then
               inHex = .true.
             end if
	   else
	     ! flat top Hex
             if (py2<(ApHexBound(1,1,crnt_Elt)
     *           *(px2-ApHexBound(2,1,crnt_Elt))
     &         +ApHexBound(3,1,crnt_Elt)).and.
     &         py2<ApHexBound(1,2,crnt_Elt).and.
     &     py2<(ApHexBound(1,3,crnt_Elt)*(px2-ApHexBound(2,3,crnt_Elt))
     &     +ApHexBound(3,3,crnt_Elt)).and.
     &     py2>(ApHexBound(1,4,crnt_Elt)*(px2-ApHexBound(2,4,crnt_Elt))
     &     +ApHexBound(3,4,crnt_Elt)).and.
     &     py2>ApHexBound(1,5,crnt_Elt).and.
     &     py2>(ApHexBound(1,6,crnt_Elt)*(px2-ApHexBound(2,6,crnt_Elt))
     &         +ApHexBound(3,6,crnt_Elt))) then
               inHex = .true.
             end if
	   end if

	   if (.not.inHex) GO TO 99
	ELSE IF (ApType.EQ.7) THEN ! Polygonal aperture
	   inPoly=.true.
	   nSides=ApVec(3)

	   !ApVec(1)=ApVec(1)
#if 0
	     print*,'PolyVtx(1:2,1)=',PolyVtx(1:2,1)
	     print*,'PolyVtx(1:2,2)=',PolyVtx(1:2,2)
	     print*,'PolyVtx(1:2,3)=',PolyVtx(1:2,3)
	     print*,'nSides= ', nSides
	     print*,'PolyEqn(1:3,1)=',PolyEqn(1:3,1)
	     print*,'PolyEqn(1:3,2)=',PolyEqn(1:3,2)
	     print*,'PolyEqn(1:3,3)=',PolyEqn(1:3,3)
#endif

	   ! ApVec(1) and ApVec(2) are element center coord
	   px2=px-ApVec(1); py2=py-ApVec(2)

	   !inPt(1:2)=ApVec(1:2)
	   inPt(1:2)=0.d0  ! ref. pt in elem-local coords

	   Do iside=1,nSides
             inPt_val=inPt(1)*PolyApEqn(1,iside,crnt_Elt)
     &                +inPt(2)*PolyApEqn(2,iside,crnt_Elt)
     &                +PolyApEqn(3,iside,crnt_Elt)
             ray_val=px2*PolyApEqn(1,iside,crnt_Elt)
     &               +py2*PolyApEqn(2,iside,crnt_Elt)
     &               +PolyApEqn(3,iside,crnt_Elt)
             if ((inPt_val>0.d0 .and. ray_val<0.d0) .or.
     &           (inPt_val<0.d0 .and. ray_val>0.d0)) then
               inPoly=.false.
               exit
             end if
           End Do
	   if (.not.inPoly) GO TO 99

	ELSE IF (ApType.EQ.8) THEN ! Tapered Polygonal aperture
           inPoly=.true.
           nSides=ApVec(3)

           ! ApVec(1) and ApVec(2) are element center coord
           px2=px-ApVec(1); py2=py-ApVec(2)

           inPt(1:2)=0.d0  ! ref. pt in elem-local coords

           Do iside=1,nSides
             inPt_val=inPt(1)*PolyApEqn(1,iside,crnt_Elt)
     &                +inPt(2)*PolyApEqn(2,iside,crnt_Elt)
     &                +PolyApEqn(3,iside,crnt_Elt)
             ray_val=px2*PolyApEqn(1,iside,crnt_Elt)
     &               +py2*PolyApEqn(2,iside,crnt_Elt)
     &               +PolyApEqn(3,iside,crnt_Elt)
             if ((inPt_val>0.d0 .and. ray_val<0.d0) .or.
     &           (inPt_val<0.d0 .and. ray_val>0.d0)) then
               inPoly=.false.
               exit
             end if
           End Do
           if (.not.inPoly) GO TO 99

	   taperFac = 1.0  ! default
           Do iside=1,nSides
              ray_val=px2*PolyApEqn(1,iside,crnt_Elt)
     &                +py2*PolyApEqn(2,iside,crnt_Elt)
     &                +PolyApEqn(3,iside,crnt_Elt)
            rd2=PolyApEqn(1,iside,crnt_Elt)*PolyApEqn(1,iside,crnt_Elt)
     &         +PolyApEqn(2,iside,crnt_Elt)*PolyApEqn(2,iside,crnt_Elt)
              dd=abs(ray_val)/sqrt(rd2)
              If (dd .LE.ApVec(4)) Then
                ! ray falls in tapered region of this side of polygon,
                ! so define appropriate tapering.
                if (iApodType==stepApod) then
                  tf=ApVec(5)
                else if (iApodType==linApod) then
                  tf=dd/ApVec(4)  ! linear tapering
                else if (iApodType==gaussApod) then
                  ! Gaussian tapering
                  sig=ApVec(5)/1d0; dd=ApVec(4)-dd;
                  tf=exp(-1d0*dd*dd/(sig*sig))
                end if
                if (tf<taperFac) taperFac=tf
              End If
           End Do
	   !print*,'*** taperFac =',taperFac

        ELSE IF (ApType.EQ.9) THEN
C  ApType= SunBlades
C  ApVec=  xc yc R1 R2 N  % R1=innerRadius,R2=outerRadius,N=numberOfBlades
          px2=px-ApVec(1); py2=py-ApVec(2)
          r1=ApVec(3); r2=ApVec(4)
          nBlades=ApVec(5)

          dth=2d0*pi/nBlades
          th1=0d0; th2=th1+dth

	  If (firstEntry_sunBlades) Then
            ! Precompute Asb(3,nBlades),Bsb(3,nBlades),Csb(3,nBlades)
            Do ib=1,nBlades
              xsb(1)=r1*cos(th1); ysb(1)=r1*sin(th1)
              xsb(2)=r1*cos(th2); ysb(2)=r1*sin(th2)
              xsb(3)=r2*cos((th1+th2)/2d0)
              ysb(3)=r2*sin((th1+th2)/2d0)
              Do j=1,3
                if (j==1) then
                  i1=1; i2=3
                else if (j==2) then
                  i1=3; i2=2
                else
                  i1=2; i2=1
                end if
                if (abs(xsb(i1)-xsb(i2))<1d-20) then
                  ! horizontal line
                  Asb(j,ib)=1;Bsb(j,ib)=0d0;Csb(j,ib)=-xsb(i1)
                else
                  Asb(j,ib)=(ysb(i2)-ysb(i1))/(xsb(i2)-xsb(i1))
                  Bsb(j,ib)=-1d0
                  Csb(j,ib)=ysb(i1)-xsb(i1)*Asb(j,ib)
                end if
              End Do
              ! xcm,ycm store center of mass of triangle
              xcm(ib)=(xsb(1)+xsb(2)+xsb(3))/3d0
              ycm(ib)=(ysb(1)+ysb(2)+ysb(3))/3d0
	      th1=th2; th2=th1+dth
            End Do
            firstEntry_sunBlades=.false.
          End If

	  inSB=.false.
          Do ib=1,nBlades
            r2=px2**2+py2**2
	    If (r2.LE.r1**2) Then
              inSB=.true.; exit
            Else If (inTrig(Asb(1:3,ib),Bsb(1:3,ib),Csb(1:3,ib),
     &               xcm(ib),ycm(ib),px2,py2)) Then
              inSB=.true.; exit
            End If
          End Do
	  if (.not.inSB) GO TO 99

	ELSE IF (ApType.EQ.-4) THEN
	  ! This section is disabled for now. 07/04, -jzlou
	  ! still rectangular aperture, but with a interior margin
          ! specified by 'rc'.
	  ! This secion is not supported in macosio.F and smacosio.F,
	  ! so it's never used currently? -jzlou
	  IF ((px.LT.ApVec(1)).OR.(px.GT.ApVec(2)).OR.(py.LT.ApVec(3))
     &	      .OR.(py.GT.ApVec(4))) GO TO 99
	  rc=ApVec(5)
	  xcmin=ApVec(1)+rc; xcmax=ApVec(2)-rc
	  ycmin=ApVec(3)+rc; ycmax=ApVec(4)-rc
	  IF ((px.LT.xcmin).AND.(py.LT.ycmin)) THEN
	    yt=ycmin-SQRT(rc*rc-(px-xcmin)**2)
	    IF (py.LT.yt) GO TO 99
	  ELSE IF ((px.GT.xcmax).AND.(py.LT.ycmin)) THEN
	    yt=ycmin-SQRT(rc*rc-(px-xcmax)**2)
	    IF (py.LT.yt) GO TO 99
	  ELSE IF ((px.LT.xcmin).AND.(py.GT.ycmax)) THEN
	    yt=ycmax+SQRT(rc*rc-(px-xcmin)**2)
	    IF (py.GT.yt) GO TO 99
	  ELSE IF ((px.GT.xcmax).AND.(py.GT.ycmax)) THEN
	    yt=ycmax+SQRT(rc*rc-(px-xcmax)**2)
	    IF (py.GT.yt) GO TO 99
	  END IF
	END IF

C  Check if ray hits an obscuration
	DO 10 iObs=1,nObs

#if 0
	  ! debugging only
	   print*,'In Obs: nObs,ObsType(iObs) =',nObs,ObsType(iObs)
	   print*,'In Obs: xgrid =',xgrid
	   print*,'In Obs: ygrid =',ygrid
	   print*,'ObsVec(1,iObs) =',ObsVec(1,iObs)
	   print*,'ObsVec(2,iObs) =',ObsVec(2,iObs)
	   print*,'ObsVec(3,iObs) =',ObsVec(3,iObs)
	   print*,'px,py =',px,py
	   if (nC==100) stop
	   nC=nC+1
#endif

	  IF (ABS(ObsType(iObs)).EQ.1) THEN
	    ! Circular
	    r2=(ObsVec(2,iObs)-px)**2+(ObsVec(3,iObs)-py)**2
	    A2=ObsVec(1,iObs)**2
	    IF (r2.LT.A2) THEN
	      IF (ObsType(iObs).GT.0) THEN
	        ifHitsPosObs=.TRUE.
	      ELSE
	        ifHitsNegObs=.TRUE.
	      END IF
	    END IF
	  ELSE IF (ABS(ObsType(iObs)).EQ.2) THEN
	    IF ((px.GT.ObsVec(1,iObs)).AND.(px.LT.ObsVec(2,iObs))
     &	    .AND.(py.GT.ObsVec(3,iObs)).AND.(py.LT.ObsVec(4,iObs))) THEN
	      IF (ObsType(iObs).GT.0) THEN
	        ifHitsPosObs=.TRUE.
	      ELSE
	        ifHitsNegObs=.TRUE.
	      END IF
 	    END IF
	  ELSE IF (ABS(ObsType(iObs)).EQ.3) THEN
	    r2=(ObsVec(2,iObs)-px)**2+(ObsVec(3,iObs)-py)**2
	    A2=ObsVec(1,iObs)**2
	    A3=ObsVec(4,iObs)**2
	    IF ((r2.LT.A2).AND.(r2.GT.A3)) THEN
	      IF (ObsType(iObs).GT.0) THEN
	        ifHitsPosObs=.TRUE.
	      ELSE
	        ifHitsNegObs=.TRUE.
	      END IF
 	    END IF
	  ELSE IF (ABS(ObsType(iObs)).EQ.4) THEN
	    ax=ObsVec(1,iObs)
	    ay=ObsVec(2,iObs)
	    bx=ObsVec(3,iObs)
	    by=ObsVec(4,iObs)
	    cx=ObsVec(5,iObs)
	    cy=ObsVec(6,iObs)
	    IF (.NOT.ChkLineSide(px,py,ax,ay,bx,by,cx,cy)) THEN
	      GO TO 10
	    ELSE
	      IF (.NOT.ChkLineSide(px,py,cx,cy,ax,ay,bx,by)) THEN
	        GO TO 10
	      ELSE
	        IF (.NOT.ChkLineSide(px,py,bx,by,cx,cy,ax,ay)) THEN
	          GO TO 10
	        ELSE
	          IF (ObsType(iObs).GT.0) THEN
	            ifHitsPosObs=.TRUE.
	          ELSE
	            ifHitsNegObs=.TRUE.
	          END IF
	        END IF
	      END IF
	    END IF
	  ELSE IF (ABS(ObsType(iObs)).EQ.5) THEN
	    r2=((ObsVec(3,iObs)-px)/ObsVec(1,iObs))**2
     &      +((ObsVec(4,iObs)-py)/ObsVec(2,iObs))**2
	    IF (r2.LT.1.d0) THEN
	      IF (ObsType(iObs).GT.0) THEN
	        ifHitsPosObs=.TRUE.
	      ELSE
	        ifHitsNegObs=.TRUE.
	      END IF
	    END IF
	  ELSE IF (ObsType(iObs).EQ.7) THEN
	    ! tapered elliptical obscuration
	    taperFac = 1.0  ! default
	    rd1 = 1-ObsVec(5,iObs)  ! percent of non-tapered
	    rd1 = rd1**2
	    r2=((px-ObsVec(3,iObs))/ObsVec(1,iObs))**2
     &         +((py-ObsVec(4,iObs))/ObsVec(2,iObs))**2
	    if (r2 < rd1) then
	      ! Ray hit obs and beyond tapered region,
	      ! turn on Obs flag
              ifHitsPosObs=.TRUE.
	    else if (r2 < 1.0 .and. r2 > rd1) then
	      ! Ray hit obs, and fell into tapered region.
	      ! In this case, set ifHitsPosObs to false (default)
	      ! tapering factor for this ray on this element only
	      r1=sqrt(r2)
              taperFac = (r1-rd1)/(1.0-rd1)
	    end if
	  ELSE IF (ABS(ObsType(iObs)).EQ.8) THEN
	    ! polygon obscuration
	    inPoly=.true.
            nSides=ObsVec(3,iObs)

	    ! ObsVec(1,iObs) and ObsVec(2,iObs) are element center coord
            px2=px-ObsVec(1,iObs); py2=py-ObsVec(2,iObs)

           inPt(1:2)=0.d0  ! ref. pt in elem-local coords

	   if (.false. .and. crnt_Elt==18) then
             do iside=1,nSides
               print*,'PolyObsEqn(1:3,iside,iObs,crnt_Elt) =',
     &                PolyObsEqn(1:3,iside,iObs,crnt_Elt)
             end do
	     stop
	   end if

	   Do iside=1,nSides
             inPt_val=inPt(1)*PolyObsEqn(1,iside,iObs,crnt_Elt)
     &            +inPt(2)*PolyObsEqn(2,iside,iObs,crnt_Elt)
     &            +PolyObsEqn(3,iside,iObs,crnt_Elt)
             ray_val=px2*PolyObsEqn(1,iside,iObs,crnt_Elt)
     &            +py2*PolyObsEqn(2,iside,iObs,crnt_Elt)
     &            +PolyObsEqn(3,iside,iObs,crnt_Elt)
             if ((inPt_val>0.d0 .and. ray_val<0.d0) .or.
     &           (inPt_val<0.d0 .and. ray_val>0.d0)) then
               inPoly=.false.
               exit
             end if
           End Do
           IF (ObsType(iObs).GT.0 .and. inPoly) THEN
	     ifHitsPosObs=.TRUE.
           ELSE
             ifHitsNegObs=.TRUE.
           END IF

	  ELSE IF (ABS(ObsType(iObs)).EQ.9) THEN
	    ! RotRectangle or RotNegRectangle
	    ! ObsVec(1:2,iObs) = W,H, ObsVec(3:4,iObs) = dx,dy
	    ! ObsVec(5,iObs) = rot_angle

	    ObsRot=ObsVec(5,iObs)   ! input in radian
            coss=cos(ObsRot); sins=sin(ObsRot)

            tmpx=px-ObsVec(3,iObs); tmpy=py-ObsVec(4,iObs)
	    px2=tmpx*coss+tmpy*sins; py2=-tmpx*sins+tmpy*coss

	    IF ((px2.GT.(-0.5*ObsVec(1,iObs))).AND.
     &          (px2.LT.(0.5*ObsVec(1,iObs)))
     &          .AND.(py2.GT.(-0.5*ObsVec(2,iObs))).AND.
     &          (py2.LT.(0.5*ObsVec(2,iObs)))) THEN
              IF (ObsType(iObs).GT.0) THEN
                ifHitsPosObs=.TRUE.
              ELSE
                ifHitsNegObs=.TRUE.
              END IF
            END IF
	  ELSE IF (ABS(ObsType(iObs)).EQ.10) THEN
	    ! RotEllipse or RotNegEllipse
	    ! ObsVec(1:2,iObs) = a,b (half axises), ObsVec(3:4,iObs) = dx,dy
	    ! ObsVec(5,iObs) = rot_angle

	    ObsRot=ObsVec(5,iObs)   ! input in radian
            coss=cos(ObsRot); sins=sin(ObsRot)
	    tmpx=px-ObsVec(3,iObs); tmpy=py-ObsVec(4,iObs)
            px2=tmpx*coss+tmpy*sins; py2=-tmpx*sins+tmpy*coss

            r2=(px2/ObsVec(1,iObs))**2+(py2/ObsVec(2,iObs))**2
            IF (r2.LT.1.d0) THEN
              IF (ObsType(iObs).GT.0) THEN
                ifHitsPosObs=.TRUE.
              ELSE
                ifHitsNegObs=.TRUE.
              END IF
            END IF
	  END IF
  10	CONTINUE

C  Check if ray is obscured
	IF (ifHitsPosObs .AND. .NOT.ifHitsNegObs) THEN
	  ChkRayTrans=.FALSE.
	ELSE
	  ChkRayTrans=.TRUE.
	END IF

#if 1
        ! Ray amplitude adjustment if needed, currently use a table
        ! lookup (nearest neighbor) scheme. A spline interpolation
        ! could be added.
        If (ChkRayTrans .AND. ifAmplSrf(crnt_Elt)) Then
          px=DDOTC(xMon(1:3,crnt_Elt),rho)
          py=DDOTC(yMon(1:3,crnt_Elt),rho)
          pOff(1:3)=pMon(1:3,crnt_Elt)-VptElt(1:3,crnt_Elt)
          pOffx=DDOTC(xMon(1:3,crnt_Elt),pOff(1:3))
          pOffy=DDOTC(yMon(1:3,crnt_Elt),pOff(1:3))
          px2=px-pOffx; py2=py-pOffy

          iamp=px2/AmplSrfdx(crnt_Elt)+nAmplMat(crnt_Elt)/2d0
          if (iamp<0) then
            iamp=1
          else if (iamp>nAmplMat(crnt_Elt)) then
            iamp=nAmplMat(crnt_Elt)
          end if
          jamp=py2/AmplSrfdx(crnt_Elt)+nAmplMat(crnt_Elt)/2d0
          if (jamp<0) then
            jamp=1
          else if (jamp>nAmplMat(crnt_Elt)) then
            jamp=nAmplMat(crnt_Elt)
          end if
          jAmpSrf=iEltToAmplSrf(crnt_Elt)
          taperFac=AmplMat(iamp,jamp,jAmpSrf)
        End If
#endif

	RETURN

C  Here if ray misses aperture
  99	CONTINUE
	ChkRayTrans=.FALSE.
	RETURN
	END FUNCTION ChkRayTrans

C***********************************************************************

	FUNCTION ChkLineSide(px,py,ax,ay,bx,by,cx,cy)

	IMPLICIT NONE

	LOGICAL ChkLineSide
	REAL*8 px,py,ax,ay,bx,by,cx,cy,vbax,vbay,vba2,vbap,hpbax,hpbay,
     &	vbac,hcbax,hcbay,sba,vpax,vpay,vcax,vcay

	vbax=bx-ax
	vbay=by-ay
	vba2=vbax*vbax+vbay*vbay
	IF (vba2.LE.0d0) THEN
	  ChkLineSide=.FALSE.
	  RETURN
	END IF

	vpax=px-ax
	vpay=py-ay

	vbap=(vbax*vpax+vbay*vpay)/vba2
	hpbax=vpax-vbap*vbax
	hpbay=vpay-vbap*vbay

	vcax=cx-ax
	vcay=cy-ay

	vbac=(vbax*vcax+vbay*vcay)/vba2
	hcbax=vcax-vbac*vbax
	hcbay=vcay-vbac*vbay

	sba=hpbax*hcbax+hpbay*hcbay

	IF (sba.GE.0d0) THEN
	  ChkLineSide=.TRUE.
	ELSE
	  ChkLineSide=.FALSE.
	END IF

	RETURN
	END FUNCTION ChkLineSide

C***********************************************************************

	SUBROUTINE HOE(ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,UDSrfCoef,
     &	 UDSrfType,UDSrfParam,psi,pv,prot,pin,ihat,Evec,
     &	 rad,pmon,xhat,yhat,zhat,na,nb,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &	 hi,hr,Order,NomWavelength,
     &	 XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

	use traceutil_mod
        use elt_mod, only: mMonCoef,mAnaCoef

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifLNsrf,ifDMBuild,ifPol,ifObs,!ChkRayTrans,
     &	LRT

	INTEGER i,j,ip3,jp3,UDSrfType,SrfType,NDP,
     &	 IWORK(*),IERROR,mDP,
     &	 ApType,nObs,ObsType(nObs)

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &   pr(3),rhat(3),L,Nhat(3),pout(3),rout(3),MonCoef(mMonCoef),
     & 	 rad,pmon(3),xhat(3),yhat(3),zhat(3),
     &	 na,nb,kxa,kxb,lambda,
     &	 AnaCoef(mAnaCoef),UDSrfCoef(*),XYZSRF(3,NDP),DZSRF(2,NDP),DWORK(*),
     &	 xGrid(3),yGrid(3),zGrid(3),ApVec(6),
     &   ObsVec(6,nObs),UDSrfParam(6)
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 Nvec(3),Nmag,dNdp(3,3),dpdd(3,3),dNdd(3,3),dNdi(3,3),
     &   Ndoti,Niinv,Nouti(3,3),ioutN(3,3),drdN(3,3),
     &   dprot(3),dddth(3,3),Kc,Kr,
     &   NoidNdi(3,3),dNda(3,3),
     &	 dpdi(3,3),dpda(3,3),dpdth(3,3),Crn,inmrn(3),smag,
     &   svec(3),dsvecdN(3,3),dsdp(3,3),dGrdp(3),imagn,
     &	 rmagn,D4(3,3),drGdp(3,3),drsdp(3,3),drds(3,3),drdGr(3)

	REAL*8 drdi(3,3),drda(3,3),drdd(3,3),drdth(3,3),
     &   dgdi(3,3),dgda(3,3),dgdth(3,3),dgdd(3,3),
     &	 dLdi(3),dLda(3),dLdth(3),dLdd(3)

	REAL*8 dshatdsvec(3,3),dLivecdth(3,3),dsdN(3,3),dsNdp(3,3),
     &	 dindp(3,3),dindi(3,3),dinda(3,3),dindth(3,3),dindd(3,3),
     &	 drndp(3,3),drndi(3,3),drnda(3,3),drndth(3,3),drndd(3,3),
     &	 Perpin(3,3),Perprn(3,3),dsdrnmin(3,3),
     &	 dsdi(3,3),dsda(3,3),dsdth(3,3),dsdd(3,3),
     &	 dGrdi(3),dGrda(3),dGrdth(3),dGrdd(3),
     &	 dDLidi(3),dDLida(3),dDLidth(3),dDLidd(3),
     &	 dDLrdi(3),dDLrda(3),dDLrdth(3),dDLrdd(3),rhoi(3),rhor(3)

	REAL*8 Nx(3,3),Perpr(3,3),
     &	 D1(9),D2(9),D3(3,3),S1,S2,hi(3),hr(3),Order,
     &	 NomWavelength,shat(3),uhat(3),taur,uti,
     &	 prel(3),rho(3),ihatn(3),rhatn(3),psig(3),
     &	 dLr,dLi,Pnhat(3,3),sti,str,Gr,Lgeom

	COMPLEX*16 Evec(3),Eout(3)

cdebug
	COMMON/DEBUG/ shat,Gr,ihatn,rhatn,dsdi,dsda,dsdth,dsdd,
     &	 dGrdi,dGrda,dGrdth,dGrdd,drdi,drda,drdd,drdth,
     &	 dgdi,dgda,dgdd,dgdth,dLdi,dLda,dLdd,dLdth,
     &	 dindi,dinda,dindth,dindd,drndi,drnda,drndth,drndd
cdebug

C  Compute incidence point

	CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &	       prot,pin,ihat,pr,Lgeom,Nvec,Nmag,Nhat,dNdp)


C  Return if ray missed surface

	IF (.NOT.LROK) RETURN

	Nsrfm=Nvec/Nmag

	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.GT.0d0) THEN
	  CALL DNEGATE(dNdp,dNdp,9)
	END IF

C  Compute grating parameters

	S2=-SIGN(1d0,DDOTC(ihat,psi))
	DO 5 i=1,3
	  prel(i)=pin(i)-pv(i)
	  rho(i)=prel(i)+Lgeom*ihat(i)
	  psig(i)=S2*psi(i)
	  ihatn(i)=pr(i)-hi(i)
	  rhatn(i)=hr(i)-pr(i)
 5	CONTINUE
        CALL DUNITIZE(psig)

	IF ((hi(1)*hi(1)+hi(2)*hi(2)+hi(3)*hi(3)).LT.1d20) THEN
	  DO 6 i=1,3
	    D1(i)=pr(i)-hi(i)
	    D2(i)=pv(i)-hi(i)
 6	  CONTINUE
	  CALL DUNITVEC(ihatn,D1,imagn,3)
	  dLi=-(imagn-SQRT(D2(1)*D2(1)+D2(2)*D2(2)+D2(3)*D2(3)))
	ELSE
	  DO 7 i=1,3
	    D1(i)=pv(i)-hi(i)
 7	  CONTINUE
	  CALL DUNITVEC(ihatn,D1,imagn,3)
	  dLi=-DDOTC(rho,ihatn)
	END IF

	IF ((hr(1)*hr(1)+hr(2)*hr(2)+hr(3)*hr(3)).LT.1d20) THEN
	  DO 8 i=1,3
	    D1(i)=hr(i)-pr(i)
	    D2(i)=hr(i)-pv(i)
 8	  CONTINUE
	  CALL DUNITVEC(rhatn,D1,rmagn,3)
	  dLr=-(rmagn-SQRT(D2(1)*D2(1)+D2(2)*D2(2)+D2(3)*D2(3)))
	ELSE
	  DO 9 i=1,3
	    D1(i)=hr(i)-pv(i)
 9	  CONTINUE
	  CALL DUNITVEC(rhatn,D1,rmagn,3)
	  dLr=DDOTC(rho,rhatn)
	END IF

	CALL DPROJECT(Pnhat,Nhat)
	CALL DSUB(D1,rhatn,ihatn,3)
	CALL DMPROD(svec,Pnhat,D1,3,3,1)
	S1=svec(1)*svec(1)+svec(2)*svec(2)+svec(3)*svec(3)
	IF (S1.GE.1D-22) THEN
	  CALL DUNITVEC(shat,svec,smag,3)
	  CALL DXPROD(uhat,Nhat,shat)
          CALL DUNITIZE(uhat)
          sti=DDOTC(shat,ihatn)
          str=DDOTC(shat,rhatn)
          Gr=sti-str
c	  spacing=Order*NomWavelength/Gr
        ELSE
          Gr=0D0
          D1(1)=Nhat(3); D1(2)=Nhat(1); D1(3)=Nhat(2)
	  CALL DMPROD(svec,Pnhat,D1,3,3,1)
	  CALL DUNITVEC(shat,svec,smag,3)
	  smag=1d22
	  CALL DXPROD(uhat,Nhat,shat)
          CALL DUNITIZE(uhat)
        END IF

C  Compute diffracted ray

	L=Lgeom+dLr+dLi

	taur=DDOTC(shat,ihat)-Gr
	uti=DDOTC(uhat,ihat)
	S1=1d0-uti*uti-taur*taur
	IF (S1.LT.0D0) THEN
	  LROK=.FALSE.
	  RETURN
	END IF
	S1=SQRT(S1)
	DO 1 i=1,3
	  rhat(i)=taur*shat(i)+uti*uhat(i)+S1*Nhat(i)
   1	CONTINUE

C  Compute partials

	IF (ifBuild) THEN

C    Rotation-to-translation partials

	  CALL DSUB(dprot,prot,pr,3)
	  CALL DXMAT(dddth,dprot)

C    Partials of the incidence point

	  CALL DOUTER(ioutN,ihat,Nhat)
	  Ndoti=DDOTC(Nhat,ihat)
	  Niinv=1d0/Ndoti
	  CALL DSMPROD(NoidNdi,ioutN,Niinv,9)
	  CALL DEQUATE(dpdd,NoidNdi,9)
	  dpdd(1,1)=dpdd(1,1)-1d0
	  dpdd(2,2)=dpdd(2,2)-1d0
	  dpdd(3,3)=dpdd(3,3)-1d0

	  CALL DNEGATE(dpda,dpdd,9)

	  S1=Lgeom
	  CALL DSMPROD(dpdi,dpda,S1,9)

	  CALL DMPROD(dpdth,dpdd,dddth,3,3,3)

	  CALL DMPROD(dLivecdth,NoidNdi,dddth,3,3,3)

C    Partials of the normal

	  CALL DMPROD(dNdd,dNdp,dpdd,3,3,3)

	  S1=-Lgeom
	  CALL DSMPROD(dNdi,dNdd,S1,9)

	  CALL DNEGATE(dNda,dNdd,9)

C    Grating partials

	  CALL DPROJECT(Pnhat,Nhat)
	  CALL DXMAT(Nx,Nhat)
	  CALL DSUB(inmrn,ihatn,rhatn,3)

	  CALL DPROJECT(Perprn,rhatn)
	  S1=-1d0/rmagn
	  CALL DSMPROD(drndp,Perprn,S1,9)
	  CALL DMPROD(drndi,drndp,dpdi,3,3,3)

	  CALL DMPROD(drnda,drndp,dpda,3,3,3)

	  CALL DMPROD(drndd,drndp,dpdd,3,3,3)

	  CALL DSUB(D1,hr,prot,3)
	  CALL DXMAT(D2,D1)
	  CALL DADD(D3,dLivecdth,D2,9)
	  CALL DMPROD(drndth,drndp,D3,3,3,3)

	  CALL DPROJECT(Perpin,ihatn)
	  S1=1d0/imagn
	  CALL DSMPROD(dindp,Perpin,S1,9)
	  CALL DMPROD(dindi,dindp,dpdi,3,3,3)

	  CALL DMPROD(dinda,dindp,dpda,3,3,3)

	  CALL DMPROD(dindd,dindp,dpdd,3,3,3)
c-FIXME!!!
	  CALL DSUB(D1,hi,prot,3)
	  CALL DXMAT(D2,D1)
	  CALL DMPROD(D1,dindp,D2,3,3,3)
	  CALL DMPROD(D2,dindd,dddth,3,3,3)
	  CALL DADD(dindth,D1,D2,9)
c-end FIXME
	  CALL DPROJECT(D1,shat)
	  S1=1d0/smag
	  CALL DSMPROD(dshatdsvec,D1,S1,9)

	  CALL DOUTER(D3,Nhat,inmrn)
	  S1=DDOTC(Nhat,inmrn)
	  D3(1,1)=D3(1,1)+S1
	  D3(2,2)=D3(2,2)+S1
	  D3(3,3)=D3(3,3)+S1
	  CALL DEQUATE(dsvecdN,D3,9)

	  CALL DMPROD(dsdN,dshatdsvec,dsvecdN,3,3,3)

	  CALL DMPROD(dsNdp,dsdN,dNdp,3,3,3)
	  CALL DMPROD(dsdrnmin,dshatdsvec,Pnhat,3,3,3)

	  CALL DSUB(D2,drndp,dindp,9)
	  CALL DMPROD(D3,dsdrnmin,D2,3,3,3)
	  CALL DADD(dsdp,dsNdp,D3,9)

	  CALL DMPROD(D1,dsNdp,dpdi,3,3,3)
	  CALL DSUB(D2,drndi,dindi,9)
	  CALL DMPROD(D3,dsdrnmin,D2,3,3,3)
	  CALL DADD(dsdi,D1,D3,9)

	  CALL DMPROD(D1,dsNdp,dpda,3,3,3)
	  CALL DSUB(D2,drnda,dinda,9)
	  CALL DMPROD(D3,dsdrnmin,D2,3,3,3)
	  CALL DADD(dsda,D1,D3,9)

	  CALL DMPROD(D2,dsNdp,dpdth,3,3,3)
	  CALL DMPROD(D3,dsdN,Nx,3,3,3)
	  CALL DSUB(D1,D2,D3,9)
	  CALL DSUB(D2,drndth,dindth,9)
	  CALL DMPROD(D3,dsdrnmin,D2,3,3,3)
	  CALL DADD(dsdth,D1,D3,9)

	  CALL DMPROD(D1,dsNdp,dpdd,3,3,3)
	  CALL DSUB(D2,drndd,dindd,9)
	  CALL DMPROD(D3,dsdrnmin,D2,3,3,3)
	  CALL DADD(dsdd,D1,D3,9)

	  CALL DSUB(D1,dindp,drndp,9)
	  CALL DMPROD(D2,shat,D1,1,3,3)
	  CALL DMPROD(D1,inmrn,dsdp,1,3,3)
	  CALL DADD(dGrdp,D1,D2,3)

	  CALL DSUB(D1,dindi,drndi,9)
	  CALL DMPROD(D2,shat,D1,1,3,3)
	  CALL DMPROD(D1,inmrn,dsdi,1,3,3)
	  CALL DADD(dGrdi,D1,D2,3)

	  CALL DSUB(D1,dinda,drnda,9)
	  CALL DMPROD(D2,shat,D1,1,3,3)
	  CALL DMPROD(D1,inmrn,dsda,1,3,3)
	  CALL DADD(dGrda,D1,D2,3)

	  CALL DSUB(D1,dindth,drndth,9)
	  CALL DMPROD(D2,shat,D1,1,3,3)
	  CALL DMPROD(D1,inmrn,dsdth,1,3,3)
	  CALL DADD(dGrdth,D1,D2,3)

	  CALL DSUB(D1,dindd,drndd,9)
	  CALL DMPROD(D2,shat,D1,1,3,3)
	  CALL DMPROD(D1,inmrn,dsdd,1,3,3)
	  CALL DADD(dGrdd,D1,D2,3)

C    Ray direction partials

	  Crn=SQRT(Ndoti*Ndoti+2d0*Gr*DDOTC(shat,ihat)-Gr*Gr)
	  CALL DOUTER(Nouti,Nhat,ihat)
	  S1=Ndoti/Crn-1d0
	  CALL DSMPROD(drdN,Nouti,S1,9)
	  S1=Crn-Ndoti
	  drdN(1,1)=drdN(1,1)+S1
	  drdN(2,2)=drdN(2,2)+S1
	  drdN(3,3)=drdN(3,3)+S1

	  S1=Ndoti/Crn-1d0
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DOUTER(D2,Nhat,D1)
	  S1=Gr/Crn
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DOUTER(D3,D1,shat)
	  CALL DADD(drdi,D2,D3,9)
	  drdi(1,1)=drdi(1,1)+1d0
	  drdi(2,2)=drdi(2,2)+1d0
	  drdi(3,3)=drdi(3,3)+1d0

	  S1=(DDOTC(shat,ihat)-Gr)/Crn
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DSUB(drdGr,D1,shat,3)

	  S1=Gr/Crn
	  CALL DSMPROD(drds,Nouti,S1,9)
	  drds(1,1)=drds(1,1)-Gr
	  drds(2,2)=drds(2,2)-Gr
	  drds(3,3)=drds(3,3)-Gr

	  CALL DMPROD(drsdp,drds,dsdp,3,3,3)
	  CALL DMPROD(D2,drsdp,dpdi,3,3,3)
	  CALL DMPROD(D3,drdN,dNdi,3,3,3)
	  CALL DADD(D1,D2,D3,9)
	  CALL DADD(D3,drdi,D1,9)
	  CALL DMPROD(drGdp,drdGr,dGrdp,3,1,3)
	  CALL DMPROD(D2,drGdp,dpdi,3,3,3)
cshouldbe	  CALL DADD(drdi,D2,D3,9)
	  CALL DADD(D1,D2,D3,9)
	  CALL DPROJECT(D2,ihat)
	  CALL DMPROD(drdi,D1,D2,3,3,3)

	  CALL DMPROD(D2,drsdp,dpda,3,3,3)
	  CALL DMPROD(D3,drdN,dNda,3,3,3)
	  CALL DADD(D1,D2,D3,9)
	  CALL DMPROD(D2,drGdp,dpda,3,3,3)
	  CALL DADD(drda,D1,D2,9)

	  CALL DMPROD(D1,drsdp,dpdd,3,3,3)
	  CALL DMPROD(D2,drGdp,dpdd,3,3,3)
	  CALL DADD(D3,D2,D1,9)
	  CALL DMPROD(D1,drdN,dNdp,3,3,3)
	  CALL DMPROD(D2,D1,dpdd,3,3,3)
	  CALL DADD(drdd,D2,D3,9)

	  CALL DMPROD(D1,drdd,dddth,3,3,3)
	  CALL DMPROD(D2,drdN,Nx,3,3,3)
	  CALL DSUB(D3,D1,D2,9)
	  CALL DXMAT(D1,shat)
	  CALL DMPROD(D2,drds,D1,3,3,3)
	  CALL DSUB(drdth,D3,D2,9)

C    Beamwalk partials

	  CALL DPROJECT(Perpr,rhat)
	  CALL DMPROD(D1,Perpr,dpdd,3,3,3)
	  CALL DNEGATE(dgda,D1,9)
	  CALL DSMPROD(dgdi,dgda,Lgeom,9)
	  S1=1d0/Ndoti
	  CALL DSMPROD(D1,Nhat,S1,3)
	  CALL DMPROD(D2,Perpr,ihat,3,3,1)
	  CALL DMPROD(dgdd,D2,D1,3,1,3)
	  CALL DMPROD(dgdth,dgdd,dddth,3,3,3)

C    Pathlength partials

	  CALL DNEGATE(D1,ihatn,3)
	  CALL DMPROD(dDLidi,D1,dpdi,1,3,3)
	  CALL DMPROD(dDLida,D1,dpda,1,3,3)
	  CALL DMPROD(dDLidd,D1,dpdd,1,3,3)
	  CALL DSUB(D1,hi,prot,3)
	  CALL DXMAT(D2,D1)
	  CALL DADD(D1,D2,dLivecdth,9)
	  CALL DMPROD(D3,ihatn,D1,1,3,3)
	  CALL DSUB(rhoi,hi,pv,3)
	  CALL DXMAT(D1,rhoi)
	  CALL DNEGATE(D2,rhoi,3)
	  CALL DUNITIZE(D2)
	  CALL DMPROD(D4,D2,D1,1,3,3)
	  CALL DSUB(dDLidth,D4,D3,3)

	  CALL DMPROD(dDLrdi,rhatn,dpdi,1,3,3)
	  CALL DMPROD(dDLrda,rhatn,dpda,1,3,3)
	  CALL DMPROD(dDLrdd,rhatn,dpdd,1,3,3)
	  CALL DSUB(D1,hr,prot,3)
	  CALL DXMAT(D2,D1)
	  CALL DADD(D1,D2,dLivecdth,9)
	  CALL DMPROD(D3,rhatn,D1,1,3,3)
	  CALL DSUB(rhor,hr,pv,3)
	  CALL DXMAT(D1,rhor)
	  CALL DEQUATE(D2,rhor,3)
	  CALL DUNITIZE(D2)
	  CALL DMPROD(D4,D2,D1,1,3,3)
	  CALL DSUB(dDLrdth,D3,D4,3)

	  CALL DSUB(D3,rhat,ihat,3)
	  CALL DMPROD(D1,D3,dpdd,1,3,3)
	  CALL DADD(D2,D1,dDLida,3)
	  CALL DADD(dLda,D2,dDLrda,3)

	  CALL DSMPROD(D2,D1,Lgeom,3)
	  CALL DADD(D2,D2,dDLidi,3)
	  CALL DADD(dLdi,D2,dDLrdi,3)

	  S1=-DDOTC(D3,ihat)/Ndoti
	  CALL DSMPROD(dLdd,Nhat,S1,3)
	  CALL DADD(dLdd,dLdd,dDLidd,3)
	  CALL DADD(dLdd,dLdd,dDLrdd,3)

	  CALL DMPROD(dLdth,dLdd,dddth,1,3,3)

C    Combined partial matrices

	  DO 4 i=1,3
	    ip3=i+3
	    DO 3 j=1,3
	      jp3=j+3
	      dxidxim1(i,j)=drdi(i,j)
	      dxidxim1(i,jp3)=drda(i,j)
	      dxidxim1(ip3,j)=dgdi(i,j)
	      dxidxim1(ip3,jp3)=dgda(i,j)
	      dxidui(i,j)=drdth(i,j)
	      dxidui(i,jp3)=drdd(i,j)
	      dxidui(ip3,j)=dgdth(i,j)
	      dxidui(ip3,jp3)=dgdd(i,j)
    3	    CONTINUE
	    dxidxim1(i,7)=0e0
	    dxidxim1(ip3,7)=0e0
	    dxidxim1(7,i)=0E0
	    dxidxim1(7,ip3)=0E0
	    dxidui(7,i)=dLdth(i)
	    dxidui(7,ip3)=dLdd(i)
    4	  CONTINUE
	  dxidxim1(7,7)=1d0
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	         psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) THEN
	    Eout(1)=0d0
	    Eout(2)=0d0
	    Eout(3)=0d0
	  END IF
	ELSE
	  LRT=.TRUE.
	END IF

C  Set output values for ray position, direction
        rout(1:3) = rhat(1:3)
        pout(1:3) = pr(1:3)
	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE HOE


C***********************************************************************
C For reflective grating - jzlou

	SUBROUTINE Grating(ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,UDSrfCoef,UDSrfType,
     &	 UDSrfParam,psi,pv,prot,pin,ihat,Evec,
     &	 rad,pmon,xhat,yhat,zhat,
     &	 na,nb,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &	 Order,RuleWidth,RuleDir, ! RuleDir = h1HOE(1,iElt)
     &	 XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

	use didesub, only: IntSrf
        use traceutil_mod
        use elt_mod, only: mMonCoef,mAnaCoef

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,LRT !ChkRayTrans

	INTEGER i,UDSrfType,SrfType,NDP,IWORK(*),IERROR,
     &	  mDP,ApType,nObs,ObsType(nObs)

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &         pr(3),rhat(3),L,Nhat(3),pout(3),rout(3),lambda,
     &	       na,nb,kxa,kxb,AnaCoef(mAnaCoef),MonCoef(mMonCoef),
     &	       RuleWidth,RuleDir(3),rad,pmon(3),xhat(3),yhat(3),zhat(3),
     &	       xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs),
     &	       UDSrfCoef(*),UDSrfParam(6)
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 Nvec(3),Nmag,dNdp(3,3),its,
     &         dSig,strho,ihatn(3),rhatn(3),
     &         XYZSRF(3,NDP),DZSRF(2,NDP),wVec(3),
     &	       DWORK(*),Kc,Kr

	REAL*8 drdi(3,3),drda(3,3),drdd(3,3),drdth(3,3),
     &         dgdi(3,3),dgda(3,3),dgdth(3,3),dgdd(3,3),
     &	       dLdi(3),dLda(3),dLdth(3),dLdd(3)

	REAL*8 dindi(3,3),dinda(3,3),dindth(3,3),dindd(3,3),
     &	       drndi(3,3),drnda(3,3),drndth(3,3),drndd(3,3),
     &	       dsdi(3,3),dsda(3,3),dsdth(3,3),dsdd(3,3),
     &	       dGrdi(3),dGrda(3),dGrdth(3),dGrdd(3)

	REAL*8 S1,Order,shat(3),uhat(3),rts,uti,
     &	       prel(3),rho(3),Gr

	COMPLEX*16 Evec(3),Eout(3)

cdebug
	COMMON/DEBUG/ shat,Gr,ihatn,rhatn,dsdi,dsda,dsdth,dsdd,
     &	 dGrdi,dGrda,dGrdth,dGrdd,drdi,drda,drdd,drdth,
     &	 dgdi,dgda,dgdd,dgdth,dLdi,dLda,dLdd,dLdth,
     &	 dindi,dinda,dindth,dindd,drndi,drnda,drndth,drndd
cdebug

        !print*,'***** in Reflective Grating ...'

C  Compute incidence point of ray

	IF (SrfType.EQ.1) THEN
	  CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &	  	pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.2) THEN
	  CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &	  	prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.3) THEN
	  CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &	  	prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
	ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
	  CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &	  	prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &	  	Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.5) THEN
	  CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &	  	XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	  	wVec,pr,L,pmon,xhat,yhat,zhat,
     &	  	Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.6) THEN
	  CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  	prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.7) THEN
	  CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &	  	UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,pmon,
     &          xhat,yhat,zhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.10) THEN
	  CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  	prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE
	  WRITE(*,*)'4: Unknown surface type. Aborting trace.'
	END IF

#if 0
! Test code
	pin(1)=RuleWidth/3d0
	!pin(1)=RuleWidth/2d0
	!pin(1)=RuleWidth
	print*,'--Grating: pv =',pv(1:3)
	print*,'--Grating: pin =',pin(1:3)
	print*,'--Grating: ihat =',ihat(1:3)
	print*,'--Grating: Nvec =',Nvec(1:3)
	print*,'--Grating: L =',L
#endif

C  Return if ray missed surface

	IF (.NOT.LROK) THEN
          RETURN
        END IF

	Nsrfm=Nvec/Nmag

	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.GT.0d0) THEN
	  CALL DNEGATE(dNdp,dNdp,9)
	END IF

C  Compute grating parameters -- see Siegman pg. 698 et. seq.

	prel(1:3)=pin(1:3)-pv(1:3)
	rho(1:3)=prel(1:3)+L*ihat(1:3)
          ! rho = (pin+L*ihat)-pv = vector from elt vertex to ray
          ! intersection point at current surface.

        !print*,'**** Order =',Order
        !print*,'**** lambda =',lambda
        !print*,'****  RuleWidth =',RuleWidth
        !stop

C  for fixed rule width: -dcr
        Gr=Order*lambda/RuleWidth
          ! Gr is the quantity pre-computed for applying grating eqn later.

C  this implements the grating equation: -dcr
        CALL DPERPUVEC(RuleDir,RuleDir,psi)

	CALL DXPROD(uhat,Nhat,RuleDir)
          ! Nhat = (negative) local srf normal
          ! uhat is in direction along rule line
        CALL DUNITIZE(uhat)
	CALL DXPROD(shat,uhat,Nhat)
          ! shat is the direction along grating (substrate) surface,
          ! and perpendicular to the rule lines -jzlou
        CALL DUNITIZE(shat)

	its=DDOTC(shat,ihat) ! cos of angle between incident ray and shat, which
                             ! equals to sin of angle betwen incident ray and
                             ! (substrate) surface normal -jzlou
	rts=its-Gr  ! **** this is the use of grating equation; ****
                    ! rts is the sin of angle between reflected ray and
                    ! (substrate) surface normal -jzlou
	uti=DDOTC(uhat,ihat)
	S1=1d0-uti*uti-rts*rts
	IF (S1.LT.0D0) THEN  ! original
	!IF (S1.LT.-1D-8) THEN
	  LROK=.FALSE.
          print*,'Gr =', Gr
          print*,'uhat =', uhat
          print*,'shat =', shat
          print*,'ihat =', ihat
          print*,'its =',its
          print*,'rts =',rts
          print*,'** S1 =',S1
          print*,'** S1 negativ, Ray does NOT pass grating ...'
	  RETURN
	END IF
	DO i=1,3
	  rhat(i)=rts*shat(i)+uti*uhat(i)+DSQRT(S1)*Nhat(i)
    	END DO

        If (.false.) Then
         !print*,'*** RuleDir =',RuleDir
         print*,' ********** Gr =', Gr
         print*,' ******* ihat=',ihat
         print*,' ******* shat=',shat
         print*,' ******* uhat=',uhat
         print*,' ******* Nhat=',Nhat
         print*,' ***** its = ',its
         print*,' ***** rts = ',rts
         print*,'***** S1 = ',S1
         !print*,'***** rhat = ',rhat(1:3)
         !stop
        End If

	strho=DDOTC(RuleDir,rho) ! strho = distance of ray to grating vertex projected onto
                                 ! direction of RuleDir
	dSig=MOD(strho,RuleWidth) ! dSig is strho moduled by grating step size
	IF (dSig.LT.0d0) dSig=dSig+RuleWidth
	L=L-Gr*dSig ! 'Gr*dSig' is the 'extra part' needed to be
                    ! subtracted from path-length 'L' to make this
                    ! ray coherent with the ray hitting grating vertex.

        ! Derivation of "L=L-Gr*dSig"   -jzlou
        ! Path difference = a*(sin(I)+sin(D));
        !   a = the separation distance between the two parallel rays;
        !   I = angle between incident rays and grating (base surface) normal;
        !   D = angle between reflected rays and grating (base surface) normal;
        ! Now:
        ! a*(sin(I)+sin(D)) = a*(sin(I)+sin(D))*RuleWidth/RuleWidth
        !  = a*(Order*lambda)/RuleWidth = a*Gr = Gr*dSig = path difference of the two rays

        ! Implication of "L=L-Gr*dsig"   -jzlou
        ! This makes sure any ray 'dSig' away from a grove starting point is conherent with
        ! the ray hitting the grove starting point. For rays hitting the grove starting
        ! or end point, the path differnce between any two of them would be
        ! k*RuleWith*(sin(I)+sin(D)) = k*Order*lambda, so the two rays are phase conherent,
        ! where assuming the two rays are k groves apart.

#if 0
	print*,'--Grating: dSig =',dSig
	print*,'--Grating: Gr*dSig =',Gr*dSig
	print*,'--Grating: L_final =',L
	print*,'--Grating: rhat =',rhat(1:3)
	!stop
#endif

C  Compute partials

	IF (ifBuild) THEN
C         Zero partial matrices
          dxidxim1(1:7,1:7)=0d0
          dxidui(1:7,1:6)=0d0
        END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	         psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) Eout(1:3)=0d0
	END IF
	IF ((SrfType.EQ.5).AND.(IERROR.NE.0)) THEN
	  LRT=.FALSE.
	ELSE
	  LRT=.TRUE.
	END IF

C  Set output values for ray position, direction
	rout(1:3)=rhat(1:3)
	pout(1:3)=pr(1:3)
	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE Grating

C***********************************************************************

! Transmission grating (used in DFS) -- jzlou

	SUBROUTINE TrGrating(ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,
     &	 LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,UDSrfCoef,UDSrfType,
     &	 UDSrfParam,psi,pv,prot,pin,ihat,Evec,
     &	 rad,pmon,xhat,yhat,zhat,
     &	 na,nb,kxa,kxb,lambda,
     &	 ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &	 pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &	 Order,RuleWidth,RuleDir,
     &	 XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

! h1HOE(1:3,iElt) = RuleDir(1:3) = groove srf tangent, -jzlou

	use Constants, only: PI
        use traceutil_mod
	use didesub, only: IntSrf
        use elt_mod, only: mMonCoef,mAnaCoef

	IMPLICIT NONE

	LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &          LRT

	INTEGER i,j,UDSrfType,SrfType,
     &	  NDP,IWORK(*),IERROR,mDP,
     &	  ApType,nObs,ObsType(nObs)

	SREAL :: dxidxim1(7,7),dxidui(7,6)

	REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &   pr(3),rhat(3),L,Nhat(3),pout(3),rout(3),lambda,
     &	 na,nb,kxa,kxb,AnaCoef(mAnaCoef),MonCoef(mMonCoef),
     &	 RuleWidth,RuleDir(3),rad,pmon(3),xhat(3),yhat(3),zhat(3),
     &	 xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs),
     &	 UDSrfCoef(*),UDSrfParam(6)
	REAL*8, dimension(:) :: AsphCoef

	REAL*8 Nvec(3),Nmag,dNdp(3,3),its,
     &   dSig,strho,ihatn(3),rhatn(3),XYZSRF(3,NDP),DZSRF(2,NDP),
     &	 wVec(3),DWORK(*),Kc,Kr

	REAL*8 drdi(3,3),drda(3,3),drdd(3,3),drdth(3,3),
     &         dgdi(3,3),dgda(3,3),dgdth(3,3),dgdd(3,3),
     &	       dLdi(3),dLda(3),dLdth(3),dLdd(3)

	REAL*8 dindi(3,3),dinda(3,3),dindth(3,3),dindd(3,3),
     &	       drndi(3,3),drnda(3,3),drndth(3,3),drndd(3,3),
     &	       dsdi(3,3),dsda(3,3),dsdth(3,3),dsdd(3,3),
     &	       dGrdi(3),dGrda(3),dGrdth(3),dGrdd(3)

	REAL*8 Ndi,cvec(3),muNdi,cfb,mu,dhat(3),
     &	 S1,S2,Order,shat(3),uhat(3),rts,uti,
     &	 prel(3),rho(3),Gr

	COMPLEX*16 Evec(3),Eout(3)

cdebug
	COMMON/DEBUG/ shat,Gr,ihatn,rhatn,dsdi,dsda,dsdth,dsdd,
     &	 dGrdi,dGrda,dGrdth,dGrdd,drdi,drda,drdd,drdth,
     &	 dgdi,dgda,dgdd,dgdth,dLdi,dLda,dLdd,dLdth,
     &	 dindi,dinda,dindth,dindd,drndi,drnda,drndth,drndd
cdebug

C  Compute incidence point of ray

	IF (SrfType.EQ.1) THEN
	  CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &	  pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.2) THEN
	  CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &	  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.3) THEN
	  CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
	ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
	  CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.5) THEN
	  CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &	  XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &	  wVec,pr,L,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.6) THEN
	  CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.7) THEN
	  CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &	  UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &	  Nvec,Nmag,Nhat,dNdp)
	ELSE IF (SrfType.EQ.10) THEN
	  CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &	  prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
	ELSE
	  WRITE(*,*)'5: Unknown surface type. Aborting trace.'
	END IF


C  Return if ray missed surface

	IF (.NOT.LROK) RETURN

	Nsrfm=Nvec/Nmag

	S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
	IF (S1.GT.0d0) THEN
	  CALL DNEGATE(dNdp,dNdp,9)
	END IF


C  Compute grating parameters

	DO i=1,3
	  prel(i)=pin(i)-pv(i)
	  rho(i)=prel(i)+L*ihat(i) ! rho = (pin+L*ihat)-pv is vector from elt vertex to ray
                                   ! intersection point at current surface.
  	END DO

C  for fixed rule width: -dcr
        Gr=Order*lambda/RuleWidth  ! = m*lambda/a, where a is grating period

        CALL DPERPUVEC(RuleDir,RuleDir,psi)  ! When RuleDir and psi are othogonal to each other,
				             ! the routine doesn't change RuleDir.
	!print*,' **TrGrating: aft RuleDir =',RuleDir

	CALL DXPROD(uhat,Nhat,RuleDir)    ! Nhat = (negative) local srf normal
					  ! uhat will be 'out of paper' direction

	CALL DUNITIZE(uhat)
	CALL DXPROD(shat,uhat,Nhat)
	CALL DUNITIZE(shat)   ! shat othogonal to Nhat and uhat (in paper)

	its=DDOTC(shat,ihat)  ! cos(th), th = angle betw. ihat and shat
	rts=its-Gr
	uti=DDOTC(uhat,ihat)
	S1=1d0-uti*uti-rts*rts

	IF (S1.LT.0D0) THEN
	  LROK=.FALSE.
	  RETURN
	END IF
	DO i=1,3
	  dhat(i)=rts*shat(i)+uti*uhat(i)-SQRT(S1)*Nhat(i)
        END DO

C	strho=DDOTC(RuleDir,rho)
	strho=DDOTC(shat,rho)
	dSig=MOD(strho,RuleWidth)
	IF (dSig.LT.0d0) dSig=dSig+RuleWidth
C	L=L-Gr*dSig
	L=L+(rts-its)*dSig


C  Compute refracted ray

	mu=na/nb
	Ndi=dhat(1)*Nhat(1)+dhat(2)*Nhat(2)+dhat(3)*Nhat(3)
	muNdi=mu*Ndi
	DO i=1,3
	  cvec(i)=Ndi*Nhat(i)-dhat(i)
        END DO

C  Check for total internal reflection

	S2=1d0-mu*mu+muNdi*muNdi
#if 0
	print*,' **TrGrating: S2 =',S2
	stop
#endif

	IF (S2.LT.0d0) THEN
	  LROK=.FALSE.
	  RETURN
	END IF


	cfb=DSQRT(S2)
	S1=0d0
	DO i=1,3
	  rhat(i)=-mu*cvec(i)-cfb*Nhat(i)
	  S1=S1+rhat(i)*rhat(i)
        END DO
	S2=DSQRT(S1)
	DO i=1,3
	  rhat(i)=rhat(i)/S2
        END DO

	! The direction of refracted ray (rhat) obeys grating equation in
	! the case of normal incident ray, a sanity check by jzlou.
	! The specific direction of h1HOE (along each groove surface) does
	! not matter in determining the direction of refracted ray emerging
	! from the transmission grating. -jzlou

#if 0
	! Pure test
	print*,' **TrGrating: rhat =',rhat(1:3)
        S2=DDOTC(rhat,ihat)
	S2=sqrt(1d0-S2*S2)
	print*,' **TrGrating: sin(beta) = ',S2
	print*,' **TrGrating: Gr =',Gr
	stop
#endif

C  Compute partials

	IF (ifBuild) THEN

C    Zero partial matrices

	  DO 7 j=1,7
	    DO 6 i=1,7
	      dxidxim1(i,j)=0d0
    6	    CONTINUE
    7	  CONTINUE
	  DO 9 j=1,6
	    DO 8 i=1,7
	      dxidui(i,j)=0d0
    8	    CONTINUE
    9	  CONTINUE
	END IF

C  Check if ray is obscured
	IF (ifObs) THEN
	  LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &	         psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
	  IF (.NOT.LRT) THEN
	    Eout(1:3)=0d0
	  END IF
	END IF
	IF ((SrfType.EQ.5).AND.(IERROR.NE.0)) THEN
	  LRT=.FALSE.
	ELSE
	  LRT=.TRUE.
	END IF

C  Set output values for ray position, direction
	rout(1:3)=rhat(1:3)
	pout(1:3)=pr(1:3)

	NhatP(1:3)=Nhat(1:3)

	RETURN
	END SUBROUTINE TrGrating


C***********************************************************************
! Currently not used.

	SUBROUTINE FzpTrGrating_(ifBuild,ifDMBuild,ifPol,
     &   ifLNsrf,ifObs,LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,
     &   UDSrfCoef,UDSrfType,UDSrfParam,psi,pv,prot,pin,ihat,Evec,
     &   rad,pmon,xhat,yhat,zhat,
     &   na,nb,kxa,kxb,lambda,
     &   ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &   pout,rout,L,Eout,dxidxim1,dxidui,SrfType,
     &   Order,RuleWidth,RuleDir,
     &   XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

! h1HOE(1:3,iElt) = RuleDir(1:3) = groove srf tangent, -jzlou

	use Constants, only: PI
        use elt_mod, only : DoeFl,mMonCoef,mAnaCoef
        use traceutil_mod
	use didesub, only: IntSrf

        IMPLICIT NONE

        LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &          LRT

        INTEGER i,j,UDSrfType,SrfType,
     &    NDP,IWORK(*),IERROR,mDP,
     &    ApType,nObs,ObsType(nObs)

        SREAL :: dxidxim1(7,7),dxidui(7,6)

        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &   pr(3),rhat(3),L,Nhat(3),pout(3),rout(3),lambda,
     &   na,nb,kxa,kxb,AnaCoef(mAnaCoef),MonCoef(mMonCoef),
     &   RuleWidth,RuleDir(3),rad,pmon(3),xhat(3),yhat(3),zhat(3),
     &   xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs),
     &   UDSrfCoef(*),UDSrfParam(6)
        REAL*8, dimension(:) :: AsphCoef

        REAL*8 Nvec(3),Nmag,dNdp(3,3),its,
     &   dSig,strho,ihatn(3),rhatn(3),XYZSRF(3,NDP),DZSRF(2,NDP),
     &   wVec(3),DWORK(*),Kc,Kr

        REAL*8 drdi(3,3),drda(3,3),drdd(3,3),drdth(3,3),
     &   dgdi(3,3),dgda(3,3),dgdth(3,3),dgdd(3,3),resid,d1_loc,
     &   dLdi(3),dLda(3),dLdth(3),dLdd(3)

        REAL*8 dindi(3,3),dinda(3,3),dindth(3,3),dindd(3,3),
     &         drndi(3,3),drnda(3,3),drndth(3,3),drndd(3,3),
     &         dsdi(3,3),dsda(3,3),dsdth(3,3),dsdd(3,3),
     &         dGrdi(3),dGrda(3),dGrdth(3),dGrdd(3)

        REAL*8 Ndi,cvec(3),muNdi,cfb,mu,dhat(3),pdd,
     &   S1,S2,DoeFocLen,focalOrder,
     &   Order,shat(3),uhat(3),rts,uti,CurRuleWidth,CurRuleDir(3),
     &   prel(3),rho(3),Gr

        COMPLEX*16 :: Evec(3),Eout(3),
     &                Calc_RayDielectricTC,Calc_RayMetalTC

cdebug
        COMMON/DEBUG/ shat,Gr,ihatn,rhatn,dsdi,dsda,dsdth,dsdd,
     &    dGrdi,dGrda,dGrdth,dGrdd,drdi,drda,drdd,drdth,
     &    dgdi,dgda,dgdd,dgdth,dLdi,dLda,dLdd,dLdth,
     &    dindi,dinda,dindth,dindd,drndi,drnda,drndth,drndd
cdebug

C  Compute incidence point of ray

        IF (SrfType.EQ.1) THEN
          CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &    pr,L,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.2) THEN
          CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &    prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.3) THEN
          CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &    prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
        ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
          CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &    prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &    Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.5) THEN
          CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &    XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &    wVec,pr,L,pmon,xhat,yhat,zhat,
     &    Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.6) THEN
          CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &    prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.7) THEN
          CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &    UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &    Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.10) THEN
          CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &    prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
        ELSE
          WRITE(*,*)'5: Unknown surface type. Aborting trace.'
        END IF

C  Return if ray missed surface
        IF (.NOT.LROK) RETURN

        Nsrfm=Nvec/Nmag

        S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
        IF (S1.GT.0d0) THEN
          CALL DNEGATE(dNdp,dNdp,9)
        END IF

C  Compute grating parameters

        DO i=1,3
          prel(i)=pin(i)-pv(i)
          rho(i)=prel(i)+L*ihat(i) ! rho = (pin+L*ihat)-pv is vector from elt vertex to ray
                                   ! intersection point at current surface.
        END DO

!
!  To implement a circularily ruled surface, need to make RuleDir in the direction from
!  element vertex to ray intersection point. -jzlou
        CurRuleDir(1:3)=rho(1:3)
        CALL DUNITIZE(CurRuleDir)

	! Fresnel zone plate focal length
	!DoeFocLen=5d3 ! 5 meters
	DoeFocLen=DoeFl(crnt_Elt)

!  RuleWidth value is computed to generate a ray direction so that the ray will focus
!  'DoeFocLen' distance from the Fresnel plate center -jzlou

	!Order=2  ! order probably needs to be positive to have rays converging to focus
	!focalOrder=1
	focalOrder=Order

	pdd=abs(DDOTC(rho,CurRuleDir)) ! distance from ray to plate center
        If (pdd<1d-18) then
	  ! needed to avoid divide by zero if is chief ray
	  CurRuleWidth=lambda*1d8 ! so sin(theta_m) ~= 0
	Else
	  CurRuleWidth=abs(focalOrder*lambda/sin(atan(pdd/DoeFocLen)))
        End If
!	print*,'**abs(DDOTC(rho,CurRuleDir)),DoeFocLen =',
!    &         abs(DDOTC(rho,CurRuleDir)),DoeFocLen
!       print*,'**CurRuleWidth =',CurRuleWidth

        Gr=Order*lambda/CurRuleWidth  ! = N*lambda/a = sin(theta_N)
	!print*,'**Gr =',Gr

        CALL DPERPUVEC(CurRuleDir,CurRuleDir,psi) ! When RuleDir and psi are othogonal to
                                         ! each other, the routine doesn't change RuleDir.
        !print*,' **TrGrating: aft CurRuleDir =',CurRuleDir

        CALL DXPROD(uhat,Nhat,CurRuleDir)    ! Nhat = (negative) local srf normal
                                             ! uhat will be 'out of paper' direction

        CALL DUNITIZE(uhat)
        CALL DXPROD(shat,uhat,Nhat)
        CALL DUNITIZE(shat)   ! shat othogonal to Nhat and uhat (in paper)

        its=DDOTC(shat,ihat)  ! cos(th), th = angle betw. ihat and shat
        rts=its-Gr
        uti=DDOTC(uhat,ihat)
        S1=1d0-uti*uti-rts*rts

        IF (S1.LT.0D0) THEN
          LROK=.FALSE.
          RETURN
        END IF
        DO i=1,3
          dhat(i)=rts*shat(i)+uti*uhat(i)-SQRT(S1)*Nhat(i)
        END DO

C       strho=DDOTC(RuleDir,rho)
        strho=DDOTC(shat,rho)
        dSig=MOD(strho,CurRuleWidth)
        IF (dSig.LT.0d0) dSig=dSig+CurRuleWidth
C       L=L-Gr*dSig
        L=L+(rts-its)*dSig ! Ray optical path changes to make sure the cross-section of
                           ! diffracted beam has coherent phase (subtract 'extra' path
                           ! length for each ray) -jzlou

#if 1
	! Now make ray path length conherent with the chief ray - difference of them
	! needs be multiple of wavelength
	d1_loc = dsqrt(pdd**2+DoeFocLen**2)+(rts-its)*dSig
	!d1_loc = dsqrt(pdd**2+DoeFocLen**2)
          ! d1_loc = distance from grating-adjusted L to FZP focal pt
        If (d1_loc > DoeFocLen) Then
          resid=mod(d1_loc-DoeFocLen,lambda)
          L=L-resid
        Else
	  resid=mod(DoeFocLen-d1_loc,lambda)
          L=L+resid
	End If
#endif



C  Compute refracted ray

        mu=na/nb
        Ndi=dhat(1)*Nhat(1)+dhat(2)*Nhat(2)+dhat(3)*Nhat(3)
        muNdi=mu*Ndi
        DO i=1,3
          cvec(i)=Ndi*Nhat(i)-dhat(i)
        END DO

C  Check for total internal reflection

        S2=1d0-mu*mu+muNdi*muNdi

        IF (S2.LT.0d0) THEN
          LROK=.FALSE.
          RETURN
        END IF

        cfb=DSQRT(S2)
        S1=0d0
        DO i=1,3
          rhat(i)=-mu*cvec(i)-cfb*Nhat(i)
          S1=S1+rhat(i)*rhat(i)
        END DO
        S2=DSQRT(S1)
        DO i=1,3
          rhat(i)=rhat(i)/S2
        END DO

        ! The direction of refracted ray (rhat) obeys grating equation in
        ! the case of normal incident ray, a sanity check by jzlou.
        ! The specific direction of h1HOE (along each groove surface) does
        ! not matter in determining the direction of refracted ray emerging
        ! from the transmission grating. -jzlou

#if 1
	If (GratingMaterial(crnt_Elt)==Dielectric) Then
	  ! Need to compute transimission coefficient for this ray
	  rayTC(crnt_traced_ray) = Calc_RayDielectricTC(lambda)
	Else If (GratingMaterial(crnt_Elt)==Metal) Then
	  ! Need to compute transimission coefficient for this ray
	  rayTC(crnt_traced_ray) = Calc_RayMetalTC(lambda)
	End If
#endif


C  Compute partials

        IF (ifBuild) THEN

C    Zero partial matrices

          DO 7 j=1,7
            DO 6 i=1,7
              dxidxim1(i,j)=0d0
    6       CONTINUE
    7     CONTINUE
          DO 9 j=1,6
            DO 8 i=1,7
              dxidui(i,j)=0d0
    8       CONTINUE
    9     CONTINUE
        END IF

C  Check if ray is obscured
        IF (ifObs) THEN
          LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &           psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
          IF (.NOT.LRT) THEN
            Eout(1:3)=0d0
          END IF
        END IF
        IF ((SrfType.EQ.5).AND.(IERROR.NE.0)) THEN
          LRT=.FALSE.
        ELSE
          LRT=.TRUE.
        END IF

C  Set output values for ray position, direction
        rout(1:3)=rhat(1:3)
        pout(1:3)=pr(1:3)
        NhatP(1:3)=Nhat(1:3)
!       print*,'** End of FzpTrGrating: pout(1:3) =',pout(1:3),
!    &    ', L =',L

        RETURN
        END SUBROUTINE FzpTrGrating_


C***********************************************************************
C
C DOE element modeled after Code V, with an input phase(actually OPD) profile
C that determines direction of diffracted rays and phase compensations thru
C the DOE element. -jzlou, 07/2011

	SUBROUTINE DoeTrGrating(ifBuild,ifDMBuild,ifPol,
     &   ifLNsrf,ifObs,LROK,LRT,Kc,Kr,AsphCoef,MonCoef,AnaCoef,
     &   UDSrfCoef,UDSrfType,UDSrfParam,psi,pv,prot,pin,ihat,Evec,
     &   rad,pmon,xhat,yhat,zhat,na,nb,kxa,kxb,lambda,
     &   ApType,ApVec,nObs,ObsType,ObsVec,xGrid,yGrid,zGrid,
     &   pout,rout,L,Eout,dxidxim1,dxidui,SrfType,Order,
     &   XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP)

! h1HOE(1:3,iElt) = RuleDir(1:3) = direction across rulings, -jzlou

        use elt_mod, only : DoeId,DoeWL,DoePhProf,mMonCoef,mAnaCoef
        use didesub, only: IntSrf
        use traceutil_mod

        IMPLICIT NONE

        LOGICAL LROK,ifBuild,ifDMBuild,ifPol,ifLNsrf,ifObs,!ChkRayTrans,
     &          LRT

        INTEGER i,UDSrfType,SrfType,
     &    NDP,IWORK(*),IERROR,mDP,
     &    ApType,nObs,ObsType(nObs)

        SREAL :: dxidxim1(7,7),dxidui(7,6)

        REAL*8 pv(3),prot(3),pin(3),ihat(3),psi(3),
     &   pr(3),rhat(3),L,Nhat(3),pout(3),rout(3),lambda,
     &   na,nb,kxa,kxb,AnaCoef(mAnaCoef),MonCoef(mMonCoef),
     &   rad,pmon(3),xhat(3),yhat(3),zhat(3),
     &   xGrid(3),yGrid(3),zGrid(3),ApVec(6),ObsVec(6,nObs),
     &   UDSrfCoef(*),UDSrfParam(6),thisDoeWL
        REAL*8, dimension(:) :: AsphCoef

        REAL*8 Nvec(3),Nmag,dNdp(3,3),
     &   its,
     &   dSig,strho,ihatn(3),rhatn(3),XYZSRF(3,NDP),DZSRF(2,NDP),
     &   wVec(3),DWORK(*),
     &   Kc,Kr

        REAL*8 drdi(3,3),drda(3,3),drdd(3,3),drdth(3,3),dgdi(3,3),
     &   dgda(3,3),dgdth(3,3),dgdd(3,3),dLdi(3),dLda(3),dLdth(3),dLdd(3)

        REAL*8
     &   dindi(3,3),dinda(3,3),dindth(3,3),dindd(3,3),
     &   drndi(3,3),drnda(3,3),drndth(3,3),drndd(3,3),
     &   dPhidR,
     &   dsdi(3,3),dsda(3,3),dsdth(3,3),dsdd(3,3),
     &   dGrdi(3),dGrda(3),dGrdth(3),dGrdd(3)

        REAL*8 Ndi,cvec(3),muNdi,cfb,mu,dhat(3),pdd,
     &   S1,S2,
     &   Order,shat(3),uhat(3),rts,uti,CurRuleWidth,CurRuleDir(3),
     &   prel(3),rho(3),Gr

	Real*8, save :: C1,C2,C3,C4,C5,C6

        COMPLEX*16 :: Evec(3),Eout(3)

cdebug
        COMMON/DEBUG/ shat,Gr,ihatn,rhatn,dsdi,dsda,dsdth,dsdd,
     &    dGrdi,dGrda,dGrdth,dGrdd,drdi,drda,drdd,drdth,
     &    dgdi,dgda,dgdd,dgdth,dLdi,dLda,dLdd,dLdth,
     &    dindi,dinda,dindth,dindd,drndi,drnda,drndth,drndd
cdebug

C  Compute incidence point of ray

        IF (SrfType.EQ.1) THEN
          CALL FlatSrf(ifBuild,ifLNsrf,LROK,psi,pv,prot,pin,ihat,
     &    pr,L,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.2) THEN
          CALL ConSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,
     &    prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.3) THEN
          CALL AsphSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,AsphCoef,psi,pv,
     &    prot,pin,ihat,pr,L,Nvec,Nmag,Nhat,dNdp,tAsphCoef)
        ELSE IF ((SrfType.EQ.4).OR.(SrfType.EQ.8)) THEN
          CALL MonSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,MonCoef,psi,pv,
     &    prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &    Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.5) THEN
          CALL IntSrf(ifBuild,ifLNsrf,LROK,Kc,Kr,psi,pv,prot,pin,ihat,
     &    XYZSRF,DZSRF,IERROR,IWORK,DWORK,NDP,mDP,
     &    wVec,pr,L,pmon,xhat,yhat,zhat,
     &    Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.6) THEN
          CALL AnaSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &    prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.7) THEN
          CALL UDSrf(ifBuild,ifLNsrf,LROK,UDSrfType,Kc,Kr,UDSrfCoef,
     &    UDSrfParam,psi,pv,prot,pin,ihat,pr,L,rad,pmon,xhat,yhat,zhat,
     &    Nvec,Nmag,Nhat,dNdp)
        ELSE IF (SrfType.EQ.10) THEN
          CALL ToricSrf(ifBuild,ifLNsrf,LROK,AnaCoef,psi,pv,
     &    prot,pin,ihat,pr,L,xhat,Nvec,Nmag,Nhat,dNdp)
        ELSE
          WRITE(*,*)'5: Unknown surface type. Aborting trace.'
        END IF

C  Return if ray missed surface
        IF (.NOT.LROK) RETURN

	C1=DoePhProf(1,DoeId(crnt_Elt)); C2=DoePhProf(2,DoeId(crnt_Elt))
	C3=DoePhProf(3,DoeId(crnt_Elt)); C4=DoePhProf(4,DoeId(crnt_Elt))
	C5=DoePhProf(5,DoeId(crnt_Elt)); C6=DoePhProf(6,DoeId(crnt_Elt))
	thisDoeWL=DoeWL(DoeId(crnt_Elt))

        Nsrfm=Nvec/Nmag

        S1=Nvec(1)*ihat(1)+Nvec(2)*ihat(2)+Nvec(3)*ihat(3)
        IF (S1.GT.0d0) THEN
          CALL DNEGATE(dNdp,dNdp,9)
        END IF

C  Compute grating parameters
        DO i=1,3
          prel(i)=pin(i)-pv(i)
          rho(i)=prel(i)+L*ihat(i) ! rho = (pin+L*ihat)-pv is vector from elt vertex to ray
                                   ! intersection point at current surface.
        END DO

!
!  To implement a circularily ruled surface, need to make RuleDir in the direction from
!  element vertex to ray intersection point. -jzlou
        CurRuleDir(1:3)=rho(1:3)
        CALL DUNITIZE(CurRuleDir)


! CurRuleWidth (local period) value is computed given diffraction order, dPhidR and
! lambda, where Phi(R) is input phase profile as a function of distance to DOE center.

        pdd=abs(DDOTC(rho,CurRuleDir)) ! distance from ray to DOE element center
        If (pdd<1d-18) then
          ! needed to avoid divide by zero if is chief ray
          !CurRuleWidth=lambda*1d8 ! so sin(theta_m) ~= 0
          CurRuleWidth=thisDoeWL*1d12 ! so sin(theta_m) ~= 0
        Else
	  ! dPhidR = sin(diffraction_angle), according to Code V
	  dPhidR=2*C1*pdd+4*C2*pdd**3+6*C3*pdd**5+8*C4*pdd**7+10*C5*pdd**9+
     &           12*C6*pdd**11
	  !CurRuleWidth=Order*lambda/dPhidR
	  CurRuleWidth=abs(thisDoeWL/dPhidR)
        End If

	! Use grating equation
        Gr=Order*lambda/CurRuleWidth  ! = N*lambda/a = sin(theta_N)

        CALL DPERPUVEC(CurRuleDir,CurRuleDir,psi) ! When RuleDir and psi are othogonal to
                                         ! each other, the routine doesn't change RuleDir.
        !print*,' **DoeTrGrating: aft CurRuleDir =',CurRuleDir

        CALL DXPROD(uhat,Nhat,CurRuleDir)    ! Nhat = (negative) local srf normal
                                             ! uhat will be in 'x' direction

        CALL DUNITIZE(uhat)  ! uhat in positive X direction
        CALL DXPROD(shat,uhat,Nhat)
        CALL DUNITIZE(shat)   ! shat in Y (up) direction

        !print*,'uhat =',uhat
	!print*,'shat =',shat
	!if (Order==-1) print*,'ihat =',ihat

        its=DDOTC(shat,ihat)  ! = sin(alpha)
        rts=its-Gr  ! = Gr+sin(alpha) = sin(beta), Gr can be positive or negative
		    ! depending on Order.
        !rts=-Gr-its    = Gr+sin(alpha) = sin(beta)

        uti=DDOTC(uhat,ihat)  ! ihat projected to X direction
        S1=1d0-uti*uti-rts*rts ! = 1 - sin(beta)^2 = cos(beta)^2 when uti=0d0

        IF (S1.LT.0D0) THEN
          LROK=.FALSE.
          RETURN
        END IF
        DO i=1,3
          dhat(i)=rts*shat(i)+uti*uhat(i)-SQRT(S1)*Nhat(i)
        END DO

C       strho=DDOTC(CurRuleDir,rho)
        strho=DDOTC(shat,rho)
        dSig=MOD(strho,CurRuleWidth)
        IF (dSig.LT.0d0) dSig=dSig+CurRuleWidth

	! Don't update L here; L is now updated by phase profile -jzlou
        !L=L+(rts-its)*dSig ! Ray optical path changes to make sure the cross-section of
                            ! diffracted beam has coherent phase (subtract 'extra' path
                            ! length for each ray) -jzlou

        ! Now update ray path length using the phase profile value Phi(R)
	L=L-lambda/thisDoeWL
     &    *(C1*pdd**2+C2*pdd**4+C3*pdd**6
     &      +C4*pdd**8+C5*pdd**10+C6*pdd**12)


C  Compute refracted ray
        mu=na/nb
        Ndi=dhat(1)*Nhat(1)+dhat(2)*Nhat(2)+dhat(3)*Nhat(3)
        muNdi=mu*Ndi
        DO i=1,3
          cvec(i)=Ndi*Nhat(i)-dhat(i)
        END DO

C  Check for total internal reflection

        S2=1d0-mu*mu+muNdi*muNdi

        IF (S2.LT.0d0) THEN
          LROK=.FALSE.
          RETURN
        END IF

        cfb=DSQRT(S2)
        S1=0d0
        DO i=1,3
          rhat(i)=-mu*cvec(i)-cfb*Nhat(i)
          S1=S1+rhat(i)*rhat(i)
        END DO
        S2=DSQRT(S1)
        DO i=1,3
          rhat(i)=rhat(i)/S2
        END DO


        ! The direction of refracted ray (rhat) obeys grating equation in
        ! the case of normal incident ray, a sanity check by jzlou.
        ! The specific direction of h1HOE (along each groove surface) does
        ! not matter in determining the direction of refracted ray emerging
        ! from the transmission grating. -jzlou

#if 0
        If (GratingMaterial(crnt_Elt)==Dielectric) Then
          ! Need to compute transimission coefficient for this ray
          rayTC(crnt_traced_ray) = Calc_RayDielectricTC(lambda)
        Else If (GratingMaterial(crnt_Elt)==Metal) Then
          ! Need to compute transimission coefficient for this ray
          rayTC(crnt_traced_ray) = Calc_RayMetalTC(lambda)
        End If
#endif


C  Compute partials
        IF (ifBuild) THEN
          ! Zero partial matrices
          dxidxim1(1:7,1:7)=0e0
          dxidui(1:7,1:6)=0e0
        END IF

C  Check if ray is obscured
        IF (ifObs) THEN
          LRT=ChkRayTrans(ifBuild,ifPol,ifLNsrf,
     &           psi,rho,ApType,ApVec,nObs,ObsType,ObsVec,
     &           xGrid,yGrid,zGrid,taperFactor(crnt_traced_ray))
          IF (.NOT.LRT) THEN
            Eout(1:3)=0d0
          END IF
        END IF
        IF ((SrfType.EQ.5).AND.(IERROR.NE.0)) THEN
          LRT=.FALSE.
        ELSE
          LRT=.TRUE.
        END IF

C  Set output values for ray position, direction
        rout(1:3)=rhat(1:3)
        pout(1:3)=pr(1:3)
        NhatP(1:3)=Nhat(1:3)
        RETURN
        END SUBROUTINE DoeTrGrating

	End module elemsub_mod

C***********************************************************************
C****** End file elemsub.F *********************************************
C***********************************************************************

