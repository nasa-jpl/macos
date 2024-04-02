C***********************************************************************
C       Begin file smacos_vars_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module smacos_vars_mod
          use Constants, only: TWOPI, BIG, INF
        IMPLICIT NONE

	INCLUDE 'fits1.inc'

	! Local variables

        Character(len=MacosCharLen) :: command,ans,glassFile
        Character*100 :: macosHome

        LOGICAL ifInit, ifEcho
        CHARACTER*1000 STACK
        INTEGER ls
        COMMON /SCCOM/ls, STACK
        COMMON /SCIO/ifEcho

C       These are other MACOS variables

        LOGICAL EXIST,IsNSElt,ifReturn1,ifReturn2

        CHARACTER*2 CNUM
        CHARACTER*8 CInteger,XLabel,YLabel
        CHARACTER*72 CTITLE

        INTEGER i,j,k,npts,iRay,iStartRay,iEndRay,iElt,iEm1,
     &   iSkip,iStartElt,iEndElt,nBeams,iPloTyp,iWF,mpt2,iaxis,
     &   iCoord,iPertElt,iCurRayElt,nPix,iCurWFElt,iLognRay,nAcross(6),
     &   iLinElt,iSlv1,iSlv2,iStop,iSpot,nPixAdd,
     &   iComElt,nBadRays,imax,jmax,nDrawRay,iDrawRay,iLocal,jLocal,
     &   iWavelen,iStarType,NoiseSeed(2),StretchType,ie,
     &   nAtmPts,i1,i0,j1,j0,NN(2)

        Integer, save :: iObsOpt

        SREAL Cout(49),Tcompact(5,7),WFEPV,
     &   StartTime,EndTime,DTime,Celt(49),
     &   udia,du,dx,xdia,X0PFP,Y0PFP,dxpix,dxWF,MinL,Ptype,
     &   dxPixAdd,drx,dry,expval,r0,pTilt,atmWavelen,kval

        REAL*8 invec7(7),DDEF(7),S1,S2,StartIndRef,Tspot(7,7),
     &   XMIN,YMIN,DXO,DYO,xcent,ycent,Dicr,Djcr,dxSource,
     &   TotalAvgOpd,TotalOpdInt,DeltaOpd,zp,fp,psip(3),Vptp(3),
     &   signpsi,LSlv1,cr1dir(3),cr2dir(3),cr1pos(3),cr2pos(3),
     &   th(3),cr2indir(3),cr2inpos(3),CrossDist,CrossPt(3),
     &   StartExtinc,Power,MaxInt,xDraw(3),yDraw(3),xPt(3),yPt(3),
     &   xInt(3),xTarg(3),crPos(3),dxPt(2),sqErr,GPt(2,2),GPtInv(2,2),
     &   thx,thy,ddxPt(2),cr1indir(3),cr1inpos(3),xInt0(3),dthPt,
     &   baseFlux,baseWavelen,WFE,DAvgL,WinCenPix(2),sizPix,EltPix(2),
     &   POISSON,xLocal(3),yLocal(3),crLocalPos(3),CentroidSpot(3),
     &   dCRCentroid(3),CRSpot(3),x,y,rAtm,tiltx,tilty,dL,u00,u01,
     &   u10,u11,zLocal(3),sigx,sigy,c1,c2,c3,xdum(7),xnom(7),
     &   DTcompact(5,7),RefSpot(2),CntrSpot(3),AvgL, !INF,BIG,
     &   pixMax,pixMin,DAvgCumL

        ! PARAMETER (TWOPI=6.283185307179586476925286766559005768394)
        ! PARAMETER (BIG=1d18,INF=1d22)

C       User COMMON blocks can be added here...

        COMMON /UserCommon4/npts,NoiseSeed,iSpot
        COMMON /UserCommon8/xcent,ycent,xLocal,yLocal,zLocal,
     &                      CntrSpot,RefSpot

C       Scratch arrays
        LOGICAL, dimension(:), allocatable, save :: L1
        INTEGER, dimension(:,:), allocatable, save :: DrawEltVec
        SREAL, dimension(:,:), allocatable, save ::
     &          R1, R2, PertVec
        SREAL, dimension(:,:,:), allocatable, save ::
     &         DrawRayVec
        REAL*8, dimension(:,:), allocatable, save :: D2
        REAL*8, dimension(:), allocatable, save ::  DV1,DV2
        COMPLEX*8, dimension(:,:), allocatable, save ::  CD1,CD2
        COMPLEX*16, dimension(:,:), allocatable, save ::  DWF

        Logical, save :: first_enter = .true.
        Integer :: m_err_smacos

	Type srcMapStr
          integer :: size
	  Real*8, dimension(:,:), allocatable :: rayMap,segMap
	End Type
        Type(srcMapStr) :: srcMap


! Variables in utility routines coded by jzlou
#include "util_vars.inc"

	End Module smacos_vars_mod
