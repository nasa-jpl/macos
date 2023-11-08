
#include "realtype.h"

	Module macos_vars_mod
	IMPLICIT NONE

	INCLUDE 'fits1.inc'

        CHARACTER(len=MacosCharLen):: ANS,macosHome,glassFile,Rx1,Rx2
        Integer, parameter :: max_cmd_hist=10
        CHARACTER(len=MacosCharLen) :: command,cmd_history(max_cmd_hist)

        LOGICAL :: EXIST,IsNSElt,ifReturn1,ifReturn2,ifFID

        !CHARACTER*2 CNUM
        CHARACTER*8 CInteger,XLabel,YLabel
        CHARACTER*72 CTITLE

        INTEGER i,j,k,iRay,iStartRay,iEndRay,iElt,iEm1,iSkip,
     &    iStartElt,iEndElt,nBeams,iPloTyp,saviPloTyp,iWF,mpt2,iaxis,
     &    iCoord,iPertElt,iCurRayElt,nPix,iCurWFElt,iLognRay,nAcross(6),
     &    iLinElt,iSlv1,iSlv2,iStop,iObsOpt,iSpot,nPixAdd,
     &    iComElt,nBadRays,imax,jmax,nDrawRay,iDrawRay,iLocal,jLocal,
     &    iWavelen,iStarType,StretchType,curFov,curwl,ie,
     &    nAtmPts,i1,i0,j1,j0,NN(2)

        SREAL Cout(49),Tcompact(5,7),
     &    StartTime,EndTime,DTime,RMSWFE,WFEPV,Celt(49),
     &    udia,du,dx,xdia,X0PFP,Y0PFP,dxpix,dxWF,MinL,Ptype,
     &    dxPixAdd,drx,dry,expval,r0,pTilt,atmWavelen,kval

        REAL*8 invec7(7),DDEF(9),S1,S2,StartIndRef,Tspot(7,7),
     &    XMIN,YMIN,DXO,DYO,Dicr,Djcr,dxSource,
     &    TotalAvgOpd,TotalOpdInt,DeltaOpd,zp,fp,psip(3),Vptp(3),
     &    signpsi,LSlv1,cr1dir(3),cr2dir(3),cr1pos(3),cr2pos(3),
     &    th(3),cr2indir(3),cr2inpos(3),CrossDist,CrossPt(3),
     &    StartExtinc,Power,MaxInt,xDraw(3),yDraw(3),xPt(3),yPt(3),
     &    xInt(3),xTarg(3),crPos(3),dxPt(2),sqErr,GPt(2,2),GPtInv(2,2),
     &    thx,thy,ddxPt(2),cr1indir(3),cr1inpos(3),xInt0(3),dthPt,
     &    baseFlux,baseWavelen,WFE,DAvgL,WinCenPix(2),sizPix,EltPix(2),
     &    crLocalPos(3),CentroidSpot(3),DAvgCumL,
     &    dCRCentroid(3),CRSpot(3),x,y,rAtm,tiltx,tilty,dL,u00,u01,
     &    u10,u11,TWOPI,sigx,sigy,c1,c2,c3,xdum(7),xnom(7),
     &    DTcompact(5,7), AvgL,INF,BIG

        PARAMETER (TWOPI=6.283185307179586476925286766559005768394)
        PARAMETER (BIG=1d18,INF=1d22)

        ! Scratch (dynamic) arrays
        LOGICAL, dimension(:), allocatable, save :: L1
        INTEGER, dimension(:,:), allocatable, save :: DrawEltVec
        SREAL, dimension(:,:), allocatable, save ::
     &          R1, R2, PertVec
        SREAL, dimension(:,:,:), allocatable, save ::
     &          DrawRayVec
        REAL*8, dimension(:,:), allocatable, save :: D2,
     &			        CRWalk,CRWalkProjLoc
        REAL*8, dimension(:), allocatable, save ::  DV1, DV2
        COMPLEX*8, dimension(:,:), allocatable, save ::  CD1, CD2
        COMPLEX*16, dimension(:,:), allocatable, save ::  DWF

        Integer :: npts,NoiseSeed(2)
        Real*8 :: xcent, ycent, xLocal(3), yLocal(3), zLocal(3),
     &            CntrSpot(3), RefSpot(2), pixMax,pixMin

c***loh Declares for journal and batch files
        INTEGER row_sz,col_sz
cc***loh

        Logical, save :: macos_init=.false.


! Variables in utility routines coded by jzlou
#include "util_vars.inc"

	End Module macos_vars_mod
