C***********************************************************************
C       Begin file src_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module src_mod
	use param_mod

	Implicit none
	Save

	Logical :: UDBeam_FLG, DrawMode, FSEG_CCAT_Flg
	Logical, dimension(:), allocatable :: ObsType_FLG,
     &		 ObsVec_FLG

	Integer, parameter :: mRing=20  ! max number flower rings allowed
	Integer :: nGridPts,GridType,BeamType,nTestPts,nMinPts,
     &		   nSpotPts,nWavelen,nStarType,nPetals(0:mRing),StopElt,
     &	           nRings,npts_Rx,opt_npts	   

	Integer, dimension(:,:), allocatable :: RayIndex,RayMap,SegMap

	Real*8 :: Aperture,Obscratn,xGrid(3),yGrid(3),zGrid(3),
     &		  ChfRayDir(3),ChfRayPos(3),rxBeam,ryBeam,width,gap,
     &		  CosPower,SegXgrid(3),StopPos(3),srcPos(3),radCtr,
     &            xGrid_sv(3),yGrid_sv(3),RingRads(0:mRing)

	Real*8, dimension(:,:), allocatable :: FiltFlux,ds1,ds2
	Real*8, dimension(:), allocatable :: FiltWavelen

	COMPLEX*16 :: Ex0,Ey0

	Character(len=MacosValLen) :: UDSrcProf

	Integer :: m_err_src


	Contains

	Subroutine src_mod_init

	DrawMode=.FALSE.  ! default

	! Deallocate if arrays already allocated	
	if (allocated(RayIndex)) then
	  deallocate(ObsType_FLG, ObsVec_FLG,
     &		     RayIndex, RayMap, SegMap,
     &		     FiltFlux, FiltWavelen,
     &	             stat=m_err_src)
	  if (m_err_src /= 0)
     &      call macos_memory_failure
     &               ('src_mod_init: deallocate failed!') 
	end if

	! Logical arrays
	allocate(ObsType_FLG(mObs), ObsVec_FLG(mObs),
     &		 stat=m_err_src)
        if (m_err_src /= 0)
     &    call macos_memory_failure
     &               ('src_mod_init: allocate failed!')
	ObsType_FLG = .false.; ObsVec_FLG = .false.

	! Integer arrays
	allocate(RayIndex(2,mRay), RayMap(mpts,mpts),
     &		 SegMap(mpts,mpts), stat=m_err_src)
	if (m_err_src /= 0)
     &    call macos_memory_failure
     &               ('src_mod_init: allocate failed!')
	RayIndex = 0; RayMap = 0; SegMap = 0

	! Real*8 arrays
	allocate(FiltFlux(mWavelen,mStarType), FiltWavelen(mWavelen),
     &		 stat=m_err_src)
        if (m_err_src /= 0) 
     &    call macos_memory_failure
     &               ('src_mod_init: allocate failed!')
	FiltFlux = dzero_; FiltWavelen = dzero_

	End Subroutine src_mod_init

	End Module src_mod
