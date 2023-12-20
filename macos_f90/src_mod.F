C***********************************************************************
C       Begin file src_mod.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

#include "realtype.h"

	Module src_mod
        use Kinds
	use param_mod

	Implicit none

        !   integer, parameter :: SrcGrid_


        ! From "prop_defs.inc"
        ! GridTypeName(1)='Circular'
        ! GridTypeName(2)='Square'
        ! GridTypeName(3)='Hex'
        ! GridTypeName(4)='Pie'
        ! GridTypeName(5)='Flower'

        Character*32, PARAMETER :: GridTypeName(6) = ['Circular',
     &                                                'Square',
     &                                                'Hex',
     &                                                'Pie',
     &                                                'Flower',
     &                                                'UDGrid']


        Save

	Logical :: UDBeam_FLG, DrawMode
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

   	  Subroutine src_mod_init_vars()

            m_err_src    = 0
            UDSrcProf    = ''

            UDBeam_FLG   = .FALSE.
            DrawMode     = .FALSE.

            nGridPts     = 0
            GridType     = 0
            BeamType     = 0
            nTestPts     = 0
            nMinPts      = 0
            nSpotPts     = 0
            nWavelen     = 0
            nStarType    = 0
            nPetals(:)   = 0
            nRings       = 0
            StopElt      = 0
            npts_Rx      = 0
            opt_npts     = 0

            Aperture     = 0e0_pr
            Obscratn     = 0e0_pr
            xGrid(:)     = 0e0_pr
            yGrid(:)     = 0e0_pr
            zGrid(:)     = 0e0_pr
            ChfRayDir(:) = 0e0_pr
            ChfRayPos(:) = 0e0_pr
            rxBeam       = 0e0_pr
            ryBeam       = 0e0_pr
            width        = 0e0_pr
            gap          = 0e0_pr
            CosPower     = 0e0_pr
            SegXgrid(:)  = 0e0_pr
            StopPos(:)   = 0e0_pr
            srcPos(:)    = 0e0_pr
            radCtr       = 0e0_pr
            RingRads(:)  = 0e0_pr

            xGrid_sv(:)  = 0e0_pr
            yGrid_sv(:)  = 0e0_pr

            Ex0          = 0e0_prc
            Ey0          = 0e0_prc

            ObsType_FLG(:) = .FALSE.
            ObsVec_FLG(:)  = .FALSE.

            RayIndex(:, :) = 0
            RayMap(:, :)   = 0
            SegMap(:, :)   = 0

            FiltFlux(:, :) = 0e0_pr
            ds1(:, :)      = 0e0_pr
            ds2(:, :)      = 0e0_pr
            FiltWavelen(:) = 0e0_pr

          end subroutine src_mod_init_vars


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

	! Integer arrays
	allocate(RayIndex(2,mRay), RayMap(mpts,mpts),
     &		 SegMap(mpts,mpts), stat=m_err_src)
	if (m_err_src /= 0)
     &    call macos_memory_failure
     &               ('src_mod_init: allocate failed!')

	! Real*8 arrays
	allocate(FiltFlux(mWavelen,mStarType), FiltWavelen(mWavelen),
     &		 stat=m_err_src)
        if (m_err_src /= 0)
     &    call macos_memory_failure
     &               ('src_mod_init: allocate failed!')

          ! initialise variables
          CALL src_mod_init_vars()

	End Subroutine src_mod_init

	End Module src_mod
