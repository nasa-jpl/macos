!
! Implement Plate Scale Estimate
!

#include "realtype.h"

        Module PSE_mod

        Contains

        Real*8 Function PSE
        use param_mod 
#ifdef CMACOS
        use macos_mod, only : RaySpot
#endif
        use macos_ops_mod
        
        CHARACTER(len=MacosCharLen) :: cmd,CARG(9)
        REAL*8 DARG(9),WFErms
        INTEGER IARG(9),idof,m_err_rtn
        LOGICAL LARG
        SREAL RARG(9)
        REAL*8, dimension(:,:), allocatable, save :: OPDMat,SpotMat
        SREAL, dimension(:,:), allocatable, save :: PixMat 

        Integer, save :: entry_count=0

        If (entry_count == 0) Then
          deallocate(OPDMat,PixMat,SpotMat)
          allocate(OPDMat(mpts,mpts),PIXMat(mPix,mPix),
     &             SpotMat(mRay,2),stat=m_err_rtn)
          If (m_err_rtn /= 0)
     &      call macos_memory_failure
     &                ('PSE: allocate failed!')
          entry_count = entry_count + 1
        End If

        ! Compute nominal SPOT at nElt
        cmd='SPOT'
        IARG(1)=nElt
        CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                 OPDMat,SpotMat,WFErms,PixMat)
        ! output SpotMat is nRay x 2 array
        !SpotMat = TRANSPOSE(RaySpot(1:nRay,1:2))

        ! Apply system perturbation: whole system rotates by 1urad as a group
        ! relative to the source 
        cmd='GPERTURB'
        CARG='GLOBAL'
        IARG(1)= 1;  ! group element id
        CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                 OPDMat,SpotMat,WFErms,PixMat)

        ! Compute SPOT of perturbd system at nElt
        cmd='SPOT'
        IARG(1)=nElt
        CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                 OPDMat,SpotMat,WFErms,PixMat)
        ! output SpotMat is nRay x 2 array

        PSE = 0d0  ! placeholder for now
        End Function PSE

        End Module PSE_mod

