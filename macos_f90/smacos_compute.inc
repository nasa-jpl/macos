C***********************************************************************
C       Begin file smacos_compute.inc
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C SMACOS wrapper - Given DOFs of a set of elements and the total
C perturbations for them, apply the perturbations, and compute
C the corresponding WF. 
C

#ifdef NLS_OPT
	Subroutine smacos_compute_perturbed
#else
	Subroutine smacos_compute_perturbed_np
#define mVarDOF_nls mVarDOF_np
#endif	
     &                   (OPDarr,ZCArr,
     &                    SPOTarr,SpotSize,PIXarr,
     &                    WFErms,n_varElt,eltArr,dofArr,
     &                    n_optZernArr,varZernArr,
     &                    n_optAsphArr,varAsphArr,
     &                    fovArr,wl,ptbArr)
        ! n_varElt = number of elements to perturb
        ! eltArr = element ids
        ! dofArr = mask array of 8 DOFs for each element perburted,
        !          1 for being perturbed, 0 for fixed.

        use macos_ops_mod
        Save

        ! Arguments
        REAL*8, dimension(:) :: OPDarr,ZCArr,PIXarr,ptbArr
        REAL*8, dimension(:,:) :: SPOTarr,fovArr
        Integer, dimension(:) :: eltArr,dofArr,
     &                           n_optZernArr,varZernArr
        Integer, dimension(:) :: n_optAsphArr,varAsphArr
        REAL*8 :: SpotSize,wl,WFErms
        Integer :: n_varElt

        ! Local dynamic arrays, arguments to SMACOS
        SREAL, dimension(:,:), allocatable, save :: PixMat
        REAL*8, dimension(:,:), allocatable, save :: OPDMat,SpotMat

        REAL*8 :: ChfRayPosNom(3),ChfRayDirNom(3),WavelenRx,ss,
     &            xGridNom(3),yGridNom(3),zGridNom(3),ToutNom(7,7)
        REAL*8, dimension(:), allocatable, save :: KrNom, KcNom,
     &                            nObsNom, IndRefNom
        REAL*8, dimension(:,:), allocatable, save :: psiNom,
     &          vptNom,rptNom,pNom,x_Nom,y_Nom,z_Nom,ObsTypeNom,
     &          ZernCoefNom,AsphCoefNom
        REAL*8, dimension(:,:,:), allocatable, save ::
     &          		  ObsVecNom,TEltNom

	CHARACTER(len=MacosCharLen) :: cmd,CARG(9)
        REAL*8 DARG(9),tmp_v
        INTEGER IARG(9),idof,ie,je,ke,id,ia,iz,ias,i,j,k,m_err_rtn
        LOGICAL LARG
        SREAL RARG(9)
	Integer, save :: entry_count=0

#ifndef NLS_OPT
	! aspheric coeffs not supported in constrained optim
        n_optAsphArr=0  
#endif

        IF (start_curr_optim) THEN
	  entry_count=entry_count+1
          if (allocated(KrNom)) then
            deallocate(OPDMat,PixMat,SpotMat)
            deallocate(KrNom,KcNom,nObsNom,IndRefNom,psiNom,
     &                 vptNom,rptNom,pNom,x_Nom,y_Nom,z_Nom,
     &                 ObsTypeNom,ObsVecNom,ZernCoefNom,
     &                 AsphCoefNom,TEltNom,
     &                 stat=m_err_rtn)
	    if (m_err_rtn /= 0)
     &        call macos_memory_failure
     &                ('computer_perturbed: deallocate failed!')
          end if

          !allocate(OPDMat(nGridPts,nGridPts),PIXMat(mPix,mPix),
          allocate(OPDMat(mpts,mpts),PIXMat(mPix,mPix),
     &             SpotMat(mRay,2),stat=m_err_rtn)
          if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &                ('smacos_compute: allocate failed!')

          allocate(psiNom(3,mElt),vptNom(3,mElt),rptNom(3,mElt),
     &             pNom(3,mElt),x_Nom(3,mElt),y_Nom(3,mElt),
     &             z_Nom(3,mElt),KrNom(mElt),KcNom(mElt),
     &             IndRefNom(mElt),nObsNom(mElt),
     &             ObsTypeNom(1,mElt),ObsVecNom(6,mObs,mElt),
     &		   ZernCoefNom(mZern,mElt),AsphCoefNom(mAsph,mElt),
     &             TEltNom(6,6,mElt), stat=m_err_rtn)
          if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &                ('smacos_compute: allocate failed!')

	  ! Save initial optical state for current optimization
          ChfRayPosNom(1:3)=ChfRayPos(1:3)
          ChfRayDirNom(1:3)=ChfRayDir(1:3)
          xGridNom(1:3)=xGrid(1:3)
          yGridNom(1:3)=yGrid(1:3)
          zGridNom(1:3)=zGrid(1:3)
          WavelenRx=Wavelen
          !
          Do ie=1, n_varElt
            id=eltArr(ie)
	    If (id>0 .and. EltGrp(0,id)==0) Then  ! not source nor EltGrp
              KrNom(id)=KrElt(id)
              KcNom(id)=KcElt(id)
              nObsNom(id)=nObs(id)
              psiNom(1:3,id)=psiElt(1:3,id)
              vptNom(1:3,id)=vptElt(1:3,id)
              rptNom(1:3,id)=rptElt(1:3,id)
              pNom(1:3,id)=pMon(1:3,id)
              x_Nom(1:3,id)=xMon(1:3,id)
              y_Nom(1:3,id)=yMon(1:3,id)
              z_Nom(1:3,id)=zMon(1:3,id)
	      if (nObs(id)>0) then
                ObsTypeNom(1:nObs(id),id)=ObsType(1:nObs(id),id)
                ObsVecNom(1:3,1:nObs(id),id)
     &            =ObsVec(1:3,1:nObs(id),id)
	      end if
              ZernCoefNom(1:mZern,id)=ZernCoef(1:mZern,id)
	      AsphCoefNom(1:mAsph,id)=AsphCoef(1:mAsph,id)
              CALL DEQUATE(TEltNom(1,1,id),TElt(1,1,id),36)
	    End If

	    ! Support double-pass configuration
            Do je=1, nElt
              if (LnkElt(je)==id) then
                KrNom(je)=KrElt(je)
                KcNom(je)=KcElt(je)
                nObsNom(je)=nObs(je)
                psiNom(1:3,je)=psiElt(1:3,je)
                vptNom(1:3,je)=vptElt(1:3,je)
                rptNom(1:3,je)=rptElt(1:3,je)
                pNom(1:3,je)=pMon(1:3,je)
                x_Nom(1:3,je)=xMon(1:3,je)
                y_Nom(1:3,je)=yMon(1:3,je)
                z_Nom(1:3,je)=zMon(1:3,je)
 		if (nObs(je)>0) then
                  ObsTypeNom(1:nObs(je),je)=ObsType(1:nObs(je),je)
                  ObsVecNom(1:3,1:nObs(je),je)
     &              =ObsVec(1:3,1:nObs(je),je)
	        end if
                ZernCoefNom(1:mZern,je)=ZernCoef(1:mZern,je)
	        AsphCoefNom(1:mAsph,je)=AsphCoef(1:mAsph,je)
                CALL DEQUATE(TEltNom(1,1,je),TElt(1,1,je),36)
              end if
            End Do

	    ! Support sub-assembly (group of elements) optimization
            If (id>0 .and. EltGrp(0,id)>0) Then
              Do k=1,EltGrp(0,id)
                je=EltGrp(k,id)
                if (je>0) then ! exclude source
                  KrNom(je)=KrElt(je)
                  KcNom(je)=KcElt(je)
                  nObsNom(je)=nObs(je)
                  psiNom(1:3,je)=psiElt(1:3,je)
                  vptNom(1:3,je)=vptElt(1:3,je)
                  rptNom(1:3,je)=rptElt(1:3,je)
                  pNom(1:3,je)=pMon(1:3,je)
                  x_Nom(1:3,je)=xMon(1:3,je)
                  y_Nom(1:3,je)=yMon(1:3,je)
                  z_Nom(1:3,je)=zMon(1:3,je)
	          if (nObs(je)>0) then
                    ObsTypeNom(1:nObs(je),je)
     &		      =ObsType(1:nObs(je),je)
                    ObsVecNom(1:3,1:nObs(je),je)
     &		      =ObsVec(1:3,1:nObs(je),je)
	  	  end if
                  ZernCoefNom(1:mZern,je)=ZernCoef(1:mZern,je)
	          AsphCoefNom(1:mAsph,je)=AsphCoef(1:mAsph,je) 
                  CALL DEQUATE(TEltNom(1,1,je),TElt(1,1,je),36)
                end if
              End Do
            End If
          End Do
          CALL DEQUATE(ToutNom,Tout,49)
          !
        ELSE
          ! Restore initial optical state for current optimization
          ChfRayPos(1:3)=ChfRayPosNom(1:3)
          ChfRayDir(1:3)=ChfRayDirNom(1:3)
          xGrid(1:3)=xGridNom(1:3)
          yGrid(1:3)=yGridNom(1:3)
          zGrid(1:3)=zGridNom(1:3)
          Do ie=1,n_varElt
            id=eltArr(ie)
	    If (id>0 .and. EltGrp(0,id)==0) Then ! not source nor EltGrp
              KrElt(id)=KrNom(id)
              KcElt(id)=KcNom(id)
              nObs(id)=nObsNom(id)
              psiElt(1:3,id)=psiNom(1:3,id)
              vptElt(1:3,id)=vptNom(1:3,id)
              rptElt(1:3,id)=rptNom(1:3,id)
              pMon(1:3,id)=pNom(1:3,id)
              xMon(1:3,id)=x_Nom(1:3,id)
              yMon(1:3,id)=y_Nom(1:3,id)
              zMon(1:3,id)=z_Nom(1:3,id)
	      if (nObs(id)>0) then
                ObsType(1:nObs(id),id)=ObsTypeNom(1:nObs(id),id)
                ObsVec(1:3,1:nObs(id),id)
     &    	  =ObsVecNom(1:3,1:nObs(id),id)
	      end if
              ZernCoef(1:mZern,id)=ZernCoefNom(1:mZern,id)
	      AsphCoef(1:mAsph,id)=AsphCoefNom(1:mAsph,id)
              CALL DEQUATE(TElt(1,1,id),TEltNom(1,1,id),36)
	    End If

            ! Support double-pass configuration
            Do je=1,nElt
              if (LnkElt(je)==id) then
                KrElt(je)=KrNom(je)
                KcElt(je)=KcNom(je)
                nObs(je)=nObsNom(je)
                psiElt(1:3,je)=psiNom(1:3,je)
                vptElt(1:3,je)=vptNom(1:3,je)
                rptElt(1:3,je)=rptNom(1:3,je)
                pMon(1:3,je)=pNom(1:3,je)
                xMon(1:3,je)=x_Nom(1:3,je)
                yMon(1:3,je)=y_Nom(1:3,je)
                zMon(1:3,je)=z_Nom(1:3,je)
	        if (nObs(je)>0) then
	          ObsType(1:nObs(je),je)=ObsTypeNom(1:nObs(je),je)
                  ObsVec(1:3,1:nObs(je),je)
     &		    =ObsVecNom(1:3,1:nObs(je),je)
	        end if
                ZernCoef(1:mZern,je)=ZernCoefNom(1:mZern,je)
	        AsphCoef(1:mAsph,je)=AsphCoefNom(1:mAsph,je)
                CALL DEQUATE(TElt(1,1,je),TEltNom(1,1,je),36)
              end if
            End Do

	    ! Support sub-assembly optimization
            If (id>0 .and. EltGrp(0,id)>0) Then
              Do k=1,EltGrp(0,id)
                je=EltGrp(k,id)
                if (je>0) then ! exclude source
                  KrElt(je)=KrNom(je)
                  KcElt(je)=KcNom(je)
                  nObs(je)=nObsNom(je)
                  psiElt(1:3,je)=psiNom(1:3,je)
                  vptElt(1:3,je)=vptNom(1:3,je)
                  rptElt(1:3,je)=rptNom(1:3,je)
                  pMon(1:3,je)=pNom(1:3,je)
                  xMon(1:3,je)=x_Nom(1:3,je)
                  yMon(1:3,je)=y_Nom(1:3,je)
                  zMon(1:3,je)=z_Nom(1:3,je)
	   	  if (nObs(je)>0) then
                    ObsType(1:nObs(je),je)
     &			=ObsTypeNom(1:nObs(je),je)
                    ObsVec(1:3,1:nObs(je),je)
     &			=ObsVecNom(1:3,1:nObs(je),je)
		  end if
                  ZernCoef(1:mZern,je)=ZernCoefNom(1:mZern,je)
	          AsphCoef(1:mAsph,je)=AsphCoefNom(1:mAsph,je)
                  CALL DEQUATE(TElt(1,1,je),TEltNom(1,1,je),36)
                end if
              End Do
            End If
          End Do
        END IF ! if start_curr_optim
        CALL DEQUATE(Tout,ToutNom,49)

        cmd='RESET'
        CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                 OPDMat,SpotMat,WFErms,PixMat)

C Set wavelength
          Wavelen=wl

C  Set stop if need FEX
          If (ifFEX_m) Then
            cmd='STOP'
            CARG(1)='OBJECT'
            CARG(2)='YES'
            DARG(1:3)=StopVec(1:3)
	    !print*,'smacos_compute: StopVec =',StopVec(1:3)
            CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat)
          End If

C Set FOV
          ChfRayDir(1:3)=fovArr(1:3,1)
          ChfRayPos(1:3)=fovArr(1:3,2)

C Find exit pupil
C         cmd='FEX'
C         IARG(1)=nElt-1
C         CARG(1)='YES'
C         CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
C    &                   OPDMat,SpotMat,WFErms,PixMat)

C Apply perturbations in ptbArr to a selected set of DOFs
        CARG(1)='ELEMENT'
        ia=0  ! index into ptbArr
        iz=0  ! index into varZernArr
	ias=0  ! index into varAsphArr
        Do ie=1, n_varElt
          IARG(1)=eltArr(ie)  ! element ID
          DARG(1:6)=0d0
          Do idof=1,6
            if (dofArr((ie-1)*mVarDOF_nls+idof)>0) then
              ia=ia+1
              ! apply perturbation for iElt=id, dof=idof
              DARG(idof)=ptbArr(ia)
            end if
          End Do

          If (EltGrp(0,eltArr(ie))==0) Then
            ! rigid body perturb
            cmd='PERTURB'
            CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat)

            if (dofArr((ie-1)*mVarDOF_nls+7)>0) then
              ! perturb ROC
              cmd='ROC_PERTURB'
              ia=ia+1
              DARG(1)=ptbArr(ia)
              CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                       OPDMat,SpotMat,WFErms,PixMat)
            end if
            if (dofArr(ie*mVarDOF_nls)>0) then
              ! perturb Conic
              cmd='CONIC_PERTURB'
              ia=ia+1
              DARG(1)=ptbArr(ia)
              CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                       OPDMat,SpotMat,WFErms,PixMat)
            end if
          Else
            ! Here an element group is defined at this element,
            ! so optimization on this element is for adjusting this group
            ! as a single rigid-body.
            cmd='GPERTURB'

	    !print*,' EltGrp(0,ie)=', EltGrp(0,ie)
            !print*,' Should NOT get in here !!!!!'
            !stop

            CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat)
            ! Note that ROC and CC DOFs have no effect for EltGrp optimization.
          End If

          If (n_optZernArr(ie)>0) Then
            ! perturb Zernike coeffs
            cmd='ZERN_PERTURB'
            IARG(2)=n_optZernArr(ie) ! no. Zernike terms for this element
            CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat,
     &                     varZernArr(iz+1:iz+n_optZernArr(ie)),
     &                     ptbArr(ia+1:ia+n_optZernArr(ie)),
     &                     varAsphArr(1:1), ! place holder
     &                     ptbArr(1:1))  ! place holder
            iz=iz+n_optZernArr(ie)
            ia=ia+n_optZernArr(ie)

            !Do i=1,n_optZernArr(ie)
            ! iz=iz+1; ia=ia+1
            ! j=varZernArr(iz)  ! Zernike term id
            ! ZernCoef(j,ie)=ZernCoef(j,ie)+ptbArr(ia)
            !End Do
          End If

	  If (n_optAsphArr(ie)>0) Then
	    ! perturb Aspheric coeffs
	    cmd='ASPH_PERTURB'
	    IARG(2)=n_optAsphArr(ie) ! no. aspheric coeffs for this element
	    CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat,
     &                     varZernArr(1:1), ! place holder
     &                     ptbArr(1:1), ! place holder
     &                     varAsphArr(ias+1:ias+n_optAsphArr(ie)),
     &                     ptbArr(ia+1:ia+n_optZernArr(ie)))
	    ias=ias+n_optAsphArr(ie)
	    ia=ia+n_optAsphArr(ie)
	  End If
        End Do  ! n_varElt

        If (ifCalcOPD) Then

C Find exit pupil
          if (ifFEX_m) then
            cmd='FEX'
            IARG(1)=nElt-1
            CARG(1)='YES'
            CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                   OPDMat,SpotMat,WFErms,PixMat)
          end if

C Compute OPD
	  If ((opt_tgt_m==WFE_TARGET) .OR. 
     &        (opt_tgt_m==WFE_ZMODE_TARGET)) Then
            cmd='OPD'
            IARG(1)=optElt_m
            CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat)
            !
            ! Copy OPDMat to OPDarr
            k=0; OPDarr=0d0
!           DO j=1,mdttl/2
!             DO i=1,mdttl/2
	    DO j=1,nGridPts
              DO i=1,nGridPts
                k=k+1
                OPDarr(k)=OPDMat(i,j)	
              END DO
            END DO
            OPDMat_m=OPDMat(1:nGridPts,1:nGridPts)  ! saved to module level
	  End If

C Trace chief ray or whole beam when BEAM_TARGET 
	  If (opt_tgt_m==BEAM_TARGET) Then
	    if (OptBeam_m%ifOptBeamDir .OR. OptBeam_m%ifOptBeamPos) then
	      cmd='CRT'
	      IARG(1)=optElt_m
              CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                       OPDMat,SpotMat,WFErms,PixMat)
	    else if (OptBeam_m%ifOptBeamSize) then
	      ! must trace whole beam
	      cmd='GBS'
              !IARG(1)=optElt_m
              IARG(1)=OptBeam_m%beamSizeElt
              CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                       OPDMat,SpotMat,WFErms,PixMat)
	    else if (OptBeam_m%ifOptBeamRefRayDir) then
	      cmd='RefRAY'
	      IARG(1)=OptBeam_m%beamRefRayElt
	      CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                       OPDMat,SpotMat,WFErms,PixMat) 
	    end if
	  End If 

C Compute SPOT diagram
          If (opt_tgt_m==SPOT_TARGET) Then
	    cmd='SPOT'
	    IARG(1)=optElt_m 
	    CALL MACOS_OPS(cmd,CARG,DARG,IARG,LARG,RARG,
     &                     OPDMat,SpotMat,WFErms,PixMat)
	    Do j=1,nSpotRays
	      Do i=1,2
	        SPOTarr(i,j)=SpotMat(j,i) 
	      End Do
	    End Do
	    nSpotRays_m=nSpotRays
	    SpotSize=-1d0 ! "Spot radius", or maximum distance from a ray
	                  ! to the chief ray
            Do j=1,nSpotRays
	      ss=sqrt(SPOTarr(1,j)*SPOTarr(1,j)+SPOTarr(2,j)*SPOTarr(2,j))
	      if (SpotSize<ss) SpotSize=ss 
	    End Do
	    !print*,'nSpotRays  =',nSpotRays
	    !print*,'SPOT Radius =',SpotSize
	  End If

C Compute Zernike coefficients if needed
          if (opt_tgt_m==WFE_ZMODE_TARGET) then
            Call zern_modes_coef(OPDMat,wf_zern_mode_m(1:n_wf_zern_m),
     &                           ZCArr(1:n_wf_zern_m),3)
            !print*,'smacos_: ZCArr=',ZCArr(1:n_wf_zern_m)
          end if

C Get chief ray direction at beamDirElt
          if (OptBeam_m%ifOptBeamDir) then
            beamDir_m(1:3)=ChfRayDirAtElt(1:3)
          end if

C Get chief ray position at beamPosElt
          if (OptBeam_m%ifOptBeamPos) then
            beamPos_m(1:3)=ChfRayPosAtElt(1:3)
          end if

C Get beam size at beamSize at beamSizeElt
	  if (OptBeam_m%ifOptBeamSize) then
	    beamSize_m=BeamSizeAtElt
	  end if

C Get reference ray directions
	  if (OptBeam_m%ifOptBeamRefRayDir) then
	    beamRefRayDir_m(1:3,1:nRefRay)=RefRayDir(1:3,1:nRefRay)
	  end if

          if (.false.) then
            ! Test only - save 2D OPD array
            Open(Unit=11,File='tgt_wf.txt',status='new')
            Do j=1,mdttl/2
              WRITE(11,*)(OPDMat(i,j),i=1,mdttl/2)
            End Do
            close(11)
            stop
          end if
        End If

C Restore Rx wavelength
        Wavelen=WavelenRx

        start_curr_optim = .false.
#ifdef NLS_OPT
        End Subroutine smacos_compute_perturbed
#else
        End Subroutine smacos_compute_perturbed_np
#endif
	
