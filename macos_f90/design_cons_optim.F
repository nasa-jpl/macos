C***********************************************************************
C       File design_cons_optim.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2007, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************
C
C  NPSOL-based constrained optical systems optimization
C  John Z. Lou, Jet Propulsion Laboratory
C  Last modified: 04/2007
C

#include "realtype.h" 

 	Module design_cons_optim_mod
 	use param_mod
        use src_mod
	use elt_mod
	use macos_mod
	use math_mod
 	use dopt_mod, only : WFE_TARGET,WFE_ZMODE_TARGET,SPOT_TARGET,
     &                       BEAM_TARGET,OptAlg,NonLin,isVarElt
	use zern_wf_mod

	Implicit none
	Save

	Integer, parameter :: MAXVAR=100, ! max # variables constrained optimizer allows 
     &            mZern=max_zern_modes, ! defined in zern_wf_mod
     &            mAsph=mAsphCoef       ! defined in elt_mod
        Integer :: opd_size, obj_size, tot_objfun_size,
     &             n_elts, nfov_m,nwls_m, nSpotRays_m,
     &             mVarDOF_np, mOptZern_np,tgt_wf_funit, optElt_m,
     &             mBeamOff, opt_tgt_m, n_wf_zern_m
        Integer, dimension(:), allocatable :: elt_arr_m, dof_arr_m,
     &                         n_optZern_m,varZern_arr_m,
     &                         n_optAsph_m,varAsph_arr_m
        Real*8, dimension(:), allocatable :: OPDm,ZCoefm,PIXm
        Real*8, dimension(:,:), allocatable :: SPOTm
        Real*8 :: StopVec(3), fovs_m(3,2,12),wls_m(6),
     &            beamDir_m(3),beamPos_m(3),beamSize_m,
     &            beamRefRayDir_m(3,mRefRay),SpotSizem
        Integer :: wf_zern_mode_m(max_zern_modes)
        Real*8, dimension(:), allocatable :: objfun_nom,da
        Logical :: start_curr_optim,ifCalcOPD,ifFEX_m,lmlsq_success
        Type(OptBeamStr) :: OptBeam_m
        Real*8, dimension(:,:), allocatable :: OPDMat_m  ! 2D array

C ----------------------------------------------------------------------
#include "setbeam_vars.inc"
C ----------------------------------------------------------------------

 	Contains

	Subroutine np_optim_dvr(start_optim,
     &                 opt_target,  ! controls whether do (what type) optimization
     &                 n_varElt,elts,dof_arr,
     &                 cons_arr,  ! cons_arr(2,:) provides constraints for each dof
     &                 n_wf_zern,wf_zern_mode, ! for optimizing WF Zernike modes
     &                 n_optZern,varZern_arr,  ! for adjusting element Zernike modes,
     &                                         ! has no effect when n_optZern(:)=0
     &                 fovs,nfov,fov_wt,
     &                 wls,nwls,
     &                 opt_tol,tot_itrs,
     &                 old_wfe,new_wfe,rtn_flg,
     &                 tgt_wf_flag,tgt_wf_filename,tgt_wf_format,
     &                 optElt,ifFEX,OptBeam)

	Save

	! Subroutine arguments
        Integer :: opt_target  ! for entire WFE (WFE_TARGET), or
                               ! for WFE Zernike modes (WFE_ZMODE_TARGET)
	Integer :: n_wf_zern
        Integer, dimension(:) :: wf_zern_mode
        Logical :: start_optim
        Integer :: n_varElt   ! number of elements to optimize
        Integer :: tot_itrs, rtn_flg, nfov, nwls
        Integer, dimension(:) :: elts  ! ids of elements to optimize
        Integer, dimension(:) :: dof_arr  ! DOF to optimize, each element
                                   ! has 8 DOFs, 0 for fixed, 1 for variable
        Integer, dimension(:) :: n_optZern  ! no. of Zern terms of for each varElt
        Integer, dimension(:) :: varZern_arr ! Zern term ids for each varElt
	Real*8, dimension(:,:) :: cons_arr 
        Real*8, dimension(:,:,:) :: fovs
        Real*8 :: fov_wt(:), wls(:), old_wfe(:,:), new_wfe(:,:)
        Real*8 :: opt_tol
        Integer :: optElt  ! Element id at which WF to be computed
        Character(len=*) :: tgt_wf_filename,tgt_wf_format
        Logical :: tgt_wf_flag,ifFEX,LCMP
        Type(OptBeamStr) :: OptBeam

	! npsol related local variables 
	Integer, parameter :: nclin=0,ncnln=0,nrowA=1,nrowj=1,nrowR=MAXVAR,
     &                        nctotl=MAXVAR+nclin+ncnln,
     &                        leniw=3*MAXVAR,lenw=20*MAXVAR 
	Integer :: istate(nctotl)
	Integer :: len_a,inform,iw(leniw),C2MM
	Real*8 :: mtxA(nrowA,MAXVAR),bl(nctotl),bu(nctotl),
     &            cnp(1+ncnln),cjac(nrowj,MAXVAR),clamda(nctotl),
     &            grad(MAXVAR),rnp(nrowR,MAXVAR),work(lenw)
        Real*8 :: objf

	! Other local variables
	Real*8, save :: dtt,dx,dp,dc,drc,dcc,dz,wfe_rms,wfe_rms_old(12,6)
	Real*8, dimension(:), allocatable :: aparams,tmpda,sig
	Integer :: i,j,k,iElt,idof,off,ios,ifov,iwl,m_err_rtn
	Logical, save :: first_entry = .true.

	! ---------  Execution starts  -----------
	If (n_varElt .le. 0) Then
          WRITE(*,*)'** Error in np_optim_dvr: n_varElt=',n_varElt
          stop
        End If

	rtn_flg=0  ! default value, indicating normal return
        tgt_wf_funit=10  ! default WF file unit

        opt_tgt_m=opt_target
        if (opt_tgt_m==WFE_ZMODE_TARGET) then
          n_wf_zern_m=n_wf_zern
          wf_zern_mode_m(1:n_wf_zern)=wf_zern_mode(1:n_wf_zern)
        end if
        nfov_m=nfov; fovs_m= fovs
        nwls_m=nwls; wls_m(1:6)=wls(1:6)
        optElt_m=optElt; ifFEX_m=ifFEX

	If (ifFEX) Then
          if (ifStopSet) then
            StopVec = StopPos
          else
            ! Must compute system stop before design optimization
            WRITE(*,*)' '
            WRITE(*,*) ' ** Must set System Stop before optimization.'
            rtn_flg=1
            return
          end if
        End If

	! At the beginning of each design optimization, the current
        ! optical state will be saved as sort of "nominal" state, and
        ! the calibrated parameter array "aparams" (which in fact is the
        ! deviation amounts from initial optical state) will be initialized
        ! to zero. The following statement makes sure the initial optical
        ! state is saved in smacos_compute_perturbed routine.
        start_curr_optim=start_optim

	If (first_entry) Then
	  C2MM=CBM*1d3  ! C2MM=conversion_factor to mm, CBM defined in elt_mod.F
          dtt = 1e-9       ! 1 nrad, differential for tip and tilt
          dx = 1d-5/C2MM   ! 10 nm, differential for dx and dy
          dp = 2d-5/C2MM   ! 20 nm, differential for piston
          dc = 1d-4        ! 0.1 mrad, differential for clock
          drc = 1d-4/C2MM  ! 0.1 um, differential for ROC
          dcc = 1d-04      ! differential for conic constant
          dz = 2d-6/C2MM   ! 2 nm, differential for Zernike coeff
          first_entry=.false.
        End If  ! first_entry

        opd_size = (mdttl/2)*(mdttl/2)

	If (opt_tgt_m==WFE_TARGET) Then
          obj_size=opd_size
          tot_objfun_size = (opd_size+mBeamOff)*nfov*nwls
        Else If (opt_tgt_m==WFE_ZMODE_TARGET) Then
          obj_size=n_wf_zern_m
          tot_objfun_size=(n_wf_zern_m+mBeamOff)*nfov*nwls
        Else
          obj_size=0
          tot_objfun_size=mBeamOff*nfov*nwls
          !print*,'np_: tot_objfun_size=',tot_objfun_size
        End If

	allocate(objfun_nom(tot_objfun_size),
     &           elt_arr_m(size(elts)),dof_arr_m(size(dof_arr)),
     &           OPDm(opd_size),ZCoefm(n_wf_zern),PIXm(mPix*mPix),
!    &           SPOTm(2,nRay),OPDMat_m(mpts/2,mpts/2),
     &           SPOTm(2,nRay),OPDMat_m(mpts,mpts),
     &           stat=m_err_rtn)
        if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('np_optim_dvr-1: allocate failed!')
        n_elts = n_varElt; elt_arr_m = elts; dof_arr_m = dof_arr

	allocate(n_optZern_m(n_elts))
        n_optZern_m=n_optZern
        j=0
        Do iElt=1,n_elts
          j=j+n_optZern(iElt)
        End Do
        if (j>0) then
          allocate(varZern_arr_m(j),stat=m_err_rtn)
          if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('nls_optim_dvr-2: allocate failed!')
        end if
        varZern_arr_m=varZern_arr

        ! Find DOFs to calibrate and define da
        ! len_a = number of DOFs to calibrate
        allocate(tmpda(n_varElt*(mVarDOF_np+mOptZern_np)))
        len_a=0; k=1
        !
        Do iElt=1, n_elts
          Do idof=1,mVarDOF_np
            if (dof_arr(k)>0) then
              len_a=len_a+1
              if (idof==1 .or. idof==2) then
                tmpda(len_a) = dtt
              else if (idof==3) then
                tmpda(len_a) = dc
              else if (idof==4 .or. idof==5) then
                tmpda(len_a) = dx
              else if (idof==6) then
                tmpda(len_a) = dp
              else if (idof==7) then
                tmpda(len_a) = drc  ! ROC
              else
                tmpda(len_a) = dcc  ! Conic
              end if
            end if

	    ! Set bl, bu from cons_arr 
	    bl(k)=cons_arr(1,k); bu(k)=cons_arr(2,k)
            k=k+1
          End Do
          ! add zern for this varElt
          Do idof=1,n_optZern(iElt)
            len_a=len_a+1
            tmpda(len_a)=dz
          End Do
        End Do
	If (len_a>MAXVAR) Then
	  WRITE(*,*)'**np_optim_dvr: Number of variables exceeds limit (100)'
	  WRITE(*,*)'  Optimizer quiting ...'
    	  rtn_flg=1
          return	  
	End If
	!
	allocate(aparams(len_a),da(len_a),sig(tot_objfun_size),
     &           stat=m_err_rtn)
        if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('np_optim_dvr-3: allocate failed!')
	da=tmpda(1:len_a)

	! Define nominal target values
        If (.not. tgt_wf_flag) Then
          objfun_nom = 0d0 ! for design optimization, target WF default to 0

          ! OptBeam currently only for single FOV and wavelength
          off=1
          if (OptBeam_m%ifOptBeamDir) then
            ! initialize nominal beam direction for beamDirElt
            objfun_nom(obj_size+off:obj_size+off+2)=OptBeam_m%nomBeamDir
            off=off+3
          end if
          if (OptBeam_m%ifOptBeamPos) then
            ! initialize nominal beam position for beamPosElt
            objfun_nom(obj_size+off:obj_size+off+2)=OptBeam_m%nomBeamPos
!           print*,'nls_: OptBeam_m%nomBeamPos(1:3)=',
!    &          OptBeam_m%nomBeamPos(1:3)
!           print*,'objfun_nom(1:3)=',objfun_nom(1:3)
          end if
        Else
	  if (nwls>1) then
            ! tgt_wf currently not supported for multiple wavelengths
            WRITE(*,*) ' User-specified target WF not supported'//
     &                 ' yet for multiple wavelengths'
            return
          end if
          ! Use user provided target WF
          If (LCMP(tgt_wf_format,'B',1)) Then
            ! binary format
            Open(Unit=tgt_wf_funit,
     &           File=tgt_wf_filename,form='unformatted',
     &           status='OLD',iostat=ios)
            if (ios /= 0) then
              WRITE(*,*)' Open binary WF file failed, filename= ',
     &                  tgt_wf_filename
              rtn_flg=1
              go to 100
            end if
            Do ifov=1,nfov
              Do j=1,mdttl/2
                READ(tgt_wf_funit)
     &            (objfun_nom((ifov-1)*opd_size+(j-1)*mdttl/2+i),
     &             i=1,mdttl/2)
              End Do
            End Do
          Else
            ! ascii format
            Open(Unit=tgt_wf_funit,
     &           File=tgt_wf_filename,form='formatted',
     &           status='OLD',iostat=ios)
            if (ios /= 0) then
              WRITE(*,*)' Open ascii WF file failed, filename= ',
     &                  tgt_wf_filename
              rtn_flg=1
              go to 100
            end if
            Do ifov=1,nfov
              Do j=1,mdttl/2
                READ(tgt_wf_funit,*)
     &            (objfun_nom((ifov-1)*opd_size+(j-1)*mdttl/2+i),
     &             i=1,mdttl/2)
              End Do
            End Do
          End If
	  close(tgt_wf_funit)
        End If  ! target WF

	! Define weighting for fov: sig=1/sqrt(fov_weight)
        j=1  ! sig=1d0
        Do iwl=1,nwls
          Do ifov=1,nfov
            if (abs(fov_wt(ifov))<1d-10) then
              sig(j:j+(obj_size+mBeamOff)-1)=1d0/sqrt(1d-10)
            else
              sig(j:j+(obj_size+mBeamOff)-1)=1d0/sqrt(abs(fov_wt(ifov)))
            end if
            j=j+obj_size+mBeamOff
          End Do
        End Do

	aparams=0d0   ! default initial values for calibrated paramtersa

	If (tot_itrs.ge.0) then
          ifCalcOPD=.true.
        Else
          ifCalcOPD=.false.  ! restoring pre-optimization state only
          ifov=1; iwl=1
          Call smacos_compute_perturbed_np(OPDm,ZCoefm,SPOTm,
     &      SpotSizem,PIXm,
     &      wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &      n_optZern,varZern_arr_m,n_optAsph_m,varAsph_arr_m,
     &      fovs(:,:,ifov),wls(iwl),aparams)
          goto 100
          ! tot_itrs negative means this round is for restoring the
          ! state prior to optimization, by calling smacos_compute_perturbed()
          ! with a=0.
        End If

        WRITE(*,*)' '
        WRITE(*,*)'Starting MACOS Constrained System Optimization  ...'
        WRITE(*,*)
     &    'Number of rays used for optimization =',opt_npts*opt_npts
        WRITE(*,*)' '

	! Compute WFE with initial state
        old_wfe=0d0
        Do iwl=1,nwls
          Do ifov=1,nfov
            Call smacos_compute_perturbed_np(OPDm,ZCoefm,SPOTm,
     &        SpotSizem,PIXm,
     &        wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &        n_optZern,varZern_arr_m,n_optAsph_m,varAsph_arr_m,
     &        fovs(:,:,ifov),wls(iwl),aparams)
            old_wfe(ifov,iwl)=wfe_rms
            wfe_rms_old(ifov,iwl)=wfe_rms
          End Do
          WRITE(*,*)'Initial RMS WFE =', old_wfe(1:nfov,iwl)
          WRITE(*,*)' '
        End Do

	! Debugging only
	if (.false.) then
  	  print*,'len_a =',len_a
	  print*,'bl =',bl(1:len_a)
	  print*,'bu =',bu(1:len_a)
	  print*,'leniw, lenw =',leniw,lenw
	  stop
	end if

        ! --------------------- Call npsol ----------------------------
	! Set option parameters for npsol
	Call npoptn('Verify level = -1')
        Call npoptn('Major print level = -1')
        Call npoptn('Linesearch Tolerance = 0.01') ! linear search step

	Call npsol(len_a,nclin,ncnln,nrowA,nrowj,nrowR,
     &             mtxA,bl,bu,
     &	           npsol_funcon,npsol_funcobj,
     &	           inform,tot_itrs,istate,cnp,cjac,clamda(nctotl),
     & 	           objf,grad,rnp,aparams,
     &             iw,leniw,work(lenw),lenw) 
	WRITE(*,*)' '
        ! -------------------------------------------------------------

	! Apply optimized state to current system and get OPD
        Do iwl=1,nwls
          Do ifov=1,nFov
            Call smacos_compute_perturbed_np(OPDm,ZCoefm,SPOTm,
     &              SpotSizem,PIXm,
     &              wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &              n_optZern,varZern_arr_m,n_optAsph_m,varAsph_arr_m,
     &              fovs(:,:,ifov),wls(iwl),aparams)
            new_wfe(ifov,iwl)=wfe_rms
            !WRITE(*,*)'ifov=',ifov,' Optimized WFE = ', new_wfe(ifov,iwl)
          End Do
        End Do

	If (opt_tgt_m==WFE_ZMode_TARGET) Then
          WRITE(*,*)' '
          WRITE(*,*)'Final Zernike modes:'
          Call zern_modes_coef(OPDMat_m,wf_zern_mode(1:n_wf_zern),
     &                         ZCoefm,4)
        End If

 100	continue

	deallocate(objfun_nom,elt_arr_m,dof_arr_m,aparams,da,tmpda,
     &             OPDm,ZCoefm,PIXm,SPOTm,OPDMat_m)
        if (allocated(n_optZern_m)) deallocate(n_optZern_m)
        if (allocated(varZern_arr_m)) deallocate(varZern_arr_m)
	End Subroutine np_optim_dvr

C *********************************************************************************

	! Calc objective function and gradient vector for npsol
	Subroutine npsol_funcobj(mode,n,a,obf,grd,nstat)
	Integer :: mode,n,nstat
	Real*8 :: a(n),obf,grd(n) ! n=size(a)

! obf = objective function
!       For WFE_TARGET, its value is the 2-norm of difference between
!       current WF and target WF. 
!       For WFE_ZMODE_TARGET, its value is the 2-norm of Zernike coefs  
!       being optimized. 
!       For BEAM_TARGET, its value is the 2-norm of difference between
!       current beam dir/pos and target beam dir/pos values.
!       When optimizing on both WF and BEAM, obf is combination of each target obf.
!
! grd = gradient vector of objective function 
!       Its dimension is n by 1, where n is size(a), a = vector of calibrated 
!       parameters. The i-th element of grd = d(obf)/d(a_i), computed with
!       finite difference.

	! Local variables
	Real*8, dimension(size(a)) :: a2
        Real*8 :: wfe_rms,yfit(size(objfun_nom)),obfp
        Integer :: size_a,ifov,iwl,off,off2,i,j,k

	size_a = size(a)

	! Compute obf - objective function
	obf=0d0; off=1
	Do iwl=1,nwls_m
          Do ifov=1,nfov_m
            Call smacos_compute_perturbed_np(OPDm,ZCoefm,SPOTm,
     &             SpotSizem,PIXm,
     &             wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &             n_optZern_m,varZern_arr_m,n_optAsph_m,varAsph_arr_m,
     &             fovs_m(:,1:2,ifov),wls_m(iwl),a)

	    If (opt_tgt_m==WFE_TARGET) Then
              yfit(off:off+opd_size-1)=OPDm(1:opd_size)
	      if (.not. OptBeam_m%ifOptBeamDir) then
                off=off+opd_size
              else
                off=off+opd_size; off2=1
                if (OptBeam_m%ifOptBeamDir) then
                  yfit(off+off2:off+off2+2)=beamDir_m(1:3)
                  off2=off2+3
                end if
                if (OptBeam_m%ifOptBeamPos) then
                  yfit(off+off2:off+off2+2)=beamPos_m(1:3)
                end if
                off=off+mBeamOff
              end if
	    Else If (opt_tgt_m==WFE_ZMODE_TARGET) Then
              yfit(off:off+n_wf_zern_m-1)=ZCoefm(1:n_wf_zern_m)
              if (.not. OptBeam_m%ifOptBeamDir) then
                off=off+obj_size
              else
                off=off+obj_size; off2=1
                if (OptBeam_m%ifOptBeamDir) then
                  yfit(off+off2:off+off2+2)=beamDir_m(1:3)
                  off2=off2+3
                end if
                if (OptBeam_m%ifOptBeamPos) then
                  yfit(off+off2:off+off2+2)=beamPos_m(1:3)
                end if
                off=off+mBeamOff
              end if
	    Else
	      off=1
              if (OptBeam_m%ifOptBeamDir) then
                yfit(off:off+2)=beamDir_m(1:3)
                off=off+3
              end if
              if (OptBeam_m%ifOptBeamPos) then
                yfit(off:off+2)=beamPos_m(1:3)
                off=off+3
              end if
	    End If
	  End Do
	End Do
	!
	! Calc scalar obf
        yfit=yfit-objfun_nom
        obf=DDOTN(yfit,yfit,size(yfit))

	! Now compute grd ( = d(obf)/da )
        Do i=1, size_a
	  a2=a; a2(i)=a2(i)+da(i)
          off=1
          Do iwl=1,nwls_m
            Do ifov=1,nfov_m
	      Call smacos_compute_perturbed_np(OPDm,ZCoefm,SPOTm,
     &          SpotSizem,PIXm,
     &          wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &          n_optZern_m,varZern_arr_m,n_optAsph_m,varAsph_arr_m,
     &          fovs_m(:,1:2,ifov),wls_m(iwl),a2)

	      If (opt_tgt_m==WFE_TARGET) Then
                yfit(off:off+opd_size-1)=OPDm(1:opd_size)
                if (.not. OptBeam_m%ifOptBeamDir) then
                  off=off+opd_size
                else
                  off=off+opd_size; off2=1
                  if (OptBeam_m%ifOptBeamDir) then
                    yfit(off+off2:off+off2+2)=beamDir_m(1:3)
                    off2=off2+3
                  end if
                  if (OptBeam_m%ifOptBeamPos) then
                    yfit(off+off2:off+off2+2)=beamPos_m(1:3)
                  end if
                  off=off+mBeamOff
                end if
	      Else If (opt_tgt_m==WFE_ZMODE_TARGET) Then
	        yfit(off:off+obj_size-1)=ZCoefm(1:obj_size-1)
                if (.not. OptBeam_m%ifOptBeamDir) then
                  off=off+obj_size
                else
                  off=off+opd_size; off2=1
                  if (OptBeam_m%ifOptBeamDir) then
                    yfit(off+off2:off+off2+2)=beamDir_m(1:3)
                    off2=off2+3
                  end if
                  if (OptBeam_m%ifOptBeamPos) then
                    yfit(off+off2:off+off2+2)=beamPos_m(1:3)
                  end if
                  off=off+mBeamOff
                end if
	      Else
		off=1
                if (OptBeam_m%ifOptBeamDir) then
                  yfit(off:off+2)=beamDir_m(1:3)
                  off=off+3
                end if
                if (OptBeam_m%ifOptBeamPos) then
                  yfit(off:off+2)=beamPos_m(1:3)
                  off=off+3
                end if
	      End If
	    End Do
	  End Do

	  ! Calc grd
          yfit=yfit-objfun_nom
          obfp=DDOTN(yfit,yfit,size(yfit)) ! perturbed obf
          grd(i)=(obfp-obf)/da(i)
	End Do  ! do size_a 
	End Subroutine npsol_funcobj


	! Function computing nonlinear constraints for npsol, a dummy
	! function since no nonlinear constraints in our optimization 
	Subroutine npsol_funcon
        End Subroutine npsol_funcon


C *********************************************************************************

#undef NLS_OPT
#define NP_OPT
#include "smacos_compute.inc"

C *********************************************************************************
#include "setbeam.inc"

 	End Module design_cons_optim_mod
