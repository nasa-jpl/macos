C***********************************************************************
C       File design_optim.F
C    +----------------------------------------------------------------+
C    |  Copyright (C) 1995-2008, California Institute of Technology.  |
C    |  U.S. Government Sponsorship Is Acknowledged.                  |
C    +----------------------------------------------------------------+
C***********************************************************************

C  Subroutines for performing MACOS Optical System Optimizations
C  John Z. Lou, Jet Propulsion Laboratory
C  Last updated: 08/2008
C
C  This MACOS system optimizer can optimize WFE, WF Zernike modes, and
C  optical beam pointing (position and/or direction at an optical surface)
C  in terms of (by adjusting) eight DOF of any optical element (six DOF
C  rigid-body motion, ROC and conic constant) as well as Zernike modes
C  and aspherical coefficients for element surface figures, at any designated
C  elements including the light source (rotations and shifts of source plane).
C  Also supported are
C   1) Double-pass configuration, add a 'Link' field in the return-path
C      element section in MACOS Rx.
C   2) Subassembly (group) based optimizations.
C   3) Multi-field (Field of View) and multi-wavelength optimizations.
C   4) Choice between Nonlinear (default) and Linear algorithms.

#include "realtype.h"

	Module design_optim_mod
	use param_mod
        use src_mod
        use elt_mod
        use macos_mod, only : ifStopSet
	use math_mod
	use nls
	use dopt_mod, only : WFE_TARGET,WFE_ZMODE_TARGET,SPOT_TARGET,
     &                       BEAM_TARGET,OptAlg,NonLin,isVarElt,
     &                       OptSpotSize,ifOptPrint,SvdSvCut,
     &                       use_saved_pinv,save_pinv,PinvFile
	use zern_wf_mod

	Implicit none

	Integer, parameter :: mZern=max_zern_modes, ! defined in zern_wf_mod
     &                        mAsph=mAsphCoef    ! defined in elt_mod
        Integer :: opd_size,obj_size,tot_objfun_size,
     &             n_elts,nfov_m,nwls_m,nSpotRays_m,
     &             mVarDOF_nls,mOptZern_nls,mOptAsph_nls,
     &             tgt_wf_funit,optElt_m,mBeamOff,
     &             opt_tgt_m,n_wf_zern_m
        Integer, dimension(:), allocatable :: elt_arr_m,dof_arr_m,
     &                         n_optZern_m,varZern_arr_m,
     &                         n_optAsph_m,varAsph_arr_m
        Real*8, dimension(:), allocatable :: OPDm,ZCoefm,PIXm
        Real*8, dimension(:,:), allocatable :: SPOTm
        Real*8 :: StopVec(3), fovs_m(3,2,12),wls_m(6),
     &            beamDir_m(3),beamPos_m(3),beamSize_m,
     &	          beamRefRayDir_m(3,mRefRay),
     &            SpotSizem,SpotSize_old,SpotSize_new
        Integer :: wf_zern_mode_m(max_zern_modes)
        LREAL, dimension(:), allocatable :: objfun_nom,da,aparams_nls_m
        Logical :: start_curr_optim,ifCalcOPD,ifFEX_m,lmlsq_success,
     &		   NomSensOnly,LinNomCtrl
        Type(OptBeamStr) :: OptBeam_m
        Real*8, dimension(:,:), allocatable :: OPDMat_m


	Contains

        ! This driver routine performs optimization on an optical system
        ! based on the Levenberg-Marquardt nonlinear least-square method,
        ! a unconstrained, nonlinear iterative minimization algorithm.
	! When in this routine, it is assumed that
	!  1) Rx is already loaded.
	!  2) Aperture Stop has been set if need FEX in optimization process.

        Subroutine nls_optim_dvr(start_optim,
     &                 opt_target, ! controls whether do (what type) optimization
     &                 n_varElt,elts,dof_arr,
     &                 n_wf_zern,wf_zern_mode, ! for optimizing WF Zernike modes.
     &                 n_optZern,varZern_arr,  ! for adjusting element Zernike modes,
     &                                         ! has no effect when n_optZern(:)=0
     &                 n_optAsph,varAsph_arr,  ! for adjusting element aspheric coeffs,
     &                                         ! has no effect when n_optAsph(:)=0
     &                 fovs,nfov,fov_wt,
     &                 wls,nwls,
     &                 opt_tol,tot_itrs,
     &                 old_wfe,new_wfe,rtn_flg,
     &                 tgt_wf_flag,tgt_wf_filename,tgt_wf_format,
     &                 optElt,ifFEX,OptBeam)

	Save
	Integer :: opt_target  ! WFE (WFE_TARGET), or SPOT diagram (SPOT_TARGET),
                               ! or WFE Zernike modes (WFE_ZMODE_TARGET)
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
        Integer, dimension(:) :: n_optAsph  ! no. of Asph terms of for each varElt
        Integer, dimension(:) :: varAsph_arr ! Asph term ids for each varElt
	Real*8, dimension(:,:,:) :: fovs
	Real*8 :: fov_wt(:), wls(:),old_wfe(:,:),new_wfe(:,:)
	Real*8 :: opt_tol
	Integer :: optElt  ! Element id at which WF to be computed
	Character(len=*) :: tgt_wf_filename,tgt_wf_format
	Logical :: tgt_wf_flag,ifFEX,LCMP
	Type(OptBeamStr) :: OptBeam

	! Local variables - LM algorithm related
	LREAL, parameter :: das_scale_factor=1d-05
	LREAL, parameter :: das_rel=1d-3  ! asphere step = das_rel*|coef|
        LREAL, dimension(:), allocatable :: aparams,tmpda,sig
        LREAL, dimension(:,:), allocatable :: covar,alpha
        Logical, dimension(:), allocatable :: maska
        LREAL :: alamda, chisq, chisq_arr(tot_itrs+1),
     &           wfe_rms, wfe_rms_old(12,6)
        Integer :: len_a, idof,ifov,iwl,k,tot_itrs_done,m_err_rtn
	LREAL, save :: dtt,dx,dp,dc,drc,dcc,dz,das

	! General local variables
	Integer :: iElt,itr,i,j,ios,off,status_flag
	Integer :: ias,jterm,ie_asph
	LREAL :: cval,hap
        Logical, save :: itr_done,first_entry = .true.
	Integer, save :: entry_cnt=0  ! debug only

	entry_cnt=entry_cnt+1

	If (n_varElt .le. 0) Then
	  print*,'** Error in nls_optim_dvr: n_varElt=',n_varElt
	  stop
	Else
	  lmlsq_success=.true.
	End If

	rtn_flg=0  ! default value, indicating normal return
	tgt_wf_funit=10  ! default WF file unit

	opt_tgt_m=opt_target
	If (opt_tgt_m==WFE_ZMODE_TARGET) Then
          n_wf_zern_m=n_wf_zern
	  wf_zern_mode_m(1:n_wf_zern)=wf_zern_mode(1:n_wf_zern)
	End If
	nfov_m = nfov;
        fovs_m = fovs

        ! Wavelength
	nwls_m       = nwls;
        wls_m(:)     = 0d0
        wls_m(:nwls) = wls(:nwls)

        optElt_m=optElt; ifFEX_m=ifFEX

!-->  Beam rows (chief direction / position, beam size, reference rays)
!     ride on ANY target since 2026-10-03: blur + telecentricity, or a
!     spot target + per-field centroid positions (dyson5 smile/keystone).
	mBeamOff=0
	OptBeam_m=OptBeam
	If (.TRUE.) Then
	  if (OptBeam%ifOptBeamDir)  mBeamOff=mBeamOff+3
	  if (OptBeam%ifOptBeamPos)  mBeamOff=mBeamOff+3
	  if (OptBeam%ifOptBeamSize) mBeamOff=mBeamOff+1
	  if (OptBeam%ifOptBeamRefRayDir) mBeamOff=mBeamOff+3*nRefRay
	  if (OptBeam%ifOptBeamRefRayDir) then
	    print*,' **design_optim: nRefRay =',nRefRay
	    print*,' **design_optim: OptBeam%nomBeamRefRayDir(1:3)=',
     &                   OptBeam%nomBeamRefRayDir(1:3)
	  end if
	End If

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
	start_curr_optim = start_optim

	If (first_entry) Then
	  dtt = 1e-9   ! 1 nrad, differential for tip and tilt
	  dx  = 1d-5   ! 10 nm, differential for dx and dy (if base unit = mm)
	  dp  = 2d-5   ! 20 nm, differential for piston
	  dc  = 1d-4   ! 0.1 mrad, differential for clock
	  drc = 1d-4   ! 0.1 um, differential for ROC
	  dcc = 1d-04  ! differential for conic constant
	  dz  = 2d-6   ! 2 nm, differential for Zernike coeff
	  das = 1d-10  ! for aspherical coeff
	  first_entry=.false.
	End If  ! first_entry

	!opd_size = mRay; !opd_size = mdttl*mdttl;
	!opd_size = nGridPts*nGridPts
        opd_size = mpts*mpts
	!print*,' ** mpts,nGridPts =',mpts,nGridPts

	If (opt_tgt_m==WFE_TARGET) Then
	  obj_size=opd_size
          tot_objfun_size = (obj_size+mBeamOff)*nfov*nwls
	Else if (opt_tgt_m==WFE_ZMODE_TARGET) Then
	  obj_size=n_wf_zern_m
	  tot_objfun_size=(n_wf_zern_m+mBeamOff)*nfov*nwls
	Else if (opt_tgt_m==SPOT_TARGET) Then
	  obj_size=1  ! target=Spot diagram radius
	  tot_objfun_size=(obj_size+mBeamOff)*nfov*nwls
	Else
	  ! BEAM_TARGET, etc.
	  obj_size=0
          tot_objfun_size=mBeamOff*nfov*nwls
	  !print*,'nls_: tot_objfun_size=',tot_objfun_size
	End If

	!print*,' **tot_objfun_size =',tot_objfun_size

	! Allocate local dynamic arrays
	allocate(objfun_nom(tot_objfun_size),
     &           elt_arr_m(size(elts)),dof_arr_m(size(dof_arr)),
     &		 OPDm(opd_size),ZCoefm(n_wf_zern),PIXm(mPix*mPix),
     &           SPOTm(2,mRay),OPDMat_m(mpts,mpts),
!    &           SPOTm(2,mRay),OPDMat_m(nGridPts,nGridPts),
     &           stat=m_err_rtn)
	if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('nls_optim_dvr(1): allocate failed!')
        n_elts = n_varElt; elt_arr_m = elts; dof_arr_m = dof_arr

	allocate(n_optZern_m(n_elts))
	n_optZern_m=n_optZern
	j=0
        Do iElt=1,n_elts
          j=j+n_optZern(iElt)
        End Do
        If (j>0) Then
	  allocate(varZern_arr_m(j),stat=m_err_rtn)
	  if (m_err_rtn /= 0)
     &	    call macos_memory_failure
     &             ('nls_optim_dvr(2): allocate failed!')
          varZern_arr_m=varZern_arr
	Else
	  allocate(varZern_arr_m(1),stat=m_err_rtn)
	  if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('nls_optim_dvr(3): allocate failed!')
	End If

	allocate(n_optAsph_m(n_elts))
        n_optAsph_m=n_optAsph
        j=0
        Do iElt=1,n_elts
          j=j+n_optAsph(iElt)
        End Do
        If (j>0) Then
          allocate(varAsph_arr_m(j),stat=m_err_rtn)
          if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('nls_optim_dvr(2): allocate failed!')
          varAsph_arr_m=varAsph_arr
        Else
          allocate(varAsph_arr_m(1),stat=m_err_rtn)
          if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('nls_optim_dvr(3.1): allocate failed!')
        End If

	! Find DOFs to calibrate and define da
	! len_a = number of DOFs to calibrate
	allocate(tmpda(n_varElt*(mVarDOF_nls
     &                 +mOptZern_nls+mOptAsph_nls)))
	len_a=0; k=1; ias=0
	!
	Do iElt=1, n_elts
	  Do idof=1,mVarDOF_nls
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
	    k=k+1
	  End Do
	  ! Add zern for this varElt
	  Do idof=1,n_optZern(iElt)
	    len_a=len_a+1
	    tmpda(len_a)=dz
	  End Do
	  ! Add asph for this varElt
!-->  The differential step for an aspheric coefficient scales with
!     the coefficient (das_rel = 1e-3 of it), not the fixed das=1e-10:
!     on a deck in metres an h^4 coefficient of 1e-3 probed by 1e-10
!     moves the sag by 1e-15 m at a 70 mm aperture -- round-off, a zero
!     derivative column and "gaussj: singular matrix" at the first LM
!     step (TO, dyson5 beat 4c sec. 3.5, 2026-10-01).  A ZERO coefficient
!     (an asphere started from a conic) takes the step that moves the sag
!     at the element's circular aperture radius by 1e-7 of |Kr|; with no
!     circular aperture declared the legacy absolute step remains, and a
!     line says so.  The order scaling das_scale_factor is kept for that
!     legacy branch only, applied per order of the term.
	  ie_asph=elts(iElt)
	  Do idof=1,n_optAsph(iElt)
            len_a=len_a+1
            jterm=varAsph_arr(ias+idof)
            cval=0d0
            if (ie_asph>0 .and. jterm>=1 .and. jterm<=mAsph)
     &        cval=AsphCoef(jterm,ie_asph)
            if (cval /= 0d0) then
              tmpda(len_a)=das_rel*abs(cval)
            else if (ie_asph>0 .and. ApType(ie_asph)==1 .and.
     &               ApVec(1,ie_asph)>0d0 .and.
     &               abs(KrElt(ie_asph))<1d10) then
              hap=ApVec(1,ie_asph)
              tmpda(len_a)=1d-7*abs(KrElt(ie_asph))/hap**(2*jterm+2)
            else
              ! legacy absolute step, scaled per ORDER of the term (the
              ! old code scaled by the term's position in the list)
              tmpda(len_a)=das*das_scale_factor**(jterm-1)
              WRITE(*,'(A,I2,A,I4,A,ES9.2)')
     &          ' ** CALIB: asphere term',jterm,' of element',ie_asph,
     &          ' has a zero coefficient and no circular aperture:'//
     &          ' legacy step',tmpda(len_a)
            end if
          End Do
          ias=ias+n_optAsph(iElt)
!<--
	End Do
	!
	if (allocated(aparams_nls_m)) deallocate(aparams_nls_m)
	allocate(aparams(len_a),da(len_a),maska(len_a),
     &           sig(tot_objfun_size),covar(len_a,len_a),
     &		 alpha(len_a,len_a),aparams_nls_m(len_a),
     &           stat=m_err_rtn)
          if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &             ('nls_optim_dvr(4): allocate failed!')
	da=tmpda(1:len_a)

	!print*,'**** len_a = ', len_a
	!print*,'**** da = ', da

	! Define nominal target values
	If ((opt_tgt_m/=SPOT_TARGET) .and. (.not. tgt_wf_flag)) Then
	  objfun_nom=0d0  ! target WF default to 0

	  ! beam-row targets are filled per (field, wavelength) block below
	  close(tgt_wf_funit)
 115      continue
	Else If (opt_tgt_m==SPOT_TARGET) Then
	  objfun_nom=OptSpotSize
	End If  ! target WF
!-->  Beam-row targets, one block per (wavelength, field): the chief
!     direction / position (per-field position when OptBeamPosFov rows
!     were given), beam size, reference-ray directions.  The 2026-06 code
!     filled the FIRST block only, so every other field was driven to a
!     zero vector (its comment said "single FOV only").
	If (mBeamOff>0) Then
	  off=0
	  Do iwl=1,nwls
	    Do ifov=1,nfov
	      j=off+obj_size
	      if (OptBeam_m%ifOptBeamDir) then
	        objfun_nom(j+1:j+3)=OptBeam_m%nomBeamDir; j=j+3
	      end if
	      if (OptBeam_m%ifOptBeamPos) then
	        if (OptBeam_m%nBeamPosFov>=ifov) then
	          objfun_nom(j+1:j+3)=OptBeam_m%nomBeamPosFov(1:3,ifov)
	        else
	          objfun_nom(j+1:j+3)=OptBeam_m%nomBeamPos
	        end if
	        j=j+3
	      end if
	      if (OptBeam_m%ifOptBeamSize) then
	        objfun_nom(j+1)=OptBeam_m%nomBeamSize; j=j+1
	      end if
	      if (OptBeam_m%ifOptBeamRefRayDir) then
	        Do k=1,nRefRay
	          objfun_nom(j+1:j+3)=OptBeam_m%nomBeamRefRayDir; j=j+3
	        End Do
	      end if
	      off=off+obj_size+mBeamOff
	    End Do
	  End Do
	End If
!<--

	! Define weighting for fov: sig=1/sqrt(fov_weight)
	j=1  ! sig=1d0
	Do iwl=1,nwls
	  Do ifov=1,nfov
	    if (abs(fov_wt(ifov))<1d-10) then
	      sig(j:j+(obj_size+mBeamOff)-1)=1d0/sqrt(1d-10)
	    else
	      sig(j:j+(obj_size+mBeamOff)-1)=1d0/sqrt(abs(fov_wt(ifov)))
	    end if
	    if (mBeamOff>0) sig(j+obj_size:j+obj_size+mBeamOff-1)=
     &      sig(j+obj_size:j+obj_size+mBeamOff-1)/sqrt(OptBeam_m%beamWt)
	    j=j+obj_size+mBeamOff
	  End Do
	End Do

        aparams=0d0   ! default initial values for calibrated paramters
	maska=.true.  ! everything in aparams will be calibrated
	alamda=-1.0

	If (tot_itrs.ge.0) Then
	  ifCalcOPD=.true.
	Else
	  ifCalcOPD=.false.  ! restoring pre-optimization state only
	  ifov=1; iwl=1
	  Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &           SpotSizem,PIXm,
     &           wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &           n_optZern,varZern_arr_m,
     &           n_optAsph,varAsph_arr_m,
     &           fovs(:,:,ifov),wls(iwl),aparams)
	  goto 100
	  ! tot_itrs negative means this round is for restoring the
          ! state prior to optimization, by calling smacos_compute_perturbed()
          ! with a=0.
	End If

	If (NomSensOnly) Then
          status_flag = -2
          WRITE(*,*)' Computing nominal system sensitivity ...'
          Call svdls_itr(objfun_nom,fov_wt,aparams,status_flag)
	  WRITE(*,*)' Done!'
          go to 100  ! was 110: a jump INTO the If (allocated) block
        End If

	If (ifOptPrint) Then
	  WRITE(*,*)' '
          WRITE(*,*)'Starting MACOS System Optimization ...'
	  WRITE(*,*)'Number of rays used for optimization =',opt_npts*opt_npts
          WRITE(*,*)' '
	End If

	! Compute WFE with initial state
	old_wfe=0d0
	Do iwl=1,nwls
	  Do ifov=1,nfov
            Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &             SpotSizem,PIXm,
     &             wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &		   n_optZern,varZern_arr_m,
     &             n_optAsph,varAsph_arr_m,
     &             fovs(:,:,ifov),wls(iwl),aparams)
 	    old_wfe(ifov,iwl)=wfe_rms
	    wfe_rms_old(ifov,iwl)=wfe_rms
	    SpotSize_old=SpotSizem
	  End Do
	  if (opt_tgt_m==WFE_TARGET) then
	    WRITE(*,*)'Initial RMS WFE =', old_wfe(1:nfov,iwl)
	  else if (opt_tgt_m==SPOT_TARGET) then
	    WRITE(*,*)'Initial spot size =', SpotSize_old
	  end if
	  !WRITE(*,*)' '
	End Do

	If (opt_tgt_m==WFE_ZMode_TARGET) Then
	  WRITE(*,*)'Initial Zernike modes:'
	  Call zern_modes_coef(OPDMat_m,wf_zern_mode(1:n_wf_zern),
     &                         ZCoefm,4)
	End If

	status_flag = 0 ! for svdls_itr() only

	! Apply LM (or SVD) optimization iterations
	Do itr=1,tot_itrs

	  If (OptAlg==NonLin) Then
	    Call lmlsq(objfun_nom,sig,aparams,maska,covar,
     &                 alpha,chisq,funcs_app,alamda,lmlsq_success)
	    if (.not. lmlsq_success) then
	      WRITE(*,*) ' ** Optimization (lmlsq) failed!'
	      WRITE(*,*) ' ** Possible reasons: '
	      WRITE(*,*) ' **   Incorrect setup for optimizer'
	      WRITE(*,*) ' **   DOFs for optimization are correlated'
	      WRITE(*,*) ' **   A differential step too small to move',
     &                   ' the merit (a zero derivative column)'
!-->  No STOP here: it ended the HOST process when the engine is a mex /
!     .so (TO, dyson5 beat 4c sec. 3.6) and left the optics at whatever
!     probe the last derivative column had applied.  Label 105 restores
!     the last accepted state and returns rtn_flg=1 through the normal
!     cleanup, which is what the CLI's CALIB handler and the bindings'
!     calib_run already expect (2026-10-01).
	      rtn_flg=1; go to 105
	    end if
	    chisq_arr(itr) = chisq
          Else
	    ! Linear algorithm
	    Call svdls_itr(objfun_nom,fov_wt,aparams,status_flag)
	  End If

	  If (tot_itrs>4 .and.
     &        (itr==int(0.21*tot_itrs) .or. itr==int(0.5*tot_itrs) .or.
     &         itr==int(0.71*tot_itrs))) Then
	    if (ifOptPrint) WRITE(*,*)itr,
     &                 'out of maximum ',tot_itrs,'iterations done ...'

	    ! Convergence check
            itr_done=.true.
	    Do iwl=1,nwls
	      Do ifov=1,nfov
	        Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &                 SpotSizem,PIXm,
     &                 wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &		       n_optZern,varZern_arr_m,
     &                 n_optAsph,varAsph_arr_m,
     &                 fovs(:,:,ifov),wls(iwl),aparams)
	        !WRITE(*,*)' wfe_rms(ifov,iwl)=', wfe_rms; stop
	        If (abs(wfe_rms-wfe_rms_old(ifov,iwl))>opt_tol .or.
     &              (wfe_rms>old_wfe(ifov,iwl))) Then
                  itr_done=.false. ! not done if wfe of any FOV exceeds tol
	        End If
	        wfe_rms_old(ifov,iwl)=wfe_rms
	      End Do
	    End Do
	  End If

	  if (itr_done .and. itr>(0.7*tot_itrs)) then
	      tot_itrs_done=itr; exit
	  end if
	End Do  ! LM iteration loop
	if (.not. itr_done) tot_itrs_done=tot_itrs
	!
	! Last run
	If (OptAlg==NonLin) Then
	  alamda = 0.0
	  Call lmlsq(objfun_nom,sig,aparams,maska,covar,
     &               alpha,chisq,funcs_app,alamda,lmlsq_success)
	  if (.not. lmlsq_success) then
            rtn_flg=1; go to 105
	  end if
	  chisq_arr(tot_itrs_done+1) = chisq
	Else
	  ! Linear algorithm
	  if (.not.LinNomCtrl) then
	    ! If not using nominal system sensitivity, clean up
	    ! senstitivity array storage.
	    status_flag=-1
            Call svdls_itr(objfun_nom,fov_wt,aparams,status_flag)
	  end if
	End If

	! Apply optimized state to current system and get OPD
	Do iwl=1,nwls
	  Do ifov=1,nFov
	    Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &              SpotSizem,PIXm,
     &              wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &		    n_optZern,varZern_arr_m,
     &              n_optAsph,varAsph_arr_m,
     & 	 	    fovs(:,:,ifov),wls(iwl),aparams)
	    new_wfe(ifov,iwl)=wfe_rms
	    SpotSize_new=SpotSizem
	    !WRITE(*,*)'ifov=',ifov,' Optimized WFE = ', new_wfe(ifov,iwl)
	  End Do
	End Do

	If (opt_tgt_m==WFE_ZMode_TARGET) Then
	  WRITE(*,*)' '
          WRITE(*,*)'Final Zernike modes:'
          Call zern_modes_coef(OPDMat_m,wf_zern_mode(1:n_wf_zern),
     &                         ZCoefm,4)
        End If
	aparams_nls_m=aparams ! for external output


	!WRITE(*,*)' len_a = ', len_a
	!WRITE(*,*)' final a = ',aparams(1:len_a)
	If (ifOptPrint) Then
	  WRITE(*,*)' '
 	  WRITE(*,*)'End of MACOS system optimization'
	  WRITE(*,*)'Optimization iterations =',tot_itrs_done
	  WRITE(*,*)' '
	End If
	go to 100

!-->  Failure exit: put the optics back at the last ACCEPTED parameter
!     vector (zero = the pre-optimization state when the first LM step
!     failed), as the success path above applies the final one, then
!     fall into the normal cleanup.
 105    continue
	Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &              SpotSizem,PIXm,
     &              wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &              n_optZern,varZern_arr_m,
     &              n_optAsph,varAsph_arr_m,
     &              fovs(:,:,1),wls(1),aparams)
	aparams_nls_m=aparams
	WRITE(*,*) ' ** Optimization aborted; optics restored to the',
     &             ' last accepted state.'
!<--

 100    continue

        If (allocated(objfun_nom)) Then
 110      deallocate(objfun_nom,elt_arr_m,dof_arr_m,aparams,da,maska,
     &               tmpda,OPDm,ZCoefm,PIXm,SPOTm,OPDMat_m,sig,covar,
     &               alpha,stat=m_err_rtn)
	  if (m_err_rtn /= 0) then
	    call macos_memory_failure
     &             ('nls_optim_dvr(1): deallocate failed!')
	    stop
	  end if
        End If
	If (allocated(n_optZern_m)) Then
	  deallocate(n_optZern_m,stat=m_err_rtn)
	  if (m_err_rtn /= 0) then
	    call macos_memory_failure
     &             ('nls_optim_dvr(2): deallocate failed!')
	  end if
	End If
        If (allocated(varZern_arr_m)) Then
	  deallocate(varZern_arr_m,stat=m_err_rtn)
	  if (m_err_rtn /= 0) then
	    call macos_memory_failure
     &             ('nls_optim_dvr(3): deallocate failed!')
	  end if
	End If
	If (allocated(n_optAsph_m)) Then
          deallocate(n_optAsph_m,stat=m_err_rtn)
          if (m_err_rtn /= 0) then
            call macos_memory_failure
     &             ('nls_optim_dvr(4): deallocate failed!')
          end if
        End If
	If (allocated(varAsph_arr_m)) Then
          deallocate(varAsph_arr_m,stat=m_err_rtn)
          if (m_err_rtn /= 0) then
            call macos_memory_failure
     &             ('nls_optim_dvr(5): deallocate failed!')
          end if
        End If
	End Subroutine nls_optim_dvr

C ***************************************************************************

C
C  Generate y(a) and dyda
C  Arguments yfit and dyda are automatic arrays in calling routine.
C
	Subroutine funcs_app(a,yfit,dyda)
	use traceutil_mod, only : LUseChfRayIfOK

        LREAL, dimension(:), intent(in) :: a
        LREAL, dimension(:), intent(out) :: yfit
        LREAL, dimension(:,:), intent(out) :: dyda

        Real*8, dimension(size(a)) :: a2
	Real*8 :: wfe_rms
        Real*8, dimension(max(1,mBeamOff)) :: bv
        Integer :: size_a, ifov,iwl,off,i
	Logical :: LUseChfRayIfOK_sv

! y = objfun_nom, yfit = OPD(a)
! call smacos_compute_perturbed to compute OPD and its derivatives
! w.r.t. to a.

	LUseChfRayIfOK_sv=LUseChfRayIfOK

	! Set to false so that rms WFE will be minimized.
        ! In OPD() routine, ray average will be subtracted from
        ! each OPL.
        LUseChfRayIfOK=.false.

        size_a = size(a)

        ! Compute yfit = OPD(a)
	off=1
	Do iwl=1,nwls_m
	  Do ifov=1,nfov_m
	    Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &                SpotSizem,PIXm,
     &      	      wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &                n_optZern_m,varZern_arr_m,
     &                n_optAsph_m,varAsph_arr_m,
     & 		      fovs_m(:,1:2,ifov),wls_m(iwl),a)

	    !print*,' ** funcs_app: wfe_rms =',wfe_rms

!-->  One block per (wavelength, field): the target's own rows, then the
!     beam rows (on ANY target since 2026-10-03), written by beam_rows_ so
!     the value, derivative and linear paths cannot drift apart again.
!     The 2026-06 code wrote the beam rows one row late (off+1), stepped
!     the Zernike-mode derivative by opd_size, used obj_size-1 modes there,
!     and reset the beam-only rows to off=1 for every field (so only the
!     last field was scored).  The SPOT stride fix (obj_size, 2026-10-01)
!     is kept.
	    if (opt_tgt_m==WFE_TARGET) then
	      yfit(off:off+obj_size-1)=OPDm(1:obj_size)
	    else if (opt_tgt_m==SPOT_TARGET) then
	      yfit(off:off+obj_size-1)=SpotSizem
	    else if (opt_tgt_m==WFE_ZMODE_TARGET) then
	      yfit(off:off+obj_size-1)=ZCoefm(1:obj_size)
	    end if
	    off=off+obj_size
	    if (mBeamOff>0) then
	      call beam_rows_(yfit,off-1)
	      off=off+mBeamOff
	    end if
!<--
	  End Do
	End Do


        ! Now compute dyda
        Do i=1,size_a
          a2=a; a2(i)=a2(i)+da(i)
	  off=1
	  Do iwl=1,nwls_m
	    Do ifov=1,nfov_m
	      Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &               SpotSizem,PIXm,
     &	             wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &               n_optZern_m,varZern_arr_m,
     &               n_optAsph_m,varAsph_arr_m,
     &               fovs_m(:,1:2,ifov),wls_m(iwl),a2)

	      if (opt_tgt_m==WFE_TARGET) then
	        dyda(off:off+obj_size-1,i)=
     &            (OPDm(1:obj_size)-yfit(off:off+obj_size-1))/da(i)
	      else if (opt_tgt_m==SPOT_TARGET) then
	        dyda(off:off+obj_size-1,i)=
     &            (SpotSizem-yfit(off:off+obj_size-1))/da(i)
	      else if (opt_tgt_m==WFE_ZMODE_TARGET) then
	        dyda(off:off+obj_size-1,i)=
     &            (ZCoefm(1:obj_size)-yfit(off:off+obj_size-1))/da(i)
	      end if
	      off=off+obj_size
	      if (mBeamOff>0) then
	        call beam_rows_(bv,0)
	        dyda(off:off+mBeamOff-1,i)=
     &            (bv(1:mBeamOff)-yfit(off:off+mBeamOff-1))/da(i)
	        off=off+mBeamOff
	      end if
 	      !print*,'funcs_app for dyda: wfe_rms =', wfe_rms
	    End Do
	  End Do
        End Do

	LUseChfRayIfOK=LUseChfRayIfOK_sv ! restore
        End Subroutine funcs_app

!-->  The beam rows of one (wavelength, field) block, in the order the
!     sizing (mBeamOff) and the target fill use: direction (3), position
!     (3), size (1), reference-ray directions (3 x nRefRay); v(off+1..).
	Subroutine beam_rows_(v,off)
	Real*8, dimension(:) :: v
	Integer :: off,j,k
	j=off
	if (OptBeam_m%ifOptBeamDir) then
	  v(j+1:j+3)=beamDir_m(1:3); j=j+3
	end if
	if (OptBeam_m%ifOptBeamPos) then
	  v(j+1:j+3)=beamPos_m(1:3); j=j+3
	end if
	if (OptBeam_m%ifOptBeamSize) then
	  v(j+1)=beamSize_m; j=j+1
	end if
	if (OptBeam_m%ifOptBeamRefRayDir) then
	  Do k=1,nRefRay
	    v(j+1:j+3)=beamRefRayDir_m(1:3,k); j=j+3
	  End Do
	end if
	End Subroutine beam_rows_
!<--


C ***************************************************************************

C
C  Implement a SVD-based linear least-square system optimizer
C

	! This routine performs a single SVD-pseudo inverse
	! iteration, which updates the params vector
 	!
	Subroutine svdls_itr(objfunNom,fovWt,params,statFlg)
	use lsq_mod

	LREAL, dimension(:) :: objfunNom,fovWt,params
	Integer :: statFlg

	LREAL, dimension(:), allocatable, save :: y0,y1
	LREAL, dimension(:,:), allocatable, save :: dy0da,dady0
        LREAL :: dparams(size(params)),wfe_rms
	Integer :: iwl,ifov,off,j,ICLEN,m_err_rtn
        Logical :: exist,LCMP
	Logical, save :: first_entry = .TRUE.

	! case statFlg==-2 for computing nominal sensitivity only
	! case statFlg==-1 for clean up
	! case statFlg==0 for performing one SVD iteration

	If (first_entry .and. (statFlg==0 .OR. statFlg==-2)) Then
	  first_entry=.FALSE.
	  allocate(y0(size(objfunNom)),y1(size(objfunNom)),
     &	           dy0da(size(objfunNom),size(params)),
     &             dady0(size(params),size(objfunNom)),
     &		   stat=m_err_rtn)
         if (m_err_rtn /= 0)
     &     call macos_memory_failure('svdls_itr: allocate failed!')

	  IF (.NOT. use_saved_pinv) THEN
	    ! Calc system sensitivity around initial state
	    Call funcs_app(params,y0,dy0da)

	    ! Calc pseudo inverse of dy0da
	    Call pinv(dy0da,dady0,SvdSvCut)
	    !
            If (save_pinv) Then
         If (LCMP(PinvFile(ICLEN(PinvFile)-2:ICLEN(PinvFile)),'bin',3))
     &        Then
                ! default filename is 'dady0_saved.bin' defined in dopt_mod.F
                Open(unit=12,file=PinvFile,
     &               form='unformatted',status='replace')
                Do j=1,size(dady0,2)
                  Write(12) dady0(1:size(dady0,1),j)
                End Do
              Else
                ! ascii format
                PinvFile='dady0_saved.txt'
                Open(unit=12,file=PinvFile,status='replace')
                Do j=1,size(dady0,2)
                  Write(12,*) dady0(1:size(dady0,1),j)
                End Do
              End If
              close(12)
            End If
	    if (statFlg == -2) return ! case for calc nominal sens only
	    statFlg = 1  ! this says argument 'statFlg' can't be constant
		          ! in actual argument, otherwise would crash.
	  ELSE
            INQUIRE (FILE=PinvFile,EXIST=exist)
            IF (.not.exist) THEN
              WRITE(*,*)'*** svdls_itr: saved pinv file not found, bye!'
              STOP
            END IF
            ! use saved pinv
            If (LCMP(PinvFile(ICLEN(PinvFile)-2:ICLEN(PinvFile)),
     &          'bin',3))Then
              Open(unit=12,file=PinvFile,
     &             form='unformatted',status='old')
              Do j=1,size(dady0,2)
                Read(12) dady0(1:size(dady0,1),j)
              End Do
            Else
              ! ascii format
              Open(unit=12,file=PinvFile,status='old')
              Do j=1,size(dady0,2)
                Read(12,*) dady0(1:size(dady0,1),j)
              End Do
            End If
            close(12)
	  END IF
          print*,' ** size(dy0da,1) =',size(dy0da,1)
          print*,' ** size(dy0da,2) =',size(dy0da,2)
          !
	Else If (statFlg == -1) Then
	  goto 10  ! storage clean up
	End If

	y0 = objfunNom

	! Compute current optical state
	off=1
        Do iwl=1,nwls_m
          Do ifov=1,nfov_m
            Call smacos_compute_perturbed(OPDm,ZCoefm,SPOTm,
     &                SpotSizem,PIXm,
     &                wfe_rms,n_elts,elt_arr_m,dof_arr_m,
     &                n_optZern_m,varZern_arr_m,
     &                n_optAsph_m,varAsph_arr_m,
     &                fovs_m(:,1:2,ifov),wls_m(iwl),params)

            if (opt_tgt_m==WFE_TARGET) then
              y1(off:off+obj_size-1)=OPDm(1:obj_size)
            else if (opt_tgt_m==WFE_ZMODE_TARGET) then
              y1(off:off+obj_size-1)=ZCoefm(1:obj_size)
            end if
            off=off+obj_size
            if (mBeamOff>0) then
              call beam_rows_(y1,off-1)
              off=off+mBeamOff
            end if
          End Do
        End Do

	y1=y0-y1
 	Call DMPROD(dparams,dady0,y1,size(dady0,1),size(dady0,2),1)
	params=params+dparams
	return

	! Final clean up
 10     If (allocated(y0)) Then
          deallocate(y0,y1,dy0da,dady0,stat=m_err_rtn)
	  if (m_err_rtn /= 0)
     &      call macos_memory_failure
     &                 ('svdls_itr: deallocate failed!')
	End If

	! After clean up, first_entry must be reset to true, so when next
	! time the linear optimzer is invoked, the arrays will be allocated
	! again (and sensitivity will be re-computed or loaded)
	first_entry=.TRUE.
	End Subroutine svdls_itr


C ***************************************************************************

#undef NP_OPT
#define NLS_OPT
#include "smacos_compute.inc"

C ***************************************************************************


#if 1
	Subroutine glsss_coef_estim(indRefArr,wavlenArr,wavlen,
     &                              glsType,coef)
	Implicit none
	LREAL, dimension(:) :: indRefArr,wavlenArr,coef
	LREAL :: idxfit(3),wlfit(3),wavlen,alamda,alph(3,3),sigh(3),
     &           covarh(3,3),chisq,chisq_arr(10)
	Integer :: glsType,idx(6),itr
	Logical :: maskh(3)
	Common /wlfitArr/ wlfit

	If (glsType==HartmannGlass) Then
	  if (size(indRefArr)<3 .OR. size(wavlenArr)<3) then
	    coef=0d0
	    WRITE(*,*)
     &       ' ** glsss_coef_estim: Can not estimate glass coefs'
	    return
	  end if
	  Call get_nearest_ngbs(wavlenArr,wavlen,idx(1:3))

	  wlfit(1:3)=wavlenArr(idx(1:3))
	  idxfit(1:3)=indRefArr(idx(1:3))
	  coef=0d0; maskh=.true.; alamda=-1d0

	  Do itr=1,3
	    Call lmlsq(idxfit,sigh,coef,maskh,covarh,alph,chisq,hmod,
     &                 alamda,lmlsq_success)
	    If (.not.lmlsq_success) Then
	      Write(*,*)' Glass coef fitting failed ...'
	      coef=0d0; return
	    chisq_arr(itr)=chisq
	    End If
	  End Do

	Else If (glsType==SellmeierGlass) Then
	  ! To be implemented
	End If
	End Subroutine glsss_coef_estim

C ********************************************************

	Subroutine get_nearest_ngbs(wlarr,wl,ida)
	LREAL, dimension(:) :: wlarr
	LREAL :: wl,wlarr_loc(size(wlarr)),dd,dmin
	Integer, dimension(:) :: ida
	Integer :: nwl,i,j,imin

	wlarr_loc=wlarr; nwl=size(wlarr_loc)
	Do j=1,3
	  dmin=1d22; imin=0
	  Do i=1,nwl
	    dd=abs(wl-wlarr_loc(i))
	    if (dd<dmin) then
	      dmin=dd; imin=i
	    end if
	  End Do
	  ida(j)=imin
	  wlarr_loc(imin)=1d22
	End Do
	End Subroutine get_nearest_ngbs

C ********************************************************

	Subroutine hmod(a,yfit,dyda)
        Implicit none
	LREAL, dimension(:), intent(in) :: a
        LREAL, dimension(:), intent(out) :: yfit
        LREAL, dimension(:,:), intent(out) :: dyda
	LREAL :: wla(3)
	Integer :: i
	Common /wlfitArr/ wla

	Do i=1,3
          yfit(i)=a(1)+a(2)/(a(3)-wla(i))**1.2
	  dyda(i,1)=a(2)/(a(3)-wla(i))**1.2;
          dyda(i,2)=1d0/(a(3)-wla(i))**1.2
	  dyda(i,3)=-1.2d0*a(2)/(a(3)-wla(i))**2.2
	End Do
	End Subroutine hmod

#endif
	End Module design_optim_mod
