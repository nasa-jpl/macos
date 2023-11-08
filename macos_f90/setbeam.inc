C
C setbeam.inc
C
	! Set beam size at a given surface by adjusting the
        ! source aperture
        ! 1) Define funobj, which computes the 'distance'
        !    between current beam size and target beam size and
        !    computes gradient of beam size wrt current aperture
        !    size.
        ! 2) Define simple constraints to be [0,aperture_nom], and
        !    define # of linear and nonlinear constraints to be zero.
        !
        Subroutine SetBeam(ifCurBeamSize,target_beam,
     &               npts,ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &               ifBuild,ifDMBuild,ifTrace,ifRegrid,
     &               ifIntSrfSet,ifRayPrint,
     &               iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &               Dicr,Djcr,dxSource,nBadRays,
     &               StartIndRef,StartExtinc,OPDMat,WFE,DAvgL,
     &               L1,Cmatrix,DrawRayVec,DrawEltVec)

        Implicit none

	! SetBeam arguments - same as OPD() for now
        Real*8 :: target_beam ! desired beam size at given surface
        Logical :: ifPol,ifGrid,ifPropagate,ifDraw,ifLNegOK,
     &             ifBuild,ifDMBuild,ifTrace,ifRegrid(mElt),
     &             ifIntSrfSet(mElt),ifRayPrint,L1(mRay),
     &             ifCurBeamSize
        Integer :: npts,iObsOpt,iEndElt,iCurWFElt,iCurRayElt,
     &             nBadRays,DrawEltVec(mDrawElt,mDrawRay)
        Real*8 :: Dicr,Djcr,dxSource,StartIndRef,StartExtinc,WFE,
     &            DAvgL,OPDMat(mpts,mpts)
        SREAL :: Cmatrix(7,mCm,bRay),
     &            DrawRayVec(2,mDrawElt,mDrawRay)
        Logical, save :: first_entry=.true.

        ! npsol variables
        Real*8, parameter :: inf=1d20
        Integer, parameter :: MAXVAR=2
        Integer, parameter :: nclin=0,ncnln=0,nrstate_vec=1,
     &           nrowa=1,nrowj=1,nrowr= MAXVAR,
     &           leniw=3*MAXVAR,lenw=20*MAXVAR
        Integer :: istate(MAXVAR+nclin+ncnln), npiwork(leniw)
        Integer :: inform,iter
        Real*8 :: anp(nrowa,MAXVAR),
     &            bl(MAXVAR+nclin+ncnln),bu(MAXVAR+nclin+ncnln),
     &            cnp(1+ncnln), cjac(nrowj,MAXVAR),
     &            clamda(MAXVAR+nclin+ncnln),grad(MAXVAR),
     &            rnp(nrowr,MAXVAR), work(lenw), rstate_vec(MAXVAR)
        Real*8 :: objf

        Real*8 :: Aperture_nom
        Integer :: m_err_setbeam

        if (first_entry .and. (.not.m_vars_set)) then
          ! Allocate module array variables, and initialize all module
          ! variables
          allocate(ifRegrid_m(mElt),ifIntSrfSet_m(mElt),L1_m(mRay),
     &             DrawEltVec_m(mDrawElt,mDrawRay),
     &             OPDMat_m(mpts,mpts),Cmatrix_m(7,mCm,bRay),
     &             DrawRayVec_m(2,mDrawElt,mDrawRay),
     &             stat=m_err_setbeam)
          if (m_err_setbeam /= 0)
     &    call macos_memory_failure
     &               ('SetBeam: allocate failed!')

          ifPol_m=ifPol; ifGrid_m=ifGrid
          ifPropagate_m=ifPropagate; ifDraw_m=ifDraw
          ifLNegOK_m=ifLNegOK; ifBuild_m=ifBuild
          ifDMBuild_m=ifDMBuild; ifTrace_m=ifTrace
          ifRegrid_m=ifRegrid; ifIntSrfSet_m=ifIntSrfSet
          ifRayPrint_m=ifRayPrint; L1_m=L1
          npts_m=npts; iObsOpt_m=iObsOpt; iEndElt_m=iEndElt
          iCurWFElt_m=iCurWFElt; iCurRayElt_m=iCurRayElt
          nBadRays_m=nBadRays; DrawEltVec_m=DrawEltVec
          Dicr_m=Dicr; Djcr_m=Djcr; dxSource_m=dxSource
          StartIndRef_m=StartIndRef; StartExtinc_m=StartExtinc
          WFE_m=WFE; DAvgL_m=DAvgL
          OPDMat_m=OPDMat; Cmatrix_m=Cmatrix
          DrawRayVec_m=DrawRayVec;target_beam_m=target_beam

          Aperture_nom=Aperture

          m_vars_set=.true.
          first_entry=.false.
        end if ! first_entry

	if (ifCurBeamSize) then
          Call GetBeamSize(iEndElt_m,target_beam)
          init_beamsize=target_beam
          return
        else
          target_beam_m=target_beam
        end if

        ! ************************************************
        ! Set up a test problem here: Rx=n18sd.in, SM = 23
        !iEndElt=23; iEndElt_m=iEndElt
        !target_beam=285.d0; target_beam_m=target_beam
        ! ************************************************

        Call npoptn('Verify level = -1')
        Call npoptn('Major print level = -1')
        Call npoptn('Linesearch Tolerance = 0.01') ! linear search step

        ! Define bl and bu arrays - bounds for all constraints
        bl(1)= 0d0; bu(1)=Aperture_nom  ! simple constraint x=[0,Aperture_nom]

        ! Initialize rstate_vec to current src aperture
        rstate_vec(1)=Aperture; init_aperture=Aperture

        Call npsol (nrstate_vec,nclin,ncnln,nrowa,nrowj,nrowr,
     &              anp,bl,bu,
     &              setbeam_funcon, setbeam_funobj,
     &              inform,iter,istate,
     &              cnp,cjac,clamda,objf,grad,rnp,rstate_vec,
     &              npiwork,leniw,work,lenw)

        if (prtmsg) then
          print*,'**** End of setbeam: inform = ', inform
        endif

	Call GetBeamSize(iEndElt_m,final_beamsize)
        print*,' '
        print*,'SetBeam command execution summary: '
        print*,'  => Set beam size at element', iEndElt_m
        print*,'  => Target beam size = ', target_beam
        print*,'  => Initial beam size = ', init_beamsize
        print*,'  => Initial source aperture = ', init_aperture
        print*,'  => Final beam size = ', final_beamsize
        print*,'  => Final source aperture = ', Aperture
        print*,' '

        Return
        End Subroutine SetBeam

C *********************************************************************** 

        ! A) Steps for computing objective function value:
        !   1) For given element, call OPD to that surface
        !   2) RayPos(1,iRay) stores ray positions on given element surface;
        !      nRay indicates total number of rays.
        !   3) Compute geometrical center of all rays, and 'radius' of ray
        !      beam (which is furthest ray from geometrical center).
        !   4) obj_fun_value=(radius-target_radius)**2
        !
        ! B) Steps for computing gradient of objective function wrt to
        !    source aperture:
        !   1) Set Aperture=Aperture+da
        !   2) Redo steps 1) to 4) in A) for new aperture, and get
        !      new_obj_fun_value.
        !   3) grad_obj_fun=(new_obj_fun_value-obj_fun_value)/da
        !
        Subroutine setbeam_funobj(mode,n,rstate_vec,objf,objgrd,nstate)
        Implicit none

        ! Subroutine arguments
        Integer :: mode,n,nstate
        Real*8 :: objf
        Real*8 :: rstate_vec(n),objgrd(n)

        Real*8 :: da,beam_size,tmp,new_objf
        Integer,save :: entry_count=0

        entry_count=entry_count+1
        !print*,'**** In setbeam_funobj(): mode = ', mode
        !print*,'**** In setbeam_funobj(): n = ', n
        if (prtmsg) then
          print*,'**** In setbeam_funobj(): rstate_vec(1)=',
     &         rstate_vec(1)
        end if

	! da is step size in computing gradient of beam size
        ! at the given surface wrt source change. Need to
        ! distinguish betw collmated and point source.
        if (abs(zSource) .ge. 1d10) then
          ! Collmated source.
          da=1d0  ! 1 mm = min resolution for collimated aperture size
        else
          ! Assume point source
          da=1d-05 ! 10 urad
        end if

        ! Use current source aperture to trace rays
        Aperture=rstate_vec(1)
        Call GetBeamSize(iEndElt_m,beam_size)

        objf=(beam_size-target_beam_m)**2

        if (prtmsg) then
          print*,'**** In setbeam_funobj: beam_size = ',beam_size
        end if

        ! Compute funobj gradient at current source aperture
        Aperture=Aperture+da
        Call GetBeamSize(iEndElt_m,beam_size)

        if (prtmsg) then
          print*,'**** Aperture+da = ', Aperture
          print*,'**** At aperture+da, beam size = ', beam_size
        end if

        new_objf=(beam_size-target_beam_m)**2
        objgrd=(new_objf-objf)/da
        if (prtmsg) then
          print*,'**** objgrd = ', objgrd; print*,' '
        end if

	! Restore aperture value to current state
        Aperture=Aperture-da

        ! npsol stops when mode < 0
        if (entry_count>max_funobj) mode=-1
        Return
        End Subroutine setbeam_funobj

C *********************************************************************** 

        Subroutine GetBeamSize(iEndEltArg,bm_size)
        use tracesub_mod
        Implicit none

        Integer :: iEndEltArg,iRay
        Real*8 :: bm_size,g_ctr(3),dist(3),tmp,WFEPV_

        CALL OPD(npts_m,ifPol_m,ifGrid_m,ifPropagate_m,ifDraw_m,
     &    ifLNegOK_m,ifBuild_m,ifDMBuild_m,ifTrace_m,ifRegrid_m,
     &    ifIntSrfSet_m,ifRayPrint_m,iObsOpt_m,iEndEltArg,
     &    iCurWFElt_m,iCurRayElt_m,Dicr_m,Djcr_m,dxSource_m,
     &    nBadRays_m,StartIndRef_m,StartExtinc_m,OPDMat_m,WFE_m,
     &    WFEPV_,DAvgL_m,L1_m,Cmatrix_m,DrawRayVec_m,DrawEltVec_m)

        ! Reset, so next OPD() will trace ray from source
        ifTrace_m=.false.; iCurRayElt_m=0

        ! Compute beam size at iEndEltArg
        g_ctr(1:3)=0d0
        Do iRay=1,nRay
          g_ctr(1:3)=g_ctr(1:3)+RayPos(1:3,iRay)
        End Do
        g_ctr=g_ctr/nRay; bm_size=0d0
        Do iRay=1,nRay
          dist(1:3)=RayPos(1:3,iRay)-g_ctr(1:3)
          tmp=dist(1)**2+dist(2)**2+dist(3)**2
          if (tmp>bm_size) bm_size=tmp
        End Do
        bm_size=sqrt(tmp)
        End Subroutine GetBeamSize

C *********************************************************************** 

#if 0
	Subroutine GetBeamSizeCmd(bmSize,npts_,
     &    ifPol_,ifGrid_,ifPropagate_,ifDraw_,
     &    ifLNegOK_,ifBuild_,ifDMBuild_,ifTrace_,ifRegrid_,
     &    ifIntSrfSet_,ifRayPrint_,iObsOpt_,iElt_,
     &    iCurWFElt_,iCurRayElt_,Dicr_,Djcr_,dxSource_,
     &    nBadRays_,StartIndRef_,StartExtinc_,OPDMat_,WFE_,
     &    DAvgL_,L1_,Cmatrix_,DrawRayVec_,DrawEltVec_)

        use tracesub_mod
        Implicit none

        Logical :: ifPol_,ifGrid_,ifPropagate_,ifDraw_,ifLNegOK_,
     &             ifBuild_,ifDMBuild_,ifTrace_,ifRayPrint_
        Logical :: ifRegrid_(mElt),ifIntSrfSet_(mElt),L1_(mRay)
        Integer :: npts_,iObsOpt_,iEndElt_,iCurWFElt_,iCurRayElt_,
     &             iElt_,nBadRays_
        Real*8 :: Dicr_,Djcr_,dxSource_,StartIndRef_,StartExtinc_,
     &            WFE_,DAvgL_,WFEPV_
        Real*8 :: OPDMat_(mpts,mpts)
        Integer :: DrawEltVec_(mDrawElt,mDrawRay)
        SREAL :: Cmatrix_(7,mCm,bRay),
     &           DrawRayVec_(2,mDrawElt,mDrawRay)

        Integer :: iRay,jRay
        Real*8 :: bmSize,g_ctr(3),dist(3),tmp_bs

	Call OPD(npts_,ifPol_,ifGrid_,ifPropagate_,ifDraw_,
     &    ifLNegOK_,ifBuild_,ifDMBuild_,ifTrace_,ifRegrid_,
     &    ifIntSrfSet_,ifRayPrint_,iObsOpt_,iElt_,
     &    iCurWFElt_,iCurRayElt_,Dicr_,Djcr_,dxSource_,
     &    nBadRays_,StartIndRef_,StartExtinc_,OPDMat_,WFE_,
     &    WFEPV_DAvgL_,L1_,Cmatrix_,DrawRayVec_,DrawEltVec_)

        ! Reset, so next OPD() will trace ray from source
        ifTrace_=.false.; iCurRayElt_=0

        ! Compute beam size at iElt_
        g_ctr(1:3)=0d0

	jRay=0; bmSize=0d0
        Do iRay=1,nRay
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
        Do iRay=1,nRay
          If (LRayOK(iRay) .and. LRayPass(iRay)) Then
            dist(1:3)=RayPos(1:3,iRay)-g_ctr(1:3)
            tmp_bs=dist(1)**2+dist(2)**2+dist(3)**2
            if (tmp_bs>bmSize) bmSize=tmp_bs
          End If
        End Do
        bmSize=sqrt(bmSize)
        End Subroutine GetBeamSizeCmd
#endif

C *********************************************************************** 

        ! Dummy routine since no nonlinear constraints
        Subroutine setbeam_funcon
        End Subroutine setbeam_funcon

C *********************************************************************** 

        !
        ! Set chief ray position and direction on given (one or
        ! two) element surfaces, by adjusting source position (in
        ! X and Y only) and direction (in three direction cosines)
        !
        Subroutine SetChfRay

        Logical, save :: first_entry=.true.

        if (first_entry .and. (.not.m_vars_set)) then

        end if
        End Subroutine SetChfRay
