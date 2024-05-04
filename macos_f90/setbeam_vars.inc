C setbeam_vars.inc

	! Module variables - OPD() arguments, used to call OPD in funobj
        ! in npsol
        Logical :: ifPol_m,ifGrid_m,ifPropagate_m,ifDraw_m,ifLNegOK_m,
     &             ifBuild_m,ifDMBuild_m,ifTrace_m,ifRayPrint_m
        Logical,dimension(:),allocatable:: ifRegrid_m,ifIntSrfSet_m,L1_m
        Integer :: npts_m,iObsOpt_m,iEndElt_m,iCurWFElt_m,iCurRayElt_m,
     &             nBadRays_m
        Integer, dimension(:,:), allocatable :: DrawEltVec_m
        Real*8 :: Dicr_m,Djcr_m,dxSource_m,StartIndRef_m,StartExtinc_m,
     &            WFE_m,DAvgL_m,target_beam_m
        SREAL, dimension(:,:,:), allocatable :: Cmatrix_m,DrawRayVec_m

        Real*8 :: init_beamsize, init_aperture,final_beamsize
        Logical, save :: m_vars_set=.false., prtmsg=.false.
        Integer :: max_funobj=100 ! max iterations allowed in npsol
