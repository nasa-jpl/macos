C  sysprop_mod.F -- first-order optical system properties (SYSPROP).
C
C  Packages the system's first-order / paraxial + diffraction properties
C  into one array, computed from existing module state AFTER the EFL
C  command has set focalRatio (efl = focalRatio*Aperture is reconstructed
C  here -- the command-local 'efl' is not module scope).  Shared by the
C  interactive/SMACOS 'SYSPROP' command and macos_api_mod's sysprop()
C  struct getter so the two never drift.
C
C  Inputs are all module variables: focalRatio,Wavelen,CWM,CBM,dxElt
C  (elt_mod), Aperture,Obscratn (src_mod), nElt (elt_mod), mdttl
C  (param_mod).  No new wrapper / no efl.inc refactor needed.

      MODULE sysprop_mod
        USE param_mod, ONLY: mdttl
        USE elt_mod,   ONLY: nElt, focalRatio, Wavelen, CWM, CBM, dxElt
        USE src_mod,   ONLY: Aperture, Obscratn
        IMPLICIT NONE

        INTEGER, PARAMETER :: nSysProp = 14
        INTEGER, PARAMETER :: SP_EFL_BU=1, SP_FNO=2, SP_DPUP_M=3,
     &    SP_DPUP_BU=4, SP_OBSC=5, SP_LAMBDA_M=6, SP_LAMBDA_WU=7,
     &    SP_LAMD_RAD=8, SP_LAMD_ARCSEC=9, SP_LAMD_PX=10,
     &    SP_PLATE_ARCSEC_PX=11, SP_PLATE_PX_RAD=12, SP_NYQUIST_BU=13,
     &    SP_DX_FOCAL_BU=14

        REAL*8  :: sp_vals(nSysProp)
        LOGICAL :: sp_valid = .FALSE.   ! angular quantities set (aper+wvl)
        LOGICAL :: sp_hasefl = .FALSE.  ! efl/F#/plate set (EFL succeeded)
        LOGICAL :: sp_haspx = .FALSE.   ! px quantities set (post-INT)

        REAL*8, PARAMETER :: ARCSEC_PER_RAD = 206264.806247096d0

      CONTAINS

C  Compute the first-order properties at focal-plane element iElt.
C  Requires focalRatio already set (run the EFL command first).
C  The pixel-based quantities (lambda/D in px, plate scales) need a
C  prior propagation to iElt to have set dxElt(iElt); they are left 0
C  otherwise.  Returns ok=.TRUE. iff the core (focal-length-based)
C  quantities were computed.
        SUBROUTINE DeriveSysProps(iElt, ok)
          IMPLICIT NONE
          INTEGER, INTENT(IN)            :: iElt
          LOGICAL, INTENT(OUT), OPTIONAL :: ok
          REAL*8 :: dpup_bu, efl_bu, dpup_m, lam_m, lam_bu, dxf_bu

          IF (PRESENT(ok)) ok = .FALSE.
          sp_valid  = .FALSE.
          sp_hasefl = .FALSE.
          sp_haspx  = .FALSE.
          sp_vals(:) = 0d0

          IF ((iElt < 1) .OR. (iElt > nElt)) RETURN
          IF ((ABS(Aperture) <= 0d0) .OR. (CBM <= 0d0)) RETURN

          dpup_bu = ABS(Aperture)
          dpup_m  = dpup_bu * CBM
          lam_m   = Wavelen * CWM
          lam_bu  = lam_m / CBM                 ! lambda in BaseUnits
          IF (dpup_m <= 0d0) RETURN

C  Angular / aperture quantities: need only the aperture + wavelength
C  (NO EFL marginal-ray trace), so always available -- this is what the
C  planet-offset tilt (lamD_rad) and most design targets use.
          sp_vals(SP_DPUP_M)      = dpup_m
          sp_vals(SP_DPUP_BU)     = dpup_bu
          sp_vals(SP_OBSC)        = Obscratn
          sp_vals(SP_LAMBDA_M)    = lam_m
          sp_vals(SP_LAMBDA_WU)   = Wavelen
          sp_vals(SP_LAMD_RAD)    = lam_m / dpup_m
          sp_vals(SP_LAMD_ARCSEC) = (lam_m / dpup_m) * ARCSEC_PER_RAD
          sp_valid = .TRUE.

C  EFL-dependent quantities: need focalRatio from a successful EFL
C  (chief + marginal ray; can fail / vignette on complex relays like a
C  coronagraph).  Left 0 (with sp_hasefl=.FALSE.) when EFL didn't run.
          IF (focalRatio > 0d0) THEN
            efl_bu = focalRatio * dpup_bu
            sp_vals(SP_EFL_BU)     = efl_bu
            sp_vals(SP_FNO)        = focalRatio
            sp_vals(SP_NYQUIST_BU) = lam_bu * focalRatio / 2d0
            sp_hasefl = .TRUE.

C  Pixel-based plate scale: exact value from dxElt(iElt) after a
C  diffraction propagation (OPD/INT/PIX) has set it; needs efl too.
            dxf_bu = 0d0
            IF (ALLOCATED(dxElt)) THEN
              IF (iElt <= SIZE(dxElt)) dxf_bu = ABS(dxElt(iElt))
            END IF
            IF (dxf_bu > 0d0) THEN
              sp_vals(SP_DX_FOCAL_BU)  = dxf_bu
              sp_vals(SP_PLATE_PX_RAD) = efl_bu / dxf_bu
              sp_vals(SP_PLATE_ARCSEC_PX) =
     &          (dxf_bu / efl_bu) * ARCSEC_PER_RAD
              sp_vals(SP_LAMD_PX) =
     &          (lam_m / dpup_m) * (efl_bu / dxf_bu)
              sp_haspx = .TRUE.
            END IF
          END IF

          IF (PRESENT(ok)) ok = .TRUE.
        END SUBROUTINE DeriveSysProps

C  Pretty-print the most recently derived properties (CLI use).
        SUBROUTINE PrintSysProps()
          IMPLICIT NONE
          IF (.NOT. sp_valid) THEN
            WRITE(*,*) ' SYSPROP: properties unavailable (load Rx,',
     &                 ' need source at infinity; EFL must succeed).'
            RETURN
          END IF
          WRITE(*,*) ' '
          WRITE(*,*) ' First-order system properties:'
          WRITE(*,900)'  Entrance pupil diameter (m)        =',
     &                sp_vals(SP_DPUP_M)
          WRITE(*,900)'  Central obscuration ratio          =',
     &                sp_vals(SP_OBSC)
          WRITE(*,900)'  Wavelength (m)                     =',
     &                sp_vals(SP_LAMBDA_M)
          WRITE(*,900)'  lambda/D (arcsec)                  =',
     &                sp_vals(SP_LAMD_ARCSEC)
          WRITE(*,900)'  lambda/D (rad)                     =',
     &                sp_vals(SP_LAMD_RAD)
          IF (sp_hasefl) THEN
            WRITE(*,900)'  Effective focal length (BaseUnits) =',
     &                  sp_vals(SP_EFL_BU)
            WRITE(*,900)'  F-number                           =',
     &                  sp_vals(SP_FNO)
            WRITE(*,900)'  Nyquist focal sampling (BaseUnits) =',
     &                  sp_vals(SP_NYQUIST_BU)
          ELSE
            WRITE(*,*)'  (EFL / F-number unavailable: the EFL chief+',
     &                'marginal-ray trace failed -- source at infinity?'
            WRITE(*,*)'   marginal ray vignetted on this system?)'
          END IF
          IF (sp_haspx) THEN
            WRITE(*,900)'  lambda/D (detector px)             =',
     &                  sp_vals(SP_LAMD_PX)
            WRITE(*,900)'  Plate scale (arcsec/px)            =',
     &                  sp_vals(SP_PLATE_ARCSEC_PX)
            WRITE(*,900)'  Source tilt -> shift (px/rad)      =',
     &                  sp_vals(SP_PLATE_PX_RAD)
            WRITE(*,900)'  Detector pixel pitch (BaseUnits)   =',
     &                  sp_vals(SP_DX_FOCAL_BU)
          ELSE
            WRITE(*,*)'  (pixel-based plate scale needs a diffraction',
     &                ' propagation -- run OPD, INT, or PIX'
            WRITE(*,*)'   at this element first, then re-run SYSPROP)'
          END IF
          WRITE(*,*) ' '
 900      FORMAT(1P,A,1X,E15.7)
        END SUBROUTINE PrintSysProps

      END MODULE sysprop_mod
