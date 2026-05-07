MODULE Constants

  USE Kinds

  REAL(pr), PARAMETER :: &
      EPS     = EPSILON(1.0_pr), & ! ~2.2204e-16 for REAL*8
      EPS_def = 2.2204e-18_pr,   & ! for now use this EPS to be consistent (was defined in module "nn_type")
      INF     = 1e22_pr,         &
      BIG     = 1e18_pr

  REAL(pr), PARAMETER :: &
      DEG2RAD = 1.745329251994329576923690768488612713442871888541725456097191440e-2_pr, &
      PI      = 3.141592653589793238462643383279502884197169399375105820974944592_pr,    &
      TWOPI   = 6.283185307179586476925286766559005768394338798750211641949889185_pr,    &
      HALFPI  = 1.570796326794896619231321691639751442098584699687552910487472296_pr

  REAL(pr), PARAMETER :: &
      Tan30   = 0.5773502691896257645091487805019574556476017512701268760186023265_pr, &  ! sqrt(3)/3
      COS30   = 0.8660254037844386467637231707529361834714026269051903140279034897_pr, &  ! sqrt(3)/2
      SIN30   = 0.5_pr

  REAL(pr), PARAMETER :: &
      SQRT_2  = 1.414213562373095048801688724209698078569671875376948073176679738_pr, &
      SQRT_3  = 1.732050807568877293527446341505872366942805253810380628055806979_pr, &
      SQRT_5  = 2.236067977499789696409173668731276235440618359611525724270897245_pr, &
      SQRT_6  = 2.449489742783178098197284074705891391965947480656670128432692567_pr, &
      SQRT_7  = 2.645751311064590590501615753639260425710259183082450180368334459_pr, &
      SQRT_8  = 2.828427124746190097603377448419396157139343750753896146353359476_pr, &
      SQRT_10 = 3.162277660168379331998893544432718533719555139325216826857504853_pr, &
      SQRT_11 = 3.316624790355399849114932736670686683927088545589353597058682146_pr, &
      SQRT_12 = 3.464101615137754587054892683011744733885610507620761256111613959_pr, &
      SQRT_14 = 3.741657386773941385583748732316549301756019807778726946303745467_pr, &
      SQRT_18 = 4.242640687119285146405066172629094235709015626130844219530039214_pr, &
      SQRT_20 = 4.472135954999579392818347337462552470881236719223051448541794491_pr, &
      SQRT_22 = 4.690415759823429554565630113544466280588228353411737153605701891_pr

  ! Epsilon-based derived constants for numerical analysis
  ! These scale automatically with the precision (pr)
  REAL(pr), PARAMETER :: &
      EPS_SQRT      = SQRT(EPS),              & ! ~1.49d-8  (for sqrt operations)
      EPS_CBRT      = EPS**(1.0_pr/3.0_pr),   & ! ~6.04d-6  (for cubic operations)
      EPS_2_3       = EPS**(2.0_pr/3.0_pr),   & ! ~1.64d-11 (for accumulated errors)
      EPS_SQ        = EPS*EPS                   ! ~4.93d-32 (for double precision)

  ! Ray tracing numerical tolerances (based on machine epsilon)
  ! These provide consistent, precision-aware thresholds for surface intersection
  REAL(pr), PARAMETER :: &
      TOL_LINEAR    = EPS_CBRT,               & ! ~6d-6    Linear case threshold
      TOL_TANGENT   = EPS_SQRT * 1e-6_pr,     & ! ~1d-14   Near-tangent detection
      TOL_CANCEL    = 1e-4_pr,                & ! 0.01%    Cancellation detection
      TOL_GEOM      = EPS_2_3,                & ! ~1d-11   Geometric validation
      TOL_GEOM2     = TOL_GEOM * TOL_GEOM,    & ! ~2.7e-22 Geometric tolerance squared
      TOL_NORMAL    = EPS_SQ * 1e2_pr,        & ! ~1d-30   Normal magnitude threshold
      TOL_ZERO      = EPS_CBRT * 1e-14_pr,    & ! ~6d-20   Near-zero threshold
      TOL_VERTEX    = 1e-10_pr,               & ! ~1d-10   Vertex region (paraxial approximation)
      TOL_BRENT     = 1e-14_pr                  ! ~1d-14   Brent solver tolerance

END MODULE Constants
