*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  CMSUBS FORTRAN
*
*     CMALF1   CMALF    CMCHK    CMPERM   CMPRT    CMQMUL   CMR1MD
*     CMRSWP   CMTSOL
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMALF1( FIRSTV, NEGSTP, BIGALF, BIGBND, PNORM,
     $                   JADD1 , JADD2 , PALFA1, PALFA2,
     $                   ISTATE, N, NROWA, NCTOTL,
     $                   ANORM, AP, AX, BL, BU, FEATOL, P, X )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            FIRSTV, NEGSTP
      INTEGER            ISTATE(NCTOTL)
      DOUBLE PRECISION   ANORM(*), AP(*), AX(*)
      DOUBLE PRECISION   BL(NCTOTL), BU(NCTOTL), FEATOL(NCTOTL),
     $                   P(N), X(N)
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      LOGICAL            CMDBG
      INTEGER            LCMDBG
      PARAMETER         (LCMDBG = 5)
      COMMON    /CMDEBG/ ICMDBG(LCMDBG), CMDBG
 
************************************************************************
*     CMALF1  finds steps PALFA1, PALFA2 such that
*        X + PALFA1*P  reaches a linear constraint that is currently not
*                      in the working set but is satisfied.
*        X + PALFA2*P  reaches a linear constraint that is currently not
*                      in the working set but is violated.
*     The constraints are perturbed by an amount FEATOL, so that PALFA1
*     is slightly larger than it should be,  and PALFA2 is slightly
*     smaller than it should be.  This gives some leeway later when the
*     exact steps are computed by CMALF.
*
*     Constraints in the working set are ignored  (ISTATE(j) .GE. 1).
*
*     If NEGSTP is true, the search direction will be taken to be  - P.
*
*
*     Values of ISTATE(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     The values  -2  and  -1  do not occur once a feasible point has
*     been found.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written  May 1980.
*     This version of CMALF1 dated 26-June-1986.
************************************************************************
      LOGICAL            LASTV
      INTRINSIC          ABS
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      IF (CMDBG  .AND.  ICMDBG(3) .GT. 0) WRITE (NOUT, 1100)
      LASTV  = .NOT. FIRSTV
      JADD1  = 0
      JADD2  = 0
      PALFA1 = BIGALF
 
      PALFA2 = ZERO
      IF (FIRSTV) PALFA2 = BIGALF
 
      DO 200 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS .LE. 0) THEN
            IF (J .LE. N) THEN
               ATX    = X(J)
               ATP    = P(J)
               ROWNRM = ONE
            ELSE
               I      = J - N
               ATX    = AX(I)
               ATP    = AP(I)
               ROWNRM = ONE  +  ANORM(I)
            END IF
            IF (NEGSTP) ATP = - ATP
 
            IF ( ABS( ATP ) .LE. EPSPT9*ROWNRM*PNORM) THEN
 
*              This constraint appears to be constant along P.  It is
*              not used to compute the step.  Give the residual a value
*              that can be spotted in the debug output.
 
               RES = - ONE
            ELSE IF (ATP .LE. ZERO  .AND.  JS .NE. -2) THEN
*              ---------------------------------------------------------
*              a'x  is decreasing and the lower bound is not violated.
*              ---------------------------------------------------------
*              First test for smaller PALFA1.
 
               ABSATP = - ATP
               IF (BL(J) .GT. (-BIGBND)) THEN
                  RES    = ATX - BL(J) + FEATOL(J)
                  IF (BIGALF*ABSATP .GT. ABS( RES )) THEN
                     IF (PALFA1*ABSATP .GT. RES)  THEN
                        PALFA1 = RES / ABSATP
                        JADD1  = J
                     END IF
                  END IF
               END IF
 
               IF (JS .EQ. -1) THEN
 
*                 The upper bound is violated.  Test for either larger
*                 or smaller PALFA2, depending on the value of FIRSTV.
 
                  RES    = ATX - BU(J) - FEATOL(J)
                  IF (BIGALF*ABSATP .GT. ABS( RES )) THEN
                     IF (FIRSTV  .AND.  PALFA2*ABSATP .GT. RES  .OR.
     $                    LASTV  .AND.  PALFA2*ABSATP .LT. RES) THEN
                        PALFA2 = RES / ABSATP
                        JADD2  = J
                     END IF
                  END IF
               END IF
            ELSE IF (ATP .GT. ZERO  .AND.  JS .NE. -1) THEN
*              ---------------------------------------------------------
*              a'x  is increasing and the upper bound is not violated.
*              ---------------------------------------------------------
*              Test for smaller PALFA1.
 
               IF (BU(J) .LT. BIGBND) THEN
                  RES = BU(J) - ATX + FEATOL(J)
                  IF (BIGALF*ATP .GT. ABS( RES )) THEN
                     IF (PALFA1*ATP .GT. RES) THEN
                        PALFA1 = RES / ATP
                        JADD1  = J
                     END IF
                  END IF
               END IF
 
               IF (JS .EQ. -2) THEN
 
*                 The lower bound is violated.  Test for a new PALFA2.
 
                  RES  = BL(J) - ATX - FEATOL(J)
                  IF (BIGALF*ATP .GT. ABS( RES )) THEN
                     IF (FIRSTV  .AND.  PALFA2*ATP .GT. RES  .OR.
     $                    LASTV  .AND.  PALFA2*ATP .LT. RES) THEN
                        PALFA2 = RES / ATP
                        JADD2  = J
                     END IF
                  END IF
               END IF
            END IF
 
            IF (CMDBG  .AND.  ICMDBG(3) .GT. 0)
     $         WRITE (NOUT, 1200) J, JS, FEATOL(J), RES,
     $                            ATP, JADD1, PALFA1, JADD2, PALFA2
         END IF
  200 CONTINUE
 
      RETURN
 
 1100 FORMAT(/ '    J  JS         FEATOL        RES             AP',
     $         '     JADD1       PALFA1     JADD2       PALFA2' /)
 1200 FORMAT(I5, I4, 3G15.5, 2(I6, G17.7))
 
*     End of  CMALF1.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMALF ( FIRSTV, HITLOW, ISTATE, INFORM, JADD,
     $                   N, NROWA, NCLIN, NCTOTL, NUMINF,
     $                   ALFA, PALFA, ATPHIT, BIGALF, BIGBND, PNORM,
     $                   ANORM, AP, AX, BL, BU, FEATOL, P, X )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            ISTATE(NCTOTL)
      DOUBLE PRECISION   ANORM(*), AP(*), AX(*),
     $                   BL(NCTOTL), BU(NCTOTL), FEATOL(NCTOTL),
     $                   P(N), X(N)
      LOGICAL            FIRSTV, HITLOW
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      LOGICAL            CMDBG
      INTEGER            LCMDBG
      PARAMETER         (LCMDBG = 5)
      COMMON    /CMDEBG/ ICMDBG(LCMDBG), CMDBG
 
************************************************************************
*  CMALF   finds a step ALFA such that the point x + ALFA*P reaches one
*  of the linear constraints (including bounds).  Two possible steps are
*  defined as follows...
*
*  ALFA1   is the maximum step that can be taken without violating
*          one of the linear constraints that is currently satisfied.
*  ALFA2   reaches a linear constraint that is currently violated.
*          Usually this will be the furthest such constraint along P,
*          but if FIRSTV = .TRUE. it will be the first one along P.
*          This is used only when the problem has been determined to be
*          infeasible, and the sum of infeasibilities are being
*          minimized.  (ALFA2  is not defined if NUMINF = 0.)
*
*  ALFA will usually be the minimum of ALFA1 and ALFA2.
*  ALFA could be negative (since we allow inactive constraints
*  to be violated by as much as FEATOL).  In such cases, a
*  third possible step is computed, to find the nearest satisfied
*  constraint (perturbed by FEATOL) along the direction  - P.
*  ALFA  will be reset to this step if it is shorter.  This is the
*  only case for which the final step  ALFA  does not move X exactly
*  onto a constraint (the one denoted by JADD).
*
*  Constraints in the working set are ignored  (ISTATE(j) ge 1).
*
*  JADD    denotes which linear constraint is reached.
*
*  HITLOW  indicates whether it is the lower or upper bound that
*          has restricted ALFA.
*
*  Values of ISTATE(j)....
*
*     - 2         - 1         0           1          2         3
*  a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*  The values -2 and -1 do not occur once a feasible point has been
*  found.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 66 version written  May 1980.
*  This version of  CMALF  dated  10-June-1986.
************************************************************************
      LOGICAL            HLOW1, HLOW2, LASTV, NEGSTP, STEP2
      INTRINSIC          ABS, MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      INFORM = 0
 
*     ------------------------------------------------------------------
*     First pass -- find steps to perturbed constraints, so that
*     PALFA1 will be slightly larger than the true step, and
*     PALFA2 will be slightly smaller than it should be.
*     In degenerate cases, this strategy gives us some freedom in the
*     second pass.  The general idea follows that described by P.M.J.
*     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
*     ------------------------------------------------------------------
 
      NEGSTP = .FALSE.
      CALL CMALF1( FIRSTV, NEGSTP, BIGALF, BIGBND, PNORM,
     $             JADD1, JADD2, PALFA1, PALFA2,
     $             ISTATE, N, NROWA, NCTOTL,
     $             ANORM, AP, AX, BL, BU, FEATOL, P, X )
 
      JSAVE1 = JADD1
      JSAVE2 = JADD2
 
*     ------------------------------------------------------------------
*     Second pass -- recompute step-lengths without perturbation.
*     Amongst constraints that are less than the perturbed steps,
*     choose the one (of each type) that makes the largest angle
*     with the search direction.
*     ------------------------------------------------------------------
      IF (CMDBG  .AND.  ICMDBG(3) .GT. 0) WRITE (NOUT, 1000)
      ALFA1  = BIGALF
      ALFA2  = ZERO
      IF (FIRSTV) ALFA2 = BIGALF
 
      APMAX1 = ZERO
      APMAX2 = ZERO
      ATP1   = ZERO
      ATP2   = ZERO
      HLOW1  = .FALSE.
      HLOW2  = .FALSE.
      LASTV  = .NOT. FIRSTV
 
      DO 400 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS .LE. 0) THEN
            IF (J  .LE. N)  THEN
               ATX    = X(J)
               ATP    = P(J)
               ROWNRM = ONE
            ELSE
               I      = J - N
               ATX    = AX(I)
               ATP    = AP(I)
               ROWNRM = ANORM(I) + ONE
            END IF
 
            IF ( ABS( ATP ) .LE. EPSPT9*ROWNRM*PNORM) THEN
 
*              This constraint appears to be constant along P.  It is
*              not used to compute the step.  Give the residual a value
*              that can be spotted in the debug output.
 
               RES = - ONE
            ELSE IF (ATP .LE. ZERO  .AND.  JS .NE. -2) THEN
*              ---------------------------------------------------------
*              a'x  is decreasing.
*              ---------------------------------------------------------
*              The lower bound is satisfied.  Test for smaller ALFA1.
 
               ABSATP = - ATP
               IF (BL(J) .GT. (-BIGBND)) THEN
                  RES    = ATX - BL(J)
                  IF (PALFA1*ABSATP .GE. RES  .OR.  J .EQ. JSAVE1) THEN
                     IF (APMAX1*ROWNRM*PNORM .LT. ABSATP) THEN
                        APMAX1 = ABSATP / (ROWNRM*PNORM)
                        ALFA1  = RES / ABSATP
                        JADD1  = J
                        ATP1   = ATP
                        HLOW1  = .TRUE.
                     END IF
                  END IF
               END IF
 
               IF (JS. EQ. -1)  THEN
 
*                 The upper bound is violated.  Test for either a bigger
*                 or smaller ALFA2,  depending on the value of FIRSTV.
 
                  RES    = ATX - BU(J)
                  IF (     (FIRSTV  .AND.  PALFA2*ABSATP .GE. RES
     $                 .OR.  LASTV  .AND.  PALFA2*ABSATP .LE. RES)
     $                 .OR.  J .EQ.  JSAVE2) THEN
                     IF (APMAX2*ROWNRM*PNORM .LT. ABSATP) THEN
                        APMAX2 = ABSATP / (ROWNRM*PNORM)
                        IF      (ABSATP .GE. ONE          ) THEN
                           ALFA2 = RES / ABSATP
                        ELSE IF (RES    .LT. BIGALF*ABSATP) THEN
                           ALFA2 = RES / ABSATP
                        ELSE
                           ALFA2 = BIGALF
                        END IF
                        JADD2  = J
                        ATP2   = ATP
                        HLOW2  = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF (ATP .GT. ZERO  .AND.  JS .NE.  -1)  THEN
*              ---------------------------------------------------------
*              a'x  is increasing and the upper bound is not violated.
*              ---------------------------------------------------------
*              Test for smaller ALFA1.
 
               IF (BU(J) .LT. BIGBND) THEN
                  RES = BU(J) - ATX
                  IF (PALFA1*ATP .GE. RES  .OR.  J .EQ. JSAVE1) THEN
                     IF (APMAX1*ROWNRM*PNORM .LT. ATP) THEN
                        APMAX1 = ATP / (ROWNRM*PNORM)
                        ALFA1  = RES / ATP
                        JADD1  = J
                        ATP1   = ATP
                        HLOW1  = .FALSE.
                     END IF
                  END IF
               END IF
 
               IF (JS .EQ. -2)  THEN
 
*                 The lower bound is violated.  Test for a new ALFA2.
 
                  RES    = BL(J) - ATX
                  IF (     (FIRSTV  .AND.  PALFA2*ATP .GE. RES
     $                 .OR.  LASTV  .AND.  PALFA2*ATP .LE. RES)
     $                 .OR.  J .EQ.  JSAVE2) THEN
                     IF (APMAX2*ROWNRM*PNORM .LT. ATP) THEN
                        APMAX2 = ATP / (ROWNRM*PNORM)
                        IF      (ATP .GE. ONE       ) THEN
                           ALFA2 = RES / ATP
                        ELSE IF (RES .LT. BIGALF*ATP) THEN
                           ALFA2 = RES / ATP
                        ELSE
                           ALFA2 = BIGALF
                        END IF
                        JADD2  = J
                        ATP2   = ATP
                        HLOW2  = .TRUE.
                     END IF
                  END IF
               END IF
            END IF
 
            IF (CMDBG  .AND.  ICMDBG(3) .GT. 0)
     $      WRITE (NOUT, 1200) J, JS, FEATOL(J), RES, ATP, JADD1,
     $                         ALFA1, JADD2, ALFA2
         END IF
  400 CONTINUE
 
*     ==================================================================
*     Determine ALFA, the step to be taken.
*     ==================================================================
*     In the infeasible case, check whether to take the step ALFA2
*     rather than ALFA1...
 
      STEP2 = NUMINF .GT. 0  .AND.  JADD2 .GT. 0
 
*     We do so if ALFA2 is less than ALFA1 or (if FIRSTV is false)
*     lies in the range  (ALFA1, PALFA1)  and has a smaller value of
*     ATP.
 
      STEP2 = STEP2 .AND. (ALFA2 .LT. ALFA1   .OR.   LASTV  .AND.
     $                     ALFA2 .LE. PALFA1  .AND.  APMAX2 .GE. APMAX1)
 
      IF (STEP2) THEN
         ALFA   = ALFA2
         PALFA  = PALFA2
         JADD   = JADD2
         ATPHIT = ATP2
         HITLOW = HLOW2
      ELSE
         ALFA   = ALFA1
         PALFA  = PALFA1
         JADD   = JADD1
         ATPHIT = ATP1
         HITLOW = HLOW1
 
*        If ALFA1 is negative, the constraint to be added (JADD)
*        remains unchanged, but ALFA may be shortened to the step
*        to the nearest perturbed satisfied constraint along  - P.
 
         NEGSTP = ALFA .LT. ZERO
         IF (NEGSTP) THEN
            CALL CMALF1( FIRSTV, NEGSTP, BIGALF, BIGBND, PNORM,
     $                   JADD1, JADD2, PALFA1, PALFA2,
     $                   ISTATE, N, NROWA, NCTOTL,
     $                   ANORM, AP, AX, BL, BU, FEATOL, P, X )
 
            IF (CMDBG  .AND.  ICMDBG(1) .GT. 0)
     $         WRITE (NOUT, 9000) ALFA, PALFA1
 
            ALFA = - MIN( ABS( ALFA ), PALFA1 )
         END IF
      END IF
 
*     Test for undefined or infinite step.
 
      IF (JADD .EQ. 0) THEN
         ALFA   = BIGALF
         PALFA  = BIGALF
      END IF
 
      IF (ALFA .GE. BIGALF) INFORM = 3
      IF (CMDBG  .AND.  ICMDBG(1) .GT. 0  .AND.  INFORM .GT. 0)
     $   WRITE (NOUT, 9010) JADD, ALFA
      RETURN
 
 1000 FORMAT(/ ' CMALF  entered'
     $       / '    J  JS         FEATOL        RES             AP',
     $         '     JADD1        ALFA1     JADD2        ALFA2 '/)
 1200 FORMAT( I5, I4, 3G15.5, 2(I6, G17.7) )
 9000 FORMAT(/ ' //CMALF //  Negative step',
     $       / ' //CMALF //           ALFA          PALFA'
     $       / ' //CMALF //', 2G15.4 )
 9010 FORMAT(/ ' //CMALF //  Unbounded step.'
     $       / ' //CMALF //  JADD           ALFA'
     $       / ' //CMALF //  ', I4, G15.4 )
 
*     End of  CMALF .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMCHK ( NERROR, MSGLVL, COLD, USERKX,
     $                   LIWORK, LWORK, LITOTL, LWTOTL,
     $                   N, NCLIN, NCNLN,
     $                   ISTATE, KX, NAMED, NAMES, LENNAM,
     $                   BL, BU, X )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*8        NAMES(*)
      LOGICAL            COLD, NAMED, USERKX
      INTEGER            ISTATE(N+NCLIN+NCNLN), KX(N)
      DOUBLE PRECISION   BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN), X(N)
 
      COMMON    /SOL1CM/ NOUT
 
************************************************************************
*  CMCHK   checks the data input to various optimizers.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 66 version written 10-May-1980.
*  Fortran 77 version written  5-October-1984.
*  This version of CMCHK dated  23-January-1987.
************************************************************************
      LOGICAL            OK
      INTRINSIC          ABS
      PARAMETER        ( ZERO   =  0.0D+0 , ONE    =  1.0D+0 )
 
      CHARACTER*5        ID(3)
      DATA                ID(1)   ,  ID(2)   ,  ID(3)
     $                 / 'VARBL'  , 'LNCON'  , 'NLCON'   /
 
      NERROR = 0
 
*     ------------------------------------------------------------------
*     Check that there is enough workspace to solve the problem.
*     ------------------------------------------------------------------
      OK     = LITOTL .LE. LIWORK  .AND.  LWTOTL .LE. LWORK
      IF (.NOT. OK)  THEN
         WRITE (NOUT, 1100) LIWORK, LWORK, LITOTL, LWTOTL
         NERROR = NERROR + 1
         WRITE (NOUT, 1110)
      ELSE IF (MSGLVL .GT. 0)  THEN
         WRITE (NOUT, 1100) LIWORK, LWORK, LITOTL, LWTOTL
      END IF
 
      IF (USERKX) THEN
*        ---------------------------------------------------------------
*        Check for a valid KX.
*        ---------------------------------------------------------------
         IFAIL = 1
         CALL CMPERM( KX, 1, N, IFAIL )
         IF (IFAIL .NE. 0) THEN
            WRITE (NOUT, 1300)
            NERROR = NERROR + 1
         END IF
      END IF
 
*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      DO 200 J = 1, N+NCLIN+NCNLN
         B1     = BL(J)
         B2     = BU(J)
         OK     = B1 .LE. B2
         IF (.NOT. OK)  THEN
            NERROR = NERROR + 1
            IF (J .GT. N+NCLIN)  THEN
               K  = J - N - NCLIN
               L  = 3
            ELSE IF (J .GT. N)  THEN
               K  = J - N
               L  = 2
            ELSE
               K = J
               L = 1
            END IF
            IF (.NOT. NAMED) WRITE (NOUT, 1200) ID(L), K, B1, B2
            IF (      NAMED) WRITE (NOUT, 1210) NAMES(J), B1, B2
         END IF
  200 CONTINUE
 
*     ------------------------------------------------------------------
*     If warm start, check  ISTATE.
*     ------------------------------------------------------------------
      IF (.NOT. COLD) THEN
         DO 420 J = 1, N+NCLIN+NCNLN
            IS     = ISTATE(J)
            OK     = IS .GE. (- 2)   .AND.   IS .LE. 4
            IF (.NOT. OK)  THEN
               NERROR = NERROR + 1
               WRITE (NOUT, 1500) J, IS
            END IF
  420    CONTINUE
      END IF
 
      RETURN
 
 1100 FORMAT(/ ' Workspace provided is     IW(', I6,
     $         '),  W(', I6, ').' /
     $         ' To solve problem we need  IW(', I6,
     $         '),  W(', I6, ').')
 1110 FORMAT(/ ' XXX  Not enough workspace to solve problem.')
 1200 FORMAT(/ ' XXX  The bounds on  ', A5, I3,
     $         '  are inconsistent.   BL =', G16.7, '   BU =', G16.7)
 1210 FORMAT(/ ' XXX  The bounds on  ', A8,
     $         '  are inconsistent.   BL =', G16.7, '   BU =', G16.7)
 1300 FORMAT(/ ' XXX  KX has not been supplied as a valid',
     $         '  permutation.' )
 1500 FORMAT(/ ' XXX  Component', I5, '  of  ISTATE  is out of',
     $         ' range...', I10)
 
*     End of  CMCHK .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMPERM( KX, M1, M2, IFAIL )
 
      INTEGER            IFAIL, M1, M2
      INTEGER            KX(M2)
 
      COMMON    /SOL1CM/ NOUT
 
************************************************************************
*     CMPERM checks that elements M1 to M2 of KX contain a valid
*     permutation of the integers M1 to M2. The contents of KX are
*     unchanged on exit.
*
*     SOL version of NAG Library routine M01ZBF.
*     Written by N.N.Maclaren, University of Cambridge.
*     This version of CMPERM dated 18-June-1986.
************************************************************************
 
      LOGICAL            CMDBG
      INTEGER            LCMDBG
      PARAMETER         (LCMDBG = 5)
      COMMON    /CMDEBG/ ICMDBG(LCMDBG), CMDBG
 
      INTEGER            I, IERR, J, K
      INTRINSIC          ABS
 
*     Check the parameters.
 
      IF (M2 .LT. 1  .OR.  M1 .LT. 1  .OR.  M1 .GT. M2) THEN
         IERR = 1
         IF (CMDBG  .AND.  ICMDBG(3) .GT. 0)
     $      WRITE (NOUT, FMT=1100) M1, M2
      ELSE
         IERR = 0
 
*        Check that KX is within range.
 
         DO 20 I = M1, M2
            J = KX(I)
            IF ((J .LT. M1) .OR. (J .GT. M2)) GO TO 100
            IF (I .NE. J) KX(I) = -J
   20    CONTINUE
 
*        Check that no value is repeated.
 
         DO 60 I = M1, M2
            K = - KX(I)
            IF (K .GE. 0) THEN
               J     = I
   40          KX(J) = K
               J     = K
               K     = - KX(J)
               IF (K .GT. 0) GO TO 40
               IF (J .NE. I) GO TO 120
            END IF
   60    CONTINUE
      END IF
 
*     Return
 
   80 IF (IERR .NE. 0) THEN
         IFAIL = IERR
      ELSE
         IFAIL = 0
      END IF
      RETURN
  100 IERR = 2
      WRITE (NOUT, FMT=1200) I, J
      GO TO 140
  120 IERR = 3
      WRITE (NOUT, FMT=1300) J
 
*     Restore KX.
 
  140 DO 160 I = M1, M2
         KX(I) = ABS(KX(I))
  160 CONTINUE
      GO TO 80
 
 1100 FORMAT(/ ' //CMPERM//  Illegal parameter values,'
     $       / ' //CMPERM//    M1    M1'
     $       / ' //CMPERM//', 2I6 )
 1200 FORMAT(/ ' XXX  KX(',I6,') contains an out-of-range value =', I16)
 1300 FORMAT(/ ' XXX  KX contains a duplicate value =',             I16)
 
*     End of CMPERM.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMPRT ( MSGLVL, NFREE, NROWA,
     $                   N, NCLIN, NCNLN, NCTOTL, BIGBND,
     $                   NAMED, NAMES, LENNAM,
     $                   NACTIV, ISTATE, KACTIV, KX,
     $                   A, BL, BU, C, CLAMDA, RLAMDA, X )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*8        NAMES(*)
      LOGICAL            NAMED
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), BL(NCTOTL), BU(NCTOTL), C(*),
     $                   CLAMDA(NCTOTL), RLAMDA(N), X(N)
 
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            CMDBG
      INTEGER            LCMDBG
      PARAMETER         (LCMDBG = 5)
      COMMON    /CMDEBG/ ICMDBG(LCMDBG), CMDBG
 
***********************************************************************
*  CMPRT   creates the expanded Lagrange multiplier vector CLAMDA.
*  If MSGLVL .EQ 1 or MSGLVL .GE. 10,  CMPRT prints  x,  A*x,
*  c(x),  their bounds, the multipliers, and the residuals (distance
*  to the nearer bound).
*  CMPRT is called by LSCORE and NPCORE just before exiting.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 77 version written  October 1984.
*  This version of  CMPRT  dated  10-June-1986.
***********************************************************************
      CHARACTER*2        LS, LSTATE(7)
      CHARACTER*5        ID(3), ID3
      CHARACTER*8        ID4
      EXTERNAL           DDOT
      INTRINSIC          ABS
 
      PARAMETER        ( ZERO  = 0.0D+0 )
      DATA               ID(1) / 'VARBL' /
      DATA               ID(2) / 'LNCON' /
      DATA               ID(3) / 'NLCON' /
      DATA               LSTATE(1) / '--' /, LSTATE(2) / '++' /
      DATA               LSTATE(3) / 'FR' /, LSTATE(4) / 'LL' /
      DATA               LSTATE(5) / 'UL' /, LSTATE(6) / 'EQ' /
      DATA               LSTATE(7) / 'TB' /
 
 
      NPLIN  = N     + NCLIN
      NZ     = NFREE - NACTIV
 
*     Expand multipliers for bounds, linear and nonlinear constraints
*     into the  CLAMDA  array.
 
      CALL DLOAD ( NCTOTL, ZERO, CLAMDA, 1 )
      NFIXED = N - NFREE
      DO 150 K = 1, NACTIV+NFIXED
         IF (K .LE. NACTIV) J = KACTIV(K) + N
         IF (K .GT. NACTIV) J = KX(NZ+K)
         CLAMDA(J) = RLAMDA(K)
  150 CONTINUE
 
      IF (MSGLVL .LT. 10  .AND.  MSGLVL .NE. 1) RETURN
 
      WRITE (NOUT, 1100)
      ID3 = ID(1)
 
      DO 500 J = 1, NCTOTL
         B1     = BL(J)
         B2     = BU(J)
         WLAM   = CLAMDA(J)
         IS     = ISTATE(J)
         LS     = LSTATE(IS + 3)
         IF (J .LE. N) THEN
 
*           Section 1 -- the variables  x.
*           ------------------------------
            K      = J
            V      = X(J)
 
         ELSE IF (J .LE. NPLIN) THEN
 
*           Section 2 -- the linear constraints  A*x.
*           -----------------------------------------
            IF (J .EQ. N + 1) THEN
               WRITE (NOUT, 1200)
               ID3 = ID(2)
            END IF
 
            K      = J - N
            V      = DDOT  ( N, A(K,1), NROWA, X, 1 )
         ELSE
 
*           Section 3 -- the nonlinear constraints  c(x).
*           ---------------------------------------------
 
            IF (J .EQ. NPLIN + 1) THEN
               WRITE (NOUT, 1300)
               ID3 = ID(3)
            END IF
 
            K      = J - NPLIN
            V      = C(K)
         END IF
 
*        Print a line for the j-th variable or constraint.
*        -------------------------------------------------
         RES    = V - B1
         RES2   = B2 - V
         IF (ABS(RES) .GT. ABS(RES2)) RES = RES2
         IP     = 1
         IF (B1 .LE. ( - BIGBND )) IP = 2
         IF (B2 .GE.     BIGBND  ) IP = IP + 2
         IF (NAMED) THEN
 
            ID4 = NAMES(J)
            IF (IP .EQ. 1) THEN
               WRITE (NOUT, 2100) ID4,    LS, V, B1, B2, WLAM, RES
            ELSE IF (IP .EQ. 2) THEN
               WRITE (NOUT, 2200) ID4,    LS, V,     B2, WLAM, RES
            ELSE IF (IP .EQ. 3) THEN
               WRITE (NOUT, 2300) ID4,    LS, V, B1,     WLAM, RES
            ELSE
               WRITE (NOUT, 2400) ID4,    LS, V,         WLAM, RES
           END IF
 
         ELSE
 
            IF (IP .EQ. 1) THEN
               WRITE (NOUT, 3100) ID3, K, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP .EQ. 2) THEN
               WRITE (NOUT, 3200) ID3, K, LS, V,     B2, WLAM, RES
            ELSE IF (IP .EQ. 3) THEN
               WRITE (NOUT, 3300) ID3, K, LS, V, B1,     WLAM, RES
            ELSE
               WRITE (NOUT, 3400) ID3, K, LS, V,         WLAM, RES
           END IF
         END IF
  500 CONTINUE
      RETURN
 
 1100 FORMAT(// ' Variable        State', 5X, ' Value',
     $   6X, ' Lower bound', 4X, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 1200 FORMAT(// ' Linear constr   State', 5X, ' Value',
     $   6X, ' Lower bound', 4X, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 1300 FORMAT(// ' Nonlnr constr   State', 5X, ' Value',
     $   6X, ' Lower bound', 4X, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 2100 FORMAT(1X, A8, 10X, A2, 3G16.7, G16.7, G16.4)
 2200 FORMAT(1X, A8, 10X, A2, G16.7, 5X, ' None', 6X, G16.7,
     $   G16.7, G16.4)
 2300 FORMAT(1X, A8, 10X, A2, 2G16.7, 5X, ' None', 6X, G16.7, G16.4)
 2400 FORMAT(1X, A8, 10X, A2,  G16.7, 5X, ' None', 11X, ' None',
     $   6X, G16.7, G16.4)
 3100 FORMAT(1X, A5, I3, 10X, A2, 3G16.7, G16.7, G16.4)
 3200 FORMAT(1X, A5, I3, 10X, A2,  G16.7,
     $   5X, ' None', 6X, G16.7, G16.7, G16.4)
 3300 FORMAT(1X, A5, I3, 10X, A2, 2G16.7, 5X, ' None', 6X,
     $   G16.7, G16.4)
 3400 FORMAT(1X, A5, I3, 10X, A2,  G16.7,
     $   5X, ' None', 11X, ' None', 6X, G16.7, G16.4)
 
*     End of  CMPRT
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMR1MD( N, NU, NRANK, NROWR, LENV, LENW,
     $                   R, U, V, W )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            N, NU, NRANK, NROWR, LENV, LENW
      DOUBLE PRECISION   R(NROWR,*), U(N,*), V(N), W(N)
************************************************************************
*     CMR1MD  modifies the  nrank*n  upper-triangular matrix  R  so that
*     Q*(R + v*w')  is upper triangular,  where  Q  is orthogonal,
*     v  and  w  are vectors, and the modified  R  overwrites the old.
*     Q  is the product of two sweeps of plane rotations (not stored).
*     If required,  the rotations are applied to the NU columns of
*     the matrix  U.
*
*     The matrix V*W' is an (LENV) by (LENW) matrix.
*     The vector V is overwritten.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version   October  1984.
*     This version of  CMR1MD  dated 18-September-1985.
************************************************************************
      INTRINSIC          MIN
 
      J = MIN( LENV, NRANK )
      IF (NRANK .GT. 0) THEN
 
*        ===============================================================
*        Reduce components  1  thru  (J-1)  of  V  to zero,  using a
*        backward sweep of rotations.  The rotations create a horizontal
*        spike in the  j-th  row of  R.  This row is stored in  V.
*        (Note that  DROT3G  sets  V(K) = 0  below as required.)
*        ===============================================================
         LROWJ  = N - J + 1
         VJ     = V(J)
         CALL DCOPY ( LROWJ, R(J,J), NROWR, V(J), 1 )
         LROWK  = LROWJ
         DO 400 K = J-1, 1, -1
            LROWK  = LROWK + 1
            CALL DROT3G( VJ, V(K), CS, SN )
            CALL DROT3 ( LROWK, V(K)  , 1, R(K,K), NROWR, CS, SN )
 
            IF (NU .GT. 0)
     $      CALL DROT3 ( NU   , U(J,1), N, U(K,1), N    , CS, SN )
  400    CONTINUE
 
*        ===============================================================
*        Add a multiple of elements  1  thru  LENW  of  W  to the row
*        spike of  R  (stored in elements  1  thru  N  of  V).
*        ===============================================================
         CALL DAXPY ( LENW, VJ, W, 1, V, 1 )
 
*        ===============================================================
*        Eliminate the row spike  (held in  V)  using a forward sweep
*        of rotations.
*        ===============================================================
         DO 600 K = 1, J-1
            LROWK  = LROWK - 1
            L      = K     + 1
            CALL DROT3G( R(K,K), V(K), CS, SN )
            CALL DROT3 ( LROWK, R(K,L), NROWR, V(L)  , 1, CS, SN )
 
            IF (NU .GT. 0)
     $      CALL DROT3 ( NU   , U(K,1), N    , U(J,1), N, CS, SN )
  600    CONTINUE
         CALL DCOPY ( LROWJ, V(J), 1, R(J,J), NROWR )
      END IF
 
      RETURN
 
*     End of  CMR1MD
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMRSWP( N, NU, NRANK, NROWR, I, J, R, U, V )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            N, NU, NRANK, NROWR, I, J
      DOUBLE PRECISION   R(NROWR,*), U(N,*), V(N)
 
************************************************************************
*     CMRSWP  interchanges the  I-th  and  J-th  (I .LT. J)  columns of
*     an  NRANK*N  upper-triangular matrix  R   and restores the
*     resulting matrix to upper-triangular form.  The final matrix  R
*     is equal to Q(R + VW')  where  V  and  W  are defined as
*         V   =  Rj  -  Ri      and    W  =  Ei  -  Ej
*     with  Ri  and  Rj  the Ith and Jth columns of  R,  Ei  and  Ej
*     unit vectors.
*
*     The vector V is used as workspace.  R is overwritten.  Q is the
*     product of two sweeps of plane rotations (not stored).
*     If required,  the rotations are applied to the  nu  columns of
*     the matrix  U.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of  CMRSWP  dated  18-September-1985.
************************************************************************
      INTRINSIC          MIN
      INTEGER            K, L, LENI1, LENJ, LROWJ, LROWK
      DOUBLE PRECISION   CS, SN, VJ
 
      LENJ   = MIN( J, NRANK )
      IF (LENJ .GT. 0) THEN
         CALL DCOPY ( LENJ, R(1,J), 1, V, 1 )
         IF (I .LE. NRANK) V(I)  = V(I) - R(I,I)
         LENI1 = MIN( I-1, NRANK )
         IF (LENI1 .GT. 0) THEN
            CALL DCOPY ( LENI1, R(1,I), 1, R(1,J), 1 )
            CALL DCOPY ( LENI1, V     , 1, R(1,I), 1 )
         END IF
      END IF
      IF (I .LE. NRANK) THEN
 
*        ===============================================================
*        Reduce components I thru  (LENJ-1) of V to zero,  using a
*        backward sweep of rotations.  The rotations create a horizontal
*        spike in the LENJ-th row of  R.  This row is stored in V.
*        (Note that  DROT3G  sets  V(K) = 0  below as required.)
*        ===============================================================
         LROWJ  = N - LENJ + 1
         VJ     = V(LENJ)
         CALL DCOPY ( LROWJ, R(LENJ,LENJ), NROWR, V(LENJ), 1 )
         LROWK  = LROWJ
         DO 400 K = LENJ-1, I, -1
            LROWK  = LROWK + 1
            CALL DROT3G( VJ, V(K), CS, SN )
            CALL DROT3 ( LROWK, V(K)     , 1, R(K,K), NROWR, CS, SN )
 
            IF (NU .GT. 0)
     $      CALL DROT3 ( NU   , U(LENJ,1), N, U(K,1), N    , CS, SN )
  400    CONTINUE
 
*        ===============================================================
*        Add a multiple of elements I thru J of W to the
*        horizontal spike of  R  (held in elements I thru J of V).
*        ===============================================================
         V(I) = V(I) + VJ
         V(J) = V(J) - VJ
 
*        ===============================================================
*        Eliminate the row spike  (held in V)  using a forward sweep
*        of rotations.
*        ===============================================================
         DO 600 K = I, LENJ-1
            LROWK  = LROWK - 1
            L      = K     + 1
            CALL DROT3G( R(K,K), V(K), CS, SN )
            CALL DROT3 ( LROWK, R(K,L), NROWR, V(L)     , 1, CS, SN )
 
            IF (NU .GT. 0)
     $      CALL DROT3 ( NU   , U(K,1), N    , U(LENJ,1), N, CS, SN )
  600    CONTINUE
         CALL DCOPY ( LROWJ, V(LENJ), 1, R(LENJ,LENJ), NROWR )
      END IF
 
      RETURN
 
*     End of  CMRSWP
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMTSOL( MODE, NROWT, N, T, Y )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            MODE, NROWT, N
      DOUBLE PRECISION   T(NROWT,*), Y(N)
 
************************************************************************
*     CMTSOL  solves equations involving a reverse-triangular matrix  T
*     and a right-hand-side vector  y,  returning the solution in  y.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written February-1985.
************************************************************************
      PARAMETER        ( ZERO = 0.0D+0 )
 
      N1 = N + 1
      IF (MODE .EQ. 1) THEN
 
*        Mode = 1  ---  Solve  T * y(new) = y(old).
 
         DO 100 J = 1, N
            JJ = N1 - J
            YJ = Y(J)/T(J,JJ)
            Y(J) = YJ
            L  = JJ - 1
            IF (L .GT. 0  .AND.  YJ .NE. ZERO)
     $      CALL DAXPY( L, (-YJ), T(J+1,JJ), 1, Y(J+1), 1 )
  100    CONTINUE
      ELSE
 
*        Mode = 2  ---  Solve  T' y(new) = y(old).
 
         DO 500 J = 1, N
            JJ = N1 - J
            YJ = Y(J)/T(JJ,J)
            Y(J) = YJ
            L  = JJ - 1
            IF (L .GT. 0  .AND.  YJ .NE. ZERO)
     $      CALL DAXPY( L, (-YJ), T(JJ,J+1), NROWT, Y(J+1), 1 )
  500    CONTINUE
      END IF
 
*     Reverse the solution vector.
 
      IF (N .GT. 1) THEN
         L = N/2
         DO 800 J = 1, L
            JJ    = N1 - J
            YJ    = Y(J)
            Y(J)  = Y(JJ)
            Y(JJ) = YJ
  800    CONTINUE
      END IF
 
      RETURN
 
*     End of  CMTSOL.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CMQMUL( MODE, N, NZ, NFREE, NQ, UNITQ,
     $                   KX, V, ZY, WRK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            KX(N)
      DOUBLE PRECISION   V(N), ZY(NQ,*), WRK(N)
 
************************************************************************
*     CMQMUL  transforms the vector  v  in various ways using the
*     matrix  Q = ( Z  Y )  defined by the input parameters.
*
*        MODE               result
*        ----               ------
*
*          1                v = Z v
*          2                v = Y v
*          3                v = Q v
*
*     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
*     on output, v  is a full n-vector.
*
*
*          4                v = Z'v
*          5                v = Y'v
*          6                v = Q'v
*
*     On input,  v  is a full n-vector.
*     On output, v  is ordered as  ( v(free)  v(fixed) ).
*
*          7                v = Y'v
*          8                v = Q'v
*
*     On input,  v  is a full n-vector.
*     On output, v  is as in modes 5 and 6 except that v(fixed) is not
*     set.
*
*     Modes  1, 4, 7 and 8  do not involve  v(fixed).
*     Original F66 version  April 1983.
*     Fortran 77 version written  9-February-1985.
*     Level 2 BLAS added 10-June-1986.
*     This version of CMQMUL dated 10-June-1986.
************************************************************************
      EXTERNAL           DDOT
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      NFIXED = N - NFREE
      J1     = 1
      J2     = NFREE
      IF (MODE .EQ. 1  .OR.  MODE .EQ. 4) J2 = NZ
      IF (MODE .EQ. 2  .OR.  MODE .EQ. 5  .OR.  MODE .EQ. 7) J1 = NZ + 1
      LENV   = J2 - J1 + 1
      IF (MODE .LE. 3) THEN
*        ===============================================================
*        Mode = 1, 2  or  3.
*        ===============================================================
 
         IF (NFREE .GT. 0) CALL DLOAD ( NFREE, ZERO, WRK, 1 )
 
*        Copy  v(fixed)  into the end of  wrk.
 
         IF (MODE .GE. 2  .AND.  NFIXED .GT. 0)
     $      CALL DCOPY ( NFIXED, V(NFREE+1), 1, WRK(NFREE+1), 1 )
 
*        Set  WRK  =  relevant part of  ZY * V.
 
         IF (LENV .GT. 0)  THEN
            IF (UNITQ) THEN
               CALL DCOPY ( LENV, V(J1), 1, WRK(J1), 1 )
            ELSE
               CALL DGEMV ( 'N', NFREE, J2-J1+1, ONE, ZY(1,J1), NQ,
     $                      V(J1), 1, ONE, WRK, 1 )
            END IF
         END IF
 
*        Expand  WRK  into  V  as a full n-vector.
 
         CALL DLOAD ( N, ZERO, V, 1 )
         DO 220 K = 1, NFREE
            J    = KX(K)
            V(J) = WRK(K)
  220    CONTINUE
 
*        Copy  WRK(fixed)  into the appropriate parts of  V.
 
         IF (MODE .GT. 1)  THEN
            DO 320 L = 1, NFIXED
               J       = KX(NFREE+L)
               V(J)    = WRK(NFREE+L)
  320       CONTINUE
         END IF
 
      ELSE
*        ===============================================================
*        Mode = 4, 5, 6, 7  or  8.
*        ===============================================================
*        Put the fixed components of  V  into the end of  WRK.
 
         IF (MODE .EQ. 5  .OR.  MODE .EQ. 6)  THEN
            DO 420 L = 1, NFIXED
               J            = KX(NFREE+L)
               WRK(NFREE+L) = V(J)
  420       CONTINUE
         END IF
 
*        Put the free  components of  V  into the beginning of  WRK.
 
         IF (NFREE .GT. 0)  THEN
            DO 520 K = 1, NFREE
               J      = KX(K)
               WRK(K) = V(J)
  520       CONTINUE
 
*           Set  V  =  relevant part of  ZY' * WRK.
 
            IF (LENV .GT. 0)  THEN
               IF (UNITQ) THEN
                  CALL DCOPY ( LENV, WRK(J1), 1, V(J1), 1 )
               ELSE
                  CALL DGEMV ( 'T', NFREE, J2-J1+1, ONE, ZY(1,J1), NQ,
     $                         WRK, 1, ZERO, V(J1), 1 )
               END IF
            END IF
         END IF
 
*        Copy the fixed components of  WRK  into the end of  V.
 
         IF (NFIXED .GT. 0  .AND.  (MODE .EQ. 5  .OR.  MODE .EQ. 6))
     $      CALL DCOPY ( NFIXED, WRK(NFREE+1), 1, V(NFREE+1), 1 )
      END IF
 
      RETURN
 
*     End of  CMQMUL.
 
      END
