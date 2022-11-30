*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  NPSUBS FORTRAN
*
*     NPALF    NPCHKD   NPCORE   NPCRSH   NPDFLT   NPFD     NPFEAS
*     NPFILE   NPIQP    NPKEY    NPLOC    NPMRT    NPOPTN   NPPRT
*     NPRSET   NPSETX   NPSRCH   NPUPDT   NPSOL
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPALF ( INFORM, N, NCLIN, NCNLN,
     $                   ALFA, ALFMIN, ALFMAX, BIGBND, DXNORM,
     $                   ANORM, ADX, AX, BL, BU,
     $                   DSLK, DX, SLK, X )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION   ANORM(*), ADX(*), AX(*), BL(*), BU(*),
     $                   DSLK(*), DX(N), SLK(*), X(N)
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      LOGICAL            CMDBG
      INTEGER            LCMDBG
      PARAMETER         (LCMDBG = 5)
      COMMON    /CMDEBG/ ICMDBG(LCMDBG), CMDBG
 
************************************************************************
*  NPALF   finds a step ALFA such that the point x + ALFA*P reaches one
*  of the slacks or linear constraints.  The step ALFA is the maximum
*  step that can be taken without violating one of the slacks or linear
*  constraints that is currently satisfied.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 77 version written  June 1986.
*  This version of NPALF dated  27-June-1986.
************************************************************************
      INTRINSIC          ABS, MAX, MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      IF (CMDBG  .AND.  ICMDBG(3) .GT. 0) WRITE (NOUT, 1000)
 
      ALFA   = ALFMAX
      J      = 1
 
*+    WHILE (J .LE. N+NCLIN+NCNLN .AND. ALFA .GT. ALFMIN) DO
  100 IF    (J .LE. N+NCLIN+NCNLN .AND. ALFA .GT. ALFMIN) THEN
 
         IF      (J .LE. N      ) THEN
            AXI    =  X(J)
            ADXI   = DX(J)
            ROWNRM = ONE
         ELSE IF (J .LE. N+NCLIN) THEN
            I      = J - N
            AXI    = AX(I)
            ADXI   = ADX(I)
            ROWNRM = ANORM(I) + ONE
         ELSE
            I      = J - N - NCLIN
            AXI    = SLK(I)
            ADXI   = DSLK(I)
            ROWNRM = ONE
         END IF
 
         RES = - ONE
         IF (ADXI .LE. - EPSPT9*ROWNRM*DXNORM) THEN
 
*           Constraint decreasing.
 
            ADXI = - ADXI
            IF (BL(J) .GT. -BIGBND) RES = AXI   - BL(J)
 
         ELSE IF (ADXI .GT.   EPSPT9*ROWNRM*DXNORM) THEN
 
*           Constraint increasing.
 
            IF (BU(J) .LT.  BIGBND) RES = BU(J) - AXI
 
         END IF
 
         IF (RES .GT. ZERO  .AND.  ALFA*ADXI .GT. RES)
     $      ALFA  = RES / ADXI
 
         IF (CMDBG  .AND.  ICMDBG(3) .GT. 0)
     $      WRITE (NOUT, 1200) J, RES, ADXI, ALFA
 
         J = J + 1
         GO TO 100
*+    END WHILE
      END IF
 
*     ==================================================================
*     Determine ALFA, the bound on the step to be taken.
*     ==================================================================
      ALFA   = MAX( ALFA, ALFMIN )
 
      INFORM = 0
      IF (ALFA .GE. ALFMAX) INFORM = 1
 
      IF (CMDBG  .AND.  ICMDBG(1) .GT. 0  .AND.  INFORM .GT. 0)
     $   WRITE (NOUT, 9010) ALFA
 
      RETURN
 
 1000 FORMAT(/ ' NPALF  entered'
     $       / '    J            RES             AP           ALFA '/)
 1200 FORMAT( I5, 3G15.5 )
 9010 FORMAT(/ ' //NPALF //  No finite step.'
     $       / ' //NPALF //             ALFA'
     $       / ' //NPALF //  ', G15.4 )
 
*     End of  NPALF .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPCHKD( INFORM, MSGNP, NSTATE, LVLDER, NFUN, NGRAD,
     $                   NROWJ, NROWUJ, N, NCNLN,
     $                   CONFUN, OBJFUN, NEEDC,
     $                   BIGBND, EPSRF, CDINT, FDINT,
     $                   FDCHK, FDNORM, OBJF, XNORM,
     $                   BL, BU, C, C1, CJAC, UJAC, CJDX,
     $                   DX, GRAD, UGRAD, HFORWD, HCNTRL,
     $                   X, WRK1, WRK2, W, LENW )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   C(*), C1(*), CJAC(NROWJ,*), UJAC(NROWUJ,*),
     $                   CJDX(*)
      DOUBLE PRECISION   BL(N), BU(N), DX(N), GRAD(N), UGRAD(N), X(N)
      DOUBLE PRECISION   HFORWD(*), HCNTRL(*)
      DOUBLE PRECISION   WRK1(N+NCNLN), WRK2(N+NCNLN), W(LENW)
      EXTERNAL           CONFUN, OBJFUN
 
************************************************************************
*  NPCHKD  performs the following...
*  (1)  Computes the objective and constraint values OBJF and C.
*  (2)  Evaluates the user-provided gradients in UJAC and UGRAD.
*  (3)  Counts the missing gradients.
*  (4)  Loads the known gradients into GRAD and CJAC.
*  (5)  Checks that the known gradients are programmed correctly.
*  (6)  Computes the missing gradient elements.
*
*  Systems Optimization Laboratory, Stanford University, California.
*  Original version written 4-September-1985.
*  This version of NPCHKD dated  14-July-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)
 
      LOGICAL            CENTRL, NEEDFD
      PARAMETER        ( RDUMMY =-11111.0)
 
      INFORM = 0
      MODE   = 2
      NFDIFF = 0
      NCDIFF = 0
      NCSET  = N*NCNLN
 
      IF (NCNLN .GT. 0) THEN
*        ===============================================================
*        Compute the constraints and Jacobian matrix.
*        ===============================================================
*        If some derivatives are missing, load the Jacobian with dummy
*        values.  Any elements left unaltered after the call to CONFUN
*        must be estimated.  A record of the missing Jacobian elements
*        is stored in  UJAC.
 
         NEEDFD = LVLDER .EQ. 0  .OR.  LVLDER .EQ. 1
 
         IF (NEEDFD) THEN
            DO 100 J = 1, N
               CALL DLOAD ( NCNLN, RDUMMY, UJAC(1,J), 1 )
  100       CONTINUE
         END IF
 
         CALL ILOAD ( NCNLN, (1), NEEDC, 1 )
 
         CALL CONFUN( MODE, NCNLN, N, NROWUJ,
     $                NEEDC, X, C, UJAC, NSTATE )
         IF (MODE .LT. 0) GO TO 999
 
         DO 110 J = 1, N
            CALL DCOPY ( NCNLN, UJAC(1,J), 1, CJAC(1,J), 1 )
  110    CONTINUE
 
         IF (NEEDFD) THEN
 
*           Count the number of missing Jacobian elements.
 
            DO 220 J = 1, N
               DO 210 I = 1, NCNLN
                  IF (UJAC(I,J) .EQ. RDUMMY) NCDIFF = NCDIFF + 1
  210          CONTINUE
  220       CONTINUE
 
            NCSET = NCSET - NCDIFF
            IF (NSTATE .EQ. 1) THEN
               IF (NCDIFF .EQ. 0) THEN
                  IF (LVLDER .EQ. 0) LVLDER = 2
                  IF (LVLDER .EQ. 1) LVLDER = 3
                  WRITE (NOUT, 1000) LVLDER
               ELSE
                  WRITE (NOUT, 1100) NCSET, N*NCNLN, NCDIFF
               END IF
            END IF
         END IF
      END IF
 
*     ==================================================================
*     Repeat the procedure above for the objective function.
*     ==================================================================
      NEEDFD = LVLDER .EQ. 0  .OR.  LVLDER .EQ. 2
 
      IF (NEEDFD)
     $   CALL DLOAD ( N, RDUMMY, UGRAD, 1 )
 
      CALL OBJFUN( MODE, N, X, OBJF, UGRAD, NSTATE )
      IF (MODE .LT. 0) GO TO 999
 
      CALL DCOPY ( N, UGRAD, 1, GRAD, 1 )
 
      IF (NEEDFD) THEN
 
*        Count the number of missing gradient elements.
 
         DO 300 J = 1, N
            IF (UGRAD(J) .EQ. RDUMMY) NFDIFF = NFDIFF + 1
  300    CONTINUE
 
         IF (NSTATE .EQ. 1) THEN
            IF (NFDIFF .EQ. 0) THEN
               IF (LVLDER .EQ. 0) LVLDER = 1
               IF (LVLDER .EQ. 2) LVLDER = 3
               WRITE (NOUT, 2000) LVLDER
            ELSE
               WRITE (NOUT, 2100) N - NFDIFF, N, NFDIFF
            END IF
         END IF
      END IF
 
      NFUN  = NFUN  + 1
      NGRAD = NGRAD + 1
 
*     ==================================================================
*     Check whatever gradient elements have been provided.
*     ==================================================================
      IF (LVRFYC .GE. 0) THEN
         IF (NCSET .GT. 0) THEN
            CALL CHKJAC( INFORM, LVLDER, MSGNP,
     $                   NCSET, N, NCNLN, NROWJ, NROWUJ,
     $                   BIGBND, EPSRF, EPSPT3, FDCHK, XNORM,
     $                   CONFUN, NEEDC,
     $                   BL, BU, C, C1, CJAC, UJAC, CJDX,
     $                   DX, WRK1, X, WRK2, W, LENW )
            IF (INFORM .LT. 0) GO TO 800
         END IF
 
         IF (NFDIFF .LT. N) THEN
            CALL CHKGRD( INFORM, MSGNP, N,
     $                   BIGBND, EPSRF, EPSPT3, FDCHK, OBJF, XNORM,
     $                   OBJFUN,
     $                   BL, BU, GRAD, UGRAD, DX, X, WRK1, W, LENW )
            IF (INFORM .LT. 0) GO TO 800
         END IF
      END IF
 
      NEEDFD = NCDIFF .GT. 0  .OR.  NFDIFF .GT. 0
      IF (NEEDFD) THEN
*        ===============================================================
*        Compute the missing gradient elements.
*        ===============================================================
         CALL CHFD  ( INFORM, MSGNP, LVLDER,
     $                N, NCNLN, NROWJ, NROWUJ,
     $                BIGBND, EPSRF, FDNORM, OBJF,
     $                OBJFUN, CONFUN, NEEDC,
     $                BL, BU, C, C1, CJDX, CJAC, UJAC,
     $                GRAD, UGRAD, HFORWD, HCNTRL, X,
     $                DX, W, LENW )
 
         IF (INFORM .LT. 0) GO TO 800
 
         IF (LFDSET .GT. 0) THEN
            CENTRL = .FALSE.
            CALL NPFD  ( CENTRL, INFORM,
     $                   NROWJ, NROWUJ, N, NCNLN,
     $                   BIGBND, CDINT, FDINT, FDNORM, OBJF,
     $                   CONFUN, OBJFUN, NEEDC,
     $                   BL, BU, C, C1, CJDX, CJAC, UJAC,
     $                   GRAD, UGRAD, HFORWD, HCNTRL, X,
     $                   W, LENW )
 
            IF (INFORM .LT. 0) GO TO 800
         END IF
      END IF
 
  800 RETURN
 
*     The user requested termination.
 
  999 INFORM = MODE
      RETURN
 
 1000 FORMAT(//' All Jacobian elements have been set.  ',
     $         ' Derivative level increased to ', I4 )
 1100 FORMAT(//' The user sets ', I6, '   out of', I6,
     $         '   Jacobian elements.'
     $       / ' Each iteration, ', I6,
     $         '   Jacobian elements will be estimated numerically.' )
 2000 FORMAT(//' All objective gradient elements have been set.  ',
     $         ' Derivative level increased to ', I4 )
 2100 FORMAT(//' The user sets ', I6, '   out of', I6,
     $         '   objective gradient elements.'
     $       / ' Each iteration, ', I6,
     $         '   gradient elements will be estimated numerically.' )
 
*     End of  NPCHKD.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPCORE( NAMED, NAMES, UNITQ, INFORM, MAJITS,
     $                   N, NCLIN, NCNLN, NCTOTL, NACTIV, NFREE, NZ,
     $                   NROWA, NROWJ, NROWUJ, NROWQP, NROWR,
     $                   NFUN, NGRAD, ISTATE, KACTIV, KX,
     $                   OBJF, FDNORM, XNORM, OBJFUN, CONFUN,
     $                   AQP, AX, BL, BU, C, CJAC, UJAC, CLAMDA,
     $                   FEATOL, GRAD, UGRAD, R, X, IW, W, LENW )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            NAMED
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      INTEGER            IW(*)
      DOUBLE PRECISION   AQP(NROWQP,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     $                   C(*), CJAC(NROWJ,*), UJAC(NROWUJ,*)
      DOUBLE PRECISION   CLAMDA(NCTOTL), FEATOL(NCTOTL), GRAD(N),
     $                   UGRAD(N), R(NROWR,*), X(N)
      DOUBLE PRECISION   W(LENW)
      EXTERNAL           OBJFUN, CONFUN
 
      DOUBLE PRECISION   ASIZE, DTMAX, DTMIN
      CHARACTER*8        NAMES(*)
 
************************************************************************
*  NPCORE  is the core routine for  NPSOL,  a sequential quadratic
*  programming (SQP) method for nonlinearly constrained optimization.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version      February-1982.
*  This version of NPCORE dated  4-August-1986.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT , NCOLT , NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE , DTMAX , DTMIN
      COMMON    /SOL6CM/ RCNDBD, RFROBN, DRMAX, DRMIN
 
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      PARAMETER         (LENNP = 35)
      COMMON    /SOL1NP/ LOCNP(LENNP)
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      PARAMETER         (LDBG = 5)
      LOGICAL            CMDBG, NPDBG
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG
 
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS
 
      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)
 
      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)
 
      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)
 
      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP
 
      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, IPADNP(15)
 
      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , EPSRF , ETA   , FDINT , FTOL  ,
     $                   RPADNP(24)
 
      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)
 
      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )
 
      LOGICAL            GOODGQ, NEWGQ
      LOGICAL            CENTRL, CONVRG, CONVPT, DONE  , ERROR , FEASQP
      LOGICAL            INFEAS, NEEDFD, OPTIML, OVERFL, UNITQ
      LOGICAL            KTCOND(2)
      INTRINSIC          ABS   , MAX   , MIN   , MOD   , REAL  , SQRT
      EXTERNAL           DDIV  , DDOT  , DNRM2
 
      CHARACTER*4        LSUMRY
      CHARACTER*2        JOB
      PARAMETER        ( JOB  = 'NP' )
      PARAMETER        ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      PARAMETER        ( GROWTH=1.0D+2                              )
 
*     Specify machine-dependent parameters.
 
      EPSMCH = WMACH(3)
      FLMAX  = WMACH(7)
      RTMAX  = WMACH(8)
 
      LANORM = LOCLS( 2)
      LRPQ   = LOCLS( 5)
      LQRWRK = LOCLS( 6)
      LHPQ   = LOCLS( 8)
      LGQ    = LOCLS( 9)
      LRLAM  = LOCLS(10)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK1  = LOCLS(14)
      LQPTOL = LOCLS(15)
 
      LIPERM = LOCNP( 2)
      LAQP   = LOCNP( 3)
      LADX   = LOCNP( 4)
      LBL    = LOCNP( 5)
      LBU    = LOCNP( 6)
      LDX    = LOCNP( 7)
      LGQ1   = LOCNP( 8)
      LX1    = LOCNP(11)
      LWRK2  = LOCNP(12)
      LCS1   = LOCNP(13)
      LCS2   = LOCNP(14)
      LC1MUL = LOCNP(15)
      LCMUL  = LOCNP(16)
      LCJDX1 = LOCNP(17)
      LDLAM  = LOCNP(18)
      LDSLK  = LOCNP(19)
      LRHO   = LOCNP(20)
      LWRK3  = LOCNP(21)
      LSLK1  = LOCNP(22)
      LSLK   = LOCNP(23)
      LNEEDC = LOCNP(24)
      LHFRWD = LOCNP(25)
      LHCTRL = LOCNP(26)
 
      LCJAC1 = LAQP   + NCLIN
      LCJDX  = LADX   + NCLIN
      LVIOLN = LWRK3
 
*     Initialize
 
      LSUMRY = '    '
      NQPINF = 0
 
      NPLIN  = N     + NCLIN
      NCQP   = NCLIN + NCNLN
      NL     = MIN( NPLIN + 1, NCTOTL )
 
      NROWJ1 = MAX( NCQP , 1 )
 
      NEEDFD = LVLDER .EQ. 0  .OR.  LVLDER .EQ. 2
     $                        .OR. (LVLDER .EQ. 1  .AND.  NCNLN .GT. 0)
 
      ALFA   = ZERO
      ALFDX  = ZERO
      RTFTOL = SQRT( FTOL )
      ROOTN  = SQRT( REAL(N) )
 
*     If debug printing is required,  turn off any extensive printing
*     until iteration  IDBG.
 
      MSGSV1 = MSGNP
      MSGSV2 = MSGQP
      IF (IDBG .LE. NMAJOR  .AND.  IDBG .GT. 0) THEN
         MSGNP = 0
         IF (MSGSV1 .GE. 5) MSGNP = 5
         MSGQP = 0
         IF (MSGSV2 .GE. 5) MSGQP = 5
      END IF
 
*     ------------------------------------------------------------------
*     Information from the feasibility phase will be used to generate a
*     hot start for the first QP subproblem.
*     ------------------------------------------------------------------
      CALL DCOPY ( NCTOTL, FEATOL, 1, W(LQPTOL), 1 )
 
      MAJITS = 0
      NSTATE = 0
 
      LVLDIF = 0
      IF (NEEDFD) LVLDIF = 1
 
      OBJALF = OBJF
      IF (NCNLN .GT. 0) THEN
         OBJALF = OBJALF - DDOT  ( NCNLN, W(LCMUL), 1, C, 1 )
 
         INCRUN = .TRUE.
         RHONRM = ZERO
         RHODMP = ONE
         SCALE  = ONE
         CALL DLOAD ( NCNLN, (ZERO), W(LRHO), 1 )
      END IF
 
      NEWGQ  = .FALSE.
 
*+    REPEAT
*+       REPEAT
 
  100       CENTRL = LVLDIF .EQ. 2
 
            IF (NEWGQ) THEN
               IF (NEEDFD) THEN
*                 ------------------------------------------------------
*                 Compute any missing gradient elements and the
*                 transformed gradient of the objective.
*                 ------------------------------------------------------
                  CALL NPFD  ( CENTRL, MODE,
     $                         NROWJ, NROWUJ, N, NCNLN,
     $                         BIGBND, CDINT, FDINT, FDNORM, OBJF,
     $                         CONFUN, OBJFUN, IW(LNEEDC),
     $                         BL, BU, C, W(LWRK2), W(LWRK3),CJAC,UJAC,
     $                         GRAD, UGRAD, W(LHFRWD), W(LHCTRL), X,
     $                         W, LENW )
                  INFORM = MODE
                  IF (MODE .LT. 0) GO TO 800
 
               END IF
 
               CALL DCOPY ( N, GRAD, 1, W(LGQ), 1 )
               CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ,
     $                      KX, W(LGQ), W(LZY), W(LWRK1) )
 
               NEWGQ  = .FALSE.
            END IF
 
*           ============================================================
*           (1) Solve an inequality quadratic program (IQP) for the
*               search direction and multiplier estimates.
*           (2) For each nonlinear inequality constraint,  compute
*               the slack variable for which the merit function is
*               minimized.
*           (3) Compute the search direction for the slack variables
*               and multipliers.
*
*           Note that the array VIOLN is WRK3.
*           ============================================================
            CALL NPIQP ( FEASQP, UNITQ, NQPERR, MINITS,
     $                   N, NCLIN, NCNLN, NROWA, NROWJ, NROWQP,NROWR,
     $                   LINACT, NLNACT, NACTIV, NFREE, NZ, NUMINF,
     $                   ISTATE, KACTIV, KX,
     $                   DXNORM, GDX, QPCURV,
     $                   AQP, W(LADX), W(LANORM), AX, BL, BU,
     $                   C, CJAC, CLAMDA, W(LCMUL), W(LCS1),
     $                   W(LDLAM), W(LDSLK), W(LDX), W(LBL), W(LBU),
     $                   W(LQPTOL), R, W(LRHO), W(LSLK), W(LVIOLN), X,
     $                   W(LWTINF), IW, W )
 
            IF (FEASQP) THEN
               NQPINF = 0
            ELSE
               NQPINF = NQPINF + 1
               LSUMRY(2:2) = 'Infeasible subproblem'
            END IF
 
*           ============================================================
*           Compute quantities needed for the convergence test.
*           ============================================================
*           Compute the norms of the projected gradient and the
*           gradient with respect to the free variables.
 
            GZNORM = ZERO
            IF (NZ .GT. 0)
     $         GZNORM = DNRM2 ( NZ   , W(LGQ), 1 )
            GFNORM = GZNORM
            IF (NFREE .GT. 0  .AND.  NACTIV .GT. 0)
     $         GFNORM = DNRM2 ( NFREE, W(LGQ), 1 )
 
*           If the forward-difference estimate of the transformed
*           gradient of the Lagrangian function is small,  switch to
*           central differences, recompute the derivatives and re-solve
*           the QP.
 
            GOODGQ = .TRUE.
            IF (NEEDFD  .AND.  .NOT. CENTRL) THEN
 
               GLNORM = DNRM2 ( N, W(LHPQ), 1 )
               IF (NCNLN .EQ. 0) THEN
                  CNORM = ZERO
               ELSE
                  CNORM = DNRM2 ( NCNLN, C, 1 )
               END IF
 
               GLTEST = (ONE + ABS(OBJF) + ABS(CNORM))*EPSRF/FDNORM
               IF (GLNORM .LE. GLTEST) THEN
                  GOODGQ      = .FALSE.
                  LSUMRY(3:3) = 'Central differences'
                  LVLDIF      = 2
                  NEWGQ       = .TRUE.
               END IF
 
            END IF
 
*+       UNTIL     (GOODGQ)
         IF (.NOT.  GOODGQ ) GO TO 100
 
*        ===============================================================
*        (1) Compute the number of constraints that are violated by more
*            than FEATOL.
*        (2) Compute the 2-norm of the residuals of the constraints in
*            the QP working set.
*        ===============================================================
         CALL NPFEAS( N, NCLIN, NCNLN, ISTATE,
     $                BIGBND, CVNORM, ERRMAX, JMAX, NVIOL,
     $                AX, BL, BU, C, FEATOL, X, W(LWRK2) )
 
*        Define small quantities that reflect the magnitude of OBJF and
*        the norm of GRAD(free).
 
         OBJSIZ = ONE + ABS( OBJF )
         XSIZE  = ONE +  XNORM
         GTEST  = MAX( OBJSIZ, GFNORM )
         DINKY  = RTFTOL * GTEST
 
         IF (NACTIV .EQ. 0) THEN
            CONDT = ZERO
         ELSE IF (NACTIV .EQ. 1) THEN
            CONDT = DTMIN
         ELSE
            CONDT = DDIV  ( DTMAX, DTMIN, OVERFL )
         END IF
 
         CALL DCOND ( N, R, NROWR+1, DRMAX, DRMIN )
 
         CONDH = DDIV  ( DRMAX, DRMIN, OVERFL )
         IF (CONDH .LT. RTMAX) THEN
            CONDH = CONDH*CONDH
         ELSE
            CONDH = FLMAX
         END IF
 
         IF (NZ .EQ. 0) THEN
            CONDHZ = ONE
         ELSE IF (NZ .EQ. N) THEN
            CONDHZ = CONDH
         ELSE
            CALL DCOND ( NZ, R, NROWR+1, DRZMAX, DRZMIN )
            CONDHZ = DDIV  ( DRZMAX, DRZMIN, OVERFL )
            IF (CONDHZ .LT. RTMAX) THEN
               CONDHZ = CONDHZ*CONDHZ
            ELSE
               CONDHZ = FLMAX
            END IF
         END IF
 
*        ---------------------------------------------------------------
*        Test for convergence.
*        The point test CONVPT checks for a K-T point at the initial
*        point or after a large change in X.
*        ---------------------------------------------------------------
         CONVPT    = GZNORM .LE. EPSPT8*GTEST  .AND.  NVIOL  .EQ. 0
 
         KTCOND(1) = GZNORM .LT. DINKY
         KTCOND(2) = NVIOL  .EQ. 0
         OPTIML    = KTCOND(1)  .AND.  KTCOND(2)
 
         CONVRG    = MAJITS .GT. 0  .AND.  ALFDX .LE. RTFTOL*XSIZE
 
         INFEAS    =       CONVRG         .AND.  .NOT. FEASQP
     $               .OR.  NQPINF .GT. 7
 
         DONE      = CONVPT  .OR.  (CONVRG  .AND. OPTIML)
     $                       .OR.   INFEAS
 
         OBJALF = OBJF
         GRDALF = GDX
         GLF1   = GDX
         IF (NCNLN .GT. 0) THEN
            GLF1   = GLF1
     $                 - DDOT( NCNLN, W(LCJDX), 1, CLAMDA(NL), 1 )
 
*           Compute the value and directional derivative of the
*           augmented Lagrangian merit function.
*           The penalty parameters may be increased or decreased.
 
            CALL NPMRT ( FEASQP, N, NCLIN, NCNLN,
     $                   OBJALF, GRDALF, QPCURV,
     $                   ISTATE,
     $                   W(LCJDX), W(LCMUL), W(LCS1),
     $                   W(LDLAM), W(LRHO), W(LVIOLN),
     $                   W(LWRK1), W(LWRK2) )
         END IF
 
*        ===============================================================
*        Print the details of this iteration.
*        ===============================================================
         CALL NPPRT ( KTCOND, CONVRG, LSUMRY, MSGNP, MSGQP,
     $                NROWR, NROWT, N, NCLIN, NCNLN,
     $                NCTOTL, NACTIV, LINACT, NLNACT, NZ, NFREE,
     $                MAJITS, MINITS, ISTATE, ALFA, NFUN,
     $                CONDHZ, CONDH, CONDT, OBJALF, OBJF,
     $                GFNORM, GZNORM, CVNORM,
     $                AX, C, R, W(LT), W(LVIOLN), X, W(LWRK1) )
 
         ALFA  = ZERO
         ERROR = MAJITS .GE. NMAJOR
 
         IF (.NOT. (DONE  .OR.  ERROR)) THEN
            MAJITS = MAJITS + 1
 
            IF (MAJITS .EQ. IDBG) THEN
               NPDBG = .TRUE.
               CMDBG =  NPDBG
               MSGNP =  MSGSV1
               MSGQP =  MSGSV2
            END IF
 
*           Make copies of information needed for the BFGS update.
 
            CALL DCOPY ( N, X     , 1, W(LX1) , 1 )
            CALL DCOPY ( N, W(LGQ), 1, W(LGQ1), 1 )
 
            IF (NCNLN .GT. 0) THEN
               CALL DCOPY ( NCNLN, W(LCJDX), 1, W(LCJDX1), 1 )
               CALL DCOPY ( NCNLN, W(LCMUL), 1, W(LC1MUL), 1 )
               CALL DCOPY ( NCNLN, W(LSLK) , 1, W(LSLK1) , 1 )
            END IF
 
*           ============================================================
*           Compute the parameters for the linesearch.
*           ============================================================
*           Compute ALFMAX, the largest feasible step.  Also compute
*           ALFBND,  a tentative upper bound on the step.  If the
*           merit function is decreasing at ALFBND and certain
*           conditions hold,  ALFBND will be increased in multiples
*           of two (subject to not being greater than ALFMAX).
 
            ALFMAX = DDIV  ( BIGDX, DXNORM, OVERFL )
            ALFMIN = ONE
            IF (.NOT. FEASQP) ALFMIN = ZERO
 
            CALL NPALF ( INFO, N, NCLIN, NCNLN,
     $                   ALFA, ALFMIN, ALFMAX, BIGBND, DXNORM,
     $                   W(LANORM), W(LADX), AX, BL, BU,
     $                   W(LDSLK), W(LDX), W(LSLK), X )
 
            ALFMAX = ALFA
            IF (ALFMAX .LT. ONE + EPSPT3  .AND.  FEASQP)
     $         ALFMAX = ONE
 
            IF (NCNLN .EQ. 0) THEN
               ALFBND = ALFMAX
            ELSE
               IF (NEEDFD) ALFMAX = ONE
               ALFBND = MIN( ONE, ALFMAX )
            END IF
            ALFA   = ONE
 
            ALFSML = ZERO
            IF (NEEDFD  .AND. .NOT. CENTRL) THEN
               ALFSML = DDIV  ( FDNORM, DXNORM, OVERFL )
               ALFSML = MIN   ( ALFSML, ALFMAX )
            END IF
 
*           ============================================================
*           Compute the steplength using safeguarded interpolation.
*           ============================================================
            CALL NPSRCH( NEEDFD, NLSERR, N, NCNLN,
     $                   NROWJ, NROWUJ, NFUN, NGRAD,
     $                   IW(LNEEDC), CONFUN, OBJFUN,
     $                   ALFA, ALFBND, ALFMAX, ALFSML, DXNORM,
     $                   EPSRF, ETA, GDX, GRDALF, GLF1, GLF2,
     $                   OBJF, OBJALF, QPCURV, XNORM,
     $                   C, CJAC, UJAC, W(LCJDX),
     $                   W(LC1MUL), W(LCMUL), W(LCS1),
     $                   W(LCS2), W(LDX), W(LDLAM), W(LDSLK), GRAD,
     $                   UGRAD, CLAMDA(NL), W(LRHO),
     $                   W(LSLK1), W(LSLK), W(LX1), X, W, LENW )
 
*           ------------------------------------------------------------
*           NPSRCH  sets NLSERR to the following values...
*
*           NLSERR will be negative if the user set MODE LT 0.
*
*           Values of NLSERR occurring with a nonzero value of ALFA.
*           1 -- if the search was successful and ALFA LT ALFMAX.
*           2 -- if the search was successful and ALFA  = ALFMAX.
*           3 -- if the search ended after MFSRCH iterations.
*
*           Values of NLSERR occurring with a zero value of ALFA....
*           4 -- if ALFMAX was too small.
*           6 -- if no improved point could be found.
*           7 -- if the input value of GDX is non-negative.
*           ------------------------------------------------------------
            IF (NLSERR .LT. 0) THEN
               INFORM = NLSERR
               GO TO 800
            END IF
 
            ERROR  = NLSERR .GE. 4
            IF (ERROR) THEN
*              ---------------------------------------------------------
*              The linesearch failed to find a better point.
*              If exact gradients or central differences are being used,
*              or the KT conditions are satisfied, stop.  Otherwise,
*              switch to central differences and re-solve the QP.
*              ---------------------------------------------------------
               IF (NEEDFD  .AND.  .NOT. CENTRL) THEN
                  IF (.NOT. OPTIML) THEN
                     ERROR       = .FALSE.
                     LSUMRY(3:3) = 'Central differences'
                     LVLDIF      = 2
                     NEWGQ       = .TRUE.
                  END IF
               END IF
            ELSE
               IF (NEEDFD) THEN
*                 ======================================================
*                 Compute the missing gradients.
*                 ======================================================
                  MODE  = 1
                  NGRAD = NGRAD + 1
 
                  IF (NCNLN .GT. 0) THEN
                     CALL ILOAD ( NCNLN, (1), IW(LNEEDC), 1 )
 
                     CALL CONFUN( MODE, NCNLN, N, NROWUJ, IW(LNEEDC),
     $                            X, W(LWRK1), UJAC, NSTATE )
                     INFORM = MODE
                     IF (MODE .LT. 0) GO TO 800
 
                     DO 410 J = 1, N
                        CALL DCOPY (NCNLN, UJAC(1,J), 1, CJAC(1,J), 1 )
  410                CONTINUE
                  END IF
 
                  CALL OBJFUN( MODE, N, X, OBJ, UGRAD, NSTATE )
                  INFORM = MODE
                  IF (MODE .LT. 0) GO TO 800
 
                  CALL DCOPY ( N, UGRAD, 1, GRAD, 1 )
 
                  CALL NPFD  ( CENTRL, MODE,
     $                         NROWJ, NROWUJ, N, NCNLN,
     $                         BIGBND, CDINT, FDINT, FDNORM, OBJF,
     $                         CONFUN, OBJFUN, IW(LNEEDC),
     $                         BL, BU, C, W(LWRK2), W(LWRK3),CJAC,UJAC,
     $                         GRAD, UGRAD, W(LHFRWD), W(LHCTRL), X,
     $                         W, LENW )
 
                  INFORM = MODE
                  IF (MODE .LT. 0) GO TO 800
 
                  GDX  =  DDOT( N, GRAD, 1, W(LDX), 1 )
                  GLF2 =  GDX
                  IF (NCNLN .GT. 0) THEN
                     CALL DGEMV ( 'N', NCNLN, N, ONE, CJAC, NROWJ,
     $                            W(LDX), 1, ZERO, W(LCJDX), 1 )
                     GLF2 = GLF2 -
     $                      DDOT( NCNLN, W(LCJDX), 1, CLAMDA(NL), 1 )
                  END IF
               END IF
 
               CALL DCOPY ( N, GRAD, 1, W(LGQ), 1 )
               CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ,
     $                      KX, W(LGQ), W(LZY), W(LWRK1) )
 
               XNORM  = DNRM2 ( N, X, 1 )
 
               IF (NCNLN .GT. 0  .AND.  ALFA .GE. ONE)
     $            CALL DCOPY ( NCNLN, CLAMDA(NL), 1, W(LCMUL), 1 )
 
               IF (NCLIN .GT. 0)
     $            CALL DAXPY ( NCLIN, ALFA, W(LADX), 1, AX, 1 )
               ALFDX   = ALFA * DXNORM
 
*              =========================================================
*              Update the factors of the approximate Hessian of the
*              Lagrangian function.
*              =========================================================
               CALL NPUPDT( LSUMRY, UNITQ,
     $                      N, NCNLN, NFREE, NZ,
     $                      NROWJ1, NROWJ, NQ, NROWR, KX,
     $                      ALFA, GLF1, GLF2, QPCURV,
     $                      W(LCJAC1), CJAC, W(LCJDX1), W(LCJDX),
     $                      W(LCS1), W(LCS2), W(LGQ1), W(LGQ),
     $                      W(LHPQ), W(LRPQ), CLAMDA(NL), R,
     $                      W(LWRK3), W(LZY), W(LWRK2), W(LWRK1) )
 
               CALL DCOND ( N, R, NROWR+1, DRMAX, DRMIN )
               COND   = DDIV  ( DRMAX, DRMIN, OVERFL )
 
               IF (      COND   .GT. RCNDBD
     $             .OR.  RFROBN .GT. ROOTN*GROWTH*DRMAX) THEN
*                 ------------------------------------------------------
*                 Reset the condition estimator and range-space
*                 partition of Q'HQ.
*                 ------------------------------------------------------
                  IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $               WRITE (NOUT, 9000) RFROBN, DRMAX, DRMIN,COND,RCNDBD
 
                  LSUMRY(4:4) = 'Refactorize Hessian'
 
                  CALL NPRSET( UNITQ,
     $                         N, NFREE, NZ, NQ, NROWR,
     $                         IW(LIPERM), KX,
     $                         W(LGQ), R, W(LZY), W(LWRK1), W(LQRWRK) )
               END IF
            END IF
         END IF
 
*+    UNTIL     (DONE  .OR.  ERROR)
      IF (.NOT. (DONE  .OR.  ERROR) ) GO TO 100
 
*     ======================end of main loop============================
 
      IF (DONE) THEN
         IF (CONVPT  .OR.  OPTIML) THEN
            INFORM = 0
         ELSE IF (INFEAS) THEN
            INFORM = 3
         END IF
      ELSE IF (ERROR) THEN
         IF (MAJITS .GE. NMAJOR) THEN
            INFORM = 4
         ELSE IF (OPTIML) THEN
            INFORM = 1
         ELSE
            INFORM = 6
         END IF
      END IF
 
*     ------------------------------------------------------------------
*     Set  CLAMDA.  Print the full solution.
*     ------------------------------------------------------------------
  800 MSGNP = MSGSV1
      MSGQP = MSGSV2
      IF (MSGNP .GT. 0)
     $   WRITE (NOUT, 2100) INFORM, MAJITS, NFUN, NGRAD
 
      CALL CMPRT ( MSGNP, NFREE, NROWQP,
     $             N, NCLIN, NCNLN, NCTOTL, BIGBND,
     $             NAMED, NAMES, LENNAM,
     $             NACTIV, ISTATE, KACTIV, KX,
     $             AQP, BL, BU, C, CLAMDA, W(LRLAM), X )
 
      RETURN
 
 2100 FORMAT(/ ' Exit  NP phase.  INFORM = ', I2, ' MAJITS = ', I5,
     $         '   NFUN = ', I5, '   NGRAD = ', I5 )
 
 9000 FORMAT(/ ' //NPCORE//        RFROBN         DRMAX         DRMIN'
     $       / ' //NPCORE//', 1P3E14.2
     $       / ' //NPCORE//          COND        RCNDBD'
     $       / ' //NPCORE//', 1P2E14.2 )
 
*     End of  NPCORE.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPCRSH( COLD, N, NCLIN, NCNLN,
     $                   NCTOTL, NACTIV, NFREE, NZ,
     $                   ISTATE, KACTIV, BIGBND, TOLACT,
     $                   BL, BU, C )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            COLD
      INTEGER            ISTATE(NCTOTL), KACTIV(N)
      DOUBLE PRECISION   C(*), BL(NCTOTL), BU(NCTOTL)
************************************************************************
*  NPCRSH  adds indices of nonlinear constraints to the initial working
*  set.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version   14-February 1985.
*  This version of  NPCRSH  dated 14-November-1985.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            NPDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      INTRINSIC          ABS, MIN
      PARAMETER        ( ONE = 1.0D+0 )
 
      NFIXED = N      - NFREE
      LINACT = NACTIV
      NPLIN  = N      + NCLIN
 
*     If a cold start is being made, initialize the status of the QP
*     working set.  First,  if  BL(j) = BU(j),  set ISTATE(j)=3.
 
      IF (COLD) THEN
         DO  130 J = NPLIN+1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J) .EQ. BU(J)) ISTATE(J) = 3
  130    CONTINUE
      END IF
 
*     Increment NACTIV and KACTIV.
*     Ensure that the number of bounds and general constraints in the
*     QP  working set does not exceed N.
 
      DO 200 J = NPLIN+1, NCTOTL
         IF (NFIXED + NACTIV .EQ. N) ISTATE(J) = 0
         IF (ISTATE(J) .GT. 0) THEN
            NACTIV = NACTIV + 1
            KACTIV(NACTIV) = J - N
         END IF
  200 CONTINUE
 
      IF (COLD) THEN
 
*        ---------------------------------------------------------------
*        If a cold start is required, an attempt is made to add as many
*        nonlinear constraints as possible to the working set.
*        ---------------------------------------------------------------
*        The following loop finds the most violated constraint.  If
*        there is room in KACTIV, it will be added to the working set
*        and the process will be repeated.
 
 
         IS     =   1
         BIGLOW = - BIGBND
         BIGUPP =   BIGBND
         TOOBIG =   TOLACT + TOLACT
 
*        while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
  500    IF    (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) THEN
            IS   = 0
            CMIN = TOLACT
 
            DO 520 I = 1, NCNLN
               J      = NPLIN + I
               IF (ISTATE(J) .EQ. 0) THEN
                  B1     = BL(J)
                  B2     = BU(J)
                  RESL   = TOOBIG
                  RESU   = TOOBIG
                  IF (B1 .GT. BIGLOW)
     $            RESL   = ABS( C(I) - B1 ) / (ONE + ABS( B1 ))
                  IF (B2 .LT. BIGUPP)
     $            RESU   = ABS( C(I) - B2 ) / (ONE + ABS( B2 ))
                  RES    = MIN( RESL, RESU )
                  IF (RES .LT. CMIN) THEN
                     CMIN = RES
                     IMIN = I
                     IS   = 1
                     IF (RESL .GT. RESU) IS = 2
                  END IF
               END IF
  520       CONTINUE
 
            IF (IS .GT. 0) THEN
               NACTIV         = NACTIV + 1
               KACTIV(NACTIV) = NCLIN  + IMIN
               J              = NPLIN  + IMIN
               ISTATE(J)      = IS
            END IF
            GO TO 500
*        end while
         END IF
      END IF
 
*     ------------------------------------------------------------------
*     An initial working set has now been selected.
*     ------------------------------------------------------------------
      NLNACT = NACTIV - LINACT
      NZ     = NFREE  - NACTIV
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) NFIXED, LINACT, NLNACT
 
      RETURN
 
 1000 FORMAT(/ ' //NPCRSH//  Working set selected....'
     $       / ' //NPCRSH// NFIXED LINACT NLNACT     '
     $       / ' //NPCRSH//', 3I7 )
 
*     End of  NPCRSH.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPDFLT( N, NCLIN, NCNLN, LENIW, LENW, TITLE )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
      CHARACTER*(*)      TITLE
 
************************************************************************
*  NPDFLT  loads the default values of parameters not set in the options
*  file.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 77 version written 10-September-1985.
*  This version of NPDFLT dated  14-July-1986.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)
 
      LOGICAL            CMDBG, LSDBG, NPDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG
 
      LOGICAL            NEWOPT
      COMMON    /SOL7NP/ NEWOPT
      SAVE      /SOL7NP/
 
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS
 
      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)
 
      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)
 
      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)
 
      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP
 
      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, IPADNP(15)
 
      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , EPSRF , ETA   , FDINT , FTOL  ,
     $                   RPADNP(24)
 
      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)
 
      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )
 
      INTRINSIC          ABS    , LEN    , MOD
      PARAMETER        ( ZERO   =  0.0D+0, ONE    =  1.0D+0 )
      PARAMETER        ( POINT3 =  3.3D-1, POINT8 =  0.8D+0 )
      PARAMETER        ( POINT9 =  0.9D+0                   )
      PARAMETER        ( RDUMMY = -11111., IDUMMY = -11111  )
      PARAMETER        ( GIGANT =  1.0D+10*.99999           )
      PARAMETER        ( WRKTOL =  1.0D-2                   )
 
      CHARACTER*4        ICRSH(0:2)
      CHARACTER*16       KEY
      DATA                ICRSH(0),  ICRSH(1),  ICRSH(2)
     $                 / 'COLD'   , 'WARM'   , 'HOT '    /
 
      EPSMCH = WMACH( 3)
      NOUT   = WMACH(11)
      NCQP   = NCLIN + NCNLN
      NPLIN  = N     + NCLIN
      NCTOTL = NPLIN + NCNLN
 
*     Make a dummy call NPKEY to ensure that the defaults are set.
 
      CALL NPKEY ( NOUT, '*', KEY )
      NEWOPT = .TRUE.
 
*     Save the optional parameters set by the user.  The values in
*     IPRMLS, RPRMLS, IPRMNP and RPRMNP may be changed to their
*     default values.
 
      CALL ICOPY ( MXPARM, IPRMLS, 1, IPSVLS, 1 )
      CALL DCOPY ( MXPARM, RPRMLS, 1, RPSVLS, 1 )
      CALL ICOPY ( MXPARM, IPRMNP, 1, IPSVNP, 1 )
      CALL DCOPY ( MXPARM, RPRMNP, 1, RPSVNP, 1 )
 
      IF (          LCRASH .LT. 0
     $    .OR.      LCRASH .GT. 2     )   LCRASH  =  0
      IF (          LVLDER .LT. 0
     $    .OR.      LVLDER .GT. 3     )   LVLDER  =  3
      IF (          LFORMH .LT. 0
     $    .OR.      LFORMH .GT. 1     )   LFORMH  =  0
      IF (          NMAJOR .LT. 0     )   NMAJOR  = MAX(50, 3*NPLIN+
     $                                                     10*NCNLN )
      IF (          NMINOR .LT. 1     )   NMINOR  = MAX(50, 3*NCTOTL)
      IF (          MJRDBG .LT. 0     )   MJRDBG  =  0
      IF (          MNRDBG .LT. 0     )   MNRDBG  =  0
      IF (          IDBG   .LT. 0
     $    .OR.      IDBG   .GT. NMAJOR)   IDBG    =  0
      IF (          MJRDBG .EQ. 0
     $    .AND.     MNRDBG .EQ. 0     )   IDBG    = NMAJOR + 1
      IF (          MSGNP  .EQ. IDUMMY)   MSGNP   = 10
      IF (          MSGQP  .EQ. IDUMMY)   MSGQP   =  0
                                          NLNF    =  N
                                          NLNJ    =  N
                                          NLNX    =  N
      IF (          JVRFY2 .LT. 0
     $    .OR.      JVRFY2 .GT. N     )   JVRFY2  =  N
      IF (          JVRFY1 .LT. 0
     $    .OR.      JVRFY1 .GT. JVRFY2)   JVRFY1  =  1
      IF (          JVRFY4 .LT. 0
     $    .OR.      JVRFY4 .GT. N     )   JVRFY4  =  N
      IF (          JVRFY3 .LT. 0
     $    .OR.      JVRFY3 .GT. JVRFY4)   JVRFY3  =  1
      IF (          LVERFY .EQ. IDUMMY
     $    .OR.      LVERFY .GT. 13    )   LVERFY  =  0
      IF (          TOLACT .LT. ZERO
     $    .OR.      TOLACT .GE. ONE   )   TOLACT  =  WRKTOL
      IF (          TOLFEA .LT. EPSMCH
     $    .OR.      TOLFEA .GE. ONE   )   TOLFEA  =  EPSPT5
      IF (          EPSRF  .LT. EPSMCH
     $    .OR.      EPSRF  .GE. ONE   )   EPSRF   =  EPSPT9
                                          LFDSET  =  0
      IF (          FDINT  .LT. ZERO  )   LFDSET  =  2
      IF (          FDINT  .EQ. RDUMMY)   LFDSET  =  0
      IF (          FDINT  .GE. EPSMCH
     $    .AND.     FDINT  .LT. ONE   )   LFDSET  =  1
      IF (          LFDSET .EQ. 1
     $    .AND.    (CDINT  .LT. EPSMCH
     $    .OR.      CDINT  .GE. ONE  ))   CDINT   = EPSRF**POINT3
      IF (          BIGBND .LE. ZERO  )   BIGBND  = GIGANT
      IF (          BIGDX  .LE. ZERO  )   BIGDX   = MAX( GIGANT,BIGBND )
      IF (          ETA    .LT. ZERO
     $    .OR.      ETA    .GE. ONE   )   ETA     = POINT9
      IF (          FTOL   .LT. EPSRF
     $    .OR.      FTOL   .GE. ONE   )   FTOL    = EPSRF**POINT8
 
                                          DCTOL   = EPSPT5
      IF (          LVLDER .LT. 2     )   DCTOL   = EPSPT3
      IF (          CTOL   .LT. EPSMCH
     $    .OR.      CTOL   .GE. ONE   )   CTOL    = DCTOL
 
      ITMAX1    = MAX( 50, 3*(N + NCLIN + NCNLN) )
      JVERFY(1) = JVRFY1
      JVERFY(2) = JVRFY2
      JVERFY(3) = JVRFY3
      JVERFY(4) = JVRFY4
 
      NPDBG = IDBG .EQ. 0
      CMDBG = NPDBG
 
      K     = 1
      MSG1  = MJRDBG
      MSG2  = MNRDBG
      DO 200 I = 1, LDBG
         INPDBG(I) = MOD( MSG1/K, 10 )
         ICMDBG(I) = INPDBG(I)
         ILSDBG(I) = MOD( MSG2/K, 10 )
         K = K*10
  200 CONTINUE
 
      IF (MSGNP .GT. 0) THEN
 
*        Print the title.
 
         LENT = LEN( TITLE )
         IF (LENT .GT. 0) THEN
            NSPACE = (81 - LENT)/2 + 1
            WRITE (NOUT, '(///// (80A1) )')
     $         (' ', J=1, NSPACE), (TITLE(J:J), J=1,LENT)
            WRITE (NOUT, '(80A1 //)')
     $         (' ', J=1, NSPACE), ('='       , J=1,LENT)
         END IF
 
         WRITE (NOUT, 2000)
         WRITE (NOUT, 2100) NCLIN , TOLFEA, ICRSH(LCRASH) ,
     $                      N     , BIGBND, TOLACT,
     $                              BIGDX
         WRITE (NOUT, 2200) NCNLN , FTOL  , EPSRF ,
     $                      NLNJ  , CTOL  ,
     $                      NLNF  , ETA   ,
     $                      EPSMCH,
     $                      LVLDER, LVERFY
         WRITE (NOUT, 2300) NMAJOR, MSGNP,
     $                      NMINOR, MSGQP
 
         IF (LVLDER .LT. 3) THEN
            IF      (LFDSET .EQ. 0) THEN
               WRITE (NOUT, 2400)
            ELSE IF (LFDSET .EQ. 1) THEN
               WRITE (NOUT, 2401) FDINT, CDINT
            ELSE IF (LFDSET .EQ. 2) THEN
               WRITE (NOUT, 2402)
            END IF
         END IF
 
      END IF
 
      RETURN
 
 2000 FORMAT(
     $//' Parameters'
     $/ ' ----------' )
 2100 FORMAT(
     $/ ' Linear constraints.....', I10,     6X,
     $  ' Linear feasibility.....', 1PE10.2, 6X,
     $  1X, A4, ' start.............'
     $/ ' Variables..............', I10,     6X,
     $  ' Infinite bound size....', 1PE10.2, 6X,
     $  ' Crash tolerance........', 1PE10.2
     $/   24X,                      16X,
     $  ' Infinite step size.....', 1PE10.2  )
 2200 FORMAT(
     $/ ' Nonlinear constraints..', I10,     6X,
     $  ' Optimality tolerance...', 1PE10.2, 6X,
     $  ' Function precision.....', 1PE10.2
     $/ ' Nonlinear Jacobian vars', I10,     6X,
     $  ' Nonlinear feasibility..', 1PE10.2
     $/ ' Nonlinear objectiv vars', I10,     6X,
     $  ' Linesearch tolerance...', 1PE10.2
     $/ ' EPS (machine precision)', 1PE10.2, 6X,
     $  ' Derivative level.......', I10,     6X,
     $  ' Verify level...........', I10)
 2300 FORMAT(
     $/ ' Major iterations limit.', I10, 6X,
     $  ' Major print level......', I10
     $/ ' Minor iterations limit.', I10, 6X,
     $  ' Minor print level......', I10 )
 2400 FORMAT(/ ' Difference intervals to be computed.' )
 2401 FORMAT(/ ' Difference interval....', 1PE10.2, 6X,
     $         ' Central diffce interval', 1PE10.2 )
 2402 FORMAT(/ ' User-supplied difference intervals.' )
 
*     End of  NPDFLT.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPFD  ( CENTRL, INFORM,
     $                   NROWJ, NROWUJ, N, NCNLN,
     $                   BIGBND, CDINT, FDINT, FDNORM, OBJF,
     $                   CONFUN, OBJFUN, NEEDC,
     $                   BL, BU, C, C1, C2, CJAC, UJAC,
     $                   GRAD, UGRAD, HFORWD, HCNTRL, X,
     $                   W, LENW )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            CENTRL
      INTEGER            NEEDC(*)
 
      DOUBLE PRECISION   BL(N), BU(N), C(*), C1(*), C2(*),
     $                   CJAC(NROWJ,*), UJAC(NROWUJ,*)
      DOUBLE PRECISION   GRAD(N), UGRAD(N), HFORWD(N), HCNTRL(N), X(N)
      DOUBLE PRECISION   W(LENW)
      EXTERNAL           CONFUN, OBJFUN
 
************************************************************************
*  NPFD   evaluates any missing gradients.
*
*  Systems Optimization Laboratory, Stanford University, California.
*  Original version written 3-July-1986.
*  This version of NPFD   dated 14-July-1986.
************************************************************************
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
 
      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      INTRINSIC          ABS   , MAX
 
      PARAMETER         (RDUMMY=-11111.0)
      PARAMETER         (ZERO  = 0.0D+0, HALF  = 0.5D+0, ONE   = 1.0D+0)
      PARAMETER         (TWO   = 2.0D+0, THREE = 3.0D+0, FOUR  = 4.0D+0)
 
      INFORM = 0
 
*     ==================================================================
*     Use the pre-assigned difference intervals to approximate the
*     derivatives.
*     ==================================================================
*     Use either the same interval for each component (LFDSET = 1),
*     or the intervals already in HFORWD or HCNTRL (LFDSET = 0 or 2).
 
      NSTATE =   0
      MODE   =   0
 
      BIGLOW = - BIGBND
      BIGUPP =   BIGBND
 
      FDNORM =   ZERO
 
      DO 340 J  = 1, N
 
         XJ     = X(J)
         NFOUND = 0
         IF (NCDIFF .GT. 0) THEN
            DO 310 I = 1, NCNLN
               IF (UJAC(I,J) .EQ. RDUMMY) THEN
                  NEEDC(I) = 1
                  NFOUND   = NFOUND + 1
               ELSE
                  NEEDC(I) = 0
               END IF
  310       CONTINUE
         END IF
 
         IF (NFOUND .GT. 0  .OR.  UGRAD(J) .EQ. RDUMMY) THEN
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J) .GT. BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J) .LT. BIGUPP) STEPBU = BU(J) - XJ
 
            IF (CENTRL) THEN
               IF (LFDSET .EQ. 1) THEN
                  DELTA = CDINT
               ELSE
                  DELTA = HCNTRL(J)
               END IF
            ELSE
               IF (LFDSET .EQ. 1) THEN
                  DELTA = FDINT
               ELSE
                  DELTA = HFORWD(J)
               END IF
            END IF
 
            DELTA  = DELTA*(ONE + ABS(XJ))
            FDNORM = MAX (FDNORM, DELTA)
            IF (HALF*(STEPBL + STEPBU) .LT. ZERO) DELTA =  - DELTA
 
            X(J) = XJ + DELTA
            IF (NFOUND .GT. 0) THEN
               CALL CONFUN( MODE, NCNLN, N, NROWUJ,
     $                      NEEDC, X, C1, UJAC, NSTATE )
               IF (MODE .LT. 0) GO TO 999
            END IF
 
            IF (UGRAD(J) .EQ. RDUMMY) THEN
               CALL OBJFUN( MODE, N, X, OBJF1, UGRAD, NSTATE )
               IF (MODE .LT. 0) GO TO 999
            END IF
 
            IF (CENTRL) THEN
*              ---------------------------------------------------------
*              Central differences.
*              ---------------------------------------------------------
               X(J)  = XJ + DELTA + DELTA
 
               IF (NFOUND .GT. 0) THEN
                  CALL CONFUN( MODE, NCNLN, N, NROWUJ,
     $                         NEEDC, X, C2, UJAC, NSTATE )
                  IF (MODE .LT. 0) GO TO 999
 
                  DO 320 I = 1, NCNLN
                     IF (NEEDC(I) .EQ. 1)
     $                  CJAC(I,J) = (FOUR*C1(I) - THREE*C(I) - C2(I))
     $                                  / (DELTA + DELTA)
  320             CONTINUE
               END IF
 
               IF (UGRAD(J) .EQ. RDUMMY) THEN
                  CALL OBJFUN( MODE, N, X, OBJF2, UGRAD, NSTATE )
                  IF (MODE .LT. 0) GO TO 999
 
                  GRAD(J) = (FOUR*OBJF1 - THREE*OBJF - OBJF2)
     $                                  / (DELTA + DELTA)
 
               END IF
            ELSE
*              ---------------------------------------------------------
*              Forward Differences.
*              ---------------------------------------------------------
               IF (NFOUND .GT. 0) THEN
                  DO 330 I = 1, NCNLN
                     IF (NEEDC(I) .EQ. 1)
     $                  CJAC(I,J) = (C1(I) -  C(I))/  DELTA
  330             CONTINUE
               END IF
 
               IF (UGRAD(J) .EQ. RDUMMY)
     $            GRAD(J) = (OBJF1 - OBJF) /  DELTA
 
            END IF
         END IF
         X(J) = XJ
 
  340 CONTINUE
 
      RETURN
 
  999 INFORM = MODE
      RETURN
 
*     End of  NPFD  .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPFEAS( N, NCLIN, NCNLN, ISTATE,
     $                   BIGBND, CVNORM, ERRMAX, JMAX, NVIOL,
     $                   AX, BL, BU, C, FEATOL, X, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            ISTATE(N+NCLIN+NCNLN)
      DOUBLE PRECISION   AX(*), BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN)
      DOUBLE PRECISION   C(*), FEATOL(N+NCLIN+NCNLN), X(N)
      DOUBLE PRECISION   WORK(N+NCLIN+NCNLN)
************************************************************************
*  NPFEAS  computes the following...
*  (1)  The number of constraints that are violated by more
*       than  FEATOL  and the 2-norm of the constraint violations.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version      April    1984.
*  This version of  NPFEAS  dated  16-October-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            NPDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      EXTERNAL           IDAMAX, DNRM2
      INTRINSIC          ABS
      PARAMETER        ( ZERO = 0.0D+0 )
 
      BIGLOW = - BIGBND
      BIGUPP =   BIGBND
 
*     ==================================================================
*     Compute NVIOL, the number of constraints violated by more than
*     FEATOL,  and CVNORM,  the 2-norm of the constraint
*     violations and residuals of the constraints in the QP working set.
*     ==================================================================
      NVIOL  = 0
 
      DO 200 J = 1, N+NCLIN+NCNLN
         FEASJ  = FEATOL(J)
         RES    = ZERO
 
         IF (J .LE. N + NCLIN) THEN
 
*           Bound or general linear constraint.
 
            IF (J .LE. N) THEN
               CON =  X(J)
            ELSE
               CON = AX(J-N)
            END IF
 
            TOLJ   = FEASJ
         ELSE
 
*           Nonlinear constraint.
 
            CON    = C(J-N-NCLIN)
            TOLJ   = ZERO
         END IF
 
*        Check for constraint violations.
 
         IF (BL(J) .GT. BIGLOW) THEN
            RES    = BL(J) - CON
            IF (RES .GT.   FEASJ ) NVIOL = NVIOL + 1
            IF (RES .GT.    TOLJ ) GO TO 190
         END IF
 
         IF (BU(J) .LT. BIGUPP) THEN
            RES    = BU(J) - CON
            IF (RES .LT. (-FEASJ)) NVIOL = NVIOL + 1
            IF (RES .LT.  (-TOLJ)) GO TO 190
         END IF
 
*        This constraint is satisfied,  but count the residual as a
*        violation if the constraint is in the working set.
 
         IS     = ISTATE(J)
 
         IF (IS .EQ. 0) THEN
            RES = ZERO
         ELSE IF (IS .EQ. 1  .OR.  IS .LE. -2) THEN
            RES = BL(J) - CON
         ELSE IF (IS .GE. 2  .OR.  IS .EQ. -1) THEN
            RES = BU(J) - CON
         END IF
 
         IF (ABS( RES ) .GT. FEASJ) NVIOL = NVIOL + 1
 
*        Set the array of violations.
 
  190    WORK(J) = RES
  200 CONTINUE
 
      JMAX   = IDAMAX( N+NCLIN+NCNLN, WORK, 1 )
      ERRMAX = ABS ( WORK(JMAX) )
 
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) ERRMAX, JMAX
 
      CVNORM = DNRM2 ( N+NCLIN+NCNLN, WORK, 1 )
 
      RETURN
 
 1000 FORMAT(/ ' //NPFEAS//  The maximum violation is ', 1PE14.2,
     $                     ' in constraint', I5 )
 
*     End of  NPFEAS.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPFILE( IOPTNS, INFORM )
      INTEGER            IOPTNS, INFORM
 
************************************************************************
*     NPFILE  reads the options file from unit  IOPTNS  and loads the
*     options into the relevant elements of  IPRMNP  and  RPRMNP.
*
*     If  IOPTNS .lt. 0  or  IOPTNS .gt. 99  then no file is read,
*     otherwise the file associated with unit  IOPTNS  is read.
*
*     Output:
*
*         INFORM = 0  if a complete  OPTIONS  file was found
*                     (starting with  BEGIN  and ending with  END);
*                  1  if  IOPTNS .lt. 0  or  IOPTNS .gt. 99;
*                  2  if  BEGIN  was found, but end-of-file
*                     occurred before  END  was found;
*                  3  if end-of-file occurred before  BEGIN  or
*                     ENDRUN  were found;
*                  4  if  ENDRUN  was found before  BEGIN.
************************************************************************
      LOGICAL             NEWOPT
      COMMON     /SOL7NP/ NEWOPT
      SAVE       /SOL7NP/
 
      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/
 
      EXTERNAL            MCHPAR, NPKEY
      LOGICAL             FIRST
      SAVE                FIRST , NOUT
      DATA                FIRST /.TRUE./
 
*     If first time in, set  NOUT.
*     NEWOPT is true first time into NPFILE or NPOPTN
*     and just after a call to NPSOL.
 
      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT = WMACH(11)
      END IF
 
      CALL OPFILE( IOPTNS, NOUT, INFORM, NPKEY )
 
      RETURN
 
*     End of  NPFILE.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPIQP ( FEASQP, UNITQ, NQPERR, MINITS,
     $                   N, NCLIN, NCNLN, NROWA, NROWJ, NROWQP, NROWR,
     $                   LINACT, NLNACT, NACTIV, NFREE, NZ, NUMINF,
     $                   ISTATE, KACTIV, KX,
     $                   DXNORM, GDX, QPCURV,
     $                   AQP, ADX, ANORM, AX, BL, BU,
     $                   C, CJAC, CLAMDA, CMUL, CS,
     $                   DLAM, DSLK, DX, QPBL, QPBU, QPTOL,
     $                   R, RHO, SLK, VIOLN, X,
     $                   WTINF, IW, W )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            FEASQP, UNITQ
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      INTEGER            IW(*)
      DOUBLE PRECISION   AQP(NROWQP,*), ADX(*), ANORM(*), AX(*),
     $                   BL(*), BU(*),
     $                   C(*), CJAC(NROWJ,*), CLAMDA(*), CMUL(*), CS(*)
      DOUBLE PRECISION   DLAM(*), DSLK(*), DX(N)
      DOUBLE PRECISION   QPBL(*), QPBU(*),
     $                   QPTOL(*), R(NROWR,*), RHO(*), SLK(*),
     $                   VIOLN(*), X(N), WTINF(*)
      DOUBLE PRECISION   W(*)
 
************************************************************************
*     NPIQP   does the following:
*
*     (1)  Generate the upper and lower bounds for the QP  subproblem.
*
*     (2)  Compute the  TQ  factors of the rows of  AQP  specified by
*          the array  ISTATE.  The part of the factorization defined by
*          the first contiguous group of linear constraints does not
*          need to be recomputed.  The remaining rows (which could be
*          comprised of both linear and nonlinear constraints) are
*          included as new rows of the  TQ  factorization stored in
*          T and ZY.  Note that if there are no nonlinear constraints,
*          no factorization is required.
*
*     (3)  Solve the  QP  subproblem.
*                 minimize     1/2 (W p - d)'(Wp - d) + g'p
*
*                 subject to   qpbl .le. (  p ) .le. qpbu,
*                                        ( Ap )
*
*          where  W  is a matrix (not stored) such that  W'W = H  and
*          WQ = R,  d  is the zero vector,  and  g  is the gradient.
*          If the subproblem is infeasible, compute the point which
*          minimizes the sum of infeasibilities.
*
*    (4)   Find the value of each slack variable for which the merit
*          function is minimized.
*
*    (5)   Compute  DSLK,  DLAM  and  DX,  the search directions for
*          the slack variables, the multipliers and the variables.
*
*  Systems Optimization Laboratory, Stanford University.
*  Fortran 66 version written 10-January-1983.
*  This version of NPIQP dated 31-July-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT , NCOLT , NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE , DTMAX , DTMIN
      COMMON    /SOL6CM/ RCNDBD, RFROBN, DRMAX , DRMIN
 
      INTEGER            LOCLS
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      LOGICAL            CMDBG, LSDBG, NPDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG
 
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS
 
      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)
 
      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)
 
      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)
 
      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP
 
      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, IPADNP(15)
 
      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , EPSRF , ETA   , FDINT , FTOL  ,
     $                   RPADNP(24)
 
      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)
 
      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )
 
      CHARACTER*8        NAMES(1)
      LOGICAL            LINOBJ, OVERFL, QPNAMD, VERTEX
      INTRINSIC          ABS   , MIN   , MAX
      EXTERNAL           DDIV  , DDOT  , DNRM2
      PARAMETER        ( QPNAMD =.FALSE.,VERTEX =.FALSE. )
      PARAMETER        ( ZERO   =0.0D+0, ONE    =1.0D+0, TWO   =2.0D+0 )
      PARAMETER        ( HUNDRD =1.0D+2                                )
 
      IDBGSV = IDBG
      IF (NPDBG) THEN
         IDBG   = 0
      ELSE
         IDBG = NMINOR + 1
      END IF
      LSDBG  = NPDBG
      CMDBG  = NPDBG
      CALL ICOPY ( LDBG, ILSDBG, 1, ICMDBG, 1 )
 
      LRPQ   = LOCLS( 5)
      LRPQ0  = LOCLS( 6)
      LHPQ   = LOCLS( 8)
      LGQ    = LOCLS( 9)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWRK1  = LOCLS(14)
 
      NRPQ   = 0
      NGQ    = 1
 
      FEASQP =  .TRUE.
      LINOBJ =  .TRUE.
 
      BIGLOW = - BIGBND
      BIGUPP =   BIGBND
      SSQ1   =   ZERO
 
      NPLIN  = N     + NCLIN
      NCTOTL = NPLIN + NCNLN
      NCQP   = NCLIN + NCNLN
      NRANK  = N
      NREJTD = 0
 
*     ==================================================================
*     Generate the upper and lower bounds upon the search direction, the
*     weights on the sum of infeasibilities and the nonlinear constraint
*     violations.
*     ==================================================================
      WSCALE = - ONE
      DO 170 J = 1, NCTOTL
 
         IF (J .LE. N) THEN
            CON = X(J)
         ELSE IF (J .LE. NPLIN) THEN
            CON = AX(J-N)
         ELSE
            CON = C(J-NPLIN)
         END IF
 
         BLJ = BL(J)
         BUJ = BU(J)
         IF (BLJ .GT. BIGLOW) BLJ = BLJ - CON
         IF (BUJ .LT. BIGUPP) BUJ = BUJ - CON
 
         WEIGHT = ONE
         IF (J .LE. NPLIN) THEN
            IF (ABS(BLJ) .LE. QPTOL(J)) BLJ = ZERO
            IF (ABS(BUJ) .LE. QPTOL(J)) BUJ = ZERO
         ELSE
            I    = J - NPLIN
            VIOL = ZERO
            IF (BL(J) .GT. BIGLOW) THEN
               IF (BLJ .GT. ZERO) THEN
                  VIOL   = BLJ
                  WEIGHT = BLJ
                  IF (RHO(I) .GT. ZERO) WEIGHT = VIOL*RHO(I)
                  WSCALE = MAX( WSCALE,   WEIGHT )
                  GO TO 160
               END IF
            END IF
 
            IF (BU(J) .LT. BIGUPP) THEN
               IF (BUJ .LT. ZERO) THEN
                  VIOL   =   BUJ
                  WEIGHT = - BUJ
                  IF (RHO(I) .GT. ZERO) WEIGHT = - VIOL*RHO(I)
                  WSCALE = MAX( WSCALE, - WEIGHT )
               END IF
            END IF
 
*           Set the vector of nonlinear constraint violations.
 
  160       VIOLN(I) = VIOL
         END IF
 
         WTINF(J) = WEIGHT
         QPBL(J)  = BLJ
         QPBU(J)  = BUJ
 
  170 CONTINUE
 
      IF (WSCALE .GT. ZERO) THEN
         WSCALE = ONE/WSCALE
         CALL DSCAL ( NCTOTL, (WSCALE), WTINF, 1 )
      END IF
 
*     Set the maximum allowable condition estimator of the constraints
*     in the working set.  Note that a relatively well-conditioned
*     working set is used to start the QP iterations.
 
      CONDMX = MAX( ONE/EPSPT3, HUNDRD )
 
      IF (NCNLN .GT. 0) THEN
*        ===============================================================
*        Refactorize part of the  QP  constraint matrix.
*        ===============================================================
*        Load the new Jacobian into the  QP  matrix  A.  Compute the
*        2-norms of the rows of the Jacobian.
 
         DO 180 J = 1, N
            CALL DCOPY ( NCNLN, CJAC(1,J), 1, AQP(NCLIN+1,J), 1 )
  180    CONTINUE
 
         DO 190 J = NCLIN+1, NCQP
            ANORM(J) = DNRM2 ( N, AQP(J,1), NROWQP )
  190    CONTINUE
 
*        Count the number of linear constraints in the working set and
*        move them to the front of KACTIV.  Compute the norm of the
*        matrix of constraints in the working set.
*        Let K1  point to the first nonlinear constraint.  Constraints
*        with indices KACTIV(K1),..., KACTIV(NACTIV)  must be
*        refactorized.
 
         ASIZE  = ZERO
         LINACT = 0
         K1     = NACTIV + 1
         DO 200 K = 1, NACTIV
            I     = KACTIV(K)
            ASIZE = MAX( ASIZE, ANORM(I) )
 
            IF (I .LE. NCLIN) THEN
               LINACT = LINACT + 1
               IF (LINACT .NE. K) THEN
                  ISWAP  = KACTIV(LINACT)
                  KACTIV(LINACT) = I
                  KACTIV(K)      = ISWAP
               END IF
            ELSE
 
*              Record the old position of the 1st. nonlinear constraint.
 
               IF (K1 .GT. NACTIV) K1 = K
            END IF
  200    CONTINUE
 
         IF (NACTIV .LE. 1 )
     $      CALL DCOND ( NCQP, ANORM, 1, ASIZE, AMIN )
 
*        Compute the absolute values of the nonlinear constraints in
*        the working set.  Use DX as workspace.
 
         DO 210 K = LINACT+1, NACTIV
            J = N + KACTIV(K)
            IF (ISTATE(J) .EQ. 1) DX(K) = ABS( QPBL(J) )
            IF (ISTATE(J) .GE. 2) DX(K) = ABS( QPBU(J) )
  210    CONTINUE
 
*        Sort the elements of KACTIV corresponding to nonlinear
*        constraints in descending order of violation (i.e.,
*        the first element of KACTIV for a nonlinear constraint
*        is associated with the most violated constraint.)
*        In this way, the rows of the Jacobian corresponding
*        to the more violated constraints tend to be included
*        in the  TQ  factorization.
 
*        The sorting procedure is taken from the simple insertion
*        sort in D. Knuth, ACP Volume 3, Sorting and Searching,
*        Page 81.  It should be replaced by a faster sort if the
*        number of active nonlinear constraints becomes large.
 
         DO 230 K = LINACT+2, NACTIV
            L     = K
            VIOL  = DX(L)
            KVIOL = KACTIV(L)
*           WHILE (L .GT. LINACT+1  .AND.  DX(L-1) .LT. VIOL) DO
  220       IF    (L .GT. LINACT+1                          ) THEN
               IF (                        DX(L-1) .LT. VIOL) THEN
                  DX(L)     = DX(L-1)
                  KACTIV(L) = KACTIV(L-1)
                  L         = L - 1
                  GO TO 220
               END IF
*           END WHILE
            END IF
            DX(L)     = VIOL
            KACTIV(L) = KVIOL
  230    CONTINUE
 
         K2     = NACTIV
         NACTIV = K1     - 1
         NZ     = NFREE  - NACTIV
 
*        Update the factors  R,  T  and  Q  to include constraints
*        K1  through  K2.
 
         IF (K1 .LE. K2)
     $      CALL LSADDS( UNITQ, VERTEX,
     $                   INFORM, K1, K2, NACTIV, NARTIF, NZ, NFREE,
     $                   NRANK, NREJTD, NRPQ, NGQ,
     $                   N, NQ, NROWQP, NROWR, NROWT,
     $                   ISTATE, KACTIV, KX,
     $                   CONDMX,
     $                   AQP, R, W(LT), W(LRPQ), W(LGQ),
     $                   W(LZY), W(LWRK1), DX )
      END IF
 
*     ==================================================================
*     Solve for DX, the vector of minimum two-norm that satisfies the
*     constraints in the working set.
*     ==================================================================
      CALL NPSETX( UNITQ,
     $             NCQP, NACTIV, NFREE, NZ,
     $             N, NLNX, NCTOTL, NQ, NROWQP, NROWR, NROWT,
     $             ISTATE, KACTIV, KX,
     $             DXNORM, GDX,
     $             AQP, ADX, QPBL, QPBU, W(LRPQ), W(LRPQ0), DX, W(LGQ),
     $             R, W(LT), W(LZY), W(LWRK1) )
 
*     ==================================================================
*     Solve a quadratic program for the search direction  DX  and
*     multiplier estimates  CLAMDA.
*     ==================================================================
*     If there is no feasible point for the subproblem,  the sum of
*     infeasibilities is minimized subject to the linear constraints
*     (1  thru  JINF)  being satisfied.
 
      JINF  = N + NCLIN
 
      NTRY  = 1
*+    REPEAT
  450    CALL LSCORE( 'QP subproblem', QPNAMD, NAMES, LINOBJ, UNITQ,
     $                NQPERR, MINITS, JINF, NCQP, NCTOTL,
     $                NACTIV, NFREE, NRANK, NZ, NZ1,
     $                N, NROWQP, NROWR,
     $                ISTATE, KACTIV, KX,
     $                GDX, SSQ, SSQ1, SUMINF, NUMINF, DXNORM,
     $                QPBL, QPBU, AQP, CLAMDA, ADX,
     $                QPTOL, R, DX, IW, W )
 
         IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $      WRITE (NOUT, 8000) NQPERR
 
         NVIOL = 0
         IF (NUMINF .GT. 0) THEN
 
*           Count the violated linear constraints.
 
            DO 460 J = 1, NPLIN
               IF (ISTATE(J) .LT. 0) NVIOL = NVIOL + 1
  460       CONTINUE
 
            IF (NVIOL .GT. 0) THEN
               NTRY   = NTRY + 1
               UNITQ  = .TRUE.
               NACTIV = 0
               NFREE  = N
               NZ     = N
               CALL ILOAD ( NCTOTL, (0), ISTATE, 1 )
 
               CALL NPSETX( UNITQ,
     $                      NCQP, NACTIV, NFREE, NZ,
     $                      N, NLNX, NCTOTL, NQ, NROWQP, NROWR, NROWT,
     $                      ISTATE, KACTIV, KX,
     $                      DXNORM, GDX,
     $                      AQP, ADX, QPBL, QPBU, W(LRPQ), W(LRPQ0),
     $                      DX, W(LGQ), R, W(LT), W(LZY), W(LWRK1) )
            END IF
         END IF
      IF (.NOT. (NVIOL .EQ. 0  .OR.  NTRY .GT. 2)) GO TO 450
*+    UNTIL (    NVIOL .EQ. 0  .OR.  NTRY .GT. 2)
 
*     ==================================================================
*     Count the number of nonlinear constraint gradients in the  QP
*     working set.  Make sure that all small  QP  multipliers associated
*     with nonlinear inequality constraints have the correct sign.
*     ==================================================================
      NLNACT  = 0
      IF (NACTIV .GT. 0  .AND.  NCNLN .GT. 0) THEN
         DO 500 K = 1, NACTIV
            L     = KACTIV(K)
            IF (L .GT. NCLIN) THEN
               NLNACT = NLNACT + 1
               J      = N      + L
               IF (ISTATE(J) .EQ. 1) CLAMDA(J) = MAX( ZERO, CLAMDA(J) )
               IF (ISTATE(J) .EQ. 2) CLAMDA(J) = MIN( ZERO, CLAMDA(J) )
            END IF
  500    CONTINUE
      END IF
 
      LINACT = NACTIV - NLNACT
 
*     ------------------------------------------------------------------
*     Extract various useful quantities from the QP solution.
*     ------------------------------------------------------------------
*     Compute  HPQ = R'R(pq)  from the transformed gradient of the QP
*     objective function and  R(pq)  from the transformed residual.
 
      CALL DSCAL ( N, (-ONE), W(LRPQ), 1 )
      CALL DAXPY ( N, (-ONE), W(LGQ) , 1, W(LHPQ), 1 )
      QPCURV = TWO*SSQ
 
      IF (NCNLN .GT. 0) THEN
         IF (NUMINF .GT. 0) THEN
            FEASQP = .FALSE.
            CALL DLOAD ( NCTOTL, (ZERO), CLAMDA, 1 )
 
            IF (NZ .GT. 0) THEN
*              ---------------------------------------------------------
*              Compute a null space component for the search direction
*              as the solution of  Z'HZ(pz) = -Z'g - Z'HY(py).
*              ---------------------------------------------------------
*              Overwrite DX with the transformed search direction
*              Q'(dx).  The first NZ components of DX are zero.
 
               CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ,
     $                      KX, DX, W(LZY), W(LWRK1) )
 
*              Overwrite the first NZ components of DX with the solution
*              of  (Rz)u = -(v + w),  where  (Rz)'w = Z'g  and  v  is
*              vector of first NZ components of  R(pq).
 
               CALL DCOPY ( NZ, W(LGQ), 1, DX, 1 )
               CALL DTRSV ( 'U', 'T', 'N', NZ, R, NROWR, DX, 1 )
 
               CALL DAXPY ( NZ, (ONE), W(LRPQ), 1, DX, 1 )
 
               CALL DTRSV ( 'U', 'N', 'N', NZ, R, NROWR, DX, 1 )
               CALL DSCAL ( NZ, (-ONE), DX, 1 )
 
*              Recompute RPQ, HPQ, GDX and QPCURV.
 
               CALL DCOPY ( NLNX, DX, 1, W(LRPQ), 1 )
               CALL DTRMV ( 'U', 'N', 'N', NLNX, R, NROWR, W(LRPQ), 1 )
               IF (NLNX .LT. N)
     $            CALL DGEMV( 'N', NLNX, N-NLNX, ONE, R(1,NLNX+1),NROWR,
     $                        DX(NLNX+1), 1, ONE, W(LRPQ), 1 )
 
               GDX    = DDOT  ( N, W(LGQ) , 1, DX     , 1 )
               QPCURV = DDOT  ( N, W(LRPQ), 1, W(LRPQ), 1 )
 
               CALL CMQMUL( 3, N, NZ, NFREE, NQ, UNITQ,
     $                      KX, DX, W(LZY), W(LWRK1) )
 
*              ---------------------------------------------------------
*              Recompute ADX and the 2-norm of DX.
*              ---------------------------------------------------------
               DXNORM  = DNRM2 ( N, DX, 1 )
               IF (NCQP .GT. 0)
     $            CALL DGEMV ( 'N', NCQP, N, ONE, AQP, NROWQP,
     $                         DX, 1, ZERO, ADX, 1 )
 
               IF (NPDBG  .AND.  INPDBG(2) .GT. 0)
     $            WRITE (NOUT, 8100) (DX(J), J = 1, N)
            END IF
 
            CALL DCOPY ( NLNX, W(LRPQ), 1, W(LHPQ), 1 )
            CALL DTRMV ( 'U', 'T', 'N', NLNX, R, NROWR, W(LHPQ), 1 )
            IF (NLNX .LT. N)
     $         CALL DGEMV ( 'T', NLNX, N-NLNX, ONE, R(1,NLNX+1), NROWR,
     $                      W(LRPQ), 1, ZERO, W(LHPQ+NLNX), 1 )
         END IF
 
*        ===============================================================
*        For given values of the objective function and constraints,
*        attempt to minimize the merit function with respect to each
*        slack variable.
*        ===============================================================
         DO 600 I = 1, NCNLN
            J      = NPLIN + I
            CON    = C(I)
 
            IF (      .NOT. FEASQP  .AND.
     $          VIOLN(I) .NE. ZERO  .AND.  RHO(I) .LE. ZERO )
     $         RHO(I) = ONE
 
            QUOTNT = DDIV  ( CMUL(I), SCALE*RHO(I), OVERFL )
 
*           Define the slack variable to be  CON - MULT / RHO.
*           Force each slack to lie within its upper and lower bounds.
 
            IF (BL(J) .GT. BIGLOW) THEN
               IF (QPBL(J) .GE. - QUOTNT) THEN
                  SLK(I) = BL(J)
                  GO TO 550
               END IF
            END IF
 
            IF (BU(J) .LT. BIGUPP) THEN
               IF (QPBU(J) .LE. - QUOTNT) THEN
                  SLK(I) = BU(J)
                  GO TO 550
               END IF
            END IF
 
            SLK(I) = CON - QUOTNT
 
*           The slack has been set within its bounds.
 
  550       CS(I)  = CON - SLK(I)
 
*           ------------------------------------------------------------
*           Compute the search direction for the slacks and multipliers.
*           ------------------------------------------------------------
            DSLK(I) = ADX(NCLIN+I) + CS(I)
 
            IF (FEASQP) THEN
*
*              If any constraint is such that  (DLAM)*(C - S)  is
*              positive,  the merit function may be reduced immediately
*              by substituting the QP multiplier.
*
               DLAM(I)  = CLAMDA(J) - CMUL(I)
               IF (DLAM(I) * CS(I) .GE. ZERO) THEN
                  CMUL(I) = CLAMDA(J)
                  DLAM(I) = ZERO
               END IF
            ELSE
 
*              The  QP  subproblem was infeasible.
 
               DLAM(I) = ZERO
 
               IF (ISTATE(J) .LT. 0  .OR.  VIOLN(I) .NE. ZERO)
     $            DSLK(I)  = ZERO
 
            END IF
  600    CONTINUE
 
         IF (.NOT. FEASQP)
     $      RHONRM = DNRM2 ( NCNLN, RHO, 1 )
 
         IF (NPDBG  .AND.  INPDBG(2) .GT. 0) THEN
            WRITE (NOUT, 8200) (VIOLN(I), I=1,NCNLN)
            WRITE (NOUT, 8300) (SLK(I)  , I=1,NCNLN)
         END IF
      END IF
 
      CALL ICOPY ( LDBG, INPDBG, 1, ICMDBG, 1 )
      IDBG   = IDBGSV
 
      RETURN
 
 8000 FORMAT(/ ' //NPIQP // NQPERR'
     $       / ' //NPIQP // ',  I6 )
 8100 FORMAT(/ ' //NPIQP // DX recomputed with null space portion...'
     $       / (5G12.3))
 8200 FORMAT(/ ' //NPIQP // Violations = '/ (1P5E15.6))
 8300 FORMAT(/ ' //NPIQP // Slacks     = '/ (1P5E15.6))
 
*     End of  NPIQP .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPKEY ( NOUT, BUFFER, KEY )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*)      BUFFER
 
************************************************************************
*     NPKEY   decodes the option contained in  BUFFER  in order to set
*     a parameter value in the relevant element of the parameter arrays.
*
*
*     Input:
*
*     NOUT   A unit number for printing error messages.
*            NOUT  must be a valid unit.
*
*     Output:
*
*     KEY    The first keyword contained in BUFFER.
*
*
*     NPKEY  calls OPNUMB and the subprograms
*                 LOOKUP, SCANNR, TOKENS, UPCASE
*     (now called OPLOOK, OPSCAN, OPTOKN, OPUPPR)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of NPKEY  dated 12-July-1986.
************************************************************************
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS
 
      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)
 
      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)
 
      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)
 
      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP
 
      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, IPADNP(15)
 
      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , EPSRF , ETA   , FDINT , FTOL  ,
     $                   RPADNP(24)
 
      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)
 
      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )
 
      EXTERNAL           OPNUMB
      LOGICAL            FIRST , MORE  , NUMBER, OPNUMB, SORTED
      SAVE               FIRST
 
      PARAMETER         (     MAXKEY = 38,  MAXTIE = 19,   MAXTOK = 10)
      CHARACTER*16       KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK)
      CHARACTER*16       KEY, KEY2, KEY3, VALUE
 
      PARAMETER         (IDUMMY = -11111,  RDUMMY = -11111.0,
     $                   SORTED = .TRUE.,  ZERO   =  0.0     )
 
      DATA                FIRST
     $                  /.TRUE./
      DATA   KEYS
     $ / 'BEGIN           ',
     $   'CENTRAL         ', 'COLD            ', 'CONSTRAINTS     ',
     $   'CRASH           ', 'DEBUG           ', 'DEFAULTS        ',
     $   'DERIVATIVE      ', 'DIFFERENCE      ',
     $   'END             ', 'FEASIBILITY     ', 'FUNCTION        ',
     $   'HESSIAN         ', 'HOT             ', 'INFINITE        ',
     $   'IPRMLS          ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'LINEAR          ', 'LINESEARCH      ',
     $   'LIST            ', 'LOWER           ',
     $   'MAJOR           ', 'MINOR           ',
     $   'NOLIST          ',
     $   'NONLINEAR       ', 'OPTIMALITY      ', 'PRINT           ',
     $   'PROBLEM         ', 'ROW             ', 'RPRMLS          ',
     $   'START           ', 'STOP            ', 'UPPER           ',
     $   'VARIABLES       ', 'VERIFY          ', 'WARM            '/
 
      DATA   TIES
     $ / 'BOUND           ', 'CONSTRAINTS     ', 'DEBUG           ',
     $   'FEASIBILITY     ', 'GRADIENTS       ',
     $   'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'JACOBIAN        ', 'LEVEL           ',
     $   'NO              ',
     $   'NO.      :NUMBER',
     $   'NUMBER          ', 'OBJECTIVE       ', 'PRINT           ',
     $   'STEP            ', 'TOLERANCE       ',
     $   'VARIABLES       ', 'YES             '/
*-----------------------------------------------------------------------
 
      IF (FIRST) THEN
         FIRST  = .FALSE.
         DO 10 I = 1, MXPARM
            RPRMLS(I) = RDUMMY
            IPRMLS(I) = IDUMMY
            RPRMNP(I) = RDUMMY
            IPRMNP(I) = IDUMMY
   10    CONTINUE
      END IF
 
*     Eliminate comments and empty lines.
*     A '*' appearing anywhere in BUFFER terminates the string.
 
      I      = INDEX( BUFFER, '*' )
      IF (I .EQ. 0) THEN
         LENBUF = LEN( BUFFER )
      ELSE
         LENBUF = I - 1
      END IF
      IF (LENBUF .LE. 0) THEN
         KEY = '*'
         GO TO 900
      END IF
 
*     ------------------------------------------------------------------
*     Extract up to MAXTOK tokens from the record.
*     NTOKEN returns how many were actually found.
*     KEY, KEY2, KEY3 are the first tokens if any, otherwise blank.
*     ------------------------------------------------------------------
      NTOKEN = MAXTOK
      CALL OPTOKN( BUFFER(1:LENBUF), NTOKEN, TOKEN )
      KEY    = TOKEN(1)
      KEY2   = TOKEN(2)
      KEY3   = TOKEN(3)
 
*     Certain keywords require no action.
 
      IF (KEY .EQ. ' '     .OR.  KEY .EQ. 'BEGIN' ) GO TO 900
      IF (KEY .EQ. 'LIST'  .OR.  KEY .EQ. 'NOLIST') GO TO 900
      IF (KEY .EQ. 'END'                          ) GO TO 900
 
*     Most keywords will have an associated integer or real value,
*     so look for it no matter what the keyword.
 
      I      = 1
      NUMBER = .FALSE.
 
   50 IF (I .LT. NTOKEN  .AND.  .NOT. NUMBER) THEN
         I      = I + 1
         VALUE  = TOKEN(I)
         NUMBER = OPNUMB( VALUE )
         GO TO 50
      END IF
 
      IF (NUMBER) THEN
         READ (VALUE, '(BN, E16.0)') RVALUE
      ELSE
         RVALUE = ZERO
      END IF
 
*     Convert the keywords to their most fundamental form
*     (upper case, no abbreviations).
*     SORTED says whether the dictionaries are in alphabetic order.
*     LOCi   says where the keywords are in the dictionaries.
*     LOCi = 0 signals that the keyword wasn't there.
 
      CALL OPLOOK( MAXKEY, KEYS, SORTED, KEY , LOC1 )
      CALL OPLOOK( MAXTIE, TIES, SORTED, KEY2, LOC2 )
 
*     ------------------------------------------------------------------
*     Decide what to do about each keyword.
*     The second keyword (if any) might be needed to break ties.
*     Some seemingly redundant testing of MORE is used
*     to avoid compiler limits on the number of consecutive ELSE IFs.
*     ------------------------------------------------------------------
      MORE   = .TRUE.
      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'CENTRAL     ') THEN
            CDINT  = RVALUE
         ELSE IF (KEY .EQ. 'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY .EQ. 'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY .EQ. 'DEBUG       ') THEN
            IDBG   = RVALUE
         ELSE IF (KEY .EQ. 'DEFAULTS    ') THEN
            DO 20 I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
               IPRMNP(I) = IDUMMY
               RPRMNP(I) = RDUMMY
   20       CONTINUE
         ELSE IF (KEY .EQ. 'DERIVATIVE  ') THEN
            LVLDER = RVALUE
         ELSE IF (KEY .EQ. 'DIFFERENCE  ') THEN
            FDINT  = RVALUE
         ELSE IF (KEY .EQ. 'FEASIBILITY ') THEN
            TOLFEA = RVALUE
            CTOL   = RVALUE
         ELSE IF (KEY .EQ. 'FUNCTION    ') THEN
            EPSRF  = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF
 
      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'HESSIAN     ') THEN
            LFORMH = 1
            IF   (KEY2.EQ. 'NO          ') LFORMH = 0
         ELSE IF (KEY .EQ. 'HOT         ') THEN
            LCRASH = 2
         ELSE IF (KEY .EQ. 'INFINITE    ') THEN
              IF (KEY2.EQ. 'BOUND       ') BIGBND = RVALUE * 0.99999
              IF (KEY2.EQ. 'STEP        ') BIGDX  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'IPRMLS      ') THEN
*           Allow things like  IPRMLS 21 = 100  to set IPRMLS(21) = 100
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, I16)') IPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'ITERATIONS  ') THEN
            NMAJOR = RVALUE
         ELSE IF (KEY .EQ. 'LINEAR      ') THEN
            IF (KEY2  .EQ. 'CONSTRAINTS ') NNCLIN = RVALUE
            IF (KEY2  .EQ. 'FEASIBILITY ') TOLFEA = RVALUE
            IF (LOC2 .EQ.  0             ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'LINESEARCH  ') THEN
            ETA    = RVALUE
         ELSE IF (KEY .EQ. 'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF
 
      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'MAJOR       ') THEN
              IF (KEY2.EQ. 'DEBUG       ') MJRDBG = RVALUE
              IF (KEY2.EQ. 'ITERATIONS  ') NMAJOR = RVALUE
              IF (KEY2.EQ. 'PRINT       ') MSGNP  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'MINOR       ') THEN
              IF (KEY2.EQ. 'DEBUG       ') MNRDBG = RVALUE
              IF (KEY2.EQ. 'ITERATIONS  ') NMINOR = RVALUE
              IF (KEY2.EQ. 'PRINT       ') MSGQP  = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'NONLINEAR   ') THEN
              IF (KEY2.EQ. 'CONSTRAINTS ') NNCNLN = RVALUE
              IF (KEY2.EQ. 'FEASIBILITY ') CTOL   = RVALUE
              IF (KEY2.EQ. 'JACOBIAN    ') NLNJ   = RVALUE
              IF (KEY2.EQ. 'OBJECTIVE   ') NLNF   = RVALUE
              IF (KEY2.EQ. 'VARIABLES   ') NLNX   = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'OPTIMALITY  ') THEN
            FTOL   = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF
 
      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'PRINT       ') THEN
            MSGNP  = RVALUE
         ELSE IF (KEY .EQ. 'PROBLEM     ') THEN
              IF (KEY2.EQ. 'NUMBER      ') NPROB  = RVALUE
         ELSE IF (KEY .EQ. 'ROW         ') THEN
              IF (KEY2.EQ. 'TOLERANCE   ') CTOL   = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'RPRMLS      ') THEN
*           Allow things like  RPRMLS 21 = 2  to set RPRMLS(21) = 2.0
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, E16.0)') RPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'START       ') THEN
              IF (KEY2.EQ. 'CONSTRAINTS ') JVRFY3 = RVALUE
              IF (KEY2.EQ. 'OBJECTIVE   ') JVRFY1 = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'STOP        ') THEN
              IF (KEY2.EQ. 'CONSTRAINTS ') JVRFY4 = RVALUE
              IF (KEY2.EQ. 'OBJECTIVE   ') JVRFY2 = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE IF (KEY .EQ. 'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY .EQ. 'VARIABLES   ') THEN
            NN     = RVALUE
         ELSE IF (KEY .EQ. 'VERIFY      ') THEN
              IF (KEY2.EQ. 'OBJECTIVE   ') LVERFY =  1
              IF (KEY2.EQ. 'CONSTRAINTS ') LVERFY =  2
              IF (KEY2.EQ. 'NO          ') LVERFY = -1
              IF (KEY2.EQ. 'YES         ') LVERFY =  3
              IF (KEY2.EQ. 'GRADIENTS   ') LVERFY =  3
              IF (KEY2.EQ. 'LEVEL       ') LVERFY =  RVALUE
              IF (LOC2.EQ.  0            ) LVERFY =  3
         ELSE IF (KEY .EQ. 'WARM        ') THEN
            LCRASH = 1
         ELSE
            WRITE(NOUT, 2300) KEY
         END IF
      END IF
 
  900 RETURN
 
 2300 FORMAT(' XXX  Keyword not recognized:         ', A)
 2320 FORMAT(' XXX  Second keyword not recognized:  ', A)
 2330 FORMAT(' XXX  Third  keyword not recognized:  ', A)
 2400 FORMAT(' XXX  The PARM subscript is out of range:', I10)
 
*     End of NPKEY
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPLOC ( N, NCLIN, NCNLN, NCTOTL, LITOTL, LWTOTL)
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
************************************************************************
*     NPLOC   allocates the addresses of the work arrays for NPCORE and
*     LSCORE.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version   14-February-1985.
*     This version of  NPLOC  dated 12-July-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ
 
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      PARAMETER         (LENNP = 35)
      COMMON    /SOL1NP/ LOCNP(LENNP)
 
      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      MINIW     = LITOTL + 1
      MINW      = LWTOTL + 1
 
*     Assign array lengths that depend upon the problem dimensions.
 
      IF (NCLIN + NCNLN .EQ. 0) THEN
         LENT      = 0
         LENZY     = 0
      ELSE
         LENT  = NROWT*NCOLT
         LENZY = NQ   *NQ
      END IF
 
      IF (NCNLN .EQ. 0) THEN
         LENAQP = 0
      ELSE
         LENAQP = (NCLIN + NCNLN)*N
      END IF
 
      LKACTV    = MINIW
      LKX       = LKACTV + N
      LNEEDC    = LKX    + N
      LIPERM    = LNEEDC + NCNLN
      MINIW     = LIPERM + NCTOTL
 
      LHFRWD    = MINW
      LHCTRL    = LHFRWD + N
      LANORM    = LHCTRL + N
      LQPGQ     = LANORM + NCLIN + NCNLN
      LGQ       = LQPGQ  + N
      LRLAM     = LGQ    + N
      LT        = LRLAM  + N
      LZY       = LT     + LENT
      MINW      = LZY    + LENZY
 
      LOCLS( 1) = LKACTV
      LOCLS( 2) = LANORM
      LOCLS( 8) = LQPGQ
      LOCLS( 9) = LGQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
 
*     Assign the addresses for the workspace arrays used by  NPIQP.
 
      LQPADX    = MINW
      LQPDX     = LQPADX + NCLIN + NCNLN
      LRPQ      = LQPDX  + N
      LRPQ0     = LRPQ   + N
      LQPHZ     = LRPQ0  + N
      LWTINF    = LQPHZ  + N
      LWRK1     = LWTINF + NCTOTL
      LQPTOL    = LWRK1  + NCTOTL
      MINW      = LQPTOL + NCTOTL
 
      LOCLS( 3) = LQPADX
      LOCLS( 4) = LQPDX
      LOCLS( 5) = LRPQ
      LOCLS( 6) = LRPQ0
      LOCLS( 7) = LQPHZ
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK1
      LOCLS(15) = LQPTOL
 
*     Assign the addresses for arrays used in NPCORE.
 
      LAQP      = MINW
      LADX      = LAQP   + LENAQP
      LBL       = LADX   + NCLIN  + NCNLN
      LBU       = LBL    + NCTOTL
      LDX       = LBU    + NCTOTL
      LGQ1      = LDX    + N
      LFEATL    = LGQ1   + N
      LX1       = LFEATL + NCTOTL
      LWRK2     = LX1    + N
      MINW      = LWRK2  + NCTOTL
 
      LOCNP( 1) = LKX
      LOCNP( 2) = LIPERM
      LOCNP( 3) = LAQP
      LOCNP( 4) = LADX
      LOCNP( 5) = LBL
      LOCNP( 6) = LBU
      LOCNP( 7) = LDX
      LOCNP( 8) = LGQ1
      LOCNP(10) = LFEATL
      LOCNP(11) = LX1
      LOCNP(12) = LWRK2
 
      LCS1      = MINW
      LCS2      = LCS1   + NCNLN
      LC1MUL    = LCS2   + NCNLN
      LCMUL     = LC1MUL + NCNLN
      LCJDX     = LCMUL  + NCNLN
      LDLAM     = LCJDX  + NCNLN
      LDSLK     = LDLAM  + NCNLN
      LRHO      = LDSLK  + NCNLN
      LWRK3     = LRHO   + NCNLN
      LSLK1     = LWRK3  + NCNLN
      LSLK      = LSLK1  + NCNLN
      MINW      = LSLK   + NCNLN
 
      LOCNP(13) = LCS1
      LOCNP(14) = LCS2
      LOCNP(15) = LC1MUL
      LOCNP(16) = LCMUL
      LOCNP(17) = LCJDX
      LOCNP(18) = LDLAM
      LOCNP(19) = LDSLK
      LOCNP(20) = LRHO
      LOCNP(21) = LWRK3
      LOCNP(22) = LSLK1
      LOCNP(23) = LSLK
      LOCNP(24) = LNEEDC
 
      LCJAC     = MINW
      LGRAD     = LCJAC  + NCNLN*N
      MINW      = LGRAD  + N
 
      LOCNP(25) = LHFRWD
      LOCNP(26) = LHCTRL
      LOCNP(27) = LCJAC
      LOCNP(28) = LGRAD
 
      LITOTL    = MINIW - 1
      LWTOTL    = MINW  - 1
 
      RETURN
 
*     End of  NPLOC .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPMRT ( FEASQP, N, NCLIN, NCNLN,
     $                   OBJALF, GRDALF, QPCURV,
     $                   ISTATE,
     $                   CJDX, CMUL, CS,
     $                   DLAM, RHO, VIOLN,
     $                   WORK1, WORK2 )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
      LOGICAL            FEASQP
 
      INTEGER            ISTATE(*)
 
      DOUBLE PRECISION   CJDX(*), CMUL(*), CS(*),
     $                   DLAM(*), RHO(*), VIOLN(*)
      DOUBLE PRECISION   WORK1(*), WORK2(*)
 
************************************************************************
*  NPMRT   computes the value and directional derivative of the
*  augmented Lagrangian merit function.  The penalty parameters RHO(j)
*  are boosted if the directional derivative of the resulting augmented
*  Lagrangian function is not sufficiently negative.  If RHO needs to
*  be increased,  the perturbation with minimum two-norm is found that
*  gives a directional derivative equal to  - p'Hp.
*
*  Systems Optimization Laboratory, Stanford University, California.
*  Original version written  27-May-1985.
*  This version of  NPMRT  dated 14-November-1985.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      LOGICAL            BOOST , OVERFL
      EXTERNAL           DDIV  , DDOT  , DNRM2
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      PARAMETER        ( ZERO   = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      PARAMETER        ( TWO    = 2.0D+0                              )
 
      IF (NCNLN .EQ. 0) RETURN
 
      RTMIN  = WMACH(6)
 
      OBJALF = OBJALF - DDOT  ( NCNLN, CMUL, 1, CS, 1 )
      GRDALF = GRDALF - DDOT  ( NCNLN, DLAM, 1, CS, 1 )
 
      CALL DCOPY ( NCNLN, CS, 1, WORK1, 1 )
 
      IF (.NOT. FEASQP) THEN
         NPLIN  = N + NCLIN
 
         DO 100 I = 1, NCNLN
            IF (ISTATE(NPLIN+I) .LT. 0  .OR.  VIOLN(I) .NE. ZERO)
     $         WORK1(I) = - CJDX(I)
  100    CONTINUE
      END IF
 
      GRDALF = GRDALF + DDOT  ( NCNLN, WORK1, 1, CMUL, 1 )
 
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) QPCURV, GRDALF
 
      IF (FEASQP) THEN
 
*        Find the quantities that define  rhomin, the vector of minimum
*        two-norm such that the directional derivative is one half of
*        approximate curvature   - (dx)'H(dx).
 
         DO 350 I = 1, NCNLN
            IF (ABS( CS(I) ) .LE. RTMIN) THEN
               WORK2(I) = ZERO
            ELSE
               WORK2(I) = CS(I)**2
            END IF
  350    CONTINUE
 
         QNORM  = DNRM2 ( NCNLN, WORK2, 1 )
         TSCL   = DDIV  ( GRDALF + HALF*QPCURV, QNORM, OVERFL )
         IF (ABS( TSCL ) .LE. RHOMAX  .AND.  .NOT. OVERFL) THEN
*           ------------------------------------------------------------
*           Bounded  RHOMIN  found.  The final value of  RHO(J)  will
*           never be less than  RHOMIN(j).  If the  QP  was feasible,  a
*           trial value  RHONEW  is computed that is equal to the
*           geometric mean of the previous  RHO  and a damped value of
*           RHOMIN.  The new  RHO  is defined as  RHONEW  if it is less
*           than half the previous  RHO  and greater than  RHOMIN.
*           ------------------------------------------------------------
            SCALE  = ONE
            DO 400 I = 1, NCNLN
               RHOMIN = MAX(  (WORK2(I)/QNORM)*TSCL, ZERO )
               RHOI   = RHO(I)
 
               RHONEW = SQRT( RHOI*(RHODMP + RHOMIN) )
               IF (RHONEW .LT. HALF*RHOI  ) RHOI = RHONEW
               IF (RHOI   .LT.      RHOMIN) RHOI = RHOMIN
               RHO(I) = RHOI
  400       CONTINUE
 
            RHO1   = RHONRM
            RHONRM = DNRM2 ( NCNLN, RHO, 1 )
 
*           ------------------------------------------------------------
*           If  INCRUN = true,  there has been a run of iterations in
*           which the norm of  RHO  has not decreased.  Conversely,
*           INCRUN = false  implies that there has been a run of
*           iterations in which the norm of RHO has not increased.  If
*           INCRUN changes during this iteration the damping parameter
*           RHODMP is increased by a factor of two.  This ensures that
*           RHO(j) will oscillate only a finite number of times.
*           ------------------------------------------------------------
            BOOST  = .FALSE.
            IF (      INCRUN  .AND.  RHONRM .LT. RHO1) BOOST = .TRUE.
            IF (.NOT. INCRUN  .AND.  RHONRM .GT. RHO1) BOOST = .TRUE.
            IF (BOOST) THEN
               RHODMP = TWO*RHODMP
               INCRUN = .NOT. INCRUN
            END IF
         END IF
 
         IF (NPDBG  .AND.  INPDBG(2) .GT. 0)
     $      WRITE (NOUT, 1200) (RHO(L), L=1,NCNLN)
 
      ELSE
 
*        The  QP  was infeasible.  Do not alter the penalty parameters,
*        but compute the scale factor so that the constraint violations
*        are reduced.
 
         CALL DDSCL ( NCNLN, RHO, 1, WORK1, 1 )
         PTERM2 = DDOT  ( NCNLN, WORK1, 1, CS, 1 )
 
         SCALE  = RHOMAX
         TSCL   = DDIV  ( GRDALF, PTERM2, OVERFL )
         IF (TSCL .GT. SCALE  .AND.  TSCL .LE. RHOMAX/(ONE+RHONRM)
     $                        .AND.  .NOT. OVERFL)
     $      SCALE = TSCL
 
         CALL DCOPY ( NCNLN, CS, 1, WORK1, 1 )
      END IF
 
*     ------------------------------------------------------------------
*     Compute the new value and directional derivative of the
*     merit function.
*     ------------------------------------------------------------------
      CALL DDSCL ( NCNLN, RHO, 1, WORK1, 1 )
 
      PTERM  = DDOT  ( NCNLN, WORK1, 1, CS, 1 )
      OBJALF = OBJALF + HALF*SCALE*PTERM
 
      IF (FEASQP)
     $  PTERM2 = PTERM
 
      GRDALF = GRDALF -      SCALE*PTERM2
 
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1100) SCALE, RHONRM, GRDALF
 
      RETURN
 
 1000 FORMAT(/ ' //NPMRT //        QPCURV        GRDALF '
     $       / ' //NPMRT //', 1P2E14.2 )
 1100 FORMAT(/ ' //NPMRT //         SCALE        RHONRM        GRDALF '
     $       / ' //NPMRT //', 1P3E14.2 )
 1200 FORMAT(/ ' //NPMRT //  Penalty parameters =       '/ (1P5E15.6))
 
*     End of  NPMRT .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPOPTN( STRING )
      CHARACTER*(*)      STRING
 
************************************************************************
*     NPOPTN  loads the option supplied in STRING into the relevant
*     element of IPRMLS, RPRMLS, IPRMNP or RPRMNP.
************************************************************************
 
      LOGICAL             NEWOPT
      COMMON     /SOL7NP/ NEWOPT
      SAVE       /SOL7NP/
 
      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/
 
      EXTERNAL            MCHPAR
      CHARACTER*16        KEY
      CHARACTER*72        BUFFER
      LOGICAL             FIRST , PRNT
      SAVE                FIRST , NOUT  , PRNT
      DATA                FIRST /.TRUE./
 
*     If first time in, set NOUT.
*     NEWOPT is true first time into NPFILE or NPOPTN
*     and just after a call to an optimization routine.
*     PRNT is set to true whenever NEWOPT is true.
 
      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT   =  WMACH(11)
      END IF
      BUFFER = STRING
 
*     Call NPKEY to decode the option and set the parameter value.
*     If NEWOPT is true, reset PRNT and test specially for NOLIST.
 
      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT   = .TRUE.
         CALL NPKEY ( NOUT, BUFFER, KEY )
 
         IF (KEY .EQ. 'NOLIST') THEN
            PRNT   = .FALSE.
         ELSE
            WRITE (NOUT, '(// A / A /)')
     $         ' Calls to NPOPTN',
     $         ' ---------------'
            WRITE (NOUT, '( 6X, A )') BUFFER
         END IF
      ELSE
         IF (PRNT)
     $      WRITE (NOUT, '( 6X, A )') BUFFER
         CALL NPKEY ( NOUT, BUFFER, KEY )
 
         IF (KEY .EQ.   'LIST') PRNT = .TRUE.
         IF (KEY .EQ. 'NOLIST') PRNT = .FALSE.
      END IF
 
      RETURN
 
*     End of  NPOPTN.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPPRT ( KTCOND, CONVRG, LSUMRY, MSGNP, MSGQP,
     $                   NROWR, NROWT, N, NCLIN, NCNLN,
     $                   NCTOTL, NACTIV, LINACT, NLNACT, NZ, NFREE,
     $                   MAJITS, MINITS, ISTATE, ALFA, NFUN,
     $                   CONDHZ, CONDH, CONDT, OBJALF, OBJF,
     $                   GFNORM, GZNORM, CVNORM,
     $                   AX, C, R, T, VIOLN, X, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*4        LSUMRY
      LOGICAL            KTCOND(2), CONVRG
      INTEGER            ISTATE(NCTOTL)
      DOUBLE PRECISION   AX(*), C(*), R(NROWR,*), T(NROWT,*), VIOLN(*)
      DOUBLE PRECISION   X(N)
      DOUBLE PRECISION   WORK(N)
************************************************************************
*  NPPRT  prints various levels of output for NPCORE.
*
*           Msg        Cumulative result
*           ---        -----------------
*
*        le   0        no output.
*
*        eq   1        nothing now (but full output later).
*
*        eq   5        one terse line of output.
*
*        ge  10        same as 5 (but full output later).
*
*        ge  20        objective function,  x,  Ax  and  c.
*
*        ge  30        diagonals of  T  and  R.
*
*  Debug print is performed depending on the logical variable NPDBG.
*  NPDBG is set true when IDBG major iterations have been performed.
*  At this point,  printing is done according to a string of binary
*  digits of the form CLSVT (stored in the integer array INPDBG).
*
*  C  set 'on'  gives detailed information from the checking routines.
*  L  set 'on'  gives information from the linesearch.
*  S  set 'on'  gives information from the maximum step routine NPALF.
*  V  set 'on'  gives various vectors in  NPCORE  and its auxiliaries.
*  T  set 'on'  gives a trace of which routine was called and an
*               indication of the progress of the run.
*
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 66 version written November-1982.
*  This version of  NPPRT  dated  14-November-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      LOGICAL            PHEAD
      EXTERNAL           DNRM2
 
      IF (MSGNP .GE. 20) WRITE (NOUT, 1000) MAJITS
 
      IF (MSGNP  .GE. 5) THEN
*        ---------------------------------------------------------------
*        Print heading and terse line.
*        ---------------------------------------------------------------
         PHEAD = MSGQP .GT. 0  .OR.  MAJITS .EQ. 0
 
         IF (NCNLN .EQ. 0) THEN
            IF (PHEAD) WRITE (NOUT, 1100)
            WRITE (NOUT, 1300) MAJITS, MINITS, ALFA, NFUN, OBJALF,
     $                         N-NFREE, LINACT, NZ,
     $                         GFNORM, GZNORM, CONDH, CONDHZ, CONDT,
     $                         CONVRG, KTCOND(1), KTCOND(2), LSUMRY
 
         ELSE
            IF (PHEAD) WRITE (NOUT, 1110)
            WRITE (NOUT, 1310) MAJITS, MINITS, ALFA, NFUN, OBJALF,
     $                         N-NFREE, LINACT, NLNACT, NZ,
     $                         GFNORM, GZNORM, CONDH, CONDHZ, CONDT,
     $                         CVNORM, SCALE*RHONRM,
     $                         CONVRG, KTCOND(1), KTCOND(2), LSUMRY
         END IF
 
         IF (MSGNP .GE. 20) THEN
            IF (NCNLN .EQ. 0) THEN
               WRITE (NOUT, 1400) OBJF
            ELSE
               CVIOLS = DNRM2 ( NCNLN, VIOLN, 1 )
               WRITE (NOUT, 1410) OBJF, CVIOLS
            END IF
 
*           ------------------------------------------------------------
*           Print the constraint values.
*           ------------------------------------------------------------
            WRITE (NOUT, 2000)
            WRITE (NOUT, 2100) (X(J), ISTATE(J), J=1,N)
            IF (NCLIN .GT. 0)
     $         WRITE (NOUT, 2200) (AX(K), ISTATE(N+K),       K=1,NCLIN )
            IF (NCNLN .GT. 0)
     $         WRITE (NOUT, 2300) (C(K) , ISTATE(N+NCLIN+K), K=1,NCNLN )
 
            IF (MSGNP .GE. 30) THEN
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               INCT   = NROWT - 1
               IF (NACTIV .GT. 0) THEN
                  CALL DCOPY( NACTIV, T(NACTIV,NZ+1), INCT, WORK, 1 )
                  WRITE (NOUT, 3000) (WORK(J), J=1,NACTIV)
               END IF
               WRITE (NOUT, 3100) (R(J,J), J=1,N)
            END IF
         END IF
      END IF
 
      IF (MSGNP .GE. 20) WRITE (NOUT, 5000)
 
      LSUMRY(1:2) = '  '
      LSUMRY(4:4) = ' '
 
      RETURN
 
 1000 FORMAT(/// ' Major iteration', I5
     $       /   ' ====================' )
 1100 FORMAT(//  '  Itn', ' ItQP', '     Step',
     $           '  Nfun', '     Objective', ' Bnd', ' Lin', '  Nz',
     $           '  Norm Gf', '  Norm Gz', '  Cond H', ' Cond Hz',
     $           '  Cond T', ' Conv' )
 1110 FORMAT(//  '  Itn', ' ItQP', '     Step',
     $           '  Nfun', '         Merit', ' Bnd', ' Lin',
     $           ' Nln', '  Nz',
     $           '  Norm Gf', '  Norm Gz', '  Cond H', ' Cond Hz',
     $           '  Cond T' , '   Norm C', '  Penalty', ' Conv' )
 1300 FORMAT(2I5, 1PE9.1, I6, 1PE14.6, 3I4, 1P2E9.1, 1P3E8.0,
     $                        1X, L1, 1X, 2L1, A4 )
 1310 FORMAT(2I5, 1PE9.1, I6, 1PE14.6, 4I4, 1P2E9.1, 1P3E8.0,
     $            1P2E9.1,    1X, L1, 1X, 2L1, A4 )
 1400 FORMAT(/ ' Nonlinear objective value = ', 1PE15.6 )
 1410 FORMAT(/ ' Nonlinear objective value = ', 1PE15.6, '   Norm of',
     $         ' the nonlinear constraint violations = ', 1PE15.6 )
 2000 FORMAT(/ ' Values of the constraints and their predicted status'
     $       / ' ----------------------------------------------------')
 2100 FORMAT(/ ' Variables                  '/ (1X, 5(1PE15.6, I4)))
 2200 FORMAT(/ ' General linear constraints '/ (1X, 5(1PE15.6, I4)))
 2300 FORMAT(/ ' Nonlinear constraints      '/ (1X, 5(1PE15.6, I4)))
 3000 FORMAT(/ ' Diagonals of  T  =         '/       (1P5E15.6))
 3100 FORMAT(/ ' Diagonals of  R  =         '/       (1P5E15.6))
 5000 FORMAT(  ' ==================================================',
     $         '======================================='///)
 
*     End of  NPPRT .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPRSET( UNITQ,
     $                   N, NFREE, NZ, NQ, NROWR,
     $                   IPERM, KX,
     $                   GQ, R, ZY, WORK, QRWORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            IPERM(N), KX(N)
      DOUBLE PRECISION   GQ(N), R(NROWR,*), ZY(NQ,*)
      DOUBLE PRECISION   WORK(N), QRWORK(2*N)
 
************************************************************************
*  NPRSET  bounds the condition estimator of the transformed Hessian.
*  On exit, R is of the form
*               ( DRz   0     )
*               (  0  sigma*I )
*  where D is a diagonal matrix such that DRz has a bounded condition
*  number,  I is the identity matrix and sigma  is the geometric mean
*  of the largest and smallest elements of DRz. The QR factorization
*  with interchanges is used to give diagonals of DRz that are
*  decreasing in modulus.
*
*  Systems Optimization Laboratory, Stanford University.
*  This version of NPRSET dated  4-August-1986.
************************************************************************
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL6CM/ RCNDBD, RFROBN, DRMAX, DRMIN
 
      LOGICAL            NPDBG
      PARAMETER         (LDBG   = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      LOGICAL            OVERFL
      INTRINSIC          MAX   , MIN   , LOG   , REAL  , SQRT
      EXTERNAL           DDIV  , DDOT  , DNORM , DNRM2
      PARAMETER        ( ZERO   =0.0D+0, HALF =0.5D+0, ONE    =1.0D+0 )
 
*     ==================================================================
*     Bound the condition estimator of Q'HQ.
*     The scheme used here reduces the modulus of the larger
*     diagonals while increasing the modulus of the smaller ones.
*     ==================================================================
      IF (NZ .GT. 1) THEN
*        ---------------------------------------------------------------
*        Refactorize Rz.  Interchanges are used to give diagonals
*        of decreasing magnitude.
*        ---------------------------------------------------------------
         CALL DGEQRP( 'Column iterchanges', NZ, NZ, R, NROWR,
     $                WORK, IPERM, QRWORK, INFO )
 
         DO 110 J = 1, NZ
            JMAX = IPERM(J)
            IF (JMAX .GT. J) THEN
               IF (UNITQ) THEN
                  JSAVE    = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J)    = JSAVE
               ELSE
                  CALL DSWAP ( NFREE, ZY(1,JMAX), 1, ZY(1,J), 1 )
               END IF
 
               GJMAX    = GQ(JMAX)
               GQ(JMAX) = GQ(J)
               GQ(J)    = GJMAX
            END IF
  110    CONTINUE
      END IF
 
      IF (NZ .EQ. 0) THEN
         DRGM  = ONE
      ELSE
         COND  = DDIV  ( ABS(R(1,1)), ABS(R(NZ,NZ)), OVERFL )
 
         IF (COND .GT. RCNDBD) THEN
            IF (N .GT. 1) THEN
               PWR = LOG(RCNDBD)/LOG(COND) - ONE
               DO 120 K = 1, NZ
                  ROWSCL = ABS( R(K,K) )**PWR
                  CALL DSCAL ( NZ-K+1, ROWSCL, R(K,K), NROWR )
  120          CONTINUE
            END IF
         END IF
         DRGM  = HALF*SQRT(ABS( R(1,1)*R(NZ,NZ) ))
      END IF
 
*     Reset the range-space partition of the Hessian.
 
      IF (NZ .LT. N) THEN
         DO 130 J = NZ+1, N
            CALL DLOAD ( J, ZERO, R(1,J), 1 )
  130    CONTINUE
         CALL DLOAD ( N-NZ, DRGM, R(NZ+1,NZ+1), NROWR+1 )
      END IF
 
*     Recompute the Frobenius norm of R.
 
      SCLE  = SQRT(REAL(N - NZ))*DRGM
      SUMSQ = ONE
      DO 140 J = 1, NZ
         CALL DSSQ  ( J, R(1,J), 1, SCLE, SUMSQ )
  140 CONTINUE
      RFROBN = DNORM( SCLE, SUMSQ )
 
      RETURN
 
*     End of  NPRSET.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPSETX( UNITQ,
     $                   NCQP, NACTIV, NFREE, NZ,
     $                   N, NLNX, NCTOTL, NQ, NROWQP, NROWR, NROWT,
     $                   ISTATE, KACTIV, KX,
     $                   DXNORM, GDX,
     $                   AQP, ADX, BL, BU, RPQ, RPQ0, DX, GQ,
     $                   R, T, ZY, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      DOUBLE PRECISION   AQP(NROWQP,*), ADX(*), BL(NCTOTL), BU(NCTOTL),
     $                   RPQ(NLNX), RPQ0(NLNX), GQ(N), R(NROWR,*),
     $                   T(NROWT,*), ZY(NQ,*), DX(N), WORK(N)
************************************************************************
*  NPSETX   defines a point which lies on the initial working set for
*  the QP subproblem.  This routine is a similar to LSSETX except that
*  advantage is taken of the fact that the initial estimate of the
*  solution of the least-squares subproblem is zero.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  Level 2 BLAS added 12-June-1986.
*  This version of NPSETX dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      EXTERNAL           DDOT, DNRM2
      INTRINSIC          ABS , MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      NFIXED = N - NFREE
 
      GDX    = ZERO
      CALL DLOAD ( N   , ZERO, DX  , 1 )
      CALL DLOAD ( NLNX, ZERO, RPQ , 1 )
      CALL DLOAD ( NLNX, ZERO, RPQ0, 1 )
 
      IF (NACTIV + NFIXED .GT. 0) THEN
 
*        Set  work = residuals for constraints in the working set.
*        Solve for  dx,  the smallest correction to  x  that gives a
*        point on the constraints in the working set.
*        Set the fixed variables on their bounds,  solve the triangular
*        system  T*(dxy) = residuals,  and define  dx = Y*(dxy).
*        Use  (dxy)  to update  d(=Pr)  as  d = d - R'(  0  ).
*                                                     ( dxy )
 
         DO 100 I = 1, NFIXED
            J   = KX(NFREE+I)
            IF (ISTATE(J) .LE. 3) THEN
               BND   = BL(J)
               IF (ISTATE(J) .EQ. 2) BND = BU(J)
               DX(J) = BND
               WORK(NFREE+I) = BND
            ELSE
               WORK(NFREE+I) = ZERO
            END IF
  100    CONTINUE
 
         DO 110 I = 1, NACTIV
            K   = KACTIV(I)
            J   = N + K
            BND = BL(J)
            IF (ISTATE(J) .EQ. 2) BND = BU(J)
            WORK(NZ+I) = BND - DDOT  ( N, AQP(K,1), NROWQP, DX, 1 )
  110    CONTINUE
 
         IF (NACTIV .GT. 0)
     $      CALL CMTSOL( 1, NROWT, NACTIV, T(1,NZ+1), WORK(NZ+1) )
         CALL DCOPY ( NACTIV+NFIXED, WORK(NZ+1), 1, DX(NZ+1), 1 )
         IF (NZ .GT. 0)
     $      CALL DLOAD ( NZ, ZERO, DX, 1 )
 
         GDX  = DDOT  ( NACTIV+NFIXED, GQ(NZ+1), 1, DX(NZ+1), 1 )
 
         IF (NZ .LT. N) THEN
            CALL DGEMV ('N', NZ, N-NZ, -ONE, R(1,NZ+1), NROWR,
     $                  DX(NZ+1), 1, ONE, RPQ, 1 )
            IF (NZ .LT. NLNX) THEN
               NR  = NROWR
               IF (NZ+1 .EQ. N) NR = 1
               CALL DCOPY ( NLNX-NZ, DX(NZ+1), 1, RPQ(NZ+1), 1 )
               CALL DSCAL ( NLNX-NZ, (-ONE),      RPQ(NZ+1), 1 )
               CALL DTRMV ( 'U', 'N', 'N', NLNX-NZ, R(NZ+1,NZ+1), NR,
     $                      RPQ(NZ+1), 1 )
               IF (NLNX .LT. N) THEN
                  NR = NROWR
                  IF (NLNX+1 .EQ. N) NR = N - NZ
                  CALL DGEMV( 'N', NLNX-NZ, N-NLNX, -ONE,R(NZ+1,NLNX+1),
     $                        NR, DX(NLNX+1), 1, ONE, RPQ(NZ+1), 1 )
               END IF
            END IF
         END IF
 
         CALL CMQMUL( 2, N, NZ, NFREE, NQ, UNITQ, KX, DX, ZY, WORK )
      END IF
 
*     ------------------------------------------------------------------
*     Compute the 2-norm of  DX.
*     Initialize  A*DX.
*     ------------------------------------------------------------------
      DXNORM  = DNRM2 ( N, DX, 1 )
      IF (NCQP .GT. 0)
     $   CALL DGEMV ( 'N', NCQP, N, ONE, AQP, NROWQP, DX, 1, ZERO,ADX,1)
 
      IF (NPDBG  .AND.  INPDBG(2) .GT. 0)
     $   WRITE (NOUT, 1200) (DX(J), J = 1, N)
 
      RETURN
 
 1200 FORMAT(/ ' //NPSETX// Variables after NPSETX ... '/ (5G12.3))
 
*     End of  NPSETX.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPSRCH( NEEDFD, INFORM, N, NCNLN,
     $                   NROWJ, NROWUJ, NFUN, NGRAD,
     $                   NEEDC, CONFUN, OBJFUN,
     $                   ALFA, ALFBND, ALFMAX, ALFSML, DXNORM,
     $                   EPSRF, ETA, GDX, GRDALF, GLF1, GLF,
     $                   OBJF, OBJALF, QPCURV, XNORM,
     $                   C, CJAC, UJAC, CJDX, CMUL1, CMUL, CS1, CS,
     $                   DX, DLAM, DSLK, GRAD, UGRAD, QPMUL, RHO,
     $                   SLK1, SLK, X1, X, W, LENW )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            NEEDFD
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   DX(N), GRAD(N), UGRAD(N), X1(N), X(N)
      DOUBLE PRECISION   C(*), CJAC(NROWJ,*), UJAC(NROWUJ,*), CJDX(*),
     $                   CMUL1(*), CMUL(*), CS1(*), CS(*)
      DOUBLE PRECISION   DLAM(*), DSLK(*), QPMUL(*),
     $                   RHO(*), SLK1(*), SLK(*)
      DOUBLE PRECISION   W(LENW)
      EXTERNAL           OBJFUN, CONFUN
 
************************************************************************
*  NPSRCH finds the steplength ALFA that gives sufficient decrease in
*  the augmented Lagrangian merit function.
*
*  On exit,  if INFORM = 1, 2 or 3,  ALFA will be a nonzero steplength
*  with an associated merit function value  OBJALF  which is lower than
*  that at the base point. If  INFORM = 4, 5, 6 or 7,  ALFA  is zero
*  and  OBJALF will be the merit value at the base point.
*
*  Systems Optimization Laboratory, Stanford University, California.
*  Original version written  27-May-1985.
*  Level 2 BLAS added 12-June-1986.
*  This version of NPSRCH dated 12-July-1986.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      PARAMETER         (LENNP = 35)
      COMMON    /SOL1NP/ LOCNP(LENNP)
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      LOGICAL            NPDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      LOGICAL            DEBUG , DONE  , FIRST , IMPRVD
      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      PARAMETER        ( ZERO   =0.0D+0, HALF   =0.5D+0, ONE   =1.0D+0 )
      PARAMETER        ( TWO    =2.0D+0                                )
      PARAMETER        ( TOLG   =1.0D-1                                )
 
      EPSMCH = WMACH(3)
 
      LC     = LOCLS(14)
      LWORK  = LOCNP(12)
      LCJDX  = LOCNP(21)
 
      IF (.NOT. NEEDFD  .AND.  NCNLN .GT. 0)
     $   CS1JDX = DDOT( NCNLN, CS1, 1, CJDX, 1 )
 
*     ------------------------------------------------------------------
*     Set the input parameters and tolerances for SRCHC and SRCHQ.
*
*     TOLRX   is the tolerance on relative changes in DX resulting from
*             changes in ALFA.
*
*     TOLAX   is the tolerance on absolute changes in DX resulting from
*             changes in ALFA.
*
*     TOLABS  is the tolerance on absolute changes in ALFA.
*
*     TOLREL  is the tolerance on relative changes in ALFA.
*
*     TOLTNY  is the magnitude of the smallest allowable value of ALFA.
*             If  M(TOLABS) - M(0) .gt. EPSAF,  the linesearch tries
*             steps in the range  TOLTNY .LE. ALFA .LE. TOLABS.
*     ------------------------------------------------------------------
      NSTATE = 0
      DEBUG  = NPDBG  .AND.  INPDBG(4) .GT. 0
 
      EPSAF  = EPSRF*(ONE + ABS( OBJALF ))
 
      TOLAX  = EPSPT8
      TOLRX  = EPSPT8
 
      TOLABS = ALFMAX
      IF (TOLRX*XNORM + TOLAX .LT. DXNORM*ALFBND)
     $   TOLABS = (TOLRX*XNORM + TOLAX) /  DXNORM
      TOLREL = MAX( TOLRX , EPSMCH )
 
      T      = ZERO
      DO 10 J = 1, N
         S = ABS( DX(J) )
         Q = ABS( X(J) )*TOLRX + TOLAX
         IF (S .GT. T*Q) T = S / Q
   10 CONTINUE
 
      TOLTNY = TOLABS
      IF (T*TOLABS .GT. ONE) TOLTNY = ONE / T
 
      OLDF   = OBJALF
      OLDG   = GRDALF
 
      IF (NCNLN .GT. 0) CALL ILOAD ( NCNLN, (1), NEEDC, 1 )
 
      MODE  = 2
      IF (NEEDFD) MODE = 0
 
      FIRST  = .TRUE.
 
*     ------------------------------------------------------------------
*     Commence main loop, entering SRCHC or SRCHQ two or more times.
*     FIRST = true for the first entry,  false for subsequent entries.
*     DONE  = true indicates termination, in which case the value of
*     INFORM gives the result of the search.
*     ------------------------------------------------------------------
*+    REPEAT
  100    IF (NEEDFD) THEN
            CALL SRCHQ ( DEBUG, DONE, FIRST, IMPRVD, INFORM,
     $                   ALFMAX, ALFSML, EPSAF, ETA,
     $                   XTRY, FTRY,       OLDF, OLDG,
     $                   TOLABS, TOLREL, TOLTNY,
     $                   ALFA, ALFBST, FBEST        )
         ELSE
            CALL SRCHC ( DEBUG, DONE, FIRST, IMPRVD, INFORM,
     $                   ALFMAX,         EPSAF, ETA,
     $                   XTRY, FTRY, GTRY, OLDF, OLDG,
     $                   TOLABS, TOLREL, TOLTNY,
     $                   ALFA, ALFBST, FBEST, GBEST )
         END IF
 
         IF (IMPRVD) THEN
            OBJF   = TOBJ
            OBJALF = FTRY
 
            IF (NCNLN .GT. 0)
     $         CALL DCOPY ( NCNLN, W(LC), 1, C, 1 )
 
            IF (.NOT. NEEDFD) THEN
               CALL DCOPY ( N, UGRAD, 1, GRAD, 1 )
               GDX    = TGDX
               GLF    = TGLF
 
               IF (NCNLN .GT. 0) THEN
                  CALL DCOPY ( NCNLN, W(LCJDX), 1, CJDX, 1 )
                  DO 120  J = 1, N
                     CALL DCOPY ( NCNLN, UJAC(1,J), 1, CJAC(1,J), 1 )
  120             CONTINUE
               END IF
            END IF
         END IF
 
*        ---------------------------------------------------------------
*        If DONE = FALSE,  the problem functions must be computed for
*        the next entry to SRCHC or SRCHQ.
*        If DONE = TRUE,   this is the last time through.
*        ---------------------------------------------------------------
         IF (.NOT. DONE) THEN
 
            NFUN  = NFUN  + 1
            IF (.NOT. NEEDFD) NGRAD = NGRAD + 1
 
            CALL DCOPY ( N,       X1, 1, X, 1 )
            CALL DAXPY ( N, ALFA, DX, 1, X, 1 )
            IF (NCNLN .GT. 0) THEN
 
*              Compute the new estimates of the multipliers and slacks.
*              If the step length is greater than one,  the multipliers
*              are fixed as the QP-multipliers.
 
               IF (ALFA .LE. ONE) THEN
                  CALL DCOPY ( NCNLN,       CMUL1, 1, CMUL, 1 )
                  CALL DAXPY ( NCNLN, ALFA, DLAM , 1, CMUL, 1 )
               END IF
               CALL DCOPY ( NCNLN,       SLK1, 1, SLK, 1 )
               CALL DAXPY ( NCNLN, ALFA, DSLK, 1, SLK, 1 )
 
*              ---------------------------------------------------------
*              Compute the new constraint vector and Jacobian.
*              ---------------------------------------------------------
               CALL CONFUN( MODE, NCNLN, N, NROWUJ,
     $                      NEEDC, X, W(LC), UJAC, NSTATE )
               IF (MODE .LT. 0) GO TO 999
 
               CALL DCOPY ( NCNLN,         W(LC), 1, CS, 1 )
               CALL DAXPY ( NCNLN, (-ONE), SLK  , 1, CS, 1 )
 
               CALL DCOPY ( NCNLN, CS , 1, W(LWORK), 1 )
               CALL DDSCL ( NCNLN, RHO, 1, W(LWORK), 1 )
 
               FTERM  =            DDOT( NCNLN, CMUL    , 1, CS, 1 ) -
     $                  HALF*SCALE*DDOT( NCNLN, W(LWORK), 1, CS, 1 )
 
            END IF
 
*           ------------------------------------------------------------
*           Compute the value and gradient of the objective function.
*           ------------------------------------------------------------
            CALL OBJFUN( MODE, N, X, TOBJ, UGRAD, NSTATE )
            IF (MODE .LT. 0) GO TO 999
 
            FTRY   = TOBJ
            IF (NCNLN .GT. 0) FTRY = TOBJ  - FTERM
 
*           ------------------------------------------------------------
*           Compute auxiliary gradient information.
*           ------------------------------------------------------------
            IF (.NOT. NEEDFD) THEN
               GTRY   = DDOT( N, UGRAD, 1, DX, 1 )
               TGDX   = GTRY
               TGLF   = GTRY
               IF (NCNLN .GT. 0) THEN
 
*                 Compute the Jacobian times the search direction.
 
                  CALL DGEMV ( 'N', NCNLN, N, ONE, UJAC, NROWUJ, DX, 1,
     $                         ZERO, W(LCJDX), 1 )
 
                  CALL DCOPY ( NCNLN,         W(LCJDX), 1, W(LWORK), 1 )
                  CALL DAXPY ( NCNLN, (-ONE), DSLK    , 1, W(LWORK), 1 )
 
                  GTRY   = GTRY - DDOT( NCNLN, CMUL, 1, W(LWORK), 1 )
                  IF (ALFA .LE. ONE)
     $               GTRY   = GTRY - DDOT( NCNLN, DLAM, 1, CS      , 1 )
 
                  CALL DDSCL ( NCNLN, RHO , 1, W(LWORK), 1 )
                  GTRY = GTRY  +
     $                     SCALE*DDOT( NCNLN, W(LWORK), 1, CS   , 1 )
 
                  TGLF = TGDX  - DDOT( NCNLN, W(LCJDX), 1, QPMUL, 1 )
 
*                 ------------------------------------------------------
*                 If ALFBND .LE. ALFA .LT. ALFMAX and the norm of the
*                 quasi-Newton update is bounded, set ALFMAX to be ALFA.
*                 This will cause the linesearch to stop if the merit
*                 function is decreasing at the boundary.
*                 ------------------------------------------------------
                  IF (ALFBND .LE. ALFA  .AND.  ALFA .LT. ALFMAX) THEN
 
                     CSJDX  = DDOT   ( NCNLN, CS, 1, W(LCJDX), 1 )
 
                     IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $                  WRITE (NOUT, 1400) CSJDX, CS1JDX, CURVLF
 
                     CURVLF = TGLF  - GLF1
                     CURVC  = ABS( CSJDX - CS1JDX )
                     RHOBFS = MAX( QPCURV*TOLG - CURVLF, ZERO )
                     IF (RHOBFS .LE. CURVC*RHOMAX) THEN
                        ALFMAX = ALFA
                     ELSE
                        ALFBND = MIN( TWO*ALFA, ALFMAX )
                     END IF
                     IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $                  WRITE(NOUT,1300) ALFBND, ALFA, ALFMAX
                  END IF
               END IF
            END IF
         END IF
*+    UNTIL (      DONE)
      IF    (.NOT. DONE) GO TO 100
 
      ALFA = ALFBST
      IF (.NOT. IMPRVD) THEN
         CALL DCOPY ( N,       X1, 1, X, 1 )
         CALL DAXPY ( N, ALFA, DX, 1, X, 1 )
         IF (NCNLN .GT. 0) THEN
            IF (ALFA .LE. ONE) THEN
               CALL DCOPY ( NCNLN,       CMUL1, 1, CMUL, 1 )
               CALL DAXPY ( NCNLN, ALFA, DLAM , 1, CMUL, 1 )
            END IF
            CALL DCOPY ( NCNLN,         SLK1 , 1, SLK, 1 )
            CALL DAXPY ( NCNLN,   ALFA, DSLK , 1, SLK, 1 )
            CALL DCOPY ( NCNLN,         C    , 1, CS , 1 )
            CALL DAXPY ( NCNLN, (-ONE), SLK  , 1, CS , 1 )
         END IF
      END IF
 
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1200) INFORM
 
      RETURN
 
*     The user wants to stop.
 
  999 INFORM = MODE
      RETURN
 
 1200 FORMAT(/ ' //NPSRCH// INFORM  = ', I4 )
 1300 FORMAT(/ ' //NPSRCH//        ALFBND          ALFA        ALFMAX'
     $       / ' //NPSRCH//', 1P3E14.2 )
 1400 FORMAT(/ ' //NPSRCH//         CSJDX        CS1JDX        CURVLF'
     $       / ' //NPSRCH//', 1P3E14.2 )
 
*     End of  NPSRCH.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPUPDT( LSUMRY, UNITQ,
     $                   N, NCNLN, NFREE, NZ,
     $                   NROWJ1, NROWJ2, NQ, NROWR, KX,
     $                   ALFA, GLF1, GLF2, QPCURV,
     $                   CJAC1, CJAC2, CJDX1, CJDX2,
     $                   CS1, CS2, GQ1, GQ2, HPQ, RPQ,
     $                   QPMUL, R, OMEGA, ZY, WRK1, WRK2 )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*4        LSUMRY
      LOGICAL            UNITQ
      INTEGER            KX(N)
      DOUBLE PRECISION   CJAC1(NROWJ1,*), CJAC2(NROWJ2,*),
     $                   CJDX1(*), CJDX2(*), CS1(*), CS2(*),
     $                   GQ1(N), GQ2(N), HPQ(N), RPQ(N), QPMUL(*),
     $                   R(NROWR,*), OMEGA(*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N+NCNLN), WRK2(N)
 
************************************************************************
*  NPUPDT  computes the BFGS update for the approximate Hessian of the
*  Lagrangian.  If the approximate curvature of the Lagrangian function
*  is negative,  a nonnegative penalty vector OMEGA(i) of minimum two
*  norm is computed such that the approximate curvature of the augmented
*  Lagrangian will be positive. If no finite penalty vector exists,  the
*  BFGS update is performed with the approximate curvature modified to
*  be a small positive value.
*
*  On entry,  GQ1 and GQ2 contain the transformed objective gradients at
*  X1 and X2,  HPQ contains  R'R(pq), the transformed Hessian times the
*  transformed search direction.  The vectors GQ1 and HPQ are not saved.
*  If the regular BFGS quasi-Newton update could not be performed, the
*  first character of LSUMRY is loaded with 'M'.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 66 version written April 1984.
*  Level 2 BLAS added 12-June-1986.
*  This version of NPUPTD dated  4-August-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL6CM/ RCNDBD, RFROBN, DRMAX, DRMIN
 
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      LOGICAL            NPDBG
      PARAMETER        ( LDBG   = 5 )
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
 
      LOGICAL            OVERFL, SSBFGS
      INTRINSIC          MAX   , MIN   , SQRT
      EXTERNAL           IDAMAX, DDIV  , DDOT  , DNRM2
      PARAMETER        ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
      PARAMETER        ( TOLG   = 1.0D-1                  )
 
      IF (NCNLN .GT. 0) CALL DLOAD ( NCNLN, ZERO, OMEGA, 1 )
 
*     ------------------------------------------------------------------
*     Set CURVL = (G2 - G1)'DX,  the approximate curvature along DX of
*     the (augmented) Lagrangian.  At first, the curvature is not scaled
*     by the steplength ALFA.
*     ------------------------------------------------------------------
      CURVL  = GLF2 -   GLF1
      TINYCL =        QPCURV * TOLG
      SSBFGS = CURVL .LE. ALFA*TINYCL
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) SSBFGS, TINYCL, CURVL
 
*     ------------------------------------------------------------------
*     Test if CURVL is sufficiently positive.  If there are no nonlinear
*     constraints,  no update can be performed.
*     ------------------------------------------------------------------
      IF (CURVL  .LT. TINYCL) THEN
         LSUMRY(1:1) = 'Modified BFGS'
         IF (NCNLN .GT. 0) THEN
            QMAX = ZERO
            DO 200 I = 1, NCNLN
               QI  = CJDX2(I)*CS2(I) - CJDX1(I)*CS1(I)
               QMAX = MAX( QMAX, QI )
               IF (QI .LE. ZERO) WRK1(I) = ZERO
               IF (QI .GT. ZERO) WRK1(I) = QI
  200       CONTINUE
 
            QNORM = DNRM2 ( NCNLN, WRK1, 1 )
 
            TEST  = MAX( TINYCL - CURVL, ZERO )
            BETA  = DDIV  ( QMAX*TEST, QNORM*QNORM, OVERFL )
            IF (BETA .LT. RHOMAX  .AND.  .NOT. OVERFL) THEN
               LSUMRY(1:1) = ' '
               BETA  = TEST/(QNORM*QNORM)
               DO 210 I = 1, NCNLN
                  QI       = WRK1(I)
                  OMEGA(I) =            BETA*QI
                  CURVL    = CURVL    + BETA*QI*QI
  210          CONTINUE
 
               IF (NPDBG) THEN
                  IMAX = IDAMAX( NCNLN, OMEGA, 1 )
                  IF (INPDBG(1) .GT. 0)
     $               WRITE (NOUT, 1250) OMEGA(IMAX)
 
                  IF (INPDBG(2) .GT. 0)
     $               WRITE (NOUT, 1300) (OMEGA(I), I=1,NCNLN)
               END IF
            END IF
         END IF
      END IF
 
*     ------------------------------------------------------------------
*     Compute the difference in the augmented Lagrangian gradient.
*     ------------------------------------------------------------------
*     Update GQ1 to include the augmented Lagrangian terms.
 
      IF (NCNLN .GT. 0) THEN
 
         DO 310 I = 1, NCNLN
            WRK1(I) = - QPMUL(I) + OMEGA(I) * CS1(I)
  310    CONTINUE
         CALL DGEMV ( 'T', NCNLN, N, ONE, CJAC1, NROWJ1, WRK1, 1,
     $                ZERO, WRK2, 1 )
 
         DO 320 I = 1, NCNLN
            WRK1(I) =   QPMUL(I) - OMEGA(I) * CS2(I)
  320    CONTINUE
         CALL DGEMV ( 'T', NCNLN, N, ONE, CJAC2, NROWJ2, WRK1, 1,
     $                ONE, WRK2, 1 )
 
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ, KX, WRK2, ZY, WRK1 )
         CALL DAXPY ( N, ONE, WRK2, 1, GQ1, 1 )
      END IF
 
      IF (NPDBG  .AND.  INPDBG(1) .GT. 0)
     $   WRITE (NOUT, 1100) ALFA  , CURVL
 
      IF (CURVL .LT. TINYCL) CURVL  = TINYCL
 
      DO 330 J = 1, N
         WRK2(J) = GQ2(J) - GQ1(J)
  330 CONTINUE
 
      RTGTP  = SQRT(QPCURV)
      RTYTS  = SQRT(ALFA*CURVL)
      ETA    = ONE
      IF (SSBFGS)
     $   ETA = RTYTS / (RTGTP*ALFA)
 
      TRACE1 = DNRM2 ( N,  HPQ, 1 ) /  RTGTP
      TRACE2 = DNRM2 ( N, WRK2, 1 ) / (RTYTS*ETA)
      RFROBN = ETA*SQRT( ABS(  (RFROBN - TRACE1)*(RFROBN + TRACE1)
     $                                 + TRACE2**2) )
 
*     ==================================================================
*     Update the Cholesky factor of  Q'HQ.
*     ==================================================================
*     Normalize the vector  RPQ ( = R(pq) ).
 
      CALL DSCAL ( N, (ONE / RTGTP), RPQ, 1 )
 
*     Do the self-scaled or regular BFGS update.
*     Form the vector WRK1 = gamma * (GQ2 - GQ1) - beta * R'R*PQ,
*     where  gamma = 1/SQRT( CURV ) = 1/SQRT( (GQ2 - GQ1)'SQ )
 
      CALL DSCAL ( N, (ONE / RTGTP), HPQ, 1 )
 
      IF (SSBFGS) THEN
         DO 410 J   = 1, N
            CALL DSCAL ( J, ETA, R(1,J), 1 )
            WRK1(J) = WRK2(J)/RTYTS  -  ETA * HPQ(J)
  410    CONTINUE
      ELSE
         DO 420 J   = 1, N
            WRK1(J) = WRK2(J)/RTYTS  -        HPQ(J)
  420    CONTINUE
      END IF
 
*     Perform the update to  R = R + RPQ*WRK1'.
*     RPQ is overwritten and HPQ is used as workspace.
 
      CALL CMR1MD( N, 0, N, NROWR, N, N, R, HPQ, RPQ, WRK1 )
 
      RETURN
 
 1000 FORMAT(/ ' //NPUPDT// SSBFGS    min. CURVL         CURVL '
     $       / ' //NPUPDT//   ', L4, 1P2E14.2 )
 1100 FORMAT(/ ' //NPUPDT//          ALFA         CURVL '
     $       / ' //NPUPDT//', 1P2E14.2 )
 1250 FORMAT(/ ' //NPUPDT//   OMEGA(IMAX)'
     $       / ' //NPUPDT//', 1PE14.2 )
 1300 FORMAT(/ ' //NPUPDT//  Penalty parameters = '  / (1P5E15.6))
 
*     End of  NPUPDT.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE NPSOL ( N, NCLIN, NCNLN, NROWA, NROWUJ, NROWR,
     $                   A, BL, BU,
     $                   CONFUN, OBJFUN,
     $                   INFORM, ITER, ISTATE,
     $                   C, UJAC, CLAMDA, OBJF, UGRAD, R, X,
     $                   IW, LENIW, W, LENW )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL           CONFUN, OBJFUN
      INTEGER            ISTATE(N+NCLIN+NCNLN)
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   A(NROWA,*), BL(N+NCLIN+NCNLN),
     $                   BU(N+NCLIN+NCNLN)
      DOUBLE PRECISION   C(*), UJAC(NROWUJ,*), CLAMDA(N+NCLIN+NCNLN)
      DOUBLE PRECISION   UGRAD(N), R(NROWR,*), X(N)
      DOUBLE PRECISION   W(LENW)
 
*-----------------------------------------------------------------------
*
*  NPSOL   solves the nonlinear programming problem
*
*            minimize                   F(x)
*
*                                    (    x  )
*            subject to    bl  .le.  (  A*x  )  .le.  bu
*                                    (  c(x) )
*
*  where  F(x)  is a smooth scalar function,  A  is a constant matrix
*  and  c(x)  is a vector of smooth nonlinear functions.  The feasible
*  region is defined by a mixture of linear and nonlinear equality or
*  inequality constraints on  x.
*
*  The dimensions of the problem are...
*
*  N        the number of variables (dimension of  x),
*
*  NCLIN    the number of linear constraints (rows of the matrix  A),
*
*  NCNLN    the number of nonlinear constraints (dimension of  c(x)),
*
*
*  NPSOL   uses a sequential quadratic programming algorithm, with a
*  positive-definite quasi-Newton approximation to the transformed
*  Hessian  Q'HQ  of the Lagrangian function (which will be stored in
*  the array  R).
*
*
*  Complete documentation for  NPSOL  is contained in Report
*  SOL 86-2, Users guide for NPSOL (Version 4.0), by P.E. Gill,
*  W. Murray, M.A. Saunders and M.H. Wright, Department of Operations
*  Research,  Stanford University, Stanford, California 94305.
*
*  Systems Optimization Laboratory, Stanford University.
*  Version 1.1,  April     12, 1983. (The less said about this one.....)
*  Version 2.0,  April     30, 1984.
*  Version 3.0,  March     20, 1985. (First Fortran 77 version).
*  Version 3.2,  August    20, 1985.
*  Version 4.0,  April     16, 1986. (First version with differences).
*  Version 4.01, June      30, 1986. (Level 2 Blas + F77 linesearch).
*  Version 4.02, August     5, 1986. (Reset SSBFGS. One call to CHFD).
*
*  Copyright  1983  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR Clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under Grants MCS-7926009 and ECS-8312142; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; the Army Research Office Contract DAA29-84-K-0156;
*  and the Office of Naval Research Grant N00014-75-C-0267.
*  ---------------------------------------------------------------------
 
*  Common blocks.
 
*-----------------------------------------------------------------------
      PARAMETER         (MXPARM = 30)
      INTEGER            IPRMLS(MXPARM), IPSVLS
      DOUBLE PRECISION   RPRMLS(MXPARM), RPSVLS
 
      COMMON    /LSPAR1/ IPSVLS(MXPARM),
     $                   IDBGLS, ITMAX1, ITMAX2, LCRASH, LDBGLS, LPROB ,
     $                   MSGLS , NN    , NNCLIN, NPROB , IPADLS(20)
 
      COMMON    /LSPAR2/ RPSVLS(MXPARM),
     $                   BIGBND, BIGDX , BNDLOW, BNDUPP, TOLACT, TOLFEA,
     $                   TOLRNK, RPADLS(23)
 
      EQUIVALENCE       (IPRMLS(1), IDBGLS), (RPRMLS(1), BIGBND)
 
      SAVE      /LSPAR1/, /LSPAR2/
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      INTEGER            IPRMNP(MXPARM), IPSVNP
      DOUBLE PRECISION   RPRMNP(MXPARM), RPSVNP
 
      COMMON    /NPPAR1/ IPSVNP(MXPARM),
     $                   IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     $                   LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     $                   NLNJ  , NLNX  , NNCNLN, IPADNP(15)
 
      COMMON    /NPPAR2/ RPSVNP(MXPARM),
     $                   CDINT , CTOL  , EPSRF , ETA   , FDINT , FTOL  ,
     $                   RPADNP(24)
 
      EQUIVALENCE       (IPRMNP(1), IDBGNP), (RPRMNP(1), CDINT)
 
      SAVE      /NPPAR1/, /NPPAR2/
*-----------------------------------------------------------------------
      EQUIVALENCE  (IDBGNP, IDBG  ), (ITMXNP, NMAJOR), (ITMAX2, NMINOR)
      EQUIVALENCE  (LDBGLS, MNRDBG), (LDBGNP, MJRDBG), (MSGLS , MSGQP )
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT , NCOLT , NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE , DTMAX , DTMIN
      COMMON    /SOL6CM/ RCNDBD, RFROBN, DRMAX , DRMIN
 
      LOGICAL            UNITQ
      COMMON    /SOL1SV/ NACTIV, NFREE , NZ   , UNITQ
      SAVE      /SOL1SV/
 
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      PARAMETER         (LENNP = 35)
      COMMON    /SOL1NP/ LOCNP(LENNP)
      COMMON    /SOL4NP/ LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON    /SOL5NP/ LVRFYC, JVERFY(4)
      LOGICAL            INCRUN
      COMMON    /SOL6NP/ RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
 
      LOGICAL            CMDBG, LSDBG, NPDBG
      PARAMETER         (LDBG = 5)
      COMMON    /NPDEBG/ INPDBG(LDBG), NPDBG
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG
 
      INTRINSIC          ABS   , MAX   , MIN   , MOD   , SQRT  , REAL
 
*     Local variables.
 
      EXTERNAL           DDIV  , DDOT  , DNORM , DNRM2
      CHARACTER*8        NAMES(1)
      LOGICAL            COLD  , LINOBJ, NAMED , OVERFL, ROWERR, VERTEX
      PARAMETER         (ZERO   =0.0D+0, POINT1 =0.1D+0, POINT3 =3.3D-1)
      PARAMETER         (POINT8 =0.8D+0, POINT9 =0.9D+0, ONE    =1.0D+0)
      PARAMETER         (GROWTH =1.0D+2                                )
 
      CHARACTER*40       TITLE
      DATA               TITLE
     $                 / 'SOL/NPSOL  ---  Version 4.02   Aug  1986' /
 
*     Set the machine-dependent constants.
 
      CALL MCHPAR()
 
      EPSMCH = WMACH( 3)
      RTEPS  = WMACH( 4)
      NOUT   = WMACH(11)
 
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
 
      RHOMAX = ONE/EPSMCH
      ROOTN  = SQRT(REAL(N))
 
*     Default names will be provided for variables during printing.
 
      NAMED  = .FALSE.
      INFORM = 0
      ITER   = 0
 
*     Set the default values for the parameters.
 
      CALL NPDFLT( N, NCLIN, NCNLN, LENIW, LENW, TITLE )
 
      COLD   = LCRASH .EQ. 0
 
      NPLIN  = N     + NCLIN
      NCTOTL = NPLIN + NCNLN
 
*     Assign the dimensions of arrays in the parameter list of NPCORE.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter MINACT and MINFXD
*     accordingly.
 
      MINACT = 0
      MINFXD = 0
 
      MXFREE = N - MINFXD
      MAXACT = MAX( 1, MIN( N, NCLIN ) )
      MAXNZ  = N - ( MINFXD + MINACT )
 
      IF (NCLIN + NCNLN .EQ. 0) THEN
         NQ    = 1
         NROWT = 1
         NCOLT = 1
      ELSE
         NQ    = MAX( 1, MXFREE )
         NROWT = MAX( MAXNZ, MAXACT )
         NCOLT = MXFREE
      END IF
 
      LENNAM = 1
 
      NROWQP = MAX( NCLIN+NCNLN, 1 )
      IF (NCNLN .EQ. 0  .AND.  NCLIN .GT. 0) NROWQP = NROWA
 
*     NPLOC  defines the arrays that contain the locations of various
*     work arrays within  W  and  IW.
 
      LITOTL = 0
      LWTOTL = 0
      CALL NPLOC( N, NCLIN, NCNLN, NCTOTL, LITOTL, LWTOTL)
 
*     Allocate certain addresses that are not allocated in NPLOC.
 
      LAX    = LWTOTL + 1
      LWTOTL = LAX    + NCLIN - 1
      LAX    = MIN( LAX, LWTOTL )
 
*     Check input parameters and storage limits.
 
      CALL CMCHK ( NERROR, MSGNP, COLD, .FALSE.,
     $             LENIW, LENW, LITOTL, LWTOTL,
     $             N, NCLIN, NCNLN,
     $             ISTATE, IW, NAMED, NAMES, LENNAM,
     $             BL, BU, X )
 
      IF (NERROR .GT. 0) THEN
         INFORM = 9
         GO TO 800
      END IF
 
      LKACTV = LOCLS( 1)
      LANORM = LOCLS( 2)
      LCJDX  = LOCLS( 3)
      LRES   = LOCLS( 5)
      LRES0  = LOCLS( 6)
      LGQ    = LOCLS( 9)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK1  = LOCLS(14)
 
      LKX    = LOCNP( 1)
      LIPERM = LOCNP( 2)
      LAQP   = LOCNP( 3)
      LDX    = LOCNP( 7)
      LFEATL = LOCNP(10)
      LWRK2  = LOCNP(12)
 
      LCMUL  = LOCNP(16)
      LWRK3  = LOCNP(21)
      LNEEDC = LOCNP(24)
      LHFRWD = LOCNP(25)
      LHCTRL = LOCNP(26)
      LCJAC  = LOCNP(27)
      LGRAD  = LOCNP(28)
 
      NROWJ  = MAX ( NCNLN, 1 )
 
      TOLRNK = ZERO
      RCNDBD = ONE/SQRT(EPSPT5)
 
      IF (TOLFEA .GT. ZERO)
     $   CALL DLOAD ( NPLIN, TOLFEA, W(LFEATL), 1 )
 
      IF (NCNLN .GT. 0  .AND.  CTOL .GT. ZERO)
     $   CALL DLOAD ( NCNLN, CTOL, W(LFEATL+NPLIN), 1 )
 
      IF (LFDSET .EQ. 0) THEN
         FDCHK = SQRT( EPSRF )
      ELSE IF (LFDSET .EQ. 1) THEN
         FDCHK = FDINT
      ELSE
         FDCHK = W(LHFRWD)
      END IF
 
      NFUN   = 0
      NGRAD  = 0
      NSTATE = 1
 
*     ------------------------------------------------------------------
*     If required,  compute the problem functions.
*     If the constraints are nonlinear,  the first call of CONFUN
*     sets up any constant elements in the Jacobian matrix.  A copy of
*     the Jacobian (with constant elements set) is placed in  UJAC.
*     ------------------------------------------------------------------
      IF (LVERFY .GE. 10) THEN
         XNORM  = DNRM2 ( N, X, 1 )
         LVRFYC = LVERFY - 10
 
         CALL NPCHKD( INFO, MSGNP, NSTATE, LVLDER, NFUN, NGRAD,
     $                NROWJ, NROWUJ, N, NCNLN,
     $                CONFUN, OBJFUN, IW(LNEEDC),
     $                BIGBND, EPSRF, CDINT, FDINT,
     $                FDCHK, FDNORM, OBJF, XNORM,
     $                BL, BU, C, W(LWRK3), W(LCJAC), UJAC, W(LCJDX),
     $                W(LDX), W(LGRAD), UGRAD, W(LHFRWD), W(LHCTRL),
     $                X, W(LWRK1), W(LWRK2), W, LENW )
 
         IF (INFO .NE. 0) THEN
            IF (INFO .GT. 0) INFORM = 7
            IF (INFO .LT. 0) INFORM = INFO
            GO TO 800
         END IF
         NSTATE = 0
      END IF
 
      IF (LCRASH .LT. 2) THEN
*        ===============================================================
*        Cold or warm start.  Use  LSCORE  to obtain a point that
*        satisfies the linear constraints.
*        ===============================================================
         CALL ICOPY ( LDBG, ILSDBG, 1, ICMDBG, 1 )
 
         IF (NCLIN .GT. 0) THEN
            IANRMJ = LANORM
            DO 110 J = 1, NCLIN
               W(IANRMJ) = DNRM2 ( N, A(J,1), NROWA )
               IANRMJ    = IANRMJ + 1
  110       CONTINUE
            CALL DCOND ( NCLIN, W(LANORM), 1, ASIZE, AMIN )
         END IF
 
         CALL DCOND ( NPLIN, W(LFEATL), 1, FEAMAX, FEAMIN )
         CALL DCOPY ( NPLIN, W(LFEATL), 1, W(LWTINF), 1 )
         CALL DSCAL ( NPLIN, (ONE/FEAMIN), W(LWTINF), 1 )
 
*        ===============================================================
*        The input values of X and (optionally)  ISTATE are used by
*        LSCRSH  to define an initial working set.
*        ===============================================================
         VERTEX = .FALSE.
         CALL LSCRSH( COLD, VERTEX,
     $                NCLIN, NPLIN, NACTIV, NARTIF,
     $                NFREE, N, NROWA,
     $                ISTATE, IW(LKACTV),
     $                BIGBND, TOLACT,
     $                A, W(LAX), BL, BU, X, W(LWRK1), W(LWRK2) )
 
         UNITQ  = .TRUE.
         NRES   = 0
         NGQ    = 0
         CONDMX = ONE / EPSPT5
 
         IKX    = LKX
         DO 120 I = 1, N
            IW(IKX) = I
            IKX     = IKX + 1
  120    CONTINUE
 
         IF (COLD) THEN
            NRANK  = 0
         ELSE
            NRANK  = NLNX
            CALL DLOAD ( NLNX, (ZERO), W(LRES0), 1 )
         END IF
 
*        ---------------------------------------------------------------
*        Re-order KX so that the free variables come first.
*        If a warm start is required, NRANK will be nonzero and the
*        factor R will be updated.
*        ---------------------------------------------------------------
         CALL LSBNDS( UNITQ,
     $                INFORM, NZ, NFREE, NRANK, NRES, NGQ,
     $                N, NQ, NROWA, NROWR, NROWT,
     $                ISTATE, IW(LKX),
     $                CONDMX,
     $                A, R, W(LT), W(LRES0), W(LGQ),
     $                W(LZY), W(LWRK1), W(LWRK2) )
 
*        ---------------------------------------------------------------
*        Factorize the initial working set.
*        ---------------------------------------------------------------
         IF (NACTIV .GT. 0) THEN
            NACT1  = NACTIV
            NACTIV = 0
 
            CALL LSADDS( UNITQ, VERTEX,
     $                   INFORM, 1, NACT1, NACTIV, NARTIF, NZ, NFREE,
     $                   NRANK, NREJTD, NRES, NGQ,
     $                   N, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, IW(LKACTV), IW(LKX),
     $                   CONDMX,
     $                   A, R, W(LT), W(LRES0), W(LGQ),
     $                   W(LZY), W(LWRK1), W(LWRK2) )
         END IF
 
         SSQ1 = ZERO
 
         LINOBJ = .FALSE.
         CALL LSSETX( LINOBJ, ROWERR, UNITQ,
     $                NCLIN, NACTIV, NFREE, NRANK, NZ,
     $                N, NPLIN, NQ, NROWA, NROWR, NROWT,
     $                ISTATE, IW(LKACTV), IW(LKX),
     $                JMAX, ERRMAX, CTX, XNORM,
     $                A, W(LAX), BL, BU, W(LGQ), W(LRES), W(LRES0),
     $                W(LFEATL), R, W(LT), X, W(LZY),W(LWRK1),W(LWRK2) )
 
*        ---------------------------------------------------------------
*        Call  LSCORE  to find a feasible  x.
*        ---------------------------------------------------------------
*        Use  WORK2  as the multiplier vector.
 
         JINF   = 0
         LCLAM  = LWRK2
 
         IDBGSV = IDBG
         IF (IDBG .GT. 0) THEN
            IDBG = NMINOR + 1
         END IF
 
         CALL LSCORE( 'FP problem', NAMED, NAMES, LINOBJ, UNITQ,
     $                NLPERR, ITER, JINF, NCLIN, NPLIN,
     $                NACTIV, NFREE, NRANK, NZ, NZ1,
     $                N, NROWA, NROWR,
     $                ISTATE, IW(LKACTV), IW(LKX),
     $                CTX, OBJ, SSQ1, SUMINF, NUMINF, XNORM,
     $                BL, BU, A, W(LCLAM), W(LAX),
     $                W(LFEATL), R, X, IW, W )
 
         IF (NLPERR .GT. 0) THEN
            INFORM = 2
            GO TO 800
         END IF
      END IF
 
      IDBG  = IDBGSV
      CALL ICOPY ( LDBG, INPDBG, 1, ICMDBG, 1 )
 
      LVRFYC = LVERFY
      IF (LVERFY .GE. 10) LVRFYC = -1
 
      CALL NPCHKD( INFO, MSGNP, NSTATE, LVLDER, NFUN, NGRAD,
     $             NROWJ, NROWUJ, N, NCNLN,
     $             CONFUN, OBJFUN, IW(LNEEDC),
     $             BIGBND, EPSRF, CDINT, FDINT,
     $             FDCHK, FDNORM, OBJF, XNORM,
     $             BL, BU, C, W(LWRK3), W(LCJAC), UJAC, W(LCJDX),
     $             W(LDX), W(LGRAD), UGRAD, W(LHFRWD), W(LHCTRL),
     $             X, W(LWRK1), W(LWRK2), W, LENW )
 
      IF (INFO .NE. 0) THEN
         IF (INFO .GT. 0) INFORM = 7
         IF (INFO .LT. 0) INFORM = INFO
         GO TO 800
      END IF
 
      CALL DCOPY ( N, W(LGRAD), 1, W(LGQ), 1 )
      CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ,
     $             IW(LKX), W(LGQ), W(LZY), W(LWRK1) )
 
      IF (COLD) THEN
*        ---------------------------------------------------------------
*        Cold start.  Initialize  R  as the identity matrix.
*        ---------------------------------------------------------------
         DO 210 J = 1, N
            CALL DLOAD ( N, ZERO, R(1,J), 1 )
  210    CONTINUE
         CALL DLOAD ( N, ONE, R, NROWR+1 )
         RFROBN = ROOTN
 
         IF (NCNLN .GT. 0) CALL DLOAD ( NCNLN, (ZERO), W(LCMUL), 1 )
      ELSE
*        ---------------------------------------------------------------
*        Warm start.
*        Set the multipliers for the nonlinear constraints.
*        Check the condition of the initial factor R.
*        ---------------------------------------------------------------
         IF (NCNLN .GT. 0)
     $      CALL DCOPY ( NCNLN, CLAMDA(NPLIN+1), 1, W(LCMUL), 1 )
 
         SCLE  = ZERO
         SUMSQ = ONE
         DO 220 J = 1, N
            CALL DSSQ  ( J, R(1,J), 1, SCLE, SUMSQ )
  220    CONTINUE
         RFROBN = DNORM( SCLE, SUMSQ )
 
         CALL DCOND ( N, R, NROWR+1, DRMAX, DRMIN )
         COND   = DDIV  ( DRMAX, DRMIN, OVERFL )
 
         IF (      COND   .GT. RCNDBD
     $       .OR.  RFROBN .GT. ROOTN*GROWTH*DRMAX) THEN
*           ------------------------------------------------------------
*           Refactorize the Hessian and bound the condition estimator.
*           ------------------------------------------------------------
            CALL NPRSET( UNITQ,
     $                   N, NFREE, NZ, NQ, NROWR,
     $                   IW(LIPERM), IW(LKX),
     $                   W(LGQ), R, W(LZY), W(LWRK1), W(LRES0) )
         END IF
      END IF
 
*     ==================================================================
*     Solve the problem.
*     ==================================================================
      IF (NCNLN .EQ. 0) THEN
*        ---------------------------------------------------------------
*        The problem has only linear constraints and bounds.
*        ---------------------------------------------------------------
         CALL NPCORE( NAMED, NAMES, UNITQ, INFORM, ITER,
     $                N, NCLIN, NCNLN, NCTOTL, NACTIV, NFREE, NZ,
     $                NROWA, NROWJ, NROWUJ, NROWQP, NROWR,
     $                NFUN, NGRAD, ISTATE, IW(LKACTV), IW(LKX),
     $                OBJF, FDNORM, XNORM, OBJFUN, CONFUN,
     $                A, W(LAX), BL, BU, C, W(LCJAC), UJAC, CLAMDA,
     $                W(LFEATL), W(LGRAD), UGRAD, R, X, IW, W, LENW )
      ELSE
*        ---------------------------------------------------------------
*        The problem has some nonlinear constraints.
*        ---------------------------------------------------------------
         IF (NCLIN .GT. 0) THEN
            LA1J = LAQP
            DO 520 J = 1, N
               CALL DCOPY ( NCLIN, A(1,J), 1, W(LA1J), 1 )
               LA1J = LA1J + NROWQP
  520       CONTINUE
         END IF
 
*        Try and add some nonlinear constraint indices to KACTIV.
*
         CALL NPCRSH( COLD, N, NCLIN, NCNLN,
     $                NCTOTL, NACTIV, NFREE, NZ,
     $                ISTATE, IW(LKACTV), BIGBND, TOLACT,
     $                BL, BU, C )
 
         CALL NPCORE( NAMED, NAMES, UNITQ, INFORM, ITER,
     $                N, NCLIN, NCNLN, NCTOTL, NACTIV, NFREE, NZ,
     $                NROWA, NROWJ, NROWUJ, NROWQP, NROWR,
     $                NFUN, NGRAD, ISTATE, IW(LKACTV),IW(LKX),
     $                OBJF, FDNORM, XNORM, OBJFUN, CONFUN,
     $                W(LAQP), W(LAX), BL, BU, C, W(LCJAC),UJAC,CLAMDA,
     $                W(LFEATL), W(LGRAD), UGRAD, R, X, IW, W, LENW )
 
      END IF
 
*     ------------------------------------------------------------------
*     If required, form the triangular factor of the Hessian.
*     ------------------------------------------------------------------
*     First,  form the square matrix  R  such that  H = R'R.
*     Compute the  QR  factorization of  R.
 
      IF (LFORMH .GT. 0) THEN
         LV     = LWRK2
         DO 400 J = 1, N
            IF (J .GT. 1)
     $         CALL DLOAD ( J-1, ZERO, W(LV), 1 )
 
            LVJ = LV + J - 1
            CALL DCOPY ( N-J+1, R(J,J), NROWR, W(LVJ), 1     )
            CALL CMQMUL( 3, N, NZ, NFREE, NQ, UNITQ,
     $                   IW(LKX), W(LV), W(LZY), W(LWRK1) )
            CALL DCOPY ( N    , W(LV) , 1    , R(J,1), NROWR )
  400    CONTINUE
 
         CALL DGEQR ( N, N, R, NROWR, W(LWRK1), INFO )
      END IF
 
*     Print messages if required.
 
  800 IF (MSGNP .GT.   0) THEN
         IF (INFORM .LT.   0) WRITE (NOUT, 3000)
         IF (INFORM .EQ.   0) WRITE (NOUT, 4000)
         IF (INFORM .EQ.   1) WRITE (NOUT, 4100)
         IF (INFORM .EQ.   2) WRITE (NOUT, 4200)
         IF (INFORM .EQ.   3) WRITE (NOUT, 4300)
         IF (INFORM .EQ.   4) WRITE (NOUT, 4400)
         IF (INFORM .EQ.   5) WRITE (NOUT, 4500)
         IF (INFORM .EQ.   6) WRITE (NOUT, 4600)
         IF (INFORM .EQ.   7) WRITE (NOUT, 4700)
         IF (INFORM .EQ.   9) WRITE (NOUT, 4900) NERROR
 
         IF (INFORM .GE. 0  .AND.  INFORM .NE. 9) THEN
            IF (NLPERR .EQ. 0) THEN
               WRITE (NOUT, 5000) OBJF
            ELSE
               IF (NLPERR .EQ. 3) THEN
                  WRITE (NOUT, 5010) SUMINF
               ELSE
                  WRITE (NOUT, 5020) SUMINF
               END IF
            END IF
         END IF
      END IF
 
*     Recover the optional parameters set by the user.
 
      CALL ICOPY ( MXPARM, IPSVLS, 1, IPRMLS, 1 )
      CALL DCOPY ( MXPARM, RPSVLS, 1, RPRMLS, 1 )
      CALL ICOPY ( MXPARM, IPSVNP, 1, IPRMNP, 1 )
      CALL DCOPY ( MXPARM, RPSVNP, 1, RPRMNP, 1 )
 
      RETURN
 
 3000 FORMAT(/ ' Exit NPSOL - User requested termination.'          )
 4000 FORMAT(/ ' Exit NPSOL - Optimal solution found.'              )
 4100 FORMAT(/ ' Exit NPSOL - Optimal solution found, ',
     $         ' but the requested accuracy could not be achieved.' )
 4200 FORMAT(/ ' Exit NPSOL - No feasible point for the linear',
     $         ' constraints.')
 4300 FORMAT(/ ' Exit NPSOL - No feasible point for the nonlinear',
     $         ' constraints.')
 4400 FORMAT(/ ' Exit NPSOL - Too many major iterations.             ')
 4500 FORMAT(/ ' Exit NPSOL - Problem is unbounded (or badly scaled).')
 4600 FORMAT(/ ' Exit NPSOL - Current point cannot be improved upon. ')
 4700 FORMAT(/ ' Exit NPSOL - Large errors found in the derivatives. ')
 
 4900 FORMAT(/ ' Exit NPSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.')
 5000 FORMAT(/ ' Final nonlinear objective value =', G16.7 )
 5010 FORMAT(/ ' Minimum sum of infeasibilities =',  G16.7 )
 5020 FORMAT(/ ' Final sum of infeasibilities =',    G16.7 )
 
*     End of  NPSOL .
 
      END
