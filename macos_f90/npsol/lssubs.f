*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  LSSUBS FORTRAN
*
*     LSADD    LSADDS   LSBNDS   LSCHOL   LSCORE   LSCRSH   LSDEL
*     LSDFLT   LSFEAS   LSFILE   LSGETP   LSGSET   LSKEY    LSLOC
*     LSMOVE   LSMULS   LSOPTN   LSPRT    LSSETX   LSSOL
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSADD ( UNITQ,
     $                   INFORM, IFIX, IADD, JADD,
     $                   NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, NROWA, NQ, NROWR, NROWT,
     $                   KX, CONDMX,
     $                   A, R, T, RES, GQ, ZY, WRK1, WRK2 )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            KX(N)
      DOUBLE PRECISION   A(NROWA,*), R(NROWR,*), T(NROWT,*),
     $                   RES(N,*), GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N), WRK2(N)
************************************************************************
*  LSADD   updates the factorization,  A(free) * (Z Y) = (0 T),  when a
*  constraint is added to the working set.  If  NRANK .gt. 0, the
*  factorization  ( R ) = PWQ  is also updated,  where  W  is the
*                 ( 0 )
*  least squares matrix,  R  is upper-triangular,  and  P  is an
*  orthogonal matrix.  The matrices  W  and  P  are not stored.
*
*  There are three separate cases to consider (although each case
*  shares code with another)...
*
*  (1) A free variable becomes fixed on one of its bounds when there
*      are already some general constraints in the working set.
*
*  (2) A free variable becomes fixed on one of its bounds when there
*      are only bound constraints in the working set.
*
*  (3) A general constraint (corresponding to row  IADD  of  A) is
*      added to the working set.
*
*  In cases (1) and (2), we assume that  KX(IFIX) = JADD.
*  In all cases,  JADD  is the index of the constraint being added.
*
*  If there are no general constraints in the working set,  the
*  matrix  Q = (Z Y)  is the identity and will not be touched.
*
*  If  NRES .GT. 0,  the row transformations are applied to the rows of
*  the  (N by NRES)  matrix  RES.
*  If  NGQ .GT. 0,  the column transformations are applied to the
*  columns of the  (NGQ by N)  matrix  GQ'.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October--1984.
*  This version of LSADD dated 29-December-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      LOGICAL            BOUND , OVERFL
      EXTERNAL           DDOT  , DDIV  , DNRM2
      INTRINSIC          MAX   , MIN
      PARAMETER         (ZERO = 0.0D+0, ONE = 1.0D+0)
 
*     If the condition estimator of the updated factors is greater than
*     CONDBD,  a warning message is printed.
 
      CONDBD = ONE / EPSPT9
 
      OVERFL = .FALSE.
      BOUND  = JADD .LE. N
      IF (BOUND) THEN
*        ===============================================================
*        A simple bound has entered the working set.  IADD  is not used.
*        ===============================================================
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1010) NACTIV, NZ, NFREE, IFIX, JADD, UNITQ
         NANEW = NACTIV
 
         IF (UNITQ) THEN
 
*           Q  is not stored, but KX defines an ordering of the columns
*           of the identity matrix that implicitly define  Q.
*           Reorder KX so that variable IFIX is moved to position
*           NFREE+1 and variables IFIX+1,...,NFREE+1 are moved one
*           position to the left.
 
            CALL DLOAD ( NFREE, (ZERO), WRK1, 1 )
            WRK1(IFIX) = ONE
 
            DO 100 I = IFIX, NFREE-1
               KX(I) = KX(I+1)
  100       CONTINUE
         ELSE
*           ------------------------------------------------------------
*           Q  is stored explicitly.
*           ------------------------------------------------------------
*           Set  WRK1 = the  (IFIX)-th  row of  Q.
*           Move the  (NFREE)-th  row of  Q  to position  IFIX.
 
            CALL DCOPY ( NFREE, ZY(IFIX,1), NQ, WRK1, 1 )
            IF (IFIX .LT. NFREE) THEN
               CALL DCOPY ( NFREE, ZY(NFREE,1), NQ, ZY(IFIX,1), NQ )
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
*        ===============================================================
*        A general constraint has entered the working set.
*        IFIX  is not used.
*        ===============================================================
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1020) NACTIV, NZ, NFREE, IADD, JADD, UNITQ
 
         NANEW  = NACTIV + 1
 
*        Transform the incoming row of  A  by  Q'.
 
         CALL DCOPY ( N, A(IADD,1), NROWA, WRK1, 1 )
         CALL CMQMUL( 8, N, NZ, NFREE, NQ, UNITQ, KX, WRK1, ZY, WRK2)
 
*        Check that the incoming row is not dependent upon those
*        already in the working set.
 
         DTNEW  = DNRM2 ( NZ, WRK1, 1 )
         IF (NACTIV .EQ. 0) THEN
 
*           This is the only general constraint in the working set.
 
            COND   = DDIV  ( ASIZE, DTNEW, OVERFL )
            TDTMAX = DTNEW
            TDTMIN = DTNEW
         ELSE
 
*           There are already some general constraints in the working
*           set. Update the estimate of the condition number.
 
            TDTMAX = MAX( DTNEW, DTMAX )
            TDTMIN = MIN( DTNEW, DTMIN )
            COND   = DDIV  ( TDTMAX, TDTMIN, OVERFL )
         END IF
 
         IF (COND .GT. CONDMX  .OR.  OVERFL) GO TO 900
 
         IF (UNITQ) THEN
 
*           This is the first general constraint to be added.
*           Set  Q = I.
 
            DO 200 J = 1, NFREE
               CALL DLOAD ( NFREE, (ZERO), ZY(1,J), 1 )
               ZY(J,J) = ONE
  200       CONTINUE
            UNITQ  = .FALSE.
         END IF
      END IF
 
      NZERO  = NZ - 1
      IF (BOUND) NZERO = NFREE - 1
 
*     ------------------------------------------------------------------
*     Use a sequence of 2*2 column transformations to reduce the
*     first NZERO elements of WRK1 to zero.  This affects ZY, except
*     when UNITQ is true.  The transformations may also be applied
*     to R, T and GQ'.
*     ------------------------------------------------------------------
      LROWR  = N
      NELM   = 1
      IROWT  = NACTIV
 
      DO 300 K = 1, NZERO
 
*        Compute the transformation that reduces WRK1(K) to zero,
*        then apply it to the relevant columns of  Z  and  GQ'.
 
         CALL DROT3G( WRK1(K+1), WRK1(K), CS, SN )
         IF (.NOT. UNITQ)
     $      CALL DROT3 ( NFREE, ZY(1,K+1), 1, ZY(1,K), 1, CS, SN )
         IF (NGQ .GT. 0)
     $      CALL DROT3 ( NGQ  , GQ(K+1,1), N, GQ(K,1), N, CS, SN )
 
         IF (K .GE. NZ  .AND.  NACTIV .GT. 0) THEN
 
*           Apply the rotation to  T.
 
            T(IROWT,K) = ZERO
            CALL DROT3 ( NELM, T(IROWT,K+1), 1, T(IROWT,K), 1, CS, SN )
            NELM  = NELM  + 1
            IROWT = IROWT - 1
         END IF
 
         IF (NRANK .GT. 0) THEN
 
*           Apply the same transformation to the columns of R.
*           This generates a subdiagonal element in R that must be
*           eliminated by a row rotation.
 
            IF (K .LT. NRANK) R(K+1,K) = ZERO
            LCOL   = MIN( K+1, NRANK )
 
            CALL DROT3 ( LCOL, R(1,K+1), 1, R(1,K), 1, CS, SN )
            IF (K .LT. NRANK) THEN
               CALL DROT3G( R(K,K), R(K+1,K), CS, SN )
               LROWR  = LROWR - 1
               CALL DROT3 ( LROWR,   R(K,K+1)  , NROWR,
     $                               R(K+1,K+1), NROWR, CS, SN )
 
               IF (NRES .GT. 0)
     $            CALL DROT3 ( NRES, RES(K,1)  , N    ,
     $                               RES(K+1,1), N    , CS, SN )
            END IF
         END IF
  300 CONTINUE
 
      IF (BOUND) THEN
 
*        The last row and column of ZY has been transformed to plus
*        or minus the unit vector E(NFREE).  We can reconstitute the
*        columns of GQ and R corresponding to the new fixed variable.
 
         IF (WRK1(NFREE) .LT. ZERO) THEN
            NFMIN = MIN( NRANK, NFREE )
            IF (NFMIN .GT. 0) CALL DSCAL ( NFMIN, -ONE, R(1,NFREE) , 1 )
            IF (NGQ   .GT. 0) CALL DSCAL ( NGQ  , -ONE, GQ(NFREE,1), N )
         END IF
 
*        ---------------------------------------------------------------
*        The diagonals of T have been altered.  Recompute the
*        largest and smallest values.
*        ---------------------------------------------------------------
         IF (NACTIV .GT. 0) THEN
            CALL DCOND( NACTIV, T(NACTIV,NZ), NROWT-1, TDTMAX, TDTMIN )
            COND   = DDIV  ( TDTMAX, TDTMIN, OVERFL )
         END IF
      ELSE
*        ---------------------------------------------------------------
*        General constraint.  Install the new row of T.
*        ---------------------------------------------------------------
         CALL DCOPY ( NANEW, WRK1(NZ), 1, T(NANEW,NZ), NROWT )
      END IF
 
*     ==================================================================
*     Prepare to exit.  Check the magnitude of the condition estimator.
*     ==================================================================
  900 IF (NANEW .GT. 0) THEN
         IF (COND .LT. CONDMX  .AND.  .NOT. OVERFL) THEN
 
*           The factorization has been successfully updated.
 
            INFORM = 0
            DTMAX  = TDTMAX
            DTMIN  = TDTMIN
            IF (COND .GE. CONDBD) WRITE (NOUT, 2000) JADD
         ELSE
 
*           The proposed working set appears to be linearly dependent.
 
            INFORM = 1
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0) THEN
               WRITE( NOUT, 3000 )
               IF (BOUND) THEN
                  WRITE (NOUT, 3010) ASIZE, DTMAX, DTMIN
               ELSE
                  IF (NACTIV .GT. 0) THEN
                     WRITE (NOUT, 3020) ASIZE, DTMAX, DTMIN, DTNEW
                  ELSE
                     WRITE (NOUT, 3030) ASIZE, DTNEW
                  END IF
               END IF
            END IF
         END IF
      END IF
 
      RETURN
 
 1010 FORMAT(/ ' //LSADD //  Simple bound added.'
     $       / ' //LSADD //  NACTIV    NZ NFREE  IFIX  JADD UNITQ'
     $       / ' //LSADD //  ', 5I6, L6 )
 1020 FORMAT(/ ' //LSADD //  General constraint added.           '
     $       / ' //LSADD //  NACTIV    NZ NFREE  IADD  JADD UNITQ'
     $       / ' //LSADD //  ', 5I6, L6 )
 2000 FORMAT(/ ' XXX  Serious ill-conditioning in the working set',
     $         ' after adding constraint ',  I5
     $       / ' XXX  Overflow may occur in subsequent iterations.'//)
 3000 FORMAT(/ ' //LSADD //  Dependent constraint rejected.' )
 3010 FORMAT(/ ' //LSADD //     ASIZE     DTMAX     DTMIN        '
     $       / ' //LSADD //', 1P3E10.2 )
 3020 FORMAT(/ ' //LSADD //     ASIZE     DTMAX     DTMIN     DTNEW'
     $       / ' //LSADD //', 1P4E10.2 )
 3030 FORMAT(/ ' //LSADD //     ASIZE     DTNEW'
     $       / ' //LSADD //', 1P2E10.2 )
 
*     End of  LSADD .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSADDS( UNITQ, VERTEX,
     $                   INFORM, K1, K2, NACTIV, NARTIF, NZ, NFREE,
     $                   NRANK, NREJTD, NRES, NGQ,
     $                   N, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, KACTIV, KX,
     $                   CONDMX,
     $                   A, R, T, RES, GQ,
     $                   ZY, WRK1, WRK2 )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ, VERTEX
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      DOUBLE PRECISION   CONDMX
      DOUBLE PRECISION   A(NROWA,*), R(NROWR,*),
     $                   T(NROWT,*), RES(N,*), GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N), WRK2(N)
 
************************************************************************
*     LSADDS  includes general constraints K1 thru K2 as new rows of
*     the TQ factorization stored in T, ZY.  If NRANK is nonzero, the
*     changes in Q are reflected in NRANK by N triangular factor R such
*     that
*                         W  =  P ( R ) Q,
*                                 ( 0 )
*     where  P  is orthogonal.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  October-31-1984.
*     This version of LSADDS dated 30-December-1985.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN
 
      EXTERNAL           DNRM2
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      RTMAX  = WMACH(8)
 
*     Estimate the condition number of the constraints that are not
*     to be refactorized.
 
      IF (NACTIV .EQ. 0) THEN
         DTMAX = ZERO
         DTMIN = ONE
      ELSE
         CALL DCOND ( NACTIV, T(NACTIV,NZ+1), NROWT-1, DTMAX, DTMIN )
      END IF
 
      DO 200 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV .LT. NFREE) THEN
 
            CALL LSADD ( UNITQ,
     $                   INFORM, IFIX, IADD, JADD,
     $                   NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, NROWA, NQ, NROWR, NROWT,
     $                   KX, CONDMX,
     $                   A, R, T, RES, GQ, ZY,
     $                   WRK1, WRK2 )
 
            IF (INFORM .EQ. 0) THEN
               NACTIV = NACTIV + 1
               NZ     = NZ     - 1
            ELSE
               ISTATE(JADD) =   0
               KACTIV(K)    = - KACTIV(K)
            END IF
         END IF
  200 CONTINUE
 
      IF (NACTIV .LT. K2) THEN
 
*        Some of the constraints were classed as dependent and not
*        included in the factorization.  Re-order the part of  KACTIV
*        that holds the indices of the general constraints in the
*        working set.  Move accepted indices to the front and shift
*        rejected indices (with negative values) to the end.
 
         L      = K1 - 1
         DO 300 K = K1, K2
            I         = KACTIV(K)
            IF (I .GE. 0) THEN
               L      = L + 1
               IF (L .NE. K) THEN
                  ISWAP     = KACTIV(L)
                  KACTIV(L) = I
                  KACTIV(K) = ISWAP
               END IF
            END IF
  300    CONTINUE
 
*        If a vertex is required, add some temporary bounds.
*        We must accept the resulting condition number of the working
*        set.
 
         IF (VERTEX) THEN
            CNDMAX = RTMAX
            NZADD  = NZ
            DO 320 IARTIF = 1, NZADD
               ROWMAX = ZERO
               DO 310 I = 1, NFREE
                  RNORM = DNRM2 ( NZ, ZY(I,1), NQ )
                  IF (ROWMAX .LT. RNORM) THEN
                     ROWMAX = RNORM
                     IFIX   = I
                  END IF
  310          CONTINUE
               JADD = KX(IFIX)
 
               CALL LSADD ( UNITQ,
     $                      INFORM, IFIX, IADD, JADD,
     $                      NACTIV, NZ, NFREE, NRANK, NRES, NGQ,
     $                      N, NROWA, NQ, NROWR, NROWT,
     $                      KX, CNDMAX,
     $                      A, R, T, RES, GQ, ZY,
     $                      WRK1, WRK2 )
 
               NFREE  = NFREE  - 1
               NZ     = NZ     - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
  320       CONTINUE
         END IF
      END IF
 
      NREJTD = K2 - NACTIV
 
      RETURN
 
*     End of  LSADDS.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSBNDS( UNITQ,
     $                   INFORM, NZ, NFREE, NRANK, NRES, NGQ,
     $                   N, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, KX,
     $                   CONDMX,
     $                   A, R, T, RES, GQ,
     $                   ZY, WRK1, WRK2 )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            ISTATE(*), KX(N)
      DOUBLE PRECISION   CONDMX
      DOUBLE PRECISION   A(NROWA,*), R(NROWR,*),
     $                   T(NROWT,*), RES(N,*), GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WRK1(N), WRK2(N)
 
************************************************************************
*     LSBNDS updates the factor R as KX is reordered to reflect the
*     status of the bound constraints given by ISTATE.  KX is reordered
*     so that the fixed variables come last.  One of two alternative
*     are used to reorder KX. One method needs fewer accesses to KX, the
*     other gives a matrix Rz with more rows and columns.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  30-December-1985.
*     This version dated 30-December-1985.
************************************************************************
 
      NFIXED = N - NFREE
 
      IF (NRANK .LT. N  .AND.  NRANK .GT. 0) THEN
*        ---------------------------------------------------------------
*        R is specified but singular.  Try and keep the dimension of Rz
*        as large as possible.
*        ---------------------------------------------------------------
         NACTV = 0
         NFREE = N
         NZ    = N
 
         J     = N
*+       WHILE (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) DO
  100    IF    (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) THEN
            IF (ISTATE(J) .GT. 0) THEN
               JADD = J
               DO 110 IFIX = NFREE, 1, -1
                  IF (KX(IFIX) .EQ. JADD) GO TO 120
  110          CONTINUE
 
*              Add bound JADD.
 
  120          CALL LSADD ( UNITQ,
     $                      INFORM, IFIX, IADD, JADD,
     $                      NACTV, NZ, NFREE, NRANK, NRES, NGQ,
     $                      N, NROWA, NQ, NROWR, NROWT,
     $                      KX, CONDMX,
     $                      A, R, T, RES, GQ, ZY,
     $                      WRK1, WRK2 )
 
               NFREE = NFREE - 1
               NZ    = NZ    - 1
            END IF
            J = J - 1
            GO TO 100
*+       END WHILE
         END IF
      ELSE
*        ---------------------------------------------------------------
*        R is of full rank,  or is not specified.
*        ---------------------------------------------------------------
         IF (NFIXED .GT. 0) THEN
 
*           Order KX so that the free variables come first.
 
            LSTART = NFREE + 1
            DO 250 K = 1, NFREE
               J = KX(K)
               IF (ISTATE(J) .GT. 0) THEN
                  DO 220 L = LSTART, N
                     J2 = KX(L)
                     IF (ISTATE(J2) .EQ. 0) GO TO 230
  220             CONTINUE
 
  230             KX(K)  = J2
                  KX(L)  = J
                  LSTART = L + 1
 
                  IF (NRANK .GT. 0)
     $               CALL CMRSWP( N, NRES, NRANK, NROWR, K, L,
     $                            R, RES, WRK1 )
               END IF
  250       CONTINUE
 
         END IF
         NZ = NFREE
      END IF
 
      RETURN
 
*     End of  LSBNDS.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSCHOL( NROWH, N, NRANK, TOLRNK, KX, H, INFORM )
 
      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      INTEGER            KX(*)
      DOUBLE PRECISION   H(NROWH,*)
 
************************************************************************
*     LSCHOL  forms the Cholesky factorization of the positive
*     semi-definite matrix H such that
*                   PHP'  =  R'R
*     where  P  is a permutation matrix and  R  is upper triangular.
*     The permutation P is chosen to maximize the diagonal of R at each
*     stage.  Only the diagonal and super-diagonal elements of H are
*     used.
*
*     Output:
*
*         INFORM = 0   the factorization was computed successfully,
*                      with the Cholesky factor written in the upper
*                      triangular part of H and P stored in KX.
*                  1   the matrix H was indefinite.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version of LSCHOL dated  2-February-1981.
*     Level 2 Blas added 29-June-1986.
*     This version of LSCHOL dated  30-June-1986.
************************************************************************
 
      COMMON    /SOL1CM/ NOUT
      INTRINSIC          ABS   , MAX   , SQRT
      EXTERNAL           IDAMAX
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      INFORM = 0
      NRANK  = 0
 
*     Main loop for computing rows of  R.
 
      DO 200 J = 1, N
 
*        Find maximum available diagonal.
 
         KMAX = J - 1 + IDAMAX( N-J+1, H(J,J), NROWH+1 )
         DMAX = H(KMAX,KMAX)
 
         IF (DMAX .LE. TOLRNK*ABS(H(1,1))) GO TO 300
 
*        Perform a symmetric interchange if necessary.
 
         IF (KMAX .NE. J) THEN
            K        = KX(KMAX)
            KX(KMAX) = KX(J)
            KX(J)    = K
 
            CALL DSWAP ( J       , H(1,J)   , 1, H(1,KMAX), 1     )
            CALL DSWAP ( KMAX-J+1, H(J,KMAX), 1, H(J,J)   , NROWH )
            CALL DSWAP ( N-KMAX+1, H(KMAX,KMAX), NROWH,
     $                             H(J,KMAX)   , NROWH )
 
            H(KMAX,KMAX) = H(J,J)
         END IF
 
*        Set the diagonal of  R.
 
         D      = SQRT( DMAX )
         H(J,J) = D
         NRANK  = NRANK + 1
 
         IF (J .LT. N) THEN
 
*           Set the super-diagonal elements of this row of R and update
*           the elements of the block that is yet to be factorized.
 
            CALL DSCAL ( N-J,   (ONE/D), H(J  ,J+1), NROWH )
            CALL DSYR  ( 'U', N-J, -ONE, H(J  ,J+1), NROWH,
     $                                   H(J+1,J+1), NROWH )
         END IF
 
  200 CONTINUE
*     ------------------------------------------------------------------
*     Check for the semi-definite case.
*     ------------------------------------------------------------------
  300 IF (NRANK .LT. N) THEN
 
*        Find the largest element in the unfactorized block.
 
         SUPMAX = ZERO
         DO 310 I = J, N-1
            K      = I + IDAMAX( N-I, H(I,I+1), NROWH )
            SUPMAX = MAX( SUPMAX, ABS(H(I,K)) )
  310    CONTINUE
 
         IF (SUPMAX .GT. TOLRNK*ABS(H(1,1))) THEN
            WRITE (NOUT, 1000) DMAX, SUPMAX
            INFORM = 1
         END IF
      END IF
 
      RETURN
 
 1000 FORMAT(' XXX  Hessian appears to be indefinite.'
     $      /' XXX  Maximum diagonal and off-diagonal ignored',
     $             ' in the Cholesky factorization:', 1P2E22.14 )
 
*     End of LSCHOL.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSCORE( PRBTYP, NAMED, NAMES, LINOBJ, UNITQ,
     $                   INFORM, ITER, JINF, NCLIN, NCTOTL,
     $                   NACTIV, NFREE, NRANK, NZ, NZ1,
     $                   N, NROWA, NROWR,
     $                   ISTATE, KACTIV, KX,
     $                   CTX, SSQ, SSQ1, SUMINF, NUMINF, XNORM,
     $                   BL, BU, A, CLAMDA, AX,
     $                   FEATOL, R, X, IW, W )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      CHARACTER*8        NAMES(*)
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      INTEGER            IW(*)
      DOUBLE PRECISION   BL(NCTOTL), BU(NCTOTL), A(NROWA,*),
     $                   CLAMDA(NCTOTL), AX(*),
     $                   FEATOL(NCTOTL), R(NROWR,*), X(N)
      DOUBLE PRECISION   W(*)
      LOGICAL            NAMED, LINOBJ, UNITQ
 
************************************************************************
*     LSCORE  is a subroutine for linearly constrained linear-least
*     squares.  On entry, it is assumed that an initial working set of
*     linear constraints and bounds is available.
*     The arrays ISTATE, KACTIV and KX will have been set accordingly
*     and the arrays T and ZY will contain the TQ factorization of
*     the matrix whose rows are the gradients of the active linear
*     constraints with the columns corresponding to the active bounds
*     removed.  the TQ factorization of the resulting (NACTIV by NFREE)
*     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
*     is reverse-triangular.
*
*     Values of ISTATE(J) for the linear constraints.......
*
*     ISTATE(J)
*     ---------
*          0    constraint J is not in the working set.
*          1    constraint J is in the working set at its lower bound.
*          2    constraint J is in the working set at its upper bound.
*          3    constraint J is in the working set as an equality.
*
*     Constraint J may be violated by as much as FEATOL(J).
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  LSCORE  dated  1-August-1986.
*
*     Copyright  1984  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under grants MCS-7926009 and ECS-8012974; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; and the Army Research Office Contract DAA29-79-C-0110.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN
 
      INTEGER            LOCLS
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      LOGICAL            CMDBG, LSDBG
      PARAMETER         (LDBG = 5)
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
      EQUIVALENCE   (MSGLS , MSGLVL), (IDBGLS, IDBG), (LDBGLS, MSGDBG)
 
      EXTERNAL           DDIV  , DDOT  , DNRM2
      INTRINSIC          ABS   , MAX   , SQRT
      LOGICAL            CONVRG, CYCLIN, ERROR , FIRSTV, HITCON,
     $                   HITLOW, NEEDFG, OVERFL, PRNT  , PRNT1 , ROWERR
      LOGICAL            SINGLR, STALL , STATPT, UNBNDD, UNCON , UNITGZ,
     $                   WEAK
      PARAMETER        ( ZERO   =0.0D+0, HALF   =0.5D+0, ONE   =1.0D+0 )
      PARAMETER        ( MREFN  =1     , MSTALL =50                    )
 
*     Specify the machine-dependent parameters.
 
      EPSMCH = WMACH(3)
      FLMAX  = WMACH(7)
      RTMAX  = WMACH(8)
 
      LANORM = LOCLS( 2)
      LAP    = LOCLS( 3)
      LPX    = LOCLS( 4)
      LRES   = LOCLS( 5)
      LRES0  = LOCLS( 6)
      LHZ    = LOCLS( 7)
      LGQ    = LOCLS( 8)
      LCQ    = LOCLS( 9)
      LRLAM  = LOCLS(10)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK   = LOCLS(14)
 
*     Set up the adresses of the contiguous arrays  ( RES0, RES )
*     and  ( GQ, CQ ).
 
      NRES   = 0
      IF (NRANK .GT. 0) NRES = 2
      NGQ    = 1
      IF (LINOBJ) NGQ = 2
 
*     Initialize.
 
      IREFN  =   0
      ITER   =   0
      ITMAX  =   ITMAX1
      JADD   =   0
      JDEL   =   0
      NCNLN  =   0
      NPHASE =   1
      NSTALL =   0
      NUMINF = - 1
      NZ1    =   0
 
      ALFA   = ZERO
      CONDMX = FLMAX
      DRZMAX = ONE
      DRZMIN = ONE
      SSQ    = ZERO
 
      CYCLIN = .FALSE.
      ERROR  = .FALSE.
      FIRSTV = .FALSE.
      PRNT   = .TRUE.
      PRNT1  = .TRUE.
      NEEDFG = .TRUE.
      STALL  = .TRUE.
      UNCON  = .FALSE.
      UNBNDD = .FALSE.
 
*     If debug output is required,  print nothing until iteration IDBG.
 
      MSGSVD = MSGLVL
      IF (IDBG .GT. 0  .AND.  IDBG .LE. ITMAX) THEN
         MSGLVL = 0
      END IF
 
*======================== start of the main loop =======================
*
*      cyclin = false
*      unbndd = false
*      error  = false
*      k      = 0
*
*      repeat
*            repeat
*                  compute Z'g,  print details of this iteration
*                  stat pt = (Z'g .eq. 0)
*                  if (not stat pt) then
*                     error =  k .ge. itmax
*                     if (not error) then
*                        compute p, alfa
*                        error = unbndd  or  cyclin
*                        if (not error) then
*                           k = k + 1
*                           x = x + alfa p
*                           if (feasible) update Z'g
*                           if necessary, add a constraint
*                        end if
*                     end if
*                  end if
*            until  stat pt  or  error
*
*            compute lam1, lam2, smllst
*            optmul =  smllst .gt. 0
*            if ( not (optmul .or. error) ) then
*                  delete an artificial or regular constraint
*            end if
*      until optmul  or  error
*
*=======================================================================
 
*     REPEAT
*        REPEAT
  100       IF (NEEDFG) THEN
               IF (NRANK .GT. 0) THEN
                  RESNRM = DNRM2 ( NRANK, W(LRES), 1 )
                  SSQ    = HALF*(SSQ1**2 + RESNRM**2 )
               END IF
 
               IF (NUMINF .NE. 0) THEN
 
*                 Compute the transformed gradient of either the sum of
*                 of infeasibilities or the objective.  Initialize
*                 SINGLR and UNITGZ.
 
                  CALL LSGSET( PRBTYP, LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                         N, NCLIN, NFREE,
     $                         NROWA, NQ, NROWR, NRANK, NZ, NZ1,
     $                         ISTATE, KX,
     $                         BIGBND, TOLRNK, NUMINF, SUMINF,
     $                         BL, BU, A, W(LRES), FEATOL,
     $                         W(LGQ), W(LCQ), R, X, W(LWTINF),
     $                         W(LZY), W(LWRK) )
 
                  IF (PRBTYP .NE. 'FP'  .AND.  NUMINF .EQ. 0
     $                                  .AND.  NPHASE .EQ. 1) THEN
                     ITMAX  = ITER + ITMAX2
                     NPHASE = 2
                  END IF
               END IF
            END IF
 
            GZNORM = ZERO
            IF (NZ  .GT. 0 ) GZNORM = DNRM2 ( NZ, W(LGQ), 1 )
 
            IF (NZ1 .EQ. NZ) THEN
               GZ1NRM = GZNORM
            ELSE
               GZ1NRM = ZERO
               IF (NZ1 .GT. 0) GZ1NRM = DNRM2 ( NZ1, W(LGQ), 1 )
            END IF
 
            GFNORM = GZNORM
            IF (NFREE .GT. 0  .AND.  NACTIV .GT. 0)
     $         GFNORM = DNRM2 ( NFREE, W(LGQ), 1 )
 
*           ------------------------------------------------------------
*           Print the details of this iteration.
*           ------------------------------------------------------------
*           Define small quantities that reflect the size of X, R and
*           the constraints in the working set.  If feasible,  estimate
*           the rank and condition number of Rz1.
*           Note that NZ1 .LE. NRANK + 1.
 
            IF (NZ1 .EQ. 0) THEN
               SINGLR = .FALSE.
            ELSE
               IF (NUMINF .GT. 0  .OR.  NZ1 .GT. NRANK) THEN
                  ABSRZZ = ZERO
               ELSE
                  CALL DCOND ( NZ1, R, NROWR+1, DRZMAX, DRZMIN )
                  ABSRZZ = ABS( R(NZ1,NZ1) )
               END IF
               SINGLR = ABSRZZ .LE. DRZMAX*TOLRNK
 
               IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $            WRITE (NOUT, 9100) SINGLR, ABSRZZ, DRZMAX, DRZMIN
 
            END IF
 
            CONDRZ = DDIV  ( DRZMAX, DRZMIN, OVERFL )
            CONDT  = ONE
            IF (NACTIV .GT. 0)
     $         CONDT  = DDIV  ( DTMAX , DTMIN , OVERFL )
 
            IF (PRNT) THEN
               CALL LSPRT ( PRBTYP, PRNT1, ISDEL, ITER, JADD, JDEL,
     $                      MSGLVL, NACTIV, NFREE, N, NCLIN,
     $                      NRANK, NROWR, NROWT, NZ, NZ1,
     $                      ISTATE,
     $                      ALFA, CONDRZ, CONDT, GFNORM, GZNORM, GZ1NRM,
     $                      NUMINF, SUMINF, CTX, SSQ,
     $                      AX, R, W(LT), X, W(LWRK) )
 
               JDEL  = 0
               JADD  = 0
               ALFA  = ZERO
            END IF
 
            IF (NUMINF .GT. 0) THEN
               DINKY  = ZERO
            ELSE
               OBJSIZ = ONE  + ABS( SSQ + CTX )
               WSSIZE = ZERO
               IF (NACTIV .GT. 0) WSSIZE = DTMAX
               DINKY  = EPSPT8 * MAX( WSSIZE, OBJSIZ, GFNORM )
               IF (UNCON) THEN
                  UNITGZ = GZ1NRM .LE. DINKY
               END IF
            END IF
 
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $         WRITE (NOUT, 9000) UNITGZ, IREFN, GZ1NRM, DINKY
 
*           If the projected gradient  Z'g  is small and Rz is of full
*           rank, X is a minimum on the working set.  An additional
*           refinement step is allowed to take care of an inaccurate
*           value of DINKY.
 
            STATPT = .NOT. SINGLR  .AND.  GZ1NRM .LE. DINKY
     $                             .OR.   IREFN  .GT. MREFN
 
            IF (.NOT. STATPT) THEN
*              ---------------------------------------------------------
*              Compute a search direction.
*              ---------------------------------------------------------
               PRNT  = .TRUE.
 
               ERROR = ITER .GE. ITMAX
               IF (.NOT. ERROR) THEN
 
                  IREFN = IREFN + 1
                  ITER  = ITER  + 1
 
                  IF (ITER .EQ. IDBG) THEN
                     LSDBG  = .TRUE.
                     CMDBG  =  LSDBG
                     MSGLVL =  MSGSVD
                  END IF
 
                  CALL LSGETP( LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                         N, NCLIN, NFREE,
     $                         NROWA, NQ, NROWR, NRANK, NUMINF, NZ1,
     $                         ISTATE, KX, CTP, PNORM,
     $                         A, W(LAP), W(LRES), W(LHZ), W(LPX),
     $                         W(LGQ), W(LCQ), R, W(LZY), W(LWRK) )
 
*                 ------------------------------------------------------
*                 Find the constraint we bump into along P.
*                 Update X and AX if the step ALFA is nonzero.
*                 ------------------------------------------------------
*                 ALFHIT is initialized to BIGALF.  If it remains
*                 that way after the call to CMALF, it will be
*                 regarded as infinite.
 
                  BIGALF = DDIV  ( BIGDX, PNORM, OVERFL )
 
                  CALL CMALF ( FIRSTV, HITLOW,
     $                         ISTATE, INFORM, JADD, N, NROWA,
     $                         NCLIN, NCTOTL, NUMINF,
     $                         ALFHIT, PALFA, ATPHIT,
     $                         BIGALF, BIGBND, PNORM,
     $                         W(LANORM), W(LAP), AX,
     $                         BL, BU, FEATOL, W(LPX), X )
 
*                 If  Rz1  is nonsingular,  ALFA = 1.0  will be the
*                 step to the least-squares minimizer on the
*                 current subspace. If the unit step does not violate
*                 the nearest constraint by more than FEATOL,  the
*                 constraint is not added to the working set.
 
                  HITCON = SINGLR  .OR.  PALFA  .LE. ONE
                  UNCON  = .NOT. HITCON
 
                  IF (HITCON) THEN
                     ALFA = ALFHIT
                  ELSE
                     JADD   = 0
                     ALFA   = ONE
                  END IF
 
*                 Check for an unbounded solution or negligible step.
 
                  UNBNDD =  ALFA .GE. BIGALF
                  STALL  = ABS( ALFA*PNORM ) .LE. EPSPT9*XNORM
                  IF (STALL) THEN
                     NSTALL = NSTALL + 1
                     CYCLIN = NSTALL .GT. MSTALL
                  ELSE
                     NSTALL = 0
                  END IF
 
                  ERROR = UNBNDD  .OR.  CYCLIN
                  IF (.NOT.  ERROR) THEN
*                    ---------------------------------------------------
*                    Set X = X + ALFA*P.  Update AX, GQ, RES and CTX.
*                    ---------------------------------------------------
                     IF (ALFA .NE. ZERO)
     $                  CALL LSMOVE( HITCON, HITLOW, LINOBJ, UNITGZ,
     $                               NCLIN, NRANK, NZ1,
     $                               N, NROWR, JADD, NUMINF,
     $                               ALFA, CTP, CTX, XNORM,
     $                               W(LAP), AX, BL, BU, W(LGQ),
     $                               W(LHZ), W(LPX), W(LRES),
     $                               R, X, W(LWRK) )
 
                     IF (HITCON) THEN
*                       ------------------------------------------------
*                       Add a constraint to the working set.
*                       Update the TQ factors of the working set.
*                       Use P as temporary work space.
*                       ------------------------------------------------
*                       Update  ISTATE.
 
                        IF (BL(JADD) .EQ. BU(JADD)) THEN
                           ISTATE(JADD) = 3
                        ELSE IF (HITLOW) THEN
                           ISTATE(JADD) = 1
                        ELSE
                           ISTATE(JADD) = 2
                        END IF
                        IADD = JADD - N
                        IF (JADD .LE. N) THEN
 
                           DO 510 IFIX = 1, NFREE
                              IF (KX(IFIX) .EQ. JADD) GO TO 520
  510                      CONTINUE
  520                   END IF
 
                        CALL LSADD ( UNITQ,
     $                               INFORM, IFIX, IADD, JADD,
     $                               NACTIV, NZ, NFREE, NRANK, NRES,NGQ,
     $                               N, NROWA, NQ, NROWR, NROWT,
     $                               KX, CONDMX,
     $                               A, R, W(LT), W(LRES), W(LGQ),
     $                               W(LZY), W(LWRK), W(LRLAM) )
 
                        NZ1    = NZ1 - 1
                        NZ     = NZ  - 1
 
                        IF (JADD .LE. N) THEN
 
*                          A simple bound has been added.
 
                           NFREE  = NFREE  - 1
                        ELSE
 
*                          A general constraint has been added.
 
                           NACTIV = NACTIV + 1
                           KACTIV(NACTIV) = IADD
                        END IF
 
                        IREFN  = 0
 
                     END IF
 
*                    ---------------------------------------------------
*                    Check the feasibility of constraints with non-
*                    negative ISTATE values.  If some violations have
*                    occurred.  Refine the current X and set INFORM so
*                    that feasibility is checked in LSGSET.
*                    ---------------------------------------------------
                     CALL LSFEAS( N, NCLIN, ISTATE,
     $                            BIGBND, CNORM, ERR1, JMAX1, NVIOL,
     $                            AX, BL, BU, FEATOL, X, W(LWRK) )
 
                     IF (ERR1 .GT. FEATOL(JMAX1)) THEN
                        CALL LSSETX( LINOBJ, ROWERR, UNITQ,
     $                               NCLIN, NACTIV, NFREE, NRANK, NZ,
     $                               N, NCTOTL, NQ, NROWA, NROWR, NROWT,
     $                               ISTATE, KACTIV, KX,
     $                               JMAX1, ERR2, CTX, XNORM,
     $                               A, AX, BL, BU, W(LCQ),
     $                               W(LRES), W(LRES0), FEATOL, R,
     $                               W(LT), X, W(LZY), W(LPX), W(LWRK) )
 
                        IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $                     WRITE (NOUT, 2100) ERR1, ERR2
                        IF (ROWERR)       WRITE (NOUT, 2200)
 
                        UNCON  =   .FALSE.
                        IREFN  =   0
                        NUMINF = - 1
                     END IF
                     NEEDFG = ALFA .NE. ZERO
                  END IF
               END IF
            END IF
 
*        UNTIL      STATPT  .OR.  ERROR
         IF (.NOT. (STATPT  .OR.  ERROR) ) GO TO 100
 
*        ===============================================================
*        Try and find the index JDEL of a constraint to drop from
*        the working set.
*        ===============================================================
         JDEL   = 0
 
         IF (NUMINF .EQ. 0  .AND.  PRBTYP .EQ. 'FP') THEN
            IF (N .GT. NZ)
     $         CALL DLOAD ( N-NZ, (ZERO), W(LRLAM), 1 )
            JTINY  = 0
            JSMLST = 0
            JBIGST = 0
         ELSE
 
            CALL LSMULS( PRBTYP,
     $                   MSGLVL, N, NACTIV, NFREE,
     $                   NROWA, NROWT, NUMINF, NZ, NZ1,
     $                   ISTATE, KACTIV, KX, DINKY,
     $                   JSMLST, KSMLST, JINF, JTINY,
     $                   JBIGST, KBIGST, TRULAM,
     $                   A, W(LANORM), W(LGQ), W(LRLAM),
     $                   W(LT), W(LWTINF) )
 
         END IF
 
         IF (.NOT. ERROR) THEN
            IF (     JSMLST .GT. 0) THEN
 
*              LSMULS found a regular constraint with multiplier less
*              than (-DINKY).
 
               JDEL   = JSMLST
               KDEL   = KSMLST
               ISDEL  = ISTATE(JDEL)
               ISTATE(JDEL) = 0
 
            ELSE IF (JSMLST .LT. 0) THEN
 
               JDEL   = JSMLST
 
            ELSE IF (NUMINF .GT. 0  .AND.  JBIGST .GT. 0) THEN
 
*              No feasible point exists for the constraints but the
*              sum of the constraint violations may be reduced by
*              moving off constraints with multipliers greater than 1.
 
               JDEL   = JBIGST
               KDEL   = KBIGST
               ISDEL  = ISTATE(JDEL)
               IF (TRULAM .LE. ZERO) IS = - 1
               IF (TRULAM .GT. ZERO) IS = - 2
               ISTATE(JDEL) = IS
               FIRSTV = .TRUE.
               NUMINF = NUMINF + 1
            END IF
 
            IF      (JDEL .NE. 0  .AND.  SINGLR) THEN
 
*              Cannot delete a constraint when Rz is singular.
*              Probably a weak minimum.
 
               JDEL = 0
            ELSE IF (JDEL .NE. 0               ) THEN
 
*              Constraint JDEL has been deleted.
*              Update the matrix factorizations.
 
               CALL LSDEL ( UNITQ,
     $                      N, NACTIV, NFREE, NRES, NGQ, NZ, NZ1,
     $                      NROWA, NQ, NROWR, NROWT, NRANK,
     $                      JDEL, KDEL, KACTIV, KX,
     $                      A, W(LRES), R, W(LT), W(LGQ),W(LZY),W(LWRK))
 
            END IF
         END IF
 
         IREFN  =  0
         CONVRG =  JDEL .EQ. 0
         PRNT   = .FALSE.
         UNCON  = .FALSE.
         NEEDFG = .FALSE.
 
*     until       convrg  .or.  error
      IF (.NOT.  (CONVRG  .OR.  ERROR)) GO TO 100
 
*  .........................End of main loop............................
 
      WEAK = JTINY .GT. 0  .OR.  SINGLR
 
      IF (ERROR) THEN
         IF (UNBNDD) THEN
            INFORM = 2
            IF (NUMINF .GT. 0) INFORM = 3
         ELSE IF (ITER .GE. ITMAX) THEN
            INFORM = 4
         ELSE IF (CYCLIN) THEN
            INFORM = 5
         END IF
      ELSE IF (CONVRG) THEN
         INFORM = 0
         IF (NUMINF .GT. 0) THEN
            INFORM = 3
         ELSE IF (PRBTYP .NE. 'FP'  .AND.  WEAK) THEN
            INFORM = 1
         END IF
      END IF
 
*     ------------------------------------------------------------------
*     Set   CLAMDA.  Print the full solution.
*     ------------------------------------------------------------------
      MSGLVL = MSGSVD
      IF (MSGLVL .GT. 0) WRITE (NOUT, 2000) PRBTYP, ITER, INFORM
 
      CALL CMPRT ( MSGLVL, NFREE, NROWA,
     $             N, NCLIN, NCNLN, NCTOTL, BIGBND,
     $             NAMED, NAMES, LENNAM,
     $             NACTIV, ISTATE, KACTIV, KX,
     $             A, BL, BU, X, CLAMDA, W(LRLAM), X )
 
      RETURN
 
 2000 FORMAT(/ ' Exit from ', A2, ' problem after ', I4, ' iterations.',
     $         '  INFORM =', I3 )
 2100 FORMAT(  ' XXX  Iterative refinement.  Maximum errors before and',
     $         ' after refinement are ',  1P2E14.2 )
 2200 FORMAT(  ' XXX  Warning.  Cannot satisfy the constraints to the',
     $         ' accuracy requested.')
 9000 FORMAT(/ ' //LSCORE//  UNITGZ IREFN     GZ1NRM      DINKY'
     $       / ' //LSCORE//  ', L6, I6, 1P2E11.2 )
 9100 FORMAT(/ ' //LSCORE//  SINGLR   ABS(RZZ1)      DRZMAX      DRZMIN'
     $       / ' //LSCORE//  ', L6,     1P3E12.4 )
 
*     End of  LSCORE.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSCRSH( COLD, VERTEX,
     $                   NCLIN, NCTOTL, NACTIV, NARTIF,
     $                   NFREE, N, NROWA,
     $                   ISTATE, KACTIV,
     $                   BIGBND, TOLACT,
     $                   A, AX, BL, BU, X, WX, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            COLD, VERTEX
      INTEGER            ISTATE(NCTOTL), KACTIV(N)
      DOUBLE PRECISION   A(NROWA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     $                   X(N), WX(N), WORK(N)
 
************************************************************************
*     LSCRSH  computes the quantities  ISTATE (optionally), KACTIV,
*     NACTIV, NZ and NFREE  associated with the working set at X.
*     The computation depends upon the value of the input parameter
*     COLD,  as follows...
*
*     COLD = TRUE.  An initial working set will be selected. First,
*                   nearly-satisfied or violated bounds are added.
*                   Next,  general linear constraints are added that
*                   have small residuals.
*
*     COLD = FALSE. The quantities KACTIV, NACTIV, NZ and NFREE are
*                   computed from ISTATE,  specified by the user.
*
*     Values of ISTATE(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSCRSH dated 27-December-1985.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      EXTERNAL           DDOT
      INTRINSIC          ABS, MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      FLMAX  = WMACH(7)
      CALL DCOPY ( N, X, 1, WX, 1 )
 
      IF (LSDBG) THEN
         IF (ILSDBG(1) .GT. 0)
     $      WRITE (NOUT, 1000) COLD, NCLIN, NCTOTL
         IF (ILSDBG(2) .GT. 0)
     $      WRITE (NOUT, 1100) (WX(J), J = 1, N)
      END IF
 
      NFIXED = 0
      NACTIV = 0
      NARTIF = 0
 
*     If a cold start is being made, initialize  ISTATE.
*     If  BL(j) = BU(j),  set  ISTATE(j)=3  for all variables and linear
*     constraints.
 
      IF (COLD) THEN
         DO 100 J = 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J) .EQ. BU(J)) ISTATE(J) = 3
  100    CONTINUE
      ELSE
         DO 110 J = 1, NCTOTL
            IF (ISTATE(J) .GT. 3  .OR.  ISTATE(J) .LT. 0) ISTATE(J) = 0
  110    CONTINUE
      END IF
 
*     Initialize NFIXED, NFREE and KACTIV.
*     Ensure that the number of bounds and general constraints in the
*     working set does not exceed N.
 
      DO 200 J = 1, NCTOTL
         IF (NFIXED + NACTIV .EQ. N) ISTATE(J) = 0
         IF (ISTATE(J) .GT. 0) THEN
            IF (J .LE. N) THEN
               NFIXED = NFIXED + 1
               IF (ISTATE(J) .EQ. 1) WX(J) = BL(J)
               IF (ISTATE(J) .GE. 2) WX(J) = BU(J)
            ELSE
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = J - N
            END IF
         END IF
  200 CONTINUE
 
*     ------------------------------------------------------------------
*     If a cold start is required,  attempt to add as many
*     constraints as possible to the working set.
*     ------------------------------------------------------------------
      IF (COLD) THEN
         BIGLOW = - BIGBND
         BIGUPP =   BIGBND
 
*        See if any bounds are violated or nearly satisfied.
*        If so,  add these bounds to the working set and set the
*        variables exactly on their bounds.
 
         J = N
*+       WHILE (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) DO
  300    IF    (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) THEN
            IF (ISTATE(J) .EQ. 0) THEN
               B1     = BL(J)
               B2     = BU(J)
               IS     = 0
               IF (B1 .GT. BIGLOW) THEN
                  IF (WX(J) - B1 .LE. (ONE + ABS( B1 ))*TOLACT) IS = 1
               END IF
               IF (B2 .LT. BIGUPP) THEN
                  IF (B2 - WX(J) .LE. (ONE + ABS( B2 ))*TOLACT) IS = 2
               END IF
               IF (IS .GT. 0) THEN
                  ISTATE(J) = IS
                  IF (IS .EQ. 1) WX(J) = B1
                  IF (IS .EQ. 2) WX(J) = B2
                  NFIXED = NFIXED + 1
               END IF
            END IF
            J = J - 1
            GO TO 300
*+       END WHILE
         END IF
 
*        ---------------------------------------------------------------
*        The following loop finds the linear constraint (if any) with
*        smallest residual less than or equal to TOLACT  and adds it
*        to the working set.  This is repeated until the working set
*        is complete or all the remaining residuals are too large.
*        ---------------------------------------------------------------
*        First, compute the residuals for all the constraints not in the
*        working set.
 
         IF (NCLIN .GT. 0  .AND.  NACTIV+NFIXED .LT. N) THEN
            DO 410 I = 1, NCLIN
               IF (ISTATE(N+I) .LE. 0)
     $         AX(I) = DDOT  (N, A(I,1), NROWA, WX, 1 )
  410       CONTINUE
 
            IS     = 1
            TOOBIG = TOLACT + TOLACT
 
*+          WHILE (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) DO
  500       IF    (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) THEN
               IS     = 0
               RESMIN = TOLACT
 
               DO 520 I = 1, NCLIN
                  J      = N + I
                  IF (ISTATE(J) .EQ. 0) THEN
                     B1     = BL(J)
                     B2     = BU(J)
                     RESL   = TOOBIG
                     RESU   = TOOBIG
                     IF (B1 .GT. BIGLOW)
     $                  RESL  = ABS( AX(I) - B1 ) / (ONE + ABS( B1 ))
                     IF (B2 .LT. BIGUPP)
     $                  RESU  = ABS( AX(I) - B2 ) / (ONE + ABS( B2 ))
                     RESIDL   = MIN( RESL, RESU )
                     IF(RESIDL .LT. RESMIN) THEN
                        RESMIN = RESIDL
                        IMIN   = I
                        IS     = 1
                        IF (RESL .GT. RESU) IS = 2
                     END IF
                  END IF
  520          CONTINUE
 
               IF (IS .GT. 0) THEN
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IMIN
                  J         = N + IMIN
                  ISTATE(J) = IS
               END IF
               GO TO 500
*+          END WHILE
            END IF
         END IF
 
*        ---------------------------------------------------------------
*        If required, add temporary bounds to make a vertex.
*        ---------------------------------------------------------------
         IF (VERTEX  .AND.  NACTIV+NFIXED .LT. N) THEN
 
*           Compute lengths of columns of selected linear constraints
*           (just the ones corresponding to free variables).
 
            DO 630 J = 1, N
               IF (ISTATE(J) .EQ. 0) THEN
                  COLSIZ = ZERO
                  DO 620 K = 1, NCLIN
                     IF (ISTATE(N+K) .GT. 0)
     $               COLSIZ = COLSIZ + ABS( A(K,J) )
  620             CONTINUE
                  WORK(J) = COLSIZ
               END IF
  630       CONTINUE
 
*           Find the  NARTIF  smallest such columns.
*           This is an expensive loop.  Later we can replace it by a
*           4-pass process (say), accepting the first col that is within
*           T  of  COLMIN, where  T = 0.0, 0.001, 0.01, 0.1 (say).
*           (This comment written in 1980).
 
*+          WHILE (NFIXED + NACTIV .LT. N) DO
  640       IF    (NFIXED + NACTIV .LT. N) THEN
               COLMIN = FLMAX
               DO 650 J = 1, N
                  IF (ISTATE(J) .EQ. 0) THEN
                     IF (NCLIN .EQ. 0) GO TO 660
                     COLSIZ = WORK(J)
                     IF (COLMIN .GT. COLSIZ) THEN
                        COLMIN = COLSIZ
                        JMIN   = J
                     END IF
                  END IF
  650          CONTINUE
               J      = JMIN
  660          ISTATE(J) = 4
               NARTIF = NARTIF + 1
               NFIXED = NFIXED + 1
               GO TO 640
*+          END WHILE
            END IF
         END IF
      END IF
 
      NFREE = N - NFIXED
 
      IF (LSDBG) THEN
         IF (ILSDBG(1) .GT. 0)
     $       WRITE (NOUT, 1300) NFIXED, NACTIV, NARTIF
         IF (ILSDBG(2) .GT. 0)
     $       WRITE (NOUT, 1200) (WX(J), J = 1, N)
      END IF
 
      RETURN
 
 1000 FORMAT(/ ' //LSCRSH// COLD NCLIN NCTOTL'
     $       / ' //LSCRSH// ', L4, I6, I7 )
 1100 FORMAT(/ ' //LSCRSH// Variables before crash... '/ (5G12.3))
 1200 FORMAT(/ ' //LSCRSH// Variables after  crash... '/ (5G12.3))
 1300 FORMAT(/ ' //LSCRSH// Working set selected ...             '
     $       / ' //LSCRSH// NFIXED NACTIV NARTIF      '
     $       / ' //LSCRSH// ', I6, 2I7 )
 
*     End of  LSCRSH.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSDEL ( UNITQ,
     $                   N, NACTIV, NFREE, NRES, NGQ, NZ, NZ1,
     $                   NROWA, NQ, NROWR, NROWT, NRANK,
     $                   JDEL, KDEL, KACTIV, KX,
     $                   A, RES, R, T, GQ, ZY, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            UNITQ
      INTEGER            KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), RES(N,*), R(NROWR,*), T(NROWT,*),
     $                   GQ(N,*), ZY(NQ,*)
      DOUBLE PRECISION   WORK(N)
 
************************************************************************
*     LSDEL   updates the least-squares factor R and the factorization
*     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
*     constraint is deleted from the working set.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSDEL dated 10-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      EXTERNAL           IDAMAX
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      IF (JDEL .GT. 0) THEN
 
*        Regular constraint or temporary bound deleted.
 
         IF (JDEL .LE. N) THEN
 
*           Case 1.  A simple bound has been deleted.
*           =======  Columns NFREE+1 and IR of R must be swapped.
 
            IR     = NZ    + KDEL
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $         WRITE (NOUT, 1100) NACTIV, NZ, NFREE, IR, JDEL, UNITQ
 
            IBEGIN = 1
            NFREE  = NFREE + 1
            IF (NFREE .LT. IR) THEN
               KX(IR)    = KX(NFREE)
               KX(NFREE) = JDEL
               IF (NRANK .GT. 0)
     $            CALL CMRSWP( N, NRES, NRANK, NROWR, NFREE, IR,
     $                         R, RES, WORK )
               CALL DSWAP ( NGQ, GQ(NFREE,1), N, GQ(IR,1), N )
            END IF
 
            IF (.NOT. UNITQ) THEN
 
*              Copy the incoming column of  A(free)  into the end of T.
 
               DO 130 KA = 1, NACTIV
                  I = KACTIV(KA)
                  T(KA,NFREE) = A(I,JDEL)
  130          CONTINUE
 
*              Expand Q by adding a unit row and column.
 
               IF (NFREE .GT. 1) THEN
                  CALL DLOAD ( NFREE-1, ZERO, ZY(NFREE,1), NQ )
                  CALL DLOAD ( NFREE-1, ZERO, ZY(1,NFREE), 1  )
               END IF
               ZY(NFREE,NFREE) = ONE
            END IF
         ELSE
 
*           Case 2.  A general constraint has been deleted.
*           =======
 
            IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $         WRITE (NOUT, 1200) NACTIV, NZ, NFREE, KDEL, JDEL, UNITQ
 
            IBEGIN = KDEL
            NACTIV = NACTIV - 1
 
*           Delete a row of T and move the ones below it up.
 
            DO 220 I = KDEL, NACTIV
               KACTIV(I) = KACTIV(I+1)
               LD        = NFREE - I
               CALL DCOPY ( I+1, T(I+1,LD), NROWT, T(I,LD), NROWT )
  220       CONTINUE
         END IF
 
*        ---------------------------------------------------------------
*        Eliminate the super-diagonal elements of  T,
*        using a backward sweep of 2*2 transformations.
*        ---------------------------------------------------------------
         K     = NFREE  - IBEGIN
         L     = NACTIV - IBEGIN
         LROWR = N      - K
 
         DO 420 I = IBEGIN, NACTIV
            CALL DROT3G( T(I,K+1), T(I,K), CS, SN )
 
            IF (L .GT. 0)
     $      CALL DROT3 ( L    , T(I+1,K+1), 1, T(I+1,K ), 1, CS, SN )
            CALL DROT3 ( NFREE, ZY(1,K+1) , 1, ZY(1,K  ), 1, CS, SN )
            CALL DROT3 ( NGQ  , GQ(K+1,1) , N, GQ(K,1)  , N, CS, SN )
 
*           Apply the column transformations to  R.  The non-zero
*           sub-diagonal that is generated must be eliminated by a row
*           rotation.
 
            IF (K .LT. NRANK) R(K+1,K) = ZERO
            LCOL   = MIN( K+1, NRANK )
            IF (LCOL .GT. 0)
     $         CALL DROT3 ( LCOL, R(1,K+1), 1, R(1,K), 1, CS, SN )
 
            IF (K .LT. NRANK) THEN
               CALL DROT3G( R(K,K), R(K+1,K), CS, SN )
 
               CALL DROT3 ( LROWR, R(K,K+1)    , NROWR,
     $                             R(K+1,K+1)  , NROWR, CS, SN )
               CALL DROT3 ( NRES , RES(K,1)    , N    ,
     $                             RES(K+1,1)  , N    , CS, SN )
            END IF
            K     = K     - 1
            L     = L     - 1
            LROWR = LROWR + 1
  420    CONTINUE
 
         NZ  = NZ  + 1
 
*        ---------------------------------------------------------------
*        Estimate the condition number of  T.
*        ---------------------------------------------------------------
         IF (NACTIV .EQ. 0) THEN
            DTMAX = ONE
            DTMIN = ONE
         ELSE
            CALL DCOND ( NACTIV, T(NACTIV,NZ+1), NROWT-1, DTMAX, DTMIN )
         END IF
 
      END IF
 
      NZ1 = NZ1 + 1
 
      IF (NZ .GT. NZ1) THEN
         IF (JDEL .GT. 0) THEN
            JART =   NZ1 - 1 + IDAMAX( NZ-NZ1+1, GQ(NZ1,1), 1 )
         ELSE
            JART = - JDEL
         END IF
 
         IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $      WRITE( NOUT, 1000 ) NZ, NZ1, JART
 
         IF (JART .GT. NZ1) THEN
 
*           Swap columns NZ1 and JART of R.
 
            IF (UNITQ) THEN
               K        = KX(NZ1)
               KX(NZ1)  = KX(JART)
               KX(JART) = K
            ELSE
               CALL DSWAP ( NFREE, ZY(1,NZ1), 1, ZY(1,JART), 1 )
            END IF
 
            CALL DSWAP ( NGQ, GQ(NZ1,1), N, GQ(JART,1), N )
            IF (NRANK .GT. 0)
     $         CALL CMRSWP( N, NRES, NRANK, NROWR, NZ1, JART,
     $                      R, RES, WORK )
         END IF
      END IF
 
      RETURN
 
 1000 FORMAT(/ ' //LSDEL //  Artificial constraint deleted.      '
     $       / ' //LSDEL //      NZ   NZ1   JART                 '
     $       / ' //LSDEL //  ', 3I6 )
 1100 FORMAT(/ ' //LSDEL //  Simple bound deleted.               '
     $       / ' //LSDEL //  NACTIV    NZ NFREE    IR  JDEL UNITQ'
     $       / ' //LSDEL //  ', 5I6, L6 )
 1200 FORMAT(/ ' //LSDEL //  General constraint deleted.         '
     $       / ' //LSDEL //  NACTIV    NZ NFREE  KDEL  JDEL UNITQ'
     $       / ' //LSDEL //  ', 5I6, L6 )
 
*     End of  LSDEL .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSDFLT( M, N, NCLIN, TITLE )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
      CHARACTER*(*)      TITLE
 
************************************************************************
*  LSDFLT  loads the default values of parameters not set by the user.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original Fortran 77 version written 17-September-1985.
*  This version of LSDFLT dated   9-September-1986.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
 
      LOGICAL            CMDBG, LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
      COMMON    /CMDEBG/ ICMDBG(LDBG), CMDBG
 
      LOGICAL            NEWOPT
      COMMON    /SOL3LS/ NEWOPT
      SAVE      /SOL3LS/
 
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
      EQUIVALENCE   (MSGLS , MSGLVL), (IDBGLS, IDBG), (LDBGLS, MSGDBG)
 
      LOGICAL            CDEFND
      CHARACTER*4        ICRSH(0:2)
      CHARACTER*3        LSTYPE(1:10)
      CHARACTER*16       KEY
      INTRINSIC          LEN    ,  MAX   , MOD
      PARAMETER        ( ZERO   =  0.0D+0, TEN    = 10.0D+0)
      PARAMETER        ( RDUMMY = -11111., IDUMMY = -11111 )
      PARAMETER        ( GIGANT = 1.0D+10*.99999           )
      PARAMETER        ( WRKTOL = 1.0D-2                   )
      DATA               ICRSH(0), ICRSH(1), ICRSH(2)
     $                 /'COLD'   ,'WARM'   ,'HOT '   /
      DATA               LSTYPE(1), LSTYPE(2)
     $                 /' FP'     ,' LP'     /
      DATA               LSTYPE(3), LSTYPE(4), LSTYPE(5), LSTYPE(6)
     $                 /'QP1'     ,'QP2'     ,'QP3'     ,'QP4'     /
      DATA               LSTYPE(7), LSTYPE(8), LSTYPE(9), LSTYPE(10)
     $                 /'LS1'     ,'LS2'     ,'LS3'     ,'LS4'     /
 
      EPSMCH = WMACH( 3)
 
*     Make a dummy call to LSKEY to ensure that the defaults are set.
 
      CALL LSKEY ( NOUT, '*', KEY )
      NEWOPT = .TRUE.
 
*     Save the optional parameters set by the user.  The values in
*     RPRMLS and IPRMLS may be changed to their default values.
 
      CALL ICOPY ( MXPARM, IPRMLS, 1, IPSVLS, 1 )
      CALL DCOPY ( MXPARM, RPRMLS, 1, RPSVLS, 1 )
 
      IF (       LPROB  .LT. 0      )  LPROB   = 7
                                       CDEFND  = LPROB .EQ. 2*(LPROB/2)
      IF (       LCRASH .LT. 0
     $    .OR.   LCRASH .GT. 2      )  LCRASH  = 0
      IF (       ITMAX1 .LT. 0      )  ITMAX1  = MAX(50, 5*(N+NCLIN))
      IF (       ITMAX2 .LT. 0      )  ITMAX2  = MAX(50, 5*(N+NCLIN))
      IF (       MSGLVL .EQ. IDUMMY )  MSGLVL  = 10
      IF (       IDBG   .LT. 0
     $    .OR.   IDBG   .GT. ITMAX1 + ITMAX2
     $                              )  IDBG    = 0
      IF (       MSGDBG .LT. 0      )  MSGDBG  = 0
      IF (       MSGDBG .EQ. 0      )  IDBG    = ITMAX1 + ITMAX2 + 1
      IF (       TOLACT .LT. ZERO   )  TOLACT  = WRKTOL
      IF (       TOLFEA .EQ. RDUMMY
     $    .OR.  (TOLFEA .GE. ZERO
     $    .AND.  TOLFEA .LT. EPSMCH))  TOLFEA  = EPSPT5
      IF (       TOLRNK .LE. ZERO
     $    .AND.  CDEFND             )  TOLRNK  = EPSPT5
      IF (       TOLRNK .LE. ZERO   )  TOLRNK  = TEN*EPSMCH
      IF (       BIGBND .LE. ZERO   )  BIGBND  = GIGANT
      IF (       BIGDX  .LE. ZERO   )  BIGDX   = MAX(GIGANT, BIGBND)
 
      LSDBG = IDBG .EQ. 0
      CMDBG = LSDBG
      K     = 1
      MSG   = MSGDBG
      DO 200 I = 1, LDBG
         ILSDBG(I) = MOD( MSG/K, 10 )
         ICMDBG(I) = ILSDBG(I)
         K = K*10
  200 CONTINUE
 
      IF (MSGLVL .GT. 0) THEN
 
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
         WRITE (NOUT, 2100) LSTYPE(LPROB),
     $                      NCLIN , TOLFEA, ICRSH(LCRASH),
     $                      N     , BIGBND, TOLACT,
     $                      M     , BIGDX , TOLRNK
         WRITE (NOUT, 2200) EPSMCH, ITMAX1, MSGLVL,
     $                              ITMAX2
      END IF
 
      RETURN
 
 2000 FORMAT(
     $//' Parameters'
     $/ ' ----------' )
 2100 FORMAT(
     $/ ' Problem type...........', 7X, A3
     $/ ' Linear constraints.....', I10,     6X,
     $  ' Feasibility tolerance..', 1PE10.2, 6X,
     $  1X, A4, ' start.............'
     $/ ' Variables..............', I10,     6X,
     $  ' Infinite bound size....', 1PE10.2, 6X,
     $  ' Crash tolerance........', 1PE10.2
     $/ ' Objective matrix rows..', I10,     6X,
     $  ' Infinite step size.....', 1PE10.2, 6X,
     $  ' Rank tolerance.........', 1PE10.2 )
 2200 FORMAT(
     $/ ' EPS (machine precision)', 1PE10.2, 6X,
     $  ' Feasibility phase itns.', I10, 6X,
     $  ' Print level............', I10
     $/ 40X,
     $  ' Optimality  phase itns.', I10 )
 
*     End of  LSDFLT.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSFEAS( N, NCLIN, ISTATE,
     $                   BIGBND, CVNORM, ERRMAX, JMAX, NVIOL,
     $                   AX, BL, BU, FEATOL, X, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            ISTATE(N+NCLIN)
      DOUBLE PRECISION   AX(*), BL(N+NCLIN), BU(N+NCLIN)
      DOUBLE PRECISION   FEATOL(N+NCLIN), X(N)
      DOUBLE PRECISION   WORK(N+NCLIN)
 
************************************************************************
*  LSFEAS  computes the following...
*  (1)  The number of constraints that are violated by more
*       than  FEATOL  and the 2-norm of the constraint violations.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version      April    1984.
*  This version of  LSFEAS  dated  17-October-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      EXTERNAL           IDAMAX, DNRM2
      INTRINSIC          ABS
      PARAMETER        ( ZERO = 0.0D+0 )
 
      BIGLOW = - BIGBND
      BIGUPP =   BIGBND
 
*     ==================================================================
*     Compute NVIOL,  the number of constraints violated by more than
*     FEATOL,  and CVNORM,  the 2-norm of the constraint violations and
*     residuals of the constraints in the working set.
*     ==================================================================
      NVIOL  = 0
 
      DO 200 J = 1, N+NCLIN
         FEASJ  = FEATOL(J)
         IS     = ISTATE(J)
         RES    = ZERO
 
         IF (IS .GE. 0  .AND.  IS .LT. 4) THEN
            IF (J .LE. N) THEN
               CON =  X(J)
            ELSE
               I   = J - N
               CON = AX(I)
            END IF
 
            TOLJ   = FEASJ
 
*           Check for constraint violations.
 
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
 
*           This constraint is satisfied,  but count the residual as a
*           violation if the constraint is in the working set.
 
            IF (IS .LE. 0) RES = ZERO
            IF (IS .EQ. 1) RES = BL(J) - CON
            IF (IS .GE. 2) RES = BU(J) - CON
            IF (ABS( RES ) .GT. FEASJ) NVIOL = NVIOL + 1
         END IF
  190    WORK(J) = RES
  200 CONTINUE
 
      JMAX   = IDAMAX( N+NCLIN, WORK, 1 )
      ERRMAX = ABS ( WORK(JMAX) )
 
      IF (LSDBG  .AND.  ILSDBG(1) .GT. 0)
     $   WRITE (NOUT, 1000) ERRMAX, JMAX
 
      CVNORM  = DNRM2 ( N+NCLIN, WORK, 1 )
 
      RETURN
 
 1000 FORMAT(/ ' //LSFEAS//  The maximum violation is ', 1PE14.2,
     $                     ' in constraint', I5 )
 
*     End of  LSFEAS.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSFILE( IOPTNS, INFORM )
      INTEGER            IOPTNS, INFORM
 
************************************************************************
*     LSFILE  reads the options file from unit  IOPTNS  and loads the
*     options into the relevant elements of  IPRMLS  and  RPRMLS.
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
      COMMON     /SOL3LS/ NEWOPT
      SAVE       /SOL3LS/
 
      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/
 
      EXTERNAL            MCHPAR, LSKEY
      LOGICAL             FIRST
      SAVE                FIRST , NOUT
      DATA                FIRST /.TRUE./
 
*     If first time in, set NOUT.
*     NEWOPT is true first time into LSFILE or LSOPTN
*     and just after a call to LSSOL.
 
      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT = WMACH(11)
      END IF
 
      CALL OPFILE( IOPTNS, NOUT, INFORM, LSKEY )
 
      RETURN
 
*     End of  LSFILE.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSGETP( LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                   N, NCLIN, NFREE,
     $                   NROWA, NQ, NROWR, NRANK, NUMINF, NZ1,
     $                   ISTATE, KX, CTP, PNORM,
     $                   A, AP, RES, HZ, P,
     $                   GQ, CQ, R, ZY, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      LOGICAL            LINOBJ, SINGLR, UNITGZ, UNITQ
      INTEGER            ISTATE(N+NCLIN), KX(N)
      DOUBLE PRECISION   A(NROWA,*), AP(*), RES(*), HZ(*), P(N),
     $                   GQ(N), CQ(*), R(NROWR,*), ZY(NQ,*)
      DOUBLE PRECISION   WORK(N)
 
************************************************************************
*     LSGETP  computes the following quantities for  LSCORE.
*     (1) The vector  (hz1) = (Rz1)(pz1).
*         If X is not yet feasible,  the product is computed directly.
*         If  Rz1 is singular,  hz1  is zero.  Otherwise  hz1  satisfies
*         the equations
*                        Rz1'hz1 = -gz1,
*         where  g  is the total gradient.  If there is no linear term
*         in the objective,  hz1  is set to  dz1  directly.
*     (2) The search direction P (and its 2-norm).  The vector P is
*         defined as  Z*(pz1), where  (pz1)  depends upon whether or
*         not X is feasible and the nonsingularity of  (Rz1).
*         If  NUMINF .GT. 0,  (pz1)  is the steepest-descent direction.
*         Otherwise,  x  is the solution of the  NZ1*NZ1  triangular
*         system   (Rz1)*(pz1) = (hz1).
*     (3) The vector Ap,  where A is the matrix of linear constraints.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of LSGETP dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          MIN
      PARAMETER        ( ZERO = 0.0D+0, ONE  = 1.0D+0 )
 
      IF (SINGLR) THEN
*        ---------------------------------------------------------------
*        The triangular factor for the current objective function is
*        singular,  i.e., the objective is linear along the last column
*        of Z1.  This can only occur when UNITGZ is TRUE.
*        ---------------------------------------------------------------
         IF (NZ1 .GT. 1) THEN
            CALL DCOPY ( NZ1-1, R(1,NZ1), 1, P, 1 )
            CALL DTRSV ( 'U', 'N', 'N', NZ1-1, R, NROWR, P, 1 )
         END IF
         P(NZ1) = - ONE
 
         GTP = DDOT  ( NZ1, GQ, 1, P, 1 )
         IF (GTP .GT. ZERO) CALL DSCAL ( NZ1, (-ONE), P, 1 )
 
         IF (NZ1 .LE. NRANK) THEN
            IF (NUMINF .EQ. 0) THEN
               IF (UNITGZ) THEN
                  HZ(NZ1) = R(NZ1,NZ1)*P(NZ1)
               ELSE
                  CALL DLOAD ( NZ1, (ZERO), HZ, 1 )
               END IF
            ELSE
               HZ(1)   = R(1,1)*P(1)
            END IF
         END IF
      ELSE
*        ---------------------------------------------------------------
*        The objective is quadratic in the space spanned by Z1.
*        ---------------------------------------------------------------
         IF (LINOBJ) THEN
            IF (UNITGZ) THEN
               IF (NZ1 .GT. 1)
     $            CALL DLOAD ( NZ1-1, (ZERO), HZ, 1 )
               HZ(NZ1) = - GQ(NZ1)/R(NZ1,NZ1)
            ELSE
               CALL DCOPY ( NZ1, GQ  , 1, HZ, 1 )
               CALL DSCAL ( NZ1, (-ONE), HZ, 1 )
               CALL DTRSV ( 'U', 'T', 'N', NZ1, R, NROWR, HZ, 1 )
            END IF
         ELSE
            CALL DCOPY ( NZ1, RES, 1, HZ, 1 )
         END IF
 
*        Solve  Rz1*pz1 = hz1.
 
         CALL DCOPY ( NZ1, HZ, 1, P, 1 )
         CALL DTRSV ( 'U', 'N', 'N', NZ1, R, NROWR, P, 1 )
      END IF
 
*     Compute  p = Z1*pz1  and its norm.
 
      IF (LINOBJ)
     $   CTP = DDOT  ( NZ1, CQ, 1, P, 1 )
      PNORM  = DNRM2 ( NZ1, P, 1 )
 
      CALL CMQMUL( 1, N, NZ1, NFREE, NQ, UNITQ, KX, P, ZY, WORK )
 
      IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 1000) (P(J), J = 1, N)
 
*     Compute  Ap.
 
      IF (NCLIN .GT. 0) THEN
         CALL DLOAD ( NCLIN, ZERO, AP, 1 )
         DO 410 J = 1, N
            IF (ISTATE(J) .LE. 0)
     $         CALL DAXPY( NCLIN, P(J), A(1,J), 1, AP, 1 )
  410    CONTINUE
         IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 1100) (AP(I), I = 1, NCLIN)
      END IF
 
      RETURN
 
 1000 FORMAT(/ ' //LSGETP//   P ... ' / (1P5E15.5))
 1100 FORMAT(/ ' //LSGETP//  AP ... ' / (1P5E15.5))
 
*     End of  LSGETP.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSGSET( PRBTYP, LINOBJ, SINGLR, UNITGZ, UNITQ,
     $                   N, NCLIN, NFREE,
     $                   NROWA, NQ, NROWR, NRANK, NZ, NZ1,
     $                   ISTATE, KX,
     $                   BIGBND, TOLRNK, NUMINF, SUMINF,
     $                   BL, BU, A, RES, FEATOL,
     $                   GQ, CQ, R, X, WTINF, ZY, WRK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      LOGICAL            LINOBJ, SINGLR, UNITGZ, UNITQ
      INTEGER            ISTATE(*), KX(N)
      DOUBLE PRECISION   BL(*), BU(*), A(NROWA,*),
     $                   RES(*), FEATOL(*)
      DOUBLE PRECISION   GQ(N), CQ(*), R(NROWR,*), X(N), WTINF(*),
     $                   ZY(NQ,*)
      DOUBLE PRECISION   WRK(N)
 
************************************************************************
*     LSGSET  finds the number and weighted sum of infeasibilities for
*     the bounds and linear constraints.   An appropriate transformed
*     gradient vector is returned in  GQ.
*
*     Positive values of  ISTATE(j)  will not be altered.  These mean
*     the following...
*
*               1             2           3
*           a'x = bl      a'x = bu     bl = bu
*
*     Other values of  ISTATE(j)  will be reset as follows...
*           a'x lt bl     a'x gt bu     a'x free
*              - 2           - 1           0
*
*     If  x  is feasible,  LSGSET computes the vector Q(free)'g(free),
*     where  g  is the gradient of the the sum of squares plus the
*     linear term.  The matrix Q is of the form
*                    ( Q(free)  0       ),
*                    (   0      I(fixed))
*     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
*     the matrix of constraints in the working set.  The transformed
*     gradients are stored in GQ.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of LSGSET dated 24-June-1986.
************************************************************************
      EXTERNAL           DDOT  , IDRANK
      INTRINSIC          ABS   , MAX   , MIN
      PARAMETER        ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
 
      BIGUPP =   BIGBND
      BIGLOW = - BIGBND
 
      NUMINF =   0
      SUMINF =   ZERO
      CALL DLOAD ( N, ZERO, GQ, 1 )
 
      DO 200 J = 1, N+NCLIN
         IF (ISTATE(J) .LE. 0) THEN
            FEASJ  = FEATOL(J)
            IF (J .LE. N) THEN
               CTX = X(J)
            ELSE
               K   = J - N
               CTX = DDOT  ( N, A(K,1), NROWA, X, 1 )
            END IF
            ISTATE(J) = 0
 
*           See if the lower bound is violated.
 
            IF (BL(J) .GT. BIGLOW) THEN
               S = BL(J) - CTX
               IF (S     .GT. FEASJ ) THEN
                  ISTATE(J) = - 2
                  WEIGHT    = - WTINF(J)
                  GO TO 160
               END IF
            END IF
 
*           See if the upper bound is violated.
 
            IF (BU(J) .GE. BIGUPP) GO TO 200
            S = CTX - BU(J)
            IF (S     .LE. FEASJ ) GO TO 200
            ISTATE(J) = - 1
            WEIGHT    =   WTINF(J)
 
*           Add the infeasibility.
 
  160       NUMINF = NUMINF + 1
            SUMINF = SUMINF + ABS( WEIGHT ) * S
            IF (J .LE. N) THEN
               GQ(J) = WEIGHT
            ELSE
               CALL DAXPY ( N, WEIGHT, A(K,1), NROWA, GQ, 1 )
            END IF
         END IF
  200 CONTINUE
 
*     ------------------------------------------------------------------
*     Install  GQ,  the transformed gradient.
*     ------------------------------------------------------------------
      SINGLR = .FALSE.
      UNITGZ = .TRUE.
 
      IF (NUMINF .GT. 0) THEN
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ, KX, GQ, ZY, WRK )
      ELSE IF (NUMINF .EQ. 0  .AND.  PRBTYP .EQ. 'FP') THEN
         CALL DLOAD ( N, ZERO, GQ, 1 )
      ELSE
 
*        Ready for the Optimality Phase.
*        Set NZ1 so that Rz1 is nonsingular.
 
         IF (NRANK .EQ. 0) THEN
            IF (LINOBJ) THEN
               CALL DCOPY ( N, CQ, 1, GQ, 1 )
            ELSE
               CALL DLOAD ( N, ZERO, GQ, 1 )
            END IF
            NZ1    = 0
         ELSE
 
*           Compute  GQ = - R' * (transformed residual)
 
            CALL DCOPY ( NRANK, RES, 1, GQ, 1 )
            CALL DSCAL ( NRANK, (-ONE), GQ, 1 )
            CALL DTRMV ( 'U', 'T', 'N', NRANK, R, NROWR, GQ, 1 )
            IF (NRANK .LT. N)
     $         CALL DGEMV( 'T', NRANK, N-NRANK, -ONE,R(1,NRANK+1),NROWR,
     $                      RES, 1, ZERO, GQ(NRANK+1), 1 )
 
            IF (LINOBJ) CALL DAXPY ( N, ONE, CQ, 1, GQ, 1 )
            UNITGZ = .FALSE.
            NZ1    = IDRANK( MIN(NRANK, NZ), R, NROWR+1, TOLRNK )
         END IF
      END IF
 
      RETURN
 
*     End of  LSGSET.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSKEY ( NOUT, BUFFER, KEY )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*)      BUFFER
 
************************************************************************
*     LSKEY   decodes the option contained in  BUFFER  in order to set
*     a parameter value in the relevant element of  IPRMLS  or  RPRMLS.
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
*     LSKEY  calls OPNUMB and the subprograms
*                 LOOKUP, SCANNR, TOKENS, UPCASE
*     (now called OPLOOK, OPSCAN, OPTOKN, OPUPPR)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version dated Jan 22, 1986.
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
 
      EXTERNAL           OPNUMB
      LOGICAL            FIRST , MORE  , NUMBER, OPNUMB, SORTED
      SAVE               FIRST
 
      PARAMETER         (     MAXKEY = 27,  MAXTIE = 10,   MAXTOK = 10,
     $                        MAXTYP = 16)
      CHARACTER*16       KEYS(MAXKEY), TIES(MAXTIE), TOKEN(MAXTOK),
     $                   TYPE(MAXTYP)
      CHARACTER*16       KEY, KEY2, KEY3, VALUE
 
      PARAMETER         (IDUMMY = -11111,  RDUMMY = -11111.0,
     $                   SORTED = .TRUE.,  ZERO   =  0.0     )
 
      DATA                FIRST
     $                  /.TRUE./
      DATA   KEYS
     $ / 'BEGIN           ',
     $   'COLD            ', 'CONSTRAINTS     ', 'CRASH           ',
     $   'DEBUG           ', 'DEFAULTS        ', 'END             ',
     $   'FEASIBILITY     ', 'HOT             ', 'INFINITE        ',
     $   'IPRMLS          ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'LINEAR          ', 'LIST            ',
     $   'LOWER           ', 'NOLIST          ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'RANK            ',
     $   'RPRMLS          ', 'START           ', 'UPPER           ',
     $   'VARIABLES       ', 'WARM            '/
 
      DATA   TIES
     $ / 'BOUND           ', 'CONSTRAINTS     ',
     $   'NO              ', 'NO.      :NUMBER', 'NUMBER          ',
     $   'PHASE           ', 'STEP            ',
     $   'TOLERANCE       ', 'TYPE            ', 'YES             '/
 
      DATA   TYPE
     $ / 'FP              ',
     $   'LEAST       :LS1', 'LINEAR       :LP', 'LP              ',
     $   'LS          :LS1', 'LS1             ', 'LS2             ',
     $   'LS3             ', 'LS4             ', 'LSQ         :LS1',
     $   'QP          :QP2', 'QP1             ', 'QP2             ',
     $   'QP3             ', 'QP4             ', 'QUADRATIC   :QP2'/
*-----------------------------------------------------------------------
 
      IF (FIRST) THEN
         FIRST  = .FALSE.
         DO 10 I = 1, MXPARM
            IPRMLS(I) = IDUMMY
            RPRMLS(I) = RDUMMY
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
         IF (KEY .EQ. 'COLD        ') THEN
            LCRASH = 0
         ELSE IF (KEY .EQ. 'CONSTRAINTS ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'CRASH       ') THEN
            TOLACT = RVALUE
         ELSE IF (KEY .EQ. 'DEBUG       ') THEN
            LDBGLS = RVALUE
         ELSE IF (KEY .EQ. 'DEFAULTS    ') THEN
            DO 20 I = 1, MXPARM
               IPRMLS(I) = IDUMMY
               RPRMLS(I) = RDUMMY
   20       CONTINUE
         ELSE IF (KEY .EQ. 'FEASIBILITY ') THEN
              IF (KEY2.EQ. 'PHASE       ') ITMAX1 = RVALUE
              IF (KEY2.EQ. 'TOLERANCE   ') TOLFEA = RVALUE
              IF (LOC2.EQ.  0            ) WRITE(NOUT, 2320) KEY2
         ELSE
            MORE   = .TRUE.
         END IF
      END IF
 
      IF (MORE) THEN
         MORE   = .FALSE.
         IF (KEY .EQ. 'HOT         ') THEN
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
            ITMAX2 = RVALUE
         ELSE IF (KEY .EQ. 'LINEAR      ') THEN
            NNCLIN = RVALUE
         ELSE IF (KEY .EQ. 'LOWER       ') THEN
            BNDLOW = RVALUE
         ELSE
            MORE   = .TRUE.
         END IF
      END IF
 
      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'OPTIMALITY  ') THEN
            ITMAX2 = RVALUE
         ELSE IF (KEY .EQ. 'PROBLEM     ') THEN
            IF      (KEY2 .EQ. 'NUMBER') THEN
               NPROB  = RVALUE
            ELSE IF (KEY2 .EQ. 'TYPE  ') THEN
 
*              Recognize     Problem type = LP     etc.
 
               CALL OPLOOK( MAXTYP, TYPE, SORTED, KEY3, LOC3 )
               IF (KEY3 .EQ. 'FP' ) LPROB = 1
               IF (KEY3 .EQ. 'LP' ) LPROB = 2
               IF (KEY3 .EQ. 'QP1') LPROB = 3
               IF (KEY3 .EQ. 'QP2') LPROB = 4
               IF (KEY3 .EQ. 'QP3') LPROB = 5
               IF (KEY3 .EQ. 'QP4') LPROB = 6
               IF (KEY3 .EQ. 'LS1') LPROB = 7
               IF (KEY3 .EQ. 'LS2') LPROB = 8
               IF (KEY3 .EQ. 'LS3') LPROB = 9
               IF (KEY3 .EQ. 'LS4') LPROB = 10
               IF (LOC3 .EQ.  0  ) WRITE(NOUT, 2330) KEY3
            ELSE
               WRITE(NOUT, 2320) KEY2
            END IF
         ELSE
            MORE   = .TRUE.
         END IF
      END IF
 
      IF (MORE) THEN
         MORE   = .FALSE.
         IF      (KEY .EQ. 'PRINT       ') THEN
            MSGLS  = RVALUE
         ELSE IF (KEY .EQ. 'RANK        ') THEN
            TOLRNK = RVALUE
         ELSE IF (KEY .EQ. 'RPRMLS      ') THEN
*           Allow things like  RPRMLS 21 = 2  to set RPRMLS(21) = 2.0
            IVALUE = RVALUE
            IF (IVALUE .GE. 1  .AND. IVALUE .LE. MXPARM) THEN
               READ (KEY3, '(BN, E16.0)') RPRMLS(IVALUE)
            ELSE
               WRITE(NOUT, 2400) IVALUE
            END IF
         ELSE IF (KEY .EQ. 'START       ') THEN
            IDBGLS = RVALUE
         ELSE IF (KEY .EQ. 'UPPER       ') THEN
            BNDUPP = RVALUE
         ELSE IF (KEY .EQ. 'VARIABLES   ') THEN
            NN     = RVALUE
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
 
*     End of LSKEY
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSLOC ( LPROB, N, NCLIN, LITOTL, LWTOTL )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
************************************************************************
*     LSLOC   allocates the addresses of the work arrays for  LSCORE.
*
*     Note that the arrays  ( GQ, CQ )  and  ( RES, RES0, HZ )  lie in
*     contiguous areas of workspace.
*     RES, RES0 and HZ are not needed for LP.
*     CQ is defined when the objective has an explicit linear term.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  29-October-1984.
*     This version of LSLOC dated 16-February-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ
 
      PARAMETER        ( LENLS = 20 )
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      LOGICAL            LSDBG
      PARAMETER        ( LDBG = 5 )
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      MINIW     = LITOTL + 1
      MINW      = LWTOTL + 1
 
 
*     Assign array lengths that depend upon the problem dimensions.
 
      IF (NCLIN .EQ. 0) THEN
         LENT      = 0
         LENZY     = 0
      ELSE
         LENT  = NROWT*NCOLT
         LENZY = NQ   *NQ
      END IF
 
      LENCQ  = 0
      IF (LPROB .EQ. 2*(LPROB/2)) LENCQ  = N
      LENRES = 0
      IF (LPROB .GT. 2          ) LENRES = N
 
      LKACTV    = MINIW
      MINIW     = LKACTV + N
 
      LANORM    = MINW
      LAP       = LANORM + NCLIN
      LPX       = LAP    + NCLIN
      LGQ       = LPX    + N
      LCQ       = LGQ    + N
      LRES      = LCQ    + LENCQ
      LRES0     = LRES   + LENRES
      LHZ       = LRES0  + LENRES
      LRLAM     = LHZ    + LENRES
      LT        = LRLAM  + N
      LZY       = LT     + LENT
      LWTINF    = LZY    + LENZY
      LWRK      = LWTINF + N  + NCLIN
      LFEATL    = LWRK   + N  + NCLIN
      MINW      = LFEATL + N  + NCLIN
 
      LOCLS( 1) = LKACTV
      LOCLS( 2) = LANORM
      LOCLS( 3) = LAP
      LOCLS( 4) = LPX
      LOCLS( 5) = LRES
      LOCLS( 6) = LRES0
      LOCLS( 7) = LHZ
      LOCLS( 8) = LGQ
      LOCLS( 9) = LCQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK
      LOCLS(15) = LFEATL
 
      LITOTL    = MINIW - 1
      LWTOTL    = MINW  - 1
 
      RETURN
 
*     End of  LSLOC .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSMOVE( HITCON, HITLOW, LINOBJ, UNITGZ,
     $                   NCLIN, NRANK, NZ1,
     $                   N, NROWR, JADD, NUMINF,
     $                   ALFA, CTP, CTX, XNORM,
     $                   AP, AX, BL, BU, GQ, HZ, P, RES,
     $                   R, X, WORK )
 
      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            HITCON, HITLOW, LINOBJ, UNITGZ
      DOUBLE PRECISION   AP(*), AX(*), BL(*), BU(*), GQ(*), HZ(*),
     $                   P(N), RES(*), R(NROWR,*), X(N)
      DOUBLE PRECISION   WORK(*)
 
************************************************************************
*     LSMOVE  changes X to X + ALFA*P and updates CTX, AX, RES and GQ
*     accordingly.
*
*     If a bound was added to the working set,  move X exactly on to it,
*     except when a negative step was taken (CMALF may have had to move
*     to some other closer constraint.)
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 27-December-1985.
*     Level 2 BLAS added 11-June-1986.
*     This version of LSMOVE dated 11-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      EXTERNAL           DDOT  , DNRM2
      INTRINSIC          ABS   , MIN
      PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )
 
      CALL DAXPY ( N, ALFA, P, 1, X, 1 )
      IF (LINOBJ) CTX = CTX + ALFA*CTP
 
      IF (HITCON  .AND.  JADD .LE. N) THEN
         BND = BU(JADD)
         IF (HITLOW) BND = BL(JADD)
         IF (ALFA .GE. ZERO) X(JADD) = BND
      END IF
      XNORM  = DNRM2 ( N, X, 1 )
 
      IF (NCLIN .GT. 0)
     $   CALL DAXPY ( NCLIN, ALFA, AP, 1, AX, 1 )
 
      IF (NZ1 .LE. NRANK) THEN
         IF (UNITGZ) THEN
            RES(NZ1) = RES(NZ1) - ALFA*HZ(NZ1)
         ELSE
            CALL DAXPY ( NZ1, (-ALFA), HZ, 1, RES, 1  )
         END IF
 
         IF (NUMINF .EQ. 0) THEN
 
*           Update the transformed gradient GQ so that
*           GQ = GQ + ALFA*R'( HZ ).
*                            ( 0  )
 
            IF (UNITGZ) THEN
               CALL DAXPY ( N-NZ1+1, ALFA*HZ(NZ1), R(NZ1,NZ1), NROWR,
     $                                             GQ(NZ1)   , 1      )
            ELSE
               CALL DCOPY ( NZ1, HZ, 1, WORK, 1 )
               CALL DTRMV ( 'U', 'T', 'N', NZ1, R, NROWR, WORK, 1 )
               IF (NZ1 .LT. N)
     $            CALL DGEMV ( 'T', NZ1, N-NZ1, ONE, R(1,NZ1+1), NROWR,
     $                         HZ, 1, ZERO, WORK(NZ1+1), 1 )
               CALL DAXPY ( N, ALFA, WORK, 1, GQ, 1 )
            END IF
         END IF
      END IF
 
      RETURN
 
*     End of  LSMOVE.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSMULS( PRBTYP,
     $                   MSGLVL, N, NACTIV, NFREE,
     $                   NROWA, NROWT, NUMINF, NZ, NZ1,
     $                   ISTATE, KACTIV, KX, DINKY,
     $                   JSMLST, KSMLST, JINF, JTINY,
     $                   JBIGST, KBIGST, TRULAM,
     $                   A, ANORMS, GQ, RLAMDA, T, WTINF )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      INTEGER            ISTATE(*), KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), ANORMS(*),
     $                   GQ(N), RLAMDA(N), T(NROWT,*), WTINF(*)
 
************************************************************************
*     LSMULS  first computes the Lagrange multiplier estimates for the
*     given working set.  It then determines the values and indices of
*     certain significant multipliers.  In this process, the multipliers
*     for inequalities at their upper bounds are adjusted so that a
*     negative multiplier for an inequality constraint indicates non-
*     optimality.  All adjusted multipliers are scaled by the 2-norm
*     of the associated constraint row.  In the following, the term
*     minimum refers to the ordering of numbers on the real line,  and
*     not to their magnitude.
*
*     JSMLST  is the index of the minimum of the set of adjusted
*             multipliers with values less than  - DINKY.  A negative
*             JSMLST defines the index in Q'g of the artificial
*             constraint to be deleted.
*     KSMLST  marks the position of general constraint JSMLST in KACTIV.
*
*     JBIGST  is the index of the largest of the set of adjusted
*             multipliers with values greater than (1 + DINKY).
*     KBIGST  marks its position in KACTIV.
*
*     On exit,  elements 1 thru NACTIV of RLAMDA contain the unadjusted
*     multipliers for the general constraints.  Elements NACTIV onwards
*     of RLAMDA contain the unadjusted multipliers for the bounds.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of LSMULS dated  30-June-1986.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      INTRINSIC          ABS, MIN
      PARAMETER        ( ZERO   =0.0D+0,ONE    =1.0D+0 )
 
      NFIXED =   N - NFREE
 
      JSMLST =   0
      KSMLST =   0
      SMLLST = - DINKY
 
      TINYLM =   DINKY
      JTINY  =   0
 
      JBIGST =   0
      KBIGST =   0
      BIGGST =   ONE + DINKY
 
      IF (NZ1 .LT. NZ) THEN
*        ---------------------------------------------------------------
*        Compute JSMLST for the artificial constraints.
*        ---------------------------------------------------------------
         DO 100 J = NZ1+1, NZ
            RLAM = - ABS( GQ(J) )
            IF (RLAM .LT. SMLLST) THEN
               SMLLST =   RLAM
               JSMLST = - J
            ELSE IF (RLAM .LT. TINYLM) THEN
               TINYLM =   RLAM
               JTINY  =   J
            END IF
  100    CONTINUE
 
         IF (MSGLVL .GE. 20)
     $      WRITE (NOUT, 1000) (GQ(K), K=NZ1+1,NZ)
 
      END IF
 
*     ------------------------------------------------------------------
*     Compute JSMLST for regular constraints and temporary bounds.
*     ------------------------------------------------------------------
*     First, compute the Lagrange multipliers for the general
*     constraints in the working set, by solving  T'*lamda = Y'g.
 
      IF (N .GT. NZ)
     $   CALL DCOPY ( N-NZ, GQ(NZ+1), 1, RLAMDA, 1 )
      IF (NACTIV .GT. 0)
     $   CALL CMTSOL( 2, NROWT, NACTIV, T(1,NZ+1), RLAMDA )
 
*     -----------------------------------------------------------------
*     Now set elements NACTIV, NACTIV+1,... of  RLAMDA  equal to
*     the multipliers for the bound constraints.
*     -----------------------------------------------------------------
      DO 190 L = 1, NFIXED
         J     = KX(NFREE+L)
         BLAM  = RLAMDA(NACTIV+L)
         DO 170 K = 1, NACTIV
            I    = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(K)
  170    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
  190 CONTINUE
 
*     -----------------------------------------------------------------
*     Find JSMLST and KSMLST.
*     -----------------------------------------------------------------
      DO 330 K = 1, N - NZ
         IF (K .GT. NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(K) + N
         END IF
 
         IS   = ISTATE(J)
 
         I    = J - N
         IF (J .LE. N) ANORMJ = ONE
         IF (J .GT. N) ANORMJ = ANORMS(I)
 
         RLAM = RLAMDA(K)
 
*        Change the sign of the estimate if the constraint is in
*        the working set at its upper bound.
 
         IF (IS .EQ. 2) RLAM =      - RLAM
         IF (IS .EQ. 3) RLAM =   ABS( RLAM )
         IF (IS .EQ. 4) RLAM = - ABS( RLAM )
 
         IF (IS .NE. 3) THEN
            SCDLAM = RLAM * ANORMJ
            IF      (SCDLAM .LT. SMLLST) THEN
               SMLLST = SCDLAM
               JSMLST = J
               KSMLST = K
            ELSE IF (SCDLAM .LT. TINYLM) THEN
               TINYLM = SCDLAM
               JTINY  = J
            END IF
         END IF
 
         IF (NUMINF .GT. 0  .AND.  J .GT. JINF) THEN
            SCDLAM = RLAM/WTINF(J)
            IF (SCDLAM .GT. BIGGST) THEN
               BIGGST = SCDLAM
               TRULAM = RLAMDA(K)
               JBIGST = J
               KBIGST = K
            END IF
         END IF
  330 CONTINUE
 
*     -----------------------------------------------------------------
*     If required, print the multipliers.
*     -----------------------------------------------------------------
      IF (MSGLVL .GE. 20) THEN
         IF (NFIXED .GT. 0)
     $      WRITE (NOUT, 1100) PRBTYP, (KX(NFREE+K),
     $                         RLAMDA(NACTIV+K), K=1,NFIXED)
         IF (NACTIV .GT. 0)
     $      WRITE (NOUT, 1200) PRBTYP, (KACTIV(K),
     $                         RLAMDA(K), K=1,NACTIV)
      END IF
 
      IF (LSDBG  .AND.  ILSDBG(1) .GT. 0) THEN
         WRITE (NOUT, 9000) JSMLST, SMLLST, KSMLST
         WRITE (NOUT, 9100) JBIGST, BIGGST, KBIGST
         WRITE (NOUT, 9200) JTINY , TINYLM
      END IF
 
      RETURN
 
 1000 FORMAT(/ ' Multipliers for the artificial constraints        '
     $       / 4(5X, 1PE11.2))
 1100 FORMAT(/ ' Multipliers for the ', A2, ' bound  constraints   '
     $       / 4(I5, 1PE11.2))
 1200 FORMAT(/ ' Multipliers for the ', A2, ' linear constraints   '
     $       / 4(I5, 1PE11.2))
 9000 FORMAT(/ ' //LSMULS//  JSMLST     SMLLST     KSMLST (Scaled) '
     $       / ' //LSMULS//  ', I6, 1PE11.2, 5X, I6 )
 9100 FORMAT(  ' //LSMULS//  JBIGST     BIGGST     KBIGST (Scaled) '
     $       / ' //LSMULS//  ', I6, 1PE11.2, 5X, I6 )
 9200 FORMAT(  ' //LSMULS//   JTINY     TINYLM                     '
     $       / ' //LSMULS//  ', I6, 1PE11.2)
 
*     End of  LSMULS.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSOPTN( STRING )
      CHARACTER*(*)      STRING
 
************************************************************************
*     LSOPTN  loads the option supplied in  STRING  into the relevant
*     element of  IPRMLS  or  RPRMLS.
************************************************************************
 
      LOGICAL             NEWOPT
      COMMON     /SOL3LS/ NEWOPT
      SAVE       /SOL3LS/
 
      DOUBLE PRECISION    WMACH(15)
      COMMON     /SOLMCH/ WMACH
      SAVE       /SOLMCH/
 
      EXTERNAL            MCHPAR
      CHARACTER*16        KEY
      CHARACTER*72        BUFFER
      LOGICAL             FIRST , PRNT
      SAVE                FIRST , NOUT  , PRNT
      DATA                FIRST /.TRUE./
 
*     If first time in, set  NOUT.
*     NEWOPT  is true first time into  LSFILE  or  LSOPTN
*     and just after a call to  LSSOL.
*     PRNT    is set to true whenever  NEWOPT  is true.
 
      IF (FIRST) THEN
         FIRST  = .FALSE.
         NEWOPT = .TRUE.
         CALL MCHPAR()
         NOUT   =  WMACH(11)
      END IF
      BUFFER = STRING
 
*     Call  LSKEY   to decode the option and set the parameter value.
*     If NEWOPT is true, reset PRNT and test specially for NOLIST.
 
      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT   = .TRUE.
         CALL LSKEY ( NOUT, BUFFER, KEY )
 
         IF (KEY .EQ. 'NOLIST') THEN
            PRNT   = .FALSE.
         ELSE
            WRITE (NOUT, '(// A / A /)')
     $         ' Calls to LSOPTN',
     $         ' ---------------'
            WRITE (NOUT, '( 6X, A )') BUFFER
         END IF
      ELSE
         IF (PRNT)
     $      WRITE (NOUT, '( 6X, A )') BUFFER
         CALL LSKEY ( NOUT, BUFFER, KEY )
 
         IF (KEY .EQ.   'LIST') PRNT = .TRUE.
         IF (KEY .EQ. 'NOLIST') PRNT = .FALSE.
      END IF
 
      RETURN
 
*     End of  LSOPTN.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSPRT ( PRBTYP, PRNT1, ISDEL, ITER, JADD, JDEL,
     $                   MSGLVL, NACTIV, NFREE, N, NCLIN,
     $                   NRANK, NROWR, NROWT, NZ, NZ1, ISTATE,
     $                   ALFA, CONDRZ, CONDT, GFNORM, GZNORM, GZ1NRM,
     $                   NUMINF, SUMINF, CTX, SSQ,
     $                   AX, R, T, X, WORK )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*2        PRBTYP
      LOGICAL            PRNT1
      INTEGER            ISTATE(*)
      DOUBLE PRECISION   AX(*), R(NROWR,*), T(NROWT,*), X(N)
      DOUBLE PRECISION   WORK(N)
 
************************************************************************
*  LSPRT  prints various levels of output for  LSCORE.
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
*        ge  20        constraint status,  x  and  Ax.
*
*        ge  30        diagonals of  T  and  R.
*
*
*  Debug printing is performed depending on the logical variable  LSDBG.
*  LSDBG  is set true when  IDBG  major iterations have been performed.
*  At this point,  printing is done according to a string of binary
*  digits of the form  SVT  (stored in the integer array  ILSDBG).
*
*  S  set 'on'  gives information from the maximum step routine  CMALF.
*  V  set 'on'  gives various vectors in  LSCORE  and its auxiliaries.
*  T  set 'on'  gives a trace of which routine was called and an
*               indication of the progress of the run.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  This version of LSPRT dated 14-January-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      CHARACTER*2        LADD, LDEL
      CHARACTER*2        LSTATE(0:5)
      DATA               LSTATE(0), LSTATE(1), LSTATE(2)
     $                  /'  '     , 'L '     , 'U '     /
      DATA               LSTATE(3), LSTATE(4), LSTATE(5)
     $                  /'E '     , 'T '     , 'Z '     /
 
      IF (MSGLVL .GE. 15) WRITE (NOUT, 1000) PRBTYP, ITER
 
      IF (MSGLVL .GE. 5) THEN
         IF      (JDEL .GT. 0) THEN
            KDEL =   ISDEL
         ELSE IF (JDEL .LT. 0) THEN
            JDEL = - JDEL
            KDEL =   5
         ELSE
            KDEL =   0
         END IF
 
         IF (JADD .GT. 0) THEN
            KADD = ISTATE(JADD)
         ELSE
            KADD = 0
         END IF
 
         LDEL   = LSTATE(KDEL)
         LADD   = LSTATE(KADD)
 
         IF (NUMINF .GT. 0) THEN
            OBJ    = SUMINF
         ELSE
            OBJ    = SSQ + CTX
         END IF
 
*        ---------------------------------------------------------------
*        Print the terse line.
*        ---------------------------------------------------------------
         IF (NRANK .EQ. 0) THEN
            IF (PRNT1  .OR.  MSGLVL .GE. 15) WRITE (NOUT, 1100)
            WRITE (NOUT, 1200) ITER, JDEL, LDEL, JADD, LADD,
     $                         ALFA, NUMINF, OBJ, N-NFREE, NACTIV,
     $                         NZ, NZ1, GFNORM, GZ1NRM, CONDT
         ELSE
            IF (PRNT1  .OR.  MSGLVL .GE. 15) WRITE (NOUT, 1110)
            WRITE (NOUT, 1200) ITER, JDEL, LDEL, JADD, LADD,
     $                         ALFA, NUMINF, OBJ, N-NFREE, NACTIV,
     $                         NZ, NZ1, GFNORM, GZ1NRM, CONDT, CONDRZ
         END IF
 
         IF (MSGLVL .GE. 20) THEN
            WRITE (NOUT, 2000) PRBTYP
            WRITE (NOUT, 2100) (X(J) , ISTATE(J)  ,  J=1,N)
            IF (NCLIN .GT. 0)
     $      WRITE (NOUT, 2200) (AX(K), ISTATE(N+K), K=1,NCLIN )
 
            IF (MSGLVL .GE. 30) THEN
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               IF (NACTIV .GT. 0) THEN
                  CALL DCOPY ( NACTIV, T(NACTIV,NZ+1), NROWT-1, WORK,1 )
                  WRITE (NOUT, 3000) PRBTYP, (WORK(J), J=1,NACTIV)
               END IF
               IF (NRANK  .GT. 0)
     $            WRITE (NOUT, 3100) PRBTYP, (R(J,J) , J=1,NRANK )
            END IF
            WRITE (NOUT, 5000)
         END IF
      END IF
 
      PRNT1 = .FALSE.
 
      RETURN
 
 1000 FORMAT(/// ' ', A2, ' iteration', I5
     $         / ' =================' )
 1100 FORMAT(// '  Itn Jdel  Jadd      Step',
     $          ' Ninf  Sinf/Objective', '  Bnd', '  Lin', '    Nz',
     $          '   Nz1   Norm Gf  Norm Gz1   Cond T' )
 1110 FORMAT(// '  Itn Jdel  Jadd      Step',
     $          ' Ninf  Sinf/Objective', '  Bnd', '  Lin', '    Nz',
     $          '   Nz1   Norm Gf  Norm Gz1   Cond T Cond Rz1' )
 1200 FORMAT(I5, I5, A1, I5, A1, 1PE9.1, I5, 1X, 1PE15.6, 2I5,
     $       2I6, 1P2E10.2, 1P2E9.1 )
 2000 FORMAT(/ ' Values and status of the ', A2, ' constraints'
     $       / ' ---------------------------------------' )
 2100 FORMAT(/ ' Variables...'                 /   (1X, 5(1PE15.6, I5)))
 2200 FORMAT(/ ' General linear constraints...'/   (1X, 5(1PE15.6, I5)))
 3000 FORMAT(/ ' Diagonals of ' , A2,' working set factor T'/(1P5E15.6))
 3100 FORMAT(/ ' Diagonals of ' , A2, ' triangle R         '/(1P5E15.6))
 5000 FORMAT(/// ' ---------------------------------------------------',
     $           '--------------------------------------------' )
 
*     End of  LSPRT .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSSETX( LINOBJ, ROWERR, UNITQ,
     $                   NCLIN, NACTIV, NFREE, NRANK, NZ,
     $                   N, NCTOTL, NQ, NROWA, NROWR, NROWT,
     $                   ISTATE, KACTIV, KX,
     $                   JMAX, ERRMAX, CTX, XNORM,
     $                   A, AX, BL, BU, CQ, RES, RES0, FEATOL,
     $                   R, T, X, ZY, P, WORK )
 
      IMPLICIT           DOUBLE PRECISION (A-H,O-Z)
      LOGICAL            LINOBJ, ROWERR, UNITQ
      INTEGER            ISTATE(NCTOTL), KACTIV(N), KX(N)
      DOUBLE PRECISION   A(NROWA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     $                   CQ(*), RES(*), RES0(*), FEATOL(NCTOTL), P(N),
     $                   R(NROWR,*), T(NROWT,*), ZY(NQ,*), X(N)
      DOUBLE PRECISION   WORK(NCTOTL)
 
************************************************************************
*  LSSETX  computes the point on a working set that is closest to the
*  input vector  x  (in the least-squares sense).  The norm of  x, the
*  transformed residual vector  Pr - RQ'x,  and the constraint values
*  Ax  are also initialized.
*
*  If the computed point gives a row error of more than the feasibility
*  tolerance, an extra step of iterative refinement is used.  If  x  is
*  still infeasible,  the logical variable  ROWERR  is set.
*
*  Systems Optimization Laboratory, Stanford University.
*  Original version written 31-October-1984.
*  This version dated 29-December-1985.
************************************************************************
      COMMON    /SOL1CM/ NOUT
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
 
      EXTERNAL           IDAMAX, DDOT
      INTRINSIC          ABS, MIN
      PARAMETER        ( NTRY  = 2 )
      PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )
 
*     ------------------------------------------------------------------
*     Move  x  onto the simple bounds in the working set.
*     ------------------------------------------------------------------
      DO 100 K = NFREE+1, N
          J   = KX(K)
          IS  = ISTATE(J)
          BND = BL(J)
          IF (IS .GE. 2) BND  = BU(J)
          IF (IS .NE. 4) X(J) = BND
  100 CONTINUE
 
*     ------------------------------------------------------------------
*     Move  x  onto the general constraints in the working set.
*     We shall make  ntry  tries at getting acceptable row errors.
*     ------------------------------------------------------------------
      KTRY   = 1
      JMAX   = 1
      ERRMAX = ZERO
 
*     REPEAT
  200    IF (NACTIV .GT. 0) THEN
 
*           Set  work = residuals for constraints in the working set.
*           Solve for p, the smallest correction to x that gives a point
*           on the constraints in the working set.  Define  p = Y*(py),
*           where  py  solves the triangular system  T*(py) = residuals.
 
            DO 220 I = 1, NACTIV
               K   = KACTIV(I)
               J   = N + K
               BND = BL(J)
               IF (ISTATE(J) .EQ. 2) BND = BU(J)
               WORK(I) = BND - DDOT  ( N, A(K,1), NROWA, X, 1 )
  220       CONTINUE
 
            CALL CMTSOL( 1, NROWT, NACTIV, T(1,NZ+1), WORK )
            CALL DLOAD ( N, ZERO, P, 1 )
            CALL DCOPY ( NACTIV, WORK, 1, P(NZ+1), 1 )
 
            CALL CMQMUL( 2, N, NZ, NFREE, NQ, UNITQ, KX, P, ZY, WORK )
            CALL DAXPY ( N, ONE, P, 1, X, 1 )
         END IF
 
*        ---------------------------------------------------------------
*        Compute the 2-norm of  x.
*        Initialize  Ax  for all the general constraints.
*        ---------------------------------------------------------------
         XNORM  = DNRM2 ( N, X, 1 )
         IF (NCLIN .GT. 0)
     $      CALL DGEMV ( 'N', NCLIN, N, ONE, A, NROWA,
     $                   X, 1, ZERO, AX, 1 )
 
*        ---------------------------------------------------------------
*        Check the row residuals.
*        ---------------------------------------------------------------
         IF (NACTIV .GT. 0) THEN
            DO 300 K = 1, NACTIV
               I   = KACTIV(K)
               J   = N + I
               IS  = ISTATE(J)
               IF (IS .EQ. 1) WORK(K) = BL(J) - AX(I)
               IF (IS .GE. 2) WORK(K) = BU(J) - AX(I)
  300       CONTINUE
 
            JMAX   = IDAMAX( NACTIV, WORK, 1 )
            ERRMAX = ABS( WORK(JMAX) )
         END IF
 
         KTRY = KTRY + 1
*     UNTIL    (ERRMAX .LE. FEATOL(JMAX) .OR. KTRY .GT. NTRY
      IF (.NOT.(ERRMAX .LE. FEATOL(JMAX) .OR. KTRY .GT. NTRY)) GO TO 200
 
      ROWERR = ERRMAX .GT. FEATOL(JMAX)
 
*     ==================================================================
*     Compute the linear objective value  c'x  and the transformed
*     residual  Pr  -  RQ'x = RES0  -  RQ'x.
*     ==================================================================
      IF (NRANK .GT. 0  .OR.  LINOBJ) THEN
         CALL DCOPY ( N, X, 1, P, 1 )
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ, KX, P, ZY, WORK )
      END IF
 
      CTX = ZERO
      IF (LINOBJ)
     $   CTX = DDOT  ( N, CQ, 1, P, 1 )
 
      IF (NRANK .GT. 0) THEN
 
         CALL DTRMV ( 'U', 'N', 'N', NRANK, R, NROWR, P, 1 )
         IF (NRANK .LT. N)
     $      CALL DGEMV ( 'N', NRANK, N-NRANK, ONE, R(1,NRANK+1), NROWR,
     $                   P(NRANK+1), 1, ONE, P, 1 )
 
         CALL DCOPY ( NRANK,       RES0, 1, RES, 1 )
         CALL DAXPY ( NRANK, -ONE, P   , 1, RES, 1 )
 
      END IF
 
      IF (LSDBG  .AND.  ILSDBG(2) .GT. 0)
     $   WRITE (NOUT, 2200) (X(J), J = 1, N)
 
      RETURN
 
 2200 FORMAT(/ ' //LSSETX// Variables after refinement ... '/ (5G12.3))
 
*     End of  LSSETX.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE LSSOL ( MM, N,
     $                   NCLIN, NROWA, NROWR,
     $                   A, BL, BU, CVEC,
     $                   ISTATE, KX, X, R, B,
     $                   INFORM, ITER, OBJ, CLAMDA,
     $                   IW, LENIW, W, LENW )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            LENIW, LENW
      INTEGER            ISTATE(N+NCLIN), KX(N)
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   BL(N+NCLIN), BU(N+NCLIN), A(NROWA,*)
      DOUBLE PRECISION   CLAMDA(N+NCLIN), CVEC(*)
      DOUBLE PRECISION   R(NROWR,*), X(N), B(*)
      DOUBLE PRECISION   W(LENW)
 
************************************************************************
*  LSSOL  solves problems of the form
*
*           Minimize               F(x)
*              x
*                                 (  x )
*           subject to    bl  .le.(    ).ge.  bu,
*                                 ( Ax )
*
*  where  '  denotes the transpose of a column vector,  x  denotes the
*  n-vector of parameters and  F(x) is one of the following functions..
*
*  FP =  None                         (find a feasible point).
*  LP =  c'x
*  QP1=        1/2 x'Rx                R  n times n, symmetric pos. def.
*  QP2=  c'x + 1/2 x'Rx                .  .   ..        ..       ..  ..
*  QP3=        1/2 x'R'Rx              R  m times n, upper triangular.
*  QP4=  c'x + 1/2 x'R'Rx              .  .   ..  .   ..      ...
*  LS1=        1/2 (b - Rx)'(b - Rx)   R  m times n, rectangular.
*  LS2=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .     ...
*  LS3=        1/2 (b - Rx)'(b - Rx)   R  m times n, upper triangular.
*  LS4=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .   ..      ...
*
*  The matrix  R  is entered as the two-dimensional array  R  (of row
*  dimension  NROWR).  If  NROWR = 0,  R  is not accessed.
*
*  The vector  c  is entered in the one-dimensional array  CVEC.
*
*  NCLIN  is the number of general linear constraints (rows of  A).
*  (NCLIN may be zero.)
*
*  The first  N  components of  BL  and   BU  are lower and upper
*  bounds on the variables.  The next  NCLIN  components are
*  lower and upper bounds on the general linear constraints.
*
*  The matrix  A  of coefficients in the general linear constraints
*  is entered as the two-dimensional array  A  (of dimension
*  NROWA by N).  If NCLIN = 0, A is not accessed.
*
*  The vector  x  must contain an initial estimate of the solution,
*  and will contain the computed solution on output.
*
*
*  Complete documentation for  LSSOL  is contained in Report SOL 86-1,
*  Users Guide for LSSOL (Version 1.0), by P.E. Gill, S. J. Hammarling,
*  W. Murray, M.A. Saunders and M.H. Wright, Department of
*  Operations Research, Stanford University, Stanford, California 94305.
*
*  Systems Optimization Laboratory, Stanford University.
*  Version 1.01 Dated  30-June-1986.
*
*  Copyright  1984  Stanford University.
*
*  This material may be reproduced by or for the U.S. Government pursu-
*  ant to the copyright license under DAR clause 7-104.9(a) (1979 Mar).
*
*  This material is based upon work partially supported by the National
*  Science Foundation under Grants MCS-7926009 and ECS-8312142; the
*  Department of Energy Contract AM03-76SF00326, PA No. DE-AT03-
*  76ER72018; the Army Research Office Contract DAA29-84-K-0156;
*  and the Office of Naval Research Grant N00014-75-C-0267.
************************************************************************
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
      COMMON    /SOL1CM/ NOUT
      COMMON    /SOL3CM/ LENNAM, NROWT, NCOLT, NQ
      COMMON    /SOL4CM/ EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON    /SOL5CM/ ASIZE, DTMAX, DTMIN
 
      PARAMETER         (LENLS = 20)
      COMMON    /SOL1LS/ LOCLS(LENLS)
 
      LOGICAL            LSDBG
      PARAMETER         (LDBG = 5)
      COMMON    /LSDEBG/ ILSDBG(LDBG), LSDBG
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
      EQUIVALENCE   (MSGLS , MSGLVL), (IDBGLS, IDBG), (LDBGLS, MSGDBG)
 
      INTRINSIC          MAX, MIN
 
*     Local variables.
 
      LOGICAL            COLD  , FACTRZ, LINOBJ, NAMED , ROWERR,
     $                   UNITQ , VERTEX
      CHARACTER*2        PRBTYP
      CHARACTER*8        NAMES(1)
      PARAMETER        ( ZERO   =0.0D+0, POINT1 =0.1D+0, POINT3 =3.3D-1)
      PARAMETER        ( POINT8 =0.8D+0, POINT9 =0.9D+0, ONE    =1.0D+0)
 
      CHARACTER*40       TITLE
      DATA               TITLE
     $                 / 'SOL/LSSOL  ---  Version 1.01   June 1986' /
 
*     Set the machine-dependent constants.
 
      CALL MCHPAR()
 
      EPSMCH = WMACH( 3)
      RTEPS  = WMACH( 4)
      NOUT   = WMACH(11)
 
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
 
      NAMED  = .FALSE.
 
      INFORM = 0
      ITER   = 0
 
      CONDMX = ONE / EPSPT5
 
      NCTOTL = N + NCLIN
 
*     Set the default values of the parameters.
 
      CALL LSDFLT( MM, N, NCLIN, TITLE )
 
*     Set all parameters determined by the problem type.
 
      IF      (LPROB .EQ. 1 ) THEN
         PRBTYP    = 'FP'
         M      = 0
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 2 ) THEN
         PRBTYP    = 'LP'
         M      = 0
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 3 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 4 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 5 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB .EQ. 6 ) THEN
         PRBTYP    = 'QP'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB .EQ. 7 ) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 8 ) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB .EQ. 9 ) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .FALSE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB .EQ. 10) THEN
         PRBTYP    = 'LS'
         M      = MM
         LINOBJ = .TRUE.
         FACTRZ = .FALSE.
      END IF
 
*     Assign the dimensions of arrays in the parameter list of LSCORE.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter MINACT and MINFXD
*     accordingly.
*     If a linear program is being solved and the matrix of general
*     constraints is fat,  i.e.,  NCLIN .LT. N,  a non-zero value is
*     known for MINFXD.  Note that in this case, VERTEX must be
*     set  .TRUE..
 
      MINACT = 0
      MINFXD = 0
 
      VERTEX = .FALSE.
      IF (      (PRBTYP .EQ. 'LP'  .OR.  PRBTYP .EQ. 'FP')
     $    .AND.  NCLIN  .LT. N   ) THEN
         MINFXD = N - NCLIN - 1
         VERTEX = .TRUE.
      END IF
 
      MXFREE = N - MINFXD
      MAXACT = MAX( 1, MIN( N, NCLIN ) )
      MAXNZ  = N - ( MINFXD + MINACT )
 
      IF (NCLIN .EQ. 0) THEN
         NQ     = 1
         NROWT  = 1
         NCOLT  = 1
         VERTEX = .FALSE.
      ELSE
         NQ     = MAX( 1, MXFREE )
         NROWT  = MAX( MAXNZ, MAXACT )
         NCOLT  = MXFREE
      END IF
 
      NCNLN  = 0
      LENNAM = 1
 
*     Allocate certain arrays that are not done in LSLOC.
 
      LITOTL = 0
 
      LAX    = 1
      LWTOTL = LAX + NCLIN  - 1
 
*     Allocate remaining work arrays.
 
      CALL LSLOC ( LPROB, N, NCLIN, LITOTL, LWTOTL )
 
      COLD  = LCRASH .EQ. 0
 
*     Check input parameters and storage limits.
 
      CALL CMCHK ( NERROR, MSGLVL, COLD, (.NOT.FACTRZ),
     $             LENIW, LENW, LITOTL, LWTOTL,
     $             N, NCLIN, NCNLN,
     $             ISTATE, KX, NAMED, NAMES, LENNAM,
     $             BL, BU, X )
 
      IF (NERROR .GT. 0) THEN
         INFORM = 6
         GO TO 800
      END IF
 
      LKACTV = LOCLS( 1)
 
      LANORM = LOCLS( 2)
      LPX    = LOCLS( 4)
      LRES   = LOCLS( 5)
      LRES0  = LOCLS( 6)
      LGQ    = LOCLS( 8)
      LCQ    = LOCLS( 9)
      LRLAM  = LOCLS(10)
      LT     = LOCLS(11)
      LZY    = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK   = LOCLS(14)
      LFEATL = LOCLS(15)
 
      IF (TOLFEA .GT. ZERO)
     $   CALL DLOAD ( N+NCLIN, (TOLFEA), W(LFEATL), 1 )
 
      IANRMJ = LANORM
      DO 200 J = 1, NCLIN
         W(IANRMJ) = DNRM2 ( N, A(J,1), NROWA )
         IANRMJ    = IANRMJ + 1
  200 CONTINUE
      IF (NCLIN .GT. 0)
     $   CALL DCOND ( NCLIN, W(LANORM), 1, ASIZE, AMIN )
 
      CALL DCOND ( NCTOTL, W(LFEATL), 1, FEAMAX, FEAMIN )
      CALL DCOPY ( NCTOTL, W(LFEATL), 1, W(LWTINF), 1 )
      CALL DSCAL ( NCTOTL, (ONE/FEAMIN), W(LWTINF), 1 )
 
      SSQ1   = ZERO
 
      IF (FACTRZ) THEN
*        ===============================================================
*        Factorize R using QR or Cholesky.  KX must be initialized.
*        ===============================================================
         DO 210 I = 1, N
            KX(I) = I
  210    CONTINUE
 
         IF      (PRBTYP .EQ. 'LP'  .OR.  PRBTYP .EQ. 'FP') THEN
            NRANK = 0
         ELSE IF (PRBTYP .EQ. 'QP') THEN
*           ------------------------------------------------------------
*           Compute the Cholesky factorization of R.  The Hessian is
*           M times M and resides in the upper left-hand corner of R.
*           ------------------------------------------------------------
            DO 220 J = M+1, N
               CALL DLOAD ( M, (ZERO), R(1,J), 1 )
  220       CONTINUE
 
            CALL LSCHOL( NROWR, M, NRANK, TOLRNK, KX, R, INFO )
 
            IF (NRANK .GT. 0)
     $         CALL DLOAD ( NRANK, (ZERO), W(LRES0), 1 )
 
         ELSE IF (PRBTYP .EQ. 'LS') THEN
*           ------------------------------------------------------------
*           Compute the orthogonal factorization PRQ = ( U ),  where P
*                                                      ( 0 )
*           is an orthogonal matrix and Q is a permutation matrix.
*           Overwrite R with the upper-triangle U.  The orthogonal
*           matrix P is applied to the residual and discarded.  The
*           permutation is stored in the array KX.  Once U has been
*           computed we need only work with vectors of length N within
*           LSCORE.  However, it is necessary to store the sum of
*           squares of the terms  B(NRANK+1),...,B(M),  where B = Pr.
*           ------------------------------------------------------------
            CALL DGEQRP( 'Column iterchanges', M, N, R, NROWR,
     $                   W(LWRK), IW(LKACTV), W(LGQ), INFO )
 
            LJ  = LKACTV
            DO 230 J = 1, N
               JMAX = IW(LJ)
               IF (JMAX .GT. J) THEN
                  JSAVE    = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J)    = JSAVE
               END IF
               LJ = LJ + 1
  230       CONTINUE
 
            CALL DGEAPQ( 'Transpose', 'Separate', M, N, R, NROWR,
     $                   W(LWRK), 1, B, M, W(LGQ), INFO )
 
            NRANK = IDRANK( MIN(N, M), R, NROWR+1, TOLRNK )
 
            IF (M .GT. NRANK) SSQ1 = DNRM2 ( M-NRANK, B(NRANK+1), 1 )
 
            IF (NRANK .GT. 0)
     $         CALL DCOPY ( NRANK, B, 1, W(LRES0), 1 )
         END IF
      ELSE
*        ===============================================================
*        R is input as an upper-triangular matrix with M rows.
*        ===============================================================
         NRANK = M
         IF (NRANK .GT. 0) THEN
            IF      (PRBTYP .EQ. 'QP') THEN
               CALL DLOAD ( NRANK, (ZERO), W(LRES0), 1 )
            ELSE IF (PRBTYP .EQ. 'LS') THEN
               CALL DCOPY ( NRANK, B, 1, W(LRES0), 1 )
            END IF
         END IF
      END IF
 
      IF (       MSGLVL .GT. 0     .AND.  NRANK  .LT. N
     $    .AND.  PRBTYP .NE. 'LP'  .AND.  PRBTYP .NE. 'FP')
     $   WRITE (NOUT, 9000) NRANK
 
*     ------------------------------------------------------------------
*     Find an initial working set.
*     ------------------------------------------------------------------
      CALL LSCRSH( COLD, VERTEX,
     $             NCLIN, NCTOTL, NACTIV, NARTIF,
     $             NFREE, N, NROWA,
     $             ISTATE, IW(LKACTV),
     $             BIGBND, TOLACT,
     $             A, W(LAX), BL, BU, X, W(LGQ), W(LWRK) )
 
*     ------------------------------------------------------------------
*     Compute the TQ factorization of the constraints while keeping R in
*     upper-triangular form.  Transformations associated with Q are
*     applied to CQ.  Transformations associated with P are applied to
*     RES0.  If some simple bounds are in the working set,  KX is
*     re-ordered so that the free variables come first.
*     ------------------------------------------------------------------
*     First, add the bounds. To save a bit of work, CQ is not loaded
*     until after KX has been re-ordered.
 
      NGQ   = 0
      NRES  = 0
      IF (NRANK .GT. 0) NRES = 1
      UNITQ = .TRUE.
 
      CALL LSBNDS( UNITQ,
     $             INFORM, NZ, NFREE, NRANK, NRES, NGQ,
     $             N, NQ, NROWA, NROWR, NROWT,
     $             ISTATE, KX,
     $             CONDMX,
     $             A, R, W(LT), W(LRES0), W(LCQ),
     $             W(LZY), W(LGQ), W(LWRK) )
 
      IF (LINOBJ) THEN
 
*        Install the transformed linear term in CQ.
*        CMQMUL applies the permutations in KX to CVEC.
 
         NGQ = 1
         CALL DCOPY ( N, CVEC, 1, W(LCQ), 1 )
         CALL CMQMUL( 6, N, NZ, NFREE, NQ, UNITQ,
     $                KX, W(LCQ), W(LZY), W(LWRK) )
      END IF
 
      IF (NACTIV .GT. 0) THEN
         NACT1  = NACTIV
         NACTIV = 0
 
         CALL LSADDS( UNITQ, VERTEX,
     $                INFORM, 1, NACT1, NACTIV, NARTIF, NZ, NFREE,
     $                NRANK, NREJTD, NRES, NGQ,
     $                N, NQ, NROWA, NROWR, NROWT,
     $                ISTATE, IW(LKACTV), KX,
     $                CONDMX,
     $                A, R, W(LT), W(LRES0), W(LCQ),
     $                W(LZY), W(LGQ), W(LWRK) )
      END IF
 
*     ------------------------------------------------------------------
*     Move the initial  x  onto the constraints in the working set.
*     Compute the transformed residual vector  Pr = Pb - RQ'x.
*     ------------------------------------------------------------------
      CALL LSSETX( LINOBJ, ROWERR, UNITQ,
     $             NCLIN, NACTIV, NFREE, NRANK, NZ,
     $             N, NCTOTL, NQ, NROWA, NROWR, NROWT,
     $             ISTATE, IW(LKACTV), KX,
     $             JMAX, ERRMAX, CTX, XNORM,
     $             A, W(LAX), BL, BU, W(LCQ), W(LRES), W(LRES0),
     $             W(LFEATL), R, W(LT), X, W(LZY), W(LPX), W(LWRK) )
 
      JINF = 0
 
      CALL LSCORE( PRBTYP, NAMED, NAMES, LINOBJ, UNITQ,
     $             INFORM, ITER, JINF, NCLIN, NCTOTL,
     $             NACTIV, NFREE, NRANK, NZ, NZ1,
     $             N, NROWA, NROWR,
     $             ISTATE, IW(LKACTV), KX,
     $             CTX, OBJ, SSQ1,
     $             SUMINF, NUMINF, XNORM,
     $             BL, BU, A, CLAMDA, W(LAX),
     $             W(LFEATL), R, X, IW, W )
 
      OBJ    = OBJ    + CTX
      IF (PRBTYP .EQ. 'LS'  .AND.  NRANK .GT. 0)
     $   CALL DCOPY ( NRANK, W(LRES), 1, B, 1 )
 
*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 IF (MSGLVL .GT.   0) THEN
         IF (INFORM .EQ.   0) THEN
            IF (PRBTYP .EQ. 'FP') THEN
               WRITE (NOUT, 2001)
            ELSE
               WRITE (NOUT, 2002) PRBTYP
            END IF
         END IF
         IF (INFORM .EQ.   1) WRITE (NOUT, 2010) PRBTYP
         IF (INFORM .EQ.   2) WRITE (NOUT, 2020) PRBTYP
         IF (INFORM .EQ.   3) WRITE (NOUT, 2030)
         IF (INFORM .EQ.   4) WRITE (NOUT, 2040)
         IF (INFORM .EQ.   5) WRITE (NOUT, 2050)
         IF (INFORM .EQ.   6) WRITE (NOUT, 2060) NERROR
 
         IF (INFORM .LT.   6) THEN
            IF      (NUMINF .EQ. 0) THEN
                IF (PRBTYP .NE. 'FP') WRITE (NOUT, 3000) PRBTYP, OBJ
            ELSE IF (INFORM .EQ. 3) THEN
               WRITE (NOUT, 3010) SUMINF
            ELSE
               WRITE (NOUT, 3020) SUMINF
            END IF
            IF (NUMINF .GT. 0) OBJ = SUMINF
         END IF
      END IF
 
*     Recover the optional parameters set by the user.
 
      CALL ICOPY ( MXPARM, IPSVLS, 1, IPRMLS, 1 )
      CALL DCOPY ( MXPARM, RPSVLS, 1, RPRMLS, 1 )
 
      RETURN
 
 2001 FORMAT(/ ' Exit LSSOL - Feasible point found.     ')
 2002 FORMAT(/ ' Exit LSSOL - Optimal ', A2, ' solution.')
 2010 FORMAT(/ ' Exit LSSOL - Weak ',    A2, ' solution.')
 2020 FORMAT(/ ' Exit LSSOL - ', A2,         ' solution is unbounded.' )
 2030 FORMAT(/ ' Exit LSSOL - Cannot satisfy the linear constraints. ' )
 2040 FORMAT(/ ' Exit LSSOL - Too many iterations.')
 2050 FORMAT(/ ' Exit LSSOL - Too many iterations without changing X.' )
 2060 FORMAT(/ ' Exit LSSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.'         )
 3000 FORMAT(/ ' Final ', A2, ' objective value =', G16.7 )
 3010 FORMAT(/ ' Minimum sum of infeasibilities =', G16.7 )
 3020 FORMAT(/ ' Final sum of infeasibilities =',   G16.7 )
 
 9000 FORMAT(/ ' Rank of the objective function data matrix = ', I5 )
 
*     End of  LSSOL .
 
      END
