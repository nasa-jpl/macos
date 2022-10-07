*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     FILE NPMAIN FORTRAN
*
*     Sample program for NPSOL Version 4.02  August 1986.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
*     Set the declared array dimensions.
*     NROWA  = the declared row dimension of  A.
*     NROWJ  = the declared row dimension of  CJAC.
*     NROWR  = the declared row dimension of  R.
*     MAXN   = maximum no. of variables allowed for.
*     MAXBND = maximum no. of variables + linear & nonlinear constrnts.
*     LIWORK = the length of the integer work array.
*     LWORK  = the length of the double precision work array.
 
      PARAMETER         (NROWA  =  5, NROWJ  =  20, NROWR =   10,
     $                   MAXN   =  9, LIWORK =  70, LWORK = 1000,
     $                   MAXBND =  MAXN + NROWA + NROWJ)
 
      INTEGER            ISTATE(MAXBND)
      INTEGER            IWORK(LIWORK)
      DOUBLE PRECISION   A(NROWA,MAXN)
      DOUBLE PRECISION   BL(MAXBND), BU(MAXBND)
      DOUBLE PRECISION   C(NROWJ), CJAC(NROWJ,MAXN), CLAMDA(MAXBND)
      DOUBLE PRECISION   OBJGRD(MAXN), R(NROWR,MAXN), X(MAXN)
      DOUBLE PRECISION   WORK(LWORK)
      EXTERNAL           OBJFN1, OBJFN2, CONFN1, CONFN2
 
      PARAMETER         (ZERO = 0.0, ONE = 1.0)
 
*     Set the actual problem dimensions.
*     N      = the number of variables.
*     NCLIN  = the number of general linear constraints (may be 0).
*     NCNLN  = the number of nonlinear constraints (may be 0).
 
      N      = 9
      NCLIN  = 4
      NCNLN  = 14
      NBND   = N + NCLIN + NCNLN
 
*     ------------------------------------------------------------------
*     Assign file numbers and the data arrays.
*     NOUT   = the unit number for printing.
*     IOPTNS = the unit number for reading the options file.
*     Bounds  .ge.    BIGBND  will be treated as plus  infinity.
*     Bounds  .le.  - BIGBND  will be treated as minus infinity.
*     A      = the linear constraint matrix.
*     BL     = the lower bounds on  x,  a'x  and  c(x).
*     BU     = the upper bounds on  x,  a'x  and  c(x).
*     X      = the initial estimate of the solution.
*     ------------------------------------------------------------------
      NOUT   = 6
      IOPTNS = 5
      BIGBND = 1.0D+15
 
*     Set the matrix  A.
 
      DO 40 J = 1, N
         DO 30 I = 1, NCLIN
            A(I,J) = ZERO
   30    CONTINUE
   40 CONTINUE
      A(1,1) = -ONE
      A(1,2) =  ONE
      A(2,2) = -ONE
      A(2,3) =  ONE
      A(3,3) =  ONE
      A(3,4) = -ONE
      A(4,4) =  ONE
      A(4,5) = -ONE
 
*     Set the bounds.
 
      DO 50 J  =  1, NBND
         BL(J) = -BIGBND
         BU(J) =  BIGBND
   50 CONTINUE
      BL(1)  =  ZERO
      BL(3)  = -ONE
      BL(5)  =  ZERO
      BL(6)  =  ZERO
      BL(7)  =  ZERO
 
      BU(3)  =  ONE
      BU(8)  =  ZERO
      BU(9)  =  ZERO
 
*     Set lower bounds of zero for all four linear constraints.
 
      DO 60 J  =  N+1,  N+NCLIN
         BL(J) =  ZERO
   60 CONTINUE
 
*     Set upper bounds of one for all 14 nonlinear constraints.
 
      DO 70 J  =  N + NCLIN + 1,  NBND
         BU(J) =  ONE
   70 CONTINUE
 
*     Set the initial estimate of  X.
 
      X(1)   =  .1
      X(2)   =  .125
      X(3)   =  .666666
      X(4)   =  .142857
      X(5)   =  .111111
      X(6)   =  .2
      X(7)   =  .25
      X(8)   = -.2
      X(9)   = -.25
 
 
*     ------------------------------------------------------------------
*     Read the options file.
*     ------------------------------------------------------------------
 
      CALL NPFILE( IOPTNS, INFORM )
      IF (INFORM .NE. 0) THEN
         WRITE (NOUT, 3000) INFORM
         STOP
      END IF
 
*     ------------------------------------------------------------------
*     Solve the problem.
*     ------------------------------------------------------------------
 
      CALL NPSOL ( N, NCLIN, NCNLN, NROWA, NROWJ, NROWR,
     $             A, BL, BU,
     $             CONFN1, OBJFN1,
     $             INFORM, ITER, ISTATE,
     $             C, CJAC, CLAMDA, OBJF, OBJGRD, R, X,
     $             IWORK, LIWORK, WORK, LWORK )
 
	if (.true.) stop
      IF (INFORM .GT. 0) GO TO 900
 
*     ------------------------------------------------------------------
*     The following is for illustrative purposes only.
*     A second run solves the same problem,  but defines the objective
*     and constraints via the subroutines OBJFN2 and CONFN2.  Some
*     objective derivatives and the constant Jacobian elements are not
*     supplied.
*     We do a warm start using
*              ISTATE    (the working set)
*              CLAMDA    (the Lagrange multipliers)
*              R         (the Hessian approximation)
*     from the previous run, but with a slightly perturbed starting
*     point.  The previous option file must have specified
*              Hessian    Yes
*     for R to be a useful approximation.
*     ------------------------------------------------------------------
 
      DO 100 J = 1, N
         X(J)  = X(J) + 0.1
  100 CONTINUE
 
*     The previous parameters are retained and updated.
 
      CALL NPOPTN( '   Derivative level               0')
      CALL NPOPTN( '   Verify                        No')
      CALL NPOPTN( '   Warm Start')
      CALL NPOPTN( '   Major iterations              20')
 
      CALL NPOPTN( '   Major print level             10')
 
      CALL NPSOL ( N, NCLIN, NCNLN, NROWA, NROWJ, NROWR,
     $             A, BL, BU,
     $             CONFN2, OBJFN2,
     $             INFORM, ITER, ISTATE,
     $             C, CJAC, CLAMDA, OBJF, OBJGRD, R, X,
     $             IWORK, LIWORK, WORK, LWORK )
 
      IF (INFORM .GT. 0) GO TO 900
      STOP
 
*     -----------
*     Error exit.
*     -----------
 
  900 WRITE (NOUT, 3010) INFORM
      STOP
 
 3000 FORMAT(/ ' NPFILE terminated with  INFORM =', I3)
 3010 FORMAT(/ ' NPSOL  terminated with  INFORM =', I3)
 
*     End of the example program for NPSOL.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE OBJFN1( MODE, N, X, OBJF, OBJGRD, NSTATE )
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION   X(N), OBJGRD(N)
 
*-----------------------------------------------------------------------
*     OBJFN1  computes the value and first derivatives of the nonlinear
*     objective function.
*-----------------------------------------------------------------------
      OBJF   = - X(2)*X(6) + X(1)*X(7) - X(3)*X(7) - X(5)*X(8)
     $         + X(4)*X(9) + X(3)*X(8)
 
      OBJGRD(1) =   X(7)
      OBJGRD(2) = - X(6)
      OBJGRD(3) = - X(7) + X(8)
      OBJGRD(4) =   X(9)
      OBJGRD(5) = - X(8)
      OBJGRD(6) = - X(2)
      OBJGRD(7) = - X(3) + X(1)
      OBJGRD(8) = - X(5) + X(3)
      OBJGRD(9) =   X(4)
 
      RETURN
 
*     End of  OBJFN1.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CONFN1( MODE, NCNLN, N, NROWJ,
     $                   NEEDC, X, C, CJAC, NSTATE )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   X(N), C(*), CJAC(NROWJ,*)
 
*-----------------------------------------------------------------------
*     CONFN1  computes the values and first derivatives of the nonlinear
*     constraints.
*
*     The zero elements of Jacobian matrix are set only once.  This
*     occurs during the first call to CONFN1  (NSTATE = 1).
*-----------------------------------------------------------------------
      PARAMETER         (ZERO = 0.0, TWO = 2.0)
 
      IF (NSTATE .EQ. 1) THEN
 
*        First call to CONFN1.  Set all Jacobian elements to zero.
*        N.B.  This will only work with `Derivative Level = 3'.
 
         DO 120 J = 1, N
            DO 110 I = 1, NCNLN
               CJAC(I,J) = ZERO
  110       CONTINUE
  120    CONTINUE
 
      END IF
 
      IF (NEEDC(1) .GT. 0) THEN
         C(1)       =   X(1)**2  +  X(6)**2
         CJAC(1,1)  =   TWO*X(1)
         CJAC(1,6)  =   TWO*X(6)
      END IF
 
      IF (NEEDC(2) .GT. 0) THEN
         C(2)       =   (X(2) - X(1))**2  +  (X(7) - X(6))**2
         CJAC(2,1)  = - TWO*(X(2) - X(1))
         CJAC(2,2)  =   TWO*(X(2) - X(1))
         CJAC(2,6)  = - TWO*(X(7) - X(6))
         CJAC(2,7)  =   TWO*(X(7) - X(6))
      END IF
 
      IF (NEEDC(3) .GT. 0) THEN
         C(3)       =   (X(3) - X(1))**2  +  X(6)**2
         CJAC(3,1)  = - TWO*(X(3) - X(1))
         CJAC(3,3)  =   TWO*(X(3) - X(1))
         CJAC(3,6)  =   TWO*X(6)
      END IF
 
      IF (NEEDC(4) .GT. 0) THEN
         C(4)       =   (X(1) - X(4))**2  +  (X(6) - X(8))**2
         CJAC(4,1)  =   TWO*(X(1) - X(4))
         CJAC(4,4)  = - TWO*(X(1) - X(4))
         CJAC(4,6)  =   TWO*(X(6) - X(8))
         CJAC(4,8)  = - TWO*(X(6) - X(8))
      END IF
 
      IF (NEEDC(5) .GT. 0) THEN
         C(5)       =   (X(1) - X(5))**2  +  (X(6) - X(9))**2
         CJAC(5,1)  =   TWO*(X(1) - X(5))
         CJAC(5,5)  = - TWO*(X(1) - X(5))
         CJAC(5,6)  =   TWO*(X(6) - X(9))
         CJAC(5,9)  = - TWO*(X(6) - X(9))
      END IF
 
      IF (NEEDC(6) .GT. 0) THEN
         C(6)       =   X(2)**2  +  X(7)**2
         CJAC(6,2)  =   TWO*X(2)
         CJAC(6,7)  =   TWO*X(7)
      END IF
 
      IF (NEEDC(7) .GT. 0) THEN
         C(7)       =   (X(3) - X(2))**2  +  X(7)**2
         CJAC(7,2)  = - TWO*(X(3) - X(2))
         CJAC(7,3)  =   TWO*(X(3) - X(2))
         CJAC(7,7)  =   TWO*X(7)
      END IF
 
      IF (NEEDC(8) .GT. 0) THEN
         C(8)       =   (X(4) - X(2))**2  +  (X(8) - X(7))**2
         CJAC(8,2)  = - TWO*(X(4) - X(2))
         CJAC(8,4)  =   TWO*(X(4) - X(2))
         CJAC(8,7)  = - TWO*(X(8) - X(7))
         CJAC(8,8)  =   TWO*(X(8) - X(7))
      END IF
 
      IF (NEEDC(9) .GT. 0) THEN
         C(9)       =   (X(2) - X(5))**2  +  (X(7) - X(9))**2
         CJAC(9,2)  =   TWO*(X(2) - X(5))
         CJAC(9,5)  = - TWO*(X(2) - X(5))
         CJAC(9,7)  =   TWO*(X(7) - X(9))
         CJAC(9,9)  = - TWO*(X(7) - X(9))
      END IF
 
      IF (NEEDC(10) .GT. 0) THEN
         C(10)      =   (X(4) - X(3))**2  +  X(8)**2
         CJAC(10,3) = - TWO*(X(4) - X(3))
         CJAC(10,4) =   TWO*(X(4) - X(3))
         CJAC(10,8) =   TWO*X(8)
      END IF
 
      IF (NEEDC(11) .GT. 0) THEN
         C(11)      =   (X(5) - X(3))**2  +  X(9)**2
         CJAC(11,3) = - TWO*(X(5) - X(3))
         CJAC(11,5) =   TWO*(X(5) - X(3))
         CJAC(11,9) =   TWO*X(9)
      END IF
 
      IF (NEEDC(12) .GT. 0) THEN
         C(12)      =   X(4)**2  +  X(8)**2
         CJAC(12,4) =   TWO*X(4)
         CJAC(12,8) =   TWO*X(8)
      END IF
 
      IF (NEEDC(13) .GT. 0) THEN
         C(13)      =   (X(4) - X(5))**2  +  (X(9) - X(8))**2
         CJAC(13,4) =   TWO*(X(4) - X(5))
         CJAC(13,5) = - TWO*(X(4) - X(5))
         CJAC(13,8) = - TWO*(X(9) - X(8))
         CJAC(13,9) =   TWO*(X(9) - X(8))
      END IF
 
      IF (NEEDC(14) .GT. 0) THEN
         C(14)      =   X(5)**2  +  X(9)**2
         CJAC(14,5) =   TWO*X(5)
         CJAC(14,9) =   TWO*X(9)
      END IF
 
      RETURN
 
*     End of  CONFN1.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE OBJFN2( MODE, N, X, OBJF, OBJGRD, NSTATE )
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION   X(N), OBJGRD(N)
 
*-----------------------------------------------------------------------
*     OBJFN2  computes the value and some first derivatives of the
*     nonlinear objective function.
*-----------------------------------------------------------------------
 
      OBJF   = - X(2)*X(6) + X(1)*X(7) - X(3)*X(7) - X(5)*X(8)
     $         + X(4)*X(9) + X(3)*X(8)
 
      OBJGRD(3) = - X(7) + X(8)
      OBJGRD(7) = - X(3) + X(1)
      OBJGRD(8) = - X(5) + X(3)
 
      RETURN
 
*     End of  OBJFN2.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE CONFN2( MODE, NCNLN, N, NROWJ,
     $                   NEEDC, X, C, CJAC, NSTATE )
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
      INTEGER            NEEDC(*)
      DOUBLE PRECISION   X(N), C(*), CJAC(NROWJ,*)
 
*-----------------------------------------------------------------------
*     CONFN2  computes the values and the non-constant derivatives of
*     the nonlinear constraints.
*-----------------------------------------------------------------------
      PARAMETER         (TWO = 2.0)
 
      IF (NEEDC(1) .GT. 0) THEN
         C(1)       =   X(1)**2  +  X(6)**2
         CJAC(1,1)  =   TWO*X(1)
         CJAC(1,6)  =   TWO*X(6)
      END IF
 
      IF (NEEDC(2) .GT. 0) THEN
         C(2)       =   (X(2) - X(1))**2  +  (X(7) - X(6))**2
         CJAC(2,1)  = - TWO*(X(2) - X(1))
         CJAC(2,2)  =   TWO*(X(2) - X(1))
         CJAC(2,6)  = - TWO*(X(7) - X(6))
         CJAC(2,7)  =   TWO*(X(7) - X(6))
      END IF
 
      IF (NEEDC(3) .GT. 0) THEN
         C(3)       =   (X(3) - X(1))**2  +  X(6)**2
         CJAC(3,1)  = - TWO*(X(3) - X(1))
         CJAC(3,3)  =   TWO*(X(3) - X(1))
         CJAC(3,6)  =   TWO*X(6)
      END IF
 
      IF (NEEDC(4) .GT. 0) THEN
         C(4)       =   (X(1) - X(4))**2  +  (X(6) - X(8))**2
         CJAC(4,1)  =   TWO*(X(1) - X(4))
         CJAC(4,4)  = - TWO*(X(1) - X(4))
         CJAC(4,6)  =   TWO*(X(6) - X(8))
         CJAC(4,8)  = - TWO*(X(6) - X(8))
      END IF
 
      IF (NEEDC(5) .GT. 0) THEN
         C(5)       =   (X(1) - X(5))**2  +  (X(6) - X(9))**2
         CJAC(5,1)  =   TWO*(X(1) - X(5))
         CJAC(5,5)  = - TWO*(X(1) - X(5))
         CJAC(5,6)  =   TWO*(X(6) - X(9))
         CJAC(5,9)  = - TWO*(X(6) - X(9))
      END IF
 
      IF (NEEDC(6) .GT. 0) THEN
         C(6)       =   X(2)**2  +  X(7)**2
         CJAC(6,2)  =   TWO*X(2)
         CJAC(6,7)  =   TWO*X(7)
      END IF
 
      IF (NEEDC(7) .GT. 0) THEN
         C(7)       =   (X(3) - X(2))**2  +  X(7)**2
         CJAC(7,2)  = - TWO*(X(3) - X(2))
         CJAC(7,3)  =   TWO*(X(3) - X(2))
         CJAC(7,7)  =   TWO*X(7)
      END IF
 
      IF (NEEDC(8) .GT. 0) THEN
         C(8)       =   (X(4) - X(2))**2  +  (X(8) - X(7))**2
         CJAC(8,2)  = - TWO*(X(4) - X(2))
         CJAC(8,4)  =   TWO*(X(4) - X(2))
         CJAC(8,7)  = - TWO*(X(8) - X(7))
         CJAC(8,8)  =   TWO*(X(8) - X(7))
      END IF
 
      IF (NEEDC(9) .GT. 0) THEN
         C(9)       =   (X(2) - X(5))**2  +  (X(7) - X(9))**2
         CJAC(9,2)  =   TWO*(X(2) - X(5))
         CJAC(9,5)  = - TWO*(X(2) - X(5))
         CJAC(9,7)  =   TWO*(X(7) - X(9))
         CJAC(9,9)  = - TWO*(X(7) - X(9))
      END IF
 
      IF (NEEDC(10) .GT. 0) THEN
         C(10)      =   (X(4) - X(3))**2  +  X(8)**2
         CJAC(10,3) = - TWO*(X(4) - X(3))
         CJAC(10,4) =   TWO*(X(4) - X(3))
         CJAC(10,8) =   TWO*X(8)
      END IF
 
      IF (NEEDC(11) .GT. 0) THEN
         C(11)      =   (X(5) - X(3))**2  +  X(9)**2
         CJAC(11,3) = - TWO*(X(5) - X(3))
         CJAC(11,5) =   TWO*(X(5) - X(3))
         CJAC(11,9) =   TWO*X(9)
      END IF
 
      IF (NEEDC(12) .GT. 0) THEN
         C(12)      =   X(4)**2  +  X(8)**2
         CJAC(12,4) =   TWO*X(4)
         CJAC(12,8) =   TWO*X(8)
      END IF
 
      IF (NEEDC(13) .GT. 0) THEN
         C(13)      =   (X(4) - X(5))**2  +  (X(9) - X(8))**2
         CJAC(13,4) =   TWO*(X(4) - X(5))
         CJAC(13,5) = - TWO*(X(4) - X(5))
         CJAC(13,8) = - TWO*(X(9) - X(8))
         CJAC(13,9) =   TWO*(X(9) - X(8))
      END IF
 
      IF (NEEDC(14) .GT. 0) THEN
         C(14)      =   X(5)**2  +  X(9)**2
         CJAC(14,5) =   TWO*X(5)
         CJAC(14,9) =   TWO*X(9)
      END IF
 
      RETURN
 
*     End of  CONFN2.
 
      END
