*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     FILE LSMAIN FORTRAN
*
*     Sample program for LSSOL  Version 1.01 June 1986.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      IMPLICIT           DOUBLE PRECISION(A-H,O-Z)
 
*     Set the declared array dimensions.
*     NROWC  = the declared row dimension of  C.
*     NROWA  = the declared row dimension of  A.
*     MAXN   = maximum no. of variables allowed for.
*     MAXM   = maximum no. of observations allowed for.
*     MAXBND = maximum no. of variables + linear constraints.
*     LIWORK = the length of the integer work array.
*     LWORK  = the length of the double precision work array.
 
      PARAMETER         (NROWC  =   3, NROWA =  10,
     $                   MAXN   =   9, MAXM  =  10,
     $                   LIWORK =  60, LWORK = 900,
     $                   MAXBND = MAXN + NROWC )
 
      INTEGER            KX(MAXN), ISTATE(MAXBND)
      INTEGER            IWORK(LIWORK)
      DOUBLE PRECISION   C(NROWC,MAXN), B(MAXM)
      DOUBLE PRECISION   BL(MAXBND), BU(MAXBND), CLAMDA(MAXBND)
      DOUBLE PRECISION   CVEC(MAXN)
      DOUBLE PRECISION   A(NROWA,MAXN), X(MAXN)
      DOUBLE PRECISION   WORK(LWORK)
 
      DOUBLE PRECISION   BIGBND
      CHARACTER*10       CBGBND
 
      INTRINSIC          FLOAT
 
      PARAMETER        ( POINT1=0.1D+0, POINT3=0.3D+0, ONEPT5=1.5D+0 )
      PARAMETER        ( ZERO  =0.0D+0, ONE   =1.0D+0, TWO   =2.0D+0 )
      PARAMETER        ( THREE =3.0D+0, FOUR  =4.0D+0, FIVE  =5.0D+0 )
      PARAMETER        ( SIX   =6.0D+0                               )
 
      BIGBND =  1.0D+15
      CBGBND = '1.0D+15'
 
*     ==========================================
*     Example 1.  A linear least-squares problem.
*     ==========================================
*     Set the actual problem dimensions.
*     M      = the number of observations (rows of A)   (may be 0).
*     N      = the number of variables.
*     NCLIN  = the number of general linear constraints (may be 0).
 
      M      = 10
      N      = 9
      NCLIN  = 3
      NBND   = N + NCLIN
 
*     ------------------------------------------------------------------
*     Assign file numbers and problem data.
*     NOUT   = the unit number for printing.
*     IOPTNS = the unit number for reading the options file.
*     A      = the least-squares matrix.
*     B      = the vector of observations.
*     C      = the general constraint matrix.
*     BL     = the lower bounds on  x  and  C*x.
*     BU     = the upper bounds on  x  and  C*x.
*     X      = the initial estimate of the solution.
*     ------------------------------------------------------------------
      IOPTNS = 5
      NOUT   = 6
 
      DO 120 J = 1, N
         DO 110 I = 1, M
            A(I,J) = ONE
            B(I)   = ONE
  110    CONTINUE
  120 CONTINUE
 
      A(2 ,2) =   TWO
      A(10,2) =   ZERO
 
      A(3,3) =    THREE
      A(6,3) =    TWO
      A(9,3) =    ZERO
 
      A(4,4) =    FOUR
      A(5,4) =    THREE
      A(8,4) =    ZERO
 
      A(7,5) =    ZERO
 
      A(6,6) =    ZERO
 
      A(2 ,7) =   TWO
      A(3 ,7) = - ONE
      A(6 ,7) =   ZERO
      A(9 ,7) =   TWO
      A(10,7) =   ZERO
 
      A(2 ,8) =   ZERO
      A(3 ,8) = - ONE
      A(6 ,8) =   ZERO
      A(9 ,8) =   TWO
      A(10,8) =   TWO
 
      A(2 ,9) =   ZERO
      A(3 ,9) = - THREE
      A(6 ,9) = - ONE
      A(9 ,9) =   THREE
      A(10,9) =   TWO
 
      DO 140 J = 1, N
         DO 130 I = 1, NCLIN
            C(I,J) = ONE
  130    CONTINUE
  140 CONTINUE
 
      C(1,9) =   FOUR
 
      C(2,2) =   TWO
      C(2,3) =   THREE
      C(2,4) =   FOUR
      C(2,5) = - TWO
 
      C(3,2) = - ONE
      C(3,4) = - ONE
 
      DO 150 J = 1, N
         BL(J) = - TWO
         BU(J) =   TWO
  150 CONTINUE
      BL( 3) = - BIGBND
 
*     Set the ranges for the general constraints.
 
      BL(N+1) =   TWO
      BU(N+1) =   BIGBND
      BL(N+2) = - BIGBND
      BU(N+2) = - TWO
      BL(N+3) = - FOUR
      BU(N+3) = - TWO
 
      DO 170 J = 1, N
         X(J) = ONE / FLOAT(J)
  170 CONTINUE
 
 
*     ------------------------------------------------------------------
*     Read the options file.
*     Add a single option using a call to LSOPTN.
*     ------------------------------------------------------------------
 
      CALL LSFILE( IOPTNS, INFORM )
      IF (INFORM .NE. 0) THEN
         WRITE (NOUT, 3000) INFORM
         STOP
      END IF
 
      CALL LSOPTN( 'Infinite Bound size ='//CBGBND )
 
*     ------------------------------------------------------------------
*     Solve the problem.
*     ------------------------------------------------------------------
 
      CALL LSSOL ( M, N,
     $             NCLIN, NROWC, NROWA,
     $             C, BL, BU, CVEC,
     $             ISTATE, KX, X, A, B,
     $             INFORM, ITER, OBJ, CLAMDA,
     $             IWORK, LIWORK, WORK, LWORK )
 
*     Test for an error condition.
 
      IF (INFORM .GT. 1) GO TO 999
 
*     ================================================
*     Example 2.  A QP with Hessian bordered by zeros.
*     ================================================
*     Set the new problem dimensions.
*     M      = the number of rows (and columns) of A  (may be 0).
*     N      = the number of variables.
*     NCLIN  = the number of general linear constraints (may be 0).
*     CVEC   = the linear part of the objective function.
 
         M      = 5
         N      = 9
         NCLIN  = 3
         NBND   = N + NCLIN
 
         DO 220 J = 1, M
            DO 210 I = 1, J-1
               A(I,J) = ONE
  210       CONTINUE
  220    CONTINUE
 
         DO 230 I = 1, M
            A(I,I) = TWO
  230    CONTINUE
 
         DO 260 J = 1, N
            BL(J) = - TWO
            BU(J) =   TWO
  260    CONTINUE
 
         BL(N+1) = - TWO
         BU(N+1) =   ONEPT5
         BL(N+2) = - TWO
         BU(N+2) =   ONEPT5
         BL(N+3) = - TWO
         BU(N+3) =   FOUR
 
         DO 270 J = 1, N
            CVEC(J) = - ONE
  270    CONTINUE
         CVEC(1) = - FOUR
         CVEC(8) = - POINT1
         CVEC(9) = - POINT3
 
         DO 280 J = 1, N
            X(J) = ZERO
  280    CONTINUE
 
 
*     ------------------------------------------------------------------
*     Assign some new options.
*     ------------------------------------------------------------------
 
      CALL LSOPTN( 'Defaults         ' )
      CALL LSOPTN( 'Problem type  QP2' )
      CALL LSOPTN( 'Rank tolerance  =  1.0E-10' )
      CALL LSOPTN( 'Feasibility tolerance = 1.0E-10' )
 
*     ------------------------------------------------------------------
*     Solve the QP problem.
*     ------------------------------------------------------------------
 
      CALL LSSOL ( M, N,
     $             NCLIN, NROWC, NROWA,
     $             C, BL, BU, CVEC,
     $             ISTATE, KX, X, A, B,
     $             INFORM, ITER, OBJ, CLAMDA,
     $             IWORK, LIWORK, WORK, LWORK )
 
*     Test for an error condition.
 
      IF (INFORM .GT. 1) GO TO 999
      STOP
 
*     Error condition.
 
  999 WRITE (NOUT, 3010) INFORM
      STOP
 
 3000 FORMAT(/ '  LSFILE terminated with  INFORM =', I3)
 3010 FORMAT(/ '  LSSOL  terminated with  INFORM =', I3)
 
*     End of the example program for LSSOL.
 
      END
