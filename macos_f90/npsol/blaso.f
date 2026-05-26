c file blaso.f
c extracted from npsol distribution blas.f
c just the "other" and "qr" routines from blas.f
c
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  BLAS FORTRAN
*
*                         Level 1  BLAS
*                         -------  ----
*     DAXPY    DCOPY    DDOT     DNRM2    DSCAL    DSWAP    IDAMAX
*
*                         Level 2  BLAS
*                         -------  ----
*     DGEMV    DGER     DSYMV    DSYR     DTRMV    DTRSV    LSAME
*     XERBLA
*                         Others
*                         ------
*     DCOND*   DDDIV*   DDIV     DDSCL    DGRFG    DLOAD    DNORM
*     DROT3*   DROT3G*  DSSQ     ICOPY*   ILOAD    IDRANK+
*
*    *Not in the Nag Blas.
*    +Differs from the Nag Blas.
*
*                         QR Routines
*                         -- --------
*     DGEQR    DGEQRP   DGEAP    DGEAPQ
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DCOND ( N, X, INCX, AXMAX, AXMIN )
 
      INTEGER            N, INCX
      DOUBLE PRECISION   AXMAX, AXMIN
      DOUBLE PRECISION   X( (N-1)*INCX+1 )
C
C     DCOND   finds the elements in  x  that are largest and smallest
C     in magnitude.
C
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      INTEGER            I, IX
      INTRINSIC          ABS, MAX, MIN
 
      IF (N .EQ. 0) THEN
         AXMAX = ZERO
         AXMIN = ZERO
      ELSE
         AXMAX = ABS( X(1) )
         AXMIN = AXMAX
         IX    = 1
         DO 100 I = 2, N
            IX    = IX + INCX
            AXMAX = MAX( AXMAX, ABS( X(IX) ) )
            AXMIN = MIN( AXMIN, ABS( X(IX) ) )
  100    CONTINUE
      END IF
 
      RETURN
 
*     End of  DCOND
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION FUNCTION DDIV  ( A, B, FAIL )
      DOUBLE PRECISION                  A, B
      LOGICAL                           FAIL
C
C  DDIV   returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0 and a/b would overflow,
C
C  where flmax is a large value, via the function name. In addition if
C  a/b would overflow then fail is returned as true, otherwise fail is
C  returned as false.
C
C  Note that when a and b are both zero, fail is returned as true,
C  but div is returned as 0.0. in all other cases of overflow div is
C  such that abs( div ) = flmax.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 26-October-1982.
C     Sven Hammarling, Nag Central Office.
C
      INTRINSIC           ABS   , SIGN
      LOGICAL             FIRST
      DOUBLE PRECISION    ABSB  , FLMAX , FLMIN
      DOUBLE PRECISION    ONE   ,         ZERO
      PARAMETER         ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      SAVE                FIRST , FLMIN , FLMAX
      DATA                FIRST / .TRUE. /
 
      IF( A.EQ.ZERO )THEN
         DDIV   = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
         RETURN
      END IF
 
      IF( FIRST )THEN
         FIRST  = .FALSE.
         FLMIN  = WMACH( 5 )
         FLMAX  = WMACH( 7 )
      END IF
 
      IF( B.EQ.ZERO )THEN
         DDIV   = SIGN( FLMAX, A )
         FAIL   = .TRUE.
      ELSE
         ABSB   = ABS( B )
         IF( ABSB.GE.ONE )THEN
            FAIL = .FALSE.
            IF( ABS( A ).GE.ABSB*FLMIN )THEN
               DDIV   = A/B
            ELSE
               DDIV   = ZERO
            END IF
         ELSE
            IF( ABS( A ).LE.ABSB*FLMAX )THEN
               FAIL   = .FALSE.
               DDIV   = A/B
            ELSE
               FAIL   = .TRUE.
               DDIV   = FLMAX
               IF( ( ( A.LT.ZERO ).AND.( B.GT.ZERO ) ).OR.
     $             ( ( A.GT.ZERO ).AND.( B.LT.ZERO ) )     )
     $            DDIV   = -DDIV
            END IF
         END IF
      END IF
 
      RETURN
 
*     End of DDIV  .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DDDIV ( N, D, INCD, X, INCX )
      INTEGER            N, INCD, INCX
      DOUBLE PRECISION   D( * ), X( * )
C
C     DDDIV  performs the operation
C
C     x := diag( d )(inverse)*x
C
      PARAMETER        ( ONE = 1.0 )
      EXTERNAL           DSCAL
      INTEGER            I     , ID    , IX
 
      IF( N.GE.1 )THEN
         IF( INCD.EQ.0 )THEN
 
            CALL DSCAL ( N, (ONE/D( 1 )), X, INCX )
 
         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = X( ID )/D( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = X( IX )/D( ID )
                  IX      = IX + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = X( IX )/D( ID )
                  ID      = ID + INCD
                  IX      = IX + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
 
      RETURN
 
*     End of DDDIV .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DDSCL ( N, D, INCD, X, INCX )
      INTEGER            N, INCD, INCX
      DOUBLE PRECISION   D( * ), X( * )
C
C  DDSCL  performs the operation
C
C     x := diag( d )*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DSCAL
      INTEGER            I     , ID    , IX
 
      IF( N.GE.1 )THEN
         IF( INCD.EQ.0 )THEN
 
            CALL DSCAL ( N, D( 1 ), X, INCX )
 
         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = D( ID )*X( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = D( ID )*X( IX )
                  IX      = IX + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = D( ID )*X( IX )
                  ID      = ID + INCD
                  IX      = IX + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
 
      RETURN
 
*     End of DDSCL .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DGRFG ( N, ALPHA, X, INCX, TOL, ZETA )
      INTEGER            N, INCX
      DOUBLE PRECISION   ALPHA, X( * ), TOL, ZETA
C
C  DGRFG  generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
      EXTERNAL           DSSQ  , DSCAL
      INTRINSIC          ABS   , MAX   , SIGN  , SQRT
      LOGICAL            FIRST
      DOUBLE PRECISION   BETA  , EPS   , SCALE , SSQ
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE
 
         EPS    =  WMACH( 3 )
 
*        Treat case where P is a 2 by 2 matrix specially.
 
         IF( N.EQ.1 )THEN
 
*           Deal with cases where  ALPHA = zero  and
*           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
 
            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  ABS( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( ABS( X( 1 ) ).LE.MAX( EPS*ABS( ALPHA ),
     $                                     TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( ABS( ALPHA ).GE.ABS( X( 1 ) ) )THEN
                  BETA = ABS ( ALPHA  )*
     $                   SQRT( ONE + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = ABS ( X( 1 ) )*
     $                   SQRT( ONE + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA   =  SQRT( ( ABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  =  BETA
            END IF
         ELSE
 
*           Now P is larger than 2 by 2.
 
            SSQ   = ONE
            SCALE = ZERO
 
            CALL DSSQ  ( N, X, INCX, SCALE, SSQ )
 
*           Treat cases where  SCALE = zero,
*           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
*           ALPHA = zero  specially.
*           Note that  SCALE = max( abs( X( i ) ) ).
 
            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.MAX( EPS*ABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*SQRT( SSQ )
 
               CALL DSCAL ( N, -ONE/ALPHA, X, INCX )
 
            ELSE
               IF( SCALE.LT.ABS( ALPHA ) )THEN
                  BETA = ABS ( ALPHA )*
     $                   SQRT( ONE + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE*
     $                   SQRT( SSQ +     ( ALPHA/SCALE )**2 )
               END IF
               ZETA = SQRT( ( BETA + ABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )BETA = -BETA
 
               CALL DSCAL( N, -ONE/( ZETA*BETA ), X, INCX )
 
               ALPHA = BETA
            END IF
         END IF
      END IF
      RETURN
 
*     End of DGRFG .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DLOAD ( N, CONST, X, INCX )
      INTEGER            N, INCX
      DOUBLE PRECISION   CONST
      DOUBLE PRECISION   X( * )
C
C  DLOAD  performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
      INTEGER            IX
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
 
      IF( N.LT.1 )RETURN
 
      IF( CONST.NE.ZERO )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            X( IX ) = CONST
   10    CONTINUE
      ELSE
         DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
            X( IX ) = ZERO
   20    CONTINUE
      END IF
 
      RETURN
 
*     End of DLOAD .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION FUNCTION DNORM ( SCALE, SSQ )
      DOUBLE PRECISION                  SCALE, SSQ
C
C  DNORM  returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
      INTRINSIC           SQRT
      LOGICAL             FIRST
      DOUBLE PRECISION    FLMAX , SQT
      DOUBLE PRECISION    ONE
      PARAMETER         ( ONE   = 1.0D+0 )
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      SAVE                FIRST , FLMAX
      DATA                FIRST / .TRUE. /
 
      IF( FIRST )THEN
         FIRST = .FALSE.
         FLMAX = WMACH( 7 )
      END IF
 
      SQT = SQRT( SSQ )
      IF( SCALE.LT.FLMAX/SQT )THEN
         DNORM  = SCALE*SQT
      ELSE
         DNORM  = FLMAX
      END IF
 
      RETURN
 
*     End of DNORM .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DROT3 ( N, X, INCX, Y, INCY, CS, SN )
 
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   CS, SN
      DOUBLE PRECISION   X(*), Y(*)
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
C
C  DROT3   applies the plane rotation defined by CS and SN to the
C  columns of a 2 by N matrix held in X and Y.  The method used requires
C  3 multiplications and 3 additions per column, as described in Gill,
C  Golub, Murray and Saunders, Mathematics of Computation 28 (1974) 505-
C  -535 (see page 508).
C
C  DROT3   guards against underflow, and overflow is extremely unlikely.
C  It is assumed that CS and SN have been generated by DROT3G, ensuring
C  that CS lies in the closed interval (0, 1),  and that the absolute
C  value of CS and SN (if nonzero) is no less than the machine precision
C  EPS.  It is also assumed that  RTMIN .lt. EPS.  Note that the magic
C  number Z is therefore no less than 0.5*EPS in absolute value, so it
C  is safe to use TOL = 2*RTMIN in the underflow test involving Z*A.
C  For efficiency we use the same TOL in the previous two tests.
C
C  Systems Optimization Laboratory, Stanford University.
C  Original version dated January 1982.
C  F77 version dated 28-June-1986.
C  This version of DROT3 dated 28-June-1986.
C
      INTEGER            I, IX, IY
      DOUBLE PRECISION   A, B, ONE, RTMIN, TOL, W, Z, ZERO
      INTRINSIC          ABS
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      IF (N .LT. 1  .OR.  SN .EQ. ZERO) RETURN
      IX = 1
      IY = 1
      IF (CS .EQ. ZERO) THEN
 
*        Just swap  x  and  y.
 
         DO 10 I = 1, N
            A     = X(IX)
            X(IX) = Y(IY)
            Y(IY) = A
            IX    = IX + INCX
            IY    = IY + INCY
   10    CONTINUE
 
      ELSE
 
         RTMIN  = WMACH(6)
         TOL    = RTMIN + RTMIN
         Z      = SN/(ONE + CS)
 
         DO 20 I = 1, N
            A     = X(IX)
            B     = Y(IY)
            W     = ZERO
            IF (ABS(A) .GT. TOL) W = CS*A
            IF (ABS(B) .GT. TOL) W = W + SN*B
            X(IX) = W
            A     = A + W
            IF (ABS(A) .GT. TOL) B = B - Z*A
            Y(IY) = - B
            IX    =   IX + INCX
            IY    =   IY + INCY
   20    CONTINUE
 
      END IF
 
      RETURN
 
*     End of  DROT3
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DROT3G( X, Y, CS, SN )
 
      DOUBLE PRECISION   X, Y, CS, SN
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
C
C  DROT3G  generates a plane rotation that reduces the vector (X, Y) to
C  the vector (A, 0),  where A is defined as follows...
C
C     If both X and Y are negligibly small, or
C     if Y is negligible relative to Y,
C     then  A = X,  and the identity rotation is returned.
C
C     If X is negligible relative to Y,
C     then  A = Y,  and the swap rotation is returned.
C
C     Otherwise,  A = sign(X) * sqrt( X**2 + Y**2 ).
C
C  In all cases,  X and Y are overwritten by A and 0,  and CS will lie
C  in the closed interval (0, 1).  Also,  the absolute value of CS and
C  SN (if nonzero) will be no less than the machine precision,  EPS.
C
C  DROT3G  guards against overflow and underflow.
C  It is assumed that  FLMIN .lt. EPS**2  (i.e.  RTMIN .lt. EPS).
C
C  Systems Optimization Laboratory, Stanford University.
C  Original version dated January 1982.
C  F77 version dated 28-June-1986.
C  This version of DROT3G dated 28-June-1986.
C
      DOUBLE PRECISION   A, B, EPS, ONE, RTMIN, ZERO
      LOGICAL            FIRST
      INTRINSIC          ABS, MAX, SQRT
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
 
      SAVE               FIRST , EPS   , RTMIN
      DATA               FIRST / .TRUE. /
 
      IF( FIRST )THEN
         FIRST = .FALSE.
         EPS    = WMACH(3)
         RTMIN  = WMACH(6)
      END IF
 
      IF (Y .EQ. ZERO) THEN
 
         CS = ONE
         SN = ZERO
 
      ELSE IF (X .EQ. ZERO) THEN
 
         CS = ZERO
         SN = ONE
         X  = Y
 
      ELSE
 
         A      = ABS(X)
         B      = ABS(Y)
         IF (MAX(A,B) .LE. RTMIN) THEN
            CS = ONE
            SN = ZERO
         ELSE
            IF (A .GE. B) THEN
               IF (B .LE. EPS*A) THEN
                  CS = ONE
                  SN = ZERO
                  GO TO 900
               ELSE
                  A  = A * SQRT( ONE + (B/A)**2 )
               END IF
            ELSE
               IF (A .LE. EPS*B) THEN
                  CS = ZERO
                  SN = ONE
                  X  = Y
                  GO TO 900
               ELSE
                  A  = B * SQRT( ONE + (A/B)**2 )
               END IF
            END IF
            IF (X .LT. ZERO) A = - A
            CS = X/A
            SN = Y/A
            X  = A
         END IF
      END IF
 
  900 Y  = ZERO
 
      RETURN
 
*     End of  DROT3G
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DSSQ  ( N, X, INCX, SCALE, SUMSQ )
      INTEGER            N, INCX
      DOUBLE PRECISION   X( * )
      DOUBLE PRECISION   SCALE, SUMSQ
C
C  DSSQ   returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = y( 1 )**2 +...+ y( n )**2 + ( scale**2 )*sumsq,
C
C  where y( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
      INTRINSIC          ABS
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
      IF( N.GE.1 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = ONE   + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
 
*     End of DSSQ  .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE ICOPY ( N, IX, INCIX, IY, INCIY )
 
      INTEGER            N, INCIX, INCIY
      INTEGER            IX(*), IY(*)
 
C
C  Copy the first N elements of IX into IY.
C
 
      INTEGER            J, JX, JY
 
      IF (N .GE. 1) THEN
         IF (INCIX .EQ. 1  .AND.  INCIY .EQ. 1) THEN
 
            DO 10 J = 1, N
               IY(J) = IX(J)
   10       CONTINUE
 
         ELSE
 
            JX = 1
            JY = 1
            DO 20 J = 1, N
               IY(JY) = IX(JX)
               JX = JX + INCIX
               JY = JY + INCIY
   20       CONTINUE
 
         END IF
      END IF
 
      RETURN
 
*     End of  ICOPY
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE ILOAD ( N, ICONST, IX, INCIX )
      INTEGER            N, INCIX
      INTEGER            ICONST
      INTEGER            IX( * )
C
C  ILOAD   performs the operation
C
C     ix = iconst*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
      INTEGER            JX
 
      IF( N.LT.1 )RETURN
 
      IF( ICONST.NE.0 )THEN
         DO 10, JX = 1, 1 + ( N - 1 )*INCIX, INCIX
            IX( JX ) = ICONST
   10    CONTINUE
      ELSE
         DO 20, JX = 1, 1 + ( N - 1 )*INCIX, INCIX
            IX( JX ) = 0
   20    CONTINUE
      END IF
 
      RETURN
 
*     End of ILOAD .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      INTEGER           FUNCTION IDRANK( N, X, INCX, TOL )
      INTEGER                            N, INCX
      DOUBLE PRECISION                   X( * ), TOL
 
C  IDRANK finds the first element of the n element vector x for which
C
C     abs( x( k ) ).le.( tol*max ( abs(x(1)), ..., abs(x(k-1)) )
C
C  and returns the value ( k - 1 ) in the function name IDRANK. If no
C  such k exists then IDRANK is returned as n.
C
C  If TOL is supplied as less than zero then the value EPSMCH, where
C  EPSMCH is the relative machine precision, is used in place of TOL.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-January-1985.
C     Sven Hammarling, Nag Central Office.
C     Modified by PEG, 19-December-1985.
 
      INTRINSIC                          ABS   , MAX
      INTEGER                            IX    , K
      DOUBLE PRECISION                   TOLRNK, XMAX  , ZERO
      PARAMETER                        ( ZERO  = 0.0D+0 )
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      K = 0
      IF (N .GE. 1) THEN
         TOLRNK = TOL
         IF (TOL .LT. ZERO) TOLRNK = WMACH(3)
 
         IF( INCX .GT. 0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
 
         XMAX = ABS( X(IX) )
 
*+       WHILE (K .LT. N) LOOP
   10    IF    (K .LT. N) THEN
            IF (ABS( X(IX) ) .LE. XMAX*TOLRNK) GO TO 20
            XMAX = MAX( XMAX, ABS( X(IX) ) )
            K    = K  + 1
            IX   = IX + INCX
            GO TO 10
         END IF
*+       END WHILE
 
      END IF
   20 IDRANK = K
      RETURN
 
*     End of IDRANK.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DGEQR ( M, N, A, LDA, ZETA, INFORM )
      INTEGER            M, N, LDA, INFORM
      DOUBLE PRECISION   A( LDA, * ), ZETA( * )
C
C  1. Purpose
C     =======
C
C  DGEQR  reduces the  m by n, m.ge.n, matrix A to upper triangular form
C  by means of orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m  orthogonal matrix and  R  is an n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used to introduce zeros  into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of z( k ) are in a( k + 1, k ), ..., a( m, k ). The elements of R are
C  returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( n, m - 1 ).
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the  upper  triangular  matrix  R  and the  M by N  strictly
C           lower triangular part of  A  will  contain  details  of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least ( n ).
C
C           On  exit, ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then   ZETA( k ) = 0.0
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and is always in the range ( 1.0, sqrt( 2.0 ) ).
C
C  INFORM - INTEGER.
C
C           On successful  exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear  algebra routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Operations involving the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRMV  and DTRSV . Note that no test for near
C  singularity of R is incorporated in this routine or in routine  DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular  then  the  auxiliary linear algebra  routine  DUTSV
C  can  be used to  determine  the  singular value decomposition  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DGEMV , DGER  , DGRFG
      INTRINSIC          MIN
      INTEGER            J     , K     , LA
      DOUBLE PRECISION   TEMP
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
*     Check the input parameters.
 
      IF( N.EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( M.LT.N ).OR.( N.LT.0 ).OR.( LDA.LT.M ) )THEN
         INFORM = 1
         RETURN
      END IF
 
*     Perform the factorization.
 
      LA = LDA
      DO 20, K = 1, MIN( M - 1, N )
 
*        Use a Householder reflection to zero the kth column of A.
*        First set up the reflection.
 
         CALL DGRFG ( M - K, A( K, K ), A( K + 1, K ), 1, ZERO,
     $                ZETA( K ) )
         IF( ( ZETA( K ).GT.ZERO ).AND.( K.LT.N ) )THEN
            IF( ( K + 1 ).EQ.N )
     $         LA = M - K + 1
            TEMP      = A( K, K )
            A( K, K ) = ZETA( K )
 
*           We now perform the operation  A := Q( k )*A.
 
*           Let B denote the bottom ( m - k + 1 ) by ( n - k ) part
*           of A.
 
*           First form  work = B'*u. ( work is stored in the elements
*           ZETA( k + 1 ), ..., ZETA( n ). )
 
            CALL DGEMV ( 'Transpose', M - K + 1, N - K,
     $                   ONE, A( K, K + 1 ), LA, A( K, K ), 1,
     $                   ZERO, ZETA( K + 1 ), 1 )
 
*           Now form  B := B - u*work'.
 
            CALL DGER  ( M - K + 1, N - K, -ONE, A( K, K ), 1,
     $                   ZETA( K + 1 ), 1, A( K, K + 1 ), LA )
 
*           Restore beta.
 
            A( K, K ) = TEMP
         END IF
   20 CONTINUE
 
*     Store the final zeta when m.eq.n.
 
      IF( M.EQ.N )
     $   ZETA( N ) = ZERO
 
      INFORM = 0
      RETURN
 
*     End of DGEQR .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DGEQRP( PIVOT, M, N, A, LDA, ZETA, PERM, WORK, INFORM )
      CHARACTER*1        PIVOT
      INTEGER            M, N, LDA, INFORM
      INTEGER            PERM( * )
      DOUBLE PRECISION   A( LDA, * ), ZETA( * ), WORK( * )
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
C  1. Purpose
C     =======
C
C  DGEQRP reduces the  m by n matrix A to upper triangular form by means
C  of orthogonal transformations and column permutations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )*P'      when   m.gt.n,
C           ( 0 )
C
C     A = Q*R*P'          when   m = n,
C
C     A = Q*( R  X )*P'   when   m.lt.n,
C
C  where  Q  is  an  m by m  orthogonal matrix, R  is a  min( m, n )  by
C  min( m, n )  upper triangular matrix and  P is an  n by n permutation
C  matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ),  which is used to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of z( k ) are in a( k + 1, k ), ..., a( m, k ). The elements of R are
C  returned in the upper triangular part of A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( m - 1, n ).
C
C  Two options are available for the column permutations. In either case
C  the column for which the  sub-diagonal elements are to be annihilated
C  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
C  The  particular column chosen as the pivot column is either that  for
C  which  the  unreduced  part  ( elements k onwards )  has the  largest
C  Euclidean  length, or  is that for  which the ratio of the  Euclidean
C  length  of the  unreduced part  to the  Euclidean length of the whole
C  column is a maximum.
C
C  3. Parameters
C     ==========
C
C  PIVOT  - CHARACTER*1.
C
C           On  entry, PIVOT  specifies  the  pivoting  strategy  to  be
C           performed as follows.
C
C           PIVOT = 'C' or 'c'
C
C              Column  interchanges  are  to be  incorporated  into  the
C              factorization, such that the  column whose unreduced part
C              has  maximum  Euclidean  length  is chosen  as the  pivot
C              column at each step.
C
C           PIVOT = 'S' or 's'
C
C              Scaled  column interchanges  are to be  incorporated into
C              the  factorization, such  that the  column for which  the
C              ratio  of the  Euclidean  length of the unreduced part of
C              the column to the original Euclidean length of the column
C              is a maximum is chosen as the  pivot column at each step.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of A will contain the upper triangular matrix R and the
C           M by min( M, N )  strictly lower triangular part of  A  will
C           contain details  of the  factorization  as  described above.
C           When m.lt.n then the remaining M by ( N - M ) part of A will
C           contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must  specify  the leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least ( n ).
C
C           On exit, ZETA( k )  contains the scalar  zeta  for  the  kth
C           transformation. If T( k ) = I then ZETA( k) = 0.0, otherwise
C           ZETA( k )  contains the scalar  zeta( k ) as described above
C           and  is  always  in  the  range  ( 1.0, sqrt( 2.0 ) ).  When
C           n .gt. m  the  elements  ZETA( m + 1 ),  ZETA( m + 2 ), ...,
C           ZETA( n )  are used as internal workspace.
C
C  PERM   - INTEGER array of DIMENSION at least min( m, n ).
C
C           On exit, PERM  contains details of the permutation matrix P,
C           such  that  PERM( k ) = k  if no  column interchange occured
C           at  the  kth  step  and  PERM( k ) = j, ( k .lt. j .le. n ),
C           if  columns  k and j  were  interchanged at  the  kth  step.
C           Note  that, although  there are  min( m - 1, n )  orthogonal
C           transformations, there are min( m, n ) permutations.
C
C  WORK   - 'real' array of DIMENSION at least ( 2*n ).
C
C           Used as internal workspace.
C
C           On exit, WORK( j ), j = 1, 2, ..., n, contains the Euclidean
C           length  of the  jth  column  of the  permuted  matrix  A*P'.
C
C  INFORM - INTEGER.
C
C           On  successful exit, INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly supplied. See the next section for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        PIVOT .ne. 'C' or 'c' or 'S' or 's'
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear algebra  routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Also following the use of this routine the operations
C
C     B := P'*B   and   B := P*B,
C
C  where B is an n by k matrix, and the operations
C
C     B := B*P    and   B := B*P',
C
C  where  B is a k by n  matrix, can  be performed by calls to the basic
C  linear  algebra  routine  DGEAP .  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL DGEAP ( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL DGEAP ( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  N  in the above two calls.
C  The operation  B := B*P  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine DGEAP for further details.
C
C  Operations involving  the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV.  Note that no test for near
C  singularity of  R is incorporated in this routine or in routine DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular then  the  auxiliary  linear algebra  routine  DUTSV
C  can  be  used  to  determine  the  singular value decomposition of R.
C  Operations  involving  the  matrix  X  can also be  performed  by the
C  Level 2  BLAS  routines.  Matrices  of  the  form   ( R  X )  can  be
C  factorized as
C
C     ( R  X ) = ( T  0 )*S',
C
C  where  T is upper triangular and S is orthogonal, using the auxiliary
C  linear algebra routine  DUTRQ .
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           MCHPAR, DGEMV , DGER  , DGRFG , DNRM2 , DSWAP
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      INTEGER            J     , JMAX  , K     , LA
      DOUBLE PRECISION   EPS   , MAXNRM, NORM  , DNRM2 , TEMP  , TOL
      DOUBLE PRECISION   LAMDA
      PARAMETER        ( LAMDA = 1.0D-2 )
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
*     Check the input parameters.
 
      IF( MIN( M, N ).EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( ( PIVOT.NE.'C' ).AND.( PIVOT.NE.'c' ).AND.
     $      ( PIVOT.NE.'S' ).AND.( PIVOT.NE.'s' )      ).OR.
     $    ( M.LT.0 ).OR.( N.LT.0 ).OR.( LDA.LT.M )           )THEN
         INFORM = 1
         RETURN
      END IF
 
*     Compute eps and the initial column norms.
 
      CALL MCHPAR()
      EPS = WMACH( 3 )
      DO 10, J = 1, N
         WORK( J )     = DNRM2 ( M, A( 1, J ), 1 )
         WORK( J + N ) = WORK( J )
   10 CONTINUE
 
*     Perform the factorization. TOL is the tolerance for DGRFG .
 
      LA = LDA
      DO 50, K = 1, MIN( M, N )
 
*        Find the pivot column.
 
         MAXNRM = ZERO
         JMAX   = K
         DO 20, J = K, N
            IF( ( PIVOT.EQ.'C' ).OR.( PIVOT.EQ.'c' ) )THEN
               IF( WORK( J + N  ).GT.MAXNRM )THEN
                  MAXNRM = WORK( J + N )
                  JMAX   = J
               END IF
            ELSE IF( WORK( J ).GT.ZERO )THEN
               IF( ( WORK( J + N )/WORK( J ) ).GT.MAXNRM )THEN
                  MAXNRM = WORK( J + N )/WORK( J )
                  JMAX   = J
               END IF
            END IF
   20    CONTINUE
         PERM( K ) = JMAX
         IF( JMAX.GT.K )THEN
            CALL DSWAP ( M, A( 1, K ), 1, A( 1, JMAX ), 1 )
            TEMP             = WORK( K )
            WORK( K )        = WORK( JMAX )
            WORK( JMAX )     = TEMP
            WORK( JMAX + N ) = WORK( K + N )
            PERM( K )        = JMAX
         END IF
         TOL = EPS*WORK( K )
         IF( K.LT.M )THEN
 
*           Use a Householder reflection to zero the kth column of A.
*           First set up the reflection.
 
            CALL DGRFG ( M - K, A( K, K ), A( K + 1, K ), 1, TOL,
     $                   ZETA( K ) )
            IF( K.LT.N )THEN
               IF( ZETA( K ).GT.ZERO )THEN
                  IF( ( K + 1 ).EQ.N )
     $               LA = M - K + 1
                  TEMP      = A( K, K )
                  A( K, K ) = ZETA( K )
 
*                 We now perform the operation  A := Q( k )*A.
 
*                 Let B denote the bottom ( m - k + 1 ) by ( n - k )
*                 part of A.
 
*                 First form  work = B'*u. ( work is stored in the
*                 elements ZETA( k + 1 ), ..., ZETA( n ). )
 
                  CALL DGEMV ( 'Transpose', M - K + 1, N - K,
     $                         ONE, A( K, K + 1 ), LA, A( K, K ), 1,
     $                         ZERO, ZETA( K + 1 ), 1 )
 
*                 Now form  B := B - u*work'.
 
                  CALL DGER  ( M - K + 1, N - K, -ONE, A( K, K ), 1,
     $                         ZETA( K + 1 ), 1, A( K, K + 1 ), LA )
 
*                 Restore beta.
 
                  A( K, K ) = TEMP
               END IF
 
*              Update the unreduced column norms. Use the Linpack
*              criterion for when to recompute the norms, except that
*              we retain the original column lengths throughout and use
*              a smaller lamda.
 
               DO 40, J = K + 1, N
                  IF( WORK( J + N ).GT.ZERO )THEN
                     TEMP = ABS( A( K, J ) )/WORK( J + N )
                     TEMP = MAX( ( ONE + TEMP )*( ONE - TEMP ), ZERO )
                     NORM = TEMP
                     TEMP = ONE +
     $                      LAMDA*TEMP*( WORK( J + N )/WORK( J ) )**2
                     IF( TEMP.GT.ONE )THEN
                        WORK( J + N ) = WORK( J + N )*SQRT( NORM )
                     ELSE
                        WORK( J + N ) = DNRM2 ( M - K,
     $                                          A( K + 1, J ), 1 )
                     END IF
                  END IF
   40          CONTINUE
            END IF
         END IF
   50 CONTINUE
 
*     Store the final zeta when m.le.n.
 
      IF( M.LE.N )
     $   ZETA( M ) = ZERO
 
      INFORM = 0
      RETURN
 
*     End of DGEQRP.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DGEAP ( SIDE, TRANS, M, N, PERM, K, B, LDB )
*     .. Scalar Arguments ..
      INTEGER            K, LDB, M, N
      CHARACTER*1        SIDE, TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * )
      INTEGER            PERM( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEAP  performs one of the transformations
*
*     B := P'*B   or   B := P*B,   where B is an m by k matrix,
*
*  or
*
*     B := B*P'   or   B := B*P,   where B is a k by m matrix,
*
*  P being an m by m permutation matrix of the form
*
*     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
*
*  where  P( i, index( i ) ) is the permutation matrix that interchanges
*  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
*  with rows and columns  i and index( i )  interchanged.  Of course, if
*  index( i ) = i  then  P( i, index( i ) ) = I.
*
*  This routine  is intended for use in  conjunction with  Nag auxiliary
*  routines that  perform  interchange  operations,  such  as  pivoting.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*  TRANS
*           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
*           TRANS  ( Transpose, or No transpose )  specify the operation
*           to be performed as follows.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := P'*B.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := P*B.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := B*P'.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := B*P.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*
*           On entry, M must specify the order of the permutation matrix
*           P.  M must be at least zero.  When  M = 0  then an immediate
*           return is effected.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*
*           On entry,  N must specify the value of n. N must be at least
*           zero.  When  N = 0  then an  immediate  return is  effected.
*
*           Unchanged on exit.
*
*  PERM   - INTEGER array of DIMENSION at least ( n ).
*
*           Before  entry,  PERM  must  contain the  n  indices  for the
*           permutation matrices. index( i ) must satisfy
*
*              1 .le. index( i ) .le. m.
*
*           It is usual for index( i ) to be at least i, but this is not
*           necessary for this routine.
*
*           Unchanged on exit.
*
*  K      - INTEGER.
*
*           On entry with  SIDE = 'L' or 'l',  K must specify the number
*           of columns of B and on entry with  SIDE = 'R' or 'r', K must
*           specify the  number of rows of B.  K must be at least  zero.
*           When  K = 0  then an immediate return is effected.
*
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of  DIMENSION  ( LDB, ncolb ),  where
*           ncolb = k   when   SIDE = 'L' or 'l'  and   ncolb = m   when
*           SIDE = 'R' or 'r'.
*
*           Before entry  with  SIDE = 'L' or 'l',  the  leading  M by K
*           part  of  the  array   B  must  contain  the  matrix  to  be
*           transformed  and  before entry with  SIDE = 'R' or 'r',  the
*           leading  K by M part of the array  B must contain the matrix
*           to  be  transformed.  On  exit,  B  is  overwritten  by  the
*           transformed matrix.
*
*  LDB    - INTEGER.
*
*           On entry,  LDB  must specify  the  leading dimension  of the
*           array  B  as declared  in the  calling  (sub) program.  When
*           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
*           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
*           Unchanged on exit.
*
*
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 13-January-1986.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, J, L
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
      IF( MIN( M, N, K ).EQ.0 )
     $   RETURN
      LEFT  = ( SIDE .EQ.'L' ).OR.( SIDE .EQ.'l' )
      RIGHT = ( SIDE .EQ.'R' ).OR.( SIDE .EQ.'r' )
      NULL  = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20, I = 1, N
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 10, J = 1, K
                     TEMP      = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40, I = N, 1, -1
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 30, J = 1, K
                     TEMP      = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60, J = 1, N
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 50, I = 1, K
                     TEMP      = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80, J = N, 1, -1
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 70, I = 1, K
                     TEMP      = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEAP . ( F06QJF )
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DGEAPQ( TRANS, WHEREZ, M, N, A, LDA, ZETA,
     $                   NCOLB, B, LDB, WORK, INFORM )
      CHARACTER*1        TRANS, WHEREZ
      INTEGER            M, N, LDA, NCOLB, LDB, INFORM
      DOUBLE PRECISION   A( LDA, * ), ZETA( * ), B( LDB, * ), WORK( * )
C
C  1. Purpose
C     =======
C
C  DGEAPQ performs one of the transformations
C
C     B := Q'*B   or   B := Q*B,
C
C  where B is an m by ncolb matrix and Q is an m by m orthogonal matrix,
C  given as the product of  Householder transformation matrices, details
C  of  which are stored in the  m by n ( m.ge.n )  array  A  and, if the
C  parameter  WHEREZ = 'S' or 's', in the array ZETA.
C
C  This  routine is  intended for use following auxiliary linear algebra
C  routines such as  DGEQR , DGEHES and DSLTRI. ( See those routines for
C  example calls. )
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k )  and  zeta( k ) must be supplied either
C  in  a( k, k )  or in  zeta( k ), depending upon the parameter WHEREZ.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining Q'.
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = ' ' or 'N' or 'n'
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'T' or 't' or 'C' or 'c'
C
C              Perform the operation  B := Q'*B.
C
C           Unchanged on exit.
C
C  WHEREZ - CHARACTER*1.
C
C           On entry, WHEREZ specifies where the elements of zeta are to
C           be found as follows.
C
C           WHEREZ = 'I' or 'i'
C
C              The elements of zeta are in A.
C
C           WHEREZ = 'S' or 's'
C
C              The elements of zeta are separate from A, in ZETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHEREZ = 'I' or 'i'  then  the  diagonal
C           elements of A must contain the elements of zeta.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least min( m - 1, n ).
C
C           Before entry with  WHEREZ = 'S' or 's', the array  ZETA must
C           contain the elements of the vector  zeta.
C
C           When  WHEREZ = 'I' or 'i', the array ZETA is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - 'real' array of DIMENSION ( LDB, ncolb ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - 'real' array of DIMENSION at least ( ncolb ).
C
C           Used as internal workspace.
C
C  INFORM - INTEGER.
C
C           On  successful exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        TRANS  .ne. ' ' or 'N' or 'n' or 'T' or 't' or 'C' or 'c'
C        WHEREZ .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DGEMV , DGER
      INTRINSIC          MIN
      INTEGER            J     , K     , KK    , LB
      DOUBLE PRECISION   TEMP
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
*     Check the input parameters.
 
      IF( MIN( N, NCOLB ).EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( ( TRANS .NE.' ' ).AND.
     $      ( TRANS .NE.'N' ).AND.( TRANS .NE.'n' ).AND.
     $      ( TRANS .NE.'T' ).AND.( TRANS .NE.'t' ).AND.
     $      ( TRANS .NE.'C' ).AND.( TRANS .NE.'c' )      ).OR.
     $    ( ( WHEREZ.NE.'I' ).AND.( WHEREZ.NE.'i' ).AND.
     $      ( WHEREZ.NE.'S' ).AND.( WHEREZ.NE.'s' )      ).OR.
     $    ( M.LT.N ).OR.( N.LT.0 ).OR.( LDA.LT.M ).OR.
     $    ( NCOLB.LT.0 ).OR.( LDB.LT.M )                      )THEN
         INFORM = 1
         RETURN
      END IF
 
*     Perform the transformation.
 
      LB = LDB
      DO 20, KK = 1, MIN( M - 1, N )
         IF( ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
     $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )THEN
 
*           Q'*B = Q( p )*...*Q( 2 )*Q( 1 )*B,     p = min( m - 1, n ).
 
            K = KK
         ELSE
 
*           Q*B  = Q( 1 )'*Q( 2 )'*...*Q( p )'*B,  p = min( m - 1, n ).
*           Note that  Q( k )' = Q( k ).
 
            K = MIN( N, M - 1 ) + 1 - KK
         END IF
         IF( ( WHEREZ.EQ.'S' ).OR.( WHEREZ.EQ.'s' ) )THEN
            TEMP      = A( K, K )
            A( K, K ) = ZETA( K )
         END IF
 
*        If ZETA( k ) is zero then Q( k ) = I and we can skip the kth
*        transformation.
 
         IF( A( K, K ).GT.ZERO )THEN
            IF( NCOLB.EQ.1 )
     $         LB = M - K + 1
 
*           Let C denote the bottom ( m - k + 1 ) by ncolb part of B.
 
*           First form  work = C'*u.
 
            DO 10, J = 1, NCOLB
               WORK( J ) = ZERO
   10       CONTINUE
            CALL DGEMV ( 'Transpose', M - K + 1, NCOLB,
     $                   ONE, B( K, 1 ), LB, A( K, K ), 1,
     $                   ZERO, WORK, 1 )
 
*           Now form  C := C - u*work'.
 
            CALL DGER  ( M - K + 1, NCOLB, -ONE, A( K, K ), 1,
     $                   WORK, 1, B( K, 1 ), LB )
         END IF
 
*        Restore the diagonal element of A.
 
         IF( ( WHEREZ.EQ.'S' ).OR.( WHEREZ.EQ.'s' ) )
     $      A( K, K ) = TEMP
   20 CONTINUE
 
      INFORM = 0
      RETURN
 
*     End of DGEAPQ.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c end of blaso.f

