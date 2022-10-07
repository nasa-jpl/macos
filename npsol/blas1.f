c file blas1.f
c extracted from npsol distribution blas.f
c just the level 1 blas
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
 
      SUBROUTINE DAXPY ( N, ALPHA, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * ), Y( * )
 
C  DAXPY  performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.
 
      INTEGER            I     , IX    , IY
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO  = 0.0D+0 )
 
      IF( N    .LT.1    )RETURN
      IF( ALPHA.EQ.ZERO )RETURN
 
      IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10    CONTINUE
      ELSE
         IF( INCY.GE.0 )THEN
            IY = 1
         ELSE
            IY = 1 - ( N - 1 )*INCY
         END IF
         IF( INCX.GT.0 )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               Y( IY ) = ALPHA*X( IX ) + Y( IY )
               IY      = IY + INCY
   20       CONTINUE
         ELSE
            IX = 1 - ( N - 1 )*INCX
            DO 30, I = 1, N
               Y( IY ) = ALPHA*X( IX ) + Y( IY )
               IX      = IX + INCX
               IY      = IY + INCY
   30       CONTINUE
         END IF
      END IF
      RETURN
 
*     End of DAXPY .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DCOPY ( N, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   X( * ), Y( * )
 
C  DCOPY  performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine DCOPY .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
 
      INTEGER            I     , IX    , IY
 
      IF( N.LT.1 )RETURN
 
      IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
         DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
            Y( IY ) = X( IY )
   10    CONTINUE
      ELSE
         IF( INCX.GE.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IX )
               IX      = IX + INCX
   20       CONTINUE
         ELSE
            IY = 1 - ( N - 1 )*INCY
            DO 30, I = 1, N
               Y( IY ) = X( IX )
               IY      = IY + INCY
               IX      = IX + INCX
   30       CONTINUE
         END IF
      END IF
 
      RETURN
 
*     End of DCOPY .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION FUNCTION DDOT  ( N, X, INCX, Y, INCY )
      INTEGER                           N, INCX, INCY
      DOUBLE PRECISION                  X( * ), Y( * )
 
C  DDOT   returns the value
C
C     DDOT   = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT  .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
 
      INTEGER             I     , IX    , IY
      DOUBLE PRECISION    SUM
      DOUBLE PRECISION    ZERO
      PARAMETER         ( ZERO  = 0.0D+0 )
 
      SUM = ZERO
      IF( N.GE.1 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
 
      DDOT   = SUM
      RETURN
 
*     End of DDOT  .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
      INTEGER                           N, INCX
      DOUBLE PRECISION                  X( * )
 
C  DNRM2  returns the Euclidean norm of a vector via the function
C  name, so that
C
C     DNRM2  := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2 .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
 
      EXTERNAL            DNORM , DSSQ
      INTRINSIC           ABS
      DOUBLE PRECISION    SCALE , DNORM , SSQ
      DOUBLE PRECISION    ONE   ,         ZERO
      PARAMETER         ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
      IF( N.LT.1 )THEN
         DNRM2  = ZERO
      ELSE IF( N.EQ.1 )THEN
         DNRM2  = ABS( X( 1 ) )
      ELSE
         SCALE  = ZERO
         SSQ    = ONE
 
         CALL DSSQ  ( N, X, INCX, SCALE, SSQ )
 
         DNRM2  = DNORM ( SCALE, SSQ )
 
      END IF
      RETURN
 
*     End of DNRM2 .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DSCAL ( N, ALPHA, X, INCX )
      INTEGER            N, INCX
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * )
 
C  DSCAL  performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
 
      INTEGER            IX
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
 
      IF( N.GE.1 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
 
      RETURN
 
*     End of DSCAL .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE DSWAP ( N, X, INCX, Y, INCY )
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   X( * ), Y( * )
 
C  DSWAP  performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP .
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
 
      INTEGER            I     , IX    , IY
      DOUBLE PRECISION   TEMP
 
      IF( N.LT.1 )RETURN
 
      IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
         DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
            TEMP    = X( IY )
            X( IY ) = Y( IY )
            Y( IY ) = TEMP
   10    CONTINUE
      ELSE
         IF( INCX.GE.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IX )
               X( IX ) = Y( IY )
               Y( IY ) = TEMP
               IX      = IX + INCX
   20       CONTINUE
         ELSE
            IY = 1 - ( N - 1 )*INCY
            DO 30, I = 1, N
               TEMP    = X( IX )
               X( IX ) = Y( IY )
               Y( IY ) = TEMP
               IY      = IY + INCY
               IX      = IX + INCX
   30       CONTINUE
         END IF
      END IF
 
      RETURN
 
*     End of DSWAP .
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      INTEGER FUNCTION IDAMAX( N, X, INCX )
      INTEGER                  N, INCX
      DOUBLE PRECISION         X( * )
 
C  IDAMAX returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
 
      INTRINSIC           ABS
      INTEGER             I     , IMAX  , IX
      DOUBLE PRECISION    XMAX
 
      IF( N.LT.1 )THEN
         IDAMAX = 0
         RETURN
      END IF
 
      IMAX = 1
      IF( N.GT.1 )THEN
         XMAX = ABS( X( 1 ) )
         IX   = 1
         DO 10, I = 2, N
            IX = IX + INCX
            IF( XMAX.LT.ABS( X( IX ) ) )THEN
               XMAX = ABS( X( IX ) )
               IMAX = I
            END IF
   10    CONTINUE
      END IF
 
      IDAMAX = IMAX
      RETURN
 
*     End of IDAMAX.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c end of blas1.f

