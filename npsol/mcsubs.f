*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  MCSUBS FORTRAN
*
*     MCHPAR   MCEPS    MCENV1   MCENV2   MCSTOR   MCSMAL   MCMIN
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE MCHPAR()
 
************************************************************************
*  MCHPAR  must define certain machine parameters as follows:
*     wmach(1)  = NBASE  = base of floating-point arithmetic.
*     wmach(2)  = NDIGIT = no. of base wmach(1) digits of precision.
*     wmach(3)  = EPS    = floating-point precision.
*     wmach(4)  = RTEPS  = sqrt(EPS).
*     wmach(5)  = RMIN   = smallest positive normalized floating-point
*                          number.
*     wmach(6)  = RTRMIN = sqrt(RMIN).
*     wmach(7)  = RMAX   = largest positive floating-point number.
*     wmach(8)  = RTRMAX = sqrt(RMAX).
*     wmach(9)  = UNDFLW = 0 if underflow is not fatal, +ve otherwise.
*     wmach(10) = NIN    = standard file number of the input stream.
*     wmach(11) = NOUT   = standard file number of the output stream.
************************************************************************
 
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
 
      EXTERNAL           MCENV2, MCEPS , MCSMAL
      INTRINSIC          SQRT
      LOGICAL            FIRST , HDWIRE
      INTEGER            EMIN  , NBASE , NDIGIT, NIN   , NOUT
      DOUBLE PRECISION   BASE  , EPS   , MCEPS , MCSMAL, RMAX  , RMIN
      DOUBLE PRECISION   RTEPS , RTMAX , RTMIN , SMALL , UNDFLW
      SAVE               FIRST
      DATA               FIRST / .TRUE. /
 
      IF (FIRST) THEN
         FIRST = .FALSE.
 
*        ---------------------------------------------------------------
*        Machine-dependent code.
*        1. Set UNDFLW, NIN, NOUT, HDWIRE as desired.
*        2. If  HDWIRE = .TRUE.  set the machine constants
*               NBASE, NDIGIT, EPS, RMIN, RMAX
*           in-line.  Otherwise, they will be computed by MCENV2.
*           A call to MCENV2 will cause eight underflows.
*        ---------------------------------------------------------------
 
         UNDFLW = 0
         NIN    = 5
         NOUT   = 6
         HDWIRE = .FALSE.
 
         IF (HDWIRE) THEN
            NBASE  = 16
            NDIGIT = 14
            BASE   = NBASE
            EPS    = BASE**(1 - NDIGIT)
            RMIN   = BASE**(-62)
            RMAX   = BASE**(+62)
         ELSE
            CALL MCENV2( NBASE, NDIGIT, EPS, EMIN, RMIN )
 
            EPS    = MCEPS ()
            SMALL  = MCSMAL()
            RMAX   = 1/SMALL
         END IF
 
         WMACH(1)  = NBASE
         WMACH(2)  = NDIGIT
         WMACH(3)  = EPS
         WMACH(4)  = SQRT( EPS )
         WMACH(5)  = RMIN
         WMACH(6)  = SQRT( RMIN )
         WMACH(7)  = RMAX
         WMACH(8)  = SQRT( RMAX )
         WMACH(9)  = UNDFLW
         WMACH(10) = NIN
         WMACH(11) = NOUT
      END IF
 
      RETURN
 
*     End of  MCHPAR.
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION FUNCTION MCEPS ()
 
*     MCEPS  returns approximately the relative machine precision via
*     the function name.
*
*     The value returned is given by
*
*        MCEPS  = (1/2)*beta**( 1 - t )   when   rnd = true
*
*        MCEPS  =       beta**( 1 - t )   when   rnd = false,
*
*     where beta is the base of the machine, t is the number of ( beta )
*     digits in the mantissa and rnd is true when rounding occurs and is
*     false when chopping occurs. This is the Wilkinson unit rounding
*     error.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine (EPSLON).
*
*  -- Written on 26-November-1984.
*     Sven Hammarling, Nag Central Office.
 
      EXTERNAL           MCENV1
      LOGICAL            FIRST , RND
      INTEGER            BETA  , T
      DOUBLE PRECISION   BASE  , EPS
 
      SAVE               EPS   , FIRST
      DATA               FIRST / .TRUE. /
 
      IF( FIRST )THEN
         FIRST = .FALSE.
 
         CALL MCENV1( BETA, T, RND )
 
         BASE = BETA
         IF( RND )THEN
            EPS = ( BASE**( 1 - T ) )/2
         ELSE
            EPS =   BASE**( 1 - T )
         END IF
      END IF
 
      MCEPS  = EPS
      RETURN
 
*     End of MCEPS  (EPSLON).
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE MCENV1( BETA, T, RND )
      LOGICAL            RND
      INTEGER            BETA, T
 
*     MCENV1 returns the machine parameters given by:
*
*        BETA - INTEGER.
*               The base of the machine.
*
*        T    - INTEGER.
*               The number of ( BETA ) digits in the mantissa.
*
*        RND  - LOGICAL.
*               Whether proper rounding ( RND = .TRUE. ) or chopping
*               ( RND = .FALSE. ) occurs in addition. This may not be a
*               reliable guide to the way in which the machine perfoms
*               its arithmetic.
*
*     The routine is based on the routine ENVRON by Malcolm
*     and incorporates suggestions by Gentleman and Marovich. See
*
*        Malcolm M. A. (1972) Algorithms to reveal properties of
*           floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*        Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*           that reveal properties of floating point arithmetic units.
*           Comms. of the ACM, 17, 276-277.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine (ENVRON).
*
*  -- Written on 26-November-1984.
*     Sven Hammarling and Mick Pont, Nag Central Office.
 
      EXTERNAL           MCSTOR
      LOGICAL            FIRST , LRND
      INTEGER            LBETA , LT
      DOUBLE PRECISION   A     , B     , C     , F     , ONE   , QTR
      DOUBLE PRECISION   MCSTOR
 
      SAVE               FIRST , LBETA , LRND  , LT
      DATA               FIRST / .TRUE. /
 
      IF( FIRST )THEN
         FIRST = .FALSE.
         ONE   = 1
 
*        LBETA, LT and LRND are the local values of BETA, T and RND.
*
*        Throughout this routine we use the function MCSTOR to ensure
*        that relevant values are stored and not held in registers, or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
 
         A = 1
         C = 1
 
*+       WHILE( C.EQ.ONE )LOOP
   10    IF   ( C.EQ.ONE )THEN
            A = 2*A
            C = MCSTOR( A,  ONE )
            C = MCSTOR( C, -A   )
            GO TO 10
         END IF
*+       END WHILE
 
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
 
         B = 1
         C = MCSTOR( A, B )
 
*+       WHILE( C.EQ.A )LOOP
   20    IF   ( C.EQ.A )THEN
            B = 2*B
            C = MCSTOR( A, B )
            GO TO 20
         END IF
*+       END WHILE
 
*        Now compute the base. a and b are neighbouring floating point
*        numbers in the interval ( beta**t, beta**( t + 1 ) ) and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
 
 
         QTR   = ONE/4
         C     = MCSTOR( C, -A )
         LBETA = C + QTR
 
*        Now determine whether rounding or chopping occurs, by adding
*        a bit less than beta/2 and a bit more than beta/2 to a.
 
         B = LBETA
         F = MCSTOR( B/2, -B/100 )
         C = MCSTOR( F, A )
         IF( C.EQ.A) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = MCSTOR( B/2,  B/100 )
         C = MCSTOR( F, A )
         IF( ( LRND ).AND.( C.EQ.A ) )
     $      LRND = .FALSE.
 
*        Now find the mantissa, t. It should be the integer part of
*        log to the base beta of a, however it is safer to determine t
*        by powering. So we find t as the smallest positive integer
*        for which
*
*           fl( beta**t + 1.0 ) = 1.0.
 
         LT = 0
         A  = 1
         C  = 1
 
*+       WHILE( C.EQ.ONE )LOOP
   30    IF   ( C.EQ.ONE )THEN
            LT = LT + 1
            A  = A*LBETA
            C  = MCSTOR( A,  ONE )
            C  = MCSTOR( C, -A   )
            GO TO 30
         END IF
*+       END WHILE
 
      END IF
 
      BETA = LBETA
      T    = LT
      RND  = LRND
      RETURN
 
*     End of MCENV1 (ENVRON).
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE MCENV2( BETA, T, EPS, EMIN, RMIN )
      INTEGER            BETA, T, EMIN
      DOUBLE PRECISION   EPS, RMIN
 
*     MCENV2 returns the machine parameters given by:
*
*        BETA - INTEGER.
*               The base of the machine.
*
*        T    - INTEGER.
*               The number of ( BETA ) digits in the mantissa.
*
*        EPS  - REAL.
*               The smallest positive number such that
*
*                  fl( 1.0 - EPS ) .lt. 1.0,
*
*               where fl denotes the computed value.
*
*        EMIN - INTEGER.
*               The minimum exponent before (gradual) underflow occurs.
*
*        RMIN - REAL.
*               The smallest normalized number for the machine given by
*               BASE**( EMIN - 1 ), where BASE is the floating point
*               value of BETA.
*
*
*     The computation of EPS, EMIN and RMIN is based on a routine,
*     PARANOIA by W. Kahan of the University of California at Berkeley.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine (ENVIRN).
*
*  -- Written on 6-January-1986.
*     Sven Hammarling, Mick Pont and Janet Welding, Nag Central Office.
 
      EXTERNAL           MCENV1, MCMIN , MCSTOR
      INTRINSIC          ABS   , MAX
      LOGICAL            FIRST , IWARN , LRND
      INTEGER            GNMIN , GPMIN , I     , LBETA , LEMIN , LEMIN1
      INTEGER            LEMIN2, LT    , NGNMIN, NGPMIN
      DOUBLE PRECISION   A     , B     , C     , HALF  , LEPS  , LRMIN
      DOUBLE PRECISION   ONE   , RBASE , SIXTH , SMALL , MCSTOR, THIRD
      DOUBLE PRECISION   TWO   , XBASE , ZERO
 
      COMMON    /SOL1CM/ NOUT
 
      SAVE               FIRST , IWARN , LBETA , LEMIN , LEPS  , LRMIN
      SAVE               LT
      DATA               FIRST / .TRUE. /      , IWARN / .FALSE. /
 
      IF( FIRST )THEN
         FIRST = .FALSE.
         ZERO  = 0
         ONE   = 1
         TWO   = 2
 
*        LBETA, LT, LEPS, LEMIN and LRMIN are the local values of BETA,
*        T, EPS, EMIN and RMIN.
*
*        Throughout this routine we use the function MCSTOR to ensure
*        that relevant values are stored and not held in registers, or
*        are not affected by optimizers.
*
*        MCENV1 returns the parameters LBETA and LT. ( LRND is not used
*        here. )
 
         CALL MCENV1( LBETA, LT, LRND )
 
*        Start to find EPS.
 
         B    = LBETA
         A    = B**( -LT )
         LEPS = A
 
*        Try some tricks to see whether or not this is the correct EPS.
 
         B     = TWO/3
         HALF  = ONE/2
         SIXTH = MCSTOR( B    , -HALF  )
         THIRD = MCSTOR( SIXTH,  SIXTH )
         B     = MCSTOR( THIRD, -HALF  )
         B     = MCSTOR( B    ,  SIXTH )
         B     = ABS   ( B )
         IF( B.LT.LEPS )
     $      B = LEPS
 
         LEPS = 1
 
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    IF   ( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )THEN
            LEPS = B
            C    = MCSTOR( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C    = MCSTOR( HALF     , -C                     )
            B    = MCSTOR( HALF     ,  C                     )
            C    = MCSTOR( HALF     , -B                     )
            B    = MCSTOR( HALF     ,  C                     )
            GO TO 10
         END IF
*+       END WHILE
 
         IF( A.LT.LEPS )
     $      LEPS = A
 
*        Computation of EPS complete. Now find EMIN.
*        First compute the next floating point value below 1.0, a, and
*        keep dividing a by BETA until (gradual) underflow occurs.
*        This is detected when we cannot recover the previous a.
 
         XBASE = LBETA
         RBASE = 1/XBASE
         SMALL = ONE
         DO 20, I = 1, LT - 1
            SMALL = MCSTOR( SMALL/LBETA, ZERO )
   20    CONTINUE
         A     = MCSTOR( ONE, SMALL )
         CALL MCMIN ( NGPMIN,  ONE, XBASE, RBASE, LBETA )
         CALL MCMIN ( NGNMIN, -ONE, XBASE, RBASE, LBETA )
         CALL MCMIN (  GPMIN,    A, XBASE, RBASE, LBETA )
         CALL MCMIN (  GNMIN,   -A, XBASE, RBASE, LBETA )
         LEMIN = 0
         IF( ( NGPMIN.EQ.NGNMIN ).AND.( GPMIN.EQ.GNMIN ) )THEN
            IF( NGPMIN.EQ.GPMIN )THEN
               LEMIN = NGPMIN
            ELSE IF( NGPMIN.LT.GPMIN )THEN
               IF( ABS( GPMIN - NGPMIN - LT ).LT.3 )THEN
                  LEMIN =  GPMIN
               ELSE
                  LEMIN =  NGPMIN
                  IWARN = .TRUE.
               END IF
            ELSE
               WRITE( NOUT, 9999 )
               STOP
            END IF
         ELSE
            IF( NGPMIN.EQ.GPMIN )THEN
               LEMIN1 = NGPMIN
            ELSE IF( NGPMIN.LT.GPMIN )THEN
               IF( ABS( GPMIN - NGPMIN - LT ).LT.3 )THEN
                  LEMIN1 =  GPMIN
               ELSE
                  LEMIN1 =  NGPMIN
                  IWARN  = .TRUE.
               END IF
            ELSE
               WRITE( NOUT, 9999 )
               STOP
            END IF
            IF( NGNMIN.EQ.GNMIN )THEN
               LEMIN2 = NGNMIN
            ELSE IF( NGNMIN.LT.GNMIN )THEN
               IF( ABS( GNMIN - NGNMIN - LT ).LT.3 )THEN
                  LEMIN2 =  GNMIN
               ELSE
                  LEMIN2 =  NGNMIN
                  IWARN  = .TRUE.
               END IF
            ELSE
               WRITE( NOUT, 9999 )
               STOP
            END IF
            LEMIN = MAX( LEMIN1, LEMIN2 )
         END IF
***
* Comment out this IF block if Emin is ok
         IF( IWARN )THEN
            FIRST = .TRUE.
            WRITE( NOUT, 9998 )LEMIN
         END IF
***
 
*        Finally compute RMIN by successive division by BETA.
*        We could compute RMIN as base**( EMIN - 1 ), but some machines
*        underflow during this computation.
 
         LRMIN = 1
         DO 30, I = 1, 1 - LEMIN
            LRMIN = LRMIN/LBETA
   30    CONTINUE
      END IF
 
      BETA = LBETA
      T    = LT
      EPS  = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      RETURN
 
 9999 FORMAT( // ' ** ERROR . No reliable value for Emin could be',
     $           ' found.' / ' Please contact Stanford University.'// )
 9998 FORMAT( // ' WARNING. The value Emin may be incorrect:-  Emin = ',
     $           I8 / ' If after inspection the value Emin looks',
     $           ' acceptable please comment out ' / ' the IF block',
     $           ' as marked within the code of routine mcenv2,' /
     $           ' otherwise contact Stanford University. ' / )
 
*     End of MCENV2 (ENVIRN).
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION  FUNCTION MCSTOR( A, B )
      DOUBLE PRECISION                   A, B
 
*     MCSTOR is intended to force A and B to be stored prior to doing
*     the addition of A and B. For use in situations where optimizers
*     might hold one of these in a register.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine (STORE).
*
*  -- Written on 28-November-1984.
*     Sven Hammarling, Nag Central Office.
 
      MCSTOR  = A + B
 
      RETURN
 
*     End of MCSTOR (STORE).
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      DOUBLE PRECISION FUNCTION MCSMAL()
 
*     MCSMAL is intended to return a small positive value such that the
*     reciprocal of MCSMAL does not overflow.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine (SMALL).
*
*  -- Written on 28-November-1984.
*     Sven Hammarling, Nag Central Office.
 
      EXTERNAL                  MCENV2
      LOGICAL                   FIRST
      INTEGER                   BETA  , EMIN  , T
      DOUBLE PRECISION          BASE  , EPS   , FLMIN , RMIN
 
      SAVE                      FIRST , FLMIN
      DATA                      FIRST / .TRUE. /
 
      IF( FIRST )THEN
         FIRST = .FALSE.
         CALL MCENV2( BETA, T, EPS, EMIN, RMIN )
         BASE  =  BETA
         FLMIN =  RMIN*BASE**4
      END IF
 
      MCSMAL = FLMIN
      RETURN
 
*     End of MCSMAL (SMALL).
 
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE MCMIN ( EMIN, START, XBASE, RBASE, BASE )
      INTEGER            EMIN, BASE
      DOUBLE PRECISION   START, XBASE, RBASE
 
*     Service routine for MCENV2.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine (GETMIN).
*
*  -- Written on 6-January-1986.
*     Sven Hammarling and Mick Pont, Nag Central Office.
 
      EXTERNAL           MCSTOR
      INTEGER            I
      DOUBLE PRECISION   A     , B1    , B2    , C1    , C2    , D1
      DOUBLE PRECISION   D2    , MCSTOR, ZERO
 
      A    = START
      ZERO = 0
      EMIN = 1
      B1   = MCSTOR( A/XBASE, ZERO )
      C1   = A
      C2   = A
      D1   = A
      D2   = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*+   $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 IF   ( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
     $       ( D1.EQ.A ).AND.( D2.EQ.A )      )THEN
         EMIN = EMIN - 1
         A    = B1
         B1   = MCSTOR( A /XBASE, ZERO )
         C1   = MCSTOR( B1*XBASE, ZERO )
         D1   = ZERO
         DO 20, I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2   = MCSTOR( A *RBASE, ZERO )
         C2   = MCSTOR( B2/RBASE, ZERO )
         D2   = ZERO
         DO 30, I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
      RETURN
 
*     End of MCMIN (GETMIN).
 
      END
