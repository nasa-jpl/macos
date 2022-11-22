        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPSOL__genmod
          INTERFACE 
            SUBROUTINE NPSOL(N,NCLIN,NCNLN,NROWA,NROWUJ,NROWR,A,BL,BU,  &
     &CONFUN,OBJFUN,INFORM,ITER,ISTATE,C,UJAC,CLAMDA,OBJF,UGRAD,R,X,IW, &
     &LENIW,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: LENIW
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWUJ
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: BL(N+NCLIN+NCNLN)
              REAL(KIND=8) :: BU(N+NCLIN+NCNLN)
              EXTERNAL CONFUN
              EXTERNAL OBJFUN
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: ITER
              INTEGER(KIND=4) :: ISTATE(N+NCLIN+NCNLN)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: UJAC(NROWUJ,*)
              REAL(KIND=8) :: CLAMDA(N+NCLIN+NCNLN)
              REAL(KIND=8) :: OBJF
              REAL(KIND=8) :: UGRAD(N)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: IW(LENIW)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE NPSOL
          END INTERFACE 
        END MODULE NPSOL__genmod
