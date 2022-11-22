        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSSOL__genmod
          INTERFACE 
            SUBROUTINE LSSOL(MM,N,NCLIN,NROWA,NROWR,A,BL,BU,CVEC,ISTATE,&
     &KX,X,R,B,INFORM,ITER,OBJ,CLAMDA,IW,LENIW,W,LENW)
              INTEGER(KIND=4) :: LENW
              INTEGER(KIND=4) :: LENIW
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: MM
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: BL(N+NCLIN)
              REAL(KIND=8) :: BU(N+NCLIN)
              REAL(KIND=8) :: CVEC(*)
              INTEGER(KIND=4) :: ISTATE(N+NCLIN)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: B(*)
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: ITER
              REAL(KIND=8) :: OBJ
              REAL(KIND=8) :: CLAMDA(N+NCLIN)
              INTEGER(KIND=4) :: IW(LENIW)
              REAL(KIND=8) :: W(LENW)
            END SUBROUTINE LSSOL
          END INTERFACE 
        END MODULE LSSOL__genmod
