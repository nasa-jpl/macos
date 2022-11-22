        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CMQMUL__genmod
          INTERFACE 
            SUBROUTINE CMQMUL(MODE,N,NZ,NFREE,NQ,UNITQ,KX,V,ZY,WRK)
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: MODE
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NFREE
              LOGICAL(KIND=4) :: UNITQ
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WRK(N)
            END SUBROUTINE CMQMUL
          END INTERFACE 
        END MODULE CMQMUL__genmod
