        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:25 2022
        MODULE DLARFT__genmod
          INTERFACE 
            SUBROUTINE DLARFT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: DIRECT
              CHARACTER(LEN=1) :: STOREV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: T(LDT,*)
            END SUBROUTINE DLARFT
          END INTERFACE 
        END MODULE DLARFT__genmod
