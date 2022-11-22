        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:13 2022
        MODULE DBDSQR__genmod
          INTERFACE 
            SUBROUTINE DBDSQR(UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,  &
     &LDC,WORK,INFO)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDU
              INTEGER(KIND=4) :: LDVT
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NCVT
              INTEGER(KIND=4) :: NRU
              INTEGER(KIND=4) :: NCC
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: VT(LDVT,*)
              REAL(KIND=8) :: U(LDU,*)
              REAL(KIND=8) :: C(LDC,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DBDSQR
          END INTERFACE 
        END MODULE DBDSQR__genmod
