        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:40:09 2022
        MODULE DGEMM__genmod
          INTERFACE 
            SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,&
     &C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: TRANSB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: C(LDC,*)
            END SUBROUTINE DGEMM
          END INTERFACE 
        END MODULE DGEMM__genmod
