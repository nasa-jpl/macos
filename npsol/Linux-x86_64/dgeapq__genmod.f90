        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE DGEAPQ__genmod
          INTERFACE 
            SUBROUTINE DGEAPQ(TRANS,WHEREZ,M,N,A,LDA,ZETA,NCOLB,B,LDB,  &
     &WORK,INFORM)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: WHEREZ
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: ZETA(*)
              INTEGER(KIND=4) :: NCOLB
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFORM
            END SUBROUTINE DGEAPQ
          END INTERFACE 
        END MODULE DGEAPQ__genmod
