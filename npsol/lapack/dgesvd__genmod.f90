        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:14 2022
        MODULE DGESVD__genmod
          INTERFACE 
            SUBROUTINE DGESVD(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,&
     &LWORK,INFO)
              INTEGER(KIND=4) :: LDVT
              INTEGER(KIND=4) :: LDU
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: JOBU
              CHARACTER(LEN=1) :: JOBVT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: U(LDU,*)
              REAL(KIND=8) :: VT(LDVT,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGESVD
          END INTERFACE 
        END MODULE DGESVD__genmod
