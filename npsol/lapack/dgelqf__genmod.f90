        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:13 2022
        MODULE DGELQF__genmod
          INTERFACE 
            SUBROUTINE DGELQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGELQF
          END INTERFACE 
        END MODULE DGELQF__genmod