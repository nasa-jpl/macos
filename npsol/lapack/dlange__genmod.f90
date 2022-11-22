        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:41:15 2022
        MODULE DLANGE__genmod
          INTERFACE 
            FUNCTION DLANGE(NORM,M,N,A,LDA,WORK)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: NORM
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: WORK(*)
              REAL(KIND=8) :: DLANGE
            END FUNCTION DLANGE
          END INTERFACE 
        END MODULE DLANGE__genmod
