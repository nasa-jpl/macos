        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPFEAS__genmod
          INTERFACE 
            SUBROUTINE NPFEAS(N,NCLIN,NCNLN,ISTATE,BIGBND,CVNORM,ERRMAX,&
     &JMAX,NVIOL,AX,BL,BU,C,FEATOL,X,WORK)
              INTEGER(KIND=4) :: NCNLN
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ISTATE(N+NCLIN+NCNLN)
              REAL(KIND=8) :: BIGBND
              REAL(KIND=8) :: CVNORM
              REAL(KIND=8) :: ERRMAX
              INTEGER(KIND=4) :: JMAX
              INTEGER(KIND=4) :: NVIOL
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: BL(N+NCLIN+NCNLN)
              REAL(KIND=8) :: BU(N+NCLIN+NCNLN)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: FEATOL(N+NCLIN+NCNLN)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WORK(N+NCLIN+NCNLN)
            END SUBROUTINE NPFEAS
          END INTERFACE 
        END MODULE NPFEAS__genmod
