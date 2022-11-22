        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSPRT__genmod
          INTERFACE 
            SUBROUTINE LSPRT(PRBTYP,PRNT1,ISDEL,ITER,JADD,JDEL,MSGLVL,  &
     &NACTIV,NFREE,N,NCLIN,NRANK,NROWR,NROWT,NZ,NZ1,ISTATE,ALFA,CONDRZ, &
     &CONDT,GFNORM,GZNORM,GZ1NRM,NUMINF,SUMINF,CTX,SSQ,AX,R,T,X,WORK)
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: PRBTYP
              LOGICAL(KIND=4) :: PRNT1
              INTEGER(KIND=4) :: ISDEL
              INTEGER(KIND=4) :: ITER
              INTEGER(KIND=4) :: JADD
              INTEGER(KIND=4) :: JDEL
              INTEGER(KIND=4) :: MSGLVL
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: ISTATE(*)
              REAL(KIND=8) :: ALFA
              REAL(KIND=8) :: CONDRZ
              REAL(KIND=8) :: CONDT
              REAL(KIND=8) :: GFNORM
              REAL(KIND=8) :: GZNORM
              REAL(KIND=8) :: GZ1NRM
              INTEGER(KIND=4) :: NUMINF
              REAL(KIND=8) :: SUMINF
              REAL(KIND=8) :: CTX
              REAL(KIND=8) :: SSQ
              REAL(KIND=8) :: AX(*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: WORK(N)
            END SUBROUTINE LSPRT
          END INTERFACE 
        END MODULE LSPRT__genmod
