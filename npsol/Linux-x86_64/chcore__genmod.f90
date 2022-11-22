        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:00 2022
        MODULE CHCORE__genmod
          INTERFACE 
            SUBROUTINE CHCORE(DEBUG,DONE,FIRST,EPSA,EPSR,FX,X,INFORM,   &
     &ITER,ITMAX,CDEST,FDEST,SDEST,ERRBND,F1,F2,H,HOPT,HPHI)
              LOGICAL(KIND=4) :: DEBUG
              LOGICAL(KIND=4) :: DONE
              LOGICAL(KIND=4) :: FIRST
              REAL(KIND=8) :: EPSA
              REAL(KIND=8) :: EPSR
              REAL(KIND=8) :: FX
              REAL(KIND=8) :: X
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: ITER
              INTEGER(KIND=4) :: ITMAX
              REAL(KIND=8) :: CDEST
              REAL(KIND=8) :: FDEST
              REAL(KIND=8) :: SDEST
              REAL(KIND=8) :: ERRBND
              REAL(KIND=8) :: F1
              REAL(KIND=8) :: F2
              REAL(KIND=8) :: H
              REAL(KIND=8) :: HOPT
              REAL(KIND=8) :: HPHI
            END SUBROUTINE CHCORE
          END INTERFACE 
        END MODULE CHCORE__genmod
