        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:03 2022
        MODULE SRCHC__genmod
          INTERFACE 
            SUBROUTINE SRCHC(DEBUG,DONE,FIRST,IMPRVD,INFORM,ALFMAX,EPSAF&
     &,ETA,XTRY,FTRY,GTRY,OLDF,OLDG,TOLABS,TOLREL,TOLTNY,ALFA,ALFBST,   &
     &FBEST,GBEST)
              LOGICAL(KIND=4) :: DEBUG
              LOGICAL(KIND=4) :: DONE
              LOGICAL(KIND=4) :: FIRST
              LOGICAL(KIND=4) :: IMPRVD
              INTEGER(KIND=4) :: INFORM
              REAL(KIND=8) :: ALFMAX
              REAL(KIND=8) :: EPSAF
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: XTRY
              REAL(KIND=8) :: FTRY
              REAL(KIND=8) :: GTRY
              REAL(KIND=8) :: OLDF
              REAL(KIND=8) :: OLDG
              REAL(KIND=8) :: TOLABS
              REAL(KIND=8) :: TOLREL
              REAL(KIND=8) :: TOLTNY
              REAL(KIND=8) :: ALFA
              REAL(KIND=8) :: ALFBST
              REAL(KIND=8) :: FBEST
              REAL(KIND=8) :: GBEST
            END SUBROUTINE SRCHC
          END INTERFACE 
        END MODULE SRCHC__genmod
