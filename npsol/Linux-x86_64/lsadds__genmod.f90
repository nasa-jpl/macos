        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSADDS__genmod
          INTERFACE 
            SUBROUTINE LSADDS(UNITQ,VERTEX,INFORM,K1,K2,NACTIV,NARTIF,NZ&
     &,NFREE,NRANK,NREJTD,NRES,NGQ,N,NQ,NROWA,NROWR,NROWT,ISTATE,KACTIV,&
     &KX,CONDMX,A,R,T,RES,GQ,ZY,WRK1,WRK2)
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWR
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: NQ
              INTEGER(KIND=4) :: N
              LOGICAL(KIND=4) :: UNITQ
              LOGICAL(KIND=4) :: VERTEX
              INTEGER(KIND=4) :: INFORM
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: K2
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NARTIF
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NRANK
              INTEGER(KIND=4) :: NREJTD
              INTEGER(KIND=4) :: NRES
              INTEGER(KIND=4) :: NGQ
              INTEGER(KIND=4) :: ISTATE(*)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: CONDMX
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: R(NROWR,*)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: RES(N,*)
              REAL(KIND=8) :: GQ(N,*)
              REAL(KIND=8) :: ZY(NQ,*)
              REAL(KIND=8) :: WRK1(N)
              REAL(KIND=8) :: WRK2(N)
            END SUBROUTINE LSADDS
          END INTERFACE 
        END MODULE LSADDS__genmod
