        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:01 2022
        MODULE LSMULS__genmod
          INTERFACE 
            SUBROUTINE LSMULS(PRBTYP,MSGLVL,N,NACTIV,NFREE,NROWA,NROWT, &
     &NUMINF,NZ,NZ1,ISTATE,KACTIV,KX,DINKY,JSMLST,KSMLST,JINF,JTINY,    &
     &JBIGST,KBIGST,TRULAM,A,ANORMS,GQ,RLAMDA,T,WTINF)
              INTEGER(KIND=4) :: NROWT
              INTEGER(KIND=4) :: NROWA
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=2) :: PRBTYP
              INTEGER(KIND=4) :: MSGLVL
              INTEGER(KIND=4) :: NACTIV
              INTEGER(KIND=4) :: NFREE
              INTEGER(KIND=4) :: NUMINF
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NZ1
              INTEGER(KIND=4) :: ISTATE(*)
              INTEGER(KIND=4) :: KACTIV(N)
              INTEGER(KIND=4) :: KX(N)
              REAL(KIND=8) :: DINKY
              INTEGER(KIND=4) :: JSMLST
              INTEGER(KIND=4) :: KSMLST
              INTEGER(KIND=4) :: JINF
              INTEGER(KIND=4) :: JTINY
              INTEGER(KIND=4) :: JBIGST
              INTEGER(KIND=4) :: KBIGST
              REAL(KIND=8) :: TRULAM
              REAL(KIND=8) :: A(NROWA,*)
              REAL(KIND=8) :: ANORMS(*)
              REAL(KIND=8) :: GQ(N)
              REAL(KIND=8) :: RLAMDA(N)
              REAL(KIND=8) :: T(NROWT,*)
              REAL(KIND=8) :: WTINF(*)
            END SUBROUTINE LSMULS
          END INTERFACE 
        END MODULE LSMULS__genmod
