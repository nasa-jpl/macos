        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 21 18:43:02 2022
        MODULE NPMRT__genmod
          INTERFACE 
            SUBROUTINE NPMRT(FEASQP,N,NCLIN,NCNLN,OBJALF,GRDALF,QPCURV, &
     &ISTATE,CJDX,CMUL,CS,DLAM,RHO,VIOLN,WORK1,WORK2)
              LOGICAL(KIND=4) :: FEASQP
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NCLIN
              INTEGER(KIND=4) :: NCNLN
              REAL(KIND=8) :: OBJALF
              REAL(KIND=8) :: GRDALF
              REAL(KIND=8) :: QPCURV
              INTEGER(KIND=4) :: ISTATE(*)
              REAL(KIND=8) :: CJDX(*)
              REAL(KIND=8) :: CMUL(*)
              REAL(KIND=8) :: CS(*)
              REAL(KIND=8) :: DLAM(*)
              REAL(KIND=8) :: RHO(*)
              REAL(KIND=8) :: VIOLN(*)
              REAL(KIND=8) :: WORK1(*)
              REAL(KIND=8) :: WORK2(*)
            END SUBROUTINE NPMRT
          END INTERFACE 
        END MODULE NPMRT__genmod
