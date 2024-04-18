
! #include "realtype.h"

MODULE macos_IO
  USE Kinds



  INTEGER(1), PARAMETER :: PrtCmdIndent = 17      ! Format output style for *.in type output (Rx)
  CHARACTER(LEN=510)    :: PrtMsgStr

  INTEGER :: PrtToFile = 0    ! if not set, print all information to screen (with specific exceptions)
                                 ! otherwise, if set, it defines the File ID.
                                 ! affected routines:


  INTERFACE PrintFmtArray
    MODULE PROCEDURE PrtFmtArray_dbl, &
                     PrtFmtArray_int
  END INTERFACE

  INTERFACE PrintFmt
    MODULE PROCEDURE PrtFmtScalar_dbl, &
                     PrtFmtScalar_int, &
                     PrtFmtScalar_str
  END INTERFACE


  CHARACTER(LEN=40), PRIVATE :: FmtStr
  CHARACTER(LEN=20), PRIVATE :: cmd, cmdE

  CONTAINS

    !
    ! ------------------------------------------------------------------
    !
    LOGICAL FUNCTION IsPrtToFile()

      IsPrtToFile = (PrtToFile /= 0)

    END FUNCTION IsPrtToFile
    !
    ! ------------------------------------------------------------------
    !
    SUBROUTINE SetPrtToFileID(ID)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ID
      ! - - - - - - - - - - - - - - - - - - -
      PrtToFile = ID

    END SUBROUTINE
    !
    ! ------------------------------------------------------------------
    !
    ! This routine prints a string on the screen
    !
    ! Matlab:    This routine prints a string on the screen and in the diary (if the diary
    !            is in use). It provides a callback to the standard C printf routine already
    !            linked inside MATLAB software
    !
    ! IMPORTANT: In MATLAB, if you want the literal % in your message, use %% in the message
    !            string since % has special meaning to printf. Failing to do so causes
    !            unpredictable results.
    !
    SUBROUTINE PrintMsg(MsgStr)

      IMPLICIT NONE
      CHARACTER*(*),   INTENT(IN):: MsgStr
      ! - - - - - - - - - - - - - - - - - - -
      IF (IsPrtToFile()) THEN
        WRITE(PrtToFile,'(A)') TRIM(MsgStr)
      ELSE
        WRITE(*,'(A)') TRIM(MsgStr)
      END IF

    END SUBROUTINE PrintMsg
    ! !
    ! ! ------------------------------------------------------------------
    ! !
    ! SUBROUTINE PrintMsg_Level(MsgStr,VerboseLevel)

    !   IMPLICIT NONE
    !   CHARACTER*(*),INTENT(IN) :: MsgStr
    !   INTEGER(1),   INTENT(IN) :: VerboseLevel
    !   ! - - - - - - - - - - - - - - - - - - -

    !   IF (VerboseLevel<=IO_PrintLevel_Set) THEN

    !     CALL PrintMsg_On(MsgStr)

    !   END IF

    ! END SUBROUTINE PrintMsg_Level


    ! ------------------------------------------------------------------
    !
    SUBROUTINE PrtFmtScalar_dbl(NameStr,DigitFmtStr,PrtScalar)
      USE Kinds

      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN):: DigitFmtStr, NameStr
      REAL(pr),         INTENT(IN):: PrtScalar
      ! - - - - - - - - - - - - - - - - - - -
      WRITE(cmd, "('(A',i0,A,',')") PrtCmdIndent,",'=',"    ! "(Axx,'=',"

      WRITE(PrtMsgStr,cmd(1:LEN_TRIM(cmd))//TRIM(DigitFmtStr)//")") TRIM(NameStr), PrtScalar
      CALL PrintMsg(PrtMsgStr)

    END SUBROUTINE PrtFmtScalar_dbl
    !
    ! ------------------------------------------------------------------
    !
    SUBROUTINE PrtFmtScalar_int(NameStr,DigitFmtStr,PrtScalar)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN):: DigitFmtStr, NameStr
      INTEGER,          INTENT(IN):: PrtScalar
      ! - - - - - - - - - - - - - - - - - - -
      WRITE(cmd, "('(A',i0,A,',')") PrtCmdIndent,",'=',"    ! "(Axx,'=',"

      WRITE(PrtMsgStr,cmd(1:LEN_TRIM(cmd))//TRIM(DigitFmtStr)//")") TRIM(NameStr), PrtScalar
      CALL PrintMsg(PrtMsgStr)

    END SUBROUTINE PrtFmtScalar_int
    !
    ! ------------------------------------------------------------------
    !
    SUBROUTINE PrtFmtScalar_str(NameStr,DigitFmtStr,PrtStr)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN):: DigitFmtStr, NameStr
        CHARACTER*(*),    INTENT(IN):: PrtStr
        ! - - - - - - - - - - - - - - - - - - -
        WRITE(cmd, "('(A',i0,A,',')") PrtCmdIndent,",'=',"    ! "(Axx,'=',"

        WRITE(PrtMsgStr,cmd(1:LEN_TRIM(cmd))//TRIM(DigitFmtStr)//")") TRIM(NameStr), TRIM(PrtStr)
        CALL PrintMsg(PrtMsgStr)

      END SUBROUTINE PrtFmtScalar_str

    ! ------------------------------------------------------------------
    !
    ! Print formated vector  -- e.g.:CALL PrintFmtArray('  AnaCoef','D23.15',3,AnaCoef(:))
    !
    SUBROUTINE PrtFmtArray_dbl(NameStr,DigitFmtStr,nCols,PrtVec)
        USE Kinds

        IMPLICIT NONE
        CHARACTER(LEN=*),         INTENT(IN):: DigitFmtStr, NameStr
        INTEGER,                  INTENT(IN):: nCols
        REAL(pr),                 INTENT(IN):: PrtVec(:)

        CHARACTER(LEN=255) :: PrtMsgStr
        INTEGER            :: nVec,k,j,u
        CHARACTER(LEN=40)  :: FmtStr
        CHARACTER(LEN=20)  :: cmd, cmdE

        INTEGER            :: iZ, nZ
        ! - - - - - - - - - - - - - - - - - - -

        nVec  = SIZE(PrtVec)
        IF (nVec==0) RETURN


        if (LEN_TRIM(NameStr)==0) then
          WRITE(cmd,"('(A',i0,A,',')")PrtCmdIndent,",' ',"    ! "(Axx,' ',"
        else
          WRITE(cmd, "('(A',i0,A,',')")PrtCmdIndent,",'=',"   ! "(Axx,'=',"
        end if
        u = PrtCmdIndent+2


        IF (nVec<=nCols) THEN
          WRITE(FmtStr,"(i0,a)") nVec,'('//TRIM(DigitFmtStr)//'))'   ! eg: "(3(1PD23.15))"
          WRITE(PrtMsgStr,cmd(1:u)//FmtStr) TRIM(NameStr),PrtVec(:)  !     "   cmd = (3(1PD23.15))"
          CALL PrintMsg(PrtMsgStr)

        ELSE
          ! write 1st row with 'nCols' columns
          WRITE(FmtStr,"(i0,a)") nCols,'('//TRIM(DigitFmtStr)//'))'   ! eg: "(3(1PD23.15))"
          WRITE(PrtMsgStr,cmd(1:u)//FmtStr) TRIM(NameStr),PrtVec(:nCols)   !     "   cmd = (3(1PD23.15))"
          CALL PrintMsg(PrtMsgStr)

          ! write complete rows with 'nCols' columns
          WRITE(FmtStr,"(a,i0,a,i0,a)") '(',PrtCmdIndent+1,'x,',nCols,'('//TRIM(DigitFmtStr)//'))'
          DO j=2, nVec/nCols
            WRITE(PrtMsgStr,FmtStr)PrtVec((j-1)*nCols+1:j*nCols)
            CALL PrintMsg(PrtMsgStr)
          END DO

          ! write incomplete rows
          nZ = MOD(nVec, nCols)
          iZ = nVec/nCols
          IF (nZ>0) THEN
            WRITE(FmtStr,"(a,i0,a,i0,a)") '(',PrtCmdIndent+1,'x,',nZ,'('//TRIM(DigitFmtStr)//'))'
            WRITE(PrtMsgStr,FmtStr)PrtVec(iZ*nCols+1:nVec)
            CALL PrintMsg(PrtMsgStr)
          END IF

        END IF

      END SUBROUTINE PrtFmtArray_dbl
      !
      ! ------------------------------------------------------------------
      !
      ! Print formated vector  -- e.g.:CALL PrtFmtArray('  AnaCoef','1PI4',3,AnaCoef(:))
      !
      SUBROUTINE PrtFmtArray_int(NameStr,DigitFmtStr,nCols,PrtVec)

        IMPLICIT NONE
        CHARACTER(LEN=*),         INTENT(IN):: DigitFmtStr, NameStr
        INTEGER,                  INTENT(IN):: nCols
        INTEGER,                  INTENT(IN):: PrtVec(:)

        CHARACTER(LEN=255) :: PrtMsgStr
        INTEGER            :: nVec,k,j,u
        CHARACTER(LEN=40)  :: FmtStr
        CHARACTER(LEN=20)  :: cmd

        INTEGER            :: iZ, nZ
        ! - - - - - - - - - - - - - - - - - - -

        nVec  = SIZE(PrtVec)
        IF (nVec==0) RETURN

        WRITE(cmd, "('(A',i0,A,',')")PrtCmdIndent,",'=',";    u = LEN_TRIM(cmd)  ! "(Axx,'=',"

        IF (nVec<=nCols) THEN
            WRITE(FmtStr,"(i0,a)") nVec,'('//TRIM(DigitFmtStr)//'))'   ! eg: "(3(1PD23.15))"
            WRITE(PrtMsgStr,cmd(1:u)//FmtStr) TRIM(NameStr),PrtVec(:)  !     "   cmd = (3(1PD23.15))"
            CALL PrintMsg(PrtMsgStr)

          ELSE
            ! write 1st row with 'nCols' columns
            WRITE(FmtStr,"(i0,a)") nCols,'('//TRIM(DigitFmtStr)//'))'   ! eg: "(3(1PD23.15))"
            WRITE(PrtMsgStr,cmd(1:u)//FmtStr) TRIM(NameStr),PrtVec(:nCols)   !     "   cmd = (3(1PD23.15))"
            CALL PrintMsg(PrtMsgStr)

            ! write complete rows with 'nCols' columns
            WRITE(FmtStr,"(a,i0,a,i0,a)") '(',PrtCmdIndent+1,'x,',nCols,'('//TRIM(DigitFmtStr)//'))'
            DO j=2, nVec/nCols
              WRITE(PrtMsgStr,FmtStr)PrtVec((j-1)*nCols+1:j*nCols)
              CALL PrintMsg(PrtMsgStr)
            END DO

            ! write incomplete rows
            nZ = MOD(nVec, nCols)
            iZ = nVec/nCols
            IF (nZ>0) THEN
              WRITE(FmtStr,"(a,i0,a,i0,a)") '(',PrtCmdIndent+1,'x,',nZ,'('//TRIM(DigitFmtStr)//'))'
              WRITE(PrtMsgStr,FmtStr)PrtVec(iZ*nCols+1:nVec)
              CALL PrintMsg(PrtMsgStr)
            END IF

          END IF

        ! IF (nRows<=1) THEN

        !  !WRITE(FmtStr,"(a,a,i0,a)") "('"//TRIM(NameStr),"=',",nCols,TRIM(DigitFmtStr)//')'
        !  !WRITE(PrtMsgStr,FmtStr)PrtVec(:)

        !   WRITE(FmtStr,"(i0,a)") nCols, '('//TRIM(DigitFmtStr)//'))'
        !   WRITE(PrtMsgStr,cmd(1:u)//FmtStr) TRIM(NameStr),PrtVec(:)
        !   CALL PrintMsg(PrtMsgStr)

        ! ELSE

        !  !WRITE(FmtStr,"(a,a,i0,a)") "('"//TRIM(NameStr),"=',",nCols,TRIM(DigitFmtStr)//')'
        !  !WRITE(PrtMsgStr,FmtStr)PrtVec(1:nCols)

        !   WRITE(FmtStr,"(i0,a)") nCols,TRIM(DigitFmtStr)//')'
        !   WRITE(PrtMsgStr,cmd(1:u)//FmtStr) TRIM(NameStr),PrtVec(1:nCols)
        !   CALL PrintMsg(PrtMsgStr)

        !   WRITE(FmtStr,"(a,i0,a,i0,a)") '(',PrtCmdIndent+1,'x,',nCols,TRIM(DigitFmtStr)//')'
        !   DO j=2,nRows
        !     WRITE(PrtMsgStr,FmtStr)PrtVec((j-1)*nCols+1:j*nCols)
        !     CALL PrintMsg(PrtMsgStr)
        !   END DO

        !   k=MOD(nVec,nCols)
        !   IF (nCols/=1 .AND. k/=0) THEN
        !     WRITE(FmtStr,"(a,i0,a,i0,a)") '(',PrtCmdIndent+1,'x,',k,TRIM(DigitFmtStr)//')'
        !     WRITE(PrtMsgStr,FmtStr)PrtVec((j-1)*nCols+1:(j-1)*nCols+k)
        !     CALL PrintMsg(PrtMsgStr)
        !   END IF

        ! END IF

      END SUBROUTINE PrtFmtArray_int


END MODULE macos_IO
