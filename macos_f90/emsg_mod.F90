!***********************************************************************
!  emsg_mod -- the engine's message buffer (PLAN_CONSOLIDATION 1.2 (ee)).
!
!  The engine's warnings, notes and refusals are Fortran unit-6 text, and
!  unit 6 does not reach file descriptor 1 across the FFI on macOS: the
!  bindings (a nested MATLAB, pytest's capfd) read nothing there.  So the
!  load path, CALIB's cause lines and the SMACOS loop's 'Unknown command'
!  print through this module: each line goes to unit 6 exactly as before
!  AND is kept here, and the api (messages_count / messages_get / messages_clear,
!  macos_api_mod) hands the kept lines to the bindings.
!
!  Two ways in:
!    - a converted WRITE(*,fmt) list  is  WRITE(emsgRec,fmt)
!      list  followed by  CALL emsg_put : the internal write fills one or
!      more records of emsgRec (a '/' in the format, a long list-directed
!      list), emsg_put prints and keeps each record written and re-arms
!      them.  The records start as CHAR(0), so a CALL emsg_put after a
!      WRITE that did not run (a one-line IF) does nothing.
!    - new code calls  emsg(text)  with the line itself.
!
!  The buffer keeps the last mEMsg lines (lEMsg characters each) since the
!  last emsg_clear; older ones are counted as dropped.  Unit-6 output is
!  unchanged except that a list-directed line is no longer wrapped at 80
!  columns by ifx (an internal record is as wide as emsgRec).
!***********************************************************************
      MODULE emsg_mod
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: emsgRec, emsg_put, emsg, emsg_clear, &
                  emsg_stats, emsg_line

        INTEGER, PARAMETER :: lEMsgRec = 1024, nEMsgRec = 16
        CHARACTER(len=lEMsgRec) :: emsgRec(nEMsgRec) = CHAR(0)

        INTEGER, PARAMETER :: mEMsg = 500, lEMsg = 256
        CHARACTER(len=lEMsg) :: emsgBuf(0:mEMsg-1)
        ! lines kept since the last clear (the ring holds the last mEMsg)
        INTEGER :: nEMsgTot = 0

      CONTAINS

        SUBROUTINE emsg_put()
          INTEGER :: i
          DO i = 1, nEMsgRec
            IF (emsgRec(i)(1:1) == CHAR(0)) EXIT
            WRITE(*,'(A)') TRIM(emsgRec(i))
            CALL emsg_keep(emsgRec(i))
            emsgRec(i) = CHAR(0)
          END DO
        END SUBROUTINE emsg_put

        SUBROUTINE emsg(text)
          CHARACTER(len=*), INTENT(in) :: text
          WRITE(*,'(A)') TRIM(text)
          CALL emsg_keep(text)
        END SUBROUTINE emsg

        SUBROUTINE emsg_keep(text)
          CHARACTER(len=*), INTENT(in) :: text
          emsgBuf(MOD(nEMsgTot, mEMsg)) = text
          nEMsgTot = nEMsgTot + 1
        END SUBROUTINE emsg_keep

        SUBROUTINE emsg_clear()
          nEMsgTot = 0
        END SUBROUTINE emsg_clear

        ! nLines kept, nChars of the kept lines joined by LF (one LF after
        ! each, trailing blanks trimmed), nDropped older lines not kept
        SUBROUTINE emsg_stats(nLines, nChars, nDropped)
          INTEGER, INTENT(out) :: nLines, nChars, nDropped
          INTEGER :: k
          nLines   = MIN(nEMsgTot, mEMsg)
          nDropped = nEMsgTot - nLines
          nChars   = 0
          DO k = 1, nLines
            nChars = nChars + LEN_TRIM(emsg_line(k)) + 1
          END DO
        END SUBROUTINE emsg_stats

        ! the k-th kept line, oldest first (1 <= k <= nLines)
        FUNCTION emsg_line(k) RESULT(line)
          INTEGER, INTENT(in) :: k
          CHARACTER(len=lEMsg) :: line
          INTEGER :: nLines
          nLines = MIN(nEMsgTot, mEMsg)
          line = ' '
          IF (k < 1 .OR. k > nLines) RETURN
          line = emsgBuf(MOD(nEMsgTot - nLines + k - 1, mEMsg))
        END FUNCTION emsg_line

      END MODULE emsg_mod
