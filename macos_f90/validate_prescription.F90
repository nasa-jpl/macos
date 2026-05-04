! ----------------------------------------------------------------------
! validate_prescription.F90
!
! Pre-validate a MACOS .in prescription file using pure character I/O
! before the parser (msmacosio.inc) gets a crack at it.
!
! Phase 1 catches the most common breakage: a 'Key=' line with no value
! (on the same line or as a continuation).  File-open errors are also
! reported.  More elaborate checks (enum-value validation, expected
! array lengths) belong in Phase 2.
!
! Public: validate_prescription_mod%ValidatePrescription
! ----------------------------------------------------------------------

      MODULE validate_prescription_mod
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: ValidatePrescription
      CONTAINS

        SUBROUTINE ValidatePrescription(filename, ios, msg)
          CHARACTER(*), INTENT(IN)  :: filename
          INTEGER,      INTENT(OUT) :: ios   ! 0 ok, /=0 bad
          CHARACTER(*), INTENT(OUT) :: msg

          INTEGER, PARAMETER :: vunit = 79
          CHARACTER(LEN=2048) :: line
          CHARACTER(LEN=64)   :: pending_key
          INTEGER :: lineno, pending_lineno, eqpos, lerr, i, n
          LOGICAL :: pending, fexist

          ios = 0
          msg = ''
          pending = .FALSE.
          pending_key = ''
          pending_lineno = 0
          lineno = 0

          INQUIRE (FILE=filename, EXIST=fexist)
          IF (.NOT. fexist) THEN
            ios = 1
            msg = 'file not found'
            RETURN
          END IF

          OPEN (UNIT=vunit, FILE=filename, STATUS='OLD', &
                ACTION='READ', IOSTAT=lerr)
          IF (lerr /= 0) THEN
            ios = 1
            msg = 'cannot open file'
            RETURN
          END IF

          DO
            READ (vunit, '(A)', IOSTAT=lerr) line
            IF (lerr /= 0) EXIT
            lineno = lineno + 1

            n = LEN_TRIM(line)
            IF (n == 0) CYCLE
            DO i = 1, n
              IF (line(i:i) /= ' ' .AND. ICHAR(line(i:i)) /= 9) EXIT
            END DO
            IF (i > n) CYCLE
            IF (line(i:i) == '%' .OR. line(i:i) == '!') CYCLE

            eqpos = INDEX(line, '=')
            IF (eqpos > 0) THEN
              ! New 'Key=...' line.  If the previous Key= was awaiting a
              ! continuation value, this new key means the previous one
              ! never got one.
              IF (pending) THEN
                ios = 1
                CALL FormatMissingValue(msg, pending_lineno, pending_key)
                CLOSE(vunit)
                RETURN
              END IF

              pending_key = ADJUSTL(line(1:eqpos-1))
              IF (LEN_TRIM(pending_key) == 0) THEN
                ios = 1
                WRITE(msg, '(A,I0,A)') &
                    'line ', lineno, ': "=" with no key'
                CLOSE(vunit)
                RETURN
              END IF

              IF (LEN_TRIM(line(eqpos+1:)) == 0) THEN
                pending = .TRUE.
                pending_lineno = lineno
              ELSE
                pending = .FALSE.
              END IF
            ELSE
              ! Continuation line (no '=').  If a key was pending, this
              ! line provides its value.
              IF (pending) pending = .FALSE.
            END IF
          END DO

          CLOSE(vunit)

          IF (pending) THEN
            ios = 1
            CALL FormatMissingValue(msg, pending_lineno, pending_key)
            RETURN
          END IF

          ios = 0
          msg = ''
          RETURN
        END SUBROUTINE ValidatePrescription

        SUBROUTINE FormatMissingValue(msg, lineno, key)
          CHARACTER(*), INTENT(OUT) :: msg
          INTEGER,      INTENT(IN)  :: lineno
          CHARACTER(*), INTENT(IN)  :: key
          WRITE(msg, '(A,I0,A,A,A)') &
              'line ', lineno, ': key "', TRIM(key), '" has no value'
        END SUBROUTINE FormatMissingValue

      END MODULE validate_prescription_mod
