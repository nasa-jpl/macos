! ----------------------------------------------------------------------
! validate_prescription.F90
!
! Pre-validate a MACOS .in prescription file using pure character I/O
! before the parser (msmacosio.inc) gets a crack at it.
!
! Phase 1 catches the most common breakages:
!   - a 'Key=' line with no value (on the same line or as a continuation)
!   - a blank line breaking a multi-row continuation block (e.g. the
!     middle of a Tout matrix), which trips the parser's fixed-format
!     READ even though every individual line is well-formed
!   - file-not-found / cannot-open
! More elaborate checks (enum-value validation, expected array lengths)
! belong in Phase 2.
!
! Per-key exceptions:
!   - EltName is allowed to be empty (Norbert request 2026-05-08).
!     Some prescriptions have unnamed elements; the parser tolerates
!     this so the validator should too.
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
          CHARACTER(LEN=64)   :: pending_key, cont_key
          INTEGER :: lineno, pending_lineno, eqpos, lerr, i, n
          INTEGER :: cont_key_lineno, blank_lineno
          LOGICAL :: pending, fexist, prev_was_cont, saw_blank

          ios = 0
          msg = ''
          pending = .FALSE.
          pending_key = ''
          pending_lineno = 0
          lineno = 0
          ! State for "blank line inside multi-row continuation" check
          prev_was_cont = .FALSE.
          saw_blank = .FALSE.
          cont_key = ''
          cont_key_lineno = 0
          blank_lineno = 0

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
            IF (n == 0) THEN
              ! Blank line.  Note it for the "blank inside multi-row
              ! continuation" check; final verdict comes when we see
              ! the next non-blank line.
              IF (prev_was_cont .AND. .NOT. saw_blank) THEN
                saw_blank = .TRUE.
                blank_lineno = lineno
              END IF
              CYCLE
            END IF
            DO i = 1, n
              IF (line(i:i) /= ' ' .AND. ICHAR(line(i:i)) /= 9) EXIT
            END DO
            IF (i > n) THEN
              ! All-whitespace line (treat like blank).
              IF (prev_was_cont .AND. .NOT. saw_blank) THEN
                saw_blank = .TRUE.
                blank_lineno = lineno
              END IF
              CYCLE
            END IF
            IF (line(i:i) == '%' .OR. line(i:i) == '!') CYCLE

            eqpos = INDEX(line, '=')
            IF (eqpos > 0) THEN
              ! New 'Key=...' line.  Reset the multi-row tracking state:
              ! a blank line just before a new key is fine -- it just
              ! separates element blocks.
              prev_was_cont = .FALSE.
              saw_blank     = .FALSE.

              ! If the previous Key= was awaiting a continuation value,
              ! this new key means the previous one never got one.
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
                ! Empty value on this line.  Most keys must get a
                ! value via continuation (or by EOF this is an error),
                ! but a few prescription keys are intentionally allowed
                ! to be empty -- handle those here without setting
                ! `pending`.
                IF (KeyEq(pending_key, 'EltName')) THEN
                  pending = .FALSE.
                ELSE
                  pending = .TRUE.
                  pending_lineno = lineno
                END IF
              ELSE
                pending = .FALSE.
                ! Inline value present: this key is now eligible to
                ! own subsequent continuation lines for the multi-row
                ! check.
                cont_key        = pending_key
                cont_key_lineno = lineno
              END IF
            ELSE
              ! Continuation line (no '=').
              !   If we saw a blank since the previous continuation,
              !   the parser's fixed-format multi-row READ for the
              !   originating key will be off by one or hit EOF.
              IF (prev_was_cont .AND. saw_blank) THEN
                ios = 1
                WRITE(msg, '(A,I0,A,A,A,I0,A)') &
                    'line ', blank_lineno, &
                    ': blank line inside multi-row block for key "', &
                    TRIM(cont_key), '" (started at line ', &
                    cont_key_lineno, ')'
                CLOSE(vunit)
                RETURN
              END IF
              IF (pending) pending = .FALSE.
              prev_was_cont = .TRUE.
              saw_blank     = .FALSE.
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

        ! Case-insensitive equality between two trimmed key strings.
        ! Used to recognize specific keyword exceptions regardless of
        ! how the user cased them in the .in file.
        LOGICAL FUNCTION KeyEq(a, b)
          CHARACTER(*), INTENT(IN) :: a, b
          INTEGER :: i, na, nb, ca, cb
          na = LEN_TRIM(a)
          nb = LEN_TRIM(b)
          IF (na /= nb) THEN
            KeyEq = .FALSE.
            RETURN
          END IF
          DO i = 1, na
            ca = IACHAR(a(i:i))
            cb = IACHAR(b(i:i))
            IF (ca >= IACHAR('a') .AND. ca <= IACHAR('z'))   &
              ca = ca - 32
            IF (cb >= IACHAR('a') .AND. cb <= IACHAR('z'))   &
              cb = cb - 32
            IF (ca /= cb) THEN
              KeyEq = .FALSE.
              RETURN
            END IF
          END DO
          KeyEq = .TRUE.
        END FUNCTION KeyEq

      END MODULE validate_prescription_mod
