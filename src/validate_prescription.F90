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
! Per-key exceptions (empty value accepted):
!   - EltName    -- unnamed element (Norbert 2026-05-08)
!   - BaseUnits  -- unspecified base units, backward compat
!   - WaveUnits  -- unspecified wavelength units, backward compat
!   These keys are routinely left blank in older prescription files
!   and the parser tolerates it, so the validator should too.
!
! Public: validate_prescription_mod%ValidatePrescription
!
! Line endings (2026-10-09, PLAN_CONSOLIDATION item 5; Dave's ruling 4a):
!   ValidatePrescription reports (optional has_cr) whether any record holds
!   a CR -- found in the pass it already makes, so an LF deck is read no
!   more than before.  MBFile6 then loads an LF copy made by RxNormalizeEOL
!   (CRLF and lone CR -> LF) in the system temp dir and deletes it on every
!   exit (RxTmpCleanup).  Why: ifx reads a CR-only (classic Mac) deck as ONE
!   record and fails or crashes on it, gfortran splits on CR; CRLF is
!   already handled by both runtimes.  SAVE writes LF and is untouched.
! ----------------------------------------------------------------------

      MODULE validate_prescription_mod
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: ValidatePrescription, RxNormalizeEOL, RxTmpCleanup
        PUBLIC :: RxTmpFile, nRxNormalized
        CHARACTER(LEN=1024), SAVE :: RxTmpFile = ''   ! the live LF copy
        INTEGER,             SAVE :: nRxNormalized = 0, nRxTmp = 0
      CONTAINS

        SUBROUTINE ValidatePrescription(filename, ios, msg, has_cr)
          CHARACTER(*), INTENT(IN)  :: filename
          INTEGER,      INTENT(OUT) :: ios   ! 0 ok, /=0 bad
          CHARACTER(*), INTENT(OUT) :: msg
          LOGICAL, OPTIONAL, INTENT(OUT) :: has_cr  ! any record holds a CR

          INTEGER, PARAMETER :: vunit = 79
          CHARACTER(LEN=2048) :: line
          CHARACTER(LEN=64)   :: pending_key, cont_key
          INTEGER :: lineno, pending_lineno, eqpos, lerr, i, n
          INTEGER :: cont_key_lineno, blank_lineno
          LOGICAL :: pending, fexist, prev_was_cont, saw_blank
          ! Block comments: '/*' or 'CommentBegin' as a line's first
          ! token opens one, '*/' or 'CommentEnd' closes it -- the same
          ! rule GET_EQ applies (iosub.inc).  Everything inside is
          ! skipped by the parser, so it must be skipped here too, or a
          ! deck the parser reads fine is refused before it gets there
          ! (the eac5mono.in failure: commented-out TElt rows read as a
          ! broken multi-row block).
          LOGICAL :: in_block
          INTEGER :: block_lineno, t2, pct

          ios = 0
          msg = ''
          IF (PRESENT(has_cr)) has_cr = .FALSE.
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
          in_block     = .FALSE.
          block_lineno = 0

          INQUIRE (FILE=filename, EXIST=fexist)
          IF (.NOT. fexist) THEN
            ios = 1
            msg = 'file not found'
            RETURN
          END IF

          ! The CR test peeks at the first 4 KB as BYTES: gfortran's
          ! formatted READ takes a CR as the record end, so a CR never
          ! reaches LINE there (ifx keeps it); the per-record check below
          ! stays for a deck whose CRs start later.
          IF (PRESENT(has_cr)) CALL PeekCR(filename, has_cr)
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
            IF (PRESENT(has_cr)) THEN
              IF (INDEX(line, CHAR(13)) > 0) has_cr = .TRUE.
            END IF

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
              IF (line(i:i) /= ' ' .AND. ICHAR(line(i:i)) /= 9 .AND. &
                  ICHAR(line(i:i)) /= 13) EXIT
            END DO
            IF (i > n) THEN
              ! All-whitespace line (treat like blank).
              IF (prev_was_cont .AND. .NOT. saw_blank) THEN
                saw_blank = .TRUE.
                blank_lineno = lineno
              END IF
              CYCLE
            END IF
            ! First token = up to the next blank / tab / CR.
            t2 = i
            DO WHILE (t2 < n)
              IF (line(t2+1:t2+1) == ' ' .OR. ICHAR(line(t2+1:t2+1)) == 9 &
                  .OR. ICHAR(line(t2+1:t2+1)) == 13) EXIT
              t2 = t2 + 1
            END DO
            IF (line(i:MIN(i+1,n)) == '/*' .OR. &
                KeyEq(line(i:t2), 'CommentBegin')) THEN
              in_block     = .TRUE.
              block_lineno = lineno
              CYCLE
            ELSE IF (line(i:MIN(i+1,n)) == '*/' .OR. &
                     KeyEq(line(i:t2), 'CommentEnd')) THEN
              in_block = .FALSE.
              CYCLE
            END IF
            IF (in_block) CYCLE
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

              pct = INDEX(line(eqpos+1:), '%')
              IF (pct == 0) pct = LEN(line) - eqpos + 1
              IF (LEN_TRIM(line(eqpos+1:eqpos+pct-1)) == 0) THEN
                ! Empty value on this line.  Most keys must get a
                ! value via continuation (or by EOF this is an error),
                ! but a few prescription keys are intentionally allowed
                ! to be empty for backward compatibility -- handle
                ! those here without setting `pending`.
                IF (KeyEq(pending_key, 'EltName')   .OR. &
                    KeyEq(pending_key, 'BaseUnits') .OR. &
                    KeyEq(pending_key, 'WaveUnits')) THEN
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

          IF (in_block) THEN
            ios = 1
            WRITE(msg, '(A,I0,A)') 'line ', block_lineno, &
                ': comment block ("/*" or CommentBegin) is never closed'
            RETURN
          END IF

          IF (pending) THEN
            ios = 1
            CALL FormatMissingValue(msg, pending_lineno, pending_key)
            RETURN
          END IF

          ios = 0
          msg = ''
          RETURN
        END SUBROUTINE ValidatePrescription

        SUBROUTINE PeekCR(filename, has_cr)
          CHARACTER(*), INTENT(IN)    :: filename
          LOGICAL,      INTENT(INOUT) :: has_cr
          CHARACTER(LEN=4096) :: b
          INTEGER :: u, e, nb
          INQUIRE(FILE=filename, SIZE=nb)
          IF (nb <= 0) RETURN
          OPEN(NEWUNIT=u, FILE=filename, ACCESS='STREAM', &
               FORM='UNFORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=e)
          IF (e /= 0) RETURN
          READ(u, IOSTAT=e) b(1:MIN(nb, LEN(b)))
          CLOSE(u)
          IF (e == 0) has_cr = has_cr .OR. &
               (INDEX(b(1:MIN(nb, LEN(b))), CHAR(13)) > 0)
        END SUBROUTINE PeekCR

        ! An LF copy of SRC (CRLF -> LF, lone CR -> LF) in the system temp
        ! dir ($TMPDIR, else TEMP, TMP, /tmp), unique per process + call,
        ! never beside the deck (its directory may be read-only or shared).
        ! DST = its name, also kept in RxTmpFile for RxTmpCleanup.
        SUBROUTINE RxNormalizeEOL(src, dst, ok)
#ifdef __INTEL_COMPILER
          USE IFPORT, ONLY: GETPID
#endif
          CHARACTER(*), INTENT(IN)  :: src
          CHARACTER(*), INTENT(OUT) :: dst
          LOGICAL,      INTENT(OUT) :: ok
          CHARACTER(LEN=:), ALLOCATABLE :: buf, out
          CHARACTER(LEN=1024) :: dir
          INTEGER :: u, e, nb, i, k, ld
          ok = .FALSE.
          dst = ''
          INQUIRE(FILE=src, SIZE=nb)
          IF (nb <= 0) RETURN
          ALLOCATE(CHARACTER(LEN=nb) :: buf, out)
          OPEN(NEWUNIT=u, FILE=src, ACCESS='STREAM', FORM='UNFORMATTED', &
               STATUS='OLD', ACTION='READ', IOSTAT=e)
          IF (e /= 0) RETURN
          READ(u, IOSTAT=e) buf
          CLOSE(u)
          IF (e /= 0) RETURN
          k = 0
          i = 1
          DO WHILE (i <= nb)
            k = k + 1
            IF (buf(i:i) == CHAR(13)) THEN
              out(k:k) = CHAR(10)
              IF (i < nb) THEN
                IF (buf(i+1:i+1) == CHAR(10)) i = i + 1
              END IF
            ELSE
              out(k:k) = buf(i:i)
            END IF
            i = i + 1
          END DO
          dir = ''
          CALL GET_ENVIRONMENT_VARIABLE('TMPDIR', dir, ld)
          IF (ld == 0) CALL GET_ENVIRONMENT_VARIABLE('TEMP', dir, ld)
          IF (ld == 0) CALL GET_ENVIRONMENT_VARIABLE('TMP', dir, ld)
          IF (ld == 0) dir = '/tmp'
          nRxTmp = nRxTmp + 1
          WRITE(dst, '(A,A,I0,A,I0,A)') TRIM(dir), '/macos_rx_', &
               GETPID(), '_', nRxTmp, '.in'
          OPEN(NEWUNIT=u, FILE=TRIM(dst), ACCESS='STREAM', &
               FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE', &
               IOSTAT=e)
          IF (e /= 0) RETURN
          WRITE(u, IOSTAT=e) out(1:k)
          CLOSE(u)
          IF (e /= 0) RETURN
          RxTmpFile = dst
          nRxNormalized = nRxNormalized + 1
          ok = .TRUE.
        END SUBROUTINE RxNormalizeEOL

        ! Delete the LF copy, if any.  Called after every CLOSE of the
        ! deck, on the validator-refusal paths, and at the start of every
        ! load (a backstop for any exit missed).
        SUBROUTINE RxTmpCleanup()
          INTEGER :: u, e
          IF (LEN_TRIM(RxTmpFile) == 0) RETURN
          OPEN(NEWUNIT=u, FILE=TRIM(RxTmpFile), STATUS='OLD', IOSTAT=e)
          IF (e == 0) CLOSE(u, STATUS='DELETE')
          RxTmpFile = ''
        END SUBROUTINE RxTmpCleanup

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
