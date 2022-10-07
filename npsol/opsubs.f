*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  OPSUBS FORTRAN
*
*     OPFILE   OPLOOK   OPNUMB   OPSCAN   OPTOKN   OPUPPR
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      SUBROUTINE OPFILE( IOPTNS, NOUT, INFORM, OPKEY )
      INTEGER            IOPTNS, NOUT, INFORM
      EXTERNAL           OPKEY
 
************************************************************************
*     OPFILE  reads the options file from unit  IOPTNS  and loads the
*     options into the relevant elements of the integer and real
*     parameter arrays.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version dated December 18, 1985.
************************************************************************
      LOGICAL             PRNT
      CHARACTER*16        KEY   , TOKEN(1)
      CHARACTER*72        BUFFER, OLDBUF
 
      PRNT   = .TRUE.
 
*     Return if the unit number is out of range.
 
      IF (IOPTNS .LT. 0  .OR.  IOPTNS .GT. 99) THEN
         INFORM = 1
         RETURN
      END IF
 
*     ------------------------------------------------------------------
*     Look for  BEGIN, ENDRUN  or  SKIP.
*     ------------------------------------------------------------------
      NREAD  = 0
   50    READ (IOPTNS, '(A)', END = 930) BUFFER
         NREAD = NREAD + 1
         NKEY  = 1
         CALL OPTOKN( BUFFER, NKEY, TOKEN )
         KEY   = TOKEN(1)
         IF (KEY .EQ. 'ENDRUN') GO TO 940
         IF (KEY .NE. 'BEGIN' ) THEN
            IF (NREAD .EQ. 1  .AND.  KEY .NE. 'SKIP') THEN
               WRITE (NOUT, 2000) IOPTNS, BUFFER
            END IF
            GO TO 50
         END IF
 
*     ------------------------------------------------------------------
*     BEGIN found.
*     This is taken to be the first line of an OPTIONS file.
*     Read the second line to see if it is NOLIST.
*     ------------------------------------------------------------------
      OLDBUF = BUFFER
      READ (IOPTNS, '(A)', END = 920) BUFFER
 
      CALL OPKEY ( NOUT, BUFFER, KEY )
 
      IF (KEY .EQ. 'NOLIST') THEN
         PRNT   = .FALSE.
      END IF
 
      IF (PRNT) THEN
         WRITE (NOUT, '(// A / A /)')
     $      ' OPTIONS file',
     $      ' ------------'
         WRITE (NOUT, '(6X, A )') OLDBUF, BUFFER
      END IF
 
*     ------------------------------------------------------------------
*     Read the rest of the file.
*     ------------------------------------------------------------------
*+    while (key .ne. 'end') loop
  100 IF    (KEY .NE. 'END') THEN
         READ (IOPTNS, '(A)', END = 920) BUFFER
         IF (PRNT)
     $      WRITE (NOUT, '( 6X, A )') BUFFER
 
         CALL OPKEY ( NOUT, BUFFER, KEY )
 
         IF (KEY .EQ.   'LIST') PRNT = .TRUE.
         IF (KEY .EQ. 'NOLIST') PRNT = .FALSE.
         GO TO 100
      END IF
*+    end while
 
      INFORM =  0
      RETURN
 
  920 WRITE (NOUT, 2200) IOPTNS
      INFORM = 2
      RETURN
 
  930 WRITE (NOUT, 2300) IOPTNS
      INFORM = 3
      RETURN
 
  940 WRITE (NOUT, '(// 6X, A)') BUFFER
      INFORM = 4
      RETURN
 
 2000 FORMAT(
     $ //' XXX  Error while looking for an OPTIONS file on unit', I7
     $ / ' XXX  The file should start with BEGIN, SKIP or ENDRUN'
     $ / ' XXX  but the first record found was the following:'
     $ //' ---->', A
     $ //' XXX  Continuing to look for OPTIONS file...')
 2200 FORMAT(//' XXX  End-of-file encountered while processing',
     $         ' an OPTIONS file on unit', I6)
 2300 FORMAT(//' XXX  End-of-file encountered while looking for',
     $         ' an OPTIONS file on unit', I6)
 
*     End of  OPFILE.
 
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPLOOK (NDICT, DICTRY, ALPHA, KEY, ENTRY)
C
C
C Description and usage:
C
C       Performs dictionary lookups.  A pointer is returned if a
C    match is found between the input key and the corresponding
C    initial characters of one of the elements of the dictionary.
C    If a "synonym" has been provided for an entry, the search is
C    continued until a match to a primary dictionary entry is found.
C    Cases of no match, or multiple matches, are also provided for.
C
C     Dictionary entries must be left-justified, and may be alphabetized
C    for faster searches.  Secondary entries, if any, are composed of
C    two words separated by one or more characters such as blank, tab,
C    comma, colon, or equal sign which are treated as non-significant
C    by OPSCAN.  The first entry of each such pair serves as a synonym
C    for the second, more fundamental keyword.
C
C       The ordered search stops after the section of the dictionary
C    having the same first letters as the key has been checked, or
C    after a specified number of entries have been examined.  A special
C    dictionary entry, the vertical bar '|', will also terminate the
C    search.  This will speed things up if an appropriate dictionary
C    length parameter cannot be determined.  Both types of search are
C    sequential.  See "Notes" below for some suggestions if efficiency
C    is an issue.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    NDICT               I    I      Number of dictionary entries to be
C                                    examined.
C    DICTRY  NDICT       C    I      Array of dictionary entries,
C                                    left-justified in their fields.
C                                    May be alphabetized for efficiency,
C                                    in which case ALPHA should be
C                                    .TRUE.  Entries with synonyms are
C                                    of the form
C                                    'ENTRY : SYNONYM', where 'SYNONYM'
C                                    is a more fundamental entry in the
C                                    same dictionary.  NOTE: Don't build
C                                    "circular" dictionaries!
C    ALPHA               L    I      Indicates whether the dictionary
C                                    is in alphabetical order, in which
C                                    case the search can be terminated
C                                    sooner.
C    KEY                 C    I/O    String to be compared against the
C                                    dictionary.  Abbreviations are OK
C                                    if they correspond to a unique
C                                    entry in the dictionary.  KEY is
C                                    replaced on termination by its most
C                                    fundamental equivalent dictionary
C                                    entry (uppercase, left-justified)
C                                    if a match was found.
C    ENTRY               I      O    Dictionary pointer.  If > 0, it
C                                    indicates which entry matched KEY.
C                                    In case of trouble, a negative
C                                    value means that a UNIQUE match
C                                    was not found - the absolute value
C                                    of ENTRY points to the second
C                                    dictionary entry that matched KEY.
C                                    Zero means that NO match could be
C                                    found.  ENTRY always refers to the
C                                    last search performed -
C                                    in searching a chain of synonyms,
C                                    a non-positive value will be
C                                    returned if there is any break,
C                                    even if the original input key
C                                    was found.
C
C
C External references:
C
C    Name    Description
C    OPSCAN  Finds first and last significant characters.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C    (2)  We have assumed that the dictionary is not too big.  If
C         many searches are to be done or if the dictionary has more
C         than a dozen or so entries, it may be advantageous to build
C         an index array of pointers to the beginning of the section
C         of the dictionary containing each letter, then pass in the
C         portion of the dictionary beginning with DICTRY (INDEX).
C         (This won't generally work for dictionaries with synonyms.)
C         For very large problems, a completely different approach may
C         be advisable, e.g. a binary search for ordered dictionaries.
C
C    (3)  OPLOOK is case sensitive.  In most applications it will be
C         necessary to use an uppercase dictionary, and to convert the
C         input key to uppercase before calling OPLOOK.  Companion
C         routines OPTOKN and PAIRS, available from the author, already
C         take care of this.
C
C    (4)  The key need not be left-justified.  Any leading (or
C         trailing) characters which are "non-significant" to OPSCAN
C         will be ignored.  These include blanks, horizontal tabs,
C         commas, colons, and equal signs.  See OPSCAN for details.
C
C    (5)  The ASCII collating sequence for character data is assumed.
C         (N.B. This means the numerals precede the alphabet, unlike
C         common practice!)  This should not cause trouble on EBCDIC
C         machines if DICTRY just contains alphabetic keywords.
C         Otherwise it may be necessary to use the FORTRAN lexical
C         library routines to force use of the ASCII sequence.
C
C    (6)  Parameter NUMSIG sets a limit on the length of significant
C         dictionary entries.  Special applications may require that
C         this be increased.  (It is 16 in the present version.)
C
C    (7)  No protection against "circular" dictionaries is provided:
C         don't claim that A is B, and that B is A.  All synonym chains
C         must terminate!  Other potential errors not checked for
C         include duplicate or mis-ordered entries.
C
C    (8)  The handling of ambiguities introduces some ambiguity:
C
C            ALPHA = .TRUE.  A potential problem, when one entry
C                            looks like an abbreviation for another
C                            (eg. does 'A' match 'A' or 'AB'?) was
C                            resolved by dropping out of the search
C                            immediately when an "exact" match is found.
C
C            ALPHA = .FALSE. The programmer must ensure that the above
C                            situation does not arise: each dictionary
C                            entry must be recognizable, at least when
C                            specified to full length.  Otherwise, the
C                            result of a search will depend on the
C                            order of entries.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    24 Feb. 1984  RAK/DAS  Initial design and coding.
C    25 Feb. 1984    RAK    Combined the two searches by suitable
C                           choice of terminator FLAG.
C    28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                           longer update KEY.
C    29 Mar. 1984    RAK    Put back replacement of KEY by its
C                           corresponding entry.
C    21 June 1984    RAK    Corrected bug in error handling for cases
C                           where no match was found.
C    23 Apr. 1985    RAK    Introduced test for exact matches, which
C                           permits use of dictionary entries which
C                           would appear to be ambiguous (for ordered
C                           case).  Return -I to point to the entry
C                           which appeared ambiguous (had been -1).
C                           Repaired loop termination - had to use
C                           equal length strings or risk quitting too
C                           soon when one entry is an abbreviation
C                           for another.  Eliminated HIT, reduced
C                           NUMSIG to 16.
C    15 Nov. 1985    MAS    Loop 20 now tests .LT. FLAG, not .LE. FLAG.
C                           If ALPHA is false, FLAG is now '|', not '{'.
C    26 Jan. 1986    PEG    Declaration of FLAG and TARGET modified to
C                           conform to ANSI-77 standard.
C-----------------------------------------------------------------------
 
 
C     Variable declarations.
C     ----------------------
 
*     IMPLICIT NONE
 
C     Parameters.
 
      INTEGER
     $   NUMSIG
      CHARACTER
     $   BLANK, VBAR
      PARAMETER
     $   (BLANK = ' ', VBAR = '|', NUMSIG = 16)
 
C     Variables.
 
      LOGICAL
     $   ALPHA
      INTEGER
     $   ENTRY, FIRST, I, LAST, LENGTH, MARK, NDICT
*     CHARACTER
*    $   DICTRY (NDICT) * (*), FLAG * (NUMSIG),
*    $   KEY * (*), TARGET * (NUMSIG)
      CHARACTER
     $   DICTRY (NDICT) * (*), FLAG * 16,
     $   KEY * (*), TARGET * 16
 
C     Procedures.
 
      EXTERNAL
     $   OPSCAN
 
 
C     Executable statements.
C     ----------------------
 
      ENTRY = 0
 
C     Isolate the significant portion of the input key (if any).
 
      FIRST = 1
      LAST  = MIN( LEN(KEY), NUMSIG )
      CALL OPSCAN (KEY, FIRST, LAST, MARK)
 
      IF (MARK .GT. 0) THEN
         TARGET = KEY (FIRST:MARK)
 
C        Look up TARGET in the dictionary.
 
   10    CONTINUE
            LENGTH = MARK - FIRST + 1
 
C           Select search strategy by cunning choice of termination test
C           flag.  The vertical bar is just about last in both the
C           ASCII and EBCDIC collating sequences.
 
            IF (ALPHA) THEN
               FLAG = TARGET
            ELSE
               FLAG = VBAR
            END IF
 
 
C           Perform search.
C           ---------------
 
            I = 0
   20       CONTINUE
               I = I + 1
               IF (TARGET (1:LENGTH) .EQ. DICTRY (I) (1:LENGTH)) THEN
                  IF (ENTRY .EQ. 0) THEN
 
C                    First "hit" - must still guard against ambiguities
C                    by searching until we've gone beyond the key
C                    (ordered dictionary) or until the end-of-dictionary
C                    mark is reached (exhaustive search).
 
                     ENTRY = I
 
C                    Special handling if match is exact - terminate
C                    search.  We thus avoid confusion if one dictionary
C                    entry looks like an abbreviation of another.
C                    This fix won't generally work for un-ordered
C                    dictionaries!
 
                     FIRST = 1
                     LAST = NUMSIG
                     CALL OPSCAN (DICTRY (ENTRY), FIRST, LAST, MARK)
                     IF (MARK .EQ. LENGTH) I = NDICT
                  ELSE
 
 
C                    Oops - two hits!  Abnormal termination.
C                    ---------------------------------------
 
                     ENTRY = -I
                     RETURN
                  END IF
               END IF
 
C           Check whether we've gone past the appropriate section of the
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.
 
            IF (DICTRY (I) (1:LENGTH) .LT. FLAG  .AND.  I .LT. NDICT)
     $         GO TO 20
 
 
C           Check for a synonym.
C           --------------------
 
            IF (ENTRY .GT. 0) THEN
 
C              Look for a second entry "behind" the first entry.  FIRST
C              and MARK were determined above when the hit was detected.
 
               FIRST = MARK + 2
               CALL OPSCAN (DICTRY (ENTRY), FIRST, LAST, MARK)
               IF (MARK .GT. 0) THEN
 
C                 Re-set target and dictionary pointer, then repeat the
C                 search for the synonym instead of the original key.
 
                  TARGET = DICTRY (ENTRY) (FIRST:MARK)
                  ENTRY = 0
                  GO TO 10
 
               END IF
            END IF
 
      END IF
      IF (ENTRY .GT. 0) KEY = DICTRY (ENTRY)
 
 
C     Normal termination.
C     -------------------
 
      RETURN
 
C     End of OPLOOK
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      FUNCTION OPNUMB( STRING )
 
      LOGICAL          OPNUMB
      CHARACTER*(*)    STRING
 
************************************************************************
*     Description and usage:
*
*        A simple(-minded) test for numeric data is implemented by
*        searching an input string for legitimate characters:
*                digits 0 to 9, D, E, -, + and .
*        Insurance is provided by requiring that a numeric string
*        have at least one digit, at most one D, E or .
*        and at most two -s or +s.  Note that a few ambiguities remain:
*
*           (a)  A string might have the form of numeric data but be
*                intended as text.  No general test can hope to detect
*                such cases.
*
*           (b)  There is no check for correctness of the data format.
*                For example a meaningless string such as 'E1.+2-'
*                will be accepted as numeric.
*
*        Despite these weaknesses, the method should work in the
*        majority of cases.
*
*
*     Parameters:
*
*        Name    Dimension  Type  I/O/S  Description
*        OPNUMB              L      O    Set .TRUE. if STRING appears
*                                        to be numerical data.
*        STRING              C    I      Input data to be tested.
*
*
*     Environment:  ANSI FORTRAN 77.
*
*
*     Notes:
*
*        (1)  It is assumed that STRING is a token extracted by
*             OPTOKN, which will have converted any lower-case
*             characters to upper-case.
*
*        (2)  OPTOKN pads STRING with blanks, so that a genuine
*             number is of the form  '1234        '.
*             Hence, the scan of STRING stops at the first blank.
*
*        (3)  COMPLEX data with parentheses will not look numeric.
*
*
*     Systems Optimization Laboratory, Stanford University.
*     12 Nov  1985    Initial design and coding, starting from the
*                     routine ALPHA from Informatics General, Inc.
************************************************************************
 
      LOGICAL         NUMBER
      INTEGER         J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS, NPOINT
      CHARACTER*1     ATOM
 
      NDIGIT = 0
      NEXP   = 0
      NMINUS = 0
      NPLUS  = 0
      NPOINT = 0
      NUMBER = .TRUE.
      LENGTH = LEN (STRING)
      J      = 0
 
   10    J    = J + 1
         ATOM = STRING (J:J)
         IF      (ATOM .GE. '0'  .AND.  ATOM .LE. '9') THEN
            NDIGIT = NDIGIT + 1
         ELSE IF (ATOM .EQ. 'D'  .OR.   ATOM .EQ. 'E') THEN
            NEXP   = NEXP   + 1
         ELSE IF (ATOM .EQ. '-') THEN
            NMINUS = NMINUS + 1
         ELSE IF (ATOM .EQ. '+') THEN
            NPLUS  = NPLUS  + 1
         ELSE IF (ATOM .EQ. '.') THEN
            NPOINT = NPOINT + 1
         ELSE IF (ATOM .EQ. ' ') THEN
            J      = LENGTH
         ELSE
            NUMBER = .FALSE.
         END IF
 
         IF (NUMBER  .AND.  J .LT. LENGTH) GO TO 10
 
      OPNUMB = NUMBER
     $         .AND.  NDIGIT .GE. 1
     $         .AND.  NEXP   .LE. 1
     $         .AND.  NMINUS .LE. 2
     $         .AND.  NPLUS  .LE. 2
     $         .AND.  NPOINT .LE. 1
 
      RETURN
 
*     End of OPNUMB
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPSCAN (STRING, FIRST, LAST, MARK)
C
C
C Description and usage:
C
C       Looks for non-blank fields ("tokens") in a string, where the
C    fields are of arbitrary length, separated by blanks, tabs, commas,
C    colons, or equal signs.  The position of the end of the 1st token
C    is also returned, so this routine may be conveniently used within
C    a loop to process an entire line of text.
C
C       The procedure examines a substring, STRING (FIRST : LAST), which
C    may of course be the entire string (in which case just call OPSCAN
C    with FIRST <= 1 and LAST >= LEN (STRING) ).  The indices returned
C    are relative to STRING itself, not the substring.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    STRING              C    I      Text string containing data to be
C                                    scanned.
C    FIRST               I    I/O    Index of beginning of substring.
C                                    If <= 1, the search begins with 1.
C                                    Output is index of beginning of
C                                    first non-blank field, or 0 if no
C                                    token was found.
C    LAST                I    I/O    Index of end of substring.
C                                    If >= LEN (STRING), the search
C                                    begins with LEN (STRING).  Output
C                                    is index of end of last non-blank
C                                    field, or 0 if no token was found.
C    MARK                I      O    Points to end of first non-blank
C                                    field in the specified substring.
C                                    Set to 0 if no token was found.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               ANSI Fortran 77, except for the tab character HT.
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C         in a non-standard way:  the CHAR function is not permitted
C         in a PARAMETER declaration (OK on VAX, though).  For Absoft
C         FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it
C         may be best to declare HT as a variable and assign
C         HT = CHAR(9) on ASCII machines, or CHAR(5) for EBCDIC.
C
C    (2)  The pseudo-recursive structure was chosen for fun.  It is
C         equivalent to three DO loops with embedded GO TOs in sequence.
C
C    (3)  The variety of separators recognized limits the usefulness of
C         this routine somewhat.  The intent is to facilitate handling
C         such tokens as keywords or numerical values.  In other
C         applications, it may be necessary for ALL printing characters
C         to be significant.  A simple modification to statement
C         function SOLID will do the trick.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                           based on SCAN_STRING by Ralph Carmichael.
C    25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C    16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                           (previous re-use of STRING was ambiguous).
C
C-----------------------------------------------------------------------
 
 
C     Variable declarations.
C     ----------------------
 
*     IMPLICIT NONE
 
C     Parameters.
 
      CHARACTER
     $   BLANK, EQUAL, COLON, COMMA, HT
      PARAMETER
     $   (BLANK = ' ', EQUAL = '=', COLON = ':', COMMA = ',')
 
C     Variables.
 
      LOGICAL
     $   SOLID
      INTEGER
     $   BEGIN, END, FIRST, LAST, LENGTH, MARK
      CHARACTER
     $   DUMMY, STRING * (*)
 
C     Statement functions.
 
      SOLID (DUMMY) = (DUMMY .NE. BLANK) .AND.
     $                (DUMMY .NE. COLON) .AND.
     $                (DUMMY .NE. COMMA) .AND.
     $                (DUMMY .NE. EQUAL) .AND.
     $                (DUMMY .NE. HT)
 
 
C     Executable statements.
C     ----------------------
 
****  HT     = CHAR(9) for ASCII machines, CHAR(5) for EBCDIC.
      HT     = CHAR(9)
      MARK   = 0
      LENGTH = LEN (STRING)
      BEGIN  = MAX (FIRST, 1)
      END    = MIN (LENGTH, LAST)
 
C     Find the first significant character ...
 
      DO 30 FIRST = BEGIN, END, +1
         IF (SOLID (STRING (FIRST : FIRST))) THEN
 
C           ... then the end of the first token ...
 
            DO 20 MARK = FIRST, END - 1, +1
               IF (.NOT.SOLID (STRING (MARK + 1 : MARK + 1))) THEN
 
C                 ... and finally the last significant character.
 
                  DO 10 LAST = END, MARK, -1
                     IF (SOLID (STRING (LAST : LAST))) THEN
                        RETURN
                     END IF
   10             CONTINUE
 
C                 Everything past the first token was a separator.
 
                  LAST = LAST + 1
                  RETURN
               END IF
   20       CONTINUE
 
C           There was nothing past the first token.
 
            LAST = MARK
            RETURN
         END IF
   30 CONTINUE
 
C     Whoops - the entire substring STRING (BEGIN : END) was composed of
C     separators !
 
      FIRST = 0
      MARK = 0
      LAST = 0
      RETURN
 
C     End of OPSCAN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPTOKN (STRING, NUMBER, LIST)
C
C
C Description and usage:
C
C       An aid to parsing input data.  The individual "tokens" in a
C    character string are isolated, converted to uppercase, and stored
C    in an array.  Here, a token is a group of significant, contiguous
C    characters.  The following are NON-significant, and hence may
C    serve as separators:  blanks, horizontal tabs, commas, colons,
C    and equal signs.  See OPSCAN for details.  Processing continues
C    until the requested number of tokens have been found or the end
C    of the input string is reached.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    STRING              C    I      Input string to be analyzed.
C    NUMBER              I    I/O    Number of tokens requested (input)
C                                    and found (output).
C    LIST    NUMBER      C      O    Array of tokens, changed to upper
C                                    case.
C
C
C External references:
C
C    Name    Description
C    OPSCAN  Finds positions of first and last significant characters.
C    OPUPPR  Converts a string to uppercase.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    16 Jan. 1984    RAK    Initial design and coding.
C    16 Mar. 1984    RAK    Revised header to reflect full list of
C                           separators, repaired faulty WHILE clause
C                           in "10" loop.
C    18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                           at a time, leaving STRING unchanged.
C
C-----------------------------------------------------------------------
 
 
C     Variable declarations.
C     ----------------------
 
*     IMPLICIT NONE
 
C     Parameters.
 
      CHARACTER
     $   BLANK
      PARAMETER
     $   (BLANK = ' ')
 
C     Variables.
 
      INTEGER
     $   COUNT, FIRST, I, LAST, MARK, NUMBER
      CHARACTER
     $   STRING * (*), LIST (NUMBER) * (*)
 
C     Procedures.
 
      EXTERNAL
     $   OPUPPR, OPSCAN
 
 
C     Executable statements.
C     ----------------------
 
C     WHILE there are tokens to find, loop UNTIL enough have been found.
 
      FIRST = 1
      LAST = LEN (STRING)
 
      COUNT = 0
   10 CONTINUE
 
C        Get delimiting indices of next token, if any.
 
         CALL OPSCAN (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN
            COUNT = COUNT + 1
 
C           Pass token to output string array, then change case.
 
            LIST (COUNT) = STRING (FIRST : MARK)
            CALL OPUPPR (LIST (COUNT))
            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 10
 
         END IF
 
 
C     Fill the rest of LIST with blanks and set NUMBER for output.
 
      DO 20 I = COUNT + 1, NUMBER
         LIST (I) = BLANK
   20 CONTINUE
 
      NUMBER = COUNT
 
 
C     Termination.
C     ------------
 
      RETURN
 
C     End of OPTOKN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPUPPR(STRING)
C
C ACRONYM:  UPper CASE
C
C PURPOSE:  This subroutine changes all lower case letters in the
C           character string to upper case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           be different on ASCII and EBCDIC machines.
C
C ARGUMENTS
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C  STRING       *       C   I/O   Character string possibly containing
C                                 some lower-case letters on input;
C                                 strictly upper-case letters on output
C                                 with no change to any non-alphabetic
C                                 characters.
C
C EXTERNAL REFERENCES:
C  LEN    - Returns the declared length of a CHARACTER variable.
C  INDEX  - Returns the position of second string within first.
C
C ENVIRONMENT:  ANSI FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   06/28/83   CLH    Initial design.
C   01/03/84   RAK    Eliminated NCHAR input.
C   06/14/84   RAK    Used integer PARAMETERs in comparison.
C   04/21/85   RAK    Eliminated DO/END DO in favor of standard code.
C   09/10/85   MAS    Eliminated CHAR,ICHAR in favor of LOW, UPP, INDEX.
C
C AUTHOR: Charles Hooper, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------
 
      CHARACTER      STRING * (*)
      INTEGER        I, J
      CHARACTER      C*1, LOW*26, UPP*26
      DATA           LOW /'abcdefghijklmnopqrstuvwxyz'/,
     $               UPP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
 
      DO 10 J = 1, LEN(STRING)
         C    = STRING(J:J)
         IF (C .GE. 'a'  .AND.  C .LE. 'z') THEN
            I           = INDEX( LOW, C )
            IF (I .GT. 0) STRING(J:J) = UPP(I:I)
         END IF
   10 CONTINUE
      RETURN
 
*     End of OPUPPR
      END
