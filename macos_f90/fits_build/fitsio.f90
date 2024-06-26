!
! main source code for writing/reading FITS (Flexible Image Transport System) files
! downloaded from http://orion.math.iastate.edu/burkardt/g_src/fitsio/fitsio.html
!

subroutine ftadef(ounit,lenrow,nfield,bcol,tform,nrows,status)
!
!*******************************************************************************
!
!! FTADEF defines the structure of the ASCII table data unit.
!
!       Ascii table data DEFinition
!       define the structure of the ASCII table data unit
!
!       ounit   i  Fortran I/O unit number
!       lenrow  i  length of a row, in characters
!       nfield  i  number of fields in the table
!       bcol    i  starting position of each column, (starting with 1)
!       tform   C  the data format of the column
!       nrows   i  number of rows in the table
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,lenrow,nfield,bcol(*),nrows,status
    character ( len = * ) tform(*)

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer ibuff,i,j,clen,c2
    character ctemp*24, cnum*3,cbcol*10,caxis1*10

    if (status > 0)return

    ibuff=bufnum(ounit)

    if (dtstrt(ibuff) < 0)then
!               freeze the header at its current size
            call fthdef(ounit,0,status)
            if (status > 0)return
    end if

    hdutyp(ibuff)=1
    tfield(ibuff)=nfield

    if (nxtfld + nfield > nf)then
!               too many columns open at one time; exceeded array dimensions
            status=111
            return
    end if

    tstart(ibuff)=nxtfld
    nxtfld=nxtfld+nfield

    if (nfield == 0)then
!           no data; the next HDU begins in the next logical block
        hdstrt(ibuff,chdu(ibuff)+1)=dtstrt(ibuff)
        heapsz(ibuff)=0
        theap(ibuff)=0
    else
!           initialize the table column parameters
        clen=len(tform(1))
        do 20 i=1,nfield
            tscale(i+tstart(ibuff))=1.
            tzero(i+tstart(ibuff))=0.
!               choose special value to indicate null values are not defined
            cnull(i+tstart(ibuff))=char(1)
            cform(i+tstart(ibuff))=tform(i)
            tbcol(i+tstart(ibuff))=bcol(i)-1
            tdtype(i+tstart(ibuff))=16
!               the repeat count is always one for ASCII tables
            trept(i+tstart(ibuff))=1
!               store the width of the field in TNULL
            c2=0
            do j=2,clen
                    if (tform(i)(j:j) >= '0' .and. &
                       tform(i)(j:j) <= '9')then
                            c2=j
                    else
                            exit
                    end if
            end do

            if (c2 == 0)then
!                       no explicit width, so assume width of 1 character
                    tnull(i+tstart(ibuff))=1
            else
                call ftc2ii(tform(i)(2:c2),tnull(i+tstart(ibuff)) &
                            ,status)
                if (status > 0)then
!                        error parsing TFORM to determine field width
                     status=261
                     ctemp=tform(i)
                     call ftpmsg('Error parsing TFORM to get field' &
                      //' width: '//ctemp)
                     return
                end if
            end if

!               check that column fits within the table
            if (tbcol(i+tstart(ibuff))+tnull(i+tstart(ibuff)) &
              > lenrow .and. lenrow /= 0)then
                status=236
                write(cnum,1000)i
                write(cbcol,1001)bcol(i)
                write(caxis1,1001)lenrow
1000                format(i3)
1001                format(i10)
                call ftpmsg('Column '//cnum//' will not fit '// &
               'within the specified width of the ASCII table.')

                call ftpmsg('TFORM='//cform(i+tstart(ibuff))// &
                '  TBCOL='//cbcol//'  NAXIS1='//caxis1)
                return
             end if
20           continue

!           calculate the start of the next header unit, based on the
!           size of the data unit
        rowlen(ibuff)=lenrow
        hdstrt(ibuff,chdu(ibuff)+1)= &
            dtstrt(ibuff)+(lenrow*nrows+2879)/2880*2880
!
!  initialize the fictitious heap starting address (immediately following
!  the table data) and a zero length heap.  This is used to find the
!  end of the table data when checking the fill values in the last block.
!  ASCII tables have no special data area.
!
        heapsz(ibuff)=0
        theap(ibuff)=rowlen(ibuff)*nrows
    end if
end
subroutine ftaini(iunit,status)
!
!*******************************************************************************
!
!! FTAINI initializes the parameters defining the structure of an ASCII table.
!
!       iunit   i  Fortran I/O unit number
!       OUTPUT PARAMETERS:
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer nrows,tfld,nkey,ibuff,i,nblank
    character keynam*8,value*70,comm*72,rec*80
    character cnum*3,cbcol*10,caxis1*10

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       store the type of HDU (1 = ASCII table extension)
    hdutyp(ibuff)=1

!       temporarily set the location of the end of the header to a huge number
    hdend(ibuff)=2000000000
    hdstrt(ibuff,chdu(ibuff)+1)=2000000000

!       check that this is a valid ASCII table, and get parameters
    call ftgttb(iunit,rowlen(ibuff),nrows,tfld,status)
    if (status > 0)go to 900

    if (tfld > nf)then
!               arrays not dimensioned large enough for this many fields
            status=111
    call ftpmsg('This ASCII table has too many fields '// &
   'to be read with FITSIO (FTAINI).')
            go to 900
     end if

!       store the number of fields in the common block
    tfield(ibuff)=tfld

    if (nxtfld + tfld > nf)then
!               too many columns open at one time; exceeded array dimensions
            status=111
            return
    end if

    tstart(ibuff)=nxtfld
    nxtfld=nxtfld+tfld

!       initialize the table field parameters
    do i=1,tfld
            tscale(i+tstart(ibuff))=1.
            tzero(i+tstart(ibuff))=0.
!               choose special value to indicate that null value is not defined
            cnull(i+tstart(ibuff))=char(1)
!               pre-set required keyword values to a null value
            tbcol(i+tstart(ibuff))=-1
            tdtype(i+tstart(ibuff))=-9999
    end do

!       initialize the fictitious heap starting address (immediately following
!       the table data) and a zero length heap.  This is used to find the
!   end of the table data when checking the fill values in the last block.
!       there is no special data following an ASCII table
    heapsz(ibuff)=0
    theap(ibuff)=rowlen(ibuff)*nrows

!       now read through the rest of the header looking for table column
!       definition keywords, and the END keyword.

    nkey=8
8       nblank=0
10      nkey=nkey+1
    call ftgrec(iunit,nkey,rec,status)
    if (status == 107)then
!               if we hit the end of file, then set status = no END card found
            status=210
    call ftpmsg('Required END keyword not found in ASCII table'// &
    ' header (FTAINI).')
            go to 900
    else if (status > 0)then
            go to 900
    end if
    keynam=rec(1:8)
    comm=rec(9:80)

    if (keynam(1:1) == 'T')then
!               get the ASCII table parameter (if it is one)
            call ftpsvc(rec,value,comm,status)
            call ftgatp(ibuff,keynam,value,status)
    else if (keynam == ' ' .and. comm == ' ')then
            nblank=nblank+1
            go to 10
    else if (keynam == 'END')then
            go to 20
    end if
    go to 8

20      continue

!       test that all the required keywords were found
    do 25 i=1,tfld
        if (tbcol(i+tstart(ibuff)) == -1)then
            status=231
            call ftkeyn('TBCOL',i,keynam,status)
            call ftpmsg('Required '//keynam// &
                        ' keyword not found (FTAINI).')
            return
        else if (tbcol(i+tstart(ibuff)) < 0 .or. &
                 tbcol(i+tstart(ibuff)) >= rowlen(ibuff) &
                 .and. rowlen(ibuff) /= 0)then
            status=234
            call ftkeyn('TBCOL',i,keynam,status)
            call ftpmsg('Value of the '//keynam// &
                        ' keyword is out of range (FTAINI).')
            return

!               check that column fits within the table
        else if (tbcol(i+tstart(ibuff))+tnull(i+tstart(ibuff)) > &
                 rowlen(ibuff) .and. rowlen(ibuff) /= 0)then
                status=236
                write(cnum,1000)i
                write(cbcol,1001)tbcol(i+tstart(ibuff))+1
                write(caxis1,1001)rowlen(ibuff)
1000                format(i3)
1001                format(i10)
                call ftpmsg('Column '//cnum//' will not fit '// &
               'within the specified width of the ASCII table.')

                call ftpmsg('TFORM='//cform(i+tstart(ibuff))// &
                '  TBCOL='//cbcol//'  NAXIS1='//caxis1)
                return
        else if (tdtype(i+tstart(ibuff)) == -9999)then
            status=232
            call ftkeyn('TFORM',i,keynam,status)
            call ftpmsg('Required '//keynam// &
                        ' keyword not found (FTAINI).')
            return
        end if
25      continue

!       now we know everything about the table; just fill in the parameters:
!       the 'END' record begins 80 bytes before the current position,
!       ignoring any trailing blank keywords just before the END keyword
    hdend(ibuff)=nxthdr(ibuff)-80*(nblank+1)

!       the data unit begins at the beginning of the next logical block
    dtstrt(ibuff)=((nxthdr(ibuff)-80)/2880+1)*2880

!       reset header pointer to the first keyword
    nxthdr(ibuff)=hdstrt(ibuff,chdu(ibuff))

!       the next HDU begins in the next logical block after the data
    hdstrt(ibuff,chdu(ibuff)+1)= &
    dtstrt(ibuff)+(rowlen(ibuff)*nrows+2879)/2880*2880

900     continue
end
subroutine ftarch(iword,jword,compid)
!
!*******************************************************************************
!
!! FTARCH figures out what kind of machine it is running on.
!
!              compid =  0  -  Big Endian (SUN, Mac, Next, SGI)
!                        1  -  Little Endian (Dec Ultrix, OSF/1, PC)
!                        2  -  Vax VMS
!                        3  -  Alpha VMS
!                        4  -  IBM mainframe
!                       -1  -  SUN F() compiler (maps I*2 variables into I*4)
!       (large neg number)  -  Cray supercomputer

    integer compid
    integer*2 iword(2)
    integer jword(2)

!       Look at the equivalent integer, to distinquish the machine type.
!       The machine type is needed when testing for NaNs.

    if (iword(1) == 16270)then
!               looks like a SUN workstation (uses IEEE word format)
            compid=0
    else if (iword(1) == 14564)then
!               looks like a Decstation, alpha OSF/1, or IBM PC (byte swapped)
            compid=1
    else if (iword(1) == 16526)then
            if (jword(1) == 954417294)then
!                   looks like a VAX VMS system
                compid=2
            else
!                   looks like ALPHA VMS system
                compid=3
            end if
    else if (iword(1) == 16657)then
!               an IBM main frame (the test for NaNs is the same as on SUNs)
            compid=4
    else if (iword(1) == 1066285284)then
!               SUN F90 compiler maps I*2 variables into I*4
            compid= (-1)
    else
!               unknown machine
            compid=0
    end if
end
subroutine ftasfm(form,dtype,width,decims,status)
!
!*******************************************************************************
!
!! FTASFM parses the ASCII table column format.
!
!  The routine parses the ASCII table TFORM column format to determine the data
!       type, the field width, and number of decimal places (if relevant)
!
!       form    c  TFORM format string
!       OUTPUT PARAMETERS:
!       dattyp  i  datatype code
!       width   i  width of the field
!       decims  i  number of decimal places
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, November 1994

    character ( len = * ) form
    integer dtype,width,decims,status
    character dattyp*1,cform*16
    integer nc,c1,i,nw

    if (status > 0)return

    cform=form

!       find first non-blank character
    nc=len(form)
    do i=1,nc
      if (form(i:i) /= ' ')then
        c1=i
        go to 10
      end if
    end do

!       error: TFORM is a blank string
    status=261
    call ftpmsg('The TFORM keyword has a blank value.')
    return

10      continue

!       now the chararcter at position c1 should be the data type code
    dattyp=form(c1:c1)

!       set the numeric datatype code
    if (dattyp == 'I')then
            dtype=41
    else if (dattyp == 'E')then
            dtype=42
    else if (dattyp == 'F')then
            dtype=42
    else if (dattyp == 'D')then
            dtype=82
    else if (dattyp == 'A')then
            dtype=16
    else
!               unknown tform datatype code
            status=262
            call ftpmsg('Unknown ASCII table TFORMn keyword '// &
                        'datatype: '//cform)
            return
    end if

!       determine the field width
    c1=c1+1
    nw=0

    do i=c1,nc
      if (form(i:i) >= '0' .and. form(i:i)<='9')then
        nw=nw+1
      else
        exit
      end if
    end do

    if (nw == 0)then
!               error, no width specified
            go to 990
    else
            call ftc2ii(form(c1:c1+nw-1),width,status)
            if (status > 0 .or. width == 0)then
!                      unrecognized characters following the type code
                   go to 990
            end if
    end if

!       determine the number of decimal places (if any)
    decims=-1
    c1=c1+nw
    if (form(c1:c1) == '.')then
        c1=c1+1
        nw=0
        do 60 i=c1,nc
            if (form(i:i) >= '0' .and. form(i:i)<='9')then
                nw=nw+1
            else
                go to 70
        end if
60          continue
70          continue

        if (nw == 0)then
!               error, no decimals specified
            go to 990
        else
            call ftc2ii(form(c1:c1+nw-1),decims,status)
            if (status > 0)then
!                   unrecognized characters
                go to 990
            end if
        end if
    else if (form(c1:c1) /= ' ')then
        go to 990
    end if

!       consistency checks
    if (dattyp == 'A' .or. dattyp == 'I')then
        if (decims == -1)then
            decims=0
        else
            go to 990
        end if
    else if (decims == -1)then
!           number of decmal places must be specified for D, E, or F fields
        go to 990
    else if (decims >= width)then
!           number of decimals must be less than the width
        go to 990
    end if

    if (dattyp == 'I')then
!           set datatype to SHORT integer if 4 digits or less
        if (width <= 4)dtype=21
    else if (dattyp == 'F')then
!           set datatype to DOUBLE if 8 digits or more
        if (width >= 8)dtype=82
    end if

    return

990     continue
    status=261
    call ftpmsg('Illegal ASCII table TFORMn keyword: '//cform)
end
subroutine ftbdef(ounit,nfield,tform,pcount,nrows,status)
!
!*******************************************************************************
!
!! FTBDEF defines the structure of the binary table data unit.
!
!       Binary table data DEFinition
!       define the structure of the binary table data unit
!
!       ounit   i  Fortran I/O unit number
!       nfield  i  number of fields in the table
!       tform   C  the data format of the column
!       nrows   i  number of rows in the table
!       pcount  i  size in bytes of the special data block following the table
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,nfield,nrows,pcount,status
    character ( len = * ) tform(*)

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer ibuff,i,j,width

    if (status > 0)return

    ibuff=bufnum(ounit)

    if (dtstrt(ibuff) < 0)then
!               freeze the header at its current size
            call fthdef(ounit,0,status)
            if (status > 0)return
    end if

    hdutyp(ibuff)=2
    tfield(ibuff)=nfield

    if (nxtfld + nfield > nf)then
!               too many columns open at one time; exceeded array dimensions
            status=111
            return
    end if

    tstart(ibuff)=nxtfld
    nxtfld=nxtfld+nfield

    if (nfield == 0)then
!           no data; the next HDU begins in the next logical block
        hdstrt(ibuff,chdu(ibuff)+1)=dtstrt(ibuff)
        heapsz(ibuff)=0
        theap(ibuff)=0
    else
!           initialize the table column parameters
        do 5 i=1,nfield
            tscale(i+tstart(ibuff))=1.
            tzero(i+tstart(ibuff))=0.
!               choose special value to indicate that null value is not defined
            tnull(i+tstart(ibuff))=123454321
!               reset character NUL string, in case it has been
!               previously defined from an ASCII table extension
            cnull(i+tstart(ibuff))=char(0)

!               parse the tform strings to get the data type and repeat count
            call ftbnfm(tform(i),tdtype(i+tstart(ibuff)), &
                        trept(i+tstart(ibuff)),width,status)
            if (tdtype(i+tstart(ibuff)) == 1)then
!                  treat Bit datatype as if it were a Byte datatype
               tdtype(i+tstart(ibuff))=11
               trept(i+tstart(ibuff))=(trept(i+tstart(ibuff))+7)/8
            else if (tdtype(i+tstart(ibuff)) == 16)then
!                       store ASCII unit string length in TNULL parameter
                    tnull(i+tstart(ibuff))=width
            end if
            if (status > 0)return
5           continue

!           determine byte offset of the beginning of each field and row length
        call ftgtbc(nfield,tdtype(1+tstart(ibuff)),trept(1+ &
             tstart(ibuff)),tbcol(1+tstart(ibuff)),rowlen(ibuff), &
                    status)

!           FITSIO deals with ASCII columns as arrays of strings, not
!           arrays of characters, so need to change the repeat count
!           to indicate the number of strings in the field, not the
!           total number of characters in the field.
!
        do i=1,nfield
            if (tdtype(i+tstart(ibuff)) == 16)then
                j=trept(i+tstart(ibuff))/tnull(i+tstart(ibuff))
                trept(i+tstart(ibuff))=max(j,1)
            end if
        end do

!           initialize the heap offset (=nrows x ncolumns)
!           set initial size of the special data area = 0;
!           update keyword with the correct final value when the HDU is closed
        heapsz(ibuff)=0
        theap(ibuff)=nrows*rowlen(ibuff)

!           calculate the start of the next header unit, based on the
!           size of the data unit (table + special data)
        hdstrt(ibuff,chdu(ibuff)+1)= &
         dtstrt(ibuff)+(rowlen(ibuff)*nrows+pcount+2879)/2880*2880
    end if
end
subroutine ftbini(iunit,status)
!
!*******************************************************************************
!
!! FTBINI initializes the parameters defining the structure of a binary table.
!
!       iunit   i  Fortran I/O unit number
!       OUTPUT PARAMETERS:
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer lenrow,nrows,pcnt,tfld,nkey,ibuff,i,j,nblank
    character keynam*8,value*70,comm*72,cnaxis*8,clen*8,rec*80
    character nulchr*16

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       store the type of HDU (2 = Binary table extension)
    hdutyp(ibuff)=2

!       temporarily set the location of the end of the header to a huge number
    hdend(ibuff)=2000000000
    hdstrt(ibuff,chdu(ibuff)+1)=2000000000

!       check that this is a valid binary table, and get parameters
    call ftgtbn(iunit,rowlen(ibuff),nrows,pcnt,tfld,status)
    if (status > 0)go to 900

    if (tfld > nf)then
!               arrays not dimensioned large enough for this many fields
            status=111
    call ftpmsg('This Binary table has too many fields '// &
   'to be read with FITSIO (FTBINI).')
            go to 900
     end if

!       store the number of fields in the common block
    tfield(ibuff)=tfld

    if (nxtfld + tfld > nf)then
!               too many columns open at one time; exceeded array dimensions
            status=111
            return
    end if

    tstart(ibuff)=nxtfld
    nxtfld=nxtfld+tfld

    do i=1,16
            nulchr(i:i) = char(0)
    end do

!       initialize the table field parameters
    do 5 i=1,tfld
            tscale(i+tstart(ibuff))=1.
            tzero(i+tstart(ibuff))=0.
            tnull(i+tstart(ibuff))=123454321
            tdtype(i+tstart(ibuff))=-9999
            trept(i+tstart(ibuff))=0
!               reset character NUL string, in case it has been previously
!               defined from an ASCII table extension
            cnull(i+tstart(ibuff))=nulchr
5       continue

!       initialize the default heap starting address (immediately following
!       the table data) and set the next empty heap address
!       PCOUNT specifies the amount of special data following the table
    heapsz(ibuff)=pcnt
    theap(ibuff)=rowlen(ibuff)*nrows

!       now read through the rest of the header looking for table column
!       definition keywords, and the END keyword.

    nkey=8
8       nblank=0
10      nkey=nkey+1
    call ftgrec(iunit,nkey,rec,status)
    if (status == 107)then
!               if we hit the end of file, then set status = no END card found
            status=210
    call ftpmsg('Required END keyword not found in Binary table'// &
    ' header (FTBINI).')
            go to 900
    else if (status > 0)then
            go to 900
    end if
    keynam=rec(1:8)
    comm=rec(9:80)

    if (keynam(1:1) == 'T')then
!               get the binary table parameter (if it is one)
            call ftpsvc(rec,value,comm,status)
            call ftgbtp(ibuff,keynam,value,status)
    else if (keynam == ' ' .and. comm == ' ')then
            nblank=nblank+1
            go to 10
    else if (keynam == 'END')then
            go to 20
    end if
    go to 8

20      continue

!       test that all the required keywords were found
    do 25 i=1,tfld
        if (tdtype(i+tstart(ibuff)) == -9999)then
            status=232
            call ftkeyn('TFORM',i,keynam,status)
            call ftpmsg('Required '//keynam// &
                        ' keyword not found (FTAINI).')
            return
        end if
25      continue

!       now we know everything about the table; just fill in the parameters:
!       the 'END' record begins 80 bytes before the current position, ignoring
!       any trailing blank keywords just before the END keyword
    hdend(ibuff)=nxthdr(ibuff)-80*(nblank+1)

!       the data unit begins at the beginning of the next logical block
    dtstrt(ibuff)=((nxthdr(ibuff)-80)/2880+1)*2880

!       reset header pointer to the first keyword
    nxthdr(ibuff)=hdstrt(ibuff,chdu(ibuff))

!       the next HDU begins in the next logical block after the data
    hdstrt(ibuff,chdu(ibuff)+1)= &
    dtstrt(ibuff)+(rowlen(ibuff)*nrows+pcnt+2879)/2880*2880

!       determine the byte offset of the beginning of each field and row length
    if (tfld > 0)then
       call ftgtbc(tfld,tdtype(1+tstart(ibuff)), &
       trept(1+tstart(ibuff)),tbcol(1+tstart(ibuff)),lenrow,status)

!          FITSIO deals with ASCII columns as arrays of strings, not
!          arrays of characters, so need to change the repeat count
!          to indicate the number of strings in the field, not the
!          total number of characters in the field.
       do 30 i=1,tfld
          if (tdtype(i+tstart(ibuff)) == 16)then
!                avoid 'divide by zero' in case TFORMn = '0A'
             if (tnull(i+tstart(ibuff)) /= 0)then
                j=trept(i+tstart(ibuff))/tnull(i+tstart(ibuff))
                trept(i+tstart(ibuff))=max(j,1)
             end if
          end if
30         continue
       if (status > 0)go to 900

!          check that the sum of the column widths = NAXIS2 value
       if (rowlen(ibuff) /= lenrow)then
            status=241
            write(cnaxis,1001)rowlen(ibuff)
            write(clen,1001)lenrow
1001            format(i8)
       call ftpmsg('NAXIS1 ='//cnaxis//' not equal'// &
       ' to the sum of the column widths ='//clen//' (FTBINI).')
       end if
    end if

900     continue
end
subroutine ftbnfm(form,dtype,rcount,width,status)
!
!*******************************************************************************
!
!! FTBNFM parses a binary table column format.
!
!
!  the routine parses the binary table column format to determine the data
!       type and the repeat count (and string width, if it is an ASCII field)
!
!       form    c  format string
!       OUTPUT PARAMETERS:
!       dattyp  i  datatype code
!       rcount  i  repeat count
!       width   i  if ASCII field, this is the width of the unit string
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) form
    integer dtype,rcount,width,status,tstat
    character dattyp
    character ( len = 16 ) cform
    integer point,nc,c1,i,nw

    if (status > 0)return

    cform=form
!
!  find first non-blank character
!
    nc=len(form)
    do 5 i=1,nc
            if (form(i:i) /= ' ')then
                    c1=i
                    go to 10
            end if
5       continue

!       error: TFORM is a blank string
    status=261
    call ftpmsg('The TFORM keyword has a blank value.')
    return

10      continue

!       find the size of the field repeat count, if present
    nw=0
    do 20 i=c1,nc
            if (form(i:i) >= '0' .and. form(i:i) <= '9')then
                    nw=nw+1
            else
                    go to 30
            end if
20      continue
30      continue
    if (nw == 0)then
!               no explicit repeat count, so assume a value of 1
            rcount=1
    else
            call ftc2ii(form(c1:c1+nw-1),rcount,status)
            if (status > 0)then
               call ftpmsg('Error in FTBNFM evaluating TFORM' &
               //' repeat value: '//cform)
               return
             end if
    end if

    c1=c1+nw

!       see if this is a variable length pointer column (e.g., 'rPt'); if so,
!       then add 1 to the starting search position in the TFORM string
    if (form(c1:c1) == 'P')then
            point=-1
            c1=c1+1
            rcount=1
    else
            point=1
    end if

!       now the chararcter at position c1 should be the data type code
    dattyp=form(c1:c1)

!       set the numeric datatype code
    if (dattyp == 'I')then
            dtype=21
            width=2
    else if (dattyp == 'J')then
            dtype=41
            width=4
    else if (dattyp == 'E')then
            dtype=42
            width=4
    else if (dattyp == 'D')then
            dtype=82
            width=8
    else if (dattyp == 'A')then
            dtype=16
    else if (dattyp == 'L')then
            dtype=14
            width=1
    else if (dattyp == 'X')then
            dtype=1
            width=1
    else if (dattyp == 'B')then
            dtype=11
            width=1
    else if (dattyp == 'C')then
            dtype=83
            width=8
    else if (dattyp == 'M')then
            dtype=163
            width=16
    else
!               unknown tform datatype code
            status=262
            call ftpmsg('Unknown Binary table TFORMn keyword '// &
                        'datatype: '//cform)
            return
    end if

!       set dtype negative if this is a variable length field ('P')
    dtype=dtype*point

!       if this is an ASCII field, determine its width
    if (dtype == 16)then
            c1=c1+1
            nw=0
            do 40 i=c1,nc
                if (form(i:i) >= '0' .and. form(i:i)<='9')then
                    nw=nw+1
            else
                    go to 50
            end if
40              continue
50              continue
            if (nw == 0)then
!                       no explicit width field, so assume that the
!                       width is the same as the repeat count
                    width=rcount
            else
                    tstat=status
                    call ftc2ii(form(c1:c1+nw-1),width,status)
                    if (status > 0)then
!                       unrecognized characters following the 'A', so ignore it
                           width=rcount
                           status=tstat
                    end if
            end if
    end if
end
subroutine ftc2d(cval,dval,status)
!
!*******************************************************************************
!
!! FTC2D converts a character string to a double precision value.
!
!       perform datatype conversion, if required

    character ( len = * ) cval
    integer ival,status
    character dtype
    logical lval
    character*16 sval
    double precision dval

    if (status > 0)return

    if (cval == ' ')then
!               null value string
            status = 204
            return
    end if

!       convert string to its intrinsic data type
    call ftc2x(cval,dtype,ival,lval,sval,dval,status)
    if (status > 0)return

    if (dtype == 'F')then
!               no datatype conversion required, so just return
    else if (dtype == 'I')then
!               convert from integer to double precision
            dval=ival
    else if (dtype == 'L')then
!               need to convert from logical to double precision
            if (lval)then
                    dval=1.
            else
                    dval=0.
            end if
    else if (dtype == 'C')then
!               can't convert a string to double precision, so return error
            dval=0
            status=406
            sval=cval
            call ftpmsg('Error in FTC2D evaluating this string '// &
            'as a double value: '//sval)
    end if
end
subroutine ftc2dd(cval,val,status)
!
!*******************************************************************************
!
!! FTC2DD converts a character string to double precision.
!
!       (assumes that the input string is left justified)
!       cval    c  input character string to be converted
!       val     d  output value
!       status  i  output error status (0 = OK)

    character ( len = * ) cval
    double precision val
    integer status,nleng
    character iform*8,sval*16

    if (status > 0)return

!       find length of the input double character string
    nleng=index(cval,' ')-1
    if (nleng == -1)nleng=len(cval)

!       construct the format statement to read the character string
    if (nleng <= 9)then
            write(iform,1000)nleng
1000            format('(F',I1,'.0)')
    else
            write(iform,1001)nleng
1001            format('(F',I2,'.0)')
    end if

    read(cval,iform,err=900)val
    return

900     status=409
    sval=cval
    call ftpmsg('Error in FTC2DD evaluating this string '// &
         'as a double: '//sval)
end
subroutine ftc2i(cval,ival,status)
!
!*******************************************************************************
!
!! FTC2I converts a character string to an integer.
!
!       perform datatype conversion, if required

    integer ival,status
    character ( len = * ) cval
    character dtype
    logical lval
    character sval*16
    double precision dval

    if (status > 0)return

    if (cval == ' ')then
!               null value string
            status = 204
            return
    end if

!       convert string to its intrinsic data type
    call ftc2x(cval,dtype,ival,lval,sval,dval,status)
    if (status > 0)return

    if (dtype == 'I')then
!               no datatype conversion required, so just return
    else if (dtype == 'F')then
!               need to convert from floating point to integer
            ival=dval
    else if (dtype == 'L')then
!               need to convert from logical to integer
            if (lval)then
                    ival=1
            else
                    ival=0
            end if
    else if (dtype == 'C')then
!               can't convert a string to an integer, so return error
            ival=0
            status=403
            sval=cval
    call ftpmsg('Error in FTC2I evaluating this string as an ' &
    //'integer: '//sval)
    end if
end
subroutine ftc2ii(cval,ival,status)
!
!*******************************************************************************
!
!! FTC2II converts a character string to an integer.
!
!       (assumes that the input string is left justified)

    integer ival,status,nleng
    character ( len = * ) cval
    character*8 iform

    if (status > 0)return

    if (cval == ' ')go to 900

!       find length of the input integer character string
    nleng=index(cval,' ')-1
    if (nleng == -1)nleng=len(cval)

!       construct the format statement to read the character string
    if (nleng <= 9)then
            write(iform,1000)nleng
1000            format('(I',I1,')')
    else
            write(iform,1001)nleng
1001            format('(I',I2,')')
    end if

    read(cval,iform,err=900)ival
    return

900     continue
!       work around for bug in the DEC Alpha VMS compiler
    if (cval(1:nleng) == '-2147483648')then
             ival=-2147483647 - 1
    else
             status=407
    end if
end
subroutine ftc2l(cval,lval,status)
!
!*******************************************************************************
!
!! FTC2L converts a character string to a logical value.
!
!       perform datatype conversion, if required

    logical lval
    integer ival,status
    character ( len = * ) cval
    character dtype
    character sval*16
    double precision dval

    if (status > 0)return

    if (cval == ' ')then
!               null value string
            status = 204
            return
    end if

!       convert string to its intrinsic data type
    call ftc2x(cval,dtype,ival,lval,sval,dval,status)
    if (status > 0)return

    if (dtype /= 'L')then
!              this is not a logical keyword, so return error
           status=404
           sval=cval
           call ftpmsg('Error in FTC2L evaluating this string '// &
            'as a logical value: '//sval)
    end if
end
subroutine ftc2ll(cval,lval,status)
!
!*******************************************************************************
!
!! FTC2LL converts a character string to a logical value.
!
!       (assumes that the input string is left justified)
    integer status
    logical lval
    character ( len = * ) cval

    if (status > 0)return

!       convert character string to logical
    if (cval(1:1) =='T')then
            lval=.true.
    else
!               any other character is considered false
            lval=.false.
    end if
end
subroutine ftc2r(cval,rval,status)
!
!*******************************************************************************
!
!! FTC2R converts a character string to a real value.
!
!       perform datatype conversion, if required

    character ( len = * ) cval
    real rval
    integer ival,status
    character dtype
    logical lval
    character*16 sval
    double precision dval

    if (status > 0)return

    if (cval == ' ')then
!               null value string
            status = 204
            return
    end if

!       convert string to its intrinsic data type
    call ftc2x(cval,dtype,ival,lval,sval,dval,status)
    if (status > 0)return

    if (dtype == 'F')then
!               convert from double to single precision
            rval=dval
    else if (dtype == 'I')then
!               convert from integer to real
            rval=ival
    else if (dtype == 'L')then
!               need to convert from logical to real
            if (lval)then
                    rval=1.
            else
                    rval=0.
            end if
    else if (dtype == 'C')then
!               can't convert a string to a real, so return error
            rval=0
            status=405
            sval=cval
            call ftpmsg('Error in FTC2R evaluating this string '// &
            'as a real value: '//sval)
    end if
end
subroutine ftc2rr(cval,val,status)
!
!*******************************************************************************
!
!! FTC2RR converts a character string to a real value.
!
!       (assumes that the input string is left justified)
!       cval    c  input character string to be converted
!       val     r  output value
!       status  i  output error status (0 = OK)

    character ( len = * ) cval
    real val
    integer status,nleng
    character iform*8,sval*16

    if (status > 0)return

    if (cval == ' ')go to 900

!       find length of the input real character string
    nleng=index(cval,' ')-1
    if (nleng == -1)nleng=len(cval)

!       construct the format statement to read the character string
    if (nleng <= 9)then
            write(iform,1000)nleng
1000            format('(F',I1,'.0)')
    else
            write(iform,1001)nleng
1001            format('(F',I2,'.0)')
    end if

    read(cval,iform,err=900)val
    return

900     status=408
    sval=cval
    call ftpmsg('Error in FTC2RR evaluating this string '// &
         'as a real: '//sval)
end
subroutine ftc2s(in,cval,status)
!
!*******************************************************************************
!
!! FTC2S converts an input quoted string to an unquoted string.
!
!       The first character of the input string must be a quote character (')
!       and at least one additional quote character must also be present in the
!       input string. This routine then simply outputs all the characters
!       between the first and last quote characters in the input string.
!
!       in      c  input quoted string
!       cval    c  output unquoted string
!       status  i  output error status (0=ok, 1=first quote missing,
!                  2=second quote character missing.

    character ( len = * ) in,cval
    integer length,i,j,i2,status
    character dtype

!       test for datatype
    call ftdtyp(in,dtype,status)
    if (status > 0)return
    if (dtype /= 'C')then
!               do no conversion and just return the raw character string
            cval=in
    else
!               convert character string to unquoted string

!               find closing quote character
            length=len(in)
            i2=length-1

            do i=length,2,-1
                    if (in(i:i) == '''') then
                      exit
                    end if
                    i2=i2-1
            end do

            if (i2 == 0)then
!                       there was no closing quote character
                    status=205
        call ftpmsg('The following keyword value string has no ' &
        //'closing quote:')
        call ftpmsg(in)
            else if (i2 == 1)then
!                       null string
                    cval=' '
            else
                    cval=in(2:i2)

!                       test for double single quote characters; if found,
!                       then  delete one of the quotes (FITS uses 2 single
!                       quote characters to represent a single quote)
                    i2=i2-2
                    do 30  i=1,i2
                        if (cval(i:i) == '''')then
                            if (cval(i+1:i+1) == '''')then
                               do j=i+1,i2
                                     cval(j:j)=cval(j+1:j+1)
                               end do
                               cval(i2:i2)=' '
                            end if
                        end if
30                      continue
            end if
    end if
end
subroutine ftc2x(cval,dtype,ival,lval,sval,dval,status)
!
!*******************************************************************************
!
!! FTC2X converts a character string into it intrinsic data type.
!
!       cval  c  input character string to be converted
!       dtype c  returned intrinsic datatype of the string (I,L,C,F)
!
!       one of  the following values is returned, corresponding to the
!       value of dtype:
!               ival i integer value
!               lval l logical value
!               sval c string value
!               dval d double precision value
!       statue i returned error status
!
  character ( len = * ) cval
  character dtype
  integer ival,status
  logical lval
  character ( len = * ) sval
  double precision dval
!
!  determine intrinsic datatype.
!
  call ftdtyp(cval,dtype,status)
!
!  convert string into its intrinsic datatype
!
  if (dtype == 'I')then
    call ftc2ii(cval,ival,status)
  else if (dtype == 'F')then
    call ftc2dd(cval,dval,status)
  else if (dtype == 'L')then
    call ftc2ll(cval,lval,status)
  else if (dtype == 'C')then
    call ftc2s(cval,sval,status)
  end if

  return
end
subroutine ftcdel(iunit,naxis1,naxis2,delbyt,fstbyt,status)
!
!*******************************************************************************
!
!! FTCDEL deletes a specified column by shifting the rows.
!
!       iunit   i  Fortran I/O unit number
!       naxis1  i  width in bytes of existing table
!       naxis2  i  number of rows in the table
!       delbyt  i  how many bytes to delete in each row
!       fstbyt  i  byte position in the row to delete the bytes (0=row start)
!       status  i  returned error status (0=ok)

    integer iunit,naxis1,naxis2,delbyt,fstbyt,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*5760 buff
    character xdummy(26240)
    common/ftheap/buff,xdummy
!

    integer ibuff,i,i1,i2,irow,newlen,nseg,nbytes,remain

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

    newlen=naxis1-delbyt

    if (newlen <= 5760)then
!
!       CASE #1: optimal case where whole new row fits in the work buffer
!
        i1=fstbyt+1
        i2=i1+delbyt
        do 10 irow=1,naxis2-1
!               read the row to be shifted
            call ftgtbs(iunit,irow,i2,newlen,buff,status)

!               set row length to its new value
            rowlen(ibuff)=newlen

!               write the row in the new place
            call ftptbs(iunit,irow,i1,newlen,buff,status)

!               reset row length to its original value
            rowlen(ibuff)=naxis1
10          continue

!           now do the last row
        remain=naxis1-(fstbyt+delbyt)
        if (remain > 0)then
!               read the row to be shifted
            call ftgtbs(iunit,naxis2,i2,remain,buff,status)

!               set row length to its new value
            rowlen(ibuff)=newlen

!               write the row in the new place
            call ftptbs(iunit,naxis2,i1,remain,buff,status)

!               reset row length to its original value
            rowlen(ibuff)=naxis1
        end if
    else
!
!       CASE #2:  whole row doesn't fit in work buffer; move row in pieces
!
        nseg=(newlen+5759)/5760

        do 40 irow=1,naxis2-1
            i1=fstbyt+1
            i2=i1+delbyt
            nbytes=newlen-(nseg-1)*5760

            do 30 i=1,nseg
!                   read the row to be shifted
                call ftgtbs(iunit,irow,i2,nbytes,buff,status)

!                   set row length to its new value
                rowlen(ibuff)=newlen

!                   write the row in the new place
                call ftptbs(iunit,irow,i1,nbytes,buff,status)

!                   reset row length to its original value
                rowlen(ibuff)=naxis1

                i1=i1+nbytes
                i2=i2+nbytes
                nbytes=5760
30              continue
40          continue

!           now do the last row
        remain=naxis1-(fstbyt+delbyt)
        if (remain > 0)then
            nseg=(remain+5759)/5760
            i1=fstbyt+1
            i2=i1+delbyt
            nbytes=remain-(nseg-1)*5760

            do 50 i=1,nseg
!                   read the row to be shifted
                call ftgtbs(iunit,naxis2,i2,nbytes,buff,status)

!                   set row length to its new value
                rowlen(ibuff)=newlen

!                   write the row in the new place
                call ftptbs(iunit,naxis2,i1,nbytes,buff,status)

!                   reset row length to its original value
                rowlen(ibuff)=naxis1

                i1=i1+nbytes
                i2=i2+nbytes
                nbytes=5760
50              continue
        end if
    end if
end
subroutine ftcdfl(iunit,status)
!
!*******************************************************************************
!
!! FTCDFL checks data unit fill values.
!
!
!       Check that the data unit is correctly filled with zeros or blanks
!       from the end of the data to the end of the current FITS 2880 byte block

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1994

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nf = 3000)
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*2880 chbuff
    character chfill,xdummy(29119)
    common/ftheap/chbuff,chfill,xdummy
!

    integer ibuff,filpos,nfill,i

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check if the data unit is null
    if (theap(ibuff) == 0)return

!       move to the beginning of the fill bytes
    filpos=dtstrt(ibuff)+theap(ibuff)+heapsz(ibuff)
    call ftmbyt(iunit,filpos,.true.,status)

!       get all the fill bytes
    nfill=(filpos+2879)/2880*2880-filpos
    if (nfill == 0)return

    call ftgcbf(iunit,nfill,chbuff,status)
    if (status > 0)then
       call ftpmsg('Error reading data unit fill bytes (FTCDFL).')
       return
    end if

!       set the correct fill value to be checked
    if (hdutyp(ibuff) == 1)then
!              this is an ASCII table; should be filled with blanks
           chfill=char(32)
    else
           chfill=char(0)
    end if

!       check for all zeros or blanks
    do 10 i=1,nfill
        if (chbuff(i:i) /= chfill)then
            status=255
            if (hdutyp(ibuff) == 1)then
                call ftpmsg('Warning: remaining bytes following'// &
                ' ASCII table data are not filled with blanks.')
            else
                call ftpmsg('Warning: remaining bytes following'// &
                ' data are not filled with zeros.')
            end if
            return
        end if
10      continue
end
subroutine ftchdu(iunit,status)
!
!*******************************************************************************
!
!! FTCHDU closes the Header Data Unit.
!
!       If we have write access to the file, then close the current HDU by:
!                 -padding remaining space in the header with blanks
!                 -writing the END keyword in the CHU
!                 -check the data fill values, and rewrite them if not correct
!                 -flushing the current buffer to disk
!                 -recover common block space containing column descriptors

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991

    integer iunit,status

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,pcount
    character*8 comm

!       ignore input status and close HDU regardless of input status value

    ibuff=bufnum(iunit)
!       check that unit number is valid (that file is actually opened)
    if (ibuff == 0)then
       if (status <= 0)status=101
       return
    end if

!       see if we have write access to this file
    if (wrmode(ibuff))then

!           if data has been written to heap, update the PCOUNT keyword
        if (heapsz(ibuff) > 0)then
           call ftgkyj(iunit,'PCOUNT',pcount,comm,status)
           if (heapsz(ibuff) > pcount)then
             call ftmkyj(iunit,'PCOUNT',heapsz(ibuff),'&',status)
           end if

!              update the variable length TFORM values if necessary
           call ftuptf(iunit, status)
        end if

!           rewrite the header END card and the following blank fill, and
!           insure that the internal data structure matches the keywords
        call ftrdef(iunit,status)

!           write the correct data fill values, if they are not already correct
        call ftpdfl(iunit,status)
    end if

!       set current column name buffer as undefined
    call ftrsnm

!       flush the buffers holding data for this HDU
    call ftflsh(ibuff,status)

!       recover common block space containing column descriptors for this HDU
    call ftfrcl(iunit,status)

    if (status > 0)then
        call ftpmsg('Error while closing current HDU (FTCHDU).')
    end if
end
subroutine ftchfl(iunit,status)
!
!*******************************************************************************
!
!! FTCHFL checks header fill values.
!
!       Check that the header unit is correctly filled with blanks from the
!       END card to the end of the current FITS 2880-byte block

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1994

    integer iunit,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,nblank,i,endpos
    character*80 rec
    logical gotend

    if (status > 0)return

    ibuff=bufnum(iunit)

!       calculate the number of blank keyword slots in the header
endpos=hdend(ibuff)
    nblank=(dtstrt(ibuff)-endpos)/80
!       move the i/o pointer to the end of the header keywords
    call ftmbyt(iunit,endpos,.true.,status)
!       find the END card (there may be blank keywords perceeding it)

    gotend=.false.
    do 10 i=1,nblank
            call ftgcbf(iunit,80,rec,status)
            if (rec(1:8) == 'END     ')then
                   if (gotend)then
!                          there is a duplicate END record
                       status=254
         call ftpmsg('Warning: Header fill area contains '// &
         'duplicate END card:')
                   end if
                   gotend=.true.
                   if (rec(9:80) /= ' ')then
!                          END keyword has extra characters
                       status=253
        call ftpmsg('Warning: END keyword contains '// &
        'extraneous non-blank characters:')
                   end if
             else if (gotend)then
                   if (rec /= ' ')then
!                          The fill area contains extraneous characters
                       status=254
         call ftpmsg('Warning: Header fill area contains '// &
         'extraneous non-blank characters:')
                    end if
            end if

            if (status > 0)then
                       call ftpmsg(rec)
                       return
            end if
10      continue
end
subroutine ftcins(iunit,naxis1,naxis2,delbyt,fstbyt,status)
!
!*******************************************************************************
!
!! FTCINS inserts DELBYT bytes after byte fstbyt in every row of the table.
!
!       iunit   i  Fortran I/O unit number
!       naxis1  i  width in bytes of existing table
!       naxis2  i  number of rows in the table
!       delbyt  i  how many bytes to insert in each row
!       fstbyt  i  byte position in the row to insert the bytes (0=row start)
!       status  i  returned error status (0=ok)

    integer iunit,naxis1,naxis2,delbyt,fstbyt,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*5760 buff
    character xdummy(26240)
    common/ftheap/buff,xdummy
!

    integer ibuff,i,i1,irow,newlen,fbyte,nseg,nbytes
    character cfill*1

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       select appropriate fill value
    if (hdutyp(ibuff) == 1)then
!           fill  header or ASCII table with space
        cfill=char(32)
    else
!           fill image or bintable data area with Null (0)
        cfill=char(0)
    end if

    newlen=naxis1+delbyt

    if (newlen <= 5760)then
!
!       CASE #1: optimal case where whole new row fits in the work buffer
!
!           write the correct fill value into the buffer
        do i=1,delbyt
            buff(i:i)=cfill
        end do
        i1=delbyt+1

!           first move the trailing bytes (if any) in the last row
        fbyte=fstbyt+1
        nbytes=naxis1-fstbyt
        call ftgtbs(iunit,naxis2,fbyte,nbytes,buff(i1:),status)

!           set row length to its new value
        rowlen(ibuff)=newlen

!           write the row (with leading fill bytes) in the new place
        nbytes=nbytes+delbyt
        call ftptbs(iunit,naxis2,fbyte,nbytes,buff,status)

!           reset row length to its original value
        rowlen(ibuff)=naxis1

!           now move the rest of the rows
        do 20 irow=naxis2-1,1,-1
!               read the row to be shifted (work backwards through the table)
            call ftgtbs(iunit,irow,fbyte,naxis1,buff(i1:),status)

!               set row length to its new value
            rowlen(ibuff)=newlen

!               write the row (with the leading fill bytes) in the new place
            call ftptbs(iunit,irow,fbyte,newlen,buff,status)

!               reset row length to its original value
            rowlen(ibuff)=naxis1
20          continue

    else
!
!       CASE #2:  whole row doesn't fit in work buffer; move row in pieces
!
!           first copy the data, then go back and write fill into the new column
!           start by copying the trailing bytes (if any) in the last row

        nbytes=naxis1-fstbyt
        nseg=(nbytes+5759)/5760
        fbyte=(nseg-1)*5760+fstbyt+1
        nbytes=naxis1-fbyte+1

        do 25 i=1,nseg
            call ftgtbs(iunit,naxis2,fbyte,nbytes,buff,status)

!               set row length to its new value
            rowlen(ibuff)=newlen

!               write the row in the new place
            call ftptbs(iunit,naxis2,fbyte+delbyt,nbytes, &
                        buff,status)

!               reset row length to its original value
            rowlen(ibuff)=naxis1

            fbyte=fbyte-5760
            nbytes=5760
25          continue

!           now move the rest of the rows
        nseg=(naxis1+5759)/5760

        do 40 irow=naxis2-1,1,-1
            fbyte=(nseg-1)*5760+fstbyt+1
            nbytes=naxis1-(nseg-1)*5760
            do 30 i=1,nseg
!                   read the row to be shifted (work backwards thru the table)
                call ftgtbs(iunit,irow,fbyte,nbytes,buff,status)

!                   set row length to its new value
                rowlen(ibuff)=newlen

!                   write the row in the new place
                call ftptbs(iunit,irow,fbyte+delbyt,nbytes, &
                            buff,status)

!                   reset row length to its original value
                rowlen(ibuff)=naxis1

                fbyte=fbyte-5760
                nbytes=5760
30              continue
40          continue

!           now write the fill values into the new column
        nbytes=min(delbyt,5760)
        do 50 i=1,nbytes
                buff(i:i)=cfill
50          continue

        nseg=(delbyt+5759)/5760

!           set row length to its new value
        rowlen(ibuff)=newlen

        do 70 irow=1,naxis2
            fbyte=fstbyt+1
            nbytes=delbyt-((nseg-1)*5760)
            do 60 i=1,nseg
!                   write the fill
                call ftptbs(iunit,irow,fbyte,nbytes,buff,status)
                fbyte=fbyte+nbytes
                nbytes=5760
60              continue
70          continue

!           reset the rowlength
        rowlen(ibuff)=naxis1
    end if
end
subroutine ftclos(iunit,status)
!
!*******************************************************************************
!
!! FTCLOS closes a FITS file that was previously opened with ftopen or ftinit.
!
!       iunit   i  Fortran I/O unit number
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,status

    logical keep

!       close the current HDU and pad the header with blanks
    call ftchdu(iunit,status)

!       don't attempt to close file if unit number is invalid
    if (status /= 101)then
!           close the file
        keep=.true.
        call ftclsx(iunit,keep,status)
    end if
end
subroutine ftcmps(templt,string,casesn,match,exact)
!
!*******************************************************************************
!
!! FTCMPS compares the template to the string and test if they match.
!
!       The strings are limited to 68 characters or less (the max. length
!       of a FITS string keyword value.  This routine reports whether
!       the two strings match and whether the match is exact or
!       involves wildcards.

!       this algorithm is very similar to the way unix filename wildcards
!       work except that this first treats a wild card as a literal character
!       when looking for a match.  If there is no literal match, then
!       it interpretes it as a wild card.  So the template 'AB*DE'
!       is considered to be an exact rather than a wild card match to
!       the string 'AB*DE'.  The '#' wild card in the template string will
!       match any consecutive string of decimal digits in the colname.

!       templt    C input template (may include ? or * wild cards)
!       string    C input string to be compared to template
!       casesn    L should comparison be case sensitive?
!       match     L (output) does the template match the string?
!       exact     L (output) are the strings an exact match (true) or
!                            is it a wildcard match (false)

!       written by Wm Pence, HEASARC/GSFC, December 1994
!       modified December 1995 to fix 2 bugs
!       modified Jan 1997 to support the # wild card

    character ( len = * ) templt,string
    logical casesn,match,exact
    character*68 temp,str
    integer tlen,slen,t1,s1

    tlen=len(templt)
    slen=len(string)
    tlen=min(tlen,68)
    slen=min(slen,68)

    match=.false.
    exact=.true.
    temp=templt
    str=string
    if (.not. casesn)then
        call ftupch(temp)
        call ftupch(str)
    end if

!       check for exact match
    if (temp == str)then
        match=.true.
        return
    end if

!       the strings are not identical, any match cannot be exact
    exact=.false.

    t1=1
    s1=1
10      continue
    if (t1 > tlen .or. s1 > slen)then
!           completely scanned one or both strings, so it must be a match
        match=.true.
        return
    end if

!       see if the characters in the 2 strings are an exact match
    if (temp(t1:t1) == str(s1:s1) .or. &
       (temp(t1:t1) == '?' .and. str(s1:s1) /= ' ') )then
!           The '?' wild card matches anything except a blank
        s1=s1+1
        t1=t1+1

    else if (temp(t1:t1) == '#' .and. (str(s1:s1) <= '9' &
        .and. str(s1:s1) >= '0' ))then
!           The '#' wild card matches any string of digits
        t1=t1+1
!           find the end of consecutive digits in the string
15          s1=s1+1
        if (str(s1:s1) <= '9' .and. str(s1:s1) >= '0')go to 15

    else if (temp(t1:t1) == '*')then
!           get next character from template and look for it in the string
        t1=t1+1
        if (t1 > tlen .or. (temp(t1:t1) == ' '))then
!               * is followed by a space, so a match is guaranteed
            match=.true.
            return
        end if

20          continue
        if (temp(t1:t1) == str(s1:s1))then
!               found a matching character
            t1=t1+1
            s1=s1+1
        else
!               increment the string pointer and try again
            s1=s1+1

!               return if hit end of string and failed to find a match
            if (s1 > slen)return

            go to 20
        end if

    else
!           match failed
        return
    end if
    go to 10
end
subroutine ftcmsg
!
!*******************************************************************************
!
!! FTCMSG clears the error message stack
    call ftxmsg(0,'dummy')
end
subroutine ftcopy (iunit,ounit,moreky,status)
!
!*******************************************************************************
!
!! FTCOPY copies the CHDU from IUNIT to the CHDU of OUNIT.
!
!       This will also reserve space in the header for MOREKY keywords
!       if MOREKY > 0.

!       iunit   i  fortran unit number of the input file to be copied
!       ounit   i  fortran unit number of the output file to be copied to
!       moreky  i  create space in header for this many more keywords
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Jan, 1992

    integer iunit,ounit,moreky,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,obuff,i,nkeys,nadd
    integer bitpix,naxis,naxes(99),pcount,gcount
    character hrec*80
    logical simple,extend

    if (status > 0)return

    if (iunit == ounit)then
            status=101
            return
    end if

    ibuff=bufnum(iunit)
    obuff=bufnum(ounit)

!       check that the output CHDU is empty
    call ftghsp(ounit,nkeys,nadd,status)
    if (nkeys /= 0)then
         call ftpmsg('Cannot copy HDU to a non-empty HDU')
         status = 201
         return
    end if

!       find out the number of keywords which exist in the input CHDU
    call ftghsp(iunit,nkeys,nadd,status)

!       copy the keywords one at a time to the output CHDU
    if ( (chdu(ibuff) == 1 .and. chdu(obuff) /= 1) .or. &
       (chdu(ibuff) /= 1 .and. chdu(obuff) == 1) )then
!               copy primary array to image extension, or vise versa

!               copy the required keywords:
            simple=.true.
            call ftghpr(iunit,99,simple,bitpix,naxis, &
            naxes,pcount,gcount,extend,status)
            if (status > 0)return
            extend=.true.
            call ftphpr(ounit,simple,bitpix,naxis, &
            naxes,pcount,gcount,extend,status)
            if (status > 0)return

!               copy remaining keywords, excluding pcount, gcount and extend
            do 10 i=naxis+4,nkeys
                call ftgrec(iunit,i,hrec,status)
                if (hrec(1:8) /= 'PCOUNT  ' .and. &
                    hrec(1:8) /= 'GCOUNT  ' .and. &
                    hrec(1:8) /= 'EXTEND  ')then
                       call ftprec(ounit,hrec,status)
                end if
10              continue
    else
!               just copy all the keys exactly from the input file to the output
            do 20 i=1,nkeys
                call ftgrec(iunit,i,hrec,status)
                call ftprec(ounit,hrec,status)
20              continue
    end if

!       reserve space for more keywords (if moreky > 0)
    call fthdef(ounit,moreky,status)

!       now ccopy the data from the input CHDU to the output CHDU
    call ftcpdt(iunit,ounit,status)

end
subroutine ftcpdt(iunit,ounit,status)
!
!*******************************************************************************
!
!! FTCPDT copies the data from the IUNIT CHDU to the data of the OUNIT CHDU.
!
!       This will overwrite any data already in the OUNIT CHDU.

!       iunit   i  fortran unit number of the input file to be copied
!       ounit   i  fortran unit number of the output file to be copied to
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Aug 1993

    integer iunit,ounit,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    character*2880 cbuff
    character xdummy(29120)
    common/ftheap/cbuff,xdummy
!

    integer ibuff,obuff,nblock,i

    if (status > 0)return

    if (iunit == ounit)then
            status=101
            return
    end if

    ibuff=bufnum(iunit)
    obuff=bufnum(ounit)

!       determine HDU structure as defined by keywords in output file
    call ftrdef(ounit,status)

!       Calculate the number of bytes to be copied.  By definition there
!       will be an integral number of 2880-byte logical blocks to be copied
    nblock=(hdstrt(ibuff,chdu(ibuff)+1)-dtstrt(ibuff))/2880

    if (nblock > 0)then
!           move to the beginning of the data in the input and output files
        call ftmbyt(iunit,dtstrt(ibuff),.false.,status)
        call ftmbyt(ounit,dtstrt(obuff),.true.,status)

!           now copy the data one block at a time
        do 30 i=1,nblock
            call ftgcbf(iunit,2880,cbuff,status)
            call ftpcbf(ounit,2880,cbuff,status)
30          continue
    end if
end
subroutine ftcrep(comm,comm1,repeat)
!
!*******************************************************************************
!
!! FTCREP checks if the first comment string is to be repeated for all keywords.
!
!       (if the last non-blank character is '&', then it is to be repeated)

!       comm    c  input comment string
!       OUTPUT PARAMETERS:
!       comm1   c  output comment string, = COMM minus the last '&' character
!       repeat  l  true if the last character of COMM was the '&" character
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) comm,comm1
    logical repeat
    integer i,j

    repeat=.false.
    j=len(comm)
    do i=j,1,-1
            if (comm(i:i) /= ' ')then
                    if (comm(i:i) == '&')then
                            comm1=comm(1:i-1)
                            repeat=.true.
                    end if
                    return
            end if
    end do
end
subroutine ftcrhd(iunit,status)
!
!*******************************************************************************
!
!! FTCRHD creates a header data unit.
!
!       'CReate Header Data unit'
!       create, initialize, and move the i/o pointer to a new extension at
!       the end of the FITS file.

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991

    integer iunit,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff

    if (status > 0)return

!       close the current HDU
    call ftchdu(iunit,status)
    if (status > 0)return

    ibuff=bufnum(iunit)

!       check that we haven't exceeded the maximum allowed number of extensions
    if (maxhdu(ibuff)+1 >= ne)then
            status=301
            return
    end if

!       move to the end of the highest known extension
    call ftmbyt(iunit,hdstrt(ibuff,maxhdu(ibuff)+1),.true.,status)

!       initialize various parameters about the CHDU
    maxhdu(ibuff)=maxhdu(ibuff)+1
    chdu(ibuff)=maxhdu(ibuff)
    nxthdr(ibuff)=hdstrt(ibuff,chdu(ibuff))
!       the logical location of the END record at the start of the header
    hdend(ibuff)=nxthdr(ibuff)
!       the data start location is undefined
    dtstrt(ibuff)=-2000000000
end
subroutine ftcsum(iunit,nrec,sum,status)
!
!*******************************************************************************
!
!! FTCSUM calculates a 32-bit 1's complement checksum of the FITS 2880-byte blocks.
!
!
!       This Fortran algorithm is based on the C algorithm developed by Rob
!       Seaman at NOAO that was presented at the 1994 ADASS conference, to be
!       published in the Astronomical Society of the Pacific Conference Series.

!       This uses a 32-bit 1's complement checksum in which the overflow bits
!       are permuted back into the sum and therefore all bit positions are
!       sampled evenly.  In this Fortran version of the original C algorithm,
!       a double precision value (which has at least 48 bits of precision)
!       is used to accumulate the checksum because standard Fortran does not
!       support an unsigned integer datatype.

!       iunit   i  fortran unit number
!       nrec    i  number of FITS 2880-byte blocks to be summed
!       sum     d  check sum value (initialize to zero before first call)
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Sept, 1994

    integer iunit,nrec,status,i,j,hibits,i4vals(720)
    double precision sum,word32
    parameter (word32=4.294967296D+09)
!       word32 is equal to 2**32

    if (status > 0)return

!       Sum the specified number of FITS 2880-byte records.  This assumes that
!       the FITSIO file pointer points to the start of the records to be summed.
    do 30 j=1,nrec

!           read the record as 720 pixel I*4 vector (do byte swapping if needed)
        call ftgi4b(iunit,720,4,i4vals,status)

        do i=1,720
            if (i4vals(i) >= 0)then
                    sum=sum+i4vals(i)
            else
!                       sign bit is set, so add the equalvalent unsigned value
                    sum=sum+(word32+i4vals(i))
            end if
        end do

!           fold any overflow bits beyond 32 back into the word
20          hibits=sum/word32
        if (hibits > 0)then
            sum=sum-(hibits*word32)+hibits
            go to 20
        end if
30      continue
end
subroutine ftd2e(val,dec,cval,vlen,status)
!
!*******************************************************************************
!
!! FTD2E converts a double precision value to an E format character string.
!
!       If it will fit, the value field will be 20 characters wide;
!       otherwise it will be expanded to up to 35 characters, left
!       justified.
!
!       val     d  input value to be converted
!       dec     i  number of decimal places to display in output string
!       cval    c  output character string
!       vlen    i  length of output string
!       status  i  output error status (0 = OK)

    double precision val
    integer dec,vlen,status
    character*35 cval,form*10

    vlen = 1
    if (status > 0)return

    if (dec >= 1 .and. dec <= 9)then
            vlen=20
            write(form,2000)dec
2000            format('(1pe20.',i1,')')
    else if (dec >= 10 .and. dec <= 28)then
            if (val < 0.)then
                vlen=max(20,dec+7)
            else
                vlen=max(20,dec+6)
            end if
            write(form,2001)vlen,dec
2001            format('(1pe',i2,'.',i2,')')
    else
!               illegal number of decimal places were specified
            status=411
            call ftpmsg('Error in FTR2E: number of decimal places ' &
                        //'is less than 1 or greater than 28.')
            return
    end if

    write(cval,form,err=900)val
    if (cval(1:1) == '*')go to 900
    return

900     status=402
    call ftpmsg('Error in FTD2E converting double to En.m string.')
end
subroutine ftd2f(val,dec,cval,status)
!
!*******************************************************************************
!
!! FTD2F converts a double precision value to F20.* format character string.
!
!       NOTE: some precision may be lost
!       val     d  input value to be converted
!       dec     i  number of decimal places to display in output string
!       cval    c  output character string
!       status  i  output error status (0 = OK)

    double precision val
    integer dec,status
    character*20 cval,form*8

    if (status > 0)return

    if (dec >= 0 .and. dec <= 9)then
            write(form,2000)dec
2000            format('(f20.',i1,')')
    else if (dec >= 10 .and. dec <18)then
            write(form,2001)dec
2001            format('(f20.',i2,')')
    else
!               illegal number of decimal places were specified
            status=411
            call ftpmsg('Error in FTD2F: number of decimal places ' &
                        //'is less than 0 or greater than 18.')
            return
    end if

    write(cval,form,err=900)val
    if (cval(1:1) == '*')go to 900
    return
900     status=402
    call ftpmsg('Error in FTD2F converting double to F20. string.')
end
subroutine ftdblk(ounit,nblock,hdrdat,status)
!
!*******************************************************************************
!
!! FTDBLK deletes 2880-byte FITS blocks at the end of the current header or data.
!
!       ounit   i  fortran output unit number
!       nblock  i  number of 2880-byte blocks to be deleted
!       hdrdat  i  delete space at end of header (0) or data (1)
!       status  i  returned error status (0=ok)

    integer ounit,nblock,hdrdat,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    character*2880 buff
    character xdummy(29120)
    common/ftheap/buff,xdummy
!

    integer ibuff,jpoint,i,tstat

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

!       get address of first block to be deleted/overwritten
    if (hdrdat == 0)then
        jpoint=dtstrt(ibuff)-2880*nblock
    else
        jpoint=hdstrt(ibuff,chdu(ibuff)+1)-2880*nblock
    end if

!       move each block up, until we reach the end of file
10      continue
!           move to the read start position
        tstat=status
        call ftmbyt(ounit,jpoint+nblock*2880,.false.,status)

!           read one 2880-byte FITS logical record
        call ftgcbf(ounit,2880,buff,status)

!           check for end of file
        if (status == 107)then
            status=tstat
            go to 20
        end if

!           move back to the write start postion
        call ftmbyt(ounit,jpoint,.false.,status)

!           write the 2880-byte FITS logical record
        call ftpcbf(ounit,2880,buff,status)

!           check for error
        if (status > 0)then
            call ftpmsg('Error deleting FITS blocks (FTDBLK)')
            return
        end if

!           increment pointer to next block and loop back
        jpoint=jpoint+2880
        go to 10
20      continue

!       now fill the last nblock blocks with zeros;  initialize the  buffer
    do 30 i=1,2880
        buff(i:i)=char(0)
30      continue

!       move back to the write start postion
    call ftmbyt(ounit,jpoint,.false.,status)

!       write the 2880-byte block NBLOCK times.
    do 40 i=1,nblock
        call ftpcbf(ounit,2880,buff,status)
40      continue

    if (hdrdat == 0)then
!           recalculate the starting location of the current data unit, if moved
        dtstrt(ibuff)=dtstrt(ibuff)-2880*nblock
    end if

!       recalculate the starting location of all subsequent HDUs
    do 50 i=chdu(ibuff)+1,maxhdu(ibuff)+1
        hdstrt(ibuff,i)=hdstrt(ibuff,i)-2880*nblock
50      continue

    if (status > 0)then
        call ftpmsg('Error deleting FITS block(s) (FTDBLK)')
    end if
end
subroutine ftdcol(iunit,colnum,status)
!
!*******************************************************************************
!
!! FTDCOL deletes a column from a table.
!
!       iunit   i  Fortran I/O unit number
!       colnum  i  number of of the column to be deleted
!       status  i  returned error status (0=ok)

    integer iunit,colnum,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,typhdu,delbyt,fstbyt,sp,tflds,i
    integer naxis1,naxis2,size,freesp,nblock,tbc
    character comm*70,keynam*8

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       test that the CHDU is an ASCII table or BINTABLE
    typhdu=hdutyp(ibuff)
    if (typhdu /= 1 .and. typhdu /= 2)then
            status=235
            call ftpmsg('Can only delete column from TABLE '// &
            'or BINTABLE extension (FTDCOL)')
            return
    end if

!       check if column number exists in the table
    tflds=tfield(ibuff)
    if (colnum < 1 .or. colnum > tflds)then
        status=302
        return
    end if

!       get the starting byte position of the column (=zero for first column)
    fstbyt=tbcol(colnum+tstart(ibuff))

!       find the width of the column
    if (typhdu == 1)then
!           tnull is used to store the width of the ASCII column field
!           NOTE: ASCII columns may not be in physical order, or may overlap.

        delbyt=tnull(colnum+tstart(ibuff))

!           delete the space(s) between the columns, if there are any.
        if (colnum < tflds)then
!               check for spaces between following column
            sp=tbcol(colnum+1+tstart(ibuff))-tbcol(colnum+ &
               tstart(ibuff))-delbyt
            if (sp > 0)then
                delbyt=delbyt+1
            end if
        else if (colnum > 1)then
!               check for space between the last and next to last columns
            sp=tbcol(colnum+tstart(ibuff))-tbcol(colnum-1+ &
               tstart(ibuff))-tnull(colnum-1+tstart(ibuff))
            if (sp > 0)then
               delbyt=delbyt+1
               fstbyt=fstbyt-1
            end if
        end if
    else
        if (colnum < tflds)then
            delbyt=tbcol(colnum+1+tstart(ibuff))- &
                   tbcol(colnum+tstart(ibuff))
        else
            delbyt=rowlen(ibuff)-tbcol(colnum+tstart(ibuff))
        end if
    end if

!       get current size of the table
    naxis1=rowlen(ibuff)
    call ftgkyj(iunit,'NAXIS2',naxis2,comm,status)

!       Calculate how many FITS blocks (2880 bytes) need to be deleted
    size=theap(ibuff)+heapsz(ibuff)
    freesp=(delbyt*naxis2) + ((size+2879)/2880)*2880 - size
    nblock=freesp/2880

!       shift each row up, deleting the desired column
    call ftcdel(iunit,naxis1,naxis2,delbyt,fstbyt,status)

!       shift the heap up and update pointer to start of heap
    size=delbyt*naxis2
    call fthpup(iunit,size,status)

!       delete the needed number of new FITS blocks at the end of the HDU
    if (nblock > 0)call ftdblk(iunit,nblock,1,status)

    if (typhdu == 1)then
!           adjust the TBCOL values of the remaining columns
        do 10 i=1,tflds
            call ftkeyn('TBCOL',i,keynam,status)
            call ftgkyj(iunit,keynam,tbc,comm,status)
            if (tbc > fstbyt)then
                 tbc=tbc-delbyt
                 call ftmkyj(iunit,keynam,tbc,'&',status)
            end if
10          continue
    end if

!       update the mandatory keywords
    call ftmkyj(iunit,'TFIELDS',tflds-1,'&',status)
    call ftmkyj(iunit,'NAXIS1',naxis1-delbyt,'&',status)

!       delete the index keywords starting with 'T' associated with the
!       deleted column and subtract 1 from index of all higher keywords
    call ftkshf(iunit,colnum,tflds,-1,status)

!       parse the header to initialize the new table structure
    call ftrdef(iunit,status)
end
subroutine ftddef(ounit,bytlen,status)
!
!*******************************************************************************
!
!! FTDDEF redefines the length of the data unit.
!
!       Data DEFinition
!       re-define the length of the data unit
!       this simply redefines the start of the next HDU
!
!       ounit   i  Fortran I/O unit number
!       bytlen  i  new length of the data unit, in bytes
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,bytlen,status

!
    integer nb,ne,nf
    parameter (nf = 3000)
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff

    if (status > 0)return

    ibuff=bufnum(ounit)

    if (dtstrt(ibuff) < 0)then
!               freeze the header at its current size
            call fthdef(ounit,0,status)
    end if

    hdstrt(ibuff,chdu(ibuff)+1)= &
            dtstrt(ibuff)+(bytlen+2879)/2880*2880

!       initialize the fictitious heap starting address (immediately following
!       the array data) and a zero length heap.  This is used to find the
!   end of the data when checking the fill values in the last block.
    heapsz(ibuff)=0
    theap(ibuff)=bytlen
end
subroutine ftdelt(iunit,status)
!
!*******************************************************************************
!
!! FTDELT deletes a FITS file that was previously opened with ftopen or ftinit.
!
!       iunit   i  Fortran I/O unit number
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1994

    integer iunit,status,ibuff

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

!       ignore input status, and delete file regardless of status value

    ibuff=bufnum(iunit)

!       set current column name buffer as undefined
    call ftrsnm

!       flush the buffers holding data for this HDU
    call ftflsh(ibuff,status)

!       recover common block space containing column descriptors for this HDU
    call ftfrcl(iunit,status)

!       delete the file
    call ftclsx(iunit,.false.,status)
end
subroutine ftdhdu(ounit,typhdu,status)
!
!*******************************************************************************
!
!! FTDHDU deletes the current HDU (as long as it is not the primary array).
!
!       ounit   i  fortran output unit number
!       typhdu  i  type of the new CHDU, after deleting the old CHDU
!       status  i  returned error status (0=ok)

    integer ounit,typhdu,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,ibuff,nhdu,nblock

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

    nhdu=chdu(ibuff)
    if (nhdu == 1)then
!            cannot delete the primary array
         status=301
         return
    end if

!       close the CHDU first, to flush buffers and free memory
    call ftchdu(ounit,status)

!       how many blocks to delete?
    nblock=(hdstrt(ibuff,nhdu+1)-hdstrt(ibuff,nhdu))/2880
    if (nblock < 1)return

!       delete the blocks
    call ftdblk(ounit,nblock,1,status)
    if (status > 0)return

!       decrement the number of HDUs in the file and their starting address
    do 10 i=nhdu+1,maxhdu(ibuff)
            hdstrt(ibuff,i)=hdstrt(ibuff,i+1)
10      continue
    maxhdu(ibuff)=maxhdu(ibuff)-1

!       try reinitializing the CHDU, if there is one
    call ftrhdu(ounit,typhdu,status)
    if (status > 0)then
!            there is no HDU after the one we just deleted so move back one HDU
         status=0
         call ftcmsg
         call ftgext(ounit,nhdu-1,typhdu,status)
    end if
end
subroutine ftdkey(iunit,keynam,status)
!
!*******************************************************************************
!
!! FTDKEY deletes a header keyword.
!
!       iunit   i  fortran output unit number
!       keynam  c  keyword name    ( 8 characters, cols.  1- 8)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    character ( len = * ) keynam
    integer iunit,status,tstat,i,lenval,nkeys,keypos
    character keybuf*80,strval*70,comm*8,value*70,bslash*1,kname*8

    if (status > 0)return

!       have to use 2 \\'s because the SUN compiler treats 1 \ as an escape
    bslash='\\'

!       find the keyword to be deleted
    call ftgcrd(iunit,keynam,keybuf,status)
    if (status == 202)then
        kname=keynam
        call ftpmsg('FTDKEY could not find the '//kname// &
        ' keyword to be deleted.')
        return
    end if

!       get the position of the keyword in the header
    call ftghps(iunit,nkeys,keypos,status)
    keypos=keypos-1

!       get position of last character in value string to see if it is a \ or &
    if (status > 0)return
    tstat=status
    call ftpsvc(keybuf,strval,comm,status)
    call ftc2s(strval,value,status)
    if (status > 0)status=tstat

    lenval=1
    do i=70,1,-1
            if (value(i:i) /= ' ')then
                    lenval=i
                    exit
            end if
    end do

20      continue
!
!       now delete this keyword
!
    call ftdrec(iunit,keypos,status)
    if (status > 0)return

!       test if this keyword was also continued
    if (value(lenval:lenval) == bslash .or. &
            value(lenval:lenval) == '&')then
            call ftgnst(iunit,value,lenval,comm,status)
            if (lenval > 0)go to 20
    end if
end
subroutine ftdrec(ounit,pos,status)
!
!*******************************************************************************
!
!! FTDREC deletes a keyword record at position POS from header.
!
!       ounit   i  fortran output unit number
!       pos     i  position of keyword to be deleted (1 = first keyword)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Jan 1995

    integer ounit,pos,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    character*80 keybuf,keytmp
    integer ibuff,i,j,nshift

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

    if (pos < 1 .or. pos > &
       (hdend(ibuff)-hdstrt(ibuff,chdu(ibuff)))/80)then
            status=203
            return
    end if

    nxthdr(ibuff)=hdstrt(ibuff,chdu(ibuff))+(pos-1)*80

!       calculate number of header records following the deleted record
    nshift=(hdend(ibuff)-nxthdr(ibuff))/80

!       go through header shifting each 80 byte record up one place to
!       fill in the gap created by the deleted keyword
    j=hdend(ibuff)
    keybuf=' '
    do i=1,nshift
            j=j-80
!               read current record contents
            call ftmbyt(ounit,j,.false.,status)
            call ftgcbf(ounit,80,keytmp,status)
!               overwrite with new contents
            call ftmbyt(ounit,j,.false.,status)
            call ftpcbf(ounit,80,keybuf,status)
            keybuf=keytmp
    end do

!       update end-of-header pointer
    hdend(ibuff)=hdend(ibuff)-80

100     continue
end
subroutine ftdrow(iunit,frow,nrows,status)
!
!*******************************************************************************
!
!! FTDROW deletes NROWS rows from a table, beginning with row FROW.
!
!       iunit   i  Fortran I/O unit number
!       frow    i  first row number to be delete
!       nrows   i  number of rows to be deleted
!       status  i  returned error status (0=ok)

    integer iunit,frow,nrows,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,naxis1,naxis2,size,freesp,nblock,row
    character comm*8

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       test that the CHDU is an ASCII table or BINTABLE
    if (hdutyp(ibuff) /= 1 .and. hdutyp(ibuff) /= 2)then
            status=235
            call ftpmsg('Can only delete rows from TABLE or '// &
            'BINTABLE extension (FTDROW)')
            return
    end if

!       get current size of the table
    call ftgkyj(iunit,'NAXIS1',naxis1,comm,status)
    call ftgkyj(iunit,'NAXIS2',naxis2,comm,status)

    if (nrows < 0)then
             status=306
             call ftpmsg('Cannot delete negative number of ' // &
             'rows in the table (FTDROW)')
             return
    else if (frow+nrows-1 > naxis2)then
             status=307
             call ftpmsg('Specified number of rows to delete ' &
             //'exceeds number of rows in table (FTDROW)')
             return
    else if (nrows == 0)then
             return
    else if (frow > naxis2)then
            status=307
            call ftpmsg('First row to delete is greater'// &
              ' than the number of rows in the table (FTDROW)')
            return
    else if (frow <= 0)then
            status=307
            call ftpmsg('Delete starting row number is less ' &
            //'than 1 (FTDROW)')
            return
    end if

!       Calculate how many FITS blocks (2880 bytes) need to be deleted
    size=theap(ibuff)+heapsz(ibuff)
    freesp=((size+2879)/2880)*2880 - size + naxis1*nrows
    nblock=freesp/2880

!       shift the rows up
    row=frow+nrows
    call ftrwup(iunit,row,naxis2,nrows,status)

!       shift the heap up
    size=naxis1*nrows
    call fthpup(iunit,size,status)

    if (nblock > 0)call ftdblk(iunit,nblock,1,status)

!       update the NAXIS2 keyword
    naxis2=naxis2-nrows
    call ftmkyj(iunit,'NAXIS2',naxis2,'&',status)
end
subroutine ftdsum(string,complm,sum)
!
!*******************************************************************************
!
!! FTDSUM decodes the 32 bit checksum.

!       If complm=.true., then the complement of the sum will be decoded.

!       This Fortran algorithm is based on the C algorithm developed by Rob
!       Seaman at NOAO that was presented at the 1994 ADASS conference, to be
!       published in the Astronomical Society of the Pacific Conference Series.
!
!       sum     d  checksum value
!       complm  l  encode the complement of the sum?
!       string  c  output ASCII encoded check sum
!       sum     d  checksum value
!
!       written by Wm Pence, HEASARC/GSFC, May, 1995

    double precision sum,all32,word32,factor(4)
    character*16 string,tmpstr
    integer offset,i,j,k,temp,hibits
    logical complm

!       all32 equals a 32 bit unsigned integer with all bits set
!       word32 is equal to 2**32
    parameter (all32=4.294967295D+09)
    parameter (word32=4.294967296D+09)

!       ASCII 0 is the offset value
    parameter (offset=48)

    data factor/16777216.0D+00,65536.0D+00,256.0D+00,1.0D+00/

    sum=0

!       shift the characters 1 place to the left, since the FITS character
!       string value starts in column 12, which is not word aligned
    tmpstr(1:15)=string(2:16)
    tmpstr(16:16)=string(1:1)

!       convert characters from machine's native character coding sequence
!       to ASCII codes.   This only affects IBM mainframe computers
!       that do not use ASCII for the internal character representation.
!        call ftc2as(tmpstr,16)

!       subtract the offset from each byte and interpret each 4 character
!       string as a 4-byte unsigned integer; sum the 4 integers
    k=0
    do i=1,4
      do j=1,4
        k=k+1
        temp=ichar(tmpstr(k:k))-offset
        sum=sum+temp*factor(j)
      end do
    end do

!       fold any overflow bits beyond 32 back into the word
30      hibits=sum/word32
    if (hibits > 0)then
            sum=sum-(hibits*word32)+hibits
            go to 30
     end if

    if (complm)then
!           complement the 32-bit unsigned integer equivalent (flip every bit)
        sum=all32-sum
    end if
end
subroutine ftdtyp(value,dtype,status)
!
!*******************************************************************************
!
!! FTDTYP determines datatype of a FITS value field,
!
!       This assumes value field conforms to FITS standards and may not
!          detect all invalid formats.
!       value   c  input value field from FITS header record only,
!                  (usually the value field is in columns 11-30 of record)
!                  The value string is left justified.
!       dtype   c  output type (C,L,I,F) for Character string, Logical,
!                    Integer, Floating point, respectively
!
!       written by Wm Pence, HEASARC/GSFC, February 1991

    character ( len = * )value,dtype
    integer status

    if (status > 0)return

    dtype=' '

    if (value(1:1) == '''')then
!               character string
            dtype='C'
    else if (value(1:1)=='T' .or. value(1:1)=='F')then
!               logical
            dtype='L'
    else if (index(value,'.') > 0)then
!               floating point
            dtype='F'
    else
!               assume it must be an integer, since it isn't anything else
            dtype='I'
    end if
end
subroutine ftesum(sum,complm,string)
!
!*******************************************************************************
!
!! FTESUM encodes the 32 bit checksum.
!
!  It does this by converting every
!       2 bits of each byte into an ASCII character (32 bit word encoded
!       as 16 character string).   Only ASCII letters and digits are used
!       to encode the values (no ASCII punctuation characters).

!       If complm=.true., then the complement of the sum will be encoded.

!       This Fortran algorithm is based on the C algorithm developed by Rob
!       Seaman at NOAO that was presented at the 1994 ADASS conference, to be
!       published in the Astronomical Society of the Pacific Conference Series.
!
!       sum     d  checksum value
!       complm  l  encode the complement of the sum?
!       string  c  output ASCII encoded check sum
!
!       written by Wm Pence, HEASARC/GSFC, Sept, 1994

    double precision sum,tmpsum,all32
    character ( len = * ) string
    character tmpstr*16
    integer offset,exclud(13),nbyte(4),ch(4),i,j,k
    integer quot,remain,check,nc
    logical complm

!       all32 equals a 32 bit unsigned integer with all bits set
    parameter (all32=4.294967295D+09)

!       ASCII 0 is the offset value
    parameter (offset=48)

!       this is the list of ASCII punctutation characters to be excluded
    data exclud/58,59,60,61,62,63,64,91,92,93,94,95,96/

!       initialize input string (in case it is greater than 16 chars long)
    string = ' '

    if (complm)then
!           complement the 32-bit unsigned integer equivalent (flip every bit)
        tmpsum=all32-sum
    else
        tmpsum=sum
    end if

!       separate each 8-bit byte into separate integers
    nbyte(1)=tmpsum/16777216.
    tmpsum=tmpsum-nbyte(1)*16777216.
    nbyte(2)=tmpsum/65536.
    tmpsum=tmpsum-nbyte(2)*65536.
    nbyte(3)=tmpsum/256.
    nbyte(4)=tmpsum-nbyte(3)*256.

!       encode each 8-bit integer as 4-characters
    do i=1,4
            quot=nbyte(i)/4+offset
            remain=nbyte(i) - (nbyte(i)/4*4)
            ch(1)=quot+remain
            ch(2)=quot
            ch(3)=quot
            ch(4)=quot

!               avoid ASCII punctuation characters by incrementing and
!               decrementing adjacent characters thus preserving checksum value
10              check=0
                do k=1,13
                    do j=1,4,2
                       if (ch(j)   == exclud(k) .or. &
                           ch(j+1) == exclud(k))then
                           ch(j)=ch(j)+1
                           ch(j+1)=ch(j+1)-1
                           check=1
                       end if
                    end do
                 end do

!               keep repeating, until all punctuation character are removed
            if (check /= 0)go to 10

!               convert the byte values to the equivalent ASCII characters
            do j=0,3
                nc=4*j+i
                tmpstr(nc:nc)=char(ch(j+1))
            end do
    end do

!       shift the characters 1 place to the right, since the FITS character
!       string value starts in column 12, which is not word aligned
    string(1:1) =tmpstr(16:16)
    string(2:16)=tmpstr(1:15)

!       convert characters from ASCII codes to machine's native character
!       coding sequence.  (The string gets converted back to ASCII when it
!       is written to the FITS file). This only affects IBM mainframe computers
!       that do not use ASCII for the internal character representation.
!        call ftas2c(string,16)
end
subroutine ftfiou(iounit,status)
!
!*******************************************************************************
!
!! FTFIOU frees a specified logical unit number; if iounit=-1, free all units.
!
    integer iounit,status

    if (status > 0)return

    call ftxiou(iounit,status)
end
subroutine ftflus(iunit,status)
!
!*******************************************************************************
!
!! FTFLUS flushes all the data in the current FITS file to disk.
!
!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, March, 1996

    integer iunit,extno,xtend,status

    if (status > 0)return

!       get the current HDU number
    call ftghdn(iunit, extno)

!       close out the current HDU
    call ftchdu(iunit,status)
    if (status > 0)then
        call ftpmsg('FTFLUS could not close the current HDU.')
        return
    end if

!       reopen the same HDU
    call ftgext(iunit,extno,xtend,status)
    if (status > 0)then
        call ftpmsg('FTFLUS could not reopen the current HDU.')
        return
    end if
end
subroutine ftfrcl(iunit,status)
!
!*******************************************************************************
!
!! FTFRCL frees space used by column descriptors in the HDU being closed.
!
!
!  The various parameters
!       describing each table column (e.g., starting byte address, datatype,
!       tscale, tzero, etc.) are stored in 1-D arrays, and the tstart
!       parameter gives the starting element number in the arrays
!       for each unit number.  If a table is closed, then all the
!       descriptors for that table columns must be overwritten by
!       shifting any descriptors that follow it in the 1-D arrays to the left.

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC,May, 1995

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer ibuff,n2shft,i,j1,j2

!       ignore input status and flush columns regardless of input status value

    ibuff=bufnum(iunit)

    if (status == -999)then
!           just initialize the descriptors as undefined
        tstart(ibuff)=-1
    else if (tstart(ibuff) < 0)then
!           descriptors are already undefined; just return
    else if (tfield(ibuff) == 0)then
!           table had no columns so just reset pointers as undefined
        tstart(ibuff)=-1
        dtstrt(ibuff)=-2000000000
    else
!           calc number of descriptors to be shifted over the recovered space
        n2shft=nxtfld-(tstart(ibuff)+tfield(ibuff))

        if (n2shft > 0)then
            j1=tstart(ibuff)
            j2=j1+tfield(ibuff)
            do 10 i=1,n2shft
!                   shift the descriptors
                j1=j1+1
                j2=j2+1
                tbcol(j1)=tbcol(j2)
                tdtype(j1)=tdtype(j2)
                trept(j1)=trept(j2)
                tscale(j1)=tscale(j2)
                tzero(j1)=tzero(j2)
                tnull(j1)=tnull(j2)
                cnull(j1)=cnull(j2)
                cform(j1)=cform(j2)
10              continue
        end if

!           update pointer to next vacant column discriptor location
        nxtfld=nxtfld-tfield(ibuff)

!           update starting pointer for other opened files
        do 20 i=1,nb
            if (tstart(i) > tstart(ibuff))then
                tstart(i)=tstart(i)-tfield(ibuff)
            end if
20          continue

!           set pointers for this unit as undefined
        tstart(ibuff)=-1
        dtstrt(ibuff)=-2000000000
    end if
end
subroutine ftg2db(ounit,group,nulval,dim1,nx,ny, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG2DB reads a 2-d image of byte values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  c*1  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   c*1  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    character array(dim1,*),nulval
    logical anyflg,ltemp
    integer fpixel,row

    anyflg=.false.
    fpixel=1
    do 10 row = 1,ny
            call ftgpvb(ounit,group,fpixel,nx,nulval, &
                array(1,row),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue

end
subroutine ftg2dd(ounit,group,nulval,dim1,nx,ny, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG2DD reads a 2-d image of r*8 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  d  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   d  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    double precision array(dim1,*),nulval
    logical anyflg,ltemp
    integer fpixel,row

    anyflg=.false.
    fpixel=1
    do 10 row = 1,ny
            call ftgpvd(ounit,group,fpixel,nx,nulval, &
                array(1,row),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue

end
subroutine ftg2de(ounit,group,nulval,dim1,nx,ny, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG2DE reads a 2-d image of real values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  r  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   r  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    real array(dim1,*),nulval
    logical anyflg,ltemp
    integer fpixel,row

    anyflg=.false.
    fpixel=1
    do 10 row = 1,ny
            call ftgpve(ounit,group,fpixel,nx,nulval, &
                array(1,row),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue

end
subroutine ftg2di(ounit,group,nulval,dim1,nx,ny, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG2DI reads a 2-d image of i*2 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  i*2  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   i*2  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    integer*2 array(dim1,*),nulval
    logical anyflg,ltemp
    integer fpixel,row

    anyflg=.false.
    fpixel=1
    do 10 row = 1,ny
            call ftgpvi(ounit,group,fpixel,nx,nulval, &
                array(1,row),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue

end
subroutine ftg2dj(ounit,group,nulval,dim1,nx,ny, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG2DJ reads a 2-d image of i*4 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  i  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   i  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    integer array(dim1,*),nulval
    logical anyflg,ltemp
    integer fpixel,row

    anyflg=.false.
    fpixel=1
    do 10 row = 1,ny
            call ftgpvj(ounit,group,fpixel,nx,nulval, &
                array(1,row),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue

end
subroutine ftg3db(ounit,group,nulval,dim1,dim2,nx,ny,nz, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG3DB reads a 3-d cube of byte values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  c*1  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   c*1  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    character array(dim1,dim2,*),nulval
    logical anyflg,ltemp
    integer fpixel,row,band

    anyflg=.false.
    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
            call ftgpvb(ounit,group,fpixel,nx,nulval, &
                array(1,row,band),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue
20      continue
end
subroutine ftg3dd(ounit,group,nulval,dim1,dim2,nx,ny,nz, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG3DD reads a 3-d cube of byte values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  d  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   d  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    double precision array(dim1,dim2,*),nulval
    logical anyflg,ltemp
    integer fpixel,row,band

    anyflg=.false.
    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
            call ftgpvd(ounit,group,fpixel,nx,nulval, &
                array(1,row,band),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue
20      continue
end
subroutine ftg3de(ounit,group,nulval,dim1,dim2,nx,ny,nz, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG3DE reads a 3-d cube of real values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  r  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   r  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    real array(dim1,dim2,*),nulval
    logical anyflg,ltemp
    integer fpixel,row,band

    anyflg=.false.
    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
            call ftgpve(ounit,group,fpixel,nx,nulval, &
                array(1,row,band),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue
20      continue
end
subroutine ftg3di(ounit,group,nulval,dim1,dim2,nx,ny,nz, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG3DI reads a 3-d cube of i*2 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  i*2  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   i*2  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    integer*2 array(dim1,dim2,*),nulval
    logical anyflg,ltemp
    integer fpixel,row,band

    anyflg=.false.
    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
            call ftgpvi(ounit,group,fpixel,nx,nulval, &
                array(1,row,band),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue
20      continue
end
subroutine ftg3dj(ounit,group,nulval,dim1,dim2,nx,ny,nz, &
                      array,anyflg,status)
!
!*******************************************************************************
!
!! FTG3DJ reads a 3-d cube of byte values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       nulval  i  undefined pixels will be set to this value (unless = 0)
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   i  the array of values to be read
!       anyflg  l  set to true if any of the image pixels were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    integer array(dim1,dim2,*),nulval
    logical anyflg,ltemp
    integer fpixel,row,band

    anyflg=.false.
    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
            call ftgpvj(ounit,group,fpixel,nx,nulval, &
                array(1,row,band),ltemp,status)
            if (ltemp)anyflg=.true.
            fpixel=fpixel+nx
10      continue
20      continue
end
subroutine ftgabc(nfield,tform,space, rowlen,tbcol,status)
!
!*******************************************************************************
!
!! FTGABC "Gets ASCII table Beginning Columns".
!
!       determine the byte offset of the beginning of each field of a
!       ASCII table, and the total width of the table

!       nfield i  number of fields in the binary table
!       tform  c  array of FITS datatype codes of each column.
!                 must be left justified in the string variable
!       space  i  number of blank spaces to insert between each column
!       OUTPUT PARAMETERS:
!       rowlen i  total width of the table, in bytes
!       tbcol  i  beginning position of each column (first column begins at 1)
!       status i  returned error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1992

    integer nfield,space,rowlen,tbcol(*),status
    character ( len = * ) tform(*)
    integer i,j,ival

    if (status > 0)return

    rowlen=0
    do 100 i=1,nfield
            if (tform(i)(2:2) == ' ')then
!                       no explicit width; assume width=1
                    ival=1
            else
!                       find the field width characters
                    j=2
10                      j=j+1
                    if (tform(i)(j:j) == ' ' .or. &
                        tform(i)(j:j) == '.')then
!                           read the width
                        call ftc2ii(tform(i)(2:j-1),ival,status)
                    else
!                           keep looking for the end of the width field
                        go to 10
                    end if
                    tbcol(i)=rowlen+1
                    rowlen=rowlen+ival+space
            end if
100     continue

!       don't add space after the last field
    rowlen=rowlen-space
end
subroutine ftgacl(iunit,colnum,xtype,xbcol,xunit,xform, &
          xscal,xzero,xnull,xdisp,status)
!
!*******************************************************************************
!
!! FTGACL returns the parameters which define the column.
!
!       iunit   i  Fortran i/o unit number
!       colnum  i  number of the column (first column = 1)
!       xtype   c  name of the column
!       xbcol   i  starting character in the row of the column
!       xunit   c  physical units of the column
!       xform   c  Fortran-77 format of the column
!       xscal   d  scaling factor for the column values
!       xzero   d  scaling zero point for the column values
!       xnull   c  value used to represent undefined values in the column
!       xdisp   c  display format for the column (if different from xform
!       status  i  returned error status

    integer iunit,colnum,xbcol,status
    double precision xscal,xzero
    character ( len = * ) xtype,xunit,xform,xnull,xdisp

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer ibuff,nfound

    if (status > 0)return

    if (colnum < 1 .or. colnum > 999)then
!               illegal column number
            status=302
            return
    end if

    ibuff=bufnum(iunit)

!       get the parameters which are stored in the common block
    xbcol=tbcol(colnum+tstart(ibuff))+1
    xform=cform(colnum+tstart(ibuff))
    xscal=tscale(colnum+tstart(ibuff))
    xzero=tzero(colnum+tstart(ibuff))
    xnull=cnull(colnum+tstart(ibuff))

!       read remaining values from the header keywords
    xtype=' '
    call ftgkns(iunit,'TTYPE',colnum,1,xtype,nfound,status)
    xunit=' '
    call ftgkns(iunit,'TUNIT',colnum,1,xunit,nfound,status)
    xdisp=' '
    call ftgkns(iunit,'TDISP',colnum,1,xdisp,nfound,status)
end
subroutine ftgatp(ibuff,keyin,valin,status)
!
!*******************************************************************************
!
!! FTGATP "Gets ASCII Table Parameter".
!
!       test if the keyword is one of the table column definition keywords
!       of an ASCII table. If so, decode it and update the value in the common
!       block

!       ibuff   i sequence number of the data buffer
!       keynam  c name of the keyword
!       valin   c value of the keyword
!       status  i returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ibuff,status
    character ( len = * ) keyin,valin

!
!       nb = number of file buffers = max. number of FITS file opened at once
!       nf = maximum number of fields allowed in a table
    integer nf,nb
    parameter (nb = 20)
    parameter (nf = 3000)

!       tfield = number of fields in the table
!       tbcol = byte offset in the row of the beginning of the column
!       rowlen = length of one row of the table, in bytes
!       tdtype = integer code representing the datatype of the column
!       trept = the repeat count = number of data values/element in the column
!       tnull = the value used to represent an undefined value in the column
!       tscale = the scale factor for the column
!       tzero = the scaling zero point for the column
!       heapsz = the total size of the binary table heap (+ gap if any)
!       theap = the starting byte offset for the binary table heap, relative
!               to the start of the binary table data
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)

!       cnull = character string representing nulls in character columns
!       cform = the Fortran format of the column
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer nfield,i,c2,bcol,tstat
    character tform*16,keynam*8,value*70

    if (status > 0)return

    keynam=keyin
    value=valin
    tstat=status

    if (keynam(1:5) == 'TFORM')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TFORMn keyword
                status=tstat
            else
!                   get the TFORM character string, without quotes
                call ftc2s(value,tform,status)
                if (status > 0)return
                if  (tform(1:1) /= 'A' .and. tform(1:1) /= 'I' &
                .and. tform(1:1) /= 'F' .and. tform(1:1) /= 'E' &
                .and. tform(1:1) /= 'D')then
                    status=311
                 call ftpmsg('Illegal '//keynam//' format code: ' &
                             //tform)
                    return
                end if

                cform(nfield+tstart(ibuff))=tform
!                   set numeric data type code to indicate an ASCII table field
                tdtype(nfield+tstart(ibuff))=16
!                   set the repeat count to 1
                trept(nfield+tstart(ibuff))=1
!                   set the TNULL parameter to the width of the field:
                c2=0
                do 10 i=2,8
                    if (tform(i:i) >= '0' .and. tform(i:i) &
                       <= '9')then
                            c2=i
                    else
                            go to 20
                    end if
10                  continue
20                  continue

                if (status > 0)return
                if (c2 == 0)then
!                       no explicit field width, so assume width=1 character
                    tnull(nfield+tstart(ibuff))=1
                else
                    call ftc2ii(tform(2:c2),tnull(nfield+ &
                                tstart(ibuff)),status)
                    if (status > 0)then
!                               error parsing the TFORM value string
                            status=261
       call ftpmsg('Error parsing '//keynam//' field width: ' &
                    //tform)
                    end if
                end if
            end if
    else if (keynam(1:5) == 'TBCOL')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TBCOLn keyword
                status=tstat
            else
!                   get the beginning column number
                call ftc2ii(value,bcol,status)
                 if (status > 0)then
                    call ftpmsg('Error reading value of '//keynam &
                    //' as an integer: '//value)
                 else
                    tbcol(nfield+tstart(ibuff))=bcol-1
                 end if
            end if
    else if (keynam(1:5) == 'TSCAL')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TSCALn keyword
                status=tstat
            else
!                   get the scale factor
                call ftc2dd(value,tscale(nfield+tstart(ibuff)), &
                            status)
                if (status > 0)then
                     call ftpmsg('Error reading value of'//keynam &
                  //' as a Double: '//value)
                end if
            end if
    else if (keynam(1:5) == 'TZERO')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TZEROn keyword
                status=tstat
            else
!                   get the scaling zero point
                call ftc2dd(value,tzero(nfield+tstart(ibuff)), &
                            status)
                if (status > 0)then
                     call ftpmsg('Error reading value of'//keynam &
                  //' as a Double: '//value)
                end if
            end if
    else if (keynam(1:5) == 'TNULL')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TNULLn keyword
                status=tstat
            else
!                   get the Null value flag (character)
                call ftc2s(value,cnull(nfield+tstart(ibuff)),status)
                if (status > 0)then
                     call ftpmsg('Error reading value of'//keynam &
                  //' as a character string: '//value)
                end if
            end if
    end if
end
subroutine ftgbcl(iunit,colnum,xtype,xunit,dtype,rcount, &
          xscal,xzero,xnull,xdisp,status)
!
!*******************************************************************************
!
!! FTGBCL returns the parameters which define the binary column.
!
!       iunit   i  Fortran i/o unit number
!       colnum  i  number of the column (first column = 1)
!       xtype   c  name of the column
!       xunit   c  physical units of the column
!       dtype   c  datatype of the column
!       rcount  i  repeat count of the column
!       xscal   d  scaling factor for the column values
!       xzero   d  scaling zero point for the column values
!       xnull   i  value used to represent undefined values in integer column
!       xdisp   c  display format for the column
!       status  i  returned error status

    integer iunit,colnum,rcount,xnull,status
    double precision xscal,xzero
    character ( len = * ) xtype,xunit,dtype,xdisp

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,nfound,tcode
    logical descrp
    character ctemp*2,fwide*4

    if (status > 0)return

    if (colnum < 1 .or. colnum > 999)then
!               illegal column number
            status=302
            return
    end if

    ibuff=bufnum(iunit)

!       get the parameters which are stored in the common block
    rcount=trept(colnum+tstart(ibuff))
    xscal=tscale(colnum+tstart(ibuff))
    xzero=tzero(colnum+tstart(ibuff))
    xnull=tnull(colnum+tstart(ibuff))

!       translate the numeric data type code
    dtype=' '
    tcode=tdtype(colnum+tstart(ibuff))
    if (tcode < 0)then
            descrp=.true.
            tcode=-tcode
    else
            descrp=.false.
    end if

    if (tcode == 21)then
            dtype='I'
    else if (tcode == 41)then
            dtype='J'
    else if (tcode == 42)then
            dtype='E'
    else if (tcode == 82)then
            dtype='D'
    else if (tcode == 16)then
!               this is an ASCII field; width of field is stored in TNULL
            write(fwide,1000)tnull(colnum+tstart(ibuff))
1000            format(i4)
            if (tnull(colnum+tstart(ibuff)) > 999)then
                dtype='A'//fwide
            else if (tnull(colnum+tstart(ibuff)) > 99)then
                dtype='A'//fwide(2:4)
            else if (tnull(colnum+tstart(ibuff)) > 9)then
                dtype='A'//fwide(3:4)
            else if (tnull(colnum+tstart(ibuff)) > 0)then
                dtype='A'//fwide(4:4)
            else
                dtype='A'
            end if
!               ASCII column don't have an integer null value
            xnull=0
    else if (tcode == 14)then
            dtype='L'
    else if (tcode == 1)then
            dtype='X'
    else if (tcode == 11)then
            dtype='B'
    else if (tcode == 83)then
            dtype='C'
    else if (tcode == 163)then
            dtype='M'
    end if

    if (descrp)then
            ctemp='P'//dtype(1:1)
            dtype=ctemp
    end if

!       read remaining values from the header keywords
    xtype=' '
    call ftgkns(iunit,'TTYPE',colnum,1,xtype,nfound,status)
    xunit=' '
    call ftgkns(iunit,'TUNIT',colnum,1,xunit,nfound,status)
    xdisp=' '
    call ftgkns(iunit,'TDISP',colnum,1,xdisp,nfound,status)
end
subroutine ftgbit(buffer,log8)
!
!*******************************************************************************
!
!! FTGBIT decodes bits within the byte into an array of logical values.
!
!       The corresponding logical value is set to
!       true if the bit is set to 1.

!       buffer  i  input integer containing the byte to be decoded
!       log8    l  output array of logical data values corresponding
!                  to the bits in the input buffer
!
!       written by Wm Pence, HEASARC/GSFC, May 1992

    integer buffer,tbuff
    logical log8(8)

    log8(1)=.false.
    log8(2)=.false.
    log8(3)=.false.
    log8(4)=.false.
    log8(5)=.false.
    log8(6)=.false.
    log8(7)=.false.
    log8(8)=.false.

!       test for special case: no bits are set
    if (buffer == 0)return

!       This algorithm tests to see if each bit is set by testing
!       the numerical value of the byte, starting with the most significant
!       bit.  If the bit is set, then it is reset to zero before testing
!       the next most significant bit, and so on.

    tbuff=buffer

!       now decode the least significant byte
    if (tbuff > 127)then
            log8(1)=.true.
            tbuff=tbuff-128
    end if
    if (tbuff > 63)then
            log8(2)=.true.
            tbuff=tbuff-64
    end if
    if (tbuff > 31)then
            log8(3)=.true.
            tbuff=tbuff-32
    end if
    if (tbuff > 15)then
            log8(4)=.true.
            tbuff=tbuff-16
    end if
    if (tbuff > 7)then
            log8(5)=.true.
            tbuff=tbuff-8
    end if
    if (tbuff > 3)then
            log8(6)=.true.
            tbuff=tbuff-4
    end if
    if (tbuff > 1)then
            log8(7)=.true.
            tbuff=tbuff-2
    end if
    if (tbuff == 1)then
            log8(8)=.true.
    end if
end
subroutine ftgbnh(iunit,nrows,nfield,ttype,tform,tunit, &
                      extnam,pcount,status)
!
!*******************************************************************************
!
!! FTGBNH is obsolete.  Call FTGHBN instead.
!

    integer iunit,nrows,nfield,pcount,status
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam

    call ftghbn(iunit,-1,nrows,nfield,ttype,tform, &
                      tunit,extnam,pcount,status)
end
subroutine ftgbtp(ibuff,keyin,valin,status)
!
!*******************************************************************************
!
!! FTGBTP "Gets Binary Table Parameter"
!
!       test if the keyword is one of the table column definition keywords
!       of a binary table. If so, decode it and update the values in the common
!       block

!       ibuff   i sequence number of the data buffer
!       keynam  c name of the keyword
!       valout  c value of the keyword
!       OUTPUT PARAMETERS:
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ibuff,status,width
    character ( len = * ) keyin,valin

!
!       nb = number of file buffers = max. number of FITS file opened at once
!       nf = maximum number of fields allowed in a table
    integer nf,nb
    parameter (nb = 20)
    parameter (nf = 3000)
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer nfield,tstat
    character tform*16,keynam*8,value*70

    if (status > 0)return
    keynam=keyin
    value=valin
    tstat=status

    if (keynam(1:5) == 'TFORM')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TFORMn keyword
                status=tstat
            else
!                   get the TFORM character string, without quotes
                call ftc2s(value,tform,status)
!                   get the datatype code and repeat count
                call ftbnfm(tform,tdtype(nfield+tstart(ibuff)), &
                   trept(nfield+tstart(ibuff)),width,status)
                if (tdtype(nfield+tstart(ibuff)) == 1)then
!                       treat Bit datatype as if it were a Byte datatype
                    tdtype(nfield+tstart(ibuff))=11
                    trept(nfield+tstart(ibuff))=(trept(nfield+ &
                    tstart(ibuff))+7)/8
                else if (tdtype(nfield+tstart(ibuff)) == 16)then
!                      store the width of the ASCII field in the TNULL parameter
                    tnull(nfield+tstart(ibuff))=width
               end if
            end if
    else if (keynam(1:5) == 'TSCAL')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TSCALn keyword
                status=tstat
            else
!                   get the scale factor
                call ftc2dd(value,tscale(nfield+tstart(ibuff)), &
                            status)
                if (status > 0)then
                     call ftpmsg('Error reading value of'//keynam &
                  //' as a Double: '//value)
                end if
            end if
    else if (keynam(1:5) == 'TZERO')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TZEROn keyword
                status=tstat
            else
!                   get the scaling zero point
                call ftc2dd(value,tzero(nfield+tstart(ibuff)), &
                            status)
                if (status > 0)then
                     call ftpmsg('Error reading value of'//keynam &
                  //' as a Double: '//value)
                end if
            end if
    else if (keynam(1:5) == 'TNULL')then
!               get the field number
            call ftc2ii(keynam(6:8),nfield,status)
            if (status > 0)then
!                   this must not have been a TNULLn keyword
                status=tstat
            else
!                   make sure this is not an ASCII column (the tnull
!                   variable is use to store the ASCII column width)
                if (tdtype(nfield+tstart(ibuff)) /= 16)then
!                       get the Null value flag (Integer)
                    call ftc2ii(value,tnull(nfield+tstart(ibuff)), &
                                status)
                    if (status > 0)then
                        call ftpmsg('Error reading value of '// &
                        keynam//' as an integer: '//value)
                    end if
                end if
            end if
    else if (keynam(1:8) == 'THEAP   ')then
!               get the heap offset value
            call ftc2ii(value,theap(ibuff),status)
            if (status > 0)then
                    call ftpmsg('Error reading value of '//keynam &
                    //' as an integer: '//value)
            end if
    end if
end
subroutine ftgcfb(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFB reads an array of byte values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   b  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    character array(*),dummy
    integer i

    do i=1,nelem
            flgval(i)=.false.
    end do

    call ftgclb(iunit,colnum,frow,felem,nelem,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfc(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFC reads an array of complex values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   cmp  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    real array(*),dummy
    integer i
    integer felemx, nelemx

!       a complex value is interpreted as a pair of float values, thus
!       need to multiply the first element and number of elements by 2

    felemx = (felem - 1) * 2 + 1
    nelemx = nelem * 2

    do 10 i=1,nelemx
            flgval(i)=.false.
10      continue

    call ftgcle(iunit,colnum,frow,felemx,nelemx,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfd(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFD reads an array of r*8 values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   d  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    double precision array(*),dummy
    integer i

    do i=1,nelem
            flgval(i)=.false.
    end do

    call ftgcld(iunit,colnum,frow,felem,nelem,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfe(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFE reads an array of R*4 values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   r  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    real array(*),dummy
    integer i

    do i=1,nelem
            flgval(i)=.false.
    end do

    call ftgcle(iunit,colnum,frow,felem,nelem,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfi(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFI reads an array of I*2 values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   i*2 returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    integer*2 array(*),dummy
    integer i

    do i=1,nelem
            flgval(i)=.false.
    end do

    call ftgcli(iunit,colnum,frow,felem,nelem,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfj(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFJ reads an array of I*4 values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   i  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    integer array(*),dummy,i

    do i=1,nelem
            flgval(i)=.false.
    end do

    call ftgclj(iunit,colnum,frow,felem,nelem,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfl(iunit,colnum,frow,felem,nelem,lray, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFL reads logical values from a specified column of the table.
!
!       The binary table column being read from must have datatype 'L'
!       and no datatype conversion will be perform if it is not.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       lray    l  returned array of data values that is read
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical lray(*),flgval(*),anynul

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer bstart,maxpix,tcode,offset
    integer ibuff,i,i1,ntodo,itodo,repeat,rstart,estart
    character buffer(80)
    logical descrp
    character messge*80

    if (status > 0)return

    ibuff=bufnum(iunit)
    tcode=tdtype(colnum+tstart(ibuff))

!       Do sanity check of input parameters
    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    else if (nelem == 0)then
      return
    end if

!       initialize the null flag array
    do 5 i=1,nelem
            flgval(i)=.false.
5       continue
    anynul=.false.

    i1=0
    ntodo=nelem
    rstart=frow-1
    estart=felem-1
    maxpix=80

    if (tcode == 14)then
            repeat=trept(colnum+tstart(ibuff))
            if (felem > repeat)then
!                   illegal element number
                write(messge,1005)felem
1005                format( &
         'Starting element number is greater than repeat: ',i10)
                call ftpmsg(messge)
                status = 308
                return
            end if
            descrp=.false.
    else if (tcode == -14)then
!               this is a variable length descriptor column
            descrp=.true.
!               read the number of elements and the starting offset:
            call ftgdes(iunit,colnum,frow,repeat, &
                                offset,status)
            if (repeat == 0)then
!                       error: null length vector
                    status=318
                    return
            else if (estart+ntodo > repeat)then
!                       error: trying to read beyond end of record
                    status=319
                    return
            end if
!               move the i/o pointer to the start of the pixel sequence
            bstart=dtstrt(ibuff)+offset+ &
                            theap(ibuff)+estart
            call ftmbyt(iunit,bstart,.true.,status)
    else
!               column must be logical data type
            status=312
            return
    end if

!       process as many contiguous pixels as possible
20      itodo=min(ntodo,repeat-estart,maxpix)

    if (.not. descrp)then
!           move the i/o pointer to the start of the sequence of pixels
        bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
               tbcol(colnum+tstart(ibuff))+estart
        call ftmbyt(iunit,bstart,.false.,status)
    end if

!       get the array of logical bytes
    call ftgcbf(iunit,itodo,buffer,status)
    if (status > 0)return

!       decode the 'T' and 'F' characters, and look for nulls (0)
    do 10 i=1,itodo
            if (buffer(i) == 'T')then
                    lray(i1+i)=.true.
            else if (buffer(i) == 'F')then
                    lray(i1+i)=.false.
            else if (ichar(buffer(i)) == 0)then
                    flgval(i1+i)=.true.
                    anynul=.true.
            else
                    status=316
                    return
            end if
10      continue

    if (status > 0)then
        write(messge,1006)i1+1,i1+itodo
1006        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCFL).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                    estart=0
                    rstart=rstart+1
            end if
            go to 20
    end if
end
subroutine ftgcfm(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFM reads an double precision complex values from a column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       array   dcmp  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    double precision array(*),dummy
    integer i
    integer felemx, nelemx

!       a complex value is interpreted as a pair of float values, thus
!       need to multiply the first element and number of elements by 2

    felemx = (felem - 1) * 2 + 1
    nelemx = nelem * 2

    do 10 i=1,nelemx
            flgval(i)=.false.
10      continue

    call ftgcld(iunit,colnum,frow,felemx,nelemx,1,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcfs(iunit,colnum,frow,felem,nelem,array, &
            flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCFS reads an array of string values from a specified column of the table.
!
!       Any undefined pixels will be have the corresponding value of FLGVAL
!       set equal to .true., and ANYNUL will be set equal to .true. if
!       any pixels are undefined.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element in the row to read
!       nelem   i  number of elements to read
!       array   c  returned array of data values that was read from FITS file
!       flgval  l  set .true. if corresponding element undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval(*),anynul
    character ( len = * ) array(*)
    character*8 dummy
    integer i

    do 10 i=1,nelem
            flgval(i)=.false.
10      continue

    call ftgcls(iunit,colnum,frow,felem,nelem,2,dummy, &
        array,flgval,anynul,status)
end
subroutine ftgcks(iunit,datsum,chksum,status)
!
!*******************************************************************************
!
!! FTGCKS calculates and encodes the checksums of the data unit and the total HDU

!       iunit   i  fortran unit number
!       datsum  d  output  checksum for the data
!       chksum  d  output  checksum for the entire HDU
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Sept, 1994

    integer iunit,status
    double precision datsum,chksum

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

    integer ibuff,nrec

    if (status > 0)return

!       calculate number of data records
    ibuff=bufnum(iunit)
    nrec=(hdstrt(ibuff,chdu(ibuff)+1)-dtstrt(ibuff))/2880

    datsum=0.
    if (nrec > 0)then

!           move to the start of the data
        call ftmbyt(iunit,dtstrt(ibuff),.true.,status)

!           accumulate the 32-bit 1's complement checksum
        call ftcsum(iunit,nrec,datsum,status)
    end if

!       move to the start of the header
    call ftmbyt(iunit,hdstrt(ibuff,chdu(ibuff)),.true.,status)

!       calculate number of FITS blocks in the header
    nrec=(dtstrt(ibuff)-hdstrt(ibuff,chdu(ibuff)))/2880

!       accumulate the header into the checksum
    chksum=datsum
    call ftcsum(iunit,nrec,chksum,status)
end
subroutine ftgcl(iunit,colnum,frow,felem,nelem,lray,status)
!
!*******************************************************************************
!
!! FTGCL reads an array of logical values from a specified column of the table.
!
!       The binary table column being read from must have datatype 'L'
!       and no datatype conversion will be perform if it is not.
!       This routine ignores any undefined values in the logical array.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       lray    l  returned array of data values that is read
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical lray(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character buffer(32000)
    common/ftheap/buffer
!

    integer bstart,maxpix,offset,tcode
    integer ibuff,i,i1,ntodo,itodo,repeat,rstart,estart
    logical descrp
    character messge*80

    if (status > 0)return

    ibuff=bufnum(iunit)
    tcode=tdtype(colnum+tstart(ibuff))

!       Do sanity check of input parameters
    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    else if (nelem == 0)then
      return
    end if

    i1=0
    ntodo=nelem
    rstart=frow-1
    estart=felem-1
    maxpix=32000

    if (tcode == 14)then
            repeat=trept(colnum+tstart(ibuff))
            if (felem > repeat)then
!                   illegal element number
                write(messge,1005)felem
1005                format( &
         'Starting element number is greater than repeat: ',i10)
                call ftpmsg(messge)
                status = 308
                return
            end if
            descrp=.false.
    else if (tcode == -14)then
!               this is a variable length descriptor column
            descrp=.true.
!               read the number of elements and the starting offset:
            call ftgdes(iunit,colnum,frow,repeat, &
                                offset,status)
            if (repeat == 0)then
!                       error: null length vector
                    status=318
                    return
            else if (estart+ntodo > repeat)then
!                       error: trying to read beyond end of record
                    status=319
                    return
            end if
!               move the i/o pointer to the start of the pixel sequence
            bstart=dtstrt(ibuff)+offset+ &
                            theap(ibuff)+estart
            call ftmbyt(iunit,bstart,.true.,status)
    else
!               column must be logical data type
            status=312
            return
    end if

!       process as many contiguous pixels as possible
20      itodo=min(ntodo,repeat-estart,maxpix)

    if (.not. descrp)then
!           move the i/o pointer to the start of the sequence of pixels
        bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
        tbcol(colnum+tstart(ibuff))+estart
        call ftmbyt(iunit,bstart,.false.,status)
    end if

!       get the array of logical bytes
    call ftgcbf(iunit,itodo,buffer,status)

!       decode the 'T' and 'F' characters,
    do 10 i=1,itodo
            if (buffer(i) == 'T')then
                    lray(i1+i)=.true.
            else if (buffer(i) == 'F')then
                    lray(i1+i)=.false.
            else if (ichar(buffer(i)) == 0)then
!                       ignore null values; leave input logical value unchanged
            else
!                       illegal logical value
                    status=316
                    return
            end if
10      continue

    if (status > 0)then
        write(messge,1006)i1+1,i1+itodo
1006        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCL).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                    estart=0
                    rstart=rstart+1
            end if
            go to 20
    end if
end
subroutine ftgclb(iunit,colnum,frow,felem,nelem,eincr, &
     nultyp,nulval,array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCLB reads byte data values from the specified column of the table.
!
!       This general purpose routine will handle null values in one
!       of two ways: if nultyp=1, then undefined array elements will be
!       set equal to the input value of NULVAL.  Else if nultyp=2, then
!       undefined array elements will have the corresponding FLGVAL element
!       set equal to .TRUE.  If NULTYP=1 and NULVAL=0, then no checks for
!       undefined values will be made, for maximum efficiency.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       eincr   i  element increment
!       nultyp  i  input code indicating how to handle undefined values
!       nulval  b  value that undefined pixels will be set to (if nultyp=1)
!       array   b  array of data values that are read from the FITS file
!       flgval  l  set .true. if corresponding element undefined (if nultyp=2)
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,eincr,nultyp,status
    character array(*),nulval
    logical flgval(*),anynul

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer nulchk,i4null,rskip
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    logical tofits,trans
    integer*2 i2null
    character sval*30,sform*13,snull*16,i1null*1,messge*80
    integer buffer(8000)
    common/fttemp/buffer

    if (status > 0)return

    call ftgcpr(iunit,colnum,frow,felem,nelem,0, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,i4null,snull,status)

    if (status > 0 .or. nelem == 0)return

!       multiply incre to just get every nth pixel
    incre = incre * eincr

!       determine if we have to check for null values
    nulchk = nultyp
    if (nultyp == 1 .and. ichar(nulval) == 0)then
!           user doesn't want to check for nulls
        nulchk=0
    else
!           user does want to check for null values
        if (tcode <= 41)then
!               check if null value is defined for integer column
            if (i4null == 123454321)then
                nulchk=0
            else
                if (tcode == 11)then
                    i1null=char(i4null)
                else if (tcode == 21)then
                    i2null=i4null
                end if
            end if
        end if
    end if

!       check for important special case: no datatype conversion required
    if (tcode == 11 .and. nulchk == 0 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
          trans=.false.
    else
          trans=.true.
    end if

    sval=' '
    i1=1
    ntodo=nelem
    rstart=0
    anynul=.false.
!       the data are being scaled from FITS to internal format
    tofits=.false.

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,(repeat-estart-1)/eincr+1,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp+(rstart * lenrow) + (estart * incre / eincr)
    call ftmbyt(iunit,bstart,.false.,status)

!       read the data from FITS file, doing datatype conversion and scaling
    if (tcode == 21)then
!               column data type is I (I*2)
!               read the data and do any machine dependent data conversion
            call ftgi2b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti2i1(buffer,itodo,scale,zero,tofits, &
            nulchk,i2null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 41)then
!               column data type is J (I*4)
!               read the data and do any machine dependent data conversion
            call ftgi4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti4i1(buffer,itodo,scale,zero,tofits, &
            nulchk,i4null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 42)then
!               column data type is E (R*4)
!               read the data and do any machine dependent data conversion
            call ftgr4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr4i1(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 82)then
!               column data type is D (R*8)
!               read the data and do any machine dependent data conversion
            call ftgr8b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr8i1(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 11)then
!               column data type is B (byte)
!               read the data and do any machine dependent data conversion
!               note that we can use the input array directly
            call ftgi1b(iunit,itodo,incre,array(i1),status)
!               check for null values, and do scaling and datatype conversion
            if (trans)then
              call fti1i1(array(i1),itodo,scale,zero,tofits,nulchk, &
              i1null,nulval,flgval(i1),anynul,array(i1),status)
            end if
    else
!               this is an ASCII table column; get the character string
            call ftgcbf(iunit,twidth,sval,status)
            if (status > 0)return

!               check for null value
            if (sval(1:16) == snull)then
                anynul=.true.
                if (nultyp == 1)then
                    array(i1)=nulval
                else if (nultyp == 2)then
                    flgval(i1)=.true.
                end if
            else
!                   read the value, then do scaling and datatype conversion
                if (sform(5:5) == 'I')then
                    read(sval,sform,err=900)ival
                    dval=ival*scale+zero
                else
                    read(sval,sform,err=900)dval
                    dval=dval*scale+zero
                end if

!                   trap any values that overflow the I*1 range
                if (dval < 255.49 .and. dval > -.49)then
                    array(i1)=char(int(dval))
                else if (dval >= 255.49)then
                    status=-11
                    array(i1)=char(255)
                else
                    status=-11
                    array(i1)=char(0)
                end if
            end if
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCLB).')
        call ftpmsg(messge)
        return
    end if

    if (ntodo > 0)then
!           increment the pointers
        i1=i1+itodo
        estart=estart+itodo*eincr
        if (estart >= repeat)then
            rskip=estart/repeat
            rstart=rstart+rskip
            estart=estart-rskip*repeat
        end if
        go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while reading FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error reading formatted data value from ASCII table
    write(messge,1002)colnum,rstart+frow
1002    format('Error reading column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to read value with format '//sform)
    status=315
end
subroutine ftgcld(iunit,colnum,frow,felem,nelem,eincr, &
     nultyp,nulval,array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCLD reads real*8 data values from the specified column of the table.
!
!       This general purpose routine will handle null values in one
!       of two ways: if nultyp=1, then undefined array elements will be
!       set equal to the input value of NULVAL.  Else if nultyp=2, then
!       undefined array elements will have the corresponding FLGVAL element
!       set equal to .TRUE.  If NULTYP=1 and NULVAL=0, then no checks for
!       undefined values will be made, for maximum efficiency.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       eincr   i  element increment
!       nultyp  i  input code indicating how to handle undefined values
!       nulval  d  value that undefined pixels will be set to (if nultyp=1)
!       array   d  array of data values that are read from the FITS file
!       flgval  l  set .true. if corresponding element undefined (if nultyp=2)
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,eincr,nultyp,status
    double precision array(*),nulval
    logical flgval(*),anynul

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer nulchk,i4null,rskip
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    logical tofits,trans
    integer*2 i2null
    character sval*30,sform*13,snull*16,i1null*1,messge*80
    character chbuff(32000)
    common/ftheap/chbuff
    integer buffer(8000)
    common/fttemp/buffer

    if (status > 0)return

    call ftgcpr(iunit,colnum,frow,felem,nelem,0, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,i4null,snull,status)

    if (status > 0 .or. nelem == 0)return

!       multiply incre to just get every nth pixel
    incre = incre * eincr

!       determine if we have to check for null values
    nulchk = nultyp
    if (nultyp == 1 .and. nulval == 0)then
!           user doesn't want to check for nulls
        nulchk=0
    else
!           user does want to check for null values
        if (tcode <= 41)then
!               check if null value is defined for integer column
            if (i4null == 123454321)then
                nulchk=0
            else
                if (tcode == 11)then
                    i1null=char(i4null)
                else if (tcode == 21)then
                    i2null=i4null
                end if
            end if
        end if
    end if

!       check for important special case: no datatype conversion required
    if (tcode == 82 .and. nulchk == 0 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
          trans=.false.
    else
          trans=.true.
    end if

    sval=' '
    i1=1
    ntodo=nelem
    rstart=0
    anynul=.false.
!       the data are being scaled from FITS to internal format
    tofits=.false.

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,(repeat-estart-1)/eincr+1,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp+(rstart * lenrow) + (estart * incre / eincr)
    call ftmbyt(iunit,bstart,.false.,status)

!       read the data from FITS file, doing datatype conversion and scaling
    if (tcode == 21)then
!               column data type is I (I*2)
!               read the data and do any machine dependent data conversion
            call ftgi2b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti2r8(buffer,itodo,scale,zero,tofits, &
            nulchk,i2null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 41)then
!               column data type is J (I*4)
!               read the data and do any machine dependent data conversion
            call ftgi4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti4r8(buffer,itodo,scale,zero,tofits, &
            nulchk,i4null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 42)then
!               column data type is E (R*4)
!               read the data and do any machine dependent data conversion
            call ftgr4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr4r8(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 82)then
!               column data type is D (R*8)
!               read the data and do any machine dependent data conversion
!               note that we can use the input array directly
            call ftgr8b(iunit,itodo,incre,array(i1),status)
!               check for null values, and do scaling and datatype conversion
            if (trans)then
              call ftr8r8(array(i1),itodo,scale,zero,tofits, &
              nulchk,nulval,flgval(i1),anynul,array(i1),status)
            end if
    else if (tcode == 11)then
!               column data type is B (byte)
!               read the data and do any machine dependent data conversion
            call ftgi1b(iunit,itodo,incre,chbuff,status)
!               check for null values, and do scaling and datatype conversion
            call fti1r8(chbuff,itodo,scale,zero,tofits, &
            nulchk,i1null,nulval,flgval(i1),anynul,array(i1),status)
    else
!               this is an ASCII table column; get the character string
            call ftgcbf(iunit,twidth,sval,status)
            if (status > 0)return

!               check for null
            if (sval(1:16) == snull)then
                    anynul=.true.
                    if (nultyp == 1)then
                            array(i1)=nulval
                    else if (nultyp == 2)then
                            flgval(i1)=.true.
                    end if

!               now read the value, then do scaling and datatype conversion
            else if (sform(5:5) == 'I')then
                    read(sval,sform,err=900)ival
                    array(i1)=ival*scale+zero
            else
                    read(sval,sform,err=900)dval
                    array(i1)=dval*scale+zero
            end if
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCLD).')
        call ftpmsg(messge)
        return
    end if

    if (ntodo > 0)then
!           increment the pointers
        i1=i1+itodo
        estart=estart+itodo*eincr
        if (estart >= repeat)then
            rskip=estart/repeat
            rstart=rstart+rskip
            estart=estart-rskip*repeat
        end if
        go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while reading FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error reading formatted data value from ASCII table
    write(messge,1002)colnum,rstart+frow
1002    format('Error reading column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to read value with format '//sform)
    status=315
end
subroutine ftgcle(iunit,colnum,frow,felem,nelem,eincr, &
     nultyp,nulval,array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCLE reads real*4 data values from the specified column of the table.
!
!       This general purpose routine will handle null values in one
!       of two ways: if nultyp=1, then undefined array elements will be
!       set equal to the input value of NULVAL.  Else if nultyp=2, then
!       undefined array elements will have the corresponding FLGVAL element
!       set equal to .TRUE.  If NULTYP=1 and NULVAL=0, then no checks for
!       undefined values will be made, for maximum efficiency.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       eincr   i  element increment
!       nultyp  i  input code indicating how to handle undefined values
!       nulval  r  value that undefined pixels will be set to (if nultyp=1)
!       array   r  array of data values that are read from the FITS file
!       flgval  l  set .true. if corresponding element undefined (if nultyp=2)
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,eincr,nultyp,status
    real array(*),nulval
    logical flgval(*),anynul

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer nulchk,i4null,rskip
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    logical tofits,trans
    integer*2 i2null
    character sval*30,sform*13,snull*16,i1null*1,messge*80
    character chbuff(32000)
    common/ftheap/chbuff
    integer buffer(8000)
    common/fttemp/buffer

    if (status > 0)return

    call ftgcpr(iunit,colnum,frow,felem,nelem,0, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,i4null,snull,status)

    if (status > 0 .or. nelem == 0)return

!       multiply incre to just get every nth pixel
    incre = incre * eincr

!       determine if we have to check for null values
    nulchk = nultyp
    if (nultyp == 1 .and. nulval == 0)then
!           user doesn't want to check for nulls
        nulchk=0
    else
!           user does want to check for null values
        if (tcode <= 41)then
!               check if null value is defined for integer column
            if (i4null == 123454321)then
                nulchk=0
            else
                if (tcode == 11)then
                    i1null=char(i4null)
                else if (tcode == 21)then
                    i2null=i4null
                end if
            end if
        end if
    end if

!       check for important special case: no datatype conversion required
    if (tcode == 42 .and. nulchk == 0 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
          trans=.false.
    else
          trans=.true.
    end if

    sval=' '
    i1=1
    ntodo=nelem
    rstart=0
    anynul=.false.
!       the data are being scaled from FITS to internal format
    tofits=.false.

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,(repeat-estart-1)/eincr+1,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp+(rstart * lenrow) + (estart * incre / eincr)
    call ftmbyt(iunit,bstart,.false.,status)

!       read the data from FITS file, doing datatype conversion and scaling
    if (tcode == 42)then
!               column data type is E (R*4)
!               read the data and do any machine dependent data conversion
!               note that we can use the input array directly
            call ftgr4b(iunit,itodo,incre,array(i1),status)
!               check for null values, and do scaling and datatype conversion
            if (trans)then
              call ftr4r4(array(i1),itodo,scale,zero,tofits,nulchk, &
              nulval,flgval(i1),anynul,array(i1),status)
            end if
    else if (tcode == 21)then
!               column data type is I (I*2)
!               read the data and do any machine dependent data conversion
            call ftgi2b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti2r4(buffer,itodo,scale,zero,tofits, &
            nulchk,i2null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 41)then
!               column data type is J (I*4)
!               read the data and do any machine dependent data conversion
            call ftgi4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti4r4(buffer,itodo,scale,zero,tofits, &
            nulchk,i4null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 82)then
!               column data type is D (R*8)
!               read the data and do any machine dependent data conversion
            call ftgr8b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr8r4(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 11)then
!               column data type is B (byte)
!               read the data and do any machine dependent data conversion
            call ftgi1b(iunit,itodo,incre,chbuff,status)
!               check for null values, and do scaling and datatype conversion
            call fti1r4(chbuff,itodo,scale,zero,tofits, &
            nulchk,i1null,nulval,flgval(i1),anynul,array(i1),status)
    else
!               this is an ASCII table column; get the character string
            call ftgcbf(iunit,twidth,sval,status)
            if (status > 0)return

!               check for null
            if (sval(1:16) == snull)then
                    anynul=.true.
                    if (nultyp == 1)then
                            array(i1)=nulval
                    else if (nultyp == 2)then
                            flgval(i1)=.true.
                    end if

!               now read the value, then do scaling and datatype conversion
            else if (sform(5:5) == 'I')then
                    read(sval,sform,err=900)ival
                    array(i1)=ival*scale+zero
            else
                    read(sval,sform,err=900)dval
                    array(i1)=dval*scale+zero
            end if
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCLE).')
        call ftpmsg(messge)
        return
    end if

    if (ntodo > 0)then
!           increment the pointers
        i1=i1+itodo
        estart=estart+itodo*eincr
        if (estart >= repeat)then
            rskip=estart/repeat
            rstart=rstart+rskip
            estart=estart-rskip*repeat
        end if
        go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while reading FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error reading formatted data value from ASCII table
    write(messge,1002)colnum,rstart+frow
1002    format('Error reading column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to read value with format '//sform)
    status=315
end
subroutine ftgcli(iunit,colnum,frow,felem,nelem,eincr, &
     nultyp,nulval,array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCLI reads integer*2 values from the specified column of the table.
!
!       This general purpose routine will handle null values in one
!       of two ways: if nultyp=1, then undefined array elements will be
!       set equal to the input value of NULVAL.  Else if nultyp=2, then
!       undefined array elements will have the corresponding FLGVAL element
!       set equal to .TRUE.  If NULTYP=1 and NULVAL=0, then no checks for
!       undefined values will be made, for maximum efficiency.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       eincr   i  element increment
!       nultyp  i  input code indicating how to handle undefined values
!       nulval  i*2  value that undefined pixels will be set to (if nultyp=1)
!       array   i*2  array of data values that are read from the FITS file
!       flgval  l  set .true. if corresponding element undefined (if nultyp=2)
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,eincr,nultyp,status
    integer*2 array(*),nulval
    logical flgval(*),anynul

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer nulchk,i4null,rskip
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    logical tofits,trans
    integer*2 i2null
    character sval*30,sform*13,snull*16,i1null*1,messge*80
    integer maxi2,mini2
    double precision i2max,i2min
    parameter (i2max=3.276749D+04)
    parameter (i2min=-3.276849D+04)
    parameter (maxi2=32767)
    parameter (mini2=-32768)
    character chbuff(32000)
    common/ftheap/chbuff
    integer buffer(8000)
    common/fttemp/buffer

    if (status > 0)return

    call ftgcpr(iunit,colnum,frow,felem,nelem,0, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,i4null,snull,status)

    if (status > 0 .or. nelem == 0)return

!       multiply incre to just get every nth pixel
    incre = incre * eincr

!       determine if we have to check for null values
    nulchk = nultyp
    if (nultyp == 1 .and. nulval == 0)then
!           user doesn't want to check for nulls
        nulchk=0
    else
!           user does want to check for null values
        if (tcode <= 41)then
!               check if null value is defined for integer column
            if (i4null == 123454321)then
                nulchk=0
            else
                if (tcode == 11)then
                    i1null=char(i4null)
                else if (tcode == 21)then
                    i2null=i4null
                end if
            end if
        end if
    end if

!       check for important special case: no datatype conversion required
    if (tcode == 21 .and. nulchk == 0 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
          trans=.false.
    else
          trans=.true.
    end if

    sval=' '
    i1=1
    ntodo=nelem
    rstart=0
    anynul=.false.
!       the data are being scaled from FITS to internal format
    tofits=.false.

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,(repeat-estart-1)/eincr+1,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp+(rstart * lenrow) + (estart * incre / eincr)
    call ftmbyt(iunit,bstart,.false.,status)

!       read the data from FITS file, doing datatype conversion and scaling
    if (tcode == 21)then
!               column data type is I (I*2)
!               read the data and do any machine dependent data conversion
!               note that we can use the input array directly
            call ftgi2b(iunit,itodo,incre,array(i1),status)
!               check for null values, and do scaling and datatype conversion
            if (trans)then
              call fti2i2(array(i1),itodo,scale,zero,tofits,nulchk, &
              i2null,nulval,flgval(i1),anynul,array(i1),status)
            end if
    else if (tcode == 41)then
!               column data type is J (I*4)
!               read the data and do any machine dependent data conversion
            call ftgi4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti4i2(buffer,itodo,scale,zero,tofits, &
            nulchk,i4null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 42)then
!               column data type is E (R*4)
!               read the data and do any machine dependent data conversion
            call ftgr4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr4i2(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 82)then
!               column data type is D (R*8)
!               read the data and do any machine dependent data conversion
            call ftgr8b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr8i2(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 11)then
!               column data type is B (byte)
!               read the data and do any machine dependent data conversion
            call ftgi1b(iunit,itodo,incre,chbuff,status)
!               check for null values, and do scaling and datatype conversion
            call fti1i2(chbuff,itodo,scale,zero,tofits, &
            nulchk,i1null,nulval,flgval(i1),anynul,array(i1),status)
    else
!               this is an ASCII table column; get the character string
            call ftgcbf(iunit,twidth,sval,status)
            if (status > 0)return

!               check for null value
            if (sval(1:16) == snull)then
                anynul=.true.
                if (nultyp == 1)then
                    array(i1)=nulval
                else if (nultyp == 2)then
                    flgval(i1)=.true.
                end if
            else
!                   read the value, then do scaling and datatype conversion
                if (sform(5:5) == 'I')then
                    read(sval,sform,err=900)ival
                    dval=ival*scale+zero
                else
                    read(sval,sform,err=900)dval
                    dval=dval*scale+zero
                end if

!                   trap any values that overflow the I*2 range
                if (dval < i2max .and. dval > i2min)then
                    array(i1)=dval
                else if (dval >= i2max)then
                    status=-11
                    array(i1)=maxi2
                else
                    status=-11
                    array(i1)=mini2
                end if
            end if
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCLI).')
        call ftpmsg(messge)
        return
    end if

    if (ntodo > 0)then
!           increment the pointers
        i1=i1+itodo
        estart=estart+itodo*eincr
        if (estart >= repeat)then
            rskip=estart/repeat
            rstart=rstart+rskip
            estart=estart-rskip*repeat
        end if
        go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while reading FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error reading formatted data value from ASCII table
    write(messge,1002)colnum,rstart+frow
1002    format('Error reading column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to read value with format '//sform)
    status=315
end
subroutine ftgclj(iunit,colnum,frow,felem,nelem,eincr, &
     nultyp,nulval,array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCLJ reads integer*4 data from the specified column of the table.
!
!       This general purpose routine will handle null values in one
!       of two ways: if nultyp=1, then undefined array elements will be
!       set equal to the input value of NULVAL.  Else if nultyp=2, then
!       undefined array elements will have the corresponding FLGVAL element
!       set equal to .TRUE.  If NULTYP=1 and NULVAL=0, then no checks for
!       undefined values will be made, for maximum efficiency.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       eincr   i  element increment
!       nultyp  i  input code indicating how to handle undefined values
!       nulval  i  value that undefined pixels will be set to (if nultyp=1)
!       array   i  array of data values that are read from the FITS file
!       flgval  l  set .true. if corresponding element undefined (if nultyp=2)
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,eincr,nultyp,status
    integer array(*),nulval
    logical flgval(*),anynul

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer nulchk,i4null,rskip
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    logical tofits,trans
    integer*2 i2null
    character sval*30,sform*13,snull*16,i1null*1,messge*80
    character chbuff(32000)
    double precision i4max,i4min
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
    common/ftheap/chbuff
    integer buffer(8000)
    common/fttemp/buffer

!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    call ftgcpr(iunit,colnum,frow,felem,nelem,0, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,i4null,snull,status)

    if (status > 0 .or. nelem == 0)return

!       multiply incre to just get every nth pixel
    incre = incre * eincr

!       determine if we have to check for null values
    nulchk = nultyp
    if (nultyp == 1 .and. nulval == 0)then
!           user doesn't want to check for nulls
        nulchk=0
    else
!           user does want to check for null values
        if (tcode <= 41)then
!               check if null value is defined for integer column
            if (i4null == 123454321)then
                nulchk=0
            else
                if (tcode == 11)then
                    i1null=char(i4null)
                else if (tcode == 21)then
                    i2null=i4null
                end if
            end if
        end if
    end if

!       check for important special case: no datatype conversion required
    if (tcode == 41 .and. nulchk == 0 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
          trans=.false.
    else
          trans=.true.
    end if

    sval=' '
    i1=1
    ntodo=nelem
    rstart=0
    anynul=.false.
!       the data are being scaled from FITS to internal format
    tofits=.false.

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,(repeat-estart-1)/eincr+1,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp+(rstart * lenrow) + (estart * incre / eincr)
    call ftmbyt(iunit,bstart,.false.,status)

!       read the data from FITS file, doing datatype conversion and scaling
    if (tcode == 41)then
!               column data type is J (I*4)
!               read the data and do any machine dependent data conversion
!               note that we can use the input array directly
            call ftgi4b(iunit,itodo,incre,array(i1),status)
!               check for null values, and do scaling and datatype conversion
            if (trans)then
              call fti4i4(array(i1),itodo,scale,zero,tofits,nulchk, &
              i4null,nulval,flgval(i1),anynul,array(i1),status)
            end if
    else if (tcode == 21)then
!               column data type is I (I*2)
!               read the data and do any machine dependent data conversion
            call ftgi2b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call fti2i4(buffer,itodo,scale,zero,tofits, &
            nulchk,i2null,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 42)then
!               column data type is E (R*4)
!               read the data and do any machine dependent data conversion
            call ftgr4b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr4i4(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 82)then
!               column data type is D (R*8)
!               read the data and do any machine dependent data conversion
            call ftgr8b(iunit,itodo,incre,buffer,status)
!               check for null values, and do scaling and datatype conversion
            call ftr8i4(buffer,itodo,scale,zero,tofits, &
            nulchk,nulval,flgval(i1),anynul,array(i1),status)
    else if (tcode == 11)then
!               column data type is B (byte)
!               read the data and do any machine dependent data conversion
            call ftgi1b(iunit,itodo,incre,chbuff,status)
!               check for null values, and do scaling and datatype conversion
            call fti1i4(chbuff,itodo,scale,zero,tofits, &
            nulchk,i1null,nulval,flgval(i1),anynul,array(i1),status)
    else
!               this is an ASCII table column; get the character string
            call ftgcbf(iunit,twidth,sval,status)
            if (status > 0)return

!               check for null value
            if (sval(1:16) == snull)then
                anynul=.true.
                if (nultyp == 1)then
                    array(i1)=nulval
                else if (nultyp == 2)then
                    flgval(i1)=.true.
                end if
            else
!                   read the value, then do scaling and datatype conversion
                if (sform(5:5) == 'I')then
                    read(sval,sform,err=900)ival
                    dval=ival*scale+zero
                else
                    read(sval,sform,err=900)dval
                    dval=dval*scale+zero
                end if

!                   trap any values that overflow the I*4 range
                if (dval < i4max .and. dval > i4min)then
                    array(i1)=dval
                else if (dval >= i4max)then
                    status=-11
                    array(i1)=maxi4
                else
                    status=-11
                    array(i1)=mini4
                end if
            end if
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error reading elements',i9,' thru',i9, &
           ' of data array (FTGCLJ).')
        call ftpmsg(messge)
        return
    end if

    if (ntodo > 0)then
!           increment the pointers
        i1=i1+itodo
        estart=estart+itodo*eincr
        if (estart >= repeat)then
            rskip=estart/repeat
            rstart=rstart+rskip
            estart=estart-rskip*repeat
        end if
        go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while reading FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error reading formatted data value from ASCII table
    write(messge,1002)colnum,rstart+frow
1002    format('Error reading column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to read value with format '//sform)
    status=315
end
subroutine ftgcls(iunit,colnum,frow,felem,nelem,nultyp,nulval, &
      sray,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGCLS reads character strings from the specified column of the table.
!
!       The binary or ASCII table column being read must have datatype 'A'
!       This general purpose routine will handle null values in one
!       of two ways: if nultyp=1, then undefined array elements will be
!       set equal to the input value of NULVAL.  Else if nultyp=2, then
!       undefined array elements will have the corresponding FLGVAL element
!       set equal to .TRUE.  If NULTYP=1 and NULVAL=0, then no checks for
!       undefined values will be made, for maximum efficiency.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       frow    i  first row to read
!       felem   i  first element within row to read
!       nelem   i  number of elements to read
!       nultyp  i  input code indicating how to handle undefined values
!       nulval  c  value that undefined pixels will be set to (if nultyp=1)
!       sray    c  array of data values to be read
!       flgval  l  set .true. if corresponding element undefined (if nultyp=2)
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,nultyp,status
    logical flgval(*),anynul
    character ( len = * ) sray(*),nulval

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer bstart,nulchk,twidth,tread,tcode,offset,repeat
    integer ibuff,i1,ntodo,rstart,estart,lennul,strlen,nulfil
    character snull*16, messge*80

    if (status > 0)return
    ibuff=bufnum(iunit)

!       Do sanity check of input parameters
    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (hdutyp(ibuff) /= 1 .and. felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    else if (nelem == 0)then
      return
    end if

    anynul=.false.
    i1=1

!       column must be character string data type

    tcode=tdtype(colnum+tstart(ibuff))
    if (tcode == 16)then
!               for ASCII columns, TNULL actually stores the field width
            twidth=tnull(colnum+tstart(ibuff))
            ntodo=nelem
            rstart=frow-1
            repeat=trept(colnum+tstart(ibuff))
            if (felem > repeat)then
!                   illegal element number
                write(messge,1005)felem
1005                format( &
         'Starting element number is greater than repeat: ',i10)
                call ftpmsg(messge)
                status = 308
                return
            end if
            estart=felem-1
            bstart=dtstrt(ibuff)+rstart*rowlen(ibuff) &
                   +tbcol(colnum+tstart(ibuff))+estart*twidth
    else if (tcode == -16)then
!               this is a variable length descriptor field
            ntodo=1
!               read the string length and the starting offset:
            call ftgdes(iunit,colnum,frow,twidth,offset,status)
!               calc the i/o pointer position for the start of the string
            bstart=dtstrt(ibuff)+offset+theap(ibuff)
    else
!               error: not a character string column
            status=309
            call ftpmsg('Cannot to read character string'// &
            ' from a non-character column of a table (FTGCLS).')
            return
    end if

!       define the max. number of charcters to be read: either
!       the length of the variable length field, or the length
!       of the character string variable, which ever is smaller
    strlen=len(sray(1))
    tread=min(twidth,strlen)

!       move the i/o pointer to the start of the sequence of pixels
    call ftmbyt(iunit,bstart,.false.,status)

    lennul=0
!       determine if we have to check for null values
    if (nultyp == 1 .and. nulval == ' ')then
!               user doesn't want to check for nulls
            nulchk=0
    else
            nulchk=nultyp
            snull=cnull(colnum+tstart(ibuff))
!               lennul = length of the string to check for null values
            lennul=min(len(sray(1)),8)
    end if

!       process one string at a time
20      continue
!       get the string of characters
    sray(i1)=' '
    call ftgcbf(iunit,tread,sray(i1),status)
    if (status > 0)return

!       check for null value, if required
    if (nulchk /= 0)then
            if (ichar(sray(i1)(1:1)) == 0 .or. &
                sray(i1)(1:lennul) == snull(1:lennul))then
                    if (nulchk == 1)then
                            sray(i1)=nulval
                            anynul=.true.
                    else
                            flgval(i1)=.true.
                            anynul=.true.
                    end if
            end if
    end if

!       check for null terminated string; pad out with blanks if found
    nulfil=index(sray(i1),char(0))
    if (nulfil > 1)then
            sray(i1)(nulfil:len(sray(1)))=' '
    end if

    if (status > 0)then
        write(messge,1006)i1
1006        format('Error reading string for element',i9, &
           ' of data array (FTGCLS).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-1
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+1
            estart=estart+1
            if (estart == repeat)then
                    rstart=rstart+1
                    estart=0
            end if
!               move to the start of the next string; need to do
!               this every time in case we didn't read all the characters
!               from the previous string.
            bstart=dtstrt(ibuff)+rstart*rowlen(ibuff) &
                   +tbcol(colnum+tstart(ibuff))+estart*twidth
!               move the i/o pointer
            call ftmbyt(iunit,bstart,.false.,status)
            go to 20
    end if
end
subroutine ftgcnn(iunit,casesn,templt,colnam,colnum,status)
!
!*******************************************************************************
!
!! FTGCNN determines column name and number from a column name template string.
!
!       The template may contain the * and ?
!       wildcards.  Status = 237 is returned if match is not unique.
!       One may call this routine again with input status=237  to
!       get the next match.

!       iunit   i  Fortran i/o unit number
!       casesn  l  true if an exact case match of the names is required
!       templt  c  templt for column name
!       colnam  c  name of (first) column that matchs the template
!       colnum  i  number of the column (first column = 1)
!                       (a value of 0 is returned if the column is not found)
!       status  i  returned error status

!       written by Wm Pence, HEASARC/GSFC, December 1994

    integer iunit,colnum,status
    character ( len = * ) templt,colnam
    logical casesn

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    integer colpnt,untpnt
    common/ftname/colpnt,untpnt
!

    integer ibuff,i,nfound,tstat,ival
    logical match,exact,founde,foundw,unique
    character*80 errmsg
    character*68 tname(999)
    save tname

    ibuff=bufnum(iunit)

!       load the common block with names, if not already defined
    if (colpnt == -999 .or. iunit /= untpnt)then
        do 10 i=1,tfield(ibuff)
            tname(i)=' '
10          continue
        call ftgkns(iunit,'TTYPE',1,nf,tname,nfound,status)
        if (status > 0)return
        untpnt=iunit
        colpnt=1
    end if

    if (status <= 0)then
        tstat=0
        colpnt=1
    else if (status == 237)then
!           search for next non-unique match, starting from the previous match
        tstat=237
        status=0
    else
        return
    end if

    colnam=' '
    colnum=0

!       set the 'found exact' and 'found wildcard' flags to false
    founde=.false.
    foundw=.false.

    do 100 i=colpnt,tfield(ibuff)
!               test for match between template and column name
            call ftcmps(templt,tname(i),casesn,match,exact)

            if (match)then
                if (founde .and. exact)then
!                       warning: this is the second exact match we've found
!                       reset pointer to first match so next search starts there
                    colpnt=colnum+1
                    status=237
                    return
                else if (founde)then
!                       already found exact match so ignore this non-exact match
                else if (exact)then
!                       this is the first exact match we have found, so save it.
                    colnam=tname(i)
                    colnum=i
                    founde=.true.
                else if (foundw)then
!                       we have already found a wild card match, so not unique
!                       continue searching for other matches
                    unique=.false.
                else
!                       this is the first wild card match we've found. save it
                    colnam=tname(i)
                    colnum=i
                    foundw=.true.
                    unique=.true.
                end if
            end if
100     continue

!       OK, we've checked all the names now see if we got any matches
    if (founde)then
!           we did find 1 exact match
        if (tstat == 237)status=237
    else if (foundw)then
!           we found one or more wildcard matches
!           report error if not unique
        if (.not. unique .or. tstat == 237)status=237
    else
!           didn't find a match; check if template is a simple positive integer
        call ftc2ii(templt,ival,tstat)
        if (tstat == 0 .and. ival <= tfield(ibuff) &
            .and. ival > 0)then
            colnum=ival
            colnam=tname(ival)
        else
            status=219
            if (tstat /= 237)then
              errmsg='FTGCNN: Could not find column: '//templt
              call ftpmsg(errmsg)
            end if
        end if
    end if

!       reset pointer so next search starts here if input status=237
    colpnt=colnum+1
end
subroutine ftgcno(iunit,casesn,templt,colnum,status)
!
!*******************************************************************************
!
!! FTGCNO determines the column number corresponding to an input column name.
!
!       This supports the * and ? wild cards in the input template.

!       iunit   i  Fortran i/o unit number
!       casesn  l  true if an exact case match of the names is required
!       templt  c  name of column as specified in a TTYPE keyword
!       colnum  i  number of the column (first column = 1)
!                       (a value of 0 is returned if the column is not found)
!       status  i  returned error status

!       modified by Wm Pence, HEASARC/GSFC, December 1994

    integer iunit,colnum,status
    character ( len = * ) templt
    logical casesn
    character*8 dummy

    call ftgcnn(iunit,casesn,templt,dummy,colnum,status)
end
subroutine ftgcpr(iunit,colnum,frow,felem,nelem,rwmode, &
   ibuff,scale,zero,tform,twidth,tcode,maxelm,startp, &
   elnum,incre,repeat,lenrow,hdtype,inull,snull,status)
!
!*******************************************************************************
!
!! FTGCPR gets column parameters.
!
!  It also tests starting row and element numbers for validity.

!       iunit   I - fortran unit number
!       colnum  I - column number (1 = 1st column of table)
!       frow    I - first row (1 = 1st row of table)
!       felem   I - first element within vector (1 = 1st)
!       nelem   I - number of elements to read or write
!       rwmode  I - = 1 if writing data, = 0 if reading data
!       ibuff   O - buffer associated with this file
!       scale   O - FITS scaling factor (TSCALn keyword value)
!       zero    O - FITS scaling zero pt (TZEROn keyword value)
!       tform   O - ASCII column format: value of TFORMn keyword
!       twidth  O - width of ASCII column (characters)
!       tcode   O - column datatype code: I*4=41, R*4=42, etc
!       maxelm  O - max number of elements that fit in buffer
!       startp  O - offset in file to starting row & column
!       elnum   O - starting element number ( 0 = 1st element)
!       incre   O - byte offset between elements within a row
!       repeat  O - number of elements in a row (vector column)
!       lenrow  O - length of a row, in bytes
!       hdtype  O - HDU type: 0, 1, 2 = primary, table, bintable
!       inull   O - null value for integer columns
!       snull   O - null value for ASCII table columns
!       status IO - error status

!       written by Wm Pence, HEASARC/GSFC, November 1996

    integer iunit,colnum,frow,felem,nelem
    integer rwmode,ibuff,twidth,tcode,maxelm,startp
    integer elnum,incre,repeat,lenrow,hdtype,inull
    integer status
    character ( len = * ) snull, tform
    double precision scale,zero

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
    integer compid
    common/ftcpid/compid
!

    integer datast, xtbcol,acode
    character*80 messge
    integer bufdim
    parameter (bufdim = 32000)

    ibuff=bufnum(iunit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(iunit,status)

!       Do sanity check of input parameters
    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (hdutyp(ibuff) /= 1 .and. felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    else if (nelem == 0)then
!         not reading or writing any pixels, so just return
      return
    end if

!       copy relevant parameters from the common block
    hdtype = hdutyp(ibuff)
    lenrow = rowlen(ibuff)
    datast = dtstrt(ibuff)
    tcode  = tdtype(colnum+tstart(ibuff))
    tform  ='(           )'
    tform(5:12)=cform(colnum+tstart(ibuff))

    acode = abs(tcode)
    if ((hdtype == 1 .and. tform(5:5) == 'A') .or. &
      (hdtype == 2 .and. acode == 16) .or. &
       acode == 14)then
!          error: illegal table format code
       status=311
       write(messge,1005)colnum,cform(colnum+tstart(ibuff))
1005       format('Cannot read or write numerical values in column', &
       i4,' with TFORM = ',a8)
       call ftpmsg(messge)
       return
    end if

    if (hdtype == 1 .and. rwmode == 1)then
       if (tform(5:5) == 'E')then
           tform(2:4)='1P,'
       else if (tform(5:5) == 'D')then
           tform(2:5)='1P,E'
       end if
    else if (hdtype == 1)then
       tform(2:4)='BN,'
    end if

    snull =  cnull(colnum+tstart(ibuff))
    scale=  tscale(colnum+tstart(ibuff))
    zero=    tzero(colnum+tstart(ibuff))
    inull=   tnull(colnum+tstart(ibuff))
    xtbcol=  tbcol(colnum+tstart(ibuff))
    repeat=  trept(colnum+tstart(ibuff))

    if (tcode /= 16)then
      twidth=max(acode/10,1)
    else
      twidth = tnull(colnum+tstart(ibuff))
    end if

!       Special case: interprete 'X' column as 'B'
    if (acode == 1)then
       tcode  = tcode * 11
       repeat = (repeat + 7) / 8
    end if

!       Special case: support the 'rAw' format in BINTABLEs
    if (hdtype == 2 .and. tcode == 16)then
      repeat =  repeat /  twidth
    end if

    if (hdtype == 1)then
!         ASCII tables don't have vector elements
      elnum = 0
    else
      elnum = felem - 1
    end if

!       interprete complex and double complex as pairs of floats or doubles
    if (abs(tcode) > 82)then
      if (tcode > 0)then
         tcode = (tcode + 1) / 2
      else
         tcode = (tcode - 1) / 2
      end if

      repeat  = repeat * 2
      twidth  = twidth / 2
    end if

    incre= twidth

!       calculate no. of pixels that fit in buffer
    if (hdtype == 1)then
!           in ASCII tables, can only process 1 value at a time
        maxelm = 1
    else
        maxelm = bufdim / twidth
    end if

!       special case for the SUN F90 compiler where integer*2
!       variables are stored in 4-byte integers
    if (compid == -1 .and. abs(tcode) == 21)then
        maxelm = bufdim / 4
    end if

!       calc starting byte position to 1st element of col
!       (this does not apply to variable length columns)
    startp = datast + ((frow - 1) * lenrow) + xtbcol

    if (hdtype == 0 .and. rwmode == 1)then

!       When writing primary arrays, set the repeat count greater than the
!       total number of pixels to be written.  This prevents an out-of-range
!       error message in cases where the final image array size is not
!       yet known or defined.

      repeat = elnum + nelem

    else if (tcode > 0)then
!         Fixed length table column

      if (elnum >= repeat)then
!                       illegal element number
         write(messge,1006)felem
1006         format( &
         'Starting element number is greater than repeat: ',i10)
         call ftpmsg(messge)
         status = 308
      else if (repeat == 1 .and. nelem > 1)then

!            When accessing a scalar column, fool the calling routine into
!            thinking that this is a vector column with very big elements.
!            This allows multiple values (up to the maxelem number of elements
!            that will fit in the buffer) to be read or written with a single
!            routine call, which increases the efficiency.

         incre = lenrow
         repeat = nelem
      end if
    else
!         Variable length Binary Table column

      tcode = tcode * (-1)

      if (rwmode ==  1)then
!           return next empty heap address for writing

!           total no. of elements in the field
        repeat = nelem + elnum

!           calculate starting position (for writing new data) in the heap
        startp = datast + heapsz(ibuff)+theap(ibuff)

!           write the descriptor into the fixed length part of table
        call ftpdes(iunit, colnum, frow, repeat, heapsz(ibuff), &
                    status)

!           increment the address to the next empty heap position
        heapsz(ibuff) = heapsz(ibuff) + (repeat * incre)
      else
!           get the read start position in the heap

        call ftgdes(iunit, colnum, frow, repeat, startp, status)

        if (tdtype(colnum+tstart(ibuff)) == -1)then
!               Special case: interprete 'X' column as 'B'
            repeat = (repeat + 7) / 8
        end if

        if (elnum >= repeat)then
!                       illegal element number
           write(messge,1006)felem
           call ftpmsg(messge)
           status = 308
        end if

        startp=datast + startp + theap(ibuff)
      end if
    end if
end
subroutine ftgcrd(iunit,keynam,card,status)
!
!*******************************************************************************
!
!! FTGCRD reads the 80 character card image of a header keyword record.
!
!    If the input name contains wild cards ('?' matches any single char
!    and '*' matches any sequence of chars, # matches any string of decimal
!    digits) then the search ends once the end of header is reached and does
!    not automatically resume from the top of the header.

!       iunit   i  Fortran I/O unit number
!       keynam  c  name of keyword to be read
!       OUTPUT PARAMETERS:
!       card    c  80 character card image that was read
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991
!       modified January 1997 to support wildcards

    integer iunit,status
    character ( len = * ) keynam,card

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

    integer i,j,ibuff,maxkey,start
    character kname*9
    character*80 keybuf
    logical wild,casesn,match,exact

    card=' '
    if (status > 0)go to 100
    casesn=.true.

!       get the number of the data buffer used for this unit
    ibuff=bufnum(iunit)

!       make sure keyword name is in uppercase
    kname=keynam
    call ftupch(kname)

!       test if input name contains wild card characters
    wild=.false.
    do i=1,9
      if (kname(i:i) == '?' .or. kname(i:i) == '*' &
     .or. kname(i:i) == '#')wild=.true.
    end do

!       Start by searching for keyword from current pointer position to the end.
!       Calculate the maximum number of keywords to be searched:
    start=nxthdr(ibuff)
    maxkey=(hdend(ibuff)-start)/80

    do 20 j=1,2
!           position I/O pointer to the next header keyword
        if (maxkey > 0)then
            call ftmbyt(iunit,start,.false.,status)
        end if

        do i=1,maxkey
            call ftgcbf(iunit,80,keybuf,status)
            if (status > 0)go to 100
            if (wild)then
              call ftcmps(kname(1:8),keybuf(1:8),casesn,match,exact)
              if (match)then
!                     setheader pointer to the following keyword
                  nxthdr(ibuff)=start+i*80
                  card=keybuf
                  return
              end if
            else if (keybuf(1:8) == kname(1:8))then
!                       setheader pointer to the following keyword
                    nxthdr(ibuff)=start+i*80
                    card=keybuf
                    return
            end if
        end do

!       end search at end of header if input name contains wildcards
        if (wild .or. (j == 2))go to 30

!           didn't find keyword yet, so now search from top down to starting pt.
!           calculate max number of keywords to be searched and reset nxthdr
        maxkey=(start-hdstrt(ibuff,chdu(ibuff)))/80
        start=hdstrt(ibuff,chdu(ibuff))
20      continue

!       keyword was not found
30      status=202

!       don't write to error stack because this innoculous error happens a lot
!       call ftpmsg('Could not find the '//kname//' keyword to read.')

100     continue
end
subroutine ftgcvb(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVB reads an array of byte values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  b  value that undefined pixels will be set to
!       array   b  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul

    character array(*),nulval

    call ftgclb(iunit,colnum,frow,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcvc(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVC reads an array of complex values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  cmp  value that undefined pixels will be set to
!       array   cmp  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    real array(*),nulval(2)
    integer felemx, nelemx

!       a complex value is interpreted as a pair of float values, thus
!       need to multiply the first element and number of elements by 2

    felemx = (felem - 1) * 2 + 1
    nelemx = nelem * 2

    call ftgcle(iunit,colnum,frow,felemx,nelemx,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcvd(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVD reads an array of r*8 values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  d  value that undefined pixels will be set to
!       array   d  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    double precision array(*),nulval

    call ftgcld(iunit,colnum,frow,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcve(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVE reads an array of R*4 values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  r  value that undefined pixels will be set to
!       array   r  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    real array(*),nulval

    call ftgcle(iunit,colnum,frow,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcvi(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVI reads an array of I*2 values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  i*2  value that undefined pixels will be set to
!       array   i*2 returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    integer*2 array(*),nulval

    call ftgcli(iunit,colnum,frow,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcvj(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVJ reads an array of I*4 values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  i  value that undefined pixels will be set to
!       array   i  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    integer array(*),nulval

    call ftgclj(iunit,colnum,frow,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcvm(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVM reads double precision complex values from a column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=0, in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element within the row to read
!       nelem   i  number of elements to read
!       nulval  dcmp  value that undefined pixels will be set to
!       array   dcmp  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    double precision array(*),nulval(2)
    integer felemx, nelemx

!       a complex value is interpreted as a pair of float values, thus
!       need to multiply the first element and number of elements by 2

    felemx = (felem - 1) * 2 + 1
    nelemx = nelem * 2

    call ftgcld(iunit,colnum,frow,felemx,nelemx,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcvs(iunit,colnum,frow,felem,nelem,nulval,array, &
            anynul,status)
!
!*******************************************************************************
!
!! FTGCVS reads an array of string values from a specified column of the table.
!
!       Any undefined pixels will be set equal to the value of NULVAL,
!       unless NULVAL=' ', in which case no checks for undefined pixels
!       will be made.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       felem   i  first element in the row to read
!       nelem   i  number of elements to read
!       nulval  c  value that undefined pixels will be set to
!       array   c  returned array of data values that was read from FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,colnum,frow,felem,nelem,status
    logical flgval,anynul
    character ( len = * ) array(*),nulval

    call ftgcls(iunit,colnum,frow,felem,nelem,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgcx(iunit,colnum,frow,fbit,nbit,lray,status)
!
!*******************************************************************************
!
!! FTGCX reads logical values from a bit or byte column of the binary table.
!
!       A logical .true. value is returned
!       if the corresponding bit is 1, and a logical .false. value is
!       returned if the bit is 0.
!       The binary table column being read from must have datatype 'B'
!       or 'X'. This routine ignores any undefined values in the 'B' array.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       fbit    i  first bit within the row to read
!       nbit    i  number of bits to read
!       lray    l  returned array of logical data values that is read
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Mar 1992

    integer iunit,colnum,frow,fbit,nbit,status
    logical lray(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer bstart,offset,tcode,fbyte,bitloc,ndone
    integer ibuff,i,ntodo,repeat,rstart,estart,buffer
    logical descrp,log8(8)
    character cbuff

    if (status > 0)return

    ibuff=bufnum(iunit)
    tcode=tdtype(colnum+tstart(ibuff))

!       check input parameters
    if (nbit <= 0)then
            return
    else if (frow < 1)then
!               error: illegal first row number
            status=307
            return
    else if (fbit < 1)then
!               illegal element number
            status=308
            return
    end if

    fbyte=(fbit+7)/8
    bitloc=fbit-(fbit-1)/8*8
    ndone=0
    ntodo=nbit
    rstart=frow-1
    estart=fbyte-1

    if (tcode == 11)then
            repeat=trept(colnum+tstart(ibuff))
            if (fbyte > repeat)then
!                       illegal element number
                    status=308
                    return
            end if
            descrp=.false.
!               move the i/o pointer to the start of the sequence of pixels
            bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
            tbcol(colnum+tstart(ibuff))+estart
    else if (tcode == -11)then
!               this is a variable length descriptor column
            descrp=.true.
!               read the number of elements and the starting offset:
            call ftgdes(iunit,colnum,frow,repeat, &
                                offset,status)
            repeat=(repeat+7)/8
            if (repeat == 0)then
!                       error: null length vector
                    status=318
                    return
            else if ((fbit+nbit+6)/8 > repeat)then
!                       error: trying to read beyond end of record
                    status=319
                    return
            end if
            bstart=dtstrt(ibuff)+offset+ &
                            theap(ibuff)+estart
    else
!               column must be byte or bit data type
            status=312
            return
    end if

!       move the i/o pointer to the start of the pixel sequence
    call ftmbyt(iunit,bstart,.false.,status)

!       get the next byte
20      call ftgcbf(iunit,1,cbuff,status)
    buffer=ichar(cbuff)
    if (buffer < 0)buffer=buffer+256

!       decode the bits within the byte into an array of logical values
    call ftgbit(buffer,log8)

    do i=bitloc,8
            ndone=ndone+1
            lray(ndone)=log8(i)
            if (ndone == ntodo)go to 100
    end do

!       not done, so get the next byte
    if (.not. descrp)then
            estart=estart+1
            if (estart == repeat)then
!                       move the i/o pointer to the next row of pixels
                    estart=0
                    rstart=rstart+1
                    bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
                           tbcol(colnum+tstart(ibuff))+estart
                    call ftmbyt(iunit,bstart,.false.,status)
            end if
    end if
    bitloc=1
    go to 20

100     continue
end
subroutine ftgcxd(iunit,colnum,frow,nrow,fbit,nbit, &
               dvalue,status)
!
!*******************************************************************************
!
!! FTGCXD reads bits from an 'X' or 'B' column as an unsigned n-bit integer.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       nrow    i  number of rows to read
!       fbit    i  first bit within the row to read
!       nbit    i  number of bits to read
!       dvalue  d  returned value(s)
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Nov 1994

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer iunit,colnum,fbit,nbit,frow,nrow,status
    integer i,k,istart,itodo,ntodo,row,ibuff
    double precision dvalue(*),power,dval
    logical lray(64)

    if (status > 0)return

    ibuff=bufnum(iunit)
    if ((fbit+nbit+6)/8 > trept(colnum+tstart(ibuff)))then
        call ftpmsg('Asked to read more bits than exist in'// &
        ' the column (ftgcxd)')
        status=308
        return
    end if

    row=frow-1
    do 30 k=1,nrow
        row=row+1
        dval=0.
        power=1.0D+00
        istart=fbit+nbit
        ntodo=nbit

10          itodo=min(ntodo,64)
        istart=istart-itodo

!           read up to 64 bits at a time
!           get the individual bits
        call ftgcx(iunit,colnum,row,istart,itodo,lray,status)
        if (status > 0)return

!           reconstruct the positive integer value
        do 20 i=itodo,1,-1
            if (lray(i))dval=dval+power
            power=power*2.0D+00
20          continue

        ntodo=ntodo-itodo
        if (itodo > 0)go to 10
        dvalue(k)=dval
30      continue
end
subroutine ftgcxi(iunit,colnum,frow,nrow,fbit,nbit, &
                      ivalue,status)
!
!*******************************************************************************
!
!! FTGCXI reads bits from an 'X' or 'B' column as an unsigned n-bit integer.
!
!  This is the case unless nbits=16 in which case the 16 bits
!       are interpreted as a 16-bit signed 2s complement word

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       nrow    i  number of rows to read
!       fbit    i  first bit within the row to read
!       nbit    i  number of bits to read
!       ivalue  i*2  returned integer value(s)
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Nov 1994

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer iunit,colnum,fbit,nbit,frow,nrow,status,i,j,k,row,ibuff
    integer*2 ivalue(*),ival,power2(16)
    logical lray(16)
    save power2
    data power2/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192, &
    16384,0/

    if (status > 0)return

    ibuff=bufnum(iunit)

    if (nbit > 16)then
        call ftpmsg('Cannot read more than 16 bits (ftgcxi)')
        status=308
        return
    else if ((fbit+nbit+6)/8 > trept(colnum+tstart(ibuff)))then
        call ftpmsg('Asked to read more bits than exist in'// &
        ' the column (ftgcxi)')
        status=308
        return
    end if


    row=frow-1
    do 30 k=1,nrow
        row=row+1
!           get the individual bits
        call ftgcx(iunit,colnum,row,fbit,nbit,lray,status)
        if (status > 0)return
        ival=0
        j=0
        if (nbit == 16 .and. lray(1))then
!               interprete this as a 16 bit negative integer
            do 10 i=16,2,-1
                j=j+1
                if (.not. lray(i))ival=ival+power2(j)
10              continue
!               make 2's complement
            ivalue(k)=-ival-1
        else
!               reconstruct the positive integer value
            do 20 i=nbit,1,-1
                j=j+1
                if (lray(i))ival=ival+power2(j)
20              continue
            ivalue(k)=ival
        end if
30      continue
end
subroutine ftgcxj(iunit,colnum,frow,nrow,fbit,nbit, &
               jvalue,status)
!
!*******************************************************************************
!
!! FTGCXJ reads bits from an 'X' or 'B' column as an unsigned n-bit integer.
!
!  This is the case unless nbits=32 in which case the 32 bits
!       are interpreted as a 32-bit signed 2s complement word

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       frow    i  first row to read
!       nrow    i  number of rows to read
!       fbit    i  first bit within the row to read
!       nbit    i  number of bits to read
!       jvalue  i  returned integer value(s)
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Nov 1994

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer iunit,colnum,fbit,nbit,frow,nrow,status,i,j,k,row,jval
    integer jvalue(*),power2(32),ibuff
    logical lray(32)
    save power2
    data power2/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192, &
    16384,32768,65536,131072,262144,524288,1048576,2097152,4194304, &
    8388608,16777216,33554432,67108864,134217728,268435456,536870912 &
    ,1073741824,0/

    if (status > 0)return

    ibuff=bufnum(iunit)

    if (nbit > 32)then
        call ftpmsg('Cannot read more than 32 bits (ftgcxj)')
        status=308
        return
    else if ((fbit+nbit+6)/8 > trept(colnum+tstart(ibuff)))then
        call ftpmsg('Asked to read more bits than exist in'// &
        ' the column (ftgcxj)')
        status=308
        return
    end if

    row=frow-1
    do 30 k=1,nrow
        row=row+1
!           get the individual bits
        call ftgcx(iunit,colnum,row,fbit,nbit,lray,status)
        if (status > 0)return

        jval=0
        j=0
        if (nbit == 32 .and. lray(1))then
!               interprete this as a 32 bit negative integer
            do 10 i=32,2,-1
                j=j+1
                if (.not. lray(i))jval=jval+power2(j)
10              continue
!               make 2's complement
            jvalue(k)=-jval-1
        else
!               reconstruct the positive integer value
            do 20 i=nbit,1,-1
                j=j+1
                if (lray(i))jval=jval+power2(j)
20              continue
            jvalue(k)=jval
        end if
30      continue
end
subroutine ftgdes(iunit,colnum,rownum,nelem,offset,status)
!
!*******************************************************************************
!
!! FTGDES reads the descriptor values from a binary table.
!
!  This is only
!       used for column which have TFORMn = 'P', i.e., for variable
!       length arrays.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to read
!       rownum  i  number of the row to read
!       nelem   i  output number of elements
!       offset  i  output byte offset of the first element
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Nov 1991

    integer iunit,colnum,rownum,nelem,offset,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bstart,iray(2)

    if (status > 0)return
    if (rownum < 1)then
!               error: illegal row number
            status=307
            return
    end if

    ibuff=bufnum(iunit)

!       check that this is really a 'P' type column
    if (tdtype(colnum+tstart(ibuff)) >= 0)then
            status=317
            return
    end if

!       move to the specified column and row:
    bstart=dtstrt(ibuff)+(rownum-1)*rowlen(ibuff) &
           +tbcol(colnum+tstart(ibuff))
    call ftmbyt(iunit,bstart,.true.,status)

!       now read the number of elements and the offset to the table:
    call ftgi4b(iunit,2,0,iray,status)
    nelem=iray(1)
    offset=iray(2)
end
subroutine ftgerr(errnum,text)
!
!*******************************************************************************
!
!! FTGERR returns a descriptive error message corresponding to the error number
!
!       errnum i  input symbolic error code presumably returned by another
!                 FITSIO subroutine
!       text   C*30  Descriptive error message

    integer errnum
    character ( len = * ) text

!       nerror specifies the maxinum number of different error messages
    integer nerror
    parameter (nerror=100)
    character ( len = 30 ) errors(nerror)
    character ( len = 30 ) er1(10),er2(10),er3(10)
    character ( len = 30 ) er4(10),er5(10),er6(10)
    character ( len = 30 ) er7(10),er8(10),er9(10),er10(10)
    integer i,errcod(nerror)
    save errors

!       we equivalence the big array to several smaller ones, so that
!       the DATA statements will not have too many continuation lines.
    equivalence (errors(1), er1(1))
    equivalence (errors(11),er2(1))
    equivalence (errors(21),er3(1))
    equivalence (errors(31),er4(1))
    equivalence (errors(41),er5(1))
    equivalence (errors(51),er6(1))
    equivalence (errors(61),er7(1))
    equivalence (errors(71),er8(1))
    equivalence (errors(81),er9(1))
    equivalence (errors(91),er10(1))

    data errcod/0,101,102,103,104,105,106,107,108,109,110,111, &
    201,202,203,204,205,206,207,208,209,211,212,213,214,215,216, &
    217,218,221,222,223,224,225,226,227,228,229,230,231,232, &
    241,251,252,261,262, &
    302,303,304,305,306,307,308,309,310,311,312,313,314,315,316, &
    317,318,319,    401,402,403,404,405,406,407,408,409,411,112, &
    210,233,220,219,301,320,321,322,263,323,113,114,234,253,254, &
    255,412,235,236,501,502,503,504,505,237/

    data er1/ &
   'OK, no error', &
   'Bad logical unit number', &
   'Too many FITS files opened', &
   'File not found; not opened', &
   'Error opening existing file', &
   'Error creating new FITS file', &
   'Error writing to FITS file', &
   'EOF while reading FITS file', &
   'Error reading FITS file', &
   'Bad blocking factor (1-28800)'/

    data er2/ &
   'Error closing FITS file', &
   'Too many columns in table', &
   'Header is not empty', &
   'Specified keyword not found', &
   'Bad keyword record number', &
   'Keyword value is undefined', &
   'Missing quote in string value', &
   'Could not construct NAMEnnn', &
   'Bad character in header record', &
   'Keywords out of order?'/

    data er3/ &
   'Bad nnn value in NAMEnnn', &
   'Illegal BITPIX keyword value', &
   'Illegal NAXIS keyword value', &
   'Illegal NAXISnnn keyword value', &
   'Illegal PCOUNT keyword value', &
   'Illegal GCOUNT keyword value', &
   'Illegal TFIELDS keyword value', &
   'Illegal NAXIS1 keyword value', &
   'Illegal NAXIS2 keyword value', &
   'SIMPLE keyword not found'/

    data er4/ &
   'BITPIX keyword not found', &
   'NAXIS  keyword not found', &
   'NAXISnnn keyword(s) not found', &
   'XTENSION keyword not found', &
   'CHDU is not an ASCII table', &
   'CHDU is not a binary table', &
   'PCOUNT keyword not found', &
   'GCOUNT keyword not found', &
   'TFIELDS keyword not found', &
   'TBCOLnnn keywords not found'/

    data er5/ &
   'TFORMnnn keywords not found', &
   'Row width not = field widths', &
   'Unknown extension type', &
   'Unknown FITS record type', &
   'Cannot parse TFORM keyword', &
   'Unknown TFORM datatype code', &
   'Column number out of range', &
   'Data structure not defined', &
   'Negative file record number', &
   'HDU start location is unknown'/

    data er6/ &
   'Requested no. of bytes < 0', &
   'Illegal first row number', &
   'Illegal first element number', &
   'Bad TFORM for Character I/O', &
   'Bad TFORM for Logical I/O', &
   'Invalid ASCII table TFORM code', &
   'Invalid BINTABLE TFORM code', &
   'Error making formated string', &
   'Null value is undefined', &
   'Internal read error of string'/

    data er7/ &
   'Illegal logical column value', &
   'Bad TFORM for descriptor I/O', &
   'Variable array has 0 length', &
   'End-of-rec in var. len. array', &
   'Int to Char conversion error', &
   'Real to Char conversion error', &
   'Illegal Char to Int conversion', &
   'Illegal Logical keyword value', &
   'Illegal Char to R*4 conversion', &
   'Illegal Char to R*8 conversion'/

    data er8/ &
   'Char to Int conversion error', &
   'Char to Real conversion error', &
   'Char to R*8 conversion error', &
   'Illegal no. of decimal places', &
   'Cannot modify a READONLY file', &
   'END header keyword not found', &
   'CHDU is not an IMAGE extension', &
   'Illegal SIMPLE keyword value', &
   'Column name (TTYPE) not found', &
   'Out of bounds HDU number'/

    data er9/ &
   'Bad no. of array dimensions', &
   'Max pixel less than min pixel', &
   'Illegal BSCALE or TSCALn = 0', &
   'Could not parse TDIMn keyword', &
   'Axis length less than 1', &
   'Incompatible FITSIO version', &
   'All LUNs have been allocated', &
   'TBCOLn value out of range', &
   'END keyword value not blank ', &
   'Header fill area not blank'/

    data er10/ &
   'Data fill area invalid', &
   'Data type conversion overflow', &
   'CHDU must be a table/bintable', &
   'Column is too wide for table', &
   'celestial angle too large', &
   'bad celestial coordinate', &
   'error in celestial coord calc', &
   'unsupported projection', &
   'missing celestial coord keywrd', &
   'column name not unique'/

!       find the matching error code number
    do i=1,nerror
            if (errnum == errcod(i))then
                    text=errors(i)
                    return
            end if
    end do

    text='Unknown FITSIO status code'
end
subroutine ftgext(iunit,extno,xtend,status)
!
!*******************************************************************************
!
!! FTGEXT moves the IO pointer to another extension.
!
!       'Get Extension'
!       move i/o pointer to another extension (or the primary HDU) and
!       initialize all the common block parameters which describe the
!       extension

!       iunit   i  fortran unit number
!       extno   i  number of the extension to point to.
!       xtend   i  type of extension:   0 = the primary HDU
!                                       1 = an ASCII table
!                                       2 = a binary table
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991

    integer iunit,extno,xtend,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,xchdu,xhdend,xmaxhd

    if (status > 0)return

    ibuff=bufnum(iunit)

!       move to the beginning of the desired extension
    call ftmbyt(iunit,hdstrt(ibuff,extno),.false.,status)
    if (status <= 0)then

!               temporarily save parameters
            xchdu=chdu(ibuff)
            xmaxhd=maxhdu(ibuff)
            xhdend=hdend(ibuff)

!               initialize various parameters about the CHDU
            chdu(ibuff)=extno
            maxhdu(ibuff)=max(extno,maxhdu(ibuff))
!               the location of the END record is currently unknown, so
!               temporarily just set it to a very large number
            hdend(ibuff)=2000000000

!               determine the structure of the CHDU
            call ftrhdu(iunit,xtend,status)
            if (status > 0)then
!                       couldn't read the extension so restore previous state
                    chdu(ibuff)= xchdu
                    maxhdu(ibuff)=xmaxhd
                    hdend(ibuff)= xhdend
            end if
    end if
end
subroutine ftggpb(iunit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTGGPB reads an array of group parameter values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter be read (starting with 1)
!       nparm   i  number of group parameters to be read
!       array   b  returned array of values that were read
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,fparm,nparm,status,row
    character nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
!       set nulval to blank to inhibit checking for undefined values
    nulval=' '
    row=max(1,group)
    call ftgclb(iunit,1,row,fparm,nparm,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftggpd(iunit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTGGPD reads an array of group parameter values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter be read (starting with 1)
!       nparm   i  number of group parameters to be read
!       array   d  returned array of values that were read
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,fparm,nparm,status,row
    double precision nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
!       set nulval to blank to inhibit checking for undefined values
    nulval=0
    row=max(1,group)
    call ftgcld(iunit,1,row,fparm,nparm,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftggpe(iunit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTGGPE reads an array of group parameter values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter be read (starting with 1)
!       nparm   i  number of group parameters to be read
!       array   r  returned array of values that were read
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,fparm,nparm,status,row
    real nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
!       set nulval to blank to inhibit checking for undefined values
    nulval=0
    row=max(1,group)
    call ftgcle(iunit,1,row,fparm,nparm,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftggpi(iunit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTGGPI reads an array of group parameter values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter be read (starting with 1)
!       nparm   i  number of group parameters to be read
!       array   i*2  returned array of values that were read
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,fparm,nparm,status,row
    integer*2 nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
!       set nulval to blank to inhibit checking for undefined values
    nulval=0
    row=max(1,group)
    call ftgcli(iunit,1,row,fparm,nparm,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftggpj(iunit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTGGPJ reads an array of group parameter values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter be read (starting with 1)
!       nparm   i  number of group parameters to be read
!       array   i  returned array of values that were read
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,fparm,nparm,status,row
    integer nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
!       set nulval to blank to inhibit checking for undefined values
    nulval=0
    row=max(1,group)
    call ftgclj(iunit,1,row,fparm,nparm,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftghad(iunit,curhdu,nxthdu)
!
!*******************************************************************************
!
!! FTGHAD returns the starting byte address of the CHDU and the next HDU.
!
!       curhdu  i  starting address of the CHDU
!       nxthdu  i  starting address of the next HDU

!       written by Wm Pence, HEASARC/GSFC, May, 1995

    integer iunit,curhdu,nxthdu

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,hdunum

    ibuff=bufnum(iunit)
    hdunum=chdu(ibuff)
    curhdu=hdstrt(ibuff,hdunum)
    nxthdu=hdstrt(ibuff,hdunum+1)
end
subroutine ftghbn(iunit,maxfld,nrows,nfield,ttype,tform, &
                      tunit,extnam,pcount,status)
!
!*******************************************************************************
!
!! FTGHBN reads required standard header keywords from a binary table extension.
!
!       iunit   i  Fortran i/o unit number
!       maxfld  i  maximum no. of fields to read; size of ttype array
!       OUTPUT PARAMETERS:
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       ttype   c  name of each field (array)
!       tform   c  format of each field (array)
!       tunit   c  units of each field (array)
!       extnam  c  name of table (optional)
!       pcount  i  size of special data area following the table (usually = 0)
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,maxfld,ncols,nrows,nfield,pcount,status,tstat
    integer maxf,i,nfind
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam
    character comm*72

!       check that this is a valid binary table and get parameters
    call ftgtbn(iunit,ncols,nrows,pcount,nfield,status)
    if (status > 0)return

    if (maxfld < 0)then
          maxf=nfield
    else if (maxfld == 0)then
          go to 20
    else
          maxf=min(maxfld,nfield)
    end if
!       initialize optional keywords
    do 10 i=1,maxf
            ttype(i)=' '
            tunit(i)=' '
10      continue

    call ftgkns(iunit,'TTYPE',1,maxf,ttype,nfind,status)
    call ftgkns(iunit,'TUNIT',1,maxf,tunit,nfind,status)

    if (status > 0)return

    call ftgkns(iunit,'TFORM',1,maxf,tform,nfind,status)
    if (status > 0 .or. nfind /= maxf)then
            status=232
            return
    end if

20      extnam=' '
    tstat=status
    call ftgkys(iunit,'EXTNAME',extnam,comm,status)
!       this keyword is not required, so ignore status
    if (status == 202)status =tstat
end
subroutine ftghdn(iunit,hdunum)
!
!*******************************************************************************
!
!! FTGHDN returns the number of the current header data unit.
!
!  The first HDU (the primary array) is number 1.

!       iunit   i  fortran unit number
!       hdunum  i  returned number of the current HDU
!
!       written by Wm Pence, HEASARC/GSFC, March, 1993

    integer iunit,hdunum

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    hdunum=chdu(bufnum(iunit))
end
subroutine ftghpr(iunit,maxdim,simple,bitpix,naxis,naxes, &
                      pcount,gcount,extend,status)
!
!*******************************************************************************
!
!! FTGHPR gets the required primary header or image extension keywords.
!
!       iunit   i  fortran unit number to use for reading
!       maxdim  i  maximum no. of dimensions to read; dimension of naxes
!       OUTPUT PARAMETERS:
!       simple  l  does file conform to FITS standard?
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       pcount  i  number of group parameters (usually 0)
!       gcount  i  number of random groups (usually 1 or 0)
!       extend  l  may extensions be present in the FITS file?
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,bitpix,naxis,naxes(*),pcount,gcount,blank,status
    integer maxdim,nblank
    logical simple,extend
    double precision fill

    call ftgphx(iunit,maxdim,simple,bitpix,naxis,naxes, &
          pcount,gcount,extend,fill,fill,blank,nblank,status)
end
subroutine ftghps(iunit,nkeys,pos,status)
!
!*******************************************************************************
!
!! FTGHPS gets the number of keywords and current position in the header.
!
!       Get Header Position
!       get the number of keywords in the header and the current position
!       in the header, i.e.,  the number of the next keyword record that
!       would be read.
!
!       iunit   i  Fortran I/O unit number
!       pos     i  current position in header (1 = beginning of header)
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Jan 1995

    integer iunit,nkeys,pos,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff

    if (status > 0)return

    ibuff=bufnum(iunit)
    nkeys=(hdend(ibuff)-hdstrt(ibuff,chdu(ibuff)))/80
    pos=(nxthdr(ibuff)-hdstrt(ibuff,chdu(ibuff)))/80+1
end
subroutine ftghsp(ounit,nexist,nmore,status)
!
!*******************************************************************************
!
!! FTGHSP returns the number of additional keywords that can fit in the header.
!
!       Get Header SPace
!       return the number of additional keywords that will fit in the header
!
!       ounit   i  Fortran I/O unit number
!       nexist  i  number of keywords already present in the CHU
!       nmore   i  number of additional keywords that will fit in header
!                 -1 indicates that there is no limit to the number of keywords
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,nexist,nmore,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff
    if (status > 0)return
    ibuff=bufnum(ounit)

    nexist=(hdend(ibuff)-hdstrt(ibuff,chdu(ibuff)))/80
    if (dtstrt(ibuff) < 0)then
!               the max size of the header has not been defined, so there
!               is no limit to the number of keywords which may be written.
            nmore=-1
    else
            nmore=(dtstrt(ibuff)-hdend(ibuff))/80-1
    end if
end
subroutine ftghtb(iunit,maxfld,ncols,nrows,nfield,ttype, &
                      tbcol,tform,tunit,extnam,status)
!
!*******************************************************************************
!
!! FTGHTB reads required standard header keywords from an ASCII table extension.
!
!       iunit   i  Fortran i/o unit number
!       maxfld  i  maximum no. of fields to read; dimension of ttype
!       OUTPUT PARAMETERS:
!       ncols   i  number of columns in the table
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       ttype   c  name of each field (array)
!       tbcol   i  beginning column of each field (array)
!       tform   c  Fortran-77 format of each field (array)
!       tunit   c  units of each field (array)
!       extnam  c  name of table (optional)
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,maxfld,ncols,nrows,nfield,status,tbcol(*)
    integer i,nfind,maxf,tstat
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam
    character comm*72

    call ftgttb(iunit,ncols,nrows,nfield,status)
    if (status > 0)return

    if (maxfld <= 0)then
            maxf=nfield
    else
            maxf=min(maxfld,nfield)
    end if

!       initialize optional keywords
    do i=1,maxf
            ttype(i)=' '
            tunit(i)=' '
    end do

    call ftgkns(iunit,'TTYPE',1,maxf,ttype,nfind,status)
    call ftgkns(iunit,'TUNIT',1,maxf,tunit,nfind,status)

    if (status > 0)return

    call ftgknj(iunit,'TBCOL',1,maxf,tbcol,nfind,status)
    if (status > 0 .or. nfind /= maxf)then
!               couldn't find the required TBCOL keywords
            status=231
    call ftpmsg('Required TBCOL keyword(s) not found in ASCII'// &
    ' table header (FTGHTB).')
            return
    end if

    call ftgkns(iunit,'TFORM',1,maxf,tform,nfind,status)
    if (status > 0 .or. nfind /= maxf)then
!               couldn't find the required TFORM keywords
            status=232
    call ftpmsg('Required TFORM keyword(s) not found in ASCII'// &
    ' table header (FTGHTB).')
            return
    end if

    extnam=' '
    tstat=status
    call ftgkys(iunit,'EXTNAME',extnam,comm,status)
!       this keyword is not required, so ignore 'keyword not found' status
    if (status == 202)status=tstat
end
subroutine ftgi1b(iunit,nvals,incre,chbuff,status)
!
!*******************************************************************************
!
!! FTGI1B reads an array of Integer*1 bytes from the input FITS file.
!
    integer nvals,incre,iunit,status,offset
    character chbuff(nvals)

!       iunit   i  fortran unit number
!       nvals   i  number of pixels in the i2vals array
!       incre   i  byte increment between values
!       chbuff  c*1 array of input byte values
!       status  i  output error status

    if (incre <= 1)then
            call ftgcbf(iunit,nvals,chbuff,status)
    else
!               offset is the number of bytes to move between each value
            offset=incre-1
            call ftgcbo(iunit,1,nvals,offset,chbuff,status)
    end if
end
subroutine ftgics(iunit,xrval,yrval,xrpix,yrpix,xinc,yinc,rot, &
                     type,status)
!
!*******************************************************************************
!
!! FTGICS reads the values of the celestial coordinate system keywords.
!
!       These values may be used as input to the subroutines that
!       calculate celestial coordinates. (FTXYPX, FTWLDP)

!       This routine assumes that the CHDU contains an image
!       with the RA type coordinate running along the first axis
!       and the DEC type coordinate running along the 2nd axis.

    double precision xrval,yrval,xrpix,yrpix,xinc,yinc,rot
    integer iunit,status,tstat
    character ( len = * ) type
    character comm*20,ctype*8

    if (status > 0)return

    call ftgkyd(iunit,'CRVAL1',xrval,comm,status)
    call ftgkyd(iunit,'CRVAL2',yrval,comm,status)

    call ftgkyd(iunit,'CRPIX1',xrpix,comm,status)
    call ftgkyd(iunit,'CRPIX2',yrpix,comm,status)

    call ftgkyd(iunit,'CDELT1',xinc,comm,status)
    call ftgkyd(iunit,'CDELT2',yinc,comm,status)

    call ftgkys(iunit,'CTYPE1',ctype,comm,status)

    if (status > 0)then
        call ftpmsg('FTGICS could not find all the required'// &
                    'celestial coordinate Keywords.')
        status=505
        return
    end if

    type=ctype(5:8)

    tstat=status
    call ftgkyd(iunit,'CROTA2',rot,comm,status)
    if (status > 0)then
!           CROTA2 is assumed to = 0 if keyword is not present
        status=tstat
        rot=0.
    end if
end
subroutine ftgiou(iounit,status)
!
!*******************************************************************************
!
!! FTGIOU gets an unallocated logical unit number.
!
    integer iounit,status

    if (status > 0)return
    iounit=0
    call ftxiou(iounit,status)
end
subroutine ftgkey(iunit,keynam,value,comm,status)
!
!*******************************************************************************
!
!! FTGKEY reads value and comment of a header keyword from the keyword buffer.
!
!       iunit   i  Fortran I/O unit number
!       keynam  c  name of keyword to be read
!       OUTPUT PARAMETERS:
!       value   c  output value of the keyword, if any
!       comm    c  output comment string, if any, of the keyword
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991

    integer iunit,status
    character ( len = * ) keynam,value,comm
    character*80 keybuf

    call ftgcrd(iunit,keynam,keybuf,status)
    if (status <= 0)then
!               parse the record to find value and comment strings
            call ftpsvc(keybuf,value,comm,status)
    end if
end
subroutine ftgknd(iunit,keywrd,nstart,nmax, &
                      dval,nfound,status)
!
!*******************************************************************************
!
!! FTGKND reads an array of real*8 values from header records.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       nstart  i  starting sequence number (usually 1)
!       nmax    i  number of keywords to read
!       OUTPUT PARAMETERS:
!       dval    d  array of output keyword values
!       nfound  i  number of keywords found
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd
    double precision dval(*)
    integer iunit,nstart,nmax,nfound,status,tstat
    integer nkeys,mkeys,i,ival,nend,namlen,indval
    logical vnull
    character inname*8,keynam*8
    character*80 rec,value,comm

    if (status > 0)return

!       for efficiency, we want to search just once through the header
!       for all the keywords which match the root.

    nfound=0
    nend=nstart+nmax-1
    inname=keywrd
    call ftupch(inname)

!       find the length of the root name
    namlen=0
    do 5 i=8,1,-1
            if (inname(i:i) /= ' ')then
                    namlen=i
                    go to 6
            end if
5       continue
6       if (namlen == 0)return

!       get the number of keywords in the header
    call ftghsp(iunit,nkeys,mkeys,status)

    vnull = .false.
    do i=3,nkeys
            call ftgrec(iunit,i,rec,status)
            if (status > 0)return
            keynam=rec(1:8)
            if (keynam(1:namlen) == inname(1:namlen))then

!                   try to interpret the remainder of the name as an integer
                tstat=status
                call ftc2ii(keynam(namlen+1:8),ival,status)
                if (status <= 0)then
                    if (ival <= nend .and. ival >= nstart)then
                        call ftpsvc(rec,value,comm,status)
                        indval=ival-nstart+1
                        call ftc2d(value,dval(indval),status)

                        if (status == 204)then
!                             value is undefined
                          status=0
                          vnull = .true.
                        end if

                        if (status > 0)then
         call ftpmsg('Error in FTGKND evaluating '//keynam// &
         ' as a Double: '//value)
                           return
                         else
                           nfound=max(nfound,indval)
                         end if
                    end if
                else
                    if (status == 407)then
                            status=tstat
                    else
                            return
                    end if
                end if
            end if
    end do

    if (status <= 0 .and. vnull)then
!           one or more values were undefined
        status = 204
    end if
end
subroutine ftgkne(iunit,keywrd,nstart,nmax, &
                      rval,nfound,status)
!
!*******************************************************************************
!
!! FTGKNE reads an array of real*4 values from header records.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       nstart  i  starting sequence number (usually 1)
!       nmax    i  number of keywords to read
!       OUTPUT PARAMETERS:
!       rval    r  array of output keyword values
!       nfound  i  number of keywords found
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd
    real rval(*)
    integer iunit,nstart,nmax,nfound,status,tstat
    integer nkeys,mkeys,i,ival,nend,namlen,indval
    logical vnull
    character inname*8,keynam*8
    character*80 rec,value,comm

    if (status > 0)return

!       for efficiency, we want to search just once through the header
!       for all the keywords which match the root.

    nfound=0
    nend=nstart+nmax-1
    inname=keywrd
    call ftupch(inname)

!       find the length of the root name
    namlen=0
    do 5 i=8,1,-1
            if (inname(i:i) /= ' ')then
                    namlen=i
                    go to 6
            end if
5       continue
6       if (namlen == 0)return

!       get the number of keywords in the header
    call ftghsp(iunit,nkeys,mkeys,status)

    vnull = .false.
    do i=3,nkeys
            call ftgrec(iunit,i,rec,status)
            if (status > 0)return
            keynam=rec(1:8)
            if (keynam(1:namlen) == inname(1:namlen))then

!                   try to interpret the remainder of the name as an integer
                tstat=status
                call ftc2ii(keynam(namlen+1:8),ival,status)
                if (status <= 0)then
                    if (ival <= nend .and. ival >= nstart)then
                        call ftpsvc(rec,value,comm,status)
                        indval=ival-nstart+1
                        call ftc2r(value,rval(indval),status)

                        if (status == 204)then
!                             value is undefined
                          status=0
                          vnull = .true.
                        end if

                        if (status > 0)then
         call ftpmsg('Error in FTGKNE evaluating '//keynam// &
         ' as a Real: '//value)
                           return
                         else
                           nfound=max(nfound,indval)
                         end if
                    end if
                else
                    if (status == 407)then
                            status=tstat
                    else
                            return
                    end if
                end if
            end if
    end do

    if (status <= 0 .and. vnull)then
!           one or more values were undefined
        status = 204
    end if
end
subroutine ftgknj(iunit,keywrd,nstart,nmax,intval, &
                      nfound,status)
!
!*******************************************************************************
!
!! FTGKNJ reads an array of integer values from  header records.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       nstart  i  starting sequence number (usually 1)
!       nmax    i  number of keywords to read
!       OUTPUT PARAMETERS:
!       intval  i  array of output keyword values
!       nfound  i  number of keywords found
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd
    integer intval(*)
    integer iunit,nstart,nmax,nfound,status,tstat
    integer nkeys,mkeys,i,ival,nend,namlen,indval
    logical vnull
    character inname*8,keynam*8
    character*80 rec,value,comm

    if (status > 0)return

!       for efficiency, we want to search just once through the header
!       for all the keywords which match the root.

    nfound=0
    nend=nstart+nmax-1
    inname=keywrd
    call ftupch(inname)

!       find the length of the root name
    namlen=0
    do 5 i=8,1,-1
            if (inname(i:i) /= ' ')then
                    namlen=i
                    go to 6
            end if
5       continue
6       if (namlen == 0)return

!       get the number of keywords in the header
    call ftghsp(iunit,nkeys,mkeys,status)

    vnull = .false.
    do i=3,nkeys
            call ftgrec(iunit,i,rec,status)
            if (status > 0)return
            keynam=rec(1:8)
            if (keynam(1:namlen) == inname(1:namlen))then

!                   try to interpret the remainder of the name as an integer
                tstat=status
                call ftc2ii(keynam(namlen+1:8),ival,status)
                if (status <= 0)then
                    if (ival <= nend .and. ival >= nstart)then
                        call ftpsvc(rec,value,comm,status)
                        indval=ival-nstart+1
                        call ftc2i(value,intval(indval),status)

                        if (status == 204)then
!                             value is undefined
                          status=0
                          vnull = .true.
                        end if

                        if (status > 0)then
         call ftpmsg('Error in FTGKNJ evaluating '//keynam// &
         ' as an integer: '//value)
                           return
                        else
                           nfound=max(nfound,indval)
                        end if
                    end if
                else
                    if (status == 407)then
                            status=tstat
                    else
                            return
                    end if
                end if
            end if
    end do

    if (status <= 0 .and. vnull)then
!           one or more values were undefined
        status = 204
    end if
end
subroutine ftgknl(iunit,keywrd,nstart,nmax,logval, &
                      nfound,status)
!
!*******************************************************************************
!
!! FTGKNL reads an array of logical values from  header records.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       nstart  i  starting sequence number (usually 1)
!       nmax    i  number of keywords to read
!       OUTPUT PARAMETERS:
!       logval  l  array of output keyword values
!       nfound  i  number of keywords found
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd
    logical logval(*), vnull
    integer iunit,nstart,nmax,nfound,status,tstat
    integer nkeys,mkeys,i,ival,nend,namlen,indval
    character inname*8,keynam*8
    character*80 rec,value,comm

    if (status > 0)return

!       for efficiency, we want to search just once through the header
!       for all the keywords which match the root.

    nfound=0
    nend=nstart+nmax-1
    inname=keywrd
    call ftupch(inname)

!       find the length of the root name
    namlen=0
    do 5 i=8,1,-1
            if (inname(i:i) /= ' ')then
                    namlen=i
                    go to 6
            end if
5       continue
6       if (namlen == 0)return

!       get the number of keywords in the header
    call ftghsp(iunit,nkeys,mkeys,status)

    vnull = .false.
    do i=3,nkeys
            call ftgrec(iunit,i,rec,status)
            if (status > 0)return
            keynam=rec(1:8)
            if (keynam(1:namlen) == inname(1:namlen))then

!                   try to interpret the remainder of the name as an integer
                tstat=status
                call ftc2ii(keynam(namlen+1:8),ival,status)
                if (status <= 0)then
                    if (ival <= nend .and. ival >= nstart)then
                        call ftpsvc(rec,value,comm,status)
                        indval=ival-nstart+1
                        call ftc2ll(value,logval(indval),status)
                        nfound=max(nfound,indval)

                        if (status == 204)then
!                             value is undefined
                          status=0
                          vnull = .true.
                        end if
                    end if
                else
                    if (status == 407)then
                            status=tstat
                    else
                            return
                    end if
                end if
            end if
    end do

    if (status <= 0 .and. vnull)then
!           one or more values were undefined
        status = 204
    end if
end
subroutine ftgkns(iunit,keywrd,nstart,nmax,strval,nfound, &
                      status)
!
!*******************************************************************************
!
!! FTGKNS reads an array of character string values from  header records.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       nstart  i  starting sequence number (usually 1)
!       nmax    i  number of keywords to read
!       OUTPUT PARAMETERS:
!       strval  c  array of output keyword values
!       nfound  i  number of keywords found
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,strval(*)
    integer iunit,nstart,nmax,nfound,status,tstat
    integer nkeys,mkeys,i,ival,nend,namlen,indval,ibuff
    logical vnull
    character inname*8,keynam*8
    character*80 value,comm

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(iunit)

!       for efficiency, we want to search just once through the header
!       for all the keywords which match the root.

    nfound=0
    nend=nstart+nmax-1
    inname=keywrd
    call ftupch(inname)

!       find the length of the root name
    namlen=0
    do i=8,1,-1
            if (inname(i:i) /= ' ')then
                    namlen=i
                    exit
            end if
    end do
    if (namlen == 0)return

!       get the number of keywords in the header
    call ftghsp(iunit,nkeys,mkeys,status)

    vnull = .false.
    do i=3,nkeys
            call ftgrec(iunit,i,value,status)
            if (status > 0)return
            keynam=value(1:8)
            if (keynam(1:namlen) == inname(1:namlen))then

!                   try to interpret the remainder of the name as an integer
                tstat=status
                call ftc2ii(keynam(namlen+1:8),ival,status)
                if (status <= 0)then
                  if (ival <= nend .and. ival >= nstart)then

!                       OK, this looks like a valid keyword; Reset the
!                       next-header-keyword pointer by one record, then
!                       call ftgkys to read it. (This does  support
!                       long continued string values)

                    nxthdr(ibuff)=nxthdr(ibuff)-80
                    indval=ival-nstart+1
                    call ftgkys(iunit,keynam,strval(indval), &
                                comm,status)

                    if (status == 204)then
!                         value is undefined
                      status=0
                      vnull = .true.
                    end if

                    nfound=max(nfound,indval)
                  end if
                else
                    if (status == 407)then
                            status=tstat
                    else
                            return
                    end if
                end if
            end if
    end do

    if (status <= 0 .and. vnull)then
!           one or more values were undefined
        status = 204
    end if
end
subroutine ftgkyd(iunit,keywrd,dval,comm,status)
!
!*******************************************************************************
!
!! FTGKYD reads a double precision value and comment string from a header record.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       dval    i  output keyword value
!       comm    c  output keyword comment
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer iunit,status
    character value*35
    double precision dval

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

!       convert character string to double precision
!       datatype conversion will be performed if necessary and if possible
    call ftc2d(value,dval,status)
end
subroutine ftgkye(iunit,keywrd,rval,comm,status)
!
!*******************************************************************************
!
!! FTGKYE reads a real*4 value and the comment string from a header record.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       rval    r  output keyword value
!       comm    c  output keyword comment
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer iunit,status
    character value*35
    real rval

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

!       convert character string to real
!       datatype conversion will be performed if necessary and if possible
    call ftc2r(value,rval,status)
end
subroutine ftgkyj(iunit,keywrd,intval,comm,status)
!
!*******************************************************************************
!
!! FTGKYJ reads an integer value and the comment string from a header record.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       intval  i  output keyword value
!       comm    c  output keyword comment
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer iunit,intval,status
    character value*35

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

!       convert character string to integer
!       datatype conversion will be performed if necessary and if possible
    call ftc2i(value,intval,status)
end
subroutine ftgkyl(iunit,keywrd,logval,comm,status)
!
!*******************************************************************************
!
!! FTGKYL reads a logical value and the comment string from a header record.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       logval  l  output keyword value
!       comm    c  output keyword comment
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer iunit,status
    character value*20
    logical logval

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

!       convert character string to logical
    call ftc2l(value,logval,status)
end
subroutine ftgkyn(iunit,nkey,keynam,value,comm,status)
!
!*******************************************************************************
!
!! FTGKYN reads value and comment of the NKEYth header record.
!
!       This routine is useful for reading the entire header, one
!       record at a time.

!       iunit   i  Fortran I/O unit number
!       nkey    i  sequence number (starting with 1) of the keyword to read
!       OUTPUT PARAMETERS:
!       keynam  c  output name of the keyword
!       value   c  output value of the keyword, if any
!       comm    c  output comment string, if any, of the keyword
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,nkey,status
    character ( len = * ) keynam,value,comm
    character keybuf*80,arec*8

    if (status > 0)return

    call ftgrec(iunit,nkey,keybuf,status)
    if (status > 0)return

    keynam=keybuf(1:8)

!       parse the value and comment fields from the record
    call ftpsvc(keybuf,value,comm,status)
    if (status > 0)return

!       Test that keyword name contains only valid characters.
!       This also serves as a check in case there was no END keyword and
!       program continues to read on into the data unit
    call fttkey(keybuf(1:8),status)
    if (status > 0)then
        write(arec,1000)nkey
1000        format(i8)
        call ftpmsg('Name of header keyword number'//arec// &
       ' contains illegal character(s):')
        call ftpmsg(keybuf)

!          see if we are at the beginning of FITS logical record
       if (nkey-1 == (nkey-1)/36*36 .and. nkey > 1)then
         call ftpmsg('(This may indicate a missing END keyword).')
       end if
    end if
end
subroutine ftgkys(iunit,keywrd,strval,comm,status)
!
!*******************************************************************************
!
!! FTGKYS reads a character string value and comment from a header record.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       strval  c  output keyword value
!       comm    c  output keyword comment
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!       modified 6/93 to support long strings which are continued
!       over several keywords.  A string may be continued by putting
!       a backslash as the last non-blank character in the keyword string,
!       then continuing the string in the next keyword which must have
!       a blank keyword name.
!       Modified 9/94 to support the new OGIP continuation convention

    character ( len = * ) keywrd,comm,strval
    integer status,iunit
    character value*70, comm2*70, bslash*1
    integer clen,i,bspos,lenval

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

!       convert character string to unquoted string
    call ftc2s(value,strval,status)

    if (status > 0)return

    clen=len(strval)

!       is last character a backslash or & ?
!       have to use 2 \\'s because the SUN compiler treats 1 \ as an escape
    bslash='\\'
    do i=70,1,-1
            if (value(i:i) /= ' ' .and. value(i:i)/='''')then
                    if (value(i:i) == bslash .or. &
                        value(i:i) == '&')then
!                               have to subtract 1 due to the leading quote char
                            bspos=i-1
                            go to 20
                    end if
!                       no continuation character, so just return
                    return
            end if
    end do
!       value field was blank, so just return
    return

!       try to get the string continuation, and new comment string
20      call ftgnst(iunit,value,lenval,comm2,status)
    if (lenval == 0)return

    if (bspos <= clen)then
            strval(bspos:)=value(1:lenval)
            bspos=bspos+lenval-1
    end if

    if (comm2 /= ' ')comm=comm2

!       see if there is another continuation line
    if (value(lenval:lenval) == bslash .or. &
        value(lenval:lenval) == '&')go to 20
end
subroutine ftgkyt(iunit,keywrd,jval,dval,comm,status)
!
!*******************************************************************************
!
!! FTGKYT reads an integer value and fractional parts of a keyword value.
!
!
!  The information is read, along with the comment string, from a header record
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       jval    i  output integer part of keyword value
!       dval    d  output fractional part of keyword value
!       comm    c  output keyword comment
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, Sept 1992

    character ( len = * ) keywrd,comm
    integer iunit,jval,status,i,dot
    double precision dval
    character value*35
    logical ed

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

!       read keyword in straight forward way first:
!       just convert character string to double precision
!       datatype conversion will be performed if necessary and if possible
    call ftc2d(value,dval,status)
    jval=dval
    if (jval >= 0)then
            dval=dval-jval
    else
            dval=dval+jval
    end if

!       now see if we have to read the fractional part again, this time
!       with more precision

!       find the decimal point, if any, and look for a D or E
    dot=0
    ed=.false.
    do i=1,35
        if (value(i:i) == '.')dot=i
        if (value(i:i) == 'E' .or. value(i:i) == 'D')ed=.true.
    end do

    if (.not. ed .and. dot > 0)then
!           convert fractional part to double precision
        call ftc2d(value(dot:),dval,status)
    end if

end
subroutine ftgmsg(text)
!
!*******************************************************************************
!
!! FTGMSG gets an error message from the stack and shift the stack up.
!
    character ( len = * ) text
    call ftxmsg(-1,text)
end
subroutine ftgnst(iunit,value,lenval,comm,status)
!
!*******************************************************************************
!
!! FTGNST gets the next string keyword.
!
!       see if the next keyword in the header is the continuation
!       of a long string keyword, and if so, return the value string,
!       the number of characters in the string, and the associated comment
!       string.

!       value  c  returned value of the string continuation
!       lenval i  number of non-blank characters in the continuation string
!       comm   C  value of the comment string, if any, in this keyword.

    character ( len = * ) value,comm
    integer iunit,lenval,status

    integer i,length,tstat,nkeys,nextky
    character record*80, strval*70

    if (status > 0)return

    tstat=status
    value=' '
    comm=' '
    lenval=0

!       get current header position
    call ftghps(iunit,nkeys,nextky,status)

!       get the next keyword record
    if (nextky <= nkeys)then
        call ftgrec(iunit,nextky,record,status)
    else
!           positioned at end of header, so there is no next keyword to read
        return
    end if

!       does this appear to be a continuation keyword (=blank keyword name
!       or CONTINUE)?
    if (record(1:10) /= ' ' .and. record(1:10) /= &
       'CONTINUE  ')return

!       return if record is blank
    if (record == ' ')return

!       set a dummy keyword name
    record(1:10)='DUMMYKEY= '

!       parse the record to get the value string and comment
    call ftpsvc(record,strval,comm,status)

!       convert character string to unquoted string
    call ftc2s(strval,value,status)
    if (status > 0)then
!               this must not be a continuation card; reset status and messages
            status=tstat
            call ftcmsg
            value=' '
            comm=' '
            return
    end if

    length=len(value)
    do i=length,1,-1
            if (value(i:i) /= ' ')then
                    lenval=i
                    return
            end if
    end do

end
subroutine ftgnxk(iunit,inclst,ninc,exclst,nexc,card,status)
!
!*******************************************************************************
!
!! FTGNXK returns the next keyword in INCLIST and not in EXCLIST.
!
!  The keyword matches one of the names in inclist
!    but does not match any of the names in exclist.  The search
!    goes from the current position to the end of the header, only.
!    Wild card characters may be used in the name lists ('*', '?' and '#').

!       iunit   i  Fortran I/O unit number
!       inclist c  list of included keyword names
!       ninc    i number of names in inclist
!       exclist c list of excluded keyword names
!       nexc    i number of names in exclist
!       OUTPUT PARAMETERS:
!       card    c  first matching 80 character card image
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, January 1997

    integer iunit,ninc,nexc,status,ii,jj
    character ( len = * ) inclst(*),exclst(*),card
    character*80 keybuf
    logical casesn,match,exact

    card=' '
    if (status > 0)return
    casesn=.false.

10      call ftgcrd(iunit,'*',keybuf,status)
    if (status <= 0)then
      do 30 ii = 1, ninc
        call ftcmps(inclst(ii),keybuf(1:8),casesn,match,exact)
        if (match)then
          do 20 jj = 1,nexc
            call ftcmps(exclst(jj),keybuf(1:8),casesn,match,exact)
!               reject this card if in exclusion list
            if (match)go to 10
20            continue

!             keyword is not excluded, so return it
          card = keybuf
          return
        end if
30        continue

!         didn't match, so go back to read next keyword
      go to 10
    end if

!       failed to read next keyword (probably hit end of header)
end
subroutine ftgpfb(iunit,group,felem,nelem, &
                      array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGPFB reads an array of byte values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will have the corresponding element of
!       FLGVAL set equal to .true.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       array   b  returned array of values that were read
!       flgval  l  set to .true. if the corresponding element is undefined
!       anynul  l  set to .true. if any returned elements are undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    character nulval,array(*)
    logical anynul,flgval(*)
    integer i

    do 10 i=1,nelem
            flgval(i)=.false.
10      continue

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgclb(iunit,2,row,felem,nelem,1,2,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpfd(iunit,group,felem,nelem, &
                      array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGPFD reads an array of r*8 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will have the corresponding element of
!       FLGVAL set equal to .true.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       array   d  returned array of values that were read
!       flgval  l  set to .true. if the corresponding element is undefined
!       anynul  l  set to .true. if any returned elements are undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    double precision nulval,array(*)
    logical anynul,flgval(*)
    integer i

    do 10 i=1,nelem
            flgval(i)=.false.
10      continue

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgcld(iunit,2,row,felem,nelem,1,2,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpfe(iunit,group,felem,nelem, &
                      array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGPFE reads an array of r*4 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will have the corresponding element of
!       FLGVAL set equal to .true.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       array   r  returned array of values that were read
!       flgval  l  set to .true. if the corresponding element is undefined
!       anynul  l  set to .true. if any returned elements are undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    real nulval,array(*)
    logical anynul,flgval(*)
    integer i

    do 10 i=1,nelem
            flgval(i)=.false.
10      continue

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgcle(iunit,2,row,felem,nelem,1,2,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpfi(iunit,group,felem,nelem, &
                      array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGPFI reads an array of I*2 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will have the corresponding element of
!       FLGVAL set equal to .true.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       array   i*2  returned array of values that were read
!       flgval  l  set to .true. if the corresponding element is undefined
!       anynul  l  set to .true. if any returned elements are undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    integer*2 nulval,array(*)
    logical anynul,flgval(*)
    integer i

    do 10 i=1,nelem
            flgval(i)=.false.
10      continue

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgcli(iunit,2,row,felem,nelem,1,2,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpfj(iunit,group,felem,nelem, &
                      array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGPFJ reads an array of I*4 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will have the corresponding element of
!       FLGVAL set equal to .true.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       array   i  returned array of values that were read
!       flgval  l  set to .true. if the corresponding element is undefined
!       anynul  l  set to .true. if any returned elements are undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    integer nulval,array(*)
    logical anynul,flgval(*)
    integer i

    do i=1,nelem
            flgval(i)=.false.
    end do

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgclj(iunit,2,row,felem,nelem,1,2,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgphx(iunit,maxdim,simple,bitpix,naxis,naxes,pcount &
                 ,gcount,extend,bscale,bzero,blank,nblank,status)
!
!*******************************************************************************
!
!! FTGPHX gets the main primary header keywords that define the array structure.
!
!       iunit   i  fortran unit number to use for reading
!       maxdim  i  maximum no. of dimensions to read; dimension of naxes
!       OUTPUT PARAMETERS:
!       simple  l  does file conform to FITS standard?
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       pcount  i  number of group parameters (usually 0)
!       gcount  i  number of random groups (usually 1 or 0)
!       extend  l  may extensions be present in the FITS file?
!       bscale  d  scaling factor
!       bzero   d  scaling zero point
!       blank   i  value used to represent undefined pixels
!       nblank  i  number of trailing blank keywords immediately before the END
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,maxdim,bitpix,naxis
    integer naxes(*),pcount,gcount,blank,status,tstat
    logical simple,extend,unknow
    character keynam*8,value*20,lngval*40,comm*72,extn*4,keybuf*80
    double precision bscale,bzero
    integer nkey,nblank,i,ibuff,taxes,maxd

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check that the first keyword is valid
    call ftgrec(iunit,1,keybuf,status)

    keynam=keybuf(1:8)
!       parse the value and comment fields from the record
    call ftpsvc(keybuf,value,comm,status)

    if (status > 0)go to 900

    simple=.true.
    unknow=.false.
    if (chdu(ibuff) == 1)then
        if (keynam == 'SIMPLE')then
            if (value == 'F')then
!                       this is not a simple FITS file; try to process it anyway
                    simple=.false.
            else if (value /= 'T')then
!                       illegal value for the SIMPLE keyword
                    status=220

     if (keybuf(9:10) /= '= ')then
       call ftpmsg('The SIMPLE keyword is missing "= " in '// &
       'columns 9-10.')
     else
       call ftpmsg('The SIMPLE keyword value is illegal:'//value &
       // '.  It must equal T or F:')
     end if

                    call ftpmsg(keybuf)
            end if
        else
            status=221
    call ftpmsg('First keyword of the file is not SIMPLE: '//keynam)
            call ftpmsg(keybuf)
            go to 900
        end if
    else
         if (keynam == 'XTENSION')then
            if (value(2:9) /= 'IMAGE   ' .and. &
                value(2:9) /= 'IUEIMAGE')then
!                    I don't know what type of extension this is, but press on
                 unknow=.true.

     if (keybuf(9:10) /= '= ')then
       call ftpmsg('The XTENSION keyword is missing "= " in '// &
       'columns 9-10.')
     else
       call ftpmsg('This is not an IMAGE extension: '//value)
     end if

                 call ftpmsg(keybuf)
             end if
         else
             status=225
             write(extn,1000)chdu(ibuff)
1000             format(i4)
             call ftpmsg('First keyword in extension '//extn// &
             ' was not XTENSION: '//keynam)
             call ftpmsg(keybuf)
         end if
    end if
    if (status > 0)go to 900

!       check that BITPIX is the second keyword
    call ftgrec(iunit,2,keybuf,status)

    keynam=keybuf(1:8)
!       parse the value and comment fields from the record
    call ftpsvc(keybuf,value,comm,status)

    if (status > 0)go to 900
    if (keynam /= 'BITPIX')then
            status=222
    call ftpmsg('Second keyword was not BITPIX: '//keynam)
            call ftpmsg(keybuf)
            go to 900
    end if
!       convert character string to integer
    call ftc2ii(value,bitpix,status)
    if (status > 0)then
!         bitpix value must be an integer
      if (keybuf(9:10) /= '= ')then
         call ftpmsg('BITPIX keyword is missing "= "'// &
        ' in columns 9-10.')
      else
          call ftpmsg('Value of BITPIX is not an integer: '//value)
      end if
      call ftpmsg(keybuf)
      status=211
      go to 900
    end if

!       test that bitpix has a legal value
    call fttbit(bitpix,status)
    if (status > 0)then
            call ftpmsg(keybuf)
            go to 900
    end if

!       check that the third keyword is NAXIS
    call ftgtkn(iunit,3,'NAXIS',naxis,status)
    if (status == 208)then
!               third keyword was not NAXIS
            status=223
    else if (status == 209)then
!               NAXIS value was not an integer
            status=212
    end if
    if (status > 0)go to 900

    if (maxdim <= 0)then
            maxd=naxis
    else
            maxd=min(maxdim,naxis)
    end if

    do i=1,naxis
!               construct keyword name
            call ftkeyn('NAXIS',i,keynam,status)
!               attempt to read the keyword
            call ftgtkn(iunit,3+i,keynam,taxes,status)
            if (status > 0)then
                    status=224
                    go to 900
            else if (taxes < 0)then
!                       NAXISn keywords must not be negative
                    status=213
                    go to 900
            else if (i <= maxd)then
                    naxes(i)=taxes
            end if
    end do

!       now look for other keywords of interest: bscale, bzero, blank, and END
!       and pcount, gcount, and extend
15      bscale=1.
    bzero=0.
    pcount=0
    gcount=1
    extend=.false.
!       choose a special value to represent the absence of a blank value
    blank=123454321

    nkey=3+naxis
18      nblank=0
20      nkey=nkey+1
    tstat=status
    call ftgrec(iunit,nkey,keybuf,status)
    if (status > 0)then
!               first, check for normal end-of-header status, and reset to 0
            if (status == 203)status=tstat
!               if we hit the end of file, then set status = no END card found
            if (status == 107)then
                   status=210
                   call ftpmsg('FITS header has no END keyword!')
            end if
            go to 900
    end if
    keynam=keybuf(1:8)
    comm=keybuf(9:80)

    if (keynam == 'BSCALE')then
!               convert character string to floating pt.
            call ftpsvc(keybuf,lngval,comm,status)
            call ftc2dd(lngval,bscale,status)
            if (status > 0)then
                 call ftpmsg('Error reading BSCALE keyword value'// &
                 ' as a Double:'//lngval)
            end if
    else if (keynam == 'BZERO')then
!               convert character string to floating pt.
            call ftpsvc(keybuf,lngval,comm,status)
            call ftc2dd(lngval,bzero,status)
            if (status > 0)then
                 call ftpmsg('Error reading BZERO keyword value'// &
                 ' as a Double:'//lngval)
            end if
    else if (keynam == 'BLANK')then
!               convert character string to integer
            call ftpsvc(keybuf,value,comm,status)
            call ftc2ii(value,blank,status)
            if (status > 0)then
                 call ftpmsg('Error reading BLANK keyword value'// &
                 ' as an integer:'//value)
            end if
    else if (keynam == 'PCOUNT')then
!               convert character string to integer
            call ftpsvc(keybuf,value,comm,status)
            call ftc2ii(value,pcount,status)
            if (status > 0)then
                 call ftpmsg('Error reading PCOUNT keyword value'// &
                 ' as an integer:'//value)
            end if
    else if (keynam == 'GCOUNT')then
!               convert character string to integer
            call ftpsvc(keybuf,value,comm,status)
            call ftc2ii(value,gcount,status)
            if (status > 0)then
                 call ftpmsg('Error reading GCOUNT keyword value'// &
                 ' as an integer:'//value)
            end if
    else if (keynam == 'EXTEND')then
!               convert character string to logical
            call ftpsvc(keybuf,value,comm,status)
            call ftc2ll(value,extend,status)
            if (status > 0)then
                 call ftpmsg('Error reading EXTEND keyword value'// &
                 ' as a Logical:'//value)
             end if
    else if (keynam == ' ' .and. comm == ' ')then
!               need to ignore trailing blank records before the END card
            nblank=nblank+1
            go to 20
    else if (keynam == 'END')then
            go to 900
    end if
    if (status > 0)go to 900
    go to 18

900     continue

    if (status > 0)then
      if (chdu(ibuff) == 1)then
        call ftpmsg('Failed to parse the required keywords in '// &
         'the Primary Array header ')
      else
        call ftpmsg('Failed to parse the required keywords in '// &
         'the Image Extension header (FTGPHX).')
      end if

    else if (unknow)then
!           set status if this was an unknown type of extension
        status=233
    end if
end
subroutine ftgprh(iunit,simple,bitpix,naxis,naxes, &
                      pcount,gcount,extend,status)
!
!*******************************************************************************
!
!! FTGPRH is obsolete; call FTGHPR instead.
!

    integer iunit,bitpix,naxis,naxes(*),pcount,gcount,blank,status
    integer nblank
    logical simple,extend
    double precision fill

    call ftgphx(iunit,0,simple,bitpix,naxis,naxes, &
          pcount,gcount,extend,fill,fill,blank,nblank,status)
end
subroutine ftgpvb(iunit,group,felem,nelem,nulval, &
                      array,anynul,status)
!
!*******************************************************************************
!
!! FTGPVB reads an array of byte values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will be set equal to NULVAL, unless NULVAL=0
!       in which case no checking for undefined values will be performed.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       nulval  b  the value to be assigned to undefined pixels
!       array   b  returned array of values that were read
!       anynul  l  set to .true. if any returned elements were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    character nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgclb(iunit,2,row,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpvd(iunit,group,felem,nelem,nulval, &
                      array,anynul,status)
!
!*******************************************************************************
!
!! FTGPVD reads an array of r*8 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will be set equal to NULVAL, unless NULVAL=0
!       in which case no checking for undefined values will be performed.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       nulval  b  the value to be assigned to undefined pixels
!       array   b  returned array of values that were read
!       anynul  l  set to .true. if any returned elements were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    double precision nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgcld(iunit,2,row,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpve(iunit,group,felem,nelem,nulval, &
                      array,anynul,status)
!
!*******************************************************************************
!
!! FTGPVE reads an array of r*4 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will be set equal to NULVAL, unless NULVAL=0
!       in which case no checking for undefined values will be performed.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       nulval  r  the value to be assigned to undefined pixels
!       array   r  returned array of values that were read
!       anynul  l  set to .true. if any returned elements were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    real nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgcle(iunit,2,row,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpvi(iunit,group,felem,nelem,nulval, &
                      array,anynul,status)
!
!*******************************************************************************
!
!! FTGPVI reads an array of i*2 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will be set equal to NULVAL, unless NULVAL=0
!       in which case no checking for undefined values will be performed.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       nulval  i*2  the value to be assigned to undefined pixels
!       array   i*2  returned array of values that were read
!       anynul  l  set to .true. if any returned elements were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    integer*2 nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgcli(iunit,2,row,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgpvj(iunit,group,felem,nelem,nulval, &
                      array,anynul,status)
!
!*******************************************************************************
!
!! FTGPVJ reads an array of i*4 values from the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).
!       Undefined elements will be set equal to NULVAL, unless NULVAL=0
!       in which case no checking for undefined values will be performed.
!       ANYNUL is return with a value of .true. if any pixels were undefined.

!       iunit   i  Fortran unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be read (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be read
!       nulval  i  the value to be assigned to undefined pixels
!       array   i  returned array of values that were read
!       anynul  l  set to .true. if any returned elements were undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,group,felem,nelem,status,row
    integer nulval,array(*)
    logical anynul,flgval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(1,group)
    call ftgclj(iunit,2,row,felem,nelem,1,1,nulval, &
        array,flgval,anynul,status)
end
subroutine ftgrec(iunit,nrec,record,status)
!
!*******************************************************************************
!
!! FTGREC reads the Nth 80-byte header record.
!
!       This routine is useful for reading the entire header, one
!       record at a time.

!       iunit   i  Fortran I/O unit number
!       nrec    i  sequence number (starting with 1) of the record to read
!       OUTPUT PARAMETERS:
!       record  c  output 80-byte record
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,nrec,status
    character ( len = * ) record

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,nbyte,endhd
    character arec*8

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(iunit)

!       calculate byte location of the record, and check if it is legal
    nbyte=hdstrt(ibuff,chdu(ibuff))+(nrec-1)*80

!   endhd=(hdend(ibuff)/2880+1)*2880
!       modified this on 4 Nov 1994 to allow for blanks before the END keyword
endhd=max(hdend(ibuff),dtstrt(ibuff)-2880)

    if (nrec == 0)then
!               simply move to the beginning of the header
!               update the keyword pointer position
            nxthdr(ibuff)=nbyte+80
            record=' '
            return
    else if (nbyte > endhd .or. nrec < 0)then
!               header record number is out of bounds
            status=203
            write(arec,1000)nrec
1000            format(i8)
            call ftpmsg('Cannot get Keyword number '//arec//'.'// &
            '  It does not exist.')
            go to 100
    end if

!       position the I/O pointer to the appropriate header keyword
    call ftmbyt(iunit,nbyte,.false.,status)

!       read the 80 byte record
    call ftgcbf(iunit,80,record,status)
    if (status > 0)then
            write(arec,1000)nrec
            call ftpmsg('FTGREC could not read header keyword'// &
              ' number '//arec//'.')
            return
    end if

!       update the keyword pointer position
    nxthdr(ibuff)=nbyte+80

100     continue
end
subroutine ftgsfb(iunit,colnum,naxis,naxes,blc,trc,inc, &
    array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGSFB reads a subsection of byte data from an image or a table column.
!
!  Returns an associated array of null value flags.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       array   i  array of data values that are read from the FITS file
!       flgval  l  set to .true. if corresponding array element is undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    character array(*),nulval
    logical anynul,anyf,flgval(*)

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSFB ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=2
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSFB, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgclb(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),flgval(i1),anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsfd(iunit,colnum,naxis,naxes,blc,trc,inc, &
    array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGSFD reads a subsection of double precision data from an image or table column.
!
!  Returns an associated array of null value flags.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       array   i  array of data values that are read from the FITS file
!       flgval  l  set to .true. if corresponding array element is undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    double precision array(*),nulval
    logical anynul,anyf,flgval(*)

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSFD ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=2
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSFD, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgcld(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),flgval(i1),anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsfe(iunit,colnum,naxis,naxes,blc,trc,inc, &
    array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGSFE reads a subsection of real data from an image or table column.
!
!  Returns an associated array of null value flags.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       array   i  array of data values that are read from the FITS file
!       flgval  l  set to .true. if corresponding array element is undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    real array(*),nulval
    logical anynul,anyf,flgval(*)

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSFE ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=2
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSFE, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgcle(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),flgval(i1),anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsfi(iunit,colnum,naxis,naxes,blc,trc,inc, &
    array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGSFI reads a subsection of integer*2 data from an image or table column.
!
!  Returns an associated array of null value flags.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       array   i  array of data values that are read from the FITS file
!       flgval  l  set to .true. if corresponding array element is undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    integer*2 array(*),nulval
    logical anynul,anyf,flgval(*)

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSFI ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=2
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSFI, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgcli(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),flgval(i1),anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsfj(iunit,colnum,naxis,naxes,blc,trc,inc, &
    array,flgval,anynul,status)
!
!*******************************************************************************
!
!! FTGSFJ reads a subsection of integer*4 data from an image or table column.
!
!  Returns an associated array of null value flags.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       array   i  array of data values that are read from the FITS file
!       flgval  l  set to .true. if corresponding array element is undefined
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    integer array(*),nulval
    logical anynul,anyf,flgval(*)

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSFJ ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=2
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSFJ, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgclj(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),flgval(i1),anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsvb(iunit,colnum,naxis,naxes,blc,trc,inc, &
    nulval,array,anynul,status)
!
!*******************************************************************************
!
!! FTGSVB reads a subsection of byte data from an image or a table column.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       nulval  i  value that undefined pixels will be set to
!       array   i  array of data values that are read from the FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    character array(*),nulval
    logical anynul,anyf

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    logical ldummy
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSVB ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=1
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSVB, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgclb(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),ldummy,anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsvd(iunit,colnum,naxis,naxes,blc,trc,inc, &
    nulval,array,anynul,status)
!
!*******************************************************************************
!
!! FTGSVD reads a subsection of double precision data from an image or table column.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       nulval  i  value that undefined pixels will be set to
!       array   i  array of data values that are read from the FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    double precision array(*),nulval
    logical anynul,anyf

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    logical ldummy
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSVD ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=1
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSVD, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgcld(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),ldummy,anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsve(iunit,colnum,naxis,naxes,blc,trc,inc, &
    nulval,array,anynul,status)
!
!*******************************************************************************
!
!! FTGSVE reads a subsection of real data from an image or table column.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       nulval  i  value that undefined pixels will be set to
!       array   i  array of data values that are read from the FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    real array(*),nulval
    logical anynul,anyf

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    logical ldummy
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSVE ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=1
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSVE, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgcle(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),ldummy,anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsvi(iunit,colnum,naxis,naxes,blc,trc,inc, &
    nulval,array,anynul,status)
!
!*******************************************************************************
!
!! FTGSVI reads a subsection of integer*2 data from an image or a table column.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       nulval  i  value that undefined pixels will be set to
!       array   i  array of data values that are read from the FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    integer*2 array(*),nulval
    logical anynul,anyf

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    logical ldummy
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSVI ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=1
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSVI, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
    end do

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgcli(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),ldummy,anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgsvj(iunit,colnum,naxis,naxes,blc,trc,inc, &
    nulval,array,anynul,status)
!
!*******************************************************************************
!
!! FTGSVJ reads a subsection of integer*4 data from an image or table column.
!
!       iunit   i  fortran unit number
!       colnum  i  number of the column to read from
!       naxis   i  number of dimensions in the FITS array
!       naxes   i  size of each dimension.
!       blc     i  'bottom left corner' of the subsection to be read
!       trc     i  'top right corner' of the subsection to be read
!       inc     i  increment to be applied in each dimension
!       nulval  i  value that undefined pixels will be set to
!       array   i  array of data values that are read from the FITS file
!       anynul  l  set to .true. if any of the returned values are undefined
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1993

    integer iunit,colnum,naxis,naxes(*),blc(*),trc(*),inc(*),status
    integer array(*),nulval
    logical anynul,anyf

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer i,i1,i2,i3,i4,i5,i6,i7,i8,i9,row,rstr,rstp,rinc
    integer str(9),stp(9),incr(9),dsize(10)
    integer felem,nelem,nultyp,ninc,ibuff,numcol
    logical ldummy
    character caxis*20

!       this routine is set up to handle a maximum of nine dimensions

    if (status > 0)return

    if (naxis < 1 .or. naxis > 9)then
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTGSVJ ' &
            //'is illegal.')
            return
    end if

!       if this is a primary array, then the input COLNUM parameter should
!       be interpreted as the row number, and we will alway read the image
!       data from column 2 (any group parameters are in column 1).

    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 0)then
!               this is a primary array, or image extension
            if (colnum == 0)then
                rstr=1
                rstp=1
            else
                rstr=colnum
                rstp=colnum
            end if
            rinc=1
            numcol=2
    else
!               this is a table, so the row info is in the (naxis+1) elements
            rstr=blc(naxis+1)
            rstp=trc(naxis+1)
            rinc=inc(naxis+1)
            numcol=colnum
    end if

    nultyp=1
    anynul=.false.
    i1=1
    do 5 i=1,9
            str(i)=1
            stp(i)=1
            incr(i)=1
            dsize(i)=1
5       continue
    do 10 i=1,naxis
            if (trc(i) < blc(i))then
                    status=321
                    write(caxis,1001)i
    call ftpmsg('In FTGSVJ, the range specified for axis '// &
    caxis(19:20)//' has the start greater than the end.')
                    return
            end if
            str(i)=blc(i)
            stp(i)=trc(i)
            incr(i)=inc(i)
            dsize(i+1)=dsize(i)*naxes(i)
10      continue

    if (naxis == 1 .and. naxes(1) == 1)then
!               This is not a vector column, so read all the rows at once
            nelem=(rstp-rstr)/rinc+1
            ninc=rinc
            rstp=rstr
    else
!               have to read each row individually, in all dimensions
            nelem=(stp(1)-str(1))/inc(1)+1
            ninc=incr(1)
    end if

    do 100 row=rstr,rstp,rinc
     do 90 i9=str(9),stp(9),incr(9)
      do 80 i8=str(8),stp(8),incr(8)
       do 70 i7=str(7),stp(7),incr(7)
        do 60 i6=str(6),stp(6),incr(6)
         do 50 i5=str(5),stp(5),incr(5)
          do 40 i4=str(4),stp(4),incr(4)
           do 30 i3=str(3),stp(3),incr(3)
            do 20 i2=str(2),stp(2),incr(2)

    felem=str(1)+(i2-1)*dsize(2)+(i3-1)*dsize(3)+(i4-1)*dsize(4) &
    +(i5-1)*dsize(5)+(i6-1)*dsize(6)+(i7-1)*dsize(7) &
    +(i8-1)*dsize(8)+(i9-1)*dsize(9)

    call ftgclj(iunit,numcol,row,felem,nelem,ninc, &
    nultyp,nulval,array(i1),ldummy,anyf,status)
    if (status > 0)return
    if (anyf)anynul=.true.
    i1=i1+nelem

20              continue
30             continue
40            continue
50           continue
60          continue
70         continue
80        continue
90       continue
100     continue
end
subroutine ftgtbb(iunit,frow,fchar,nchars,value,status)
!
!*******************************************************************************
!
!! FTGTBB reads a consecutive string of bytes from an ascii or binary table.
!
!       This will span multiple rows of the table if NCHARS+FCHAR is
!       greater than the length of a row.

!       iunit   i  fortran unit number
!       frow    i  starting row number (1st row = 1)
!       fchar   i  starting character/byte in the row to read (1st character=1)
!       nchars  i  number of characters/bytes to read (can span multiple rows)
!       value   i  returned string of bytes
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1991

    integer iunit,frow,fchar,nchars,status
    integer value(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bstart

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check for errors
    if (nchars <= 0)then
!               zero or negative number of character requested
            return
    else if (frow < 1)then
!               error: illegal first row number
            status=307
            return
    else if (fchar < 1)then
!               error: illegal starting character
            status=308
            return
    end if

!       move the i/o pointer to the start of the sequence of characters
    bstart=dtstrt(ibuff)+(frow-1)*rowlen(ibuff)+fchar-1
    call ftmbyt(iunit,bstart,.false.,status)

!       get the string of bytes
    call ftgbyt(iunit,nchars,value,status)
end
subroutine ftgtbc(tfld,tdtype,trept,tbcol,lenrow,status)
!
!*******************************************************************************
!
!! FTGTBC "Gets Table Beginning Columns."
!
!       determine the byte offset of the beginning of each field of a
!       binary table

!       tfld   i  number of fields in the binary table
!       tdtype i array of numerical datatype codes of each column
!       trept  i array of repetition factors for each column
!       OUTPUT PARAMETERS:
!       tbcol  i array giving the byte offset to the start of each column
!       lenrow i total width of the table, in bytes
!       status i  returned error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!       modified 6/17/92 to deal with ASCII column trept values measured
!       in units of characters rather than in terms of number of repeated
!       strings.

    integer tfld,tdtype(*),trept(*),tbcol(*),lenrow
    integer status,i,nbytes
    character ifld*4

    if (status > 0)return

!       the first column always begins at the first byte of the row:
    tbcol(1)=0

    do 100 i=1,tfld-1
            if (tdtype(i) == 16)then
!                       ASCII field; each character is 1 byte
                    nbytes=1
            else if (tdtype(i) > 0)then
                    nbytes=tdtype(i)/10
            else if (tdtype(i) == 0)then
!                   error: data type of column not defined! (no TFORM keyword)
                    status=232
                    write(ifld,1000)i
1000                    format(i4)
                    call ftpmsg('Field'//ifld//' of the binary'// &
                    ' table has no TFORMn keyword')
                    return
            else
!                       this is a descriptor field: 2J
                    nbytes=8
            end if

            if (nbytes == 0)then
!                       this is a bit array
                    tbcol(i+1)=tbcol(i)+(trept(i)+7)/8
            else
                    tbcol(i+1)=tbcol(i)+trept(i)*nbytes
            end if
100     continue

!       determine the total row width
    if (tdtype(tfld) == 16)then
!               ASCII field; each character is 1 byte
            nbytes=1
    else if (tdtype(tfld) > 0)then
            nbytes=tdtype(tfld)/10
    else if (tdtype(i) == 0)then
!            error: data type of column not defined! (no TFORM keyword)
            status=232
            write(ifld,1000)tfld
            call ftpmsg('Field'//ifld//' of the binary'// &
                    ' table is missing required TFORMn keyword.')
            return
    else
!               this is a descriptor field: 2J
            nbytes=8
    end if
    if (nbytes == 0)then
!               this is a bit array
            lenrow=tbcol(tfld)+(trept(tfld)+7)/8
    else
            lenrow=tbcol(tfld)+trept(tfld)*nbytes
    end if

end
subroutine ftgtbh(iunit,ncols,nrows,nfield,ttype,tbcol, &
                      tform,tunit,extnam,status)
!
!*******************************************************************************
!
!! FTGTBH is obsolete.  Call FTGHTB instead.
!

    integer iunit,ncols,nrows,nfield,status,tbcol(*)
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam

    call ftghtb(iunit,0,ncols,nrows,nfield,ttype, &
                      tbcol,tform,tunit,extnam,status)
end
subroutine ftgtbn(iunit,ncols,nrows,pcount,nfield,status)
!
!*******************************************************************************
!
!! FTGTBN checks that this is a valid binary table and get parameters.
!
!       iunit   i  Fortran i/o unit number
!       ncols   i  width of each row of the table, in bytes
!       nrows   i  number of rows in the table
!       pcount  i  size of special data area following the table (usually = 0)
!       nfield  i  number of fields in the table
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,ncols,nrows,nfield,pcount,status
    character keynam*8,value*10,comm*8,rec*80

    if (status > 0)return

!       check for correct type of extension
    call ftgrec(iunit,1,rec,status)
    if (status > 0)go to 900

    keynam=rec(1:8)

    if (keynam == 'XTENSION')then
            call ftpsvc(rec,value,comm,status)
            if (status > 0)go to 900

            if (value(2:9) /= 'BINTABLE' .and. &
                value(2:9) /= 'A3DTABLE' .and. &
                value(2:9) /= '3DTABLE ')then
!                       this is not a binary table extension
                    status=227
                    go to 900
             end if
    else
             status=225
             go to 900
    end if

!       check that the second keyword is BITPIX = 8
    call fttkyn(iunit,2,'BITPIX','8',status)
    if (status == 208)then
!               BITPIX keyword not found
            status=222
    else if (status == 209)then
!               illegal value of BITPIX
            status=211
    end if
    if (status > 0)go to 900

!       check that the third keyword is NAXIS = 2
    call fttkyn(iunit,3,'NAXIS','2',status)
    if (status == 208)then
!               NAXIS keyword not found
            status=223
    else if (status == 209)then
!               illegal NAXIS value
            status=212
    end if
    if (status > 0)go to 900

!       check that the 4th keyword is NAXIS1 and get it's value
    call ftgtkn(iunit,4,'NAXIS1',ncols,status)
    if (status == 208)then
!               NAXIS1 keyword not found
            status=224
    else if (status == 209)then
!               illegal value of NAXISnnn
            status=213
    end if
    if (status > 0)go to 900

!       check that the 5th keyword is NAXIS2 and get it's value
    call ftgtkn(iunit,5,'NAXIS2',nrows,status)
    if (status == 208)then
!               NAXIS2 keyword not found
            status=224
    else if (status == 209)then
!               illegal value of NAXISnnn
            status=213
    end if
    if (status > 0)go to 900

!       check that the 6th keyword is PCOUNT and get it's value
    call ftgtkn(iunit,6,'PCOUNT',pcount,status)
    if (status == 208)then
!               PCOUNT keyword not found
            status=228
    else if (status == 209)then
!               illegal PCOUNT value
            status=214
    end if
    if (status > 0)go to 900

!       check that the 7th keyword is GCOUNT = 1
    call fttkyn(iunit,7,'GCOUNT','1',status)
    if (status == 208)then
!               GCOUNT keyword not found
            status=229
    else if (status == 209)then
!               illegal value of GCOUNT
            status=215
    end if
    if (status > 0)go to 900

!       check that the 8th keyword is TFIELDS and get it's value
    call ftgtkn(iunit,8,'TFIELDS',nfield,status)
    if (status == 208)then
!               TFIELDS keyword not found
            status=230
    else if (status == 209)then
!               illegal value of TFIELDS
            status=216
    end if

900     continue
    if (status > 0)then
        call ftpmsg('Failed to parse the required keywords in '// &
         'the binary BINTABLE header (FTGTTB).')
    end if
end
subroutine ftgtbs(iunit,frow,fchar,nchars,svalue,status)
!
!*******************************************************************************
!
!! FTGTBS reads a string of characters from an ascii or binary table.
!
!       This will span multiple rows of the table if NCHARS+FCHAR is
!       greater than the length of a row.

!       iunit   i  fortran unit number
!       frow    i  starting row number (1st row = 1)
!       fchar   i  starting character/byte in the row to read (1st character=1)
!       nchars  i  number of characters/bytes to read (can span multiple rows)
!       svalue  c  returned string of characters
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,frow,fchar,nchars,status
    character ( len = * ) svalue

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bstart,nget

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check for errors
    if (nchars <= 0)then
!               zero or negative number of character requested
            return
    else if (frow < 1)then
!               error: illegal first row number
            status=307
            return
    else if (fchar < 1)then
!               error: illegal starting character
            status=308
            return
    end if

!       move the i/o pointer to the start of the sequence of characters
    bstart=dtstrt(ibuff)+(frow-1)*rowlen(ibuff)+fchar-1
    call ftmbyt(iunit,bstart,.false.,status)

!       get the string of characters, (up to the length of the input string)
    if (len(svalue) /= 1)then
        svalue=' '
        nget=min(nchars,len(svalue))
    else
!           assume svalue was dimensioned as: character svalue(nchars)
        nget=nchars
    end if
    call ftgcbf(iunit,nget,svalue,status)
end
subroutine ftgtcl(iunit,colnum,datcod,repeat,width,status)
!
!*******************************************************************************
!
!! FTGTCL gets the datatype, repeat count and string width of a column.
!
!       get the datatype of the column, as well as the vector
!       repeat count and (if it is an ASCII character column) the
!       width of a unit string within the column.  This supports the
!       TFORMn = 'rAw' syntax for specifying arrays of substrings.


!       iunit   i  Fortran i/o unit number
!       colnum  i  number of the column (first column = 1)

!       datcod  i  returned datatype code
!       repeat  i  number of elements in the vector column
!       width   i  width of unit string in character columns
!       status  i  returned error status
!
!       written by Wm Pence, HEASARC/GSFC, November 1994

    integer iunit,colnum,datcod,repeat,width,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,dummy
    character keywrd*8,tform*24,comm*20

    if (status > 0)return

!       construct the keyword name
    call ftkeyn('TFORM',colnum,keywrd,status)

!       get the keyword value
    call ftgkys(iunit,keywrd,tform,comm,status)
    if (status > 0)then
        call ftpmsg('Could not read the '//keywrd//' keyword.')
        return
    end if

!       parse the keyword value
    ibuff=bufnum(iunit)
    if (hdutyp(ibuff) == 1)then
!           this is an ASCII table
        repeat=1
        call ftasfm(tform,datcod,width,dummy,status)

    else if (hdutyp(ibuff) == 2)then
!           this is a binary table
        call ftbnfm(tform,datcod,repeat,width,status)

    else
!           error: this HDU is not a table
        status=235
        return
    end if
end
subroutine ftgtcs(iunit,xcol,ycol,xrval,yrval,xrpix,yrpix, &
                     xinc,yinc,rot,type,status)
!
!*******************************************************************************
!
!! FTGTCS reads the values of the celestial coordinate system keywords.
!
!       The values are read from a FITS table where the X and Y or RA and
!       DEC coordinates are stored in separate column.
!
!       These values may be used as input to the subroutines that
!       calculate celestial coordinates. (FTXYPX, FTWLDP)

!       xcol (integer) number of the column containing the RA type coordinate
!       ycol (integer) number of the column containing the DEC type coordinate

    double precision xrval,yrval,xrpix,yrpix,xinc,yinc,rot
    integer iunit,xcol,ycol,status
    character ( len = * ) type
    character comm*20,ctype*8,keynam*8,xnum*3,ynum*3

    if (status > 0)return

    call ftkeyn('TCRVL',xcol,keynam,status)
    xnum=keynam(6:8)
    call ftgkyd(iunit,keynam,xrval,comm,status)

    call ftkeyn('TCRVL',ycol,keynam,status)
    ynum=keynam(6:8)
    call ftgkyd(iunit,keynam,yrval,comm,status)

    keynam='TCRPX'//xnum
    call ftgkyd(iunit,keynam,xrpix,comm,status)
    keynam='TCRPX'//ynum
    call ftgkyd(iunit,keynam,yrpix,comm,status)

    keynam='TCDLT'//xnum
    call ftgkyd(iunit,keynam,xinc,comm,status)
    keynam='TCDLT'//ynum
    call ftgkyd(iunit,keynam,yinc,comm,status)

    keynam='TCTYP'//xnum
    call ftgkys(iunit,keynam,ctype,comm,status)

    if (status > 0)then
        call ftpmsg('FTGTCS could not find all the required'// &
                    ' celestial coordinate Keywords.')
        status=505
        return
    end if

    type=ctype(5:8)

    rot=0.
end
subroutine ftgtdm(iunit,colnum,maxdim,naxis,naxes,status)
!
!*******************************************************************************
!
!! FTGTDM parses the TDIMnnn keyword to get the dimensionality of a column.
!
!       iunit   i  fortran unit number to use for reading
!       colnum  i  column number to read
!       maxdim  i  maximum no. of dimensions to read; dimension of naxes
!       OUTPUT PARAMETERS:
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, October 1993

    integer iunit,colnum,maxdim,naxis,naxes(*),status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,nfound,c1,c2,clast,dimval
    logical last
    character*120 tdim

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

    if (colnum < 1 .or. colnum > tfield(ibuff))then
!               illegal column number
            status=302
            return
    end if

    nfound=0
!       try getting the TDIM keyword value
    call ftgkns(iunit,'TDIM',colnum,1,tdim,nfound,status)

    if (nfound /= 1)then
!               no TDIM keyword found
            naxis=1
            naxes(1)=trept(colnum+tstart(ibuff))
            return
    end if

    naxis=0
!       first, find the opening ( and closing )
    c1=index(tdim,'(')+1
    c2=index(tdim,')')-1
    if (c1 == 1 .or. c2 == -1)go to 900

    last=.false.
!       find first non-blank character
10      if (tdim(c1:c1) /= ' ')go to 20
    c1=c1+1
    go to 10

!       find the comma separating the dimension sizes
20      clast=index(tdim(c1:c2),',')+c1-2
    if (clast == c1-2)then
            last=.true.
            clast=c2
    end if

!       read the string of characters as the (integer) dimension size
    call ftc2ii(tdim(c1:clast),dimval,status)
    if (status > 0)then
         call ftpmsg('Error in FTGTDM parsing dimension string: ' &
         //tdim)
         go to 900
    end if

    naxis=naxis+1
    if (naxis <= maxdim)naxes(naxis)=dimval

    if (last)return

    c1=clast+2
    go to 10

!       could not parse the tdim value
900     status=263
end
subroutine ftgthd(tmplat,card,hdtype,status)
!
!*******************************************************************************
!
!! FTGTHD parses a template header line.
!
!
!       'Get Template HeaDer'
!       parse a template header line and create a formated
!       80-character string which is suitable for appending to a FITS header

!       tmplat  c  input header template string
!       card    c  returned 80-character string = FITS header record
!       hdtype  i  type of operation that should be applied to this keyword:
!                      -2  =  modify the name of a keyword; the new name
!                             is returned in characters 41:48 of CARD.
!                      -1  =  delete this keyword
!                       0  =  append (if it doesn't already exist) or
!                             overwrite this keyword (if it does exist)
!                       1  =  append this comment keyword ('HISTORY',
!                             'COMMENT', or blank keyword name)
!                       2  =  this is an END record; do not append it
!                             to a FITS header!
!       status  i  returned error status
!               if a positive error status is returned then the first
!               80 characters of the offending input line are returned
!               by the CARD parameter

    integer hdtype,status,tstat
    character ( len = * ) tmplat,card
    integer i1,i2,com1,strend,length
    character inline*100,keynam*8,ctemp*80,qc*1
    logical number
    double precision dvalue

    if (status > 0)return
    card=' '
    hdtype=0

    inline=tmplat

!       test if columns 1-8 are blank; if so, this is a FITS comment record;
!       just copy it verbatim to the FITS header
    if (inline(1:8) == ' ')then
            card=inline(1:80)
            go to 999
    end if

!       parse the keyword name = the first token separated by a space or a '='
!       1st locate the first nonblank character (we know it is not all blank):
    i1=0
20      i1=i1+1
!       test for a leading minus sign which flags name of keywords to be deleted
    if (inline(i1:i1) == '-')then
            hdtype=-1
!               test for a blank keyword name
            if (inline(i1+1:i1+8) == '        ')then
                   card=' '
                   i2=i1+9
                   go to 35
            end if
            go to 20
    else if (inline(i1:i1) == ' ')then
            go to 20
    end if

!       now find the last character of the keyword name
    i2=i1
30      i2=i2+1
    if (inline(i2:i2) /= ' ' .and. inline(i2:i2) /= '=')go to 30

!       test for legal keyword name length (max 8 characters)
    if (i2-i1 > 8)then
            status=207
            card=inline(1:80)
            go to 999
    end if

    keynam=inline(i1:i2-1)

!       convert to upper case and test for illegal characters in keyword name
    call ftupch(keynam)
    call fttkey(keynam,status)
    if (status > 0)then
            card=inline(1:80)
            go to 999
    end if

!       if this is the 'END' then this is the end of the input file
    if (keynam == 'END     ')goto 998

!       copy the keyword name to the output record string
    card(1:8)=keynam

!       jump if this is just the name of keyword to be deleted
    if (hdtype < 0)go to 35

!       test if this is a COMMENT or HISTORY record
    if (keynam == 'COMMENT' .or. keynam == 'HISTORY')then
!               append next 72 characters from input line to output record
            card(9:80)=inline(i2:)
            hdtype=1
            go to 999
    else
!               this keyword must have a value, so append the '= ' to output
            card(9:10)='= '
    end if

!       now locate the value token in the input line.  If it includes
!       embedded spaces it must be enclosed in single quotes. The value must
!       be separated by at least one blank space from the comment string

!       find the first character of the value string
35      i1=i2-1
40      i1=i1+1
    if (i1 > 100)then
!               no value is present in the input line
            if (hdtype < 0)then
!                       this is normal; just quit
                    go to 999
            else
                    status=204
                    card=inline(1:80)
                    go to 999
            end if
    end if
    if (hdtype < 0 .and. inline(i1:i1) == '=')then
!               The leading minus sign, plus the presence of an equal sign
!               between the first 2 tokens is taken to mean that the
!               keyword with the first token name is to be deleted.
            go to 999
    else if (inline(i1:i1)== ' ' .or.inline(i1:i1)== '=')then
            go to 40
    end if

!       is the value a quoted string?
    if (inline(i1:i1) == '''')then
!               find the closing quote
            i2=i1
50              i2=i2+1
            if (i2 > 100)then
!                       error: no closing quote on value string
                    status=205
                    card=inline(1:80)
        call ftpmsg('Keyword value string has no closing quote:')
        call ftpmsg(card)
                    go to 999
            end if
            if (inline(i2:i2) == '''')then
                    if (inline(i2+1:i2+1) == '''')then
!                               ignore 2 adjacent single quotes
                            i2=i2+1
                            go to 50
                    end if
            else
                    go to 50
            end if
!               value string can't be more than 70 characters long (cols 11-80)
            length=i2-i1
            if (length > 69)then
                    status=205
                    card=inline(1:80)
        call ftpmsg('Keyword value string is too long:')
        call ftpmsg(card)
                    go to 999
            end if

!               append value string to output, left justified in column 11
            card(11:11+length)=inline(i1:i2)
!               com1 is the starting position for the comment string
            com1=max(32,13+length)

!               FITS string must be at least 8 characters long
            if (length < 9)then
                    card(11+length:11+length)=' '
                    card(20:20)=''''
            end if
    else
!               find the end of the value field
            i2=i1
60              i2=i2+1
            if (i2 > 100)then
!                       error: value string is too long
                    status=205
                    card=inline(1:80)
        call ftpmsg('Keyword value string is too long:')
        call ftpmsg(card)
                    go to 999
            end if
            if (inline(i2:i2) /= ' ')go to 60

!               test if this is a logical value
            length=i2-i1
            if (length == 1 .and. (inline(i1:i1) == 'T' &
                .or. inline(i1:i1) == 'F'))then
                    card(30:30)=inline(i1:i1)
                    com1=32
            else
!                   test if this is a numeric value; try reading it as
!                   double precision value; if it fails, it must be a string
                number=.true.
                tstat=status
                call ftc2dd(inline(i1:i2-1),dvalue,status)
                if (status > 0)then
                    status=tstat
                    number=.false.
                else
!                       check the first character to make sure this is a number
!                       since certain non-numeric character strings pass the
!                       above test on SUN machines.
                    qc=inline(i1:i1)
                    if (qc /= '+' .and. qc /= '-' .and. qc /= &
                    '.' .and. (qc < '0' .or. qc > '9'))then
!                              This really was not a number!
                           number=.false.
                    end if
                end if

                if (number)then
                    if (length <= 20)then
!                               write the value right justified in col 30
                            card(31-length:30)=inline(i1:i2-1)
                            com1=32
                    else
!                               write the long value left justified in col 11
                            card(11:10+length)=inline(i1:i2-1)
                            com1=max(32,12+length)
                    end if
                else
!                       value is a character string datatype
                    card(11:11)=''''
                    strend=11+length
                    card(12:strend)=inline(i1:i2-1)
!                       need to expand any embedded single quotes into 2 quotes
                    i1=11
70                      i1=i1+1
                    if (i1 > strend) go to 80
                    if (card(i1:i1) == '''')then
                            i1=i1+1
                            if (card(i1:i1) /= '''')then
!                                       have to insert a 2nd quote into string
                                    ctemp=card(i1:strend)
                                    card(i1:i1)=''''
                                    strend=strend+1
                                    i1=i1+1
                                    card(i1:strend)=ctemp
                            end if
                    end if
                    go to 70

80                      strend=max(20,strend+1)
                    card(strend:strend)=''''
                    com1=max(32,strend+2)
                end if
            end if
    end if

!       check if this was a request to modify a keyword name
    if (hdtype == -1)then
            hdtype = -2
!               the keyword value is really the new keyword name
!               return the new name in characters 41:48 of the output card
            keynam=card(12:19)
!               convert to upper case and test for illegal characters in name
            call ftupch(keynam)
            call fttkey(keynam,status)
            if (status > 0)then
                    card=inline(1:80)
                    go to 999
            else
                    card(9:80)=' '
                    card(41:48)=keynam
                    go to 999
            end if
    end if

!       is there room for a comment string?
    if (com1 < 79)then
!               now look for the beginning of the comment string
            i1=i2
90              i1=i1+1
!               if no comment field then just quit
            if (i1 > 100)go to 999
            if (inline(i1:i1) == ' ')go to 90

!               append the comment field
            if (inline(i1:i1) == '/')then
                    card(com1:80)=inline(i1:)
            else
                    card(com1:80)='/ '//inline(i1:)
            end if
    end if

    go to 999

!   end of input file was detected
998     hdtype=2

999     continue
end
subroutine ftgtkn(iunit,nkey,keynam,ival,status)
!
!*******************************************************************************
!
!! FTGTKN tests that a key has the right name, and gets its value.
!
!
!  The routine tests that keyword number NKEY has name = KEYNAM and get the
!       integer value of the keyword.  Return an error if the keyword
!       name does not match the input KEYNAM, or if the value of the
!       keyword is not a positive integer.
!
!       iunit   i  Fortran I/O unit number
!       nkey    i  sequence number of the keyword to test
!       keynam  c  name that the keyword is supposed to have
!       OUTPUT PARAMETERS:
!       ival    i  returned value of the integer keyword
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!
    integer iunit,nkey,status,ival
    character ( len = * ) keynam
    character kname*8,value*30,comm*48,npos*8,keybuf*80

    if (status > 0)return

!       read the name and value of the keyword
    call ftgrec(iunit,nkey,keybuf,status)

    kname=keybuf(1:8)
!       parse the value and comment fields from the record
    call ftpsvc(keybuf,value,comm,status)

    if (status > 0)go to 900

!       test if the keyword has the correct name
    if (kname /= keynam)then
            status=208
            go to 900
    end if

!       convert character string to integer
    call ftc2ii(value,ival,status)
    if (status > 0 .or. ival < 0 )then
!               keyword value must be zero or positive integer
            status=209
    end if

900     continue

    if (status > 0)then
        write(npos,1000)nkey
1000        format(i8)
        call ftpmsg('FTGTKN found unexpected keyword or value '// &
        'for header keyword number '//npos//'.')
        call ftpmsg('  Was expecting positive integer keyword '// &
        keynam(1:8))
        if (keybuf(9:10) /= '= ')then
            call ftpmsg('  but found the keyword '//kname// &
            ' with no value field (no "= " in cols. 9-10).')
        else
          call ftpmsg('  but instead found keyword = '//kname// &
          ' with value = '//value)
        end if
        call ftpmsg(keybuf)
    end if
end
subroutine ftgttb(iunit,ncols,nrows,nfield,status)
!
!*******************************************************************************
!
!! FTGTTB tests that this is a legal ASCII table, and gets some keywords.
!
!       iunit   i  Fortran i/o unit number
!       OUTPUT PARAMETERS:
!       ncols   i  number of columns in the table
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,ncols,nrows,nfield,status
    character keynam*8,value*10,comm*8,keybuf*80

    if (status > 0)return

!       check for correct type of extension
    call ftgrec(iunit,1,keybuf,status)

    keynam=keybuf(1:8)
!       parse the value and comment fields from the record
    call ftpsvc(keybuf,value,comm,status)

    if (status > 0)go to 900

    if (keynam == 'XTENSION')then
            if (value(2:9) /= 'TABLE   ')then
!                       this is not a ASCII table extension
                    status=226
    call ftpmsg('Was expecting an ASCII table; instead got '// &
    'XTENSION= '//value)
                    call ftpmsg(keybuf)
                    go to 900
             end if
    else
             status=225
    call ftpmsg('First keyword of extension was not XTENSION:'// &
             keynam)
             call ftpmsg(keybuf)
             go to 900
    end if

!       check that the second keyword is BITPIX = 8
    call fttkyn(iunit,2,'BITPIX','8',status)
    if (status == 208)then
!               BITPIX keyword not found
            status=222
    else if (status == 209)then
!               illegal value of BITPIX
            status=211
    end if
    if (status > 0)go to 900

!       check that the third keyword is NAXIS = 2
    call fttkyn(iunit,3,'NAXIS','2',status)
    if (status == 208)then
!               NAXIS keyword not found
            status=223
    else if (status == 209)then
!               illegal value of NAXIS
            status=212
    end if
    if (status > 0)go to 900

!       check that the 4th keyword is NAXIS1 and get it's value
    call ftgtkn(iunit,4,'NAXIS1',ncols,status)
    if (status == 208)then
!               NAXIS1 keyword not found
            status=224
    else if (status == 209)then
!               illegal NAXIS1 value
            status=213
    end if
    if (status > 0)go to 900

!       check that the 5th keyword is NAXIS2 and get it's value
    call ftgtkn(iunit,5,'NAXIS2',nrows,status)
    if (status == 208)then
!               NAXIS2 keyword not found
            status=224
    else if (status == 209)then
!               illegal NAXIS2 value
            status=213
    end if
    if (status > 0)go to 900

!       check that the 6th keyword is PCOUNT = 0
    call fttkyn(iunit,6,'PCOUNT','0',status)
    if (status == 208)then
!               PCOUNT keyword not found
            status=228
    else if (status == 209)then
!               illegal PCOUNT value
            status=214
    end if
    if (status > 0)go to 900

!       check that the 7th keyword is GCOUNT = 1
    call fttkyn(iunit,7,'GCOUNT','1',status)
    if (status == 208)then
!               GCOUNT keyword not found
            status=229
    else if (status == 209)then
!               illegal value of GCOUNT
            status=215
    end if
    if (status > 0)go to 900

!       check that the 8th keyword is TFIELDS
    call ftgtkn(iunit,8,'TFIELDS',nfield,status)
    if (status == 208)then
!               TFIELDS keyword not found
            status=230
    else if (status == 209)then
!               illegal value of TFIELDS
            status=216
    end if

900     continue
    if (status > 0)then
        call ftpmsg('Failed to parse the required keywords in '// &
         'the ASCII TABLE header (FTGTTB).')
    end if
end
subroutine ftgunt(iunit,keywrd,kunit,status)
!
!*******************************************************************************
!
!! FTGUNT reads the unit string from the comment string from a header record.
!
!       iunit   i  fortran input unit number
!       keywrd  c  keyword name
!       OUTPUT PARAMETERS:
!       kunit   c  output keyword units
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1997

    character ( len = * ) keywrd,kunit
    integer iunit,ii,status,ulen
    character value*35,comm*72

    if (status > 0)return

    kunit = ' '

!       find the keyword and return value and comment as character strings
    call ftgkey(iunit,keywrd,value,comm,status)

    if (status > 0)return

!       look for brackets enclosing the units string
    if (comm(1:1) == '[')then
       ulen=2
       do ii = 3,72
          if (comm(ii:ii) == ']')then
            kunit=comm(2:ulen)
            return
          end if
          ulen=ii
       end do
       return

    end if
end
subroutine fthdef(ounit,moreky,status)
!
!*******************************************************************************
!
!! FTHDEF defines the size of the current header unit.
!
!       define the size of the current header unit; this simply lets
!       us determine where the data unit will start
!
!       ounit   i  Fortran I/O unit number
!       moreky  i  number of additional keywords to reserve space for
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,moreky,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,mkeys

    if (status > 0)return

!       based on the number of keywords which have already been written,
!       plus the number of keywords to reserve space for, we then can
!       define where the data unit should start (it must start at the
!       beginning of a 2880-byte logical block).

    ibuff=bufnum(ounit)

    mkeys=max(moreky,0)
    dtstrt(ibuff)=((hdend(ibuff)+mkeys*80)/2880+1)*2880
end
subroutine fthpdn(ounit,nbytes,status)
!
!*******************************************************************************
!
!! FTHPDN shifts the binary table heap down by nbyte bytes.
!
!       ounit   i  fortran output unit number
!       nbytes  i  number of bytes by which to move the heap
!       status  i  returned error status (0=ok)

    integer ounit,nbytes,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*5760 buff
    character xdummy(26240)
    common/ftheap/buff,xdummy
!

    integer i,ibuff,ntodo,jpoint,nchar,tstat

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

    if (heapsz(ibuff) > 0)then
        ntodo=heapsz(ibuff)

!           set pointer to the end of the heap
        jpoint=dtstrt(ibuff)+theap(ibuff)+heapsz(ibuff)

10          nchar=min(ntodo,5760)
        jpoint=jpoint-nchar

!           move to the read start position
        call ftmbyt(ounit,jpoint,.false.,status)

!           read the heap
        call ftgcbf(ounit,nchar,buff,status)

!           move forward to the write start postion
        call ftmbyt(ounit,jpoint+nbytes,.true.,status)

!           write the heap
        call ftpcbf(ounit,nchar,buff,status)

!           check for error
        if (status > 0)then
           call ftpmsg('Error while moving heap down (FTDNHP)')
           return
        end if

!           check for more data in the heap
        ntodo=ntodo-nchar
        if (ntodo > 0)go to 10

!           now overwrite the old fill data with zeros
        do 20 i=1,5760
            buff(i:i)=char(0)
20          continue

        jpoint=dtstrt(ibuff)+theap(ibuff)
        call ftmbyt(ounit,jpoint,.false.,status)

        ntodo=nbytes
30          nchar=min(ntodo,5760)
        call ftpcbf(ounit,nchar,buff,status)
        ntodo=ntodo-nchar
        if (ntodo > 0)go to 30
    end if

!       update the heap starting address
    theap(ibuff)=theap(ibuff)+nbytes

!       try updating the keyword value, if it exists
    tstat=status
    call ftmkyj(ounit,'THEAP',theap(ibuff),'&',status)
    if (status == 202)status=tstat
end
subroutine fthpup(ounit,nbytes,status)
!
!*******************************************************************************
!
!! FTHPUP shifts the binary table heap up by nbytes bytes.
!
!       ounit   i  fortran output unit number
!       nbytes  i  number of bytes by which to move the heap
!       status  i  returned error status (0=ok)

    integer ounit,nbytes,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*5760 buff
    character xdummy(26240)
    common/ftheap/buff,xdummy
!

    integer i,ibuff,ntodo,jpoint,nchar,tstat

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

    if (heapsz(ibuff) > 0)then
        ntodo=heapsz(ibuff)

!           set pointer to the start of the heap
        jpoint=dtstrt(ibuff)+theap(ibuff)

10          nchar=min(ntodo,5760)

!           move to the read start position
        call ftmbyt(ounit,jpoint,.false.,status)

!           read the heap
        call ftgcbf(ounit,nchar,buff,status)

!           move back to the write start postion
        call ftmbyt(ounit,jpoint-nbytes,.false.,status)

!           write the heap
        call ftpcbf(ounit,nchar,buff,status)

!           check for error
        if (status > 0)then
           call ftpmsg('Error while moving heap up (FTUPHP)')
           return
        end if

!           check for more data in the heap
        ntodo=ntodo-nchar
        jpoint=jpoint+nchar
        if (ntodo > 0)go to 10

!           now overwrite the old fill data with zeros
        do 20 i=1,5760
            buff(i:i)=char(0)
20          continue

        jpoint=dtstrt(ibuff)+theap(ibuff)+heapsz(ibuff)-nbytes
        call ftmbyt(ounit,jpoint,.false.,status)

        ntodo=nbytes
30          nchar=min(ntodo,5760)
        call ftpcbf(ounit,nchar,buff,status)
        ntodo=ntodo-nchar
        if (ntodo > 0)go to 30
    end if

!       update the heap starting address
    theap(ibuff)=theap(ibuff)-nbytes

!       try updating the keyword value, if it exists
    tstat=status
    call ftmkyj(ounit,'THEAP',theap(ibuff),'&',status)
    if (status == 202)status=tstat
end
subroutine fti1i1(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI1I1 copies input i*1 values to output i*1 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   c*1 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  c*1 value in the input array that is used to indicated nulls
!       setval  c*1 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  c*1 returned array of values
!       status  i  output error status (0 = ok)

    character input(*),chkval
    character output(*),setval
    integer n,i,chktyp,status,itemp
    double precision scale,zero,dval
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                      itemp=ichar(input(i))
                      if (itemp < 0)itemp=itemp+256
                      dval=(itemp-zero)/scale
!                         trap any values that overflow the I*1 range
                      if (dval< 255.49 .and. dval> -.49)then
                            output(i)=char(nint(dval))
                      else if (dval >= 255.49)then
                            status=-11
                            output(i)=char(255)
                      else
                            status=-11
                            output(i)=char(0)
                      end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                   don't have to check for nulls
                if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                else
                    do 40 i=1,n
                        itemp=ichar(input(i))
                        if (itemp < 0)itemp=itemp+256
                        dval=itemp*scale+zero
!                           trap any values that overflow the I*1 range
                      if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
40                      continue
                end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                        output(i)=input(i)
                                    end if
50                              continue
                    else
                     do 60 i=1,n
                      if (input(i) == chkval)then
                                anynul=.true.
                                if (chktyp == 1)then
                                    output(i)=setval
                                else
                                    flgray(i)=.true.
                                end if
                      else
                        itemp=ichar(input(i))
                        if (itemp < 0)itemp=itemp+256
                        dval=itemp*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
                      end if
60                       continue
                    end if
            end if
    end if
end
subroutine fti1i2(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI1I2 copies input i*1 values to output i*2 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   c*1 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  c*1 value in the input array that is used to indicated nulls
!       setval  i*2 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i*2 returned array of values
!       status  i  output error status (0 = ok)

    character input(*),chkval
    integer*2 output(*),setval,mini2,maxi2
    integer n,i,chktyp,status,itemp
    double precision scale,zero,dval,i2max,i2min
    logical tofits,flgray(*),anynul,noscal

    parameter (maxi2=32767)
    parameter (mini2=-32768)
    parameter (i2max=3.276749D+04)
    parameter (i2min=-3.276849D+04)

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            itemp=ichar(input(i))
                            if (itemp < 0)itemp=itemp+256
                            output(i)=itemp
10                      continue
            else
                    do 20 i=1,n
                        itemp=ichar(input(i))
                        if (itemp < 0)itemp=itemp+256
                        dval=(itemp-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i2max .and. dval>i2min)then
                            output(i)=nint(dval)
                        else if (dval >= i2max)then
                            status=-11
                            output(i)=maxi2
                        else
                            status=-11
                            output(i)=mini2
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                  itemp=ichar(input(i))
                                  if (itemp < 0)itemp=itemp+256
                                  output(i)=itemp
30                              continue
                    else
                        do 40 i=1,n
                          itemp=ichar(input(i))
                          if (itemp < 0)itemp=itemp+256
                          dval=itemp*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=dval
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                  itemp=ichar(input(i))
                                  if (itemp < 0)itemp=itemp+256
                                  output(i)=itemp
                                    end if
50                              continue
                else
                    do 60 i=1,n
                        if (input(i) == chkval)then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          itemp=ichar(input(i))
                          if (itemp < 0)itemp=itemp+256
                          dval=itemp*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=dval
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine fti1i4(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI1I4 copies input i*1 values to output i*4 values.
!
!  The routine does optional scaling and checking for null values.
!

!       input   c*1 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  c*1 value in the input array that is used to indicated nulls
!       setval  i   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i   returned array of values
!       status  i  output error status (0 = ok)

    character input(*),chkval
    integer output(*),setval
    integer n,i,chktyp,status,itemp
    double precision scale,zero,dval,i4max,i4min
    logical tofits,flgray(*),anynul,noscal
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            itemp=ichar(input(i))
                            if (itemp < 0)itemp=itemp+256
                            output(i)=itemp
10                      continue
            else
                    do 20 i=1,n
                        itemp=ichar(input(i))
                        if (itemp < 0)itemp=itemp+256
                        dval=(itemp-zero)/scale
!                           trap any values that overflow the I*4 range
                        if (dval<i4max .and. dval>i4min)then
                            output(i)=nint(dval)
                        else if (dval >= i4max)then
                            status=-11
                            output(i)=maxi4
                        else
                            status=-11
                            output(i)=mini4
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                  itemp=ichar(input(i))
                                  if (itemp < 0)itemp=itemp+256
                                  output(i)=itemp
30                              continue
                    else
                        do 40 i=1,n
                          itemp=ichar(input(i))
                          if (itemp < 0)itemp=itemp+256
                          dval=itemp*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=dval
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                  itemp=ichar(input(i))
                                  if (itemp < 0)itemp=itemp+256
                                  output(i)=itemp
                                    end if
50                              continue
                else
                    do 60 i=1,n
                        if (input(i) == chkval)then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          itemp=ichar(input(i))
                          if (itemp < 0)itemp=itemp+256
                          dval=itemp*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=dval
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine fti1r4(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI1R4 copies input i*1 values to output r*4 values.
!
!  The routine also does optional scaling and checking for null values
!
!       input   c*1 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  c*1 value in the input array that is used to indicated nulls
!       setval  r   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  r   returned array of values

    character input(*),chkval
    real output(*),setval
    integer n,i,chktyp,status,itemp
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            itemp=ichar(input(i))
                            if (itemp < 0)itemp=itemp+256
                            output(i)=itemp
10                      continue
            else
                    do 20 i=1,n
                            itemp=ichar(input(i))
                            if (itemp < 0)itemp=itemp+256
                            output(i)=(itemp-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    itemp=ichar(input(i))
                                    if (itemp < 0)itemp=itemp+256
                                    output(i)=itemp
30                              continue
                    else
                            do 40 i=1,n
                              itemp=ichar(input(i))
                              if (itemp < 0)itemp=itemp+256
                              output(i)=itemp*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                        itemp=ichar(input(i))
                                    if (itemp < 0)itemp=itemp+256
                                        output(i)=itemp
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                              itemp=ichar(input(i))
                              if (itemp < 0)itemp=itemp+256
                              output(i)=itemp*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine fti1r8(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI1R8 copies input i*1 values to output r*8 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   c*1 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  c*1 value in the input array that is used to indicated nulls
!       setval  d   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  d   returned array of values

    character input(*),chkval
    double precision output(*),setval
    integer n,i,chktyp,status,itemp
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            itemp=ichar(input(i))
                            if (itemp < 0)itemp=itemp+256
                            output(i)=itemp
10                      continue
            else
                    do 20 i=1,n
                            itemp=ichar(input(i))
                            if (itemp < 0)itemp=itemp+256
                            output(i)=(itemp-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    itemp=ichar(input(i))
                                    if (itemp < 0)itemp=itemp+256
                                    output(i)=itemp
30                              continue
                    else
                            do 40 i=1,n
                              itemp=ichar(input(i))
                              if (itemp < 0)itemp=itemp+256
                              output(i)=itemp*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                        itemp=ichar(input(i))
                                    if (itemp < 0)itemp=itemp+256
                                        output(i)=itemp
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                              itemp=ichar(input(i))
                              if (itemp < 0)itemp=itemp+256
                              output(i)=itemp*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine fti2c(ival,cval,status)
!
!*******************************************************************************
!
!! FTI2C converts an integer value to a C*20 character string, right justified.
!
    integer ival,status
    character*20 cval

    if (status > 0)return

    write(cval,1000,err=900)ival
1000    format(i20)
    if (cval(1:1) == '*')go to 900
    return
900     status=401
    call ftpmsg('Error in FTI2C converting integer to C*20 string.')
end
subroutine fti2i1(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI2I4 copies input i*2 values to output i*1 values.
!
!  The routine does optional scaling and checking for null values.
!
!
!       input   i*2 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i*2 value in the input array that is used to indicated nulls
!       setval  c*1 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  c*1 returned array of values
!       status  i  output error status (0 = ok)

    integer*2 input(*),chkval
    character output(*),setval
    integer n,i,chktyp,itemp,status
    double precision scale,zero,dval
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
!       have to use a temporary variable because of IBM mainframe
                        itemp=input(i)
!                           trap any values that overflow the I*1 range
                        if (itemp<= 255 .and. itemp>= 0)then
                            output(i)=char(itemp)
                        else if (itemp > 255)then
                            status=-11
                            output(i)=char(255)
                        else
                            status=-11
                            output(i)=char(0)
                        end if
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                            output(i)=char(nint(dval))
                        else if (dval >= 255.49)then
                            status=-11
                            output(i)=char(255)
                        else
                            status=-11
                            output(i)=char(0)
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                        do 30 i=1,n
!       have to use a temporary variable because of IBM mainframe
                            itemp=input(i)
!                               trap any values that overflow the I*1 range
                            if (itemp<= 255 .and. itemp>= 0)then
                                output(i)=char(itemp)
                            else if (itemp > 255)then
                                status=-11
                                output(i)=char(255)
                            else
                                status=-11
                                output(i)=char(0)
                            end if
30                          continue
                    else
                      do 40 i=1,n
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
40                        continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                      do 50 i=1,n
                         if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                         else
!       have to use a temporary variable because of IBM mainframe
                            itemp=input(i)
!                               trap any values that overflow the I*1 range
                            if (itemp<= 255 .and. itemp>= 0)then
                                output(i)=char(itemp)
                            else if (itemp > 255)then
                                status=-11
                                output(i)=char(255)
                            else
                                status=-11
                                output(i)=char(0)
                            end if
                         end if
50                        continue
                  else
                      do 60 i=1,n
                        if (input(i) == chkval)then
                                anynul=.true.
                                if (chktyp == 1)then
                                    output(i)=setval
                                else
                                    flgray(i)=.true.
                                end if
                      else
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
                      end if
60                       continue
                end if
            end if
    end if
end
subroutine fti2i2(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI2I2 copies input i*2 values to output i*2 values.
!
!  The routine does optional scaling and checking for null values.
!

!       input   i*2 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i*2 value in the input array that is used to indicated nulls
!       setval  i*2 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i*2 returned array of values
!       status  i  output error status (0 = ok)

!        integer*2 j (this was only needed to workaround the Microsoft bug)

    integer*2 input(*),output(*),chkval,setval,mini2,maxi2
    integer n,i,chktyp,status
    double precision scale,zero,dval,i2max,i2min
    logical tofits,flgray(*),anynul,noscal

    parameter (maxi2=32767)
    parameter (mini2=-32768)
    parameter (i2max=3.276749D+04)
    parameter (i2min=-3.276849D+04)

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits)then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n

!      The following workaround was removed Dec 1996.  Hopefully this
!      compiler bug is fixed in later versions, but in any case, it is more
!      important to remove this workaround to make the code more efficient
!      on other machines
!                       Have to use internal variable j to work around
!                       a bug in the Microsoft v5.0 compiler on IBM PCs
!                               j=input(i)
!                               output(i)=j

                           output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i2max .and. dval>i2min)then
                            output(i)=nint(dval)
                        else if (dval >= i2max)then
                            status=-11
                            output(i)=maxi2
                        else
                            status=-11
                            output(i)=mini2
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
!                               Have to use internal variable j to work around
!                               a bug in the Microsoft v5.0 compiler on IBM PCs
!                                        j=input(i)
!                                        output(i)=j

                                    output(i)=input(i)
30                              continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=dval
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
!                               Have to use internal variable j to work around
!                               a bug in the Microsoft v5.0 compiler on IBM PCs
!                                                j=input(i)
!                                                output(i)=j

                                            output(i)=input(i)
                                    end if
50                              continue
                else
                    do 60 i=1,n
                        if (input(i) == chkval)then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=dval
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine fti2i4(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI2I4 copies input i*2 values to output i*4 values.
!
!  The routine does optional scaling and checking for null values.
!

!       input   i*2 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i*2 value in the input array that is used to indicated nulls
!       setval  i   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i   returned array of values
!       status  i  output error status (0 = ok)

    integer*2 input(*),chkval
    integer output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval,i4max,i4min
    logical tofits,flgray(*),anynul,noscal
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i4max .and. dval>i4min)then
                            output(i)=nint(dval)
                        else if (dval >= i4max)then
                            status=-11
                            output(i)=maxi4
                        else
                            status=-11
                            output(i)=mini4
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=dval
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                else
                    do 60 i=1,n
                        if (input(i) == chkval)then
                             anynul=.true.
                             if (chktyp == 1)then
                                  output(i)=setval
                              else
                                  flgray(i)=.true.
                              end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=dval
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
                        end if
60                      continue
                 end if
            end if
    end if
end
subroutine fti2r4(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI2R4 copies input i*2 values to output r*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   i*2 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i*2 value in the input array that is used to indicated nulls
!       setval  r   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  r   returned array of values

    integer*2 input(*),chkval
    real output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine fti2r8(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI2R8 copies input i*2 values to output r*8 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   i*2 input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i*2 value in the input array that is used to indicated nulls
!       setval  d   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  d   returned array of values

    integer*2 input(*),chkval
    double precision output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine fti4i1(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI4I1 copies input i*4 values to output i*1 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   i input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i value in the input array that is used to indicated nulls
!       setval  c*1 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  c*1 returned array of values
!       status  i  output error status (0 = ok)

    integer input(*),chkval
    character output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                do 10 i=1,n
!                       trap any values that overflow the I*1 range
                    if (input(i)<= 255 .and. input(i)>= 0)then
                            output(i)=char(input(i))
                    else if (input(i) > 255)then
                            status=-11
                            output(i)=char(255)
                    else
                            status=-11
                            output(i)=char(0)
                    end if
10                  continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                            output(i)=char(nint(dval))
                        else if (dval >= 255.49)then
                            status=-11
                            output(i)=char(255)
                        else
                            status=-11
                            output(i)=char(0)
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                   don't have to check for nulls
                if (noscal)then
                  do 30 i=1,n
!                       trap any values that overflow the I*1 range
                    if (input(i)<= 255 .and. input(i)>= 0)then
                            output(i)=char(input(i))
                    else if (input(i) > 255)then
                            status=-11
                            output(i)=char(255)
                    else
                            status=-11
                            output(i)=char(0)
                    end if
30                    continue
                else
                    do 40 i=1,n
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
40                      continue
                end if
            else
!                   must test for null values
                if (noscal)then
                     do 50 i=1,n
                         if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                          else
!                               trap any values that overflow the I*1 range
                            if (input(i)<= 255 .and. &
                                input(i)>= 0)then
                                output(i)=char(input(i))
                            else if (input(i) > 255)then
                                status=-11
                                output(i)=char(255)
                            else
                                status=-11
                                output(i)=char(0)
                            end if
                         end if
50                       continue
                else
                  do 60 i=1,n
                    if (input(i) == chkval)then
                                anynul=.true.
                                if (chktyp == 1)then
                                    output(i)=setval
                                else
                                    flgray(i)=.true.
                                end if
                     else
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
                     end if
60                     continue
                end if
            end if
    end if
end
subroutine fti4i2(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI4I2 copies input i*4 values to output i*2 values.
!
!  The routine does optional scaling and checking for null values.
!

!       input   i  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i  value in the input array that is used to indicated nulls
!       setval  i*2 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i*2 returned array of values
!       status  i  output error status (0 = ok)

    integer input(*),chkval
    integer*2 output(*),setval
    integer n,i,chktyp,status,maxi2,mini2
    double precision scale,zero,dval,i2max,i2min
    logical tofits,flgray(*),anynul,noscal
    parameter (i2max=3.276749D+04)
    parameter (i2min=-3.276849D+04)
    parameter (maxi2=32767)
    parameter (mini2=-32768)

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
!                           trap any values that overflow the I*2 range
                        if (input(i) <= maxi2 .and. &
                            input(i) >= mini2)then
                                output(i)=input(i)
                        else if (input(i) > maxi2)then
                                status=-11
                                output(i)=maxi2
                        else
                                status=-11
                                output(i)=mini2
                        end if
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i2max .and. dval>i2min)then
                            output(i)=nint(dval)
                        else if (dval >= i2max)then
                            status=-11
                            output(i)=maxi2
                        else
                            status=-11
                            output(i)=mini2
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                      do 30 i=1,n
!                           trap any values that overflow the I*2 range
                        if (input(i) <= maxi2 .and. &
                            input(i) >= mini2)then
                                output(i)=input(i)
                        else if (input(i) > maxi2)then
                                status=-11
                                output(i)=maxi2
                        else
                                status=-11
                                output(i)=mini2
                        end if
30                        continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=dval
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                       do 50 i=1,n
                          if (input(i) == chkval)then
                                    anynul=.true.
                                    if (chktyp == 1)then
                                            output(i)=setval
                                    else
                                            flgray(i)=.true.
                                    end if
                          else
!                               trap any values that overflow the I*2 range
                            if (input(i) <= maxi2 .and. &
                                input(i) >= mini2)then
                                    output(i)=input(i)
                            else if (input(i) > maxi2)then
                                    status=-11
                                    output(i)=maxi2
                            else
                                    status=-11
                                    output(i)=mini2
                            end if
                          end if
50                         continue
                else
                    do 60 i=1,n
                        if (input(i) == chkval)then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=dval
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine fti4i4(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI4I4 copies input i*4 values to output i*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   i  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i   value in the input array that is used to indicated nulls
!       setval  i   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i   returned array of values
!       status  i  output error status (0 = ok)

    integer input(*),chkval
    integer output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval,i4max,i4min
    logical tofits,flgray(*),anynul,noscal
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i4max .and. dval>i4min)then
                            output(i)=nint(dval)
                        else if (dval >= i4max)then
                            status=-11
                            output(i)=maxi4
                        else
                            status=-11
                            output(i)=mini4
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=dval
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                else
                    do 60 i=1,n
                        if (input(i) == chkval)then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=dval
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine fti4r4(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI4R4 copies input i*4 values to output r*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   i  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i  value in the input array that is used to indicated nulls
!       setval  r  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       output  r  returned array of values

    integer input(*),chkval
    real output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine fti4r8(input,n,scale,zero,tofits, &
            chktyp,chkval,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTI4R8 copies input i*4 values to output r*8 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   i  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       chkval  i  value in the input array that is used to indicated nulls
!       setval  d  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       output  d  returned array of values

    integer input(*),chkval
    double precision output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (input(i) == chkval)then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine ftibin(ounit,nrows,nfield,ttype,tform,tunit, &
                      extnam,pcount,status)
!
!*******************************************************************************
!
!! FTIBIN inserts a binary table extension following the current HDU.
!
!       ounit   i  fortran output unit number
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       ttype   c  name of each field (array) (optional)
!       tform   c  format of each field (array)
!       tunit   c  units of each field (array) (optional)
!       extnam  c  name of table extension (optional)
!       pcount  i  size of special data area following the table (usually = 0)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,nrows,nfield,pcount,status
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,nhdu,i,savstr,nblock,hsize,nkey

    if (status > 0)return
    ibuff=bufnum(ounit)

!       close the current HDU to make sure END and fill values are written
    call ftchdu(ounit,status)
    if (status > 0)return

!       save the starting address of the next HDU
    nhdu=chdu(ibuff)+1
    savstr=hdstrt(ibuff,nhdu)

!       count number of optional TUNITS keywords to be written
    nkey=0
    do 5 i=1,nfield
           if (tunit(i) /= ' ')nkey=nkey+1
5       continue
    if (extnam /= ' ')nkey=nkey+1

!       calc min size of header
    nblock=(9 + 2*nfield + nkey +35)/36
    hsize=nblock*2880

!       define a fake CHDU with a minimum header
    dtstrt(ibuff)=hdstrt(ibuff,chdu(ibuff))+hsize

!       define the size of the new HDU (this modifies hdstrt(ibuff,nhdu))
    call ftbdef(ounit,nfield,tform,pcount,nrows,status)

!       use start of next HDU to calc. how big this new HDU is
    nblock=(hdstrt(ibuff,nhdu)-hdstrt(ibuff,nhdu-1))/2880

!       reset the start of the next HDU back to it original value
    hdstrt(ibuff,nhdu)=savstr

!       insert the required number of blocks at the end of the real CHDU
!       (first define hdutyp so that the correct fill value will be used)
    hdutyp(ibuff)=2
    call ftiblk(ounit,nblock,1,status)
    if (status > 0)return

!       increment the number of HDUs in the file and their starting address
    maxhdu(ibuff)=maxhdu(ibuff)+1
    do 10 i=maxhdu(ibuff),nhdu,-1
            hdstrt(ibuff,i+1)=hdstrt(ibuff,i)
10      continue

!       again, reset the start of the next HDU back to it original value
    hdstrt(ibuff,nhdu)=savstr

!       flush the buffers holding data for the old HDU
    call ftflsh(ibuff,status)

!       recover common block space containing column descriptors for old HDU
    call ftfrcl(ounit,status)

!       move to the new (empty) HDU
    chdu(ibuff)=nhdu

!       set parameters describing an empty header
    hdutyp(ibuff)=2
    nxthdr(ibuff)=hdstrt(ibuff,nhdu)
    hdend(ibuff)= hdstrt(ibuff,nhdu)
    dtstrt(ibuff)=hdstrt(ibuff,nhdu)+hsize

!       write the header keywords
    call ftphbn(ounit,nrows,nfield,ttype,tform,tunit,extnam, &
                pcount,status)

!       define the structure of the new HDU
    call ftbdef(ounit,nfield,tform,pcount,nrows,status)
end
subroutine ftiblk(ounit,nblock,hdrdat,status)
!
!*******************************************************************************
!
!! FTIBLK inserts a 2880-byte block at the end of the current header or data.
!
!       ounit   i  fortran output unit number
!       nblock  i  number of blocks to insert
!       hdrdat  i  insert space in header (0) or data (1)
!       status  i  returned error status (0=ok)

    integer ounit,nblock,hdrdat,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    character*2880 buff(2)
    character xdummy(26240)
    common/ftheap/buff,xdummy
!

    integer ibuff,ipoint,jpoint,i,tstat,thdu,nshift,in,out,tin
    character cfill

    if (status > 0)return
    tstat=status

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

!       set the appropriate fill value
    if (hdrdat == 0 .or. hdutyp(ibuff) == 1)then
!               fill  header or ASCII table with space
            cfill=char(32)
    else
!               fill with Null (0) in image or bintable data area
            cfill=char(0)
    end if

!       find position in file to insert new block
    if (hdrdat == 0)then
        ipoint=dtstrt(ibuff)
    else
        ipoint=hdstrt(ibuff,chdu(ibuff)+1)
    end if


    if (nblock == 1 .and. hdrdat == 0)then
!
!  Don't use this algoritm, even though it may be faster (but initial
!  tests showed it didn't make any difference on a SUN) because it is
!  less safe than the other more general algorithm.  If there is
!  not enough disk space available for the added block, this faster
!  algorithm won't fail until it tries to move the last block, thus leaving
!  the FITS file in a corrupted state.   The other more general
!  algorithm tries to add a new empty block to the file as the
!  first step.  If this fails, it still leaves the current FITS
!  file unmodified, which is better for the user.
!
!  (Note added later:)
!  Will use this algorithm anyway when inserting one block in a FITS
!  header because the more general algorithm results in a status=252 error
!  in cases where the number of rows in a table has not yet been defined
!
!           use this more efficient algorithm if just adding a single block
!           initialize the first buffer
        do 5 i=1,2880
           buff(1)(i:i)=cfill
5           continue

        in=2
        out=1

!           move to the read start position
10          call ftmbyt(ounit,ipoint,.false.,status)

!           read one 2880-byte FITS logical record into the input buffer
        call ftgcbf(ounit,2880,buff(in),status)

!           check for End-Of-File
        if (status == 107)go to 20

!           move back to the write start postion
        call ftmbyt(ounit,ipoint,.false.,status)

!           write the 2880-byte FITS logical record stored in the output buffer
        call ftpcbf(ounit,2880,buff(out),status)

!           check for error during write (the file may not have write access)
        if (status > 0)return

!           swap the input and output buffer pointers and move to next block
        tin=in
        in=out
        out=tin
        ipoint=ipoint+2880

!           now repeat the process until we reach the End-Of-File
        go to 10

!           we have reached the end of file; now append the last block
20          status=tstat

!           move back to the write start postion
        call ftmbyt(ounit,ipoint,.true.,status)

!           write the 2880-byte FITS logical record stored in the output buffer
        call ftpcbf(ounit,2880,buff(out),status)

    else
!           use this general algorithm for adding arbitrary number of blocks

!           first, find the end of file
        thdu=chdu(ibuff)

30          call ftmahd(ounit,maxhdu(ibuff)+1,i,status)

        if (status == 107)then
            status=tstat
!               move back to the current extension
            call ftmahd(ounit,thdu,i,status)
            go to 100
        else if (status <= 0)then
            go to 30
        else
            call ftpmsg('Error while seeking End of File (FTIBLK)')
            return
        end if

!           calculate number of 2880-byte blocks that have to be shifted down
100         continue
        nshift=(hdstrt(ibuff,maxhdu(ibuff)+1)-ipoint)/2880
        jpoint=hdstrt(ibuff,maxhdu(ibuff)+1)-2880

!           move all the blocks, one at a time, starting at end of file and
!           working back to the insert position
        do 110 i=1,nshift

!               move to the read start position
            call ftmbyt(ounit,jpoint,.false.,status)

!               read one 2880-byte FITS logical record
            call ftgcbf(ounit,2880,buff,status)

!               move forward to the write start postion
            call ftmbyt(ounit,jpoint+nblock*2880,.true.,status)

!               write the 2880-byte FITS logical record
            call ftpcbf(ounit,2880,buff,status)

!               check for error
            if (status > 0)then
                call ftpmsg('Error inserting empty FITS block(s) '// &
                '(FTIBLK)')
                return
            end if
            jpoint=jpoint-2880
110         continue

        do i=1,2880
            buff(1)(i:i)=cfill
        end do

!           move back to the write start postion
        call ftmbyt(ounit,ipoint,.true.,status)

        do 130 i=1,nblock
!               write the 2880-byte FITS logical record
            call ftpcbf(ounit,2880,buff,status)
130         continue
    end if

    if (hdrdat == 0)then
!               recalculate the starting location of the current data unit
            dtstrt(ibuff)=dtstrt(ibuff)+2880*nblock
    end if

!       recalculate the starting location of all subsequent HDUs
    do 140 i=chdu(ibuff)+1,maxhdu(ibuff)+1
                hdstrt(ibuff,i)=hdstrt(ibuff,i)+2880*nblock
140     continue
    if (status > 0)then
        call ftpmsg('Error inserting FITS block(s) (FTIBLK)')
    end if
end
subroutine fticls(iunit,fstcol,ncols,ttype,tform,status)
!
!*******************************************************************************
!
!! FTICLS inserts one or more new columns into an existing table.
!
!     iunit   i  Fortran I/O unit number
!     fstcol  i  number (position) for the new column; 1 = first column
!                  any existing columns will be moved up NCOLS positions
!     ncols   I  number of columns to insert
!     ttype   c  array of column names (values for TTYPEn keyword)
!     tform   c  array of column formats (values for TFORMn keyword)
!     status  i  returned error status (0=ok)

  integer iunit,fstcol,ncols,status
  character ( len = * ) ttype(*),tform(*)

!
  integer nf,nb,ne
  parameter (nb = 20)
  parameter (nf = 3000)
  parameter (ne = 512)
  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
  integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
  integer theap
  double precision tscale,tzero
  common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

  integer ibuff,colnum,typhdu,datcod,repeat,width,decims,delbyt
  integer naxis1,naxis2,size,freesp,nblock,tflds,tbc,fstbyt,i
  character comm*70,tfm*30,keynam*8

  if (status > 0)return

!     define the number of the buffer used for this file
  ibuff=bufnum(iunit)

!     test that the CHDU is an ASCII table or BINTABLE
  typhdu=hdutyp(ibuff)
  if (typhdu /= 1 .and. typhdu /= 2)then
          status=235
          call ftpmsg('Can only append column to TABLE or '// &
          'BINTABLE extension (FTICOL)')
          return
      end if

!     check that the column number is valid
  tflds=tfield(ibuff)
  if (fstcol < 1)then
      status=302
      return
  else if (fstcol > tflds)then
      colnum=tflds+1
  else
      colnum=fstcol
      end if

!     parse the tform values and calc number of bytes to add to each row
!     make sure format characters are in upper case:
  delbyt=0
  do 5 i=1,ncols
      tfm=tform(i)
      call ftupch(tfm)

      if (typhdu == 1)then
          call ftasfm(tfm,datcod,width,decims,status)
!             add one space between the columns
          delbyt=delbyt+width+1
      else
          call ftbnfm(tfm,datcod,repeat,width,status)
          if (datcod == 1)then
!                 bit column; round up to a multiple of 8 bits
              delbyt=delbyt+(repeat+7)/8
          else if (datcod == 16)then
!                 ASCII string column
              delbyt=delbyt+repeat
           else
!                numerical data type
              delbyt=delbyt+(datcod/10)*repeat
          end if
      end if
5     continue

!     quit on error, or if column is zero byte wide (repeat=0)
  if (status > 0 .or. delbyt == 0)return

!     get current size of the table
  naxis1=rowlen(ibuff)
  call ftgkyj(iunit,'NAXIS2',naxis2,comm,status)

!     Calculate how many more FITS blocks (2880 bytes) need to be added
  size=theap(ibuff)+heapsz(ibuff)
  freesp=(delbyt*naxis2) - ((size+2879)/2880)*2880 + size
  nblock=(freesp+2879)/2880

!     insert the needed number of new FITS blocks at the end of the HDU
  if (nblock > 0)call ftiblk(iunit,nblock,1,status)

!     shift the heap down, and update pointers to start of heap
  size=delbyt*naxis2
  call fthpdn(iunit,size,status)

!     calculate byte position in the row where to insert the new column
  if (colnum > tflds)then
      fstbyt=naxis1
  else
      fstbyt=tbcol(colnum+tstart(ibuff))
      end if

!     insert DELBYT bytes in every row, at byte position FSTBYT
  call ftcins(iunit,naxis1,naxis2,delbyt,fstbyt,status)

  if (typhdu == 1)then
!         adjust the TBCOL values of the existing columns
      do 10 i=1,tflds
          call ftkeyn('TBCOL',i,keynam,status)
          call ftgkyj(iunit,keynam,tbc,comm,status)
          if (tbc > fstbyt)then
               tbc=tbc+delbyt
               call ftmkyj(iunit,keynam,tbc,'&',status)
          end if
10        continue
      end if

!     update the mandatory keywords
  call ftmkyj(iunit,'TFIELDS',tflds+ncols,'&',status)
  call ftmkyj(iunit,'NAXIS1',naxis1+delbyt,'&',status)

!     increment the index value on any existing column keywords
  call ftkshf(iunit,colnum,tflds,ncols,status)

!     add the required keywords for the new columns
  do 15 i=1,ncols
      comm='label for field'
      call ftpkns(iunit,'TTYPE',colnum,1,ttype(i),comm,status)

      comm='format of field'
      tfm=tform(i)
      call ftupch(tfm)
      call ftpkns(iunit,'TFORM',colnum,1,tfm,comm,status)

      if (typhdu == 1)then
          comm='beginning column of field '
          if (colnum == tflds+1)then
!                 allow for the space between preceding column
              tbc=fstbyt+2
!                 set tflds 0, so this branch will not be executed again
          else
              tbc=fstbyt+1
          end if
          call ftpknj(iunit,'TBCOL',colnum,1,tbc,comm,status)

!             increment the column starting position for the next column
          call ftasfm(tfm,datcod,width,decims,status)
!             add one space between the columns
          fstbyt=fstbyt+width+1
      end if

      colnum=colnum+1
15    continue

!     parse the header to initialize the new table structure
  call ftrdef(iunit,status)
end
subroutine fticol(iunit,numcol,ttype,tform,status)
!
!*******************************************************************************
!
!! FTICOL inserts a new column into an existing table.
!
!       iunit   i  Fortran I/O unit number
!       numcol  i  number (position) for the new column; 1 = first column
!                  any existing columns will be moved up one position
!       ttype   c  name of column (value for TTYPEn keyword)
!       tform   c  column format (value for TFORMn keyword)
!       status  i  returned error status (0=ok)

    integer iunit,numcol,status
    character ( len = * ) ttype,tform

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,colnum,typhdu,datcod,repeat,width,decims,delbyt
    integer naxis1,naxis2,size,freesp,nblock,tflds,tbc,fstbyt,i
    character comm*70,tfm*30,keynam*8

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       test that the CHDU is an ASCII table or BINTABLE
    typhdu=hdutyp(ibuff)
    if (typhdu /= 1 .and. typhdu /= 2)then
            status=235
            call ftpmsg('Can only append column to TABLE or '// &
            'BINTABLE extension (FTICOL)')
            return
    end if

!       check that the column number is valid
    tflds=tfield(ibuff)
    if (numcol < 1)then
        status=302
        return
    else if (numcol > tflds)then
        colnum=tflds+1
    else
        colnum=numcol
    end if

!       parse the tform value and calc number of bytes to add to each row
!       make sure format characters are in upper case:
    tfm=tform
    call ftupch(tfm)

    if (typhdu == 1)then
        call ftasfm(tfm,datcod,width,decims,status)
!           add one space between the columns
        delbyt=width+1
    else
        call ftbnfm(tfm,datcod,repeat,width,status)
        if (datcod == 1)then
!               bit column; round up to a multiple of 8 bits
            delbyt=(repeat+7)/8
        else if (datcod == 16)then
!               ASCII string column
            delbyt=repeat
        else
!               numerical data type
            delbyt=(datcod/10)*repeat
        end if
    end if

!       quit on error, or if column is zero byte wide (repeat=0)
    if (status > 0 .or. delbyt == 0)return

!       get current size of the table
    naxis1=rowlen(ibuff)
    call ftgkyj(iunit,'NAXIS2',naxis2,comm,status)

!       Calculate how many more FITS blocks (2880 bytes) need to be added
    size=theap(ibuff)+heapsz(ibuff)
    freesp=(delbyt*naxis2) - ((size+2879)/2880)*2880 + size
    nblock=(freesp+2879)/2880

!       insert the needed number of new FITS blocks at the end of the HDU
    if (nblock > 0)call ftiblk(iunit,nblock,1,status)

!       shift the heap down, and update pointers to start of heap
    size=delbyt*naxis2
    call fthpdn(iunit,size,status)

!       calculate byte position in the row where to insert the new column
    if (colnum > tflds)then
        fstbyt=naxis1
    else
        fstbyt=tbcol(colnum+tstart(ibuff))
    end if

!       insert DELBYT bytes in every row, at byte position FSTBYT
    call ftcins(iunit,naxis1,naxis2,delbyt,fstbyt,status)

    if (typhdu == 1)then
!           adjust the TBCOL values of the existing columns
        do 10 i=1,tflds
            call ftkeyn('TBCOL',i,keynam,status)
            call ftgkyj(iunit,keynam,tbc,comm,status)
            if (tbc > fstbyt)then
                 tbc=tbc+delbyt
                 call ftmkyj(iunit,keynam,tbc,'&',status)
            end if
10          continue
    end if

!       update the mandatory keywords
    call ftmkyj(iunit,'TFIELDS',tflds+1,'&',status)
    call ftmkyj(iunit,'NAXIS1',naxis1+delbyt,'&',status)

!       increment the index value on any existing column keywords
    call ftkshf(iunit,colnum,tflds,1,status)

!       add the required keywords for the new column
    comm='label for field'
    call ftpkns(iunit,'TTYPE',colnum,1,ttype,comm,status)

    comm='format of field'
    call ftpkns(iunit,'TFORM',colnum,1,tfm,comm,status)

    if (typhdu == 1)then
        comm='beginning column of field '
        if (colnum == tflds+1)then
!               allow for the space between preceding column
            tbc=fstbyt+2
        else
            tbc=fstbyt+1
        end if
        call ftpknj(iunit,'TBCOL',colnum,1,tbc,comm,status)
    end if

!       parse the header to initialize the new table structure
    call ftrdef(iunit,status)
end
subroutine ftiimg(ounit,bitpix,naxis,naxes,status)
!
!*******************************************************************************
!
!! FTIIMG inserts an IMAGE extension following the current HDU.
!
!       ounit   i  fortran output unit number
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       status  i  returned error status (0=ok)

    integer ounit,bitpix,naxis,naxes(*),status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

    integer ibuff,nhdu,i,savstr,nblock

    if (status > 0)return
    ibuff=bufnum(ounit)

    if (chdu(ibuff) == 1)then
      if ( hdend(ibuff) == hdstrt(ibuff,chdu(ibuff)) )then
!           Nothing has been written to the file yet, so write primary array
        call ftphpr(ounit,.true., bitpix,naxis,naxes,0,1, &
                    .true.,status)
       return
      end if
    end if

!       close the current HDU to make sure END and fill values are written
    call ftchdu(ounit,status)
    if (status > 0)return

!       save the starting address of the next HDU
    nhdu=chdu(ibuff)+1
    savstr=hdstrt(ibuff,nhdu)

!       define a fake CHDU with a one block header
    dtstrt(ibuff)=hdstrt(ibuff,chdu(ibuff))+2880

!       define the size of the new HDU (this modifies hdstrt(ibuff,nhdu))
    call ftpdef(ounit,bitpix,naxis,naxes,0,1,status)

!       use start of next HDU to calc. how big this new HDU is
    nblock=(hdstrt(ibuff,nhdu)-hdstrt(ibuff,nhdu-1))/2880

!       reset the start of the next HDU back to it original value
    hdstrt(ibuff,nhdu)=savstr

!       insert the required number of blocks at the end of the real CHDU
!       (first define hdutyp so that the correct fill value will be used)
    hdutyp(ibuff)=0
    call ftiblk(ounit,nblock,1,status)
    if (status > 0)return

!       increment the number of HDUs in the file and their starting address
    maxhdu(ibuff)=maxhdu(ibuff)+1
    do 10 i=maxhdu(ibuff),nhdu,-1
            hdstrt(ibuff,i+1)=hdstrt(ibuff,i)
10      continue

!       again, reset the start of the next HDU back to it original value
    hdstrt(ibuff,nhdu)=savstr

!       flush the buffers holding data for the old HDU
    call ftflsh(ibuff,status)

!       recover common block space containing column descriptors for old HDU
    call ftfrcl(ounit,status)

!       move to the new (empty) HDU
    chdu(ibuff)=nhdu

!       set parameters describing an empty 1 block header
    hdutyp(ibuff)=0
    nxthdr(ibuff)=hdstrt(ibuff,nhdu)
    hdend(ibuff)= hdstrt(ibuff,nhdu)
    dtstrt(ibuff)=hdstrt(ibuff,nhdu)+2880

!       write the header keywords
    call ftphpr(ounit,.true.,bitpix,naxis,naxes,0,1,.true.,status)

!       define the structure of the new HDU
    call ftpdef(ounit,bitpix,naxis,naxes,0,1,status)
end
subroutine ftikyd(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTIKYD inserts a double E keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim

    character value*35,key*8,com*47
    character*80 record
    integer nkeys,keypos,vlen

    if (status > 0)return

!       convert double to F format character string and construct the record
    call ftd2e(dval,decim,value,vlen,status)
    key=keywrd
    com=comm
    record=key//'= '//value(1:vlen)//' / '//com

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,record,status)
end
subroutine ftikye(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTIKYE inserts a real*4 E keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993

    character ( len = * ) keywrd,comm
    integer ounit,status,decim
    real rval

    character value*20,key*8,com*47
    character*80 record
    integer nkeys,keypos

    if (status > 0)return

!       convert real to F format character string and construct the full record
    call ftr2e(rval,decim,value,status)
    key=keywrd
    com=comm
    record=key//'= '//value//' / '//com

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,record,status)
end
subroutine ftikyf(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTIKYF inserts a real*4 F keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993

    character ( len = * ) keywrd,comm
    integer ounit,status,decim
    real rval

    character value*20,key*8,com*47
    character ( len = 80 ) record
    integer nkeys,keypos

    if (status > 0)return

!       convert real to F format character string and construct the full record
    call ftr2f(rval,decim,value,status)
    key=keywrd
    com=comm
    record=key//'= '//value//' / '//com

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,record,status)
end
subroutine ftikyg(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTIKYG inserts a double F keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993

    character ( len = * ) keywrd,comm
    integer ounit,status,decim
    double precision dval

    character value*20,key*8,com*47
    character*80 record
    integer nkeys,keypos

    if (status > 0)return

!       convert double to F format character string and construct the record
    call ftd2f(dval,decim,value,status)
    key=keywrd
    com=comm
    record=key//'= '//value//' / '//com

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,record,status)
end
subroutine ftikyj(ounit,keywrd,intval,comm,status)
!
!*******************************************************************************
!
!! FTIKYJ inserts an integer keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       intval  i  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993

    character ( len = * ) keywrd,comm
    integer ounit,status,intval

    character value*20,key*8,com*47
    character ( len = 80 ) record
    integer nkeys,keypos

    if (status > 0)return

!       convert integer to character string and construct the full record
    call fti2c(intval,value,status)
    key=keywrd
    com=comm
    record=key//'= '//value//' / '//com

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,record,status)
end
subroutine ftikyl(ounit,keywrd,logval,comm,status)
!
!*******************************************************************************
!
!! FTIKYL inserts a logical keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       logval  l  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993

    character ( len = * ) keywrd,comm
    integer ounit,status
    logical logval

    character value*20,key*8,com*47
    character ( len = 80 ) record
    integer nkeys,keypos

    if (status > 0)return

!       convert logical to character string and construct the full record
    call ftl2c(logval,value,status)
    key=keywrd
    com=comm
    record=key//'= '//value//' / '//com

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,record,status)
end
subroutine ftikys(ounit,keywrd,strval,comm,status)
!
!*******************************************************************************
!
!! FTIKYS inserts a string keyword into the header at the current position.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       strval  c  keyword value
!       comm    c  keyword comment
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, March 1993
!       Modifed 9/94 to call FTPKLS, supporting the OGIP long string convention

    character ( len = * ) keywrd,comm,strval
    integer ounit,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer lenval,length,i,nspace,ibuff,nexthd,endhd,nkeys,keypos

    if (status > 0)return

!       find how many keywords are required to write the string, in case it
!       cannot fit onto one keyword and has to be continued on multiple lines.

    lenval=len(strval)
    length=0
    do i=lenval,1,-1
            if (strval(i:i) /= ' ')then
                    length=i
                    exit
            end if
    end do
    nspace=max(1,(length-2)/67+1)

!       save current pointer values
    ibuff=bufnum(ounit)
endhd=hdend(ibuff)
    nexthd=nxthdr(ibuff)

!       insert enough spaces in the header at the current location
    call ftghps(ounit,nkeys,keypos,status)

    do 30 i=1,nspace
        call ftirec(ounit,keypos,' ',status)
30      continue

!       temporarily reset position of the end of header to force keyword
!       to be written at the current header position.
    hdend(ibuff)=nexthd

!       write the keyword (supporting the OGIP long string convention)
    call ftpkls(ounit,keywrd,strval,comm,status)

!       reset the next keyword pointer to follow the inserted keyword
    nxthdr(ibuff)=nexthd+80*nspace

!       reset the end-of-header pointer to its real location
    hdend(ibuff)=endhd+80*nspace
end
subroutine ftikyu(ounit,keywrd,comm,status)
!
!*******************************************************************************
!
!! FTIKYU inserts a null-valued keyword to a header record.
!
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1997

    character ( len = * ) keywrd,comm
    integer ounit,status
    character keynam*8,card*80
    integer nkeys,keypos

    if (status > 0)return

    keynam=keywrd
    card=keynam//'=                      / '//comm

    call ftghps(ounit,nkeys,keypos,status)
    call ftirec(ounit,keypos,card,status)
end
subroutine ftinit(funit,fname,block,status)
!
!*******************************************************************************
!
!! FTINIT opens a new FITS file with write access.
!
!       funit   i  Fortran I/O unit number
!       fname   c  name of file to be opened
!       block   i  input record length blocking factor
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer funit,status,block,strlen,i
    character ( len = * ) fname

    if (status > 0)return

!       ignore any leading blanks in the file name
    strlen=len(fname)
    do 10 i=1,strlen
        if (fname(i:i) /= ' ')then

!               call the machine dependent routine which creates the file
            call ftopnx(funit,fname(i:),1,1,block,status)
            if (status > 0)then
     call ftpmsg('FTINIT failed to create the following new file:')
     call ftpmsg(fname)
                return
            end if

!               set column descriptors as undefined
            call ftfrcl(funit,-999)

!               set current column name buffer as undefined
            call ftrsnm
            return
        end if
10      continue

!       if we got here, then the input filename was all blanks
    status=105
    call ftpmsg('FTINIT: Name of file to create is blank.')
end
subroutine ftirec(ounit,pos,record,status)
!
!*******************************************************************************
!
!! FTIREC inserts a 80 character keyword record into the header.
!
!  The keyword record is inserted at the pos-th keyword
!       position (i.e., immediately before the current keyword at position POS.
!
!       ounit   i  fortran output unit number
!       pos     i  keyword will be inserted at this position (1 = 1st keyword)
!       record  c*80  keyword record
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Jan 1995

    character ( len = * ) record
    integer ounit,pos,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    character*80 outrec, inrec
    integer ibuff, fkey, lkey, i, nexthd, nkey

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

!       calculate number of existing keywords
    nkey=(hdend(ibuff)-hdstrt(ibuff,chdu(ibuff)))/80

    if (pos == nkey+1)then
!               simply append the record to the header
            call ftprec(ounit,record,status)
            return
    else if (pos < 1 .or. pos >  nkey)then
            status=203
            return
    end if

    outrec=record

!       move to the insert position
    nexthd=hdstrt(ibuff,chdu(ibuff))+(pos-1)*80
    call ftmbyt(ounit,nexthd,.false.,status)
    nxthdr(ibuff)=nexthd

!       calculated the first and last keyword to be rewritten
    fkey=pos
    lkey=fkey + (hdend(ibuff)-nexthd)/80 - 1

!       now sequentially read each keyword and overwrite it with the previous
    do i=fkey,lkey
            call ftgrec(ounit,i,inrec,status)
            call ftmodr(ounit,outrec,status)
            outrec=inrec
    end do

!       finally, write the last keyword
    call ftprec(ounit,outrec,status)

!       reset the next keyword pointer to follow the inserted keyword
    nxthdr(ibuff)=nexthd+80
end
subroutine ftirow(iunit,frow,nrows,status)
!
!*******************************************************************************
!
!! FTIROW inserts NROWS blank rows immediately after row FROW.
!
!       iunit   i  Fortran I/O unit number
!       frow    i  row number after which the new rows will be inserted.
!                  Specify  0 to add rows to the beginning of the table.
!       nrows   i  number of rows to add to the table (must be greater than 0)
!       status  i  returned error status (0=ok)

    integer iunit,frow,nrows,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,naxis1,naxis2,size,freesp,nblock
    character comm*8

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       test that the CHDU is an ASCII table or BINTABLE
    if (hdutyp(ibuff) /= 1 .and. hdutyp(ibuff) /= 2)then
            status=235
            call ftpmsg('Can only add rows to TABLE or BINTABLE '// &
            'extension (FTIROW)')
            return
    end if

    if (nrows < 0)then
             status=306
             call ftpmsg('Cannot insert negative number of ' // &
             'rows in the table (FTIROW)')
             return
    else if (nrows == 0)then
             return
    end if

!       get current size of the table
    call ftgkyj(iunit,'NAXIS1',naxis1,comm,status)
    call ftgkyj(iunit,'NAXIS2',naxis2,comm,status)

    if (frow > naxis2)then
            status=307
            call ftpmsg('Insert position is greater than the '// &
              'number of rows in the table (FTIROW)')
            return
    else if (frow < 0)then
            status=307
            call ftpmsg('Insert starting row number is less than 0' &
            //' (FTIROW)')
            return
    end if

!       Calculate how many more FITS blocks (2880 bytes) need to be added
    size=theap(ibuff)+heapsz(ibuff)
    freesp=((size+2879)/2880)*2880 - size
    size=naxis1*nrows-freesp
    nblock=(size+2879)/2880

!       insert the needed number of new FITS blocks
    if (nblock > 0)call ftiblk(iunit,nblock,1,status)

!       shift the heap down, and update pointers to start of heap
    size=naxis1*nrows
    call fthpdn(iunit,size,status)

!       shift the rows down
    call ftrwdn(iunit,frow,naxis2,nrows,status)

!       update the NAXIS2 keyword
    naxis2=naxis2+nrows
    call ftmkyj(iunit,'NAXIS2',naxis2,'&',status)
end
subroutine ftitab(ounit,rowlen,nrows,nfield,ttype,tbcol, &
                      tform,tunit,extnam,status)
!
!*******************************************************************************
!
!! FTITAB inserts an ASCII table extension following the current HDU.
!
!       ounit   i  fortran output unit number
!       rowlen  i  width of a row, in characters
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       ttype   c  name of each field (array) (optional)
!       tform   c  format of each field (array)
!       tunit   c  units of each field (array) (optional)
!       extnam  c  name of table extension (optional)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,rowlen,nrows,nfield,tbcol(*),status
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

    integer ibuff,nhdu,i,savstr,nblock,hsize,nkey

    if (status > 0)return
    ibuff=bufnum(ounit)

!       close the current HDU to make sure END and fill values are written
    call ftchdu(ounit,status)
    if (status > 0)return

!       save the starting address of the next HDU
    nhdu=chdu(ibuff)+1
    savstr=hdstrt(ibuff,nhdu)

!       count number of optional TUNITS keywords to be written
    nkey=0
    do i=1,nfield
           if (tunit(i) /= ' ')nkey=nkey+1
    end do
    if (extnam /= ' ')nkey=nkey+1

!       calc min size of header
    nblock=(9 + 3*nfield + nkey +35)/36
    hsize=nblock*2880

!       define a fake CHDU with minimum header
    dtstrt(ibuff)=hdstrt(ibuff,chdu(ibuff))+hsize

!       define the size of the new HDU (this modifies hdstrt(ibuff,nhdu))
    call ftadef(ounit,rowlen,nfield,tbcol,tform,nrows,status)

!       use start of next HDU to calc. how big this new HDU is
    nblock=(hdstrt(ibuff,nhdu)-hdstrt(ibuff,nhdu-1))/2880

!       reset the start of the next HDU back to it original value
    hdstrt(ibuff,nhdu)=savstr

!       insert the required number of blocks at the end of the real CHDU
!       (first define hdutyp so that the correct fill value will be used)
    hdutyp(ibuff)=1
    call ftiblk(ounit,nblock,1,status)
    if (status > 0)return

!       increment the number of HDUs in the file and their starting address
    maxhdu(ibuff)=maxhdu(ibuff)+1
    do 10 i=maxhdu(ibuff),nhdu,-1
            hdstrt(ibuff,i+1)=hdstrt(ibuff,i)
10      continue

!       again, reset the start of the next HDU back to it original value
    hdstrt(ibuff,nhdu)=savstr

!       flush the buffers holding data for the old HDU
    call ftflsh(ibuff,status)

!       recover common block space containing column descriptors for old HDU
    call ftfrcl(ounit,status)

!       move to the new (empty) HDU
    chdu(ibuff)=nhdu

!       set parameters describing an empty header
    hdutyp(ibuff)=1
    nxthdr(ibuff)=hdstrt(ibuff,nhdu)
    hdend(ibuff)= hdstrt(ibuff,nhdu)
    dtstrt(ibuff)=hdstrt(ibuff,nhdu)+hsize

!       write the header keywords
    call ftphtb(ounit,rowlen,nrows,nfield,ttype,tbcol,tform,tunit, &
                extnam,status)

!       define the structure of the new HDU
    call ftadef(ounit,rowlen,nfield,tbcol,tform,nrows,status)

end
subroutine ftkeyn(keywrd,nseq,keyout,status)
!
!*******************************************************************************
!
!! FTKEYN makes a keyword name by from the root name and a sequence number.
!
!       keywrd  c  root keyword name
!       nseq    i  sequence number
!       OUTPUT PARAMETERS:
!       keyout  c  output concatinated keyword name
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, February 1991

    character ( len = * ) keywrd,keyout
    integer nseq,status,nspace,i
    character value*20,work*8

    work=keywrd

!       find end of keyword string
    nspace=1
    do i=1,8
            if (work(i:i) == ' ')then
              exit
            end if
            nspace=nspace+1
    end do

!       append sequence number to keyword root only if there is room
    if (nseq < 0)then
!               illegal value
            go to 900
    else if (nseq < 10 .and. nspace <= 8)then
            write(work(nspace:nspace),1001,err=900)nseq
    else if (nseq < 100 .and. nspace <= 7)then
            write(work(nspace:nspace+1),1002,err=900)nseq
    else if (nseq < 1000 .and. nspace <= 6)then
            write(work(nspace:nspace+2),1003,err=900)nseq
    else if (nseq < 10000 .and. nspace <= 5)then
            write(work(nspace:nspace+3),1004,err=900)nseq
    else if (nseq < 100000 .and. nspace <= 4)then
            write(work(nspace:nspace+4),1005,err=900)nseq
    else if (nseq < 1000000 .and. nspace <= 3)then
            write(work(nspace:nspace+5),1006,err=900)nseq
    else if (nseq < 10000000 .and. nspace <= 2)then
            write(work(nspace:nspace+6),1007,err=900)nseq
    else
!               number too big to fit in keyword
            go to 900
    end if

1001    format(i1)
1002    format(i2)
1003    format(i3)
1004    format(i4)
1005    format(i5)
1006    format(i6)
1007    format(i7)

    keyout=work
    return
!       come here if error concatinating the seq. no. to the root string
900     continue

    if (status > 0)return
    status=206
    write(value,1008)nseq
1008    format(i20)
    call ftpmsg('Could not concatinate the integer '//value// &
   ' to the root keyword named: '//work)
end
subroutine ftkshf(iunit,colmin,colmax,incre,status)
!
!*******************************************************************************
!
!! FTKSHF shifts the index value on any existing column keywords.
!
!       This routine will modify the name of any keyword that begins with 'T'
!       and has an index number in the range COLMIN - COLMAX, inclusive.

!       if incre is positive, then the index values will be incremented.
!       if incre is negative, then the kewords with index = COLMIN
!       will be deleted and the index of higher numbered keywords will
!       be decremented.

!       iunit   i  Fortran I/O unit number
!       colmin  i  starting column number to be incremented
!       colmax  i  maximum column number to be increment
!       incre   i  amount by which the index value should be shifted
!       status  i  returned error status (0=ok)

    integer iunit,colmin,colmax,incre,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,typhdu,tflds,nkeys,nmore,nrec,ival,tstat,i1
    character rec*80,newkey*8,q*4

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       test that the CHDU is an ASCII table or BINTABLE
    typhdu=hdutyp(ibuff)
    if (typhdu /= 1 .and. typhdu /= 2)then
            status=235
            call ftpmsg('Can only operate on TABLE or '// &
            'BINTABLE extension (FTKSHF)')
            return
    end if

!       test column number limits
    tflds=tfield(ibuff)
    if (colmin < 1 .or. colmax < 1)then
         status=302
         return
    else if (colmin > colmax .or. colmin > tflds)then
         return
    end if

!       get the number of keywords in the header
    call ftghsp(iunit,nkeys,nmore,status)

!       go thru header starting with the 9th keyword looking for 'TxxxxNNN'

    nrec=9
100     call ftgrec(iunit,nrec,rec,status)

    if (rec(1:1) == 'T')then
        q=rec(2:5)
        i1=6

!           search list of 5-character 'official' indexed keywords
        if ( q == 'BCOL' .or. q == 'FORM' .or. q == 'TYPE' &
        .or. q == 'UNIT' .or. q == 'NULL' .or. q == 'SCAL' &
        .or. q == 'ZERO' .or. q == 'DISP')go to 20

!           search list of 5-character 'local' indexed keywords
        if ( q == 'LMIN' .or. q == 'LMAX' .or. q == 'DMIN' &
        .or. q == 'DMAX' .or. q == 'CTYP' .or. q == 'CRPX' &
        .or. q == 'CRVL' .or. q == 'CDLT' .or. q == 'CROT' &
        .or. q == 'CUNI')go to 20

        q=rec(1:4)
        i1=5
!           search list of 4-character 'official' indexed keywords
        if (q == 'TDIM')go to 20

!           no match so go on to next keyword
        go to 90

20          continue
!           try reading the index number suffix
        tstat=0
        call ftc2ii(rec(i1:8),ival,tstat)
        if (tstat == 0 .and. ival >= colmin .and. &
            ival <= colmax)then
            if (incre <= 0 .and. ival == colmin)then
!                   delete keyword related to this column
                call ftdrec(iunit,nrec,status)
                nkeys=nkeys-1
                nrec=nrec-1
            else
                ival=ival+incre
                i1=i1-1
                call ftkeyn(rec(1:i1),ival,newkey,status)
                rec(1:8)=newkey
!                   modify the index number of this keyword
                call ftmrec(iunit,nrec,rec,status)
            end if
        end if
    end if

90      nrec=nrec+1
    if (nrec <= nkeys)go to 100
end
subroutine ftl2c(lval,cval,status)
!
!*******************************************************************************
!
!! FTL2C converts a logical value to a C*20 right justified character string.
!
    integer status
    logical lval
    character*20 cval

    if (status > 0)return

    if (lval)then
            cval='                   T'
    else
            cval='                   F'
    end if
end
subroutine ftmahd(iunit,extno,xtend,status)
!
!*******************************************************************************
!
!! FTMAHD moves the i/o pointer to the specified HDU.
!
!  It also initializes all
!       the common block parameters which describe the extension

!       iunit   i  fortran unit number
!       extno   i  number of the extension to point to.
!       xtend   i  returned type of extension:   0 = the primary HDU
!                                                1 = an ASCII table
!                                                2 = a binary table
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991

    integer iunit,extno,xtend,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,movto,tstat

    if (status > 0)then
        return
    else if (extno <= 0 .or. extno >= ne)then
        status=301
        return
    end if

    ibuff=bufnum(iunit)

!       check if we are already positioned to the correct HDU
    if (extno == chdu(ibuff))then
!           just return the type of extension
        xtend=hdutyp(ibuff)
    else

!           now move to the extension, or the highest one we know about
10          movto=min(extno,maxhdu(ibuff)+1)

!           before closing out the CHDU, make sure the new extension exists
        call ftmbyt(iunit,hdstrt(ibuff,movto),.false.,status)
        if (status > 0)return

!           close out the current HDU before moving to the new one
        call ftchdu(iunit,status)
        if (status > 0)then
            call ftpmsg('FTMAHD could not close the'// &
                ' current HDU before moving to the new HDU.')
            return
        end if

        call ftgext(iunit,movto,xtend,status)
        if (status > 0)then
!               failed to move to new extension, so restore previous extension
            tstat=0
            call ftrhdu(iunit,movto,tstat)
            return
        end if

!           continue reading extensions until we get to the one we want
        if (movto < extno)go to 10
    end if
end
subroutine ftmcom(ounit,keywrd,comm,status)
!
!*******************************************************************************
!
!! FTMCOM modifies a the comment string in a header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       comm    c  new keyword comment (max of 72 characters long)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    character ( len = * ) keywrd,comm
    integer ounit,status,lenval,ncomm
    character value*80,knam*8,cmnt*72

    if (status > 0)return

    knam=keywrd

!       find the old keyword + value string
    call ftgcrd(ounit,knam,value,status)
    if (status == 202)then
      call ftpmsg('FTMCOM Could not find the '//knam//' keyword.')
      return
    end if

    call ftprsv(value,lenval,status)

    cmnt=comm

!       find amount of space left for comment string (3 spaces needed for ' / ')
    ncomm=77-lenval

!       write the keyword record if there is space
    if (ncomm > 0)then
      call ftmodr(ounit, &
      value(1:lenval)//' / '//cmnt(1:ncomm),status)
    end if
end
subroutine ftmcrd(ounit,keywrd,card,status)
!
!*******************************************************************************
!
!! FTMCRD modifies (overwrites) a given header record specified by keyword name.
!
!       This can be used to overwrite the name of the keyword as well as
!       the value and comment fields.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       card    c  new 80-character card image to be written
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    character ( len = * ) keywrd,card
    integer ounit,status
    character value*80

    if (status > 0)return

!       find the old keyword string
    call ftgcrd(ounit,keywrd,value,status)

    value=card

!       make sure new keyword name is in upper case
    call ftupch(value(1:8))

!       test that keyword name contains only legal characters
    call fttkey(value(1:8),status)

!       write the new keyword record
    call ftmodr(ounit,value,status)
end
subroutine ftmkey(ounit,keywrd,value,comm,status)
!
!*******************************************************************************
!
!! FTMKEY modifies an existing simple FITS keyword record.
!
!  The format used is:
!            "KEYWORD = VALUE / COMMENT"
!               VALUE is assumed to be 20 characters long
!               COMMENT is assumed to be 47 characters long
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       value   c  keyword value   (20 characters, cols. 11-30)
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,value,comm
    integer ounit,status
    character key*8, val*20, com*47

    key=keywrd
    val=value
    com=comm

!       overwrite the preceeding 80 characters to the output buffer:
    call ftmodr(ounit,key//'= '//val//' / '//com,status)
end
subroutine ftmkyd(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTMKYD modifies a double precision value header record in E format.
!
!       If it will fit, the value field will be 20 characters wide;
!       otherwise it will be expanded to up to 35 characters, left
!       justified.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (max. 47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim,vlen
    character value*35,key*8,cmnt*48

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

    key=keywrd
!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

!       convert double precision to E format character string
    call ftd2e(dval,decim,value,vlen,status)

!       write the keyword record
    call ftmodr(ounit,key//'= '//value(1:vlen)//' / '//cmnt,status)
end
subroutine ftmkye(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTMKYE modifies a real*4 value header record in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    real rval
    integer ounit,status,decim
    character value*20,cmnt*48

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

!       convert real to E format character string
    call ftr2e(rval,decim,value,status)

!       modify the keyword record
    call ftmkey(ounit,keywrd,value,cmnt,status)
end
subroutine ftmkyf(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTMKYF modifies a real*4 value header record in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    real rval
    integer ounit,status,decim
    character value*20,cmnt*48

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

!       convert real to F format character string
    call ftr2f(rval,decim,value,status)

!       write the keyword record
    call ftmkey(ounit,keywrd,value,cmnt,status)
end
subroutine ftmkyg(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTMKYG modifies a double precision value header record in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim
    character value*20,cmnt*48

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

!       convert double precision to F format character string
    call ftd2f(dval,decim,value,status)

!       modify the keyword record
    call ftmkey(ounit,keywrd,value,cmnt,status)
end
subroutine ftmkyj(ounit,keywrd,intval,comm,status)
!
!*******************************************************************************
!
!! FTMKYJ modifies an integer value header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       intval  i  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer ounit,status,intval
    character value*20,cmnt*48

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

!       convert integer to character string
    call fti2c(intval,value,status)

!       modify the keyword record
    call ftmkey(ounit,keywrd,value,cmnt,status)
end
subroutine ftmkyl(ounit,keywrd,logval,comm,status)
!
!*******************************************************************************
!
!! FTMKYL modifies a logical value header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       logval  l  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer ounit,status
    logical logval
    character value*20,cmnt*48

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

!       convert logical to character string
    call ftl2c(logval,value,status)

!       modify the keyword record
    call ftmkey(ounit,keywrd,value,cmnt,status)
end
subroutine ftmkys(ounit,keywrd,strval,comm,status)
!
!*******************************************************************************
!
!! FTMKYS modifies a character string value header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       strval  c  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!       modifed 7/93 to support string keywords continued over multiple cards
!       modified 9/94 to support the OGIP long string convention

    character ( len = * ) keywrd,strval,comm
    integer ounit,status

    integer clen,i,nvalue,ncomm
    character keynam*8,value*70,cmnt*48,bslash

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    if (status > 0)return

!       check if the new value is too long to fit in a single 'card image'
    clen=len(strval)
    if (clen <= 68)go to 20

    do i=clen,69,-1
            if (strval(i:i) /= ' ')go to 100
    end do

!       now check that the old keyword is not continued over multiple cards
!       read the (first line) of the existing keyword

20      call ftgkey(ounit,keywrd,value,cmnt,status)
    if (status > 0)go to 900

!       is last character of the value a backslash or & ?
!       have to use 2 \\'s because the SUN compiler treats 1 \ as an escape
    bslash='\\'
    do 30 i=70,1,-1
            if (value(i:i) /= ' '.and. value(i:i)/='''')then
                if (value(i:i) == bslash .or. &
                    value(i:i) == '&')then
!                     backspace the current header pointer by one record
                  nxthdr(bufnum(ounit))=nxthdr(bufnum(ounit))-80
                  go to 100
                else
                  go to 40
                end if
            end if
30      continue

!       OK, we can simply overwrite the old keyword with the new
40      continue

!       overwrite the old comment unless user supplied the magic value
    if (comm /= '&')then
            cmnt=comm
    end if
!       convert string to quoted character string (max length = 70 characters)
    call fts2c(strval,value,clen,status)
    if (status > 0)go to 900

!       find amount of space left for comment string
!       (assume 10 char. for 'keyword = ', and 3 between value and comment)
!       which leaves 67 spaces for the value string + comment string
    nvalue=max(20,clen)
    ncomm=67-nvalue

!       write the keyword record
    keynam=keywrd
    if (ncomm > 0)then
!         there is space for a comment
      call ftmodr(ounit, &
      keynam//'= '//value(1:nvalue)//' / '//cmnt(1:ncomm),status)
    else
!         no room for a comment
      call ftmodr(ounit, &
      keynam//'= '//value(1:nvalue)//'   ',status)
    end if
    go to 900

100     continue

!       Either the old or new keyword is continued over multiple
!       header card images, so have to use a less efficient way to modify
!       the keyword by completely deleting the old and inserting the new.

!       read the old comment, if we need to preserve it
    if (comm == '&')then
            call ftgkys(ounit,keywrd,value,cmnt,status)
            if (status > 0)go to 900
!               reset the current header pointer by 2 records to make
!               it faster (usually) to find and delete the keyword
            nxthdr(bufnum(ounit))=nxthdr(bufnum(ounit))-160
    else
            cmnt=comm
    end if

!       delete the old keyword
    call ftdkey(ounit,keywrd,status)
    if (status > 0)go to 900

!       insert the new keyword
    call ftikys(ounit,keywrd,strval,cmnt,status)

900     continue
end
subroutine ftmkyu(ounit,keywrd,comm,status)
!
!*******************************************************************************
!
!! FTMKYU modifies a null-valued keyword.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1997

    character ( len = * ) keywrd,comm
    integer ounit,status
    character value*80,cmnt*80

    if (status > 0)return

!       find the old keyword
    call ftgkey(ounit,keywrd,value,cmnt,status)

!       check for special symbol indicating that comment should not be changed
    if (comm /= '&')then
          cmnt=comm
    end if

    value = ' '

!       modify the keyword record
    call ftmkey(ounit,keywrd,value,cmnt,status)
end
subroutine ftmnam(ounit,oldkey,newkey,status)
!
!*******************************************************************************
!
!! FTMNAM modifies (overwrites) the name of an existing keyword.
!
!  The routine preserves the current value and comment fields associated
!  with the keyword.
!
!       ounit   i  fortran output unit number
!       oldkey  c  old keyword name    ( 8 characters, cols.  1- 8)
!       newkey  c  new keyword name to be written
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    character ( len = * ) oldkey,newkey
    integer ounit,status
    character card*80

    if (status > 0)return

!       find the old keyword string
    call ftgcrd(ounit,oldkey,card,status)

    card(1:8)=newkey

!       make sure new keyword name is in upper case
    call ftupch(card(1:8))

!       test that keyword name contains only legal characters
    call fttkey(card(1:8),status)

!       write the new keyword record
    call ftmodr(ounit,card,status)
end
subroutine ftmodr(ounit,record,status)
!
!*******************************************************************************
!
!! FTMODR modifies the preceeding 80 character record in the FITS header.
!
!       ounit   i  fortran output unit number
!       record  c  input 80 character header record
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) record
    character*80  rec
    integer ounit,status,ibuff

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

    rec=record

!       make sure keyword name is in upper case
    call ftupch(rec(1:8))

!       test that keyword name contains only legal characters
    call fttkey(rec(1:8),status)

!       move the I/O pointer back to the beginning of the preceeding keyword
    call ftmbyt(ounit,nxthdr(ibuff)-80,.false.,status)

!       overwrite the 80 characters to the output buffer:
    call ftpcbf(ounit,80,rec,status)
end
subroutine ftmrec(ounit,nkey,record,status)
!
!*******************************************************************************
!
!! FTMREC modifies the Nth keyword in the CHU.
!
!  The keyword is replaced with the input 80 character string.
!
!       ounit   i  fortran output unit number
!       nkey    i  sequence number (starting with 1) of the keyword to read
!       record  c  80-character string to replace the record with
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,nkey,status
    character ( len = * ) record
    character rec*80

!       find the old keyword; just use REC as a temporary variable
    call ftgrec(ounit,nkey,rec,status)

    rec=record
!       overwrite the keyword with the new record
    call ftmodr(ounit,rec,status)
end
subroutine ftmrhd(iunit,extmov,xtend,status)
!
!*******************************************************************************
!
!! FTMRHD moves the i/o pointer to the specified HDU.
!
!  It also initializes all
!       the common block parameters which describe the extension

!       iunit   i  fortran unit number
!       extmov  i  number of the extension to point to, relative to the CHDU
!       xtend   i  returned type of extension:   0 = the primary HDU
!                                                1 = an ASCII table
!                                                2 = a binary table
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991

    integer iunit,extmov,xtend,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,extno

    if (status > 0)return

    ibuff=bufnum(iunit)

!       calculate the absolute HDU number, then move to it
    extno=chdu(ibuff)+extmov
    call ftmahd(iunit,extno,xtend,status)
end
subroutine ftnkey(nseq,keywrd,keyout,status)
!
!*******************************************************************************
!
!! FTNKEY makes a keyword name from a sequence number and the root name.
!
!  (Sequence number is prepended to the name)
!
!       nseq    i  sequence number
!       keywrd  c  root keyword name
!       OUTPUT PARAMETERS:
!       keyout  c  output concatinated keyword name
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Aug 1994

    character ( len = * ) keywrd,keyout
    integer nseq,status,nspace,i
    character value*20,work*8

    work=keywrd

!       find end of keyword string
    nspace=0
    do i=8,1,-1
            if (work(i:i) /= ' ')then
              exit
            end if
            nspace=nspace+1
    end do

!       prepend sequence number to keyword root only if there is room
    if (nseq < 0)then
!               illegal value
            go to 900
    else if (nseq < 10 .and. nspace >= 1)then
            write(keyout,1001,err=900)nseq,work(1:7)
    else if (nseq < 100 .and. nspace >= 2)then
            write(keyout,1002,err=900)nseq,work(1:6)
    else if (nseq < 1000 .and. nspace >= 3)then
            write(keyout,1003,err=900)nseq,work(1:5)
    else if (nseq < 10000 .and. nspace >= 4)then
            write(keyout,1004,err=900)nseq,work(1:4)
    else if (nseq < 100000 .and. nspace >= 5)then
            write(keyout,1005,err=900)nseq,work(1:3)
    else if (nseq < 1000000 .and. nspace >= 6)then
            write(keyout,1006,err=900)nseq,work(1:2)
    else if (nseq < 10000000 .and. nspace >= 7)then
            write(keyout,1007,err=900)nseq,work(1:1)
    else
!               number too big to fit in keyword
            go to 900
    end if

1001    format(i1,a7)
1002    format(i2,a6)
1003    format(i3,a5)
1004    format(i4,a4)
1005    format(i5,a3)
1006    format(i6,a2)
1007    format(i7,a1)

    return
!       come here if error concatinating the seq. no. to the root string
900     continue

    if (status > 0)return
    status=206
    write(value,1008)nseq
1008    format(i20)
    call ftpmsg('Could not concatinate the integer '//value// &
   ' and the root keyword named: '//work)
end
subroutine ftnulc(input,np,chktyp,setval,flgray,anynul, &
                      scaled,scale,zero)
!
!*******************************************************************************
!
!! FTNULC checks input complex array for nulls and applies scaling.
!
!       if chktyp=1 then set the undefined pixel = SETVAL
!       if chktyp=2 then set the corresponding FLGRAY = .true.

!       When scaling complex data values,  both the real and imaginary
!       components of the value are scaled by SCALE, but the offset
!       given by ZERO is only applied to the real part of the complex number

!       input   r  input array of values
!       np      i  number of pairs of values
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  r  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       scaled  l  does data need to be scaled?
!       scale   d  scale factor
!       zero    d  offset

    real input(*),setval(2)
    integer np,i,chktyp,j
    double precision scale,zero
    logical flgray(*),anynul,scaled
    logical fttrnn
    external fttrnn

    if (chktyp == 2)then
!               initialize the null flag values
            do 5 i=1,np
                    flgray(i)=.false.
5               continue
    end if

    j=1
    do 10 i=1,np
!               do the real part of the complex number
            if (chktyp /= 0 .and. fttrnn(input(j)))then
                anynul=.true.
                if (chktyp == 1)then
!                               set both parts of the complex number to the
!                               specified special value
                            input(j)=setval(1)
                            input(j+1)=setval(2)
                else
!                               set the corresponding flag value to true
                            flgray(i)=.true.
                end if
                j=j+2
            else if (scaled)then
                input(j)=input(j)*scale+zero
                j=j+1

!                   do the imaginary part of the complex number
                if (chktyp /= 0 .and. fttrnn(input(j)))then
                        anynul=.true.
                        if (chktyp == 1)then
!                               set both parts of the complex number to the
!                               specified special value
                            input(j-1)=setval(1)
                            input(j)=setval(2)
                        else
!                               set the corresponding flag value to true
                            flgray(i)=.true.
                        end if
                else if (scaled)then
                    input(j)=input(j)*scale
                end if
                j=j+1
            else
                j=j+2
            end if
10      continue
end
subroutine ftnulm(input,np,chktyp,setval,flgray,anynul, &
                      scaled,scale,zero)
!
!*******************************************************************************
!
!! FTNULM checks input double complex array for nulls and applies scaling.
!
!       if chktyp=1 then set the undefined pixel = SETVAL
!       if chktyp=2 then set the corresponding FLGRAY = .true.

!       When scaling complex data values,  both the real and imaginary
!       components of the value are scaled by SCALE, but the offset
!       given by ZERO is only applied to the real part of the complex number

!       input   d  input array of values
!       np      i  number of pairs of values
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  d  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       scaled  l  does data need to be scaled?
!       scale   d  scale factor
!       zero    d  offset

    double precision input(*),setval(2)
    integer np,i,chktyp,j
    double precision scale,zero
    logical flgray(*),anynul,scaled
    logical fttdnn
    external fttdnn

    if (chktyp == 2)then
!               initialize the null flag values
            do 5 i=1,np
                    flgray(i)=.false.
5               continue
    end if

    j=1
    do 10 i=1,np
!               do the real part of the complex number
            if (chktyp /= 0 .and. fttdnn(input(j)))then
                anynul=.true.
                if (chktyp == 1)then
!                               set both parts of the complex number to the
!                               specified special value
                            input(j)=setval(1)
                            input(j+1)=setval(2)
                else
!                               set the corresponding flag value to true
                            flgray(i)=.true.
                end if
                j=j+2
            else if (scaled)then
                input(j)=input(j)*scale+zero
                j=j+1

!                   do the imaginary part of the complex number
                if (chktyp /= 0 .and. fttdnn(input(j)))then
                        anynul=.true.
                        if (chktyp == 1)then
!                               set both parts of the complex number to the
!                               specified special value
                            input(j-1)=setval(1)
                            input(j)=setval(2)
                        else
!                               set the corresponding flag value to true
                            flgray(i)=.true.
                        end if
                else if (scaled)then
                    input(j)=input(j)*scale
                end if
                j=j+1
            else
                j=j+2
            end if
10      continue
end
subroutine ftopen(funit,fname,rwmode,block,status)
!
!*******************************************************************************
!
!! FTOPEN opens an existing FITS file with readonly or read/write access.
!
!       funit   i  Fortran I/O unit number
!       fname   c  name of file to be opened
!       rwmode  i  file access mode: 0 = readonly; else = read and write
!       block   i  returned record length blocking factor
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer funit,rwmode,block,status,strlen,i,xtend
    character ( len = * ) fname

    if (status > 0)return

!       ignore any leading blanks in the file name
    strlen=len(fname)
    do 10 i=1,strlen
        if (fname(i:i) /= ' ')then

!               call the machine dependent routine which opens the file
            call ftopnx(funit,fname(i:),0,rwmode,block,status)
            if (status > 0)then
                 call ftpmsg('FTOPEN failed to Find and/or Open'// &
                           ' the following file:')
                 call ftpmsg(fname)
                 return
            end if

!               set column descriptors as undefined
            call ftfrcl(funit,-999)

!               determine the structure and size of the primary HDU
            call ftrhdu(funit,xtend,status)
            if (status > 0)then
              call ftpmsg('FTOPEN could not interpret primary ' &
                //'array header keywords of file:')
              call ftpmsg(fname)
              if (status == 252)then
                  call ftpmsg('Is this a FITS file??')
              end if
            end if

!               set current column name buffer as undefined
            call ftrsnm
            return
        end if
10      continue

!       if we got here, then the input filename was all blanks
    status=104
    call ftpmsg('FTOPEN: Name of file to open is blank.')
    return

end
subroutine ftp2db(ounit,group,dim1,nx,ny,array,status)
!
!*******************************************************************************
!
!! FTP2DB writes a 2-d image of byte values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   c*1  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    character array(dim1,*)
    integer fpixel,row

    fpixel=1
    do 10 row = 1,ny
            call ftpprb(ounit,group,fpixel,nx,array(1,row),status)
            fpixel=fpixel+nx
10      continue

end
subroutine ftp2dd(ounit,group,dim1,nx,ny,array,status)
!
!*******************************************************************************
!
!! FTP2DD writes a 2-d image of r*8 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   d  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    double precision array(dim1,*)
    integer fpixel,row

    fpixel=1
    do 10 row = 1,ny
            call ftpprd(ounit,group,fpixel,nx,array(1,row),status)
            fpixel=fpixel+nx
10      continue

end
subroutine ftp2de(ounit,group,dim1,nx,ny,array,status)
!
!*******************************************************************************
!
!! FTP2DE writes a 2-d image of r*4 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   r  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    real array(dim1,*)
    integer fpixel,row

    fpixel=1
    do 10 row = 1,ny
            call ftppre(ounit,group,fpixel,nx,array(1,row),status)
            fpixel=fpixel+nx
10      continue

end
subroutine ftp2di(ounit,group,dim1,nx,ny,array,status)
!
!*******************************************************************************
!
!! FTP2DI writes a 2-d image of i*2 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   i*2  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    integer*2 array(dim1,*)
    integer fpixel,row

    fpixel=1
    do 10 row = 1,ny
            call ftppri(ounit,group,fpixel,nx,array(1,row),status)
            fpixel=fpixel+nx
10      continue

end
subroutine ftp2dj(ounit,group,dim1,nx,ny,array,status)
!
!*******************************************************************************
!
!! FTP2DJ writes a 2-d image of i*4 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       nx      i  size of the image in the x direction
!       ny      i  size of the image in the y direction
!       array   i  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,nx,ny,status
    integer array(dim1,*)
    integer fpixel,row

    fpixel=1
    do 10 row = 1,ny
            call ftpprj(ounit,group,fpixel,nx,array(1,row),status)
            fpixel=fpixel+nx
10      continue

end
subroutine ftp3db(ounit,group,dim1,dim2,nx,ny,nz,array,status)
!
!*******************************************************************************
!
!! FTP3DB writes a 3-d cube of byte values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   c*1  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    character array(dim1,dim2,*)
    integer fpixel,row,band

    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
        call ftpprb(ounit,group,fpixel,nx,array(1,row,band),status)
        fpixel=fpixel+nx
10      continue
20      continue

end
subroutine ftp3dd(ounit,group,dim1,dim2,nx,ny,nz,array,status)
!
!*******************************************************************************
!
!! FTP3DD writes a 3-d cube of r*8 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   r*8  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    double precision array(dim1,dim2,*)
    integer fpixel,row,band

    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
        call ftpprd(ounit,group,fpixel,nx,array(1,row,band),status)
        fpixel=fpixel+nx
10      continue
20      continue

end
subroutine ftp3de(ounit,group,dim1,dim2,nx,ny,nz,array,status)
!
!*******************************************************************************
!
!! FTP3DE writes a 3-d cube of r*4 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   r  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    real array(dim1,dim2,*)
    integer fpixel,row,band

    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
        call ftppre(ounit,group,fpixel,nx,array(1,row,band),status)
        fpixel=fpixel+nx
10      continue
20      continue

end
subroutine ftp3di(ounit,group,dim1,dim2,nx,ny,nz,array,status)
!
!*******************************************************************************
!
!! FTP3DI writes a 3-d cube of i*2 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   i*2  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    integer*2 array(dim1,dim2,*)
    integer fpixel,row,band

    fpixel=1
    do 20 band=1,nz
      do row = 1,ny
        call ftppri(ounit,group,fpixel,nx,array(1,row,band),status)
        fpixel=fpixel+nx
      end do
20      continue

end
subroutine ftp3dj(ounit,group,dim1,dim2,nx,ny,nz,array,status)
!
!*******************************************************************************
!
!! FTP3DJ writes a 3-d cube of i*4 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       dim1    i  actual first dimension of ARRAY
!       dim2    i  actual second dimension of ARRAY
!       nx      i  size of the cube in the x direction
!       ny      i  size of the cube in the y direction
!       nz      i  size of the cube in the z direction
!       array   i  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,dim1,dim2,nx,ny,nz,status
    integer array(dim1,dim2,*)
    integer fpixel,row,band

    fpixel=1
    do 20 band=1,nz
    do 10 row = 1,ny
        call ftpprj(ounit,group,fpixel,nx,array(1,row,band),status)
        fpixel=fpixel+nx
10      continue
20      continue

end
subroutine ftpbit(setbit,wrbit,buffer)
!
!*******************************************************************************
!
!! FTPBIT encodes bits within the byte as specified by the input logical array.
!
!       The corresponding bit is set to
!       1 if the logical array element is true.  Only the bits
!       between begbit and endbit, inclusive, are set or reset;
!       the remaining bits, if any, remain unchanged.

!       setbit  l  input array of logical data values corresponding
!                  to the bits to be set in the output buffer
!                  TRUE means corresponding bit is to be set.
!       wrbit   l  input array of logical values indicating which
!                  bits in the byte are to be modified.  If FALSE,
!                  then the corresponding bit should remain unchanged.
!       buffer  i  output integer containing the encoded byte
!
!       written by Wm Pence, HEASARC/GSFC, May 1992

    integer buffer,tbuff,outbit
    logical setbit(8),wrbit(8)

    outbit=0
    tbuff=buffer

!       test each of the 8 bits, starting with the most significant
    if (tbuff > 127)then
!           the bit is currently set in the word
        if (wrbit(1) .and. (.not.setbit(1)))then
!                only in this case do we reset the bit
        else
!               in all other cases we want the bit to be set
            outbit=outbit+128
        end if
        tbuff=tbuff-128
    else
!           bit is currently not set; set it only if requested to
        if (wrbit(1) .and. setbit(1))outbit=outbit+128
    end if

    if (tbuff > 63)then
        if (wrbit(2) .and. (.not.setbit(2)))then
        else
            outbit=outbit+64
        end if
        tbuff=tbuff-64
    else
        if (wrbit(2) .and. setbit(2))outbit=outbit+64
    end if

    if (tbuff > 31)then
        if (wrbit(3) .and. (.not.setbit(3)))then
        else
            outbit=outbit+32
        end if
        tbuff=tbuff-32
    else
        if (wrbit(3) .and. setbit(3))outbit=outbit+32
    end if

    if (tbuff > 15)then
        if (wrbit(4) .and. (.not.setbit(4)))then
        else
            outbit=outbit+16
        end if
        tbuff=tbuff-16
    else
        if (wrbit(4) .and. setbit(4))outbit=outbit+16
    end if

    if (tbuff > 7)then
        if (wrbit(5) .and. (.not.setbit(5)))then
        else
            outbit=outbit+8
        end if
        tbuff=tbuff-8
    else
        if (wrbit(5) .and. setbit(5))outbit=outbit+8
    end if

    if (tbuff > 3)then
        if (wrbit(6) .and. (.not.setbit(6)))then
        else
            outbit=outbit+4
        end if
        tbuff=tbuff-4
    else
        if (wrbit(6) .and. setbit(6))outbit=outbit+4
    end if

    if (tbuff > 1)then
        if (wrbit(7) .and. (.not.setbit(7)))then
        else
            outbit=outbit+2
        end if
        tbuff=tbuff-2
    else
        if (wrbit(7) .and. setbit(7))outbit=outbit+2
    end if

    if (tbuff == 1)then
        if (wrbit(8) .and. (.not.setbit(8)))then
        else
            outbit=outbit+1
        end if
    else
        if (wrbit(8) .and. setbit(8))outbit=outbit+1
    end if

    buffer=outbit
end
subroutine ftpbnh(ounit,nrows,nfield,ttype,tform,tunit, &
                      extnam,pcount,status)
!
!*******************************************************************************
!
!! FTPBNH is obsolete.  Call FTPHBN instead.
!

    integer ounit,nrows,nfield,pcount,status
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam

    call ftphbn(ounit,nrows,nfield,ttype,tform,tunit, &
                      extnam,pcount,status)
end
subroutine ftpcks(iunit,status)
!
!*******************************************************************************
!
!! FTPCKS creates or updates the checksum keywords in the CHU.
!
!  These keywords provide a checksum verification of the FITS HDU based on
!  the ASCII coded 1's complement checksum algorithm developed by Rob Seaman
!  at NOAO.
!
!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Sept, 1994

    integer iunit,status

    integer nf,nb,ne
    parameter (nf = 3000)
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)

    double precision sum,dsum,odsum
    integer ibuff,nrec,dd,mm,yy,dummy,i,tstat
    character datstr*8,string*16,comm*40
    character ( len = 16 ) oldcks
    character ( len = 20 ) datsum
    logical complm

    if (status > 0)return

    ibuff=bufnum(iunit)

!       generate current date string to put into the keyword comment
    call ftgsdt(dd,mm,yy,status)
    if (status > 0)return

    datstr='  /  /  '
    write(datstr(1:2),1001)dd
    write(datstr(4:5),1001)mm
    write(datstr(7:8),1001)yy
1001    format(i2)

!       replace blank with leading 0 in each field if required
    if (datstr(1:1) == ' ')datstr(1:1)='0'
    if (datstr(4:4) == ' ')datstr(4:4)='0'
    if (datstr(7:7) == ' ')datstr(7:7)='0'

!       get the checksum keyword, if it exists, otherwise initialize it
    tstat=status
    call ftgkys(iunit,'CHECKSUM',oldcks,comm,status)
    if (status == 202)then
      status=tstat
      oldcks=' '
      comm='encoded HDU checksum updated on '//datstr
      call ftpkys(iunit,'CHECKSUM','0000000000000000',comm,status)
    end if

!       get the DATASUM keyword and convert it to a double precision value
!       if it exists, otherwise initialize it
    tstat=status
    call ftgkys(iunit,'DATASUM',datsum,comm,status)
    if (status == 202)then
      status=tstat
      odsum=0.
!         set the CHECKSUM keyword as undefined
      oldcks=' '
      comm='data unit checksum updated on '//datstr
      call ftpkys(iunit,'DATASUM','         0',comm,status)
    else
!         decode the datasum into a double precision variable
      do i=1,20
        if (datsum(i:i) /= ' ')then
            call ftc2dd(datsum(i:20),odsum,status)
            if (status == 409)then
!                   couldn't read the keyword; assume it is out of date
                status=tstat
                odsum=-1.
            end if
            go to 15
        end if
      end do
      odsum=0.
    end if

!       rewrite the header END card, and following blank fill
15      call ftwend(iunit,status)
    if (status > 0)return

!       now re-read the required keywords to determine the structure
    call ftrhdu(iunit,dummy,status)

!       write the correct data fill values, if they are not already correct
    call ftpdfl(iunit,status)

!       calc. checksum of the data records; first, calc number of data records
    nrec=(hdstrt(ibuff,chdu(ibuff)+1)-dtstrt(ibuff))/2880
    dsum=0.

    if (nrec > 0)then
!           move to the start of the data
        call ftmbyt(iunit,dtstrt(ibuff),.true.,status)

!           accumulate the 32-bit 1's complement checksum
        call ftcsum(iunit,nrec,dsum,status)
    end if

    if (dsum /= odsum)then
!               modify the DATASUM keyword with the correct value
            comm='data unit checksum updated on '//datstr
!               write the datasum into an I10 integer string
            write(datsum,2000)dsum
2000            format(f11.0)
            call ftmkys(iunit,'DATASUM',datsum(1:10),comm,status)
!               set the CHECKSUM keyword as undefined
            oldcks=' '
    end if

!       if DATASUM was correct, check if CHECKSUM is still OK
    if (oldcks /= ' ')then

!           move to the start of the header
        call ftmbyt(iunit,hdstrt(ibuff,chdu(ibuff)),.true.,status)

!           accumulate the header checksum into the previous data checksum
        nrec= (dtstrt(ibuff)-hdstrt(ibuff,chdu(ibuff)))/2880
        sum=dsum
        call ftcsum(iunit,nrec,sum,status)

!           encode the COMPLEMENT of the checksum into a 16-character string
        complm=.true.
        call ftesum(sum,complm,string)

!           return if the checksum is correct
        if (string == '0000000000000000')then
            return
        else if (oldcks == '0000000000000000')then
!               update the CHECKSUM keyword value with the checksum string
            call ftmkys(iunit,'CHECKSUM',string,'&',status)
            return
        end if
    end if

!       Zero the checksum and compute the new value
    comm='encoded HDU checksum updated on '//datstr
    call ftmkys(iunit,'CHECKSUM','0000000000000000',comm,status)

!       move to the start of the header
    call ftmbyt(iunit,hdstrt(ibuff,chdu(ibuff)),.true.,status)

!       accumulate the header checksum into the previous data checksum
    nrec= (dtstrt(ibuff)-hdstrt(ibuff,chdu(ibuff)))/2880
    sum=dsum
    call ftcsum(iunit,nrec,sum,status)

!       encode the COMPLEMENT of the checksum into a 16-character string
    complm=.true.
    call ftesum(sum,complm,string)

!       update the CHECKSUM keyword value with the checksum string
    call ftmkys(iunit,'CHECKSUM',string,'&',status)
end
subroutine ftpclb(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLB writes unsigned byte data to the specified column of the table.
!
!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   i  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    character array(*)

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    real rval
    logical tofits,lval,trans
    integer*2 i2val
    character sval*30,sform*13,snull*8,i1val*1,messge*80
    double precision i4max,i4min
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
    character chbuff(32000)
    common/ftheap/chbuff
    integer buffer(8000)
    common/fttemp/buffer

!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    call ftgcpr(ounit,colnum,frow,felem,nelem,1, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,ival,snull,status)

    if (status > 0 .or. nelem == 0)return

    i1=1
    ntodo=nelem
    rstart=0
!       the data are being scaled from internal format to FITS:
    tofits=.true.

!       see if we can write the raw input bytes, or whether we have to
!       copy data to temporary array prior to byteswapping or scaling
!       (Note that byteswapping is not a factor for byte data type).
    if (abs(tcode) == 11 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
              trans=.false.
    else
              trans=.true.
    end if

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,repeat-estart,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp + rstart*lenrow + estart*incre
    call ftmbyt(ounit,bstart,.true.,status)

!       copy data to buffer, doing scaling and datatype conversion, if required
    if (tcode == 11)then
!           column data type is B (byte)
        if (trans)then
!               convert the input data into a temporary buffer
            call fti1i1(array(i1),itodo,scale,zero,tofits, &
                ival,i1val,i1val,lval,lval,chbuff,status)
!               do any machine dependent conversion and write the byte data
            call ftpi1b(ounit,itodo,incre,chbuff,status)
        else
!               directly write the input array
            call ftpi1b(ounit,itodo,incre,array(i1),status)
        end if
    else if (tcode == 21)then
!               column data type is I (I*2)
            call fti1i2(array(i1),itodo,scale,zero,tofits, &
               ival,i1val,i2val,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*2 data
            call ftpi2b(ounit,itodo,incre,buffer,status)
    else if (tcode == 41)then
!               column data type is J (I*4)
            call fti1i4(array(i1),itodo,scale,zero,tofits, &
            ival,i1val,ival,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*4 data
            call ftpi4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 42)then
!               column data type is E (R*4)
            call fti1r4(array(i1),itodo,scale,zero,tofits, &
            ival,i1val,rval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*4 data
            call ftpr4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 82)then
!               column data type is D (R*8)
            call fti1r8(array(i1),itodo,scale,zero,tofits, &
            ival,i1val,dval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*8 data
            call ftpr8b(ounit,itodo,incre,buffer,status)
    else
!               this is an ASCII table column
            ival=ichar(array(i1))
            if (ival < 0)ival=ival+256
            dval=(ival-zero)/scale

            if (sform(5:5) == 'I')then
!                   column data type is integer
!                   trap any values that overflow the I*4 range
                if (dval < i4max .and. dval > i4min)then
                    ival=nint(dval)
                else if (ival >= i4max)then
                    status=-11
                    ival=maxi4
                else
                    status=-11
                    ival=mini4
                end if

!                   create the formated character string
                write(sval,sform,err=900)ival
            else
!                   create the formated character string
                write(sval,sform,err=900)dval
            end if

!               write the character string to the FITS file
            call ftpcbf(ounit,twidth,sval,status)
    end if

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error writing elements',i9,' thru',i9, &
           ' of input data array (FTPCLB).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                  estart=0
                  rstart=rstart+1
            end if
            go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while writing FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error writing formatted data value to ASCII table
    write(messge,1002)colnum,rstart+1
1002    format('Error writing column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to write value with format '//sform)
    status=313
end
subroutine ftpclc(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLC writes complex data to the specified column of the table.
!
!       The binary table column being written to must have datatype 'C'
!       and no datatype conversion will be perform if it is not.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   cmp  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
!       the input array is really complex data type
    real array(*)
    integer felemx, nelemx

!       simply multiply the number of elements by 2, and call ftpcle
!       Technically, this is not strictly correct because the data scaling
!       (with TSCALn and TZEROn) is applied differently to complex numbers.
!       In practice, complex number will probably never be scaled so
!       this complication will be ignored.

    felemx = (felem - 1) * 2 + 1
    nelemx  = nelem * 2
    call ftpcle(ounit,colnum,frow,felemx,nelemx,array,status)

end
subroutine ftpcld(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLD writes double precision data to the specified column of the table.
!
!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   d  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    double precision array(*)
    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype,ival
    integer bstart,i1,ntodo,itodo,rstart
    double precision scale,zero,dval
    real rval
    logical tofits,lval,trans
    integer*2 i2val
    character sval*30,sform*13,snull*8,i1val*1,messge*80
    double precision i4max,i4min
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)

    character chbuff(32000)
    common/ftheap/chbuff
    double precision buffer(4000)
    common/fttemp/buffer
    integer compid
    common/ftcpid/compid

!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    call ftgcpr(ounit,colnum,frow,felem,nelem,1, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,ival,snull,status)

    if (status > 0 .or. nelem == 0)return

    i1=1
    ntodo=nelem
    rstart=0
!       the data are being scaled from internal format to FITS:
    tofits=.true.

!       see if we can write the raw input bytes, or whether we have to
!       copy data to temporary array prior to byteswapping or scaling
    if ((compid == 0 .or. compid == -1) .and. &
        abs(tcode) == 82 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
              trans=.false.
    else
              trans=.true.
    end if

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,repeat-estart,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp + rstart*lenrow + estart*incre
    call ftmbyt(ounit,bstart,.true.,status)

!       copy data to buffer, doing scaling and datatype conversion, if required
    if (tcode == 82)then
!           column data type is D (R*8)
        if (trans)then
!               convert the input data into a temporary buffer
            call ftr8r8(array(i1),itodo,scale,zero,tofits, &
                ival,dval,lval,lval,buffer,status)
!               do any machine dependent conversion and write the R*8 data
            call ftpr8b(ounit,itodo,incre,buffer,status)
        else
!               directly write the input array
            call ftpr8b(ounit,itodo,incre,array(i1),status)
        end if
    else if (tcode == 21)then
!               column data type is I (I*2)
            call ftr8i2(array(i1),itodo,scale,zero,tofits, &
            ival,i2val,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*2 data
            call ftpi2b(ounit,itodo,incre,buffer,status)
    else if (tcode == 41)then
!               column data type is J (I*4)
            call ftr8i4(array(i1),itodo,scale,zero,tofits, &
            ival,ival,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*4 data
            call ftpi4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 42)then
!               column data type is E (R*4)
            call ftr8r4(array(i1),itodo,scale,zero,tofits, &
            ival,rval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*4 data
            call ftpr4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 11)then
!               column data type is B (byte)
            call ftr8i1(array(i1),itodo,scale,zero,tofits, &
            ival,i1val,lval,lval,chbuff,status)
!               do any machine dependent data conversion and write the byte data
            call ftpi1b(ounit,itodo,incre,chbuff,status)
    else
!               this is an ASCII table column
            dval=(array(i1)-zero)/scale

            if (sform(5:5) == 'I')then
!                 column data type is integer
!                 trap any values that overflow the I*4 range
              if (dval < i4max .and. dval > i4min)then
                  ival=nint(dval)
              else if (dval >= i4max)then
                  status=-11
                  ival=maxi4
              else
                  status=-11
                  ival=mini4
              end if

!                 create the formated character string
              write(sval,sform,err=900)ival
            else
!                 create the formated character string
              write(sval,sform,err=900)dval
            end if

!               write the character string to the FITS file
            call ftpcbf(ounit,twidth,sval,status)
    end if

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error writing elements',i9,' thru',i9, &
           ' of input data array (FTPCLD).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                  estart=0
                  rstart=rstart+1
            end if
            go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while writing FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error writing formatted data value to ASCII table
    write(messge,1002)colnum,rstart+1
1002    format('Error writing column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to write value with format '//sform)
    status=313
end
subroutine ftpcle(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLE writes real data values to the specified column of the table.
!
!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   r  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    real array(*)
    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    real rval
    logical tofits,lval,trans
    integer*2 i2val
    character sval*30,sform*13,snull*8,i1val*1,messge*80
    double precision i4max,i4min
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)

    character chbuff(32000)
    common/ftheap/chbuff
    real buffer(8000)
    common/fttemp/buffer
    integer compid
    common/ftcpid/compid

!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    call ftgcpr(ounit,colnum,frow,felem,nelem,1, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,ival,snull,status)

    if (status > 0 .or. nelem == 0)return

    i1=1
    ntodo=nelem
    rstart=0
!       the data are being scaled from internal format to FITS:
    tofits=.true.

!       see if we can write the raw input bytes, or whether we have to
!       copy data to temporary array prior to byteswapping or scaling
    if ((compid == 0 .or. compid == -1) .and. &
        abs(tcode) == 42 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
              trans=.false.
    else
              trans=.true.
    end if

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,repeat-estart,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp + rstart*lenrow + estart*incre
    call ftmbyt(ounit,bstart,.true.,status)

!       copy data to buffer, doing scaling and datatype conversion, if required
    if (tcode == 42)then
!           column data type is E (R*4)
        if (trans)then
!               convert the input data into a temporary buffer
            call ftr4r4(array(i1),itodo,scale,zero,tofits, &
                ival,rval,lval,lval,buffer,status)
!               do any machine dependent conversion and write the R*4 data
            call ftpr4b(ounit,itodo,incre,buffer,status)
        else
!               directly write the input array
            call ftpr4b(ounit,itodo,incre,array(i1),status)
        end if
    else if (tcode == 21)then
!               column data type is I (I*2)
            call ftr4i2(array(i1),itodo,scale,zero,tofits, &
            ival,i2val,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*2 data
            call ftpi2b(ounit,itodo,incre,buffer,status)
    else if (tcode == 41)then
!               column data type is J (I*4)
            call ftr4i4(array(i1),itodo,scale,zero,tofits, &
            ival,ival,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*4 data
            call ftpi4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 82)then
!               column data type is D (R*8)
            call ftr4r8(array(i1),itodo,scale,zero,tofits, &
            ival,dval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*8 data
            call ftpr8b(ounit,itodo,incre,buffer,status)
    else if (tcode == 11)then
!               column data type is B (byte)
            call ftr4i1(array(i1),itodo,scale,zero,tofits, &
            ival,i1val,lval,lval,chbuff,status)
!               do any machine dependent data conversion and write the byte data
            call ftpi1b(ounit,itodo,incre,chbuff,status)
    else
!               this is an ASCII table column
            dval=(array(i1)-zero)/scale

            if (sform(5:5) == 'I')then
!                 column data type is integer
!                 trap any values that overflow the I*4 range
              if (dval < i4max .and. dval > i4min)then
                  ival=nint(dval)
              else if (dval >= i4max)then
                  status=-11
                  ival=maxi4
              else
                  status=-11
                  ival=mini4
              end if

!                 create the formated character string
              write(sval,sform,err=900)ival
            else
!                 create the formated character string
              write(sval,sform,err=900)dval
            end if

!               write the character string to the FITS file
            call ftpcbf(ounit,twidth,sval,status)
    end if

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error writing elements',i9,' thru',i9, &
           ' of input data array (FTPCLE).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                  estart=0
                  rstart=rstart+1
            end if
            go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while writing FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error writing formatted data value to ASCII table
    write(messge,1002)colnum,rstart+1
1002    format('Error writing column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to write value with format '//sform)
    status=313
end
subroutine ftpcli(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLI writes integer*2 data values to a column of the table.
!
!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   i*2  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    integer*2 array(*)

    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    real rval
    logical tofits,lval,trans
    integer*2 i2val
    character sval*30,sform*13,snull*8,i1val*1,messge*80

    double precision i4max,i4min
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
    character chbuff(32000)
    common/ftheap/chbuff
    integer*2 buffer(16000)
    common/fttemp/buffer
    integer compid
    common/ftcpid/compid

!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    call ftgcpr(ounit,colnum,frow,felem,nelem,1, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,ival,snull,status)

    if (status > 0 .or. nelem == 0)return

    i1=1
    ntodo=nelem
    rstart=0
!       the data are being scaled from internal format to FITS:
    tofits=.true.

!       see if we can write the raw input bytes, or whether we have to
!       copy data to temporary array prior to byteswapping or scaling
    if ((compid == 0) .and. &
        abs(tcode) == 21 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
              trans=.false.
    else
              trans=.true.
    end if

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,repeat-estart,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp + rstart*lenrow + estart*incre
    call ftmbyt(ounit,bstart,.true.,status)

!       copy data to buffer, doing scaling and datatype conversion, if required
    if (tcode == 21)then
!           column data type is I (I*2)
        if (trans)then
!               convert the input data into a temporary buffer
            call fti2i2(array(i1),itodo,scale,zero,tofits, &
                ival,i2val,i2val,lval,lval,buffer,status)
!               do any machine dependent conversion and write the I*2 data
            call ftpi2b(ounit,itodo,incre,buffer,status)
        else
!               directly write the input array
            call ftpi2b(ounit,itodo,incre,array(i1),status)
        end if
    else if (tcode == 41)then
!               column data type is J (I*4)
            call fti2i4(array(i1),itodo,scale,zero,tofits, &
            ival,i2val,ival,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*4 data
            call ftpi4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 42)then
!               column data type is E (R*4)
            call fti2r4(array(i1),itodo,scale,zero,tofits, &
            ival,i2val,rval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*4 data
            call ftpr4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 82)then
!               column data type is D (R*8)
            call fti2r8(array(i1),itodo,scale,zero,tofits, &
            ival,i2val,dval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*8 data
            call ftpr8b(ounit,itodo,incre,buffer,status)
    else if (tcode == 11)then
!               column data type is B (byte)
            call fti2i1(array(i1),itodo,scale,zero,tofits, &
            ival,i2val,i1val,lval,lval,chbuff,status)
!               do any machine dependent data conversion and write the byte data
            call ftpi1b(ounit,itodo,incre,chbuff,status)
    else
!               this is an ASCII table column
            dval=(array(i1)-zero)/scale

            if (sform(5:5) == 'I')then
!                 column data type is integer
!                 trap any values that overflow the I*4 range
              if (dval < i4max .and. dval > i4min)then
                  ival=nint(dval)
              else if (dval >= i4max)then
                  status=-11
                  ival=maxi4
              else
                  status=-11
                  ival=mini4
              end if

!                 create the formated character string
              write(sval,sform,err=900)ival
            else
!                 create the formated character string
              write(sval,sform,err=900)dval
            end if

!               write the character string to the FITS file
            call ftpcbf(ounit,twidth,sval,status)
    end if

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error writing elements',i9,' thru',i9, &
           ' of input data array (FTPCLI).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                  estart=0
                  rstart=rstart+1
            end if
            go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while writing FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error writing formatted data value to ASCII table
    write(messge,1002)colnum,rstart+1
1002    format('Error writing column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to write value with format '//sform)
    status=313
end
subroutine ftpclj(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLJ writes integer data values to the specified column of the table.
!
!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   i  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    integer array(*)
    integer ibuff,twidth,tcode,maxpix,startp
    integer estart,incre,repeat,lenrow,hdtype
    integer bstart,i1,ntodo,itodo,rstart,ival
    double precision scale,zero,dval
    real rval
    logical tofits,lval,trans
    integer*2 i2val
    character sval*30,sform*13,snull*8,i1val*1,messge*80
    double precision i4max,i4min
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)

    character chbuff(32000)
    common/ftheap/chbuff
    integer buffer(8000)
    common/fttemp/buffer
    integer compid
    common/ftcpid/compid

!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    call ftgcpr(ounit,colnum,frow,felem,nelem,1, &
     ibuff,scale,zero,sform,twidth,tcode,maxpix,startp, &
     estart,incre,repeat,lenrow,hdtype,ival,snull,status)

    if (status > 0 .or. nelem == 0)return

    i1=1
    ntodo=nelem
    rstart=0
!       the data are being scaled from internal format to FITS:
    tofits=.true.

!       see if we can write the raw input bytes, or whether we have to
!       copy data to temporary array prior to byteswapping or scaling
    if ((compid == 0 .or. compid == -1) .and. &
        abs(tcode) == 41 .and. &
        scale == 1.D00 .and. zero == 0.D00)then
              trans=.false.
    else
              trans=.true.
    end if

!       process as many contiguous pixels as possible, up to buffer size
20      itodo=min(ntodo,repeat-estart,maxpix)

!       move the i/o pointer to the start of the sequence of pixels
    bstart=startp + rstart*lenrow + estart*incre
    call ftmbyt(ounit,bstart,.true.,status)

!       copy data to buffer, doing scaling and datatype conversion, if required
    if (tcode == 41)then
!           column data type is J (I*4)
        if (trans)then
!               convert the input data into a temporary buffer
            call fti4i4(array(i1),itodo,scale,zero,tofits, &
                ival,ival,ival,lval,lval,buffer,status)
!               do any machine dependent conversion and write the I*4 data
            call ftpi4b(ounit,itodo,incre,buffer,status)
        else
!               directly write the input array
            call ftpi4b(ounit,itodo,incre,array(i1),status)
        end if
    else if (tcode == 21)then
!               column data type is I (I*2)
            call fti4i2(array(i1),itodo,scale,zero,tofits, &
            ival,ival,i2val,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the I*2 data
            call ftpi2b(ounit,itodo,incre,buffer,status)
    else if (tcode == 42)then
!               column data type is E (R*4)
            call fti4r4(array(i1),itodo,scale,zero,tofits, &
            ival,ival,rval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*4 data
            call ftpr4b(ounit,itodo,incre,buffer,status)
    else if (tcode == 82)then
!               column data type is D (R*8)
            call fti4r8(array(i1),itodo,scale,zero,tofits, &
            ival,ival,dval,lval,lval,buffer,status)
!               do any machine dependent data conversion and write the R*8 data
            call ftpr8b(ounit,itodo,incre,buffer,status)
    else if (tcode == 11)then
!               column data type is B (byte)
            call fti4i1(array(i1),itodo,scale,zero,tofits, &
            ival,ival,i1val,lval,lval,chbuff,status)
!               do any machine dependent data conversion and write the byte data
            call ftpi1b(ounit,itodo,incre,chbuff,status)
    else
!               this is an ASCII table column
            dval=(array(i1)-zero)/scale

            if (sform(5:5) == 'I')then
!                 column data type is integer
!                 trap any values that overflow the I*4 range
              if (dval < i4max .and. dval > i4min)then
                  ival=nint(dval)
              else if (dval >= i4max)then
                  status=-11
                  ival=maxi4
              else
                  status=-11
                  ival=mini4
              end if

!                 create the formated character string
              write(sval,sform,err=900)ival
            else
!                 create the formated character string
              write(sval,sform,err=900)dval
            end if

!               write the character string to the FITS file
            call ftpcbf(ounit,twidth,sval,status)
    end if

    if (status > 0)then
        write(messge,1001)i1,i1+itodo-1
1001        format('Error writing elements',i9,' thru',i9, &
           ' of input data array (FTPCLJ).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+itodo
            estart=estart+itodo
            if (estart == repeat)then
                  estart=0
                  rstart=rstart+1
            end if
            go to 20
    end if

!       check for any overflows
    if (status == -11)then
       status=412
       messge='Numerical overflow during type '// &
              'conversion while writing FITS data.'
       call ftpmsg(messge)
    end if
    return

900     continue
!       error writing formatted data value to ASCII table
    write(messge,1002)colnum,rstart+1
1002    format('Error writing column',i4,', row',i9, &
    ' of the ASCII Table.')
    call ftpmsg(messge)
    call ftpmsg('Tried to write value with format '//sform)
    status=313
end
subroutine ftpcll(ounit,colnum,frow,felem,nelem,lray,status)
!
!*******************************************************************************
!
!! FTPCLL writes logical values to the specified column of the table.
!
!       The binary table column being written to must have datatype 'L'
!       and no datatype conversion will be perform if it is not.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       lray    l  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    logical lray(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character buffer(32000)
    common/ftheap/buffer
!

    integer bstart,maxpix,i
    parameter (maxpix = 32000)
    character messge*80
    integer ibuff,i1,ntodo,itodo,repeat,rstart,estart,tcode
    logical descrp

    if (status > 0)return

    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (nelem == 0)then
!         just return if zero rows to write
      return
    end if

    ibuff=bufnum(ounit)
!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

    if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    end if

    i1=1
    ntodo=nelem
    rstart=frow-1
    estart=felem-1

!       column must be logical data type
    tcode=tdtype(colnum+tstart(ibuff))
    if (tcode == 14)then
            descrp=.false.
            repeat=trept(colnum+tstart(ibuff))
            if (felem > repeat)then
!                  illegal element number
               write(messge,1005)felem
1005               format( &
         'Starting element number is greater than repeat: ',i10)
               call ftpmsg(messge)
               status = 308
               return
            end if
    else if (tcode == -14)then
            descrp=.true.
            repeat=nelem+estart
!               write the number of elements and the starting offset:
            call ftpdes(ounit,colnum,frow,repeat, &
                                heapsz(ibuff),status)
!               move the i/o pointer to the start of the pixel sequence
            bstart=dtstrt(ibuff)+heapsz(ibuff)+ &
                            theap(ibuff)+estart
            call ftmbyt(ounit,bstart,.true.,status)
!               increment the empty heap starting address:
            heapsz(ibuff)=heapsz(ibuff)+repeat
    else
!               error illegal data type code
            status=310
            return
    end if

!       process as many contiguous pixels as possible
20      itodo=min(ntodo,repeat-estart,maxpix)

    if (.not. descrp)then
!           move the i/o pointer to the start of the sequence of pixels
        bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
        tbcol(colnum+tstart(ibuff))+estart
        call ftmbyt(ounit,bstart,.true.,status)
    end if

!       create the buffer of logical bytes
    do i=1,itodo
            if (lray(i1))then
                    buffer(i)='T'
            else
                    buffer(i)='F'
            end if
            i1=i1+1
    end do

!       write out the buffer
    call ftpcbf(ounit,itodo,buffer,status)

    if (status > 0)then
        write(messge,1006)i1,i1+itodo-1
1006        format('Error writing elements',i9,' thru',i9, &
           ' of input data array (FTPCLL).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    if (ntodo > 0)then
!               increment the pointers
            estart=estart+itodo
            if (estart == repeat)then
                    estart=0
                    rstart=rstart+1
            end if
            go to 20
    end if
end
subroutine ftpclm(ounit,colnum,frow,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPCLM writes double precision complex data to a column of the table.
!
!       The binary table column being written to must have datatype 'M'
!       and no datatype conversion will be perform if it is not.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   dcmp  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
!       array is really double precison complex
    double precision array(*)
    integer felemx, nelemx

!       simply multiply the number of elements by 2, and call ftpcld
!       Technically, this is not strictly correct because the data scaling
!       (with TSCALn and TZEROn) is applied differently to complex numbers.
!       In practice, complex number will probably never be scaled so
!       this complication will be ignored.

    felemx = (felem - 1) * 2 + 1
    nelemx  = nelem * 2
    call ftpcld(ounit,colnum,frow,felemx,nelemx,array,status)

end
subroutine ftpcls(ounit,colnum,frow,felem,nelem,sray,status)
!
!*******************************************************************************
!
!! FTPCLS writes character string values to the specified column of the table.
!
!       The binary or ASCII table column being written to must have datatype 'A'

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       sray    c  array of data values to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status
    character ( len = * ) sray(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer bstart,strlen,c1,c2,repeat,twidth
    integer ibuff,i1,ntodo,rstart,estart,nchars,clen,tcode
    character sbuff*80,blank*80,messge*80
    logical small,fill

    if (status > 0)return
    ibuff=bufnum(ounit)

    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (nelem == 0)then
!         just return if zero rows to write
      return
    end if

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

    if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    end if

    blank=' '
    i1=1

!       column must be character string data type
    tcode=tdtype(colnum+tstart(ibuff))
    if (tcode == 16)then
!               for ASCII columns, TNULL actually stores the field width
            twidth=tnull(colnum+tstart(ibuff))
            ntodo=nelem
            rstart=frow-1
            repeat=trept(colnum+tstart(ibuff))
            estart=felem-1
            if (estart >= repeat)then
!                  illegal element number
               write(messge,1005)felem
1005               format( &
         'Starting element number is greater than repeat: ',i10)
               call ftpmsg(messge)
               status = 308
               return
            end if
            bstart=dtstrt(ibuff)+rstart*rowlen(ibuff) &
                   +tbcol(colnum+tstart(ibuff))+estart*twidth
    else if (tcode == -16)then
!               this is a variable length descriptor field
!               the length of the output string is defined by nelem
            twidth=nelem
            ntodo=1
            repeat=1
!               write the number of string length and the starting offset:
            call ftpdes(ounit,colnum,frow,twidth, &
                                heapsz(ibuff),status)
!               calc the i/o pointer position for the start of the string
            bstart=dtstrt(ibuff)+heapsz(ibuff)+theap(ibuff)
!               increment the empty heap starting address:
            heapsz(ibuff)=heapsz(ibuff)+twidth
    else
!               error: not a character string column
            status=309
            return
    end if

!       move the i/o pointer to the start of the sequence of pixels
    call ftmbyt(ounit,bstart,.true.,status)

!       is the input string short enough to completely fit in buffer?
    strlen=len(sray(1))
    if (strlen > 80 .and. twidth > 80)then
            small=.false.
    else
            small=.true.
    end if

!       do we need to pad the FITS string field with trailing blanks?
    if (twidth > strlen)then
            fill=.true.
    else
            fill=.false.
    end if

!       process one string at a time
20      continue
    nchars=min(strlen,twidth)
    if (small)then
!               the whole input string fits in the temporary buffer
            sbuff=sray(i1)
!               output the string
            call ftpcbf(ounit,nchars,sbuff,status)
    else
!               have to write the string in several pieces
            c1=1
            c2=80
30              sbuff=sray(i1)(c1:c2)
!               output the string
            clen=c2-c1+1
            call ftpcbf(ounit,clen,sbuff,status)
            nchars=nchars-clen
            if (nchars > 0)then
                    c1=c1+80
                    c2=min(c2+80,c1+nchars-1)
                    go to 30
            end if
    end if

!       pad any remaining space in the column with blanks
    if (fill)then
            nchars=twidth-strlen
40              clen=min(nchars,80)
            call ftpcbf(ounit,clen,blank,status)
            nchars=nchars-80
            if (nchars > 0)go to 40
    end if

    if (status > 0)then
        write(messge,1006)i1
1006        format('Error writing element',i9, &
           ' of input string array (FTPCLS).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-1
    if (ntodo > 0)then
!               increment the pointers
            i1=i1+1
            estart=estart+1
            if (estart == repeat)then
                    estart=0
                    rstart=rstart+1
                    bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
                    tbcol(colnum+tstart(ibuff))
!                       move the i/o pointer
                    call ftmbyt(ounit,bstart,.true.,status)
            end if
            go to 20
    end if
end
subroutine ftpclu(ounit,colnum,frow,felem,nelem,status)
!
!*******************************************************************************
!
!! FTPCLU sets elements of a table to be undefined.
!
!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,frow,felem,nelem,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
    character snull*500
    character xdummy(31500)
    common/ftheap/snull,xdummy
!

    integer bytpix,bstart,i4null(2),tcode,nchars,i,offset,nulval
    integer ibuff,ntodo,itodo,repeat,rstart,estart
    integer*2 i2null,i1
    integer rnull(2)
    logical descrp
    character i1null
    character messge*80

    if (status > 0)return

    if (frow < 1)then
      write(messge,1001)frow
1001      format('Starting row number is out of range: ',i10)
      call ftpmsg(messge)
      status = 307
      return
    else if (felem < 1)then
      write(messge,1002)felem
1002      format('Starting element number is out of range: ',i10)
      call ftpmsg(messge)
      status = 308
      return
    else if (nelem < 0)then
      write(messge,1003)nelem
1003      format('Negative no. of elements to read or write: ',i10)
      call ftpmsg(messge)
      status = 306
      return
    else if (nelem == 0)then
      return
    end if

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

    if (colnum < 1 .or. colnum > tfield(ibuff))then
      write(messge,1004)colnum
1004      format('Specified column number is out of range: ',i10)
      call ftpmsg(messge)
      status = 302
      return
    end if

    tcode=tdtype(colnum+tstart(ibuff))
    bytpix=max(abs(tcode)/10,1)

    descrp=.false.
    ntodo=nelem
    rstart=frow-1
    estart=felem-1
    i1=1

    if (tcode == 16)then
!               this is an ASCII field
            repeat=trept(colnum+tstart(ibuff))
            if (felem > repeat)then
!                  illegal element number
               write(messge,1005)felem
1005               format( &
         'Starting element number is greater than repeat: ',i10)
               call ftpmsg(messge)
               status = 308
               return
            end if
            if (cnull(colnum+tstart(ibuff))(1:1) == char(1))then
!                       error: null value has not been defined
                    status=314
            call ftpmsg('Null value string for ASCII table'// &
            ' column has not yet been defined (FTPCLU).')
                    return
            end if
!               the TNULL parameter stores the width of the character field
            bytpix=tnull(colnum+tstart(ibuff))
    else
!               this is a binary table
            nulval=tnull(colnum+tstart(ibuff))

            if (tcode > 0)then
                    if (hdutyp(ibuff) == 0)then
!                           if this is a primary array or image extension, then
!                           set repeat as large as needed to write all
!                           the pixels.  This prevents an error message if
!                           array size is not yet known.  The actual array
!                           dimension must be defined by the NAXISn keywords
!                           before closing this HDU.
                        repeat=estart+nelem
                    else
                        repeat=trept(colnum+tstart(ibuff))
                    end if

                    if (felem > repeat)then
!                           illegal element number
                        write(messge,1004)felem
                        call ftpmsg(messge)
                        status = 308
                        return
                    end if


            else
!                       this is a variable length descriptor column
                    descrp=.true.
                    tcode=-tcode
!                       read the number of elements and the starting offset:
                    call ftgdes(ounit,colnum,frow,repeat, &
                                offset,status)
                    if (ntodo+estart > repeat)then
!                               error:  tried to write past end of record
                            status=319
                            return
                    end if

!                       move the i/o pointer to the start of the pixel sequence
                    bstart=dtstrt(ibuff)+offset+ &
                            theap(ibuff)+estart*bytpix
                    call ftmbyt(ounit,bstart,.true.,status)
            end if

            if (tcode==11 .or. tcode==21 .or. tcode==41)then
                    if (nulval == 123454321)then
!                               error: null value has not been defined
                            status=314
            call ftpmsg('Null value for integer'// &
            ' column has not yet been defined (FTPCLU).')
                            return
                    end if
            else
!                       set the floating point Not-a-Number values
                    do 10 i=1,2
                      rnull(i) = -1
10                      continue
            end if

    end if

!       process as many contiguous pixels as possible
20      itodo=min(ntodo,repeat-estart)

    if (.not. descrp)then
!           move the i/o pointer to the start of the sequence of pixels
        bstart=dtstrt(ibuff)+rstart*rowlen(ibuff) &
               +tbcol(colnum+tstart(ibuff))+estart*bytpix
        call ftmbyt(ounit,bstart,.true.,status)
    end if

!       write the appropriate null value to the pixels
    if (tcode == 21)then
!               column data type is I (I*2)
            do 5 i=1,itodo
                    i2null=nulval
                    call ftpi2b(ounit,1,0,i2null,status)
5               continue
    else if (tcode == 41)then
!               column data type is J (I*4)
            do 15 i=1,itodo
                    i4null(1)=nulval
                    call ftpi4b(ounit,1,0,i4null,status)
15              continue
    else if (tcode == 42)then
!               column data type is E (R*4)
            do 25 i=1,itodo
                    call ftpbyt(ounit,4,rnull,status)
25              continue
    else if (tcode == 82 .or. tcode == 83)then
!               column data type is D (R*8), or C complex 2 x R*4
            do 35 i=1,itodo
                    call ftpbyt(ounit,8,rnull,status)
35              continue
    else if (tcode == 16)then
!               this is an ASCII table column
            snull=cnull(colnum+tstart(ibuff))
!               write up to 500 characters in the column, remainder unchanged
!               (500 is the maximum size string allowed in IBM AIX compiler)
            nchars=min(bytpix,500)
            do 45 i=1,itodo
                    call ftpcbf(ounit,nchars,snull,status)
45              continue
    else if (tcode == 11)then
!               column data type is B (byte)
            i1null=char(nulval)
            do 55 i=1,itodo
                    call ftpcbf(ounit,1,i1null,status)
55              continue
    else if (tcode == 163)then
!               column data type is double complex (M)
            do 65 i=1,itodo*2
                    call ftpbyt(ounit,8,rnull,status)
65              continue
    else if (tcode == 14)then
!               column data type is logical (L)
            i4null(1)=0
            do 85 i=1,itodo
                    call ftpbyt(ounit,1,i4null,status)
85              continue
    end if


    if (status > 0)then
        write(messge,1006)i1,i1+itodo-1
1006        format('Error writing NULL elements',i9,' thru',i9, &
           ' (FTPCLU).')
        call ftpmsg(messge)
        return
    end if

!       find number of pixels left to do, and quit if none left
    ntodo=ntodo-itodo
    i1 = i1 + itodo
    if (ntodo > 0)then
!               increment the pointers
            estart=estart+itodo
            if (estart == repeat)then
                    estart=0
                    rstart=rstart+1
            end if
            go to 20
    end if
end
subroutine ftpclx(iunit,colnum,frow,fbit,nbit,lray,status)
!
!*******************************************************************************
!
!! FTPCLX writes logical values to a bit or byte column of the binary table.
!
!  If the LRAY parameter is .true.,
!       then the corresponding bit is set to 1, otherwise it is set
!       to 0.
!       The binary table column being written to must have datatype 'B'
!       or 'X'.

!       iunit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       fbit    i  first bit within the row to write
!       nbit    i  number of bits to write
!       lray    l  array of logical data values corresponding to the bits
!                        to be written
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Mar 1992
!       modified by Wm Pence May 1992 to remove call to system dependent
!                                     bit testing and setting routines.

    integer iunit,colnum,frow,fbit,nbit,status
    logical lray(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer bstart,offset,tcode,fbyte,bitloc,ndone,tstat
    integer ibuff,i,ntodo,repeat,rstart,estart,buffer
    logical descrp,wrbit(8),setbit(8)
    character cbuff
    character crow*9

    if (status > 0)return

    ibuff=bufnum(iunit)
    tcode=tdtype(colnum+tstart(ibuff))

!       check input parameters
    if (nbit <= 0)then
            return
    else if (frow < 1)then
!               error: illegal first row number
            status=307
            write(crow,2000)frow
2000            format(i9)
            call ftpmsg('Starting row number for table write '// &
            'request is out of range:'//crow//' (FTPCLX).')
            return
    else if (fbit < 1)then
!               illegal element number
            status=308
            write(crow,2000)fbit
            call ftpmsg('Starting element number for write '// &
            'request is out of range:'//crow//' (FTPCLX).')
            return
    end if

    fbyte=(fbit+7)/8
    bitloc=fbit-(fbit-1)/8*8
    ndone=0
    ntodo=nbit
    rstart=frow-1
    estart=fbyte-1

    if (tcode == 11)then
            descrp=.false.
!               N.B: REPEAT is the number of bytes, not number of bits
            repeat=trept(colnum+tstart(ibuff))
            if (fbyte > repeat)then
!                               illegal element number
                            status=308
                            write(crow,2000)fbit
                call ftpmsg('Starting element number for write '// &
                'request is out of range:'//crow//' (FTPCLX).')
                            return
            end if
!               calc the i/o pointer location to start of sequence of pixels
            bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
            tbcol(colnum+tstart(ibuff))+estart
    else if (tcode == -11)then
!               this is a variable length descriptor column
            descrp=.true.
!               only bit arrays (tform = 'X') are supported for variable
!               length arrays.  REPEAT is the number of BITS in the array.
            repeat=fbit+nbit-1
            offset=heapsz(ibuff)
!               write the number of elements and the starting offset:
            call ftpdes(iunit,colnum,frow,repeat, &
                                offset,status)
!               calc the i/o pointer location to start of sequence of pixels
            bstart=dtstrt(ibuff)+offset+ &
                            theap(ibuff)+estart
!               increment the empty heap starting address (in bytes):
            repeat=(repeat+7)/8
            heapsz(ibuff)=heapsz(ibuff)+repeat
    else
!               column must be byte or bit data type
            status=310
            return
    end if

!       move the i/o pointer to the start of the pixel sequence
    call ftmbyt(iunit,bstart,.true.,status)
    tstat=0

!       read the next byte (we may only be modifying some of the bits)
20      call ftgcbf(iunit,1,cbuff,status)
    if (status == 107)then
!            hit end of file trying to read the byte, so just set byte = 0
         status=tstat
         cbuff=char(0)
    end if

    buffer=ichar(cbuff)
    if (buffer < 0)buffer=buffer+256
!       move back, to be able to overwrite the byte
    call ftmbyt(iunit,bstart,.true.,status)

!       reset flags indicating which bits are to be set
    wrbit(1)=.false.
    wrbit(2)=.false.
    wrbit(3)=.false.
    wrbit(4)=.false.
    wrbit(5)=.false.
    wrbit(6)=.false.
    wrbit(7)=.false.
    wrbit(8)=.false.

!       flag the bits that are to be set
    do 10 i=bitloc,8
            wrbit(i)=.true.
            ndone=ndone+1
            if(lray(ndone))then
                    setbit(i)=.true.
            else
                    setbit(i)=.false.
            end if
            if (ndone == ntodo)go to 100
10      continue

!       set or reset the bits within the byte
    call ftpbit(setbit,wrbit,buffer)

!       write the new byte
    cbuff=char(buffer)
    call ftpcbf(iunit,1,cbuff,status)

!       not done, so get the next byte
    bstart=bstart+1
    if (.not. descrp)then
            estart=estart+1
            if (estart == repeat)then
!                       move the i/o pointer to the next row of pixels
                    estart=0
                    rstart=rstart+1
                    bstart=dtstrt(ibuff)+rstart*rowlen(ibuff)+ &
                           tbcol(colnum+tstart(ibuff))+estart
                    call ftmbyt(iunit,bstart,.true.,status)
            end if
    end if
    bitloc=1
    go to 20

100     continue
!       set or reset the bits within the byte
    call ftpbit(setbit,wrbit,buffer)

!       write the new byte
    cbuff=char(buffer)
    call ftpcbf(iunit,1,cbuff,status)
end
subroutine ftpcnb(ounit,colnum,frow,felem,nelem,array,nulval, &
                      status)
!
!*******************************************************************************
!
!! FTPCNB writes character (byte) pixels to the specified column of a table.
!
!  Any input pixels equal to the value of NULVAL will
!       be replaced by the appropriate null value in the output FITS file.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   c*1  array of data values to be written
!       nulval  c*1  pixel value used to represent an undefine pixel
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,colnum,frow,felem,nelem,status
    character array(*),nulval

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,repeat,first,ngood,nbad,i,fstelm,fstrow

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

!       get the column repeat count and calculate the absolute position within
!       the column of the first element to be written
    repeat=trept(colnum+tstart(ibuff))
    first=(frow-1)*repeat+felem-1

    ngood=0
    nbad=0
    do 10 i=1,nelem
        if (array(i) /= nulval)then
            ngood=ngood+1
            if (nbad > 0)then
!                   write the previous consecutive set of null pixels
                fstelm=i-nbad+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
                nbad=0
            end if
        else
            nbad=nbad+1
            if (ngood > 0)then
!                   write the previous consecutive set of good pixels
                fstelm=i-ngood+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclb(ounit,colnum,fstrow,fstelm,ngood, &
                            array(i-ngood),status)
                ngood=0
            end if
        end if
10      continue

!       finished;  now just write the last set of pixels
    if (nbad > 0)then
!           write the consecutive set of null pixels
        fstelm=i-nbad+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
    else
!           write the consecutive set of good pixels
        fstelm=i-ngood+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclb(ounit,colnum,fstrow,fstelm,ngood, &
                    array(i-ngood),status)
    end if
end
subroutine ftpcnd(ounit,colnum,frow,felem,nelem,array,nulval, &
                      status)
!
!*******************************************************************************
!
!! FTPCND writes double precision pixels to the specified column of a table.
!
!       Any input pixels equal to the value of NULVAL will
!       be replaced by the appropriate null value in the output FITS file.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   d  array of data values to be written
!       nulval  d  pixel value used to represent an undefine pixel
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,colnum,frow,felem,nelem,status
    double precision array(*),nulval

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,repeat,first,ngood,nbad,i,fstelm,fstrow

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

!       get the column repeat count and calculate the absolute position within
!       the column of the first element to be written
    repeat=trept(colnum+tstart(ibuff))
    first=(frow-1)*repeat+felem-1

    ngood=0
    nbad=0
    do 10 i=1,nelem
        if (array(i) /= nulval)then
            ngood=ngood+1
            if (nbad > 0)then
!                   write the previous consecutive set of null pixels
                fstelm=i-nbad+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
                nbad=0
            end if
        else
            nbad=nbad+1
            if (ngood > 0)then
!                   write the previous consecutive set of good pixels
                fstelm=i-ngood+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpcld(ounit,colnum,fstrow,fstelm,ngood, &
                            array(i-ngood),status)
                ngood=0
            end if
        end if
10      continue

!       finished;  now just write the last set of pixels
    if (nbad > 0)then
!           write the consecutive set of null pixels
        fstelm=i-nbad+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
    else
!           write the consecutive set of good pixels
        fstelm=i-ngood+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpcld(ounit,colnum,fstrow,fstelm,ngood, &
                    array(i-ngood),status)
    end if
end
subroutine ftpcne(ounit,colnum,frow,felem,nelem,array,nulval, &
                      status)
!
!*******************************************************************************
!
!! FTPCNE writes floating point pixels to the specified column of a table.
!
!  Any input pixels equal to the value of NULVAL will
!       be replaced by the appropriate null value in the output FITS file.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   r  array of data values to be written
!       nulval  r  pixel value used to represent an undefine pixel
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,colnum,frow,felem,nelem,status
    real array(*),nulval

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,repeat,first,ngood,nbad,i,fstelm,fstrow

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

!       get the column repeat count and calculate the absolute position within
!       the column of the first element to be written
    repeat=trept(colnum+tstart(ibuff))
    first=(frow-1)*repeat+felem-1

    ngood=0
    nbad=0
    do 10 i=1,nelem
        if (array(i) /= nulval)then
            ngood=ngood+1
            if (nbad > 0)then
!                   write the previous consecutive set of null pixels
                fstelm=i-nbad+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
                nbad=0
            end if
        else
            nbad=nbad+1
            if (ngood > 0)then
!                   write the previous consecutive set of good pixels
                fstelm=i-ngood+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpcle(ounit,colnum,fstrow,fstelm,ngood, &
                            array(i-ngood),status)
                ngood=0
            end if
        end if
10      continue

!       finished;  now just write the last set of pixels
    if (nbad > 0)then
!           write the consecutive set of null pixels
        fstelm=i-nbad+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
    else
!           write the consecutive set of good pixels
        fstelm=i-ngood+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpcle(ounit,colnum,fstrow,fstelm,ngood, &
                    array(i-ngood),status)
    end if
end
subroutine ftpcni(ounit,colnum,frow,felem,nelem,array,nulval, &
                      status)
!
!*******************************************************************************
!
!! FTPCNI writes integer*2 pixels to the specified column of a table.
!
!  Any input pixels equal to the value of NULVAL will
!       be replaced by the appropriate null value in the output FITS file.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   i*2  array of data values to be written
!       nulval  i*2  pixel value used to represent an undefine pixel
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,colnum,frow,felem,nelem,status
    integer*2 array(*),nulval

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,repeat,first,ngood,nbad,i,fstelm,fstrow

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

!       get the column repeat count and calculate the absolute position within
!       the column of the first element to be written
    repeat=trept(colnum+tstart(ibuff))
    first=(frow-1)*repeat+felem-1

    ngood=0
    nbad=0
    do 10 i=1,nelem
        if (array(i) /= nulval)then
            ngood=ngood+1
            if (nbad > 0)then
!                   write the previous consecutive set of null pixels
                fstelm=i-nbad+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
                nbad=0
            end if
        else
            nbad=nbad+1
            if (ngood > 0)then
!                   write the previous consecutive set of good pixels
                fstelm=i-ngood+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpcli(ounit,colnum,fstrow,fstelm,ngood, &
                            array(i-ngood),status)
                ngood=0
            end if
        end if
10      continue

!       finished;  now just write the last set of pixels
    if (nbad > 0)then
!           write the consecutive set of null pixels
        fstelm=i-nbad+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
    else
!           write the consecutive set of good pixels
        fstelm=i-ngood+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpcli(ounit,colnum,fstrow,fstelm,ngood, &
                    array(i-ngood),status)
    end if
end
subroutine ftpcnj(ounit,colnum,frow,felem,nelem,array,nulval, &
                      status)
!
!*******************************************************************************
!
!! FTPCNJ writes integer pixels to the specified column of a table.
!
!  Any input pixels equal to the value of NULVAL will
!       be replaced by the appropriate null value in the output FITS file.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       frow    i  first row to write
!       felem   i  first element within the row to write
!       nelem   i  number of elements to write
!       array   i  array of data values to be written
!       nulval  i  pixel value used to represent an undefine pixel
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,colnum,frow,felem,nelem,status
    integer array(*),nulval

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,repeat,first,ngood,nbad,i,fstelm,fstrow

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)

!       get the column repeat count and calculate the absolute position within
!       the column of the first element to be written
    repeat=trept(colnum+tstart(ibuff))
    first=(frow-1)*repeat+felem-1

    ngood=0
    nbad=0
    do 10 i=1,nelem
        if (array(i) /= nulval)then
            ngood=ngood+1
            if (nbad > 0)then
!                   write the previous consecutive set of null pixels
                fstelm=i-nbad+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
                nbad=0
            end if
        else
            nbad=nbad+1
            if (ngood > 0)then
!                   write the previous consecutive set of good pixels
                fstelm=i-ngood+first
!                   calculate the row and element of the first pixel to write
                fstrow=(fstelm-1)/repeat+1
                fstelm=fstelm-(fstrow-1)*repeat
                call ftpclj(ounit,colnum,fstrow,fstelm,ngood, &
                            array(i-ngood),status)
                ngood=0
            end if
        end if
10      continue

!       finished;  now just write the last set of pixels
    if (nbad > 0)then
!           write the consecutive set of null pixels
        fstelm=i-nbad+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclu(ounit,colnum,fstrow,fstelm,nbad,status)
    else
!           write the consecutive set of good pixels
        fstelm=i-ngood+first
        fstrow=(fstelm-1)/repeat+1
        fstelm=fstelm-(fstrow-1)*repeat
        call ftpclj(ounit,colnum,fstrow,fstelm,ngood, &
                    array(i-ngood),status)
    end if
end
subroutine ftpcom(ounit,commnt,status)
!
!*******************************************************************************
!
!! FTPCOM writes a COMMENT record to the FITS header.
!
!       ounit   i  fortran output unit number
!       commnt c  input comment string
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,status,strlen,actlen,i,nkeys,c1,c2
    character ( len = * ) commnt
    character*80  rec

    if (status > 0)return

!       find the length of the string, and write it out 70 characters at a time
    nkeys=1
    strlen=len(commnt)
    actlen=strlen
    do i=strlen,1,-1
            if (commnt(i:i) /= ' ')then
                    actlen=i
                    exit
            end if
    end do

    c1=1
    c2=min(actlen,70)
    nkeys=(actlen-1)/70+1
    do i=1,nkeys
            rec='COMMENT   '//commnt(c1:c2)
            call ftprec(ounit,rec,status)
            c1=c1+70
            c2=min(actlen,c2+70)
    end do

end
subroutine ftpdat(ounit,status)
!
!*******************************************************************************
!
!! FTPDAT writes the current date to the DATE keyword in the ounit CHU.
!
!       ounit   i  fortran output unit number
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Jan 1992

    integer ounit,status,dd,mm,yy
    character datstr*8

!       call the system dependent routine to get the current date
    call ftgsdt(dd,mm,yy,status)
    if (status > 0)return

    datstr='  /  /  '
    write(datstr(1:2),1001)dd
    write(datstr(4:5),1001)mm
    write(datstr(7:8),1001)yy
1001    format(i2)

!       replace blank with leading 0 in each field if required
    if (datstr(1:1) == ' ')datstr(1:1)='0'
    if (datstr(4:4) == ' ')datstr(4:4)='0'
    if (datstr(7:7) == ' ')datstr(7:7)='0'

!       update the DATE keyword
    call ftukys(ounit,'DATE',datstr, &
               'FITS file creation date (dd/mm/yy)',status)
end
subroutine ftpdef(ounit,bitpix,naxis,naxes,pcount,gcount, &
                      status)
!
!*******************************************************************************
!
!! FTPDEF defines the structure of the primary data unit.
!
!
!       Primary data DEFinition
!       define the structure of the primary data unit or an IMAGE extension
!
!       ounit   i  Fortran I/O unit number
!       bitpix  i  bits per pixel value
!       naxis   i  number of data axes
!       naxes   i  length of each data axis (array)
!       pcount  i  number of group parameters
!       gcount  i  number of 'random groups'
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,bitpix,naxis,naxes(*),pcount,gcount,status

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,ttype,bytlen,npix,i,pcnt,gcnt
    character caxis*20

    if (status > 0)return

    ibuff=bufnum(ounit)

    if (dtstrt(ibuff) < 0)then
!               freeze the header at its current size
            call fthdef(ounit,0,status)
            if (status > 0)return
    end if

!       check for error conditions
    if (naxis < 0)then
            status=212
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPDEF ' &
            //'is illegal.')

    else if (pcount < 0)then
            status=214
    else if (gcount < 0)then
            status=215
    else
            go to 5
    end if
    return

!       test that bitpix has a legal value and set the datatype code value
5       if (bitpix == 8)then
            ttype=11
            bytlen=1
    else if (bitpix == 16)then
            ttype=21
            bytlen=2
    else if (bitpix == 32)then
            ttype=41
            bytlen=4
    else if (bitpix == -32)then
            ttype=42
            bytlen=4
    else if (bitpix == -64)then
            ttype=82
            bytlen=8
    else
!               illegal value of bitpix
            status=211
            return
    end if

!       calculate the number of pixels in the array
    if (naxis == 0)then
!               no data
            npix=0
            gcnt=0
            pcnt=0
    else
!               make sure that the gcount is not zero
            gcnt=max(gcount,1)
            pcnt=pcount
            npix=1
            do 10 i=1,naxis
                    if (naxes(i) >= 0)then
!       The convention used by 'random groups' with NAXIS1 = 0 is not
!       directly supported here.  If one wants to write a 'random group'
!       FITS file, then one should call FTPDEF with naxes(1) = 1, but
!       then write the required header keywords (with FTPHPR) with
!       naxes(1) = 0.
                            npix=npix*naxes(i)
                    else if (naxes(i) < 0)then
                            status=213
                            return
                    end if
10              continue
    end if
!       the next HDU begins in the next logical block after the data
    hdstrt(ibuff,chdu(ibuff)+1)= &
            dtstrt(ibuff)+((pcnt+npix)*bytlen*gcnt+2879)/2880*2880

!       the primary array is actually interpreted as a binary table.  There
!       are two columns: the first column contains the
!       group parameters, if any, and the second column contains the
!       primary array of data.  Each group is a separate row in the table.
!       The scaling and null values are set to the default values.

    hdutyp(ibuff)=0
    tfield(ibuff)=2

    if (nxtfld + 2 > nf)then
!               too many columns open at one time; exceeded array dimensions
            status=111
    else
            tstart(ibuff)=nxtfld
            nxtfld=nxtfld+2
            tdtype(1+tstart(ibuff))=ttype
            tdtype(2+tstart(ibuff))=ttype
            trept(1+tstart(ibuff))=pcnt
            trept(2+tstart(ibuff))=npix
!               choose a special value to represent the absence of a blank value
            tnull(1+tstart(ibuff))=123454321
            tnull(2+tstart(ibuff))=123454321
            tscale(1+tstart(ibuff))=1.
            tscale(2+tstart(ibuff))=1.
            tzero(1+tstart(ibuff))=0.
            tzero(2+tstart(ibuff))=0.
            tbcol(1+tstart(ibuff))=0
            tbcol(2+tstart(ibuff))=pcnt*bytlen
            rowlen(ibuff)=(pcnt+npix)*bytlen
    end if

!       initialize the fictitious heap starting address (immediately following
!       the array data) and a zero length heap.  This is used to find the
!   end of the data when checking the fill values in the last block.
    heapsz(ibuff)=0
    theap(ibuff)=(pcnt+npix)*bytlen*gcnt
end
subroutine ftpdes(ounit,colnum,rownum,nelem,offset,status)
!
!*******************************************************************************
!
!! FTPDES writes the descriptor values to a binary table.
!
!  This is only
!       used for column which have TFORMn = 'P', i.e., for variable
!       length arrays.

!       ounit   i  fortran unit number
!       colnum  i  number of the column to write to
!       rownum  i  number of the row to write
!       nelem   i  input number of elements
!       offset  i  input byte offset of the first element
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Nov 1991

    integer ounit,colnum,rownum,nelem,offset,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bstart,iray(2)

    if (status > 0)return
    if (rownum < 1)then
!               error: illegal row number
            status=307
            return
    end if

    ibuff=bufnum(ounit)

!       check that this is really a 'P' type column
    if (tdtype(colnum+tstart(ibuff)) >= 0)then
            status=317
            return
    end if

!       move to the specified column and row:
    bstart=dtstrt(ibuff)+(rownum-1)*rowlen(ibuff) &
           +tbcol(colnum+tstart(ibuff))
    call ftmbyt(ounit,bstart,.true.,status)

!       now write the number of elements and the offset to the table:
    iray(1)=nelem
    iray(2)=offset
    call ftpi4b(ounit,2,0,iray,status)
end
subroutine ftpdfl(iunit,status)
!
!*******************************************************************************
!
!! FTPDFL writes the Data Unit Fill values if they are not already correct.
!
!       Fill the data unit with zeros or blanks depending on the type of HDU
!       from the end of the data to the end of the current FITS 2880 byte block

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1994

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nf = 3000)
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character *2880 chbuff
    character chfill,xdummy(29119)
    common/ftheap/chbuff,chfill,xdummy
!

    integer ibuff,filpos,nfill,i,tstat

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check if the data unit is null

    if (theap(ibuff) == 0)return

    filpos=dtstrt(ibuff)+theap(ibuff)+heapsz(ibuff)
    nfill=(filpos+2879)/2880*2880-filpos

!       return if there are no fill bytes
    if (nfill == 0)return

!       set the correct fill value to be checked
    if (hdutyp(ibuff) == 1)then
!              this is an ASCII table; should be filled with blanks
           chfill=char(32)
    else
           chfill=char(0)
    end if

!       move to the beginning of the fill bytes and read them
    tstat=status
    call ftmbyt(iunit,filpos,.true.,status)
    call ftgcbf(iunit,nfill,chbuff,status)

    if (status > 0)then
!           fill bytes probably haven't been written yet so have to write them
        status=tstat
        go to 100
    end if

!       check if all the fill values are correct
    do 10 i=1,nfill
        if (chbuff(i:i) /= chfill)go to 100
10      continue

!       fill bytes were correct, so just return
    return

100     continue

!       fill the buffer with the correct fill value
    do 20 i=1,nfill
           chbuff(i:i)=chfill
20      continue

!       move to the beginning of the fill bytes
    call ftmbyt(iunit,filpos,.true.,status)

!       write all the fill bytes
    call ftpcbf(iunit,nfill,chbuff,status)

    if (status > 0)then
       call ftpmsg('Error writing Data Unit fill bytes (FTPDFL).')
    end if
end
subroutine ftpgpb(ounit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTPGPB writes an array of group parmeters into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter to be written (starting with 1)
!       nparm   i  number of group parameters to be written
!       array   b  the array of group parameters to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,fparm,nparm,status,row

    character array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpclb(ounit,1,row,fparm,nparm,array,status)
end
subroutine ftpgpd(ounit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTPGPD writes an array of group parameters into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter to be written (starting with 1)
!       nparm   i  number of group parameters to be written
!       array   d  the array of group parameters to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,fparm,nparm,status,row
    double precision array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcld(ounit,1,row,fparm,nparm,array,status)
end
subroutine ftpgpe(ounit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTPGPE writes an array of group parmeters into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter to be written (starting with 1)
!       nparm   i  number of group parameters to be written
!       array   r  the array of group parameters to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,fparm,nparm,status,row
    real array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcle(ounit,1,row,fparm,nparm,array,status)
end
subroutine ftpgpi(ounit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTPGPI writes an array of group parmeters into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter to be written (starting with 1)
!       nparm   i  number of group parameters to be written
!       array   i*2  the array of group parameters to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,fparm,nparm,status,row
    integer*2 array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcli(ounit,1,row,fparm,nparm,array,status)
end
subroutine ftpgpj(ounit,group,fparm,nparm,array,status)
!
!*******************************************************************************
!
!! FTPGPJ writes an array of group parmeters into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       fparm   i  the first group parameter to be written (starting with 1)
!       nparm   i  number of group parameters to be written
!       array   i  the array of group parameters to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,fparm,nparm,status,row
    integer array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpclj(ounit,1,row,fparm,nparm,array,status)
end
subroutine ftphbn(ounit,nrows,nfield,ttype,tform,tunit, &
                      extnam,pcount,status)
!
!*******************************************************************************
!
!! FTPHBN writes required standard header keywords for a binary table extension.
!
!       ounit   i  fortran output unit number
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       ttype   c  name of each field (array) (optional)
!       tform   c  format of each field (array)
!       tunit   c  units of each field (array) (optional)
!       extnam  c  name of table extension (optional)
!       pcount  i  size of special data area following the table (usually = 0)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,nrows,nfield,pcount,status
    integer i,lenrow,dtype,rcount,xbcol,length,width
    integer nkeys,nmore
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam
    character comm*48,tfm*40

    if (status > 0)return

    call ftghsp(ounit,nkeys,nmore,status)
    if (nkeys /= 0)then
!            some keywords have already been written
         status=201
         return
    end if

    comm='binary table extension'
    call ftpkys(ounit,'XTENSION','BINTABLE',comm,status)

    comm='8-bit bytes'
    call ftpkyj(ounit,'BITPIX',8,comm,status)

    comm='2-dimensional binary table'
    call ftpkyj(ounit,'NAXIS',2,comm,status)

    if (status > 0)return

!       calculate the total width of each row, in bytes
    lenrow=0
    do 10 i=1,nfield
!               get the numerical datatype and repeat count of the field
            call ftbnfm(tform(i),dtype,rcount,width,status)
            if (dtype == 1)then
!                       treat Bit datatype as if it were a Byte datatype
                    dtype=11
                    rcount=(rcount+7)/8
            end if
!               get the width of the field
            call ftgtbc(1,dtype,rcount,xbcol,length,status)
            lenrow=lenrow+length
10      continue

    comm='width of table in bytes'
    call ftpkyj(ounit,'NAXIS1',lenrow,comm,status)

    if (status > 0)return

    if (nrows >= 0)then
            comm='number of rows in table'
            call ftpkyj(ounit,'NAXIS2',nrows,comm,status)
    else
            status=218
    end if

    if (status > 0)return

    if (pcount >= 0)then
            comm='size of special data area'
            call ftpkyj(ounit,'PCOUNT',pcount,comm,status)
    else
            status=214
    end if

    comm='one data group (required keyword)'
    call ftpkyj(ounit,'GCOUNT',1,comm,status)

    comm='number of fields in each row'
    call ftpkyj(ounit,'TFIELDS',nfield,comm,status)

    if (status > 0)return

    do 20 i=1,nfield
        if (ttype(i) /= ' ' .and. ichar(ttype(i)(1:1))/=0)then
            comm='label for field '
            write(comm(17:19),1000)i
1000            format(i3)
            call ftpkns(ounit,'TTYPE',i,1,ttype(i),comm,status)
        end if

        comm='data format of field'
!           make sure format characters are in upper case:
        tfm=tform(i)
        call ftupch(tfm)

!           Add datatype to the comment string:
        call ftbnfm(tfm,dtype,rcount,width,status)
        if (dtype == 21)then
            comm(21:)=': 2-byte INTEGER'
        else if(dtype == 41)then
            comm(21:)=': 4-byte INTEGER'
        else if(dtype == 42)then
            comm(21:)=': 4-byte REAL'
        else if(dtype == 82)then
            comm(21:)=': 8-byte DOUBLE'
        else if(dtype == 16)then
            comm(21:)=': ASCII Character'
        else if(dtype == 14)then
            comm(21:)=': 1-byte LOGICAL'
        else if(dtype == 11)then
            comm(21:)=': BYTE'
        else if(dtype == 1)then
            comm(21:)=': BIT'
        else if(dtype == 83)then
            comm(21:)=': COMPLEX'
        else if(dtype == 163)then
            comm(21:)=': DOUBLE COMPLEX'
        else if(dtype < 0)then
            comm(21:)=': variable length array'
        end if

        call ftpkns(ounit,'TFORM',i,1,tfm,comm,status)

        if (tunit(i) /= ' ' .and. ichar(tunit(i)(1:1))/=0)then
            comm='physical unit of field'
            call ftpkns(ounit,'TUNIT',i,1,tunit(i),comm,status)
        end if
        if (status > 0)return
20      continue

    if (extnam /= ' ' .and. ichar(extnam(1:1)) /= 0)then
            comm='name of this binary table extension'
            call ftpkys(ounit,'EXTNAME',extnam,comm,status)
    end if
end
subroutine ftphis(ounit,histry,status)
!
!*******************************************************************************
!
!! FTPHIS writes a HISTORY record to the FITS header.
!
!       ounit   i  fortran output unit number
!       histry  c  input history string
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,status,strlen,actlen,i,nkeys,c1,c2
    character ( len = * ) histry
    character*80  rec

    if (status > 0)return

!       find the length of the string, and write it out 70 characters at a time
    nkeys=1
    strlen=len(histry)
    actlen=strlen
    do i=strlen,1,-1
            if (histry(i:i) /= ' ')then
                    actlen=i
                    exit
            end if
    end do

    c1=1
    c2=min(actlen,70)
    nkeys=(actlen-1)/70+1
    do i=1,nkeys
            rec='HISTORY   '//histry(c1:c2)
            call ftprec(ounit,rec,status)
            c1=c1+70
            c2=min(actlen,c2+70)
    end do
end
subroutine ftphpr(ounit,simple,bitpix,naxis,naxes, &
                      pcount,gcount,extend,status)
!
!*******************************************************************************
!
!! FTPHPR writes required primary header keywords.
!
!       ounit   i  fortran output unit number
!       simple  l  does file conform to FITS standard?
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       pcount  i  number of group parameters
!       gcount  i  number of random groups
!       extend  l  may extensions be present in the FITS file?
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,bitpix,naxis,naxes(*),pcount,gcount,status,i,ibuff
    character comm*50,caxis*20,clen*3
    logical simple,extend

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    if (status > 0)return

    ibuff=bufnum(ounit)

    if ( hdend(ibuff) /= hdstrt(ibuff,chdu(ibuff)) )then
!            some keywords have already been written
         status=201
         return
    end if

    if (chdu(ibuff) == 1)then
        if (simple)then
            comm='file does conform to FITS standard'
        else
            comm='file does not conform to FITS standard'
        end if
        call ftpkyl(ounit,'SIMPLE',simple,comm,status)
    else
        comm='IMAGE extension'
        call ftpkys(ounit,'XTENSION','IMAGE',comm,status)
    end if

!       test for legal value of bitpix
    call fttbit(bitpix,status)
    comm='number of bits per data pixel'
    call ftpkyj(ounit,'BITPIX',bitpix,comm,status)
    if (status > 0)go to 900

    if (naxis >= 0 .and. naxis <= 999)then
            comm='number of data axes'
            call ftpkyj(ounit,'NAXIS',naxis,comm,status)
    else
!               illegal value of naxis
            status=212
            write(caxis,1000)naxis
1000            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPHPR ' &
            //'is illegal.')
            go to 900
    end if

    comm='length of data axis'
    do i=1,naxis
            if (naxes(i) >= 0)then
                    if (i <= 9)then
                      write(comm(21:21),1001)i
                    else if (i <= 99)then
                      write(comm(21:22),1002)i
                    else
                      write(comm(21:23),1003)i
                    end if
1001                    format(i1)
1002                    format(i2)
1003                    format(i3)
                    call ftpknj(ounit,'NAXIS',i,1,naxes(i),comm, &
                                status)
            else
!                       illegal NAXISnnn keyword value
                    status=213
                    write(clen,1003)i
                    write(caxis,1000)naxes(i)
    call ftpmsg('In call to FTPHPR, axis '//clen// &
    ' has illegal negative size: '//caxis)
                    go to 900
            end if
    end do

    if (chdu(ibuff) == 1)then
!               only write the EXTEND keyword to primary header if true
            if (extend)then
                    comm='FITS dataset may contain extensions'
                    call ftpkyl(ounit,'EXTEND',extend,comm,status)
            end if

!               write the PCOUNT and GCOUNT values if nonstandard
            if (pcount > 0 .or. gcount > 1)then
                comm='random group records are present'
                call ftpkyl(ounit,'GROUPS',.true.,comm,status)
                comm='number of random group parameters'
                call ftpkyj(ounit,'PCOUNT',pcount,comm,status)
                comm='number of random groups'
                call ftpkyj(ounit,'GCOUNT',gcount,comm,status)
            end if

            call ftpcom(ounit,'FITS (Flexible Image Transport '// &
   'System) format defined in Astronomy and',status)
            call ftpcom(ounit,'Astrophysics Supplement Series '// &
   'v44/p363, v44/p371, v73/p359, v73/p365.',status)
            call ftpcom(ounit,'Contact the NASA Science '// &
   'Office of Standards and Technology for the',status)
            call ftpcom(ounit,'FITS Definition document '// &
   '#100 and other FITS information.',status)

    else
            comm='required keyword; must = 0'
            call ftpkyj(ounit,'PCOUNT',pcount,comm,status)
            comm='required keyword; must = 1'
            call ftpkyj(ounit,'GCOUNT',gcount,comm,status)
    end if

900     continue
end
subroutine ftphps(ounit,bitpix,naxis,naxes,status)
!
!*******************************************************************************
!
!! FTPHPS writes required primary header keywords.
!
!       ounit   i  fortran output unit number
!       simple  l  does file conform to FITS standard?
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       pcount  i  number of group parameters
!       gcount  i  number of random groups
!       extend  l  may extensions be present in the FITS file?
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,bitpix,naxis,naxes(*),status

    call ftphpr(ounit,.true.,bitpix,naxis,naxes, &
                      0,1,.true.,status)
end
subroutine ftphtb(ounit,ncols,nrows,nfield,ttype,tbcol, &
    tform,tunit,extnam,status)
!
!*******************************************************************************
!
!! FTPHTB writes required standard header keywords for an ASCII table extension.
!
!       ounit   i  fortran output unit number
!       ncols   i  number of columns in the table
!       nrows   i  number of rows in the table
!       nfield  i  number of fields in the table
!       ttype   c  name of each field (array) (optional)
!       tbcol   i  beginning column of each field (array)
!       tform   c  Fortran-77 format of each field (array)
!       tunit   c  units of each field (array) (optional)
!       extnam  c  name of table extension (optional)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,ncols,nrows,nfield,tbcol(*),status,i
    integer nkeys,nmore
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam
    character comm*48,tfm*20

    if (status > 0)return

    call ftghsp(ounit,nkeys,nmore,status)
    if (nkeys /= 0)then
!            some keywords have already been written
         status=201
         return
    end if

    comm='ASCII table extension'
    call ftpkys(ounit,'XTENSION','TABLE',comm,status)

    comm='8-bit ASCII characters'
    call ftpkyj(ounit,'BITPIX',8,comm,status)

    comm='2-dimensional ASCII table'
    call ftpkyj(ounit,'NAXIS',2,comm,status)

    if (status > 0)return

    if (ncols >= 0)then
            comm='width of table in characters'
            call ftpkyj(ounit,'NAXIS1',ncols,comm,status)
    else
!               illegal table width
            status=217
    call ftpmsg('ASCII table has negative width (NAXIS1) in'// &
    ' call to FTPHTB')
            return
    end if

    if (status > 0)return

    if (nrows >= 0)then
            comm='number of rows in table'
            call ftpkyj(ounit,'NAXIS2',nrows,comm,status)
    else
!               illegal number of rows in table
            status=218
    call ftpmsg('ASCII table has negative number of rows in'// &
    ' call to FTPHTB')
    end if

    if (status > 0)return

    comm='no group parameters (required keyword)'
    call ftpkyj(ounit,'PCOUNT',0,comm,status)

    comm='one data group (required keyword)'
    call ftpkyj(ounit,'GCOUNT',1,comm,status)

    if (status > 0)return

    if (nfield >= 0)then
            comm='number of fields in each row'
            call ftpkyj(ounit,'TFIELDS',nfield,comm,status)
    else
!               illegal number of fields
            status=216
    call ftpmsg('ASCII table has negative number of fields in'// &
    ' call to FTPHTB')
    end if

    if (status > 0)return

    do i=1,nfield
        if (ttype(i) /= ' ' .and. ichar(ttype(i)(1:1))/=0)then
            comm='label for field '
            write(comm(17:19),1000)i
1000            format(i3)
            call ftpkns(ounit,'TTYPE',i,1,ttype(i),comm,status)
        end if

        comm='beginning column of field '
        write(comm(27:29),1000)i
        call ftpknj(ounit,'TBCOL',i,1,tbcol(i),comm,status)

        comm='Fortran-77 format of field'
!           make sure format characters are in upper case:
        tfm=tform(i)
        call ftupch(tfm)
        call ftpkns(ounit,'TFORM',i,1,tfm,comm,status)

        if (tunit(i) /= ' ' .and. ichar(tunit(i)(1:1))/=0)then
            comm='physical unit of field'
            call ftpkns(ounit,'TUNIT',i,1,tunit(i),comm,status)
        end if
    if (status > 0)return
    end do

    if (extnam /= ' ' .and. ichar(extnam(1:1)) /= 0)then
            comm='name of this ASCII table extension'
            call ftpkys(ounit,'EXTNAME',extnam,comm,status)
    end if
end
subroutine ftpi1b(ounit,nvals,incre,chbuff,status)
!
!*******************************************************************************
!
!! FTPI1B writes an array of Integer*1 bytes to the output FITS file.
!
    integer nvals,incre,ounit,status,offset
    character chbuff(nvals)

!       ounit   i  fortran unit number
!       nvals   i  number of pixels in the i2vals array
!       incre   i  byte increment between values
!       chbuff  c*1 array of input byte values
!       status  i  output error status

    if (incre <= 1)then
            call ftpcbf(ounit,nvals,chbuff,status)
    else
!               offset is the number of bytes to move between each value
            offset=incre-1
            call ftpcbo(ounit,1,nvals,offset,chbuff,status)
    end if
end
subroutine ftpini(iunit,status)
!
!*******************************************************************************
!
!! FTPINI initializes the parameters defining the structure of the primary data.
!
!       iunit   i  Fortran I/O unit number
!       OUTPUT PARAMETERS:
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bitpix,naxis,naxes(99),pcnt,gcnt,ttype
    integer blank,bytlen,npix,i,nblank,tstat
    double precision bscale,bzero
    logical simple,extend,groups
    character*8 comm

    if (status > 0)return
    groups=.false.

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       store the type of HDU (0=primary array or image extension)
    hdutyp(ibuff)=0

!       temporarily set the location of the end of the header to a huge number
    hdend(ibuff)=2000000000
    hdstrt(ibuff,chdu(ibuff)+1)=2000000000

!       get the standard header keywords
    tstat=status
    call ftgphx(iunit,99,simple,bitpix,naxis,naxes, &
          pcnt,gcnt,extend,bscale,bzero,blank,nblank,status)

    if (status == 251)then
!               ignore 'unknown extension type' error, and go on
            status=tstat
    else if (status > 0)then
            return
    end if

    if (naxis > 99)then
!               the image array has too many dimensions for me to handle
            status=111
    call ftpmsg('This FITS image has too many dimensions (FTPINI)')
            return
    end if

!       the 'END' record is 80 bytes before the current position, ignoring
!       any trailing blank keywords just before the END keyword.
    hdend(ibuff)=nxthdr(ibuff)-80*(nblank+1)

!       the data unit begins at the beginning of the next logical block
    dtstrt(ibuff)=((nxthdr(ibuff)-80)/2880+1)*2880

!       test for the presence of 'random groups' structure
    if (naxis > 0 .and. naxes(1) == 0)then
            tstat=status
            call ftgkyl(iunit,'GROUPS',groups,comm,status)
            if (status > 0)then
                    status=tstat
                    groups=.false.
            end if
    end if

!       test  bitpix and set the datatype code value
    if (bitpix == 8)then
            ttype=11
            bytlen=1
    else if (bitpix == 16)then
            ttype=21
            bytlen=2
    else if (bitpix == 32)then
            ttype=41
            bytlen=4
    else if (bitpix == -32)then
            ttype=42
            bytlen=4
    else if (bitpix == -64)then
            ttype=82
            bytlen=8
    end if

!       calculate the size of the primary array
    if (naxis == 0)then
            npix=0
    else
            if (groups)then
!                       NAXIS1 = 0 is a special flag for 'random groups'
                    npix=1
            else
                    npix=naxes(1)
            end if

            do 10 i=2,naxis
                    npix=npix*naxes(i)
10              continue
    end if

!       now we know everything about the array; just fill in the parameters:
!       the next HDU begins in the next logical block after the data
    hdstrt(ibuff,chdu(ibuff)+1)= &
    dtstrt(ibuff)+((pcnt+npix)*bytlen*gcnt+2879)/2880*2880

!       initialize the fictitious heap starting address (immediately following
!       the array data) and a zero length heap.  This is used to find the
!   end of the data when checking the fill values in the last block.
    heapsz(ibuff)=0
    theap(ibuff)=(pcnt+npix)*bytlen*gcnt

!       quit if there is no data
    if (naxis == 0)then
            tfield(ibuff)=0
            rowlen(ibuff)=0
            go to 900
    end if

!       the primary array is actually interpreted as a binary table.  There
!       are two columns: the first column contains the
!       group parameters, if any, and the second column contains the
!       primary array of data.  Each group is in a separate row of the table.

    tfield(ibuff)=2
    if (nxtfld + 2 > nf)then
!               too many columns open at one time; exceeded array dimensions
            status=111
    else
            tstart(ibuff)=nxtfld
            nxtfld=nxtfld+2
            tdtype(1+tstart(ibuff))=ttype
            tdtype(2+tstart(ibuff))=ttype
            trept(1+tstart(ibuff))=pcnt
            trept(2+tstart(ibuff))=npix
            tnull(1+tstart(ibuff))=blank
            tnull(2+tstart(ibuff))=blank
            tscale(1+tstart(ibuff))=1.
            tscale(2+tstart(ibuff))=bscale
            tzero(1+tstart(ibuff))=0.
            tzero(2+tstart(ibuff))=bzero
            tbcol(1+tstart(ibuff))=0
            tbcol(2+tstart(ibuff))=pcnt*bytlen
            rowlen(ibuff)=(pcnt+npix)*bytlen
    end if

900     continue
end
subroutine ftpkey(ounit,keywrd,value,comm,status)
!
!*******************************************************************************
!
!! FTPKEY writes a simple FITS keyword record.
!
!  The format used is:
!            "KEYWORD = VALUE / COMMENT"
!               VALUE is assumed to be 20 characters long
!               COMMENT is assumed to be 47 characters long
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       value   c  keyword value   (20 characters, cols. 11-30)
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,value,comm
    integer ounit,status
    character key*8, val*20, com*47

    key=keywrd
    val=value
    com=comm

!       append the 80 characters to the output buffer:
    call ftprec(ounit,key//'= '//val//' / '//com,status)
end
subroutine ftpkls(ounit,keywrd,strval,comm,status)
!
!*******************************************************************************
!
!! FTPKLS writes a character string value to a header record.
!
!  The routine supports the OGIP long string convention.  If the keyword
!  string value is longer than 68 characters (which is the maximum that will
!  fit on a single 80 character keyword record) then the value string will
!       be continued over multiple keywords.  This OGIP convention uses the
!       '&' character at the end of a string to indicate that it is continued
!       on the next keyword.  The name of all the continued keywords is
!       'CONTINUE'.
!
!       The FTPLSW subroutine should be called prior to using this
!       subroutine, to write a warning message in the header
!       describing how the convention works.

!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       strval  c  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Sept 1994

    character ( len = * ) keywrd,comm,strval
    integer ounit,status,lenval,ncomm,nvalue
    integer clen,i,strlen,nseg,c1,c2
    character value*70,keynam*10
    character ( len = 48 ) cmnt

    if (status > 0)return

    keynam=keywrd
    keynam(9:10)='= '
    cmnt=comm

!       find the number of characters in the input string
    clen=len(strval)
    do i=clen,1,-1
            if (strval(i:i) /= ' ')then
                    strlen=i
                    go to 20
            end if
    end do
    strlen=1

!       calculate the number of keywords needed to write the whole string
20      nseg=max(1,(strlen-2)/67+1)

    c1=1
    do i=1,nseg
            c2=min(c1+67,strlen)
!               convert string to quoted character string

!        fts2c was modified on 29 Nov 1994, so this code is no longer needed
!                (remember to declare character*70 ctemp if this code is used)
!                if (i > 1 .and. strval(c1:c1) == ' ')then
!C                   have to preserve leading blanks on continuation cards
!                    ctemp='A'//strval(c1+1:c2)
!                    call fts2c(ctemp,value,lenval,status)
!C                   now reset the first character of the string back to a blank
!                    value(2:2)=' '
!                else

                call fts2c(strval(c1:c2),value,lenval,status)

!                end if

            if (i /= nseg .and. lenval /= 70)then
!                       if the string is continued, preserve trailing blanks
                    value(lenval:69)=' '
                    value(70:70)=''''
                    lenval=70
            end if

!               overwrite last character with a '&' if string is continued.
            if (i < nseg)then
                    value(69:69)='&'
            end if

!               find amount of space left for comment string (assume
!               10 char. for 'keyword = ', and 3 between value and comment)
!               which leaves 67 spaces for the value + comment strings

            nvalue=max(20,lenval)
            ncomm=67-nvalue

!               write the keyword record
            if (ncomm > 0)then
!                       there is space for a comment
                    call ftprec(ounit,keynam// &
                    value(1:nvalue)//' / '//cmnt(1:ncomm),status)
            else
!                       no room for a comment
                    call ftprec(ounit,keynam// &
                    value(1:nvalue)//'   ',status)
            end if

!               initialize for the next segment of the string, if any
            c1=c1+67
            keynam='CONTINUE  '
    end do
end
subroutine ftpknd(ounit,keywrd,nstart,nkey,dval,decim,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKND writes an array of real*8 values to header records in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       dval    d  array of keyword values
!       decim   i  number of decimal places to display in the value field
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm(*)
    integer nstart,nkey,decim,ounit,status,i,j
    double precision dval(*)
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
!       (if the last non-blank character is '&', then it is to be repeated)
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
              call ftpkyd(ounit,keynam,dval(i),decim,comm1,status)
            else
              call ftpkyd(ounit,keynam,dval(i),decim,comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpkne(ounit,keywrd,nstart,nkey,rval,decim,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKNE writes an array of real*4 values to header records in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       rval    r  array of keyword values
!       decim   i  number of decimal places to display in the value field
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm(*)
    integer nstart,nkey,decim,ounit,status,i,j
    real rval(*)
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
!       (if the last non-blank character is '&', then it is to be repeated)
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
              call ftpkye(ounit,keynam,rval(i),decim,comm1,status)
            else
              call ftpkye(ounit,keynam,rval(i),decim,comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpknf(ounit,keywrd,nstart,nkey,rval,decim,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKNF writes an array of real*4 values to header records in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       rval    r  array of keyword values
!       decim   i  number of decimal places to display in the value field
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm(*)
    integer nstart,nkey,decim,ounit,status,i,j
    real rval(*)
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
!       (if the last non-blank character is '&', then it is to be repeated)
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
              call ftpkyf(ounit,keynam,rval(i),decim,comm1,status)
            else
              call ftpkyf(ounit,keynam,rval(i),decim,comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpkng(ounit,keywrd,nstart,nkey,dval,decim,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKNG writes an array of real*8 values to header records in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       dval    d  array of keyword values
!       decim   i  number of decimal places to display in the value field
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm(*)
    integer nstart,nkey,decim,ounit,status,i,j
    double precision dval(*)
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
!       (if the last non-blank character is '&', then it is to be repeated)
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
              call ftpkyg(ounit,keynam,dval(i),decim,comm1,status)
            else
              call ftpkyg(ounit,keynam,dval(i),decim,comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpknj(ounit,keywrd,nstart,nkey,intval,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKNJ writes an array of integer values to header records.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       intval  i  array of keyword values
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm(*)
    integer nstart,nkey,ounit,status,intval(*),i,j
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
!       (if the last non-blank character is '&', then it is to be repeated)
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
               call ftpkyj(ounit,keynam,intval(i),comm1,status)
            else
               call ftpkyj(ounit,keynam,intval(i),comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpknl(ounit,keywrd,nstart,nkey,logval,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKNL writes an array of logical values to header records.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       logval  l  array of keyword values
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm(*)
    integer nstart,nkey,ounit,status,i,j
    logical logval(*)
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
!       (if the last non-blank character is '&', then it is to be repeated)
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
              call ftpkyl(ounit,keynam,logval(i),comm1,status)
            else
              call ftpkyl(ounit,keynam,logval(i),comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpkns(ounit,keywrd,nstart,nkey,strval,comm, &
                      status)
!
!*******************************************************************************
!
!! FTPKNS writes an array of character string values to header records.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       nstart  i  starting sequence number (usually 1)
!       nkey    i  number of keywords to write
!       strval  c  array of keyword values
!       comm    c  array of keyword comments (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,strval(*),comm(*)
    integer nstart,nkey,ounit,status,i,j
    character keynam*8,comm1*48
    logical repeat

    if (status > 0)return

!       check if the first comment string is to be repeated for all keywords
    call ftcrep(comm(1),comm1,repeat)

    j=nstart
    do i=1,nkey
!               construct keyword name:
            call ftkeyn(keywrd,j,keynam,status)

!               write the keyword record
            if (repeat)then
               call ftpkys(ounit,keynam,strval(i),comm1,status)
            else
               call ftpkys(ounit,keynam,strval(i),comm(i),status)
            end if
            if (status > 0)return
            j=j+1
    end do
end
subroutine ftpkyd(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTPKYD writes a double precision value to a header record in E format.
!
!       If it will fit, the value field will be 20 characters wide;
!       otherwise it will be expanded to up to 35 characters, left
!       justified.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (max. 47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim,vlen
    character value*35,key*8,cmnt*48

    key=keywrd
    cmnt=comm

!       convert double precision to E format character string
    call ftd2e(dval,decim,value,vlen,status)

!       write the keyword record
    call ftprec(ounit,key//'= '//value(1:vlen)//' / '//cmnt,status)
end
subroutine ftpkye(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTPKYE writes a real*4 value to a header record in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    real rval
    integer ounit,status,decim
    character value*20

!       convert real to E format character string
    call ftr2e(rval,decim,value,status)

!       write the keyword record
    call ftpkey(ounit,keywrd,value,comm,status)
end
subroutine ftpkyf(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTPKYF writes a real*4 value to a header record in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    real rval
    integer ounit,status,decim
    character value*20

!       convert real to F format character string
    call ftr2f(rval,decim,value,status)

!       write the keyword record
    call ftpkey(ounit,keywrd,value,comm,status)
end
subroutine ftpkyg(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTPKYG writes a double precision value to a header record in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim
    character value*20

!       convert double precision to F format character string
    call ftd2f(dval,decim,value,status)

!       write the keyword record
    call ftpkey(ounit,keywrd,value,comm,status)
end
subroutine ftpkyj(ounit,keywrd,intval,comm,status)
!
!*******************************************************************************
!
!! FTPKYJ writes an integer value to a header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       intval  i  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer ounit,status,intval
    character value*20

!       convert integer to character string
    call fti2c(intval,value,status)

!       write the keyword record
    call ftpkey(ounit,keywrd,value,comm,status)
end
subroutine ftpkyl(ounit,keywrd,logval,comm,status)
!
!*******************************************************************************
!
!! FTPKYL writes a logical value to a header record .
!
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       logval  l  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keywrd,comm
    integer ounit,status
    logical logval
    character ( len = 20 ) value

!       convert logical to character string
    call ftl2c(logval,value,status)

!       write the keyword record
    call ftpkey(ounit,keywrd,value,comm,status)
end
subroutine ftpkys(ounit,keywrd,strval,comm,status)
!
!*******************************************************************************
!
!! FTPKYS writes a character string value to a header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       strval  c  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!       modified 6/93 to handle long string values by continuing the
!       string onto subsequent comment keywords (with a blank keyword name)

!       Modified again in 9/94 to remove support for long string values;
!       Now, one must call ftpkls to write a long string values.

    character ( len = * ) keywrd,comm,strval
    integer ounit,status,lenval,ncomm,nvalue
    character strtmp*68,value*70,keynam*8,cmnt*48

    if (status > 0)return

    strtmp=strval
    keynam=keywrd
    cmnt=comm

!       convert string to quoted character string (max length = 70 characters)
    call fts2c(strtmp,value,lenval,status)

!       find amount of space left for comment string
!       (assume 10 char. for 'keyword = ', and 3 between value and comment)
!       which leaves 67 spaces for the value string + comment string
    nvalue=max(20,lenval)
    ncomm=67-nvalue

!       write the keyword record
    if (ncomm > 0)then
!         there is space for a comment
      call ftprec(ounit, &
      keynam//'= '//value(1:nvalue)//' / '//cmnt(1:ncomm),status)
    else
!         no room for a comment
      call ftprec(ounit, &
      keynam//'= '//value(1:nvalue)//'   ',status)
    end if
end
subroutine ftpkyt(ounit,keywrd,jval,dval,comm,status)
!
!*******************************************************************************
!
!! FTPKYT writes an integer and double precision fraction to the FITS header.
!
!  The routine concatenates the integer with the double precision fraction
!       and writes it to the FITS header along with the comment string
!       The value will be displayed in F28.16 format
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       jval    i  integer part of the keyword value
!       dval    d  fractional part of the keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Sept 1992

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,jval,status,dlen,i,fchar
    character dstr*35,jstr*20,key*8,cmnt*48

    if (status > 0)return

    if (dval >= 1.0  .or. dval <  0.)then
            status = 402
    end if

    key=keywrd
    cmnt=comm

!       convert integer to C*20 character string
    call fti2c(jval,jstr,status)

!       ignore leading spaces
    fchar=10
    do i=10,20
      if (jstr(i:i) /= ' ')then
         fchar = i
         exit
      end if
    end do

!       convert double precision to E23.16 format character string
    call ftd2e(dval,15,dstr,dlen,status)

!       write the concatinated keyword record
    call ftprec(ounit,key//'= '//jstr(fchar:20)//'.'// &
     dstr(1:1)//dstr(3:17)//' / '//cmnt,status)
end
subroutine ftpkyu(ounit,keywrd,comm,status)
!
!*******************************************************************************
!
!! FTPKYU writes a null-valued keyword to a header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1997

    character ( len = * ) keywrd,comm
    integer ounit,status
    character keynam*8,card*80

    if (status > 0)return

    keynam=keywrd
    card=keynam//'=                      / '//comm

    call ftprec(ounit,card,status)
end
subroutine ftplsw(ounit,status)
!
!*******************************************************************************
!
!! FTPLSW "Puts Long String Warning".
!
!       write the LONGSTRN keyword and a few COMMENT keywords to the header
!       (if they don't already exist) to warn users that this FITS file
!       may use the OGIP long string convention.

!       This subroutine should be called whenever FTPKLS is called.

    integer ounit,status,tstat
    character value*8
    character ( len = 8 ) comm

    if (status > 0)return

    tstat=status
    call ftgkys(ounit,'LONGSTRN',value,comm,status)
    if (status == 0)then
!             The keyword already exists so just exit
          return
     end if

     status=tstat
     call ftpkys(ounit,'LONGSTRN','OGIP 1.0', &
     'The HEASARC Long String Convention may be used.',status)

     call ftpcom(ounit, &
   'This FITS file may contain long string keyword values that are' &
    ,status)
       call ftpcom(ounit, &
   'continued over multiple keywords.  The HEASARC convention uses' &
    //' the &',status)
        call ftpcom(ounit, &
   'character at the end of each substring which is then continued' &
    ,status)
        call ftpcom(ounit, &
   'on the next keyword which has the name CONTINUE.' &
    ,status)
end
subroutine ftpmsg(text)
!
!*******************************************************************************
!
!! FTPMSG puts an error message onto stack.
!
    character ( len = * ) text
    call ftxmsg(1,text)
end
subroutine ftpnul(ounit,blank,status)
!
!*******************************************************************************
!
!! FTPNUL defines the null value for an integer primary array.
!
!       ounit   i  Fortran I/O unit number
!       blank   i  the value to be use to signify undefined data
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,blank,status

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,i,ngroup

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)
    if (status > 0)return

!       test for proper HDU type
    if (hdutyp(ibuff) /= 0)then
        status=233
        return
    end if

!       the primary array is actually interpreted as a binary table.  There
!       are two columns for each group: the first column contains the
!       group parameters, if any, and the second column contains the
!       primary array of data.

    ngroup=tfield(ibuff)/2
    do 10 i=1,ngroup
            tnull(i*2+tstart(ibuff))=blank
10      continue
end
subroutine ftppnb(ounit,group,felem,nelem,array,nulval,status)
!
!*******************************************************************************
!
!! FTPPNB writes an array of c*1 (byte) values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same as the
!       array being written).  Any input pixels equal to the value of NULVAL
!       will be replaced by the appropriate null value in the output FITS file.

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   c*1  the array of values to be written
!       nulval  c*1  pixel value used to represent an undefine pixel
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,group,felem,nelem,status,row
    character array(*),nulval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcnb(ounit,2,row,felem,nelem,array,nulval,status)
end
subroutine ftppnd(ounit,group,felem,nelem,array,nulval,status)
!
!*******************************************************************************
!
!! FTPPND writes an array of double precision values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same as the
!       array being written).  Any input pixels equal to the value of NULVAL
!       will be replaced by the appropriate null value in the output FITS file.

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   d  the array of values to be written
!       nulval  d  pixel value used to represent an undefine pixel
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,group,felem,nelem,status,row
    double precision array(*),nulval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcnd(ounit,2,row,felem,nelem,array,nulval,status)
end
subroutine ftppne(ounit,group,felem,nelem,array,nulval,status)
!
!*******************************************************************************
!
!! FTPPNE writes an array of real values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same as the
!       array being written).  Any input pixels equal to the value of NULVAL
!       will be replaced by the appropriate null value in the output FITS file.

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   r  the array of values to be written
!       nulval  r  pixel value used to represent an undefine pixel
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,group,felem,nelem,status,row
    real array(*),nulval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcne(ounit,2,row,felem,nelem,array,nulval,status)
end
subroutine ftppni(ounit,group,felem,nelem,array,nulval,status)
!
!*******************************************************************************
!
!! FTPPNI writes an array of i*2 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same as the
!       array being written).  Any input pixels equal to the value of NULVAL
!       will be replaced by the appropriate null value in the output FITS file.

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   i*2  the array of values to be written
!       nulval  i*2  pixel value used to represent an undefine pixel
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,group,felem,nelem,status,row
    integer*2 array(*),nulval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcni(ounit,2,row,felem,nelem,array,nulval,status)
end
subroutine ftppnj(ounit,group,felem,nelem,array,nulval,status)
!
!*******************************************************************************
!
!! FTPPNJ writes an array of i values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same as the
!       array being written).  Any input pixels equal to the value of NULVAL
!       will be replaced by the appropriate null value in the output FITS file.

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   i  the array of values to be written
!       nulval  i  pixel value used to represent an undefine pixel
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1994

    integer ounit,group,felem,nelem,status,row
    integer array(*),nulval

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcnj(ounit,2,row,felem,nelem,array,nulval,status)
end
subroutine ftpprb(ounit,group,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPPRB writes an array of byte values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   b  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,felem,nelem,status,row

    character array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpclb(ounit,2,row,felem,nelem,array,status)
end
subroutine ftpprd(ounit,group,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPPRD writes an array of r*8 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   d  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,felem,nelem,status,row
    double precision array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcld(ounit,2,row,felem,nelem,array,status)
end
subroutine ftppre(ounit,group,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPPRE writes an array of r*4 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   r  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,felem,nelem,status,row
    real array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcle(ounit,2,row,felem,nelem,array,status)
end
subroutine ftpprh(ounit,simple,bitpix,naxis,naxes, &
                      pcount,gcount,extend,status)
!
!*******************************************************************************
!
!! FTPPRH is obsolete; call FTPHPR instead.
!

    integer ounit,bitpix,naxis,naxes(*),pcount,gcount,status
    logical simple,extend

    call ftphpr(ounit,simple,bitpix,naxis,naxes, &
                      pcount,gcount,extend,status)
end
subroutine ftppri(ounit,group,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPPRI writes an array of i*2 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   i*2  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,felem,nelem,status,row
    integer*2 array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpcli(ounit,2,row,felem,nelem,array,status)
end
subroutine ftpprj(ounit,group,felem,nelem,array,status)
!
!*******************************************************************************
!
!! FTPPRJ writes an array of i*4 values into the primary array.
!
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being written).

!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be written
!       array   i  the array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,felem,nelem,status,row
    integer array(*)

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpclj(ounit,2,row,felem,nelem,array,status)
end
subroutine ftppru(ounit,group,felem,nelem,status)
!
!*******************************************************************************
!
!! FTPPRU sets elements of the primary array equal to the undefined value.
!
!       ounit   i  Fortran output unit number
!       group   i  number of the data group, if any
!       felem   i  the first pixel to be written (this routine treats
!                  the primary array a large one dimensional array of
!                  values, regardless of the actual dimensionality).
!       nelem   i  number of data elements to be set to undefined
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,group,felem,nelem,status,row

!       the primary array is represented as a binary table:
!               each group of the primary array is a row in the table,
!               where the first column contains the group parameters
!               and the second column contains the image itself
    row=max(group,1)
    call ftpclu(ounit,2,row,felem,nelem,status)
end
subroutine ftprec(ounit,record,status)
!
!*******************************************************************************
!
!! FTPREC writes an 80 character record to the FITS header.
!
!       ounit   i  fortran output unit number
!       record  c  input 80 character header record
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) record
    character*80  rec
    integer ounit,status,ibuff

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    if (status > 0)return

!       get the number of the data buffer used for this unit
    ibuff=bufnum(ounit)

    if (dtstrt(ibuff) > 0 &
      .and.(dtstrt(ibuff)-hdend(ibuff)) <= 80)then
!               not enough room in the header for another keyword

!               try getting more header space
            call ftiblk(ounit,1,0,status)
            if (status > 0)then
                    go to 900
            end if
    end if

    rec=record

!       make sure keyword name is in upper case
    call ftupch(rec(1:8))

!       test that keyword name contains only legal characters
    call fttkey(rec(1:8),status)

!       test that the rest of the record contains only legal values
    call fttrec(rec(9:80),status)

!       position the I/O pointer to the end of the header
    call ftmbyt(ounit,hdend(ibuff),.true.,status)

!       append the 80 characters to the output buffer:
    call ftpcbf(ounit,80,rec,status)
    if (status > 0)go to 900

!       increment the pointer to the last header record
    hdend(ibuff)=hdend(ibuff)+80

!       the following statement was added in v4.00 and removed again
!       in v4.09.  There appears to be no good reason to reset the
!       'next keyword' pointer after appending a new keyword to the
!       header, since this effectively just resets the pointer to the
!       beginning of the header.
!        nxthdr(ibuff)=hdend(ibuff)

900     continue
end
subroutine ftprsv(keyin,lenval,status)
!
!*******************************************************************************
!
!! FTPRSV finds the length of the keyword+value string in a keyword record.
!
!       keyrec  c  80 column header record
!       OUTPUT PARAMETERS:
!       lenval  i  output length of keyword+value string
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keyin
    integer lenval,status,j,c1
    character*80 keyrec

    if (status > 0)return

    keyrec=keyin
    if (keyrec(1:8) =='COMMENT ' .or. keyrec(1:8)=='HISTORY ' &
    .or. keyrec(1:8)=='END     ' .or. keyrec(1:8)=='        ') &
    then
!           this is a COMMENT or HISTORY record, with no value
         lenval=8
    else if (keyrec(9:10) == '= ')then
!           this keyword has a value field; now find the first character:
        do j=10,80
            if (keyrec(j:j) /= ' ')then
                    c1=j
                    go to 15
            end if
        end do
!           error: value is blank
        status=204
        call ftpmsg('The keyword '//keyrec(1:8)// &
        ' has no value string after the equal sign:')
        call ftpmsg(keyrec)
        return

15          if (keyrec(c1:c1) == '''')then
!               This is a string value.
!               Work forward to find a single quote.  Two single quotes
!               in succession is to be interpreted as a literal single
!               quote character as part of the character string, not as
!               the end of the character string.  Everything to the right
!               of the closing quote is assumed to be the comment.
            do 20 j=c1+1,80
                if (keyrec(j:j) == '''')then
                    if (j<80 .and. keyrec(j+1:j+1)=='''')then
!                               found 2 successive quote characters; this is
!                               interpreted as a literal quote character
                    else
                            lenval=max(30,j)
                            go to 30
                    end if
                end if
20              continue
!               error: no closing quote character
            status=205
        call ftpmsg('The following Keyword value string has '// &
              'no closing quote:')
        call ftpmsg(keyrec)
            return
        else
!               This is either an integer, floating point, or logical value.
!               Extract the first token as the value; remainder = comment
            do 25 j=c1,80
                if (keyrec(j:j) == ' ')then
                    lenval=j-1
                    go to 30
                end if
25              continue
!               the first token went all the way to column 80:
            lenval=80
        end if
    else
!               illegal keyword record format; must have '= ' in columns 9-10
!                status=210
!            Modified July 1993:  this is actually not an error.  The
!            keyword should simply be interpreted as a comment.
         lenval=8
    end if
30      continue
end
subroutine ftpscl(ounit,bscale,bzero,status)
!
!*******************************************************************************
!
!! FTPSCL defines the scaling factor for the primary header data.
!
!       ounit   i  Fortran I/O unit number
!       bscale  d  scaling factor
!       bzero   d  scaling zero point
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,status
    double precision bscale,bzero

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,i,ngroup

    if (status > 0)return

    if (bscale == 0.)then
!               illegal bscale value
            status=322
            return
    end if

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)
    if (status > 0)return

!       test for proper HDU type
    if (hdutyp(ibuff) /= 0)then
        status=233
        return
    end if

!       the primary array is actually interpreted as a binary table.  There
!       are two columns for each group: the first column contains the
!       group parameters, if any, and the second column contains the
!       primary array of data.
    ngroup=tfield(ibuff)/2
    do 10 i=1,ngroup
            tscale(i*2+tstart(ibuff))=bscale
            tzero(i*2+tstart(ibuff))=bzero
10      continue
end
subroutine ftpssb(iunit,group,naxis,naxes,fpixel,lpixel, &
                      array,status)
!
!*******************************************************************************
!
!! FTPSSB writes a subsection of byte values to the primary array.
!
!       A subsection is defined to be any contiguous rectangular
!       array of pixels within the n-dimensional FITS data file.
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran input unit number
!       group   i  number of the data group to be written, if any
!       naxis   i  number of data axes in the FITS array
!       naxes   i  (array) size of each FITS axis
!       fpixel  i  (array) the first pixel in each dimension to be included
!                  in the subsection (first pixel = 1)
!       lpixel  i  (array) the last pixel in each dimension to be included
!                  in the subsection
!       array   c*1  array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    integer iunit,group,naxis,naxes(*),fpixel(*),lpixel(*),status
    character array(*)
    integer fpix(7),irange(7),dimen(7),astart,pstart
    integer off2,off3,off4,off5,off6,off7
    integer st10,st20,st30,st40,st50,st60,st70
    integer st1,st2,st3,st4,st5,st6,st7
    integer i,i1,i2,i3,i4,i5,i6,i7
    character caxis*20

    if (status > 0)return

    if (naxis < 1 .or. naxis > 7)then
!               this routine only supports up to 7 dimensions
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPSSB ' &
            //'is illegal.')
            return
    end if

!       calculate the sizes and number of loops to perform in each dimension
    do 10 i=1,7
         fpix(i)=1
         irange(i)=1
         dimen(i)=1
10      continue

    do 20 i=1,naxis
         fpix(i)=fpixel(i)
         irange(i)=lpixel(i)-fpixel(i)+1
         dimen(i)=naxes(i)
20      continue
    i1=irange(1)

!       compute the pixel offset between each dimension
    off2=     dimen(1)
    off3=off2*dimen(2)
    off4=off3*dimen(3)
    off5=off4*dimen(4)
    off6=off5*dimen(5)
    off7=off6*dimen(6)

    st10=fpix(1)
    st20=(fpix(2)-1)*off2
    st30=(fpix(3)-1)*off3
    st40=(fpix(4)-1)*off4
    st50=(fpix(5)-1)*off5
    st60=(fpix(6)-1)*off6
    st70=(fpix(7)-1)*off7

!       store the initial offset in each dimension
    st1=st10
    st2=st20
    st3=st30
    st4=st40
    st5=st50
    st6=st60
    st7=st70

    astart=1

    do 170 i7=1,irange(7)
    do 160 i6=1,irange(6)
    do 150 i5=1,irange(5)
    do 140 i4=1,irange(4)
    do 130 i3=1,irange(3)
    pstart=st1+st2+st3+st4+st5+st6+st7
    do 120 i2=1,irange(2)
            call ftpprb(iunit,group,pstart,i1, &
                array(astart),status)
            astart=astart+i1
            pstart=pstart+off2
120     continue
    st2=st20
    st3=st3+off3
130     continue
    st3=st30
    st4=st4+off4
140     continue
    st4=st40
    st5=st5+off5
150     continue
    st5=st50
    st6=st6+off6
160     continue
    st6=st60
    st7=st7+off7
170     continue
end
subroutine ftpssd(iunit,group,naxis,naxes,fpixel,lpixel, &
                      array,status)
!
!*******************************************************************************
!
!! FTPSSD writes a subsection of double precision values to the primary array.
!
!       A subsection is defined to be any contiguous rectangular
!       array of pixels within the n-dimensional FITS data file.
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran input unit number
!       group   i  number of the data group to be written, if any
!       naxis   i  number of data axes in the FITS array
!       naxes   i  (array) size of each FITS axis
!       fpixel  i  (array) the first pixel in each dimension to be included
!                  in the subsection (first pixel = 1)
!       lpixel  i  (array) the last pixel in each dimension to be included
!                  in the subsection
!       array   d  array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    integer iunit,group,naxis,naxes(*),fpixel(*),lpixel(*),status
    double precision array(*)
    integer fpix(7),irange(7),dimen(7),astart,pstart
    integer off2,off3,off4,off5,off6,off7
    integer st10,st20,st30,st40,st50,st60,st70
    integer st1,st2,st3,st4,st5,st6,st7
    integer i,i1,i2,i3,i4,i5,i6,i7
    character caxis*20

    if (status > 0)return

    if (naxis < 1 .or. naxis > 7)then
!               this routine only supports up to 7 dimensions
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPSSD ' &
            //'is illegal.')
            return
    end if

!       calculate the sizes and number of loops to perform in each dimension
    do i=1,7
         fpix(i)=1
         irange(i)=1
         dimen(i)=1
    end do

    do 20 i=1,naxis
         fpix(i)=fpixel(i)
         irange(i)=lpixel(i)-fpixel(i)+1
         dimen(i)=naxes(i)
20      continue
    i1=irange(1)

!       compute the pixel offset between each dimension
    off2=     dimen(1)
    off3=off2*dimen(2)
    off4=off3*dimen(3)
    off5=off4*dimen(4)
    off6=off5*dimen(5)
    off7=off6*dimen(6)

    st10=fpix(1)
    st20=(fpix(2)-1)*off2
    st30=(fpix(3)-1)*off3
    st40=(fpix(4)-1)*off4
    st50=(fpix(5)-1)*off5
    st60=(fpix(6)-1)*off6
    st70=(fpix(7)-1)*off7

!       store the initial offset in each dimension
    st1=st10
    st2=st20
    st3=st30
    st4=st40
    st5=st50
    st6=st60
    st7=st70

    astart=1

    do 170 i7=1,irange(7)
    do 160 i6=1,irange(6)
    do 150 i5=1,irange(5)
    do 140 i4=1,irange(4)
    do 130 i3=1,irange(3)
    pstart=st1+st2+st3+st4+st5+st6+st7
    do 120 i2=1,irange(2)
            call ftpprd(iunit,group,pstart,i1, &
                array(astart),status)
            astart=astart+i1
            pstart=pstart+off2
120     continue
    st2=st20
    st3=st3+off3
130     continue
    st3=st30
    st4=st4+off4
140     continue
    st4=st40
    st5=st5+off5
150     continue
    st5=st50
    st6=st6+off6
160     continue
    st6=st60
    st7=st7+off7
170     continue
end
subroutine ftpsse(iunit,group,naxis,naxes,fpixel,lpixel, &
                      array,status)
!
!*******************************************************************************
!
!! FTPSSE writes a subsection of real values to the primary array.
!
!       A subsection is defined to be any contiguous rectangular
!       array of pixels within the n-dimensional FITS data file.
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran input unit number
!       group   i  number of the data group to be written, if any
!       naxis   i  number of data axes in the FITS array
!       naxes   i  (array) size of each FITS axis
!       fpixel  i  (array) the first pixel in each dimension to be included
!                  in the subsection (first pixel = 1)
!       lpixel  i  (array) the last pixel in each dimension to be included
!                  in the subsection
!       array   r  array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    integer iunit,group,naxis,naxes(*),fpixel(*),lpixel(*),status
    real array(*)
    integer fpix(7),irange(7),dimen(7),astart,pstart
    integer off2,off3,off4,off5,off6,off7
    integer st10,st20,st30,st40,st50,st60,st70
    integer st1,st2,st3,st4,st5,st6,st7
    integer i,i1,i2,i3,i4,i5,i6,i7
    character caxis*20

    if (status > 0)return

    if (naxis < 1 .or. naxis > 7)then
!               this routine only supports up to 7 dimensions
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPSSE ' &
            //'is illegal.')
            return
    end if

!       calculate the sizes and number of loops to perform in each dimension
    do 10 i=1,7
         fpix(i)=1
         irange(i)=1
         dimen(i)=1
10      continue

    do 20 i=1,naxis
         fpix(i)=fpixel(i)
         irange(i)=lpixel(i)-fpixel(i)+1
         dimen(i)=naxes(i)
20      continue
    i1=irange(1)

!       compute the pixel offset between each dimension
    off2=     dimen(1)
    off3=off2*dimen(2)
    off4=off3*dimen(3)
    off5=off4*dimen(4)
    off6=off5*dimen(5)
    off7=off6*dimen(6)

    st10=fpix(1)
    st20=(fpix(2)-1)*off2
    st30=(fpix(3)-1)*off3
    st40=(fpix(4)-1)*off4
    st50=(fpix(5)-1)*off5
    st60=(fpix(6)-1)*off6
    st70=(fpix(7)-1)*off7

!       store the initial offset in each dimension
    st1=st10
    st2=st20
    st3=st30
    st4=st40
    st5=st50
    st6=st60
    st7=st70

    astart=1

    do 170 i7=1,irange(7)
    do 160 i6=1,irange(6)
    do 150 i5=1,irange(5)
    do 140 i4=1,irange(4)
    do 130 i3=1,irange(3)
    pstart=st1+st2+st3+st4+st5+st6+st7
    do 120 i2=1,irange(2)
            call ftppre(iunit,group,pstart,i1, &
                array(astart),status)
            astart=astart+i1
            pstart=pstart+off2
120     continue
    st2=st20
    st3=st3+off3
130     continue
    st3=st30
    st4=st4+off4
140     continue
    st4=st40
    st5=st5+off5
150     continue
    st5=st50
    st6=st6+off6
160     continue
    st6=st60
    st7=st7+off7
170     continue
end
subroutine ftpssi(iunit,group,naxis,naxes,fpixel,lpixel, &
                      array,status)
!
!*******************************************************************************
!
!! FTPSSI writes a subsection of integer*2 values to the primary array.
!
!       A subsection is defined to be any contiguous rectangular
!       array of pixels within the n-dimensional FITS data file.
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran input unit number
!       group   i  number of the data group to be written, if any
!       naxis   i  number of data axes in the FITS array
!       naxes   i  (array) size of each FITS axis
!       fpixel  i  (array) the first pixel in each dimension to be included
!                  in the subsection (first pixel = 1)
!       lpixel  i  (array) the last pixel in each dimension to be included
!                  in the subsection
!       array   i*2  array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    integer iunit,group,naxis,naxes(*),fpixel(*),lpixel(*),status
    integer*2 array(*)
    integer fpix(7),irange(7),dimen(7),astart,pstart
    integer off2,off3,off4,off5,off6,off7
    integer st10,st20,st30,st40,st50,st60,st70
    integer st1,st2,st3,st4,st5,st6,st7
    integer i,i1,i2,i3,i4,i5,i6,i7
    character caxis*20

    if (status > 0)return

    if (naxis < 1 .or. naxis > 7)then
!               this routine only supports up to 7 dimensions
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPSSI ' &
            //'is illegal.')
            return
    end if

!       calculate the sizes and number of loops to perform in each dimension
    do 10 i=1,7
         fpix(i)=1
         irange(i)=1
         dimen(i)=1
10      continue

    do 20 i=1,naxis
         fpix(i)=fpixel(i)
         irange(i)=lpixel(i)-fpixel(i)+1
         dimen(i)=naxes(i)
20      continue
    i1=irange(1)

!       compute the pixel offset between each dimension
    off2=     dimen(1)
    off3=off2*dimen(2)
    off4=off3*dimen(3)
    off5=off4*dimen(4)
    off6=off5*dimen(5)
    off7=off6*dimen(6)

    st10=fpix(1)
    st20=(fpix(2)-1)*off2
    st30=(fpix(3)-1)*off3
    st40=(fpix(4)-1)*off4
    st50=(fpix(5)-1)*off5
    st60=(fpix(6)-1)*off6
    st70=(fpix(7)-1)*off7

!       store the initial offset in each dimension
    st1=st10
    st2=st20
    st3=st30
    st4=st40
    st5=st50
    st6=st60
    st7=st70

    astart=1

    do 170 i7=1,irange(7)
    do 160 i6=1,irange(6)
    do 150 i5=1,irange(5)
    do 140 i4=1,irange(4)
    do 130 i3=1,irange(3)
    pstart=st1+st2+st3+st4+st5+st6+st7
    do 120 i2=1,irange(2)
            call ftppri(iunit,group,pstart,i1, &
                array(astart),status)
            astart=astart+i1
            pstart=pstart+off2
120     continue
    st2=st20
    st3=st3+off3
130     continue
    st3=st30
    st4=st4+off4
140     continue
    st4=st40
    st5=st5+off5
150     continue
    st5=st50
    st6=st6+off6
160     continue
    st6=st60
    st7=st7+off7
170     continue
end
subroutine ftpssj(iunit,group,naxis,naxes,fpixel,lpixel, &
                      array,status)
!
!*******************************************************************************
!
!! FTPSSJ writes a subsection of integer values to the primary array.
!
!       A subsection is defined to be any contiguous rectangular
!       array of pixels within the n-dimensional FITS data file.
!       Data conversion and scaling will be performed if necessary
!       (e.g, if the datatype of the FITS array is not the same
!       as the array being read).

!       iunit   i  Fortran input unit number
!       group   i  number of the data group to be written, if any
!       naxis   i  number of data axes in the FITS array
!       naxes   i  (array) size of each FITS axis
!       fpixel  i  (array) the first pixel in each dimension to be included
!                  in the subsection (first pixel = 1)
!       lpixel  i  (array) the last pixel in each dimension to be included
!                  in the subsection
!       array   i  array of values to be written
!       status  i  returned error stataus

!       written by Wm Pence, HEASARC/GSFC, Feb 1992

    integer iunit,group,naxis,naxes(*),fpixel(*),lpixel(*),status
    integer array(*)
    integer fpix(7),irange(7),dimen(7),astart,pstart
    integer off2,off3,off4,off5,off6,off7
    integer st10,st20,st30,st40,st50,st60,st70
    integer st1,st2,st3,st4,st5,st6,st7
    integer i,i1,i2,i3,i4,i5,i6,i7
    character caxis*20

    if (status > 0)return

    if (naxis < 1 .or. naxis > 7)then
!               this routine only supports up to 7 dimensions
            status=320
            write(caxis,1001)naxis
1001            format(i20)
            call ftpmsg('NAXIS ='//caxis//' in the call to FTPSSJ ' &
            //'is illegal.')
            return
    end if

!       calculate the sizes and number of loops to perform in each dimension
    do 10 i=1,7
         fpix(i)=1
         irange(i)=1
         dimen(i)=1
10      continue

    do 20 i=1,naxis
         fpix(i)=fpixel(i)
         irange(i)=lpixel(i)-fpixel(i)+1
         dimen(i)=naxes(i)
20      continue
    i1=irange(1)

!       compute the pixel offset between each dimension
    off2=     dimen(1)
    off3=off2*dimen(2)
    off4=off3*dimen(3)
    off5=off4*dimen(4)
    off6=off5*dimen(5)
    off7=off6*dimen(6)

    st10=fpix(1)
    st20=(fpix(2)-1)*off2
    st30=(fpix(3)-1)*off3
    st40=(fpix(4)-1)*off4
    st50=(fpix(5)-1)*off5
    st60=(fpix(6)-1)*off6
    st70=(fpix(7)-1)*off7

!       store the initial offset in each dimension
    st1=st10
    st2=st20
    st3=st30
    st4=st40
    st5=st50
    st6=st60
    st7=st70

    astart=1

    do 170 i7=1,irange(7)
    do 160 i6=1,irange(6)
    do 150 i5=1,irange(5)
    do 140 i4=1,irange(4)
    do 130 i3=1,irange(3)
    pstart=st1+st2+st3+st4+st5+st6+st7
    do 120 i2=1,irange(2)
            call ftpprj(iunit,group,pstart,i1, &
                array(astart),status)
            astart=astart+i1
            pstart=pstart+off2
120     continue
    st2=st20
    st3=st3+off3
130     continue
    st3=st30
    st4=st4+off4
140     continue
    st4=st40
    st5=st5+off5
150     continue
    st5=st50
    st6=st6+off6
160     continue
    st6=st60
    st7=st7+off7
170     continue
end
subroutine ftpsvc(keyin,value,comm,status)
!
!*******************************************************************************
!
!! FTPSVC parses the header record to find value and comment strings.
!
!       keyrec  c  80 column header record
!       OUTPUT PARAMETERS:
!       value   c  output keyword value string
!       comm    c  output keyword comment string
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    character ( len = * ) keyin,value,comm
    character*80 keyrec,keytmp,ctemp
    integer status,j,k,c1

    if (status > 0)return

    keyrec=keyin
    if (keyrec(1:8) =='COMMENT ' .or. keyrec(1:8)=='HISTORY ' &
    .or. keyrec(1:8)=='END     ' .or. keyrec(1:8)=='        ') &
    then
!           this is a COMMENT or HISTORY record, with no value
        value=' '
        comm=keyrec(9:80)
    else if (keyrec(9:10) == '= ')then
!           this keyword has a value field; now find the first character:
        do j=10,80
            if (keyrec(j:j) /= ' ')then
                    c1=j
                    go to 15
            end if
        end do

!       the absence of a value string is legal, and simply indicates
!       that the keyword value is undefined.  Don't write an error
!       message in this case.

!            status=204
!            call ftpmsg('The keyword '//keyrec(1:8)//
!     &      ' has no value string after the equal sign:')
!            call ftpmsg(keyrec)

        value=' '
        comm=' '
        return

15          if (keyrec(c1:c1) == '/')then
!               keyword has no defined value (has a null value)
            value=' '
            ctemp=keyrec(c1:80)
        else if (keyrec(c1:c1) == '''')then
!               This is a string value.
!               Work forward to find a single quote.  Two single quotes
!               in succession is to be interpreted as a literal single
!               quote character as part of the character string, not as
!               the end of the character string.  Everything to the right
!               of the closing quote is assumed to be the comment.
!               First, copy input to temporary string variable
            keytmp=keyrec
            do 20 j=c1+1,80
                if (keytmp(j:j) == '''')then
                    if (j<80 .and. keytmp(j+1:j+1)=='''')then
!                               found 2 successive quote characters; this is
!                               interpreted as a literal quote character; remove
!                               one of the quotes from the string, and continue
!                               searching for the closing quote character:
                            do 18 k=j+2,80
                                keytmp(k-1:k-1)=keytmp(k:k)
18                              continue
                            keytmp(80:80)=' '
                    else
                            value=keytmp(c1:j)
                            if (j < 80)then
                                    ctemp=keytmp(j+1:80)
                            else
                                    ctemp=' '
                            end if
                            go to 30
                    end if
                end if
20              continue
!               error: no closing quote character
            status=205
        call ftpmsg('The following Keyword value string has '// &
              'no closing quote:')
        call ftpmsg(keyrec)
            return
        else
!               This is either an integer, floating point, or logical value.
!               Extract the first token as the value; remainder = comment
            do 25 j=c1,80
                if (keyrec(j:j) == ' ')then
                    value=keyrec(c1:j-1)
                    ctemp=keyrec(j+1:80)
                    go to 30
                end if
25              continue
!               the first token went all the way to column 80:
            value=keyrec(c1:80)
            ctemp=' '
        end if

30          comm=' '
!           look for first character in the comment string
        do 40 j=1,78
            if (ctemp(j:j)/=' ')then
                    if (ctemp(j:j)=='/')then
!                            ignore first space, if it exists
                         if (ctemp(j+1:j+1) == ' ')then
                            comm=ctemp(j+2:80)
                         else
                            comm=ctemp(j+1:80)
                         end if
                    else
                            comm=ctemp(j:80)
                    end if
                    go to 50
            end if
40          continue
    else
!           illegal keyword record format; must have '= ' in columns 9-10
!           status=210
!           Modified July 1993:  this is actually not an error.  The
!           keyword should simply be interpreted as a comment.
        value=' '
        comm=keyrec(9:80)
    end if
50      continue
end
subroutine ftptbb(iunit,frow,fchar,nchars,value,status)
!
!*******************************************************************************
!
!! FTPTBB writes a consecutive string of bytes to an ascii or binary table.
!
!       This will span multiple rows of the table if NCHARS+FCHAR is
!       greater than the length of a row.

!       iunit   i  fortran unit number
!       frow    i  starting row number (1st row = 1)
!       fchar   i  starting byte in the row to write (1st character=1)
!       nchars  i  number of bytes to write (can span multiple rows)
!       value   i  array of bytes to write
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1991

    integer iunit,frow,fchar,nchars,status
    integer value(*)

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bstart

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check for errors
    if (nchars <= 0)then
!               zero or negative number of character requested
            return
    else if (frow < 1)then
!               error: illegal first row number
            status=307
            return
    else if (fchar < 1)then
!               error: illegal starting character
            status=308
            return
    end if

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(iunit,status)

!       move the i/o pointer to the start of the sequence of characters
    bstart=dtstrt(ibuff)+(frow-1)*rowlen(ibuff)+fchar-1
    call ftmbyt(iunit,bstart,.true.,status)

!       put the string of bytes
    call ftpbyt(iunit,nchars,value,status)
end
subroutine ftptbh(ounit,ncols,nrows,nfield,ttype,tbcol, &
    tform,tunit,extnam,status)
!
!*******************************************************************************
!
!! FTPTBH is obsolete.  Call FTPHTB instead.
!

    integer ounit,ncols,nrows,nfield,tbcol(*),status
    character ( len = * ) ttype(*),tform(*),tunit(*),extnam

    call ftphtb(ounit,ncols,nrows,nfield,ttype,tbcol, &
    tform,tunit,extnam,status)
end
subroutine ftptbs(iunit,frow,fchar,nchars,svalue,status)
!
!*******************************************************************************
!
!! FTPTBS writes a string of characters to an ascii or binary table.
!
!       This will span multiple rows of the table if NCHARS+FCHAR is
!       greater than the length of a row.

!       iunit   i  fortran unit number
!       frow    i  starting row number (1st row = 1)
!       fchar   i  starting character/byte in the row to write (1st character=1)
!       nchars  i  number of characters/bytes to write (can span multiple rows)
!       svalue  c  string of characters to write
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1991

    integer iunit,frow,fchar,nchars,status
    character ( len = * ) svalue

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff,bstart

    if (status > 0)return

    ibuff=bufnum(iunit)

!       check for errors
    if (nchars <= 0)then
!               zero or negative number of character requested
            return
    else if (frow < 1)then
!               error: illegal first row number
            status=307
            return
    else if (fchar < 1)then
!               error: illegal starting character
            status=308
            return
    end if

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(iunit,status)

!       move the i/o pointer to the start of the sequence of characters
    bstart=dtstrt(ibuff)+(frow-1)*rowlen(ibuff)+fchar-1
    call ftmbyt(iunit,bstart,.true.,status)

!       put the string of characters
    call ftpcbf(iunit,nchars,svalue,status)
end
subroutine ftptdm(iunit,colnum,naxis,naxes,status)
!
!*******************************************************************************
!
!! FTPTDM writes the TDIMnnn keyword describing the dimensionality of a column.
!
!       iunit   i  fortran unit number to use for reading
!       colnum  i  column number to read
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       OUTPUT PARAMETERS:
!       status  i  output error status (0=OK)
!
!       written by Wm Pence, HEASARC/GSFC, October 1993

    integer iunit,colnum,naxis,naxes(*),status

    integer i,j,nextsp
    character tdim*120, cval*20

    if (status > 0)return

    if (naxis < 1 .or. naxis > 100)then
!               illegal number of axes
            status=320
            return
    else if (colnum < 1 .or. colnum > 999)then
!               illegal column number
            status=302
            return
    end if

!       construct the keyword value
    tdim='('

    nextsp=2
    do 100 i=1,naxis
            if (naxes(i) < 1)then
                    status=323
                    return
            end if

!               convert integer to right justified C*20 string
            call fti2c(naxes(i),cval,status)
            if (status > 0)return

            do 20 j=20,1,-1
                    if (cval(j:j) == ' ')then
                            tdim(nextsp:)=cval(j+1:20)
                            nextsp=nextsp+21-j
                            tdim(nextsp-1:)=','
                            go to 100
                    end if
20              continue
100     continue

    tdim(nextsp-1:)=')'

    call ftpkns(iunit,'TDIM',colnum,1,tdim, &
            'size of the multidimensional array',status)
end
subroutine ftpthp(ounit,heap,status)
!
!*******************************************************************************
!
!! FTPTHP defines the starting address for the heap for a binary table.
!
!       The default address is NAXIS1 * NAXIS2.  It is in units of
!       bytes relative to the beginning of the regular binary table data.
!       This subroutine also writes the appropriate THEAP keyword to the
!       FITS header.

!       ounit   i  Fortran I/O unit number
!       heap   i  starting address of the heap
!       OUTPUT PARAMETERS:
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, Nov 1991

    integer ounit,heap,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff

    if (status > 0)return
    ibuff=bufnum(ounit)
    theap(ibuff)=heap

!       write the keyword
    call ftukyj(ounit,'THEAP',heap,'Byte offset of heap area', &
                status)
end
subroutine ftpunt(ounit,keywrd,kunit,status)
!
!*******************************************************************************
!
!! FTPUNT writes the units string in a header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       kunit   c  keyword units string
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1997

    character ( len = * ) keywrd,kunit
    integer ounit,status,lenval,ii,clen,olen
    character card*80,value*80,knam*8,ocmnt*72,ncmnt*72

    if (status > 0)return

    knam=keywrd

!       find the old keyword
    call ftgcrd(ounit,knam,card,status)
    if (status == 202)then
      call ftpmsg('FTPUNT Could not find the '//knam//' keyword.')
      return
    end if

!       parse the record to find value and comment strings
    call ftpsvc(card,value,ocmnt,status)

!       get the length of the keyword name + value string
    call ftprsv(card,lenval,status)

    if (status > 0)return

!       write the units string, in square brackets, to the new comment

    clen=1
    if (kunit /= ' ')then
      ncmnt='['//kunit

      do ii = 72,1,-1
          if (ncmnt(ii:ii) /= ' ')then
                clen = ii+1
                ncmnt(clen:)='] '
                clen=clen+2
                exit
          end if
      end do
    end if

!       check for existing units field in the comment
    olen=1
    if (ocmnt(1:1) == '[')then
        do 30 ii = 2,72
            if (ocmnt(ii:ii) == ']')then
                olen=ii+1
                if (ocmnt(olen:olen) == ' ')olen=olen+1
                go to 40
            end if
30          continue
    end if
40      continue

!       concatinate the old comment string to the new string
    ncmnt(clen:)=ocmnt(olen:)

!       construct the whole new card
    card(lenval+1:)=' / '//ncmnt

!       modify the keyword record
    call ftmodr(ounit,card,status)
end
subroutine ftr2e(val,dec,cval,status)
!
!*******************************************************************************
!
!! FTR2E converts real value to E20.* format character string.
!
!       val     r  input value to be converted
!       dec     i  number of decimal places to display in output string
!       cval    c  output character string
!       status  i  output error status (0 = OK)

    real val
    integer dec,status
    character*20 cval,form*10

    if (status > 0)return

    if (dec >= 1 .and. dec <= 9)then
            write(form,2000)dec
2000            format('(1pe20.',i1,')')
    else if (dec >= 10 .and. dec <= 13)then
            write(form,2001)dec
2001            format('(1pe20.',i2,')')
    else
!               illegal number of decimal places were specified
            status=411
            call ftpmsg('Error in FTR2E: number of decimal places ' &
                        //'is less than 1 or greater than 13.')
            return
    end if

    write(cval,form,err=900)val
    if (cval(1:1) == '*')go to 900
    return

900     status=402
    call ftpmsg('Error in FTR2E converting real to E20. string.')
end
subroutine ftr2f(val,dec,cval,status)
!
!*******************************************************************************
!
!! FTR2F converts real value to F20.* format character string.
!
!       val     r  input value to be converted
!       dec     i  number of decimal places to display in output string
!       cval    c  output character string
!       status  i  output error status (0 = OK)

    real val
    integer dec,status
    character*20 cval,form*8

    if (status > 0)return

    if (dec >= 0 .and. dec <= 9)then
            write(form,2000)dec
2000            format('(f20.',i1,')')
    else if (dec >= 10 .and. dec <18)then
            write(form,2001)dec
2001            format('(f20.',i2,')')
    else
            status=411
            call ftpmsg('Error in FTR2F: number of decimal places ' &
                        //'is less than 0 or greater than 18.')
            return
    end if

    write(cval,form,err=900)val
    if (cval(1:1) == '*')go to 900
    return
900     status=402
    call ftpmsg('Error in FTR2F converting real to F20. string.')
end
subroutine ftr4i1(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR4I1 copies input r*4 values to output i*1 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   r input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  c*1 value to set  array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  c*1 returned array of values
!       status  i  output error status (0 = ok)

    real input(*)
    character output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval
    logical tofits,flgray(*),anynul,noscal
    logical fttrnn
    external fttrnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                do 10 i=1,n
!                       trap any values that overflow the I*1 range
                    if (input(i)< 255.49 .and. &
                        input(i)> -.49)then
                            output(i)=char(nint(input(i)))
                    else if (input(i) >= 255.49)then
                            status=-11
                            output(i)=char(255)
                    else
                            status=-11
                            output(i)=char(0)
                    end if
10                  continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                            output(i)=char(nint(dval))
                        else if (dval >= 255.49)then
                            status=-11
                            output(i)=char(255)
                        else
                            status=-11
                            output(i)=char(0)
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                   don't have to check for nulls
                if (noscal)then
                  do 30 i=1,n
!                       trap any values that overflow the I*1 range
                    if (input(i)< 255.49 .and. &
                        input(i)> -.49)then
                            output(i)=char(int(input(i)))
                    else if (input(i) >= 255.49)then
                            status=-11
                            output(i)=char(255)
                    else
                            status=-11
                            output(i)=char(0)
                    end if
30                    continue
                else
                    do 40 i=1,n
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
40                      continue
                end if
            else
!                   must test for null values
                if (noscal)then
                     do 50 i=1,n
                         if (fttrnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                          else
!                               trap any values that overflow the I*1 range
                            if (input(i)< 255.49 .and. &
                               input(i)> -.49)then
                               output(i)=char(int(input(i)))
                            else if (input(i) >= 255.49)then
                                status=-11
                                output(i)=char(255)
                            else
                                status=-11
                                output(i)=char(0)
                            end if
                          end if
50                       continue
                else
                   do 60 i=1,n
                      if (fttrnn(input(i)))then
                                anynul=.true.
                                if (chktyp == 1)then
                                    output(i)=setval
                                else
                                    flgray(i)=.true.
                                end if
                      else
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
                      end if
60                     continue
                end if
            end if
    end if
end
subroutine ftr4i2(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR4I2 copies input r*4 values to output i*2 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   r  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  i*2 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i*2 returned array of values
!       status  i  output error status (0 = ok)

    real input(*)
    integer*2 output(*),setval,mmini2,mmaxi2
    integer n,i,chktyp,status
    double precision scale,zero,dval,i2max,i2min
    logical tofits,flgray(*),anynul,noscal
    logical fttrnn
    parameter (i2max=3.276749D+04)
    parameter (i2min=-3.276849D+04)
    real mini2,maxi2
    parameter (maxi2=32767.49)
    parameter (mini2=-32768.49)
    parameter (mmaxi2=32767)
    parameter (mmini2=-32768)
    external fttrnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
!                           trap any values that overflow the I*2 range
                        if (input(i) <= maxi2 .and. &
                            input(i) >= mini2)then
                                output(i)=nint(input(i))
                        else if (input(i) > maxi2)then
                                status=-11
                                output(i)=mmaxi2
                        else
                                status=-11
                                output(i)=mmini2
                        end if
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i2max .and. dval>i2min)then
                            output(i)=nint(dval)
                        else if (dval >= i2max)then
                            status=-11
                            output(i)=mmaxi2
                        else
                            status=-11
                            output(i)=mmini2
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                      do 30 i=1,n
!                           trap any values that overflow the I*2 range
                        if (input(i) <= maxi2 .and. &
                            input(i) >= mini2)then
                                output(i)=int(input(i))
                        else if (input(i) > maxi2)then
                                status=-11
                                output(i)=mmaxi2
                        else
                                status=-11
                                output(i)=mmini2
                        end if
30                        continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=int(dval)
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=mmaxi2
                          else
                              status=-11
                              output(i)=mmini2
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                    do 50 i=1,n
                        if (fttrnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
!                               trap any values that overflow the I*2 range
                            if (input(i) <= maxi2 .and. &
                                input(i) >= mini2)then
                                    output(i)=int(input(i))
                            else if (input(i) > maxi2)then
                                    status=-11
                                    output(i)=mmaxi2
                            else
                                    status=-11
                                    output(i)=mmini2
                            end if
                        end if
50                      continue
                else
                    do 60 i=1,n
                        if (fttrnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=int(dval)
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=mmaxi2
                          else
                              status=-11
                              output(i)=mmini2
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine ftr4i4(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR4I4 copies input r*4 values to output i*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   r  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  i   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i   returned array of values
!       status  i  output error status (0 = ok)

    real input(*)
    integer output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval,i4min,i4max
    logical tofits,flgray(*),anynul,noscal
    logical fttrnn
    parameter (i4max= 2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    real mini4,maxi4
!       Warning: only have about 7 digits of precision, so don't try
!       to set the maxi4 and mini4 limits any closer to the I*4 range.
    parameter (maxi4= 2.1474835E+09)
    parameter (mini4=-2.1474835E+09)
    integer mmaxi4,mmini4
    parameter (mmaxi4=2147483647)
    external fttrnn
!       work around for bug in the DEC Alpha VMS compiler
    mmini4=-2147483647 - 1

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
!                           trap any values that overflow the I*4 range
                        if (input(i) <= maxi4 .and. &
                            input(i) >= mini4)then
                                output(i)=nint(input(i))
                        else if (input(i) > maxi4)then
                                status=-11
                                output(i)=mmaxi4
                        else
                                status=-11
                                output(i)=mmini4
                        end if
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*4 range
                        if (dval<i4max .and. dval>i4min)then
                            output(i)=nint(dval)
                        else if (dval >= i4max)then
                            status=-11
                            output(i)=mmaxi4
                        else
                            status=-11
                            output(i)=mmini4
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                      do 30 i=1,n
!                           trap any values that overflow the I*4 range
                        if (input(i) <= maxi4 .and. &
                            input(i) >= mini4)then
                                output(i)=int(input(i))
                        else if (input(i) > maxi4)then
                                status=-11
                                output(i)=mmaxi4
                        else
                                status=-11
                                output(i)=mmini4
                        end if
30                        continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=int(dval)
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=mmaxi4
                          else
                              status=-11
                              output(i)=mmini4
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                    do 50 i=1,n
                        if (fttrnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                    output(i)=setval
                            else
                                    flgray(i)=.true.
                            end if
                        else
!                               trap any values that overflow the I*4 range
                            if (input(i) <= maxi4 .and. &
                                input(i) >= mini4)then
                                    output(i)=int(input(i))
                            else if (input(i) > maxi4)then
                                    status=-11
                                    output(i)=mmaxi4
                            else
                                    status=-11
                                    output(i)=mmini4
                            end if
                        end if
50                      continue
                else
                    do 60 i=1,n
                        if (fttrnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                 output(i)=setval
                            else
                                 flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=int(dval)
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=mmaxi4
                          else
                              status=-11
                              output(i)=mmini4
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine ftr4r4(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR4R4 copies input r*4 values to output r*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!
!       input   r  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  r  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       output  r  returned array of values

    real input(*)
    real output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal
    logical fttrnn
    external fttrnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (fttrnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (fttrnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine ftr4r8(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR4R8 copies input r*4 values to output r*8 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   r  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  d  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       output  d  returned array of values

    real input(*)
    double precision output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal
    logical fttrnn
    external fttrnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (fttrnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (fttrnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine ftr8i1(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR8I1 copies input r*8 values to output i*1 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   d input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  c*1 value to set  array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  c*1 returned array of values
!       status  i  output error status (0 = ok)

    double precision input(*)
    character output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval
    logical tofits,flgray(*),anynul,noscal
    logical fttdnn
    external fttdnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                do 10 i=1,n
!                       trap any values that overflow the I*1 range
                    if (input(i)< 255.49 .and. &
                        input(i)> -.49)then
                            output(i)=char(nint(input(i)))
                    else if (input(i) >= 255.49)then
                            status=-11
                            output(i)=char(255)
                    else
                            status=-11
                            output(i)=char(0)
                    end if
10                  continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                            output(i)=char(nint(dval))
                        else if (dval >= 255.49)then
                            status=-11
                            output(i)=char(255)
                        else
                            status=-11
                            output(i)=char(0)
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                   don't have to check for nulls
                if (noscal)then
                  do 30 i=1,n
!                       trap any values that overflow the I*1 range
                    if (input(i)< 255.49 .and. &
                        input(i)> -.49)then
                            output(i)=char(int(input(i)))
                    else if (input(i) >= 255.49)then
                            status=-11
                            output(i)=char(255)
                    else
                            status=-11
                            output(i)=char(0)
                    end if
30                    continue
                else
                    do 40 i=1,n
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
40                      continue
                end if
            else
!                   must test for null values
                if (noscal)then
                     do 50 i=1,n
                         if (fttdnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                          else
!                               trap any values that overflow the I*1 range
                            if (input(i)< 255.49 .and. &
                               input(i)> -.49)then
                               output(i)=char(int(input(i)))
                            else if (input(i) >= 255.49)then
                                status=-11
                                output(i)=char(255)
                            else
                                status=-11
                                output(i)=char(0)
                            end if
                          end if
50                       continue
                 else
                    do 60 i=1,n
                      if (fttdnn(input(i)))then
                                anynul=.true.
                                if (chktyp == 1)then
                                    output(i)=setval
                                else
                                    flgray(i)=.true.
                                end if
                      else
                        dval=input(i)*scale+zero
!                           trap any values that overflow the I*1 range
                        if (dval< 255.49 .and. dval> -.49)then
                                output(i)=char(int(dval))
                        else if (dval >= 255.49)then
                                status=-11
                                output(i)=char(255)
                        else
                                status=-11
                                output(i)=char(0)
                        end if
                      end if
60                     continue
                end if
            end if
    end if
end
subroutine ftr8i2(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR8I2 copies input r*8 values to output i*2 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   d  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  i*2 value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i*2 returned array of values
!       status  i  output error status (0 = ok)

    double precision input(*)
    integer*2 output(*),setval,maxi2,mini2
    integer n,i,chktyp,status
    double precision scale,zero,dval,i2max,i2min
    logical tofits,flgray(*),anynul,noscal
    logical fttdnn
    parameter (i2max=3.276749D+04)
    parameter (i2min=-3.276849D+04)

    parameter (maxi2=32767)
    parameter (mini2=-32768)
    external fttdnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
!                           trap any values that overflow the I*2 range
                        if (input(i) <= i2max .and. &
                            input(i) >= i2min)then
                                output(i)=nint(input(i))
                        else if (input(i) > i2max)then
                                status=-11
                                output(i)=maxi2
                        else
                                status=-11
                                output(i)=mini2
                        end if
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*2 range
                        if (dval<i2max .and. dval>i2min)then
                            output(i)=nint(dval)
                        else if (dval >= i2max)then
                            status=-11
                            output(i)=maxi2
                        else
                            status=-11
                            output(i)=mini2
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                      do 30 i=1,n
!                           trap any values that overflow the I*2 range
                        if (input(i) <= i2max .and. &
                            input(i) >= i2min)then
                                output(i)=int(input(i))
                        else if (input(i) > i2max)then
                                status=-11
                                output(i)=maxi2
                        else
                                status=-11
                                output(i)=mini2
                        end if
30                        continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=int(dval)
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                    do 50 i=1,n
                        if (fttdnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                 output(i)=setval
                            else
                                 flgray(i)=.true.
                            end if
                        else
!                               trap any values that overflow the I*2 range
                            if (input(i) <= i2max .and. &
                                input(i) >= i2min)then
                                    output(i)=int(input(i))
                            else if (input(i) > i2max)then
                                    status=-11
                                    output(i)=maxi2
                            else
                                    status=-11
                                    output(i)=mini2
                            end if
                        end if
50                      continue
                    else
                      do 60 i=1,n
                        if (fttdnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*2 range
                          if (dval<i2max .and. dval>i2min)then
                              output(i)=int(dval)
                          else if (dval >= i2max)then
                              status=-11
                              output(i)=maxi2
                          else
                              status=-11
                              output(i)=mini2
                          end if
                        end if
60                        continue
                    end if
            end if
    end if
end
subroutine ftr8i4(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR8I4 copies input r*8 values to output i*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   d  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  i   value to set output array to if value is undefined
!       flgray  l   array of logicals indicating if corresponding value is null
!       anynul  l   set to true if any nulls were set in the output array
!       output  i   returned array of values
!       status  i  output error status (0 = ok)

    double precision input(*)
    integer output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero,dval,i4min,i4max
    logical tofits,flgray(*),anynul,noscal
    logical fttdnn
    parameter (i4max=2.14748364749D+09)
    parameter (i4min=-2.14748364849D+09)
    integer maxi4,mini4
    parameter (maxi4=2147483647)
    external fttdnn
!       work around for bug in the DEC Alpha VMS compiler
    mini4=-2147483647 - 1

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
!                           trap any values that overflow the I*4 range
                        if (input(i) <= i4max .and. &
                            input(i) >= i4min)then
                                output(i)=nint(input(i))
                        else if (input(i) > i4max)then
                                status=-11
                                output(i)=maxi4
                        else
                                status=-11
                                output(i)=mini4
                        end if
10                      continue
            else
                    do 20 i=1,n
                        dval=(input(i)-zero)/scale
!                           trap any values that overflow the I*4 range
                        if (dval<i4max .and. dval>i4min)then
                            output(i)=nint(dval)
                        else if (dval >= i4max)then
                            status=-11
                            output(i)=maxi4
                        else
                            status=-11
                            output(i)=mini4
                        end if
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                      do 30 i=1,n
!                           trap any values that overflow the I*4 range
                        if (input(i) <= i4max .and. &
                            input(i) >= i4min)then
                                output(i)=int(input(i))
                        else if (input(i) > i4max)then
                                status=-11
                                output(i)=maxi4
                        else
                                status=-11
                                output(i)=mini4
                        end if
30                        continue
                    else
                        do 40 i=1,n
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=int(dval)
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
40                          continue
                    end if
            else
!                   must test for null values
                if (noscal)then
                    do 50 i=1,n
                        if (fttdnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                output(i)=setval
                            else
                                flgray(i)=.true.
                            end if
                        else
!                               trap any values that overflow the I*4 range
                            if (input(i) <= i4max .and. &
                                input(i) >= i4min)then
                                    output(i)=int(input(i))
                            else if (input(i) > i4max)then
                                    status=-11
                                    output(i)=maxi4
                            else
                                    status=-11
                                    output(i)=mini4
                            end if
                        end if
50                      continue
                else
                    do 60 i=1,n
                        if (fttdnn(input(i)))then
                            anynul=.true.
                            if (chktyp == 1)then
                                 output(i)=setval
                            else
                                 flgray(i)=.true.
                            end if
                        else
                          dval=input(i)*scale+zero
!                             trap any values that overflow the I*4 range
                          if (dval<i4max .and. dval>i4min)then
                              output(i)=int(dval)
                          else if (dval >= i4max)then
                              status=-11
                              output(i)=maxi4
                          else
                              status=-11
                              output(i)=mini4
                          end if
                        end if
60                      continue
                end if
            end if
    end if
end
subroutine ftr8r4(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR8R4 copies input r*8 values to output r*4 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   d  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  r  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       output  r  returned array of values

    double precision input(*)
    real output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal
    logical fttdnn
    external fttdnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (fttdnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (fttdnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine ftr8r8(input,n,scale,zero,tofits, &
            chktyp,setval,flgray,anynul,output,status)
!
!*******************************************************************************
!
!! FTR8R8 copies input r*8 values to output r*8 values.
!
!  The routine does optional scaling and checking for null values.
!
!       input   d  input array of values
!       n       i  number of values
!       scale   d  scaling factor to be applied
!       zero    d  scaling zero point to be applied
!       tofits  l  true if converting from internal format to FITS
!       chktyp  i  type of null value checking to be done if TOFITS=.false.
!                       =0  no checking for null values
!                       =1  set null values = SETVAL
!                       =2  set corresponding FLGRAY value = .true.
!       setval  d  value to set output array to if value is undefined
!       flgray  l  array of logicals indicating if corresponding value is null
!       anynul  l  set to true if any nulls were set in the output array
!       output  d  returned array of values

    double precision input(*)
    double precision output(*),setval
    integer n,i,chktyp,status
    double precision scale,zero
    logical tofits,flgray(*),anynul,noscal
    logical fttdnn
    external fttdnn

    if (status > 0)return

    if (scale == 1. .and. zero == 0)then
            noscal=.true.
    else
            noscal=.false.
    end if

    if (tofits) then
!               we don't have to worry about null values when writing to FITS
            if (noscal)then
                    do 10 i=1,n
                            output(i)=input(i)
10                      continue
            else
                    do 20 i=1,n
                            output(i)=(input(i)-zero)/scale
20                      continue
            end if
    else
!               converting from FITS to internal format; may have to check nulls
            if (chktyp == 0)then
!                       don't have to check for nulls
                    if (noscal)then
                            do 30 i=1,n
                                    output(i)=input(i)
30                              continue
                    else
                            do 40 i=1,n
                                    output(i)=input(i)*scale+zero
40                              continue
                    end if
            else
!                       must test for null values
                    if (noscal)then
                            do 50 i=1,n
                                    if (fttdnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                            output(i)=input(i)
                                    end if
50                              continue
                    else
                            do 60 i=1,n
                                    if (fttdnn(input(i)))then
                                        anynul=.true.
                                        if (chktyp == 1)then
                                            output(i)=setval
                                        else
                                            flgray(i)=.true.
                                        end if
                                    else
                                     output(i)=input(i)*scale+zero
                                    end if
60                              continue
                    end if
            end if
    end if
end
subroutine ftrdef(ounit,status)
!
!*******************************************************************************
!
!! FTRDEF rereads the CHDU header keywords to determine the structure.
!
!
!       ReDEFine the structure of a data unit.  This routine re-reads
!       the CHDU header keywords to determine the structure and length of the
!       current data unit.  This redefines the start of the next HDU.
!
!       ounit   i  Fortran I/O unit number
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1993

    integer ounit,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,dummy

    if (status > 0)return

    ibuff=bufnum(ounit)

!       see if we have write access to this file (no need to go on, if not)
    if (wrmode(ibuff))then
!           rewrite the header END card, and following blank fill
        call ftwend(ounit,status)
        if (status > 0)return

!           now re-read the required keywords to determine the structure
        call ftrhdu(ounit,dummy,status)
    end if
end
subroutine ftrhdu(iunit,xtend,status)
!
!*******************************************************************************
!
!! FTRHDU reads the CHDU structure.
!
!  It does this by reading the header keywords which define
!  the size and structure of the header and data units.
!
!       iunit   i  Fortran I/O unit number
!       OUTPUT PARAMETERS:
!       xtend   i  returned type of extension:   0 = the primary HDU
!                                                1 = an ASCII table
!                                                2 = a binary table
!                                               -1 = unknown
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer iunit,xtend,status,i,ic,tstat
    character keynam*8,exttyp*10,comm*30,keybuf*80
    logical endof

    if (status > 0)return

!       read first keyword to determine the type of the CHDU
    call ftgrec(iunit,1,keybuf,status)

    if (status > 0)then
      call ftpmsg('Cannot read first keyword in header (FTRHDU)')
            return
    end if

!       release any current column descriptors for this unit
    call ftfrcl(iunit,status)

    keynam=keybuf(1:8)
!       parse the value and comment fields from the record
    call ftpsvc(keybuf,exttyp,comm,status)

    if (status > 0)then
!               unknown type of FITS record; can't read it
      call ftpmsg('Cannot parse value of first keyword; unknown ' &
        //'type of FITS record (FTRHDU):')

    else if (keynam == 'SIMPLE')then
!               initialize the parameters describing the primay HDU
            call ftpini(iunit,status)
            xtend=0
    else if (keynam=='XTENSION')then
            if (exttyp(1:1) /= '''')then
!                       value of XTENSION is not a quoted character string!
                    if (keybuf(9:10) /= '= ')then
                        call ftpmsg('XTENSION keyword does not ' &
                       //'have "= " in cols 9-10.')
                    else
                    call ftpmsg('Unknown type of extension; value' &
                 //' of XTENSION keyword is not a quoted string:')
                    end if
                    status=251
                    call ftpmsg(keybuf)
            else if (exttyp(2:9) == 'TABLE   ')then
!                       initialize the parameters for the ASCII table extension
                    call ftaini(iunit,status)
                    xtend=1
            else if (exttyp(2:9) == 'BINTABLE' .or. exttyp(2:9) &
              == 'A3DTABLE' .or. exttyp(2:9) == '3DTABLE ')then
!                       initialize the parameters for the binary table extension
                    call ftbini(iunit,status)
                    xtend=2
            else
!                       try to initialize the parameters describing extension
                    tstat=status
                    call ftpini(iunit,status)
                    xtend=0
                    if (status == 251)then
!                           unknown type of extension
                        xtend=-1
                        status=tstat
                    end if
            end if
    else
!               unknown record
!               If file is created on a VAX with 512-byte records, then
!               the FITS file may have fill bytes (ASCII NULs) at the end.
!               Also, if file has been editted on a SUN, an extra ASCII 10
!               character may appear at the end of the file.  Finally, if
!               file is not a multiple of the record length long, then
!               the last truncated record may be filled with ASCII blanks.
!               So, if the record only contains NULS, LF, and blanks, then
!               assume we found the end of file.  Otherwise report an error.

        endof=.true.
            do 10 i=1,80
                ic=ichar(keybuf(i:i))
                if (ic /= 0 .and. ic /= 10 .and. ic /= 32) &
               endof=.false.
10              continue
            if (endof)then
                 status=107
                 call ftpmsg('ASCII 0s, 10s, or 32s at start of ' &
               //'extension are treated as EOF (FTRHDU):')
            else
                 status=252
                 call ftpmsg('Extension does not start with SIMPLE' &
                 //' or XTENSION keyword (FTRHDU):')
            end if
            xtend=-1
            call ftpmsg(keybuf)
    end if
end
subroutine ftrsim(ounit,bitpix,naxis,naxes,status)
!
!*******************************************************************************
!
!! FTRSIM resizes an existing primary array or IMAGE extension.
!
!       ounit   i  fortran output unit number
!       bitpix  i  number of bits per data value
!       naxis   i  number of axes in the data array
!       naxes   i  array giving the length of each data axis
!       status  i  returned error status (0=ok)

!       written by Wm Pence, HEASARC/GSFC, July 1997

    integer ounit,bitpix,naxis,naxes(*),status
    integer i,bytlen,nblock,minax
    integer nsize,osize,obitpx,onaxis,onaxes(99),pcount,gcount
    logical simple,extend
    character*8 keynm

    if (status > 0)return

    call ftghpr(ounit,99,simple,obitpx,onaxis,onaxes, &
                      pcount,gcount,extend,status)
    if (status > 0)return

!       check for error conditions
    if (naxis < 0 .or. naxis > 999)then
            status=212
           return
    end if

!       test that bitpix has a legal value and set the datatype code value
5       if (bitpix == 8)then
            bytlen=1
    else if (bitpix == 16)then
            bytlen=2
    else if (bitpix == 32)then
            bytlen=4
    else if (bitpix == -32)then
            bytlen=4
    else if (bitpix == -64)then
            bytlen=8
    else
!               illegal value of bitpix
            status=211
            return
    end if

!       calculate the number of pixels in the new image
    if (naxis == 0)then
!               no data
            nsize=0
    else
            nsize=1
            do 10 i=1,naxis
                    if (naxes(i) >= 0)then
                            nsize=nsize*naxes(i)
                    else
                            status=213
                            return
                    end if
10              continue
    end if

!       calculate the number of pixels in the old image
    if (onaxis == 0)then
!               no data
            osize=0
    else
            osize=1
            do i=1,onaxis
                    if (onaxes(i) >= 0)then
                            osize=osize*onaxes(i)
                    else
                            status=213
                            return
                    end if
            end do
    end if

!       sizes of old and new images, in bytes
    osize=(osize+pcount) * gcount * abs(obitpx)/8
    nsize=(nsize+pcount) * gcount * bytlen

!       sizes of old and new images, in blocks
    osize=(osize+2879)/2880
    nsize=(nsize+2879)/2880

!       insert or delete blocks, as necessary
    if (nsize > osize)then
         nblock=nsize-osize
         call ftiblk(ounit,nblock,1,status)
    else if (osize > nsize)then
         nblock=osize-nsize
         call ftdblk(ounit,nblock,1,status)
    end if
    if (status > 0)return

!       update the header keywords

    if (bitpix /= obitpx)then
        call ftmkyj(ounit,'BITPIX',bitpix,'&',status)
    end if

    if (naxis /= onaxis)then
        call ftmkyj(ounit,'NAXIS',naxis,'&',status)
    end if

!       update all the existing keywords
    minax=min(naxis,onaxis)
    do 20 i=1,minax
        call ftkeyn('NAXIS',i,keynm,status)
        call ftmkyj(ounit,keynm,naxes(i),'&',status)
20      continue

    if (naxis > onaxis)then
!           insert more NAXISn keywords
        do 25 i=onaxis+1,naxis
            call ftkeyn('NAXIS',i,keynm,status)
            call ftikyj(ounit,keynm,naxes(i), &
                        'length of data axis',status)
25          continue
    else if (onaxis > naxis)then
!           delete old NAXISn keywords
        do 30 i=naxis+1,onaxis
            call ftkeyn('NAXIS',i,keynm,status)
            call ftdkey(ounit,keynm,status)
30          continue
    end if

!       re-read the header, to make sure structures are updated
    call ftrdef(ounit,status)
end
subroutine ftrsnm
!
!*******************************************************************************
!
!! FTRSNM resets the column names as undefined.
!
!       this will force ftgcnn to read the column names from the
!       file the next time it is called

!       written by Wm Pence, HEASARC/GSFC, Feb 1995

    integer colpnt,untpnt
    common/ftname/colpnt,untpnt

    colpnt= -999
    untpnt=0
end
subroutine ftrwdn(iunit,frow,lrow,nshift,status)
!
!*******************************************************************************
!
!! FTRWDN shifts rows in a table down by NROWS rows, inserting blank rows.
!
!       iunit   i  Fortran I/O unit number
!       frow    i  rows *AFTER* this one are to be moved down
!       lrow    i  last row to be moved down (last row of the table)
!       nshift  i  how far to shift the rows
!       status  i  returned error status (0=ok)

    integer iunit,frow,lrow,nshift,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*5760 buff(2)
    character xdummy(20480)
    common/ftheap/buff,xdummy
!

    integer ibuff,kshift,nchar,fchar,in,out,i,j,irow,tin,jrow
    integer lstptr,inptr,outptr,nseg
    character cfill*1

    if (status > 0)return

!       don't have to do anything if inserting blank rows at end of the table
    if (frow == lrow)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       select appropriate fill value
    if (hdutyp(ibuff) == 1)then
!           fill  header or ASCII table with space
        cfill=char(32)
    else
!           fill image or bintable data area with Null (0)
        cfill=char(0)
    end if

!       how many rows will fit in the single buffer?
    kshift=2880/rowlen(ibuff)

!
!       CASE #1: optimal case where the NSHIFT number of rows will all
!       fit in the 2880-byte work buffer simultaneously.  The rows can
!       be shifted down in one efficient pass through the table.
!
    if (kshift >= nshift)then

!    Note: the f77 compiler with the -O flag on a linux PC gives
!    incorrect results with the following 2 lines:
!       kshift=nshift
!       nchar=kshift*rowlen(ibuff)
!    Apparently the compiler simply ignores the first statement
!    so kshift is left with it's old value when multipying times rowlen

    nchar=nshift*rowlen(ibuff)
    fchar=1

!       initialize the first buffer
    in=2
    out=1

    do 5 i=1,2880
        buff(1)(i:i)=cfill
5       continue

    do 10 irow=frow+1,lrow,nshift

!           read the row(s) to be shifted
        call ftgtbs(iunit,irow,fchar,nchar,buff(in),status)

!           overwrite these row(s) with the previous row(s)
        call ftptbs(iunit,irow,fchar,nchar,buff(out),status)

!           swap the input and output buffer pointers and move to next rows
        tin=in
        in=out
        out=tin
        jrow=irow
10      continue

!       write the last row(s) out
    irow=jrow+nshift
    nchar=(lrow-jrow+1)*rowlen(ibuff)

    call ftptbs(iunit,irow,fchar,nchar,buff(out),status)
    return

!
!       CASE #2: One or more rows of the table will fit in the work buffer,
!       but cannot fit all NSHIFT rows in the buffer at once.  Note that
!       since we do not need 2 buffers, as in the previous case, we can
!       combine both buffers into one single 2880*2 byte buffer, to handle
!       wider tables.  This algorithm copies then moves blocks of contiguous
!       rows at one time, working upwards from the bottom of the table.
!
    else if (rowlen(ibuff) <= 5760)then

!       how many rows can we move at one time?
    kshift=5760/rowlen(ibuff)
    fchar=1

!       initialize pointers
    lstptr=lrow
    inptr=lrow-kshift+1

20      if (inptr <= frow)inptr=frow+1
    nchar=(lstptr-inptr+1)*rowlen(ibuff)
    outptr=inptr+nshift

!       read the row(s) to be shifted
    call ftgtbs(iunit,inptr,fchar,nchar,buff,status)

!       write the row(s) to the new location
    call ftptbs(iunit,outptr,fchar,nchar,buff,status)

!       If there are more rows, update pointers and repeat
    if (inptr > frow+1)then
        lstptr=lstptr-kshift
        inptr =inptr -kshift
        go to 20
    end if

!       initialize the buffer with the fill value
    do 25 i=1,5760
        buff(1)(i:i)=cfill
25      continue

!       fill the empty rows with blanks or nulls
    nchar=rowlen(ibuff)
    do 30 i=1,nshift
        outptr=frow+i
        call ftptbs(iunit,outptr,fchar,nchar,buff,status)
30      continue
    return

!
!       CASE #3:  Cannot fit a whole row into the work buffer, so have
!       to move each row in pieces.
!
    else

    nseg=(rowlen(ibuff)+5759)/5760
    nchar=5760

    do 60 j=1,nseg
        fchar=(j-1)*5760+1
        if (j == nseg)nchar=rowlen(ibuff)-(nseg-1)*5760

        do 40 i=lrow,frow+1,-1
!               read the row to be shifted
            call ftgtbs(iunit,i,fchar,nchar,buff,status)

!               write the row(s) to the new location
            call ftptbs(iunit,i+nshift,fchar,nchar,buff,status)
40          continue

!           initialize the buffer with the fill value
        do 45 i=1,5760
            buff(1)(i:i)=cfill
45          continue

!           fill the empty rows with blanks or nulls
        do 50 i=1,nshift
            outptr=frow+i
            call ftptbs(iunit,outptr,fchar,nchar,buff,status)
50          continue
60      continue

    end if
end
subroutine ftrwup(iunit,frow,lrow,nshift,status)
!
!*******************************************************************************
!
!! FTRWUP shifts rows in a table up by NROWS rows, overwriting the rows above.
!
!       iunit   i  Fortran I/O unit number
!       frow    i  first row to be moved up
!       lrow    i  last row to be moved up (last row of the table)
!       nshift  i  how far to shift the rows (number of rows)
!       status  i  returned error status (0=ok)

    integer iunit,frow,lrow,nshift,status

!
    integer nf,nb,ne
    parameter (nb = 20)
    parameter (nf = 3000)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character*5760 buff
    character xdummy(26240)
    common/ftheap/buff,xdummy
!

    integer ibuff,kshift,nchar,fchar,i,j
    integer lstptr,inptr,outptr,nseg
    character cfill*1

    if (status > 0)return

!       define the number of the buffer used for this file
    ibuff=bufnum(iunit)

!       select appropriate fill value
    if (hdutyp(ibuff) == 1)then
!           fill  header or ASCII table with space
        cfill=char(32)
    else
!           fill image or bintable data area with Null (0)
        cfill=char(0)
    end if

!
!       CASE #1: One or more rows of the table will fit in the work buffer,
!
    if (rowlen(ibuff) <= 5760)then

!       how many rows can we move at one time?
    kshift=5760/rowlen(ibuff)
    fchar=1

!       check if we just need to clear the last NSHIFT rows of the table
    if (frow == lrow+1)go to 25

!       initialize pointers
    inptr=frow
    lstptr=inptr+kshift-1

20      if (lstptr > lrow)lstptr=lrow
    nchar=(lstptr-inptr+1)*rowlen(ibuff)
    outptr=inptr-nshift

!       read the row(s) to be shifted
    call ftgtbs(iunit,inptr,fchar,nchar,buff,status)

!       write the row(s) to the new location
    call ftptbs(iunit,outptr,fchar,nchar,buff,status)

!       If there are more rows, update pointers and repeat
    if (lstptr < lrow)then
        inptr =inptr +kshift
        lstptr=lstptr+kshift
        go to 20
    end if

!       initialize the buffer with the fill value
25      continue
    do 30 i=1,5760
        buff(i:i)=cfill
30      continue

!       fill the empty rows at the bottom of the table with blanks or nulls
    nchar=rowlen(ibuff)
    do 35 i=1,nshift
        outptr=lrow-nshift+i
        call ftptbs(iunit,outptr,fchar,nchar,buff,status)
35      continue
    return

!
!       CASE #2:  Cannot fit a whole row into the work buffer, so have
!       to move each row in pieces.
!
    else

    nseg=(rowlen(ibuff)+5759)/5760
    nchar=5760

    do 60 j=1,nseg
        fchar=(j-1)*5760+1
        if (j == nseg)nchar=rowlen(ibuff)-(nseg-1)*5760

!           check if we just need to clear the last NSHIFT rows of the table
        if (frow == lrow+1)go to 45

        do 40 i=frow,lrow
!               read the row to be shifted
            call ftgtbs(iunit,i,fchar,nchar,buff,status)

!               write the row(s) to the new location
            call ftptbs(iunit,i-nshift,fchar,nchar,buff,status)
40          continue

!           initialize the buffer with the fill value
45          continue
        do 50 i=1,5760
            buff(i:i)=cfill
50          continue

!           fill the empty rows with blanks or nulls
        do 55 i=1,nshift
            outptr=lrow-nshift+i
            call ftptbs(iunit,outptr,fchar,nchar,buff,status)
55          continue
60      continue
    end if
end
subroutine fts2c(in,cval,lenval,status)
!
!*******************************************************************************
!
!! FTS2C converts an input string to a left justified quoted string.
!
!               The minimum length FITS string is 8 characters, so
!               pad the quoted string with spaces if necessary.
!       cval = returned quoted string
!       lenval = length of the cval string, including the 2 quote characters
    character ( len = * ) in,cval
    integer length,i,j,i1,i2,lenval,status

    if (status > 0)return

    i1=1
    i2=1
!       test for blank input string
    if (in == ' ')then
            cval='''        '''
            lenval=10
            return
    end if

    length=len(in)
!       find first and last non-blank characters

!       modified 29 Nov 1994 to treat leading spaces as significant
!        do 5 i=1,length
!                i1=i
!                if (in(i:i) /= ' ')go to 10
!5       continue
!10      continue

    do 15 i=length,1,-1
            i2=i
            if (in(i:i) /= ' ')go to 20
15      continue
20      continue

    cval=''''//in(i1:i2)

!       test if there are any single quotes in the string;  if so, replace
!       them with two successive single quotes
    lenval=i2-i1+2
    do 30 i=lenval,2,-1
            if (cval(i:i) == '''')then
!                  shift all the characters over 1 space
               do 40 j=len(cval),i+1,-1
                  cval(j:j)=cval(j-1:j-1)
40                 continue
               i2=i2+1
            end if
30      continue

!       find location of closing quote
    lenval=max(10,i2-i1+3)
    lenval=min(lenval,len(cval))
    if (lenval == 70 .and. cval(69:70) == '''''')then
!           this occurs if the string ends with a literal appostrophy
        cval(70:70) = ' '
    else
        cval(lenval:lenval)=''''
    end if
end
subroutine ftsdnn(value)
!
!*******************************************************************************
!
!! FTSDNN sets a 64-bit pattern equal to an IEEE Not-a-Number value.
!
!       A NaN has all the exponent bits=1, and the fractional part
!       not=0.
!
!       written by Wm Pence, HEASARC/GSFC, February 1991

    integer value(2)

!       there are many NaN values;  choose a simple one in which all bits=1
    value(1)=-1
    value(2)=-1
end
subroutine ftsnul(ounit,colnum,nulval,status)
!
!*******************************************************************************
!
!! FTSNUL defines the null value for an ASCII table column.
!
!       ounit   i  Fortran I/O unit number
!       colnum  i  number of the column to be defined
!       nulval  c  the string to be use to signify undefined data
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,status
    character ( len = * ) nulval

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
    character cnull*16, cform*8
    common/ft0003/cnull(nf),cform(nf)
!

    integer ibuff

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)
    if (status > 0)return

!       test for proper HDU type
    if (hdutyp(ibuff) /= 1)then
        status=226
        return
    end if

    if (colnum > tfield(ibuff) .or. colnum < 1)then
         status=302
         return
    end if

    cnull(colnum+tstart(ibuff))=nulval
end
subroutine ftsrnn(value)
!
!*******************************************************************************
!
!! FTSRNN sets a 32-bit pattern equal to an IEEE Not-a-Number value.
!
!       A NaN has all the exponent bits=1, and the fractional part
!       not=0.
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer value

!       there are many NaN values;  choose a simple one in which all bits=1
    value=-1
end
subroutine fttbit(bitpix,status)
!
!*******************************************************************************
!
!! FTTBIT tests that bitpix has a legal value.
!
    integer bitpix,status
    character value*20

    if (status > 0)return

    if (bitpix /= 8 .and. bitpix /= 16 .and. bitpix /= 32 &
        .and. bitpix /= -32 .and. bitpix /= -64)then
            status=211
            write(value,1000)bitpix
1000            format(i20)
            call ftpmsg('Illegal BITPIX value: '//value)
    end if
end
    function fttdnn(value)
!
!*******************************************************************************
!
!! FTTDNN tests if a R*8 value has a IEEE Not-a-Number value.
!
!       A NaN has all the exponent bits=1, and the fractional part
!       not=0.
!       Exponent field is in bits 20-30 in the most significant 4-byte word
!       Mantissa field is in bits 0-19 of most sig. word and entire 2nd word
!
!       written by Wm Pence, HEASARC/GSFC, May 1992
!       modified Aug 1994 to handle all IEEE special values.

    integer value(2)
!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    logical fttdnn
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer compid
    common/ftcpid/compid
!
    integer word1,word2
!       COMPID specifies what type of floating point word structure
!       is used on this machine, and determines how to test for NaNs.

!       COMPID value:
!          2 or 3    generic machine: simply test for NaNs with all bits set
!           1        like a decstation or alpha OSF/1, or IBM PC
!           0        SUN workstation, or IBM mainframe
!          -2305843009213693952   Cray (64-bit) machine

    fttdnn=.false.
    if (compid >= 2)then
!           on the VAX we can assume that all NaNs will be set to all bits on
!           (which is equivalent to an integer with a value of -1) because
!           this is what the IEEE to VAX conversion MACRO program returns
        if (value(1) == -1 .and. value(2) == -1)fttdnn=.true.
    else if (compid >= -1)then
        if (compid <= 0)then
!               this is for SUN-like machines, or IBM main frames
            word1=value(1)
            word2=value(2)
        else
!               this is for DECstation and IBM PCs.  The 2 32 bit integer words
!               are reversed from what you get on the SUN.
            word1=value(2)
            word2=value(1)
        end if

!           efficiently search the number space for NaNs and underflows
        if (word2 == -1)then
            if ((word1 >= -1048577 .and. word1 <= -1) &
             .or. (word1 >= 2146435071))then
                  fttdnn=.true.
            else if ((word1 < -2146435072) .or. &
            (word1 >= 0 .and. word1 < 1048576))then
                  value(1)=0
                  value(2)=0
            end if
         else if (word2 == 0)then
            if ((word1 > -1048577 .and. word1 <= -1) &
             .or. (word1 > 2146435071))then
                  fttdnn=.true.
            else if ((word1 <= -2146435072) .or. &
            (word1 >= 0 .and. word1 <= 1048576))then
                  value(1)=0
                  value(2)=0
            end if
         else
            if ((word1 > -1048577 .and. word1 <= -1) &
             .or. (word1 > 2146435071))then
                  fttdnn=.true.
            else if ((word1 < -2146435072) .or. &
            (word1 >= 0 .and. word1 < 1048576))then
                  value(1)=0
                  value(2)=0
            end if
         end if
    else
!           branch for the Cray:  COMPID stores the negative integer
!           which corresponds to the 3 most sig digits set to 1.   If these
!           3 bits are set in a floating point number, then it represents
!           a reserved value (i.e., a NaN)
        if (value(1)< 0 .and. value(1) >= compid)fttdnn=.true.
    end if
end
subroutine fttkey(keynam,status)
!
!*******************************************************************************
!
!! FTTKEY tests that keyword name contains only legal characters.
!
!  Legal characters are uppercase letters, numbers, hyphen, underscore, or space
!         (but no embedded spaces)

!       keynam  c*8  keyword name
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)

    character keynam*(*)
    integer status,i
    character c1,pos
    logical spaces

    if (status > 0)return

    spaces=.false.
    do 20 i=1,8
        c1=keynam(i:i)
        if ((c1 >= 'A' .and. c1 <= 'Z') .or. &
            (c1 >= '0' .and. c1 <= '9') .or. &
             c1 == '-' .or. c1 == '_')then
             if (spaces)then
!                   error: name contains embedded space
                status=207
                call ftpmsg('Keyword name contains embedded '// &
                'space(s): '//keynam(1:8))
                return
             end if
        else if (c1 == ' ')then
             spaces=.true.
        else
!                illegal character found
             status=207
             write(pos,1000)i
1000             format(i1)
             call ftpmsg('Character '//pos//' in this keyword name' &
             //' is illegal: "'//keynam(1:8)//'"')
!                explicitly test for the 2 most common cases:
             if (ichar(c1) == 0)then
               call ftpmsg('(This is an ASCII NUL (0) character).')
             else if (ichar(c1) == 9)then
               call ftpmsg('(This is an ASCII TAB (9) character).')
             end if
             return
        end if
20      continue
end
subroutine fttkyn(iunit,nkey,keynam,keyval,status)
!
!*******************************************************************************
!
!! FTTKYN tests that a keyword has a given name and value.
!
!  The routine tests that the keyword number NKEY has name = KEYNAM
!       and has value = KEYVAL
!
!       iunit   i  Fortran I/O unit number
!       nkey    i  sequence number of the keyword to test
!       keynam  c  name that the keyword is supposed to have
!       keyval  c  value that the keyword is supposed to have
!       OUTPUT PARAMETERS:
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!
    integer iunit,nkey,status
    character ( len = * ) keynam,keyval
    character kname*8,value*30,comm*48,npos*8,keybuf*80
    character errmsg*80

    if (status > 0)return

!       read the name and value of the keyword

!       get the whole record
    call ftgrec(iunit,nkey,keybuf,status)

    kname=keybuf(1:8)
!       parse the value and comment fields from the record
    call ftpsvc(keybuf,value,comm,status)
    if (status > 0)go to 900

!       test if the keyword has the correct name
    if (kname /= keynam)then
            status=208
            go to 900
    end if

!       check that the keyword has the correct value
    if (value /= keyval)then
            status=209
    end if

900     continue
    if (status > 0)then

        write(npos,1000)nkey
1000        format(i8)
        errmsg='FTTKYN found unexpected keyword or value '// &
        'for header keyword number '//npos//'.'
        call ftpmsg(errmsg)
        errmsg='  Was expecting keyword '//keynam// &
        ' with value = '//keyval
        call ftpmsg(errmsg)
        if (keybuf(9:10) /= '= ')then
          errmsg='      but found keyword '//kname// &
        ' with no "= " in cols. 9-10.'
        else
          errmsg='      but found keyword '//kname// &
        ' with value = '//value
        end if
        call ftpmsg(errmsg)
        call ftpmsg(keybuf)
    end if
end
subroutine fttnul(ounit,colnum,inull,status)
!
!*******************************************************************************
!
!! FTTNUL defines the null value for a table column.
!
!       ounit   i  Fortran I/O unit number
!       colnum  i  number of the column to be defined
!       inull   i  the value to be use to signify undefined data
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,inull,status

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff

    if (status > 0)return

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)
    if (status > 0)return

!       test for proper HDU type
    if (hdutyp(ibuff) == 0)then
        status=235
        return
    end if

    if (colnum > tfield(ibuff) .or. colnum < 1)then
         status=302
         return
    end if

    tnull(colnum+tstart(ibuff))=inull
end
subroutine fttrec(string,status)
!
!*******************************************************************************
!
!! FTTREC makes sure the characters in a header record are printable.
!
!       i.e., with ASCII codes greater than or equal to 32 (a blank)
!       and less than or equal to 126 (tilda).

!       Note: this will not detect the delete character (ASCII 127)
!       because of the difficulties in also supporting this check
!       on IBM mainframes, where the collating sequence is entirely
!       different.

!       Dec 1996:  since support for non-ASCII character sets has
!       been dropped, the test for characters greater than 126
!       has been restated.

!       string  c*72 keyword record
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)

!       optimized in 7/93 to compare "ichar(string(i:i)) < space"
!       rather than                       "(string(i:i)) < ' ' "
!       This is much faster on SUNs and DECstations,
!       and decreases the time needed to write a keyword (ftprec) by 10%.
!       This change made no difference on a VAX

    integer space,tilda
!       the following 2 lines are only correct for machines that use ASCII
    parameter (space = 32)
    parameter (tilda = 126)
    character string*(*)
    integer status,i
    character pos*2

    if (status > 0)return

    do 20 i=1,72
        if (ichar(string(i:i)) < space .or. &
            ichar(string(i:i)) > tilda) then
!                 illegal character found
              status=207
              write(pos,1000)i
1000              format(i2)
    call ftpmsg('Character #'//pos//' in this keyword value or '// &
    'comment string is illegal:')
    call ftpmsg(string)
              return
        end if
20      continue
end
    function fttrnn(value)
!
!*******************************************************************************
!
!! FTTRNN tests if a R*4 value has a IEEE Not-a-Number (NaN) value.
!
!       A NaN has all the exponent bits=1, and the fractional part not=0.
!       The exponent field occupies bits 23-30,  (least significant bit = 0)
!       The mantissa field occupies bits 0-22

!       This routine also sets any underflow values to zero.

!       written by Wm Pence, HEASARC/GSFC, May 1992
!       modified Aug 1994 to handle all IEEE special values.

    integer value

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    logical fttrnn
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer compid
    common/ftcpid/compid
!

!       COMPID specifies what type of floating point word structure
!       is used on this machine, and determines how to test for NaNs.

!       COMPID value:
!         2 or 3 VAX or generic machine: simply test for NaNs with all bits set
!           1   like a decstation or alpha OSF/1, or IBM PC
!           0   SUN workstation, or IBM mainframe
!          -2305843009213693952   Cray (64-bit) machine

    fttrnn=.false.
    if (compid >= 2)then
!           on the VAX we can assume that all NaNs will be set to all bits on
!           (which is equivalent to an integer with a value of -1) because
!           this is what the IEEE to VAX conversion MACRO program returns
        if (value == -1)fttrnn=.true.
    else if (compid >= -1)then
!           the following test works on all other machines (except Cray)
!           the sign bit may be either 1 or 0 so have to test both possibilites.
!           Note: overflows and infinities are also flagged as NaNs.
        if (value >= 2139095039 .or. (value < 0 .and. &
               value >= -8388609))then
               fttrnn=.true.
        else if ((value > 0 .and. value <= 8388608) .or. &
               value <= -2139095040)then
!                  set underflows and denormalized values to zero
               value=0
        end if
    else
!           branch for the Cray:  COMPID stores the negative integer
!           which corresponds to the 3 most sig digits set to 1.   If these
!           3 bits are set in a floating point number, then it represents
!           a reserved value (i.e., a NaN)
        if (value < 0 .and. value >= compid)fttrnn=.true.
    end if
end
subroutine fttscl(ounit,colnum,bscale,bzero,status)
!
!*******************************************************************************
!
!! FTTSCL defines the scaling factor for a table column.
!
!       ounit   i  Fortran I/O unit number
!       colnum  i  number of the column to be defined
!       bscale  d  scaling factor
!       bzero   d  scaling zero point
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991

    integer ounit,colnum,status
    double precision bscale,bzero

!
    integer nb,ne,nf
    parameter (nb = 20)
    parameter (ne = 512)
    parameter (nf = 3000)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    integer ibuff

    if (status > 0)return

    if (bscale == 0.)then
!               illegal bscale value
            status=322
            return
    end if

    ibuff=bufnum(ounit)

!       if HDU structure is not defined then scan the header keywords
    if (dtstrt(ibuff) < 0)call ftrdef(ounit,status)
    if (status > 0)return

!       test for proper HDU type
    if (hdutyp(ibuff) == 0)then
        status=235
        return
    end if

    if (colnum > tfield(ibuff) .or. colnum < 1)then
         status=302
         return
    end if

    tscale(colnum+tstart(ibuff))=bscale
    tzero(colnum+tstart(ibuff))=bzero
end
subroutine ftucks(iunit,status)
!
!*******************************************************************************
!
!! FTUCKS updates the CHECKSUM keyword value.
!
!  This assumes that the DATASUM
!       keyword exists and has the correct value.

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, May, 1995

    integer iunit,status

!
    integer nf,nb,ne
    parameter (nf = 3000)
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
    integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
    integer theap
    double precision tscale,tzero
    common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb), &
    tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb) &
    ,theap(nb)
!

    double precision sum,dsum
    integer ibuff,nrec,dd,mm,yy,i,tstat
    character datstr*8,string*16,comm*40,datsum*20,oldcks*16
    logical complm

    if (status > 0)return

    ibuff=bufnum(iunit)

!       get the DATASUM keyword value
    call ftgkys(iunit,'DATASUM',datsum,comm,status)
    if (status == 202)then
         call ftpmsg('DATASUM keyword not found (FTUCKS)')
         return
    end if
!
!       decode the datasum string into a double precision variable
!
    do i=1,20
        if (datsum(i:i) /= ' ')then
            call ftc2dd(datsum(i:20),dsum,status)
            go to 15
        end if
    end do
    dsum=0.

!       generate current date string to put into the keyword comment
15      call ftgsdt(dd,mm,yy,status)
    if (status > 0)return

    datstr='  /  /  '
    write(datstr(1:2),1001)dd
    write(datstr(4:5),1001)mm
    write(datstr(7:8),1001)yy
1001    format(i2)

!       replace blank with leading 0 in each field if required
    if (datstr(1:1) == ' ')datstr(1:1)='0'
    if (datstr(4:4) == ' ')datstr(4:4)='0'
    if (datstr(7:7) == ' ')datstr(7:7)='0'

!       get the CHECKSUM keyword value if it exists
    tstat=status
    call ftgkys(iunit,'CHECKSUM',oldcks,comm,status)
    if (status == 202)then
      status=tstat
      oldcks='0000000000000000'
      comm='encoded HDU checksum updated on '//datstr
      call ftpkys(iunit,'CHECKSUM','0000000000000000',comm,status)
    end if

!       rewrite the header END card, and following blank fill
    call ftwend(iunit,status)
    if (status > 0)return

!       move to the start of the header
    call ftmbyt(iunit,hdstrt(ibuff,chdu(ibuff)),.true.,status)

!       accumulate the header checksum into the previous data checksum
    nrec= (dtstrt(ibuff)-hdstrt(ibuff,chdu(ibuff)))/2880
    sum=dsum
    call ftcsum(iunit,nrec,sum,status)

!       encode the COMPLEMENT of the checksum into a 16-character string
    complm=.true.
    call ftesum(sum,complm,string)

!       return if the checksum is correct
    if (string == '0000000000000000')return

    if (oldcks == '0000000000000000')then
!               update the CHECKSUM keyword value with the checksum string
            call ftmkys(iunit,'CHECKSUM',string,'&',status)
    else

!           Zero the checksum and compute the new value
        comm='encoded HDU checksum updated on '//datstr
        call ftmkys(iunit,'CHECKSUM','0000000000000000',comm,status)

!           move to the start of the header
        call ftmbyt(iunit,hdstrt(ibuff,chdu(ibuff)),.true.,status)

!           accumulate the header checksum into the previous data checksum
        sum=dsum
        call ftcsum(iunit,nrec,sum,status)

!           encode the COMPLEMENT of the checksum into a 16-character string
        complm=.true.
        call ftesum(sum,complm,string)

!           update the CHECKSUM keyword value with the checksum string
        call ftmkys(iunit,'CHECKSUM',string,'&',status)
    end if
end
subroutine ftucrd(ounit,keywrd,card,status)
!
!*******************************************************************************
!
!! FTUCRD updates a 80-character FITS header card/record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       card    c  80-character FITS card image
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, May 1995

    character ( len = * ) keywrd,card
    integer ounit,status,tstat

    if (status > 0)return
    tstat=status

!       try modifying the card, if it exists
    call ftmcrd(ounit,keywrd,card,status)

    if (status == 202)then
!               card doesn't exist, so create it
            status=tstat
            call ftprec(ounit,card,status)
    end if
end
subroutine ftukyd(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTUKYD updates a double precision value header record in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (max. 47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim,tstat

    if (status > 0)return
    tstat=status

!       try modifying the keyword, if it exists
    call ftmkyd(ounit,keywrd,dval,decim,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkyd(ounit,keywrd,dval,decim,comm,status)
    end if
end
subroutine ftukye(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTUKYE updates a real*4 value header record in E format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,comm
    real rval
    integer ounit,status,decim,tstat

    if (status > 0)return
    tstat=status

!       try modifying the keyword, if it exists
    call ftmkye(ounit,keywrd,rval,decim,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkye(ounit,keywrd,rval,decim,comm,status)
    end if
end
subroutine ftukyf(ounit,keywrd,rval,decim,comm,status)
!
!*******************************************************************************
!
!! FTUKYF updates a real*4 value header record in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       rval    r  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,comm
    real rval
    integer ounit,status,decim,tstat

    if (status > 0)return
    tstat=status

!       try modifying the keyword, if it exists
    call ftmkyf(ounit,keywrd,rval,decim,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkyf(ounit,keywrd,rval,decim,comm,status)
    end if
end
subroutine ftukyg(ounit,keywrd,dval,decim,comm,status)
!
!*******************************************************************************
!
!! FTUKYG updates a double precision value header record in F format.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       dval    d  keyword value
!       decim   i  number of decimal places to display in value field
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,comm
    double precision dval
    integer ounit,status,decim,tstat

    if (status > 0)return
    tstat=status

!       try modifying the keyword, if it exists
    call ftmkyg(ounit,keywrd,dval,decim,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkyg(ounit,keywrd,dval,decim,comm,status)
    end if
end
subroutine ftukyj(ounit,keywrd,intval,comm,status)
!
!*******************************************************************************
!
!! FTUKYJ updates an integer value header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       intval  i  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,comm
    integer ounit,status,intval,tstat

    if (status > 0)return
    tstat=status

!       try modifying the keyword, if it exists
    call ftmkyj(ounit,keywrd,intval,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkyj(ounit,keywrd,intval,comm,status)
    end if
end
subroutine ftukyl(ounit,keywrd,logval,comm,status)
!
!*******************************************************************************
!
!! FTUKYL updates a logical value header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       logval  l  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,comm
    integer ounit,status,tstat
    logical logval

    if (status > 0)return
    tstat=status

!       try modifying the keyword, if it exists
    call ftmkyl(ounit,keywrd,logval,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkyl(ounit,keywrd,logval,comm,status)
    end if
end
subroutine ftukys(ounit,keywrd,strval,comm,status)
!
!*******************************************************************************
!
!! FTUKYS updates a character string value header record.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       strval  c  keyword value
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS:
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Oct 1994

    character ( len = * ) keywrd,strval,comm
    integer ounit,status,tstat

    if (status > 0)return

    tstat=status
!       try modifying the keyword, if it exists
    call ftmkys(ounit,keywrd,strval,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
!               note that this supports the HEASARC long-string conventions
            call ftpkls(ounit,keywrd,strval,comm,status)
    end if
end
subroutine ftukyu(ounit,keywrd,comm,status)
!
!*******************************************************************************
!
!! FTUKYU updates a null-valued keyword.
!
!       ounit   i  fortran output unit number
!       keywrd  c  keyword name    ( 8 characters, cols.  1- 8)
!       comm    c  keyword comment (47 characters, cols. 34-80)
!       OUTPUT PARAMETERS
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, July 1997

    character ( len = * ) keywrd,comm
    integer ounit,status,tstat

    if (status > 0)return

    tstat=status

!       try modifying the keyword, if it exists
    call ftmkyu(ounit,keywrd,comm,status)

    if (status == 202)then
!               keyword doesn't exist, so create it
            status=tstat
            call ftpkyu(ounit,keywrd,comm,status)
    end if
end
subroutine ftupch(string)
!
!*******************************************************************************
!
!! FTUPCH converts input string to upper case.
!
!       written by Wm Pence, HEASARC/GSFC, February 1991
!       modified 7/93 to use ichar comparisons, to improve performance

    character ( len = * ) string
    integer i,length
    integer a,z

    a=ichar('a')
    z=ichar('z')

    length=len(string)
    do i=1,length
            if   (ichar(string(i:i)) >= a &
            .and. ichar(string(i:i)) <= z)then
                    string(i:i)=char(ichar(string(i:i))-32)
            end if
      end do
end
subroutine ftuptf(iunit,status)
!
!*******************************************************************************
!
!! FTUPTF updates the TFORM keywords for the variable length array columns.
!
!  This is to make sure they all have the form 1Pt(len) or Pt(len)
!       where 'len' is the maximum length of the vector in the table (e.g.,
!       '1PE(400)')

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Jan, 1997

    integer iunit,status

    integer ii,tflds,naxis2,maxlen,jj,length,addr,endpos,cpos
    character comment*80, keynam*8,tform*40,newfrm*40
    character message*80,lenstr*20

    call ftgkyj(iunit,'TFIELDS', tflds, comment, status)
    call ftgkyj(iunit,'NAXIS2', naxis2, comment, status)

    do 100 ii = 1,tflds
       call ftkeyn('TFORM',ii,keynam,status)
       call ftgkys(iunit,keynam,tform,comment,status)

       if (status > 0)then
           message='Error while updating variable length TFORMn '// &
           'values (ftuptf)'
           call ftpmsg(message)
       end if

!          test if this is a variable length array column
       if (tform(1:1) == 'P' .or. tform(2:2) == 'P')then

!             test if the length field is missing
          if (tform(5:5) == ' ')then

             maxlen = 0
             do 50 jj=1,naxis2
                call ftgdes(iunit,ii,jj,length,addr,status)
                maxlen = max(maxlen,length)
50               continue

             if (tform(1:1) == 'P')then
            endpos=3
             else
            endpos=4
             end if

!                convert integer to C*20 string, and find first digit
             call fti2c(maxlen,lenstr,status)
             do 60 jj = 1, 20
                cpos = jj
                if (lenstr(cpos:cpos) /= ' ')go to 70
60               continue

!                construct new keyword value
70               newfrm=tform
             newfrm(endpos:)='('//lenstr(cpos:20)//')'

!                now modify the old TFORMn keyword
             call ftmkys(iunit,keynam,newfrm,comment,status)
          end if
       end if
100     continue
end
subroutine ftuscc(input,np,scaled,scale,zero,output)
!
!*******************************************************************************
!
!! FTUSCC unscales a complex array prior to writing to the FITS file.
!
!       input  r  array of complex numbers (pairs of real/imaginay numbers)
!       np     i  total number of values to scale (no. of pairs times 2)
!       scaled l  is the data scaled?
!       scale  d  scale factor
!       zero   d  offset
!       output r  output array

    integer np,i,j
    logical scaled
    real input(np),output(np)
    double precision scale,zero

    j=1
    if (scaled)then
        do 10 i=1,np/2
            output(j)=(input(j)-zero)/scale
            j=j+1
!               the imaginary part of the number is not offset!!
            output(j)=input(j)/scale
            j=j+1
10          continue
    else
        do 20 i=1,np
            output(i)=input(i)
20          continue
    end if
end
subroutine ftuscm(input,np,scaled,scale,zero,output)
!
!*******************************************************************************
!
!! FTUSCM unscales a complex array prior to writing to the FITS file.
!
!       input  d  array of complex numbers (pairs of real/imaginay numbers)
!       np     i  total number of values to scale (no. of pairs times 2)
!       scaled l  is the data scaled?
!       scale  d  scale factor
!       zero   d  offset
!       output d  output array

    integer np,i,j
    logical scaled
    double precision input(np),output(np)
    double precision scale,zero

    j=1
    if (scaled)then
        do 10 i=1,np/2
            output(j)=(input(j)-zero)/scale
            j=j+1
!               the imaginary part of the number is not offset!!
            output(j)=input(j)/scale
            j=j+1
10          continue
    else
        do 20 i=1,np
            output(i)=input(i)
20          continue
    end if
end
subroutine ftvcks(iunit,dataok,hduok,status)
!
!*******************************************************************************
!
!! FTVCKS verifies the HDU by comparing checksums.
!
!  The routine compares the value of the computed checksums against
!       the values of the DATASUM and CHECKSUM keywords if they are present.

!       iunit   i  fortran unit number
!       dataok  i  output verification code for the data unit alone
!       hduok   i  output verification code for the entire HDU
!                  the code values = 1  verification is correct
!                                  = 0  checksum keyword is not present
!                                  = -1 verification not correct
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Dec, 1994

    integer iunit,dataok,hduok,status,tstat,i
    double precision datsum,chksum,dsum
    character keyval*20,comm*8
    logical cexist,dexist

    if (status > 0)return

!       check if the CHECKSUM keyword exists
    tstat=status
    call ftgkys(iunit,'CHECKSUM',keyval,comm,status)
    if (status <= 0)then
        cexist=.true.
    else
        hduok=0
        cexist=.false.
        status=tstat
    end if

!       check if the DATASUM keyword exists and get its value
    call ftgkys(iunit,'DATASUM',keyval,comm,status)
    if (status <= 0)then
        dexist=.true.
    else
        dataok=0
        dexist=.false.
        status=tstat
    end if

!       return if neither keyword exists
    if (.not. cexist .and. .not. dexist)return

!       calculate the data checksum and the HDU checksum
    call ftgcks(iunit,datsum,chksum,status)
    if (status > 0)return

    if (dexist)then

!           decode the datasum string into a double precision variable
        do i=1,20
            if (keyval(i:i) /= ' ')then
                call ftc2dd(keyval(i:20),dsum,status)
                if (status == 409)then
!                       couldn't read the keyword; assume it is out of date
                    status=tstat
                    dsum=-1.
                end if
                go to 15
            end if
        end do
        dsum=0.
15          continue

        if (dsum == datsum)then
            dataok=1
        else
            dataok=-1
        end if
    end if

    if (cexist)then
        if (chksum == 0 .or. chksum == 4.294967295D+09)then
            hduok=1
        else
            hduok=-1
        end if
    end if
end
subroutine ftvers(vernum)
!
!*******************************************************************************
!
!! FTVERS returns the current revision number of the FITSIO package.
!
!
!       The revision number will be incremented whenever any modifications,
!       bug fixes, or enhancements are made to the package

    real vernum

    vernum=5.03

!       version 5.03    1 Aug   1997   null keywords; keyword units
!       version 5.02   11 Apr   1997   F90 portability modifications
!       version 5.01    8 Apr   1997   OK if nelems = 0 when accessing tables
!       version 5.00   21 Mar   1997   Major overhaul:
!                                      more efficient; F90-compatible
!       version 4.14   13 Aug   1996   initialize lenval in ftdkey; check for
!                                      valid unit number in ftclos
!       version 4.13   22 Mar   1996   add ftflus; prevent duplicate header
!                                      keyword in ftphpr, ftphtb, and ftphbn
!       version 4.12   28 Feb   1996   added fticls subroutine
!       version 4.11    8 Feb   1996   add calls to ftrdef in ftptbb and ftptbs
!                                      bug in ftdrow
!       version 4.10    1 Dec   1995   fixed pattern matching bugs in ftcmps;
!       version 4.09    8 Nov   1995   don't update header pointer in ftprec;
!                                      open blocksize optimized in fitsvax.f
!       version 4.08    3 Oct   1995   bug in ftiimg: data offset by 8 bytes
!       version 4.07    7 Sept  1995   fticol failed on ASCII columns
!       version 4.06   18 Aug   1995   ftdelt bug; ftpmsg saves latest errors
!       version 4.05    2 Aug   1995   another bug in ftfrcl in reseting tstart
!       version 4.04   12 Jul   1995   bug in ftfrcl in resetting tstart
!       version 4.03    3 Jul   1995   bug in restoring CHDU when moving to EOF
!       version 4.02   20 Jun   1995   modified checksum algorithm
!       version 4.01   30 May   1995   many changes
!       version 3.711  30 Jan   1995   ftgphx was cutting BSCALE to 20 chars
!       version 3.710  27 Jan   1995   fix ftgcnn, fitsmac; add ftirec, ftdrec
!       version 3.700  29 Dec   1994   public release
!       version 3.623   8 Nov   1994   ftgkys, ftgnst, checksum
!       version 3.622   7 Nov   1994   ftgclj R*8 alignment; I*2 overflow fti4i2
!       version 3.621   4 Nov   1994   fixed endhd position in ftgrec
!       version 3.62    2 Nov   1994   ftgcx[ijd] routines added
!       version 3.612  31 Oct   1994   restored previous FTIBLK algorithm
!       version 3.61   26 Oct   1994   ftirow and ftdrow to modify tables
!       version 3.6    18 Oct   1994   ftukyX, range checking, new EOF checks
!       version 3.512  20 Sep   1994   fixed writing header fill in FTWEND
!       version 3.511  20 Sep   1994   removed '=' from CONTINUE on long strings
!       version 3.51   14 Sep   1994   long string convention and IEEE support
!       version 3.504  22 Aug   1994   fixed bug in ftcopy making files too big
!       version 3.503   8 Aug   1994   fixed bug in ftcopy making files too big
!       version 3.502  26 Jul   1994   explicitly write data fill bytes
!       version 3.501  19 Jul   1994   minor changes for FTOOLS release
!       version 3.500  29 Jun   1994   added error message stack
!       version 3.415  07 Jun   1994   fixed ftmahd and ftgrec
!       version 3.414  18 May   1994   modify ftmoff and ftpbyt for status 112
!       version 3.413  18 Mar   1994   Cray port added
!       version 3.412  01 Mar   1994   SUN internal read problem in ftgthd
!       version 3.411  25 Feb   1994   fixed 107 error when reading byte column
!       version 3.410  21 Jan   1994   bug fixes in Alpha VMS version
!       version 3.409  21 Dec   1993   long string bug; HP support
!       version 3.408  09 Nov   1993   Alpha VMS open; ftgthd -; 210 status
!       version 3.407  02 Nov   1993   initialize TABLEs with blanks; ftrdef
!       version 3.406  26 Oct   1993   ftgtdm bug - last not initialized
!                                      modified to read unknown extenstions
!       version 3.405  21 Oct   1993   ftpini bug with GROUP format files
!       version 3.404   7 Oct   1993   new TDIM subroutines, new error status
!       version 3.403   1 Sept  1993   initialize strlen in ftpkys
!       version 3.402  23 Aug   1993   bug in ftgcno
!       version 3.401  20 Aug   1993   minor change to ftpi1b
!       version 3.4  - 11 Aug   1993
!       version 3.31 -  2 Feb   1993
!       version 3.3  - 28 Oct   1992
!       version 3.21 -  8 July  1992
!       version 3.20 - 30 Mar   1992
!       version 3.10 -  4 Nov   1991
!       version 3.01 - 27 Sept  1991
!       version 3.00 - 12 Sept  1991
!       version 2.99 - 24 July  1991
!       version 2.0  -  1 May   1991
!       version 1.3  -  2 April 1991
!       version 1.22 - 22 March 1991
!       version 1.21 - 20 March 1991

end
subroutine ftwend(iunit,status)
!
!*******************************************************************************
!
!! FTWEND writes the END card.
!
!       write the END card, and following fill values in the CHDU

!       iunit   i  fortran unit number
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, Aug 1994

    integer iunit,status

!
    integer nb,ne
    parameter (nb = 20)
    parameter (ne = 512)
    integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
    integer nxtfld
    logical wrmode
    common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
    wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
!

    integer ibuff,nblank,i,endpos
    character*80 rec

    if (status > 0)return

    ibuff=bufnum(iunit)

!       calc the data starting position if not currently defined
    if (dtstrt(ibuff) < 0)then
            dtstrt(ibuff)=(hdend(ibuff)/2880 + 1)*2880
    end if

!       calculate the number of blank keyword slots in the header
endpos=hdend(ibuff)
    nblank=(dtstrt(ibuff)-endpos)/80
!       move the i/o pointer to the end of the header keywords
    call ftmbyt(iunit,endpos,.true.,status)

!       fill all the slots with blanks
    rec=' '
    do 10 i=1,nblank
            call ftpcbf(iunit,80,rec,status)
10      continue

!               The END keyword must either be placed
!               immediately after the last keyword that was written
!               (as indicated by the HDEND parameter), or must be in the
!               first 80 bytes of the FITS record immediately preceeding
!               the data unit, whichever is further in the file.
!               The latter will occur if the user reserved room for more
!               header keywords which have not (yet) been filled.

!       move pointer to where the END card should be
endpos=max(endpos,dtstrt(ibuff)-2880)
    call ftmbyt(iunit,endpos,.true.,status)

!       write the END record to the output buffer:
    rec='END'
    call ftpcbf(iunit,80,rec,status)

    if (status > 0)then
        call ftpmsg('Error while writing END card (FTWEND).')
    end if
end
subroutine ftwldp(xpix,ypix,xref,yref,xrefpix,yrefpix, &
                      xinc,yinc,rot,type,xpos,ypos,status)
!
!*******************************************************************************
!
!! FTWLDP determines world position from pixel coordinates.
!
!       Fortran version of worldpos.c -- WCS Algorithms from Classic AIPS
!       Translated by James Kent Blackburn -- HEASARC/GSFC/NASA -- November 1994
!       routine to determine accurate position from pixel coordinates
!       does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections
!       returns   0 = good,
!               501 = angle too large for projection;
! Input:
! dbl   xpix    x pixel number  (RA or long without rotation)
! dbl   ypiy    y pixel number  (dec or lat without rotation)
! dbl   xref    x reference coordinate value (deg)
! dbl   yref    y reference coordinate value (deg)
! dbl   xrefpix x reference pixel
! dbl   yrefpix y reference pixel
! dbl   xinc    x coordinate increment (deg)
! dbl   yinc    y coordinate increment (deg)
! dbl   rot     rotation (deg)  (from N through E)
! chr   type    projection type code e.g. "-SIN"
! Output:
! dbl    xpos   x (RA) coordinate (deg)
! dbl    ypos   y (dec) coordinate (deg)
! int   status  error status flag, zero

    integer           status
    double precision  xpix,ypix,xref,yref,xrefpix,yrefpix
    double precision  xinc,yinc,rot,xpos,ypos
    character ( len = * )     type
    integer           error1,error4
    parameter         (error1=501)
    parameter         (error4=504)

    double precision  cosr,sinr,dx,dy,dz,temp
    double precision  sins,coss,dect,rat,dt,l,m,mg,da,dd,cos0,sin0
    double precision  dec0,ra0,decout,raout
    double precision  geo1,geo2,geo3
    double precision  cond2r
    parameter         (cond2r=1.745329252d-2)
    double precision  twopi,deps
    parameter         (twopi = 6.28318530717959)
    parameter         (deps = 1.0d-5)
    integer           i,itype
    character*4       ctypes(8)
    data ctypes/ '-SIN', '-TAN', '-ARC', '-NCP', &
                 '-GLS', '-MER', '-AIT', '-STG' /

    if (status > 0) return
!  Offset from ref pixel
    dx = (xpix-xrefpix) * xinc
    dy = (ypix-yrefpix) * yinc
!  Take out rotation
    cosr = dcos(rot*cond2r)
    sinr = dsin(rot*cond2r)
    if (rot /= 0.0) then
      temp = dx * cosr - dy * sinr
      dy = dy * cosr + dx * sinr
      dx = temp
    end if
!  Find type of coordinate transformation (0 is linear)
!   WDP, 1/97: removed support for default type, to give better error ck
!       itype = 0
    itype = -1
    do 10 i = 1, 8
      if (ctypes(i) == type) itype = i
  10    continue
!  default, linear result for error return
    xpos = xref + dx
    ypos = yref + dy
!  Convert to radians
    ra0 = xref * cond2r
    dec0 = yref * cond2r
    l = dx * cond2r
    m = dy * cond2r
    sins = l*l + m*m
    decout = 0.0
    raout = 0.0
    cos0 = dcos(dec0)
    sin0 = dsin(dec0)
!  Process by case
    if (itype == 0) then
!  LINEAR
      rat =  ra0 + l
      dect = dec0 + m
    else if (itype == 1) then
!  SINE from '-SIN' type
      if (sins > 1.0) then
        status = error1
        goto 30
      end if
      coss = dsqrt(1.0 - sins)
      dt = sin0 * coss + cos0 * m
      if ((dt > 1.0) .or. (dt < -1.0)) then
        status = error1
        goto 30
      end if
      dect = dasin(dt)
      rat = cos0 * coss - sin0 * m
      if ((rat == 0.0) .and. (l == 0.0)) then
        status = error1
        goto 30
      end if
      rat = datan2 (l, rat) + ra0
    else if (itype == 2) then
!  TANGENT from '-TAN' type
      if (sins > 1.0) then
        status = error1
        goto 30
      end if
      dect = cos0 - m * sin0
      if (dect == 0.0) then
        status = error1
        goto 30
      end if
      rat = ra0 + datan2(l, dect)
      dect = datan(dcos(rat-ra0) * (m * cos0 + sin0) / dect)
    else if (itype == 3) then
!  Arc from '-ARC' type
      if (sins >= twopi * twopi / 4.0) then
        status = error1
        goto 30
      end if
      sins = dsqrt(sins)
      coss = dcos(sins)
      if (sins /= 0.0) then
        sins = dsin(sins) / sins
      else
        sins = 1.0
      end if
      dt = m * cos0 * sins + sin0 * coss
      if ((dt > 1.0) .or. (dt < -1.0)) then
        status = error1
        goto 30
      end if
      dect = dasin(dt)
      da = coss - dt * sin0
      dt = l * sins * cos0
      if ((da == 0.0) .and. (dt == 0.0)) then
        status = error1
        goto 30
      end if
      rat = ra0 + datan2(dt, da)
    else if (itype == 4) then
!  North Celestial Pole from '-NCP' type
      dect = cos0 - m * sin0
      if (dect == 0.0) then
        status = error1
        goto 30
      end if
      rat = ra0 + datan2(l, dect)
      dt = dcos(rat-ra0)
      if (dt == 0.0) then
        status = error1
        goto 30
      end if
      dect = dect / dt
      if ((dect > 1.0) .or. (dect < -1.0)) then
        status = error1
        goto 30
      end if
      dect = dacos(dect)
      if (dec0 < 0.0) dect = -dect
    else if (itype == 5) then
!  Global Sinusoid from '-GLS' type
      dect = dec0 + m
      if (dabs(dect) > twopi/4.0) then
        status = error1
        goto 30
      end if
      coss = dcos(dect)
      if (dabs(l) > twopi*coss/2.0) then
        status = error1
        goto 30
      end if
      rat = ra0
      if (coss > deps) rat = rat + l / coss
    else if (itype == 6) then
!  Mercator from '-MER' type
      dt = yinc * cosr + xinc * sinr
      if (dt == 0.0) dt = 1.0
      dy = (yref/2.0 + 45.0) * cond2r
      dx = dy + dt / 2.0 * cond2r
      dy = dlog(dtan(dy))
      dx = dlog(dtan(dx))
      geo2 = dt * cond2r / (dx - dy)
      geo3 = geo2 * dy
      geo1 = dcos(yref * cond2r)
      if (geo1 <= 0.0) geo1 = 1.0
      rat = l / geo1 + ra0
      if (dabs(rat - ra0) > twopi) then
        status = error1
        goto 30
      end if
      dt = 0.0
      if (geo2 /= 0.0) dt = (m + geo3) / geo2
      dt = dexp(dt)
      dect = 2.0 * datan(dt) - twopi / 4.0
    else if (itype == 7) then
!  Aitoff from '-AIT' type
      dt = yinc * cosr + xinc * sinr
      if (dt == 0.0) dt = 1.0
      dt = dt * cond2r
      dy = yref * cond2r
      dx = dsin(dy+dt)/dsqrt((1.0+dcos(dy+dt))/2.0) - &
           dsin(dy)/dsqrt((1.0+dcos(dy))/2.0)
      if (dx == 0.0) dx = 1.0
      geo2 = dt / dx
      dt = xinc * cosr - yinc * sinr
      if (dt == 0.0) dt = 1.0
      dt = dt * cond2r
      dx = 2.0 * dcos(dy) * dsin(dt/2.0)
      if (dx == 0.0) dx = 1.0
      geo1 = dt * dsqrt((1.0+dcos(dy)*dcos(dt/2.0))/2.0) / dx
      geo3 = geo2 * dsin(dy) / dsqrt((1.0+dcos(dy))/2.0)
      rat = ra0
      dect = dec0
      if ((l == 0.0) .and. (m == 0.0)) goto 20
      dz = 4.0-l*l/(4.0*geo1*geo1)-((m+geo3)/geo2)*((m+geo3)/geo2)
      if ((dz > 4.0) .or. (dz < 2.0)) then
        status = error1
        goto 30
      end if
      dz = 0.5 * dsqrt(dz)
      dd = (m+geo3) * dz / geo2
      if (dabs(dd) > 1.0) then
        status = error1
        goto 30
      end if
      dd = dasin(dd)
      if (dabs(dcos(dd)) < deps) then
        status = error1
        goto 30
      end if
      da = l * dz / (2.0 * geo1 * dcos(dd))
      if (dabs(da) > 1.0) then
        status = error1
        goto 30
      end if
      da = dasin(da)
      rat = ra0 + 2.0 * da
      dect = dd
    else if (itype == 8) then
!  Stereographic from '-STG' type
      dz = (4.0 - sins) / (4.0 + sins)
      if (dabs(dz) > 1.0) then
        status = error1
        goto 30
      end if
      dect = dz * sin0 + m * cos0 * (1.0+dz) / 2.0
      if (dabs(dect) > 1.0) then
        status = error1
        goto 30
      end if
      dect = dasin(dect)
      rat = dcos(dect)
      if (dabs(rat) < deps) then
        status = error1
        goto 30
      end if
      rat = l * (1.0+dz) / (2.0 * rat)
      if (dabs(rat) > 1.0) then
        status = error1
        goto 30
      end if
      rat = dasin(rat)
      mg = 1.0 + dsin(dect)*sin0 + dcos(dect)*cos0*dcos(rat)
      if (dabs(mg) < deps) then
        status = error1
        goto 30
      end if
      mg = 2.0 * (dsin(dect)*cos0 - dcos(dect)*sin0*dcos(rat)) / mg
      if (dabs(mg-m) > deps) rat = twopi/2.0 - rat
      rat = ra0 + rat
    else
!  Unsupported Projection
      status = error4
      goto 30
    end if
  20    continue
!  Return RA in range
    raout = rat
    decout = dect
    if (raout-ra0 > twopi/2.0) raout = raout - twopi
    if (raout-ra0 < (-twopi)/2.0) raout = raout + twopi
    if (raout < 0.0) raout = raout + twopi
!  Correct units back to degrees
    xpos = raout / cond2r
    ypos = decout / cond2r
  30    continue
end
subroutine ftxiou(iounit,status)
!
!*******************************************************************************
!
!! FTXIOU manages logical unit numbers in the range 50-99.
!
    integer iounit,status,i
    integer*2 array(50)
    save array
    data array/50*0/

    if (iounit == 0)then
!           get an unused logical unit number
        do i=50,1,-1

!        The following would be a more robust way of testing for
!        an available unit number, however, this cannot work
!        when building FITSIO using the IRAF/SPP version, because
!        IRAF does not use Fortran I/O.
!
!                inquire(unit=iounit, exist=exists, opened=open)
!                if(exists .and. .not. open)then
!                    array(iounit-49)=1
!                    return
!                end if

             if (array(i) == 0)then
                 array(i)=1
                 iounit=i+49
                 return
             end if
        end do
!           error: all units are allocated
        iounit=-1
        status=114
        call ftpmsg('FTGIOU has no more available unit numbers.')

    else if (iounit == -1)then
!           deallocate all the unit numbers
        do i=1,50
             array(i)=0
        end do

    else
!            deallocate a specific unit number
         if (iounit >= 50 .and. iounit <= 99)then
             array(iounit-49)=0
         end if
    end if

    return
end
subroutine ftxmsg(action,text)
!
!*******************************************************************************
!
!! FTXMSG gets, puts, or clears the error message stack.
!

    integer action
    character ( len = * ) text

    integer nbuff,i
    parameter (nbuff=50)
    character ( len = 80 ) txbuff(nbuff)
    save txbuff
    data txbuff/nbuff * ' '/

    if (action == -1)then

!           get error message from top of stack and shift the stack up one
        text=txbuff(1)
        do i=1,nbuff-1
            txbuff(i) = txbuff(i+1)
        end do
        txbuff(nbuff)=' '

    else if (action == 1)then

!           put error message onto stack.
        do i=1,nbuff
            if (txbuff(i) == ' ')then
               txbuff(i)=text
               return
            end if
        end do
!           stack is full so discard oldest message
        do i=1,nbuff-1
            txbuff(i) = txbuff(i+1)
        end do
        txbuff(nbuff)=text

    else if (action == 0)then

!           clear the error message stack
        do i=1,nbuff
            txbuff(i) = ' '
        end do

    end if
end
subroutine ftxypx(xpos,ypos,xref,yref,xrefpix,yrefpix, &
                      xinc,yinc,rot,type,xpix,ypix,status)
!
!*******************************************************************************
!
!! FTXYPX determines pixel coordinates from right ascension and declination.
!
!       Fortran version of worldpos.c -- WCS Algorithms from Classic AIPS
!       Translated by James Kent Blackburn -- HEASARC/GSFC/NASA -- November 1994
!       routine to determine accurate pixel coordinates from an RA and Dec
!       does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections
!       returns   0 = good,
!               501 = angle too large for projection;
!               502 = bad values
!               503 = ???undocumented error - looks like an underflow???
! Input:
! dbl   xpos   x (RA) coordinate (deg)
! dbl   ypos   y (dec) coordinate (deg)
! dbl   xref    x reference coordinate value (deg)
! dbl   yref    y reference coordinate value (deg)
! dbl   xrefpix x reference pixel
! dbl   yrefpix y reference pixel
! dbl   xinc    x coordinate increment (deg)
! dbl   yinc    y coordinate increment (deg)
! dbl   rot     rotation (deg)  (from N through E)
! chr   type    projection type code e.g. "-SIN"
! Output:
! dbl   xpix    x pixel number  (RA or long without rotation)
! dbl   ypiy    y pixel number  (dec or lat without rotation)
! int   status  error status flag, zero

    integer           status
    double precision  xpos,ypos,xref,yref,xrefpix,yrefpix
    double precision  xinc,yinc,rot,xpix,ypix
    character ( len = * )     type
    integer           error1,error2,error3,error4
    parameter         (error1=501)
    parameter         (error2=502)
    parameter         (error3=503)
    parameter         (error4=504)
    double precision  dx,dy,dz,r,ra0,dec0,ra,dec
    double precision  coss,sins,dt,da,dd,sint,oldxpos
    double precision  l,m,geo1,geo2,geo3,sinr,cosr
    double precision  cond2r
    parameter         (cond2r=1.745329252d-2)
    double precision  twopi,deps
    parameter         (twopi = 6.28318530717959)
    parameter         (deps = 1.0d-5)
    integer           i,itype
    character*4       ctypes(8)
    data ctypes/ '-SIN', '-TAN', '-ARC', '-NCP', &
                 '-GLS', '-MER', '-AIT', '-STG' /

    if (status > 0) return
!  0 hour wrap around test
    oldxpos = xpos
    dt = (xpos - xref)
    if (dt > +180) xpos = xpos - 360
    if (dt < -180) xpos = xpos + 360
!  Default values - Linear
    dx = xpos - xref
    dy = ypos - yref
    dz = 0.0
!  Correct for rotation
    r = rot * cond2r
    cosr = dcos(r)
    sinr = dsin(r)
    dz = dx * cosr + dy * sinr
    dy = dy * cosr - dx * sinr
    dx = dz
!  Check axis increments - bail out if either 0
    if ((xinc == 0.0) .or. (yinc == 0.0)) then
      xpix = 0.0
      ypix = 0.0
      status = error2
      goto 30
    end if
    xpix = dx / xinc + xrefpix
    ypix = dy / yinc + yrefpix
!  Find type of coordinate transformation (0 is linear)
!   WDP, 1/97: removed support for default type, to give better error ck
!       itype = 0
    itype = -1
    do 10 i = 1, 8
      if (ctypes(i) == type) itype = i
  10    continue
!  Done if linear
    if (itype == 0) goto 30
!  Non-Linear position
    ra0 = xref * cond2r
    dec0 = yref * cond2r
    ra = xpos * cond2r
    dec = ypos * cond2r
!  Compute directional cosine
    coss = dcos(dec)
    sins = dsin(dec)
    l = dsin(ra-ra0) * coss
    sint = sins * dsin(dec0) + coss * dcos(dec0) * dcos(ra-ra0)
!  Process by case
    if (itype == 1) then
!  SINE from '-SIN' type
      if (sint < 0.0) then
        status = error1
        goto 30
      end if
      m = sins * dcos(dec0) - coss * dsin(dec0) * dcos(ra-ra0)
    else if (itype == 2) then
!  TANGENT from '-TAN' type
      if (sint <= 0.0) then
        status = error1
        goto 30
      end if
      m = sins * dsin(dec0) + coss * dcos(dec0) * dcos(ra-ra0)
      l = l / m
      m = (sins*dcos(dec0) - coss*dsin(dec0)*dcos(ra-ra0)) / m
    else if (itype == 3) then
!  Arc from '-ARC' type
      m = sins*dsin(dec0) + coss*dcos(dec0)*dcos(ra-ra0)
      if (m < -1.0) m = -1.0
      if (m > 1.0) m = 1.0
      m = dacos(m)
      if (m /= 0) then
        m = m / dsin(m)
      else
        m = 1.0
      end if
      l = l * m
      m = (sins*dcos(dec0) - coss*dsin(dec0)*dcos(ra-ra0)) * m
    else if (itype == 4) then
!  North Celestial Pole from '-NCP' type
      if (dec0 == 0.0) then
        status = error1
        goto 30
      else
        m = (dcos(dec0) - coss * dcos(ra-ra0)) / dsin(dec0)
      end if
    else if (itype == 5) then
!  Global Sinusoid from '-GLS' type
      dt = ra - ra0
      if (dabs(dec) > twopi/4.0) then
        status = error1
        goto 30
      end if
      if (dabs(dec0) > twopi/4.0) then
        status = error1
        goto 30
      end if
      m = dec - dec0
      l = dt * coss
    else if (itype == 6) then
!  Mercator from '-MER' type
      dt = yinc * cosr + xinc * sinr
      if (dt == 0.0) dt = 1.0
      dy = (yref/2.0 + 45.0) * cond2r
      dx = dy + dt / 2.0 * cond2r
      dy = dlog(dtan(dy))
      dx = dlog(dtan (dx))
      geo2 = dt * cond2r / (dx - dy)
      geo3 = geo2 * dy
      geo1 = cos (yref * cond2r)
      if (geo1 <= 0.0) geo1 = 1.0
      dt = ra - ra0
      l = geo1 * dt
      dt = dec / 2.0 + twopi / 8.0
      dt = dtan(dt)
      if (dt < deps) then
        status = error2
        goto 30
      end if
      m = geo2 * dlog(dt) - geo3
    else if (itype == 7) then
! Aitoff from '-AIT' type
      l = 0.0
      m = 0.0
      da = (ra - ra0) / 2.0
      if (dabs(da) > twopi/4.0) then
        status = error1
        goto 30
      end if
      dt = yinc * cosr + xinc * sinr
      if (dt == 0.0) dt = 1.0
      dt = dt * cond2r
      dy = yref * cond2r
      dx = dsin(dy+dt)/dsqrt((1.0+dcos(dy+dt))/2.0) - &
           dsin(dy)/dsqrt((1.0+dcos(dy))/2.0)
      if (dx == 0.0) dx = 1.0
      geo2 = dt / dx
      dt = xinc * cosr - yinc * sinr
      if (dt == 0.0) dt = 1.0
      dt = dt * cond2r
      dx = 2.0 * dcos(dy) * dsin(dt/2.0)
      if (dx == 0.0) dx = 1.0
      geo1 = dt*dsqrt((1.0+dcos(dy)*dcos(dt/2.0))/2.0)/dx
      geo3 = geo2 * dsin(dy) / dsqrt((1.0+dcos(dy))/2.0)
      dt = dsqrt ((1.0 + dcos(dec) * dcos(da))/2.0)
      if (dabs(dt) < deps) then
        status = error3
        goto 30
      end if
      l = 2.0 * geo1 * dcos(dec) * dsin(da) / dt
      m = geo2 * dsin(dec) / dt - geo3
    else if (itype == 8) then
!  Stereographic from '-STG' type
      da = ra - ra0
      if (dabs(dec) > twopi/4.0) then
        status = error1
        goto 30
      end if
      dd = 1.0 + sins*dsin(dec0) + coss*dcos(dec0)*dcos(da)
      if (dabs(dd) < deps) then
        status = error1
        goto 30
      end if
      dd = 2.0 / dd
      l = l * dd
      m = dd * (sins*dcos(dec0) - coss*dsin(dec0)*dcos(da))
    else
!  Unsupported Projection
      status = error4
      goto 30
    end if
!  Convert back to degrees
    dx = l / cond2r
    dy = m / cond2r
!  Correct for rotation
    dz = dx * cosr + dy * sinr
    dy = dy * cosr - dx * sinr
    dx = dz
!  Convert to PIXELS ... yeah!
    xpix = dx / xinc + xrefpix
    ypix = dy / yinc + yrefpix
  30    continue
!  reset xpos to correct for in place modification
    xpos = oldxpos
end


