!
! fitsfor.F
! A support routine for writing FITS file,
! downloaded from http://orion.math.iastate.edu/burkardt/g_src/fitsio/fitsfort.f90
!

subroutine ftopnx(funit,fname,oldnew,rwmode,block,status)
!
!*******************************************************************************
!
!! FTOPNX is a low-level routine to create or open a new file.
!
!
!  Author:
!
!    William Pence, HEASARC/GSFC
!
!  Parameters:
!
!       funit   i  Fortran I/O unit number
!       fname   c  name of file to be opened
!       oldnew  i  file status: 0 = open old/existing file; else open new file
!       rwmode  i  file access mode: 0 = readonly; else = read/write
!       block   i  FITS record blocking factor
!       status  i  returned error status (0=ok)
!
  integer funit,oldnew,rwmode,block,status,i,ibuff,inital,size
  character ( len = * ) fname
  logical igneof,found

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
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

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer pb
  parameter (pb = 20)
  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)

  integer compid
  common/ftcpid/compid
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  real rword
  double precision dword

  save inital
  data inital/0/

  if (status > 0)return

  if (inital == 0)then
!           first time through need to initialize pointers
      nxtfld=0
      maxbuf=pb
      do i=1,nb
          buflun(i)=0
      end do
      do i=1,pb
          logbuf(i)=0
          recnum(i)=0
          modify(i)=.false.
          pindex(i)=i
      end do
      inital=1

!           Determine at run time what type of machine we are running on.
!           Initialize a real and double value to arbitrary values.
      rword=1.1111111111
      dword=1.1111111111D+00

!           ftarch looks at the equivalent integer value
      call ftarch(rword,dword,compid)
  end if
!
!       check for valid unit number.
!
  if (funit < 1 .or. funit > 199)then
          status=101
          return
  end if
!
!       find available logical buffer slot for this file
!
  do i=1,nb
          if (buflun(i) == 0)then
                  ibuff=i
                  go to 20
          end if
  end do

!       error: no vacant logical buffer slots left
  status=102
  return

20      continue

  if (oldnew == 0)then
      igneof = .false.
!           test if file exists
      inquire(file=fname,exist=found)
      if (.not. found)then
!               error: file doesn't exist??
          status=103
          return
      end if
  else
      igneof = .true.
  end if

  call ftopnf(funit,fname,oldnew,rwmode,block,size,status)

!       initialize the HDU parameters
  maxrec(ibuff)=size

  if (oldnew == 1 .or. block <= 1)then
!           new files always have a record length of 2880 bytes
      reclen(ibuff)=2880
  else
      reclen(ibuff)=block
  end if

  bufnum(funit)=ibuff
  chdu(ibuff)=1
  hdutyp(ibuff)=0
  maxhdu(ibuff)=1
  hdstrt(ibuff,1)=0
  hdend(ibuff)=0
  nxthdr(ibuff)=0
!
!  data start location is undefined
!
  dtstrt(ibuff)=-1000000000

  heapsz(ibuff)=0
  theap(ibuff)=0
  tfield(ibuff)=0
  rowlen(ibuff)=0
!
!  initialize the logical buffer parameters
!
  buflun(ibuff)=funit
  currnt(ibuff)=0

  if (rwmode == 0)then
          wrmode(ibuff)=.false.
  else
          wrmode(ibuff)=.true.
  end if

!       load the first record of the file
  call ftldrc(funit,1,igneof,status)

  return
end
subroutine ftclsx(iunit,keep,status)
!
!*******************************************************************************
!
!! FTCLSX is a low level routine to close a file.
!
!
!  Parameters:
!
!       iunit   i  Fortran I/O unit number
!       keep    l  keep the file? (else delete it)
!       status  i  returned error status (0=ok)
!
!       written by Wm Pence, HEASARC/GSFC, Aug 1992

  integer iunit,status
  logical keep

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne
  parameter (nb = 20)
  parameter (ne = 512)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer ibuff

  ibuff=bufnum(iunit)

  if (ibuff == 0)return
!
!       reset file common block parameters.
!
  bufnum(iunit)=0
  buflun(ibuff)=0
  wrmode(ibuff)=.false.
  currnt(ibuff)=0
  reclen(ibuff)=0
  bytnum(ibuff)=0

  if (keep)then
    close(iunit,err=900)
  else
    close(iunit,status='DELETE',err=900)
  end if

  return

900     continue
!       set error code, if it has not previous been set
  if (status <= 0)status=110

  return
end
subroutine ftmbyt(iunit,bytno,igneof,status)
!
!*******************************************************************************
!
!! FTMBYT resets the i/o pointer to a certain byte number in a FITS file.
!
!
!  Discussion:
!
!    Subsequent read or write operations will begin at this point.
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       bytno   i  number of the byte to point to.
!       igneof  l  ignore end-of-file (107) error?
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, June, 1991
!       rewritten Feb, 1995

  integer iunit,bytno,status
  logical igneof

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne,pb
  parameter (nb = 20)
  parameter (ne = 512)
  parameter (pb = 20)
  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer lbuff,record,byten

  if (status > 0)then
          return
  else if (bytno < 0)then
!               error: negative byte number
          status=304
  else
          lbuff=bufnum(iunit)

!               calculate the record number and byte offset to move to
          record=bytno/reclen(lbuff)+1
          byten=mod(bytno,reclen(lbuff))

          if (record /= recnum(currnt(lbuff)))then
!                       not the current record, so load the new record;
                  call ftldrc(iunit,record,igneof,status)
          end if
          bytnum(lbuff)=byten
  end if

  return
end
subroutine ftpcbf(ounit,nchar,cbuff,status)
!
!*******************************************************************************
!
!! FTPCBF "Puts the Character BuFfer".
!
!
!       copy input buffer of characters to the output character buffer.
!
!  Parameters:
!
!       ounit   i  Fortran output unit number
!       nchar   i  number of characters in the string
!       cbuff   c  input character string
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!       modified Feb 1995

  integer ounit,nchar,status
  character ( len = * ) cbuff

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne,pb
  parameter (nb = 20)
  parameter (ne = 512)
  parameter (pb = 20)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)

!       have to use separate character arrays because of compiler limitations
  character ( len = 2880 ) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
  common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer lbuff,pbuff,buflen,lastb,nleft,in1,nbyt,nrec

  if (status > 0)return

  lbuff=bufnum(ounit)
  buflen=reclen(lbuff)

  if (nchar < 0)then
!               error: negative number of bytes to write
          status=306
          return
  else if (.not. wrmode(lbuff))then
!           don't have write access to this file
      status=112
      return
  end if

!       lastb   = position of last byte read from input buffer
!       nleft   = number of bytes left in the input buffer
!       in1     = position of first byte remaining in the input buffer
!       nbyt    = number of bytes to transfer from input to output

  nleft=nchar
  in1=1

!       find the number of bytes that will fit in output buffer
200     pbuff=currnt(lbuff)
  lastb=bytnum(lbuff)
  nbyt=min(nleft,buflen-lastb)
  if (nbyt > 0)then
!           append the input buffer to the output physical buffer
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
      19,20)pbuff

!               if got here, then pbuff is out of range
          status=101
          return

1               b1(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
2               b2(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
3               b3(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
4               b4(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
5               b5(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
6               b6(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
7               b7(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
8               b8(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
9               b9(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
10              b10(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
11              b11(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
12              b12(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
13              b13(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
14              b14(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
15              b15(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
16              b16(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
17              b17(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
18              b18(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
19              b19(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
20              b20(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)

100         modify(pbuff)=.true.
      bytnum(lbuff)=bytnum(lbuff)+nbyt
      in1=in1+nbyt
      nleft=nleft-nbyt
  end if

!       process more bytes, if any
  if (nleft > 0)then
    nrec=recnum(pbuff)+1

    if (nleft > buflen)then
!           first, flush any current buffers to disk
      call ftflsh(lbuff,status)

!           write whole blocks directly to the FITS file by-passing buffers
150         write(ounit,rec=nrec,err=900)cbuff(in1:in1+buflen-1)
      in1=in1+buflen
      nleft=nleft-buflen
      bytnum(lbuff)=bytnum(lbuff)+buflen
      nrec=nrec+1
      if (nleft > buflen)go to 150

!           Save maximum record written, for comparison in ftread
      maxrec(lbuff) = max(maxrec(lbuff), nrec-1)
    end if

!         load the next file record into a physical buffer
    call ftldrc(ounit,nrec,.true.,status)
    if (status > 0)return
    go to 200
  end if
  return

!       come here if there was a disk write error of some sort
900     status=106

  return
end
subroutine ftpcbo(ounit,gsize,ngroup,offset,cbuff,status)
!
!*******************************************************************************
!
!! FTPCBO "Puts the Character BuFfer with Offsets".
!
!
!       copy input buffer of characters to the output character buffer.
!
!  Parameters:
!
!       ounit   i  Fortran output unit number
!       gsize   i  size of each group of bytes
!       ngroup  i  number of groups to write
!       offset  i  size of gap between groups
!       cbuff   c  input character string
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1996

  integer ounit,gsize,ngroup,offset,status
  character ( len = * ) cbuff

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne,pb
  parameter (nb = 20)
  parameter (ne = 512)
  parameter (pb = 20)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)

!       have to use separate character arrays because of compiler limitations
  character ( len = 2880 ) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
  common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer lbuff,pbuff,buflen,lastb,nleft,in1,nbyt
  integer i,bytno,record,oldrec,incre

  if (status > 0)return

  lbuff=bufnum(ounit)

  if (.not. wrmode(lbuff))then
!           don't have write access to this file
      status=112
      return
  end if

  buflen=reclen(lbuff)
  pbuff =currnt(lbuff)
  oldrec=recnum(pbuff)
!       lastb = position of last byte read or written in FITS buffer
  lastb =bytnum(lbuff)
  bytno =(oldrec-1) * buflen + lastb
!       in1   = position of first byte remaining in the input buffer
  in1   =1
  incre =gsize+offset
  nbyt  = 0

  do i = 1,ngroup

!           nleft   = number of bytes left in the input buffer
      nleft=gsize
!           nbyt    = number of bytes to transfer from input to output
      nbyt=min(nleft,buflen-lastb)
      if (nbyt == 0)go to 300

200         continue
!           append the input buffer to the output physical buffer
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
         19,20)pbuff

!               if got here, then pbuff is out of range
          status=101
          return

1               b1(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
2               b2(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
3               b3(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
4               b4(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
5               b5(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
6               b6(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
7               b7(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
8               b8(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
9               b9(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
10              b10(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
11              b11(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
12              b12(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
13              b13(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
14              b14(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
15              b15(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
16              b16(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
17              b17(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
18              b18(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
19              b19(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
          go to 100
20              b20(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)

100         in1=in1+nbyt
      nleft=nleft-nbyt

!           process more bytes, if any
300         continue
      if (nleft > 0)then
!               entire group did not fit in the buffer
!               load the next file record into a physical buffer
          oldrec=oldrec+1
          modify(pbuff)=.true.
          call ftldrc(ounit,oldrec,.true.,status)
          if (status > 0)return
          pbuff=currnt(lbuff)
          lastb=0
          nbyt=nleft
          go to 200
      end if

      if (i /= ngroup)then
!               move to the position of the next group
          bytno=bytno+incre
          record=bytno/buflen+1
          lastb=mod(bytno,buflen)

          if (record /= oldrec)then
!                   not the current record, so load the new record;
              modify(pbuff)=.true.
              call ftldrc(ounit,record,.true.,status)
              if (status > 0)return
              oldrec=record
              pbuff=currnt(lbuff)
          end if
      end if
  end do

  modify(pbuff)=.true.
  bytnum(lbuff)=lastb+nbyt

  return
end
subroutine ftgcbf(iunit,nchar,array,status)
!
!*******************************************************************************
!
!! FTGCBF "Gets the Character BuFfer".
!
!
!       read NCHAR characters from the character buffer.
!
!  Parameters:
!
!       iunit   i  Fortran unit number for reading from disk
!       nchar   i  number of characters to read
!       array   c  output character string
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, June 1991
!       modified Feb 1995

  integer iunit,nchar,status
  character ( len = * ) array

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne,pb
  parameter (nb = 20)
  parameter (ne = 512)
  parameter (pb = 20)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)

!       have to use separate character arrays because of compiler limitations
  character ( len = 2880 ) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
  common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer nleft,nbyt,lastb,in1,lbuff,pbuff,buflen,nrec,ios,i

  if (status > 0)return

  if (nchar < 0)then
!               error: negative number of bytes to read
          status=306
          return
  end if

  lbuff=bufnum(iunit)
  buflen=reclen(lbuff)

!       lastb   = position of last byte read from input buffer
!       nleft   = number of bytes left in the input buffer
!       in1     = position of first byte remaining in the input buffer
!       nbyt    = number of bytes to transfer from input to output

  nleft=nchar
  in1=1

!       find the number of remaining bytes that can be read from buffer
200     pbuff=currnt(lbuff)
  lastb=bytnum(lbuff)
  nbyt=min(nleft,buflen-lastb)

!       get characters from the physical buffer to the output string
  if (nbyt > 0)then

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
      19,20)pbuff

!               if got here, then pbuff is out of range
          status=101
          return

1               array(in1:in1+nbyt-1)=b1(lastb+1:lastb+nbyt)
          go to 100
2               array(in1:in1+nbyt-1)=b2(lastb+1:lastb+nbyt)
          go to 100
3               array(in1:in1+nbyt-1)=b3(lastb+1:lastb+nbyt)
          go to 100
4               array(in1:in1+nbyt-1)=b4(lastb+1:lastb+nbyt)
          go to 100

!  The SUN F90 compiler gives a segmentation fault on the following
!  statement when executing testprog, while reading a complex (C) column
!  when using the Linux F90 routines (fitsf90_nag.f).
5               array(in1:in1+nbyt-1)=b5(lastb+1:lastb+nbyt)
          go to 100
6               array(in1:in1+nbyt-1)=b6(lastb+1:lastb+nbyt)
          go to 100
7               array(in1:in1+nbyt-1)=b7(lastb+1:lastb+nbyt)
          go to 100
8               array(in1:in1+nbyt-1)=b8(lastb+1:lastb+nbyt)
          go to 100
9               array(in1:in1+nbyt-1)=b9(lastb+1:lastb+nbyt)
          go to 100
10              array(in1:in1+nbyt-1)=b10(lastb+1:lastb+nbyt)
          go to 100
11              array(in1:in1+nbyt-1)=b11(lastb+1:lastb+nbyt)
          go to 100
12              array(in1:in1+nbyt-1)=b12(lastb+1:lastb+nbyt)
          go to 100
13              array(in1:in1+nbyt-1)=b13(lastb+1:lastb+nbyt)
          go to 100
14              array(in1:in1+nbyt-1)=b14(lastb+1:lastb+nbyt)
          go to 100
15              array(in1:in1+nbyt-1)=b15(lastb+1:lastb+nbyt)
          go to 100
16              array(in1:in1+nbyt-1)=b16(lastb+1:lastb+nbyt)
          go to 100
17              array(in1:in1+nbyt-1)=b17(lastb+1:lastb+nbyt)
          go to 100
18              array(in1:in1+nbyt-1)=b18(lastb+1:lastb+nbyt)
          go to 100
19              array(in1:in1+nbyt-1)=b19(lastb+1:lastb+nbyt)
          go to 100
20              array(in1:in1+nbyt-1)=b20(lastb+1:lastb+nbyt)

100         bytnum(lbuff)=bytnum(lbuff)+nbyt
      in1=in1+nbyt
      nleft=nleft-nbyt
  end if

!       process more bytes, if any
  if (nleft > 0)then
    nrec=recnum(pbuff)+1

150       continue

    if (nleft > buflen)then
!           read whole blocks directly from the FITS file by-passing buffers

!           test if desired record exists before trying to read it
      if (nrec + nleft/buflen - 1 > maxrec(lbuff)) then
!               record doesn't exist, so return EOF error
          status=107
          return
      end if

!           check if record is already loaded in one of the physical buffers
!           must read it from buffer since it may have been modified
      do i=1,maxbuf
         if (logbuf(i) == lbuff .and. recnum(i) == nrec)then
!                 found the desired record; don't have to read it
            go to 170
         end if
      end do

!           record not already loaded in buffer, so read it from disk
      read(iunit,rec=nrec,iostat=ios)array(in1:in1+buflen-1)

      if (ios /= 0)then
!               assume that this error indicates an end of file condition
          status=107
          return
      end if

      bytnum(lbuff)=bytnum(lbuff)+buflen
      in1=in1+buflen
      nleft=nleft-buflen
      nrec=nrec+1
      go to 150
    end if

!         load the next file record into a physical buffer
170       call ftldrc(iunit,nrec,.false.,status)
    if (status > 0)return
    go to 200
  end if

  return
end
subroutine ftgcbo(iunit,gsize,ngroup,offset,array,status)
!
!*******************************************************************************
!
!! FTGCBO "Gets the Character BuFfer with Offsets".
!
!
!       read characters from the character buffer.
!
!  Parameters:
!
!       iunit   i  Fortran output unit number
!       gsize   i  size of each group of bytes
!       ngroup  i  number of groups to read
!       offset  i  size of gap between groups
!       array   c  output character string
!       status  i  output error status (0 = ok)
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1996

  integer iunit,gsize,ngroup,offset,status
  character ( len = * ) array

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne,pb
  parameter (nb = 20)
  parameter (ne = 512)
  parameter (pb = 20)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)

!       have to use separate character arrays because of compiler limitations
  character ( len = 2880 ) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
  common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer lbuff,pbuff,buflen,lastb,nleft,in1,nbyt
  integer i,bytno,record,oldrec,incre

  if (status > 0)return

  lbuff =bufnum(iunit)
  buflen=reclen(lbuff)
  pbuff =currnt(lbuff)
  oldrec=recnum(pbuff)
!       lastb = position of last byte read from input buffer
  lastb =bytnum(lbuff)
  bytno =(oldrec-1) * buflen + lastb
!       in1   = position of first byte remaining in the input buffer
  in1   =1
  nbyt  =0
  incre =gsize+offset

  do i=1,ngroup

!           nleft   = number of bytes left in the input buffer
      nleft=gsize
!           nbyt    = number of bytes to transfer from input to output
      nbyt=min(nleft,buflen-lastb)
      if (nbyt == 0)go to 300

200         continue
!           get characters from the physical buffer to the output string
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
         19,20)pbuff

!               if got here, then pbuff is out of range
          status=101
          return

1               array(in1:in1+nbyt-1)=b1(lastb+1:lastb+nbyt)
          go to 100
2               array(in1:in1+nbyt-1)=b2(lastb+1:lastb+nbyt)
          go to 100
3               array(in1:in1+nbyt-1)=b3(lastb+1:lastb+nbyt)
          go to 100
4               array(in1:in1+nbyt-1)=b4(lastb+1:lastb+nbyt)
          go to 100
5               array(in1:in1+nbyt-1)=b5(lastb+1:lastb+nbyt)
          go to 100
6               array(in1:in1+nbyt-1)=b6(lastb+1:lastb+nbyt)
          go to 100
7               array(in1:in1+nbyt-1)=b7(lastb+1:lastb+nbyt)
          go to 100
8               array(in1:in1+nbyt-1)=b8(lastb+1:lastb+nbyt)
          go to 100
9               array(in1:in1+nbyt-1)=b9(lastb+1:lastb+nbyt)
          go to 100
10              array(in1:in1+nbyt-1)=b10(lastb+1:lastb+nbyt)
          go to 100
11              array(in1:in1+nbyt-1)=b11(lastb+1:lastb+nbyt)
          go to 100
12              array(in1:in1+nbyt-1)=b12(lastb+1:lastb+nbyt)
          go to 100
13              array(in1:in1+nbyt-1)=b13(lastb+1:lastb+nbyt)
          go to 100
14              array(in1:in1+nbyt-1)=b14(lastb+1:lastb+nbyt)
          go to 100
15              array(in1:in1+nbyt-1)=b15(lastb+1:lastb+nbyt)
          go to 100
16              array(in1:in1+nbyt-1)=b16(lastb+1:lastb+nbyt)
          go to 100
17              array(in1:in1+nbyt-1)=b17(lastb+1:lastb+nbyt)
          go to 100
18              array(in1:in1+nbyt-1)=b18(lastb+1:lastb+nbyt)
          go to 100
19              array(in1:in1+nbyt-1)=b19(lastb+1:lastb+nbyt)
          go to 100
20              array(in1:in1+nbyt-1)=b20(lastb+1:lastb+nbyt)

100         in1=in1+nbyt
      nleft=nleft-nbyt

!           process more bytes, if any
300         continue
      if (nleft > 0)then
!               load the next file record into a physical buffer
          oldrec=oldrec+1
          call ftldrc(iunit,oldrec,.false.,status)
          if (status > 0)return
          pbuff=currnt(lbuff)
          lastb=0
          nbyt=nleft
          go to 200
      end if

      if (i /= ngroup)then
!               move to the position of the next group
          bytno=bytno+incre
          record=bytno/buflen+1
          lastb=mod(bytno,buflen)

          if (record /= oldrec)then
!                   not the current record, so load the new record;
              call ftldrc(iunit,record,.false.,status)
              if (status > 0)return
              oldrec=record
              pbuff=currnt(lbuff)
          end if
      end if

  end do

  bytnum(lbuff)=lastb+nbyt

  return
end
subroutine ftgrsz(iunit,nrows,status)
!
!*******************************************************************************
!
!! FTGRSZ returns the optimal number of rows to move at one time.
!
!
!  This is the optimal value for the number of rows that should be
!  read or written at one time in a binary table for maximum efficiency.
!       Accessing more rows than this may cause excessive flushing and
!       rereading of buffers to/from disk.
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       nrows   i  optimal number of rows to access
!       status  i  output error status
!
!       written by Wm Pence, HEASARC/GSFC, December, 1996

  integer iunit,nrows,status

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,nf,ne,pb
  parameter (nb = 20)
  parameter (nf = 3000)
  parameter (ne = 512)
  parameter (pb = 20)

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

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer ibuff, ii, jj, unique, nfiles
!
!  There are pb internal buffers available each reclen(nb) bytes long
!
  ibuff=bufnum(iunit)
!
!  if HDU structure is not defined then scan the header keywords
!
  if (dtstrt(ibuff) < 0)call ftrdef(iunit,status)
  if (status > 0)return
!
!  Determine how many different FITS files are currently open
!
  nfiles = 0
  do ii = 1,nb
     unique = 1
     do jj = 1, ii-1
       if (buflun(ii) <= 0 .or. buflun(ii) == buflun(jj))then
         unique = 0
         exit
       end if
     end do

     if (unique == 1)nfiles=nfiles+1
  end do
!
!  one buffer (at least) is always allocated to each open file.
!  assume record size is 2880 bytes (not necessarily true on Vax)
!
  nrows = ((pb - nfiles) * 2880) / max(1,rowlen(ibuff))
  nrows = max(1, nrows)

  return
end
subroutine ftflsh(lbuff,status)
!
!*******************************************************************************
!
!! FTFLSH flushes any modified buffers associated with lbuff to disk.
!
!
!       Make the contents of the buffers undefined.
!
!  Parameters:
!
!       lbuff   i  logical buffer assocaiated with this file
!       status  i  output error status

  integer lbuff,status

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,pb
  parameter (nb = 20)
  parameter (pb = 20)

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer ounit,rlen,i

!       ignore input status and flush buffers regardless of status value

  ounit=buflun(lbuff)
  rlen=reclen(lbuff)

!       find any buffer associated with this file
  do i=1,maxbuf
      if (logbuf(i) == lbuff)then
          if (modify(i))then
!                   write the modified buffer to disk
              call ftwrit(ounit,recnum(i),rlen,i,status)
              modify(i)=.false.
          end if

!               erase the association of this buffer with the file
          logbuf(i)=0
          recnum(i)=0
      end if
  end do

  return
end
subroutine ftldrc(iunit,nrec,igneof,status)
!
!*******************************************************************************
!
!! FTLDRC loads a specified record from a file into a physical buffer.
!
!
!       The record is only loaded if it is not already loaded.  Reset all
!       pointers to make this the new current record for that file.
!       Update ages of all the physical buffers.
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       nrec    i  direct access file record number to be loaded
!       igneof  l  ignore end of file error (107)?
!       status  i  output error status

  integer iunit,nrec,status
  logical igneof

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne,pb
  parameter (nb = 20)
  parameter (ne = 512)
  parameter (pb = 20)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)
  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)
  integer compid
  common/ftcpid/compid
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer i,lbuff,pbuff,ounit,olen,orec,a1,tstat

  if (status > 0)return

  lbuff=bufnum(iunit)

!       check if record is already loaded in one of the physical buffers
  do i=1,maxbuf
      if (logbuf(i) == lbuff .and. recnum(i) == nrec)then
!               found the desired record; don't have to read it
          pbuff=i
          go to 20
      end if
  end do
!
!       the record is not already loaded, so we have to read it from disk.
!       First, decide which physical buffer into which to read it.
!
  call ftwhbf(lbuff,pbuff)

  if (modify(pbuff))then
!           old buffer has been modified, so we have to flush it to disk
      ounit=buflun(logbuf(pbuff))
      olen=reclen(logbuf(pbuff))
      orec=recnum(pbuff)
      call ftwrit(ounit,orec,olen,pbuff,status)
      modify(pbuff)=.false.
  end if

!       now read the record into the physical buffer
  olen=reclen(lbuff)
  tstat=0
  call ftread(iunit,nrec,olen,pbuff,tstat)

  if (.not. igneof .and. tstat == 107)then
!           return if hit EOF and told not to ignore it
      status=107
      return
  else if (tstat == 107)then
!           apparently hit end of file

      if (.not. wrmode(lbuff))then
!               just return if we don't have write access to the file
          return
      else
!               fill the new buffer with the desired value
          if (hdutyp(lbuff) == 1)then
!                   ASCII table: fill buffer with blanks
              call ftflbl(pbuff)
          else if (compid >= -1)then
!                   initialize buffer = 0 (except on Cray machines)
              call ftflzr(pbuff)
          else
!                   call special routine for Cray machines, since words
!                   are twice as long (integers are 8-bytes long)
              call ftzrcr(pbuff)
          end if

!               mark the new record as having been modified
          modify(pbuff)=.true.
      end if
  end if

!       define log. buffer and the record number contained in the phys. buffer
  logbuf(pbuff)=lbuff
  recnum(pbuff)=nrec

20      continue
!       this is now the current buffer for this logical buffer
  currnt(lbuff)=pbuff
  bytnum(lbuff)=0

!       find the current position of the buffer in the age index
  do i=1,maxbuf
      if (pindex(i) == pbuff)then
         a1=i
         go to 35
      end if
  end do

35      continue
!       rebuild the indices giving the chronological ordering of the buffers
  do i=a1,maxbuf-1
          pindex(i)=pindex(i+1)
  end do
!       this buffer is now the youngest (= last in the index)
  pindex(maxbuf)=pbuff

  return
end
subroutine ftflzr(pbuff)
!
!*******************************************************************************
!
!! FTFLZR initializea the common block buffer efficiently with zeros.
!
!
!  This routine should not be used on Cray computers.
!
!  Parameters:
!
!       pbuff  i  number of the physical buffer to initialize

  integer pbuff,i
  integer pb
  parameter (pb = 20)
  double precision buff
  common /ftbuff/buff(360,pb)

  do i=1,360
    buff(i,pbuff)=0.
  end do

  return
end
subroutine ftzrcr(pbuff)
!
!*******************************************************************************
!
!! FTZRCR initializes the common block buffer with zeros.
!
!
!       This routine is reserved for Cray computers.
!
!  Parameters:
!
!       pbuff  i  number of the physical buffer to initialize

  integer pbuff,i
  integer pb
  parameter (pb = 20)
  integer buff,dummy
  common /ftbuff/buff(360,pb),dummy(360,pb)

!       The dummy array was added to the common block to eliminate
!       compiler warnings on other platforms besides Cray.  On all
!       other machines, an integer is 4 bytes long, hence the dummy
!       array is needed to increase the size of the common block to
!       the same length as other routines where the array is declared
!       as character ( len = 2880 ) buff(pb).

!       This will cause compiler warnings on the CRAY, but these
!       can be ignored since the dummy array is never used.

  do i=1,360
    buff(i,pbuff)=0
  end do

  return
end
subroutine ftflbl(pbuff)
!
!*******************************************************************************
!
!! FTFLBL initializes the common block buffer efficiently with blanks.
!
!
!  Parameters:
!
!       pbuff  i  number of the physical buffer to initialize

  integer pbuff
  integer pb
  parameter (pb = 20)
  character ( len = 2880 ) cbuff
  common /ftbuff/cbuff(pb)

  cbuff(pbuff) = ' '

  return
end
subroutine ftwhbf ( ilbuff, pbuff )
!
!*******************************************************************************
!
!! FTWHBF decides which physical buffer to use to load in a new record.
!
!
!  Parameters:
!
!    Input, integer ILBUFF, the logical buffer number of the record to
!    be loaded.
!
!    Output, integer PBUFF, the physical buffer that should be used.
!
  integer, parameter :: nb = 20
  integer, parameter :: pb = 20
!
  integer i
  integer ilbuff
  integer lbuff
  integer num
  integer pbuff
!
  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

  integer maxbuf,logbuf,recnum,pindex
  logical modify
  common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb), &
  pindex(pb)
!
!  Search from the oldest to the youngest for an unlocked record.
!
!  NUM is the number of the next oldest buffer.
!
!  LBUFF is the logical buffer associated with this physical buffer.
!
  do i = 1, maxbuf

    num = pindex(i)
    lbuff = logbuf(num)
!
!  Is the physical buffer a current buffer (i.e., locked)?
!
    if ( lbuff == 0 ) then
      pbuff = num
      return
    else if ( currnt(lbuff) /= num ) then
      pbuff = num
      return
    end if

  end do
!
!  If all the buffers are locked, we have to reuse the current one.
!
  pbuff = currnt(ilbuff)

  return
end
subroutine ftwrit(ounit,nrec,length,pbuff,status)
!
!*******************************************************************************
!
!! FTWRIT writes a physical buffer to the disk file.
!
!
!  Parameters:
!
!       ounit   i  Fortran unit number to write to
!       nrec    i  number of the file record to write
!       length  i  number of bytes to write
!       pbuff   i  number of the physical buffer to write from
!       status  i  output error status

  integer ounit,nrec,length,pbuff,status

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne
  parameter (nb = 20)
  parameter (ne = 512)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

!       have to use separate character arrays because of compiler limitations
  character ( len = 2880 ) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
  common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer ibuff

!       Note: performance testing on a SUN showed that writing a
!       c*2888 string was MUCH (11x) faster than writing a C*1(2880) array
!       with write(...)(b1(i),i=1,2880).  It was also 2-3 times faster
!       than if the array was declared as a double and written with
!       write(...)(darray(i),i=1,360).  The VAX took about the same
!       time for all 3 different ways to write the bytes.

  ibuff=bufnum(ounit)
!       Save maximum record written, for comparison in ftread
  maxrec(ibuff) = max(maxrec(ibuff), nrec)

  go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
  19,20)pbuff
!
!       if got here, then pbuff is out of range.
!
  status=101
  return

1       write(ounit,rec=nrec,err=900)b1(1:length)
  return
2       write(ounit,rec=nrec,err=900)b2(1:length)
  return
3       write(ounit,rec=nrec,err=900)b3(1:length)
  return
4       write(ounit,rec=nrec,err=900)b4(1:length)
  return
5       write(ounit,rec=nrec,err=900)b5(1:length)
  return
6       write(ounit,rec=nrec,err=900)b6(1:length)
  return
7       write(ounit,rec=nrec,err=900)b7(1:length)
  return
8       write(ounit,rec=nrec,err=900)b8(1:length)
  return
9       write(ounit,rec=nrec,err=900)b9(1:length)
  return
10      write(ounit,rec=nrec,err=900)b10(1:length)
  return
11      write(ounit,rec=nrec,err=900)b11(1:length)
  return
12      write(ounit,rec=nrec,err=900)b12(1:length)
  return
13      write(ounit,rec=nrec,err=900)b13(1:length)
  return
14      write(ounit,rec=nrec,err=900)b14(1:length)
  return
15      write(ounit,rec=nrec,err=900)b15(1:length)
  return
16      write(ounit,rec=nrec,err=900)b16(1:length)
  return
17      write(ounit,rec=nrec,err=900)b17(1:length)
  return
18      write(ounit,rec=nrec,err=900)b18(1:length)
  return
19      write(ounit,rec=nrec,err=900)b19(1:length)
  return
20      write(ounit,rec=nrec,err=900)b20(1:length)
  return

!       come here if there was a disk write error of some sort
900     status=106

  return
end
subroutine ftread(iunit,nrec,length,pbuff,status)
!
!*******************************************************************************
!
!! FTREAD reads a disk file record into a physical buffer.
!
!
!  Parameters:
!
!       iunit   i  Fortran unit number to read from
!       nrec    i  number of the file record to read
!       length  i  number of bytes to read
!       pbuff   i  number of the physical buffer to read into
!       status  i  output error status

  integer iunit,nrec,length,pbuff,status

!       COMMON BLOCK DEFINITIONS:--------------------------------------------
  integer nb,ne
  parameter (nb = 20)
  parameter (ne = 512)

  integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
  integer nxtfld
  logical wrmode
  common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb), &
  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

  integer buflun,currnt,reclen,bytnum,maxrec
  common/ftlbuf/buflun(nb),currnt(nb),reclen(nb), &
  bytnum(nb),maxrec(nb)

!       have to use separate character arrays because of compiler limitations
  character ( len = 2880 ) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
  common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14, &
  b15,b16,b17,b18,b19,b20
!       END OF COMMON BLOCK DEFINITIONS-----------------------------------

  integer ibuff,ios

!       test if desired record exists before trying to read it
  ibuff=bufnum(iunit)
  if (nrec > maxrec(ibuff)) then
!             record doesn't exist, so return EOF error
        status=107
        return
  end if

  go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, &
  19,20)pbuff

!       if got here, then pbuff is out of range
  status=101
  return

1       read(iunit,rec=nrec,iostat=ios)b1(1:length)
  go to 100
2       read(iunit,rec=nrec,iostat=ios)b2(1:length)
  go to 100
3       read(iunit,rec=nrec,iostat=ios)b3(1:length)
  go to 100
4       read(iunit,rec=nrec,iostat=ios)b4(1:length)
  go to 100
5       read(iunit,rec=nrec,iostat=ios)b5(1:length)
  go to 100
6       read(iunit,rec=nrec,iostat=ios)b6(1:length)
  go to 100
7       read(iunit,rec=nrec,iostat=ios)b7(1:length)
  go to 100
8       read(iunit,rec=nrec,iostat=ios)b8(1:length)
  go to 100
9       read(iunit,rec=nrec,iostat=ios)b9(1:length)
  go to 100
10      read(iunit,rec=nrec,iostat=ios)b10(1:length)
  go to 100
11      read(iunit,rec=nrec,iostat=ios)b11(1:length)
  go to 100
12      read(iunit,rec=nrec,iostat=ios)b12(1:length)
  go to 100
13      read(iunit,rec=nrec,iostat=ios)b13(1:length)
  go to 100
14      read(iunit,rec=nrec,iostat=ios)b14(1:length)
  go to 100
15      read(iunit,rec=nrec,iostat=ios)b15(1:length)
  go to 100
16      read(iunit,rec=nrec,iostat=ios)b16(1:length)
  go to 100
17      read(iunit,rec=nrec,iostat=ios)b17(1:length)
  go to 100
18      read(iunit,rec=nrec,iostat=ios)b18(1:length)
  go to 100
19      read(iunit,rec=nrec,iostat=ios)b19(1:length)
  go to 100
20      read(iunit,rec=nrec,iostat=ios)b20(1:length)

100     continue
  if (ios /= 0)then
!               assume that this error indicates an end of file condition
          status=107
  end if

  return
end
subroutine ftpi2b(ounit,nvals,incre,i2vals,status)
!
!*******************************************************************************
!
!! FTPI2B writes an array of Integer*2 bytes to the output FITS file.
!
!
!       Does any required translation from internal machine format to FITS.
!
!  Parameters:
!
!       ounit   i  fortran unit number
!       nvals   i  number of pixels in the i2vals array
!       incre   i  byte increment between values
!       i2vals  i*2 array of input integer*2 values
!       status  i  output error status

  integer nvals,incre,ounit,status,offset
  integer*2 i2vals(nvals)


  integer*2 temp(4)
  integer ierr,cray2ieg,remain
  integer compid
  common/ftcpid/compid
  character ctemp
!
  if (compid == 0)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid == -1)then
!           SUN F90 compiler maps I*2 -> I*4; have to pack bytes
      call ftpki2(i2vals,nvals,ctemp)
  else if (compid >= 1)then
!           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
      call ftswby(i2vals,nvals)
  else
!           must be a CRAY
!           convert from cray I*8 to IEEE I*2
!           there is a bug in cray2ieg if the number of values to convert
!           is  1 less than a  multiple of 4 2-byte words.  (3, 7, 11, etc)
      remain=nvals-nvals/4*4
      if (remain == 3)then
        if (nvals > 3)then
          ierr= cray2ieg(7,nvals-3,i2vals,0,i2vals,1,' ')
        end if
        temp(3)=i2vals(nvals)
        temp(2)=i2vals(nvals-1)
        temp(1)=i2vals(nvals-2)
        ierr=cray2ieg(7,4,i2vals(nvals/4+1),0,temp,1,' ')
      else
        ierr=cray2ieg(7,nvals,i2vals,0,i2vals,1,' ')
      end if
  end if

  if (incre <= 2)then
          call ftpbyt(ounit,nvals*2,i2vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-2
          call ftpbyo(ounit,2,nvals,offset,i2vals,status)
  end if

  return
end
subroutine ftpi4b(ounit,nvals,incre,i4vals,status)
!
!*******************************************************************************
!
!! FTPI4B writes an array of Integer*4 bytes to the output FITS file.
!
!
!       Does any required translation from internal machine format to FITS.
!
!  Parameters:
!
!       ounit   i  fortran unit number
!       nvals   i  number of pixels in the i4vals array
!       incre   i  byte increment between values
!       i4vals  i  array of input integer*4 values
!       status  i  output error status

  integer nvals,incre,ounit,status,offset
  integer i4vals(nvals)


  integer compid
  common/ftcpid/compid

  integer cray2ieg,neven,ierr

  if (compid == 0 .or. compid == -1)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid >= 1)then
!           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
      call ftswi4(i4vals,nvals)
  else
!         must be a CRAY
!         there is a bug in cray2ieg if the number of values to convert
!         is not a multiple of 8 bytes.
    neven=nvals/2*2
    if (neven > 0)then
        ierr= cray2ieg(1,neven,i4vals,0,i4vals,1,' ')
    end if

    if (neven /= nvals)then
!           have to do the remaining odd word separately
      ierr= cray2ieg(1,1,i4vals(nvals/2+1),0,i4vals(nvals),1,' ')
    end if
  end if

  if (incre <= 4)then
          call ftpbyt(ounit,nvals*4,i4vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-4
          call ftpbyo(ounit,4,nvals,offset,i4vals,status)
  end if

  return
end
subroutine ftpr4b(ounit,nvals,incre,r4vals,status)
!
!*******************************************************************************
!
!! FTPR4B writes an array of Real*4 bytes to the output FITS file.
!
!
!       Does any required translation from internal machine format to FITS.
!
!  Parameters:
!
!       ounit   i  fortran unit number
!       nvals   i  number of pixels in the r4vals array
!       incre   i  byte increment between values
!       r4vals  r  array of input real*4 values
!       status  i  output error status

  integer nvals,incre,ounit,status,offset
  real r4vals(nvals)

  integer compid
  common/ftcpid/compid

  integer i,neven,ierr,cray2ieg

  if (compid == 0 .or. compid == -1)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid == 1)then
!           little endian machine (e.g. DEC or PC) must be byte swapped
      call ftswi4(r4vals,nvals)
  else if (compid >= 2)then
!           convert from VAX format to IEEE
      do i=1,nvals
              r4vals(i)=r4vals(i)*0.25
      end do
      call ftswby(r4vals,nvals*2)
  else
!         must be a CRAY
!         there is a bug in cray2ieg if the number of values to convert
!         is not a multiple of 8 bytes.
    neven=nvals/2*2
    ierr= cray2ieg(2,neven,r4vals,0,r4vals,1,' ')
    if (neven /= nvals)then
!           have to do the remaining odd word separately
      ierr= cray2ieg(2,1,r4vals(nvals/2+1),0,r4vals(nvals),1,' ')
    end if
  end if

  if (incre <= 4)then
          call ftpbyt(ounit,nvals*4,r4vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-4
          call ftpbyo(ounit,4,nvals,offset,r4vals,status)
  end if

  return
end
subroutine ftpr8b(ounit,nvals,incre,r8vals,status)
!
!*******************************************************************************
!
!! FTPR8B writes an array of Real*8 bytes to the output FITS file.
!
!
!       Does any required translation from internal machine format to FITS.
!
!  Parameters:
!
!       ounit   i  fortran unit number
!       nvals   i  number of pixels in the r4vals array
!       incre   i  byte increment between values
!       r8vals  d  array of input real*8 values
!       status  i  output error status

  integer nvals,incre,ounit,status,offset
  double precision r8vals(nvals)


  integer compid
  common/ftcpid/compid

  integer i,ierr,cray2ieg

  if (compid == 0 .or. compid == -1)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid == 1)then
!           little endian machine (e.g. DEC or PC) must be byte swapped
      call ftswi8(r8vals,nvals)
  else if (compid == 2)then
!           convert from VAX format to IEEE
      call ieevpd(r8vals,r8vals,nvals)
  else if (compid == 3)then
!           convert from Alpha VMS format to IEEE
      do 5 i=1,nvals
              r8vals(i)=r8vals(i)*0.25
5           continue
      call ftswby(r8vals,nvals*4)
  else
!           must be a CRAY
      ierr= cray2ieg(3,nvals,r8vals,0,r8vals,1,' ')
  end if

  if (incre <= 8)then
      call ftpbyt(ounit,nvals*8,r8vals,status)
  else
!           offset is the number of bytes to move between each value
      offset=incre-8
      call ftpbyo(ounit,8,nvals,offset,r8vals,status)
  end if

  return
end
subroutine ftgi2b(iunit,nvals,incre,i2vals,status)
!
!*******************************************************************************
!
!! FTGI2B reads an array of Integer*2 bytes from the input FITS file.
!
!
!       Does any required translation from FITS to internal machine format
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       nvals   i  number of pixels to read
!       incre   i  byte increment between values
!       i2vals  i*2 output array of integer*2 values
!       status  i  output error status

  integer nvals,iunit,incre,status,offset
  integer*2 i2vals(nvals)


  integer compid
  common/ftcpid/compid

  integer ierr,ieg2cray,i,nloop,fpixel,ntodo
  integer*2 temp(4)
  character ctemp

  if (incre <= 2)then
          call ftgbyt(iunit,nvals*2,i2vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-2
          call ftgbyo(iunit,2,nvals,offset,i2vals,status)
  end if

  if (compid == 0)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid == -1)then
!           SUN F90 compiler maps I*2 -> I*4; have to unpack bytes
      call ftupi2(i2vals,nvals,ctemp)
  else if (compid >= 1)then
!           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
      call ftswby(i2vals,nvals)
  else
!           must be a CRAY
!           convert from IEEE I*2 to cray I*8

!           have to use a temporary array if nvals = 2 or 3
      if (nvals <= 3)then
        ierr=ieg2cray(7,nvals,i2vals,0,temp,1,' ')
        do 5 i=1,nvals
            i2vals(i)=temp(i)
5             continue
      else

!             have to work backwards, so as to not overwrite the input data
        nloop=(nvals-1)/4+1
        fpixel = (nloop*4)-3
        ntodo=nvals-(nloop-1)*4
        do i=nloop,1,-1
          ierr=ieg2cray(7,ntodo,i2vals(i),0,i2vals(fpixel),1,' ')
          fpixel=fpixel-4
          ntodo=4
        end do
      end if
  end if

  return
end
subroutine ftgi4b(iunit,nvals,incre,i4vals,status)
!
!*******************************************************************************
!
!! FTGI4B reads an array of Integer*4 bytes from the input FITS file.
!
!
!       Does any required translation from FITS to internal machine format.
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       nvals   i  number of pixels to read
!       incre   i  byte increment between values
!       i4vals  i  output array of integer values
!       status  i  output error status

  integer nvals,iunit,incre,status,offset
  integer i4vals(nvals)


  integer ierr,ieg2cray,nloop,fpixel,ntodo,i

  integer compid
  common/ftcpid/compid

  if (incre <= 4)then
          call ftgbyt(iunit,nvals*4,i4vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-4
          call ftgbyo(iunit,4,nvals,offset,i4vals,status)
  end if

  if (compid == 0 .or. compid == -1)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid >= 1)then
!           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
      call ftswi4(i4vals,nvals)
  else
!           must be a CRAY
!           convert from IEEE I*4 to cray I*8
!           have to work backwards, so as to not overwrite the input data

      nloop=(nvals+1)/2
      fpixel = (nloop*2)-1
      ntodo=nvals-(nloop-1)*2
      do i=nloop,1,-1
          ierr=ieg2cray(1,ntodo,i4vals(i),0,i4vals(fpixel),1,' ')
          fpixel=fpixel-2
          ntodo=2
      end do
  end if

  return
end
subroutine ftgr4b(iunit,nvals,incre,r4vals,status)
!
!*******************************************************************************
!
!! FTGR4B reads an array of Real*4 bytes from the input FITS file.
!
!
!       Does any required translation from FITS to internal machine format.
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       nvals   i  number of pixels to read
!       incre   i  byte increment between values
!       r4vals  r  output array of real*4 values
!       status  i  output error status

  integer nvals,iunit,incre,status,offset
  real r4vals(nvals)


  integer ierr,ieg2cray,nloop,fpixel,ntodo,i

  integer compid
  common/ftcpid/compid

  if (incre <= 4)then
          call ftgbyt(iunit,nvals*4,r4vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-4
          call ftgbyo(iunit,4,nvals,offset,r4vals,status)
  end if

  if (compid == 0 .or. compid == -1)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid == 1)then
!           little endian machine (e.g. DEC or PC) must be byte swapped
      call ftswi4(r4vals,nvals)
  else if (compid >= 2)then
!           convert the values from IEEE format to VAX floating point format
!           first, swap the bytes
      call ftswby(r4vals,nvals*2)
!           then test for IEEE special values and multiply value by 4.0
      call ftr4vx(r4vals,r4vals,nvals)
  else
!           must be a CRAY
!           convert from IEEE R*4 to cray R*8
!           have to work backwards, so as to not overwrite the input data

      nloop=(nvals+1)/2
      fpixel = (nloop*2)-1
      ntodo=nvals-(nloop-1)*2
      do i=nloop,1,-1
          ierr=ieg2cray(2,ntodo,r4vals(i),0,r4vals(fpixel),1,' ')
          fpixel=fpixel-2
          ntodo=2
      end do
  end if

  return
end
subroutine ftgr8b(iunit,nvals,incre,r8vals,status)
!
!*******************************************************************************
!
!! FTGR8B reads an array of Real*8 bytes from the input FITS file.
!
!
!       Does any required translation from FITS to internal machine format.
!
!  Parameters:
!
!       iunit   i  fortran unit number
!       nvals   i  number of pixels to read
!       incre   i  byte increment between values
!       r8vals  d  output array of real*8 values
!       status  i  output error status

  integer nvals,iunit,incre,status,offset
  double precision r8vals(nvals)


  integer compid
  common/ftcpid/compid

  integer ierr,ieg2cray,nloop,fpixel,ntodo,i

  if (incre <= 8)then
          call ftgbyt(iunit,nvals*8,r8vals,status)
  else
!               offset is the number of bytes to move between each value
          offset=incre-8
          call ftgbyo(iunit,8,nvals,offset,r8vals,status)
  end if

  if (compid == 0 .or. compid == -1)then
!           big endian machine (e.g., SUN) doesn't need byte swapping
  else if (compid == 1)then
!           little endian machine (e.g. DEC or PC) must be byte swapped
      call ftswi8(r8vals,nvals)
  else if (compid == 2)then
!           convert the values from IEEE format to VAX double (D) format
      call ieevud(r8vals,r8vals,nvals)
  else if (compid == 3)then
!           convert the values from IEEE format to VMS double (G) format
!           first, swap the bytes
      call ftswby(r8vals,nvals*4)
!           then test for IEEE special values and multiply value by 4.0
      call ftr8vx(r8vals,r8vals,r8vals,nvals)
  else
!           must be a CRAY
!           convert from IEEE R*8 to cray R*16
!           have to work backwards, so as to not overwrite the input data

      nloop=(nvals+1)/2
      fpixel = (nloop*2)-1
      ntodo=nvals-(nloop-1)*2
      do i=nloop,1,-1
          ierr=ieg2cray(3,ntodo,r8vals(i),0,r8vals(fpixel),1,' ')
          fpixel=fpixel-2
          ntodo=2
      end do
  end if

  return
end
subroutine ftswby(buffer,npix)
!
!*******************************************************************************
!
!! FTSWBY swaps pairs of bytes to convert I*2 values to/from IEEE.
!
!
!  Parameters:
!
!       buffer  i  input buffer of 16-bit words to be swapped
!       npix    i  number of 16-bit words to be swapped
!
!       written by Wm Pence, HEASARC/GSFC, February 1991

  integer npix,k
  logical*1 buffer(npix*2),temp

  do k=2,npix*2,2
     temp = buffer(k-1)
     buffer(k-1) = buffer(k)
     buffer(k) = temp
  end do

  return
end
subroutine ftswi4(buffer,npix)
!
!*******************************************************************************
!
!! FTSWI4 swaps 4 successive bytes in values to convert to/from IEEE.
!
!
!               if the original byte order is: 1 2 3 4
!               the output order is:           4 3 2 1
!
!  Parameters:
!
!       buffer  i  input buffer of bytes to be swapped
!       npix    i  number of 4 bytes-words to be swapped
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1996

!  The equivalence statements in this algorithm cause errors with
!  the Linux NAG F90 compiler.
!        integer npix, i
!        integer buffer(npix)
!        integer*4 temp1,temp2
!        logical*1 l1(4), l2(4)
!        equivalence (temp1,l1)
!        equivalence (temp2,l2)
!
!        do i=1,npix
!            temp1=buffer(i)
!            l2(4)=l1(1)
!            l2(3)=l1(2)
!            l2(2)=l1(3)
!            l2(1)=l1(4)
!            buffer(i)=temp2
!        end do

!       alternate swapping algorithm, Wm Pence, April 97
  integer npix,i4
  logical*1 buffer(npix*4),temp
!
!  Swap bytes 1 and 4, 2 and 3.
!
  do i4=4,npix*4,4

     temp=buffer(i4-3)
     buffer(i4-3)=buffer(i4)
     buffer(i4)=temp

     temp=buffer(i4-2)
     buffer(i4-2)=buffer(i4-1)
     buffer(i4-1)=temp

  end do

  return
end
subroutine ftswi8(buffer,npix)
!
!*******************************************************************************
!
!! FTSWI8 swaps 8 successive bytes in values to convert to/from IEEE.
!
!
!               if the original byte order is: 1 2 3 4 5 6 7 8
!               the output order is:           8 7 6 5 4 3 2 1
!
!  Parameters:
!
!       buffer  i  input buffer of bytes to be swapped
!       npix    i  number of 8 bytes-words to be swapped
!
!       written by Wm Pence, HEASARC/GSFC, Dec 1996

!  The equivalence statements in this algorithm cause errors with
!  the Linux NAG F90 compiler.
!        integer npix, i
!        double precision buffer(npix),temp1,temp2
!        logical*1 l1(8), l2(8)
!        equivalence (temp1,l1)
!        equivalence (temp2,l2)
!
!        do i=1,npix
!            temp1=buffer(i)
!            l2(8)=l1(1)
!            l2(7)=l1(2)
!            l2(6)=l1(3)
!            l2(5)=l1(4)
!            l2(4)=l1(5)
!            l2(3)=l1(6)
!            l2(2)=l1(7)
!            l2(1)=l1(8)
!            buffer(i)=temp2
!        end do

!       alternate swapping algorithm, Wm Pence, April 97
  integer npix,i8
  logical*1 buffer(npix*8),temp
!
!  Swap bytes 1 and 8, 2 and 7, 3 and 6, 4 and 5.
!
  do i8=8,npix*8,8

     temp=buffer(i8-7)
     buffer(i8-7)=buffer(i8)
     buffer(i8)=temp

     temp=buffer(i8-6)
     buffer(i8-6)=buffer(i8-1)
     buffer(i8-1)=temp

     temp=buffer(i8-5)
     buffer(i8-5)=buffer(i8-2)
     buffer(i8-2)=temp

     temp=buffer(i8-4)
     buffer(i8-4)=buffer(i8-3)
     buffer(i8-3)=temp

  end do

  return
end
subroutine ftr4vx(r4vals,i2vals,nvals)
!
!*******************************************************************************
!
!! FTR4VX converts IEEE 32-bit floating point numbers to VAX floating point.
!
!
!       This routine is only called on Vax and Alpha VMS systems.
!
!  Parameters:
!
  real r4vals(*)
  integer*2 i2vals(*)
  integer nvals,i,j

  j=1
  do i=1,nvals
!           test for NaNs (treat +/- infinity the same as a NaN)
      if (i2vals(j) >= 32640 .or. (i2vals(j) < 0 .and. &
           i2vals(j) >= -128))then
           i2vals(j)=-1
           i2vals(j+1)=-1

!           set underflows and -0 (8000000 hex) = to zero
      else if (i2vals(j) <= -32641 .or. (i2vals(j) >= 0 .and. &
           i2vals(j) <= 127))then
           r4vals(i)=0.0
      else
!
!  Must be a real number, so multiply by 4.0 to convert to Vax.
!
           r4vals(i)=r4vals(i)*4.0
      end if
      j=j+2
  end do

  return
end
subroutine ftr8vx(r8vals,i4vals,i2vals,nvals)
!
!*******************************************************************************
!
!! FTR8VX converts IEEE 32-bit floating point numbers to VAX floating point.
!
!
!       This routine is only called on VAX computers.
!
!  Parameters:
!
  double precision r8vals(*)
  integer*2 i2vals(*)
  integer i4vals(*)
  integer nvals,i,j,k

  j=1
  k=1
  do i=1,nvals
!           test for NaNs (treat +/- infinity the same as a NaN)
      if (i2vals(j) >= 32752 .or. (i2vals(j) < 0 .and. &
           i2vals(j) >= -16))then
           i4vals(k)  =-1
           i4vals(k+1)=-1

!           set underflows and -0 (8000000 hex) = to zero
      else if (i2vals(j) <= -32753 .or. (i2vals(j) >= 0 .and. &
           i2vals(j) <= 15))then
           r8vals(i)=0.0
      else
!                Must be a real number, so multiply by 4.0 to convert to Vax
           r8vals(i)=r8vals(i)*4.0
      end if
      j=j+4
      k=k+2
  end do

  return
end
subroutine ftpki2(i2vals,nvals,temp)
!
!*******************************************************************************
!
!! FTPKI2 packs an array of 4-byte integers into sequence of 2-byte integers.
!
!
!       This routine is only currently used on the SUN Solaris F90
!       which does not directly support integer*2 variables and instead
!       maps them into integer*4 variables.
!
!  Parameters:
!
  integer nvals,ii,jj
  integer*2 temp
  character i2vals(nvals*4)

  jj = 2
  do ii = 4,nvals*4,4
       i2vals(jj-1) = i2vals(ii - 1)
       i2vals(jj)   = i2vals(ii)
       jj = jj +2
  end do

  return
end
subroutine ftupi2(i2vals,nvals,temp)
!
!*******************************************************************************
!
!! FTUPI2 unpacks an array of 2-byte integers into 4-integers.
!
!
!       This routine is only currently used on the SUN Solaris F90
!       which does not directly support integer*2 variables and instead
!       maps them into integer*4 variables.
!
!  Parameters:
!
  integer nvals,ii,jj
  integer*2 temp
  character i2vals(nvals*4),zero,neg,pos

  zero = char(0)
!       neg is use to extend the sign bit if the value is negative
  neg  = char(255)
  pos  = char(127)

  jj=nvals*4
  do ii = nvals*2,2,-2
       i2vals(jj)=i2vals(ii)
       i2vals(jj-1)=i2vals(ii - 1)

!            fill in the 2 most-significant bytes
       if (i2vals(jj-1) <= pos)then
           i2vals(jj-2)=zero
           i2vals(jj-3)=zero
       else
           i2vals(jj-2)=neg
           i2vals(jj-3)=neg
       end if
       jj=jj-4
  end do

  return
end


