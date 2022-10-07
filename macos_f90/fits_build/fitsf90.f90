! 
! A support routine for writing a FITS file
! downloaded from http://orion.math.iastate.edu/burkardt/g_src/fitsio/fitsio.html
!

subroutine ftopnf ( funit, fname, oldnew, rwmode, block, size, status )
!
!*******************************************************************************
!
!! FTOPNF creates or opens a new file.
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
!       size    i  min size of file, in bytes
!       status  i  returned error status (0=ok)
!
  integer block
  character ( len = 2880 ) dummy
  character ( len = 9 ) fmode
  character ( len = * ) fname
  character ( len = 3 ) fstat
  integer funit
  integer lenrec
  integer oldnew
  integer rwmode
  integer size
  integer status
!
!  Some compilers require that dummy be initialized!
!
  dummy = ' '
!
!  Simply return the default block size.
!  No easy way to determine the size of the file, 
!  so just set value to a huge number
!
  if (oldnew == 0)then

    fstat='OLD'
    block=1
    size=1000000000
!
!  Create a new file.
!
  else

    fstat='NEW'
    size = 0

  end if

  if ( rwmode == 0 ) then
    fmode = 'READ'
  else
    fmode = 'READWRITE'
  end if
!
!  determine the record length needed to read or write 2880-byte records
!  (could be 720 words or 2880 bytes, depending on the compiler)
!
  inquire(iolength=lenrec) dummy
  open(unit=funit,file=fname,status=fstat,err=900, action=fmode, &
       recl=lenrec,form='UNFORMATTED',access='DIRECT')

  return
!
!       error opening file:
!
900     status=104 + oldnew
end
subroutine ftgsdt(dd,mm,yy,status)
!
!*******************************************************************************
!
!! FTGSDT gets the current date from the system
!
!       dd      i  day of the month (1-31)
!       mm      i  month of the year (1-12)
!       yy      i  last 2 digits of the year (1992 = 92, 2001 = 01)
!
!     This routine is not valid after the year 1999; must be modified
!     before the year 2000.

  integer dd,mm,yy,status
  integer iarray(8)

  if (status > 0)return

  call date_and_time(values=iarray)

  dd=iarray(3)
  mm=iarray(2)
  yy=mod(iarray(1),100)
end
subroutine ftpbyt(ounit,nbytes,array,status)
!
!*******************************************************************************
!
!! FTPBYT writes string of data bytes to output buffer.
!
!
!  Author:
!
!    William Pence, HEASARC/GSFC
!
!       ounit   i  fortran unit number
!       nbytes  i  number of bytes
!       array   i  integer array
!       status  i  output error status
!
  integer nbytes,ounit,status
  integer array((nbytes-1) / 4 + 1)
  character temp

!       Trick the compiler into passing the integer array to a routine
!       that is expecting a character string argument.
!       This is very ugly, but I/O efficiency is of primary importance.

  call ftpxbf(ounit,nbytes,array,temp,status)
end
subroutine ftpxbf(ounit,nbytes,array,temp,status)
!
!*******************************************************************************
!
!! FTPXBF ???
!
  integer temp
  integer nbytes,ounit,status
  character ( len = nbytes) array

!       array is now declared as a string!
  call ftpcbf(ounit,nbytes,array,status) 
end
subroutine ftpbyo(ounit,gsize,ngroup,offset,array,status)
!
!*******************************************************************************
!
!! FTPBYO "Puts Bytes with Offsets".
!
!
!  Author:
!
!    William Pence, HEASARC/GSFC
!
!       copy input buffer of bytes to the output character buffer.
!
!       ounit   i  Fortran output unit number
!       gsize   i  size of each group of bytes
!       ngroup  i  number of groups to write
!       offset  i  size of gap between groups
!       array   i  input array of bytes
!       status  i  output error status (0 = ok)
!
  integer ounit,gsize,ngroup,offset,status
  integer array(*)
  character temp

!       Trick the compiler into passing the integer array to a routine
!       that is expecting a character string argument.
!       This is very ugly, but I/O efficiency is of primary importance.

  call ftpxbo(ounit,gsize,ngroup,offset,array,temp,status)
end
subroutine ftpxbo(ounit,gsize,ngroup,offset,array,temp,status)
!
!*******************************************************************************
!
!! FTPXBO ???
!
  integer temp
  integer ounit,gsize,ngroup,offset,status
  character ( len = gsize*ngroup) array

!       array is now declared as a string!
  call ftpcbo(ounit,gsize,ngroup,offset,array,status)
end
subroutine ftgbyt(iunit,nbytes,array,status)
!
!*******************************************************************************
!
!! FTGBYT reads string of data bytes from input buffer.
!
!
!  Author:
!
!    William Pence, HEASARC/GSFC
!
!       iunit   i  fortran unit number
!       nbytes  i  number of bytes
!       array   i  integer array
!       status  i  output error status
!
  integer nbytes,iunit,status
  integer array(*)
  character temp

!       Trick the compiler into passing the integer array to a routine
!       that is expecting a character string argument.
!       This is very ugly, but I/O efficiency is of primary importance.

  call ftgxbf(iunit,nbytes,array,temp,status)

end
subroutine ftgxbf(iunit,nbytes,array,temp,status)
!
!*******************************************************************************
!
!! FTGXBF ???
!
  integer temp
  integer nbytes,iunit,status
  character ( len = nbytes ) array

!       array is now declared as a string!
  call ftgcbf(iunit,nbytes,array,status) 
end
subroutine ftgbyo(iunit,gsize,ngroup,offset,array,status)
!
!*******************************************************************************
!
!! GTGBYO "Gets BYtes with Offsets".
!
!       read bytes from the character buffer.
!
!  Author:
!
!    William Pence, HEASARC/GSFC
!
!       iunit   i  Fortran output unit number
!       gsize   i  size of each group of bytes
!       ngroup  i  number of groups to read
!       offset  i  size of gap between groups
!       array   i  output array of bytes
!       status  i  output error status (0 = ok)
!
  integer iunit,gsize,ngroup,offset,status
  integer array(*)
  character temp

!       Trick the compiler into passing the integer array to a routine
!       that is expecting a character string argument.
!       This is very ugly, but I/O efficiency is of primary importance.

  call ftgxbo(iunit,gsize,ngroup,offset,array,temp,status)
end
subroutine ftgxbo(iunit,gsize,ngroup,offset,array,temp,status)
!
!*******************************************************************************
!
!! FTGXBO ???
!
  integer temp
  integer iunit,gsize,ngroup,offset,status
  character ( len = gsize*ngroup ) array

!       array is now declared as a string!
  call ftgcbo(iunit,gsize,ngroup,offset,array,status)
end
subroutine ieevud(dbl1, dbl2, int)
!
!*******************************************************************************
!
!! IEEVUD ???
!
  double precision dbl1,dbl2
  integer int
!       dummy routine; only used on Vax VMS computers
end
subroutine ieevpd(dbl1, dbl2, int)
!
!*******************************************************************************
!
!! IEEVPD ???
!
  double precision dbl1,dbl2
  integer int
!       dummy routine; only used on Vax VMS computers
end
function cray2ieg(i1,i2,i3,i4,i5,i6,s1)
!
!*******************************************************************************
!
!! CRAY2IEG ???
!
  integer cray2ieg
  integer i1,i2,i3,i4,i5,i6
  character s1
!       dummy routine; only used on Cray supercomputers
  cray2ieg=0
end
function ieg2cray(i1,i2,i3,i4,i5,i6,s1)
!
!*******************************************************************************
!
!! IEG2CRAY ???
!
  integer ieg2cray
  integer i1,i2,i3,i4,i5,i6
  character s1
!       dummy routine; only used on Cray supercomputers
  ieg2cray=0
end

