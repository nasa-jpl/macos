!
! Sample program to show fits image read and write
! downloaded from http://orion.math.iastate.edu/burkardt/g_src/fitsio/cookbook.f90
!

program cookbook
!
!*******************************************************************************
!
!! COOKBOOK calls all the example subroutines

  character(len=80) :: filname

  write ( *, * ) ' '
  write ( *, * ) 'COOKBOOK'
  write ( *, * ) '  Some simple FITSIO tests.'

  filname = 'cookbook_1.fit'
  call write_image(filname)

  if (.false.) then
  call write_ascii
  call write_bintable
  call copy_hdu
  call select_rows
  call read_header
  call read_image
  call read_table
  end if

  write ( *, * ) ' '
  write ( *, * ) 'COOKBOOK'
  write ( *, * ) '  Normal end of FITSIO tests.'

  stop
end


subroutine write_image(filename)
character (len=*) :: filename
!
!*******************************************************************************
!
!! WRITE_IMAGE creates a FITS primary array containing a 2-D image.
!
  integer status,unit,blocksize,bitpix,naxis,naxes(2)
  integer i,j,group,fpixel,nelements,array(300,200)
  !character ( len = 80 ) filename
  logical simple,extend
!
  write ( *, * ) ' '
  write ( *, * ) 'WRITE_IMAGE:'

  status=0
!
!  Name of the FITS file to be created:
!
  !filename = 'cookbook_1.fit'
!
!  Delete the file if it already exists, so we can then recreate it.
!
  call delete_file ( filename, status )
!
!  Get an unused Logical Unit Number to use to open the FITS file.
!
  call ftgiou ( unit, status )
!
!  Create the new empty FITS file.
!
  blocksize=1
  call ftinit(unit,filename,blocksize,status)
!
!  Initialize parameters about the FITS image (300 x 200 16-bit integers).
!
  simple=.true.
  !bitpix=16   ! original
  bitpix=32
  naxis=2
  naxes(1)=300
  naxes(2)=200
  extend=.true.
!
!  Write the required header keywords.
!
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!
!  Initialize the values in the image with a linear ramp function
!
  do j=1,naxes(2)
    do i=1,naxes(1)
      array(i,j)=i+j
    end do
  end do
!
!  Write the array to the FITS file.
!
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)
  !call ftpprj(unit,group,fpixel,nelements,array,status)  ! original
  call ftpprd(unit,group,fpixel,nelements,array,status)
!
!  Write another optional keyword to the header.
!
  call ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)
!
!  Close the file and free the unit number.
!
  call ftclos(unit, status)
  call ftfiou(unit, status)
!
!  Check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine write_ascii
!
!*******************************************************************************
!
!! WRITE_ASCIII creates an ASCII table containing 3 columns and 6 rows.
!
  integer status,unit,readwrite,blocksize,tfields,nrows,rowlen
  integer nspace,tbcol(3),diameter(6), colnum,frow,felem
  real density(6)
  character ( len = 40 ) filename
  character ( len = 16 ) extname
  character ( len = 16 ) ttype(3)
  character ( len = 16 ) tform(3)
  character ( len = 16 ) tunit(3)
  character ( len = 16 ) name(6)
  data ttype/'Name','Diameter','Density'/
  data tform/'A8','I6','F4.2'/
  data tunit/' ','km','g/cm'/
  data name/'Mercury','Venus','Earth','Mars','Jupiter','Saturn'/
  data diameter/4880,12112,12742,6800,143000,121000/
  data density/5.1,5.3,5.52,3.94,1.33,0.69/
!
  write ( *, * ) ' '
  write ( *, * ) 'WRITE_ASCII:'

  status=0
!
!  Name of the FITS file to append the ASCII table to:
!
  filename='cookbook_1.fit'
!
!  Get an unused Logical Unit Number to use to open the FITS file
!
  call ftgiou(unit,status)
!
!  open the FITS file, with write access
!
  readwrite=1
  call ftopen(unit,filename,readwrite,blocksize,status)
!
!  append a new empty extension onto the end of the primary array
!
  call ftcrhd(unit,status)
!
!  define parameters for the ASCII table (see the above data statements)
!
  tfields=3
  nrows=6
  extname='PLANETS_ASCII'
!     
!  calculate the starting position of each column, and the total row length
!
  nspace=1
  call ftgabc(tfields,tform,nspace,rowlen,tbcol,status)
!
!  write the required header parameters for the ASCII table
!
  call ftphtb(unit,rowlen,nrows,tfields,ttype,tbcol,tform,tunit, &
    extname,status)
!
!  write names to the first column, diameters to 2nd col., and density to 3rd
!
  frow=1
  felem=1
  colnum=1
  call ftpcls(unit,colnum,frow,felem,nrows,name,status)
  colnum=2
  call ftpclj(unit,colnum,frow,felem,nrows,diameter,status)  
  colnum=3
  call ftpcle(unit,colnum,frow,felem,nrows,density,status)  
!
!  close the FITS file and free the unit number
!
  call ftclos(unit, status)
  call ftfiou(unit, status)
!
!  check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine write_bintable
!
!*******************************************************************************
!
!! WRITE_BINTABLE creates a binary table containing 3 columns and 6 rows.
!
  integer status,unit,readwrite,blocksize,hdutype,tfields,nrows
  integer varidat,diameter(6), colnum,frow,felem
  real density(6)
  character ( len = 40 ) filename
  character ( len = 16 ) extname
  character ( len = 16 ) ttype(3)
  character ( len = 16 ) tform(3)
  character ( len = 16 ) tunit(3)
  character ( len = 16 ) name(6)
  data ttype/'Name','Diameter','Density'/
  data tform/'8A','1J','1E'/
  data tunit/' ','km','g/cm'/
  data name/'Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'/
  data diameter/6800,143000,121000,47000,45000,6000/
  data density/3.94,1.33,0.69,1.56,2.27,1.0/
!
  write ( *, * ) ' '
  write ( *, * ) 'WRITE_BINTABLE:'

  status=0
!
!  Name of the FITS file to append the ASCII table to:
!
  filename='cookbook_1.fit'
!
!  Get an unused Logical Unit Number to use to open the FITS file
!
  call ftgiou(unit,status)
!
!  open the FITS file, with write access
!
  readwrite=1
  call ftopen(unit,filename,readwrite,blocksize,status)
!
!  move to the last (2nd) HDU in the file
!
  call ftmahd(unit,2,hdutype,status)
!
!  append/create a new empty HDU onto the end of the file and move to it
!
  call ftcrhd(unit,status)
!
!  define parameters for the binary table (see the above data statements)
!
  tfields=3
  nrows=6
  extname='PLANETS_BINARY'
  varidat=0
!     
!  write the required header parameters for the binary table
!
  call ftphbn(unit,nrows,tfields,ttype,tform,tunit, &
    extname,varidat,status)
!
!  write names to the first column, diameters to 2nd col., and density to 3rd
!
  frow=1
  felem=1
  colnum=1
  call ftpcls(unit,colnum,frow,felem,nrows,name,status)
  colnum=2
  call ftpclj(unit,colnum,frow,felem,nrows,diameter,status)  
  colnum=3
  call ftpcle(unit,colnum,frow,felem,nrows,density,status)  
!
!  close the FITS file and free the unit number
!
  call ftclos(unit, status)
  call ftfiou(unit, status)
!
!  check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine copy_hdu
!
!*******************************************************************************
!
!! COPY_HDU copies the 1st and 3rd HDUs from the input file to a new FITS file.
!
  integer status,inunit,outunit,readwrite,blocksize,morekeys,hdutype
  character ( len = 40 ) infilename
  character ( len = 40 ) outfilename
!
  write ( *, * ) ' '
  write ( *, * ) 'COPY_HDU:'

  status=0
!
!  Name of the FITS files:
!
  infilename='cookbook_1.fit'
  outfilename='cookbook_2.fit'
!
!  Delete the file if it already exists, so we can then recreate it
!
  call delete_file(outfilename,status)
!
!  Get unused Logical Unit Numbers to use to open the FITS files
!
  call ftgiou(inunit,status)
  call ftgiou(outunit,status)
!
!  open the input FITS file, with readonly access
!
  readwrite=0
  call ftopen(inunit,infilename,readwrite,blocksize,status)
!
!  create the new empty FITS file with the standard block size
!
  call ftinit(outunit,outfilename,blocksize,status)
!
!  copy the primary array from the input file to the output file
!
  morekeys=0
  call ftcopy(inunit,outunit,morekeys,status)
!
!  append/create a new empty extension on the end of the output file
!
  call ftcrhd(outunit,status)
!
!  skip to the 3rd extension in the input file
!
  call ftmahd(inunit,3,hdutype,status)
!
!  copy this extension from the input file to the output file
!
  call ftcopy(inunit,outunit,morekeys,status)  
!
!  close the FITS file and free the unit numbers
!
  call ftclos(inunit, status)
  call ftclos(outunit, status)
  call ftfiou(-1, status)
!
!  check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine select_rows
!
!*******************************************************************************
!
!! SELECT_ROWS selects rows from an input table and copies to the output table.
!
  integer status,inunit,outunit,readwrite,blocksize,hdutype
  integer nkeys,nspace,naxes(2),nfound,colnum,frow,felem
  integer noutrows,irow,temp(100),i
  real nullval,density(6)
  character ( len = 40 ) infilename
  character ( len = 40 ) outfilename
  character ( len = 80 ) record
  logical exact,anynulls
!
  write ( *, * ) ' '
  write ( *, * ) 'SELECT_ROWS:'

  status=0
!
!  Names of the FITS files:
!
  infilename='cookbook_1.fit'
  outfilename='cookbook_2.fit'
!
!  Get unused Logical Unit Numbers to use to open the FITS files
!
  call ftgiou(inunit,status)
  call ftgiou(outunit,status)
!
!  open the FITS files, with the appropriate read/write access
!
  readwrite=0
  call ftopen(inunit,infilename,readwrite,blocksize,status)
  readwrite=1
  call ftopen(outunit,outfilename,readwrite,blocksize,status)
!
!  move to the 3rd HDU in the input file (a binary table in this case)
!
  call ftmahd(inunit,3,hdutype,status)
!
!  move to the last extension in the output file
!
  do while (status == 0)
      call ftmrhd(outunit,1,hdutype,status)
  end do

  if (status == 107)then
!
!  this is normal; it just means we hit the end of file
!
      status=0
      call ftcmsg
  end if
!
!  create a new empty extension in the output file
!
  call ftcrhd(outunit,status)
  call ftghsp(inunit,nkeys,nspace,status)
!
!  copy all the keywords from the input to the output extension
!
  do i=1,nkeys
      call ftgrec(inunit,i,record,status)
      call ftprec(outunit,record,status)
  end do
!
!  force FITSIO to read the output file keywords to define the data structure
!
  call ftrdef(outunit,status)
!
!  get the width of the table (in bytes) and the number of rows
!
  call ftgknj(inunit,'NAXIS',1,2,naxes,nfound,status)
!
!  find which column contains the DENSITY values
!
  exact=.false.
  call ftgcno(inunit,exact,'DENSITY',colnum,status)
!
!  read the DENSITY column values.
!
  frow=1
  felem=1
  nullval=-99.
  call ftgcve(inunit,colnum,frow,felem,naxes(2),nullval, &
    density,anynulls,status)
!
!  If the density is less than 3.0, copy the row to the output table
!
  noutrows=0
  do irow=1,naxes(2)
      if (density(irow) < 3.0)then
          noutrows=noutrows+1
          call ftgtbb(inunit,irow,1,naxes(1),temp,status)
          call ftptbb(outunit,noutrows,1,naxes(1),temp,status)
      end if
  end do
    
  call ftmkyj(outunit,'NAXIS2',noutrows,'&',status)
!
!  close the FITS file and free the unit numbers
!
  call ftclos(inunit, status)
  call ftclos(outunit, status)
  call ftfiou(-1, status)
!
!  check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine read_header
!
!*******************************************************************************
!
!! READ_HEADER prints the header keywords in all extensions of a FITS file.
!
  integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i
  character ( len = 80 ) filename
  character ( len = 80 ) record
!
  write ( *, * ) ' '
  write ( *, * ) 'READ_HEADER:'

  status=0
!
!  Get an unused Logical Unit Number to use to open the FITS file.
!
  call ftgiou(unit,status)
!
!  name of FITS file.
!
  filename='cookbook_1.fit'
!
!  open the FITS file, with read-only access.
!
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)

100   continue
!
!  Determine the number of keywords in the header
!
  call ftghsp(unit,nkeys,nspace,status)
!
!  Read each 80-character keyword record, and print it out
!
  do i = 1, nkeys
      call ftgrec(unit,i,record,status)
      print *,record
  end do
!
!  Print out and END record, and a blank line to mark the end of the header
!
  if (status == 0)then
      print *,'END'
      print *,' '
  end if
!
!  try moving to the next extension in the FITS file, if it exists.
!
  call ftmrhd(unit,1,hdutype,status)

  if (status == 0)then
!
!  success, so loop back and print out keywords in this extension
!
      go to 100

  else if (status == 107)then
!
!  hit end of file, so quit
!
      print *,'***** END OF FILE *****'
      status=0
      call ftcmsg
  end if
!
!     close the file, free the unit number, and exit
!
  call ftclos(unit, status)
  call ftfiou(unit, status)
!
!  check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine read_image
!
!*******************************************************************************
!
!! READ_IMAGE reads a FITS image and determines the pixel range.
!
  integer status,unit,readwrite,blocksize,naxes(2),nfound
  integer group,firstpix,nbuffer,npixels,i
  real datamin,datamax,nullval,buffer(100)
  logical anynull
  character ( len = 80 ) filename
!
  write ( *, * ) ' '
  write ( *, * ) 'READ_IMAGE:'

  status=0
!
!  Get an unused Logical Unit Number to use to open the FITS file.
!
  call ftgiou(unit,status)
!
!  open the FITS file previously created by WRITE_IMAGE.
!
  filename='cookbook_1.fit'
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)
!
!  determine the size of the image
!
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
!
!  check that it found both NAXIS1 and NAXIS2 keywords
!
  if (nfound /= 2)then
      print *,'READ_IMAGE failed to read the NAXISn keywords.'
      return
   end if
!
!  initialize variables.
!
  npixels=naxes(1)*naxes(2)
  group=1
  firstpix=1
  nullval=-999
  datamin=1.0E30
  datamax=-1.0E30

  do while (npixels > 0)
!
!  read up to 100 pixels at a time 
!
      nbuffer=min(100,npixels)
  
     call ftgpve(unit,group,firstpix,nbuffer,nullval,buffer,anynull,status)
!
!  find the min and max values
!
      do i=1,nbuffer
          datamin=min(datamin,buffer(i))
          datamax=max(datamax,buffer(i))
      end do
!
!  increment pointers and loop back to read the next group of pixels
!
      npixels=npixels-nbuffer
      firstpix=firstpix+nbuffer
  end do
!
!  print out the min and max values
!
  print *,'Min and max values in the image are:',datamin,datamax
!
!  close the file and free the unit number
!
  call ftclos(unit, status)
  call ftfiou(unit, status)
!
!  Check for any error, and if so print out error messages
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine read_table
!
!*******************************************************************************
!
!! READ_TABLE reads and prints data values from an ASCII or binary table.
!
  integer status,unit,readwrite,blocksize,hdutype,ntable
  integer felem,nelems,nullj,diameter,nfound,irow,colnum
  real nulle,density
  character ( len = 40 ) filename
  character nullstr
  character ( len = 8 ) name
  character ( len = 10 ) ttype(3)
  logical anynull
!
  write ( *, * ) ' '
  write ( *, * ) 'READ_TABLE:'

  status=0
!
!  Get an unused Logical Unit Number to use to open the FITS file.
!
  call ftgiou(unit,status)
!
!  Open the FITS file previously created by WRITE_IMAGE
!
  filename='cookbook_1.fit'
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)
!
!  Loop twice, first reading the ASCII table, then the binary table
!
  do ntable=1,2
!
!  Move to the next extension
!
      call ftmrhd(unit,1,hdutype,status)

      print *,' '
      if (hdutype == 1)then
          print *,'Extension ',ntable,' is an ASCII table.'
      else if (hdutype == 2)then
          print *,'Extension ',ntable,' is a binary table.'
      end if
!
!  Read the TTYPEn keywords, which give the names of the columns.
!
      call ftgkns(unit,'TTYPE',1,3,ttype,nfound,status)
      write(*,2000)ttype
2000      format(8x,3a10)
!
!  Read the data, one row at a time, and print them out
!
      felem=1
      nelems=1
      nullstr=' '
      nullj=0
      nulle=0.
      do irow=1,6
          colnum=1
          call ftgcvs(unit,colnum,irow,felem,nelems,nullstr,name, &
            anynull,status)
          colnum=2
          call ftgcvj(unit,colnum,irow,felem,nelems,nullj,diameter, &
            anynull,status)
          colnum=3
          call ftgcve(unit,colnum,irow,felem,nelems,nulle,density, &
            anynull,status)
          write(*,2001)irow,name,diameter,density
2001          format(i4,a10,i10,f10.2)
      end do
  end do
!
!  Close the file and free the unit number
!
  call ftclos(unit, status)
  call ftfiou(unit, status)
!
!  Check for any error, and if so print out error messages.
!
  if (status > 0) then
    call print_error(status)
  end if

  return
end


subroutine print_error(status)
!
!*******************************************************************************
!
!! PRINT_ERROR prints out the FITSIO error messages to the user.
!
  integer status
  character ( len = 30 ) errtext
  character ( len = 80 ) errmessage
!
!  Check if status is OK (no error); if so, simply return.
!
  if (status <= 0) then
    return
  end if
!
!  Get the text string which describes the error
!
  call ftgerr(status,errtext)
  print *,'FITSIO Error Status =',status,': ',errtext
!
!  Read and print out all the error messages on the FITSIO stack
!
  call ftgmsg(errmessage)
  do while (errmessage .ne. ' ')
      print *,errmessage
      call ftgmsg(errmessage)
  end do

  return
end


subroutine delete_file(filename,status)
!
!*******************************************************************************
!
!! DELETE_FILE deletes a FITS file.
!
  integer status,unit,blocksize
  character ( len = * ) filename
!
!  Simply return if status is greater than zero.
!
  if (status > 0) then
    return
  end if
!
!  Get an unused Logical Unit Number to use to open the FITS file
!
  call ftgiou ( unit, status )
!
!  Try to open the file, to see if it exists
!
  call ftopen ( unit, filename, 1, blocksize, status )

  if ( status == 0 ) then
!
!  File was opened;  so now delete it 
!
    call ftdelt(unit,status)

  else if (status == 103)then
!
!  File doesn't exist, so just reset status to zero and clear errors
!
    status=0
    call ftcmsg

  else
!
!  There was some other error opening the file; delete the file anyway
!
      status=0
      call ftcmsg
      call ftdelt(unit,status)
  end if
!
!  Free the unit number for later reuse.
!
  call ftfiou(unit, status)

  return
end


