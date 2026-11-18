!
! macos_fitswirte.f90
! Routines to write FITS files
!

Subroutine fitswrite(filname,image,img_sz)
  Implicit none
  Character(len=80) :: filname
  Real*8 :: image(img_sz,img_sz)  ! input image in 64-bit floating point
  Integer :: img_sz

  write ( *, * ) ' '
  write ( *, * ) '  Writing FITS image to file ', filname

  Call write_image(filname,image,img_sz)

  write ( *, * ) '  Normal end of FITS image write'
End Subroutine


Subroutine write_image(filename,array,arr_sz)
Character (len=80) :: filename
Integer :: arr_sz
Real*8 :: array(arr_sz,arr_sz)
!
!*******************************************************************************
!
!! WRITE_IMAGE dumps a FITS primary array containing a 2D image to a file.
!
  integer status,unit,blocksize,bitpix,naxis,naxes(2)
  integer i,j,group,fpixel,nelements
  logical simple,extend
 
  !write ( *, * ) ' '
  !write ( *, * ) '  WRITE_IMAGE:'
  status=0
 
!
!  Delete the file if it already exists, so we can then recreate it.
!
  call delete_file (filename,status)
!
!  Get an unused Logical Unit Number to use to open the FITS file.
!
  call ftgiou (unit,status)
!
!  Create the new empty FITS file.
!
  blocksize=1
  call ftinit(unit,filename,blocksize,status)

!
!  Initialize parameters about the FITS image (300 x 200 16-bit integers).
!  bitpix value of -64 is for double floating type.
  simple=.true.; bitpix=-64; naxis=2
  naxes(1)=arr_sz; naxes(2)=arr_sz
  !print*,'  ** naxes(1),naxes(2) =',naxes(1),naxes(2)
  extend=.true.

!
!  Write the required header keywords.
!
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!
!  Write the array to the FITS file.
!
  group=1; fpixel=1
  nelements=naxes(1)*naxes(2)
  !call ftpprj(unit,group,fpixel,nelements,array,status)
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
End Subroutine

