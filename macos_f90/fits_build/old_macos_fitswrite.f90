!
! Routine to write FITS files
!

Subroutine fitswrite(filname,image,img_sz)
  Implicit none
  Character(len=80) :: filname
  Real*8 :: image(img_sz,img_sz)  ! input image in 64-bit floating point
  Integer :: image_gs(img_sz,img_sz)
  Real :: gmin,gmax
  Integer :: img_sz,ndim(2),nscale,i,j
  Real :: glvl

  nscale=2**8-1  ! 0-255 grayscale 
 
  ! Scale image into an integer array of 16-bit grayscale
  ndim(1)=size(image,1); ndim(2)=size(image,2)
  print*,'  ndim(1),ndim(2) =',ndim(1),ndim(2)

  gmin=minval(image); gmax=maxval(image)
  glvl=(gmax-gmin)/nscale
 
  Do i=1,ndim(1)
    Do j=1,ndim(2)
      image_gs(i,j)=(image(i,j)-gmin)/glvl
    End Do
  End Do

  !Open(unit=12,file='Opd.gs',status='replace')
  !Do j=1,img_sz
  !  Write(12,*) image_gs(1:img_sz,j)
  !End Do;
  !Close(12)

  write ( *, * ) ' '
  write ( *, * ) '  Writing FITS image to file ', filname
  write ( *, * ) ' '

  call write_image(filname,image_gs,img_sz)

  write ( *, * ) ' '
  write ( *, * ) '  Normal end of FITS image write'
  return
end


subroutine write_image(filename,array,arr_sz)
character (len=80) :: filename
integer :: array(arr_sz,arr_sz)
integer :: arr_sz
!
!*******************************************************************************
!
!! WRITE_IMAGE dumps a FITS primary array containing a 2D image to a file.
!
  integer status,unit,blocksize,bitpix,naxis,naxes(2)
  integer i,j,group,fpixel,nelements
  logical simple,extend
 
  write ( *, * ) ' '
  write ( *, * ) '  WRITE_IMAGE:'

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
!
  simple=.true.; bitpix=16; naxis=2
  naxes(1)=size(array,1); naxes(2)=size(array,2)
  !print*,'  ** naxes(1),naxes(2) =',naxes(1),naxes(2)
  extend=.true.

!
!  Write the required header keywords.
!
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!
!  Write the array to the FITS file.
!
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)
  call ftpprj(unit,group,fpixel,nelements,array,status)


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
