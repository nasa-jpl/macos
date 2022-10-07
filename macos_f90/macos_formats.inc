C
C  macos_formats.inc
C

C  2. FORMAT statements
 106	FORMAT(1P,' Exit pupil distance from previous element=',d11.4/
     &	'           Closest proximity of chief rays=',d11.4)
 107	FORMAT(1P,' Exit pupil finder results:'/
     &	'   Old f=',D17.9,'   New f=',D17.9/
     &	'   Old z=',D17.9,'   New z=',D17.9/
     &	' Old psi=',D17.9,' New psi=',D17.9,2(/9x,D17.9,9x,D17.9)/
     &	' Old Vpt=',D17.9,' New Vpt=',D17.9,2(/9x,D17.9,9x,D17.9))
 109	FORMAT(1P,' Reference surface recalculation results:'/
     &	'   Old f=',D17.9,'   New f=',D17.9/
     &	'   Old z=',D17.9,'   New z=',D17.9/
     &	' Old psi=',D17.9,' New psi=',D17.9,2(/9x,D17.9,9x,D17.9)/
     &	' Old Vpt=',D17.9,' New Vpt=',D17.9,2(/9x,D17.9,9x,D17.9))
 110    FORMAT(' Partial of ray ',i4,' at Element ',i4,' (',a8,
     &  ') to ray at Element ',i4,' (',a8,')')
 111    FORMAT(' Partial of ray ',i4,' at Element ',i4,' (',a8,
     &  ') to Element ',i4,' (',a8,') perturbations')
 112    FORMAT(' Computing ray number ',i4,' out of ',i4)
 113    FORMAT(/' Storing data in output file'/)
 114    FORMAT(' Computing ray matrices')

#if 0
 115    FORMAT(1P,/' Chief Ray from Element ',i4,' (',a8,
     &  ') to Element ',i4,' (',a8,'):'/
     &  '  Starting point=',3D17.9/'       End point=',3D17.9/
     &  '       Direction=',3D17.9/'          Length=',D17.9)
#endif

 115    FORMAT(1P,/' Chief Ray from Element ',i4,' (',a8,
     &  ') to Element ',i4,' (',a8,'):'/
     &  '  Starting point=',3D22.14/'       End point=',3D22.14/
     &  '       Direction=',3D22.14/'          Length=',D22.14)

 116    FORMAT(' Compute time was',f10.4,' sec')
 117	FORMAT(' Tracing ',i6,' rays...')
 118	FORMAT(' Tracing ',i6,' rays and BUILDing linear model...')
 119	FORMAT(' Tracing ',i6,' rays but not BUILDing linear model.'/
     &	' BUILD limited to ',i6,' rays.  MOD npts and recompute.')
 120    FORMAT(1x,A,3D17.9)
 150    FORMAT('  Must load optical data first: use NEW or OLD.')
 151    FORMAT(1x,A,I3,A,D17.9)
 152    FORMAT(' Must BUILD first.')
 157    FORMAT(' Must LPERTURB before plotting linear model results.')
 158    FORMAT('  Must run OPD command first')
 159    FORMAT(' ** Computing OPD Zernike fit coeffs (Noll ordered) **')
 160    FORMAT(' ** Removing Zernike modes from OPD (Noll ordered) **')
 161    FORMAT(
     &  ' ** Computing OPD with selected Zernike modes (Noll ordered)')
 162    FORMAT(A,I3,A,3f14.9)
 163    FORMAT(A,f12.4,' ',A)
 164    FORMAT(2x,A)
 165    FORMAT(A,$)
 166    FORMAT(x,E17.9,$)
 167    FORMAT(x,A,D17.9,x,A)
 168    FORMAT(x,A,I3,A2,A10,A2,A,D17.9,x,A)
 169    FORMAT(x,A,I4,A2,A10,A2,A,$)
!169    FORMAT(x,A,I4,A,$)

 201    FORMAT(2x,A,3E17.9)
 202    FORMAT(2x,A,E17.9)
 203    FORMAT(A,i3,A)
 204    FORMAT(2x,6D17.9)

 207    FORMAT(' Status for file ',A,':'/
     &  '     Current element for OPD calculations=',i4/
     &  ' Current element for WF/Spot calculations=',i4/
     &  '         Obscuration option for ray-trace=',i4/
     &  '                     Current WF plot type=',A/
     &  '                           Composed image=',L2/
     &  '   Pixel location set by chief ray option=',L2/
     &  '                             Polarization=',L2)
 217    FORMAT(' Map of rays at input aperture:')
 218    FORMAT(2x,26i3)
 219    FORMAT(' Map of segments at input aperture:')
 221    FORMAT(2x,26i3)
 222    FORMAT(2x,20i4)
 223    FORMAT(2x,16i5)
 224    FORMAT(2x,13i6)
 225    FORMAT(2x,11i7)
 301	FORMAT(' MACOS>')
 302	FORMAT(' Command not understood.  Please reenter...')

!303	FORMAT(/' Modeling and Analysis for Controlled Optical Systems'/
!    &	' MACOS Version 3.2, 2006'//
!    &	' This run is limited to',i4,' surfaces,',i4,' segments,',

 303    FORMAT(/'  This run is limited to',i4,
     &  ' surfaces,',i4,' segments,',i7/
     &	'  rays (',i5,' for sensitivity calcs). Up to',i2,
     &	' diffraction planes are stored.'/'  Up to ',i3,
     &	' interpolated surfaces of up to ',i6,
     &	' data points are supported.'/
     &	'  Maximum pixel array size is',i6,' by',i6,'.'/
     &  '  Diffraction grid size is',i6,' by',i6,'.'/)
 312	FORMAT(
     & '  Glass table ',A,' read, with ',I3,' glasses included.'//
     & '  For a list of commands, type HELp.'/
     & '  Use RX command to show prescription files in current folder.'/
     & '  Use up/down arrow keys to access command history.')
 313	FORMAT(' Environment variable MACOS_HOME not set.'/)
 314	FORMAT('  No glass tables were found.'//
     & '  For a list of commands, type HELp.'/
     & '  Use RX command to show prescription files in current folder.'/
     & '  Use up/down arrow keys to access command history.')
 300	FORMAT(' MACOS command summary:'/
     &	'  Commands:   Quit or END MACOS'/
     &	'              NEW: build new optical system dataset'/
     &	'              OLD: load existing optical system dataset'/
     &  '              EXEcute: run a .jou journal file'/
     &  '              JOUrnal: generate a .jou journal file'/
     &	'              SUMmarize optical system data'/
     &	'              STAtus of current calculations'/
     &  '              RESet options and defaults to starting values'/
     &	'              SHOw data for a particular element'/
     &	'              MODify optical system data'/
     &	'              SAVe optical system data'/
     &	'              EXPort specified results in various formats')
 307    FORMAT('  Commands for ray-trace analysis are:'/
     &	'           MAP prints ray and segment ID maps'/
     &	'           RAYtrace a specified ray'/
     &	'           LNEG allows negative pathlengths (default)'/
     &	'           NOLNEG turns off negative pathlengths'/
     &	'           POLarized light'/
     &	'           NOPOL: no polarization calculations'/
     &	'           OBS: spot diagram obscuration option'/
     &	'           SINT: setup interpolated surface data'/
     &	'           UDSINIT: setup user-defined surface data'/
     &	'           PERturb optical system'/
     &	'           PREad: perturb using input from a file'/
     &	'           OPD plot and analysis'/
     &	'           DOPd plot differential OPD maps'/
     &	'           CENTEr beam in system stop'/
     &	'           STOP defines system stop'/
     &	'           COOrd: compute beam coordinates at an element'/
     &	'           CENTRoid sets centroid as reference pt for FEX'/
     &	'           FFP: find field point to place chief ray at a'/
     &	'           specified point on an element'/
     &	'           PFP: find field point to place chief ray at a'/
     &	'           specified pixel location on a detector'/
     &	'           SPCenter: specify spot diagram centering option'/
     &	'           SPOt diagram plot'/
     &  '           GBS: report beam sizes on selected elements'/
     &  '           FDP finds detector plane position and orientation')
 304    FORMAT('  Commands for linear model building and analysis are:'/
     &	'           BUild linear model using ray partials'/
     &	'           DMBUild including a deformable mirror'/
     &	'           PARtials: print partials of a specified ray'/
     &	'           LPErturb optical system using linear model'/
     &	'           LPRead: perturb using input from a file'/
     &	'           LREset: zero previous linear perturbations'/
     &	'           LSPot: perturbed linear system spot diagram'/
     &	'           LOPD: perturbed linear system OPD map'/
     &	'           LPIxilate: plot perturbed pixilated image'/
     &	'           LINtensity: plot perturbed wavefront')
 305	FORMAT('  Commands for diffraction analysis are:'/
     &	'           PROpagate diffraction wavefront'/
     &	'           BEAm: set beam type'/
     &	'           VECtor diffraction (with polarized light)'/
     &	'           SCAlar diffraction (with polarized light)'/
     &	'           REGrid at specified element during propagation'/
     &	'           NOREGrid: regridding off at specified element'/
     &	'           ORS to optimize reference surface'/
     &	'           SRS to slave one reference surface to another'/
     &	'           LNEGok allows negative ray lengths at ref surfs'/
     &	'           NOLNeg prohibits negative ray lengths'/
     &	'           COMpose a multiple object/color image'/
     &	'  Commands to plot wavefronts and images are:'/
     &	'              INtensity of the wavefront output'/
     &	'              PIXillated intensity output'/
     &	'              ADD current intensity to composed image'/
     &	'              DADD: display current composed image'/
     &	'              LOg intensity output')
 306	FORMAT('              AMplitude and phase output'/
     &	'              REal and imaginary wavefront output'/
     &	'              PFP: find field point to place chief ray at a'/
     &	'              location specified in pixel coordinates'/
     &	'              PLOcate: specify location of pixel array'/
     &	'              SEEd sets image noise seed value'/
     &	'              NOIse adds noise to a COMposed image'/
     &	'              STRetch specifies display data stretch'/
     &	'  Commands to set wavefront plot type are:'/
     &	'              WIRe sets wireframe surface plotting'/
     &	'              SLIce sets "slice" surface plotting'/
     &	'              GRay sets gray-scale surface plotting'/
     &	'              COLumn sets column line-plotting'/
     &	'              CONtour sets contour surface plotting'/
     &	'              TEXt print output'/
     &	'              BINary image file output'/
     &	'              FITS image file output'/
     &	'              MATlab image file output')
 321     FORMAT('  Commands for system optimizations are:'/
     &  '              SETOPT set options for optimization'/
     &  '              AVAR add a variable element'/
     &  '              MVAR modify a variable element'/
     &  '              DVAR delete a variable element'/
     &  '              VARS list all variable elements'/
     &  '              AFOV add a field of view'/
     &  '              DFOV delete a field of view'/
     &  '              FOVS list all field of view'/
     &  '              CALib perform system optimization'/
     &	'  Help repeats this message')
 308	FORMAT(1P,' Field angle finder results:'/
     &	' Did',i5,' iterations, error=',D17.9/
     &	'  Old dx=',D17.9,'  New dx=',D17.9/9x,D17.9,9x,D17.9/
     &	' Old crd=',D17.9,' New crd=',D17.9,2(/9x,D17.9,9x,D17.9)/
     &	' Old crp=',D17.9,' New crp=',D17.9,2(/9x,D17.9,9x,D17.9))
 309	FORMAT(1P,' Beam centered in',i4,' iterations, error=',D17.9)
 310	FORMAT(1P,'  Computed StopPos =',3D20.12)
 315    FORMAT(1P,'  Stop specified at Elt',I3,', offset =',2D20.12)
 311	FORMAT(1P,' Computed xLocal=',3D24.16/
     &	' Computed yLocal=',3D24.16/' Computed zLocal=',3D24.16)
 316    FORMAT(1P,' xLocal=',3D24.16/
     &            ' yLocal=',3D24.16/' zLocal=',3D24.16)
 317    FORMAT(1P,A11,3D22.12)
 318    FORMAT(1P,A12,3D22.12)
 319    FORMAT(1P,A7,I3,A2,A14,A2,A12,A2,A18,D18.10)
 320    FORMAT(5P,I3,I8,'  ',D17.9,'  ',D17.9,'  ',D17.9)
 322    FORMAT(1P,255D17.9)
 323    FORMAT(1P,511D17.9)
 324    FORMAT(1P,1023D17.9)
 325    FORMAT(1P,127D17.9)
 402    FORMAT(1P,' RMS OPD error is ',e13.6)

 442    FORMAT(1P,' RMS OPD error is ',D17.9/
     &         1P,' P-V OPD error is ',D17.9/
     &	       1P,' Average OPD is   ',D24.16)

 443    FORMAT(1P,' RMS OPD error is ',D17.9,A3/
     &         1P,' P-V OPD error is ',D17.9,A3/
     &         1P,' Average OPD is   ',D24.16,A3)

 403	FORMAT(' Enter 5 output coordinate 7-vectors:')
 407	FORMAT(1P,
     &	14x,'Chief ray location: x=',D14.7,' y=',D14.7,' z=',D14.7/
     &  15x,'Centroid location: x=',D14.7,' y=',D14.7,' z=',D14.7/
     &  '  Centroid offset from chief ray: x=',D14.7,' y=',D14.7,
     &	' z=',D14.7)
 408	FORMAT(1P,'        Chief ray location: x=',D14.7,' y=',D14.7/
     &	'  Centroid of spot diagram: x=',D14.7,' y=',D14.7)
 410	FORMAT(1P,' Average OPD is ',e13.6)
 411	FORMAT(1P,' Average total path is ',e20.13,/,
     &	' Average delta path is ',e20.13)
 412    FORMAT('   Computed system effective focal length = ',D12.5,
     &  ' mm')
 816    FORMAT(' Compute time was',f10.4,' sec')


