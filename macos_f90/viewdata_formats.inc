
 100    FORMAT(1P,'ChfRayDir=',3D23.15/'ChfRayPos=',3D23.15/
     &  ' zSource= ',D23.15/'  IndRef= ',D23.15/'  Extinc= ',D23.15/
     &  'BaseUnits= ',A16/'WaveUnits= ',A16/' Wavelen= ',D23.15/
     &  '    Flux= ',D23.15/'GridType=   ',A16/'Aperture= ',D23.15/
     &  'Obscratn= ',D23.15/'nGridpts= ',i4/'   xGrid= ',3D23.15/
     &  '   yGrid= ',3D23.15/'    nElt= ',i4)
 116    FORMAT(1P,10x,i5,5x,i5,5x,i5)
 117    FORMAT(1P,'    nSeg= ',i4/'   width= ',D23.15/
     &  '     gap= ',D23.15/'SegXgrid= ',D23.15,2(2x,D23.15)/
     &  'SegCoord= ',i5,2(5x,i5))
 118    FORMAT(1P,'    nSeg= ',i4/
     &  '  radCtr= ',D23.15/'SegXgrid= ',D23.15,2(2x,D23.15)/
     &  'SegCoord= ',i5,5x,i5)

C Basic params & optical params & conic params
  98    FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8/
     &  '  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
  99    FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8/
     &  '   KrElt= ',D23.15/'   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
 101    FORMAT(1P/'    iElt= ',i4/' EltName= ',a16/
     &  ' Element= ',a12/ ' Surface= ',a8/
     &  '    fElt= ',D23.15/'    eElt= ',D23.15/
     &  '   KrElt= ',D23.15/'   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
  95    FORMAT(1P,'  IndRef= ',D23.15/'  Extinc= ',D23.15)
  96    FORMAT(1P,'GlassElt= ',A16/'  Extinc= ',D23.15)
  97    FORMAT(1P,'GlassElt= ',A16/'GlassCoef= ',6D23.15/
     &  '  Extinc= ',D23.15)
 102    FORMAT(1P,'    nObs= ',i4)
 121    FORMAT(1P,'  ApType=   ',A16)

C  asphere params
 103    FORMAT(1P,'AsphCoef= ',4D23.15)

C  anamorph params
 104    FORMAT(1P,' AnaCoef= ',4D23.15)

C  monomial params
 105    FORMAT(1P,' MonCoef= ',6D23.15/  19(10x,6D23.15/))

C  user-defined params
 106    FORMAT(1P,'UDSrfType= ',i4/'UDSrfFile= ',A/
     &  'UDSrfParam= ',6D23.15)

C  grid surface params
 132    FORMAT(1P,'GridFile= ',A/'nGridMat= ',i4/
     &  'GridSrfdx= ',D23.15)

C  zernike params
 107    FORMAT(1P,'ZernCoef= ',6D23.15/ 6(10x,6D23.15/) (10x,3D23.15))

C  zernike type
 131    FORMAT(1P,'ZernType=   ',A16)

C  surface coord params
 108    FORMAT(1P,'    pMon= ',3D23.15/'    xMon= ',3D23.15/
     &  '    yMon= ',3D23.15/'    zMon= ',3D23.15)
 109    FORMAT(1P,'    lMon= ',D23.15)

C  HOE params & Grating params
 110    FORMAT(1P,'   h1HOE= ',3D23.15/'   h2HOE= ',3D23.15/
     &  'OrderHOE= ',D23.15/' WaveHOE= ',D23.15)
 111    FORMAT(1P,'   h1HOE= ',3D23.15/
     &  'OrderHOE= ',D23.15/'RuleWidth= ',D23.15)

C  lens array params
 112    FORMAT(1P,'LensArrayType= ',i2/'LensArrayWidth= ',D23.15)

C  obscuration params
 122    FORMAT(1P,'   ApVec= ',6D23.15)
 124    FORMAT(1P,' ObsType=   ',A16)
 125    FORMAT(1P,'  ObsVec= ',6D23.15)
 126    FORMAT(1P,'    xObs= ',3D23.15)
 127    FORMAT(1P,'   Poly3DApVec= ',I4)
 128    FORMAT(1P,'         ',3D23.15)
 129    FORMAT(1P,'   ApVec= ',2D23.15,I4)
 133    FORMAT(1P,'   EltGrp= ',$)
 136    FORMAT(1P,'         ',2D23.15)
 139    FORMAT(1P,'  ObsVec= ',2D23.15,I4)

C  coatings
 141    FORMAT(1P,'   nCoat= ',i4/'CoatIndx= ',5D23.15/3(10X,5D23.15))
 142    FORMAT(1P,'CoatThk= ',5D23.15/3(10X,5D23.15))

C  pinhole
 151    FORMAT(1P,' PinHole= ',D23.15)

C  basic params
 113    FORMAT(1P,'    zElt= ',D23.15/'PropType=   ',A16)
 114    FORMAT(1P,' nECoord= ',i4)
 115    FORMAT(1P,'    TElt= ',6D23.15)
 162    FORMAT(1P,10x,6D23.15)

C  Output Variables
 161    FORMAT(1P/'nOutCord= ',i4/'    Tout= ',7D23.15,6(/10x,7D23.15))
 171    FORMAT(1P,'Save3DApVec=  Y')

c sws
c some old stuff for showdata
c increment by 100
 230    FORMAT(1P/'    iElt= ',i4/' EltName= ',a32/
     &  ' Element= ',a12/ ' Surface= ',a8)
 233    FORMAT(1P,'    psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
 234    FORMAT(1P,'    fElt= ',D23.15/
     &  '    eElt= ',D23.15/
     &  '   KrElt= ',D23.15/'   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
 235    FORMAT(1P,'    KrElt= ',D23.15/
     &  '   KcElt= ',D23.15/'  psiElt= ',3D23.15/
     &  '  VptElt= ',3D23.15/'  RptElt= ',3D23.15)
