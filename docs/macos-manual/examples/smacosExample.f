C***********************************************************************
C       Begin file smacosExample.f

      IMPLICIT NONE

C   SMACOS variable declarations

      INCLUDE 'smacosvars.cmn'

C   These are local variables:

      INTEGER

      REAL*4

      REAL*8

C   Load prescription file

      command='OLD'
      CARG(1)='cassExitPupil'
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Set plot type

      command='GRAY'
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Set system stop

      command='STOP'
      CARG(1)='ELT'
      IARG(1)=3
      DARG(1)=0d0
      DARG(2)=0d0
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Set exit pupil

      command='FEX'
      IARG(1)=5
      CARG(1)='YES'
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Plot OPD at exit pupil after resetting ray grid size

      nGridPts=128
      command='RESET'
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

      command='OPD'
      IARG(1)=5
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Plot detector intensity after resetting ray grid size

      nGridPts=96
      command='RESET'

      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

      command='INTENSITY'
      IARG(1)=6
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Plot detector using a log10 stretch

      command='STRETCH'
      CARG(1)='LOG10'
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

      command='INTENSITY'
      IARG(1)=6
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

C   Quit!

      command='QUIT'
      CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
     & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,
     & PixArray,Cmatrix)

      STOP
      END

C       End file smacosExample.f
C***********************************************************************
