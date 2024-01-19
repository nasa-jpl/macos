function [GridOptic]=create_grid_srf(Element_dat,xMon,ngrid,diameter,lMon)
param_names=fieldnames(Element_dat.Data);

if ~isempty(find(strcmp(param_names,'iElt'),1));
    iElt=Element_dat.Data.iElt;
end;
if ~isempty(find(strcmp(param_names,'EltName'),1));
    EltName=Element_dat.Data.EltName;
end;
if ~isempty(find(strcmp(param_names,'Element'),1));
    Element=Element_dat.Data.Element;
end;
if ~isempty(find(strcmp(param_names,'Surface'),1));
    Surface=Element_dat.Data.Surface;
end;
if ~isempty(find(strcmp(param_names,'KrElt'),1));
    KrElt=Element_dat.Data.KrElt;
end;
if ~isempty(find(strcmp(param_names,'KcElt'),1));
    KcElt=Element_dat.Data.KcElt;
end
if ~isempty(find(strcmp(param_names,'psiElt'),1));
    psiElt=Element_dat.Data.psiElt;
end
if ~isempty(find(strcmp(param_names,'VptElt'),1));
    VptElt=Element_dat.Data.VptElt;
end
if ~isempty(find(strcmp(param_names,'RptElt'),1));
    RptElt=Element_dat.Data.RptElt;
end
if ~isempty(find(strcmp(param_names,'IndRef'),1));
    IndRef=Element_dat.Data.IndRef;
end
if ~isempty(find(strcmp(param_names,'Extinc'),1));
    Extinc=Element_dat.Data.Extinc;
end
if ~isempty(find(strcmp(param_names,'ZernType'),1));
    ZernType=Element_dat.Data.ZernType;
else
    ZernType='Malacara';
end
if ~isempty(find(strcmp(param_names,'ZernCoef'),1));
    ZernCoef=Element_dat.Data.ZernCoef;
else
    ZernCoef=zeros(8,6);
end
if ~isempty(find(strcmp(param_names,'GridFile'),1));
    GridFile=Element_dat.Data.GridFile;
else
    GridFile='None';
end
if ~isempty(find(strcmp(param_names,'nObs'),1));
    nObs=Element_dat.Data.nObs;
end
if ~isempty(find(strcmp(param_names,'ApType'),1));
    ApType=Element_dat.Data.ApType;
end
if ~isempty(find(strcmp(param_names,'zElt'),1));
    zElt=Element_dat.Data.zElt;
end
if ~isempty(find(strcmp(param_names,'PropType'),1));
    PropType=Element_dat.Data.PropType;
end
if ~isempty(find(strcmp(param_names,'nECoord'),1));
    nECoord=Element_dat.Data.nECoord;
end
if ~isempty(find(strcmp(param_names,'TElt'),1));
    TElt=Element_dat.Data.TElt;
end

%--------------------------------------------------------------------------
% Here we define nGridMat, GridSrfdx, pMon, xMon, yMon, zMon, and pMon
%--------------------------------------------------------------------------
zMon=psiElt;
yMon=TElt(1:3,2)';%cross(zMon,xMon);
pMon=RptElt;
nGridMat=ngrid;
GridSrfdx=diameter/(ngrid-1);

%--------------------------------------------------------------------------
% Creating a structural array with new grid surface parameters
%--------------------------------------------------------------------------
GridOptic.iElt=iElt;
GridOptic.EltName=EltName;
GridOptic.Element=Element;
GridOptic.Surface=Surface;
GridOptic.KrElt=KrElt;
GridOptic.KcElt=KcElt;
GridOptic.psiElt=psiElt;
GridOptic.VptElt=VptElt;
GridOptic.RptElt=RptElt;
GridOptic.IndRef=IndRef;
GridOptic.Extinc=Extinc;
GridOptic.ZernType=ZernType;
GridOptic.ZernCoef=ZernCoef;
GridOptic.GridFile=GridFile;
GridOptic.nGridMat=nGridMat;
GridOptic.GridSrfdx=GridSrfdx;
GridOptic.pMon=pMon;
GridOptic.xMon=xMon;
GridOptic.yMon=yMon;
GridOptic.zMon=zMon;
GridOptic.lMon=lMon;
GridOptic.nObs=nObs;
GridOptic.ApType=ApType;
GridOptic.zElt=zElt;
GridOptic.PropType=PropType;
GridOptic.nECoord=nECoord;
if ~isempty(find(strcmp(param_names,'TElt'),1));
    GridOptic.TElt=TElt;
end

return




     

     
     