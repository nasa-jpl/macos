function[SegOptic]=mono2segs(seginfo)
Kc=seginfo.seg_Kc;
Kr=seginfo.seg_Kr;
diameter=seginfo.Hd;
TElt=seginfo.seg_TElt;
RptElt=seginfo.seg_ap_cen;
numseg=seginfo.numseg;
nGridMat=seginfo.nGridMat;

for ii=1:numseg;
%--------------------------------------------------------------------------
% make PolyApVec
%--------------------------------------------------------------------------
    PolyApVec=[6 0 0;seginfo.seg_data(:,:,ii)];
%--------------------------------------------------------------------------
% generate structural array
%--------------------------------------------------------------------------
    SegOptic{ii,1}.iElt=  ii;
    SegOptic{ii,1}.EltName=  seginfo.seg_name{ii};
    SegOptic{ii,1}.Element=  'NSReflector';
    SegOptic{ii,1}.NSCount=  1;
    SegOptic{ii,1}.Surface=  'Conic';
    SegOptic{ii,1}.KrElt=  -Kr;   %radius of curvature
    SegOptic{ii,1}.KcElt=  Kc;   %conic constant
    SegOptic{ii,1}.psiElt=  [0  0 -1];   %element orientation, usually segs keep mono PM psiElt
    SegOptic{ii,1}.VptElt=  [0  0  0];   %keeps mono PM vertex
    SegOptic{ii,1}.RptElt= RptElt(ii,:);    %this is the center of the segement
    SegOptic{ii,1}.IndRef=  1.00000000000D+00;
    SegOptic{ii,1}.Extinc=  1.00000000000D+22;
    SegOptic{ii,1}.ZernType=  'Malacara';
    SegOptic{ii,1}.ZernCoef=  zeros(8,6);
    SegOptic{ii,1}.GridFile=  'None';
    SegOptic{ii,1}.nGridMat=  nGridMat;
    SegOptic{ii,1}.GridSrfdx=  (diameter/(nGridMat-1));
    SegOptic{ii,1}.pMon=  RptElt(ii,:);
    SegOptic{ii,1}.xMon=  TElt(1:3,1,ii)';
    SegOptic{ii,1}.yMon=  TElt(1:3,2,ii)';
    SegOptic{ii,1}.zMon=  TElt(1:3,3,ii)';
    SegOptic{ii,1}.lMon=  diameter/2;
    SegOptic{ii,1}.nObs=  0;
    SegOptic{ii,1}.xObs= [1.00000000000e+00  0.00000000000e+00  0.00000000000e+00];
    SegOptic{ii,1}.ApType=  'Polygonal';
    SegOptic{ii,1}.PolyApVec=  PolyApVec;
    SegOptic{ii,1}.zElt=  Kr/2;
    SegOptic{ii,1}.PropType=  'Geometric';
    SegOptic{ii,1}.nECoord=  6;
    SegOptic{ii,1}.TElt=  TElt(:,:,ii);
end
return
            