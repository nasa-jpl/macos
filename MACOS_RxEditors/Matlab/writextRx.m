%                     FOR OFFICIAL USE ONLY
%***********************************************************************
%    +----------------------------------------------------------------+
%    |  Copyright (C) 2009, California Institute of Technology.       |  
%    |  U.S. Government Sponsorship is acknowledged.		          |
%    +----------------------------------------------------------------+
%***********************************************************************
%
% Purpose: This particular file loads pre-created structural array labeled
%          most2prc_allseg.mat which contains the Source and Segment
%          parameter values from most2prc_allseg.in. It then writes the
%          Source and the segments indicated by segparam to a text file.
%
% Record:
%       Date                  Author                    Description
% ================      ========================      ====================
% May 2nd, 2008           Luis F. Marchen               Original Code
%
% Define Variables and Constants:
% filename    ---mat file containing the Rx elements with modified
%                parameters
% segparam    ---contains number array with number of elements to be
%                written to text file, this indicates which segment must be
%                written to the file

function writeRxfile=writextRx(filename,segparam);
load ([filename]); %loads the structural array wich source
%==========================================================================
% Code to write segments to the text file above 
%==========================================================================
%--------------------------------------------------------------------------
% Writing source to file specified by filename
%--------------------------------------------------------------------------
fid = fopen('most2prc_allseg.txt', 'wt'); %open text file to write

fprintf(fid, 'BaseUnits=  %s\n',most2prc_allseg.Source.BaseUnits);
fprintf(fid, 'WaveUnits=  %s\n',most2prc_allseg.Source.WaveUnits);

fprintf(fid, 'ChfRayDir=  %12.14f\t' ,most2prc_allseg.Source.ChfRayDir(1,1));
fprintf(fid, '%12.14f\t' ,most2prc_allseg.Source.ChfRayDir(1,2:3));

fprintf(fid, '\nChfRayPos=  %12.14f\t' ,most2prc_allseg.Source.ChfRayPos(1,1));
fprintf(fid, '%12.14f\t' ,most2prc_allseg.Source.ChfRayPos(1,2:3));

fprintf(fid, '\n  zSource=  %s\n',strrep(num2str(most2prc_allseg.Source.zSource,'%10.14e'),'e+0','D+'));
fprintf(fid, '   IndRef=  %12.14f\n',most2prc_allseg.Source.IndRef);
fprintf(fid, '   Extinc=  %12.14f\n',most2prc_allseg.Source.Extinc);
fprintf(fid, '  Wavelen=  %s\n', strrep(num2str(most2prc_allseg.Source.Wavelen,'%10.14e'),'e-0','D-'));
fprintf(fid, '     Flux=  %12.14f\n',most2prc_allseg.Source.Flux);
fprintf(fid, ' GridType=  %s\n',most2prc_allseg.Source.GridType);
fprintf(fid, ' Aperture=  %s\n',strrep(num2str(most2prc_allseg.Source.Aperture,'%10.14e'),'e+0','D+'));
fprintf(fid, ' Obscratn=  %12.14f\n',most2prc_allseg.Source.Obscratn);
fprintf(fid, ' nGridpts=  %1.0f\n',most2prc_allseg.Source.nGridpts);

fprintf(fid, '    xGrid=  %12.14f\t' ,most2prc_allseg.Source.xGrid(1,1));
fprintf(fid, '%12.14f\t' ,most2prc_allseg.Source.xGrid(1,2:3));

fprintf(fid, '\n    yGrid=  %12.14f\t' ,most2prc_allseg.Source.yGrid(1,1));
fprintf(fid, '%12.14f\t' ,most2prc_allseg.Source.yGrid(1,2:3));

fprintf(fid, '\n SegXGrid=  %12.14f\t' ,most2prc_allseg.Source.SegXGrid(1,1));
fprintf(fid, '%12.14f\t' ,most2prc_allseg.Source.SegXGrid(1,2:3));
fprintf(fid, '\n     nElt=  %1.0f\n',most2prc_allseg.Source.nElt);

%--------------------------------------------------------------------------
% Writing segments to most2prc_allseg.txt
%--------------------------------------------------------------------------
for j=segparam;
fprintf(fid, '\n     iElt=  %1.0f',most2prc_allseg.Segment(j).iElt); 
fprintf(fid, '\n  EltName=  %s',most2prc_allseg.Segment(j).EltName); 
fprintf(fid, '\n  Element=  %s',most2prc_allseg.Segment(j).Element); 
fprintf(fid, '\n  NSCount=  %1.0f',most2prc_allseg.Segment(j).NSCount); 
fprintf(fid, '\n  Surface=  %s',most2prc_allseg.Segment(j).Surface); 
fprintf(fid, '\n    KrElt= %s',strrep(num2str(most2prc_allseg.Segment(j).KrElt,'%10.14e'),'e+0','D+'));
fprintf(fid, '\n    KcElt= %s\n',strrep(num2str(most2prc_allseg.Segment(j).KcElt,'%10.14e'),'e+0','D+'));

fprintf(fid, '   psiElt=');
for ii=1:3;
    fprintf(fid, '  %10.14f\t',most2prc_allseg.Segment(j).psiElt(1,ii));
end

fprintf(fid, '\n   VptElt=');
for ii=1:3;
    fprintf(fid, '  %s\t',strrep(num2str(most2prc_allseg.Segment(j).VptElt(1,ii),'%10.14e'),'e+0','D+'));
end
 
fprintf(fid, '\n   RptElt=');
for ii=1:3;
    fprintf(fid, '  %s\t',strrep(num2str(most2prc_allseg.Segment(j).RptElt(1,ii),'%10.14e'),'e+0','D+'));
end

fprintf(fid, '\n   IndRef=  %s',strrep(num2str(most2prc_allseg.Segment(j).IndRef,'%10.14e'),'e+0','D+'));
fprintf(fid, '\n   Extinc=  %s',strrep(num2str(most2prc_allseg.Segment(j).Extinc,'%1.0e'),'e+0','D+'));
fprintf(fid, '\n ZernType=  %s',most2prc_allseg.Segment(j).ZernType); 

fprintf(fid, '\n ZernCoef=');
for ii=1:6;
    fprintf(fid, '  %1.1f',most2prc_allseg.Segment(j).ZernCoef(1,ii));
end
fprintf(fid, '\n          ');
for jj=2:7;
    for ii=1:6;
        fprintf(fid, '  %1.1f',most2prc_allseg.Segment(j).ZernCoef(jj,ii));
    end
    fprintf(fid, '\n          ');
end;
for ii=1:3;
    fprintf(fid, '  %1.1f',most2prc_allseg.Segment(j).ZernCoef(8,ii));
end

fprintf(fid, '\n GridFile=  %s',most2prc_allseg.Segment(j).GridFile); 
fprintf(fid, '\n nGridMat=  %1.0f',most2prc_allseg.Segment(j).nGridMat); 
fprintf(fid,'\n GridSrfdx=  %6.14f',most2prc_allseg.Segment(j).GridSrfdx); 

fprintf(fid, '\n     pMon=');
for ii=1:3;
    fprintf(fid, '  %6.14f',most2prc_allseg.Segment(j).pMon(1,ii));
end

fprintf(fid, '\n     xMon=');
for ii=1:3;
    fprintf(fid, '  %6.14f',most2prc_allseg.Segment(j).xMon(1,ii));
end

fprintf(fid, '\n     yMon=');
for ii=1:3;
    fprintf(fid, '  %6.14f',most2prc_allseg.Segment(j).yMon(1,ii));
end

fprintf(fid, '\n     zMon=');
for ii=1:3;
    fprintf(fid, '  %6.14f',most2prc_allseg.Segment(j).zMon(1,ii));
end

fprintf(fid, '\n     lMon=  %s',strrep(num2str(most2prc_allseg.Segment(j).lMon,'%10.14e'),'e+0','D+'));
fprintf(fid, '\n     nObs=  %1.0f',most2prc_allseg.Segment(j).nObs); 

fprintf(fid, '\n     xObs=');
for ii=1:3;
    fprintf(fid, '  %s',strrep(num2str(most2prc_allseg.Segment(j).xObs(1,ii),'%1.0e'),'e+00','d'));
end

fprintf(fid, '\nPolyApVec=  %1.0f',most2prc_allseg.Segment(j).PolyApVec(1,1));
fprintf(fid, '\n      ');
for jj=2:6;
    for ii=1:3;
        fprintf(fid, '     %6.14f',most2prc_allseg.Segment(j).PolyApVec(jj,ii));
    end
    fprintf(fid, '\n      ');
end;
for ii=1:3;
    fprintf(fid, '     %6.14f',most2prc_allseg.Segment(j).PolyApVec(7,ii));
end

fprintf(fid,'\n     zElt=  %6.14f',most2prc_allseg.Segment(j).zElt); 
fprintf(fid,'\n PropType=  %s',most2prc_allseg.Segment(j).PropType); 
fprintf(fid,'\n  nECoord=  %1.0f',most2prc_allseg.Segment(j).nECoord); 

fprintf(fid, '\n    TElt=');
for jj=1:6;
    for ii=1:6;
        fprintf(fid, '  %10.14f',most2prc_allseg.Segment(j).TElt(jj,ii));
    end
    fprintf(fid, '\n         ');
end;
fprintf(fid, '\n');
end;
fclose(fid) %close text file and finish writing

return





