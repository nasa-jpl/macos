%                       FOR OFFICIAL USE ONLY
%
%**************************************************************************
%   + -----------------------------------------------------------------+
%   | Copyright [2009], by the California Institute of Technology.     |
%   | ALL RIGHTS RESERVED. United States Government Sponsorship        | 
%   | acknowledged.                                                    |
%   |                                                                  |
%   | Any commercial use must be negotiated with the Office of         |
%   | Technology Transfer at the California Institute of Technology.   |
%   | This software may be subject to U.S. export control laws. By     | 
%   | accepting this software, the user agrees to comply with all      |
%   | applicable U.S. export laws and regulations. User has the        |
%   | responsibility to obtain export licenses, or other export        |
%   | authority as may be required before exporting such information   |
%   | to foreign countries or providing access to foreign persons.     |
%   *------------------------------------------------------------------+
%**************************************************************************
%
% Purpose: This file reads MACOS optical prescriptions, and creates 
%          structural arrays which contain all elements in the prescription
%          with parameters and their values. The file is able to read files
%          with comments and empty lines; however, comments at the end of
%          none empty lines should be avoided.
%
% Record:
%       Date                  Author                    Description
% ================      ========================      ====================
%  May 28th, 2008           Luis F. Marchen               Original Code
%
% Define Variables and Constants:
% fid         ---file identifier for eading lines              
% count       ---number of iteration of while loop not counting empty lines
%                or comments
function []=WriteRxFile(write_fname,ElementData,credits)
%==========================================================================
% Comments or credits to be written to top of write file, and parameters 
%==========================================================================
overwrite_wfname=1;
note=credits;
%==========================================================================
% The code below specifies the MACOS file to be read, reads the file, and
% creates a structural array which contains the source, all optical
% elements containing iElt param, and the Tout and nOutCord parameters.
%==========================================================================
OpticalData=ElementData;
num_Elt=length(OpticalData);
%==========================================================================
% This file is currently labeled "Read ..." but we will use part of it
% below to actually define some parameters to print file. We also need to
% write the code which will allow us to make updates to the parameters.
%==========================================================================
%--------------------------------------------------------------------------
% If write file exits it will be delited from directory 
%--------------------------------------------------------------------------
if exist(write_fname) == 2;	%fname already exists as file
    if overwrite_wfname;
        warning([write_fname ' already exists, and it is being deleted'])
        disp(['Deleting ' write_fname])
        delete(write_fname);
    else
        warning([write_fname ' already exists!!!'])
        return
    end
end

fname_id = fopen(write_fname,'w');
%==========================================================================
% Write file name and credits (comments) to top of write file
%==========================================================================
fprintf(fname_id,'%% %s\n',write_fname);
if exist('credits') == 1
  % allow multiline notes
  [m,n]= size(note);
  for i=1:m
    fprintf(fname_id,'%% %s\n', credits(i,:));    
  end
end
%date file
fprintf(fname_id,'%% Date Created: %s\n\n',datestr(date, 2));
%==========================================================================
% Write parameters for source and optical elements
%==========================================================================
for elt=1:num_Elt;
    field_names = fieldnames(OpticalData{elt}); %Field names of structure, or public fields of object
    param_names = strvcat(field_names);      %parameter labels (names)
    [m n]=size(param_names);
    disp(['Element Being Written is:' num2str(elt)])
    for pnum=1:m;       %loops over all parameters of current element
        param_label = param_names(pnum,:);
        param_label=strtrim(param_label);

%--------------------------------------------------------------------------
% Code below writes source parameters to write file
%--------------------------------------------------------------------------
        switch param_label;
            case 'BaseUnits'
                temp_param=OpticalData{elt}.BaseUnits;
                fprintf(fname_id, '\nBaseUnits=  %s',temp_param);
            case 'WaveUnits'
                temp_param=OpticalData{elt}.WaveUnits;
                fprintf(fname_id, '\nWaveUnits=  %s',temp_param);
            case 'ChfRayDir'
                temp_param=OpticalData{elt}.ChfRayDir;
                [row col]=size(temp_param);
                fprintf(fname_id, 'ChfRayDir=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'ChfRayPos'
                temp_param=OpticalData{elt}.ChfRayPos;
                [row col]=size(temp_param);
                fprintf(fname_id, '\nChfRayPos=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'zSource'
                temp_param=OpticalData{elt}.zSource;
                fprintf(fname_id, '\n  zSource=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'IndRef'
                temp_param=eval(['OpticalData{elt}.' param_names(pnum,:) ';']);
                fprintf(fname_id, '\n   IndRef=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'Extinc'
                temp_param=OpticalData{elt}.Extinc;
                fprintf(fname_id, '\n   Extinc=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'Wavelen'
                temp_param=OpticalData{elt}.Wavelen;
                fprintf(fname_id, '\n  Wavelen=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'Flux'
                temp_param=OpticalData{elt}.Flux;
                fprintf(fname_id, '\n     Flux=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'GridType'
                temp_param=OpticalData{elt}.GridType;
                fprintf(fname_id, '\n GridType=  %s',temp_param);
            case 'Aperture'
                temp_param=OpticalData{elt}.Aperture;
                fprintf(fname_id, '\n Aperture=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'Obscratn'
                temp_param=OpticalData{elt}.Obscratn;
                fprintf(fname_id, '\n Obscratn=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'nGridpts'
                temp_param=OpticalData{elt}.nGridpts;
                fprintf(fname_id, '\n nGridpts=');
                fprintf(fname_id, '  %1.0f',temp_param);
            case 'xGrid'
                temp_param=OpticalData{elt}.xGrid;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n    xGrid=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'yGrid'
                temp_param=OpticalData{elt}.yGrid;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n    yGrid=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'nElt'
                temp_param=eval(['OpticalData{elt}.' param_names(pnum,:) ';']);
                fprintf(fname_id, '\n     nElt=');
                fprintf(fname_id, '  %1.0f',temp_param);
            case 'nSeg'
                temp_param=OpticalData{elt}.nSeg;
                fprintf(fname_id, '\n     nSeg=  %1.0f',temp_param);
            case 'width'
                temp_param=OpticalData{elt}.width;
                fprintf(fname_id, '\n    width=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.3e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.3e'),'e+0','D+'));
                end;
            case 'gap'
                temp_param=OpticalData{elt}.gap;
                fprintf(fname_id, '\n     gap=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;                
            case 'SegXgrid'
                temp_param=OpticalData{elt}.SegXgrid;
                [row col]=size(temp_param);
                fprintf(fname_id, '\nSegXgrid=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'SegCoord'
                temp_param=OpticalData{elt}.SegCoord;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n SegCoord=');
                for ii=1:col;
                    if temp_param(1,ii)<1;
                    fprintf(fname_id, '      %.0f',temp_param(1,ii));
                    elseif temp_param(1,ii)==0
                    fprintf(fname_id,'         %.0f',temp_param(1,ii));
                    else
                    fprintf(fname_id,'       %.0f',temp_param(1,ii));   
                    end
                end
                fprintf(fname_id, '\n          ');
                for jj=2:row;
                    for ii=1:col;
                        if temp_param(jj,ii)<1;
                        fprintf(fname_id, '      %.0f',temp_param(jj,ii));
                        elseif temp_param(jj,ii)==0
                        fprintf(fname_id,'         %.0f',temp_param(jj,ii));
                        else
                        fprintf(fname_id,'       %.0f',temp_param(jj,ii));  
                        end
                    end
                    fprintf(fname_id,'\n          ');
                end;             
%--------------------------------------------------------------------------
% Code below writes optical element parameters to write file                
%--------------------------------------------------------------------------
            case 'iElt'
                temp_param=OpticalData{elt}.iElt;
                fprintf(fname_id, '\n\n     iElt=  %1.0f',temp_param);
            case 'EltName'
                temp_param=OpticalData{elt}.EltName;
                fprintf(fname_id, '\n  EltName=  %s',temp_param);
            case 'Element'
                temp_param=OpticalData{elt}.Element;
                fprintf(fname_id, '\n  Element=  %s',temp_param);
            case 'NSCount'
                temp_param=OpticalData{elt}.NSCount;
                fprintf(fname_id, '\n  NSCount=  %1.0f',temp_param);
            case 'Surface'
                temp_param=OpticalData{elt}.Surface;
                fprintf(fname_id, '\n  Surface=  %s',temp_param);
            case 'UDSrfType'
                temp_param=OpticalData{elt}.UDSrfType;
                fprintf(fname_id, '\n UDSrfType=  %1.0f',temp_param);
            case 'UDSrfParam'
                temp_param=OpticalData{elt}.UDSrfParam;
                [row col]=size(temp_param);
                fprintf(fname_id, '\nUDSrfParam=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'KrElt'
                temp_param=OpticalData{elt}.KrElt;
                fprintf(fname_id, '\n    KrElt=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'KcElt'
                temp_param=OpticalData{elt}.KcElt;
                fprintf(fname_id, '\n    KcElt=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'AsphCoef'
                temp_param=OpticalData{elt}.AsphCoef;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   AsphCoef=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'psiElt'
                temp_param=OpticalData{elt}.psiElt;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   psiElt=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'VptElt'
                temp_param=OpticalData{elt}.VptElt;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   VptElt=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'RptElt'
                temp_param=OpticalData{elt}.RptElt;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   RptElt=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'ZernType'
                temp_param=OpticalData{elt}.ZernType;
                fprintf(fname_id, '\n ZernType=  %s',temp_param);
            case 'ZernCoef'
                temp_param=OpticalData{elt}.ZernCoef;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n ZernCoef=');
                for ii=1:col;
                    fprintf(fname_id, '  %1.1f',temp_param(1,ii));
                end
                fprintf(fname_id, '\n          ');
                for jj=2:row-1;
                    for ii=1:6;
                        fprintf(fname_id, '  %1.1f',temp_param(jj,ii));
                    end
                    fprintf(fname_id, '\n          ');
                end;
                for ii=1:3;
                    fprintf(fname_id, '  %1.1f',temp_param(8,ii));
                end
            case 'GridFile'
                temp_param=OpticalData{elt}.GridFile;
                fprintf(fname_id, '\n GridFile=  %s',temp_param);              
            case 'nGridMat'
                temp_param=OpticalData{elt}.nGridMat;
                fprintf(fname_id, '\n nGridMat=  %1.0f',temp_param);
            case 'GridSrfdx'
                temp_param=OpticalData{elt}.GridSrfdx;
                fprintf(fname_id, '\nGridSrfdx=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'pMon'
                temp_param=OpticalData{elt}.pMon;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n     pMon=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'xMon'
                temp_param=OpticalData{elt}.xMon;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n     xMon=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'yMon'
                temp_param=OpticalData{elt}.yMon;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n     yMon=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'zMon'
                temp_param=OpticalData{elt}.zMon;
               [row col]=size(temp_param);
                fprintf(fname_id, '\n     zMon=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'lMon'
                temp_param=OpticalData{elt}.lMon;
                fprintf(fname_id, '\n     lMon=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'nObs'
                temp_param=OpticalData{elt}.nObs;
                fprintf(fname_id, '\n     nObs=  %1.0f',temp_param);
            case 'ObsType'
                temp_param=OpticalData{elt}.ObsType;
                fprintf(fname_id, '\n  ObsType=  %s',temp_param);
            case 'ObsType1'
                temp_param=OpticalData{elt}.ObsType1;
                fprintf(fname_id, '\n  ObsType=  %s',temp_param);
            case 'ObsType2'
                temp_param=OpticalData{elt}.ObsType2;
                fprintf(fname_id, '\n  ObsType=  %s',temp_param);  
            case 'ObsVec'
                temp_param=OpticalData{elt}.ObsVec;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   ObsVec=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'ObsVec1'
                temp_param=OpticalData{elt}.ObsVec1;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   ObsVec=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
            case 'ObsVec2'
                temp_param=OpticalData{elt}.ObsVec2;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n   ObsVec=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end   
            case 'xObs'
                temp_param=OpticalData{elt}.xObs;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n     xObs=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end       
            case 'ApType'
                temp_param=OpticalData{elt}.ApType;
                fprintf(fname_id, '\n   ApType=  %s',temp_param);
            case 'ApVec'
                temp_param=OpticalData{elt}.ApVec;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n    ApVec=');
                %for ii=1:col;
                %    fprintf(fname_id, '  %4.1f',temp_param(1,ii));
                %end
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end 
            case 'PolyApVec'
                temp_param=OpticalData{elt}.PolyApVec;
                [row col]=size(temp_param);
                fprintf(fname_id, '\nPolyApVec=');
                for ii=1;
                    fprintf(fname_id, '  %1.0f',temp_param(1,ii));
                end          
                fprintf(fname_id, '\n          ');
                for jj=2:row;
                    for ii=1:col;
                        if abs(temp_param(jj,ii))<1 && abs(temp_param(jj,ii)~=0);
                            fprintf(fname_id, '  %s',strrep(num2str(temp_param(jj,ii),'%1.15e'),'e-0','D-'));
                        else
                            fprintf(fname_id, '  %s',strrep(num2str(temp_param(jj,ii),'%1.15e'),'e+0','D+'));
                        end;
                    end
                    if jj<row
                        fprintf(fname_id, '\n          ');
                    end
                end;
            case 'zElt'
                temp_param=OpticalData{elt}.zElt;
                fprintf(fname_id, '\n     zElt=');
                if abs(temp_param)<1 && abs(temp_param~=0);
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.15e'),'e+0','D+'));
                end;
            case 'PropType'
                temp_param=OpticalData{elt}.PropType;
                fprintf(fname_id, '\n PropType=  %s',temp_param);
            case 'nECoord'
                temp_param=OpticalData{elt}.nECoord;
                fprintf(fname_id, '\n  nECoord=  %1.0f',temp_param);
                
            case 'TElt'
                temp_param=OpticalData{elt}.TElt;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n     TElt=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii)~=0);
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    elseif abs(temp_param(1,ii))>1 && abs(temp_param(1,ii)~=0);
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end
                fprintf(fname_id, '\n          ');
                for jj=2:row;
                    for ii=1:col;
                        if abs(temp_param(jj,ii))<1 && abs(temp_param(jj,ii)~=0);
                            fprintf(fname_id, '  %s',strrep(num2str(temp_param(jj,ii),'%1.15e'),'e-0','D-'));
                        else
                            fprintf(fname_id, '  %s',strrep(num2str(temp_param(jj,ii),'%1.15e'),'e+0','D+'));
                        end;
                    end
                    if jj<row
                        fprintf(fname_id, '\n          ');
                    end
                end;
%                 case 'TElt'
%                 temp_param=OpticalData{elt}.TElt;
%                 [row col]=size(temp_param);
%                 fprintf(fname_id, '\n     TElt=');
%                 for ii=1:col;
%                     if temp_param(1,ii)~=0;
%                         fprintf(fname_id, '  %2.15f',temp_param(1,ii));
%                     else
%                         fprintf(fname_id, '  %2.1f',temp_param(1,ii));
%                     end;
%                 end
%                 fprintf(fname_id, '\n          ');
%                 for jj=2:row;
%                     for ii=1:col;
%                         if temp_param(jj,ii)~=0;
%                             fprintf(fname_id, '  %2.15f',temp_param(jj,ii));
%                         else
%                             fprintf(fname_id, '  %2.1f',temp_param(jj,ii));
%                         end;
%                     end
%                     if jj<row
%                         fprintf(fname_id, '\n          ');
%                     end
%                 end;
%--------------------------------------------------------------------------
% The code below writes output coordinate parameters to write file
%--------------------------------------------------------------------------
            case 'nOutCord'
                temp_param=OpticalData{elt}.nOutCord;
                fprintf(fname_id, '\n\n nOutCord=  %1.0f',temp_param);
            case 'Tout'
                temp_param=OpticalData{elt}.Tout;
                [row col]=size(temp_param);
                fprintf(fname_id, '\n     Tout=');
                for ii=1:col;
                    if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii)~=0);
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e-0','D-'));
                    else
                        fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.15e'),'e+0','D+'));
                    end;
                end          
                fprintf(fname_id, '\n          ');
                for jj=2:row;
                    for ii=1:7;
                        if abs(temp_param(jj,ii))<1 && abs(temp_param(jj,ii)~=0);
                            fprintf(fname_id, '  %s',strrep(num2str(temp_param(jj,ii),'%1.15e'),'e-0','D-'));
                        else
                            fprintf(fname_id, '  %s',strrep(num2str(temp_param(jj,ii),'%1.15e'),'e+0','D+'));
                        end;
                    end
                    fprintf(fname_id, '\n          ');
                end;
              
        end
    end
end;
fclose(fname_id);
return

