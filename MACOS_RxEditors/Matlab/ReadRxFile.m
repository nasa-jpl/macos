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
% Purpose: Reads MACOS optical prescriptions, and creates cell arrays of 
%          structures which contain all elements in the prescription
%          with parameters and their values. The function is able to read
%          files with comments and empty lines; however, comments at the 
%          end of none empty lines should be avoided.
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
% OpticalData ---cell array which contains parameters and values of optics
% param_names ---lists the names of all parameters in MACOS prescription
% macos_fname ---name of the MACOS prescription to be read
%
function [OpticalData,numElt,List_EltNames]= ReadRxFile(macos_fname) 
%==========================================================================
% Reading source parameters and values, and making the structural array
% section for source
%==========================================================================
%--------------------------------------------------------------------------
% Opening Prescription File to be Read
%--------------------------------------------------------------------------
fid=fopen(macos_fname);
%--------------------------------------------------------------------------
% Reading Parameter Names and Values, and Creating Structural
% Array Containing Optical Element Information Including Source and
% nOutCord (Tout and nOutCord are one element in structural array).
%--------------------------------------------------------------------------
count=0;
eltcount=1;
obs2_count=0;
obs3_count=0;

while ~feof(fid)
    dat= fgetl(fid);
    [str, remain] = strtok(dat);
    n=strfind(str,'='); %find index of "="sign to eliminate from str
    
    ss=isletter(remain);
    if isempty(ss);
        ss=0;
    end
    
    if ss(1)==1
    numeric_param=[];
    else
        numeric_param=str2num(remain);
    end
    
    string_param=str(1:n-1); %parameter name without "=" sign, otherwise empty if part
                             %of previous parameter
                             
%--------------------------------------------------------------------------
% The code below uses a continue statement to ignore empty lines or 
% comments. The code also determines when there
% is a new optical element, including the source, and the nOutCord, so that 
% the parameter of each elements can be grouped together in the structural
% array
%--------------------------------------------------------------------------
    if isempty(string_param)|| strncmp(string_param,'%',1);
        %disp('The Line Is Empty and No Structure Is Created');
        continue
    end;
    count=count+1;
%--------------------------------------------------------------------------
% Record how many obscurations an element has so that when reading an
% element to the individual structural array we keep track of the
% obscurations in the element
%--------------------------------------------------------------------------
if strcmp(str,'nObs=')==1;
    obs_count=numeric_param;
    obs2_count=0; %set to zero for a new element
    obs3_count=0; %set to zero for a new element
end
    
%--------------------------------------------------------------------------
% Determines and records when we begin reading a new optical element, 
% including the source, and the nOutCord, so that the parameter of each 
% elements can be grouped together in its corresponding structure
%--------------------------------------------------------------------------
    if strcmp(str,'iElt=')==1 || strcmp(str,'nOutCord=')==1 ;
        eltcount=eltcount+1;
    end
%--------------------------------------------------------------------------
% Code below extracts all values of ZernCoef from MACOS prescription
%--------------------------------------------------------------------------
    if strcmp(string_param,'ZernCoef');
        for jj=1:7;
            dat=fgetl(fid);
            if length(str2num(dat))==3;
                numeric_param=[numeric_param;0 0 0 str2num(dat)];
            else
                numeric_param=[numeric_param;str2num(dat)];
            end
        end
    end
%--------------------------------------------------------------------------
% Code below extracts all values of PolyApVec from MACOS prescription
%--------------------------------------------------------------------------
    if strcmp(string_param,'PolyApVec');
        numeric_param=[numeric_param 0 0];
        for jj=1:6;
            dat=fgetl(fid);
            numeric_param=[numeric_param;str2num(dat)];
        end
    end  
%--------------------------------------------------------------------------
% Code below extracts all values of TElt from MACOS prescription
%--------------------------------------------------------------------------
    if strcmp(string_param,'TElt');
        for jj=1:5;
            dat=fgetl(fid);
            %disp('Extracting TElt Values');
            numeric_param=[numeric_param; str2num(dat)];
        end
    end    
%--------------------------------------------------------------------------
% Code below extracts all values of TElt from MACOS prescription
%--------------------------------------------------------------------------
    if strcmp(string_param,'SegCoord');
        for jj=1:OpticalData{1}.nSeg-1;
            dat=fgetl(fid);
            %disp('Extracting TElt Values');
            numeric_param=[numeric_param; str2num(dat)];
        end
    end 
%--------------------------------------------------------------------------
% Code below extracts all values of Tout from MACOS prescription
%--------------------------------------------------------------------------
    if strcmp(string_param,'Tout');
        for jj=1:4;
            dat=fgetl(fid);
            %disp('Extracting TElt Values');
            numeric_param=[numeric_param; str2num(dat)];
        end
    end  
%--------------------------------------------------------------------------
% Code below extracts all values of ObsVec from MACOS prescription
%--------------------------------------------------------------------------
if exist('obs_count');
    if obs_count>1
    if strcmp(string_param,'ObsType');
        obs2_count=obs2_count+1;
        string_param=sprintf('ObsType%s',num2str(obs2_count));
    elseif strcmp(string_param,'ObsVec');
        obs3_count=obs3_count+1;
        string_param=sprintf('ObsVec%s',num2str(obs3_count));
    end
    end
end
%--------------------------------------------------------------------------
% Code below creates the structural array with all elements
%--------------------------------------------------------------------------
    if isempty(numeric_param); %this is the case when param value is a string and not a number
        aa=isspace(remain);
        ii=find(aa==0);
        remain=remain(ii);
        eval(['OpticalData{eltcount,1}.' string_param '=' 'remain' ';' ]);
    else %this is the case when param value is a number
        eval(['OpticalData{eltcount,1}.' string_param '=' 'numeric_param' ';' ]);
    end
end
fclose(fid);
numElt=eltcount;
disp(['Number of Optical Elements in ' macos_fname ' is ' num2str(numElt-2) '.']);

%==========================================================================
% This code generates a list of the names of the optics
%==========================================================================
List_EltNames=cell(length(OpticalData)-2,1);
for jj=1:length(OpticalData)-2;
    List_EltNames{jj}=OpticalData{jj+1}.EltName;
end

return
