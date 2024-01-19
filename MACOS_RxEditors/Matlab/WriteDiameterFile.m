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
% Elements    ---structural array which contains parameters and values
% param_names ---lists the names of all parameters in MACOS prescription
% read_fname ---name of the MACOS prescription to be read
%
%==========================================================================
% Reading source parameters and values, and making the structural array
% section for source
%==========================================================================
function [] = WriteDiameterFile(read_fname,write_fname,param)
[rr cc]=size(param.rbSrf)
elt=param.rbSrf(2:rr,1); %the first row is 2 if param.rbSrf contains source
%--------------------------------------------------------------------------
% Opening Prescription File to be Read
%--------------------------------------------------------------------------
fid=fopen(read_fname);
%--------------------------------------------------------------------------
% Reading Parameter Names and Values, and Creating Structural
% Array Containing Optical Element Information Including Source and
% nOutCord (Tout and nOutCord are one element in structural array).
%--------------------------------------------------------------------------
eltcount=0;
while ~feof(fid)
    dat= fgetl(fid);
    n=strfind(dat,'='); %find index of "="sign to eliminate from str
    [nr nc]=size(dat);
    str=dat(n(2)+1:nc);    
    numeric_param=str2num(str);%converts remain to number, if not a number then empty
    string_param=dat(1:n(2));
    str=dat(1:n(1)-1);
    str=deblank(str); %removes trailing blank spaces
    a=isspace(str);
    b=find(a~=0);
    [nr1 nc1]=size(str);
    str=str(b+1:nc1)    
%--------------------------------------------------------------------------
% Determines and records when we begin reading a new optical element, 
% including the source, and the nOutCord, so that the parameter of each 
% elements can be grouped together in its corresponding structure
%--------------------------------------------------------------------------
    if strcmp(str,'iElt')==1;
        eltcount=eltcount+1
    end 

%--------------------------------------------------------------------------
% Code below creates the structural array with all elements
%--------------------------------------------------------------------------
    if isempty(numeric_param); %this is the case when param value is a string and not a number
        aa=isspace(remain);
        ii=find(aa==0);
        remain=remain(ii);
        eval(['diameter{eltcount,1}.' str 'eltcount' '=' 'remain' ';' ]);
    else %this is the case when param value is a number
        eval(['diameter{eltcount,1}' '=' 'numeric_param' ';' ]);
    end
end
numElt=eltcount;
fclose(fid);
disp(['Number of Optical Elements in ' read_fname ' is ' num2str(numElt) '.']);

%==========================================================================
% In this section we write the beam diameter to a file only the numeric
% data.
%==========================================================================
overwrite_wfname=1;
if exist(write_fname) == 2;	%fname already exists as file
    if overwrite_wfname;
        warning([write_fname ' already exists, and it is being deleted'])
        disp(['Deleting ' write_fname])
        delete(write_fname);
    else
        warning([write_fname ' already exists!!!'])

    end
end
fname_id = fopen(write_fname,'w');

%==========================================================================
% Write diameter for element to a text file to utilized by the Generic
% Coronagraph Errorbudget Tool
%==========================================================================
[uu kk]=size(elt)
for jj=1:uu;
    if strcmp(param.units,'m')==1;
        temp_param=diameter{elt(jj)}
    elseif strcmp(param.units,'mm')==1;
        temp_param=diameter{elt(jj)}./1000;
    end
    fprintf(fname_id, '%s',num2str(temp_param,'%1.15e'));
    if 1<jj<uu-1
        fprintf(fname_id, '\n');
    end
end
return




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    