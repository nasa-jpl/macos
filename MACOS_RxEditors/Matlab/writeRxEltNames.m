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
%
% fid         ---file identifier for eading lines              
% count       ---number of iteration of while loop not counting empty lines
%                or comments
%
%

function []=writeRxEltNames(write_fname,ElementData,wElt)
%==========================================================================
% Comments or credits to be written to top of write file, and parameters 
%==========================================================================
overwrite_wfname=1;
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
% Write parameters for source and optical elements
%==========================================================================
[tt uu]=size(wElt);
for jj=1:tt;
    elt=wElt(jj,1)+1
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
            case 'EltName'
                temp_param=OpticalData{elt}.EltName;
                fprintf(fname_id, '%s\n',temp_param);
        end
    end
end;
fclose(fname_id);
return

