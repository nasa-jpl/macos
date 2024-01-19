%==========================================================================
% We delete selected surfaces from read Rx file using the doubly linked
% list constructed in previous step. There is two options: 1)Only
% disconnect the element from the list, or 2)disconnect and then destroy
% the element (this is done with the simple command delete). The iElt
% parameters is updated at this stage and a new doubly linked list is
% created while the old one is deleted from command window.
%==========================================================================
if strcmpi(Eltdel_Flag,'On'); 
[rdel cdel]=size(Eltdel_list);
    for ii=1:rdel;
        vstr_del=genvarname(sprintf('%s', Eltdel_list(ii,:)));
        if strcmpi(Eltdel_choice,'disconnect');
            disconnect(eval(vstr_del));
        elseif strcmpi(Eltdel_choice,'Delete');
            delete(eval(vstr_del));
        end %Eltdel_choice
    end %Eltdel_list
    
return
%==========================================================================
% convert doubly linked list to cell array of structures after making all
% necessary updates. This generates a cell array of structures labeled
% ElementData
%==========================================================================
dblist2cell %command with no inputs to convert doubly linked list
            %to  cell array of structures
%--------------------------------------------------------------------------
% After converting doubly linked list back to a cell array of structures we
% use the WriteRxFile function to write prescription to write_fname
%--------------------------------------------------------------------------
write_fname='temp_Rx.in';
credits=strvcat('author: Luis Marchen','Description: File Updated by Deleting optics');
WriteRxFile(write_fname,ElementData,credits);

%--------------------------------------------------------------------------
% We clear all optical nodes, linked lists, and parameters from command
% window
%--------------------------------------------------------------------------
clear write_fname credits
clc

%==========================================================================
% Reload the temporay macos file to create optical nodes and to update
% parameter iElt for all optics
%==========================================================================
macos_fname='temp_Rx.in';
[OpticalData,num_Elt,List_EltNames]=ReadRxFile(macos_fname);

for ii=2:num_Elt-1;
    OpticalData{ii}.iElt=ii-1;
end
%Create list of element names for macos_fname prescription
%New_EltNames=new_EltNames(OpticalData);
%==========================================================================
% Here we use OpticalData to create a doubly linked list of optical
% elements for current prescription. 
%==========================================================================
create_dbl_list %uses OpticalData and creates doubly linked list 
end %delete optic on/off


