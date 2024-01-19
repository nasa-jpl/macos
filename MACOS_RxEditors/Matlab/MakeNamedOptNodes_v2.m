% Purpose: File reads MACOS optical prescriptions, creates a cell array
%          of structures which contains all elements in the prescription 
%          with their respective parameters as fields, and corresponding 
%          field values. We use a class definition script, labeled 
%          OptNodeDlList, to input the cell array and create a "doubly
%          linked list".
%
% Record:
%       Date                 Author                      Description
% =================     =======================      ====================
%  June 25th, 2008         Luis F. Marchen               Original Code
%
% Define Variables and Constants:
% OpticalData --cell array of structure containing optical data              
% num_Elt     --number of elements in read Rx includes Source and Tout
% vstr        --variable string to name each node as it is built and linked
% macos_fname --name of the MACOS prescription to be read

clear all
clc
%==========================================================================
% We read the Rx file created by converting Code V sequence file to MACOS; 
% otherwise, we choose the basic file to create the desired Rx file. After 
% reading file we build a cell array of structures containg all optics and 
% parameters in the variable OpticalData.
%==========================================================================
%macos_fname='most2prc_allseg_v3.in';
macos_fname='macosfile.in';
[OpticalData,num_Elt,List_EltNames]=ReadRxFile(macos_fname);

%Create list of element names for macos_fname prescription
%New_EltNames=new_EltNames(OpticalData);
%==========================================================================
% Here we use OpticalData to create a doubly linked list of optical
% elements for current prescription. 
%==========================================================================
create_dbl_list %uses OpticalData and creates doubly linked list 

%==========================================================================
% We delete selected surfaces from read Rx file using the doubly linked
% list constructed in previous step. There is two options: 1)Only
% disconnect the element from the list, or 2)disconnect and then destroy
% the element (this is done with the simple command delete). Notice that
% the iElt parameter will not be updated at this stage since we will be
% adding other optics to the list, and it will need to be updated as well.
%==========================================================================
Eltdel_Flag='Off'; %turn on/off delete optics
Eltdel_list=strvcat('Optic3','Optic4'); %list of optic/s to delete
Eltdel_choice='Delete'; %this must be either "delte" or "disconnect" only
%Elt_Delete %script used to delete specified optics from doubly linked list

elt_delete %m file which takes Eltdel_flag, Eltdel_list, and Eltdel_choice
           %and deletes or disconnects Eltdel_list from doubly linked list
%==========================================================================
% Read previous build and create optical nodes for all elements. Insert 
% segmented primary mirror into the list at the position where the
% monolithic primary is located in the prescription. 
%==========================================================================
%--------------------------------------------------------------------------
% Previous model prescription containing latest segment information,
% obscurations, and booms. This can help us define or insert elements which
% are not part of the new prescription such as booms, obscurations,
% segments, etc. (This is an optional part of the code which is up to the
% user but it is important if the elements such as booms and obscurations 
% have not been updated since the previous build)
%--------------------------------------------------------------------------
OldData_Flag='Off'; %flag to turn on/off reading of old build prescription
OldRx_fname='most2frnt_dbl.in';
create_old_optnodes %uses OldData_Flag and OldRx_fname to create old optical
                    %nodes of the class OpticalNamedNodes

%==========================================================================
% Creating a grid surface using create_grid_srf. To include the resulting
% node in the doubly linked list the node must be of the same class as the
% nodes in the doubly linked list (meaning node must be created from nodes 
% created by same class function), it may be renamed to be consistent with
% lists naming convention, and then inserted into the list by command
% "insert.Next" or "insert.Before"
%==========================================================================
if 0;
grid_elt_dat=RM2; %Optical Node which is being converted to a grid surface
xMon=[1 0 0];
ngrid=99;
diameter=6.52e-1;
lMon=6.52e-2;

[GridOptic]=create_grid_srf(grid_elt_dat,xMon,ngrid,diameter,lMon);                     
 vstrgridopt = genvarname(sprintf('%s','Optic100'));                      
 eval([vstrgridopt '= OpticalNamedNode(vstrgridopt,GridOptic);']);                    
                       
%Note: if the element being converted to a grid surface is already part of
%the doubly linked list we need to disconnect, not delete the element, and
%then we can use create_grid_srf.m to created the grid surface. But we must
%convert the resultant optic to an optical node of class OpticalNamedNode
%just as I did above then we may insert the element back in the same
%location as before

insertAfter(Optic100,Optic14); %Optic100 is RM2 converted to grid surface 
                               %and it is inserted after RM1 in the list
end                              
%==========================================================================
% convert doubly linked list to cell array of structures after making all
% necessary updates
%==========================================================================
dblist2cell %command with no inputs to convert doubly linked list
            %to  cell array of structures
%--------------------------------------------------------------------------
% After converting doubly linked list back to a cell array of structures we
% use the WriteRxFile function to write prescription to write_fname
%--------------------------------------------------------------------------
write_fname='pecover14_v1zmx.in';
credits=strvcat('author: Luis Marchen','Description: This files is only a test to test WriteRxFile.m');
WriteRxFile(write_fname,ElementData,credits)


