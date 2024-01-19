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
macos_fname='most3frntdbl.in';
[OpticalData,num_Elt]=ReadRxFile(macos_fname);

%Create list of element names for macos_fname prescription
New_EltNames=cell(length(OpticalData)-2,1);
for jj=1:length(OpticalData)-2;
    New_EltNames{jj}=OpticalData{jj+1}.EltName;
end

%==========================================================================
% Here we use OpticalData to create a doubly linked list of optical
% elements for current prescription. This list allows us to:
%    -Construct a node and assign parameter values to it
%    -Insert a node after a specified node
%    -Insert a node before a specified node
%    -Remove a node from the list
%    -Display data property on the command line
%    -Remove a node from the list before it is destroyed
%==========================================================================
for num=1:num_Elt;
    if num==1
        vstr = genvarname(sprintf('%s','Optic0'));
        sourcelabel='Optic0';
        eval([vstr '= OpticalNamedNode(sourcelabel(1,:),OpticalData{num});']);
    elseif num==num_Elt; %last element in structure (nOutCord)
        noutcordlabel='Opticf';
        vstr = genvarname(sprintf('%s','Opticf'));
        eval([vstr '= OpticalNamedNode(noutcordlabel(1,:),OpticalData{num});']);
    else
        vstr = genvarname(sprintf('%s%1.0f','Optic ', (num-1))); 
        %eval([vstr '= OpticalNamedNode(New_EltNames{num-1},OpticalData{num});']);
        eval([vstr '= OpticalNamedNode(vstr,OpticalData{num});']);
    end
    
    if num==2; 
        Optic1.insertAfter(Optic0)
    elseif num>2
        vstr2=genvarname(sprintf('%s%1.0f','Optic ', (num-2)));
        insertAfter(eval(vstr),eval(vstr2));
    end
end
stop
%==========================================================================
% We delete selected surfaces from read Rx file using the doubly linked
% list constructed in previous step. There is two options: 1)Only
% disconnect the element from the list, or 2)disconnect and then destroy
% the element (this is done with the simple command delete). Notice that
% the iElt parameter will not be updated at this stage since we will be
% adding other optics to the list, and it will need to be updated as well.
%==========================================================================
Eltdel_Flag='On'; %turn on/off delete optics
Eltdel_list=strvcat('Optic3','Optic4'); %list of optic/s to delete
Eltdel_choice='Delete'; %this must be either "delte" or "disconnect" only
%Elt_Delete %script used to delete specified optics from doubly linked list

if strcmpi(Eltdel_Flag,'On'); 
[rdel cdel]=size(Eltdel_list);
    for ii=1:rdel;
        vstr_del=genvarname(sprintf('%s', Eltdel_list(ii,:)));
        if strcmpi(Eltdel_choice,'disconnect');
            disconnect(eval(vstr_del));
        elseif strcmpi(Eltdel_choice,'delete');
            delete(eval(vstr_del));
        end %Eltdel_choice
    end %Eltdel_list
end %delete optic on/off

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
OldData_Flag='On'; %flag to turn on/off reading of old build prescription
OldRx_fname='most2prc_allseg_v3.in';

if strcmpi(OldData_Flag,'On');
    [Old_OpticalData,Old_num_Elt]=ReadRxFile(OldRx_fname);
    Old_EltNames=cell(length(Old_OpticalData)-2,1);

    for jj=1:length(Old_OpticalData)-2;
        Old_EltNames{jj}=Old_OpticalData{jj+1}.EltName;
    end

    %sort out specific elements from Old_EltNames
    seg_posi=strcmp(Old_EltNames,'Seg_A1');
    seg_posf=strcmp(Old_EltNames,'Seg_C6');
    segposi=find(seg_posi==1);
    segposf=find(seg_posf==1);

    for num=1:Old_num_Elt;
        if num==1
            disp('Source is ignored in constructing old build optical nodes');
        elseif num==Old_num_Elt; %last element in structure (nOutCord)
            disp('nOutCord is ignored in constructing old build optical nodes');
        else
            vstr = genvarname(sprintf('%s',Old_EltNames{num-1}));
            eval([vstr '= OpticalNamedNode(Old_EltNames{num-1},Old_OpticalData{num});']);
        end
    end
end %OldData_Flag

%--------------------------------------------------------------------------
% Here after having defined nodes for segments, booms and obscurations we
% rename those nodes and insert them in the doubly linked list
%--------------------------------------------------------------------------
Seg_InsertFlag=1; %flag to turn on/off insertion of segments
Seg_Insert_List=strvcat('Seg_A4','Seg_A5','Seg_C1','Seg_C4');


%==========================================================================
% Creating a grid surface using create_grid_srf. To include the resulting
% node in the doubly linked list the node must be of the same class as the
% nodes in the doubly linked list (meaning node must be created from nodes 
% created by same class function), it may be renamed to be consistent with
% lists naming convention, and then inserted into the list by command
% "insert.Next" or "insert.Before"
%==========================================================================
if 0;
Element_dat=RM2; %Optical Node which is being conveted to a grid surface
xMon=[1 0 0];
ngrid=99;
diameter=6.52e-1;
lMon=6.52e-2;

[GridOptic]=create_grid_srf(Element_dat,xMon,ngrid,diameter,lMon);                     
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
% After all changes have been done to the list we convert that back into a
% cell array of structures, done only to make it easy to pass that through
% a function to write the updated MACOS prescription.
%==========================================================================
ElementData{1}=Optic0.Data; %This assigns parameter and fields to a cell 
count=0;
nextElt=Optic0.Name;
ElementData{1}=Optic0.Data;

while ~strcmp(nextElt,'Opticf');%~isempty(eval([nextElt '.Next.Data']));
    count=count+1;
    ElementData{count+1}=eval([nextElt '.Next.Data']);
    nextElt=eval([nextElt '.Next.Name']);
end

%I need to make sure that the the optical node name matches the optical
%node display name in order to ran the the above 9 lines

%==========================================================================
% After converting doubly linked list back to a cell array of structures we
% use the WriteRxFile function to write prescription to write_fname
%==========================================================================
write_fname='testingwriteRxFile.in';
credits=strvcat('author: Luis Marchen','Description: This files is only a test to test WriteRxFile.m');
WriteRxFile(write_fname,ElementData,credits)


