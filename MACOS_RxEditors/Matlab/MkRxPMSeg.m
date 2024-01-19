clear all
close all
clc
ifplot=1;

% macos_fname='most3frntdbl.in';
% [OpticalData,num_Elt,List_EltNames]=ReadRxFile(macos_fname);
% 
%--------------------------------------------------------------------------
% Set up parameters for segment generator
%--------------------------------------------------------------------------
Hd=1162; %vertext to vertex distance, units mm
He=1006; %side to opposing side distance, units mm
gap=65;  %size of gab between segments, units mm
Kr=7400; %radius of curvature, units mm
Kc=-1;   %conic constant
ifSurface='Conic'; %specifies the type of primary mirror surface Conic or Aspheric

%==========================================================================
%the following script computes vertices for segments (set up for 18
%segmetns only for AMD), when the script runs it generates
%Segment_Vertices.
%==========================================================================
[seg_data,SegCen,TElt]=PMSegGenerator(Hd/1000,He/1000,Kc,Kr/1000,gap/1000,ifSurface);
return
%--------------------------------------------------------------------------
%parameters from PMSegGenerator which are required by mono2segs.m
%--------------------------------------------------------------------------
seginfo.seg_Kc=Kc;
seginfo.seg_Kr=Kr.*1000; %converts meters to mm
seginfo.Hd=Hd.*1000;     %convert meters to mm
seginfo.seg_TElt=TElt;
seginfo.seg_ap_cen=SegCen.*1000; %converts meters to mm
seginfo.numseg=6;
seginfo.nGridMat=256;
seginfo.seg_data=seg_data.*1000; %converts meters to mm
seginfo.seg_name={'A1','A2','A3','A4', 'A5','A6','B1','C1','B2','C2','B3',...
                                       'C3','B4','C4','B5','C5','B6','C6'};
%--------------------------------------------------------------------------
% plot surface: conic, or aspheric
%--------------------------------------------------------------------------
%PLOT_AMD_PARABOLA%plots the PM parabola if ifplot=1, and creates fields 
                  %for parabola parameters which will be use in MACOS Rx
%--------------------------------------------------------------------------
% the user defined function mono2segs.m generate a structural array with
% each segment parameters given segment data in AMD structural array. the
% script seg_optnodes generates the nodes for the doubly linked list.
%--------------------------------------------------------------------------
SegOpticData=mono2segs(seginfo);
seg_optnodes
%==========================================================================
% We read the Rx file created by converting Code V sequence file to MACOS; 
% otherwise, we choose the basic file to create the desired Rx file. After 
% reading file we build a cell array of structures containg all optics and 
% parameters in the variable OpticalData.
%==========================================================================
%macos_fname='most2prc_allseg_v3.in';
macos_fname='Optiix1.in';
[OpticalData,num_Elt,List_EltNames]=ReadRxFile(macos_fname);

%--------------------------------------------------------------------------
% Here we use OpticalData to create a doubly linked list of optical
% elements for current prescription. 
%--------------------------------------------------------------------------
create_dbl_list %uses OpticalData and creates doubly linked list 

%==========================================================================
% Delete specified surfaces (Eltdel_list) from doubly linked list. There is
% two options: 1)Only disconnect the element from the list, or 
% 2)disconnect and then destroy the element (this is done with the simple 
%   command delete). Notice that the iElt parameter will not be updated at 
%   this stage since we will be adding other optics to the list, and it 
%   will need to be updated as well.
%==========================================================================
if 0
    Eltdel_Flag='On'; %turn on/off delete optics
    %Eltdel_list=strvcat('Optic22','Optic24'); %list of optic/s to delete
    Eltdel_list=strvcat('Optic2'); %list of optic/s to delete
    
    Eltdel_choice='Delete'; %this must be either "delte" or "disconnect" only
    %Elt_Delete %script used to delete specified optics from doubly linked list
    
    elt_delete %m file which takes Eltdel_flag, Eltdel_list, and Eltdel_choice
    %and deletes or disconnects Eltdel_list from doubly linked list
end
%==========================================================================
% Creating a grid surface using create_grid_srf. To include the resulting
% node in the doubly linked list the node must be of the same class as the
% nodes in the doubly linked list (meaning node must be created from nodes 
% created by same class function), it may be renamed to be consistent with
% lists naming convention, and then inserted into the list by command
% "insert.Next" or "insert.Before"
%==========================================================================
nextElt=Optic1.Name; %defines location of node where we will insert segments
                      %into the doubly linked list
%--------------------------------------------------------------------------
% inserts the forward path segments to MACOS prescription in the place of
% the monolithic PM defined by nextElt above
%--------------------------------------------------------------------------
for ii=1:seginfo.numseg
    inForElt=['Seg' num2str(ii) 'Forward']; 
    inNextElt=[nextElt '.Next'];

    eval(['insertAfter(' inForElt ',' inNextElt ')'])
    nextElt=eval([nextElt '.Next.Name']);
end
                              
%--------------------------------------------------------------------------
% inserts the return path segments to MACOS prescription in the place of
% the monolithic return PM defined by last nextElt from loop above 
%--------------------------------------------------------------------------
if 0
nextElt=eval([nextElt '.Next.Name']);
for ii=1:AMD.numseg   
    inForElt=['Seg' num2str(ii) 'Return']; 
    inNextElt=[nextElt '.Next'];

    eval(['insertAfter(' inForElt ',' inNextElt ')'])
    nextElt=eval([nextElt '.Next.Name']);
end
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
write_fname='Optiix1_seg_v1.in';
credits=strvcat('author: Luis Marchen','Description: File with segmented PM');
WriteRxFile(write_fname,ElementData,credits)





