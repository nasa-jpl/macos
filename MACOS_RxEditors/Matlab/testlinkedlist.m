clear all
clc

%==========================================================================
% We read the Rx file created by converting Code V sequence file to MACOS; 
% otherwise, we choose the basic file to create the desired Rx file. After 
% reading file we build a cell array of structures containg all optics and 
% parameters in the variable OpticalData.
%==========================================================================
macos_fname='most2prc_allseg_v3.in';
[OpticalData,num_Elt]=ReadRxFile(macos_fname);

%==========================================================================
% Here we use OpticalData to create a doubly linked list of optical
% elements. This list allows us to:
%    -Construct a node and assign parameter values to it
%    -Insert a node after a specified node
%    -Insert a node before a specified node
%    -Remove a node from the list
%    -Display data property on the command line
%    -Remove a node from the list before it is destroyed
%==========================================================================
for num=1:num_Elt;
    if num==1
        vstr = genvarname(sprintf('%s','Source '));
    elseif num==num_Elt; %last element in structure (nOutCord)
        vstr = genvarname(sprintf('%s','nOutCord'));
    else
        vstr = genvarname(sprintf('%s%1.0f','Optic ', (num-1)));
    end
eval([vstr '= OptNodeDlList(OpticalData{num});']);
if num==2; %link elements to make linked list
    Optic1.insertAfter(Source)
elseif num>2
    vstr2=genvarname(sprintf('%s%1.0f','Optic ', (num-2)));
    insertAfter(eval(vstr),eval(vstr2));
end
end

