function [ElementDataMod]=MakedblpassRx(singlepass_macos_fname,Ortho_Flag,write_txtfname,write_dblpass_fname)
%==========================================================================
% This scripts takes a single pass prescription and creates a double pass
% prescription retaining all elements and parameters from single pass. If
% there are lenses in the single pass prescription we must make sure we
% look at the double pass prescription and change the index of refraction.
% Another important thing is that we must also redefine the source for
% double pass Rx before it can work properly, this is done in command line
% MACOS
%==========================================================================
%==========================================================================
% We read the single pass Rx file and create a cell array of structures
% containing all elements and parameters exactly as they appear in Rx file
%==========================================================================
[OpticalData,num_Elt,List_EltNames]=ReadRxFile(singlepass_macos_fname);

%==========================================================================
% Creating cell array to input structures for double pass prescription
% created from the single pass prescription macos_fname
%==========================================================================
[row col]=size(OpticalData);
dblpasize=row+row-2;

dblpaData=cell(dblpasize,1);

f=flipdim(1:row-1,2); %determines the size of the forward elements
fw=[1 f(1:row-2)];    %vector containing forward elements iterations
ret=[2:row];          %vector containing iteration for return part of Rx
dblpavec=[fw,ret];    %vector combines iterations for forward and return 

%==========================================================================
% Creating double pass prescription. The source will need to be redefined
% and the index of refraction for lenses must be reversed from surface 2 of
% lens to surface 1 of lens for all lenses.
%==========================================================================

for jj=1:dblpasize;
    dblpaData{jj}=OpticalData{dblpavec(jj)};
end

%--------------------------------------------------------------------------
%Here we check prescription to see if elements contain TElt, and to
%determine if components are orthogonal, and to make nECoord positive if
%TElt exists and if it is negative
%--------------------------------------------------------------------------
if strcmpi(Ortho_Flag,'On');
ElementData=dblpaData;
ElementDataMod=TEltOrthoTst(ElementData,write_txtfname); %command which takes ElementData from dblist2cell, checks TElt
             %for orthogonality, and makes nECoord positive if negative
end
%--------------------------------------------------------------------------
% After converting doubly linked list back to a cell array of structures we
% use the WriteRxFile function to write prescription to write_fname
%--------------------------------------------------------------------------
if 1
credits=strvcat('author: Luis Marchen','Description: This files is a double pass prescription');
WriteRxFile(write_dblpass_fname,ElementDataMod,credits)                               
end                            
return
