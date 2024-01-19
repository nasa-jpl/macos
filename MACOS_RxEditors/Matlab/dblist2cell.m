%==========================================================================
% After all changes have been done to the list we convert that back into a
% cell array of structures, done only to make it easy to pass that to
% a function to write the updated MACOS prescription.
%==========================================================================
clear ElementData;
ElementData{1}=Optic0.Data; %This assigns parameter and fields to a cell 
count=0;
nextElt=Optic0.Name;
%ElementData{1}=Optic0.Data;

while ~strcmp(nextElt,'Opticf');%~isempty(eval([nextElt '.Next.Data']));
    count=count+1;
    ElementData{count+1}=eval([nextElt '.Next.Data']);
    nextElt=eval([nextElt '.Next.Name']);
end

%==========================================================================
% Update the iElt parameters, that is the element number in Rx
%==========================================================================
len=length(ElementData);

for jj=1:(len-2);
   ElementData{jj+1}.iElt= jj; 
end