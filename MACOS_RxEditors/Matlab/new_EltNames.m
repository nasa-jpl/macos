%==========================================================================
% This code generates a list of the names of the optics
%==========================================================================
function [New_EltNames]=new_EltNames(OpticalData)
%Create list of element names for macos_fname prescription
New_EltNames=cell(length(OpticalData)-2,1);
for jj=1:length(OpticalData)-2;
    New_EltNames{jj}=OpticalData{jj+1}.EltName;
end
return