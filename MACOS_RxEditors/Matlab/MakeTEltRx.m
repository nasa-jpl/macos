
function[OpticalData]=MakeTEltRx(macos_fname,write_fname,credits)
%==========================================================================
% We read the Rx file created by converting Code V sequence file to MACOS; 
% otherwise, we choose the basic file to create the desired Rx file. After 
% reading file we build a cell array of structures containg all optics and 
% parameters in the variable OpticalData.
%==========================================================================
[OpticalData,num_Elt,List_EltNames]=ReadRxFile(macos_fname); 

%==========================================================================
% Reads xyzLocal.txt and creates TEltData a multidimensional array of
% TElt's
%==========================================================================
xyzLocal_fname='xyzLocal.txt';
TEltData=ReadxyzLocals(xyzLocal_fname);

for mm=1:num_Elt-2;
    TEltData(1:3,3,mm)=OpticalData{mm+1}.psiElt';
    TEltData(4:6,6,mm)=OpticalData{mm+1}.psiElt';
end

%==========================================================================
%Create TElt for elements
%==========================================================================
for jj=1:num_Elt-2;
   OpticalData{jj+1}.TElt=TEltData(:,:,jj);
end

%==========================================================================
% Check for TElt orthogonality, and if TElt is not orthogonal for a
% particular element it is orthogonalized by get_localframeg.m.
%==========================================================================
for ii=2:num_Elt-1;
    if any(strcmpi('TElt',fieldnames(OpticalData{ii})));
        TElt=OpticalData{ii}.TElt;

        mt1=TElt(1:3,1:3);
        mt2=TElt(1:3,1:3);
        test1=det(mt1'*mt1);
        test2=det(mt2'*mt2);
        residual(ii)=1-test1;
        if abs(1-test1)<1e-14 && abs(1-test2)<1e-14
            disp(['TElt Exists for Element ' num2str(OpticalData{ii}.iElt) '(' OpticalData{ii}.EltName ')' 'and the local axes are orthogonal with 1-det(TEltT*TElt)= ' num2str(residual(ii))]);
        elseif abs(1-test1)>1e-13 || abs(1-test2)>1e-13
            warning(['Local axes are not orthogonal at element ' num2str(OpticalData{ii}.iElt) '(' OpticalData{ii}.EltName ')' 'with 1-det(TEltT*TElt)= ' num2str(residual(ii))]);
            disp(['Orthogonalizing local axes for element ' num2str(OpticalData{ii}.iElt) '(' OpticalData{ii}.EltName ')']);
%--------------------------------------------------------------
%orthogonalizing x and y-components of TElt to z-component
%--------------------------------------------------------------
            v=[OpticalData{ii}.TElt(1:3,3) OpticalData{ii}.TElt(1:3,1)];
            c = get_localframeg(v,[3 1]);
            c=c';
%-----------------------------------------------------------------------
%Create new TElt for current optic with new orthogonal set of components
%-----------------------------------------------------------------------
            OpticalData{ii}.TElt=zeros(6,6);
            OpticalData{ii}.TElt(1:3,1:3)=c;
            OpticalData{ii}.TElt(4:6,4:6)=c;
        end
%--------------------------------------------------------------------------
% Here if the element contains a TElt matrix for local perturbations we
% check nECoord, if it is negative it will be made positive otherwise it
% will remain positive
%--------------------------------------------------------------------------
        if OpticalData{ii}.nECoord==-6;
            OpticalData{ii}.nECoord=OpticalData{ii}.nECoord*-1;
        end
    else
        disp(['Element' num2str(OpticalData{ii}.iElt) '(' OpticalData{ii}.EltName ')' 'does not contain TElt']);
    end
end
%==========================================================================
% Write macos prescription after adding the TElt parameter to each element
% and changing the sign of the nECoord parameter
%==========================================================================
if 1
WriteRxFile(write_fname,OpticalData,credits)                               
end
return


    