%==========================================================================
% This script takes ElementData which is created by dblist2cell, checks for
% TElt orthogonality of components for both 3x3 blocks, and also check
% nECoord sign, if current element in iteration has a TElt the nECoord
% parameter must be positive, if it is negative it is made positive
% otherwise it remains the same. (script created November 5th, 2008 by Luis 
% Marchen)
%==========================================================================
%macos_fname='most2prc_allseg_v3.in';
function[ElementDataMod]=TEltOrthoTst(ElementData,write_fname)
[row col]=size(ElementData); 
overwrite_wfname=1;
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
% Write file name and credits (comments) to top of write file
%==========================================================================
fprintf(fname_id,'%% %s\n',write_fname);
if exist('credits') == 1
  % allow multiline notes
  [m,n]= size(note);
  for i=1:m
    fprintf(fname_id,'%% %s\n', credits(i,:));    
  end
end
%date file
fprintf(fname_id,'%% Date Created: %s\n\n',datestr(date, 2))

%==========================================================================
% Check for TElt orthogonality
%==========================================================================
for ii=2:row-1;
    if any(strcmpi('TElt',fieldnames(ElementData{ii})));
        TElt=ElementData{ii}.TElt;
        mt1=TElt(1:3,1:3);
        mt2=TElt(1:3,1:3);
        test1=det(mt1'*mt1);
        test2=det(mt2'*mt2);
        residual(ii)=1-test1;
        if abs(1-test1)<1e-14 && abs(1-test2)<1e-14
            %disp(['TElt Exists for Element ' num2str(ElementData{ii}.iElt) '(' ElementData{ii}.EltName ')' 'and the local axes are orthogonal with 1-det(TEltT*TElt)= ' num2str(residual(ii))]);
            %fprintf(fname_id,'%% Date Created: %s\n\n',datestr(date, 2))
            fprintf(fname_id, '  %s\n',['TElt exists for element ' num2str(ElementData{ii}.iElt) '(' ElementData{ii}.EltName ')' ', and the local axes are orthogonal with 1-det(TEltT*TElt)= ' num2str(residual(ii))]);
        elseif abs(1-test1)>1e-13 || abs(1-test2)>1e-13
            warning(['Local axes are not orthogonal at element ' num2str(ElementData{ii}.iElt) '(' ElementData{ii}.EltName ')' 'with 1-det(TEltT*TElt)= ' num2str(residual(ii))]);
        end
%--------------------------------------------------------------------------
% Here if the element contains a TElt matrix for local perturbations we
% check nECoord, if it is negative it will be made positive otherwise it
% will remain positive
%--------------------------------------------------------------------------  
       if ElementData{ii}.nECoord==-6;
          ElementData{ii}.nECoord=ElementData{ii}.nECoord*-1;
       end
           
    else
        %disp(['Element' num2str(ElementData{ii}.iElt) '(' ElementData{ii}.EltName ')' 'does not contain TElt']);
        fprintf(fname_id, '  %s\n',['TElt does not exists for element ' num2str(ElementData{ii}.iElt) '(' ElementData{ii}.EltName ')' ]);
    end
end
fclose(fname_id);
ElementDataMod=ElementData;
return