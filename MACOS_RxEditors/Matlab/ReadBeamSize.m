function [beamDiameter,num_Elt]= ReadBeamSize(beamSize_fname) 
%==========================================================================
% Reads xyzLocals and creates:
% 
% iEltcount => number of elements
% count     => index of multidimensional array
% TElt      => multidimensional array with TElt's for each element
%==========================================================================

%==========================================================================
% Running a while loop to do an element count to create a multidimensional
% array for TElt
%==========================================================================
fid2=fopen(beamSize_fname);
iEltcount=0;
while ~feof(fid2)
    dat= fgetl(fid2);
    [str, remain] = strtok(dat);
%--------------------------------------------------------------------------
% Determines and records when we begin reading a new optical element
%--------------------------------------------------------------------------
    if strcmp(str,'iElt');
        iEltcount=iEltcount+1;
    end
end
fclose(fid2);

%==========================================================================
% Read text file with xyzLocal data and create the TElt matrices in a
% multidimensional array
%==========================================================================
fid1=fopen(beamSize_fname);

%--------------------------------------------------------------------------
% Reading beam diameter for each optic. the code identifies the equal sign
% when it occurs the second instance then uses strtok to separate number
% from string.
%--------------------------------------------------------------------------
count=0;
beamDiameter=zeros(iEltcount,1);
while ~feof(fid1)
    dat= fgetl(fid1);
%--------------------------------------------------------------------------
% Element count
%--------------------------------------------------------------------------
[str, remain] = strtok(dat);
    if strcmp(str,'iElt');
        count=count+1;
    end  
    
%--------------------------------------------------------------------------
% commence data extraction
%--------------------------------------------------------------------------
    k=strfind(dat,'='); %finds instance when equal sign occurs second time
    str1=dat(k(2):end);  %creates a string from the equal sign to end of string
    [str2, remain] = strtok(str1);
    numeric_param=str2num(remain);%converts remain to number, if not a number then empty
    beamDiameter(count,1)=numeric_param;%./1000;                                             

end
num_Elt=count+2;
fclose(fid1);

return







