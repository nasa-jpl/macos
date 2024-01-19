function [TEltData]= ReadxyzLocals(xyzLocal_fname) 
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
fid2=fopen(xyzLocal_fname);
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
fid1=fopen(xyzLocal_fname);

%--------------------------------------------------------------------------
% Reading Parameter Names and Values, and Creating Structural
% Array Containing Optical Element Information Including Source and
% nOutCord (Tout and nOutCord are one element in structural array).
%--------------------------------------------------------------------------
count=0;
TElt=zeros(6,6,iEltcount);
while ~feof(fid1)
    dat= fgetl(fid1);
    [str, remain] = strtok(dat);
    numeric_param=str2num(remain);%converts remain to number, if not a number then empty
    n=strfind(str,'='); %find index of "=" sign to eliminate from str
    string_param=str(1:n-1); %parameter name without "=" sign, otherwise empty if part
                             %of previous parameter
                             
%--------------------------------------------------------------------------
% Element count
%--------------------------------------------------------------------------
    if strcmp(str,'iElt');
        count=count+1;
    end                             
%--------------------------------------------------------------------------
% Code below extracts all values of xyzLocals from xyzLocals file
%--------------------------------------------------------------------------
    if strcmp(string_param,'xLocal') ||strcmp(string_param,'yLocal') || strcmp(string_param,'zLocal');
        for jj=1:2;
            dat=fgetl(fid1);
            [strxyzL, remainxyzL] = strtok(dat);
            numeric_param=[numeric_param;str2num(remainxyzL)];
        end
        TElt(1:3,1:3,count)=numeric_param';
        TElt(4:6,4:6,count)=numeric_param';
    end
    
end
fclose(fid1);

TEltData=TElt;
    
return