fid1=fopen(macos_fname);
Elt_Count=0;
while ~feof(fid1)
    dat= fgetl(fid1);
    [str, remain] = strtok(dat);
    numeric_param=str2num(remain);%converts remain to number, if not a number then empty
    n=strfind(str,'='); %find index of "="sign to eliminate from str
    string_param=str(1:n-1); %parameter name without "=" sign, otherwise empty if part
                             %of previous parameter
    if isempty(string_param)|| strncmp(string_param,'%',1);
        %disp('The Line Is Empty and No Structure Is Created');
        continue
    end;
    count=count+1;
    if  strcmp(string_param,'iElt')
        Elt_Count=Elt_Count+1;
    end
end
numElt=Elt_Count;