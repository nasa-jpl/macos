%--------------------------------------------------------------------------
% Previous model prescription containing latest segment information,
% obscurations, and booms. This can help us define or insert elements which
% are not part of the new prescription such as booms, obscurations,
% segments, etc. (This is an optional part of the code which is up to the
% user but it is important if the elements such as booms and obscurations 
% have not been updated since the previous build)
%--------------------------------------------------------------------------
if strcmpi(OldData_Flag,'On');
    [Old_OpticalData,Old_num_Elt]=ReadRxFile(OldRx_fname);
    Old_EltNames=cell(length(Old_OpticalData)-2,1);

    for jj=1:length(Old_OpticalData)-2;
        Old_EltNames{jj}=Old_OpticalData{jj+1}.EltName;
    end

    %sort out specific elements from Old_EltNames
    seg_posi=strcmp(Old_EltNames,'Seg_A1');
    seg_posf=strcmp(Old_EltNames,'Seg_C6');
    segposi=find(seg_posi==1);
    segposf=find(seg_posf==1);

    for num=1:Old_num_Elt;
        if num==1
            disp('Source is ignored in constructing old build optical nodes');
        elseif num==Old_num_Elt; %last element in structure (nOutCord)
            disp('nOutCord is ignored in constructing old build optical nodes');
        else
            vstr = genvarname(sprintf('%s',Old_EltNames{num-1}));
            eval([vstr '= OpticalNamedNode(Old_EltNames{num-1},Old_OpticalData{num});']);
        end
    end
end %OldData_Flag