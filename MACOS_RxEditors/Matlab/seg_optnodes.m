%--------------------------------------------------------------------------
% Previous model prescription containing latest segment information,
% obscurations, and booms. This can help us define or insert elements which
% are not part of the new prescription such as booms, obscurations,
% segments, etc. (This is an optional part of the code which is up to the
% user but it is important if the elements such as booms and obscurations 
% have not been updated since the previous build)
%--------------------------------------------------------------------------
    for jj=1:length(SegOpticData);
        Seg_EltNames{jj}=SegOpticData{jj}.EltName;
    end
    seg_numElt=length(Seg_EltNames);
  
    for num=1:seg_numElt;
        vstr_f = genvarname(['Seg' num2str(num) 'Forward']); %nodes for foward path
        vstr_r = genvarname(['Seg' num2str(num) 'Return']); %nodes for return path
        eval([vstr_f '= OpticalNamedNode('  '''Seg' num2str(num) 'Forward' ''' ,SegOpticData{num});']);
        %eval([vstr_r '= OpticalNamedNode(SegReturn' num2str(num) ',SegOpticData{num});']);
        eval([vstr_r '= OpticalNamedNode('  '''Seg' num2str(num) 'Return' ''' ,SegOpticData{num});']);
    end
