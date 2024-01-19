%This code is an example to show how to build a doubly linked list given
%the OpticalData data structure and the number of elements num_Elt.


function [Source,nOutCord,varargout]=ttest(OpticalData,num_Elt)

for num=1:num_Elt;
    if num==1
        vstr = genvarname(sprintf('%s','Source '));
    elseif num==num_Elt; %last element in structure (nOutCord)
        vstr = genvarname(sprintf('%s','nOutCord'));
    else
        vstr = genvarname(sprintf('%s%1.0f','Optic ', (num-1)));
    end
    eval([vstr '= OptNodeDlList(OpticalData{num});']);
    if num==2; %link elements to make linked list
        Optic1.insertAfter(Source)
    elseif num>2
        vstr2=genvarname(sprintf('%s%1.0f','Optic ', (num-2)));
        insertAfter(eval(vstr),eval(vstr2));
    end
end
return