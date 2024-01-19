%==========================================================================
% Here we use OpticalData to create a doubly linked list of optical
% elements for current prescription. This list allows us to:
%    -Construct a node and assign parameter values to it
%    -Insert a node after a specified node
%    -Insert a node before a specified node
%    -Remove a node from the list
%    -Display data property on the command line
%    -Remove a node from the list before it is destroyed
%==========================================================================
for num=1:num_Elt;
    if num==1
        vstr = genvarname(sprintf('%s','Optic0'));
        sourcelabel='Optic0';
        eval([vstr '= OpticalNamedNode(sourcelabel(1,:),OpticalData{num});']);
    elseif num==num_Elt; %last element in structure (nOutCord)
        noutcordlabel='Opticf';
        vstr = genvarname(sprintf('%s','Opticf'));
        eval([vstr '= OpticalNamedNode(noutcordlabel(1,:),OpticalData{num});']);
    else
        vstr = genvarname(sprintf('%s%1.0f','Optic ', (num-1))); 
        %eval([vstr '= OpticalNamedNode(New_EltNames{num-1},OpticalData{num});']);
        eval([vstr '= OpticalNamedNode(vstr,OpticalData{num});']);
    end
    
    if num==2; 
        Optic1.insertAfter(Optic0)
    elseif num>2
        vstr2=genvarname(sprintf('%s%1.0f','Optic ', (num-2)));
        insertAfter(eval(vstr),eval(vstr2));
    end
end