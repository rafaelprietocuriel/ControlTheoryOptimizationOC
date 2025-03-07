function b=isidenticalmodel(ppdeObj1,ppdeObj2)
%
% ISIDENTICALMODEL test if two models are identical
%
% ISIDENTICALMODEL(OCOBJ1,OCOBJ2) two models are identical if they are of
% the same class, denote the same problem (modelname) and the parameter
% values are the same

b=1;
if ~strcmp(class(ppdeObj1),class(ppdeObj2))
    b=0;
    return
end

if ~strcmp(modelname(ppdeObj1),modelname(ppdeObj2))
    b=0;
    return
end
[parval1 parvar1]=parametervalue(ppdeObj1);
[parval2 parvar2]=parametervalue(ppdeObj2);
if numel(parval1)~=numel(parval2)
    b=0;
    return
end
if ~all(strcmp(parvar2,parvar1))
    b=0;
    return
end