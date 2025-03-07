function b=isidenticalmodel(docObj1,docObj2)
%
% ISIDENTICALMODEL test if two models are identical
%
% ISIDENTICALMODEL(OCOBJ1,OCOBJ2) two models are identical if they are of
% the same class, denote the same problem (modelname) and the parameter
% values are the same

b=1;
if ~strcmp(class(docObj1),class(docObj2))
    b=0;
    return
end

if ~strcmp(modelname(docObj1),modelname(docObj2))
    b=0;
    return
end
[parval1 parvar1]=parametervalue(docObj1);
[parval2 parvar2]=parametervalue(docObj2);
if numel(parval1)~=numel(parval2)
    b=0;
    return
end
if ~all(strcmp(parvar2,parvar1))
    b=0;
    return
end