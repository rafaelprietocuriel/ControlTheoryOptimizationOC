function b=isidenticalmodel(ocObj1,ocObj2)
%
% ISIDENTICALMODEL test if two models are identical
%
% ISIDENTICALMODEL(OCOBJ1,OCOBJ2) two models are identical if they are of
% the same class, denote the same problem (modelname) and the parameter
% values are the same

b=1;
if ~strcmp(class(ocObj1),class(ocObj2))
    b=0;
    return
end

if ~strcmp(modelname(ocObj1),modelname(ocObj2))
    b=0;
    return
end
[parval1 parvar1]=parametervalue(ocObj1);
[parval2 parvar2]=parametervalue(ocObj2);
if numel(parval1)~=numel(parval2)
    b=0;
    return
end
if ~all(strcmp(parvar2,parvar1))
    b=0;
    return
end