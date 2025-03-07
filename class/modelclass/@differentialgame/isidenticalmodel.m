function b=isidenticalmodel(dgObj1,dgObj2)
%
% ISIDENTICALMODEL test if two models are identical
%
% ISIDENTICALMODEL(OCOBJ1,OCOBJ2) two models are identical if they are of
% the same class, denote the same problem (modelname) and the parameter
% values are the same

b=1;
if ~strcmp(class(dgObj1),class(dgObj2))
    b=0;
    return
end

if ~strcmp(modelname(dgObj1),modelname(dgObj2))
    b=0;
    return
end
[parval1 parvar1]=parametervalue(dgObj1);
[parval2 parvar2]=parametervalue(dgObj2);
if numel(parval1)~=numel(parval2)
    b=0;
    return
end
if ~all(strcmp(parvar2,parvar1))
    b=0;
    return
end