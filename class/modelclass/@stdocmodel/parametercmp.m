function b=parametercmp(ocObj1,ocObj2)
% compares parameter values for those parameter variables that appear in
% both models
b=false;
if isempty(ocObj1) || isempty(ocObj2) 
    return
end

[parval1,parvar1]=parametervalue(ocObj1);
[parval2,parvar2]=parametervalue(ocObj2);

[dum idx1 idx2]=intersect(parvar1,parvar2);

b=all(parval1(idx1)==parval2(idx2));