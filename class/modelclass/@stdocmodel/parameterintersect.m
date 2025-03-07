function [parvar idx1 idx2]=parameterintersect(ocObj1,ocObj2)
% compares parameter values for those parameter variables that appear in
% both models
idx1=[];
idx2=[];

if isempty(ocObj1) || isempty(ocObj2) 
    return
end

[parval1,parvar1]=parametervalue(ocObj1);
[parval2,parvar2]=parametervalue(ocObj2);
