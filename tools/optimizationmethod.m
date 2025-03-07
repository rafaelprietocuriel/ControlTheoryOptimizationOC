function method=optimizationmethod(ocStruct)

method='';
if isempty(ocStruct)
    return
end
method=ocStruct.optimization.method;