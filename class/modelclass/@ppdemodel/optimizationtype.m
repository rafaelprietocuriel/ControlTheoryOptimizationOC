function otype=optimizationtype(ocObj)

otype='';
if isempty(ocObj)
    return
end
otype=optimizationtype(ocObj.Model);
