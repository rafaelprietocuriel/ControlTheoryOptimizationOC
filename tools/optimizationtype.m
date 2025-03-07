function otype=optimizationtype(ocStruct)

otype='';
if isempty(ocStruct)
    return
end
otype=ocStruct.optimizationtype;