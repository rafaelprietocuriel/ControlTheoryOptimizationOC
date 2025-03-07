function subtype=dgsubtype(ocStruct)

subtype='';
if isempty(ocStruct)
    return
end
subtype=retrievedifferentialgameinformation(ocStruct,'subtype');
subtype=subtype.value;