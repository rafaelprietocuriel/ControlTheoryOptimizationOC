function mext=basicextension(ocStruct)

mext='';
if isempty(ocStruct)
    return
end
try
    mext=ocStruct.basicextension;
catch
    mext='oc';
end