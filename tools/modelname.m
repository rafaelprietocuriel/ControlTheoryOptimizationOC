function name=modelname(ocStruct)

name='';
if isempty(ocStruct)
    return
end
name=ocStruct.modelname;