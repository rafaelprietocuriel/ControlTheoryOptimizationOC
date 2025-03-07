function name=initfilename(ocStruct)

name='';
if isempty(ocStruct)
    return
end

if isfield(ocStruct,'initfilename')
    name=[ocStruct.initfilename '.' ocStruct.basicextension 'm'];
else
    name=[ocStruct.modelname '.' ocStruct.basicextension 'm'];
end