function name=submodelname(mmObj)

name='';
if isempty(mmObj)
    return
end

for ii=1:numberofmodels(mmObj)
    name{ii}=modelname(mmObj.Model{ii});
end