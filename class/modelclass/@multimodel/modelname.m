function out=modelname(mmObj)

out='';
if isempty(mmObj)
    return
end

out=mmObj.modelname;