function idx=parameterindex(mmObj,parvar)
nummod=numberofmodels(mmObj);
idx=cell(1,nummod);

if isempty(mmObj)
    return
end
for ii=1:nummod
    if iscell(parvar)
        idx{ii}=parameterindex(mmObj.Model{ii},parvar{ii});
    else
        idx{ii}=parameterindex(mmObj.Model{ii},parvar);
    end
end