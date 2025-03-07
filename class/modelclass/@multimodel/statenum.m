function out=statenum(mmObj)
nummod=numberofmodels(mmObj);
out=cell(1,nummod);
if isempty(mmObj)
    return
end
for ii=1:nummod
    out{ii}=statenum(mmObj.Model{ii});
end