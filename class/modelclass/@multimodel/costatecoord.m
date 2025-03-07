function out=costatecoord(mmObj)
nummod=numberofmodels(mmObj);
out=cell(1,nummod);
if isempty(mmObj)
    return
end
for ii=1:nummod
    out{ii}=costatecoord(mmObj.Model{ii});
end