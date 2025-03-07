function out=isautonomous(mmObj)
nummod=numberofmodels(mmObj);
out=cell(1,nummod);
if isempty(mmObj)
    return
end
for ii=1:nummod
    out{ii}=isautonomous(mmObj.Model{ii});
end