function out=constraintcombination(ocObj,arcarg)

out=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'constraintcombination',num2str(arcarg));
out=info.value;