function idx=constraintcombinationindex(ocObj,arcarg)

idx=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'constraintcombinationindex',num2str(arcarg));
idx=info.value;