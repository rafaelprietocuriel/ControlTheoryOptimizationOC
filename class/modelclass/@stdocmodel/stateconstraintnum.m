function num=stateconstraintnum(ocObj)

num=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'inequalitystateconstraintnum');
num=info.value;