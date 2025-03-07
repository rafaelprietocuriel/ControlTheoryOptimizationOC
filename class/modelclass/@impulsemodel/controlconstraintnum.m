function num=controlconstraintnum(ocObj)

num=[];
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'inequalitycontrolconstraintnum');
num=info.value;