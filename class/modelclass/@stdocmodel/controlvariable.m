function var=controlvariable(ocObj)
%
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'controlname');
var=info.value;