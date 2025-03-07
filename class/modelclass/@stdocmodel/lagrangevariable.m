function var=lagrangevariable(ocObj)
%
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'lagrangemultipliercontrolname');
var=info.value;