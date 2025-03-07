function var=dependentvariablenum(ocObj,arcid)
%
% ARCARGUMENT

if isnumeric(arcid)
    arcid=num2str(arcid);
end
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'canonicalsystemequationnum',arcid);
var=info.value;