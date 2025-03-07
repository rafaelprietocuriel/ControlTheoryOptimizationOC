function var=dependentvariablenum(dgObj,arcid)
%
% ARCARGUMENT

if isnumeric(arcid)
    arcid=num2str(arcid);
end
if isempty(dgObj)
    return
end
info=retrievemodelinformation(dgObj.Model,'canonicalsystemequationnum',arcid);
var=info.value;