function var=dependentvariablenum(docObj,arcid)
%
% ARCARGUMENT

if isnumeric(arcid)
    arcid=num2str(arcid);
end
if isempty(docObj)
    return
end
info=retrievediffmodelinformation(docObj.Model,'canonicalsystemequationnum',arcid);
var=info.value;