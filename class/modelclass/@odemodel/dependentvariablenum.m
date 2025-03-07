function var=dependentvariablenum(odeObj,arcid)
%
% ARCARGUMENT

if isnumeric(arcid)
    arcid=num2str(arcid);
end
if isempty(odeObj)
    return
end
info=retrieveodemodelinformation(odeObj.Model,'statenum');
var=info.value;