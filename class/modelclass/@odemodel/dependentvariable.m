function var=dependentvariable(odeObj,arcarg)
%
% ARCARGUMENT
if nargin==1
    arcarg=0;
end
if isnumeric(arcarg)
    arcarg=num2str(arcarg);
end
if isempty(odeObj)
    return
end
info=retrieveodemodelinformation(odeObj.Model,'equationvariablename',arcarg);
var=info.value;