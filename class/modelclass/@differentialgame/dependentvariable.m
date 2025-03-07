function var=dependentvariable(dgObj,arcarg)
%
% ARCARGUMENT
if nargin==1
    arcarg=0;
end
if isnumeric(arcarg)
    arcarg=num2str(arcarg);
end
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'equationvariablename',arcarg);
var=info.value;