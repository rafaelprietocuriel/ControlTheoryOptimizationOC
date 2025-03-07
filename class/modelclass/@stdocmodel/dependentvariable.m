function var=dependentvariable(ocObj,arcarg)
%
% ARCARGUMENT
if nargin==1
    arcarg=0;
end
if isnumeric(arcarg)
    arcarg=num2str(arcarg);
end
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'equationvariablename',arcarg);
var=info.value;