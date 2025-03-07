function arcid=arcidentifier(ocObj)
%
% ARCIDENTIFIER arcidentifiers are characters of numbers, consecutively
% numerated.

arcid='';
if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'arcidentifier');
arcid=info.value;