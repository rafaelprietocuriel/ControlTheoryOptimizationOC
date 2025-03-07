function var=discountrate(ocObj)
%
% DISCOUNTRATE value of the discountrate

if isempty(ocObj)
    return
end
info=retrievemodelinformation(ocObj.Model,'discountrate');
var=info.value;