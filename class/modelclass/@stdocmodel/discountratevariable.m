function var=discountratevariable(ocObj)
%
% DISCOUNTRATE value of the discountrate

var=[];
if isempty(ocObj)
    return
end
var=ocObj.Model.objective.integral.discountrate;
%info=retrievemodelinformation(ocObj.Model,'discountrate');
%var=info.value;