function var=discountrate(ocObj)
%
% DISCOUNTRATE value of the discountrate

if isempty(ocObj)
    return
end
var=subsparametervalue(ocObj,ocObj.Model.objective.integral.discountrate);
%info=retrievemodelinformation(ocObj.Model,'discountrate');
%var=info.value;