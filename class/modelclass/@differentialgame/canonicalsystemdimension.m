function [out,descr]=canonicalsystemdimension(dgObj)

out=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'canonicalsystemequationnum');
out=info.value;
descr=info.description;