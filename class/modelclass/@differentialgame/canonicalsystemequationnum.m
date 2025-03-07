function [out,descr]=canonicalsystemequationnum(dgObj)

out=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'canonicalsystemequationnum');
out=info.value;
descr=info.description;