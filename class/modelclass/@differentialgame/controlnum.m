function out=controlnum(dgObj)

out=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'controlnum');
out=info.value;