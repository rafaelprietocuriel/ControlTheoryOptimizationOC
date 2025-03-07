function out=statenum(dgObj)

out=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'statenum');
out=info.value;