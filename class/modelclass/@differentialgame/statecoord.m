function out=statecoord(dgObj)

out=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'statenum');
out=1:info.value;