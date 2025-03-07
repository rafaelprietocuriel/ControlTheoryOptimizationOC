function out=statecoord(odeObj)

out=[];
if isempty(odeObj)
    return
end
info=retrieveodemodelinformation(odeObj.Model,'statenum');
out=1:info.value;