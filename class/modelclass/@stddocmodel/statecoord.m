function out=statecoord(docObj)

out=[];
if isempty(docObj)
    return
end
info=retrievemodelinformation(docObj.Model,'statenum');
out=1:info.value;