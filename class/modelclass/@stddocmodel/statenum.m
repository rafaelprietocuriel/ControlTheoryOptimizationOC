function out=statenum(docObj)

out=[];
if isempty(docObj)
    return
end
info=retrievediffmodelinformation(docObj.Model,'statenum');
out=info.value;