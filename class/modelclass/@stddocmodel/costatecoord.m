function out=costatecoord(docObj)

out=[];
if isempty(docObj)
    return
end
info=retrievemodelinformation(docObj.Model,'statenum');
out=(1:info.value)+info.value;