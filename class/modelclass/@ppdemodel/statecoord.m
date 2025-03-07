function out=statecoord(ppdeObj)

out=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'statenum');
out=1:info.value;