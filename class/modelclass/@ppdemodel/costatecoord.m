function out=costatecoord(ppdeObj)

out=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'statenum');
out=(1:info.value)+info.value;