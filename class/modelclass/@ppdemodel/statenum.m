function out=statenum(ppdeObj)

out=[];
if isempty(ppdeObj)
    return
end
info=retrieveppdemodelinformation(ppdeObj.Model,'statenum');
out=info.value;