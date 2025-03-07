function ppdeObjNew=mergeresult(ppdeObjNew,ppdeObj)

if isempty(ppdeObj)
    ppdeObjNew=ppdeObj;
    return
end
if ~strcmp(modelname(ppdeObjNew),modelname(ppdeObj))
end
fn=fieldnames(ppdeObj.Result);
for ii=1:length(fn)
    if isfield(ppdeObjNew.Result,fn{ii})
        ppdeObjNew.Result.(fn{ii})=[ppdeObjNew.Result.(fn{ii}) ppdeObj.Result.(fn{ii})];
    else
        ppdeObjNew.Result.(fn{ii})=ppdeObj.Result.(fn{ii});
    end
end