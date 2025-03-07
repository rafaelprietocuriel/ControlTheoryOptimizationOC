function out=contresult(ppdeObj)

out=[];
if isempty(ppdeObj)
    return
end
resultStruct=result(ppdeObj);
if isempty(resultStruct) || ~isfield(resultStruct,'Continuation')
    if nargout==1
        out=[];
    end
    return
end
out=resultStruct.Continuation;
