function out=contresult(docObj)

out=[];
if isempty(docObj)
    return
end
resultStruct=result(docObj);
if isempty(resultStruct) || ~isfield(resultStruct,'Continuation')
    if nargout==1
        out=[];
    end
    return
end
out=resultStruct.Continuation;
