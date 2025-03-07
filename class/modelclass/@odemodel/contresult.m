function out=contresult(odeObj)

out=[];
if isempty(odeObj)
    return
end
resultStruct=result(odeObj);
if isempty(resultStruct) || ~isfield(resultStruct,'Continuation')
    if nargout==1
        out=[];
    end
    return
end
out=resultStruct.Continuation;
