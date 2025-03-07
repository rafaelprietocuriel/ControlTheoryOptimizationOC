function out=matcontresult(odeObj)

out=[];
if isempty(odeObj)
    return
end
resultStruct=result(odeObj);
if isempty(resultStruct) || ~isfield(resultStruct,'MatContContinuation')
    if nargout==1
        out=[];
    end
    return
end
out=resultStruct.MatContContinuation;
