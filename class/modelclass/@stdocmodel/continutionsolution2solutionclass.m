function out=continutionsolution2solutionclass(ocObj,contidx,solutionidx,varargin)

out=[];
if isempty(ocObj)
    return
end
resultStruct=result(ocObj);
if isempty(resultStruct) || ~isfield(resultStruct,'Continuation')
    if nargout==1
        out=[];
    end
    return
end
resultStruct.Continuation;
